      !-------------------------------------------------------------------------
      !  Subroutine   :                    ppm_rc_init_seg
      !-------------------------------------------------------------------------
      ! Copyright (c) 2016 MOSAIC Group (MPI-CBG Dresden)
      !
      !
      ! This file is part of the PPM_RC program.
      !
      ! PPM_RC is free software: you can redistribute it and/or modify
      ! it under the terms of the GNU Lesser General Public License
      ! as published by the Free Software Foundation, either
      ! version 3 of the License, or (at your option) any later
      ! version.
      !
      ! PPM_RC is distributed in the hope that it will be useful,
      ! but WITHOUT ANY WARRANTY; without even the implied warranty of
      ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
      ! GNU General Public License for more details.
      !
      ! You should have received a copy of the GNU General Public License
      ! and the GNU Lesser General Public License along with PPM_RC. If not,
      ! see <http://www.gnu.org/licenses/>.
      !-------------------------------------------------------------------------
      !  MOSAIC Group
      !  Max Planck Institute of Molecular Cell Biology and Genetics
      !  Pfotenhauerstr. 108, 01307 Dresden, Germany
      !
      !  Author           - y.afshar           April   2013
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Subroutine   :                    ppm_rc_init_seg
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Initialize the RC client and read all input files.
      !
      !  Input        :
      !
      !  Input/output :
      !
      !  Output       : info          (I) return status. 0 on success.
      !
      !  Routines     : ppm_rc_read_input
      !                 ppm_rc_defaults
      !                 ppm_rc_read_ctrl
      !
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      SUBROUTINE DTYPE(ppm_rc_init_seg)(info)

        !-------------------------------------------------------------------------
        !  Modules
        !-------------------------------------------------------------------------
        USE ppm_module_data, ONLY : ppm_error_error,             &
        & ppm_param_decomp_xpencil,ppm_param_decomp_xy_slab,     &
        & ppm_param_assign_internal,ppm_rank,ppm_param_decomp_cartesian
        ! ppm_param_assign_metis_comm,ppm_param_decomp_bisection
        USE ppm_module_io, ONLY : ppm_io_set_unit
        USE ppm_module_mktopo, ONLY : ppm_mktopo
        USE ppm_module_topo_typedef, ONLY : ppm_t_topo,ppm_topo
        USE ppm_module_mesh_typedef, ONLY : ppm_t_subpatch_
        USE ppm_module_topo_get, ONLY : ppm_topo_get
        USE ppm_module_mpi
        USE ppm_module_mapping_typedef

        USE ppm_rc_module_write, ONLY : DTYPE(ppm_rc_write_image), &
        &   DTYPE(ppm_rc_write_image_label),DTYPE(ppm_rc_read_write_image), &
        &   DTYPE(ppm_rc_read_write_unique_label)
        USE ppm_rc_module_read, ONLY : ppm_rc_read_image_info, &
        &   DTYPE(ppm_rc_read_image)
        USE ppm_rc_module_util
        USE ppm_rc_module_filter
        USE ppm_rc_module_linkedlist
        USE ppm_rc_module_fire
        USE ppm_rc_module_topologicalnumber, ONLY : FG_ConnectivityType
        IMPLICIT NONE

        !-------------------------------------------------------------------------
        !  Includes
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        INTEGER, INTENT(  OUT) :: info

        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
        CLASS(ppm_t_subpatch_), POINTER :: sbpitr

        CLASS(ppm_t_discr_info_),  POINTER :: imgdinfo
        CLASS(ppm_t_discr_info_),  POINTER :: lbldinfo
        CLASS(ppm_t_discr_info_),  POINTER :: initdinfo
        !!! temporary pointers to the data and bookkeeping information
        !!! for a mesh or particle on which fieldID has been discretized.

        TYPE(ppm_rc_list), POINTER :: seed

        TYPE(ppm_t_topo), POINTER :: topo

        REAL(ppm_kind_double), DIMENSION(__DIME)               :: Offset
        REAL(ppm_kind_double)                             :: t0
        REAL(MK)                                          :: FieldMin,FieldMax
        REAL(MK)                                          :: sx,ex,lx
        REAL(MK)                                          :: sy,ey,ly
#if   __DIME == __2D
        REAL(MK), CONTIGUOUS,  DIMENSION(:,:),    POINTER :: wpr
#elif __DIME == __3D
        REAL(MK), CONTIGUOUS,  DIMENSION(:,:,:),  POINTER :: wpr
        REAL(MK)                                          :: sz,ez,lz
#endif
#ifdef __MPI
        REAL(MK),             DIMENSION(:,:), POINTER     :: min_sub
        REAL(MK),             DIMENSION(:,:), POINTER     :: max_sub
#endif
        REAL(MK), CONTIGUOUS, DIMENSION(:,:), POINTER     :: xp
#if   __DIME == __2D
        INTEGER, CONTIGUOUS,  DIMENSION(:,:),  POINTER :: wpil
        INTEGER, CONTIGUOUS,  DIMENSION(:,:),  POINTER :: wpip
#elif __DIME == __3D
        INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER :: wpil
        INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER :: wpip
#endif
        INTEGER,             DIMENSION(:),     POINTER :: seedn
        INTEGER,             DIMENSION(:),     POINTER :: Nm
        INTEGER,             DIMENSION(:),     POINTER :: hi_a
        INTEGER,             DIMENSION(:),     POINTER :: lo_a
        INTEGER,             DIMENSION(__DIME)              :: ldu,ld,dd
        INTEGER                                        :: nproc
        INTEGER                                        :: vLabel
        INTEGER                                        :: i,j,k,l,iopt
        INTEGER                                        :: ipatch,ipart
        INTEGER                                        :: nsize,iseed
#ifdef __Linux
        INTEGER                                        :: memory
#endif
#ifdef __MPI
        INTEGER,             DIMENSION(:), ALLOCATABLE :: ineigh,nneigh
        INTEGER,             DIMENSION(:), ALLOCATABLE :: displ,flag
        INTEGER                                        :: request
#endif

        CHARACTER(LEN=ppm_char) :: caller = 'ppm_rc_init_seg'

        LOGICAL :: IsEnclosedByLabel

        !-------------------------------------------------------------------------
        !  Externals
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
        !  Initialize
        !-------------------------------------------------------------------------
        CALL substart(caller,t0,info)

!         debug=2
        !-------------------------------------------------------------------------
        !  Find the right number of processors for reading the image file.
        !
        ! TODO
        ! Should fix this part to be able to work with arbitrary number
        ! of processors. Currently for the 3D stack, the number of readers
        ! (processors who actually read the image) is limited by the number of Z.
        !-------------------------------------------------------------------------

        !-------------------------------------------------------------------------
        ! First : save the total number of processors
        !-------------------------------------------------------------------------
        nproc=ppm_nproc

        SELECT CASE (ppm_nproc)
        CASE (1)
           !-------------------------------------------------------------------------
           !  define the ghostsize for IO topology
           !  On one processor we only have one topology, one mesh, and etc.
           !  So we are using one set of ghost for all of them
           !-------------------------------------------------------------------------
           IF (nsteql.GT.0) THEN
              ghostsize=MAX(ioghostsize,inighostsize,ghostsize_run,ghostsize_equil)
           ELSE
              ghostsize=MAX(ioghostsize,inighostsize,ghostsize_run)
           ENDIF
        CASE DEFAULT
           !-------------------------------------------------------------------------
           !  define the ghostsize for IO topology
           !-------------------------------------------------------------------------
           ghostsize=ioghostsize

           !-------------------------------------------------------------------------
           !  If the user provides the number of IO processors
           !  use it as a base.
           !-------------------------------------------------------------------------
           IF (n_procs_read.GT.0.AND.n_procs_read.LT.ppm_nproc) THEN
              ppm_nproc=n_procs_read
              IF (Ngrid(ppm_rc_dim).LE.ppm_nproc) THEN
                 ppm_nproc=Ngrid(ppm_rc_dim)-1
              ENDIF
           ELSE IF (n_procs_read.GE.ppm_nproc) THEN
              IF (Ngrid(ppm_rc_dim).LE.ppm_nproc) THEN
                 ppm_nproc=Ngrid(ppm_rc_dim)-1
              ENDIF
           ELSE
              IF (ppm_nproc.GE.Ngrid(ppm_rc_dim)) THEN
                 ppm_nproc=Ngrid(ppm_rc_dim)-1
              ENDIF
           ENDIF !(n_procs_read.GT.1)
        END SELECT

        DEALLOCATE(ioghostsize,STAT=info)
        or_fail_alloc("ioghostsize")

#ifdef __Linux
        !---------------------------------------------------------------------
        ! We can check the memory consumption for debugging purposes.
        ! This one is only working on LINUX
        !---------------------------------------------------------------------
        IF (debug.GT.1) THEN
           CALL ppm_rc_mem_usage(memory,info)
           stdout("mem_usage at program initialization=",memory)
        ENDIF
#endif

        !---------------------------------------------------------------------
        ! Create new topology suitable for reading and writing images.
        !---------------------------------------------------------------------
        assig  = ppm_param_assign_internal
        decomp = MERGE(ppm_param_decomp_xpencil,ppm_param_decomp_xy_slab,ppm_rc_dim.EQ.2)

        ! creating dummy particle
        NULLIFY(xp)
        iopt=ppm_param_alloc_fit
        ldu(1)=__DIME
        ldu(2)=1
        CALL ppm_alloc(xp,ldu,iopt,info)
        or_fail_alloc("xp")
        xp=zero

        !-------------------------------------------------------------------------
        ! second : Create IO topology
        !-------------------------------------------------------------------------
        CALL ppm_mktopo(iotopoid,iomeshid,xp,0,decomp,assig,min_phys, &
        &    max_phys,bcdef,ghostsize,cost,Ngrid(1:__DIME),info)
        or_fail('Failed to create new topology.')

        ! Get rid of the dummy particle
        iopt=ppm_param_dealloc
        CALL ppm_alloc(xp,ldu,iopt,info)
        or_fail_dealloc("xp")

        IF (nproc.GT.1) THEN
           IF (ppm_nproc.EQ.1) THEN
              IF (rank.GT.0) THEN
                 topo => ppm_topo(iotopoid)%t
                 topo%nsublist=0
                 topo%sub2proc=0
              ENDIF
           ENDIF
        ENDIF

        !-------------------------------------------------------------------------
        ! third : correct the total number of processors to the real numbers
        !-------------------------------------------------------------------------
        ppm_nproc=nproc

        !-------------------------------------------------------------------------
        ! IO mesh to read and write image file
        !-------------------------------------------------------------------------
        iomesh => ppm_mesh%at(iomeshid)

        !-------------------------------------------------------------------------
        ! Create uniform mesh over the domain
        !-------------------------------------------------------------------------
        CALL iomesh%def_uniform(info)
        or_fail("Failed to create uniform iomesh.")

#ifdef __Linux
        IF (debug.GT.1) THEN
           CALL ppm_rc_mem_usage(memory,info)
           stdout("mem_usage at topo and IO mesh are created=",memory)
        ENDIF
#endif

        !---------------------------------------------------------------------
        ! Craete one image from input image files and exit the program
        ! This one is for creating one stack from some single inputs
        !---------------------------------------------------------------------
        SELECT CASE (createoneimage)
        CASE (.TRUE.)
           !---------------------------------------------------------------------
           ! Define image field variable to read the image pixels intensity
           ! This field has a real type
           !---------------------------------------------------------------------
           ALLOCATE(ppm_t_field::image,STAT=info)
           or_fail_alloc('Failed to allocate field.')

           !---------------------------------------------------------------------
           ! Create a field for image density
           ! image is a scalar field with a real type
           !---------------------------------------------------------------------
           CALL image%create(1,info,name="image_pixel_intensity")
           or_fail("Create field failed!" )

           CALL image%discretize_on(iomesh,info)
           or_fail("image discretize_on failed!")
           ! Descretize the field over the mesh

           CALL DTYPE(ppm_rc_read_write_image)(image,iomesh,inputimage,inputimage,info)
           or_fail('Failed to read image files and write them into one image file.')

           !To stop the program
           info=exit_gracefully
           GOTO 9999
        END SELECT

        !---------------------------------------------------------------------
        ! Craete one unique label from input label(s) file(s) and exit the program
        !---------------------------------------------------------------------
        SELECT CASE (createuniquelabel)
        CASE (.TRUE.)
           IF (LGT(TRIM(initimage),"")) THEN
              ninputimage=ninitimage

              !---------------------------------------------------------------------
              ! Define labels field variable to assign label
              ! to the image pixels
              !---------------------------------------------------------------------
              ALLOCATE(ppm_t_field::labels,STAT=info)
              or_fail_alloc('Failed to allocate labels field.')

              CALL labels%create(1,info,dtype=ppm_type_int,name="pixel_labels")
              or_fail("Create labels field failed!" )

              ! Descretize the field over the iomesh
              CALL labels%discretize_on(iomesh,info)
              or_fail("labels discretize_on failed!")

              CALL ppm_rc_read_image_info(initimage,ninputimage,info)
              or_fail('Failed to read initimage information.')

              CALL DTYPE(ppm_rc_read_write_unique_label)(labels,iomesh,initimage,initimage,info)
              or_fail('Failed to read label files and write them into one image file.')

              !To stop the program
              info=exit_gracefully
              GOTO 9999
           ELSE
              fail("The initimage filename is not available!",ppm_error=ppm_error_fatal)
           ENDIF
        END SELECT

        !---------------------------------------------------------------------
        ! Define image field variable to read the image pixels intensity
        ! This field has a real type
        !---------------------------------------------------------------------
        ALLOCATE(ppm_t_field::image,STAT=info)
        or_fail_alloc('Failed to allocate field.')

        !---------------------------------------------------------------------
        ! Create a field for image density
        ! image is a scalar field with a real type
        !---------------------------------------------------------------------
        CALL image%create(1,info,name="image_pixel_intensity")
        or_fail("Create field failed!" )

        NULLIFY(imgdinfo,lbldinfo,initdinfo)

        !---------------------------------------------------------------------
        ! Descretize the scalar field (image) on IO mesh
        !---------------------------------------------------------------------
        CALL image%discretize_on(iomesh,info,discr_info=imgdinfo)
        or_fail("image discretize_on failed!")
        ! Descretize the field over the mesh

#ifdef __Linux
        IF (debug.GT.1) THEN
           CALL ppm_rc_mem_usage(memory,info)
           stdout("mem_usage at image is created and descretized on IO topo=",memory)
        ENDIF
#endif

        !---------------------------------------------------------------------
        ! Read the input file
        ! Store the intensity on image field which is defined on iomesh
        ! as a real variable
        !---------------------------------------------------------------------
        CALL DTYPE(ppm_rc_read_image)(image,iomesh,inputimage,info)
        or_fail('Failed to read image.')

        !---------------------------------------------------------------------
        ! In case we have another input file for initialization
        ! We create the labels field
        !---------------------------------------------------------------------
        IF (LGT(TRIM(initimage),"")) THEN
           ninputimage=ninitimage

           !---------------------------------------------------------------------
           ! Define labels field variable to assign label
           ! to the image pixels
           !---------------------------------------------------------------------
           ALLOCATE(ppm_t_field::labels,STAT=info)
           or_fail_alloc('Failed to allocate labels field.')

           CALL labels%create(1,info,dtype=ppm_type_int,name="pixel_labels")
           or_fail("Create labels field failed!" )

           ! Descretize the field over the iomesh
           CALL labels%discretize_on(iomesh,info,discr_info=lbldinfo)
           or_fail("labels discretize_on failed!")

           CALL ppm_rc_read_image_info(initimage,ninputimage,info)
           or_fail('Failed to read initimage information.')

           CALL DTYPE(ppm_rc_read_image)(labels,iomesh,initimage,info)
           or_fail('Failed to read initimage.')
        ENDIF

        !---------------------------------------------------------------------
        ! Create new topology suitable for the initialization.
        ! Use bisection decomposition strategy and map the image.
        !---------------------------------------------------------------------
        SELECT CASE (ppm_nproc)
        CASE (1)
           initopoid=iotopoid

           imesh => iomesh

           initmeshid = iomeshid

           IF (.NOT.ASSOCIATED(lbldinfo)) THEN
              !---------------------------------------------------------------------
              ! Define labels field variable to assign label
              ! the image pixels
              !---------------------------------------------------------------------
              ALLOCATE(ppm_t_field::labels,STAT=info)
              or_fail_alloc('Failed to allocate labels field.')

              CALL labels%create(1,info,dtype=ppm_type_int,name="pixel_labels")
              or_fail("Create labels field failed!" )

              CALL labels%discretize_on(imesh,info)
              or_fail("labels discretize_on failed!")
              !discretize a label field over the new initialization mesh
           ENDIF

           NULLIFY(imgdinfo,lbldinfo)
        CASE DEFAULT
           IF (IAND(ppm_nproc,ppm_nproc-1).NE.0) THEN
              WRITE(cbuf,'(A,A)') "The current version of RC client is only ", &
              & "working with the number of processors equlas to any power of 2"
              fail(cbuf,ppm_error=ppm_error_fatal)
           ENDIF

           !!! TODO
           !!! In the new update I should get rid of this
           !!! Currently I disbaled initopo map to topo
           IF (nsteql.GT.0) THEN
              ghostsize=MAX(inighostsize,ghostsize_run,ghostsize_equil)
           ELSE
              ghostsize=MAX(inighostsize,ghostsize_run)
           ENDIF

           assig  = ppm_param_assign_internal
           decomp = ppm_param_decomp_cartesian
           !ppm_param_decomp_bisection

           !Dummy particle
           !It is only necessary for topology createaion
           NULLIFY(xp)
           iopt=ppm_param_alloc_fit
           ldu(1)=__DIME
           ldu(2)=1
           CALL ppm_alloc(xp,ldu,iopt,info)
           or_fail_alloc("xp")
           xp=zero

           CALL ppm_mktopo(initopoid,initmeshid,xp,0,decomp,assig,min_phys, &
           &    max_phys,bcdef,ghostsize,cost,Ngrid(1:__DIME),info)
           or_fail('Failed to create new topology.')

           imesh => ppm_mesh%at(initmeshid)
           ! Init mesh

           CALL imesh%def_uniform(info)
           or_fail("Failed to create uniform mesh.")
           ! Create uniform initlization mesh over the domain

#ifdef __Linux
           IF (debug.GT.1) THEN
              CALL ppm_rc_mem_usage(memory,info)
              stdout("mem_usage at Init topo and Init mesh are created=",memory)
           ENDIF
#endif
           !---------------------------------------------------------------------
           ! Descretize the scalar field (image) on Init mesh
           !---------------------------------------------------------------------
           CALL image%discretize_on(imesh,info,discr_info=initdinfo)
           or_fail("image discretize_on failed!")

#ifdef __Linux
           IF (debug.GT.1) THEN
              CALL ppm_rc_mem_usage(memory,info)
              stdout("mem_usage at image is descretized on Init topo=",memory)
           ENDIF
#endif

           !---------------------------------------------------------------------
           ! Global mapping to the new mesh in the new topology
           ! the mapping is between IO and initialization topology
           !---------------------------------------------------------------------
           CALL iomesh%map(imesh,info)
           or_fail("Failed to do global mapping from old to the new topology.")
           !---------------------------------------------------------------------
           ! Push image data in the buffer
           !---------------------------------------------------------------------
           CALL image%map_push(iomesh,info)
           or_fail("Failed to push image data in the buffer!")
#ifdef __Linux
           IF (debug.GT.1) THEN
              CALL ppm_rc_mem_usage(memory,info)
              stdout("mem_usage at global mapping to Init topo=",memory)
           ENDIF
#endif
           !---------------------------------------------------------------------
           ! Start of the Non-blocking send/receive of the buffer
           !---------------------------------------------------------------------
           CALL iomesh%map_isend(info,sendrecv=.TRUE.)
           or_fail("Failed to send image data to the new mesh!")

           !---------------------------------------------------------------------
           ! Get rid of the dummy particle
           !---------------------------------------------------------------------
           iopt=ppm_param_dealloc
           CALL ppm_alloc(xp,ldu,iopt,info)
           or_fail_dealloc("xp")
           NULLIFY(xp)

           !---------------------------------------------------------------------
           ! Get rid of the object which will not be used anymore to free memory
           !---------------------------------------------------------------------
           CALL image%discr_info%remove(info,imgdinfo)
           or_fail("image%discr_info%remove")

           DEALLOCATE(imgdinfo,STAT=info)
           or_fail_dealloc("Failed to deallocate discr_info")

           !We need to free memory which will not be used later
           !image data on IO mesh will not be needed anymore
           NULLIFY(imgdinfo)

           !---------------------------------------------------------------------
           ! Wait to make sure we have the buffer
           !---------------------------------------------------------------------
           CALL iomesh%map_isend(info,sendrecv=.FALSE.)
           or_fail("Failed to send image data to the new mesh!")
           !---------------------------------------------------------------------
           ! Pop received buffer (image) data in the new mesh
           !---------------------------------------------------------------------
           CALL image%map_pop(imesh,info)
           or_fail("Failed to pop image data in the new mesh!")

#ifdef __Linux
           IF (debug.GT.1) THEN
              CALL ppm_rc_mem_usage(memory,info)
              stdout("mem_usage at destroy the descretized image on IO mesh=",memory)
           ENDIF
#endif

           IF (ASSOCIATED(lbldinfo)) THEN
              !---------------------------------------------------------------------
              ! Global mapping to the new mesh in the new topology
              ! the mapping is between IO and initialization topology
              !---------------------------------------------------------------------
              !---------------------------------------------------------------------
              ! Push labels data in the buffer
              !---------------------------------------------------------------------
              CALL labels%map_push(iomesh,info)
              or_fail("Failed to push labels data in the buffer!")
              !---------------------------------------------------------------------
              ! Start of the Non-blocking send/receive of the buffer
              !---------------------------------------------------------------------
              CALL iomesh%map_isend(info,sendrecv=.TRUE.)
              or_fail("Failed to send image data to the new mesh!")

              CALL labels%discretize_on(imesh,info)
              or_fail("labels discretize_on failed!")
              !discretize a label field over the new initialization mesh

#ifdef __Linux
              IF (debug.GT.1) THEN
                 CALL ppm_rc_mem_usage(memory,info)
                 stdout("mem_usage at descretizing labels on Init mesh=",memory)
              ENDIF
#endif

              CALL labels%discr_info%remove(info,lbldinfo)
              or_fail("labels%discr_info%remove")

              DEALLOCATE(lbldinfo,STAT=info)
              or_fail_dealloc("Failed to deallocate discr_info")

              !We need to free memory which will not be used later
              !labels data on IO mesh will not be needed anymore
              NULLIFY(lbldinfo)

#ifdef __Linux
              IF (debug.GT.1) THEN
                 CALL ppm_rc_mem_usage(memory,info)
                 stdout("mem_usage at removing discr_info from IO mesh=",memory)
              ENDIF
#endif

              !---------------------------------------------------------------------
              ! Wait to make sure we have the buffer
              !---------------------------------------------------------------------
              CALL iomesh%map_isend(info,sendrecv=.FALSE.)
              or_fail("Failed to send image data to the new mesh!")
              !---------------------------------------------------------------------
              ! Pop receiveed buffer (image) data in the new mesh
              !---------------------------------------------------------------------
              CALL labels%map_pop(imesh,info)
              or_fail("Failed to pop labels data in the new mesh!")

              !---------------------------------------------------------------------
              ! Destroy the IO mesh
              !---------------------------------------------------------------------
              CALL iomesh%destroy(info)
              or_fail('Failed to destroy iomesh.')

              dealloc_pointer("iomesh")

#ifdef __Linux
              IF (debug.GT.1) THEN
                 CALL ppm_rc_mem_usage(memory,info)
                 stdout("mem_usage at Io mesh destroy=",memory)
              ENDIF
#endif

              !---------------------------------------------------------------------
              ! Destroy the IO topo
              !---------------------------------------------------------------------
              topo => ppm_topo(iotopoid)%t

              iopt=ppm_param_dealloc
              CALL ppm_alloc(topo,iopt,info)
              or_fail_dealloc("Failed to deallocate topology (iotopoid)!")

              NULLIFY(ppm_topo(iotopoid)%t)

#ifdef __Linux
              IF (debug.GT.1) THEN
                 CALL ppm_rc_mem_usage(memory,info)
                 stdout("mem_usage at Io topo destroy=",memory)
              ENDIF
#endif
           ELSE
              !If we do not have labels already descretized on IO mesh, we can first destroy
              !IO mesh and free memory, then descretize labels on new topo and mesh
              !---------------------------------------------------------------------
              ! Destroy the IO mesh
              !---------------------------------------------------------------------
              CALL iomesh%destroy(info)
              or_fail('Failed to destroy iomesh.')

              dealloc_pointer("iomesh")

#ifdef __Linux
              IF (debug.GT.1) THEN
                 CALL ppm_rc_mem_usage(memory,info)
                 stdout("mem_usage at Io mesh destroy=",memory)
              ENDIF
#endif
              !---------------------------------------------------------------------
              ! Destroy the IO topo
              !---------------------------------------------------------------------
              topo => ppm_topo(iotopoid)%t

              iopt=ppm_param_dealloc
              CALL ppm_alloc(topo,iopt,info)
              or_fail_dealloc("Failed to deallocate topology (iotopoid)!")

              NULLIFY(ppm_topo(iotopoid)%t)

#ifdef __Linux
              IF (debug.GT.1) THEN
                 CALL ppm_rc_mem_usage(memory,info)
                 stdout("mem_usage at Io topo destroy=",memory)
              ENDIF
#endif

              !---------------------------------------------------------------------
              ! Define labels field variable to assign label
              ! the image pixels
              !---------------------------------------------------------------------
              ALLOCATE(ppm_t_field::labels,STAT=info)
              or_fail_alloc('Failed to allocate labels field.')

              CALL labels%create(1,info,dtype=ppm_type_int,name="pixel_labels")
              or_fail("Create labels field failed!" )

              CALL labels%discretize_on(imesh,info)
              or_fail("labels discretize_on failed!")
              !discretize a label field over the new initialization mesh

#ifdef __Linux
              IF (debug.GT.1) THEN
                 CALL ppm_rc_mem_usage(memory,info)
                 stdout("mem_usage at descretizing labels on Init mesh=",memory)
              ENDIF
#endif
           ENDIF

        END SELECT

        !-------------------------------------------------------------------------
        ! Normalize and Shift
        !-------------------------------------------------------------------------
        SELECT CASE (lNormalize)
        CASE (.TRUE.)
           !-------------------------------------------------------------------------
           ! Normalize and Shift the input image to [0,1]
           !-------------------------------------------------------------------------
           CALL DTYPE(ppm_rc_CopyImageAndNormalize)(image,image,imesh,info, &
           &    Normalize=.TRUE.,FieldMinVal=FieldMin,FieldMaxVal=FieldMax)

           ImageNormalfac=FieldMax-FieldMin
           MinShiftVal=FieldMin

           IF (debug.GT.1) THEN
              IF (rank.EQ.0) THEN
                 stdout("ImageNormalfac=",ImageNormalfac)
                 stdout("MinShiftVal=",MinShiftVal)
              ENDIF
           ENDIF

           !-------------------------------------------------------------------------
           ! Correct the tolerance, when we are normalizing the image
           !-------------------------------------------------------------------------
           Tolerance_=REAL(Tolerance,MK)/ImageNormalfac

           IF (intensity_ths.GT.one) THEN
              intensity_ths=intensity_ths/ImageNormalfac
           ENDIF
        CASE DEFAULT
           Tolerance_=REAL(Tolerance,MK)
        END SELECT

        !-------------------------------------------------------------------------
        !  Initialize init related terms
        !-------------------------------------------------------------------------
        ALLOCATE(e_Init,STAT=info)
        or_fail_alloc("e_Init")

        ldu=FLOOR(init_rd)
        ld=FLOOR(init_sp)

        CALL e_Init%CreateElement(ppm_rc_dim,info,ForegroundValue=1, &
        &    BackgroundValue=0,size_=ldu,spacing_=ld)
        or_fail("e_Init%CreateElement")
        !create initialization elements

        CALL e_Init%DTYPE(GetOutput)(image,imesh,labels,info)
        or_fail("e_Init%GetOutput")
        !create initialization regions

        !-------------------------------------------------------------------------
        !  Updating the Ghost information on every processor
        !-------------------------------------------------------------------------
        IF (ppm_nproc.GT.1) THEN
           !-------------------------------------------------------------------------
           !  Ghost get
           !-------------------------------------------------------------------------
           CALL imesh%map_ghost_get(info,ghostsize=inighostsize)
           or_fail("imesh%map_ghost_get")
           !---------------------------------------------------------------------
           ! Push labels data for ghost update in the buffer
           !---------------------------------------------------------------------
           CALL labels%map_ghost_push(imesh,info)
           or_fail("labels%map_ghost_push")
           !---------------------------------------------------------------------
           ! Start of the Non-blocking send/receive of the buffer
           !---------------------------------------------------------------------
           CALL imesh%map_isend(info,sendrecv=.TRUE.)
           or_fail("imesh%map_isend")
        ENDIF

        CALL e_Init%DestroyElement(info)
        or_fail("e_Init%DestroyElement")
        ! destroy the initialization element and assigned memory

        DEALLOCATE(e_Init,inighostsize,STAT=info)
        or_fail_dealloc("e_Init & inighostsize")

        NULLIFY(e_Init)

        IF (debug.GT.1) THEN
           CALL DTYPE(ppm_rc_write_image_label)(labels,imesh,"afterput",info,idn=FORBIDDEN)
           or_fail("ppm_rc_write_image_label")
        ENDIF

        IF (ppm_nproc.GT.1) THEN
           !---------------------------------------------------------------------
           ! Wait to make sure we have the buffer
           !---------------------------------------------------------------------
           CALL imesh%map_isend(info,sendrecv=.FALSE.)
           or_fail("imesh%map_isend")
           !---------------------------------------------------------------------
           ! Pop received buffer (labels) data in the ghost layer
           !---------------------------------------------------------------------
           CALL labels%map_ghost_pop(imesh,info)
           or_fail("labels%map_ghost_pop")
        ENDIF

#ifdef __Linux
        IF (debug.GT.1) THEN
           CALL ppm_rc_mem_usage(memory,info)
           stdout("mem_usage at Init complete=",memory)
        ENDIF
#endif

        !-------------------------------------------------------------------------
        !
        !-------------------------------------------------------------------------
        ipatch=imesh%subpatch%nb
        ld(1)=SUM(ppm_rc_seeds(1:ipatch)%nb)

        ALLOCATE(nlabels(ld(1)),STAT=info)
        or_fail_alloc("nlabels")

        !Assign the new label for each seed region
        !which is correspond to one intial region
        DO iseed=1,ld(1)
           nlabels(iseed)=loc_label
           loc_label=loc_label-1
        ENDDO

        sbpitr => imesh%subpatch%begin()
        nsize=0
        DO WHILE (ASSOCIATED(sbpitr))
           !The size of ghost layer around each subpatch
           !we compute only one layer of ghost around each patch
           Nm => sbpitr%nnodes

           IF (ALL(sbpitr%istart(1:ppm_rc_dim).EQ.1.AND. &
           &       sbpitr%iend(1:ppm_rc_dim).EQ.Ngrid(1:__DIME))) THEN
              !do nothing
           ELSE
#if   __DIME == __2D
              IF (sbpitr%istart(1).NE.1) THEN
                 nsize=nsize+Nm(2)
              ENDIF
              IF (sbpitr%iend(1).NE.Ngrid(1)) THEN
                 nsize=nsize+Nm(2)
              ENDIF
              IF (sbpitr%istart(2).NE.1) THEN
                 nsize=nsize+Nm(1)+2
              ENDIF
              IF (sbpitr%iend(2).NE.Ngrid(2)) THEN
                 nsize=nsize+Nm(1)+2
              ENDIF
              !nsize=nsize+(Nm(2)+Nm(1)+2)*2
#elif __DIME == __3D
              IF (sbpitr%istart(1).NE.1) THEN
                 nsize=nsize+Nm(2)*Nm(3)
              ENDIF
              IF (sbpitr%iend(1).NE.Ngrid(1)) THEN
                 nsize=nsize+Nm(2)*Nm(3)
              ENDIF
              IF (sbpitr%istart(2).NE.1) THEN
                 nsize=nsize+Nm(3)*(Nm(1)+2)
              ENDIF
              IF (sbpitr%iend(2).NE.Ngrid(2)) THEN
                 nsize=nsize+Nm(3)*(Nm(1)+2)
              ENDIF
              IF (sbpitr%istart(3).NE.1) THEN
                 nsize=nsize+(Nm(1)+2)*(Nm(2)+2)
              ENDIF
              IF (sbpitr%iend(3).NE.Ngrid(3)) THEN
                 nsize=nsize+(Nm(1)+2)*(Nm(2)+2)
              ENDIF
              !nsize=nsize+(Nm(2)*Nm(3)+(Nm(1)+2)*Nm(3)+(Nm(1)+2)*(Nm(2)+2))*2
#endif
           ENDIF

           sbpitr => imesh%subpatch%next()
        ENDDO

        CALL DTYPE(ppm_rc_initforestfire)(imesh,nsize,info)
        or_fail("ppm_rc_forestfire")

        DEALLOCATE(nlabels,STAT=info)
        or_fail_dealloc("nlabels")

        !-------------------------------------------------------------------------
        !  Now recreate topology and map to the new topology
        !-------------------------------------------------------------------------
        IF (nsteql.GT.0) THEN
           ghostsize=MAX(ghostsize_equil,ghostsize_run)
        ELSE
           ghostsize=ghostsize_run
        ENDIF

        SELECT CASE (ppm_nproc)
        CASE (1)
           topoid = initopoid

           mesh => iomesh

           meshid = iomeshid

           NULLIFY(imesh)

           NULLIFY(wpil,wpr)

           sbpitr => mesh%subpatch%begin()
           DO WHILE (ASSOCIATED(sbpitr))
              Nm   => sbpitr%nnodes
              hi_a => sbpitr%hi_a
              lo_a => sbpitr%lo_a

              CALL sbpitr%get_field(labels,wpil,info)
              or_fail("Failed to get field wpil data.")

              CALL sbpitr%get_field(image,wpr,info)
              or_fail("Failed to get field wpr data.")

#if   __DIME == __2D
              IF (sbpitr%istart(1).EQ.1) THEN
                 DO j=lo_a(2),hi_a(2)
                    DO i=lo_a(1),1
                       wpil(i,j)=FORBIDDEN
                    ENDDO
                 ENDDO
                 DO j=lo_a(2),hi_a(2)
                    DO i=lo_a(1),0
                       wpr(i,j)=zero
                    ENDDO
                 ENDDO
              ENDIF
              IF (sbpitr%istart(2).EQ.1) THEN
                 DO j=lo_a(2),1
                    DO i=lo_a(1),hi_a(1)
                       wpil(i,j)=FORBIDDEN
                    ENDDO
                 ENDDO
                 DO j=lo_a(2),0
                    DO i=lo_a(1),hi_a(1)
                       wpr(i,j)=zero
                    ENDDO
                 ENDDO
              ENDIF
              IF (sbpitr%iend(1).EQ.Ngrid(1)) THEN
                 DO j=lo_a(2),hi_a(2)
                    DO i=Nm(1),hi_a(1)
                       wpil(i,j)=FORBIDDEN
                    ENDDO
                 ENDDO
                 DO j=lo_a(2),hi_a(2)
                    DO i=Nm(1)+1,hi_a(1)
                       wpr(i,j)=zero
                    ENDDO
                 ENDDO
              ENDIF
              IF (sbpitr%iend(2).EQ.Ngrid(2)) THEN
                 DO j=Nm(2),hi_a(2)
                    DO i=lo_a(1),hi_a(1)
                       wpil(i,j)=FORBIDDEN
                    ENDDO
                 ENDDO
                 DO j=Nm(2)+1,hi_a(2)
                    DO i=lo_a(1),hi_a(1)
                       wpr(i,j)=zero
                    ENDDO
                 ENDDO
              ENDIF
#elif __DIME == __3D
              IF (sbpitr%istart(1).EQ.1) THEN
                 DO k=lo_a(3),hi_a(3)
                    DO j=lo_a(2),hi_a(2)
                       DO i=lo_a(1),1
                          wpil(i,j,k)=FORBIDDEN
                       ENDDO
                    ENDDO
                 ENDDO
                 DO k=lo_a(3),hi_a(3)
                    DO j=lo_a(2),hi_a(2)
                       DO i=lo_a(1),0
                          wpr(i,j,k)=zero
                       ENDDO
                    ENDDO
                 ENDDO
              ENDIF
              IF (sbpitr%istart(2).EQ.1) THEN
                 DO k=lo_a(3),hi_a(3)
                    DO j=lo_a(2),1
                       DO i=lo_a(1),hi_a(1)
                          wpil(i,j,k)=FORBIDDEN
                       ENDDO
                    ENDDO
                 ENDDO
                 DO k=lo_a(3),hi_a(3)
                    DO j=lo_a(2),0
                       DO i=lo_a(1),hi_a(1)
                          wpr(i,j,k)=zero
                       ENDDO
                    ENDDO
                 ENDDO
              ENDIF
              IF (sbpitr%istart(3).EQ.1) THEN
                 DO k=lo_a(3),1
                    DO j=lo_a(2),hi_a(2)
                       DO i=lo_a(1),hi_a(1)
                          wpil(i,j,k)=FORBIDDEN
                       ENDDO
                    ENDDO
                 ENDDO
                 DO k=lo_a(3),0
                    DO j=lo_a(2),hi_a(2)
                       DO i=lo_a(1),hi_a(1)
                          wpr(i,j,k)=zero
                       ENDDO
                    ENDDO
                 ENDDO
              ENDIF
              IF (sbpitr%iend(1).EQ.Ngrid(1)) THEN
                 DO k=lo_a(3),hi_a(3)
                    DO j=lo_a(2),hi_a(2)
                       DO i=Nm(1),hi_a(1)
                          wpil(i,j,k)=FORBIDDEN
                       ENDDO
                    ENDDO
                 ENDDO
                 DO k=lo_a(3),hi_a(3)
                    DO j=lo_a(2),hi_a(2)
                       DO i=Nm(1)+1,hi_a(1)
                          wpr(i,j,k)=zero
                       ENDDO
                    ENDDO
                 ENDDO
              ENDIF
              IF (sbpitr%iend(2).EQ.Ngrid(2)) THEN
                 DO k=lo_a(3),hi_a(3)
                    DO j=Nm(2),hi_a(2)
                       DO i=lo_a(1),hi_a(1)
                          wpil(i,j,k)=FORBIDDEN
                       ENDDO
                    ENDDO
                 ENDDO
                 DO k=lo_a(3),hi_a(3)
                    DO j=Nm(2)+1,hi_a(2)
                       DO i=lo_a(1),hi_a(1)
                          wpr(i,j,k)=zero
                       ENDDO
                    ENDDO
                 ENDDO
              ENDIF
              IF (sbpitr%iend(3).EQ.Ngrid(3)) THEN
                 DO k=Nm(3),hi_a(3)
                    DO j=lo_a(2),hi_a(2)
                       DO i=lo_a(1),hi_a(1)
                          wpil(i,j,k)=FORBIDDEN
                       ENDDO
                    ENDDO
                 ENDDO
                 DO k=Nm(3)+1,hi_a(3)
                    DO j=lo_a(2),hi_a(2)
                       DO i=lo_a(1),hi_a(1)
                          wpr(i,j,k)=zero
                       ENDDO
                    ENDDO
                 ENDDO
              ENDIF
#endif

              sbpitr => mesh%subpatch%next()
           ENDDO !WHILE(ASSOCIATED(sbpitr))

           NULLIFY(wpil,wpr)

!            DO i=1,108
           !check if all went well on image
           CALL DTYPE(ppm_rc_write_image)(image,mesh,"image",info)
           or_fail("ppm_rc_write_image")
!            ENDDO

#ifdef __MPI
        CASE DEFAULT
           NULLIFY(wpil)

           sbpitr => imesh%subpatch%begin()
           ipart=0
           DO WHILE (ASSOCIATED(sbpitr))
              Nm   => sbpitr%nnodes

              CALL sbpitr%get_field(labels,wpil,info)
              or_fail("Failed to get field wpil data.")

#if   __DIME == __2D
              IF (sbpitr%istart(1).EQ.1) THEN
                 DO j=1,Nm(2)
                    wpil(1,j)=FORBIDDEN
                 ENDDO
              ENDIF
              IF (sbpitr%istart(2).EQ.1) THEN
                 DO i=1,Nm(1)
                    wpil(i,1)=FORBIDDEN
                 ENDDO
              ENDIF
              IF (sbpitr%iend(1).EQ.Ngrid(1)) THEN
                 DO j=1,Nm(2)
                    wpil(Nm(1),j)=FORBIDDEN
                 ENDDO
              ENDIF
              IF (sbpitr%iend(2).EQ.Ngrid(2)) THEN
                 DO i=1,Nm(1)
                    wpil(i,Nm(2))=FORBIDDEN
                 ENDDO
              ENDIF

!               DO j=1,Nm(2)
!                  DO i=1,Nm(1)
!                     IF (wpil(i,j).LT.0) THEN
!                        ipart=ipart+1
!                     ENDIF
!                  ENDDO
!               ENDDO
#elif __DIME == __3D
              IF (sbpitr%istart(1).EQ.1) THEN
                 DO k=1,Nm(3)
                    DO j=1,Nm(2)
                       wpil(1,j,k)=FORBIDDEN
                    ENDDO
                 ENDDO
              ENDIF
              IF (sbpitr%istart(2).EQ.1) THEN
                 DO k=1,Nm(3)
                    DO i=1,Nm(1)
                       wpil(i,1,k)=FORBIDDEN
                    ENDDO
                 ENDDO
              ENDIF
              IF (sbpitr%istart(3).EQ.1) THEN
                 DO j=1,Nm(2)
                    DO i=1,Nm(1)
                       wpil(i,j,1)=FORBIDDEN
                    ENDDO
                 ENDDO
              ENDIF
              IF (sbpitr%iend(1).EQ.Ngrid(1)) THEN
                 DO k=1,Nm(3)
                    DO j=1,Nm(2)
                       wpil(Nm(1),j,k)=FORBIDDEN
                    ENDDO
                 ENDDO
              ENDIF
              IF (sbpitr%iend(2).EQ.Ngrid(2)) THEN
                 DO k=1,Nm(3)
                    DO i=1,Nm(1)
                       wpil(i,Nm(2),k)=FORBIDDEN
                    ENDDO
                 ENDDO
              ENDIF
              IF (sbpitr%iend(3).EQ.Ngrid(3)) THEN
                 DO j=1,Nm(2)
                    DO i=1,Nm(1)
                       wpil(i,j,Nm(3))=FORBIDDEN
                    ENDDO
                 ENDDO
              ENDIF

!               DO k=1,Nm(3)
!                  DO j=1,Nm(2)
!                     DO i=1,Nm(1)
!                        IF (wpil(i,j,k).LT.0) THEN
!                           ipart=ipart+1
!                        ENDIF
!                     ENDDO !i=1,Nm(1)
!                  ENDDO !i=1,Nm(2)
!               ENDDO !i=1,Nm(3)
#endif

              sbpitr => imesh%subpatch%next()
           ENDDO !WHILE(ASSOCIATED(sbpitr))

           NULLIFY(wpil,xp)

!            iopt=ppm_param_alloc_fit
!            ldu(1)=__DIME
!            ldu(2)=ipart
!            CALL ppm_alloc(xp,ldu,iopt,info)
!            or_fail_alloc("xp")

!            sbpitr => imesh%subpatch%begin()
!            ipart=0
!            DO WHILE (ASSOCIATED(sbpitr))
!               Nm   => sbpitr%nnodes
!
!               CALL sbpitr%get_field(labels,wpil,info)
!               or_fail("Failed to get field wpil data.")
!
! #if   __DIME == __2D
!               DO j=1,Nm(2)
!                  DO i=1,Nm(1)
!                     IF (wpil(i,j).LT.0) THEN
!                        ipart=ipart+1
!                        xp(1,ipart)=REAL(i+sbpitr%istart(1)-2,MK)
!                        xp(2,ipart)=REAL(j+sbpitr%istart(2)-2,MK)
!                     ENDIF
!                  ENDDO
!               ENDDO
! #elif __DIME == __3D
!               DO k=1,Nm(3)
!                  DO j=1,Nm(2)
!                     DO i=1,Nm(1)
!                        IF (wpil(i,j,k).LT.0) THEN
!                           ipart=ipart+1
!                           xp(1,ipart)=REAL(i+sbpitr%istart(1)-2,MK)
!                           xp(2,ipart)=REAL(j+sbpitr%istart(2)-2,MK)
!                           xp(3,ipart)=REAL(k+sbpitr%istart(3)-2,MK)
!                        ENDIF
!                     ENDDO
!                  ENDDO
!               ENDDO
! #endif
!
!               sbpitr => imesh%subpatch%next()
!            ENDDO !WHILE(ASSOCIATED(sbpitr))

!            topoid = 0
!            meshid =-1
!            assig  = ppm_param_assign_internal
!            decomp = ppm_param_decomp_bisection
!
!            CALL ppm_mktopo(topoid,meshid,xp,ipart,decomp,assig,min_phys, &
!            &    max_phys,bcdef,ghostsize,cost,Ngrid(1:__DIME),info)
!            or_fail('Failed to create new topology.')

           topoid = initopoid
           topo => ppm_topo(topoid)%t

           IF (topo%nsubs.NE.ppm_nproc) THEN
              fail("This implementation is only works with one subdomain per processor", &
              & ppm_error=ppm_error_fatal)
           ENDIF

           ALLOCATE(nneigh(topo%nsubs),SOURCE=0,STAT=info)
           or_fail_alloc("nneigh")

           CALL MPI_Iallgather(topo%nneighsubs(1),1,MPI_INTEGER,nneigh,1, &
           &    MPI_INTEGER,comm,request,info)
           or_fail_MPI("MPI_Iallgather")


           Offset=zerod

           ALLOCATE(ppm_t_equi_mesh::mesh,STAT=info)
           or_fail_alloc('Failed to allocate mesh pointer.')

           CALL mesh%create(topoid,Offset,info, &
           &    Nm=Ngrid(1:__DIME),ghostsize=ghostsize)
           or_fail("Failed to create mesh.")

!            mesh => ppm_mesh%at(meshid)

           CALL mesh%def_uniform(info)
           or_fail("Failed to create uniform mesh.")

           CALL image%discretize_on(mesh,info)
           or_fail("image discretize_on failed!")

           !---------------------------------------------------------------------
           ! Global mapping to the new mesh in the new topology
           !---------------------------------------------------------------------
           CALL imesh%map(mesh,info)
           or_fail("Failed to do global mapping from old to the new topology.")

           !---------------------------------------------------------------------
           ! starting image global mapping to the new mesh in the new topology
           !---------------------------------------------------------------------
           CALL image%map_push(imesh,info)
           or_fail("Failed to push image data in the buffer!")
           CALL imesh%map_isend(info,sendrecv=.TRUE.)
           or_fail("Failed to send image data to the new mesh!")

           iopt=ppm_param_dealloc
           CALL ppm_alloc(xp,ldu,iopt,info)
           or_fail_dealloc("Failed to deallocate xp!")

           CALL labels%discretize_on(mesh,info)
           or_fail("labels discretize_on failed!")

           !---------------------------------------------------------------------
           ! complete image global mapping to the new mesh in the new topology
           !---------------------------------------------------------------------
           CALL imesh%map_isend(info,sendrecv=.FALSE.)
           or_fail("Failed to send image data to the new mesh!")
           CALL image%map_pop(mesh,info)
           or_fail("Failed to pop image data in the new mesh!")

           IF (ASSOCIATED(initdinfo)) THEN
              !---------------------------------------------------------------------
              ! Get rid of the object which will not be used anymore to free memory
              !---------------------------------------------------------------------
              CALL image%discr_info%remove(info,initdinfo)
              or_fail("image%discr_info%remove")

              DEALLOCATE(initdinfo,STAT=info)
              or_fail_dealloc("Failed to deallocate discr_info")

              !We need to free memory which will not be used later
              !image data on IO mesh will not be needed anymore
              NULLIFY(initdinfo)
           ENDIF

           !---------------------------------------------------------------------
           ! starting labels global mapping to the new mesh in the new topology
           !---------------------------------------------------------------------
           CALL labels%map_push(imesh,info)
           or_fail("Failed to push image data in the buffer!")
           CALL imesh%map_isend(info,sendrecv=.TRUE.)
           or_fail("Failed to send image data to the new mesh!")

           !wait for nneigh to be avilable on all ranks
           CALL MPI_Wait(request,MPI_STATUS_IGNORE,info)
           or_fail_MPI("MPI_Wait")

           !Now on every processor we have the nieghbor count
           nsize=SUM(nneigh)
           ALLOCATE(ineigh(nsize),displ(ppm_nproc),STAT=info)
           or_fail_alloc("ineigh")

           ineigh=0
           displ(1)=0
           DO i=2,ppm_nproc
              displ(i)=displ(i-1)+nneigh(i-1)
           ENDDO

           CALL MPI_Iallgatherv(topo%ineighsubs(1,1),topo%nneighsubs(1), &
           &    MPI_INTEGER,ineigh,nneigh,displ,MPI_INTEGER,comm,request,info)
           or_fail_MPI("MPI_Iallgather")

           !-------------------------------------------------------------------------
           ! check if all went well on image
           ! Writing the image based on the new decomposition
           !-------------------------------------------------------------------------
!            DO i=1,57
           CALL DTYPE(ppm_rc_write_image)(image,mesh,"image",info)
           or_fail("ppm_rc_write_image")
!            ENDDO

           !---------------------------------------------------------------------
           ! complete labels global mapping to the new mesh in the new topology
           !---------------------------------------------------------------------
           CALL imesh%map_isend(info,sendrecv=.FALSE.)
           or_fail("Failed to send image data to the new mesh!")
           CALL labels%map_pop(mesh,info)
           or_fail("Failed to pop image data in the new mesh!")



           !-------------------------------------------------------------------------
           ! ghost update on the new mesh at new topo
           ! ghostsize_run is set to ghostsize_equil during equilibrium steps
           !-------------------------------------------------------------------------
           CALL mesh%map_ghost_get(info,ghostsize=ghostsize_run)
           or_fail("mesh%map_ghost_get")
           !-------------------------------------------------------------------------
!---Start--! starting labels ghost update by pushing data and
           ! sending them with NONBLOCKING sends
           !-------------------------------------------------------------------------
           CALL labels%map_ghost_push(mesh,info)
           or_fail("labels%map_ghost_push")
           CALL mesh%map_isend(info,sendrecv=.TRUE.)
           or_fail("mesh%map_send")

           !---------------------------------------------------------------------
           ! Taking advantage of computation over communication
           ! taking care of stuff between sending and receiving messages
           !---------------------------------------------------------------------
#ifdef __Linux
           IF (debug.GT.1) THEN
              CALL ppm_rc_mem_usage(memory,info)
              stdout("mem_usage at ghost update on new meshid at new topoid=",memory)
           ENDIF
#endif
           !---------------------------------------------------------------------
           ! Destroy the initial imesh
           !---------------------------------------------------------------------
           CALL imesh%destroy(info)
           or_fail('Failed to destroy imesh.')

           dealloc_pointer("imesh")

!            !---------------------------------------------------------------------
!            ! Destroy the init topo
!            !---------------------------------------------------------------------
!            topo => ppm_topo(initopoid)%t
!
!            iopt=ppm_param_dealloc
!            CALL ppm_alloc(topo,iopt,info)
!            or_fail_dealloc("Failed to deallocate topology (initopoid)!")

!            NULLIFY(ppm_topo(initopoid)%t)

           !-------------------------------------------------------------------------
!---End----! complete ghost update by popping received data
           !-------------------------------------------------------------------------
           CALL mesh%map_isend(info,sendrecv=.FALSE.)
           or_fail("mesh%map_send")
           CALL labels%map_ghost_pop(mesh,info)
           or_fail("labels%map_ghost_pop")

           !-------------------------------------------------------------------------
!---Start--! starting image ghost update by pushing data and
           ! sending them with NONBLOCKING sends
           !-------------------------------------------------------------------------
           CALL image%map_ghost_push(mesh,info)
           or_fail("image%map_ghost_push")
           CALL mesh%map_isend(info,sendrecv=.TRUE.)
           or_fail("mesh%map_send")

           !---------------------------------------------------------------------
           ! Taking advantage of computation over communication
           ! taking care of stuff between sending and receiving messages
           !---------------------------------------------------------------------
           NULLIFY(wpil,wpr)

           sbpitr => mesh%subpatch%begin()
           DO WHILE (ASSOCIATED(sbpitr))
              Nm   => sbpitr%nnodes
              hi_a => sbpitr%hi_a
              lo_a => sbpitr%lo_a

              CALL sbpitr%get_field(labels,wpil,info)
              or_fail("Failed to get field i_wp data.")

              CALL sbpitr%get_field(image,wpr,info)
              or_fail("Failed to get field r_wp data.")

#if   __DIME == __2D
              IF (sbpitr%istart(1).EQ.1) THEN
                 DO j=lo_a(2),hi_a(2)
                    DO i=lo_a(1),0
                       wpil(i,j)=FORBIDDEN
                    ENDDO
                 ENDDO
                 DO j=lo_a(2),hi_a(2)
                    DO i=lo_a(1),0
                       wpr(i,j)=zero
                    ENDDO
                 ENDDO
              ENDIF
              IF (sbpitr%istart(2).EQ.1) THEN
                 DO j=lo_a(2),0
                    DO i=lo_a(1),hi_a(1)
                       wpil(i,j)=FORBIDDEN
                    ENDDO
                 ENDDO
                 DO j=lo_a(2),0
                    DO i=lo_a(1),hi_a(1)
                       wpr(i,j)=zero
                    ENDDO
                 ENDDO
              ENDIF
              IF (sbpitr%iend(1).EQ.Ngrid(1)) THEN
                 DO j=lo_a(2),hi_a(2)
                    DO i=Nm(1)+1,hi_a(1)
                       wpil(i,j)=FORBIDDEN
                    ENDDO
                 ENDDO
                 DO j=lo_a(2),hi_a(2)
                    DO i=Nm(1)+1,hi_a(1)
                       wpr(i,j)=zero
                    ENDDO
                 ENDDO
              ENDIF
              IF (sbpitr%iend(2).EQ.Ngrid(2)) THEN
                 DO j=Nm(2)+1,hi_a(2)
                    DO i=lo_a(1),hi_a(1)
                       wpil(i,j)=FORBIDDEN
                    ENDDO
                 ENDDO
                 DO j=Nm(2)+1,hi_a(2)
                    DO i=lo_a(1),hi_a(1)
                       wpr(i,j)=zero
                    ENDDO
                 ENDDO
              ENDIF
#elif __DIME == __3D
              IF (sbpitr%istart(1).EQ.1) THEN
                 DO k=lo_a(3),hi_a(3)
                    DO j=lo_a(2),hi_a(2)
                       DO i=lo_a(1),0
                          wpil(i,j,k)=FORBIDDEN
                       ENDDO
                    ENDDO
                 ENDDO
                 DO k=lo_a(3),hi_a(3)
                    DO j=lo_a(2),hi_a(2)
                       DO i=lo_a(1),0
                          wpr(i,j,k)=zero
                       ENDDO
                    ENDDO
                 ENDDO
              ENDIF
              IF (sbpitr%istart(2).EQ.1) THEN
                 DO k=lo_a(3),hi_a(3)
                    DO j=lo_a(2),0
                       DO i=lo_a(1),hi_a(1)
                          wpil(i,j,k)=FORBIDDEN
                       ENDDO
                    ENDDO
                 ENDDO
                 DO k=lo_a(3),hi_a(3)
                    DO j=lo_a(2),0
                       DO i=lo_a(1),hi_a(1)
                          wpr(i,j,k)=zero
                       ENDDO
                    ENDDO
                 ENDDO
              ENDIF
              IF (sbpitr%istart(3).EQ.1) THEN
                 DO k=lo_a(3),0
                    DO j=lo_a(2),hi_a(2)
                       DO i=lo_a(1),hi_a(1)
                          wpil(i,j,k)=FORBIDDEN
                       ENDDO
                    ENDDO
                 ENDDO
                 DO k=lo_a(3),0
                    DO j=lo_a(2),hi_a(2)
                       DO i=lo_a(1),hi_a(1)
                          wpr(i,j,k)=zero
                       ENDDO
                    ENDDO
                 ENDDO
              ENDIF
              IF (sbpitr%iend(1).EQ.Ngrid(1)) THEN
                 DO k=lo_a(3),hi_a(3)
                    DO j=lo_a(2),hi_a(2)
                       DO i=Nm(1)+1,hi_a(1)
                          wpil(i,j,k)=FORBIDDEN
                       ENDDO
                    ENDDO
                 ENDDO
                 DO k=lo_a(3),hi_a(3)
                    DO j=lo_a(2),hi_a(2)
                       DO i=Nm(1)+1,hi_a(1)
                          wpr(i,j,k)=zero
                       ENDDO
                    ENDDO
                 ENDDO
              ENDIF
              IF (sbpitr%iend(2).EQ.Ngrid(2)) THEN
                 DO k=lo_a(3),hi_a(3)
                    DO j=Nm(2)+1,hi_a(2)
                       DO i=lo_a(1),hi_a(1)
                          wpil(i,j,k)=FORBIDDEN
                       ENDDO
                    ENDDO
                 ENDDO
                 DO k=lo_a(3),hi_a(3)
                    DO j=Nm(2)+1,hi_a(2)
                       DO i=lo_a(1),hi_a(1)
                          wpr(i,j,k)=zero
                       ENDDO
                    ENDDO
                 ENDDO
              ENDIF
              IF (sbpitr%iend(3).EQ.Ngrid(3)) THEN
                 DO k=Nm(3)+1,hi_a(3)
                    DO j=lo_a(2),hi_a(2)
                       DO i=lo_a(1),hi_a(1)
                          wpil(i,j,k)=FORBIDDEN
                       ENDDO
                    ENDDO
                 ENDDO
                 DO k=Nm(3)+1,hi_a(3)
                    DO j=lo_a(2),hi_a(2)
                       DO i=lo_a(1),hi_a(1)
                          wpr(i,j,k)=zero
                       ENDDO
                    ENDDO
                 ENDDO
              ENDIF
#endif

              sbpitr => mesh%subpatch%next()
           ENDDO !WHILE(ASSOCIATED(sbpitr))

           NULLIFY(wpil,wpr)

           !-------------------------------------------------------------------------
!---End----! complete ghost update by popping received data
           !-------------------------------------------------------------------------
           CALL mesh%map_isend(info,sendrecv=.FALSE.)
           or_fail("mesh%map_send")
           CALL image%map_ghost_pop(mesh,info)
           or_fail("image%map_ghost_pop")
#endif

        END SELECT

        !TODO
        ALLOCATE(ppm_t_field::pind,STAT=info)
        or_fail_alloc('Failed to allocate pind.')

        CALL pind%create(1,info,dtype=ppm_type_int,name="particle_index")
        or_fail("Create pind field failed!" )

        CALL pind%discretize_on(mesh,info)
        or_fail("pind discretize_on failed!")

#ifdef __Linux
        IF (debug.GT.1) THEN
           CALL ppm_rc_mem_usage(memory,info)
           stdout("mem_usage at pind creation=",memory)
        ENDIF
#endif

        !---------------------------------------------------------------------
        !  Initialize the particles on the image
        !---------------------------------------------------------------------
        ipatch=mesh%subpatch%nb
        IF (ipatch.GT.0) THEN
           ALLOCATE(Candidates(ipatch),InnerContourContainer(ipatch), &
           &        CompetingRegions(ipatch),m_Seeds(ipatch),ppm_rc_seeds_to_remove(ipatch),STAT=info)
           or_fail_alloc("Candidates, InnerContourContainer, CompetingRegions & m_Seeds")
        ELSE
           fail("There is no patch defined for the mesh.")
        ENDIF

        NULLIFY(wpip,wpil)

        sbpitr => mesh%subpatch%begin()
        ipart=0
        ipatch=1
        DO WHILE (ASSOCIATED(sbpitr))
           Nm   => sbpitr%nnodes

           CALL sbpitr%get_field(labels,wpil,info)
           or_fail("Failed to get field i_wp data.")

           CALL sbpitr%get_field(pind,wpip,info)
           or_fail("Failed to get field wpip data.")

           wpip=0

#if   __DIME == __2D
           DO j=1,Nm(2)
              DO i=1,Nm(1)
                 IF (wpil(i,j).LT.0) THEN
                    vLabel=-wpil(i,j)
                    dd(1)=i
                    dd(2)=j

                    IsEnclosedByLabel=.TRUE.
                    DO l=1,FG_ConnectivityType%NumberOfNeighbors
                       ld=dd+FG_ConnectivityType%NeighborsPoints(:,l)
                       IF (ABS(wpil(ld(1),ld(2))).NE.vLabel) THEN
                          IsEnclosedByLabel=.FALSE.
                          EXIT
                       ENDIF
                    ENDDO !j=1,FG_ConnectivityType%NumberOfNeighbors

                    IF (IsEnclosedByLabel) THEN
                       wpil(i,j)=vLabel
                       CYCLE
                    ENDIF

                    ipart=ipart+1

                    wpip(i,j)=ipart

                    ALLOCATE(seed,STAT=info)
                    or_fail_alloc("seed")
                    CALL seed%add(dd)
                    CALL InnerContourContainer(ipatch)%push(seed,info)
                    or_fail("could not add new seed to the collection")
                 ENDIF
              ENDDO !i=1,Nm(1)
           ENDDO !j=1,Nm(2)
#elif __DIME == __3D
           DO k=1,Nm(3)
              DO j=1,Nm(2)
                 DO i=1,Nm(1)
                    IF (wpil(i,j,k).LT.0) THEN
                       vLabel=-wpil(i,j,k)

                       dd(1)=i
                       dd(2)=j
                       dd(3)=k

                       IsEnclosedByLabel=.TRUE.
                       DO l=1,FG_ConnectivityType%NumberOfNeighbors
                          ld=dd+FG_ConnectivityType%NeighborsPoints(:,l)
                          IF (ABS(wpil(ld(1),ld(2),ld(3))).NE.vLabel) THEN
                             IsEnclosedByLabel=.FALSE.
                             EXIT
                          ENDIF
                       ENDDO !j=1,FG_ConnectivityType%NumberOfNeighbors

                       IF (IsEnclosedByLabel) THEN
                          wpil(i,j,k)=vLabel
                          CYCLE
                       ENDIF

                       ipart=ipart+1

                       wpip(i,j,k)=ipart

                       ALLOCATE(seed,STAT=info)
                       or_fail_alloc("seed")
                       CALL seed%add(dd)
                       CALL InnerContourContainer(ipatch)%push(seed,info)
                       or_fail("could not add new seed to the collection")
                    ENDIF
                 ENDDO !i=1,Nm(1)
              ENDDO !j=1,Nm(2)
           ENDDO !k=1,Nm(3)
#endif

           sbpitr => mesh%subpatch%next()
           ipatch=ipatch+1
        ENDDO !WHILE(ASSOCIATED(sbpitr))

        NULLIFY(wpil,wpip)

        Npart=ipart

        ALLOCATE(ppm_t_particles_s::Part,STAT=info)
        or_fail_alloc("Part")
        !Allocate Part

        CALL Part%create(Npart,info,"Part")
        or_fail("Part%create")
        !Create Particles

        NULLIFY(xp)
        CALL Part%get_xp(xp,info)
        or_fail("Part%get_xp")

        sbpitr => mesh%subpatch%begin()
        ipatch=1
        DO WHILE (ASSOCIATED(sbpitr))
           CALL sbpitr%get_field(pind,wpip,info)
           or_fail("Failed to get field wpip data.")

           seed => InnerContourContainer(ipatch)%begin()
           DO WHILE (ASSOCIATED(seed))
              seedn => seed%first%getValue()

#if   __DIME == __2D
              ipart=wpip(seedn(1),seedn(2))

              !TODO change to get_pos2-or 3d
              xp(1,ipart)=REAL(seedn(1)+sbpitr%istart(1)-2,MK)
              xp(2,ipart)=REAL(seedn(2)+sbpitr%istart(2)-2,MK)
#elif __DIME == __3D
              ipart=wpip(seedn(1),seedn(2),seedn(3))

              ! correct the particle position
              xp(1,ipart)=REAL(seedn(1)+sbpitr%istart(1)-2,MK)
              xp(2,ipart)=REAL(seedn(2)+sbpitr%istart(2)-2,MK)
              xp(3,ipart)=REAL(seedn(3)+sbpitr%istart(3)-2,MK)
#endif

              seed => InnerContourContainer(ipatch)%next()
           ENDDO !ASSOCIATED(seed)

           sbpitr => mesh%subpatch%next()
           ipatch=ipatch+1
        ENDDO !WHILE (ASSOCIATED(sbpitr))

        NULLIFY(wpip)

        !set_xp with no read_only will change all the flags
        CALL Part%set_xp(xp,info,read_only=.TRUE.)
        or_fail("Part%set_xp")

        !Yaser
        !TOCHECK
        !this is a hack to make it work
        Part%active_topoid            =topoid
        !The particle ghostlayer should be chosen according
        !the required ghost layer, which is dependant on the energy
        !type
        Part%ghostlayer=1.5_MK    !MAXVAL(ghostsize)-0.5_MK
        !Particles are inside the computational domain
        Part%flags(ppm_part_areinside)=.TRUE.
        !Particles are on the right processors
        Part%flags(ppm_part_partial)  =.TRUE.
        ! Dangereous to use the ghosts
        Part%flags(ppm_part_ghosts)   =.FALSE.

        CALL Part%comp_global_index(info)
        or_fail("Part%comp_global_index")

#ifdef __Linux
        IF (debug.GT.1) THEN
           CALL ppm_rc_mem_usage(memory,info)
           stdout("mem_usage at computing particle global index property=",memory)
        ENDIF
#endif

        CALL Part%map_ghosts(info)
        or_fail("Part%map_ghosts")
        !I need more particles in case of Energy functional which
        !needs a bigger ghost of particles

        !no change on flags
        CALL Part%get_xp(xp,info,with_ghosts=.TRUE.)
        or_fail("Part%get_xp")

        sbpitr => mesh%subpatch%begin()
        ipatch=1
        DO WHILE (ASSOCIATED(sbpitr))
           Nm => sbpitr%nnodes

           CALL sbpitr%get_field(pind,wpip,info)
           or_fail("Failed to get field wpip data.")

           DO ipart=Part%Npart+1,Part%Mpart
              ld(1)=NINT(xp(1,ipart))+2-sbpitr%istart(1)
              ld(2)=NINT(xp(2,ipart))+2-sbpitr%istart(2)

#if   __DIME == __2D
              wpip(ld(1),ld(2))=ipart
#elif __DIME == __3D
              ld(3)=NINT(xp(3,ipart))+2-sbpitr%istart(3)
              wpip(ld(1),ld(2),ld(3))=ipart
#endif

              IF (ANY(ld.LT.0.OR.ld.GT.Nm+1)) CYCLE
              !We are only interested in the ghost particles which are
              !adjacent to the domain

              ALLOCATE(seed,STAT=info)
              or_fail_alloc("seed")
              CALL seed%add(ld)
              CALL InnerContourContainer(ipatch)%push(seed,info)
              or_fail("could not add new seed to the collection")
              !Add one layer of ghost particles to the part contour
              !We would need them for the energy calculation
           ENDDO

           sbpitr => mesh%subpatch%next()
           ipatch=ipatch+1
        ENDDO !WHILE (ASSOCIATED(sbpitr))

        NULLIFY(wpip)

        !set_xp with no read_only will change all the flags
        CALL Part%set_xp(xp,info,read_only=.TRUE.)
        or_fail("Part%set_xp")

        Mpart=Part%Mpart

        ALLOCATE(ppm_t_field::plabels,STAT=info)
        or_fail_alloc('Failed to allocate plabels.')

        IF (freqoutput.LT.maxiter) THEN
           CALL DTYPE(ppm_rc_write_image_label)(labels,mesh,outputfile,info,idn=0)
           or_fail("ppm_rc_write_image_label")
           !CALL ppm_vtk_particles("yaser",Part,info,with_ghosts=.TRUE.)
        ENDIF

#ifdef  __MPI
        sbpitr => mesh%subpatch%begin()
        nsize=0
        DO WHILE (ASSOCIATED(sbpitr))
           !The size of ghost layer around each subpatch
           !we compute only one layer of ghost around each patch
           Nm => sbpitr%nnodes

           IF (ALL(sbpitr%istart(1:ppm_rc_dim).EQ.1.AND. &
           &       sbpitr%iend(1:ppm_rc_dim).EQ.Ngrid(1:__DIME))) THEN
              !do nothing
           ELSE
#if   __DIME == __2D
              IF (sbpitr%istart(1).NE.1) THEN
                 nsize=nsize+Nm(2)
              ENDIF
              IF (sbpitr%iend(1).NE.Ngrid(1)) THEN
                 nsize=nsize+Nm(2)
              ENDIF
              IF (sbpitr%istart(2).NE.1) THEN
                 nsize=nsize+Nm(1)+2
              ENDIF
              IF (sbpitr%iend(2).NE.Ngrid(2)) THEN
                 nsize=nsize+Nm(1)+2
              ENDIF
#elif __DIME == __3D
              IF (sbpitr%istart(1).NE.1) THEN
                 nsize=nsize+Nm(2)*Nm(3)
              ENDIF
              IF (sbpitr%iend(1).NE.Ngrid(1)) THEN
                 nsize=nsize+Nm(2)*Nm(3)
              ENDIF
              IF (sbpitr%istart(2).NE.1) THEN
                 nsize=nsize+Nm(3)*(Nm(1)+2)
              ENDIF
              IF (sbpitr%iend(2).NE.Ngrid(2)) THEN
                 nsize=nsize+Nm(3)*(Nm(1)+2)
              ENDIF
              IF (sbpitr%istart(3).NE.1) THEN
                 nsize=nsize+(Nm(1)+2)*(Nm(2)+2)
              ENDIF
              IF (sbpitr%iend(3).NE.Ngrid(3)) THEN
                 nsize=nsize+(Nm(1)+2)*(Nm(2)+2)
              ENDIF
#endif
           ENDIF

           sbpitr => mesh%subpatch%next()
        ENDDO !WHILE (ASSOCIATED(sbpitr))

        ghost_size=nsize
#endif

        ALLOCATE(list_del_parts(Npart),STAT=info)
        or_fail_alloc("list_del_parts")

        del_parts=0

#ifdef  __MPI
        topo => ppm_topo(topoid)%t

        IF (ppm_nproc.GT.1) THEN
           max_sub => topo%max_subs
           min_sub => topo%min_subs

           !wait for ineigh to be avilable on all ranks
           CALL MPI_Wait(request,MPI_STATUS_IGNORE,info)
           or_fail_MPI("MPI_Wait")

           ALLOCATE(flag(topo%nsubs),STAT=info)
           or_fail_alloc("flag")

           flag=0

           DO i=1,topo%nsubs
              j=topo%sub2proc(i)+1

              IF (flag(i).EQ.0) flag(i)=1

              DO k=displ(j)+1,displ(j)+nneigh(j)
                 l=ineigh(k)
                 IF (flag(l).EQ.0) then
                    sx=MAX(min_sub(1,i),min_sub(1,l))
                    ex=MIN(max_sub(1,i),max_sub(1,l))
                    lx=ABS(ex-sx)
                    sy=MAX(min_sub(2,i),min_sub(2,l))
                    ey=MIN(max_sub(2,i),max_sub(2,l))
                    ly=ABS(ey-sy)
#if   __DIME == __2D
                    IF (lx.LE.lmyeps.AND.ly.LE.lmyeps) THEN
                       flag(l)=flag(i)
                    ELSE IF (lx.LE.lmyeps.OR.ly.LE.lmyeps) THEN
                       flag(l)=MERGE(1,2,flag(i).EQ.2)
                    ELSE
                       fail("St is wrong!!!",ppm_error=ppm_error_fatal)
                    ENDIF
#elif __DIME == __3D
                    sz=MAX(min_sub(3,i),min_sub(3,l))
                    ez=MIN(max_sub(3,i),max_sub(3,l))
                    lz=ABS(ez-sz)
                    IF (lx.LE.lmyeps.AND.ly.LE.lmyeps.AND.lz.LE.lmyeps) THEN
                       flag(l)=MERGE(1,2,flag(i).EQ.2)
                    ELSE IF (lx.LE.lmyeps.AND.ly.LE.lmyeps) THEN
                       flag(l)=flag(i)
                    ELSE IF (lx.LE.lmyeps.AND.lz.LE.lmyeps) THEN
                       flag(l)=flag(i)
                    ELSE IF (ly.LE.lmyeps.AND.lz.LE.lmyeps) THEN
                       flag(l)=flag(i)
                    ELSE IF (lx.LE.lmyeps) THEN
                       flag(l)=MERGE(1,2,flag(i).EQ.2)
                    ELSE IF (ly.LE.lmyeps) THEN
                       flag(l)=MERGE(1,2,flag(i).EQ.2)
                    ELSE IF (lz.LE.lmyeps) THEN
                       flag(l)=MERGE(1,2,flag(i).EQ.2)
                    ELSE
                       fail("St is wrong!!!",ppm_error=ppm_error_fatal)
                    ENDIF
#endif
                 ENDIF !(flag(l).EQ.0)
              ENDDO !k=1,nneigh(p)
           ENDDO !i=1,topo%nsubs

           DEALLOCATE(nneigh,ineigh,displ,STAT=info)
           or_fail_dealloc("nneigh,ineigh & displ")
        ENDIF !(ppm_nproc.GT.1)

        CALL MPI_Ibarrier(comm,request,info)
        or_fail_MPI("MPI_Ibarrier")

        nneighproc = topo%nneighproc

        ALLOCATE(ineighproc(1:nneighproc),SOURCE=topo%ineighproc(1:nneighproc),STAT=info)
        or_fail_alloc("ineighproc")

        !convert rank to index
        ineighproc=ineighproc+1

        ALLOCATE(procflag(nneighproc+1),STAT=info)
        or_fail_alloc("procflag")

        IF (ppm_nproc.GT.1) THEN
           DO i=1,topo%nsubs
              j=topo%sub2proc(i)+1
              IF (rank+1.EQ.j) THEN
                 procflag(1)=flag(i)
              ELSE
                 DO k=1,nneighproc
                    IF (ineighproc(k).EQ.j) THEN
                       procflag(k+1)=flag(i)
                    ENDIF
                 ENDDO
              ENDIF
           ENDDO

           DEALLOCATE(flag,STAT=info)
           or_fail_dealloc("flag")
        ELSE
           procflag(1)=1
        ENDIF
#endif

#ifdef __Linux
        IF (debug.GT.1) THEN
           CALL ppm_rc_mem_usage(memory,info)
           stdout("mem_usage at program Initialization complete=",memory)
        ENDIF
#endif

#ifdef  __MPI
        !wait for everyone to pass the barrier
        CALL MPI_Wait(request,MPI_STATUS_IGNORE,info)
        or_fail_MPI("MPI_Wait")
#endif
        !-------------------------------------------------------------------------
        !  End of init
        !-------------------------------------------------------------------------
        IF (rank.EQ.0) THEN
           stdout("RC initialization is complete.")
        ENDIF

        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
      9999 CONTINUE
        CALL substop(caller,t0,info)
        RETURN
      END SUBROUTINE DTYPE(ppm_rc_init_seg)
