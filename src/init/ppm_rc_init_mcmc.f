      !-------------------------------------------------------------------------
      !  Subroutine   :                    ppm_rc_init_mcmc
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
      !  Author           - y.afshar           April   2015
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Subroutine   :                    ppm_rc_init_mcmc
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
      SUBROUTINE DTYPE(ppm_rc_init_mcmc)(info)

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
        USE ppm_rc_module_energy, ONLY : e_data
        USE ppm_rc_module_topologicalnumber, ONLY : FG_ConnectivityType
        USE ppm_rc_module_hash, ONLY : ppm_rc_HashIndextable
        USE ppm_rc_module_mcmc, ONLY : CreateMCMClengthProposalMask,MCMCcellNm,   &
        &   MCMCFreeMemory,MCMCInsertCandidatesToContainers,MCMCcellsize,         &
        &   MCMCboundarycellIndex,MCMCboundarycellSize,MCMCRegularParticles,      &
        &   MCMCRegularParticlesInCell,MCMCFloatingParticlesInCell,               &
        &   MCMCRegularParticlesProposalNormalizer,MCMCInsertLabelInRegionLabel,  &
        &   MCMCRegularParticlesProposalNormalizerlocal,MCMCIsParticleTopoValid,  &
        &   MCMCgetRegularParticlesAtIndex,MCMCinteriorcellIndex,                 &
        &   MCMCinteriorcellDisp
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

        !-------------------------------------------------------------------------
        !  Temporary pointers to the data and bookkeeping information
        !  for a mesh on which fieldID has been discretized.
        !-------------------------------------------------------------------------
        CLASS(ppm_t_discr_info_),  POINTER :: imgdinfo
        CLASS(ppm_t_discr_info_),  POINTER :: lbldinfo
        CLASS(ppm_t_discr_info_),  POINTER :: initdinfo

        TYPE(ppm_t_topo), POINTER :: topo

        TYPE(ppm_rc_HashIndextable) :: htablei

#if   __DIME == __2D
        TYPE(MCMCParticle), DIMENSION( 9) :: MCMCParticles
#elif __DIME == __3D
        TYPE(MCMCParticle), DIMENSION(27) :: MCMCParticles
#endif

        REAL(ppm_kind_double), DIMENSION(__DIME)          :: Offset
        REAL(ppm_kind_double)                             :: t0
        REAL(MK)                                          :: FieldMin,FieldMax
#if   __DIME == __2D
        REAL(MK), CONTIGUOUS,  DIMENSION(:,:),    POINTER :: wpi
#elif __DIME == __3D
        REAL(MK), CONTIGUOUS,  DIMENSION(:,:,:),  POINTER :: wpi
#endif
        REAL(MK), CONTIGUOUS,  DIMENSION(:,:),    POINTER :: xp

#if   __DIME == __2D
        INTEGER, CONTIGUOUS,  DIMENSION(:,:),  POINTER :: wpl
        INTEGER,              DIMENSION(:,:),  POINTER :: tmplabels_
        INTEGER,              DIMENSION(:,:),  POINTER :: tmplabels
#elif __DIME == __3D
        INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER :: wpl
        INTEGER,             DIMENSION(:,:,:), POINTER :: tmplabels_
        INTEGER,             DIMENSION(:,:,:), POINTER :: tmplabels
#endif
        INTEGER,             DIMENSION(:),     POINTER :: Nm
        INTEGER,             DIMENSION(:),     POINTER :: hi_a
        INTEGER,             DIMENSION(:),     POINTER :: lo_a
        INTEGER,             DIMENSION(__DIME)              :: ldu,ld
        INTEGER                                        :: nproc
        INTEGER                                        :: vLabel
        INTEGER                                        :: nLabel
        INTEGER                                        :: i,j
#if   __DIME == __2D
        INTEGER,             DIMENSION(4)              :: ll
#endif
#if   __DIME == __3D
        INTEGER                                        :: k
        INTEGER,             DIMENSION(8)              :: ll
#endif
        INTEGER                                        :: l,iopt,m
        INTEGER                                        :: nsize,iseed
        INTEGER                                        :: nx,ny,cbox
        INTEGER                                        :: lieven,ljeven
#if   __DIME == __3D
        INTEGER                                        :: nz
        INTEGER                                        :: lkeven
#endif
#ifdef __MPI
        INTEGER                                        :: request
#endif

        CHARACTER(LEN=ppm_char) :: caller = 'ppm_rc_init_mcmc'

        LOGICAL :: Replaced
        !-------------------------------------------------------------------------
        !  Externals
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
        !  Initialize
        !-------------------------------------------------------------------------
        CALL substart(caller,t0,info)

        SELECT CASE (MCMCcontinue)
        ! In this case we start a MCMC sampling from input image and
        ! a label files (either as an input with different labels or create one)
        CASE (.FALSE.)
           ! debug=2
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
              ghostsize=MAX(ioghostsize,inighostsize,ghostsize_run)
           CASE DEFAULT
              !-------------------------------------------------------------------------
              !  define the ghostsize for IO topology
              !-------------------------------------------------------------------------
              ghostsize=ioghostsize

              !-------------------------------------------------------------------------
              !  If the user provides the number of IO processors
              !  use it as a base.
              !-------------------------------------------------------------------------
              IF (n_procs_read.GT.1.AND.n_procs_read.LT.ppm_nproc) THEN
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

           !---------------------------------------------------------------------
           ! Free memory.
           !---------------------------------------------------------------------
           DEALLOCATE(ioghostsize,ghostsize_equil,STAT=info)
           or_fail_alloc("ioghostsize & ghostsize_equil")

           !---------------------------------------------------------------------
           ! Create new topology suitable for reading and writing images.
           !---------------------------------------------------------------------
           assig  = ppm_param_assign_internal
           decomp = MERGE(ppm_param_decomp_xpencil,ppm_param_decomp_xy_slab,ppm_rc_dim.EQ.2)

           !---------------------------------------------------------------------
           ! Create dummy particle.
           !---------------------------------------------------------------------
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

           !---------------------------------------------------------------------
           ! Get rid of the dummy particle
           !---------------------------------------------------------------------
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
              or_fail("Failed to read image files and write them into one image file.")

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

           !---------------------------------------------------------------------
           ! Read the input file
           ! Store the intensity on image field which is defined & descretized
           ! on iomesh as a real variable
           !---------------------------------------------------------------------
           CALL DTYPE(ppm_rc_read_image)(image,iomesh,inputimage,info)
           or_fail('Failed to read image.')

           !---------------------------------------------------------------------
           ! In case we have another input file for the initialization
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
                 WRITE(cbuf,'(A,A)') "The current version of DRS client is only ", &
                 & "working with the number of processors equlas to any power of 2"
                 fail(cbuf,ppm_error=ppm_error_fatal)
              ENDIF

              !!! TODO
              !!! In the new update I should get rid of this
              !!! Currently I disbaled initopo map to topo
              ghostsize=MAX(inighostsize,ghostsize_run)

              assig  = ppm_param_assign_internal
              decomp = ppm_param_decomp_cartesian
              !ppm_param_decomp_bisection

              !---------------------------------------------------------------------
              ! Create dummy particle
              !---------------------------------------------------------------------
              NULLIFY(xp)
              iopt=ppm_param_alloc_fit
              ldu(1)=__DIME
              ldu(2)=1
              CALL ppm_alloc(xp,ldu,iopt,info)
              or_fail_alloc("xp")
              xp=zero

              !-------------------------------------------------------------------------
              ! Create Init topology
              !-------------------------------------------------------------------------
              CALL ppm_mktopo(initopoid,initmeshid,xp,0,decomp,assig,min_phys, &
              &    max_phys,bcdef,ghostsize,cost,Ngrid(1:__DIME),info)
              or_fail('Failed to create new topology.')

              !-------------------------------------------------------------------------
              ! Init mesh
              !-------------------------------------------------------------------------
              imesh => ppm_mesh%at(initmeshid)

              !-------------------------------------------------------------------------
              ! Create uniform initlization mesh over the domain
              !-------------------------------------------------------------------------
              CALL imesh%def_uniform(info)
              or_fail("Failed to create uniform mesh.")

              !---------------------------------------------------------------------
              ! Descretize the scalar field (image) on Init mesh
              !---------------------------------------------------------------------
              CALL image%discretize_on(imesh,info,discr_info=initdinfo)
              or_fail("image discretize_on failed!")

              !---------------------------------------------------------------------
              ! Global mapping to the new mesh in the new topology
              ! the mapping is between IO and Init topology
              !---------------------------------------------------------------------
              CALL iomesh%map(imesh,info)
              or_fail("Failed to do global mapping from old to the new topology.")
              !---------------------------------------------------------------------
              ! Push image data in the buffer
              !---------------------------------------------------------------------
              CALL image%map_push(iomesh,info)
              or_fail("Failed to push image data in the buffer!")
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

              !---------------------------------------------------------------------
              ! We need to free memory which will not be used later
              ! image data on IO mesh will not be needed anymore
              !---------------------------------------------------------------------
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

                 !---------------------------------------------------------------------
                 ! Discretize a label field over the new initialization mesh
                 ! Take advantage of computation over communication
                 !---------------------------------------------------------------------
                 CALL labels%discretize_on(imesh,info)
                 or_fail("labels discretize_on failed!")

                 !---------------------------------------------------------------------
                 ! Get rid of the object which will not be used anymore to free memory
                 !---------------------------------------------------------------------
                 CALL labels%discr_info%remove(info,lbldinfo)
                 or_fail("labels%discr_info%remove")

                 DEALLOCATE(lbldinfo,STAT=info)
                 or_fail_dealloc("Failed to deallocate discr_info")

                 !---------------------------------------------------------------------
                 ! We need to free memory which will not be used later
                 ! labels data on IO mesh will not be needed anymore
                 !---------------------------------------------------------------------
                 NULLIFY(lbldinfo)

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

                 !---------------------------------------------------------------------
                 ! Destroy the IO topo
                 !---------------------------------------------------------------------
                 topo => ppm_topo(iotopoid)%t

                 iopt=ppm_param_dealloc
                 CALL ppm_alloc(topo,iopt,info)
                 or_fail_dealloc("Failed to deallocate topology (iotopoid)!")
              ELSE
                 !If we do not have labels already descretized on IO mesh, we can first destroy
                 !IO mesh and free memory, then descretize labels on new topo and mesh
                 !---------------------------------------------------------------------
                 ! Destroy the IO mesh
                 !---------------------------------------------------------------------
                 CALL iomesh%destroy(info)
                 or_fail('Failed to destroy iomesh.')

                 dealloc_pointer("iomesh")

                 !---------------------------------------------------------------------
                 ! Destroy the IO topo
                 !---------------------------------------------------------------------
                 topo => ppm_topo(iotopoid)%t

                 iopt=ppm_param_dealloc
                 CALL ppm_alloc(topo,iopt,info)
                 or_fail_dealloc("Failed to deallocate topology (iotopoid)!")

                 !---------------------------------------------------------------------
                 ! Define labels field variable to assign label
                 ! the image pixels
                 !---------------------------------------------------------------------
                 ALLOCATE(ppm_t_field::labels,STAT=info)
                 or_fail_alloc('Failed to allocate labels field.')

                 CALL labels%create(1,info,dtype=ppm_type_int,name="pixel_labels")
                 or_fail("Create labels field failed!" )

                 !---------------------------------------------------------------------
                 ! Discretize a label field over the new initialization mesh
                 !---------------------------------------------------------------------
                 CALL labels%discretize_on(imesh,info)
                 or_fail("labels discretize_on failed!")
              ENDIF

           END SELECT

           Replaced=.FALSE.
           !-------------------------------------------------------------------------
           ! Normalize and Shift
           !-------------------------------------------------------------------------
           SELECT CASE (vInitKind)
           CASE (e_localmax)
              lNormalize=.TRUE.

              !-------------------------------------------------------------------------
              ! Normalize and Shift the input image to [0,1]
              !-------------------------------------------------------------------------
              CALL DTYPE(ppm_rc_CopyImageAndNormalize)(image,image,imesh,info, &
              &    Normalize=lNormalize,FieldMinVal=FieldMin,FieldMaxVal=FieldMax)

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
                 Replaced=.TRUE.
                 intensity_ths=intensity_ths/ImageNormalfac
              ENDIF
           CASE DEFAULT
              !-------------------------------------------------------------------------
              ! Do not use normalized imaeg for MCMC
              !-------------------------------------------------------------------------
              lNormalize=.FALSE.

              Tolerance_=REAL(Tolerance,MK)
           END SELECT

           !-------------------------------------------------------------------------
           !  Initialize init related terms
           !-------------------------------------------------------------------------
           ALLOCATE(e_Init,STAT=info)
           or_fail_alloc("e_Init")

           ldu=FLOOR(init_rd)
           ld=FLOOR(init_sp)
           !-------------------------------------------------------------------------
           !  Create initialization elements
           !-------------------------------------------------------------------------
           CALL e_Init%CreateElement(ppm_rc_dim,info,ForegroundValue=1, &
           &    BackgroundValue=0,size_=ldu,spacing_=ld)
           or_fail("e_Init%CreateElement")

           !-------------------------------------------------------------------------
           !  Create initialization regions
           !-------------------------------------------------------------------------
           CALL e_Init%DTYPE(GetOutput)(image,imesh,labels,info)
           or_fail("e_Init%GetOutput")

           !-------------------------------------------------------------------------
           !  Updating Ghost information on every processor
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

           !-------------------------------------------------------------------------
           ! Destroy the initialization element and assigned memory
           !-------------------------------------------------------------------------
           CALL e_Init%DestroyElement(info)
           or_fail("e_Init%DestroyElement")

           DEALLOCATE(e_Init,inighostsize,STAT=info)
           or_fail_dealloc("e_Init & inighostsize")

           NULLIFY(e_Init)

           IF (debug.GT.1) THEN
              CALL DTYPE(ppm_rc_write_image_label)(labels,imesh,"afterput",info,idn=FORBIDDEN)
              or_fail("ppm_rc_write_image_label")
           ENDIF

           !-------------------------------------------------------------------------
           ! Normalize and Shift
           !-------------------------------------------------------------------------
           SELECT CASE (vInitKind)
           CASE (e_localmax)
              lNormalize=.FALSE.

              !-------------------------------------------------------------------------
              ! Scale the Normalized input image to its original
              !-------------------------------------------------------------------------
              CALL DTYPE(ppm_rc_CopyImageAndNormalize)(image,image,imesh,info, &
              &    Normalize=lNormalize,Scalefac=ImageNormalfac)

              !-------------------------------------------------------------------------
              ! Correct the tolerance, when we are normalizing the image
              !-------------------------------------------------------------------------
              Tolerance_=REAL(Tolerance,MK)

              IF (Replaced) THEN
                 intensity_ths=intensity_ths*ImageNormalfac
              ENDIF

              ImageNormalfac=one
           END SELECT

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

           !-------------------------------------------------------------------------
           ! Count the seed points and allocate newlabels
           !-------------------------------------------------------------------------
           nsize=imesh%subpatch%nb
           ld(1)=SUM(ppm_rc_seeds(1:nsize)%nb)

           ALLOCATE(nlabels(ld(1)),STAT=info)
           or_fail_alloc("nlabels")

           !-------------------------------------------------------------------------
           ! Assign the new label for each seed region
           ! which corresponds to one intial region
           !-------------------------------------------------------------------------
           DO iseed=1,ld(1)
              nlabels(iseed)=loc_label
              loc_label=loc_label-1
           ENDDO

           sbpitr => imesh%subpatch%begin()
           nsize=0
           DO WHILE (ASSOCIATED(sbpitr))
              !-------------------------------------------------------------------------
              ! The size of ghost layer around each subpatch
              ! is only one row of ghost
              !-------------------------------------------------------------------------
              Nm => sbpitr%nnodes

              IF (ALL(sbpitr%istart(1:ppm_rc_dim).EQ.1.AND. &
              &       sbpitr%iend(1:ppm_rc_dim).EQ.Ngrid(1:__DIME))) THEN
                 ! do nothing
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

              sbpitr => imesh%subpatch%next()
           ENDDO

           !-------------------------------------------------------------------------
           !  Start fire from seed points and label regions
           !-------------------------------------------------------------------------
           CALL DTYPE(ppm_rc_initforestfire)(imesh,nsize,info)
           or_fail("ppm_rc_forestfire")

           DEALLOCATE(nlabels,STAT=info)
           or_fail_dealloc("nlabels")

           !-------------------------------------------------------------------------
           !  Now create new topology and map to the new topology
           !-------------------------------------------------------------------------
           ghostsize=ghostsize_run

           SELECT CASE (ppm_nproc)
           CASE (1)
              topoid = initopoid

              mesh => iomesh

              meshid = iomeshid

              NULLIFY(imesh)

              NULLIFY(wpl,wpi)

              sbpitr => mesh%subpatch%begin()
              DO WHILE (ASSOCIATED(sbpitr))
                 Nm   => sbpitr%nnodes
                 hi_a => sbpitr%hi_a
                 lo_a => sbpitr%lo_a

                 CALL sbpitr%get_field(labels,wpl,info)
                 or_fail("Failed to get field wpl data.")

                 CALL sbpitr%get_field(image,wpi,info)
                 or_fail("Failed to get field wpi data.")

#if   __DIME == __2D
                 IF (sbpitr%istart(1).EQ.1) THEN
                    DO j=lo_a(2),hi_a(2)
                       DO i=lo_a(1),1
                          wpl(i,j)=FORBIDDEN
                       ENDDO
                    ENDDO
                    DO j=lo_a(2),hi_a(2)
                       DO i=lo_a(1),0
                          wpi(i,j)=zero
                       ENDDO
                    ENDDO
                 ENDIF
                 IF (sbpitr%istart(2).EQ.1) THEN
                    DO j=lo_a(2),1
                       DO i=lo_a(1),hi_a(1)
                          wpl(i,j)=FORBIDDEN
                       ENDDO
                    ENDDO
                    DO j=lo_a(2),0
                       DO i=lo_a(1),hi_a(1)
                          wpi(i,j)=zero
                       ENDDO
                    ENDDO
                 ENDIF
                 IF (sbpitr%iend(1).EQ.Ngrid(1)) THEN
                    DO j=lo_a(2),hi_a(2)
                       DO i=Nm(1),hi_a(1)
                          wpl(i,j)=FORBIDDEN
                       ENDDO
                    ENDDO
                    DO j=lo_a(2),hi_a(2)
                       DO i=Nm(1)+1,hi_a(1)
                          wpi(i,j)=zero
                       ENDDO
                    ENDDO
                 ENDIF
                 IF (sbpitr%iend(2).EQ.Ngrid(2)) THEN
                    DO j=Nm(2),hi_a(2)
                       DO i=lo_a(1),hi_a(1)
                          wpl(i,j)=FORBIDDEN
                       ENDDO
                    ENDDO
                    DO j=Nm(2)+1,hi_a(2)
                       DO i=lo_a(1),hi_a(1)
                          wpi(i,j)=zero
                       ENDDO
                    ENDDO
                 ENDIF
#elif __DIME == __3D
                 IF (sbpitr%istart(1).EQ.1) THEN
                    DO k=lo_a(3),hi_a(3)
                       DO j=lo_a(2),hi_a(2)
                          DO i=lo_a(1),1
                             wpl(i,j,k)=FORBIDDEN
                          ENDDO
                       ENDDO
                    ENDDO
                    DO k=lo_a(3),hi_a(3)
                       DO j=lo_a(2),hi_a(2)
                          DO i=lo_a(1),0
                             wpi(i,j,k)=zero
                          ENDDO
                       ENDDO
                    ENDDO
                 ENDIF
                 IF (sbpitr%istart(2).EQ.1) THEN
                    DO k=lo_a(3),hi_a(3)
                       DO j=lo_a(2),1
                          DO i=lo_a(1),hi_a(1)
                             wpl(i,j,k)=FORBIDDEN
                          ENDDO
                       ENDDO
                    ENDDO
                    DO k=lo_a(3),hi_a(3)
                       DO j=lo_a(2),0
                          DO i=lo_a(1),hi_a(1)
                             wpi(i,j,k)=zero
                          ENDDO
                       ENDDO
                    ENDDO
                 ENDIF
                 IF (sbpitr%istart(3).EQ.1) THEN
                    DO k=lo_a(3),1
                       DO j=lo_a(2),hi_a(2)
                          DO i=lo_a(1),hi_a(1)
                             wpl(i,j,k)=FORBIDDEN
                          ENDDO
                       ENDDO
                    ENDDO
                    DO k=lo_a(3),0
                       DO j=lo_a(2),hi_a(2)
                          DO i=lo_a(1),hi_a(1)
                             wpi(i,j,k)=zero
                          ENDDO
                       ENDDO
                    ENDDO
                 ENDIF
                 IF (sbpitr%iend(1).EQ.Ngrid(1)) THEN
                    DO k=lo_a(3),hi_a(3)
                       DO j=lo_a(2),hi_a(2)
                          DO i=Nm(1),hi_a(1)
                             wpl(i,j,k)=FORBIDDEN
                          ENDDO
                       ENDDO
                    ENDDO
                    DO k=lo_a(3),hi_a(3)
                       DO j=lo_a(2),hi_a(2)
                          DO i=Nm(1)+1,hi_a(1)
                             wpi(i,j,k)=zero
                          ENDDO
                       ENDDO
                    ENDDO
                 ENDIF
                 IF (sbpitr%iend(2).EQ.Ngrid(2)) THEN
                    DO k=lo_a(3),hi_a(3)
                       DO j=Nm(2),hi_a(2)
                          DO i=lo_a(1),hi_a(1)
                             wpl(i,j,k)=FORBIDDEN
                          ENDDO
                       ENDDO
                    ENDDO
                    DO k=lo_a(3),hi_a(3)
                       DO j=Nm(2)+1,hi_a(2)
                          DO i=lo_a(1),hi_a(1)
                             wpi(i,j,k)=zero
                          ENDDO
                       ENDDO
                    ENDDO
                 ENDIF
                 IF (sbpitr%iend(3).EQ.Ngrid(3)) THEN
                    DO k=Nm(3),hi_a(3)
                       DO j=lo_a(2),hi_a(2)
                          DO i=lo_a(1),hi_a(1)
                             wpl(i,j,k)=FORBIDDEN
                          ENDDO
                       ENDDO
                    ENDDO
                    DO k=Nm(3)+1,hi_a(3)
                       DO j=lo_a(2),hi_a(2)
                          DO i=lo_a(1),hi_a(1)
                             wpi(i,j,k)=zero
                          ENDDO
                       ENDDO
                    ENDDO
                 ENDIF
#endif

                 sbpitr => mesh%subpatch%next()
              ENDDO !WHILE(ASSOCIATED(sbpitr))

              NULLIFY(wpl,wpi)

              !-------------------------------------------------------------------------
              !  Check if all went well on image
              !-------------------------------------------------------------------------
              CALL DTYPE(ppm_rc_write_image)(image,mesh,"image",info)
              or_fail("ppm_rc_write_image")

#ifdef __MPI
           CASE DEFAULT
              NULLIFY(wpl)

              sbpitr => imesh%subpatch%begin()
              DO WHILE (ASSOCIATED(sbpitr))
                 Nm   => sbpitr%nnodes

                 CALL sbpitr%get_field(labels,wpl,info)
                 or_fail("Failed to get field wpl data.")

#if   __DIME == __2D
                 IF (sbpitr%istart(1).EQ.1) THEN
                    DO j=1,Nm(2)
                       wpl(1,j)=FORBIDDEN
                    ENDDO
                 ENDIF
                 IF (sbpitr%istart(2).EQ.1) THEN
                    DO i=1,Nm(1)
                       wpl(i,1)=FORBIDDEN
                    ENDDO
                 ENDIF
                 IF (sbpitr%iend(1).EQ.Ngrid(1)) THEN
                    DO j=1,Nm(2)
                       wpl(Nm(1),j)=FORBIDDEN
                    ENDDO
                 ENDIF
                 IF (sbpitr%iend(2).EQ.Ngrid(2)) THEN
                    DO i=1,Nm(1)
                       wpl(i,Nm(2))=FORBIDDEN
                    ENDDO
                 ENDIF
#elif __DIME == __3D
                 IF (sbpitr%istart(1).EQ.1) THEN
                    DO k=1,Nm(3)
                       DO j=1,Nm(2)
                          wpl(1,j,k)=FORBIDDEN
                       ENDDO
                    ENDDO
                 ENDIF
                 IF (sbpitr%istart(2).EQ.1) THEN
                    DO k=1,Nm(3)
                       DO i=1,Nm(1)
                          wpl(i,1,k)=FORBIDDEN
                       ENDDO
                    ENDDO
                 ENDIF
                 IF (sbpitr%istart(3).EQ.1) THEN
                    DO j=1,Nm(2)
                       DO i=1,Nm(1)
                          wpl(i,j,1)=FORBIDDEN
                       ENDDO
                    ENDDO
                 ENDIF
                 IF (sbpitr%iend(1).EQ.Ngrid(1)) THEN
                    DO k=1,Nm(3)
                       DO j=1,Nm(2)
                          wpl(Nm(1),j,k)=FORBIDDEN
                       ENDDO
                    ENDDO
                 ENDIF
                 IF (sbpitr%iend(2).EQ.Ngrid(2)) THEN
                    DO k=1,Nm(3)
                       DO i=1,Nm(1)
                          wpl(i,Nm(2),k)=FORBIDDEN
                       ENDDO
                    ENDDO
                 ENDIF
                 IF (sbpitr%iend(3).EQ.Ngrid(3)) THEN
                    DO j=1,Nm(2)
                       DO i=1,Nm(1)
                          wpl(i,j,Nm(3))=FORBIDDEN
                       ENDDO
                    ENDDO
                 ENDIF
#endif

                 sbpitr => imesh%subpatch%next()
              ENDDO !WHILE(ASSOCIATED(sbpitr))

              NULLIFY(wpl,xp)

              topoid = initopoid
              topo => ppm_topo(topoid)%t

              IF (topo%nsubs.NE.ppm_nproc) THEN
                 fail("This implementation is only works with one subdomain per processor", &
                 & ppm_error=ppm_error_fatal)
              ENDIF

              Offset=zerod

              ALLOCATE(ppm_t_equi_mesh::mesh,STAT=info)
              or_fail_alloc('Failed to allocate mesh pointer.')

              CALL mesh%create(topoid,Offset,info, &
              &    Nm=Ngrid(1:__DIME),ghostsize=ghostsize)
              or_fail("Failed to create mesh.")

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

                 !-------------------------------------------------------------------------
                 ! We need to free memory which will not be used later
                 ! image data on IO mesh will not be needed anymore
                 !-------------------------------------------------------------------------
                 NULLIFY(initdinfo)
              ENDIF

              !---------------------------------------------------------------------
              ! starting labels global mapping to the new mesh in the new topology
              !---------------------------------------------------------------------
              CALL labels%map_push(imesh,info)
              or_fail("Failed to push image data in the buffer!")
              CALL imesh%map_isend(info,sendrecv=.TRUE.)
              or_fail("Failed to send image data to the new mesh!")

              !-------------------------------------------------------------------------
              ! Check if all went well on image
              ! Writing the image based on the new decomposition
              !-------------------------------------------------------------------------
              CALL DTYPE(ppm_rc_write_image)(image,mesh,"image",info)
              or_fail("ppm_rc_write_image")

              !---------------------------------------------------------------------
              ! complete labels global mapping to the new mesh in the new topology
              !---------------------------------------------------------------------
              CALL imesh%map_isend(info,sendrecv=.FALSE.)
              or_fail("Failed to send image data to the new mesh!")
              CALL labels%map_pop(mesh,info)
              or_fail("Failed to pop image data in the new mesh!")

              !-------------------------------------------------------------------------
              ! ghost update on the new mesh at new topo
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
              !---------------------------------------------------------------------
              ! Destroy the initial imesh
              !---------------------------------------------------------------------
              CALL imesh%destroy(info)
              or_fail('Failed to destroy imesh.')

              dealloc_pointer("imesh")

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
              NULLIFY(wpl,wpi)

              sbpitr => mesh%subpatch%begin()
              DO WHILE (ASSOCIATED(sbpitr))
                 Nm   => sbpitr%nnodes
                 hi_a => sbpitr%hi_a
                 lo_a => sbpitr%lo_a

                 CALL sbpitr%get_field(labels,wpl,info)
                 or_fail("Failed to get field i_wp data.")

                 CALL sbpitr%get_field(image,wpi,info)
                 or_fail("Failed to get field r_wp data.")

#if   __DIME == __2D
                 IF (sbpitr%istart(1).EQ.1) THEN
                    DO j=lo_a(2),hi_a(2)
                       DO i=lo_a(1),0
                          wpl(i,j)=FORBIDDEN
                       ENDDO
                    ENDDO
                    DO j=lo_a(2),hi_a(2)
                       DO i=lo_a(1),0
                          wpi(i,j)=zero
                       ENDDO
                    ENDDO
                 ENDIF
                 IF (sbpitr%istart(2).EQ.1) THEN
                    DO j=lo_a(2),0
                       DO i=lo_a(1),hi_a(1)
                          wpl(i,j)=FORBIDDEN
                       ENDDO
                    ENDDO
                    DO j=lo_a(2),0
                       DO i=lo_a(1),hi_a(1)
                          wpi(i,j)=zero
                       ENDDO
                    ENDDO
                 ENDIF
                 IF (sbpitr%iend(1).EQ.Ngrid(1)) THEN
                    DO j=lo_a(2),hi_a(2)
                       DO i=Nm(1)+1,hi_a(1)
                          wpl(i,j)=FORBIDDEN
                       ENDDO
                    ENDDO
                    DO j=lo_a(2),hi_a(2)
                       DO i=Nm(1)+1,hi_a(1)
                          wpi(i,j)=zero
                       ENDDO
                    ENDDO
                 ENDIF
                 IF (sbpitr%iend(2).EQ.Ngrid(2)) THEN
                    DO j=Nm(2)+1,hi_a(2)
                       DO i=lo_a(1),hi_a(1)
                          wpl(i,j)=FORBIDDEN
                       ENDDO
                    ENDDO
                    DO j=Nm(2)+1,hi_a(2)
                       DO i=lo_a(1),hi_a(1)
                          wpi(i,j)=zero
                       ENDDO
                    ENDDO
                 ENDIF
#elif __DIME == __3D
                 IF (sbpitr%istart(1).EQ.1) THEN
                    DO k=lo_a(3),hi_a(3)
                       DO j=lo_a(2),hi_a(2)
                          DO i=lo_a(1),0
                             wpl(i,j,k)=FORBIDDEN
                          ENDDO
                       ENDDO
                    ENDDO
                    DO k=lo_a(3),hi_a(3)
                       DO j=lo_a(2),hi_a(2)
                          DO i=lo_a(1),0
                             wpi(i,j,k)=zero
                          ENDDO
                       ENDDO
                    ENDDO
                 ENDIF
                 IF (sbpitr%istart(2).EQ.1) THEN
                    DO k=lo_a(3),hi_a(3)
                       DO j=lo_a(2),0
                          DO i=lo_a(1),hi_a(1)
                             wpl(i,j,k)=FORBIDDEN
                          ENDDO
                       ENDDO
                    ENDDO
                    DO k=lo_a(3),hi_a(3)
                       DO j=lo_a(2),0
                          DO i=lo_a(1),hi_a(1)
                             wpi(i,j,k)=zero
                          ENDDO
                       ENDDO
                    ENDDO
                 ENDIF
                 IF (sbpitr%istart(3).EQ.1) THEN
                    DO k=lo_a(3),0
                       DO j=lo_a(2),hi_a(2)
                          DO i=lo_a(1),hi_a(1)
                             wpl(i,j,k)=FORBIDDEN
                          ENDDO
                       ENDDO
                    ENDDO
                    DO k=lo_a(3),0
                       DO j=lo_a(2),hi_a(2)
                          DO i=lo_a(1),hi_a(1)
                             wpi(i,j,k)=zero
                          ENDDO
                       ENDDO
                    ENDDO
                 ENDIF
                 IF (sbpitr%iend(1).EQ.Ngrid(1)) THEN
                    DO k=lo_a(3),hi_a(3)
                       DO j=lo_a(2),hi_a(2)
                          DO i=Nm(1)+1,hi_a(1)
                             wpl(i,j,k)=FORBIDDEN
                          ENDDO
                       ENDDO
                    ENDDO
                    DO k=lo_a(3),hi_a(3)
                       DO j=lo_a(2),hi_a(2)
                          DO i=Nm(1)+1,hi_a(1)
                             wpi(i,j,k)=zero
                          ENDDO
                       ENDDO
                    ENDDO
                 ENDIF
                 IF (sbpitr%iend(2).EQ.Ngrid(2)) THEN
                    DO k=lo_a(3),hi_a(3)
                       DO j=Nm(2)+1,hi_a(2)
                          DO i=lo_a(1),hi_a(1)
                             wpl(i,j,k)=FORBIDDEN
                          ENDDO
                       ENDDO
                    ENDDO
                    DO k=lo_a(3),hi_a(3)
                       DO j=Nm(2)+1,hi_a(2)
                          DO i=lo_a(1),hi_a(1)
                             wpi(i,j,k)=zero
                          ENDDO
                       ENDDO
                    ENDDO
                 ENDIF
                 IF (sbpitr%iend(3).EQ.Ngrid(3)) THEN
                    DO k=Nm(3)+1,hi_a(3)
                       DO j=lo_a(2),hi_a(2)
                          DO i=lo_a(1),hi_a(1)
                             wpl(i,j,k)=FORBIDDEN
                          ENDDO
                       ENDDO
                    ENDDO
                    DO k=Nm(3)+1,hi_a(3)
                       DO j=lo_a(2),hi_a(2)
                          DO i=lo_a(1),hi_a(1)
                             wpi(i,j,k)=zero
                          ENDDO
                       ENDDO
                    ENDDO
                 ENDIF
#endif

                 sbpitr => mesh%subpatch%next()
              ENDDO !WHILE(ASSOCIATED(sbpitr))

              NULLIFY(wpl,wpi)

              !-------------------------------------------------------------------------
   !---End----! complete ghost update by popping received data
              !-------------------------------------------------------------------------
              CALL mesh%map_isend(info,sendrecv=.FALSE.)
              or_fail("mesh%map_send")
              CALL image%map_ghost_pop(mesh,info)
              or_fail("image%map_ghost_pop")
#endif

           END SELECT !ppm_nproc
        END SELECT !MCMCcontinue

        !-------------------------------------------------------------------------
        ! Free the memory which we do not need nor use anymore
        ! This does not affect the new MCMC start
        !-------------------------------------------------------------------------
        CALL MCMCFreeMemory(info)
        or_fail("MCMCFreeMemory")

        !-------------------------------------------------------------------------
        ! TODO
        ! We should find the better solution
        !
        ! Make sure that no subdomain is bigger than
        ! 2047 in 3D and 65535 in 2D. The current limitation is due to the
        ! Hash structure which is currently used. We are hashing 4 integers
        ! (X,Y,Z,label) into one long integer to save memory in 3D
        !-------------------------------------------------------------------------
        sbpitr => mesh%subpatch%begin()
        DO WHILE (ASSOCIATED(sbpitr))
           Nm => sbpitr%nnodes
           IF (ppm_rc_dim.EQ.3) THEN
              check_true(<#ALL(Nm.LE.2046)#>, &
              & "The subdomain size (X,Y & Z <= 2047) limit is violated! Please use more resources for decomposition!")
           ELSE
              check_true(<#ALL(Nm.LE.65534)#>, &
              & "The subdomain size (X & Y <= 65535) limit is violated! Please use more resources for decomposition!")
           ENDIF
           sbpitr => mesh%subpatch%next()
        ENDDO

        !-------------------------------------------------------------------------
        ! Set up the edge image and the edge image discrete distribution
        ! the output is MCMCEdgeImageDistr which contains integer random
        ! distribution with probability of EdgeImage
        !-------------------------------------------------------------------------
!         CALL DTYPE(ppm_rc_EdgeDetection)(image,mesh,info)
!         or_fail("ppm_rc_EdgeDetection")

        !---------------------------------------------------------------------
        ! Prepare a fast proposal computation
        !---------------------------------------------------------------------
        CALL CreateMCMClengthProposalMask(info)
        or_fail("CreateMCMClengthProposalMask")

        !---------------------------------------------------------------------
        !  Determine cell structure for MCMC sampling
        !  MCMCcellsize: The size of the cell
        !  MCMCcellNm:   Number of cells in each direction
        !---------------------------------------------------------------------
        ALLOCATE(MCMCcellsize(ppm_rc_dim),MCMCcellNm(ppm_rc_dim),STAT=info)
        or_fail_alloc("MCMCcellsize & MCMCcellNm")

        MCMCcellsize=REAL(ghostsize_run,MK)

        sbpitr => mesh%subpatch%begin()
        DO WHILE (ASSOCIATED(sbpitr))
           Nm => sbpitr%nnodes

           !---------------------------------------------------------------------
           !  Determine number of cell boxes and effective cell size.
           !---------------------------------------------------------------------
           DO i=1,ppm_rc_dim
              ! number of cells based on a cellsize = cutoff
              MCMCcellNm(i)=INT((REAL(Nm(i),MK)+small*ten)/MCMCcellsize(i))
              ! make at least one box
              IF (MCMCcellNm(i).LT.1) MCMCcellNm(i)=1
              MCMCcellsize(i)=(REAL(Nm(i),MK)+small*ten)/REAL(MCMCcellNm(i),MK)
           ENDDO

           sbpitr => mesh%subpatch%next()
        ENDDO !WHILE(ASSOCIATED(sbpitr))

        !---------------------------------------------------------------------
        ! Number of boxes which belong to the boundary area of the domain
        !---------------------------------------------------------------------
#if   __DIME == __2D
        !MCMCghostcellSize=2X+2(Y-2)
        MCMCboundarycellSize=2*(MCMCcellNm(1)+MCMCcellNm(2)-2)
#elif __DIME == __3D
        !MCMCghostcellSize=2XY+2Y(Z-2)+2(X-2)(Z-2)
        MCMCboundarycellSize=2*MCMCcellNm(1)*MCMCcellNm(2)+     &
        &                    2*MCMCcellNm(2)*(MCMCcellNm(3)-2)+ &
        &                    2*(MCMCcellNm(1)-2)*(MCMCcellNm(3)-2)
#endif

        ! the size is 4 or 8 for 2D and 3D respectively.
        nsize=ISHFT(1,ppm_rc_dim)
        !---------------------------------------------------------------------
        ! Array contains all the boundary cell indices
        !---------------------------------------------------------------------
        ALLOCATE(MCMCboundarycellIndex(MCMCboundarycellSize),MCMCinteriorcellDisp(0:nsize),STAT=info)
        or_fail_alloc("MCMCboundarycellIndex & MCMCinteriorcellDisp")

        MCMCinteriorcellDisp=0
        !---------------------------------------------------------------------
        ! Fill the MCMCboundarycellIndex with boundary box indices
        !---------------------------------------------------------------------
        l=0
        nsize=0
#if   __DIME == __2D
        DO j=1,MCMCcellNm(2)
           ljeven=MOD(j,2)
           DO i=1,MCMCcellNm(1)
              nsize=nsize+1
              IF (i.EQ.1.OR.i.EQ.MCMCcellNm(1).OR. &
              &   j.EQ.1.OR.j.EQ.MCMCcellNm(2)) THEN
                 l=l+1
                 MCMCboundarycellIndex(l)=nsize
              ELSE
                 !odd-odd   --> color=1
                 !odd-even  --> color=2
                 !even-odd  --> color=3
                 !even-even --> color=4
                 lieven=MOD(i,2)
                 m=cellcoloring(lieven,ljeven)
                 MCMCinteriorcellDisp(m)=MCMCinteriorcellDisp(m)+1
              ENDIF
           ENDDO !i=1,MCMCcellNm(1)
        ENDDO !j=1,MCMCcellNm(2)
#elif __DIME == __3D
        DO k=1,MCMCcellNm(3)
           lkeven=MOD(k,2)
           DO j=1,MCMCcellNm(2)
              ljeven=MOD(j,2)
              DO i=1,MCMCcellNm(1)
                 nsize=nsize+1
                 IF (i.EQ.1.OR.i.EQ.MCMCcellNm(1).OR. &
                 &   j.EQ.1.OR.j.EQ.MCMCcellNm(2).OR. &
                 &   k.EQ.1.OR.k.EQ.MCMCcellNm(3)) THEN
                    l=l+1
                    MCMCboundarycellIndex(l)=nsize
                 ELSE
                    !odd-odd-odd    --> color=1
                    !odd-even-odd   --> color=2
                    !even-odd-odd   --> color=3
                    !even-even-odd  --> color=4
                    !odd-odd-even   --> color=5
                    !odd-even-even  --> color=6
                    !even-odd-even  --> color=7
                    !even-even-even --> color=8
                    lieven=MOD(i,2)
                    m=cellcoloring(lieven,ljeven,lkeven)
                    MCMCinteriorcellDisp(m)=MCMCinteriorcellDisp(m)+1
                 ENDIF
              ENDDO !i=1,MCMCcellNm(1)
           ENDDO !j=1,MCMCcellNm(2)
        ENDDO !k=1,MCMCcellNm(3)
#endif

        nsize=SUM(MCMCinteriorcellDisp)
        ALLOCATE(MCMCinteriorcellIndex(nsize),STAT=info)
        or_fail_alloc("MCMCinteriorcellIndex")

#if   __DIME == __2D
        ll(1)=0
        ll(2)=MCMCinteriorcellDisp(1)
        ll(3)=MCMCinteriorcellDisp(2)+ll(2)
        ll(4)=MCMCinteriorcellDisp(3)+ll(3)
        DO j=2,MCMCcellNm(2)-1
           ljeven=MOD(j,2)
           DO i=2,MCMCcellNm(1)-1
              lieven=MOD(i,2)
              l=cellcoloring(lieven,ljeven)
              ll(l)=ll(l)+1
              cbox=i+(j-1)*MCMCcellNm(1)
              MCMCinteriorcellIndex(ll(l))=cbox
           ENDDO
        ENDDO
        DO i=1,4
           MCMCinteriorcellDisp(i)=ll(i)
        ENDDO
#elif __DIME == __3D
        ll(1)=0
        ll(2)=MCMCinteriorcellDisp(1)
        ll(3)=MCMCinteriorcellDisp(2)+ll(2)
        ll(4)=MCMCinteriorcellDisp(3)+ll(3)
        ll(5)=MCMCinteriorcellDisp(4)+ll(4)
        ll(6)=MCMCinteriorcellDisp(5)+ll(5)
        ll(7)=MCMCinteriorcellDisp(6)+ll(6)
        ll(8)=MCMCinteriorcellDisp(7)+ll(7)
        DO k=2,MCMCcellNm(3)-1
           lkeven=MOD(k,2)
           DO j=2,MCMCcellNm(2)-1
              ljeven=MOD(j,2)
              DO i=2,MCMCcellNm(1)-1
                 lieven=MOD(i,2)
                 l=cellcoloring(lieven,ljeven,lkeven)
                 ll(l)=ll(l)+1
                 cbox=i+MCMCcellNm(1)*(j-1+(k-1)*MCMCcellNm(2))
                 MCMCinteriorcellIndex(ll(l))=cbox
              ENDDO
           ENDDO
        ENDDO
        DO i=1,8
           MCMCinteriorcellDisp(i)=ll(i)
        ENDDO
#endif
        !---------------------------------------------------------------------
        ! How many cells
        !---------------------------------------------------------------------
        nsize=PRODUCT(MCMCcellNm)
        ALLOCATE(MCMCRegularParticlesInCell(nsize),MCMCFloatingParticlesInCell(nsize),STAT=info)
        or_fail_alloc("MCMCRegularParticlesInCell & MCMCFloatingParticlesInCell")

        !---------------------------------------------------------------------
        ! How many regions
        !---------------------------------------------------------------------
        nsize=e_data%size()
        ALLOCATE(MCMCRegularParticles(nsize),STAT=info)
        or_fail_alloc("MCMCRegularParticles")

        ALLOCATE(MCMCRegularParticlesProposalNormalizer(nsize),      &
        &        MCMCRegularParticlesProposalNormalizerlocal(nsize), &
        &        SOURCE=zero,STAT=info)
        or_fail_alloc("MCMCRegularParticlesProposalNormalizer")

!         !---------------------------------------------------------------------
!         ! Initialize MCMCRegionLabel by putting the BG label in that
!         !---------------------------------------------------------------------
!         info=MCMCInsertLabelInRegionLabel()
!         or_fail("MCMCInsertLabelInRegionLabel failed to add BG label to the array!")

        !---------------------------------------------------------------------
        ! Creating the hash to insert all checked index and avoid checking them again
        !---------------------------------------------------------------------
        CALL htablei%create(512,info)
        or_fail("htablei%create")

        !---------------------------------------------------------------------
        ! Initialize the particles on the image label
        !---------------------------------------------------------------------
        NULLIFY(wpl)

        sbpitr => mesh%subpatch%begin()
        DO WHILE (ASSOCIATED(sbpitr))
           Nm => sbpitr%nnodes

           CALL sbpitr%get_field(labels,wpl,info)
           or_fail("Failed to get field i_wp data.")

#if   __DIME == __2D
           DO j=1,Nm(2)
              DO i=1,Nm(1)
                 IF (wpl(i,j).LT.0) THEN
                    vLabel=-wpl(i,j)

!                     info=MCMCInsertLabelInRegionLabel(vLabel)
!                     or_fail("MCMCInsertLabelInRegionLabel")

                    tmplabels_ => wpl(i-2:i+2,j-2:j+2)

                    IF (.NOT.htablei%search(i,j)) THEN
                       CALL htablei%insert(i,j,info)
                       or_fail("htablei%insert")

                       tmplabels => tmplabels_(2:4,2:4)

                       nsize=MCMCgetRegularParticlesAtIndex(tmplabels,i,j,Nm,MCMCParticles)

                       IF (nsize.GT.0) THEN
                          nx=INT(REAL(i,MK)/MCMCcellsize(1))
                          IF (nx.EQ.MCMCcellNm(1)) nx=nx-1

                          ny=INT(REAL(j,MK)/MCMCcellsize(2))
                          IF (ny.EQ.MCMCcellNm(2)) ny=ny-1

                          cbox=1+nx+ny*MCMCcellNm(1)

                          DO m=1,nsize
                             !---------------------------------------------------------------------
                             ! Returns false if topology is changed when applying this particle. This
                             ! is based on the current state of the label image.
                             ! TODO: this method should be changed to achieve full topological control.
                             ! now the particle is rejected if it changes somehow the topology.
                             !---------------------------------------------------------------------
                             IF (MCMCIsParticleTopoValid(tmplabels,MCMCParticles(m)%candlabel)) THEN
                                info=MCMCInsertCandidatesToContainers(i,j,cbox, &
                                &    MCMCParticles(m),vLabel,.FALSE.,Replaced)
                                or_fail("MCMCInsertCandidatesToContainers")
                             ENDIF
                          ENDDO !m=1,nsize
                       ENDIF !(nsize.GT.0)
                    ENDIF !(.NOT.htablei%search(i,j))

                    DO l=1,FG_ConnectivityType%NumberOfNeighbors
                       ld(1)=i+FG_ConnectivityType%NeighborsPoints(1,l)
                       ld(2)=j+FG_ConnectivityType%NeighborsPoints(2,l)

                       IF (ANY(ld.LT.1.OR.ld.GT.Nm)) CYCLE

                       nLabel=ABS(wpl(ld(1),ld(2)))
                       IF (nLabel.NE.vLabel.AND.nLabel.NE.FORBIDDEN) THEN
                          IF (.NOT.htablei%search(ld(1),ld(2))) THEN
                             CALL htablei%insert(ld(1),ld(2),info)
                             or_fail("htablei%insert")

                             SELECT CASE (l)
                             CASE (1)
                                tmplabels => tmplabels_(2:4,1:3)
                             CASE (2)
                                tmplabels => tmplabels_(1:3,2:4)
                             CASE (3)
                                tmplabels => tmplabels_(3:5,2:4)
                             CASE (4)
                                tmplabels => tmplabels_(2:4,3:5)
                             END SELECT

                             nsize=MCMCgetRegularParticlesAtIndex(tmplabels,ld(1),ld(2),Nm,MCMCParticles)

                             IF (nsize.GT.0) THEN
                                nx=INT(REAL(ld(1),MK)/MCMCcellsize(1))
                                IF (nx.EQ.MCMCcellNm(1)) nx=nx-1

                                ny=INT(REAL(ld(2),MK)/MCMCcellsize(2))
                                IF (ny.EQ.MCMCcellNm(2)) ny=ny-1

                                cbox=1+nx+ny*MCMCcellNm(1)

                                DO m=1,nsize
                                   IF (MCMCIsParticleTopoValid(tmplabels,MCMCParticles(m)%candlabel)) THEN
                                      info=MCMCInsertCandidatesToContainers(ld(1),ld(2),cbox, &
                                      &    MCMCParticles(m),nLabel,.FALSE.,Replaced)
                                      or_fail("MCMCInsertCandidatesToContainers")
                                   ENDIF
                                ENDDO !m=1,nsize
                             ENDIF !(nsize.GT.0)
                          ENDIF !.NOT.htablei%search(ld(1),ld(2))
                       ENDIF !(ABS(wpl(ld(1),ld(2))).NE.vLabel)
                    ENDDO !l=1,FG_ConnectivityType%NumberOfNeighbors
                 ENDIF !(wpl(i,j).LT.0)
              ENDDO !i=1,Nm(1)
           ENDDO !j=1,Nm(2)
#elif __DIME == __3D
           DO k=1,Nm(3)
              DO j=1,Nm(2)
                 DO i=1,Nm(1)
                    IF (wpl(i,j,k).LT.0) THEN
                       vLabel=-wpl(i,j,k)

!                        info=MCMCInsertLabelInRegionLabel(vLabel)
!                        or_fail("MCMCInsertLabelInRegionLabel")

                       tmplabels_ => wpl(i-2:i+2,j-2:j+2,k-2:k+2)

                       IF (.NOT.htablei%search(i,j,k)) THEN
                          CALL htablei%insert(i,j,k,info)
                          or_fail("htablei%insert")

                          tmplabels => tmplabels_(2:4,2:4,2:4)

                          nsize=MCMCgetRegularParticlesAtIndex(tmplabels,i,j,k,Nm,MCMCParticles)

                          IF (nsize.GT.0) THEN
                             nx=INT(REAL(i,MK)/MCMCcellsize(1))
                             IF (nx.EQ.MCMCcellNm(1)) nx=nx-1

                             ny=INT(REAL(j,MK)/MCMCcellsize(2))
                             IF (ny.EQ.MCMCcellNm(2)) ny=ny-1

                             nz=INT(REAL(k,MK)/MCMCcellsize(3))
                             IF (nz.EQ.MCMCcellNm(3)) nz=nz-1

                             cbox=1+nx+MCMCcellNm(1)*(ny+nz*MCMCcellNm(2))

                             DO m=1,nsize
                                !-------------------------------------------------------------------------
                                !!! Returns false if topology is changed when applying this particle. This
                                !!! is based on the current state of the label image.
                                !!! TODO: this method should be changed to achieve full topological control.
                                !!! now the particle is rejected if it changes somehow the topology.
                                !-------------------------------------------------------------------------
                                IF (MCMCIsParticleTopoValid(tmplabels,MCMCParticles(m)%candlabel)) THEN
                                   info=MCMCInsertCandidatesToContainers(i,j,k,cbox, &
                                   &    MCMCParticles(m),vLabel,.FALSE.,Replaced)
                                   or_fail("MCMCInsertCandidatesToContainers")
                                ENDIF
                             ENDDO !m=1,nsize
                          ENDIF !(nsize.GT.0)
                       ENDIF !(.NOT.htablei%search(i,j,k))

                       DO l=1,FG_ConnectivityType%NumberOfNeighbors
                          ld(1)=i+FG_ConnectivityType%NeighborsPoints(1,l)
                          ld(2)=j+FG_ConnectivityType%NeighborsPoints(2,l)
                          ld(3)=k+FG_ConnectivityType%NeighborsPoints(3,l)

                          IF (ANY(ld.LT.1.OR.ld.GT.Nm)) CYCLE

                          nLabel=ABS(wpl(ld(1),ld(2),ld(3)))
                          IF (nLabel.NE.vLabel.AND.nLabel.NE.FORBIDDEN) THEN
                             IF (.NOT.htablei%search(ld(1),ld(2),ld(3))) THEN
                                CALL htablei%insert(ld(1),ld(2),ld(3),info)
                                or_fail("htablei%insert")

                                SELECT CASE (l)
                                CASE (1)
                                   tmplabels => tmplabels_(2:4,2:4,1:3)
                                CASE (2)
                                   tmplabels => tmplabels_(2:4,1:3,2:4)
                                CASE (3)
                                   tmplabels => tmplabels_(1:3,2:4,2:4)
                                CASE (4)
                                   tmplabels => tmplabels_(3:5,2:4,2:4)
                                CASE (5)
                                   tmplabels => tmplabels_(2:4,3:5,2:4)
                                CASE (6)
                                   tmplabels => tmplabels_(2:4,2:4,3:5)
                                END SELECT

                                nsize=MCMCgetRegularParticlesAtIndex(tmplabels,ld(1),ld(2),ld(3),Nm,MCMCParticles)

                                IF (nsize.GT.0) THEN
                                   nx=INT(REAL(ld(1),MK)/MCMCcellsize(1))
                                   IF (nx.EQ.MCMCcellNm(1)) nx=nx-1

                                   ny=INT(REAL(ld(2),MK)/MCMCcellsize(2))
                                   IF (ny.EQ.MCMCcellNm(2)) ny=ny-1

                                   nz=INT(REAL(ld(3),MK)/MCMCcellsize(3))
                                   IF (nz.EQ.MCMCcellNm(3)) nz=nz-1

                                   cbox=1+nx+MCMCcellNm(1)*(ny+nz*MCMCcellNm(2))

                                   DO m=1,nsize
                                      IF (MCMCIsParticleTopoValid(tmplabels,MCMCParticles(m)%candlabel)) THEN
                                         info=MCMCInsertCandidatesToContainers(ld(1),ld(2),ld(3),cbox,   &
                                         &    MCMCParticles(m),ABS(wpl(ld(1),ld(2),ld(3))),.FALSE.,Replaced)
                                         or_fail("MCMCInsertCandidatesToContainers")
                                      ENDIF
                                   ENDDO !m=1,nsize
                                ENDIF !(nsize.GT.0)
                             ENDIF !.NOT.htablei%search(ld(1),ld(2),ld(3))
                          ENDIF !ABS(wpl(ld(1),ld(2),ld(3))).NE.vLabel
                       ENDDO !l=1,FG_ConnectivityType%NumberOfNeighbors
                    ENDIF !wpl(i,j,k).LT.0
                 ENDDO !i=1,Nm(1)
              ENDDO !j=1,Nm(2)
           ENDDO !k=1,Nm(3)
#endif

           sbpitr => mesh%subpatch%next()
        ENDDO !WHILE(ASSOCIATED(sbpitr))
        NULLIFY(wpl)

        ! Free memory
        CALL htablei%destroy(info)
        or_fail("htablei%destroy")

        ALLOCATE(ppm_t_field::Backuplabels,STAT=info)
        or_fail_alloc('Failed to allocate Backuplabels.')

        CALL Backuplabels%create(1,info,dtype=ppm_type_int,name="particle_index")
        or_fail("Create Backuplabels field failed!" )

        CALL Backuplabels%discretize_on(mesh,info)
        or_fail("Backuplabels discretize_on failed!")

        !-------------------------------------------------------------------------
        ! We need a copy of the label image for the reconstruction of the
        ! final results (marginals). It is a backup of the inital state.
        !-------------------------------------------------------------------------
        CALL DTYPE(ppm_rc_CopyImageAndNormalize)(labels, &
        &    Backuplabels,mesh,info,withGhost=.TRUE.)
        or_fail("ppm_rc_CopyImageAndNormalize")

#ifdef  __MPI
        sbpitr => mesh%subpatch%begin()
        nsize=0
        DO WHILE (ASSOCIATED(sbpitr))
           !-------------------------------------------------------------------------
           ! The size of ghost layer around each subpatch
           ! we compute only one layer of ghost around each patch
           !-------------------------------------------------------------------------
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

        CALL MPI_Ibarrier(comm,request,info)
        or_fail_MPI("MPI_Ibarrier")
#endif

        topo => ppm_topo(topoid)%t
        !-------------------------------------------------------------------------
        ! proccolor, is this processor color
        !-------------------------------------------------------------------------
        proccolor=topo%ineighcolor(0)

        !-------------------------------------------------------------------------
        ! Number of colors which will be used in MCMC sampling
        !-------------------------------------------------------------------------
        SELECT CASE (ppm_rc_dim)
        CASE (2)
           ALLOCATE(procflag(4),SOURCE=[1,2,3,4],STAT=info)
        CASE (3)
           ALLOCATE(procflag(8),SOURCE=[1,2,3,4,5,6,7,8],STAT=info)
        END SELECT
        or_fail_alloc("procflag")
#ifdef  __MPI
        !wait for everyone to pass the barrier
        CALL MPI_Wait(request,MPI_STATUS_IGNORE,info)
        or_fail_MPI("MPI_Wait")
#endif
        !-------------------------------------------------------------------------
        !  End of init
        !-------------------------------------------------------------------------
        IF (rank.EQ.0) THEN
           stdout("MCMC initialization is complete.")
        ENDIF
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
      9999 CONTINUE
        CALL substop(caller,t0,info)
        RETURN
        CONTAINS
#if   __DIME == __2D
          INTEGER FUNCTION cellcoloring(lx,ly)
#elif __DIME == __3D
          INTEGER FUNCTION cellcoloring(lx,ly,lz)
#endif
          IMPLICIT NONE
          INTEGER, INTENT(IN   ) :: lx
          INTEGER, INTENT(IN   ) :: ly
#if   __DIME == __3D
          INTEGER, INTENT(IN   ) :: lz
#endif
#if   __DIME == __2D
          INTEGER, DIMENSION(0:3), PARAMETER :: ccl=[4,3,2,1]
          !odd-odd   --> 11 --> 3 --> color=1
          !odd-even  --> 10 --> 2 --> color=2
          !even-odd  --> 01 --> 1 --> color=3
          !even-even --> 00 --> 0 --> color=4
          cellcoloring=ccl(IOR(ISHFT(lx,1),ly))
#elif __DIME == __3D
          INTEGER, DIMENSION(0:7), PARAMETER :: ccl=[8,4,7,3,6,2,5,1]

          !odd-odd-odd    --> 111 --> 7 --> color=1
          !odd-even-odd   --> 101 --> 5 --> color=2
          !even-odd-odd   --> 011 --> 3 --> color=3
          !even-even-odd  --> 001 --> 1 --> color=4
          !odd-odd-even   --> 110 --> 6 --> color=5
          !odd-even-even  --> 100 --> 4 --> color=6
          !even-odd-even  --> 010 --> 2 --> color=7
          !even-even-even --> 000 --> 0 --> color=8
          cellcoloring=ccl(IOR(ISHFT(lx,2),IOR(ISHFT(ly,1),lz)))
#endif
          RETURN
          END FUNCTION cellcoloring
      END SUBROUTINE DTYPE(ppm_rc_init_mcmc)
