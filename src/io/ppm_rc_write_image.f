      !-------------------------------------------------------------------------
      !  Subroutine   :                    ppm_rc_write_image
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
      !  Author           - y.afshar           June   2014
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Subroutine   :                    ppm_rc_write_image
      !-------------------------------------------------------------------------
      !
      !  Purpose      :
      !
      !  Input        :
      !
      !  Output       : info          (I) return status. 0 on success.
      !
      !  Routines     :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------

      SUBROUTINE DTYPE(ppm_rc_write_image)(FieldIn,MeshIn,OutputFileName,info,MeshOut,idn)

        !-------------------------------------------------------------------------
        !  Modules
        !-------------------------------------------------------------------------
        USE ppm_module_topo_typedef, ONLY : ppm_t_topo,ppm_topo
        USE ppm_module_interfaces, ONLY : ppm_t_discr_info_

        USE ppm_rc_module_tiff, ONLY : ppm_rc_open_write_tiff,ppm_rc_open_write_bigtiff, &
        &   ppm_rc_close_tiff,ppm_rc_write_tiff_header,ppm_rc_write_tiff_strip,bitsPerSampleW,  &
        &   ppm_rc_write_tiff_scanline
        IMPLICIT NONE

        !-------------------------------------------------------------------------
        !  Includes
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
        !  Externals
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        CLASS(ppm_t_field_),               POINTER       :: FieldIn

        CLASS(ppm_t_equi_mesh_),           POINTER       :: MeshIn

        CHARACTER(LEN=*),                  INTENT(IN   ) :: OutputFileName

        INTEGER,                           INTENT(  OUT) :: info

        CLASS(ppm_t_equi_mesh_), OPTIONAL, POINTER       :: MeshOut

        INTEGER,                 OPTIONAL, INTENT(IN   ) :: idn

        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
        TYPE(ppm_t_topo), POINTER :: topo

        CLASS(ppm_t_discr_info_), POINTER :: tmp_dinfo

        CLASS(ppm_t_equi_mesh_), POINTER :: MeshOut_

#if   __DIME == __2D
        REAL(ppm_kind_single), CONTIGUOUS, DIMENSION(:,:),   POINTER :: DTYPE(wp_rs) => NULL()
        REAL(MK),              CONTIGUOUS, DIMENSION(:,:),   POINTER :: DTYPE(wp_r) => NULL()
#elif __DIME == __3D
        REAL(ppm_kind_single), CONTIGUOUS, DIMENSION(:,:,:), POINTER :: DTYPE(wp_rs) => NULL()
        REAL(MK),              CONTIGUOUS, DIMENSION(:,:,:), POINTER :: DTYPE(wp_r) => NULL()
#endif
        REAL(ppm_kind_single), CONTIGUOUS, DIMENSION(:,:),   POINTER :: bufrs => NULL()
        REAL(MK),              CONTIGUOUS, DIMENSION(:,:),   POINTER :: bufr => NULL()
        REAL(ppm_kind_double) :: t0
        REAL(ppm_kind_single) :: Normalfacs

#if   __DIME == __2D
        INTEGER, CONTIGUOUS, DIMENSION(:,:),   POINTER :: DTYPE(wp_i) => NULL()
#elif __DIME == __3D
        INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER :: DTYPE(wp_i) => NULL()
#endif
        INTEGER, CONTIGUOUS, DIMENSION(:,:),   POINTER :: bufi => NULL()
        INTEGER,             DIMENSION(:),     POINTER :: Nm
        INTEGER,             DIMENSION(:),     POINTER :: istart
        ! Lower-left coordinates on the global mesh
        INTEGER                                        :: isub,x,y,i,ipatch
        INTEGER                                        :: istrip
        INTEGER                                        :: ldu(2),iopt

        CHARACTER(LEN=ppm_char) :: filename
        CHARACTER(LEN=ppm_char) :: caller = 'ppm_rc_write_image'

        LOGICAL :: IsBigTIFF
        !-------------------------------------------------------------------------
        !  Initialize
        !-------------------------------------------------------------------------
        CALL substart(caller,t0,info)

        CALL check()

        NULLIFY(tmp_dinfo)

        IF (PRESENT(MeshOut)) THEN
           IF (.NOT.FieldIn%is_discretized_on(MeshOut)) THEN
              CALL FieldIn%discretize_on(MeshOut,info,discr_info=tmp_dinfo)
              or_fail("FieldIn discretize_on failed!")
           ENDIF

           IF (.NOT.ASSOCIATED(MeshIn,MeshOut)) THEN
              !---------------------------------------------------------------------
              ! Global mapping to the new mesh in the new topology
              !---------------------------------------------------------------------
              CALL MeshIn%map(MeshOut,info)
              or_fail("Failed to do global mapping from old to the new topology.")

              CALL FieldIn%map_push(MeshIn,info)
              or_fail("Failed to push field data in the buffer!")
              CALL MeshIn%map_isend(info)
              or_fail("Failed to send field data to the new mesh!")
              CALL FieldIn%map_pop(MeshOut,info)
              or_fail("Failed to pop field data in the new mesh!")
           ENDIF

           MeshOut_ => MeshOut
        ELSE
           MeshOut_ => MeshIn
        ENDIF

        !-------------------------------------------------------------------------
        !  Create format string
        !-------------------------------------------------------------------------
        topo => ppm_topo(MeshOut_%topoid)%t

        !---------------------------------------------------------------------
        ! First we find out which image slice to read
        ! We use the istart array for that
        !---------------------------------------------------------------------
        sub_loop: DO i=1,topo%nsublist

           isub = topo%isublist(i)

           patch_loop: DO ipatch=1,MeshOut_%subpatch_by_sub(isub)%nsubpatch

              SELECT TYPE(p => MeshOut_%subpatch_by_sub(isub)%vec(ipatch)%t)
              TYPE IS (ppm_t_subpatch)
                 Nm     => p%nnodes
                 istart => p%istart

                 IF (PRESENT(idn)) THEN
#if   __DIME == __2D
                    WRITE(filename,'(A,A1,I0,A2,I0,A1,I0,A2,I0,A1,I0,A4,A1)') &
                    & TRIM(OutputFileName),'_',rank,'__',istart(1),'_',istart(2),'__',ipatch,'_',idn,'.tif',CHAR(0)
!                     & TRIM(OutputFileName),'_',isub,'__',istart(1),'_',istart(2),'__',ipatch,'_',idn,'.tif',CHAR(0)
#elif __DIME == __3D
                    WRITE(filename,'(A,A1,I0,A2,I0,A1,I0,A1,I0,A2,I0,A1,I0,A4,A1)') &
                    & TRIM(OutputFileName),'_',rank,'__',istart(1),'_',istart(2),'_',istart(3),'__',ipatch,'_',idn,'.tif',CHAR(0)
!                     & TRIM(OutputFileName),'_',isub,'__',istart(1),'_',istart(2),'_',istart(3),'__',ipatch,'_',idn,'.tif',CHAR(0)
#endif
                 ELSE
#if   __DIME == __2D
                    WRITE(filename,'(A,A1,I0,A2,I0,A1,I0,A2,I0,A4,A1)') &
                    & TRIM(OutputFileName),'_',rank,'__',istart(1),'_',istart(2),'__',ipatch,'.tif',CHAR(0)
!                     & TRIM(OutputFileName),'_',isub,'__',istart(1),'_',istart(2),'__',ipatch,'.tif',CHAR(0)
#elif __DIME == __3D
                    WRITE(filename,'(A,A1,I0,A2,I0,A1,I0,A1,I0,A2,I0,A4,A1)') &
                    & TRIM(OutputFileName),'_',rank,'__',istart(1),'_',istart(2),'_',istart(3),'__',ipatch,'.tif',CHAR(0)
!                     & TRIM(OutputFileName),'_',isub,'__',istart(1),'_',istart(2),'_',istart(3),'__',ipatch,'.tif',CHAR(0)
#endif
                 ENDIF

                 !The TIFF file format uses 32bit offsets and, as such, is limited to 4 gigabytes
#if   __DIME == __2D
                 IsBigTIFF=REAL(Nm(1),MK)*REAL(Nm(2),MK)*REAL(bitsPerSampleW/8,MK)/1024._MK.GE.4194304._MK
#elif __DIME == __3D
                 IsBigTIFF=REAL(Nm(1),MK)*REAL(Nm(2),MK)*REAL(Nm(3),MK)*REAL(bitsPerSampleW/8,MK)/1024._MK.GE.4194304._MK
#endif

                 SELECT CASE (IsBigTIFF)
                 CASE (.TRUE.)
                    info=ppm_rc_open_write_bigtiff(filename)
                 CASE DEFAULT
                    info=ppm_rc_open_write_tiff(filename)
                 END SELECT !(IsBigTIFF)
                 or_fail('Error opening tiff file to write in.')

                 info=ppm_rc_write_tiff_header(bitsPerSampleW,1,Nm(1),Nm(2))
                 or_fail('Error writing tiff file header.')

                 SELECT CASE (FieldIn%data_type)
                 CASE (ppm_type_int)
                    CALL p%get_field(FieldIn,DTYPE(wp_i),info)
                    or_fail("Failed to get field wp_i.")

#if   __DIME == __2D
                    info=ppm_rc_write_tiff_strip(DTYPE(wp_i)(1:Nm(1),1:Nm(2)),bitsPerSampleW,Nm(1),Nm(2),0,1)
                    or_fail("ppm_rc_write_tiff_strip")
#elif __DIME == __3D
                    DO istrip=1,Nm(3)
                       info=ppm_rc_write_tiff_header(istrip-1,Nm(3))
                       or_fail('Error writing tiff file header.')

                       info=ppm_rc_write_tiff_strip(DTYPE(wp_i)(1:Nm(1),1:Nm(2),istrip),bitsPerSampleW,Nm(1),Nm(2),istrip-1,Nm(3))
                       or_fail("ppm_rc_write_tiff_strip")
                    ENDDO ! loop through slices
#endif

                 CASE (ppm_type_real)
                 !!! For the real data type we do check whether the data is been normalized
                 !!! or not, also in the case of real data, it would be the image dada itself
                 !!! so we do not need to correct the range of data
                    CALL p%get_field(FieldIn,DTYPE(wp_r),info)
                    or_fail("Failed to get field wp_r.")

                    IF (lNormalize) THEN
                       iopt=ppm_param_alloc_fit
                       ldu(1:2)=Nm(1:2)
                       SELECT CASE (bitsPerSampleW)
                       CASE (32)
                          CALL ppm_alloc(bufr,ldu,iopt,info)

                       CASE DEFAULT
                          CALL ppm_alloc(bufi,ldu,iopt,info)

                       END SELECT
                       or_fail_alloc('Failed to allocate C buffer.')
                    ENDIF !(lNormalize)

#if   __DIME == __2D
                    SELECT CASE (lNormalize)
                    CASE (.TRUE.)
                       SELECT CASE (bitsPerSampleW)
                       CASE (32)
                          FORALL (x=1:Nm(1),y=1:Nm(2))
                             bufr(x,y)=DTYPE(wp_r)(x,y)*Normalfac
                          END FORALL

                          info=ppm_rc_write_tiff_strip(bufr,bitsPerSampleW,Nm(1),Nm(2),0,1)
                          or_fail("ppm_rc_write_tiff_strip")

                       CASE DEFAULT
                          FORALL (x=1:Nm(1),y=1:Nm(2))
                             bufi(x,y)=INT(DTYPE(wp_r)(x,y)*Normalfac)
                          END FORALL

                          info=ppm_rc_write_tiff_strip(bufi,bitsPerSampleW,Nm(1),Nm(2),0,1)
                          or_fail("ppm_rc_write_tiff_strip")

                       END SELECT !(bitsPerSampleW)

                    CASE (.FALSE.)
                       info=ppm_rc_write_tiff_strip(DTYPE(wp_r)(1:Nm(1),1:Nm(2)),bitsPerSampleW,Nm(1),Nm(2),0,1)
                       or_fail("ppm_rc_write_tiff_strip")

                    END SELECT !(lNormalize)
#elif __DIME == __3D
                    DO istrip=1,Nm(3)
                       info=ppm_rc_write_tiff_header(istrip-1,Nm(3))
                       or_fail('Error writing tiff file header.')

                       SELECT CASE (lNormalize)
                       CASE (.TRUE.)
                          SELECT CASE (bitsPerSampleW)
                          CASE (32)
                             FORALL (x=1:Nm(1),y=1:Nm(2))
                                bufr(x,y)=DTYPE(wp_r)(x,y,istrip)*Normalfac
                             END FORALL

                             info=ppm_rc_write_tiff_strip(bufr,bitsPerSampleW,Nm(1),Nm(2),istrip-1,Nm(3))
                             or_fail("ppm_rc_write_tiff_strip")
!                              info=ppm_rc_write_tiff_scanline(bufr,0,bitsPerSampleW,0,Nm(1),Nm(2),istrip-1,Nm(3))
!                              or_fail("ppm_rc_write_tiff_scanline")

                          CASE DEFAULT
                             FORALL (x=1:Nm(1),y=1:Nm(2))
                                bufi(x,y)=INT(DTYPE(wp_r)(x,y,istrip)*Normalfac)
                             END FORALL

                             info=ppm_rc_write_tiff_strip(bufi,bitsPerSampleW,Nm(1),Nm(2),istrip-1,Nm(3))
                             or_fail("ppm_rc_write_tiff_strip")

                          END SELECT !(bitsPerSampleW)

                       CASE (.FALSE.)
                          info=ppm_rc_write_tiff_strip(DTYPE(wp_r)(1:Nm(1),1:Nm(2),istrip),bitsPerSampleW,Nm(1),Nm(2),istrip-1,Nm(3))
                          or_fail("ppm_rc_write_tiff_strip")

                       END SELECT !(lNormalize)
                    ENDDO ! loop through slices
#endif

                 CASE (ppm_type_real_single)
                 !!! For the real data type we do check whether the data is been normalized
                 !!! or not, also in the case of real data, it would be the image dada itself
                 !!! so we do not need to correct the range of data
                    CALL p%get_field(FieldIn,DTYPE(wp_rs),info)
                    or_fail("Failed to get field wp_rs.")

                    IF (lNormalize) THEN
                       iopt=ppm_param_alloc_fit
                       ldu(1:2)=Nm(1:2)
                       SELECT CASE (bitsPerSampleW)
                       CASE (32)
                          CALL ppm_alloc(bufrs,ldu,iopt,info)

                       CASE DEFAULT
                          CALL ppm_alloc(bufi,ldu,iopt,info)

                       END SELECT
                       or_fail_alloc('Failed to allocate C buffer.')

                       Normalfacs=REAL(Normalfac,ppm_kind_single)
                    ENDIF !(lNormalize)

#if   __DIME == __2D
                    SELECT CASE (lNormalize)
                    CASE (.TRUE.)
                       SELECT CASE (bitsPerSampleW)
                       CASE (32)
                          FORALL (x=1:Nm(1),y=1:Nm(2))
                             bufrs(x,y)=DTYPE(wp_rs)(x,y)*Normalfacs
                          END FORALL

                          info=ppm_rc_write_tiff_strip(bufrs,bitsPerSampleW,Nm(1),Nm(2),0,1)
                          or_fail("ppm_rc_write_tiff_strip")

                       CASE DEFAULT
                          FORALL (x=1:Nm(1),y=1:Nm(2))
                             bufi(x,y)=INT(DTYPE(wp_rs)(x,y)*Normalfacs)
                          END FORALL

                          info=ppm_rc_write_tiff_strip(bufi,bitsPerSampleW,Nm(1),Nm(2),0,1)
                          or_fail("ppm_rc_write_tiff_strip")

                       END SELECT !(bitsPerSampleW)

                    CASE (.FALSE.)
                       info=ppm_rc_write_tiff_strip(DTYPE(wp_rs)(1:Nm(1),1:Nm(2)),bitsPerSampleW,Nm(1),Nm(2),0,1)
                       or_fail("ppm_rc_write_tiff_strip")

                    END SELECT !(lNormalize)
#elif __DIME == __3D
                    DO istrip=1,Nm(3)
                       info=ppm_rc_write_tiff_header(istrip-1,Nm(3))
                       or_fail('Error writing tiff file header.')

                       SELECT CASE (lNormalize)
                       CASE (.TRUE.)
                          SELECT CASE (bitsPerSampleW)
                          CASE (32)
                             FORALL (x=1:Nm(1),y=1:Nm(2))
                                bufrs(x,y)=DTYPE(wp_rs)(x,y,istrip)*Normalfacs
                             END FORALL

                             info=ppm_rc_write_tiff_strip(bufrs,bitsPerSampleW,Nm(1),Nm(2),istrip-1,Nm(3))
                             or_fail("ppm_rc_write_tiff_strip")
!                              info=ppm_rc_write_tiff_scanline(bufr,0,bitsPerSampleW,0,Nm(1),Nm(2),istrip-1,Nm(3))
!                              or_fail("ppm_rc_write_tiff_scanline")

                          CASE DEFAULT
                             FORALL (x=1:Nm(1),y=1:Nm(2))
                                bufi(x,y)=INT(DTYPE(wp_rs)(x,y,istrip)*Normalfac)
                             END FORALL

                             info=ppm_rc_write_tiff_strip(bufi,bitsPerSampleW,Nm(1),Nm(2),istrip-1,Nm(3))
                             or_fail("ppm_rc_write_tiff_strip")

                          END SELECT !(bitsPerSampleW)

                       CASE (.FALSE.)
                          info=ppm_rc_write_tiff_strip(DTYPE(wp_rs)(1:Nm(1),1:Nm(2),istrip),bitsPerSampleW,Nm(1),Nm(2),istrip-1,Nm(3))
                          or_fail("ppm_rc_write_tiff_strip")

                       END SELECT !(lNormalize)
                    ENDDO ! loop through slices
#endif

                 CASE DEFAULT
                    fail("FieldIn data_type is not supported")

                 END SELECT !(FieldIn%data_type)

                 info=ppm_rc_close_tiff()
                 or_fail('Error closing tiff file.')

              END SELECT !SELECT TYPE

           ENDDO patch_loop

        ENDDO sub_loop ! loop through subs

        iopt=ppm_param_dealloc
        CALL ppm_alloc(bufi,ldu,iopt,info)
        CALL ppm_alloc(bufr,ldu,iopt,info)
        CALL ppm_alloc(bufrs,ldu,iopt,info)
        or_fail_dealloc('Could not deallocate C buffer.')

        NULLIFY(DTYPE(wp_r),DTYPE(wp_rs),DTYPE(wp_i),bufi,bufr)

        IF (ASSOCIATED(tmp_dinfo)) THEN
           CALL MeshOut%field_ptr%remove(info,FieldIn)
           or_fail("MeshOut%field_ptr%remove")

           CALL FieldIn%discr_info%remove(info,tmp_dinfo)
           or_fail("FieldIn%discr_info%remove")

           DEALLOCATE(tmp_dinfo,STAT=info)
           or_fail_dealloc("Failed to deallocate tmp_dinfo")

           NULLIFY(tmp_dinfo)
        ENDIF

        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
      9999 CONTINUE
        CALL substop(caller,t0,info)
        RETURN
      CONTAINS
        SUBROUTINE check
          IMPLICIT NONE
          check_true(<#ASSOCIATED(MeshIn)#>, &
          & "Input Mesh does not exist!",exit_point=8888)

          check_true(<#FieldIn%is_discretized_on(MeshIn)#>, &
          & "This Field is not discretized on Input Mesh!",exit_point=8888)

          IF (PRESENT(MeshOut)) THEN
             check_true(<#ASSOCIATED(MeshOut)#>, &
             & "Output Mesh does not exist!",exit_point=8888)
          ENDIF

        8888 CONTINUE
          RETURN
        END SUBROUTINE
      END SUBROUTINE DTYPE(ppm_rc_write_image)


      SUBROUTINE DTYPE(ppm_rc_write_image_label)(FieldIn,MeshIn, &
      &          OutputFileName,info,MeshOut,idn,liotopo)

        !-------------------------------------------------------------------------
        !  Modules
        !-------------------------------------------------------------------------
        USE ppm_module_data, ONLY : ppm_param_assign_internal, &
        & ppm_param_decomp_xpencil,ppm_param_decomp_xy_slab
        USE ppm_module_mktopo, ONLY : ppm_mktopo
        USE ppm_module_topo_typedef, ONLY : ppm_t_topo,ppm_topo
        USE ppm_module_interfaces, ONLY : ppm_t_discr_info_

        USE ppm_rc_module_tiff, ONLY : ppm_rc_open_write_tiff,ppm_rc_open_write_bigtiff, &
        &   ppm_rc_close_tiff,ppm_rc_write_tiff_header,ppm_rc_write_tiff_strip !,bitsPerSampleW
        IMPLICIT NONE

        !-------------------------------------------------------------------------
        !  Includes
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        CLASS(ppm_t_field_),               POINTER       :: FieldIn

        CLASS(ppm_t_equi_mesh_),           POINTER       :: MeshIn

        CHARACTER(LEN=*),                  INTENT(IN   ) :: OutputFileName

        INTEGER,                           INTENT(  OUT) :: info

        CLASS(ppm_t_equi_mesh_), OPTIONAL, POINTER       :: MeshOut

        INTEGER,                 OPTIONAL, INTENT(IN   ) :: idn
        !!!The image number which will be used for header information
        !!!The output file will use this header information
        LOGICAL,                 OPTIONAL, INTENT(IN   ) :: liotopo
        !!!This option set to TRUE will create IO topo IO mesh and write the
        !!!label on the new iomesh and new iotopo
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
        TYPE(ppm_t_topo), POINTER :: topo

        CLASS(ppm_t_discr_info_), POINTER :: tmp_dinfo

        CLASS(ppm_t_equi_mesh_), POINTER :: MeshOut_
        CLASS(ppm_t_equi_mesh_), POINTER :: tiomesh
        !!! mesh for I/O


        REAL(MK), DIMENSION(:,:), POINTER :: xp
        REAL(ppm_kind_double)             :: t0

#if   __DIME == __2D
        INTEGER, CONTIGUOUS, DIMENSION(:,:),   POINTER :: DTYPE(wp_i) => NULL()
#elif __DIME == __3D
        INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER :: DTYPE(wp_i) => NULL()
#endif
        INTEGER, CONTIGUOUS, DIMENSION(:,:),   POINTER :: buf => NULL()
        INTEGER,             DIMENSION(:),     POINTER :: Nm
        INTEGER,             DIMENSION(:),     POINTER :: istart
        INTEGER                                        :: isub,x,y,i,ipatch
        INTEGER                                        :: white
#if   __DIME == __3D
        INTEGER                                        :: istrip
#endif
        INTEGER                                        :: ldu(2),iopt
        INTEGER                                        :: nproc,ld(__DIME)
        INTEGER                                        :: tiotopoid,tiomeshid

        CHARACTER(LEN=ppm_char) :: filename
        CHARACTER(LEN=ppm_char) :: caller = 'ppm_rc_write_image_label'

        LOGICAL :: IsBigTIFF
        LOGICAL :: labelInit
        !-------------------------------------------------------------------------
        !  Initialize
        !-------------------------------------------------------------------------
        CALL substart(caller,t0,info)

        CALL check()

        !Write the label image file for using as an initialization input
        labelInit=.FALSE.

        NULLIFY(MeshOut_)

        IF (PRESENT(liotopo)) THEN
           IF (liotopo.AND.ppm_nproc.GT.1) THEN
              !-------------------------------------------------------------------------
              !  Create an IO topology for reading the image file
              !-------------------------------------------------------------------------
              nproc=ppm_nproc

              IF (n_io_procs_read.GT.1.AND.n_io_procs_read.LT.ppm_nproc) THEN
                 ppm_nproc=n_io_procs_read
                 IF (Ngrid(ppm_rc_dim).LE.ppm_nproc) THEN
                    ppm_nproc=Ngrid(ppm_rc_dim)-1
                 ENDIF
              ELSE IF (n_io_procs_read.GE.ppm_nproc) THEN
                 IF (Ngrid(ppm_rc_dim).LE.ppm_nproc) THEN
                    ppm_nproc=Ngrid(ppm_rc_dim)-1
                 ENDIF
              ELSE
                 !   ppm_nproc   n_io_procs
                 ! +------------------------+
                 ! |     2               2  |
                 ! |     4               2  |
                 ! |     8               2  |
                 ! |    16               8  |
                 ! |    32               8  |
                 ! |    64               8  |
                 ! |   128               8  |
                 ! |   256               8  |
                 ! |   512              16  |
                 ! |  1024              16  |
                 ! |  2048              16  |
                 ! +------------------------+
                 SELECT CASE (ppm_nproc)
                 CASE (2:8)
                    ppm_nproc=MAX(INT(LOG(REAL(ppm_nproc))/LOG(2.0)),2)
                    DO WHILE (IAND(ppm_nproc,ppm_nproc-1).NE.0)
                       ppm_nproc=ppm_nproc-1
                    ENDDO

                 CASE DEFAULT
                    ppm_nproc=INT(LOG(REAL(ppm_nproc))/LOG(2.0))
                    DO WHILE (IAND(ppm_nproc,ppm_nproc-1).NE.0)
                       ppm_nproc=ppm_nproc+1
                    ENDDO

                 END SELECT
              ENDIF

              ALLOCATE(xp(__DIME,1),STAT=info)
              or_fail_alloc("xp")

              xp=zero
              ld=0

              assig  = ppm_param_assign_internal
              decomp = MERGE(ppm_param_decomp_xpencil,ppm_param_decomp_xy_slab,ppm_rc_dim.EQ.2)

              tiotopoid=0
              tiomeshid=-1

              CALL ppm_mktopo(tiotopoid,tiomeshid,xp,0,decomp,assig, &
              &    min_phys,max_phys,bcdef,ld,cost,Ngrid,info)
              or_fail('Failed to create new topology.')

              DEALLOCATE(xp,STAT=info)
              or_fail_dealloc("xp")

              ppm_nproc=nproc

              tiomesh => ppm_mesh%at(tiomeshid)
              ! IO mesh to read and write image file

              CALL tiomesh%def_uniform(info)
              or_fail("Failed to create uniform tiomesh.")
              ! Create uniform mesh over the domain

              CALL FieldIn%discretize_on(tiomesh,info)
              or_fail("FieldIn discretize_on failed!")

              !---------------------------------------------------------------------
              ! Global mapping to the new mesh in the new topology
              !---------------------------------------------------------------------
              CALL MeshIn%map(tiomesh,info)
              or_fail("Failed to do global mapping from old to the new topology.")

              CALL FieldIn%map_push(MeshIn,info)
              or_fail("Failed to push field data in the buffer!")
              CALL MeshIn%map_isend(info)
              or_fail("Failed to send field data to the new mesh!")
              CALL FieldIn%map_pop(tiomesh,info)
              or_fail("Failed to pop field data in the new mesh!")

              MeshOut_ => tiomesh
           ENDIF !liotopo.AND.ppm_nproc.GT.1
        ENDIF !PRESENT(liotopo)

        NULLIFY(tmp_dinfo)

        IF (PRESENT(MeshOut)) THEN
           IF (.NOT.FieldIn%is_discretized_on(MeshOut)) THEN
              CALL FieldIn%discretize_on(MeshOut,info,discr_info=tmp_dinfo)
              or_fail("FieldIn discretize_on failed!")
           ENDIF

           IF (.NOT.ASSOCIATED(MeshIn,MeshOut)) THEN
              !---------------------------------------------------------------------
              ! Global mapping to the new mesh in the new topology
              !---------------------------------------------------------------------
              CALL MeshIn%map(MeshOut,info)
              or_fail("Failed to do global mapping from old to the new topology.")

              CALL FieldIn%map_push(MeshIn,info)
              or_fail("Failed to push field data in the buffer!")
              CALL MeshIn%map_isend(info)
              or_fail("Failed to send field data to the new mesh!")
              CALL FieldIn%map_pop(MeshOut,info)
              or_fail("Failed to pop field data in the new mesh!")
           ENDIF

           MeshOut_ => MeshOut
        ELSE
           IF (.NOT.ASSOCIATED(MeshOut_)) THEN
              MeshOut_ => MeshIn
           ENDIF
        ENDIF

        !-------------------------------------------------------------------------
        !  Create format string
        !-------------------------------------------------------------------------
        topo => ppm_topo(MeshOut_%topoid)%t

        !---------------------------------------------------------------------
        ! First we find out which image slice to read
        ! We use the istart array for that
        !---------------------------------------------------------------------
        sub_loop: DO i=1,topo%nsublist

           isub = topo%isublist(i)

           patch_loop: DO ipatch=1,MeshOut_%subpatch_by_sub(isub)%nsubpatch

              SELECT TYPE(p => MeshOut_%subpatch_by_sub(isub)%vec(ipatch)%t)
              TYPE IS (ppm_t_subpatch)
                 Nm     => p%nnodes
                 istart => p%istart

                 IF (PRESENT(idn)) THEN
                    IF (idn.EQ.FORBIDDEN) labelInit=.TRUE.

#if   __DIME == __2D
                    WRITE(filename,'(A,A1,I0,A2,I0,A1,I0,A2,I0,A1,I0,A4,A1)') &
                    & TRIM(OutputFileName),'_',rank,'__',istart(1),'_',istart(2),'__',ipatch,'_',idn,'.tif',CHAR(0)
!                     & TRIM(OutputFileName),'_',isub,'__',istart(1),'_',istart(2),'__',ipatch,'_',idn,'.tif',CHAR(0)
#elif __DIME == __3D
                    WRITE(filename,'(A,A1,I0,A2,I0,A1,I0,A1,I0,A2,I0,A1,I0,A4,A1)') &
                    & TRIM(OutputFileName),'_',rank,'__',istart(1),'_',istart(2),'_',istart(3),'__',ipatch,'_',idn,'.tif',CHAR(0)
!                     & TRIM(OutputFileName),'_',isub,'__',istart(1),'_',istart(2),'_',istart(3),'__',ipatch,'_',idn,'.tif',CHAR(0)
#endif
                 ELSE
#if   __DIME == __2D
                    WRITE(filename,'(A,A1,I0,A2,I0,A1,I0,A2,I0,A4,A1)') &
                    & TRIM(OutputFileName),'_',rank,'__',istart(1),'_',istart(2),'__',ipatch,'.tif',CHAR(0)
!                     & TRIM(OutputFileName),'_',isub,'__',istart(1),'_',istart(2),'__',ipatch,'.tif',CHAR(0)
#elif __DIME == __3D
                    WRITE(filename,'(A,A1,I0,A2,I0,A1,I0,A1,I0,A2,I0,A4,A1)') &
                    & TRIM(OutputFileName),'_',rank,'__',istart(1),'_',istart(2),'_',istart(3),'__',ipatch,'.tif',CHAR(0)
!                     & TRIM(OutputFileName),'_',isub,'__',istart(1),'_',istart(2),'_',istart(3),'__',ipatch,'.tif',CHAR(0)
#endif
                 ENDIF

                 !The TIFF file format uses 32bit offsets and, as such, is limited to 4 gigabytes
#if   __DIME == __2D
                 IF (PRESENT(idn)) THEN
                    IsBigTIFF=REAL(Nm(1),MK)*REAL(Nm(2),MK)/1024._MK.GE.4194304._MK
                 ELSE
                    IsBigTIFF=REAL(Nm(1),MK)*REAL(Nm(2),MK)*REAL(maxiter,MK)/1024._MK.GE.4194304._MK
                 ENDIF
#elif __DIME == __3D
                 IsBigTIFF=REAL(Nm(1),MK)*REAL(Nm(2),MK)*REAL(Nm(3),MK)/1024._MK.GE.4194304._MK
#endif

                 SELECT CASE (IsBigTIFF)
                 CASE (.TRUE.)
                    info=ppm_rc_open_write_bigtiff(filename)
                 CASE DEFAULT
                    info=ppm_rc_open_write_tiff(filename)
                 END SELECT
                 or_fail('Error opening tiff file to write in.')

!                  info=ppm_rc_write_tiff_header(bitsPerSampleW,1,Nm(1),Nm(2))
!                  or_fail('Error writing tiff file header.')

                 info=ppm_rc_write_tiff_header(8,1,Nm(1),Nm(2))
                 or_fail('Error writing tiff file header.')

                 SELECT CASE (FieldIn%data_type)
                 CASE (ppm_type_int)
                    CALL p%get_field(FieldIn,DTYPE(wp_i),info)
                    or_fail("Failed to get field wp_i.")

                    iopt=ppm_param_alloc_fit
                    ldu(1:2)=Nm(1:2)
                    CALL ppm_alloc(buf,ldu,iopt,info)
                    or_fail_alloc('Failed to allocate C buffer.')

                    !The white pixel value in the TIFF image
                    white=255

                    SELECT CASE (labelInit)
                    CASE (.TRUE.)
#if   __DIME == __2D
                       !Take each labeled regoin
                       DO y=1,Nm(2)
                          DO x=1,Nm(1)
                             IF (DTYPE(wp_i)(x,y).EQ.0.OR.DTYPE(wp_i)(x,y).EQ.FORBIDDEN) THEN
                                buf(x,y)=0
                             ELSE
                                white=MOD(ABS(DTYPE(wp_i)(x,y)),150)+106
                                buf(x,y)=white
                             ENDIF
                          ENDDO
                       ENDDO

                       info=ppm_rc_write_tiff_strip(buf,8,Nm(1),Nm(2),0,1)
                       or_fail("ppm_rc_write_tiff_strip")
#elif __DIME == __3D
                       DO istrip=1,Nm(3)
                          !Take each labeled regoin
                          DO y=1,Nm(2)
                             DO x=1,Nm(1)
                                IF (DTYPE(wp_i)(x,y,istrip).EQ.0.OR.DTYPE(wp_i)(x,y,istrip).EQ.FORBIDDEN) THEN
                                   buf(x,y)=0
                                ELSE
                                   white=MOD(ABS(DTYPE(wp_i)(x,y,istrip)),150)+106
                                   buf(x,y)=white
                                ENDIF
                             ENDDO
                          ENDDO

                          info=ppm_rc_write_tiff_header(istrip-1,Nm(3))
                          or_fail('Error writing tiff file header.')

                          info=ppm_rc_write_tiff_strip(buf,8,Nm(1),Nm(2),istrip-1,Nm(3))
                          or_fail("libtiff_write_tiff_strip")
                       ENDDO ! loop through slices
#endif

                    CASE (.FALSE.)
#if   __DIME == __2D
                       !Only take the boundary of each labeled regoin
                       DO y=1,Nm(2)
                          DO x=1,Nm(1)
                             IF (DTYPE(wp_i)(x,y).LT.0) THEN
!                                 white=MOD(ABS(DTYPE(wp_i)(x,y)),150)+106
                                buf(x,y)=white
                             ELSE
                                buf(x,y)=0
                             ENDIF
                          ENDDO
                       ENDDO

                       info=ppm_rc_write_tiff_strip(buf,8,Nm(1),Nm(2),0,1)
                       or_fail("ppm_rc_write_tiff_strip")
#elif __DIME == __3D
                       DO istrip=1,Nm(3)
                          !Only take the boundary of each labeled regoin
                          DO y=1,Nm(2)
                             DO x=1,Nm(1)
                                IF (DTYPE(wp_i)(x,y,istrip).LT.0) THEN
                                   white=MOD(ABS(DTYPE(wp_i)(x,y,istrip)),150)+106
                                   buf(x,y)=white
                                ELSE
                                   buf(x,y)=0
                                ENDIF
                             ENDDO
                          ENDDO

                          info=ppm_rc_write_tiff_header(istrip-1,Nm(3))
                          or_fail('Error writing tiff file header.')

                          info=ppm_rc_write_tiff_strip(buf,8,Nm(1),Nm(2),istrip-1,Nm(3))
                          or_fail("libtiff_write_tiff_strip")
                       ENDDO ! loop through slices
#endif

                    END SELECT !(labelInit)

                 CASE DEFAULT
                    fail("FieldIn data_type is not supported")

                 END SELECT !(FieldIn%data_type)

                 info=ppm_rc_close_tiff()
                 or_fail('Error closing tiff file.')

              END SELECT !SELECT TYPE

           ENDDO patch_loop

        ENDDO sub_loop ! loop through subs

        iopt=ppm_param_dealloc
        CALL ppm_alloc(buf,ldu,iopt,info)
        or_fail_dealloc('Could not deallocate C buffer.')

        NULLIFY(DTYPE(wp_i),buf)

        IF (ASSOCIATED(tmp_dinfo)) THEN
           CALL MeshOut_%field_ptr%remove(info,FieldIn)
           or_fail("MeshOut_%field_ptr%remove")

           CALL FieldIn%discr_info%remove(info,tmp_dinfo)
           or_fail("FieldIn%discr_info%remove")

           DEALLOCATE(tmp_dinfo,STAT=info)
           or_fail_dealloc("Failed to deallocate tmp_dinfo")

           NULLIFY(tmp_dinfo)
        ENDIF

        IF (PRESENT(liotopo)) THEN
           IF (liotopo.AND.ppm_nproc.GT.1) THEN
              !---------------------------------------------------------------------
              ! Destroy the IO mesh
              !---------------------------------------------------------------------
              CALL tiomesh%destroy(info)
              or_fail('Failed to destroy tiomesh.')

              dealloc_pointer("tiomesh")

              !---------------------------------------------------------------------
              ! Destroy the IO topo
              !---------------------------------------------------------------------
              topo => ppm_topo(tiotopoid)%t

              iopt=ppm_param_dealloc
              CALL ppm_alloc(topo,iopt,info)
              or_fail_dealloc("Failed to deallocate topology (iotopoid)!")
           ENDIF !liotopo.AND.ppm_nproc.GT.1
        ENDIF !PRESENT(liotopo)

        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
      9999 CONTINUE
        CALL substop(caller,t0,info)
        RETURN
      CONTAINS
        SUBROUTINE check
          IMPLICIT NONE
          check_true(<#ASSOCIATED(MeshIn)#>, &
          & "Input Mesh does not exist!",exit_point=8888)

          check_true(<#FieldIn%is_discretized_on(MeshIn)#>, &
          & "This Field is not discretized on Input Mesh!",exit_point=8888)

          IF (PRESENT(MeshOut)) THEN
             check_true(<#ASSOCIATED(MeshOut)#>, &
             & "Output Mesh does not exist!",exit_point=8888)
          ENDIF

          IF (PRESENT(MeshOut).AND.PRESENT(liotopo)) THEN
             IF (liotopo.AND.ppm_nproc.GT.1) THEN
                fail("It is not possible to have MeshOut and iomesh at the same time!", &
                & ppm_error=ppm_error_fatal,exit_point=8888)
             ENDIF
          ENDIF

        8888 CONTINUE
          RETURN
        END SUBROUTINE
      END SUBROUTINE DTYPE(ppm_rc_write_image_label)

