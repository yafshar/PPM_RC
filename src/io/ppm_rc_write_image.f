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
      !  Purpose      : Write the FieldIn data to the output TIFF file
      !
      !  Input        : FieldIn,MeshIn,OutputFileName,MeshOut,idn,Scalefac
      !
      !  Output       : info          (I) return status. 0 on success.
      !
      !  Routines     :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------

      SUBROUTINE DTYPE(ppm_rc_write_image)(FieldIn,MeshIn,OutputFileName,info, &
      &          MeshOut,idn,Scalefac)

        !-------------------------------------------------------------------------
        !  Modules
        !-------------------------------------------------------------------------
#ifdef __F2003
        USE ISO_C_BINDING, ONLY : C_INT
#endif

        USE ppm_module_topo_typedef, ONLY : ppm_t_topo,ppm_topo
        USE ppm_module_interfaces, ONLY : ppm_t_discr_info_

        USE ppm_rc_module_tiff, ONLY : ppm_rc_open_write_tiff,ppm_rc_open_write_bigtiff, &
        &   ppm_rc_close_tiff,ppm_rc_write_tiff_header,ppm_rc_write_tiff_strip, &
        &   bitsPerSampleW,ppm_rc_write_tiff_scanline
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
        !!! Input Field

        CLASS(ppm_t_equi_mesh_),           POINTER       :: MeshIn
        !!! Input Filed is discretized on this Mesh.

        CHARACTER(LEN=*),                  INTENT(IN   ) :: OutputFileName
        !!! Output file name for writing the data in.

        INTEGER,                           INTENT(  OUT) :: info

        CLASS(ppm_t_equi_mesh_), OPTIONAL, POINTER       :: MeshOut
        !!! Output Mesh which we want to write the data in
        !!! MeshOut can be different from MeshIn, As we want to work on
        !!! one topology and the corresponding mesh and want to write the
        !!! data on a different topology and its coressponding Mesh

        INTEGER,                 OPTIONAL, INTENT(IN   ) :: idn
        !!! ID number of the file.
        !!! It can be iteration index.

        REAL(MK),                OPTIONAL, INTENT(IN   ) :: Scalefac
        !!! Scaling factor for the output data
        !!! We can write the normalized image in a 8,16 or 32 BITs format
        !!! For writing the normal image in 8 BITs format Scalefac should be
        !!! 255
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
        TYPE(ppm_t_topo), POINTER :: topo

        CLASS(ppm_t_discr_info_), POINTER :: tmp_dinfo

        CLASS(ppm_t_equi_mesh_), POINTER :: MeshOut_

#if   __DIME == __2D
        REAL(ppm_kind_single), CONTIGUOUS, DIMENSION(:,:),   POINTER :: wprs
        REAL(MK),              CONTIGUOUS, DIMENSION(:,:),   POINTER :: wpr
#elif __DIME == __3D
        REAL(ppm_kind_single), CONTIGUOUS, DIMENSION(:,:,:), POINTER :: wprs
        REAL(MK),              CONTIGUOUS, DIMENSION(:,:,:), POINTER :: wpr
#endif
        REAL(ppm_kind_single),             DIMENSION(:,:),   POINTER :: wprsp
        REAL(MK),                          DIMENSION(:,:),   POINTER :: wprp
        REAL(ppm_kind_single), CONTIGUOUS, DIMENSION(:,:),   POINTER :: bufrs
        REAL(MK),              CONTIGUOUS, DIMENSION(:,:),   POINTER :: bufr
        REAL(MK)                                                     :: Scalefac_
        REAL(MK)                                                     :: MinShiftVal_
        REAL(ppm_kind_double) :: t0
        REAL(ppm_kind_single) :: Scalefacs
        REAL(ppm_kind_single) :: MinShiftVals

#if   __DIME == __2D
        INTEGER, CONTIGUOUS, DIMENSION(:,:),   POINTER :: wpi
#elif __DIME == __3D
        INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER :: wpi
#endif
        INTEGER, CONTIGUOUS, DIMENSION(:,:),   POINTER :: wpip
        INTEGER, CONTIGUOUS, DIMENSION(:,:),   POINTER :: bufi
        INTEGER,             DIMENSION(:),     POINTER :: Nm
        INTEGER,             DIMENSION(:),     POINTER :: istart
        ! Lower-left coordinates on the global mesh
        INTEGER                                        :: isub,x,y,i,ipatch
        INTEGER                                        :: istrip
        INTEGER                                        :: ldu(2),iopt
#ifdef __F2003
        INTEGER(C_INT)                                 :: bitsPerSampleW_
#else
        INTEGER                                        :: bitsPerSampleW_
#endif

        CHARACTER(LEN=ppm_char) :: filename
        CHARACTER(LEN=ppm_char) :: caller = 'ppm_rc_write_image'

        LOGICAL :: IsBigTIFF
        !-------------------------------------------------------------------------
        !  Initialize
        !-------------------------------------------------------------------------
        CALL substart(caller,t0,info)

        CALL check()

        IF (PRESENT(Scalefac)) THEN
           IF (Scalefac.GT.two) THEN
              Scalefac_=Scalefac
              MinShiftVal_=zero
              IF      (Scalefac_.LE.255._MK) THEN
                 bitsPerSampleW_=8
              ELSE IF (Scalefac_.LE.65535._MK) THEN
                 bitsPerSampleW_=16
              ELSE
                 bitsPerSampleW_=32
              ENDIF
           ELSE IF (Scalefac.GT.zero) THEN
              Scalefac_=Scalefac
              MinShiftVal_=zero
              bitsPerSampleW_=bitsPerSampleW
           ELSE
              Scalefac_=ImageNormalfac
              MinShiftVal_=MinShiftVal
              bitsPerSampleW_=bitsPerSampleW
           ENDIF
        ELSE
           Scalefac_=ImageNormalfac
           MinShiftVal_=MinShiftVal
           bitsPerSampleW_=bitsPerSampleW
        ENDIF

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

        NULLIFY(wpr,wprs,wpi,bufi,bufrs,bufr)

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
#elif __DIME == __3D
                    WRITE(filename,'(A,A1,I0,A2,I0,A1,I0,A1,I0,A2,I0,A1,I0,A4,A1)') &
                    & TRIM(OutputFileName),'_',rank,'__',istart(1),'_',istart(2),'_',istart(3),'__',ipatch,'_',idn,'.tif',CHAR(0)
#endif
                 ELSE
#if   __DIME == __2D
                    WRITE(filename,'(A,A1,I0,A2,I0,A1,I0,A2,I0,A4,A1)') &
                    & TRIM(OutputFileName),'_',rank,'__',istart(1),'_',istart(2),'__',ipatch,'.tif',CHAR(0)
#elif __DIME == __3D
                    WRITE(filename,'(A,A1,I0,A2,I0,A1,I0,A1,I0,A2,I0,A4,A1)') &
                    & TRIM(OutputFileName),'_',rank,'__',istart(1),'_',istart(2),'_',istart(3),'__',ipatch,'.tif',CHAR(0)
#endif
                 ENDIF

                 !The TIFF file format uses 32bit offsets and, as such, is limited to 4 gigabytes
#if   __DIME == __2D
                 IsBigTIFF=REAL(Nm(1),MK)*REAL(Nm(2),MK)*REAL(bitsPerSampleW_/8,MK)/1024._MK.GE.4194304._MK
#elif __DIME == __3D
                 IsBigTIFF=REAL(Nm(1),MK)*REAL(Nm(2),MK)*REAL(Nm(3),MK)*REAL(bitsPerSampleW_/8,MK)/1024._MK.GE.4194304._MK
#endif

                 SELECT CASE (IsBigTIFF)
                 CASE (.TRUE.)
                    info=ppm_rc_open_write_bigtiff(filename)
                 CASE DEFAULT
                    info=ppm_rc_open_write_tiff(filename)
                 END SELECT !(IsBigTIFF)
                 or_fail('Error opening tiff file to write in.')

                 info=ppm_rc_write_tiff_header(bitsPerSampleW_,1,Nm(1),Nm(2))
                 or_fail('Error writing tiff file header.')

                 SELECT CASE (FieldIn%data_type)
                 CASE (ppm_type_int)
                    CALL p%get_field(FieldIn,wpi,info)
                    or_fail("Failed to get field wpi.")

#if   __DIME == __2D
                    wpip => wpi(1:Nm(1),1:Nm(2))

                    info=ppm_rc_write_tiff_strip(wpip,bitsPerSampleW_,Nm(1),Nm(2),0,1)
                    or_fail("ppm_rc_write_tiff_strip")
#elif __DIME == __3D
                    DO istrip=1,Nm(3)
                       info=ppm_rc_write_tiff_header(istrip-1,Nm(3))
                       or_fail('Error writing tiff file header.')

                       wpip => wpi(1:Nm(1),1:Nm(2),istrip)

                       info=ppm_rc_write_tiff_strip(wpip,bitsPerSampleW_,Nm(1),Nm(2),istrip-1,Nm(3))
                       or_fail("ppm_rc_write_tiff_strip")
                    ENDDO ! loop through slices
#endif

                 CASE (ppm_type_real)
                 !!! For the real data type we do check whether the data is been normalized
                 !!! or not, also in the case of real data, it would be the image dada itself
                 !!! so we do not need to correct the range of data
                    CALL p%get_field(FieldIn,wpr,info)
                    or_fail("Failed to get field wpr.")

                    IF (lNormalize) THEN
                       iopt=ppm_param_alloc_fit
                       SELECT CASE (bitsPerSampleW_)
                       CASE (32)
                          CALL ppm_alloc(bufr,Nm,iopt,info)
                       CASE DEFAULT
                          CALL ppm_alloc(bufi,Nm,iopt,info)
                       END SELECT
                       or_fail_alloc('Failed to allocate C buffer.')
                    ENDIF !(lNormalize)

#if   __DIME == __2D
                    SELECT CASE (lNormalize)
                    CASE (.TRUE.)
                       SELECT CASE (bitsPerSampleW_)
                       CASE (32)
                          IF (MinShiftVal_.GT.zero.OR.MinShiftVal_.LT.zero) THEN
                             FORALL (x=1:Nm(1),y=1:Nm(2))
                                bufr(x,y)=wpr(x,y)*Scalefac_+MinShiftVal_
                             END FORALL
                          ELSE
                             FORALL (x=1:Nm(1),y=1:Nm(2))
                                bufr(x,y)=wpr(x,y)*Scalefac_
                             END FORALL
                          ENDIF

                          info=ppm_rc_write_tiff_strip(bufr,bitsPerSampleW_,Nm(1),Nm(2),0,1)
                          or_fail("ppm_rc_write_tiff_strip")

                       CASE DEFAULT
                          IF (MinShiftVal_.GT.zero.OR.MinShiftVal_.LT.zero) THEN
                             FORALL (x=1:Nm(1),y=1:Nm(2))
                                bufi(x,y)=INT(wpr(x,y)*Scalefac_+MinShiftVal_)
                             END FORALL
                          ELSE
                             FORALL (x=1:Nm(1),y=1:Nm(2))
                                bufi(x,y)=INT(wpr(x,y)*Scalefac_)
                             END FORALL
                          ENDIF

                          info=ppm_rc_write_tiff_strip(bufi,bitsPerSampleW_,Nm(1),Nm(2),0,1)
                          or_fail("ppm_rc_write_tiff_strip")

                       END SELECT !(bitsPerSampleW_)

                    CASE (.FALSE.)
                       wprp => wpr(1:Nm(1),1:Nm(2))

                       info=ppm_rc_write_tiff_strip(wprp,bitsPerSampleW_,Nm(1),Nm(2),0,1)
                       or_fail("ppm_rc_write_tiff_strip")

                    END SELECT !(lNormalize)
#elif __DIME == __3D
                    DO istrip=1,Nm(3)
                       info=ppm_rc_write_tiff_header(istrip-1,Nm(3))
                       or_fail('Error writing tiff file header.')

                       SELECT CASE (lNormalize)
                       CASE (.TRUE.)
                          SELECT CASE (bitsPerSampleW_)
                          CASE (32)
                             IF (MinShiftVal_.GT.zero.OR.MinShiftVal_.LT.zero) THEN
                                FORALL (x=1:Nm(1),y=1:Nm(2))
                                   bufr(x,y)=wpr(x,y,istrip)*Scalefac_+MinShiftVal_
                                END FORALL
                             ELSE
                                FORALL (x=1:Nm(1),y=1:Nm(2))
                                   bufr(x,y)=wpr(x,y,istrip)*Scalefac_
                                END FORALL
                             ENDIF

                             info=ppm_rc_write_tiff_strip(bufr,bitsPerSampleW_,Nm(1),Nm(2),istrip-1,Nm(3))
                             or_fail("ppm_rc_write_tiff_strip")
!                              info=ppm_rc_write_tiff_scanline(bufr,0,bitsPerSampleW_,0,Nm(1),Nm(2),istrip-1,Nm(3))
!                              or_fail("ppm_rc_write_tiff_scanline")

                          CASE DEFAULT
                             IF (MinShiftVal_.GT.zero.OR.MinShiftVal_.LT.zero) THEN
                                FORALL (x=1:Nm(1),y=1:Nm(2))
                                   bufi(x,y)=INT(wpr(x,y,istrip)*Scalefac_+MinShiftVal_)
                                END FORALL
                             ELSE
                                FORALL (x=1:Nm(1),y=1:Nm(2))
                                   bufi(x,y)=INT(wpr(x,y,istrip)*Scalefac_)
                                END FORALL
                             ENDIF

                             info=ppm_rc_write_tiff_strip(bufi,bitsPerSampleW_,Nm(1),Nm(2),istrip-1,Nm(3))
                             or_fail("ppm_rc_write_tiff_strip")

                          END SELECT !(bitsPerSampleW_)

                       CASE (.FALSE.)
                          wprp => wpr(1:Nm(1),1:Nm(2),istrip)

                          info=ppm_rc_write_tiff_strip(wprp,bitsPerSampleW_,Nm(1),Nm(2),istrip-1,Nm(3))
                          or_fail("ppm_rc_write_tiff_strip")

                       END SELECT !(lNormalize)
                    ENDDO ! loop through slices
#endif

                 !!! For the real data type we do check whether the data is been normalized
                 !!! or not, also in the case of real data, it would be the image dada itself
                 !!! so we do not need to correct the range of data
                 CASE (ppm_type_real_single)
                    CALL p%get_field(FieldIn,wprs,info)
                    or_fail("Failed to get field wprs.")

                    IF (lNormalize) THEN
                       iopt=ppm_param_alloc_fit
                       SELECT CASE (bitsPerSampleW_)
                       CASE (32)
                          CALL ppm_alloc(bufrs,Nm,iopt,info)
                       CASE DEFAULT
                          CALL ppm_alloc(bufi,Nm,iopt,info)
                       END SELECT
                       or_fail_alloc('Failed to allocate C buffer.')

                       Scalefacs=REAL(Scalefac_,ppm_kind_single)
                       MinShiftVals=REAL(MinShiftVal_,ppm_kind_single)
                    ENDIF !(lNormalize)

#if   __DIME == __2D
                    SELECT CASE (lNormalize)
                    CASE (.TRUE.)
                       SELECT CASE (bitsPerSampleW_)
                       CASE (32)
                          IF (MinShiftVal_.GT.zero.OR.MinShiftVal_.LT.zero) THEN
                             FORALL (x=1:Nm(1),y=1:Nm(2))
                                bufrs(x,y)=wprs(x,y)*Scalefacs+MinShiftVals
                             END FORALL
                          ELSE
                             FORALL (x=1:Nm(1),y=1:Nm(2))
                                bufrs(x,y)=wprs(x,y)*Scalefacs
                             END FORALL
                          ENDIF

                          info=ppm_rc_write_tiff_strip(bufrs,bitsPerSampleW_,Nm(1),Nm(2),0,1)
                          or_fail("ppm_rc_write_tiff_strip")

                       CASE DEFAULT
                          IF (MinShiftVal_.GT.zero.OR.MinShiftVal_.LT.zero) THEN
                             FORALL (x=1:Nm(1),y=1:Nm(2))
                                bufi(x,y)=INT(wprs(x,y)*Scalefacs+MinShiftVals)
                             END FORALL
                          ELSE
                             FORALL (x=1:Nm(1),y=1:Nm(2))
                                bufi(x,y)=INT(wprs(x,y)*Scalefacs)
                             END FORALL
                          ENDIF

                          info=ppm_rc_write_tiff_strip(bufi,bitsPerSampleW_,Nm(1),Nm(2),0,1)
                          or_fail("ppm_rc_write_tiff_strip")

                       END SELECT !(bitsPerSampleW_)

                    CASE (.FALSE.)
                       wprsp => wprs(1:Nm(1),1:Nm(2))

                       info=ppm_rc_write_tiff_strip(wprsp,bitsPerSampleW_,Nm(1),Nm(2),0,1)
                       or_fail("ppm_rc_write_tiff_strip")

                    END SELECT !(lNormalize)
#elif __DIME == __3D
                    DO istrip=1,Nm(3)
                       info=ppm_rc_write_tiff_header(istrip-1,Nm(3))
                       or_fail('Error writing tiff file header.')

                       SELECT CASE (lNormalize)
                       CASE (.TRUE.)
                          SELECT CASE (bitsPerSampleW_)
                          CASE (32)
                             IF (MinShiftVal_.GT.zero.OR.MinShiftVal_.LT.zero) THEN
                                FORALL (x=1:Nm(1),y=1:Nm(2))
                                   bufrs(x,y)=wprs(x,y,istrip)*Scalefacs+MinShiftVals
                                END FORALL
                             ELSE
                                FORALL (x=1:Nm(1),y=1:Nm(2))
                                   bufrs(x,y)=wprs(x,y,istrip)*Scalefacs
                                END FORALL
                             ENDIF

                             info=ppm_rc_write_tiff_strip(bufrs,bitsPerSampleW_,Nm(1),Nm(2),istrip-1,Nm(3))
                             or_fail("ppm_rc_write_tiff_strip")
!                              info=ppm_rc_write_tiff_scanline(bufr,0,bitsPerSampleW_,0,Nm(1),Nm(2),istrip-1,Nm(3))
!                              or_fail("ppm_rc_write_tiff_scanline")

                          CASE DEFAULT
                             IF (MinShiftVal_.GT.zero.OR.MinShiftVal_.LT.zero) THEN
                                FORALL (x=1:Nm(1),y=1:Nm(2))
                                   bufi(x,y)=INT(wprs(x,y,istrip)*Scalefacs+MinShiftVals)
                                END FORALL
                             ELSE
                                FORALL (x=1:Nm(1),y=1:Nm(2))
                                   bufi(x,y)=INT(wprs(x,y,istrip)*Scalefacs)
                                END FORALL
                             ENDIF

                             info=ppm_rc_write_tiff_strip(bufi,bitsPerSampleW_,Nm(1),Nm(2),istrip-1,Nm(3))
                             or_fail("ppm_rc_write_tiff_strip")

                          END SELECT !(bitsPerSampleW_)

                       CASE (.FALSE.)
                          wprsp => wprs(1:Nm(1),1:Nm(2),istrip)

                          info=ppm_rc_write_tiff_strip(wprsp,bitsPerSampleW_,Nm(1),Nm(2),istrip-1,Nm(3))
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

        NULLIFY(wpr,wprs,wpi,bufi,bufrs,bufr)

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
          check_true(<#ASSOCIATED(MeshIn)#>,"Input Mesh does not exist!",exit_point=8888)

          check_true(<#FieldIn%is_discretized_on(MeshIn)#>, &
          & "This Field is not discretized on Input Mesh!",exit_point=8888)

          IF (PRESENT(MeshOut)) THEN
             check_true(<#ASSOCIATED(MeshOut)#>,"Output Mesh does not exist!",exit_point=8888)
          ENDIF
          IF (PRESENT(Scalefac)) THEN
             IF (Scalefac.GT.zero) THEN
                IF (FieldIn%data_type.EQ.ppm_type_int) THEN
                   fail("No support to convert one type to the other! Scalefac can not be > 0 for integer data type!", &
                   & ppm_error=ppm_error_fatal,exit_point=8888)
                ENDIF
             ENDIF
          ENDIF
        8888 CONTINUE
          RETURN
        END SUBROUTINE
      END SUBROUTINE DTYPE(ppm_rc_write_image)


      SUBROUTINE DTYPE(ppm_rc_write_image_label)(FieldIn,MeshIn, &
      &          OutputFileName,info,MeshOut,idn,liotopo,        &
      &          bitsPerSampleW,lcast)

        !-------------------------------------------------------------------------
        !  Modules
        !-------------------------------------------------------------------------
#ifdef __F2003
        USE ISO_C_BINDING, ONLY : C_INT
#endif

        USE ppm_module_data, ONLY : ppm_param_assign_internal, &
        & ppm_param_decomp_xpencil,ppm_param_decomp_xy_slab
        USE ppm_module_mktopo, ONLY : ppm_mktopo
        USE ppm_module_topo_typedef, ONLY : ppm_t_topo,ppm_topo
        USE ppm_module_interfaces, ONLY : ppm_t_discr_info_

        USE ppm_rc_module_tiff, ONLY : ppm_rc_open_write_tiff,ppm_rc_open_write_bigtiff, &
        &   ppm_rc_close_tiff,ppm_rc_write_tiff_header,ppm_rc_write_tiff_strip
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

        INTEGER,                 OPTIONAL, INTENT(IN   ) :: bitsPerSampleW

        LOGICAL,                 OPTIONAL, INTENT(IN   ) :: lcast
        !!! Cast the labels to the output bitsPerSample and write the results
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
        TYPE(ppm_t_topo), POINTER :: topo

        CLASS(ppm_t_discr_info_), POINTER :: tmp_dinfo

        CLASS(ppm_t_equi_mesh_), POINTER :: MeshOut_
        CLASS(ppm_t_equi_mesh_), POINTER :: tiomesh
        !!! mesh for I/O


        REAL(MK), CONTIGUOUS, DIMENSION(:,:), POINTER :: xp
        REAL(ppm_kind_double)                         :: t0

#if   __DIME == __2D
        INTEGER, CONTIGUOUS, DIMENSION(:,:),   POINTER :: wpi
#elif __DIME == __3D
        INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER :: wpi
#endif
        INTEGER, CONTIGUOUS, DIMENSION(:,:),   POINTER :: buf
        INTEGER,             DIMENSION(:),     POINTER :: Nm
        INTEGER,             DIMENSION(:),     POINTER :: istart
        INTEGER                                        :: isub,x,y,i,ipatch
        INTEGER                                        :: white
        INTEGER                                        :: mwhite
        INTEGER                                        :: rwhite
#if   __DIME == __3D
        INTEGER                                        :: istrip
#endif
        INTEGER                                        :: ldu(2),iopt
        INTEGER                                        :: nproc,ld(__DIME)
        INTEGER                                        :: tiotopoid,tiomeshid
#ifdef __F2003
        INTEGER(C_INT)                                 :: bitsPerSampleW_
#else
        INTEGER                                        :: bitsPerSampleW_
#endif

        CHARACTER(LEN=ppm_char) :: filename
        CHARACTER(LEN=ppm_char) :: caller = 'ppm_rc_write_image_label'

        LOGICAL :: IsBigTIFF
        LOGICAL :: labelInit
        !-------------------------------------------------------------------------
        !  Initialize
        !-------------------------------------------------------------------------
        CALL substart(caller,t0,info)

        CALL check()

        bitsPerSampleW_=MERGE(bitsPerSampleW,8,PRESENT(bitsPerSampleW))

        !Write the label image file for using as an initialization input
        labelInit=.FALSE.

        NULLIFY(MeshOut_)

        IF (PRESENT(liotopo)) THEN
           IF (liotopo.AND.ppm_nproc.GT.1) THEN
              !-------------------------------------------------------------------------
              !  Create an IO topology for writing the image file
              !-------------------------------------------------------------------------
              nproc=ppm_nproc

              IF (n_procs_write.GT.0.AND.n_procs_write.LT.ppm_nproc) THEN
                 ppm_nproc=n_procs_write
                 IF (Ngrid(ppm_rc_dim).LE.ppm_nproc) THEN
                    ppm_nproc=Ngrid(ppm_rc_dim)-1
                 ENDIF
              ELSE IF (n_procs_write.GE.ppm_nproc) THEN
                 IF (Ngrid(ppm_rc_dim).LE.ppm_nproc) THEN
                    ppm_nproc=Ngrid(ppm_rc_dim)-1
                 ENDIF
              ELSE
                 IF (ppm_nproc.GE.Ngrid(ppm_rc_dim)) THEN
                    ppm_nproc=Ngrid(ppm_rc_dim)-1
                 ENDIF
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

              IF (nproc.GT.1) THEN
                 IF (ppm_nproc.EQ.1) THEN
                    IF (rank.GT.0) THEN
                       topo => ppm_topo(tiotopoid)%t
                       topo%nsublist=0
                       topo%sub2proc=0
                    ENDIF
                 ENDIF
              ENDIF

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

        NULLIFY(wpi,buf)

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
#elif __DIME == __3D
                    WRITE(filename,'(A,A1,I0,A2,I0,A1,I0,A1,I0,A2,I0,A1,I0,A4,A1)') &
                    & TRIM(OutputFileName),'_',rank,'__',istart(1),'_',istart(2),'_',istart(3),'__',ipatch,'_',idn,'.tif',CHAR(0)
#endif
                 ELSE
#if   __DIME == __2D
                    WRITE(filename,'(A,A1,I0,A2,I0,A1,I0,A2,I0,A4,A1)') &
                    & TRIM(OutputFileName),'_',rank,'__',istart(1),'_',istart(2),'__',ipatch,'.tif',CHAR(0)
#elif __DIME == __3D
                    WRITE(filename,'(A,A1,I0,A2,I0,A1,I0,A1,I0,A2,I0,A4,A1)') &
                    & TRIM(OutputFileName),'_',rank,'__',istart(1),'_',istart(2),'_',istart(3),'__',ipatch,'.tif',CHAR(0)
#endif
                 ENDIF

                 !The TIFF file format uses 32bit offsets and, as such, is limited to 4 gigabytes
#if   __DIME == __2D
                 IF (PRESENT(idn)) THEN
                    IsBigTIFF=REAL(Nm(1),MK)*REAL(Nm(2),MK)*REAL(bitsPerSampleW_/8,MK)/1024._MK.GE.4194304._MK
                 ELSE
                    IF (maxiter.GT.0) THEN
                       IsBigTIFF=REAL(Nm(1),MK)*REAL(Nm(2),MK)*REAL(maxiter,MK)*REAL(bitsPerSampleW_/8,MK)/1024._MK.GE.4194304._MK
                    ELSE
                       IsBigTIFF=REAL(Nm(1),MK)*REAL(Nm(2),MK)*REAL(bitsPerSampleW_/8,MK)/1024._MK.GE.4194304._MK
                    ENDIF
                 ENDIF
#elif __DIME == __3D
                 IsBigTIFF=REAL(Nm(1),MK)*REAL(Nm(2),MK)*REAL(Nm(3),MK)*REAL(bitsPerSampleW_/8,MK)/1024._MK.GE.4194304._MK
#endif

                 SELECT CASE (IsBigTIFF)
                 CASE (.TRUE.)
                    info=ppm_rc_open_write_bigtiff(filename)
                 CASE DEFAULT
                    info=ppm_rc_open_write_tiff(filename)
                 END SELECT
                 or_fail('Error opening tiff file to write in.')

                 info=ppm_rc_write_tiff_header(bitsPerSampleW_,1,Nm(1),Nm(2))
                 or_fail('Error writing tiff file header.')

                 SELECT CASE (FieldIn%data_type)
                 CASE (ppm_type_int)
                    CALL p%get_field(FieldIn,wpi,info)
                    or_fail("Failed to get field wpi.")

                    IF (PRESENT(lcast)) THEN
                       !The white pixel value in the TIFF image
                       SELECT CASE (bitsPerSampleW_)
                       CASE (8)
                          iopt=ppm_param_alloc_fit
                          CALL ppm_alloc(buf,Nm,iopt,info)
                          or_fail_alloc('Failed to allocate C buffer.')
#if   __DIME == __2D
                          !Take each labeled regoin
                          DO y=1,Nm(2)
                             DO x=1,Nm(1)
                                IF (wpi(x,y).EQ.0.OR.wpi(x,y).EQ.FORBIDDEN) THEN
                                   buf(x,y)=0
                                ELSE
                                   IF (ABS(wpi(x,y)).GT.255) THEN
                                      buf(x,y)=255
                                   ELSE
                                      buf(x,y)=ABS(wpi(x,y))
                                   ENDIF
                                ENDIF
                             ENDDO
                          ENDDO

                          info=ppm_rc_write_tiff_strip(buf,bitsPerSampleW_,Nm(1),Nm(2),0,1)
                          or_fail("ppm_rc_write_tiff_strip")
#elif __DIME == __3D
                          DO istrip=1,Nm(3)
                             !Take each labeled regoin
                             DO y=1,Nm(2)
                                DO x=1,Nm(1)
                                   IF (wpi(x,y,istrip).EQ.0.OR.wpi(x,y,istrip).EQ.FORBIDDEN) THEN
                                      buf(x,y)=0
                                   ELSE
                                      IF (ABS(wpi(x,y,istrip)).GT.255) THEN
                                         buf(x,y)=255
                                      ELSE
                                         buf(x,y)=ABS(wpi(x,y,istrip))
                                      ENDIF
                                   ENDIF
                                ENDDO
                             ENDDO

                             info=ppm_rc_write_tiff_header(istrip-1,Nm(3))
                             or_fail('Error writing tiff file header.')

                             info=ppm_rc_write_tiff_strip(buf,bitsPerSampleW_,Nm(1),Nm(2),istrip-1,Nm(3))
                             or_fail("libtiff_write_tiff_strip")
                          ENDDO ! loop through slices
#endif

                          iopt=ppm_param_dealloc
                          CALL ppm_alloc(buf,ldu,iopt,info)
                          or_fail_dealloc('Could not deallocate C buffer.')
                       CASE (16)
                          iopt=ppm_param_alloc_fit
                          CALL ppm_alloc(buf,Nm,iopt,info)
                          or_fail_alloc('Failed to allocate C buffer.')

#if   __DIME == __2D
                          !Take each labeled regoin
                          DO y=1,Nm(2)
                             DO x=1,Nm(1)
                                IF (wpi(x,y).EQ.0.OR.wpi(x,y).EQ.FORBIDDEN) THEN
                                   buf(x,y)=0
                                ELSE
                                   IF (ABS(wpi(x,y)).GT.65535) THEN
                                      buf(x,y)=65535
                                   ELSE
                                      buf(x,y)=ABS(wpi(x,y))
                                   ENDIF
                                ENDIF
                             ENDDO
                          ENDDO

                          info=ppm_rc_write_tiff_strip(buf,bitsPerSampleW_,Nm(1),Nm(2),0,1)
                          or_fail("ppm_rc_write_tiff_strip")
#elif __DIME == __3D
                          DO istrip=1,Nm(3)
                             !Take each labeled regoin
                             DO y=1,Nm(2)
                                DO x=1,Nm(1)
                                   IF (wpi(x,y,istrip).EQ.0.OR.wpi(x,y,istrip).EQ.FORBIDDEN) THEN
                                      buf(x,y)=0
                                   ELSE
                                      IF (ABS(wpi(x,y,istrip)).GT.65535) THEN
                                         buf(x,y)=65535
                                      ELSE
                                         buf(x,y)=ABS(wpi(x,y,istrip))
                                      ENDIF
                                   ENDIF
                                ENDDO
                             ENDDO

                             info=ppm_rc_write_tiff_header(istrip-1,Nm(3))
                             or_fail('Error writing tiff file header.')

                             info=ppm_rc_write_tiff_strip(buf,bitsPerSampleW_,Nm(1),Nm(2),istrip-1,Nm(3))
                             or_fail("libtiff_write_tiff_strip")
                          ENDDO ! loop through slices
#endif


                          iopt=ppm_param_dealloc
                          CALL ppm_alloc(buf,ldu,iopt,info)
                          or_fail_dealloc('Could not deallocate C buffer.')
                       CASE (32)
                          iopt=ppm_param_alloc_fit
                          CALL ppm_alloc(xp,Nm,iopt,info)
                          or_fail_alloc('Failed to allocate C buffer.')


#if   __DIME == __2D
                          !Take each labeled regoin
                          DO y=1,Nm(2)
                             DO x=1,Nm(1)
                                IF (wpi(x,y).EQ.0.OR.wpi(x,y).EQ.FORBIDDEN) THEN
                                   xp(x,y)=0.0_MK
                                ELSE
                                   xp(x,y)=REAL(ABS(wpi(x,y)),MK)
                                ENDIF
                             ENDDO
                          ENDDO

                          info=ppm_rc_write_tiff_strip(xp,bitsPerSampleW_,Nm(1),Nm(2),0,1)
                          or_fail("ppm_rc_write_tiff_strip")
#elif __DIME == __3D
                          DO istrip=1,Nm(3)
                             !Take each labeled regoin
                             DO y=1,Nm(2)
                                DO x=1,Nm(1)
                                   IF (wpi(x,y,istrip).EQ.0.OR.wpi(x,y,istrip).EQ.FORBIDDEN) THEN
                                      xp(x,y)=0.0_MK
                                   ELSE
                                      xp(x,y)=REAL(ABS(wpi(x,y,istrip)),MK)
                                   ENDIF
                                ENDDO
                             ENDDO

                             info=ppm_rc_write_tiff_header(istrip-1,Nm(3))
                             or_fail('Error writing tiff file header.')

                             info=ppm_rc_write_tiff_strip(xp,bitsPerSampleW_,Nm(1),Nm(2),istrip-1,Nm(3))
                             or_fail("libtiff_write_tiff_strip")
                          ENDDO ! loop through slices
#endif

                          iopt=ppm_param_dealloc
                          CALL ppm_alloc(xp,ldu,iopt,info)
                          or_fail_dealloc('Could not deallocate C buffer.')
                       END SELECT
                    ELSE !(PRESENT(lcast))
                       iopt=ppm_param_alloc_fit
                       CALL ppm_alloc(buf,Nm,iopt,info)
                       or_fail_alloc('Failed to allocate C buffer.')

                       !The white pixel value in the TIFF image
                       white=MERGE(255,65535,bitsPerSampleW_.EQ.8)
                       mwhite=MERGE(205,60035,bitsPerSampleW_.EQ.8)
                       rwhite=white-mwhite+1

                       SELECT CASE (labelInit)
                       CASE (.TRUE.)
#if   __DIME == __2D
                          !Take each labeled regoin
                          DO y=1,Nm(2)
                             DO x=1,Nm(1)
                                IF (wpi(x,y).EQ.0.OR.wpi(x,y).EQ.FORBIDDEN) THEN
                                   buf(x,y)=0
                                ELSE
                                   white=MOD(ABS(wpi(x,y)),mwhite)+rwhite
                                   buf(x,y)=white
                                ENDIF
                             ENDDO
                          ENDDO

                          info=ppm_rc_write_tiff_strip(buf,bitsPerSampleW_,Nm(1),Nm(2),0,1)
                          or_fail("ppm_rc_write_tiff_strip")
#elif __DIME == __3D
                          DO istrip=1,Nm(3)
                             !Take each labeled regoin
                             DO y=1,Nm(2)
                                DO x=1,Nm(1)
                                   IF (wpi(x,y,istrip).EQ.0.OR.wpi(x,y,istrip).EQ.FORBIDDEN) THEN
                                      buf(x,y)=0
                                   ELSE
                                      white=MOD(ABS(wpi(x,y,istrip)),mwhite)+rwhite
                                      buf(x,y)=white
                                   ENDIF
                                ENDDO
                             ENDDO

                             info=ppm_rc_write_tiff_header(istrip-1,Nm(3))
                             or_fail('Error writing tiff file header.')

                             info=ppm_rc_write_tiff_strip(buf,bitsPerSampleW_,Nm(1),Nm(2),istrip-1,Nm(3))
                             or_fail("libtiff_write_tiff_strip")
                          ENDDO ! loop through slices
#endif

                       CASE (.FALSE.)
#if   __DIME == __2D
                          !Only take the boundary of each labeled regoin
                          DO y=1,Nm(2)
                             DO x=1,Nm(1)
                                IF (wpi(x,y).LT.0) THEN
!                                  white=MOD(ABS(wpi(x,y)),mwhite)+rwhite
                                   buf(x,y)=white
                                ELSE
                                   buf(x,y)=0
                                ENDIF
                             ENDDO
                          ENDDO

                          info=ppm_rc_write_tiff_strip(buf,bitsPerSampleW_,Nm(1),Nm(2),0,1)
                          or_fail("ppm_rc_write_tiff_strip")
#elif __DIME == __3D
                          DO istrip=1,Nm(3)
                             !Only take the boundary of each labeled regoin
                             DO y=1,Nm(2)
                                DO x=1,Nm(1)
                                   IF (wpi(x,y,istrip).LT.0) THEN
   !                                    white=MOD(ABS(wpi(x,y,istrip)),mwhite)+rwhite
                                      buf(x,y)=white
                                   ELSE
                                      buf(x,y)=0
                                   ENDIF
                                ENDDO
                             ENDDO

                             info=ppm_rc_write_tiff_header(istrip-1,Nm(3))
                             or_fail('Error writing tiff file header.')

                             info=ppm_rc_write_tiff_strip(buf,bitsPerSampleW_,Nm(1),Nm(2),istrip-1,Nm(3))
                             or_fail("libtiff_write_tiff_strip")
                          ENDDO ! loop through slices
#endif

                       END SELECT !(labelInit)

                       iopt=ppm_param_dealloc
                       CALL ppm_alloc(buf,ldu,iopt,info)
                       or_fail_dealloc('Could not deallocate C buffer.')

                    ENDIF !(PRESENT(lcast))

                 CASE DEFAULT
                    fail("FieldIn data_type is not supported")

                 END SELECT !(FieldIn%data_type)

                 info=ppm_rc_close_tiff()
                 or_fail('Error closing tiff file.')

              END SELECT !SELECT TYPE

           ENDDO patch_loop

        ENDDO sub_loop ! loop through subs

        NULLIFY(wpi,buf)

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

