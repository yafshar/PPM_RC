      !-------------------------------------------------------------------------
      !  Subroutine   :                    ppm_rc_read_write_unique_label
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
      !  Author           - y.afshar           May    2016
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Subroutine   :                    ppm_rc_read_write_unique_label
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Read a bunch of TIFF slabs using Libtiff library from one
      !                 or several images.
      !                 Fill up the 3D mesh with intensity values obtained from
      !                 the image labels and write the uniquely labels resulted
      !                 read into OutputFileName
      !
      !  Input        : root name of the files to be read.
      !
      !  Output       : info          (I) return status. 0 on success.
      !
      !  Routines     : EXTERNAL C function call
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      SUBROUTINE DTYPE(ppm_rc_read_write_unique_label)(FieldIn, &
      &          MeshIn,InputFileName,OutputFileName,info)

        !-------------------------------------------------------------------------
        !  Modules
        !-------------------------------------------------------------------------
        USE ppm_module_field_typedef, ONLY : ppm_t_field_
        USE ppm_module_mesh_typedef, ONLY : ppm_t_subpatch_

        USE ppm_rc_module_tiff, ONLY : ppm_rc_open_tiff,ppm_rc_open_bigtiff,   &
        &   ppm_rc_read_tiff_scanline,ppm_rc_close_tiff,ppm_rc_read_tiff_info, &
        &   bitsPerSample
        USE ppm_rc_module_util, ONLY : ppm_rc_label_exist
        USE ppm_rc_module_fire, ONLY : ppm_rc_floodFill
        IMPLICIT NONE

        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        CLASS(ppm_t_field_),     POINTER       :: FieldIn
        !!! the source filed which the label image is supposed to be kept on that
        CLASS(ppm_t_equi_mesh_), POINTER       :: MeshIn
        !!! the FieldIn is discretized on this mesh

        CHARACTER(LEN=*),        INTENT(IN   ) :: InputFileName
        !!! Input file's name
        CHARACTER(LEN=*),        INTENT(IN   ) :: OutputFileName
        !!! Output file's name

        INTEGER,                 INTENT(  OUT) :: info
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
        CLASS(ppm_t_subpatch_), POINTER :: sbpitr

        REAL(ppm_kind_double) :: t0

#if    __DIME == __2D
        INTEGER, CONTIGUOUS, DIMENSION(:,:),   POINTER :: wpi
        INTEGER, CONTIGUOUS, DIMENSION(:),     POINTER :: wpip
#elif  __DIME == __3D
        INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER :: wpi
        INTEGER, CONTIGUOUS, DIMENSION(:,:),   POINTER :: wpip
        INTEGER                                        :: k,page
#endif

        INTEGER,             DIMENSION(:),     POINTER :: Nm
        INTEGER,             DIMENSION(__DIME)         :: Ngrid_
        INTEGER                                        :: i,j,z
        INTEGER                                        :: nlb
        INTEGER                                        :: nsize
        INTEGER                                        :: clabel
        INTEGER                                        :: vlabel
        INTEGER                                        :: istrip
        INTEGER                                        :: rownm

#if    __DIME == __3D
        LOGICAL :: IsBigTIFF
#endif

        CHARACTER(LEN=ppm_char) :: filename
        CHARACTER(LEN=ppm_char) :: caller='ppm_rc_read_write_unique_label'

        !-------------------------------------------------------------------------
        !  Initialize
        !-------------------------------------------------------------------------
        CALL substart(caller,t0,info)

        CALL check()

        NULLIFY(wpi)

        !---------------------------------------------------------------------
        ! First we find out which image slice to read
        ! We use the istart array for that
        !---------------------------------------------------------------------
        sbpitr => MeshIn%subpatch%begin()
        DO WHILE (ASSOCIATED(sbpitr))
           Nm => sbpitr%nnodes

           CALL sbpitr%get_field(FieldIn,wpi,info)
           or_fail("Failed to get field wp.")

#if    __DIME == __2D
           WRITE(filename,'(A,A4,A1)') TRIM(InputFileName),'.tif',CHAR(0)

           info=ppm_rc_read_tiff_info(filename,Ngrid_)
           or_fail('Error in reading a tiff header file!')

           info=ppm_rc_open_tiff(filename)
           or_fail('Error opening tiff file.')

           z=0
           !loop through slices
           DO istrip=sbpitr%istart(__DIME),sbpitr%iend(__DIME)
              z=z+1
              !convert the row number to C style
              rownm=istrip-1

              wpip => wpi(1:Nm(1),z)

              info=ppm_rc_read_tiff_scanline(wpip,rownm,bitsPerSample,0,Nm(1))
           ENDDO !istrip=sbpitr%istart(__DIME),sbpitr%iend(__DIME)

           info=ppm_rc_close_tiff()
           or_fail('Error closing tiff file.')
#elif  __DIME == __3D
           z=0

           SELECT CASE (inputformat)
           CASE (1)
              !loop through images
              DO i=1,ninputimage
                 SELECT CASE (ninputimage)
                 CASE (1)
                    WRITE(filename,'(A,A4,A1)') TRIM(InputFileName),'.tif',CHAR(0)
                 CASE (2:9)
                    WRITE(filename,'(A,A1,I1.1,A4,A1)') TRIM(InputFileName),'_',i,'.tif',CHAR(0)
                 CASE (10:99)
                    WRITE(filename,'(A,A1,I2.2,A4,A1)') TRIM(InputFileName),'_',i,'.tif',CHAR(0)
                 CASE (100:999)
                    WRITE(filename,'(A,A1,I3.3,A4,A1)') TRIM(InputFileName),'_',i,'.tif',CHAR(0)
                 CASE (1000:9999)
                    WRITE(filename,'(A,A1,I4.4,A4,A1)') TRIM(InputFileName),'_',i,'.tif',CHAR(0)
                 CASE (10000:99999)
                    WRITE(filename,'(A,A1,I5.5,A4,A1)') TRIM(InputFileName),'_',i,'.tif',CHAR(0)
                 END SELECT

                 info=ppm_rc_read_tiff_info(filename,Ngrid_)
                 or_fail('Error in reading a tiff header file!')

                 !The TIFF file format uses 32bit offsets and, as such, is limited to 4 gigabytes
                 IsBigTIFF=REAL(Ngrid_(1),MK)*REAL(Ngrid_(2),MK)*REAL(Ngrid_(3),MK)*REAL(bitsPerSample/8,MK)/1024._MK.GE.4194304._MK

                 SELECT CASE (IsBigTIFF)
                 CASE (.TRUE.)
                    info=ppm_rc_open_bigtiff(filename)
                 CASE DEFAULT
                    info=ppm_rc_open_tiff(filename)
                 END SELECT !(IsBigTIFF)
                 or_fail('Error opening tiff file.')

                 DO page=1,Ngrid_(3)
                    z=z+1

                    wpip => wpi(1:Nm(1),1:Nm(2),z)

                    info=ppm_rc_read_tiff_scanline(wpip,0,bitsPerSample,0,Nm(1),Nm(2),page-1,Ngrid_(3))
                    or_fail('Error in ppm_rc_read_tiff_scanline !!!')
                 ENDDO

                 info=ppm_rc_close_tiff()
                 or_fail('Error closing tiff file.')
              ENDDO !i=1,ninputimage
           CASE (2)
              !loop through images
              DO i=1,ninputimage
                 SELECT CASE (ninputimage)
                 CASE (2:9)
                    WRITE(filename,'(A,I1.1,A4,A1)') TRIM(InputFileName),i,'.tif',CHAR(0)
                 CASE (10:99)
                    WRITE(filename,'(A,I2.2,A4,A1)') TRIM(InputFileName),i,'.tif',CHAR(0)
                 CASE (100:999)
                    WRITE(filename,'(A,I3.3,A4,A1)') TRIM(InputFileName),i,'.tif',CHAR(0)
                 CASE (1000:9999)
                    WRITE(filename,'(A,I4.4,A4,A1)') TRIM(InputFileName),i,'.tif',CHAR(0)
                 CASE (10000:99999)
                    WRITE(filename,'(A,I5.5,A4,A1)') TRIM(InputFileName),i,'.tif',CHAR(0)
                 END SELECT

                 info=ppm_rc_read_tiff_info(filename,Ngrid_)
                 or_fail('Error in reading a tiff header file!')

                 info=ppm_rc_open_tiff(filename)
                 or_fail('Error opening tiff file.')

                 DO page=1,Ngrid_(3)
                    z=z+1

                    wpip => wpi(1:Nm(1),1:Nm(2),z)

                    info=ppm_rc_read_tiff_scanline(wpip,0,bitsPerSample,0,Nm(1),Nm(2),page-1,Ngrid_(3))
                    or_fail('Error in ppm_rc_read_tiff_scanline !!!')
                 ENDDO

                 info=ppm_rc_close_tiff()
                 or_fail('Error closing tiff file.')
              ENDDO !i=1,ninputimage
           CASE (3)
              !loop through images
              DO i=1,ninputimage
                 WRITE(filename,'(A,I0,A4,A1)') TRIM(InputFileName),i,'.tif',CHAR(0)

                 info=ppm_rc_read_tiff_info(filename,Ngrid_)
                 or_fail('Error in reading a tiff header file!')

                 info=ppm_rc_open_tiff(filename)
                 or_fail('Error opening tiff file.')

                 DO page=1,Ngrid_(3)
                    z=z+1

                    wpip => wpi(1:Nm(1),1:Nm(2),z)

                    info=ppm_rc_read_tiff_scanline(wpip,0,bitsPerSample,0,Nm(1),Nm(2),page-1,Ngrid_(3))
                    or_fail('Error in ppm_rc_read_tiff_scanline !!!')
                 ENDDO

                 info=ppm_rc_close_tiff()
                 or_fail('Error closing tiff file.')
              ENDDO !i=1,ninputimage
           END SELECT
#endif

           sbpitr => MeshIn%subpatch%next()
        ENDDO !(ASSOCIATED(sbpitr))
        NULLIFY(wpi)


        nsize=16
        ALLOCATE(nlabels(nsize),STAT=info)
        or_fail_alloc("nlabels")

        nlb=0
        vlabel=65536

        sbpitr => MeshIn%subpatch%begin()
        DO WHILE (ASSOCIATED(sbpitr))
           Nm => sbpitr%nnodes

           CALL sbpitr%get_field(FieldIn,wpi,info)
           or_fail("Failed to get field wp.")

#if    __DIME == __2D
           DO j=1,Nm(2)
              DO i=1,Nm(1)
                 ! Current label
                 clabel=ABS(wpi(i,j))

                 IF (clabel.NE.0) THEN
                    IF (.NOT.ppm_rc_label_exist(clabel,nlabels,nlb)) THEN
                       ! New label
                       vlabel=vlabel+1

                       IF (nlb+1.GT.nsize) THEN
                          ALLOCATE(bufi(nsize*2),STAT=info)
                          or_fail_alloc("Tmp array allocation failed!")

                          FORALL (z=1:nlb) bufi(z)=nlabels(z)
                          nsize=nsize*2

                          CALL MOVE_ALLOC(bufi,nlabels)
                       ENDIF

                       nlb=nlb+1
                       nlabels(nlb)=vlabel

                       CALL ppm_rc_floodFill(wpi,Nm,(/i,j/),clabel,vlabel,1,info)
                       or_fail("ppm_rc_floodFill")
                    ENDIF
                 ENDIF
              ENDDO !i=1,Nm(1)
           ENDDO !j=1,Nm(2)
#elif  __DIME == __3D
           DO k=1,Nm(3)
              DO j=1,Nm(2)
                 DO i=1,Nm(1)
                    ! Current label
                    clabel=ABS(wpi(i,j,k))

                    IF (clabel.NE.0) THEN
                       IF (.NOT.ppm_rc_label_exist(clabel,nlabels,nlb)) THEN
                          ! New label
                          vlabel=vlabel+1

                          IF (nlb+1.GT.nsize) THEN
                             ALLOCATE(bufi(nsize*2),STAT=info)
                             or_fail_alloc("Tmp array allocation failed!")

                             FORALL (z=1:nlb) bufi(z)=nlabels(z)
                             nsize=nsize*2

                             CALL MOVE_ALLOC(bufi,nlabels)
                          ENDIF

                          nlb=nlb+1
                          nlabels(nlb)=vlabel

                          CALL ppm_rc_floodFill(wpi,Nm,(/i,j,k/),clabel,vlabel,1,info)
                          or_fail("ppm_rc_floodFill")
                       ENDIF
                    ENDIF
                 ENDDO !i=1,Nm(1)
              ENDDO !j=1,Nm(2)
           ENDDO !k=1,Nm(3)
#endif

           DEALLOCATE(nlabels,STAT=info)
           or_fail_dealloc("nlabels")

           sbpitr => MeshIn%subpatch%next()
        ENDDO !(ASSOCIATED(sbpitr))
        NULLIFY(wpi)

        ! The first label starts from 2
        ! This is done to avoid extra communication at forestfire
        IF (nlb.LT.65535) THEN
           sbpitr => MeshIn%subpatch%begin()
           DO WHILE (ASSOCIATED(sbpitr))
              Nm => sbpitr%nnodes

              CALL sbpitr%get_field(FieldIn,wpi,info)
              or_fail("Failed to get field wp.")

#if    __DIME == __2D
              FORALL (i=1:Nm(1),j=1:Nm(2),wpi(i,j).GT.0)
                 wpi(i,j)=wpi(i,j)-65535
              END FORALL
#elif  __DIME == __3D
              FORALL (i=1:Nm(1),j=1:Nm(2),k=1:Nm(3),wpi(i,j,k).GT.0)
                 wpi(i,j,k)=wpi(i,j,k)-65535
              END FORALL
#endif

              sbpitr => MeshIn%subpatch%next()
           ENDDO !(ASSOCIATED(sbpitr))
        ENDIF
        NULLIFY(wpi)

        SELECT CASE (nlb)
        CASE (0:254)
           CALL DTYPE(ppm_rc_write_image_label)(FieldIn,MeshIn,OutputFileName,info, &
           &    bitsPerSampleW=8,lcast=.TRUE.)
        CASE (255:65534)
           CALL DTYPE(ppm_rc_write_image_label)(FieldIn,MeshIn,OutputFileName,info, &
           &    bitsPerSampleW=16,lcast=.TRUE.)
        CASE DEFAULT
           CALL DTYPE(ppm_rc_write_image_label)(FieldIn,MeshIn,OutputFileName,info, &
           &    bitsPerSampleW=32,lcast=.TRUE.)
        END SELECT
        or_fail("ppm_rc_write_image failed!")

        stdout("")
        CALL ppm_log(caller,cbuf,info)
        stdout("The algorithm found : ",nlb," unique labels.")
        CALL ppm_log(caller,cbuf,info)
        stdout("")

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
          & "Input Mesh does not exist",exit_point=8888)

          check_true(<#FieldIn%is_discretized_on(MeshIn)#>, &
          & "This Field is not discretized on Input Mesh",exit_point=8888)

          check_true(<#FieldIn%data_type.EQ.ppm_type_int#>, &
          & "This label field data type is not correct",exit_point=8888)

          check_true(<#ppm_nproc.EQ.1#>, &
          & "Creating unique label is only working on one processor",exit_point=8888)

          check_true(<#MeshIn%subpatch%nb.EQ.1#>, &
          & "This mesh should only have one subpatch",exit_point=8888)
        8888 CONTINUE
          RETURN
        END SUBROUTINE
      END SUBROUTINE DTYPE(ppm_rc_read_write_unique_label)
