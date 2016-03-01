      !-------------------------------------------------------------------------
      !  Subroutine   :                    ppm_rc_read_image
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
      !  Subroutine   :                    ppm_rc_read_image
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Read a bunch of TIFF slabs using Libtiff library.
      !                 Fill up the 3D mesh with intensity values obtained from
      !                 the images.
      !
      !  Input        : root name of the files to be read.
      !
      !  Output       : info          (I) return status. 0 on success.
      !
      !  Routines     :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      SUBROUTINE DTYPE(ppm_rc_read_image)(FieldIn,MeshIn,InputFileName,info)

        !-------------------------------------------------------------------------
        !  Modules
        !-------------------------------------------------------------------------
        USE ppm_module_field_typedef, ONLY : ppm_t_field_
        USE ppm_module_mesh_typedef, ONLY : ppm_t_subpatch_

        USE ppm_rc_module_tiff, ONLY : ppm_rc_open_tiff,ppm_rc_open_bigtiff, &
        &   ppm_rc_read_tiff_scanline,ppm_rc_close_tiff,bitsPerSample,       &
        &   samplesPerPixel
        IMPLICIT NONE

        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        CLASS(ppm_t_field_),     POINTER       :: FieldIn
        !!! the source filed which the image is supposed to be kept on that

        CLASS(ppm_t_equi_mesh_), POINTER       :: MeshIn
        !!! the FieldIn is discretized on this mesh

        CHARACTER(LEN=*),        INTENT(IN   ) :: InputFileName
        ! Input file's name

        INTEGER,                 INTENT(  OUT) :: info
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
        CLASS(ppm_t_subpatch_), POINTER :: sbpitr

#if   __DIME == __2D
        REAL(ppm_kind_single), CONTIGUOUS, DIMENSION(:,:),   POINTER :: DTYPE(wp_rs)
        REAL(MK),              CONTIGUOUS, DIMENSION(:,:),   POINTER :: DTYPE(wp_r)
#elif __DIME == __3D
        REAL(ppm_kind_single), CONTIGUOUS, DIMENSION(:,:,:), POINTER :: DTYPE(wp_rs)
        REAL(MK),              CONTIGUOUS, DIMENSION(:,:,:), POINTER :: DTYPE(wp_r)
#endif
        REAL(ppm_kind_double)                                        :: t0

#if   __DIME == __2D
        INTEGER, CONTIGUOUS, DIMENSION(:,:),   POINTER :: DTYPE(wp_i)
#elif __DIME == __3D
        INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER :: DTYPE(wp_i)
#endif
        INTEGER,             DIMENSION(:),     POINTER :: Nm
#if   __DIME == __2D
        INTEGER                                        :: rownm
#elif __DIME == __3D
        INTEGER                                        :: page,npages
#endif
        INTEGER                                        :: z,istrip

        CHARACTER(LEN=ppm_char) :: filename
        CHARACTER(LEN=ppm_char) :: caller='ppm_rc_read_image'

        LOGICAL :: IsBigTIFF
        !-------------------------------------------------------------------------
        !  Initialize
        !-------------------------------------------------------------------------
        CALL substart(caller,t0,info)

        CALL check()

        !-------------------------------------------------------------------------
        !  Create format string
        !-------------------------------------------------------------------------
        SELECT CASE (ninputimage)
        CASE (1)
           WRITE(filename,'(A,A4,A1)') TRIM(InputFileName),'.tif',CHAR(0)

           !The TIFF file format uses 32bit offsets and, as such, is limited to 4 gigabytes
#if   __DIME == __2D
           IsBigTIFF=REAL(Ngrid(1),MK)*REAL(Ngrid(2),MK)*REAL(bitsPerSample*samplesPerPixel/8,MK)/1024._MK.GE.4194304._MK
#elif __DIME == __3D
           IsBigTIFF=REAL(Ngrid(1),MK)*REAL(Ngrid(2),MK)*REAL(Ngrid(3),MK)*REAL(bitsPerSample*samplesPerPixel/8,MK)/1024._MK.GE.4194304._MK
#endif

           SELECT CASE (IsBigTIFF)
           CASE (.TRUE.)
              info=ppm_rc_open_bigtiff(filename)
           CASE DEFAULT
              info=ppm_rc_open_tiff(filename)
           END SELECT !(IsBigTIFF)
           or_fail('Error opening tiff file.')

        END SELECT !(ninputimage)

        !---------------------------------------------------------------------
        ! First we find out which image slice to read
        ! We use the istart array for that
        !---------------------------------------------------------------------
        NULLIFY(DTYPE(wp_i),DTYPE(wp_r),DTYPE(wp_rs))
        sbpitr => MeshIn%subpatch%begin()
        DO WHILE (ASSOCIATED(sbpitr))
           Nm => sbpitr%nnodes

           SELECT CASE (FieldIn%data_type)
           CASE (ppm_type_real_single)
              CALL sbpitr%get_field(FieldIn,DTYPE(wp_rs),info)
              or_fail("Failed to get field wp_rs.")

           CASE (ppm_type_real)
              CALL sbpitr%get_field(FieldIn,DTYPE(wp_r),info)
              or_fail("Failed to get field wp_r.")

           CASE (ppm_type_int)
              CALL sbpitr%get_field(FieldIn,DTYPE(wp_i),info)
              or_fail("Failed to get field wp_i.")

           END SELECT !(FieldIn%data_type)

#if   __DIME == __3D
           SELECT CASE (ninputimage)
           CASE (1)
              npages=Ngrid(__DIME)
              !sbpitr%iend(__DIME)-sbpitr%istart(__DIME)+1

           CASE DEFAULT
              page=0
              npages=1

           END SELECT
#endif

           z=0

           SELECT CASE (samplesPerPixel)
           CASE (1)
              !loop through slices
              DO istrip=sbpitr%istart(__DIME),sbpitr%iend(__DIME)
                 z=z+1

#if   __DIME == __2D
                 !convert the row number to C style
                 rownm=istrip-1

                 SELECT CASE (FieldIn%data_type)
                 CASE (ppm_type_real_single)
                    info=ppm_rc_read_tiff_scanline(DTYPE(wp_rs)(1:Nm(1),z),rownm,bitsPerSample,0,Nm(1))

                 CASE (ppm_type_real)
                    info=ppm_rc_read_tiff_scanline(DTYPE(wp_r)(1:Nm(1),z),rownm,bitsPerSample,0,Nm(1))

                 CASE (ppm_type_int)
                    info=ppm_rc_read_tiff_scanline(DTYPE(wp_i)(1:Nm(1),z),rownm,bitsPerSample,0,Nm(1))

                 END SELECT !(FieldIn%data_type)
                 or_fail('Error in ppm_rc_read_tiff_scanline !!!')
#elif __DIME == __3D
                 SELECT CASE (ninputimage)
                 CASE (1)
                    page=istrip-1

                 CASE (2:9)
                    WRITE(filename,'(A,A1,I1.1,A4,A1)') TRIM(InputFileName),'_',istrip,'.tif',CHAR(0)
                    info=ppm_rc_open_tiff(filename)
                    or_fail('Error opening tiff file.')

                 CASE (10:99)
                    WRITE(filename,'(A,A1,I2.2,A4,A1)') TRIM(InputFileName),'_',istrip,'.tif',CHAR(0)
                    info=ppm_rc_open_tiff(filename)
                    or_fail('Error opening tiff file.')

                 CASE (100:999)
                    WRITE(filename,'(A,A1,I3.3,A4,A1)') TRIM(InputFileName),'_',istrip,'.tif',CHAR(0)
                    info=ppm_rc_open_tiff(filename)
                    or_fail('Error opening tiff file.')

                 CASE (1000:9999)
                    WRITE(filename,'(A,A1,I4.4,A4,A1)') TRIM(InputFileName),'_',istrip,'.tif',CHAR(0)
                    info=ppm_rc_open_tiff(filename)
                    or_fail('Error opening tiff file.')

                 CASE (10000:99999)
                    WRITE(filename,'(A,A1,I5.5,A4,A1)') TRIM(InputFileName),'_',istrip,'.tif',CHAR(0)
                    info=ppm_rc_open_tiff(filename)
                    or_fail('Error opening tiff file.')

                 END SELECT

                 SELECT CASE (FieldIn%data_type)
                 CASE (ppm_type_real_single)
                    info=ppm_rc_read_tiff_scanline(DTYPE(wp_rs)(1:Nm(1),1:Nm(2),z),0,bitsPerSample,0,Nm(1),Nm(2),page,npages)

                 CASE (ppm_type_real)
                    info=ppm_rc_read_tiff_scanline(DTYPE(wp_r)(1:Nm(1),1:Nm(2),z),0,bitsPerSample,0,Nm(1),Nm(2),page,npages)

                 CASE (ppm_type_int)
                    info=ppm_rc_read_tiff_scanline(DTYPE(wp_i)(1:Nm(1),1:Nm(2),z),0,bitsPerSample,0,Nm(1),Nm(2),page,npages)

                 END SELECT !(FieldIn%data_type)
                 or_fail('Error in ppm_rc_read_tiff_scanline !!!')

                 SELECT CASE (ninputimage)
                 CASE (1)
                 CASE DEFAULT
                    info=ppm_rc_close_tiff()
                    or_fail('Error closing tiff file.')

                 END SELECT !(ninputimage)
#endif
              ENDDO !istrip=sbpitr%istart(__DIME),sbpitr%iend(__DIME)

           END SELECT !(samplesPerPixel)

           sbpitr => MeshIn%subpatch%next()
        ENDDO !(ASSOCIATED(sbpitr))

        SELECT CASE (ninputimage)
        CASE (1)
           info=ppm_rc_close_tiff()
           or_fail('Error closing tiff file.')

        END SELECT !(ninputimage)
        NULLIFY(DTYPE(wp_i),DTYPE(wp_r),DTYPE(wp_rs))

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

        8888 CONTINUE
          RETURN
        END SUBROUTINE
      END SUBROUTINE DTYPE(ppm_rc_read_image)
