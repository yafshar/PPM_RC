      !-------------------------------------------------------------------------
      !  Subroutine   :                    ppm_rc_read_write_image
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
      !  Subroutine   :                    ppm_rc_read_write_image
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Read a bunch of TIFF slabs using Libtiff library from several images.
      !                 Fill up the 3D mesh with intensity values obtained from
      !                 the images and write the resulted read into OutputFileName
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
      SUBROUTINE DTYPE(ppm_rc_read_write_image)(FieldIn,MeshIn,InputFileName,OutputFileName,info)

        !-------------------------------------------------------------------------
        !  Modules
        !-------------------------------------------------------------------------
        USE ppm_module_field_typedef, ONLY : ppm_t_field_
        USE ppm_module_mesh_typedef, ONLY : ppm_t_subpatch_

        USE ppm_rc_module_tiff, ONLY : open_tiff,read_tiff_scanline, &
        &   close_tiff,read_tiff_info,bitsPerSample
        IMPLICIT NONE

        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        CLASS(ppm_t_field_),     POINTER       :: FieldIn
        !!! the source filed which the image is supposed to be kept on that
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
#if  __DIME == __3D
        CLASS(ppm_t_subpatch_), POINTER :: sbpitr

        REAL(ppm_kind_single), CONTIGUOUS, DIMENSION(:,:,:), POINTER :: DTYPE(wp_rs) => NULL()
        REAL(MK),              CONTIGUOUS, DIMENSION(:,:,:), POINTER :: DTYPE(wp_r) => NULL()
#endif
        REAL(ppm_kind_double)                                        :: t0

#if   __DIME == __3D
        INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER :: DTYPE(wp_i) => NULL()
        INTEGER,             DIMENSION(:),     POINTER :: Nm
        INTEGER,             DIMENSION(3)              :: Ngrid_
        INTEGER                                        :: i,z,page

        CHARACTER(LEN=ppm_char) :: filename
#endif
        CHARACTER(LEN=ppm_char) :: caller='ppm_rc_read_write_image'

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
           GOTO 9999
        END SELECT

        !This part is only working on one processor
        SELECT CASE (ppm_nproc)
        CASE (1)
        CASE DEFAULT
           fail("Combining tiff files into one is only working on one processor", &
           & ppm_error=ppm_error_fatal)
        END SELECT

        SELECT CASE (MeshIn%subpatch%nb)
        CASE (1)
        CASE DEFAULT
           fail("This mesh should only have one subpatch", &
           & ppm_error=ppm_error_fatal)
        END SELECT

#if   __DIME == __3D
        !---------------------------------------------------------------------
        ! First we find out which image slice to read
        ! We use the istart array for that
        !---------------------------------------------------------------------
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

           END SELECT

           z=0

           !loop through images
           DO i=1,ninputimage
              SELECT CASE (ninputimage)
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

              info=read_tiff_info(filename,Ngrid_)
              or_fail('Error in reading a tiff header file!')

              info=open_tiff(filename)
              or_fail('Error opening tiff file.')

              DO page=1,Ngrid_(3)
                 z=z+1
                 SELECT CASE (FieldIn%data_type)
                 CASE (ppm_type_real_single)
                    info=read_tiff_scanline(DTYPE(wp_rs)(1:Nm(1),1:Nm(2),z), &
                    &    0,bitsPerSample,0,Nm(1),Nm(2),page-1,Ngrid_(3))

                 CASE (ppm_type_real)
                    info=read_tiff_scanline(DTYPE(wp_r)(1:Nm(1),1:Nm(2),z), &
                    &    0,bitsPerSample,0,Nm(1),Nm(2),page-1,Ngrid_(3))

                 CASE (ppm_type_int)
                    info=read_tiff_scanline(DTYPE(wp_i)(1:Nm(1),1:Nm(2),z), &
                    &    0,bitsPerSample,0,Nm(1),Nm(2),page-1,Ngrid_(3))

                 END SELECT
                 or_fail('Error in read_tiff_scanline !!!')
              ENDDO

              info=close_tiff()
              or_fail('Error closing tiff file.')
           ENDDO !i=1,ninputimage

           sbpitr => MeshIn%subpatch%next()
        ENDDO !(ASSOCIATED(sbpitr))

        CALL DTYPE(ppm_rc_write_image)(FieldIn,MeshIn,OutputFileName,info)
        or_fail("ppm_rc_write_image failed!")

        NULLIFY(DTYPE(wp_i),DTYPE(wp_rs),DTYPE(wp_r))
#endif

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
        8888 CONTINUE
          RETURN
        END SUBROUTINE
      END SUBROUTINE DTYPE(ppm_rc_read_write_image)
