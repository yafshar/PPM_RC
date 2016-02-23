      !-------------------------------------------------------------------------
      !  Subroutine   :                    ppm_rc_read_image_info
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Read a strip of a TIF image to get the size of the
      !                 'image' mesh
      !
      !  Input        : root name of the files to be read.
      !
      !  Output       : info          (I) return status. 0 on success.
      !
      !  Routines     : EXTERNAL C function call
      !
      !  Remarks      : There is no assumtion about each .TIF file
      !                 So, basically it could be a stack of 2D strips or
      !                 several .TIF files.
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  MOSAIC Group
      !  Max Planck Institute of Molecular Cell Biology and Genetics
      !  Pfotenhauerstr. 108, 01307 Dresden, Germany
      !
      !  Author           - y.afshar           June   2014
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_rc_read_image_info(InputFileName,NmInputFile,info)

        !-------------------------------------------------------------------------
        !  Modules
        !-------------------------------------------------------------------------
        USE ppm_rc_module_tiff, ONLY : read_tiff_info
        IMPLICIT NONE

        !-------------------------------------------------------------------------
        !  Includes
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        CHARACTER(LEN=*), INTENT(IN   ) :: InputFileName

        INTEGER,          INTENT(IN   ) :: NmInputFile
        INTEGER,          INTENT(  OUT) :: info
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
        REAL(ppm_kind_double) :: t0

        INTEGER               :: i,z
        INTEGER, DIMENSION(3) :: Ngrid_

        CHARACTER(LEN=ppm_char) :: filename
        CHARACTER(LEN=ppm_char) :: caller='ppm_rc_read_image_info'

        !-------------------------------------------------------------------------
        !  Initialize
        !-------------------------------------------------------------------------
        CALL substart(caller,t0,info)

        z=0
        DO i=1,NmInputFile
           SELECT CASE (NmInputFile)
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

           info=read_tiff_info(filename,Ngrid_)
           or_fail('Error in reading a tiff header file!')

           z=z+MERGE(Ngrid_(3),1,ppm_rc_dim.EQ.3)
        ENDDO !i=1,NmInputFile

        IF (z.LT.NmInputFile) THEN
           stdout_f('(A,I7,2X,A,I7)',"The number of image z slices",'z', &
           & "is less than the number of input images",NmInputFile)
           fail("Failed",ppm_error=ppm_error_fatal)
        ENDIF
        IF (ppm_rc_dim.EQ.3.AND.z.EQ.1) THEN
           stdout_f('(A,I7,2X,A,I7)',"The number of image z slices=",'z', &
           & "does not match with the selected DIMENSIONALITY of ",ppm_rc_dim)
           fail("Failed",ppm_error=ppm_error_fatal)
        ENDIF

        IF (ALL(Ngrid(1:ppm_rc_dim).EQ.0)) THEN
           Ngrid(1:2)=Ngrid_(1:2)
           Ngrid(3)=z
        ELSE
           check_true(<#(Ngrid_(1).EQ.Ngrid(1).AND.Ngrid_(2).EQ.Ngrid(2).AND.z.EQ.Ngrid(3))#>, &
           & "The input and init image sizes do not match, you can not use this init image for Initialization!!!")
        ENDIF

        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
      9999 CONTINUE

        CALL substop(caller,t0,info)
        RETURN
      END SUBROUTINE ppm_rc_read_image_info
