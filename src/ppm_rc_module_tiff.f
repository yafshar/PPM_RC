      !-------------------------------------------------------------------------
      !  Module       :                    ppm_rc_module_tiff
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
      !  Author           - y.afshar           March   2015
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Module       :                    ppm_rc_module_tiff
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Reading and Writing TIFF files
      !
      !  Remarks      : Read_tiff_strip is only working for 3D data
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      MODULE ppm_rc_module_tiff
#ifdef __F2003
#ifdef __TIFF
        USE ISO_C_BINDING

        USE ppm_module_data
        IMPLICIT NONE

        PRIVATE

        INTEGER(C_INT) :: bitsPerSample=8
        INTEGER(C_INT) :: bitsPerSampleW=0
        !!! Bits Per Sample in TIFF file
        INTEGER(C_INT) :: samplesPerPixel=1
        !!! Number of samples in the TIFF file

        !----------------------------------------------------------------------
        !  Some work memory on the heap
        !----------------------------------------------------------------------
        !----------------------------------------------------------------------
        !  Define module interfaces
        !----------------------------------------------------------------------
        INTERFACE
          FUNCTION ppm_rc_open_tiff(filename) bind(C,NAME='ppm_rc_open_tiff')
            IMPORT            :: C_INT,C_CHAR
            CHARACTER(C_CHAR) :: filename(*)
            !!! Name of the TIFF file to open
            INTEGER(C_INT)    :: ppm_rc_open_tiff
          END FUNCTION
        END INTERFACE

        INTERFACE
          FUNCTION ppm_rc_open_bigtiff(filename) bind(C,NAME='ppm_rc_open_bigtiff')
            IMPORT            :: C_INT,C_CHAR
            CHARACTER(C_CHAR) :: filename(*)
            !!! Name of the TIFF file to open
            INTEGER(C_INT)    :: ppm_rc_open_bigtiff
          END FUNCTION
        END INTERFACE

        INTERFACE
          FUNCTION ppm_rc_open_write_tiff(filename) bind(C,NAME='ppm_rc_open_write_tiff')
            IMPORT            :: C_INT,C_CHAR
            CHARACTER(C_CHAR) :: filename(*)
            INTEGER(C_INT)    :: ppm_rc_open_write_tiff
          END FUNCTION
        END INTERFACE

        INTERFACE
          FUNCTION ppm_rc_open_write_bigtiff(filename) bind(C,NAME='ppm_rc_open_write_bigtiff')
            IMPORT            :: C_INT,C_CHAR
            CHARACTER(C_CHAR) :: filename(*)
            INTEGER(C_INT)    :: ppm_rc_open_write_bigtiff
          END FUNCTION
        END INTERFACE

        INTERFACE
          FUNCTION ppm_rc_close_tiff() bind(C,NAME='ppm_rc_close_tiff')
            IMPORT         :: C_INT
            INTEGER(C_INT) :: ppm_rc_close_tiff
          END FUNCTION
        END INTERFACE

        INTERFACE
          FUNCTION read_tiff_infoC(filename,ngrid,bitsPerSample_,samplesPerPixel_,rank,comm) bind(C,NAME='read_tiff_infoC')
            IMPORT            :: C_INT,C_CHAR
            CHARACTER(C_CHAR) :: filename(*)
            INTEGER(C_INT)    :: ngrid(3)
            INTEGER(C_INT)    :: bitsPerSample_
            INTEGER(C_INT)    :: samplesPerPixel_
            INTEGER(C_INT)    :: rank
            INTEGER(C_INT)    :: comm
            INTEGER(C_INT)    :: read_tiff_infoC
          END FUNCTION
        END INTERFACE
        INTERFACE ppm_rc_write_tiff_header
          FUNCTION write_tiff_header1(bitsPerSample_,samplesPerPixel_,width,length) bind(C,NAME='write_tiff_header1')
            IMPORT         :: C_INT
            INTEGER(C_INT) :: bitsPerSample_
            INTEGER(C_INT) :: samplesPerPixel_
            INTEGER(C_INT) :: width
            INTEGER(C_INT) :: length
            INTEGER(C_INT) :: write_tiff_header1
          END FUNCTION
          FUNCTION write_tiff_header2(page,npages) bind(C,NAME='write_tiff_header2')
            IMPORT         :: C_INT
            INTEGER(C_INT) :: page
            INTEGER(C_INT) :: npages
            INTEGER(C_INT) :: write_tiff_header2
          END FUNCTION
        END INTERFACE

        INTERFACE ppm_rc_read_tiff_scanline
          FUNCTION read_tiff_scanline_int1(scanline,rownm,bitsPerSample_,Sample,width) bind(C,NAME='read_tiff_scanline_int1')
            IMPORT         :: C_INT
            INTEGER(C_INT) :: width
            INTEGER(C_INT) :: scanline(width)
            INTEGER(C_INT) :: rownm
            INTEGER(C_INT) :: bitsPerSample_
            INTEGER(C_INT) :: Sample
            INTEGER(C_INT) :: read_tiff_scanline_int1
          END FUNCTION
          FUNCTION read_tiff_scanline_long1(scanline,rownm,bitsPerSample_,Sample,width) bind(C,NAME='read_tiff_scanline_long1')
            IMPORT          :: C_INT,C_LONG
            INTEGER(C_INT)  :: width
            INTEGER(C_LONG) :: scanline(width)
            INTEGER(C_INT)  :: rownm
            INTEGER(C_INT)  :: bitsPerSample_
            INTEGER(C_INT)  :: Sample
            INTEGER(C_INT)  :: read_tiff_scanline_long1
          END FUNCTION
          FUNCTION read_tiff_scanline_float1(scanline,rownm,bitsPerSample_,Sample,width) bind(C,NAME='read_tiff_scanline_float1')
            IMPORT         :: C_INT,C_FLOAT
            INTEGER(C_INT) :: width
            REAL(C_FLOAT)  :: scanline(width)
            INTEGER(C_INT) :: rownm
            INTEGER(C_INT) :: bitsPerSample_
            INTEGER(C_INT) :: Sample
            INTEGER(C_INT) :: read_tiff_scanline_float1
          END FUNCTION
          FUNCTION read_tiff_scanline_double1(scanline,rownm,bitsPerSample_,Sample,width) bind(C,NAME='read_tiff_scanline_double1')
            IMPORT         :: C_INT,C_DOUBLE
            INTEGER(C_INT) :: width
            REAL(C_DOUBLE) :: scanline(width)
            INTEGER(C_INT) :: rownm
            INTEGER(C_INT) :: bitsPerSample_
            INTEGER(C_INT) :: Sample
            INTEGER(C_INT) :: read_tiff_scanline_double1
          END FUNCTION

          FUNCTION read_tiff_scanline_int2(scanline,rownm,bitsPerSample_,Sample,width,length,page,npages) bind(C,NAME='read_tiff_scanline_int2')
            IMPORT         :: C_INT
            INTEGER(C_INT) :: width
            INTEGER(C_INT) :: length
            INTEGER(C_INT) :: scanline(width,length)
            INTEGER(C_INT) :: rownm
            INTEGER(C_INT) :: bitsPerSample_
            INTEGER(C_INT) :: Sample
            INTEGER(C_INT) :: page
            INTEGER(C_INT) :: npages
            INTEGER(C_INT) :: read_tiff_scanline_int2
          END FUNCTION
          FUNCTION read_tiff_scanline_long2(scanline,rownm,bitsPerSample_,Sample,width,length,page,npages) bind(C,NAME='read_tiff_scanline_long2')
            IMPORT          :: C_INT,C_LONG
            INTEGER(C_INT)  :: width
            INTEGER(C_INT)  :: length
            INTEGER(C_LONG) :: scanline(width,length)
            INTEGER(C_INT)  :: rownm
            INTEGER(C_INT)  :: bitsPerSample_
            INTEGER(C_INT)  :: Sample
            INTEGER(C_INT)  :: page
            INTEGER(C_INT)  :: npages
            INTEGER(C_INT)  :: read_tiff_scanline_long2
          END FUNCTION
          FUNCTION read_tiff_scanline_float2(scanline,rownm,bitsPerSample_,Sample,width,length,page,npages) bind(C,NAME='read_tiff_scanline_float2')
            IMPORT         :: C_INT,C_FLOAT
            INTEGER(C_INT) :: width
            INTEGER(C_INT) :: length
            REAL(C_FLOAT)  :: scanline(width,length)
            INTEGER(C_INT) :: rownm
            INTEGER(C_INT) :: bitsPerSample_
            INTEGER(C_INT) :: Sample
            INTEGER(C_INT) :: page
            INTEGER(C_INT) :: npages
            INTEGER(C_INT) :: read_tiff_scanline_float2
          END FUNCTION
          FUNCTION read_tiff_scanline_double2(scanline,rownm,bitsPerSample_,Sample,width,length,page,npages) bind(C,NAME='read_tiff_scanline_double2')
            IMPORT         :: C_INT,C_DOUBLE
            INTEGER(C_INT) :: width
            INTEGER(C_INT) :: length
            REAL(C_DOUBLE) :: scanline(width,length)
            INTEGER(C_INT) :: rownm
            INTEGER(C_INT) :: bitsPerSample_
            INTEGER(C_INT) :: Sample
            INTEGER(C_INT) :: page
            INTEGER(C_INT) :: npages
            INTEGER(C_INT) :: read_tiff_scanline_double2
          END FUNCTION
       END INTERFACE

       INTERFACE ppm_rc_read_tiff_strip
          !The strip function is only working for 3D data
          FUNCTION read_tiff_strip_int1(strip,bitsPerSample_,Sample,width,length,page,npages) bind(C,NAME='read_tiff_strip_int1')
            IMPORT         :: C_INT
            INTEGER(C_INT) :: width
            INTEGER(C_INT) :: length
            INTEGER(C_INT) :: strip(width,length)
            INTEGER(C_INT) :: bitsPerSample_
            INTEGER(C_INT) :: Sample
            INTEGER(C_INT) :: page
            INTEGER(C_INT) :: npages
            INTEGER(C_INT) :: read_tiff_strip_int1
          END FUNCTION
          FUNCTION read_tiff_strip_long1(strip,bitsPerSample_,Sample,width,length,page,npages) bind(C,NAME='read_tiff_strip_long1')
            IMPORT          :: C_INT,C_LONG
            INTEGER(C_INT)  :: width
            INTEGER(C_INT)  :: length
            INTEGER(C_LONG) :: strip(width,length)
            INTEGER(C_INT)  :: bitsPerSample_
            INTEGER(C_INT)  :: Sample
            INTEGER(C_INT)  :: page
            INTEGER(C_INT)  :: npages
            INTEGER(C_INT)  :: read_tiff_strip_long1
          END FUNCTION
          FUNCTION read_tiff_strip_float1(strip,bitsPerSample_,Sample,width,length,page,npages) bind(C,NAME='read_tiff_strip_float1')
            IMPORT         :: C_INT,C_FLOAT
            INTEGER(C_INT) :: width
            INTEGER(C_INT) :: length
            REAL(C_FLOAT)  :: strip(width,length)
            INTEGER(C_INT) :: bitsPerSample_
            INTEGER(C_INT) :: Sample
            INTEGER(C_INT) :: page
            INTEGER(C_INT) :: npages
            INTEGER(C_INT) :: read_tiff_strip_float1
          END FUNCTION
          FUNCTION read_tiff_strip_double1(strip,bitsPerSample_,Sample,width,length,page,npages) bind(C,NAME='read_tiff_strip_double1')
            IMPORT         :: C_INT,C_DOUBLE
            INTEGER(C_INT) :: width
            INTEGER(C_INT) :: length
            REAL(C_DOUBLE) :: strip(width,length)
            INTEGER(C_INT) :: bitsPerSample_
            INTEGER(C_INT) :: Sample
            INTEGER(C_INT) :: page
            INTEGER(C_INT) :: npages
            INTEGER(C_INT) :: read_tiff_strip_double1
          END FUNCTION
        END INTERFACE

        INTERFACE ppm_rc_write_tiff_scanline
          FUNCTION write_tiff_scanline_int1(scanline,rownm,bitsPerSample_,Sample,width) bind(C,NAME='write_tiff_scanline_int1')
            IMPORT         :: C_INT
            INTEGER(C_INT) :: width
            INTEGER(C_INT) :: scanline(width)
            INTEGER(C_INT) :: rownm
            INTEGER(C_INT) :: bitsPerSample_
            INTEGER(C_INT) :: Sample
            INTEGER(C_INT) :: write_tiff_scanline_int1
          END FUNCTION
          FUNCTION write_tiff_scanline_long1(scanline,rownm,bitsPerSample_,Sample,width) bind(C,NAME='write_tiff_scanline_long1')
            IMPORT          :: C_INT,C_LONG
            INTEGER(C_INT)  :: width
            INTEGER(C_LONG) :: scanline(width)
            INTEGER(C_INT)  :: rownm
            INTEGER(C_INT)  :: bitsPerSample_
            INTEGER(C_INT)  :: Sample
            INTEGER(C_INT)  :: write_tiff_scanline_long1
          END FUNCTION
          FUNCTION write_tiff_scanline_float1(scanline,rownm,bitsPerSample_,Sample,width) bind(C,NAME='write_tiff_scanline_float1')
            IMPORT         :: C_INT,C_FLOAT
            INTEGER(C_INT) :: width
            REAL(C_FLOAT)  :: scanline(width)
            INTEGER(C_INT) :: rownm
            INTEGER(C_INT) :: bitsPerSample_
            INTEGER(C_INT) :: Sample
            INTEGER(C_INT) :: write_tiff_scanline_float1
          END FUNCTION
          FUNCTION write_tiff_scanline_double1(scanline,rownm,bitsPerSample_,Sample,width) bind(C,NAME='write_tiff_scanline_double1')
            IMPORT         :: C_INT,C_DOUBLE
            INTEGER(C_INT) :: width
            REAL(C_DOUBLE) :: scanline(width)
            INTEGER(C_INT) :: rownm
            INTEGER(C_INT) :: bitsPerSample_
            INTEGER(C_INT) :: Sample
            INTEGER(C_INT) :: write_tiff_scanline_double1
          END FUNCTION
          FUNCTION write_tiff_scanline_int2(scanline,rownm,bitsPerSample_,Sample,width,length,page,npages) bind(C,NAME='write_tiff_scanline_int2')
            IMPORT         :: C_INT
            INTEGER(C_INT) :: width
            INTEGER(C_INT) :: length
            INTEGER(C_INT) :: scanline(width,length)
            INTEGER(C_INT) :: rownm
            INTEGER(C_INT) :: bitsPerSample_
            INTEGER(C_INT) :: Sample
            INTEGER(C_INT) :: page
            INTEGER(C_INT) :: npages
            INTEGER(C_INT) :: write_tiff_scanline_int2
          END FUNCTION
          FUNCTION write_tiff_scanline_long2(scanline,rownm,bitsPerSample_,Sample,width,length,page,npages) bind(C,NAME='write_tiff_scanline_long2')
            IMPORT          :: C_INT,C_LONG
            INTEGER(C_INT)  :: width
            INTEGER(C_INT)  :: length
            INTEGER(C_LONG) :: scanline(width,length)
            INTEGER(C_INT)  :: rownm
            INTEGER(C_INT)  :: bitsPerSample_
            INTEGER(C_INT)  :: Sample
            INTEGER(C_INT)  :: page
            INTEGER(C_INT)  :: npages
            INTEGER(C_INT)  :: write_tiff_scanline_long2
          END FUNCTION
          FUNCTION write_tiff_scanline_float2(scanline,rownm,bitsPerSample_,Sample,width,length,page,npages) bind(C,NAME='write_tiff_scanline_float2')
            IMPORT         :: C_INT,C_FLOAT
            INTEGER(C_INT) :: width
            INTEGER(C_INT) :: length
            REAL(C_FLOAT)  :: scanline(width,length)
            INTEGER(C_INT) :: rownm
            INTEGER(C_INT) :: bitsPerSample_
            INTEGER(C_INT) :: Sample
            INTEGER(C_INT) :: page
            INTEGER(C_INT) :: npages
            INTEGER(C_INT) :: write_tiff_scanline_float2
          END FUNCTION
          FUNCTION write_tiff_scanline_double2(scanline,rownm,bitsPerSample_,Sample,width,length,page,npages) bind(C,NAME='write_tiff_scanline_double2')
            IMPORT         :: C_INT,C_DOUBLE
            INTEGER(C_INT) :: width
            INTEGER(C_INT) :: length
            REAL(C_DOUBLE) :: scanline(width,length)
            INTEGER(C_INT) :: rownm
            INTEGER(C_INT) :: bitsPerSample_
            INTEGER(C_INT) :: Sample
            INTEGER(C_INT) :: page
            INTEGER(C_INT) :: npages
            INTEGER(C_INT) :: write_tiff_scanline_double2
          END FUNCTION
       END INTERFACE

       INTERFACE ppm_rc_write_tiff_strip
          FUNCTION write_tiff_strip_int1(strip,bitsPerSample_,width,length,page,npages) bind(C,NAME='write_tiff_strip_int1')
            IMPORT         :: C_INT
            INTEGER(C_INT) :: width
            INTEGER(C_INT) :: length
            INTEGER(C_INT) :: strip(width,length)
            INTEGER(C_INT) :: bitsPerSample_
            INTEGER(C_INT) :: page
            INTEGER(C_INT) :: npages
            INTEGER(C_INT) :: write_tiff_strip_int1
          END FUNCTION
          FUNCTION write_tiff_strip_long1(strip,bitsPerSample_,width,length,page,npages) bind(C,NAME='write_tiff_strip_long1')
            IMPORT          :: C_INT,C_LONG
            INTEGER(C_INT)  :: width
            INTEGER(C_INT)  :: length
            INTEGER(C_LONG) :: strip(width,length)
            INTEGER(C_INT)  :: bitsPerSample_
            INTEGER(C_INT)  :: page
            INTEGER(C_INT)  :: npages
            INTEGER(C_INT)  :: write_tiff_strip_long1
          END FUNCTION
          FUNCTION write_tiff_strip_float1(strip,bitsPerSample_,width,length,page,npages) bind(C,NAME='write_tiff_strip_float1')
            IMPORT         :: C_INT,C_FLOAT
            INTEGER(C_INT) :: width
            INTEGER(C_INT) :: length
            REAL(C_FLOAT)  :: strip(width,length)
            INTEGER(C_INT) :: bitsPerSample_
            INTEGER(C_INT) :: page
            INTEGER(C_INT) :: npages
            INTEGER(C_INT) :: write_tiff_strip_float1
          END FUNCTION
          FUNCTION write_tiff_strip_double1(strip,bitsPerSample_,width,length,page,npages) bind(C,NAME='write_tiff_strip_double1')
            IMPORT         :: C_INT,C_DOUBLE
            INTEGER(C_INT) :: width
            INTEGER(C_INT) :: length
            REAL(C_DOUBLE) :: strip(width,length)
            INTEGER(C_INT) :: bitsPerSample_
            INTEGER(C_INT) :: page
            INTEGER(C_INT) :: npages
            INTEGER(C_INT) :: write_tiff_strip_double1
          END FUNCTION
        END INTERFACE

        PUBLIC :: bitsPerSample
        PUBLIC :: bitsPerSampleW
        PUBLIC :: samplesPerPixel

        PUBLIC :: ppm_rc_open_tiff
        PUBLIC :: ppm_rc_open_bigtiff
        PUBLIC :: ppm_rc_open_write_tiff
        PUBLIC :: ppm_rc_open_write_bigtiff
        PUBLIC :: ppm_rc_close_tiff
        PUBLIC :: ppm_rc_write_tiff_header
        PUBLIC :: ppm_rc_read_tiff_scanline
        PUBLIC :: ppm_rc_read_tiff_strip
        PUBLIC :: ppm_rc_write_tiff_scanline
        PUBLIC :: ppm_rc_write_tiff_strip

        INTERFACE ppm_rc_read_tiff_info
          MODULE PROCEDURE read_tiff_info_f1
          MODULE PROCEDURE read_tiff_info_f2
        END INTERFACE

        PUBLIC :: ppm_rc_read_tiff_info

      CONTAINS
        INTEGER FUNCTION read_tiff_info_f1(filename,ngrid)
          IMPLICIT NONE
          CHARACTER(LEN=ppm_char), INTENT(IN   ) :: filename
          INTEGER, DIMENSION(:),   INTENT(INOUT) :: ngrid
          INTEGER, DIMENSION(3) :: ngrid_
          read_tiff_info_f1=read_tiff_infoC(filename,ngrid_,bitsPerSample,samplesPerPixel,ppm_rank,ppm_comm)
          ngrid(1:ppm_dim)=ngrid_(1:ppm_dim)
          IF (bitsPerSampleW.EQ.0) bitsPerSampleW=bitsPerSample
          RETURN
        END FUNCTION read_tiff_info_f1
        INTEGER FUNCTION read_tiff_info_f2(filename)
          IMPLICIT NONE
          CHARACTER(LEN=ppm_char), INTENT(IN   ) :: filename
          INTEGER, DIMENSION(3) :: ngrid_
          read_tiff_info_f2=read_tiff_infoC(filename,ngrid_,bitsPerSample,samplesPerPixel,ppm_rank,ppm_comm)
          IF (bitsPerSampleW.EQ.0) bitsPerSampleW=bitsPerSample
          RETURN
        END FUNCTION read_tiff_info_f2
#endif
#endif
      END MODULE ppm_rc_module_tiff
