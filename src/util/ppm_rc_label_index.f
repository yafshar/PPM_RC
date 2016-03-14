      !-------------------------------------------------------------------------
      !  Subroutine   :                    ppm_rc_label_index
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
      !  Author           - y.afshar           Feb   2016
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Subroutine   :                    ppm_rc_label_index
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Finding the newlabel index in the sorted array given the
      !                 array and the newlabel
      !
      !
      !  Input        :
      !
      !  Input/output :
      !
      !  Output       : index where the newlabel can sit in the sorted array
      !
      !  Routines     :
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      INTEGER FUNCTION ppm_rc_label_index(newlabel,arrayoflabels,arraysize)

        !-------------------------------------------------------------------------
        !  Modules
        !-------------------------------------------------------------------------
        IMPLICIT NONE

        !-------------------------------------------------------------------------
        !  Includes
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        INTEGER,               INTENT(IN   ) :: newlabel
        !!! newlabel we want to find the index for this in a sorted array
        INTEGER, DIMENSION(:), INTENT(IN   ) :: arrayoflabels
        !!! Sorted array of indices
        INTEGER,               INTENT(IN   ) :: arraysize
        !!! Size of the arrayoflabels
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
        INTEGER :: lower_part
        INTEGER :: mid_part
        INTEGER :: upper_part
        !-------------------------------------------------------------------------
        !
        !-------------------------------------------------------------------------
        IF  (arraysize.LT.1) THEN
           ppm_rc_label_index=0
           RETURN
        ELSE IF (newlabel.LE.0) THEN
           ppm_rc_label_index=0
           RETURN
        ELSE IF (newlabel.LE.arrayoflabels(1)) THEN
           ppm_rc_label_index=1
           RETURN
        ELSE IF (newlabel.GE.arrayoflabels(arraysize)) THEN
           ppm_rc_label_index=arraysize+1
           RETURN
        ENDIF

        ! Limits for the search
        lower_part=1
        upper_part=arraysize

        DO WHILE (upper_part-lower_part.GT.1)
           mid_part=(lower_part+upper_part)/2
           IF      (newlabel.LT.arrayoflabels(mid_part)) THEN
              upper_part=mid_part
           ELSE IF (newlabel.GT.arrayoflabels(mid_part)) THEN
              lower_part=mid_part
           ELSE
              ppm_rc_label_index=mid_part
              RETURN
           ENDIF
        ENDDO

        mid_part=(lower_part+upper_part)/2
        ppm_rc_label_index=mid_part+1
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        RETURN
      END FUNCTION ppm_rc_label_index
