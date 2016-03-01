      !-------------------------------------------------------------------------
      !  Subroutine   :                    ppm_rc_label_exist
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
      !  Author           - y.afshar           Nov   2015
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Subroutine   :                    ppm_rc_label_exist
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Finding whether the particle exists given the array of particle
      !
      !
      !  Input        :
      !
      !  Input/output :
      !
      !  Output       : TRUE or FALSE
      !
      !  Routines     :
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      LOGICAL FUNCTION ppm_rc_label_exist_asc(label,arrayoflabels,arraysize)

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
        INTEGER,               INTENT(IN   ) :: label
        !!!
        INTEGER, DIMENSION(:), INTENT(IN   ) :: arrayoflabels
        !!! Sorted array of indices
        INTEGER,               INTENT(IN   ) :: arraysize
        !!!
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
           ppm_rc_label_exist_asc=.FALSE.
           RETURN
        ELSE IF (label.LE.0) THEN
           ppm_rc_label_exist_asc=.FALSE.
           RETURN
        ELSE IF (label.LT.arrayoflabels(1)) THEN
           ppm_rc_label_exist_asc=.FALSE.
           RETURN
        ELSE IF (label.EQ.arrayoflabels(1)) THEN
           ppm_rc_label_exist_asc=.TRUE.
           RETURN
        ELSE IF (label.EQ.arrayoflabels(arraysize)) THEN
           ppm_rc_label_exist_asc=.TRUE.
           RETURN
        ELSE IF (label.GT.arrayoflabels(arraysize)) THEN
           ppm_rc_label_exist_asc=.FALSE.
           RETURN
        ENDIF

        ! Limits for the search
        lower_part=1
        upper_part=arraysize

        DO WHILE (upper_part-lower_part.GT.1)
           mid_part=(lower_part+upper_part)/2
           IF      (label.LT.arrayoflabels(mid_part)) THEN
              upper_part=mid_part
           ELSE IF (label.GT.arrayoflabels(mid_part)) THEN
              lower_part=mid_part
           ELSE
              ppm_rc_label_exist_asc=.TRUE.
              RETURN
           ENDIF
        ENDDO

        ppm_rc_label_exist_asc=.FALSE.
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        RETURN
      END FUNCTION ppm_rc_label_exist_asc

      LOGICAL FUNCTION ppm_rc_label_exist_asc2(label,arrayoflabels,arraysize,arrayindex)

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
        INTEGER,               INTENT(IN   ) :: label
        !!!
        INTEGER, DIMENSION(:), INTENT(IN   ) :: arrayoflabels
        !!! Sorted array of indices
        INTEGER,               INTENT(IN   ) :: arraysize
        !!!
        INTEGER,               INTENT(  OUT) :: arrayindex
        !!!
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
           ppm_rc_label_exist_asc2=.FALSE.
           RETURN
        ELSE IF (label.LE.0) THEN
           ppm_rc_label_exist_asc2=.FALSE.
           RETURN
        ELSE IF (label.LT.arrayoflabels(1)) THEN
           ppm_rc_label_exist_asc2=.FALSE.
           RETURN
        ELSE IF (label.EQ.arrayoflabels(1)) THEN
           ppm_rc_label_exist_asc2=.TRUE.
           arrayindex=1
           RETURN
        ELSE IF (label.EQ.arrayoflabels(arraysize)) THEN
           ppm_rc_label_exist_asc2=.TRUE.
           arrayindex=arraysize
           RETURN
        ELSE IF (label.GT.arrayoflabels(arraysize)) THEN
           ppm_rc_label_exist_asc2=.FALSE.
           RETURN
        ENDIF

        ! Limits for the search
        lower_part=1
        upper_part=arraysize

        DO WHILE (upper_part-lower_part.GT.1)
           mid_part=(lower_part+upper_part)/2
           IF      (label.LT.arrayoflabels(mid_part)) THEN
              upper_part=mid_part
           ELSE IF (label.GT.arrayoflabels(mid_part)) THEN
              lower_part=mid_part
           ELSE
              ppm_rc_label_exist_asc2=.TRUE.
              arrayindex=mid_part
              RETURN
           ENDIF
        ENDDO

        ppm_rc_label_exist_asc2=.FALSE.
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        RETURN
      END FUNCTION ppm_rc_label_exist_asc2

      LOGICAL FUNCTION ppm_rc_label_exist_asc3(label,arrayoflabels,arrayofranks,arraysize)

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
        INTEGER,               INTENT(IN   ) :: label
        !!!
        INTEGER, DIMENSION(:), INTENT(IN   ) :: arrayoflabels
        !!! array of labels
        INTEGER, DIMENSION(:), POINTER       :: arrayofranks
        !!! array of indices of the Sorted array
        INTEGER,               INTENT(IN   ) :: arraysize
        !!!
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
        INTEGER :: lower_part
        INTEGER :: mid_part
        INTEGER :: upper_part
        INTEGER :: arrayindex
        !-------------------------------------------------------------------------
        !
        !-------------------------------------------------------------------------
        IF  (arraysize.LT.1) THEN
           ppm_rc_label_exist_asc3=.FALSE.
           RETURN
        ELSE IF (label.LE.0) THEN
           ppm_rc_label_exist_asc3=.FALSE.
           RETURN
        ELSE IF (label.LT.arrayoflabels(arrayofranks(1))) THEN
           ppm_rc_label_exist_asc3=.FALSE.
           RETURN
        ELSE IF (label.EQ.arrayoflabels(arrayofranks(1))) THEN
           ppm_rc_label_exist_asc3=.TRUE.
           RETURN
        ELSE IF (label.EQ.arrayoflabels(arrayofranks(arraysize))) THEN
           ppm_rc_label_exist_asc3=.TRUE.
           RETURN
        ELSE IF (label.GT.arrayoflabels(arrayofranks(arraysize))) THEN
           ppm_rc_label_exist_asc3=.FALSE.
           RETURN
        ENDIF

        ! Limits for the search
        lower_part=1
        upper_part=arraysize

        DO WHILE (upper_part-lower_part.GT.1)
           mid_part=(lower_part+upper_part)/2
           arrayindex=arrayofranks(mid_part)
           IF      (label.LT.arrayoflabels(arrayindex)) THEN
              upper_part=mid_part
           ELSE IF (label.GT.arrayoflabels(arrayindex)) THEN
              lower_part=mid_part
           ELSE
              ppm_rc_label_exist_asc3=.TRUE.
              RETURN
           ENDIF
        ENDDO

        ppm_rc_label_exist_asc3=.FALSE.
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        RETURN
      END FUNCTION ppm_rc_label_exist_asc3

      LOGICAL FUNCTION ppm_rc_label_exist_asc4(label,arrayoflabels,arrayofranks,arraysize,arrayindex)

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
        INTEGER,               INTENT(IN   ) :: label
        !!!
        INTEGER, DIMENSION(:), INTENT(IN   ) :: arrayoflabels
        !!! array of labels
        INTEGER, DIMENSION(:), POINTER       :: arrayofranks
        !!! array of indices of the sorted arrayoflabels
        INTEGER,               INTENT(IN   ) :: arraysize
        !!!
        INTEGER,               INTENT(  OUT) :: arrayindex
        !!!
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
           ppm_rc_label_exist_asc4=.FALSE.
           RETURN
        ELSE IF (label.LE.0) THEN
           ppm_rc_label_exist_asc4=.FALSE.
           RETURN
        ELSE IF (label.LT.arrayoflabels(arrayofranks(1))) THEN
           ppm_rc_label_exist_asc4=.FALSE.
           RETURN
        ELSE IF (label.EQ.arrayoflabels(arrayofranks(1))) THEN
           ppm_rc_label_exist_asc4=.TRUE.
           arrayindex=arrayofranks(1)
           RETURN
        ELSE IF (label.EQ.arrayoflabels(arrayofranks(arraysize))) THEN
           ppm_rc_label_exist_asc4=.TRUE.
           arrayindex=arrayofranks(arraysize)
           RETURN
        ELSE IF (label.GT.arrayoflabels(arrayofranks(arraysize))) THEN
           ppm_rc_label_exist_asc4=.FALSE.
           RETURN
        ENDIF

        ! Limits for the search
        lower_part=1
        upper_part=arraysize

        DO WHILE (upper_part-lower_part.GT.1)
           mid_part=(lower_part+upper_part)/2
           arrayindex=arrayofranks(mid_part)
           IF      (label.LT.arrayoflabels(arrayindex)) THEN
              upper_part=mid_part
           ELSE IF (label.GT.arrayoflabels(arrayindex)) THEN
              lower_part=mid_part
           ELSE
              ppm_rc_label_exist_asc4=.TRUE.
              RETURN
           ENDIF
        ENDDO

        ppm_rc_label_exist_asc4=.FALSE.
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        RETURN
      END FUNCTION ppm_rc_label_exist_asc4

      LOGICAL FUNCTION ppm_rc_label_exist_dsc(label,arrayoflabels,arraysize,descend)

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
        INTEGER,               INTENT(IN   ) :: label
        !!!
        INTEGER, DIMENSION(:), INTENT(IN   ) :: arrayoflabels
        !!! Sorted array of indices
        INTEGER,               INTENT(IN   ) :: arraysize
        !!!
        LOGICAL,               INTENT(IN   ) :: descend
        !!! array of indices is sorted in a descending order
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
           ppm_rc_label_exist_dsc=.FALSE.
           RETURN
        ELSE IF (label.LE.0) THEN
           ppm_rc_label_exist_dsc=.FALSE.
           RETURN
        ELSE IF (label.EQ.arrayoflabels(1)) THEN
           ppm_rc_label_exist_dsc=.TRUE.
           RETURN
        ELSE IF (label.EQ.arrayoflabels(arraysize)) THEN
           ppm_rc_label_exist_dsc=.TRUE.
           RETURN
        ELSE IF (.NOT.descend) THEN
           ppm_rc_label_exist_dsc=ppm_rc_label_exist(label,arrayoflabels,arraysize)
           RETURN
        ELSE IF (label.GT.arrayoflabels(1)) THEN
           ppm_rc_label_exist_dsc=.FALSE.
           RETURN
        ELSE IF (label.LT.arrayoflabels(arraysize)) THEN
           ppm_rc_label_exist_dsc=.FALSE.
           RETURN
        ENDIF

        ! Limits for the search
        lower_part=1
        upper_part=arraysize

        DO WHILE (upper_part-lower_part.GT.1)
           mid_part=(lower_part+upper_part)/2
           IF      (label.GT.arrayoflabels(mid_part)) THEN
              upper_part=mid_part
           ELSE IF (label.LT.arrayoflabels(mid_part)) THEN
              lower_part=mid_part
           ELSE
              ppm_rc_label_exist_dsc=.TRUE.
              RETURN
           ENDIF
        ENDDO

        ppm_rc_label_exist_dsc=.FALSE.
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        RETURN
      END FUNCTION ppm_rc_label_exist_dsc

      LOGICAL FUNCTION ppm_rc_label_exist_dsc2(label,arrayoflabels,arraysize,descend,arrayindex)

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
        INTEGER,               INTENT(IN   ) :: label
        !!!
        INTEGER, DIMENSION(:), INTENT(IN   ) :: arrayoflabels
        !!! Sorted array of indices
        INTEGER,               INTENT(IN   ) :: arraysize
        !!!
        LOGICAL,               INTENT(IN   ) :: descend
        !!! array of indices is sorted in a descending order
        INTEGER,               INTENT(  OUT) :: arrayindex
        !!!
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
           ppm_rc_label_exist_dsc2=.FALSE.
           RETURN
        ELSE IF (label.LE.0) THEN
           ppm_rc_label_exist_dsc2=.FALSE.
           RETURN
        ELSE IF (label.EQ.arrayoflabels(1)) THEN
           ppm_rc_label_exist_dsc2=.TRUE.
           arrayindex=1
           RETURN
        ELSE IF (label.EQ.arrayoflabels(arraysize)) THEN
           ppm_rc_label_exist_dsc2=.TRUE.
           arrayindex=arraysize
           RETURN
        ELSE IF (.NOT.descend) THEN
           ppm_rc_label_exist_dsc2=ppm_rc_label_exist(label,arrayoflabels,arraysize,arrayindex)
           RETURN
        ELSE IF (label.GT.arrayoflabels(1)) THEN
           ppm_rc_label_exist_dsc2=.FALSE.
           RETURN
        ELSE IF (label.LT.arrayoflabels(arraysize)) THEN
           ppm_rc_label_exist_dsc2=.FALSE.
           RETURN
        ENDIF

        ! Limits for the search
        lower_part=1
        upper_part=arraysize

        DO WHILE (upper_part-lower_part.GT.1)
           mid_part=(lower_part+upper_part)/2
           IF      (label.GT.arrayoflabels(mid_part)) THEN
              upper_part=mid_part
           ELSE IF (label.LT.arrayoflabels(mid_part)) THEN
              lower_part=mid_part
           ELSE
              ppm_rc_label_exist_dsc2=.TRUE.
              arrayindex=mid_part
              RETURN
           ENDIF
        ENDDO

        ppm_rc_label_exist_dsc2=.FALSE.
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        RETURN
      END FUNCTION ppm_rc_label_exist_dsc2

      LOGICAL FUNCTION ppm_rc_label_exist_dsc3(label,arrayoflabels,arrayofranks,arraysize,descend)

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
        INTEGER,               INTENT(IN   ) :: label
        !!!
        INTEGER, DIMENSION(:), INTENT(IN   ) :: arrayoflabels
        !!! Sorted array of indices
        INTEGER, DIMENSION(:), POINTER       :: arrayofranks
        !!!
        INTEGER,               INTENT(IN   ) :: arraysize
        !!!
        LOGICAL,               INTENT(IN   ) :: descend
        !!! array of indices is sorted in a descending order
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
        INTEGER :: lower_part
        INTEGER :: mid_part
        INTEGER :: upper_part
        INTEGER :: arrayindex
        !-------------------------------------------------------------------------
        !
        !-------------------------------------------------------------------------
        IF  (arraysize.LT.1) THEN
           ppm_rc_label_exist_dsc3=.FALSE.
           RETURN
        ELSE IF (label.LE.0) THEN
           ppm_rc_label_exist_dsc3=.FALSE.
           RETURN
        ELSE IF (label.EQ.arrayoflabels(arrayofranks(1))) THEN
           ppm_rc_label_exist_dsc3=.TRUE.
           RETURN
        ELSE IF (label.EQ.arrayoflabels(arrayofranks(arraysize))) THEN
           ppm_rc_label_exist_dsc3=.TRUE.
           RETURN
        ELSE IF (.NOT.descend) THEN
           ppm_rc_label_exist_dsc3=ppm_rc_label_exist(label,arrayoflabels,arrayofranks,arraysize)
           RETURN
        ELSE IF (label.GT.arrayoflabels(arrayofranks(1))) THEN
           ppm_rc_label_exist_dsc3=.FALSE.
           RETURN
        ELSE IF (label.LT.arrayoflabels(arrayofranks(arraysize))) THEN
           ppm_rc_label_exist_dsc3=.FALSE.
           RETURN
        ENDIF

        ! Limits for the search
        lower_part=1
        upper_part=arraysize

        DO WHILE (upper_part-lower_part.GT.1)
           mid_part=(lower_part+upper_part)/2
           arrayindex=arrayofranks(mid_part)
           IF      (label.GT.arrayoflabels(arrayindex)) THEN
              upper_part=mid_part
           ELSE IF (label.LT.arrayoflabels(arrayindex)) THEN
              lower_part=mid_part
           ELSE
              ppm_rc_label_exist_dsc3=.TRUE.
              RETURN
           ENDIF
        ENDDO

        ppm_rc_label_exist_dsc3=.FALSE.
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        RETURN
      END FUNCTION ppm_rc_label_exist_dsc3

      LOGICAL FUNCTION ppm_rc_label_exist_dsc4(label,arrayoflabels,arrayofranks,arraysize,descend,arrayindex)

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
        INTEGER,               INTENT(IN   ) :: label
        !!!
        INTEGER, DIMENSION(:), INTENT(IN   ) :: arrayoflabels
        !!! Sorted array of indices
        INTEGER, DIMENSION(:), POINTER       :: arrayofranks
        !!!
        INTEGER,               INTENT(IN   ) :: arraysize
        !!!
        LOGICAL,               INTENT(IN   ) :: descend
        !!! array of indices is sorted in a descending order
        INTEGER,               INTENT(  OUT) :: arrayindex
        !!!
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
           ppm_rc_label_exist_dsc4=.FALSE.
           RETURN
        ELSE IF (label.LE.0) THEN
           ppm_rc_label_exist_dsc4=.FALSE.
           RETURN
        ELSE IF (label.EQ.arrayoflabels(arrayofranks(1))) THEN
           ppm_rc_label_exist_dsc4=.TRUE.
           arrayindex=arrayofranks(1)
           RETURN
        ELSE IF (label.EQ.arrayoflabels(arrayofranks(arraysize))) THEN
           ppm_rc_label_exist_dsc4=.TRUE.
           arrayindex=arrayofranks(arraysize)
           RETURN
        ELSE IF (.NOT.descend) THEN
           ppm_rc_label_exist_dsc4=ppm_rc_label_exist(label,arrayoflabels,arrayofranks,arraysize,arrayindex)
           RETURN
        ELSE IF (label.GT.arrayoflabels(arrayofranks(1))) THEN
           ppm_rc_label_exist_dsc4=.FALSE.
           RETURN
        ELSE IF (label.LT.arrayoflabels(arrayofranks(arraysize))) THEN
           ppm_rc_label_exist_dsc4=.FALSE.
           RETURN
        ENDIF

        ! Limits for the search
        lower_part=1
        upper_part=arraysize

        DO WHILE (upper_part-lower_part.GT.1)
           mid_part=(lower_part+upper_part)/2
           arrayindex=arrayofranks(mid_part)
           IF      (label.GT.arrayoflabels(arrayindex)) THEN
              upper_part=mid_part
           ELSE IF (label.LT.arrayoflabels(arrayindex)) THEN
              lower_part=mid_part
           ELSE
              ppm_rc_label_exist_dsc4=.TRUE.
              RETURN
           ENDIF
        ENDDO

        ppm_rc_label_exist_dsc4=.FALSE.
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        RETURN
      END FUNCTION ppm_rc_label_exist_dsc4
