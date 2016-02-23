      !-------------------------------------------------------------------------
      !  Subroutine   :                    ppm_rc_index_exist
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
      !  MOSAIC Group
      !  Max Planck Institute of Molecular Cell Biology and Genetics
      !  Pfotenhauerstr. 108, 01307 Dresden, Germany
      !
      !  Author           - y.afshar           Nov   2015
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
           ppm_rc_label_exist_dsc=label_exist(label,arrayoflabels,arraysize)
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
           ppm_rc_label_exist_dsc2=label_exist(label,arrayoflabels,arraysize,arrayindex)
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


