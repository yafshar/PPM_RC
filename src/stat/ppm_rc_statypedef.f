        TYPE :: ppm_rc_stat
          PRIVATE
          REAL(ppm_kind_double), DIMENSION(:), POINTER :: value => NULL()
          !!! value stored in link
        CONTAINS
           PROCEDURE :: addValue
           PROCEDURE :: addValue3
           PROCEDURE :: addValue4

           ! add value to linked list
           PROCEDURE :: getValue
           !!! return value pointer
           PROCEDURE :: destroy
           !!! destroy the value
           GENERIC   :: add =>     &
           &            addValue,  &
           &            addValue3, &
           &            addValue4
        END TYPE ppm_rc_stat

        !collection of pointers to objects (removing an element
        !does not delete the object itself, just the pointer to the object)
        TYPE ppm_rc_t_ptr_stat
          TYPE(ppm_rc_stat), POINTER :: t => NULL()
        END TYPE ppm_rc_t_ptr_stat
        !
        TYPE, EXTENDS(ppm_t_container) :: ppm_rc_c_stat
          TYPE(ppm_rc_t_ptr_stat), DIMENSION(:), POINTER :: vec      => NULL()
!           TYPE(ppm_rc_stat),                     POINTER :: iterator => NULL()
        CONTAINS
          PROCEDURE :: begin      => ppm_rc_c_stat_begin
          PROCEDURE :: next       => ppm_rc_c_stat_next
          PROCEDURE :: last       => ppm_rc_c_stat_last
          PROCEDURE :: prev       => ppm_rc_c_stat_prev
          PROCEDURE :: destroy    => ppm_rc_c_stat_destroy
          PROCEDURE :: exists     => ppm_rc_c_stat_exists
          PROCEDURE :: has        => ppm_rc_c_stat_has
          PROCEDURE :: get_id     => ppm_rc_c_stat_get_id
          PROCEDURE :: push       => ppm_rc_c_stat_push
          PROCEDURE :: remove     => ppm_rc_c_stat_remove
          PROCEDURE :: grow_size  => ppm_rc_c_stat_grow_size
          PROCEDURE :: at         => ppm_rc_c_stat_at
        END TYPE ppm_rc_c_stat