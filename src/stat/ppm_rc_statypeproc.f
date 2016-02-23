        FUNCTION getValue(this)
          CLASS(ppm_rc_stat)                           :: this

          REAL(ppm_kind_double), DIMENSION(:), POINTER :: getValue

          getValue => this%value
        END FUNCTION getValue

        SUBROUTINE addValue(this,value)
          !!! Add a value to a link (by creating a link) and the link to the list
          CLASS(ppm_rc_stat)                                 :: this

          REAL(ppm_kind_double), DIMENSION(:), INTENT(IN   ) :: value

          IF (.NOT.ASSOCIATED(this%value)) THEN
             ALLOCATE(this%value(4),SOURCE=value)
          ELSE
             this%value(2)=this%value(2)+value(1)
             this%value(3)=this%value(3)+value(2)
             this%value(4)=this%value(4)+value(3)
          ENDIF
        END SUBROUTINE addValue

        SUBROUTINE addValue3(this,value1,value2,value3)
          !!! Add a value to a link (by creating a link) and the link to the list
          CLASS(ppm_rc_stat)                   :: this

          REAL(ppm_kind_double), INTENT(IN   ) :: value1
          REAL(ppm_kind_double), INTENT(IN   ) :: value2
          REAL(ppm_kind_double), INTENT(IN   ) :: value3

          INTEGER :: info

          IF (.NOT.ASSOCIATED(this%value)) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_argument, "addValue", &
             &    'ASSOCIATED(this%value) is not true!', __LINE__, info)
          ELSE
             this%value(2)=this%value(2)+value1
             this%value(3)=this%value(3)+value2
             this%value(4)=this%value(4)+value3
          ENDIF
        END SUBROUTINE addValue3

        SUBROUTINE addValue4(this,value1,value2,value3,value4)
          !!! Add a value to a link (by creating a link) and the link to the list
          CLASS(ppm_rc_stat)                   :: this

          REAL(ppm_kind_double), INTENT(IN   ) :: value1
          REAL(ppm_kind_double), INTENT(IN   ) :: value2
          REAL(ppm_kind_double), INTENT(IN   ) :: value3
          REAL(ppm_kind_double), INTENT(IN   ) :: value4

          INTEGER :: info

          IF (.NOT.ASSOCIATED(this%value)) THEN
             ALLOCATE(this%value(4))
             this%value(1)=value1
             this%value(2)=value2
             this%value(3)=value3
             this%value(4)=value4
          ELSE
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_argument, "addValue", &
             &    'ASSOCIATED(this%value) is true but we are adding four values!', __LINE__, info)
          ENDIF

        END SUBROUTINE addValue4

        SUBROUTINE destroy(this)
          !!! Destroy a list
          CLASS(ppm_rc_stat) :: this

          IF (ASSOCIATED(this%value)) DEALLOCATE(this%value)
          NULLIFY(this%value)

        END SUBROUTINE destroy

        !BEGIN
        FUNCTION ppm_rc_c_stat_begin(this) RESULT (iterator)

            IMPLICIT NONE
            !-------------------------------------------------------------------------
            !  Arguments
            !-------------------------------------------------------------------------
            CLASS(ppm_rc_c_stat), INTENT(INOUT) :: this

            TYPE(ppm_rc_stat),    POINTER       :: iterator

            IF (this%nb.GT.0) THEN
               this%iter_id = this%min_id
               IF (this%iter_id.LE.this%max_id) THEN
                  iterator => this%vec(this%min_id)%t
                  RETURN
               ENDIF
            ENDIF
            iterator => NULL()
        END FUNCTION ppm_rc_c_stat_begin
        !NEXT
        FUNCTION ppm_rc_c_stat_next(this) RESULT (iterator)

            IMPLICIT NONE
            !-------------------------------------------------------------------------
            !  Arguments
            !-------------------------------------------------------------------------
            CLASS(ppm_rc_c_stat), INTENT(INOUT) :: this

            TYPE(ppm_rc_stat),    POINTER       :: iterator

            IF (this%nb.GT.0) THEN
               this%iter_id = this%iter_id + 1
               IF (this%iter_id.GE.this%min_id .AND. this%iter_id.LE.this%max_id) THEN
                  iterator => this%vec(this%iter_id)%t
                  RETURN
               ENDIF
            ENDIF
            iterator => NULL()
        END FUNCTION ppm_rc_c_stat_next
        !PREVIOUS
        FUNCTION ppm_rc_c_stat_prev(this) RESULT (iterator)

            IMPLICIT NONE
            !-------------------------------------------------------------------------
            !  Arguments
            !-------------------------------------------------------------------------
            CLASS(ppm_rc_c_stat), INTENT(INOUT) :: this

            TYPE(ppm_rc_stat),    POINTER       :: iterator

            IF (this%nb.GT.0) THEN
               this%iter_id = this%iter_id - 1
               IF (this%iter_id.GE.this%min_id .AND. this%iter_id.LE.this%max_id) THEN
                  iterator => this%vec(this%iter_id)%t
                  RETURN
               ENDIF
            ENDIF
            iterator => NULL()
        END FUNCTION ppm_rc_c_stat_prev
        !LAST
        FUNCTION ppm_rc_c_stat_last(this) RESULT (iterator)

            IMPLICIT NONE
            !-------------------------------------------------------------------------
            !  Arguments
            !-------------------------------------------------------------------------
            CLASS(ppm_rc_c_stat), INTENT(INOUT) :: this

            TYPE(ppm_rc_stat),    POINTER       :: iterator

            IF (this%nb.GT.0) THEN
               this%iter_id = this%max_id
               IF (this%iter_id.GE.this%min_id) THEN
                  iterator => this%vec(this%max_id)%t
                  RETURN
               ENDIF
            ENDIF
            iterator => NULL()
        END FUNCTION ppm_rc_c_stat_last
        !DESTROY CONTAINER
        SUBROUTINE ppm_rc_c_stat_destroy(this,info)

            IMPLICIT NONE
            !-------------------------------------------------------------------------
            !  Arguments
            !-------------------------------------------------------------------------
            CLASS(ppm_rc_c_stat), INTENT(INOUT) :: this

            INTEGER,              INTENT(  OUT) :: info

            !-------------------------------------------------------------------------
            ! local variables
            !-------------------------------------------------------------------------
            TYPE(ppm_rc_stat), POINTER :: p

            start_subroutine("ppm_rc_c_stat_destroy")

            p => this%begin()
            DO WHILE(ASSOCIATED(p))
               !yaser I think this is necessary for complete destroyment
               CALL p%destroy()

               DEALLOCATE(P,STAT=info)
               or_fail_dealloc("Could not deallocate collection element")

               p => this%next()
            ENDDO

            IF (ASSOCIATED(this%vec)) DEALLOCATE(this%vec,STAT=info)
            or_fail_dealloc("Could not deallocate collection array")

            this%iter_id = 0
            this%min_id = 0
            this%max_id = 0
            this%nb = 0
            this%vec_size=0

            end_subroutine()
        END SUBROUTINE ppm_rc_c_stat_destroy

        !EXISTS
        FUNCTION ppm_rc_c_stat_exists(this,id) RESULT(exists)
            !!! Check whether an element exists and can be accessed at this id

            IMPLICIT NONE
            !-------------------------------------------------------------------------
            !  Arguments
            !-------------------------------------------------------------------------
            CLASS(ppm_rc_c_stat), INTENT(IN   ) :: this
            !!! Data structure containing the particles
            INTEGER,              INTENT(IN   ) :: id
            !!! id where the data is stored
            LOGICAL                             :: exists
            !!! Return status, on success 0.

            !-------------------------------------------------------------------------
            ! local variables
            !-------------------------------------------------------------------------

            !-------------------------------------------------------------------------
            ! Check arguments
            !-------------------------------------------------------------------------
            IF (id.LE.0 .OR. id.LT.this%min_id .OR. id.GT.this%max_id) THEN
               exists = .FALSE.
               RETURN
            ENDIF
            IF (.NOT.ASSOCIATED(this%vec)) THEN
               exists = .FALSE.
               RETURN
            ENDIF
            IF (.NOT.ASSOCIATED(this%vec(id)%t)) THEN
               exists = .FALSE.
               RETURN
            ENDIF
            exists = .TRUE.

            RETURN
        END FUNCTION ppm_rc_c_stat_exists

        !HAS
        FUNCTION ppm_rc_c_stat_has(this,element) RESULT(has)
            !!! Check whether an element is present in the collection (slow...)

            IMPLICIT NONE
            !-------------------------------------------------------------------------
            !  Arguments
            !-------------------------------------------------------------------------
            CLASS(ppm_rc_c_stat),       INTENT(INOUT) :: this
            !!! Data structure containing the particles
            TYPE(ppm_rc_stat),  TARGET, INTENT(IN   ) :: element
            !!! element which is being searched for
            LOGICAL                                   :: has
            !!! true if this element belongs to the collection

            !-------------------------------------------------------------------------
            ! local variables
            !-------------------------------------------------------------------------
            TYPE(ppm_rc_stat), POINTER :: p

            has = .TRUE.

            p => this%begin()
            DO WHILE(ASSOCIATED(p))
               IF (ASSOCIATED(p,element)) RETURN
               p => this%next()
            ENDDO

            has = .FALSE.
            RETURN

        END FUNCTION ppm_rc_c_stat_has

        !GET_ID
        FUNCTION ppm_rc_c_stat_get_id(this,element) RESULT(id)
            !!! Returns the id of an element in the collection (slow...)
            !!! Returs -1 if the element is not found.

            IMPLICIT NONE
            !-------------------------------------------------------------------------
            !  Arguments
            !-------------------------------------------------------------------------
            CLASS(ppm_rc_c_stat),      INTENT(INOUT) :: this
            !!! Data structure containing the particles
            TYPE(ppm_rc_stat), TARGET, INTENT(IN   ) :: element
            !!! Element which is being searched for
            INTEGER                                  :: id
            !!! id where the data is stored

            !-------------------------------------------------------------------------
            ! local variables
            !-------------------------------------------------------------------------
            TYPE(ppm_rc_stat), POINTER :: p

            p => this%begin()
            DO WHILE(ASSOCIATED(p))
               IF (ASSOCIATED(p,element)) THEN
                  id = this%iter_id
                  RETURN
               ENDIF
               p => this%next()
            ENDDO

            id = -1
            RETURN

        END FUNCTION ppm_rc_c_stat_get_id

        !PUSH
        SUBROUTINE ppm_rc_c_stat_push(this,element,info,id)
            !!! add an element into the collection

            IMPLICIT NONE
            !-------------------------------------------------------------------------
            !  Arguments
            !-------------------------------------------------------------------------
            CLASS(ppm_rc_c_stat),           INTENT(INOUT) :: this

            TYPE(ppm_rc_stat),    POINTER                 :: element

            INTEGER,                        INTENT(  OUT) :: info
            INTEGER,              OPTIONAL, INTENT(  OUT) :: id
            !!! index of the element in the collection

            !-------------------------------------------------------------------------
            ! local variables
            !-------------------------------------------------------------------------
            start_subroutine("ppm_rc_c_stat_push")

            !add the element at the end of the array
            this%min_id = 1
            this%max_id = this%max_id + 1
            this%nb = this%nb + 1
            IF (PRESENT(id)) id = this%max_id

            IF (this%max_id.GT.this%vec_size) THEN
               CALL this%grow_size(info)
               or_fail("could not grow ppm_rc_c_stat to a larger size")
            ENDIF

            IF (ASSOCIATED(this%vec(this%max_id)%t)) THEN
               fail("Pointer at position of new element is already associated. Something wrong in the Collection data structure")
            ENDIF

            this%vec(this%max_id)%t => element

            check_associated_noscope(<#this%vec(this%max_id)%t#>,"Pushing element into collection failed unexpectedly")

            element => NULL()

            end_subroutine()
        END SUBROUTINE ppm_rc_c_stat_push
        !REMOVE
        SUBROUTINE ppm_rc_c_stat_remove(this,info,element)
            !!! If element is present, remove it from the collection
            !!! else, remove the current element (as defined by the iterator pointer)

            IMPLICIT NONE
            !-------------------------------------------------------------------------
            !  Arguments
            !-------------------------------------------------------------------------
            CLASS(ppm_rc_c_stat),        INTENT(INOUT) :: this

            INTEGER,                     INTENT(  OUT) :: info

            TYPE(ppm_rc_stat), OPTIONAL, INTENT(INOUT) :: element

            !-------------------------------------------------------------------------
            ! local variables
            !-------------------------------------------------------------------------
            INTEGER :: del_id
            INTEGER :: iter_id_save

            iter_id_save = this%iter_id

            IF (PRESENT(element)) THEN
               del_id = this%get_id(element)
            ELSE
               del_id = this%iter_id
            ENDIF

            !deallocate the element
            CALL this%vec(del_id)%t%destroy()

            !swap with the last non-empty element of the collection
            IF (this%max_id.GT.this%min_id) THEN
               this%vec(del_id)%t => this%vec(this%max_id)%t
               this%vec(this%max_id)%t => NULL()
            ELSE
               this%vec(del_id)%t => NULL()
            ENDIF

            this%nb = this%nb - 1
            this%max_id = this%max_id - 1
            this%iter_id = iter_id_save - 1
            IF (this%nb.EQ.0 .OR. this%max_id.EQ.0) this%min_id = 0

            info=0
        END SUBROUTINE ppm_rc_c_stat_remove
        !GROW COLLECTION SIZE
        SUBROUTINE ppm_rc_c_stat_grow_size(this,info)
            !!! Reallocate collection to a larger size
            !!! (twice the size seems good)

            IMPLICIT NONE
            !-------------------------------------------------------------------------
            !  Arguments
            !-------------------------------------------------------------------------
            CLASS(ppm_rc_c_stat), INTENT(INOUT) :: this

            INTEGER,              INTENT(  OUT) :: info

            !-------------------------------------------------------------------------
            !  Local variables
            !-------------------------------------------------------------------------
            TYPE(ppm_rc_t_ptr_stat), DIMENSION(:), POINTER :: vec_temp

            INTEGER :: i

            start_subroutine("ppm_rc_c_stat_grow_size")

            IF (this%vec_size.LE.0 .OR. .NOT.ASSOCIATED(this%vec)) THEN
                !if the array is empty, allocate with a reasonable size
                this%vec_size = 10
                ALLOCATE(this%vec(this%vec_size),STAT=info)
                or_fail_alloc("could not allocate collection array")
            ELSE
                !if the array is full, double its size
                ! allocate a temporary array to store the current element
                ALLOCATE(vec_temp(this%vec_size),STAT=info)
                or_fail_alloc("could not allocate temporary array")

                ! copy the elements to the temporary array
                DO i=1,this%vec_size
                   vec_temp(i)%t => this%vec(i)%t
                ENDDO

                !reallocate the collection to a larger size
                DEALLOCATE(this%vec,STAT=info)
                or_fail_dealloc("could not deallocate collection array")

                ALLOCATE(this%vec(2*this%vec_size),STAT=info)
                or_fail_alloc("could not allocate collection array")

                !copy the elements back from the temporary array to the collection
                DO i=1,this%vec_size
                   this%vec(i)%t => vec_temp(i)%t
                ENDDO
                this%vec_size = 2 * this%vec_size

                DEALLOCATE(vec_temp,STAT=info)
                or_fail_dealloc("could not deallocate temporary array")
                NULLIFY(vec_temp)
            ENDIF

            end_subroutine()
        END SUBROUTINE ppm_rc_c_stat_grow_size
        !AT
        FUNCTION ppm_rc_c_stat_at(this,i) RESULT (element)
            !!! Access the i-th element of the collection
            !!! (the ordering is the same as that given by the iterators next()
            !!! begin(), prev() and last() )
            !!! If there is less than i elements in the collection, the function
            !!! returns a NULL pointer.

            IMPLICIT NONE
            !-------------------------------------------------------------------------
            !  Arguments
            !-------------------------------------------------------------------------
            CLASS(ppm_rc_c_stat), INTENT(INOUT) :: this

            INTEGER,              INTENT(IN   ) :: i

            TYPE(ppm_rc_stat),    POINTER       :: element

            !-------------------------------------------------------------------------
            !  Local variables
            !-------------------------------------------------------------------------
            IF (i.GE.this%min_id .AND. i.LE.this%max_id) THEN
               element => this%vec(i)%t
               this%iter_id=i
               DO WHILE (.NOT.ASSOCIATED(element))
                  IF (this%iter_id.GT.this%max_id) THEN
                     element => NULL()
                     RETURN
                  ENDIF
                  element => this%next()
               ENDDO
            ELSE
               element => NULL()
            ENDIF
            RETURN
        END FUNCTION ppm_rc_c_stat_at
