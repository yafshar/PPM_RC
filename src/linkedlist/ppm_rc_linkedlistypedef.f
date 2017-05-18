        TYPE :: ppm_rc_link
          PRIVATE
          INTEGER, DIMENSION(:), POINTER :: value => NULL()
          !!! value stored in link
          TYPE(ppm_rc_link),     POINTER :: next => NULL()
          !!! next link in the list
        CONTAINS
          PROCEDURE :: getValue_
          PROCEDURE :: getValueMCMCParticle_2d
          PROCEDURE :: getValueMCMCParticle_3d
          PROCEDURE :: getValueMCMCHistoryParticle_2d
          PROCEDURE :: getValueMCMCHistoryParticle_3d
          GENERIC   :: getValue =>                     &
          &            getValue_,                      &
          &            getValueMCMCParticle_2d,        &
          &            getValueMCMCParticle_3d,        &
          &            getValueMCMCHistoryParticle_2d, &
          &            getValueMCMCHistoryParticle_3d
          !!! return value pointer
          PROCEDURE :: nextLink
          !!! return next pointer
          PROCEDURE :: setNextLink
          !!! set next pointer
        END TYPE ppm_rc_link

        INTERFACE ppm_rc_link
          MODULE PROCEDURE constructor
          MODULE PROCEDURE constructor1
          MODULE PROCEDURE constructor2
          MODULE PROCEDURE constructor3
          MODULE PROCEDURE constructor4
          MODULE PROCEDURE constructor5
          MODULE PROCEDURE constructor6
          MODULE PROCEDURE constructor7
          MODULE PROCEDURE constructor8
          MODULE PROCEDURE constructorMCMCParticle_2d
          MODULE PROCEDURE constructorMCMCParticle_3d
          MODULE PROCEDURE constructorMCMCHistoryParticle_2d
          MODULE PROCEDURE constructorMCMCHistoryParticle_3d
          !!! construct/initialize a link
        END INTERFACE

        TYPE ppm_rc_list
           TYPE(ppm_rc_link), POINTER :: first => NULL()
           !!! first link in list
           TYPE(ppm_rc_link), POINTER :: last  => NULL()
           !!! last link in list
!            TYPE(ppm_rc_link), POINTER :: current => NULL()
!            !!! current link in list
         CONTAINS
           PROCEDURE :: addValue
           PROCEDURE :: addValue1
           PROCEDURE :: addValue2
           PROCEDURE :: addValue3
           PROCEDURE :: addValue4
           PROCEDURE :: addValue5
           PROCEDURE :: addValue6
           PROCEDURE :: addValue7
           PROCEDURE :: addValue8
           PROCEDURE :: addValueMCMCParticle_2d
           PROCEDURE :: addValueMCMCParticle_3d
           PROCEDURE :: addValueMCMCHistoryParticle_2d
           PROCEDURE :: addValueMCMCHistoryParticle_3d
           GENERIC   :: add    =>                       &
           &            addValue,                       &
           &            addValue1,                      &
           &            addValue2,                      &
           &            addValue3,                      &
           &            addValue4,                      &
           &            addValue5,                      &
           &            addValue6,                      &
           &            addValue7,                      &
           &            addValue8,                      &
           &            addValueMCMCParticle_2d,        &
           &            addValueMCMCParticle_3d,        &
           &            addValueMCMCHistoryParticle_2d, &
           &            addValueMCMCHistoryParticle_3d
           ! add value to linked list
           PROCEDURE :: remove => removeLink
           ! remove a link from the list
           PROCEDURE :: destroy
           ! destroy the linked list
           PROCEDURE :: merge => mergeLinks
           ! swap a root into the boundary
           PROCEDURE :: swap1
           PROCEDURE :: swap2
           GENERIC   :: swap => swap1,swap2
        END TYPE ppm_rc_list

        !collection of pointers to objects (removing an element
        !does not delete the object itself, just the pointer to the object)
        TYPE ppm_rc_t_ptr_list
          TYPE(ppm_rc_list), POINTER :: t => NULL()
        END TYPE ppm_rc_t_ptr_list
        !
        TYPE, EXTENDS(ppm_t_container) :: ppm_rc_c_list
          TYPE(ppm_rc_t_ptr_list), DIMENSION(:), POINTER :: vec => NULL()
        CONTAINS
          PROCEDURE :: begin      => ppm_rc_c_list_begin
          PROCEDURE :: next       => ppm_rc_c_list_next
          PROCEDURE :: last       => ppm_rc_c_list_last
          PROCEDURE :: prev       => ppm_rc_c_list_prev
          PROCEDURE :: destroy    => ppm_rc_c_list_destroy
          PROCEDURE :: exists     => ppm_rc_c_list_exists
          PROCEDURE :: has        => ppm_rc_c_list_has
          PROCEDURE :: get_id     => ppm_rc_c_list_get_id
          PROCEDURE :: push       => ppm_rc_c_list_push
          PROCEDURE :: remove     => ppm_rc_c_list_remove
          PROCEDURE :: grow_size  => ppm_rc_c_list_grow_size
          PROCEDURE :: at         => ppm_rc_c_list_at
        END TYPE ppm_rc_c_list



        !!! This type will only get one INTEGER value as an input
        !!! In case of two or three inputs, they are hashed to one value
        !!! Using utility IndexHashFunctor, which is useful for converting
        !!! x,y,z to one INTEGER to save memory
        TYPE :: ppm_rc_link_
          PRIVATE

          INTEGER                     :: value
          !!! value stored in link
          TYPE(ppm_rc_link_), POINTER :: next => NULL()
          !!! next link in the list
        CONTAINS
          PROCEDURE :: getValue    => getValue__
          !!! return value pointer
          PROCEDURE :: getXYValue  => getXYValue_
          PROCEDURE :: getXYZValue => getXYZValue_
          !!! return value pointer
          PROCEDURE :: nextLink    => nextLink_
          !!! return next pointer
          PROCEDURE :: setNextLink => setNextLink_
          !!! set next pointer
        END TYPE ppm_rc_link_

        INTERFACE ppm_rc_link_
          MODULE PROCEDURE constructor_
          MODULE PROCEDURE constructor_2d
          MODULE PROCEDURE constructor_3d
          !!! construct/initialize a link
        END INTERFACE

        TYPE ppm_rc_list_
           TYPE(ppm_rc_link_), POINTER :: first => NULL()
           !!! first link in list
           TYPE(ppm_rc_link_), POINTER :: last  => NULL()
           !!! last link in list
!            TYPE(ppm_rc_link_), POINTER :: current => NULL()
!            !!! current link in list
         CONTAINS
           PROCEDURE :: addValue_
           PROCEDURE :: addValue__
           PROCEDURE :: addValue_2d
           PROCEDURE :: addValue_3d
           GENERIC   :: add     =>   &
           &            addValue_,   &
           &            addValue__,  &
           &            addValue_2d, &
           &            addValue_3d
           ! add value to linked list
           PROCEDURE :: remove => removeLink_
           ! remove a link from the list
           PROCEDURE :: destroy => destroy_
           ! destroy the linked list
           PROCEDURE :: merge => mergeLinks_
           ! swap a root into the boundary
           PROCEDURE :: swap1_
           PROCEDURE :: swap2_
           GENERIC   :: swap    => swap1_,swap2_
        END TYPE ppm_rc_list_

        !collection of pointers to objects (removing an element
        !does not delete the object itself, just the pointer to the object)
        TYPE ppm_rc_t_ptr_list_
          TYPE(ppm_rc_list_), POINTER :: t => NULL()
        END TYPE ppm_rc_t_ptr_list_
        !
        TYPE, EXTENDS(ppm_t_container) :: ppm_rc_c_list_
          TYPE(ppm_rc_t_ptr_list_), DIMENSION(:), POINTER :: vec => NULL()
        CONTAINS
          PROCEDURE :: begin      => ppm_rc_c_list_begin_
          PROCEDURE :: next       => ppm_rc_c_list_next_
          PROCEDURE :: last       => ppm_rc_c_list_last_
          PROCEDURE :: prev       => ppm_rc_c_list_prev_
          PROCEDURE :: destroy    => ppm_rc_c_list_destroy_
          PROCEDURE :: exists     => ppm_rc_c_list_exists_
          PROCEDURE :: has        => ppm_rc_c_list_has_
          PROCEDURE :: get_id     => ppm_rc_c_list_get_id_
          PROCEDURE :: push       => ppm_rc_c_list_push_
          PROCEDURE :: remove     => ppm_rc_c_list_remove_
          PROCEDURE :: grow_size  => ppm_rc_c_list_grow_size_
          PROCEDURE :: at         => ppm_rc_c_list_at_
        END TYPE ppm_rc_c_list_

        !!! This type will only get one INTEGER value as an input
        !!! In case of two or three inputs, they are iteration number and label index
        TYPE :: ppm_rc_link_2
          PRIVATE
          INTEGER(ppm_kind_int64)      :: Iteration
          INTEGER                      :: Label
          !!! value stored in link
          TYPE(ppm_rc_link_2), POINTER :: next => NULL()
          !!! next link in the list
        CONTAINS
          PROCEDURE :: getValue    => getValue_2
          !!! return value pointer
          PROCEDURE :: nextLink    => nextLink_2
          !!! return next pointer
          PROCEDURE :: setNextLink => setNextLink_2
          !!! set next pointer
        END TYPE ppm_rc_link_2

        INTERFACE ppm_rc_link_2
          MODULE PROCEDURE constructor_2
          !!! construct/initialize a link
        END INTERFACE

        TYPE ppm_rc_list_2
           TYPE(ppm_rc_link_2), POINTER :: first => NULL()
           !!! first link in list
           TYPE(ppm_rc_link_2), POINTER :: last  => NULL()
           !!! last link in list
         CONTAINS
           PROCEDURE :: addValue_2
           ! add value to linked list
           PROCEDURE :: removeLink_2
           ! remove a link from the list
           PROCEDURE :: destroy_2
           ! destroy the linked list
           PROCEDURE :: mergeLinks_2
           ! swap a root into the boundary
           PROCEDURE :: swap1_2
           PROCEDURE :: swap2_2

           GENERIC   :: destroy => destroy_2
           GENERIC   :: add     => addValue_2
           GENERIC   :: remove  => removeLink_2
           GENERIC   :: merge   => mergeLinks_2
           GENERIC   :: swap    => swap1_2,swap2_2
        END TYPE ppm_rc_list_2

        !collection of pointers to objects (removing an element
        !does not delete the object itself, just the pointer to the object)
        TYPE ppm_rc_t_ptr_list_2
          TYPE(ppm_rc_list_2), POINTER :: t => NULL()
        END TYPE ppm_rc_t_ptr_list_2
        !
        TYPE, EXTENDS(ppm_t_container) :: ppm_rc_c_list_2
          TYPE(ppm_rc_t_ptr_list_2), DIMENSION(:), POINTER :: vec => NULL()
        CONTAINS
          PROCEDURE :: begin      => ppm_rc_c_list_2_begin
          PROCEDURE :: next       => ppm_rc_c_list_2_next
          PROCEDURE :: last       => ppm_rc_c_list_2_last
          PROCEDURE :: prev       => ppm_rc_c_list_2_prev
          PROCEDURE :: destroy    => ppm_rc_c_list_2_destroy
          PROCEDURE :: exists     => ppm_rc_c_list_2_exists
          PROCEDURE :: has        => ppm_rc_c_list_2_has
          PROCEDURE :: get_id     => ppm_rc_c_list_2_get_id
          PROCEDURE :: push       => ppm_rc_c_list_2_push
          PROCEDURE :: remove     => ppm_rc_c_list_2_remove
          PROCEDURE :: grow_size  => ppm_rc_c_list_2_grow_size
          PROCEDURE :: at         => ppm_rc_c_list_2_at
        END TYPE ppm_rc_c_list_2



