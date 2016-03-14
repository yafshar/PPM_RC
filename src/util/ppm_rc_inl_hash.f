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


      !-------------------------------------------------------------------------
      !  Subroutine   :                  ppm_neighlist_hash
      !-------------------------------------------------------------------------
      ! Copyright (c) 2011 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich),
      !                    Center for Fluid Dynamics (DTU)
      !
      !
      ! This file is part of the Parallel Particle Mesh Library (PPM).
      !
      ! PPM is free software: you can redistribute it and/or modIFy
      ! it under the terms of the GNU Lesser General Public License
      ! as published by the Free Software Foundation, either
      ! version 3 of the License, or (at your option) any later
      ! version.
      !
      ! PPM is distributed in the hope that it will be useful,
      ! but WITHOUT ANY WARRANTY; without even the implied warranty of
      ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
      ! GNU General Public License for more details.
      !
      ! You should have received a copy of the GNU General Public License
      ! and the GNU Lesser General Public License along with PPM. If not,
      ! see <http://www.gnu.org/licenses/>.
      !
      ! Parallel Particle Mesh Library (PPM)
      ! ETH Zurich
      ! CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
      SUBROUTINE create_htable(table,nelement,info)
      !!! Given number of rows of the table, creates the hash table. Number of
      !!! rows will be greater than nelement, as we use the first value that is
      !!! power of 2 and greater than nelement.
      !!!
      !!! [WARNING]
      !!! If you allocate a hashtable with more than 2^31-1 elements the hash
      !!! function will most probably produce incorrect hash keys and fail.
      USE ppm_module_data
      USE ppm_module_alloc
      USE ppm_module_error
      USE ppm_module_substart
      USE ppm_module_substop
      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      CLASS(ppm_rc_htable)   :: table
      !!! The hashtable to create.
      INTEGER, INTENT(IN   ) :: nelement
      !!! Number of desired elements.
      INTEGER, INTENT(  OUT) :: info

      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double) :: t0

      INTEGER :: i

      CHARACTER(LEN=ppm_char) :: caller='create_htable'

      CALL substart(caller,t0,info)

      IF (nelement.LE.0) THEN
         table%nrow = 1
      ELSE
         table%nrow = 2**(CEILING(LOG(REAL(nelement))/LOG(2.0)))
      ENDIF

      !---------------------------------------------------------------------
      !  Allocate array for hash table keys and array for positions on "borders" array.
      !---------------------------------------------------------------------
      ALLOCATE(table%keys(table%nrow),table%borders_pos(table%nrow),STAT=Info)
      or_fail_alloc('Failed to alocate htable_keys & htable_borders_pos!',ppm_error=ppm_error_fatal)
      !---------------------------------------------------------------------
      !  Set everything to NULL.
      !---------------------------------------------------------------------
      FORALL (i=1:table%nrow)
         table%keys(i)        = htable_null_li
         table%borders_pos(i) = -one
      END FORALL

      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      END SUBROUTINE create_htable

      SUBROUTINE destroy_htable(table,info)

      USE ppm_module_data
      USE ppm_module_alloc
      USE ppm_module_error
      USE ppm_module_substart
      USE ppm_module_substop
      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      CLASS(ppm_rc_htable)   :: table
      !!! The hashtable to create. The pointer must not be NULL
      INTEGER, INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double)   :: t0

      CHARACTER(LEN=ppm_char) :: caller='destroy_htable'

      CALL substart(caller,t0,info)

      !---------------------------------------------------------------------
      !  Deallocate array for hash table keys and array for positions on "borders" array.
      !---------------------------------------------------------------------
      DEALLOCATE(table%keys,table%borders_pos,STAT=info)
      or_fail_dealloc('Failed to deallocate htable_keys & htable_borders_pos!',ppm_error=ppm_error_fatal)

      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      END SUBROUTINE destroy_htable

      ELEMENTAL FUNCTION h_func(table, key, seed) RESULT(hash_val)
        IMPLICIT NONE
        !---------------------------------------------------------------------
        !  Arguments
        !---------------------------------------------------------------------
        CLASS(ppm_rc_htable),    INTENT(IN   ) :: table
        !!! The hashtable to create. The pointer must not be NULL

        INTEGER(ppm_kind_int64), INTENT(IN   ) :: key
        !!! Input key
        INTEGER(ppm_kind_int64), INTENT(IN   ) :: seed
        !!! Seed to be used for mixing
        INTEGER(ppm_kind_int64)               :: hash_val
        !!! Result of the hash function
        !!!
        !!! [NOTE]
        !!! we need to work with 64 bit integers because
        !!! jump*h_func might overflow, and as there is no unsigned
        !!! integer the resulting value might be negative and not
        !!! a valid array index

        !---------------------------------------------------------------------
        !  Local variables and parameters
        !---------------------------------------------------------------------
        INTEGER(ppm_kind_int64), PARAMETER  :: m = 1431374979_ppm_kind_int64
        INTEGER                             :: h
        INTEGER(ppm_kind_int64)             :: data
        INTEGER                             :: k
        INTEGER                             :: len

        len = 4

        h = IEOR(seed, len)
        data = key

        DO WHILE (len .GE. 4)
            k = IBITS(data, 0, 8) !data, pos, len. len = 1 always!
            k = IOR(k, ISHFT(IBITS(data,  8, 8),  8))
            k = IOR(k, ISHFT(IBITS(data, 16, 8), 16))
            k = IOR(k, ISHFT(IBITS(data, 24, 8), 24))

            k = k*m
            k = IEOR(k, ISHFT(k, -24))
            k = k*m

            h = h*m
            h = IEOR(h, k)

            data = data + 4
            len  = len - 1
        ENDDO

        SELECT CASE (len)
        CASE (3)
           h = IEOR(h, ISHFT(IBITS(data, 16, 8), 16))
           h = IEOR(h, ISHFT(IBITS(data,  8, 8),  8))
           h = IEOR(h, IBITS(data, 0, 8))
           h = h*m

        CASE (2)
           h = IEOR(h, ISHFT(IBITS(data,  8, 8),  8))
           h = IEOR(h, IBITS(data, 0, 8))
           h = h*m

        CASE (1)
           h = IEOR(h, IBITS(data, 0, 8))
           h = h*m

        END SELECT

        h = IEOR(h, ISHFT(h, -13))
        h = h*m
        h = IEOR(h, ISHFT(h, -15))
        hash_val = IAND(h, table%nrow - 1)
        RETURN
      END FUNCTION

      ELEMENTAL FUNCTION h_key(table, key, jump) RESULT(address)
      !!! Given the key and jump value, returns corresponding address on
      !!! "borders" array.

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_htable),    INTENT(IN   ) :: table
      !!! The hashtable to create. The pointer must not be NULL
      INTEGER(ppm_kind_int64), INTENT(IN   ) :: key
      !!! Input key, which corresponds to address requested
      INTEGER,                 INTENT(IN   ) :: jump
      !!! Jump value for double hashing
      INTEGER                                :: address
      !!! Address that corresponds to given key
      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: int_addr
      ! we need to work with 64 bit integers because
      ! jump*h_func might overflow, and as there is no unsigned
      ! integer the resulting value might be negative and not
      ! a valid array index

      int_addr = 1_ppm_kind_int64 + &
      & MOD((table%h_func(key,seed1)+jump*table%h_func(key,seed2)),table%nrow)
      address  = INT(int_addr)
      RETURN
      END FUNCTION h_key

      SUBROUTINE hash_insert(table, key, value, info)
      !!! Given the key and the value, stores both in the hash table. Info is
      !!! set to -1 if size of the hash table is not sufficient.
      !!!
      !!! [NOTE]
      !!! This routine needs to be very fast, therefor we skip the usual
      !!! chit-chat and get right to it. (-> no substart,substop unless
      !!! compiled with __DEBUG flag)
      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_htable)                   :: table
      !!! The hashtable
      INTEGER(ppm_kind_int64), INTENT(IN   ) :: key
      !!! Key to be stored
      REAL(MK),                INTENT(IN   ) :: value
      !!! Value that corresponds to given key
      INTEGER,                 INTENT(  OUT) :: info
      !!! Info for whether insertion was successful or not. 0 if SUCCESSFUL
      !!! and -1 otherwise.

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER :: jump
      INTEGER :: spot

      info = 0
      jump = 0
      ! Keep on searching withing bounds of hash table.
      DO WHILE (jump .LT. table%nrow)
         ! Get the address corresponding to given key
         spot = table%h_key(key, jump)
         ! If an empty slot found ...
         IF (table%keys(spot).EQ.htable_null_li) THEN
            ! Store the key and the corresponding value and RETURN.
            table%keys(spot) = key
            table%borders_pos(spot) = value
            RETURN
         !If the key is the same the value should be updated
         ELSE IF (table%keys(spot).EQ.key) THEN
            table%borders_pos(spot) = value
            RETURN
         ENDIF
         ! If the current slot is occupied, jump to next key that results
         ! in same hash function.
         jump = jump + 1
      ENDDO

      ! If NOT returned within the while-loop, that means our hash table is
      ! not sufficiently large. Hence, regrowth will be needed!
      CALL table%grow(info,key,value)
      RETURN
      END SUBROUTINE hash_insert

      SUBROUTINE hash_insert_(table,key_,value,info)

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_htable)    :: table
      !!! The hashtable
      INTEGER,  INTENT(IN   ) :: key_
      !!! Key to be stored
      REAL(MK), INTENT(IN   ) :: value
      !!! Value that corresponds to given key
      INTEGER,  INTENT(  OUT) :: info
      !!! Info for whether insertion was successful or not. 0 if SUCCESSFUL
      !!! and -1 otherwise.

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: key

      key=INT(key_,KIND=ppm_kind_int64)
      CALL table%hash_insert(key,value,info)
      RETURN
      END SUBROUTINE hash_insert_

      SUBROUTINE hash_insert__(table,key_,value,info)

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_htable)                 :: table
      !!! The hashtable
      INTEGER, DIMENSION(:), INTENT(IN   ) :: key_
      !!! Key to be stored
      REAL(MK),              INTENT(IN   ) :: value
      !!! Value that corresponds to given key
      INTEGER,               INTENT(  OUT) :: info
      !!! Info for whether insertion was successful or not. 0 if SUCCESSFUL
      !!! and -1 otherwise.

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: key
      INTEGER                 :: ssize

      ssize=SIZE(key_)
      SELECT CASE (ssize)
      CASE (2)
         key=IndexHashFunctor64_2d(key_)
      CASE (3)
         key=IndexHashFunctor64_3d(key_)
      CASE DEFAULT
         info=ppm_error_fatal
         RETURN
      END SELECT
      CALL table%hash_insert(key,value,info)
      RETURN
      END SUBROUTINE hash_insert__

      SUBROUTINE hash_insert_2d(table,key_1,key_2,value,info)

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_htable)    :: table
      !!! The hashtable
      INTEGER,  INTENT(IN   ) :: key_1
      INTEGER,  INTENT(IN   ) :: key_2
      !!! Key to be stored
      REAL(MK), INTENT(IN   ) :: value
      !!! Value that corresponds to given key
      INTEGER,  INTENT(  OUT) :: info
      !!! Info for whether insertion was successful or not. 0 if SUCCESSFUL
      !!! and -1 otherwise.

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: key
      key=IndexHashFunctor64_2d(key_1,key_2)
      CALL table%hash_insert(key,value,info)
      RETURN
      END SUBROUTINE hash_insert_2d

      SUBROUTINE hash_insert_3d(table,key_1,key_2,key_3,value,info)

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_htable)    :: table
      !!! The hashtable
      INTEGER,  INTENT(IN   ) :: key_1
      INTEGER,  INTENT(IN   ) :: key_2
      INTEGER,  INTENT(IN   ) :: key_3
      !!! Key to be stored
      REAL(MK), INTENT(IN   ) :: value
      !!! Value that corresponds to given key
      INTEGER,  INTENT(  OUT) :: info
      !!! Info for whether insertion was successful or not. 0 if SUCCESSFUL
      !!! and -1 otherwise.

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: key

      key=IndexHashFunctor64_3d(key_1,key_2,key_3)
      CALL table%hash_insert(key,value,info)
      RETURN
      END SUBROUTINE hash_insert_3d

      FUNCTION hash_search(table,key) RESULT(value)
      !!! Given the key, searchs the key on the hash table and returns the
      !!! corresponding value.
      !!!
      !!! [NOTE]
      !!! This routine needs to be very fast, therefor we skip the usual
      !!! chit-chat and get right to it. (-> no substart,substop unless
      !!! compiled with __DEBUG flag)
      IMPLICIT NONE
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_htable),    INTENT(IN   ) :: table
      !!! The hashtable to create. The pointer must not be NULL
      INTEGER(ppm_kind_int64), INTENT(IN   ) :: key
      !!! Input key, which the corresponding value is asked for
      REAL(MK)                               :: value
      !!! Value corresponding to the input key

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER :: jump
      INTEGER :: spot

      jump = 0
      ! Keep on searching while we don't come across a NULL value or we don't
      ! exceed bounds of hash table.
      loop: DO WHILE(jump .LT. table%nrow)
          ! Get the other key that results in same hash key as for the input
          ! key.
          spot = table%h_key(key, jump)
          ! If key matches ...
          IF (table%keys(spot).EQ.key) THEN
             ! Set the return value and return
             value = table%borders_pos(spot)
             RETURN
          ELSE IF (table%keys(spot).EQ.htable_null_li) THEN
             IF (.NOT.ANY(key.EQ.table%keys)) EXIT loop
          ENDIF
          ! Otherwise, keep on incrementing jump distance
          jump = jump + 1
      ENDDO loop
      value = -one
      RETURN
      END FUNCTION hash_search

      FUNCTION hash_search_(table,key_) RESULT(value)
      IMPLICIT NONE
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_htable), INTENT(IN   ) :: table
      !!! The hashtable
      INTEGER,              INTENT(IN   ) :: key_
      !!! Input key, which the corresponding value is asked for
      REAL(MK)                            :: value
      !!! Value corresponding to the input key

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: key

      key=INT(key_,KIND=ppm_kind_int64)
      value=table%hash_search(key)
      RETURN
      END FUNCTION hash_search_

      FUNCTION hash_search_2d(table,key_1,key_2) RESULT(value)

      IMPLICIT NONE
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_htable), INTENT(IN   ) :: table
      !!! The hashtable
      INTEGER,              INTENT(IN   ) :: key_1
      INTEGER,              INTENT(IN   ) :: key_2
      !!! Input key, which the corresponding value is asked for
      REAL(MK)                            :: value
      !!! Value corresponding to the input key

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: key

      key=IndexHashFunctor64_2d(key_1,key_2)
      value=table%hash_search(key)
      RETURN
      END FUNCTION hash_search_2d

      FUNCTION hash_search_3d(table,key_1,key_2,key_3) RESULT(value)
      IMPLICIT NONE
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_htable), INTENT(IN   ) :: table
      !!! The hashtable
      INTEGER,              INTENT(IN   ) :: key_1
      INTEGER,              INTENT(IN   ) :: key_2
      INTEGER,              INTENT(IN   ) :: key_3
      !!! Input key, which the corresponding value is asked for
      REAL(MK)                            :: value
      !!! Value corresponding to the input key

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: key

      key=IndexHashFunctor64_3d(key_1,key_2,key_3)
      value=table%hash_search(key)
      RETURN
      END FUNCTION hash_search_3d

      !TODO check this sub
      !Yaser
      !This function PROBABLY suffers from a bug!
      SUBROUTINE hash_remove(table, key, info,existed)
      !!! Given the key, removes the elements in the hash table. Info is
      !!! set to -1 if the key was NOT found.
      !!!
      !!! [NOTE]
      !!! This routine needs to be very fast, therefor we skip the usual
      !!! chit-chat and get right to it. (-> no substart,substop unless
      !!! compiled with __DEBUG flag)
      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_htable)                   :: table
      !!! The hashtable

      INTEGER(ppm_kind_int64), INTENT(IN   ) :: key
      !!! Key to be removed
      INTEGER,                 INTENT(  OUT) :: info
      !!! Info for whether removal was successful or not. 0 if SUCCESSFUL
      !!! and -1 otherwise.

      LOGICAL, OPTIONAL,       INTENT(IN   ) :: existed
      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER :: jump
      INTEGER :: spot !,spot0

      info = 0
      IF (PRESENT(existed)) THEN
         IF (existed) THEN
            jump = 0
            ! Keep on searching withing bounds of hash table.
            DO WHILE (jump.LT.table%nrow)
               ! Get the address corresponding to given key
               spot = table%h_key(key, jump)
               ! If an empty slot found ...
               IF (table%keys(spot).EQ.key) THEN
                  !Remove the key and the corresponding value and RETURN.
                  table%borders_pos(spot)=-one
                  table%keys(spot)=htable_null_li
                  RETURN
               ENDIF
               ! If the current slot is occupied, jump to next key that results
               ! in same hash function.
               jump = jump + 1
            ENDDO
            ! If NOT returned within the while-loop, that means the key was NOT found
            info = ppm_error_fatal
         ENDIf
         RETURN
      ELSE
         IF (ANY(key.EQ.table%keys)) THEN
            jump = 0
            ! Keep on searching withing bounds of hash table.
            DO WHILE (jump.LT.table%nrow)
               ! Get the address corresponding to given key
               spot = table%h_key(key, jump)
               ! If an empty slot found ...
               IF (table%keys(spot).EQ.key) THEN
                  !Remove the key and the corresponding value and RETURN.
                  table%borders_pos(spot)=-one
                  table%keys(spot)=htable_null_li
                  RETURN
               ENDIF
               ! If the current slot is occupied, jump to next key that results
               ! in same hash function.
               jump = jump + 1
            ENDDO
            ! If NOT returned within the while-loop, that means the key was NOT found
            info = ppm_error_fatal
         ENDIF
      ENDIF
      RETURN
      END SUBROUTINE hash_remove

      SUBROUTINE hash_remove_(table, key_, info,existed)

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_htable)             :: table
      !!! The hashtable
      INTEGER,           INTENT(IN   ) :: key_
      !!! Key to be removed
      INTEGER,           INTENT(  OUT) :: info
      !!! Info for whether removal was successful or not. 0 if SUCCESSFUL
      !!! and -1 otherwise.
      LOGICAL, OPTIONAL, INTENT(IN   ) :: existed

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: key

      key=INT(key_,KIND=ppm_kind_int64)
      IF (PRESENT(existed)) THEN
         CALL table%hash_remove(key,info,existed)
      ELSE
         CALL table%hash_remove(key,info)
      ENDIF
      RETURN
      END SUBROUTINE hash_remove_

      SUBROUTINE hash_remove__(table, key_, info,existed)

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_htable)                 :: table
      !!! The hashtable
      INTEGER, DIMENSION(:), INTENT(IN   ) :: key_
      !!! Key to be removed
      INTEGER,               INTENT(  OUT) :: info
      !!! Info for whether removal was successful or not. 0 if SUCCESSFUL
      !!! and -1 otherwise.
      LOGICAL, OPTIONAL,     INTENT(IN   ) :: existed

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: key
      INTEGER                 :: ssize

      ssize=SIZE(key_)

      SELECT CASE (ssize)
      CASE (2)
         key=IndexHashFunctor64_2d(key_)
      CASE (3)
         key=IndexHashFunctor64_3d(key_)
      CASE DEFAULT
         info=ppm_error_fatal
         RETURN
      END SELECT
      IF (PRESENT(existed)) THEN
         CALL table%hash_remove(key,info,existed)
      ELSE
         CALL table%hash_remove(key,info)
      ENDIF
      RETURN
      END SUBROUTINE hash_remove__

      SUBROUTINE hash_remove_2d(table,key_1,key_2,info,existed)

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_htable)             :: table
      !!! The hashtable
      INTEGER,           INTENT(IN   ) :: key_1
      INTEGER,           INTENT(IN   ) :: key_2
      !!! Key to be removed
      INTEGER,           INTENT(  OUT) :: info
      !!! Info for whether removal was successful or not. 0 if SUCCESSFUL
      !!! and -1 otherwise.
      LOGICAL, OPTIONAL, INTENT(IN   ) :: existed

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: key

      key=IndexHashFunctor64_2d(key_1,key_2)
      IF (PRESENT(existed)) THEN
         CALL table%hash_remove(key,info,existed)
      ELSE
         CALL table%hash_remove(key,info)
      ENDIF
      RETURN
      END SUBROUTINE hash_remove_2d

      SUBROUTINE hash_remove_3d(table,key_1,key_2,key_3,info,existed)

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_htable)             :: table
      !!! The hashtable
      INTEGER,           INTENT(IN   ) :: key_1
      INTEGER,           INTENT(IN   ) :: key_2
      INTEGER,           INTENT(IN   ) :: key_3
      !!! Key to be removed
      INTEGER,           INTENT(  OUT) :: info
      !!! Info for whether removal was successful or not. 0 if SUCCESSFUL
      !!! and -1 otherwise.
      LOGICAL, OPTIONAL, INTENT(IN   ) :: existed

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: key

      key=IndexHashFunctor64_3d(key_1,key_2,key_3)
      IF (PRESENT(existed)) THEN
         CALL table%hash_remove(key,info,existed)
      ELSE
         CALL table%hash_remove(key,info)
      ENDIF
      RETURN
      END SUBROUTINE hash_remove_3d

      SUBROUTINE grow_htable(table,info,key_,value_)
      !!! Based on the number of rows of the table, creates the hash table with
      !!! double size.
      !!!
      !!! [WARNING]
      !!! If you allocate a hashtable with more than 2^31-1 elements the hash
      !!! function will most probably produce incorrect hash keys and fail.
      USE ppm_module_data
      USE ppm_module_alloc
      USE ppm_module_error
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_util_qsort
      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      CLASS(ppm_rc_htable)                             :: table
      !!! The hashtable to grow.

      INTEGER,                           INTENT(  OUT) :: info
      INTEGER(ppm_kind_int64), OPTIONAL, INTENT(IN   ) :: key_
      !!! Key to be stored

      REAL(MK),                OPTIONAL, INTENT(IN   ) :: value_
      !!! Value that corresponds to given key
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double)               :: t0
      REAL(MK), DIMENSION(:), ALLOCATABLE :: borders_pos_tmp

      INTEGER(ppm_kind_int64), DIMENSION(:), ALLOCATABLE :: keys_tmp
      INTEGER,                 DIMENSION(:), POINTER     :: ranklist
      INTEGER                                            :: nsize
      INTEGER                                            :: i,j

      CHARACTER(LEN=ppm_char) :: caller='grow_htable'

      CALL substart(caller,t0,info)

      nsize=table%nrow

      ALLOCATE(keys_tmp(nsize),SOURCE=table%keys,STAT=info)
      or_fail_alloc("keys_tmp")

      ALLOCATE(borders_pos_tmp(nsize),SOURCE=table%borders_pos,STAT=info)
      or_fail_alloc("keys_tmp & borders_pos_tmp")

      NULLIFY(ranklist)
      CALL ppm_util_qsort(keys_tmp,ranklist,info,nsize)
      or_fail("ppm_util_qsort")

      nsize=table%nrow*2

      IF (nsize.GE.ppm_big_i-1) THEN
         !TOCHCECK
         fail("hashtable with more than 2^31-1 elements will fail",ppm_error=ppm_error_fatal)
      ENDIF

      CALL table%destroy(info)
      or_fail("table%destroy")

      CALL table%create(nsize,info)
      or_fail("table%create")

      DO i=1,SIZE(keys_tmp)
         j=ranklist(i)
         IF (keys_tmp(j).EQ.htable_null_li) CYCLE
         CALL table%insert(keys_tmp(j),borders_pos_tmp(j),info)
      ENDDO

      DEALLOCATE(keys_tmp,borders_pos_tmp,STAT=info)
      or_fail_dealloc("keys_tmp & borders_pos_tmp")

      !---------------------------------------------------------------------
      !  Deallocate ranklist array.
      !---------------------------------------------------------------------
      CALL ppm_alloc(ranklist,(/0/),ppm_param_dealloc,info)
      or_fail_dealloc('htable_keys',ppm_error=ppm_error_fatal)

      IF (PRESENT(key_)) THEN
         IF (PRESENT(value_)) THEN
            CALL table%insert(key_,value_,info)
         ENDIF
      ENDIF

      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      END SUBROUTINE grow_htable

      SUBROUTINE shrink_htable(table,info,shrinkage_ratio)
      !!! Based on the number of rows of the table, shrinks the hash table with
      !!! half a size or more.
      USE ppm_module_data
      USE ppm_module_alloc
      USE ppm_module_error
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_util_qsort
      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      CLASS(ppm_rc_htable)             :: table
      !!! The hashtable to shrink.

      INTEGER,           INTENT(  OUT) :: info
      INTEGER, OPTIONAL, INTENT(IN   ) :: shrinkage_ratio
      !!! OPTIONAL shrinkage_ratio (positive value).
      !!! If the size of hash table is shrinkage_ratio times bigger than the
      !!! real elements inside table, we reduce the table size
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double)               :: t0
      REAL(MK), DIMENSION(:), ALLOCATABLE :: borders_pos_tmp

      INTEGER(ppm_kind_int64), DIMENSION(:), ALLOCATABLE :: keys_tmp
      INTEGER                                            :: shrinkage_ratio_
      INTEGER,                 DIMENSION(:), POINTER     :: ranklist
      INTEGER                                            :: nsize,ssize
      INTEGER                                            :: i,j

      CHARACTER(LEN=ppm_char) :: caller='shrink_htable'

      CALL substart(caller,t0,info)

      shrinkage_ratio_=MERGE(shrinkage_ratio,4,PRESENT(shrinkage_ratio))
      shrinkage_ratio_=MERGE(4,shrinkage_ratio_,shrinkage_ratio_.LE.0)

      nsize=table%nrow
      ssize=table%size()

      IF (nsize.GE.shrinkage_ratio_*ssize) THEN
         ALLOCATE(keys_tmp(nsize),SOURCE=table%keys,STAT=info)
         or_fail_alloc("keys_tmp")

         ALLOCATE(borders_pos_tmp(nsize),SOURCE=table%borders_pos,STAT=info)
         or_fail_alloc("keys_tmp & borders_pos_tmp")

         NULLIFY(ranklist)
         CALL ppm_util_qsort(keys_tmp,ranklist,info,nsize)
         or_fail("ppm_util_qsort")

         CALL table%destroy(info)
         or_fail("table%destroy")

         CALL table%create(ssize,info)
         or_fail("table%create")

         DO i=1,nsize
            j=ranklist(i)
            IF (keys_tmp(j).EQ.htable_null_li) CYCLE
            CALL table%insert(keys_tmp(j),borders_pos_tmp(j),info)
         ENDDO

         DEALLOCATE(keys_tmp,borders_pos_tmp,STAT=info)
         or_fail_dealloc("keys_tmp & borders_pos_tmp")

         !---------------------------------------------------------------------
         !  Deallocate ranklist array.
         !---------------------------------------------------------------------
         CALL ppm_alloc(ranklist,(/0/),ppm_param_dealloc,info)
         or_fail_dealloc('htable_keys',ppm_error=ppm_error_fatal)
      ENDIF

      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      END SUBROUTINE shrink_htable

      FUNCTION hash_size(table) RESULT(value)
      IMPLICIT NONE
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_htable) :: table
      !!! The hashtable
      INTEGER              :: value
      !!!
      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      IF (table%nrow.GT.0) THEN
         value=COUNT(table%keys.NE.htable_null_li)
      ELSE
         value=0
      ENDIF
      RETURN
      END FUNCTION hash_size

      SUBROUTINE MCMCParticle_create_htable(table,nelement,info)
      !!! Given number of rows of the table, creates the hash table. Number of
      !!! rows will be greater than nelement, as we use the first value that is
      !!! power of 2 and greater than nelement.
      !!!
      !!! [WARNING]
      !!! If you allocate a hashtable with more than 2^31-1 elements the hash
      !!! function will most probably produce incorrect hash keys and fail.
      USE ppm_module_data
      USE ppm_module_alloc
      USE ppm_module_error
      USE ppm_module_substart
      USE ppm_module_substop
      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      CLASS(ppm_rc_MCMCParticlehtable) :: table
      !!! The hashtable to create.
      INTEGER,           INTENT(IN   ) :: nelement
      !!! Number of desired elements.
      INTEGER,           INTENT(  OUT) :: info

      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double) :: t0

      INTEGER :: i

      CHARACTER(LEN=ppm_char) :: caller='MCMCParticle_create_htable'

      CALL substart(caller,t0,info)

      IF (nelement.LE.0) THEN
         table%nrow = 1
      ELSE
         table%nrow = 2**(CEILING(LOG(REAL(nelement))/LOG(2.0)))
      ENDIF

      !---------------------------------------------------------------------
      !  Allocate array for hash table keys and array for positions on "borders" array.
      !---------------------------------------------------------------------
      ALLOCATE(table%keys(table%nrow),table%borders_pos(table%nrow),STAT=Info)
      or_fail_alloc('Failed to alocate htable_keys & htable_borders_pos!',ppm_error=ppm_error_fatal)
      !---------------------------------------------------------------------
      !  Set everything to NULL.
      !---------------------------------------------------------------------
      FORALL (i=1:table%nrow)
         table%keys(i)        = htable_null_li
         table%borders_pos(i) = MCMCParticle(-1,zero)
      END FORALL

      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      END SUBROUTINE MCMCParticle_create_htable

      SUBROUTINE MCMCParticle_destroy_htable(table,info)

      USE ppm_module_data
      USE ppm_module_alloc
      USE ppm_module_error
      USE ppm_module_substart
      USE ppm_module_substop
      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      CLASS(ppm_rc_MCMCParticlehtable) :: table
      !!! The hashtable to create. The pointer must not be NULL
      INTEGER,           INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double) :: t0

      CHARACTER(LEN=ppm_char) :: caller='MCMCParticle_destroy_htable'

      CALL substart(caller,t0,info)

      !---------------------------------------------------------------------
      !  Deallocate array for hash table keys and array for positions on "borders" array.
      !---------------------------------------------------------------------
      DEALLOCATE(table%keys,table%borders_pos,STAT=info)
      or_fail_dealloc('Failed to deallocate htable_keys & htable_borders_pos!',ppm_error=ppm_error_fatal)

      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      END SUBROUTINE MCMCParticle_destroy_htable

      ELEMENTAL FUNCTION MCMCParticle_h_func(table,key,seed) RESULT(hash_val)
        IMPLICIT NONE
        !---------------------------------------------------------------------
        !  Arguments
        !---------------------------------------------------------------------
        CLASS(ppm_rc_MCMCParticlehtable), INTENT(IN   ) :: table
        !!! The hashtable to create. The pointer must not be NULL
        INTEGER(ppm_kind_int64),          INTENT(IN   ) :: key
        !!! Input key
        INTEGER(ppm_kind_int64),          INTENT(IN   ) :: seed
        !!! Seed to be used for mixing
        INTEGER(ppm_kind_int64)                         :: hash_val
        !!! Result of the hash function
        !!!
        !!! [NOTE]
        !!! we need to work with 64 bit integers because
        !!! jump*h_func might overflow, and as there is no unsigned
        !!! integer the resulting value might be negative and not
        !!! a valid array index

        !---------------------------------------------------------------------
        !  Local variables and parameters
        !---------------------------------------------------------------------
        INTEGER(ppm_kind_int64), PARAMETER  :: m = 1431374979_ppm_kind_int64
        INTEGER                             :: h
        INTEGER(ppm_kind_int64)             :: data
        INTEGER                             :: k
        INTEGER                             :: len

        len = 4

        h = IEOR(seed, len)
        data = key

        DO WHILE (len .GE. 4)
            k = IBITS(data, 0, 8) !data, pos, len. len = 1 always!
            k = IOR(k, ISHFT(IBITS(data,  8, 8),  8))
            k = IOR(k, ISHFT(IBITS(data, 16, 8), 16))
            k = IOR(k, ISHFT(IBITS(data, 24, 8), 24))

            k = k*m
            k = IEOR(k, ISHFT(k, -24))
            k = k*m

            h = h*m
            h = IEOR(h, k)

            data = data + 4
            len  = len - 1
        ENDDO

        SELECT CASE (len)
        CASE (3)
           h = IEOR(h, ISHFT(IBITS(data, 16, 8), 16))
           h = IEOR(h, ISHFT(IBITS(data,  8, 8),  8))
           h = IEOR(h, IBITS(data, 0, 8))
           h = h*m

        CASE (2)
           h = IEOR(h, ISHFT(IBITS(data,  8, 8),  8))
           h = IEOR(h, IBITS(data, 0, 8))
           h = h*m

        CASE (1)
           h = IEOR(h, IBITS(data, 0, 8))
           h = h*m

        END SELECT

        h = IEOR(h, ISHFT(h, -13))
        h = h*m
        h = IEOR(h, ISHFT(h, -15))
        hash_val = IAND(h, table%nrow - 1)
        RETURN
      END FUNCTION

      ELEMENTAL FUNCTION MCMCParticle_h_key(table,key,jump) RESULT(address)
      !!! Given the key and jump value, returns corresponding address on
      !!! "borders" array.

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_MCMCParticlehtable), INTENT(IN   ) :: table
      !!! The hashtable to create. The pointer must not be NULL
      INTEGER(ppm_kind_int64),          INTENT(IN   ) :: key
      !!! Input key, which corresponds to address requested
      INTEGER,                          INTENT(IN   ) :: jump
      !!! Jump value for double hashing
      INTEGER                                         :: address
      !!! Address that corresponds to given key
      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: int_addr
      ! we need to work with 64 bit integers because
      ! jump*h_func might overflow, and as there is no unsigned
      ! integer the resulting value might be negative and not
      ! a valid array index

      int_addr = 1_ppm_kind_int64 + &
      & MOD((table%h_func(key,seed1)+jump*table%h_func(key,seed2)),table%nrow)
      address  = INT(int_addr)
      RETURN
      END FUNCTION MCMCParticle_h_key

      SUBROUTINE MCMCParticle_hash_insert(table,key,value,info)
      !!! Given the key and the value, stores both in the hash table. Info is
      !!! set to -1 if size of the hash table is not sufficient.
      !!!
      !!! [NOTE]
      !!! This routine needs to be very fast, therefor we skip the usual
      !!! chit-chat and get right to it. (-> no substart,substop unless
      !!! compiled with __DEBUG flag)
      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_MCMCParticlehtable)       :: table
      !!! The hashtable
      INTEGER(ppm_kind_int64), INTENT(IN   ) :: key
      !!! Key to be stored
      TYPE(MCMCParticle),      INTENT(IN   ) :: value
      !!! Value that corresponds to given key
      INTEGER,                 INTENT(  OUT) :: info
      !!! Info for whether insertion was successful or not. 0 if SUCCESSFUL
      !!! and -1 otherwise.

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER :: jump
      INTEGER :: spot

      info = 0
      jump = 0
      ! Keep on searching withing bounds of hash table.
      DO WHILE (jump .LT. table%nrow)
         ! Get the address corresponding to given key
         spot = table%h_key(key, jump)
         ! If an empty slot found ...
         IF (table%keys(spot).EQ.htable_null_li) THEN
            ! Store the key and the corresponding value and RETURN.
            table%keys(spot) = key
            table%borders_pos(spot) = value
            RETURN
         !If the key is the same the value should be updated
         ELSE IF (table%keys(spot).EQ.key) THEN
            table%borders_pos(spot) = value
            RETURN
         ENDIF
         ! If the current slot is occupied, jump to next key that results
         ! in same hash function.
         jump = jump + 1
      ENDDO

      ! If NOT returned within the while-loop, that means our hash table is
      ! not sufficiently large. Hence, regrowth will be needed!
      CALL table%grow(info,key,value)
      RETURN
      END SUBROUTINE MCMCParticle_hash_insert

      SUBROUTINE MCMCParticle_hash_insert_(table,key_,value,info)

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_MCMCParticlehtable)  :: table
      !!! The hashtable
      INTEGER,            INTENT(IN   ) :: key_
      !!! Key to be stored
      TYPE(MCMCParticle), INTENT(IN   ) :: value
      !!! Value that corresponds to given key
      INTEGER,            INTENT(  OUT) :: info
      !!! Info for whether insertion was successful or not. 0 if SUCCESSFUL
      !!! and -1 otherwise.

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: key

      key=INT(key_,KIND=ppm_kind_int64)
      CALL table%MCMCParticle_hash_insert(key,value,info)
      RETURN
      END SUBROUTINE MCMCParticle_hash_insert_

      SUBROUTINE MCMCParticle_hash_insert__(table,key_,value,info)

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_MCMCParticlehtable)     :: table
      !!! The hashtable
      INTEGER, DIMENSION(:), INTENT(IN   ) :: key_
      !!! Key to be stored
      TYPE(MCMCParticle),    INTENT(IN   ) :: value
      !!! Value that corresponds to given key
      INTEGER,               INTENT(  OUT) :: info
      !!! Info for whether insertion was successful or not. 0 if SUCCESSFUL
      !!! and -1 otherwise.

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: key
      INTEGER                 :: ssize

      ssize=SIZE(key_)
      SELECT CASE (ssize)
      CASE (2)
         key=IndexHashFunctor64_2d(key_)
      CASE (3)
         key=IndexHashFunctor64_3d(key_)
      CASE DEFAULT
         info=ppm_error_fatal
         RETURN
      END SELECT
      CALL table%MCMCParticle_hash_insert(key,value,info)
      RETURN
      END SUBROUTINE MCMCParticle_hash_insert__

      SUBROUTINE MCMCParticle_hash_insert_2d(table,key_1,key_2,value,info)

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_MCMCParticlehtable)  :: table
      !!! The hashtable
      INTEGER,            INTENT(IN   ) :: key_1
      INTEGER,            INTENT(IN   ) :: key_2
      !!! Key to be stored
      TYPE(MCMCParticle), INTENT(IN   ) :: value
      !!! Value that corresponds to given key
      INTEGER,            INTENT(  OUT) :: info
      !!! Info for whether insertion was successful or not. 0 if SUCCESSFUL
      !!! and -1 otherwise.

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: key
      key=IndexHashFunctor64_2d(key_1,key_2)
      CALL table%MCMCParticle_hash_insert(key,value,info)
      RETURN
      END SUBROUTINE MCMCParticle_hash_insert_2d

      SUBROUTINE MCMCParticle_hash_insert_3d(table,key_1,key_2,key_3,value,info)

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_MCMCParticlehtable)  :: table
      !!! The hashtable
      INTEGER,            INTENT(IN   ) :: key_1
      INTEGER,            INTENT(IN   ) :: key_2
      INTEGER,            INTENT(IN   ) :: key_3
      !!! Key to be stored
      TYPE(MCMCParticle), INTENT(IN   ) :: value
      !!! Value that corresponds to given key
      INTEGER,            INTENT(  OUT) :: info
      !!! Info for whether insertion was successful or not. 0 if SUCCESSFUL
      !!! and -1 otherwise.

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: key

      key=IndexHashFunctor64_3d(key_1,key_2,key_3)
      CALL table%MCMCParticle_hash_insert(key,value,info)
      RETURN
      END SUBROUTINE MCMCParticle_hash_insert_3d

      FUNCTION MCMCParticle_hash_search(table,key) RESULT(value)
      !!! Given the key, searchs the key on the hash table and returns the
      !!! corresponding value.
      !!!
      !!! [NOTE]
      !!! This routine needs to be very fast, therefor we skip the usual
      !!! chit-chat and get right to it. (-> no substart,substop unless
      !!! compiled with __DEBUG flag)
      IMPLICIT NONE
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_MCMCParticlehtable), INTENT(IN   ) :: table
      !!! The hashtable to create. The pointer must not be NULL
      INTEGER(ppm_kind_int64),          INTENT(IN   ) :: key
      !!! Input key, which the corresponding value is asked for
      TYPE(MCMCParticle)                              :: value
      !!! Value corresponding to the input key

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER :: jump
      INTEGER :: spot
      jump = 0
      ! Keep on searching while we don't come across a NULL value or we don't
      ! exceed bounds of hash table.
      loop: DO WHILE(jump .LT. table%nrow)
          ! Get the other key that results in same hash key as for the input
          ! key.
          spot = table%h_key(key,jump)
          ! If key matches ...
          IF (table%keys(spot).EQ.key) THEN
             ! Set the return value and return
             value = table%borders_pos(spot)
             RETURN
          ELSE IF (table%keys(spot).EQ.htable_null_li) THEN
             IF (.NOT.ANY(key.EQ.table%keys)) EXIT loop
          ENDIF
          ! Otherwise, keep on incrementing jump distance
          jump = jump + 1
      ENDDO loop
      value = MCMCParticle(-1,zero)
      RETURN
      END FUNCTION MCMCParticle_hash_search

      FUNCTION MCMCParticle_hash_search_(table,key_) RESULT(value)
      IMPLICIT NONE
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_MCMCParticlehtable), INTENT(IN   ) :: table
      !!! The hashtable
      INTEGER,                          INTENT(IN   ) :: key_
      !!! Input key, which the corresponding value is asked for
      TYPE(MCMCParticle)                              :: value
      !!! Value corresponding to the input key

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: key

      key=INT(key_,KIND=ppm_kind_int64)
      value=table%MCMCParticle_hash_search(key)
      RETURN
      END FUNCTION MCMCParticle_hash_search_

      FUNCTION MCMCParticle_hash_search_2d(table,key_1,key_2) RESULT(value)

      IMPLICIT NONE
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_MCMCParticlehtable), INTENT(IN   ) :: table
      !!! The hashtable
      INTEGER,                          INTENT(IN   ) :: key_1
      INTEGER,                          INTENT(IN   ) :: key_2
      !!! Input key, which the corresponding value is asked for
      TYPE(MCMCParticle)                              :: value
      !!! Value corresponding to the input key

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: key

      key=IndexHashFunctor64_2d(key_1,key_2)
      value=table%MCMCParticle_hash_search(key)
      RETURN
      END FUNCTION MCMCParticle_hash_search_2d

      FUNCTION MCMCParticle_hash_search_3d(table,key_1,key_2,key_3) RESULT(value)

      IMPLICIT NONE
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_MCMCParticlehtable), INTENT(IN   ) :: table
      !!! The hashtable
      INTEGER,                          INTENT(IN   ) :: key_1
      INTEGER,                          INTENT(IN   ) :: key_2
      INTEGER,                          INTENT(IN   ) :: key_3
      !!! Input key, which the corresponding value is asked for
      TYPE(MCMCParticle)                              :: value
      !!! Value corresponding to the input key

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: key

      key=IndexHashFunctor64_3d(key_1,key_2,key_3)
      value=table%MCMCParticle_hash_search(key)
      RETURN
      END FUNCTION MCMCParticle_hash_search_3d

      LOGICAL FUNCTION MCMCParticle_hash_contains(table,key,value)
      !!! Given the key, searchs the key on the hash table and returns the
      !!! corresponding value.
      !!!
      !!! [NOTE]
      !!! This routine needs to be very fast, therefor we skip the usual
      !!! chit-chat and get right to it. (-> no substart,substop unless
      !!! compiled with __DEBUG flag)
      IMPLICIT NONE
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_MCMCParticlehtable), INTENT(IN   ) :: table
      !!! The hashtable to create. The pointer must not be NULL
      INTEGER(ppm_kind_int64),          INTENT(IN   ) :: key
      !!! Input key, which the corresponding value is asked for
      TYPE(MCMCParticle),               INTENT(IN   ) :: value
      !!! Value corresponding to the input key

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER :: jump
      INTEGER :: spot

      jump = 0
      ! Keep on searching while we don't come across a NULL value or we don't
      ! exceed bounds of hash table.
      loop: DO WHILE(jump .LT. table%nrow)
          ! Get the other key that results in same hash key as for the input
          ! key.
          spot = table%h_key(key,jump)
          ! If key matches ...
          IF (table%keys(spot).EQ.key) THEN
             ! Set the return value and return
             ASSOCIATE (tmpParticle => table%borders_pos(spot))
                IF (value%candlabel.NE.tmpParticle%candlabel) THEN
                   MCMCParticle_hash_contains=.FALSE.
                   RETURN
                ENDIF
             END ASSOCIATE
             MCMCParticle_hash_contains=.TRUE.
             RETURN
          ELSE IF (table%keys(spot).EQ.htable_null_li) THEN
             IF (.NOT.ANY(key.EQ.table%keys)) EXIT loop
          ENDIF
          ! Otherwise, keep on incrementing jump distance
          jump = jump + 1
      ENDDO loop
      MCMCParticle_hash_contains=.FALSE.
      RETURN
      END FUNCTION MCMCParticle_hash_contains

      LOGICAL FUNCTION MCMCParticle_hash_contains_(table,key_,value)
      IMPLICIT NONE
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_MCMCParticlehtable), INTENT(IN   ) :: table
      !!! The hashtable
      INTEGER,                          INTENT(IN   ) :: key_
      !!! Input key, which the corresponding value is asked for
      TYPE(MCMCParticle),               INTENT(IN   ) :: value
      !!! Value corresponding to the input key

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: key

      key=INT(key_,KIND=ppm_kind_int64)
      MCMCParticle_hash_contains_=table%MCMCParticle_hash_contains(key,value)
      RETURN
      END FUNCTION MCMCParticle_hash_contains_

      LOGICAL FUNCTION MCMCParticle_hash_contains_2d(table,key_1,key_2,value)

      IMPLICIT NONE
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_MCMCParticlehtable), INTENT(IN   ) :: table
      !!! The hashtable
      INTEGER,                          INTENT(IN   ) :: key_1
      INTEGER,                          INTENT(IN   ) :: key_2
      !!! Input key, which the corresponding value is asked for
      TYPE(MCMCParticle),               INTENT(IN   ) :: value
      !!! Value corresponding to the input key

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: key

      key=IndexHashFunctor64_2d(key_1,key_2)
      MCMCParticle_hash_contains_2d=table%MCMCParticle_hash_contains(key,value)
      RETURN
      END FUNCTION MCMCParticle_hash_contains_2d

      LOGICAL FUNCTION MCMCParticle_hash_contains_3d(table,key_1,key_2,key_3,value)

      IMPLICIT NONE
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_MCMCParticlehtable), INTENT(IN   ) :: table
      !!! The hashtable
      INTEGER,                          INTENT(IN   ) :: key_1
      INTEGER,                          INTENT(IN   ) :: key_2
      INTEGER,                          INTENT(IN   ) :: key_3
      !!! Input key, which the corresponding value is asked for
      TYPE(MCMCParticle),               INTENT(IN   ) :: value
      !!! Value corresponding to the input key

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: key

      key=IndexHashFunctor64_3d(key_1,key_2,key_3)
      MCMCParticle_hash_contains_3d=table%MCMCParticle_hash_contains(key,value)
      RETURN
      END FUNCTION MCMCParticle_hash_contains_3d

      !TODO check this sub
      !Yaser
      !This function PROBABLY suffers from a bug!
      SUBROUTINE MCMCParticle_hash_remove(table, key, info,existed)
      !!! Given the key, removes the elements in the hash table. Info is
      !!! set to -1 if the key was NOT found.
      !!!
      !!! [NOTE]
      !!! This routine needs to be very fast, therefor we skip the usual
      !!! chit-chat and get right to it. (-> no substart,substop)
      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_MCMCParticlehtable)       :: table
      !!! The hashtable
      INTEGER(ppm_kind_int64), INTENT(IN   ) :: key
      !!! Key to be removed
      INTEGER,                 INTENT(  OUT) :: info
      !!! Info for whether removal was successful or not. 0 if SUCCESSFUL
      !!! and -1 otherwise.
      LOGICAL, OPTIONAL,       INTENT(IN   ) :: existed
      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER :: jump
      INTEGER :: spot !,spot0

      info = 0
      IF (PRESENT(existed)) THEN
         IF (existed) THEN
            jump = 0
            ! Keep on searching withing bounds of hash table.
            DO WHILE (jump.LT.table%nrow)
               ! Get the address corresponding to given key
               spot = table%h_key(key, jump)
               ! If an empty slot found ...
               IF (table%keys(spot).EQ.key) THEN
                  !Remove the key and the corresponding value and RETURN.
                  table%borders_pos(spot)=MCMCParticle(-1,zero)
                  table%keys(spot)=htable_null_li
                  RETURN
               ENDIF
               ! If the current slot is occupied, jump to next key that results
               ! in same hash function.
               jump = jump + 1
            ENDDO
            ! If NOT returned within the while-loop, that means the key was NOT found
            info = ppm_error_fatal
         ENDIF
         RETURN
      ELSE
         IF (ANY(key.EQ.table%keys)) THEN
            jump = 0
            ! Keep on searching withing bounds of hash table.
            DO WHILE (jump.LT.table%nrow)
               ! Get the address corresponding to given key
               spot = table%h_key(key, jump)
               ! If an empty slot found ...
               IF (table%keys(spot).EQ.key) THEN
                  !Remove the key and the corresponding value and RETURN.
                  table%borders_pos(spot)=MCMCParticle(-1,zero)
                  table%keys(spot)=htable_null_li
                  RETURN
               ENDIF
               ! If the current slot is occupied, jump to next key that results
               ! in same hash function.
               jump = jump + 1
            ENDDO
            ! If NOT returned within the while-loop, that means the key was NOT found
            info = ppm_error_fatal
         ENDIF
      ENDIF
      RETURN
      END SUBROUTINE MCMCParticle_hash_remove

      SUBROUTINE MCMCParticle_hash_remove_(table, key_, info,existed)

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_MCMCParticlehtable) :: table
      !!! The hashtable
      INTEGER,           INTENT(IN   ) :: key_
      !!! Key to be removed
      INTEGER,           INTENT(  OUT) :: info
      !!! Info for whether removal was successful or not. 0 if SUCCESSFUL
      !!! and -1 otherwise.
      LOGICAL, OPTIONAL, INTENT(IN   ) :: existed

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: key

      key=INT(key_,KIND=ppm_kind_int64)
      IF (PRESENT(existed)) THEN
         CALL table%MCMCParticle_hash_remove(key,info,existed)
      ELSE
         CALL table%MCMCParticle_hash_remove(key,info)
      ENDIF
      RETURN
      END SUBROUTINE MCMCParticle_hash_remove_

      SUBROUTINE MCMCParticle_hash_remove__(table, key_, info,existed)

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_MCMCParticlehtable)     :: table
      !!! The hashtable
      INTEGER, DIMENSION(:), INTENT(IN   ) :: key_
      !!! Key to be removed
      INTEGER,               INTENT(  OUT) :: info
      !!! Info for whether removal was successful or not. 0 if SUCCESSFUL
      !!! and -1 otherwise.
      LOGICAL, OPTIONAL,     INTENT(IN   ) :: existed

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: key
      INTEGER                 :: ssize

      ssize=SIZE(key_)
      SELECT CASE (ssize)
      CASE (2)
         key=IndexHashFunctor64_2d(key_)
      CASE (3)
         key=IndexHashFunctor64_3d(key_)
      CASE DEFAULT
         info=ppm_error_fatal
         RETURN
      END SELECT
      IF (PRESENT(existed)) THEN
         CALL table%MCMCParticle_hash_remove(key,info,existed)
      ELSE
         CALL table%MCMCParticle_hash_remove(key,info)
      ENDIF
      RETURN
      END SUBROUTINE MCMCParticle_hash_remove__

      SUBROUTINE MCMCParticle_hash_remove_2d(table,key_1,key_2,info,existed)

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_MCMCParticlehtable) :: table
      !!! The hashtable
      INTEGER,           INTENT(IN   ) :: key_1
      INTEGER,           INTENT(IN   ) :: key_2
      !!! Key to be removed
      INTEGER,           INTENT(  OUT) :: info
      !!! Info for whether removal was successful or not. 0 if SUCCESSFUL
      !!! and -1 otherwise.
      LOGICAL, OPTIONAL, INTENT(IN   ) :: existed

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: key

      key=IndexHashFunctor64_2d(key_1,key_2)
      IF (PRESENT(existed)) THEN
         CALL table%MCMCParticle_hash_remove(key,info,existed)
      ELSE
         CALL table%MCMCParticle_hash_remove(key,info)
      ENDIF
      RETURN
      END SUBROUTINE MCMCParticle_hash_remove_2d

      SUBROUTINE MCMCParticle_hash_remove_3d(table,key_1,key_2,key_3,info,existed)

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_MCMCParticlehtable) :: table
      !!! The hashtable
      INTEGER,           INTENT(IN   ) :: key_1
      INTEGER,           INTENT(IN   ) :: key_2
      INTEGER,           INTENT(IN   ) :: key_3
      !!! Key to be removed
      INTEGER,           INTENT(  OUT) :: info
      !!! Info for whether removal was successful or not. 0 if SUCCESSFUL
      !!! and -1 otherwise.
      LOGICAL, OPTIONAL, INTENT(IN   ) :: existed

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: key

      key=IndexHashFunctor64_3d(key_1,key_2,key_3)
      IF (PRESENT(existed)) THEN
         CALL table%MCMCParticle_hash_remove(key,info,existed)
      ELSE
         CALL table%MCMCParticle_hash_remove(key,info)
      ENDIF
      RETURN
      END SUBROUTINE MCMCParticle_hash_remove_3d

      SUBROUTINE MCMCParticle_grow_htable(table,info,key_,value_)
      !!! Based on the number of rows of the table, creates the hash table with
      !!! double size.
      !!!
      !!! [WARNING]
      !!! If you allocate a hashtable with more than 2^31-1 elements the hash
      !!! function will most probably produce incorrect hash keys and fail.
      USE ppm_module_data
      USE ppm_module_alloc
      USE ppm_module_error
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_util_qsort
      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      CLASS(ppm_rc_MCMCParticlehtable)                 :: table
      !!! The hashtable to grow.

      INTEGER,                           INTENT(  OUT) :: info
      INTEGER(ppm_kind_int64), OPTIONAL, INTENT(IN   ) :: key_
      !!! Key to be stored

      TYPE(MCMCParticle),      OPTIONAL, INTENT(IN   ) :: value_
      !!! Value that corresponds to given key
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      TYPE(MCMCParticle),      DIMENSION(:), ALLOCATABLE :: borders_pos_tmp

      REAL(ppm_kind_double) :: t0

      INTEGER(ppm_kind_int64), DIMENSION(:), ALLOCATABLE :: keys_tmp
      INTEGER,                 DIMENSION(:), POINTER     :: ranklist
      INTEGER                                            :: nsize
      INTEGER                                            :: i,j

      CHARACTER(LEN=ppm_char) :: caller='MCMCParticle_grow_htable'

      CALL substart(caller,t0,info)

      nsize=table%nrow

      ALLOCATE(keys_tmp(nsize),SOURCE=table%keys,STAT=info)
      or_fail_alloc("keys_tmp")

      ALLOCATE(borders_pos_tmp(nsize),SOURCE=table%borders_pos,STAT=info)
      or_fail_alloc("keys_tmp & borders_pos_tmp")

      NULLIFY(ranklist)
      CALL ppm_util_qsort(keys_tmp,ranklist,info,nsize)
      or_fail("ppm_util_qsort")

      nsize=table%nrow*2

      IF (nsize.GE.ppm_big_i-1) THEN
         !TOCHCECK
         fail("hashtable with more than 2^31-1 elements will fail",ppm_error=ppm_error_fatal)
      ENDIF

      CALL table%destroy(info)
      or_fail("table%destroy")

      CALL table%create(nsize,info)
      or_fail("table%create")

      DO i=1,SIZE(keys_tmp)
         j=ranklist(i)
         IF (keys_tmp(j).EQ.htable_null_li) CYCLE
         CALL table%MCMCParticle_hash_insert(keys_tmp(j),borders_pos_tmp(j),info)
      ENDDO

      DEALLOCATE(keys_tmp,borders_pos_tmp,STAT=info)
      or_fail_dealloc("keys_tmp & borders_pos_tmp")

      !---------------------------------------------------------------------
      !  Deallocate ranklist array.
      !---------------------------------------------------------------------
      CALL ppm_alloc(ranklist,(/0/),ppm_param_dealloc,info)
      or_fail_dealloc('htable_keys',ppm_error=ppm_error_fatal)

      IF (PRESENT(key_)) THEN
         IF (PRESENT(value_)) THEN
            CALL table%MCMCParticle_hash_insert(key_,value_,info)
         ENDIF
      ENDIF

      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      END SUBROUTINE MCMCParticle_grow_htable

      SUBROUTINE MCMCParticle_shrink_htable(table,info,shrinkage_ratio)
      !!! Based on the number of rows of the table, shrinks the hash table with
      !!! half a size or more.
      USE ppm_module_data
      USE ppm_module_alloc
      USE ppm_module_error
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_util_qsort
      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      CLASS(ppm_rc_MCMCParticlehtable) :: table
      !!! The hashtable to shrink.

      INTEGER,           INTENT(  OUT) :: info
      INTEGER, OPTIONAL, INTENT(IN   ) :: shrinkage_ratio
      !!! OPTIONAL shrinkage_ratio (positive value).
      !!! If the size of hash table is shrinkage_ratio times bigger than the
      !!! real elements inside table, we reduce the table size
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      TYPE(MCMCParticle),      DIMENSION(:), ALLOCATABLE :: borders_pos_tmp

      REAL(ppm_kind_double) :: t0

      INTEGER(ppm_kind_int64), DIMENSION(:), ALLOCATABLE :: keys_tmp
      INTEGER,                 DIMENSION(:), POINTER     :: ranklist
      INTEGER                                            :: shrinkage_ratio_
      INTEGER                                            :: nsize,ssize
      INTEGER                                            :: i,j

      CHARACTER(LEN=ppm_char) :: caller='MCMCParticle_shrink_htable'

      CALL substart(caller,t0,info)

      shrinkage_ratio_=MERGE(shrinkage_ratio,4,PRESENT(shrinkage_ratio))
      shrinkage_ratio_=MERGE(4,shrinkage_ratio_,shrinkage_ratio_.LE.0)

      nsize=table%nrow
      ssize=table%size()

      IF (nsize.GE.shrinkage_ratio_*ssize) THEN
         ALLOCATE(keys_tmp(nsize),SOURCE=table%keys,STAT=info)
         or_fail_alloc("keys_tmp")

         ALLOCATE(borders_pos_tmp(nsize),SOURCE=table%borders_pos,STAT=info)
         or_fail_alloc("keys_tmp & borders_pos_tmp")

         NULLIFY(ranklist)
         CALL ppm_util_qsort(keys_tmp,ranklist,info,nsize)
         or_fail("ppm_util_qsort")

         CALL table%destroy(info)
         or_fail("table%destroy")

         CALL table%create(ssize,info)
         or_fail("table%create")

         DO i=1,nsize
            j=ranklist(i)
            IF (keys_tmp(j).EQ.htable_null_li) CYCLE
            CALL table%insert(keys_tmp(j),borders_pos_tmp(j),info)
         ENDDO

         DEALLOCATE(keys_tmp,borders_pos_tmp,STAT=info)
         or_fail_dealloc("keys_tmp & borders_pos_tmp")

         !---------------------------------------------------------------------
         !  Deallocate ranklist array.
         !---------------------------------------------------------------------
         CALL ppm_alloc(ranklist,(/0/),ppm_param_dealloc,info)
         or_fail_dealloc('htable_keys',ppm_error=ppm_error_fatal)
      ENDIF

      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      END SUBROUTINE MCMCParticle_shrink_htable

      FUNCTION MCMCParticle_hash_size(table) RESULT(value)
      IMPLICIT NONE
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_MCMCParticlehtable) :: table
      !!! The hashtable. The pointer must not be NULL
      INTEGER                          :: value
      !!! Value corresponding to the number of values in table

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      IF (table%nrow.GT.0) THEN
         value=COUNT(table%keys.NE.htable_null_li)
      ELSE
         value=0
      ENDIF
      RETURN
      END FUNCTION MCMCParticle_hash_size

      FUNCTION MCMCParticle_hash_packproposal(table,proposal,proposalSize) RESULT(info)
      USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_fatal
      USE ppm_module_substart, ONLY : substart
      USE ppm_module_substop, ONLY : substop
      USE ppm_module_error, ONLY : ppm_error,ppm_err_alloc,ppm_err_dealloc
      IMPLICIT NONE
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_MCMCParticlehtable)                   :: table
      !!! The hashtable. The pointer must not be NULL
      REAL(MK), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: proposal

      INTEGER,                             INTENT(IN   ) :: proposalSize
      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      REAL(ppm_kind_double) :: t0

      INTEGER :: i,nsize
      INTEGER :: info

      CHARACTER(LEN=*), PARAMETER :: caller="MCMCParticle_hash_packproposal"

      !-------------------------------------------------------------------------
      !  Initialize
      !-------------------------------------------------------------------------
      CALL substart(caller,t0,info)

      IF (table%nrow.GT.0) THEN
         IF (proposalSize.NE.SIZE(proposal)) THEN
            DEALLOCATE(proposal,STAT=info)
            or_fail_dealloc("proposal",ppm_error=ppm_error_fatal)
            ALLOCATE(proposal(proposalSize),STAT=info)
            or_fail_alloc("proposal",ppm_error=ppm_error_fatal)
         ENDIF
         nsize=0
         DO i=1,table%nrow
            IF (table%keys(i).NE.htable_null_li) THEN
               nsize=nsize+1
               proposal(nsize)=table%borders_pos(i)%proposal
            ENDIF
         ENDDO
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      9999 CONTINUE
        CALL substop(caller,t0,info)
        RETURN
      END FUNCTION MCMCParticle_hash_packproposal

      FUNCTION MCMCParticle_hash_elementAt(table,elementIndex,elementKey) RESULT(value)
      !!! This function returns the elementIndex element in the hash table
      !!! If the elementKey is present it gets the coordinates from Hash key
      IMPLICIT NONE
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_MCMCParticlehtable)               :: table
      !!! The hashtable. The pointer must not be NULL
      INTEGER,                         INTENT(IN   ) :: elementIndex
      INTEGER, DIMENSION(:), OPTIONAL, INTENT(  OUT) :: elementKey
      !!! coordinates which has been hashed as a hash key

      TYPE(MCMCParticle)                             :: Value
      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER :: i,elementIndex_
      !-------------------------------------------------------------------------
      !  Initialize
      !-------------------------------------------------------------------------
      IF (table%nrow.GT.0) THEN
         IF (elementIndex.LE.table%nrow) THEN
            elementIndex_=0
            DO i=1,table%nrow
               IF (table%keys(i).NE.htable_null_li) THEN
                  elementIndex_=elementIndex_+1
                  IF (elementIndex_.EQ.elementIndex) THEN
                     value=table%borders_pos(i)
                     IF (PRESENT(elementKey)) THEN
                        SELECT CASE (SIZE(elementKey))
                        CASE (2)
                           CALL HashIndexFunctor_2d(table%keys(i),elementKey)
                        CASE (3)
                           CALL HashIndexFunctor_3d(table%keys(i),elementKey)
                        END SELECT
                     ENDIF
                     RETURN
                  ENDIF !elementIndex_.EQ.elementIndex
               ENDIF !table%keys(i).NE.htable_null_li
            ENDDO !i=1,table%nrow
         ENDIF !elementIndex.LE.table%nrow
      ENDIF !table%nrow.GT.0
      Value=MCMCParticle(-1,zero)
      IF (PRESENT(elementKey)) THEN
         elementKey=-1
      ENDIF
      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      RETURN
      END FUNCTION MCMCParticle_hash_elementAt

      SUBROUTINE MCMCHistoryParticle_create_htable(table,nelement,info)
      !!! Given number of rows of the table, creates the hash table. Number of
      !!! rows will be greater than nelement, as we use the first value that is
      !!! power of 2 and greater than nelement.
      !!!
      !!! [WARNING]
      !!! If you allocate a hashtable with more than 2^31-1 elements the hash
      !!! function will most probably produce incorrect hash keys and fail.
      USE ppm_module_data
      USE ppm_module_alloc
      USE ppm_module_error
      USE ppm_module_substart
      USE ppm_module_substop
      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      CLASS(ppm_rc_MCMCHistoryParticlehtable) :: table
      !!! The hashtable to create.
      INTEGER,                  INTENT(IN   ) :: nelement
      !!! Number of desired elements.
      INTEGER,                  INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double) :: t0

      INTEGER :: i

      CHARACTER(LEN=ppm_char) :: caller='MCMCHistoryParticle_create_htable'

      CALL substart(caller,t0,info)

      IF (nelement.LE.0) THEN
         table%nrow = 1
      ELSE
         table%nrow = 2**(CEILING(LOG(REAL(nelement))/LOG(2.0)))
      ENDIF

      !---------------------------------------------------------------------
      !  Allocate array for hash table keys and array for positions on "borders" array.
      !---------------------------------------------------------------------
      ALLOCATE(table%keys(table%nrow),table%borders_pos(table%nrow),STAT=Info)
      or_fail_alloc('Failed to alocate htable_keys & htable_borders_pos!',ppm_error=ppm_error_fatal)
      !---------------------------------------------------------------------
      !  Set everything to NULL.
      !---------------------------------------------------------------------
      FORALL (i=1:table%nrow)
         table%keys(i)        = htable_null_li
         table%borders_pos(i) = MCMCHistoryParticle(-1,-1,zero,.FALSE.)
      END FORALL

      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      END SUBROUTINE MCMCHistoryParticle_create_htable

      SUBROUTINE MCMCHistoryParticle_destroy_htable(table,info)

      USE ppm_module_data
      USE ppm_module_alloc
      USE ppm_module_error
      USE ppm_module_substart
      USE ppm_module_substop
      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      CLASS(ppm_rc_MCMCHistoryParticlehtable) :: table
      !!! The hashtable to create. The pointer must not be NULL
      INTEGER,                  INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double) :: t0

      CHARACTER(LEN=ppm_char) :: caller='MCMCHistoryParticle_destroy_htable'

      CALL substart(caller,t0,info)

      !---------------------------------------------------------------------
      !  Deallocate array for hash table keys and array for positions on "borders" array.
      !---------------------------------------------------------------------
      DEALLOCATE(table%keys,table%borders_pos,STAT=info)
      or_fail_dealloc('Failed to deallocate htable_keys & htable_borders_pos!',ppm_error=ppm_error_fatal)

      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      END SUBROUTINE MCMCHistoryParticle_destroy_htable

      ELEMENTAL FUNCTION MCMCHistoryParticle_h_func(table,key,seed) RESULT(hash_val)
        IMPLICIT NONE
        !---------------------------------------------------------------------
        !  Arguments
        !---------------------------------------------------------------------
        CLASS(ppm_rc_MCMCHistoryParticlehtable), INTENT(IN   ) :: table
        !!! The hashtable to create. The pointer must not be NULL
        INTEGER(ppm_kind_int64),                 INTENT(IN   ) :: key
        !!! Input key
        INTEGER(ppm_kind_int64),                 INTENT(IN   ) :: seed
        !!! Seed to be used for mixing
        INTEGER(ppm_kind_int64)                                :: hash_val
        !!! Result of the hash function
        !!!
        !!! [NOTE]
        !!! we need to work with 64 bit integers because
        !!! jump*h_func might overflow, and as there is no unsigned
        !!! integer the resulting value might be negative and not
        !!! a valid array index

        !---------------------------------------------------------------------
        !  Local variables and parameters
        !---------------------------------------------------------------------
        INTEGER(ppm_kind_int64), PARAMETER  :: m = 1431374979_ppm_kind_int64
        INTEGER                             :: h
        INTEGER(ppm_kind_int64)             :: data
        INTEGER                             :: k
        INTEGER                             :: len

        len = 4

        h = IEOR(seed, len)
        data = key

        DO WHILE (len .GE. 4)
            k = IBITS(data, 0, 8) !data, pos, len. len = 1 always!
            k = IOR(k, ISHFT(IBITS(data,  8, 8),  8))
            k = IOR(k, ISHFT(IBITS(data, 16, 8), 16))
            k = IOR(k, ISHFT(IBITS(data, 24, 8), 24))

            k = k*m
            k = IEOR(k, ISHFT(k, -24))
            k = k*m

            h = h*m
            h = IEOR(h, k)

            data = data + 4
            len  = len - 1
        ENDDO

        SELECT CASE (len)
        CASE (3)
           h = IEOR(h, ISHFT(IBITS(data, 16, 8), 16))
           h = IEOR(h, ISHFT(IBITS(data,  8, 8),  8))
           h = IEOR(h, IBITS(data, 0, 8))
           h = h*m

        CASE (2)
           h = IEOR(h, ISHFT(IBITS(data,  8, 8),  8))
           h = IEOR(h, IBITS(data, 0, 8))
           h = h*m

        CASE (1)
           h = IEOR(h, IBITS(data, 0, 8))
           h = h*m

        END SELECT

        h = IEOR(h, ISHFT(h, -13))
        h = h*m
        h = IEOR(h, ISHFT(h, -15))
        hash_val = IAND(h, table%nrow - 1)
        RETURN
      END FUNCTION

      ELEMENTAL FUNCTION MCMCHistoryParticle_h_key(table,key,jump) RESULT(address)
      !!! Given the key and jump value, returns corresponding address on
      !!! "borders" array.

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_MCMCHistoryParticlehtable), INTENT(IN   ) :: table
      !!! The hashtable to create. The pointer must not be NULL
      INTEGER(ppm_kind_int64),                 INTENT(IN   ) :: key
      !!! Input key, which corresponds to address requested
      INTEGER,                                 INTENT(IN   ) :: jump
      !!! Jump value for double hashing
      INTEGER                                                :: address
      !!! Address that corresponds to given key
      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: int_addr
      ! we need to work with 64 bit integers because
      ! jump*h_func might overflow, and as there is no unsigned
      ! integer the resulting value might be negative and not
      ! a valid array index

      int_addr = 1_ppm_kind_int64 + &
      & MOD((table%h_func(key,seed1)+jump*table%h_func(key,seed2)),table%nrow)
      address  = INT(int_addr)
      RETURN
      END FUNCTION MCMCHistoryParticle_h_key

      SUBROUTINE MCMCHistoryParticle_hash_insert(table,key,value,info)
      !!! Given the key and the value, stores both in the hash table. Info is
      !!! set to -1 if size of the hash table is not sufficient.
      !!!
      !!! [NOTE]
      !!! This routine needs to be very fast, therefor we skip the usual
      !!! chit-chat and get right to it. (-> no substart,substop unless
      !!! compiled with __DEBUG flag)
      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_MCMCHistoryParticlehtable)  :: table
      !!! The hashtable
      INTEGER(ppm_kind_int64),   INTENT(IN   ) :: key
      !!! Key to be stored
      TYPE(MCMCHistoryParticle), INTENT(IN   ) :: value
      !!! Value that corresponds to given key
      INTEGER,                   INTENT(  OUT) :: info
      !!! Info for whether insertion was successful or not. 0 if SUCCESSFUL
      !!! and -1 otherwise.

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER :: jump
      INTEGER :: spot
      info = 0
      jump = 0
      ! Keep on searching withing bounds of hash table.
      DO WHILE (jump.LT.table%nrow)
         ! Get the address corresponding to given key
         spot = table%h_key(key,jump)
         ! If an empty slot found ...
         IF (table%keys(spot).EQ.htable_null_li) THEN
            ! Store the key and the corresponding value and RETURN.
            table%keys(spot) = key
            table%borders_pos(spot) = value
            RETURN
         !If the key is the same the value should be updated
         ELSE IF (table%keys(spot).EQ.key) THEN
            table%borders_pos(spot) = value
            RETURN
         ENDIF
         ! If the current slot is occupied, jump to next key that results
         ! in same hash function.
         jump = jump + 1
      ENDDO

      ! If NOT returned within the while-loop, that means our hash table is
      ! not sufficiently large. Hence, regrowth will be needed!
      CALL table%grow(info,key,value)
      RETURN
      END SUBROUTINE MCMCHistoryParticle_hash_insert

      SUBROUTINE MCMCHistoryParticle_hash_insert_(table,key_,value,info)

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_MCMCHistoryParticlehtable)  :: table
      !!! The hashtable
      INTEGER,                   INTENT(IN   ) :: key_
      !!! Key to be stored
      TYPE(MCMCHistoryParticle), INTENT(IN   ) :: value
      !!! Value that corresponds to given key
      INTEGER,                   INTENT(  OUT) :: info
      !!! Info for whether insertion was successful or not. 0 if SUCCESSFUL
      !!! and -1 otherwise.

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: key

      key=INT(key_,KIND=ppm_kind_int64)
      CALL table%MCMCHistoryParticle_hash_insert(key,value,info)
      RETURN
      END SUBROUTINE MCMCHistoryParticle_hash_insert_

      SUBROUTINE MCMCHistoryParticle_hash_insert__(table,key_,value,info)

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_MCMCHistoryParticlehtable) :: table
      !!! The hashtable
      INTEGER, DIMENSION(:),     INTENT(IN   ) :: key_
      !!! Key to be stored
      TYPE(MCMCHistoryParticle), INTENT(IN   ) :: value
      !!! Value that corresponds to given key
      INTEGER,                   INTENT(  OUT) :: info
      !!! Info for whether insertion was successful or not. 0 if SUCCESSFUL
      !!! and -1 otherwise.

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: key
      INTEGER                 :: ssize

      ssize=SIZE(key_)
      SELECT CASE (ssize)
      CASE (2)
         key=IndexHashFunctor64_2d(key_)
      CASE (3)
         key=IndexHashFunctor64_3d(key_)
      CASE DEFAULT
         info=ppm_error_fatal
         RETURN
      END SELECT
      CALL table%MCMCHistoryParticle_hash_insert(key,value,info)
      RETURN
      END SUBROUTINE MCMCHistoryParticle_hash_insert__

      SUBROUTINE MCMCHistoryParticle_hash_insert_2d(table,key_1,key_2,value,info)

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_MCMCHistoryParticlehtable)  :: table
      !!! The hashtable
      INTEGER,                   INTENT(IN   ) :: key_1
      INTEGER,                   INTENT(IN   ) :: key_2
      !!! Key to be stored
      TYPE(MCMCHistoryParticle), INTENT(IN   ) :: value
      !!! Value that corresponds to given key
      INTEGER,                   INTENT(  OUT) :: info
      !!! Info for whether insertion was successful or not. 0 if SUCCESSFUL
      !!! and -1 otherwise.

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: key
      key=IndexHashFunctor64_2d(key_1,key_2)
      CALL table%MCMCHistoryParticle_hash_insert(key,value,info)
      RETURN
      END SUBROUTINE MCMCHistoryParticle_hash_insert_2d

      SUBROUTINE MCMCHistoryParticle_hash_insert_3d(table,key_1,key_2,key_3,value,info)

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_MCMCHistoryParticlehtable)  :: table
      !!! The hashtable
      INTEGER,                   INTENT(IN   ) :: key_1
      INTEGER,                   INTENT(IN   ) :: key_2
      INTEGER,                   INTENT(IN   ) :: key_3
      !!! Key to be stored
      TYPE(MCMCHistoryParticle), INTENT(IN   ) :: value
      !!! Value that corresponds to given key
      INTEGER,                   INTENT(  OUT) :: info
      !!! Info for whether insertion was successful or not. 0 if SUCCESSFUL
      !!! and -1 otherwise.

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: key

      key=IndexHashFunctor64_3d(key_1,key_2,key_3)
      CALL table%MCMCHistoryParticle_hash_insert(key,value,info)
      RETURN
      END SUBROUTINE MCMCHistoryParticle_hash_insert_3d

      FUNCTION MCMCHistoryParticle_hash_search(table,key) RESULT(value)
      !!! Given the key, searchs the key on the hash table and returns the
      !!! corresponding value.
      !!!
      !!! [NOTE]
      !!! This routine needs to be very fast, therefor we skip the usual
      !!! chit-chat and get right to it. (-> no substart,substop unless
      !!! compiled with __DEBUG flag)
      IMPLICIT NONE
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_MCMCHistoryParticlehtable), INTENT(IN   ) :: table
      !!! The hashtable to create. The pointer must not be NULL
      INTEGER(ppm_kind_int64),                 INTENT(IN   ) :: key
      !!! Input key, which the corresponding value is asked for
      TYPE(MCMCHistoryParticle)                              :: value
      !!! Value corresponding to the input key

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER :: jump
      INTEGER :: spot

      jump = 0
      ! Keep on searching while we don't come across a NULL value or we don't
      ! exceed bounds of hash table.
      loop: DO WHILE(jump .LT. table%nrow)
          ! Get the other key that results in same hash key as for the input
          ! key.
          spot = table%h_key(key, jump)
          ! If key matches ...
          IF (table%keys(spot).EQ.key) THEN
             ! Set the return value and return
             value = table%borders_pos(spot)
             RETURN
          ELSE IF (table%keys(spot).EQ.htable_null_li) THEN
             IF (.NOT.ANY(key.EQ.table%keys)) EXIT loop
          ENDIF
          ! Otherwise, keep on incrementing jump distance
          jump = jump + 1
      ENDDO loop
      value = MCMCHistoryParticle(-1,-1,zero,.FALSE.)
      RETURN
      END FUNCTION MCMCHistoryParticle_hash_search

      FUNCTION MCMCHistoryParticle_hash_search_(table,key_) RESULT(value)

      IMPLICIT NONE
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_MCMCHistoryParticlehtable), INTENT(IN   ) :: table
      !!! The hashtable
      INTEGER,                                 INTENT(IN   ) :: key_
      !!! Input key, which the corresponding value is asked for
      TYPE(MCMCHistoryParticle)                              :: value
      !!! Value corresponding to the input key

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: key

      key=INT(key_,KIND=ppm_kind_int64)
      value=table%MCMCHistoryParticle_hash_search(key)
      RETURN
      END FUNCTION MCMCHistoryParticle_hash_search_

      FUNCTION MCMCHistoryParticle_hash_search_2d(table,key_1,key_2) RESULT(value)

      IMPLICIT NONE
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_MCMCHistoryParticlehtable), INTENT(IN   ) :: table
      !!! The hashtable
      INTEGER,                                 INTENT(IN   ) :: key_1
      INTEGER,                                 INTENT(IN   ) :: key_2
      !!! Input key, which the corresponding value is asked for
      TYPE(MCMCHistoryParticle)                              :: value
      !!! Value corresponding to the input key

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: key

      key=IndexHashFunctor64_2d(key_1,key_2)
      value=table%MCMCHistoryParticle_hash_search(key)
      RETURN
      END FUNCTION MCMCHistoryParticle_hash_search_2d

      FUNCTION MCMCHistoryParticle_hash_search_3d(table,key_1,key_2,key_3) RESULT(value)

      IMPLICIT NONE
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_MCMCHistoryParticlehtable), INTENT(IN   ) :: table
      !!! The hashtable
      INTEGER,                                 INTENT(IN   ) :: key_1
      INTEGER,                                 INTENT(IN   ) :: key_2
      INTEGER,                                 INTENT(IN   ) :: key_3
      !!! Input key, which the corresponding value is asked for
      TYPE(MCMCHistoryParticle)                              :: value
      !!! Value corresponding to the input key

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: key

      key=IndexHashFunctor64_3d(key_1,key_2,key_3)
      value=table%MCMCHistoryParticle_hash_search(key)
      RETURN
      END FUNCTION MCMCHistoryParticle_hash_search_3d

      !TODO check this sub
      !Yaser
      !This function PROBABLY suffers from a bug!
      SUBROUTINE MCMCHistoryParticle_hash_remove(table, key, info,existed)
      !!! Given the key, removes the elements in the hash table. Info is
      !!! set to -1 if the key was NOT found.
      !!!
      !!! [NOTE]
      !!! This routine needs to be very fast, therefor we skip the usual
      !!! chit-chat and get right to it. (-> no substart,substop unless
      !!! compiled with __DEBUG flag)
      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_MCMCHistoryParticlehtable) :: table
      !!! The hashtable
      INTEGER(ppm_kind_int64),  INTENT(IN   ) :: key
      !!! Key to be removed
      INTEGER,                  INTENT(  OUT) :: info
      !!! Info for whether removal was successful or not. 0 if SUCCESSFUL
      !!! and -1 otherwise.
      LOGICAL, OPTIONAL,        INTENT(IN   ) :: existed

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER :: jump
      INTEGER :: spot !,spot0

      info = 0

      IF (PRESENT(existed)) THEN
         IF (existed) THEN
            jump = 0
            ! Keep on searching withing bounds of hash table.
            DO WHILE (jump.LT.table%nrow)
               ! Get the address corresponding to given key
               spot = table%h_key(key, jump)
               ! If an empty slot found ...
               IF (table%keys(spot).EQ.key) THEN
                  !Remove the key and the corresponding value and RETURN.
                  table%borders_pos(spot)=MCMCHistoryParticle(-1,-1,zero,.FALSE.)
                  table%keys(spot)=htable_null_li
                  RETURN
               ENDIF
               ! If the current slot is occupied, jump to next key that results
               ! in same hash function.
               jump = jump + 1
            ENDDO
            ! If NOT returned within the while-loop, that means the key was NOT found
            info = ppm_error_fatal
         ENDIF
         RETURN
      ELSE
         IF (ANY(key.EQ.table%keys)) THEN
            jump = 0
            ! Keep on searching withing bounds of hash table.
            DO WHILE (jump.LT.table%nrow)
               ! Get the address corresponding to given key
               spot = table%h_key(key, jump)
               ! If an empty slot found ...
               IF (table%keys(spot).EQ.key) THEN
                  !Remove the key and the corresponding value and RETURN.
                  table%borders_pos(spot)=MCMCHistoryParticle(-1,-1,zero,.FALSE.)
                  table%keys(spot)=htable_null_li
                  RETURN
               ENDIF
               ! If the current slot is occupied, jump to next key that results
               ! in same hash function.
               jump = jump + 1
            ENDDO
            ! If NOT returned within the while-loop, that means the key was NOT found
            info = ppm_error_fatal
         ENDIF
      ENDIF
      RETURN
      END SUBROUTINE MCMCHistoryParticle_hash_remove

      SUBROUTINE MCMCHistoryParticle_hash_remove_(table, key_, info,existed)

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_MCMCHistoryParticlehtable) :: table
      !!! The hashtable
      INTEGER,                  INTENT(IN   ) :: key_
      !!! Key to be removed
      INTEGER,                  INTENT(  OUT) :: info
      !!! Info for whether removal was successful or not. 0 if SUCCESSFUL
      !!! and -1 otherwise.
      LOGICAL, OPTIONAL,        INTENT(IN   ) :: existed

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: key

      key=INT(key_,KIND=ppm_kind_int64)
      IF (PRESENT(existed)) THEN
         CALL table%MCMCHistoryParticle_hash_remove(key,info,existed)
      ELSE
         CALL table%MCMCHistoryParticle_hash_remove(key,info)
      ENDIF
      RETURN
      END SUBROUTINE MCMCHistoryParticle_hash_remove_

      SUBROUTINE MCMCHistoryParticle_hash_remove__(table, key_, info,existed)

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_MCMCHistoryParticlehtable) :: table
      !!! The hashtable
      INTEGER, DIMENSION(:),    INTENT(IN   ) :: key_
      !!! Key to be removed
      INTEGER,                  INTENT(  OUT) :: info
      !!! Info for whether removal was successful or not. 0 if SUCCESSFUL
      !!! and -1 otherwise.
      LOGICAL, OPTIONAL,        INTENT(IN   ) :: existed

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: key
      INTEGER                 :: ssize

      ssize=SIZE(key_)

      SELECT CASE (ssize)
      CASE (2)
         key=IndexHashFunctor64_2d(key_)
      CASE (3)
         key=IndexHashFunctor64_3d(key_)
      CASE DEFAULT
         info=ppm_error_fatal
         RETURN
      END SELECT
      IF (PRESENT(existed)) THEN
         CALL table%MCMCHistoryParticle_hash_remove(key,info,existed)
      ELSE
         CALL table%MCMCHistoryParticle_hash_remove(key,info)
      ENDIF
      RETURN
      END SUBROUTINE MCMCHistoryParticle_hash_remove__

      SUBROUTINE MCMCHistoryParticle_hash_remove_2d(table,key_1,key_2,info,existed)

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_MCMCHistoryParticlehtable) :: table
      !!! The hashtable
      INTEGER,                  INTENT(IN   ) :: key_1
      INTEGER,                  INTENT(IN   ) :: key_2
      !!! Key to be removed
      INTEGER,                  INTENT(  OUT) :: info
      !!! Info for whether removal was successful or not. 0 if SUCCESSFUL
      !!! and -1 otherwise.
      LOGICAL, OPTIONAL,        INTENT(IN   ) :: existed

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: key

      key=IndexHashFunctor64_2d(key_1,key_2)
      IF (PRESENT(existed)) THEN
         CALL table%MCMCHistoryParticle_hash_remove(key,info,existed)
      ELSE
         CALL table%MCMCHistoryParticle_hash_remove(key,info)
      ENDIF
      RETURN
      END SUBROUTINE MCMCHistoryParticle_hash_remove_2d

      SUBROUTINE MCMCHistoryParticle_hash_remove_3d(table,key_1,key_2,key_3,info,existed)

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_MCMCHistoryParticlehtable) :: table
      !!! The hashtable
      INTEGER,                  INTENT(IN   ) :: key_1
      INTEGER,                  INTENT(IN   ) :: key_2
      INTEGER,                  INTENT(IN   ) :: key_3
      !!! Key to be removed
      INTEGER,                  INTENT(  OUT) :: info
      !!! Info for whether removal was successful or not. 0 if SUCCESSFUL
      !!! and -1 otherwise.
      LOGICAL, OPTIONAL,        INTENT(IN   ) :: existed

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: key

      key=IndexHashFunctor64_3d(key_1,key_2,key_3)
      IF (PRESENT(existed)) THEN
         CALL table%MCMCHistoryParticle_hash_remove(key,info,existed)
      ELSE
         CALL table%MCMCHistoryParticle_hash_remove(key,info)
      ENDIF
      RETURN
      END SUBROUTINE MCMCHistoryParticle_hash_remove_3d

      SUBROUTINE MCMCHistoryParticle_grow_htable(table,info,key_,value_)
      !!! Based on the number of rows of the table, creates the hash table with
      !!! double size.
      !!!
      !!! [WARNING]
      !!! If you allocate a hashtable with more than 2^31-1 elements the hash
      !!! function will most probably produce incorrect hash keys and fail.
      USE ppm_module_data
      USE ppm_module_alloc
      USE ppm_module_error
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_util_qsort
      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      CLASS(ppm_rc_MCMCHistoryParticlehtable)            :: table
      !!! The hashtable to grow.

      INTEGER,                             INTENT(  OUT) :: info
      INTEGER(ppm_kind_int64),   OPTIONAL, INTENT(IN   ) :: key_
      !!! Key to be stored

      TYPE(MCMCHistoryParticle), OPTIONAL, INTENT(IN   ) :: value_
      !!! Value that corresponds to given key
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      TYPE(MCMCHistoryParticle), DIMENSION(:), ALLOCATABLE :: borders_pos_tmp

      REAL(ppm_kind_double) :: t0

      INTEGER(ppm_kind_int64), DIMENSION(:), ALLOCATABLE :: keys_tmp
      INTEGER,                 DIMENSION(:), POINTER     :: ranklist
      INTEGER                                            :: nsize
      INTEGER                                            :: i,j

      CHARACTER(LEN=ppm_char) :: caller='MCMCHistoryParticle_grow_htable'

      CALL substart(caller,t0,info)

      nsize=table%nrow

      ALLOCATE(keys_tmp(nsize),SOURCE=table%keys,STAT=info)
      or_fail_alloc("keys_tmp")

      ALLOCATE(borders_pos_tmp(nsize),SOURCE=table%borders_pos,STAT=info)
      or_fail_alloc("keys_tmp & borders_pos_tmp")

      NULLIFY(ranklist)
      CALL ppm_util_qsort(keys_tmp,ranklist,info,nsize)
      or_fail("ppm_util_qsort")

      nsize=table%nrow*2

      IF (nsize.GE.ppm_big_i-1) THEN
         !TOCHCECK
         fail("hashtable with more than 2^31-1 elements will fail",ppm_error=ppm_error_fatal)
      ENDIF

      CALL table%destroy(info)
      or_fail("table%destroy")

      CALL table%create(nsize,info)
      or_fail("table%create")

      DO i=1,SIZE(keys_tmp)
         j=ranklist(i)
         IF (keys_tmp(j).EQ.htable_null_li) CYCLE
         CALL table%MCMCHistoryParticle_hash_insert(keys_tmp(j),borders_pos_tmp(j),info)
      ENDDO

      DEALLOCATE(keys_tmp,borders_pos_tmp,STAT=info)
      or_fail_dealloc("keys_tmp & borders_pos_tmp")

      !---------------------------------------------------------------------
      !  Deallocate ranklist array.
      !---------------------------------------------------------------------
      CALL ppm_alloc(ranklist,(/0/),ppm_param_dealloc,info)
      or_fail_dealloc('htable_keys',ppm_error=ppm_error_fatal)

      IF (PRESENT(key_)) THEN
         IF (PRESENT(value_)) THEN
            CALL table%insert(key_,value_,info)
         ENDIF
      ENDIF

      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      END SUBROUTINE MCMCHistoryParticle_grow_htable

      SUBROUTINE MCMCHistoryParticle_shrink_htable(table,info,shrinkage_ratio)
      !!! Based on the number of rows of the table, shrinks the hash table with
      !!! half a size or more.
      USE ppm_module_data
      USE ppm_module_alloc
      USE ppm_module_error
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_util_qsort
      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      CLASS(ppm_rc_MCMCHistoryParticlehtable) :: table
      !!! The hashtable to shrink.

      INTEGER,           INTENT(  OUT) :: info
      INTEGER, OPTIONAL, INTENT(IN   ) :: shrinkage_ratio
      !!! OPTIONAL shrinkage_ratio (positive value).
      !!! If the size of hash table is shrinkage_ratio times bigger than the
      !!! real elements inside table, we reduce the table size
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      TYPE(MCMCHistoryParticle), DIMENSION(:), ALLOCATABLE :: borders_pos_tmp

      REAL(ppm_kind_double) :: t0

      INTEGER(ppm_kind_int64), DIMENSION(:), ALLOCATABLE :: keys_tmp
      INTEGER,                 DIMENSION(:), POINTER     :: ranklist
      INTEGER                                            :: shrinkage_ratio_
      INTEGER                                            :: nsize,ssize
      INTEGER                                            :: i,j

      CHARACTER(LEN=ppm_char) :: caller='MCMCHistoryParticle_shrink_htable'

      CALL substart(caller,t0,info)

      shrinkage_ratio_=MERGE(shrinkage_ratio,4,PRESENT(shrinkage_ratio))
      shrinkage_ratio_=MERGE(4,shrinkage_ratio_,shrinkage_ratio_.LE.0)

      nsize=table%nrow
      ssize=table%size()

      IF (nsize.GE.shrinkage_ratio_*ssize) THEN
         ALLOCATE(keys_tmp(nsize),SOURCE=table%keys,STAT=info)
         or_fail_alloc("keys_tmp")

         ALLOCATE(borders_pos_tmp(nsize),SOURCE=table%borders_pos,STAT=info)
         or_fail_alloc("keys_tmp & borders_pos_tmp")

         NULLIFY(ranklist)
         CALL ppm_util_qsort(keys_tmp,ranklist,info,nsize)
         or_fail("ppm_util_qsort")

         CALL table%destroy(info)
         or_fail("table%destroy")

         CALL table%create(ssize,info)
         or_fail("table%create")

         DO i=1,nsize
            j=ranklist(i)
            IF (keys_tmp(j).EQ.htable_null_li) CYCLE
            CALL table%insert(keys_tmp(j),borders_pos_tmp(j),info)
         ENDDO

         DEALLOCATE(keys_tmp,borders_pos_tmp,STAT=info)
         or_fail_dealloc("keys_tmp & borders_pos_tmp")

         !---------------------------------------------------------------------
         !  Deallocate ranklist array.
         !---------------------------------------------------------------------
         CALL ppm_alloc(ranklist,(/0/),ppm_param_dealloc,info)
         or_fail_dealloc('htable_keys',ppm_error=ppm_error_fatal)
      ENDIF

      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      END SUBROUTINE MCMCHistoryParticle_shrink_htable

      FUNCTION MCMCHistoryParticle_hash_size(table) RESULT (value)

      IMPLICIT NONE
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_MCMCHistoryParticlehtable)  :: table
      !!! The hashtable
      INTEGER                                  :: value
      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      IF (table%nrow.GT.0) THEN
         value=COUNT(table%keys.NE.htable_null_li)
      ELSE
         value=0
      ENDIF
      RETURN
      END FUNCTION MCMCHistoryParticle_hash_size

      SUBROUTINE HashIndex_create_htable(table,nelement,info)
      !!! Given number of rows of the table, creates the hash table. Number of
      !!! rows will be greater than nelement, as we use the first value that is
      !!! power of 2 and greater than nelement.
      !!!
      !!! [WARNING]
      !!! If you allocate a hashtable with more than 2^31-1 elements the hash
      !!! function will most probably produce incorrect hash keys and fail.
      USE ppm_module_data
      USE ppm_module_alloc
      USE ppm_module_error
      USE ppm_module_substart
      USE ppm_module_substop
      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      CLASS(ppm_rc_HashIndextable) :: table
      !!! The hashtable to create.

      INTEGER,       INTENT(IN   ) :: nelement
      !!! Number of desired elements.
      INTEGER,       INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double) :: t0

      INTEGER :: i

      CHARACTER(LEN=ppm_char) :: caller='HashIndex_create_htable'

      CALL substart(caller,t0,info)

      IF (nelement.LE.0) THEN
         table%nrow = 1
      ELSE
         table%nrow = 2**(CEILING(LOG(REAL(nelement))/LOG(2.0)))
      ENDIF

      !---------------------------------------------------------------------
      !  Allocate array for hash table keys and array for positions on "borders" array.
      !---------------------------------------------------------------------
      ALLOCATE(table%keys(table%nrow),STAT=Info)
      or_fail_alloc('Failed to alocate htable_keys!',ppm_error=ppm_error_fatal)
      !---------------------------------------------------------------------
      !  Set everything to NULL.
      !---------------------------------------------------------------------
      FORALL (i=1:table%nrow) table%keys(i) = htable_null_li

      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      END SUBROUTINE HashIndex_create_htable

      SUBROUTINE HashIndex_destroy_htable(table,info)

      USE ppm_module_data
      USE ppm_module_alloc
      USE ppm_module_error
      USE ppm_module_substart
      USE ppm_module_substop
      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      CLASS(ppm_rc_HashIndextable) :: table
      !!! The hashtable to create. The pointer must not be NULL

      INTEGER,       INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double) :: t0

      CHARACTER(LEN=ppm_char) :: caller='HashIndex_destroy_htable'

      CALL substart(caller,t0,info)

      !---------------------------------------------------------------------
      !  Deallocate array for hash table keys and array for positions on "borders" array.
      !---------------------------------------------------------------------
      DEALLOCATE(table%keys,STAT=info)
      or_fail_dealloc('Failed to deallocate htable_keys!',ppm_error=ppm_error_fatal)

      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      END SUBROUTINE HashIndex_destroy_htable

      ELEMENTAL FUNCTION HashIndex_h_func(table,key,seed) RESULT(hash_val)
        IMPLICIT NONE
        !---------------------------------------------------------------------
        !  Arguments
        !---------------------------------------------------------------------
        CLASS(ppm_rc_HashIndextable), INTENT(IN   ) :: table
        !!! The hashtable to create. The pointer must not be NULL

        INTEGER(ppm_kind_int64),      INTENT(IN   ) :: key
        !!! Input key
        INTEGER(ppm_kind_int64),      INTENT(IN   ) :: seed
        !!! Seed to be used for mixing
        INTEGER(ppm_kind_int64)                     :: hash_val
        !!! Result of the hash function
        !!!
        !!! [NOTE]
        !!! we need to work with 64 bit integers because
        !!! jump*h_func might overflow, and as there is no unsigned
        !!! integer the resulting value might be negative and not
        !!! a valid array index

        !---------------------------------------------------------------------
        !  Local variables and parameters
        !---------------------------------------------------------------------
        INTEGER(ppm_kind_int64), PARAMETER  :: m = 1431374979_ppm_kind_int64
        INTEGER                             :: h
        INTEGER(ppm_kind_int64)             :: data
        INTEGER                             :: k
        INTEGER                             :: len

        len = 4

        h = IEOR(seed, len)
        data = key

        DO WHILE (len .GE. 4)
            k = IBITS(data, 0, 8) !data, pos, len. len = 1 always!
            k = IOR(k, ISHFT(IBITS(data,  8, 8),  8))
            k = IOR(k, ISHFT(IBITS(data, 16, 8), 16))
            k = IOR(k, ISHFT(IBITS(data, 24, 8), 24))

            k = k*m
            k = IEOR(k, ISHFT(k, -24))
            k = k*m

            h = h*m
            h = IEOR(h, k)

            data = data + 4
            len  = len - 1
        ENDDO

        SELECT CASE (len)
        CASE (3)
           h = IEOR(h, ISHFT(IBITS(data, 16, 8), 16))
           h = IEOR(h, ISHFT(IBITS(data,  8, 8),  8))
           h = IEOR(h, IBITS(data, 0, 8))
           h = h*m

        CASE (2)
           h = IEOR(h, ISHFT(IBITS(data,  8, 8),  8))
           h = IEOR(h, IBITS(data, 0, 8))
           h = h*m

        CASE (1)
           h = IEOR(h, IBITS(data, 0, 8))
           h = h*m

        END SELECT

        h = IEOR(h, ISHFT(h, -13))
        h = h*m
        h = IEOR(h, ISHFT(h, -15))
        hash_val = IAND(h, table%nrow - 1)
        RETURN
      END FUNCTION

      ELEMENTAL FUNCTION HashIndex_h_key(table,key,jump) RESULT(address)
      !!! Given the key and jump value, returns corresponding address on
      !!! "borders" array.

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_HashIndextable), INTENT(IN   ) :: table
      !!! The hashtable to create. The pointer must not be NULL

      INTEGER(ppm_kind_int64),      INTENT(IN   ) :: key
      !!! Input key, which corresponds to address requested
      INTEGER,                      INTENT(IN   ) :: jump
      !!! Jump value for double hashing
      INTEGER                                     :: address
      !!! Address that corresponds to given key
      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: int_addr
      ! we need to work with 64 bit integers because
      ! jump*h_func might overflow, and as there is no unsigned
      ! integer the resulting value might be negative and not
      ! a valid array index

      int_addr = 1_ppm_kind_int64 + &
      & MOD((table%h_func(key,seed1)+jump*table%h_func(key,seed2)),table%nrow)
      address  = INT(int_addr)
      RETURN
      END FUNCTION HashIndex_h_key

      SUBROUTINE HashIndex_hash_insert(table,key,info)
      !!! Given the key and the value, stores both in the hash table. Info is
      !!! set to -1 if size of the hash table is not sufficient.
      !!!
      !!! [NOTE]
      !!! This routine needs to be very fast, therefor we skip the usual
      !!! chit-chat and get right to it. (-> no substart,substop unless
      !!! compiled with __DEBUG flag)
      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_HashIndextable)           :: table
      !!! The hashtable

      INTEGER(ppm_kind_int64), INTENT(IN   ) :: key
      !!! Key to be stored
      INTEGER,                 INTENT(  OUT) :: info
      !!! Info for whether insertion was successful or not. 0 if SUCCESSFUL
      !!! and -1 otherwise.

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER :: jump
      INTEGER :: spot

      info = 0
      jump = 0
      ! Keep on searching withing bounds of hash table.
      DO WHILE (jump.LT.table%nrow)
         ! Get the address corresponding to given key
         spot = table%h_key(key,jump)
         ! If an empty slot found ...
         IF (table%keys(spot).EQ.htable_null_li) THEN
            ! Store the key and the corresponding value and RETURN.
            table%keys(spot) = key
            RETURN
         !If the key is the same the value should be updated
         ELSE IF (table%keys(spot).EQ.key) THEN
            RETURN
         ENDIF
         ! If the current slot is occupied, jump to next key that results
         ! in same hash function.
         jump = jump + 1
      ENDDO

      ! If NOT returned within the while-loop, that means our hash table is
      ! not sufficiently large. Hence, regrowth will be needed!
      CALL table%grow(info,key)
      RETURN
      END SUBROUTINE HashIndex_hash_insert

      SUBROUTINE HashIndex_hash_insert_(table,key_,info)

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_HashIndextable) :: table
      !!! The hashtable

      INTEGER,       INTENT(IN   ) :: key_
      !!! Key to be stored
      INTEGER,       INTENT(  OUT) :: info
      !!! Info for whether insertion was successful or not. 0 if SUCCESSFUL
      !!! and -1 otherwise.

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: key

      key=INT(key_,KIND=ppm_kind_int64)
      CALL table%HashIndex_hash_insert(key,info)
      RETURN
      END SUBROUTINE HashIndex_hash_insert_

      SUBROUTINE HashIndex_hash_insert__(table,key_,info)

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_HashIndextable)         :: table
      !!! The hashtable

      INTEGER, DIMENSION(:), INTENT(IN   ) :: key_
      !!! Key to be stored
      INTEGER,               INTENT(  OUT) :: info
      !!! Info for whether insertion was successful or not. 0 if SUCCESSFUL
      !!! and -1 otherwise.

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: key
      INTEGER                 :: ssize

      ssize=SIZE(key_)
      SELECT CASE (ssize)
      CASE (2)
         key=IndexHashFunctor64_2d(key_)
      CASE (3)
         key=IndexHashFunctor64_3d(key_)
      CASE DEFAULT
         info=ppm_error_fatal
         RETURN
      END SELECT
      CALL table%HashIndex_hash_insert(key,info)
      RETURN
      END SUBROUTINE HashIndex_hash_insert__

      SUBROUTINE HashIndex_hash_insert_2d(table,key_1,key_2,info)

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_HashIndextable) :: table
      !!! The hashtable

      INTEGER,       INTENT(IN   ) :: key_1
      INTEGER,       INTENT(IN   ) :: key_2
      !!! Key to be stored
      INTEGER,       INTENT(  OUT) :: info
      !!! Info for whether insertion was successful or not. 0 if SUCCESSFUL
      !!! and -1 otherwise.

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: key
      key=IndexHashFunctor64_2d(key_1,key_2)
      CALL table%HashIndex_hash_insert(key,info)
      RETURN
      END SUBROUTINE HashIndex_hash_insert_2d

      SUBROUTINE HashIndex_hash_insert_3d(table,key_1,key_2,key_3,info)

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_HashIndextable) :: table
      !!! The hashtable

      INTEGER,       INTENT(IN   ) :: key_1
      INTEGER,       INTENT(IN   ) :: key_2
      INTEGER,       INTENT(IN   ) :: key_3
      !!! Key to be stored
      INTEGER,       INTENT(  OUT) :: info
      !!! Info for whether insertion was successful or not. 0 if SUCCESSFUL
      !!! and -1 otherwise.

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: key
      key=IndexHashFunctor64_3d(key_1,key_2,key_3)
      CALL table%HashIndex_hash_insert(key,info)
      RETURN
      END SUBROUTINE HashIndex_hash_insert_3d

      LOGICAL FUNCTION HashIndex_hash_search(table,key)
      !!! Given the key, searchs the key on the hash table and returns the
      !!! corresponding value.
      !!!
      !!! [NOTE]
      !!! This routine needs to be very fast, therefor we skip the usual
      !!! chit-chat and get right to it. (-> no substart,substop unless
      !!! compiled with __DEBUG flag)
      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_HashIndextable), INTENT(IN   ) :: table
      !!! The hashtable to create. The pointer must not be NULL

      INTEGER(ppm_kind_int64),      INTENT(IN   ) :: key
      !!! Input key
      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER :: spot
      INTEGER :: jump

      jump=0
      ! Keep on searching while we don't come across a NULL value or we don't
      ! exceed bounds of hash table.
      loop: DO WHILE(jump .LT. table%nrow)
          ! Get the other key that results in same hash key as for the input
          ! key.
          spot = table%h_key(key,jump)
          ! If key matches ...
          IF (table%keys(spot).EQ.key) THEN
             ! Set the return value and return
             HashIndex_hash_search = .TRUE.
             RETURN
          ELSE IF (table%keys(spot).EQ.htable_null_li) THEN
             IF (ANY(key.EQ.table%keys)) THEN
                HashIndex_hash_search = .TRUE.
                RETURN
             ELSE
                HashIndex_hash_search = .FALSE.
                RETURN
             ENDIF
          ENDIF
          ! Otherwise, keep on incrementing jump distance
          jump = jump + 1
      ENDDO loop
      HashIndex_hash_search = .FALSE.
      RETURN
      END FUNCTION HashIndex_hash_search

      LOGICAL FUNCTION HashIndex_hash_search_(table,key_)

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_HashIndextable), INTENT(IN   ) :: table
      !!! The hashtable

      INTEGER,                      INTENT(IN   ) :: key_
      !!! Input key
      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: key

      key=INT(key_,KIND=ppm_kind_int64)
      HashIndex_hash_search_=table%HashIndex_hash_search(key)
      RETURN
      END FUNCTION HashIndex_hash_search_

      LOGICAL FUNCTION HashIndex_hash_search_2d(table,key_1,key_2)

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_HashIndextable), INTENT(IN   ) :: table
      !!! The hashtable

      INTEGER,                      INTENT(IN   ) :: key_1
      INTEGER,                      INTENT(IN   ) :: key_2
      !!! Input key
      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: key

      key=IndexHashFunctor64_2d(key_1,key_2)
      HashIndex_hash_search_2d=table%HashIndex_hash_search(key)
      RETURN
      END FUNCTION HashIndex_hash_search_2d

      LOGICAL FUNCTION HashIndex_hash_search_3d(table,key_1,key_2,key_3)

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_HashIndextable), INTENT(IN   ) :: table
      !!! The hashtable

      INTEGER,                      INTENT(IN   ) :: key_1
      INTEGER,                      INTENT(IN   ) :: key_2
      INTEGER,                      INTENT(IN   ) :: key_3
      !!! Input key
      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: key

      key=IndexHashFunctor64_3d(key_1,key_2,key_3)
      HashIndex_hash_search_3d=table%HashIndex_hash_search(key)
      RETURN
      END FUNCTION HashIndex_hash_search_3d

      !TODO check this sub
      !Yaser
      !This function PROBABLY suffers from a bug!
      SUBROUTINE HashIndex_hash_remove(table, key, info,existed)
      !!! Given the key, removes the elements in the hash table. Info is
      !!! set to -1 if the key was NOT found.
      !!!
      !!! [NOTE]
      !!! This routine needs to be very fast, therefor we skip the usual
      !!! chit-chat and get right to it. (-> no substart,substop unless
      !!! compiled with __DEBUG flag)
      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_HashIndextable)           :: table
      !!! The hashtable

      INTEGER(ppm_kind_int64), INTENT(IN   ) :: key
      !!! Key to be removed
      INTEGER,                 INTENT(  OUT) :: info
      !!! Info for whether removal was successful or not. 0 if SUCCESSFUL
      !!! and -1 otherwise.

      LOGICAL, OPTIONAL,        INTENT(IN   ) :: existed

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER :: jump
      INTEGER :: spot !,spot0

      info = 0
      IF (PRESENT(existed)) THEN
         IF (existed) THEN
            jump = 0
            ! Keep on searching withing bounds of hash table.
            DO WHILE (jump.LT.table%nrow)
               ! Get the address corresponding to given key
               spot = table%h_key(key, jump)
               ! If an empty slot found ...
               IF (table%keys(spot).EQ.key) THEN
                  !Remove the key
                  table%keys(spot)=htable_null_li
                  RETURN
               ENDIF
               ! If the current slot is occupied, jump to next key that results
               ! in same hash function.
               jump = jump + 1
            ENDDO
            ! If NOT returned within the while-loop, that means the key was NOT found
            info = ppm_error_fatal
         ENDIF
         RETURN
      ELSE
         IF (ANY(key.EQ.table%keys)) THEN
            jump = 0
            ! Keep on searching withing bounds of hash table.
            DO WHILE (jump.LT.table%nrow)
               ! Get the address corresponding to given key
               spot = table%h_key(key, jump)
               ! If an empty slot found ...
               IF (table%keys(spot).EQ.key) THEN
                  !Remove the key
                  table%keys(spot)=htable_null_li
                  RETURN
               ENDIF
               ! If the current slot is occupied, jump to next key that results
               ! in same hash function.
               jump = jump + 1
            ENDDO
            ! If NOT returned within the while-loop, that means the key was NOT found
            info = ppm_error_fatal
         ENDIF
      ENDIF
      RETURN
      END SUBROUTINE HashIndex_hash_remove

      SUBROUTINE HashIndex_hash_remove_(table, key_, info,existed)

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_HashIndextable)     :: table
      !!! The hashtable

      INTEGER,           INTENT(IN   ) :: key_
      !!! Key to be removed
      INTEGER,           INTENT(  OUT) :: info
      !!! Info for whether removal was successful or not. 0 if SUCCESSFUL
      !!! and -1 otherwise.

      LOGICAL, OPTIONAL, INTENT(IN   ) :: existed

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: key

      key=INT(key_,KIND=ppm_kind_int64)
      IF (PRESENT(existed)) THEN
         CALL table%HashIndex_hash_remove(key,info,existed)
      ELSE
         CALL table%HashIndex_hash_remove(key,info)
      ENDIF
      RETURN
      END SUBROUTINE HashIndex_hash_remove_

      SUBROUTINE HashIndex_hash_remove__(table, key_, info,existed)

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_HashIndextable)         :: table
      !!! The hashtable

      INTEGER, DIMENSION(:), INTENT(IN   ) :: key_
      !!! Key to be removed
      INTEGER,               INTENT(  OUT) :: info
      !!! Info for whether removal was successful or not. 0 if SUCCESSFUL
      !!! and -1 otherwise.

      LOGICAL, OPTIONAL,     INTENT(IN   ) :: existed

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: key
      INTEGER                 :: ssize

      ssize=SIZE(key_)

      SELECT CASE (ssize)
      CASE (2)
         key=IndexHashFunctor64_2d(key_)
      CASE (3)
         key=IndexHashFunctor64_3d(key_)
      CASE DEFAULT
         info=ppm_error_fatal
         RETURN
      END SELECT
      IF (PRESENT(existed)) THEN
         CALL table%HashIndex_hash_remove(key,info,existed)
      ELSE
         CALL table%HashIndex_hash_remove(key,info)
      ENDIF
      RETURN
      END SUBROUTINE HashIndex_hash_remove__

      SUBROUTINE HashIndex_hash_remove_2d(table,key_1,key_2,info,existed)

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_HashIndextable)     :: table
      !!! The hashtable

      INTEGER,           INTENT(IN   ) :: key_1
      INTEGER,           INTENT(IN   ) :: key_2
      !!! Key to be removed
      INTEGER,           INTENT(  OUT) :: info
      !!! Info for whether removal was successful or not. 0 if SUCCESSFUL
      !!! and -1 otherwise.

      LOGICAL, OPTIONAL, INTENT(IN   ) :: existed

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: key

      key=IndexHashFunctor64_2d(key_1,key_2)
      IF (PRESENT(existed)) THEN
         CALL table%HashIndex_hash_remove(key,info,existed)
      ELSE
         CALL table%HashIndex_hash_remove(key,info)
      ENDIF
      RETURN
      END SUBROUTINE HashIndex_hash_remove_2d

      SUBROUTINE HashIndex_hash_remove_3d(table,key_1,key_2,key_3,info,existed)

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_HashIndextable)     :: table
      !!! The hashtable

      INTEGER,           INTENT(IN   ) :: key_1
      INTEGER,           INTENT(IN   ) :: key_2
      INTEGER,           INTENT(IN   ) :: key_3
      !!! Key to be removed
      INTEGER,           INTENT(  OUT) :: info
      !!! Info for whether removal was successful or not. 0 if SUCCESSFUL
      !!! and -1 otherwise.

      LOGICAL, OPTIONAL, INTENT(IN   ) :: existed
      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: key

      key=IndexHashFunctor64_3d(key_1,key_2,key_3)
      IF (PRESENT(existed)) THEN
         CALL table%HashIndex_hash_remove(key,info,existed)
      ELSE
         CALL table%HashIndex_hash_remove(key,info)
      ENDIF
      RETURN
      END SUBROUTINE HashIndex_hash_remove_3d

      SUBROUTINE HashIndex_grow_htable(table,info,key_)
      !!! Based on the number of rows of the table, creates the hash table with
      !!! double size.
      !!!
      !!! [WARNING]
      !!! If you allocate a hashtable with more than 2^31-1 elements the hash
      !!! function will most probably produce incorrect hash keys and fail.
      USE ppm_module_data
      USE ppm_module_alloc
      USE ppm_module_error
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_util_qsort
      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      CLASS(ppm_rc_HashIndextable)                     :: table
      !!! The hashtable to grow.

      INTEGER,                           INTENT(  OUT) :: info
      INTEGER(ppm_kind_int64), OPTIONAL, INTENT(IN   ) :: key_
      !!! Key to be stored
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double) :: t0

      INTEGER(ppm_kind_int64), DIMENSION(:), ALLOCATABLE :: keys_tmp
      INTEGER,                 DIMENSION(:), POINTER     :: ranklist
      INTEGER                                            :: nsize
      INTEGER                                            :: i,j

      CHARACTER(LEN=ppm_char) :: caller='HashIndex_grow_htable'

      CALL substart(caller,t0,info)

      nsize=table%nrow

      ALLOCATE(keys_tmp(nsize),SOURCE=table%keys,STAT=info)
      or_fail_alloc("keys_tmp")

      NULLIFY(ranklist)
      CALL ppm_util_qsort(keys_tmp,ranklist,info,nsize)
      or_fail("ppm_util_qsort")

      nsize=table%nrow*2

      IF (nsize.GE.ppm_big_i-1) THEN
         !TOCHCECK
         fail("hashtable with more than 2^31-1 elements will fail",ppm_error=ppm_error_fatal)
      ENDIF

      CALL table%destroy(info)
      or_fail("table%destroy")

      CALL table%create(nsize,info)
      or_fail("table%create")

      DO i=1,SIZE(keys_tmp)
         j=ranklist(i)
         IF (keys_tmp(j).EQ.htable_null_li) CYCLE
         CALL table%HashIndex_hash_insert(keys_tmp(j),info)
      ENDDO

      DEALLOCATE(keys_tmp,STAT=info)
      or_fail_dealloc("keys_tmp")

      !---------------------------------------------------------------------
      !  Deallocate ranklist array.
      !---------------------------------------------------------------------
      CALL ppm_alloc(ranklist,(/0/),ppm_param_dealloc,info)
      or_fail_dealloc('htable_keys',ppm_error=ppm_error_fatal)

      IF (PRESENT(key_)) THEN
         CALL table%insert(key_,info)
      ENDIF

      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      END SUBROUTINE HashIndex_grow_htable

      SUBROUTINE HashIndex_shrink_htable(table,info,shrinkage_ratio)
      !!! Based on the number of rows of the table, shrinks the hash table with
      !!! half a size or more.
      USE ppm_module_data
      USE ppm_module_alloc
      USE ppm_module_error
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_util_qsort
      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      CLASS(ppm_rc_HashIndextable)     :: table
      !!! The hashtable to shrink.

      INTEGER,           INTENT(  OUT) :: info
      INTEGER, OPTIONAL, INTENT(IN   ) :: shrinkage_ratio
      !!! OPTIONAL shrinkage_ratio (positive value).
      !!! If the size of hash table is shrinkage_ratio times bigger than the
      !!! real elements inside table, we reduce the table size
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double) :: t0

      INTEGER(ppm_kind_int64), DIMENSION(:), ALLOCATABLE :: keys_tmp
      INTEGER,                 DIMENSION(:), POINTER     :: ranklist
      INTEGER                                            :: shrinkage_ratio_
      INTEGER                                            :: nsize,ssize
      INTEGER                                            :: i,j

      CHARACTER(LEN=ppm_char) :: caller='HashIndextable_shrink'

      CALL substart(caller,t0,info)

      shrinkage_ratio_=MERGE(shrinkage_ratio,4,PRESENT(shrinkage_ratio))
      shrinkage_ratio_=MERGE(4,shrinkage_ratio_,shrinkage_ratio_.LE.0)

      nsize=table%nrow
      ssize=table%size()

      IF (nsize.GE.shrinkage_ratio_*ssize) THEN
         ALLOCATE(keys_tmp(nsize),SOURCE=table%keys,STAT=info)
         or_fail_alloc("keys_tmp")

         NULLIFY(ranklist)
         CALL ppm_util_qsort(keys_tmp,ranklist,info,nsize)
         or_fail("ppm_util_qsort")

         CALL table%destroy(info)
         or_fail("table%destroy")

         CALL table%create(ssize,info)
         or_fail("table%create")

         DO i=1,nsize
            j=ranklist(i)
            IF (keys_tmp(j).EQ.htable_null_li) CYCLE
            CALL table%insert(keys_tmp(j),info)
         ENDDO

         DEALLOCATE(keys_tmp,STAT=info)
         or_fail_dealloc("keys_tmp")

         !---------------------------------------------------------------------
         !  Deallocate ranklist array.
         !---------------------------------------------------------------------
         CALL ppm_alloc(ranklist,(/0/),ppm_param_dealloc,info)
         or_fail_dealloc('htable_keys',ppm_error=ppm_error_fatal)
      ENDIF

      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      END SUBROUTINE HashIndex_shrink_htable

      FUNCTION HashIndex_hash_size(table) RESULT(value)

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_HashIndextable) :: table
      !!! The hashtable

      INTEGER                      :: value
      !!! Input key
      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      IF (table%nrow.GT.0) THEN
         value=COUNT(table%keys.NE.htable_null_li)
      ELSE
         value=0
      ENDIF
      RETURN
      END FUNCTION HashIndex_hash_size

      !----------------------------------------------------------------------
      !  realloc
      !----------------------------------------------------------------------
      SUBROUTINE MCMCParticle_realloc(Particlehtable,info,nsize)
      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      TYPE(ppm_rc_MCMCParticlehtable), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: Particlehtable

      INTEGER,                                                    INTENT(  OUT) :: info
      INTEGER, OPTIONAL,                                          INTENT(IN   ) :: nsize
      !!! new size minus one for the background region which is indexed as 0
      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      REAL(ppm_kind_double) :: t0

      INTEGER :: ssize,i,j,k,l,nsize_

      CHARACTER(LEN=*), PARAMETER :: caller="MCMCParticle_realloc"

      CALL substart(caller,t0,info)

      IF (ALLOCATED(Particlehtable)) THEN
         ssize=SIZE(Particlehtable)-1
      ELSE
         ssize=15
      ENDIF
      IF (PRESENT(nsize)) THEN
         IF (nsize.LT.ssize) THEN
            l=0
            DO i=1,ssize
               IF (Particlehtable(i)%nrow.GT.0) CYCLE
               l=l+1
            ENDDO
            IF (nsize.LT.ssize-l) THEN
               fail("Some elements will be killed!!!",ppm_error=ppm_error_fatal)
            ENDIF
         ENDIF
         nsize_=nsize
      ELSE
         nsize_=MERGE(15,ssize*2+1,ssize.EQ.15)
      ENDIF

      ALLOCATE(Particlehtabletmp(0:ssize),STAT=info)
      or_fail_alloc("Particlehtabletmp")

      l=Particlehtable(0)%nrow
      IF (l.GT.0) THEN
         CALL Particlehtabletmp(0)%create(l,info)
         or_fail("Particlehtabletmp(0)%create")

         FORALL (i=1:l)
            Particlehtabletmp(0)%keys(i)=Particlehtable(0)%keys(i)
            Particlehtabletmp(0)%borders_pos(i)=Particlehtable(0)%borders_pos(i)
         END FORALL

         CALL Particlehtable(0)%destroy(info)
         or_fail("Particlehtable(0)%destroy")
      ENDIF

      DO i=1,ssize
         l=Particlehtable(i)%nrow
         IF (l.LE.0) CYCLE
         CALL Particlehtabletmp(i)%create(l,info)
         or_fail("Particlehtabletmp(i)%create")

         FORALL (j=1:l)
            Particlehtabletmp(i)%keys(j)=Particlehtable(i)%keys(j)
            Particlehtabletmp(i)%borders_pos(j)=Particlehtable(i)%borders_pos(j)
         END FORALL

         CALL Particlehtable(i)%destroy(info)
         or_fail("Particlehtable(i)%destroy")
      ENDDO

      DEALLOCATE(Particlehtable,STAT=info)
      or_fail_dealloc("Particlehtable")

      ALLOCATE(Particlehtable(0:nsize_),STAT=info)
      or_fail_alloc("Particlehtable")

      l=Particlehtabletmp(0)%nrow
      IF (l.GT.0) THEN
         CALL Particlehtable(0)%create(l,info)
         or_fail("Particlehtable(0)%create")

         FORALL (i=1:l)
            Particlehtable(0)%keys(i)       =Particlehtabletmp(0)%keys(i)
            Particlehtable(0)%borders_pos(i)=Particlehtabletmp(0)%borders_pos(i)
         END FORALL

         CALL Particlehtabletmp(0)%destroy(info)
         or_fail("Particlehtabletmp(0)%destroy")
      ENDIF

      IF (nsize_.LT.ssize) THEN
         !We are reducing the size and get rid off all the extra regions
         k=0
         DO i=1,ssize
            l=Particlehtabletmp(i)%nrow
            IF (l.LE.0) CYCLE
            k=k+1

            CALL Particlehtable(k)%create(l,info)
            or_fail("Particlehtable(k)%create")

            FORALL (j=1:l)
               Particlehtable(k)%keys(j)       =Particlehtabletmp(i)%keys(j)
               Particlehtable(k)%borders_pos(j)=Particlehtabletmp(i)%borders_pos(j)
            END FORALL

            CALL Particlehtabletmp(i)%destroy(info)
            or_fail("Particlehtabletmp(i)%destroy")
         ENDDO
      ELSE
         DO i=1,ssize
            l=Particlehtabletmp(i)%nrow
            IF (l.LE.0) CYCLE
            CALL Particlehtable(i)%create(l,info)
            or_fail("Particlehtable(i)%create")

            FORALL (j=1:l)
               Particlehtable(i)%keys(j)       =Particlehtabletmp(i)%keys(j)
               Particlehtable(i)%borders_pos(j)=Particlehtabletmp(i)%borders_pos(j)
            END FORALL

            CALL Particlehtabletmp(i)%destroy(info)
            or_fail("Particlehtabletmp(i)%destroy")
         ENDDO
      ENDIF

      DEALLOCATE(Particlehtabletmp,STAT=info)
      or_fail_dealloc("Particlehtabletmp")

      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      END SUBROUTINE MCMCParticle_realloc

      SUBROUTINE MCMCHistoryParticle_realloc(HistoryParticlehtable,info,nsize)
      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      TYPE(ppm_rc_MCMCHistoryParticlehtable), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: HistoryParticlehtable

      INTEGER,                                                           INTENT(  OUT) :: info
      INTEGER, OPTIONAL,                                                 INTENT(IN   ) :: nsize
      !!! new size minus one for the background region which is indexed as 0
      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      REAL(ppm_kind_double) :: t0

      INTEGER :: ssize,nsize_
      INTEGER :: i,j,k,l

      CHARACTER(LEN=*), PARAMETER :: caller="MCMCHistoryParticle_realloc"

      CALL substart(caller,t0,info)

      IF (ALLOCATED(HistoryParticlehtable)) THEN
         ssize=SIZE(HistoryParticlehtable)-1
      ELSE
         ssize=15
      ENDIF
      IF (PRESENT(nsize)) THEN
         IF (nsize.LT.ssize) THEN
            l=0
            DO i=1,ssize
               IF (HistoryParticlehtable(i)%nrow.GT.0) CYCLE
               l=l+1
            ENDDO
            IF (nsize.LT.ssize-l) THEN
               fail("Some elements will be killed!!!",ppm_error=ppm_error_fatal)
            ENDIF
         ENDIF
         nsize_=nsize
      ELSE
         nsize_=MERGE(15,ssize*2+1,ssize.EQ.15)
      ENDIF

      ALLOCATE(HistoryParticlehtabletmp(0:ssize),STAT=info)
      or_fail_alloc("HistoryParticlehtabletmp")

      l=HistoryParticlehtable(0)%nrow
      IF (l.GT.0) THEN
         CALL HistoryParticlehtabletmp(0)%create(l,info)
         or_fail("HistoryParticlehtabletmp(0)%create")

         FORALL (i=1:l)
            HistoryParticlehtabletmp(0)%keys(i)=HistoryParticlehtable(0)%keys(i)
            HistoryParticlehtabletmp(0)%borders_pos(i)=HistoryParticlehtable(0)%borders_pos(i)
         END FORALL

         CALL HistoryParticlehtable(0)%destroy(info)
         or_fail("HistoryParticlehtable(0)%destroy")
      ENDIF

      DO i=1,ssize
         l=HistoryParticlehtable(i)%nrow
         IF (l.LE.0) CYCLE
         CALL HistoryParticlehtabletmp(i)%create(l,info)
         or_fail("HistoryParticlehtabletmp(i)%create")

         FORALL (j=1:l)
            HistoryParticlehtabletmp(i)%keys(j)=HistoryParticlehtable(i)%keys(j)
            HistoryParticlehtabletmp(i)%borders_pos(j)=HistoryParticlehtable(i)%borders_pos(j)
         END FORALL

         CALL HistoryParticlehtable(i)%destroy(info)
         or_fail("HistoryParticlehtable(i)%destroy")
      ENDDO

      DEALLOCATE(HistoryParticlehtable,STAT=info)
      or_fail_dealloc("HistoryParticlehtable")

      ALLOCATE(HistoryParticlehtable(0:nsize_),STAT=info)
      or_fail_alloc("HistoryParticlehtable")

      l=HistoryParticlehtabletmp(0)%nrow
      IF (l.GT.0) THEN
         CALL HistoryParticlehtable(0)%create(l,info)
         or_fail("HistoryParticlehtable(0)%create")

         FORALL (i=1:l)
            HistoryParticlehtable(0)%keys(i)       =HistoryParticlehtabletmp(0)%keys(i)
            HistoryParticlehtable(0)%borders_pos(i)=HistoryParticlehtabletmp(0)%borders_pos(i)
         END FORALL

         CALL HistoryParticlehtabletmp(0)%destroy(info)
         or_fail("HistoryParticlehtabletmp(0)%destroy")
      ENDIF

      IF (nsize_.LT.ssize) THEN
         !We are reducing the size and get rid off all the extra regions
         k=0
         DO i=1,ssize
            l=HistoryParticlehtabletmp(i)%nrow
            IF (l.LE.0) CYCLE
            k=k+1

            CALL HistoryParticlehtable(k)%create(l,info)
            or_fail("HistoryParticlehtable(k)%create")

            FORALL (j=1:l)
               HistoryParticlehtable(k)%keys(j)       =HistoryParticlehtabletmp(i)%keys(j)
               HistoryParticlehtable(k)%borders_pos(j)=HistoryParticlehtabletmp(i)%borders_pos(j)
            END FORALL

            CALL HistoryParticlehtabletmp(i)%destroy(info)
            or_fail("HistoryParticlehtabletmp(i)%destroy")
         ENDDO
      ELSE
         DO i=1,ssize
            l=HistoryParticlehtabletmp(i)%nrow
            IF (l.LE.0) CYCLE
            CALL HistoryParticlehtable(i)%create(l,info)
            or_fail("HistoryParticlehtable(i)%create")

            FORALL (j=1:l)
               HistoryParticlehtable(i)%keys(j)       =HistoryParticlehtabletmp(i)%keys(j)
               HistoryParticlehtable(i)%borders_pos(j)=HistoryParticlehtabletmp(i)%borders_pos(j)
            END FORALL

            CALL HistoryParticlehtabletmp(i)%destroy(info)
            or_fail("HistoryParticlehtabletmp(i)%destroy")
         ENDDO
      ENDIF

      DEALLOCATE(HistoryParticlehtabletmp,STAT=info)
      or_fail_dealloc("HistoryParticlehtabletmp")

      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      END SUBROUTINE MCMCHistoryParticle_realloc
