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

#if   __DIME == __3D
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

      CHARACTER(LEN=ppm_char) :: caller='create_htable'

      CALL substart(caller,t0,info)

      IF (nelement.LE.0) THEN
         table%nrow = 1_ppm_kind_int64
      ELSE
         table%nrow = 2_ppm_kind_int64**(CEILING(LOG(REAL(nelement))/LOG(2.0)))
      ENDIF

      !---------------------------------------------------------------------
      !  Allocate array for hash table keys and array for positions on "borders" array.
      !  Set everything to NULL.
      !---------------------------------------------------------------------
      ALLOCATE(table%keys(table%nrow),SOURCE=htable_null_li,STAT=Info)
      or_fail_alloc('Failed to alocate htable_keys!',ppm_error=ppm_error_fatal)

      ALLOCATE(table%borders_pos(table%nrow),SOURCE=-one,STAT=Info)
      or_fail_alloc('Failed to alocate htable_borders_pos!',ppm_error=ppm_error_fatal)

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
      IF (ALLOCATED(table%keys)) THEN
         DEALLOCATE(table%keys,table%borders_pos,STAT=info)
         or_fail_dealloc('Failed to deallocate htable_keys & htable_borders_pos!',ppm_error=ppm_error_fatal)
      ENDIF

      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      END SUBROUTINE destroy_htable

      FUNCTION h_key(table,spot1,spot2,jump) RESULT(address)
      !!! Given the key and jump value, returns corresponding address on
      !!! "borders" array.
      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_htable),    INTENT(IN   ) :: table
      !!! The hashtable to create.

      INTEGER(ppm_kind_int64), INTENT(IN   ) :: spot1
      !!! First spot value to avoid double computing
      INTEGER(ppm_kind_int64), INTENT(IN   ) :: spot2
      !!! Second spot value to avoid double computing
      INTEGER(ppm_kind_int64), INTENT(IN   ) :: jump
      !!! Jump value for double hashing
      INTEGER(ppm_kind_int64)                :: address
      !!! Address that corresponds to given key
      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      ! we need to work with 64 bit integers because
      ! jump*h_func might overflow, and as there is no unsigned
      ! integer the resulting value might be negative and not
      ! a valid array index
      address=MOD(spot1+spot2*jump,table%nrow)+1_ppm_kind_int64
      RETURN
      END FUNCTION h_key

      SUBROUTINE hash_insert(table,key,value,info)
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
      INTEGER(ppm_kind_int64) :: jump
      INTEGER(ppm_kind_int64) :: spot
      INTEGER(ppm_kind_int64) :: fspot
      INTEGER(ppm_kind_int64) :: sspot

      info = 0
      jump = 0_ppm_kind_int64

      ! Get the address corresponding to given key
      fspot=HashFunc(key,seed1,table%nrow)
      spot=fspot+1_ppm_kind_int64

      ! Keep on searching withing bounds of hash table.
      DO WHILE (jump.LT.table%nrow)
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
         jump=jump + 1_ppm_kind_int64
         IF (jump.EQ.1_ppm_kind_int64) THEN
            sspot=HashFunc(key,seed2,table%nrow)
         ENDIF

         spot=table%h_key(fspot,sspot,jump)
      ENDDO

      ! If NOT returned within the while-loop, that means our hash table is
      ! not sufficiently large. Hence, regrowth will be needed!
      CALL table%grow(info,key,value)
      RETURN
      END SUBROUTINE hash_insert

      SUBROUTINE hash_insert_32(table,key_32,value,info)

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_htable)    :: table
      !!! The hashtable
      INTEGER,  INTENT(IN   ) :: key_32
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

      key=INT(key_32,KIND=ppm_kind_int64)
      CALL table%hash_insert(key,value,info)
      RETURN
      END SUBROUTINE hash_insert_32

      SUBROUTINE hash_insert_2D3D(table,key_2D3D,value,info)

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_htable)                 :: table
      !!! The hashtable
      INTEGER, DIMENSION(:), INTENT(IN   ) :: key_2D3D
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

      SELECT CASE (SIZE(key_2D3D))
      CASE (2)
         key=IndexHashFunctor64_2d(key_2D3D)
      CASE (3)
         key=IndexHashFunctor64_3d(key_2D3D)
      CASE DEFAULT
         info=ppm_error_fatal
         RETURN
      END SELECT
      CALL table%hash_insert(key,value,info)
      RETURN
      END SUBROUTINE hash_insert_2D3D

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
      INTEGER(ppm_kind_int64) :: jump
      INTEGER(ppm_kind_int64) :: spot
      INTEGER(ppm_kind_int64) :: fspot
      INTEGER(ppm_kind_int64) :: sspot

      LOGICAL :: KeyExist

      jump=0_ppm_kind_int64
      KeyExist=.TRUE.

      ! Get the other key that results in same
      ! hash key as for the inputkey.
      fspot=HashFunc(key,seed1,table%nrow)
      spot=fspot+1_ppm_kind_int64

      ! Keep on searching while we don't come across a NULL value or we don't
      ! exceed bounds of hash table.
      loop: DO WHILE(jump.LT.table%nrow)
          ! If key matches ...
          IF (table%keys(spot).EQ.key) THEN
             ! Set the return value and return
             value = table%borders_pos(spot)
             RETURN
          ELSE IF (table%keys(spot).EQ.htable_null_li) THEN
             IF (KeyExist) THEN
                KeyExist=.FALSE.
                IF (.NOT.ANY(key.EQ.table%keys)) EXIT loop
             ENDIF
          ENDIF
          ! Otherwise, keep on incrementing jump distance
          jump=jump+1_ppm_kind_int64
          IF (jump.EQ.1_ppm_kind_int64) THEN
             sspot=HashFunc(key,seed2,table%nrow)
          ENDIF

          spot=table%h_key(fspot,sspot,jump)
      ENDDO loop
      value = -one
      RETURN
      END FUNCTION hash_search

      FUNCTION hash_search_32(table,key_32) RESULT(value)
      IMPLICIT NONE
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_htable), INTENT(IN   ) :: table
      !!! The hashtable
      INTEGER,              INTENT(IN   ) :: key_32
      !!! Input key, which the corresponding value is asked for
      REAL(MK)                            :: value
      !!! Value corresponding to the input key

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: key

      key=INT(key_32,KIND=ppm_kind_int64)
      value=table%hash_search(key)
      RETURN
      END FUNCTION hash_search_32

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
      SUBROUTINE hash_remove(table,key,info,existed)
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
      INTEGER(ppm_kind_int64) :: jump
      INTEGER(ppm_kind_int64) :: spot
      INTEGER(ppm_kind_int64) :: fspot
      INTEGER(ppm_kind_int64) :: sspot

      info = 0
      IF (PRESENT(existed)) THEN
         IF (existed) THEN
            jump=0_ppm_kind_int64

            ! Get the address corresponding to given key
            fspot=HashFunc(key,seed1,table%nrow)
            spot=fspot+1_ppm_kind_int64

            ! Keep on searching withing bounds of hash table.
            DO WHILE (jump.LT.table%nrow)
               ! If an empty slot found ...
               IF (table%keys(spot).EQ.key) THEN
                  !Remove the key and the corresponding value and RETURN.
                  table%borders_pos(spot)=-one
                  table%keys(spot)=htable_null_li
                  RETURN
               ENDIF
               ! If the current slot is occupied, jump to next key that results
               ! in same hash function.
               jump=jump+1_ppm_kind_int64
               IF (jump.EQ.1_ppm_kind_int64) THEN
                  sspot=HashFunc(key,seed2,table%nrow)
               ENDIF

               spot=table%h_key(fspot,sspot,jump)
            ENDDO
            ! If NOT returned within the while-loop, that means the key was NOT found
            info = ppm_error_fatal
         ENDIf
         RETURN
      ELSE
         IF (ANY(key.EQ.table%keys)) THEN
            jump=0_ppm_kind_int64

            ! Get the address corresponding to given key
            fspot=HashFunc(key,seed1,table%nrow)
            spot=fspot+1_ppm_kind_int64

            ! Keep on searching withing bounds of hash table.
            DO WHILE (jump.LT.table%nrow)
               ! If an empty slot found ...
               IF (table%keys(spot).EQ.key) THEN
                  !Remove the key and the corresponding value and RETURN.
                  table%borders_pos(spot)=-one
                  table%keys(spot)=htable_null_li
                  RETURN
               ENDIF
               ! If the current slot is occupied, jump to next key that results
               ! in same hash function.
               jump=jump+1_ppm_kind_int64
               IF (jump.EQ.1_ppm_kind_int64) THEN
                  sspot=HashFunc(key,seed2,table%nrow)
               ENDIF

               spot=table%h_key(fspot,sspot,jump)
            ENDDO
            ! If NOT returned within the while-loop, that means the key was NOT found
            info = ppm_error_fatal
         ENDIF
      ENDIF
      RETURN
      END SUBROUTINE hash_remove

      SUBROUTINE hash_remove_32(table,key_32,info,existed)

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_htable)             :: table
      !!! The hashtable
      INTEGER,           INTENT(IN   ) :: key_32
      !!! Key to be removed
      INTEGER,           INTENT(  OUT) :: info
      !!! Info for whether removal was successful or not. 0 if SUCCESSFUL
      !!! and -1 otherwise.
      LOGICAL, OPTIONAL, INTENT(IN   ) :: existed

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: key

      key=INT(key_32,KIND=ppm_kind_int64)
      IF (PRESENT(existed)) THEN
         CALL table%hash_remove(key,info,existed)
      ELSE
         CALL table%hash_remove(key,info)
      ENDIF
      RETURN
      END SUBROUTINE hash_remove_32

      SUBROUTINE hash_remove_2D3D(table,key_2D3D, info,existed)

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_htable)                 :: table
      !!! The hashtable
      INTEGER, DIMENSION(:), INTENT(IN   ) :: key_2D3D
      !!! Key to be removed
      INTEGER,               INTENT(  OUT) :: info
      !!! Info for whether removal was successful or not. 0 if SUCCESSFUL
      !!! and -1 otherwise.
      LOGICAL, OPTIONAL,     INTENT(IN   ) :: existed

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: key

      SELECT CASE (SIZE(key_2D3D))
      CASE (2)
         key=IndexHashFunctor64_2d(key_2D3D)
      CASE (3)
         key=IndexHashFunctor64_3d(key_2D3D)
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
      END SUBROUTINE hash_remove_2D3D

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

      IF (table%nrow*2_ppm_kind_int64.GE.INT(ppm_big_i-1,ppm_kind_int64)) THEN
         !TOCHCECK
         fail("hashtable with more than 2^31-1 elements will fail",ppm_error=ppm_error_fatal)
      ENDIF

      nsize=INT(table%nrow)

      SELECT CASE (nsize)
      CASE (0)
         CALL table%create(1,info)
         or_fail("table%create")
      CASE DEFAULT
         ALLOCATE(keys_tmp(nsize),SOURCE=table%keys,STAT=info)
         or_fail_alloc("keys_tmp")

         ALLOCATE(borders_pos_tmp(nsize),SOURCE=table%borders_pos,STAT=info)
         or_fail_alloc("keys_tmp & borders_pos_tmp")

         NULLIFY(ranklist)
         CALL ppm_util_qsort(keys_tmp,ranklist,info,nsize)
         or_fail("ppm_util_qsort")

         CALL table%destroy(info)
         or_fail("table%destroy")

         CALL table%create(nsize*2,info)
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
      END SELECT

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

      nsize=INT(table%nrow)
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

      CHARACTER(LEN=ppm_char) :: caller='MCMCParticle_create_htable'

      CALL substart(caller,t0,info)

      IF (nelement.LE.0) THEN
         table%nrow = 1_ppm_kind_int64
      ELSE
         table%nrow = 2_ppm_kind_int64**(CEILING(LOG(REAL(nelement))/LOG(2.0)))
      ENDIF

      !---------------------------------------------------------------------
      !  Allocate array for hash table keys and array for positions on "borders" array.
      ! &
      !  Set everything to NULL.
      !---------------------------------------------------------------------
      ALLOCATE(table%keys(table%nrow),SOURCE=htable_null_li,STAT=Info)
      or_fail_alloc('Failed to alocate htable_keys!',ppm_error=ppm_error_fatal)

      ALLOCATE(table%borders_pos(table%nrow),SOURCE=MCMCParticle(-1,zero),STAT=Info)
      or_fail_alloc('Failed to alocate htable_borders_pos!',ppm_error=ppm_error_fatal)

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
      IF (ALLOCATED(table%keys)) THEN
         DEALLOCATE(table%keys,table%borders_pos,STAT=info)
         or_fail_dealloc('Failed to deallocate htable_keys & htable_borders_pos!',ppm_error=ppm_error_fatal)
      ENDIF

      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      END SUBROUTINE MCMCParticle_destroy_htable

      FUNCTION MCMCParticle_h_key(table,spot1,spot2,jump) RESULT(address)
      !!! Given the spots and jump value, returns corresponding address on
      !!! "borders" array.

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_MCMCParticlehtable), INTENT(IN   ) :: table
      !!! The hashtable to create. The pointer must not be NULL

      INTEGER(ppm_kind_int64),          INTENT(IN   ) :: spot1
      !!! First spot value to avoid double computing
      INTEGER(ppm_kind_int64),          INTENT(IN   ) :: spot2
      !!! Second spot value to avoid double computing
      INTEGER(ppm_kind_int64),          INTENT(IN   ) :: jump
      !!! Jump value for double hashing
      INTEGER(ppm_kind_int64)                         :: address
      !!! Address that corresponds to given key
      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      ! we need to work with 64 bit integers because
      ! jump*h_func might overflow, and as there is no unsigned
      ! integer the resulting value might be negative and not
      ! a valid array index
      address=MOD(spot1+spot2*jump,table%nrow)+1_ppm_kind_int64
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
      INTEGER(ppm_kind_int64) :: jump
      INTEGER(ppm_kind_int64) :: spot
      INTEGER(ppm_kind_int64) :: fspot
      INTEGER(ppm_kind_int64) :: sspot

      info = 0
      jump = 0_ppm_kind_int64

      ! Get the address corresponding to given key
      fspot=HashFunc(key,seed1,table%nrow)
      spot=fspot+1_ppm_kind_int64

      ! Keep on searching withing bounds of hash table.
      DO WHILE (jump.LT.table%nrow)
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
         jump=jump + 1_ppm_kind_int64
         IF (jump.EQ.1_ppm_kind_int64) THEN
            sspot=HashFunc(key,seed2,table%nrow)
         ENDIF

         spot=table%h_key(fspot,sspot,jump)
      ENDDO

      ! If NOT returned within the while-loop, that means our hash table is
      ! not sufficiently large. Hence, regrowth will be needed!
      CALL table%grow(info,key,value)
      RETURN
      END SUBROUTINE MCMCParticle_hash_insert

      SUBROUTINE MCMCParticle_hash_insert_2D3D(table,key_2D3D,value,info)

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_MCMCParticlehtable)     :: table
      !!! The hashtable
      INTEGER, DIMENSION(:), INTENT(IN   ) :: key_2D3D
      !!! Key to be stored
      TYPE(MCMCParticle),    INTENT(IN   ) :: value
      !!! Value that corresponds to given key
      INTEGER,               INTENT(  OUT) :: info
      !!! Info for whether insertion was successful or not. 0 if SUCCESSFUL
      !!! and -1 otherwise.

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: jump
      INTEGER(ppm_kind_int64) :: spot
      INTEGER(ppm_kind_int64) :: fspot
      INTEGER(ppm_kind_int64) :: sspot
      INTEGER(ppm_kind_int64) :: key
      INTEGER                 :: ssize

      ssize=SIZE(key_2D3D)
      SELECT CASE (ssize)
      CASE (2)
         info=0

         SELECT CASE (value%candlabel)
         CASE (0)
            ! Get the address corresponding to given key
            fspot=HashFunc_XY(key_2D3D(1),key_2D3D(2),seed1,table%nrow,key)
         CASE DEFAULT
            ! Get the address corresponding to given key
            fspot=HashFunc_XYLabel(key_2D3D(1),key_2D3D(2),value%candlabel,seed1,table%nrow,key)
         END SELECT
      CASE (3)
         info=0

         SELECT CASE (value%candlabel)
         CASE (0)
            ! Get the address corresponding to given key
            fspot=HashFunc_XYZ(key_2D3D(1),key_2D3D(2),key_2D3D(3),seed1,table%nrow,key)
         CASE DEFAULT
            ! Get the address corresponding to given key
            fspot=HashFunc_XYZLabel(key_2D3D(1),key_2D3D(2),key_2D3D(3),value%candlabel,seed1,table%nrow,key)
         END SELECT
      CASE DEFAULT
         info=ppm_error_fatal
         RETURN
      END SELECT

      jump = 0_ppm_kind_int64

      spot=fspot+1_ppm_kind_int64

      ! Keep on searching withing bounds of hash table.
      DO WHILE (jump.LT.table%nrow)
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
         jump=jump + 1_ppm_kind_int64
         IF (jump.EQ.1_ppm_kind_int64) THEN
            sspot=HashFunc(key,seed2,table%nrow)
         ENDIF

         spot=table%h_key(fspot,sspot,jump)
      ENDDO

      ! If NOT returned within the while-loop, that means our hash table is
      ! not sufficiently large. Hence, regrowth will be needed!
      CALL table%grow(info,key,value)
      RETURN
      END SUBROUTINE MCMCParticle_hash_insert_2D3D

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
      INTEGER(ppm_kind_int64) :: jump
      INTEGER(ppm_kind_int64) :: spot
      INTEGER(ppm_kind_int64) :: fspot
      INTEGER(ppm_kind_int64) :: sspot
      INTEGER(ppm_kind_int64) :: key

      info = 0
      jump = 0_ppm_kind_int64

      SELECT CASE (value%candlabel)
      CASE (0)
         ! Get the address corresponding to given key
         fspot=HashFunc_XY(key_1,key_2,seed1,table%nrow,key)
      CASE DEFAULT
         ! Get the address corresponding to given key
         fspot=HashFunc_XYLabel(key_1,key_2,value%candlabel,seed1,table%nrow,key)
      END SELECT

      spot=fspot+1_ppm_kind_int64

      ! Keep on searching withing bounds of hash table.
      DO WHILE (jump.LT.table%nrow)
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
         jump=jump + 1_ppm_kind_int64
         IF (jump.EQ.1_ppm_kind_int64) THEN
            sspot=HashFunc(key,seed2,table%nrow)
         ENDIF

         spot=table%h_key(fspot,sspot,jump)
      ENDDO

      ! If NOT returned within the while-loop, that means our hash table is
      ! not sufficiently large. Hence, regrowth will be needed!
      CALL table%grow(info,key,value)
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
      INTEGER(ppm_kind_int64) :: jump
      INTEGER(ppm_kind_int64) :: spot
      INTEGER(ppm_kind_int64) :: fspot
      INTEGER(ppm_kind_int64) :: sspot
      INTEGER(ppm_kind_int64) :: key

      info = 0
      jump = 0_ppm_kind_int64

      SELECT CASE (value%candlabel)
      CASE (0)
         ! Get the address corresponding to given key
         fspot=HashFunc_XYZ(key_1,key_2,key_3,seed1,table%nrow,key)
      CASE DEFAULT
         ! Get the address corresponding to given key
         fspot=HashFunc_XYZLabel(key_1,key_2,key_3,value%candlabel,seed1,table%nrow,key)
      END SELECT

      spot=fspot+1_ppm_kind_int64

      ! Keep on searching withing bounds of hash table.
      DO WHILE (jump.LT.table%nrow)
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
         jump=jump + 1_ppm_kind_int64
         IF (jump.EQ.1_ppm_kind_int64) THEN
            sspot=HashFunc(key,seed2,table%nrow)
         ENDIF

         spot=table%h_key(fspot,sspot,jump)
      ENDDO

      ! If NOT returned within the while-loop, that means our hash table is
      ! not sufficiently large. Hence, regrowth will be needed!
      CALL table%grow(info,key,value)
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
      INTEGER(ppm_kind_int64) :: jump
      INTEGER(ppm_kind_int64) :: spot
      INTEGER(ppm_kind_int64) :: fspot
      INTEGER(ppm_kind_int64) :: sspot

      LOGICAL :: KeyExist

      jump=0_ppm_kind_int64

      KeyExist=.TRUE.

      ! Get the other key that results in same
      ! hash key as for the inputkey.
      fspot=HashFunc(key,seed1,table%nrow)
      spot=fspot+1_ppm_kind_int64

      ! Keep on searching while we don't come across a NULL value or we don't
      ! exceed bounds of hash table.
      loop: DO WHILE(jump.LT.table%nrow)
          ! If key matches ...
          IF (table%keys(spot).EQ.key) THEN
             ! Set the return value and return
             value = table%borders_pos(spot)
             RETURN
          ELSE IF (table%keys(spot).EQ.htable_null_li) THEN
             IF (KeyExist) THEN
                KeyExist=.FALSE.
                IF (.NOT.ANY(key.EQ.table%keys)) EXIT loop
             ENDIF
          ENDIF
          ! Otherwise, keep on incrementing jump distance
          jump=jump+1_ppm_kind_int64
          IF (jump.EQ.1_ppm_kind_int64) THEN
             sspot=HashFunc(key,seed2,table%nrow)
          ENDIF

          spot=table%h_key(fspot,sspot,jump)
      ENDDO loop
      value = MCMCParticle(-1,zero)
      RETURN
      END FUNCTION MCMCParticle_hash_search

      FUNCTION MCMCParticle_hash_search_2d(table,key_1,key_2,candlabel) RESULT(value)

      IMPLICIT NONE
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_MCMCParticlehtable), INTENT(IN   ) :: table
      !!! The hashtable

      INTEGER,                          INTENT(IN   ) :: key_1
      INTEGER,                          INTENT(IN   ) :: key_2
      !!! Input key, which the corresponding value is asked for
      INTEGER,                          INTENT(IN   ) :: candlabel
      !!! Input label

      TYPE(MCMCParticle)                              :: value
      !!! Value corresponding to the input key
      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: jump
      INTEGER(ppm_kind_int64) :: spot
      INTEGER(ppm_kind_int64) :: fspot
      INTEGER(ppm_kind_int64) :: sspot
      INTEGER(ppm_kind_int64) :: key

      LOGICAL :: KeyExist

      jump=0_ppm_kind_int64

      KeyExist=.TRUE.

      SELECT CASE (candlabel)
      CASE (0)
         ! Get the address corresponding to given key
         fspot=HashFunc_XY(key_1,key_2,seed1,table%nrow,key)
      CASE DEFAULT
         ! Get the address corresponding to given key
         fspot=HashFunc_XYLabel(key_1,key_2,candlabel,seed1,table%nrow,key)
      END SELECT

      spot=fspot+1_ppm_kind_int64

      ! Keep on searching while we don't come across a NULL value or we don't
      ! exceed bounds of hash table.
      loop: DO WHILE(jump.LT.table%nrow)
          ! If key matches ...
          IF (table%keys(spot).EQ.key) THEN
             ! Set the return value and return
             value = table%borders_pos(spot)
             RETURN
          ELSE IF (table%keys(spot).EQ.htable_null_li) THEN
             IF (KeyExist) THEN
                KeyExist=.FALSE.
                IF (.NOT.ANY(key.EQ.table%keys)) EXIT loop
             ENDIF
          ENDIF
          ! Otherwise, keep on incrementing jump distance
          jump=jump+1_ppm_kind_int64
          IF (jump.EQ.1_ppm_kind_int64) THEN
             sspot=HashFunc(key,seed2,table%nrow)
          ENDIF

          spot=table%h_key(fspot,sspot,jump)
      ENDDO loop
      value = MCMCParticle(-1,zero)
      RETURN
      END FUNCTION MCMCParticle_hash_search_2d

      FUNCTION MCMCParticle_hash_search_3d(table,key_1,key_2,key_3,candlabel) RESULT(value)

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
      INTEGER,                          INTENT(IN   ) :: candlabel
      !!! Input label

      TYPE(MCMCParticle)                              :: value
      !!! Value corresponding to the input key
      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: jump
      INTEGER(ppm_kind_int64) :: spot
      INTEGER(ppm_kind_int64) :: fspot
      INTEGER(ppm_kind_int64) :: sspot
      INTEGER(ppm_kind_int64) :: key

      LOGICAL :: KeyExist

      jump=0_ppm_kind_int64

      KeyExist=.TRUE.

      SELECT CASE (candlabel)
      CASE (0)
         ! Get the address corresponding to given key
         fspot=HashFunc_XYZ(key_1,key_2,key_3,seed1,table%nrow,key)
      CASE DEFAULT
         ! Get the address corresponding to given key
         fspot=HashFunc_XYZLabel(key_1,key_2,key_3,candlabel,seed1,table%nrow,key)
      END SELECT

      spot=fspot+1_ppm_kind_int64

      ! Keep on searching while we don't come across a NULL value or we don't
      ! exceed bounds of hash table.
      loop: DO WHILE(jump.LT.table%nrow)
          ! If key matches ...
          IF (table%keys(spot).EQ.key) THEN
             ! Set the return value and return
             value = table%borders_pos(spot)
             RETURN
          ELSE IF (table%keys(spot).EQ.htable_null_li) THEN
             IF (KeyExist) THEN
                KeyExist=.FALSE.
                IF (.NOT.ANY(key.EQ.table%keys)) EXIT loop
             ENDIF
          ENDIF
          ! Otherwise, keep on incrementing jump distance
          jump=jump+1_ppm_kind_int64
          IF (jump.EQ.1_ppm_kind_int64) THEN
             sspot=HashFunc(key,seed2,table%nrow)
          ENDIF

          spot=table%h_key(fspot,sspot,jump)
      ENDDO loop
      value = MCMCParticle(-1,zero)
      RETURN
      END FUNCTION MCMCParticle_hash_search_3d

      LOGICAL FUNCTION MCMCParticle_hash_contains(table,key,value)
      !!! Given the key, searchs the key on the hash table and returns the
      !!! true if the corresponding value is the same as an input value
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
      INTEGER(ppm_kind_int64) :: jump
      INTEGER(ppm_kind_int64) :: spot
      INTEGER(ppm_kind_int64) :: fspot
      INTEGER(ppm_kind_int64) :: sspot

      LOGICAL :: KeyExist

      IF (table%nrow.GT.0_ppm_kind_int64) THEN
         jump=0_ppm_kind_int64

         KeyExist=.TRUE.

         ! Get the other key that results in same
         ! hash key as for the inputkey.
         fspot=HashFunc(key,seed1,table%nrow)
         spot=fspot+1_ppm_kind_int64

         ! Keep on searching while we don't come across a NULL value or we don't
         ! exceed bounds of hash table.
         loop: DO WHILE(jump.LT.table%nrow)
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
               IF (KeyExist) THEN
                  KeyExist=.FALSE.
                  IF (.NOT.ANY(key.EQ.table%keys)) EXIT loop
               ENDIF
            ENDIF
            ! Otherwise, keep on incrementing jump distance
            jump=jump+1_ppm_kind_int64
            IF (jump.EQ.1_ppm_kind_int64) THEN
               sspot=HashFunc(key,seed2,table%nrow)
            ENDIF

            spot=table%h_key(fspot,sspot,jump)
         ENDDO loop
      ENDIF ! (table%nrow.GT.0_ppm_kind_int64)
      MCMCParticle_hash_contains=.FALSE.
      RETURN
      END FUNCTION MCMCParticle_hash_contains

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
      INTEGER(ppm_kind_int64) :: jump
      INTEGER(ppm_kind_int64) :: spot
      INTEGER(ppm_kind_int64) :: fspot
      INTEGER(ppm_kind_int64) :: sspot
      INTEGER(ppm_kind_int64) :: key

      LOGICAL :: KeyExist

      IF (table%nrow.GT.0_ppm_kind_int64) THEN
         jump=0_ppm_kind_int64

         KeyExist=.TRUE.

         SELECT CASE (value%candlabel)
         CASE (0)
            ! Get the address corresponding to given key
            fspot=HashFunc_XY(key_1,key_2,seed1,table%nrow,key)
         CASE DEFAULT
            ! Get the address corresponding to given key
            fspot=HashFunc_XYLabel(key_1,key_2,value%candlabel,seed1,table%nrow,key)
         END SELECT

         spot=fspot+1_ppm_kind_int64

         ! Keep on searching while we don't come across a NULL value or we don't
         ! exceed bounds of hash table.
         loop: DO WHILE(jump.LT.table%nrow)
            ! If key matches ...
            IF (table%keys(spot).EQ.key) THEN
               ! Set the return value and return
               ASSOCIATE (tmpParticle => table%borders_pos(spot))
                  IF (value%candlabel.NE.tmpParticle%candlabel) THEN
                     MCMCParticle_hash_contains_2d=.FALSE.
                     RETURN
                  ENDIF
               END ASSOCIATE
               MCMCParticle_hash_contains_2d=.TRUE.
               RETURN
            ELSE IF (table%keys(spot).EQ.htable_null_li) THEN
               IF (KeyExist) THEN
                  KeyExist=.FALSE.
                  IF (.NOT.ANY(key.EQ.table%keys)) EXIT loop
               ENDIF
            ENDIF
            ! Otherwise, keep on incrementing jump distance
            jump=jump+1_ppm_kind_int64
            IF (jump.EQ.1_ppm_kind_int64) THEN
               sspot=HashFunc(key,seed2,table%nrow)
            ENDIF

            spot=table%h_key(fspot,sspot,jump)
         ENDDO loop
      ENDIF ! (table%nrow.GT.0_ppm_kind_int64)
      MCMCParticle_hash_contains_2d=.FALSE.
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
      INTEGER(ppm_kind_int64) :: jump
      INTEGER(ppm_kind_int64) :: spot
      INTEGER(ppm_kind_int64) :: fspot
      INTEGER(ppm_kind_int64) :: sspot
      INTEGER(ppm_kind_int64) :: key

      LOGICAL :: KeyExist

      IF (table%nrow.GT.0_ppm_kind_int64) THEN
         jump=0_ppm_kind_int64

         KeyExist=.TRUE.

         SELECT CASE (value%candlabel)
         CASE (0)
            ! Get the address corresponding to given key
            fspot=HashFunc_XYZ(key_1,key_2,key_3,seed1,table%nrow,key)
         CASE DEFAULT
            ! Get the address corresponding to given key
            fspot=HashFunc_XYZLabel(key_1,key_2,key_3,value%candlabel,seed1,table%nrow,key)
         END SELECT

         spot=fspot+1_ppm_kind_int64

         ! Keep on searching while we don't come across a NULL value or we don't
         ! exceed bounds of hash table.
         loop: DO WHILE(jump.LT.table%nrow)
            ! If key matches ...
            IF (table%keys(spot).EQ.key) THEN
               ! Set the return value and return
               ASSOCIATE (tmpParticle => table%borders_pos(spot))
                  IF (value%candlabel.NE.tmpParticle%candlabel) THEN
                     MCMCParticle_hash_contains_3d=.FALSE.
                     RETURN
                  ENDIF
               END ASSOCIATE
               MCMCParticle_hash_contains_3d=.TRUE.
               RETURN
            ELSE IF (table%keys(spot).EQ.htable_null_li) THEN
               IF (KeyExist) THEN
                  KeyExist=.FALSE.
                  IF (.NOT.ANY(key.EQ.table%keys)) EXIT loop
               ENDIF
            ENDIF
            ! Otherwise, keep on incrementing jump distance
            jump=jump+1_ppm_kind_int64
            IF (jump.EQ.1_ppm_kind_int64) THEN
               sspot=HashFunc(key,seed2,table%nrow)
            ENDIF

            spot=table%h_key(fspot,sspot,jump)
         ENDDO loop
      ENDIF ! (table%nrow.GT.0_ppm_kind_int64)
      MCMCParticle_hash_contains_3d=.FALSE.
      RETURN
      END FUNCTION MCMCParticle_hash_contains_3d

      !TODO check this sub
      !Yaser
      !This function PROBABLY suffers from a bug!
      SUBROUTINE MCMCParticle_hash_remove(table,key,info,existed)
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
      INTEGER(ppm_kind_int64) :: jump
      INTEGER(ppm_kind_int64) :: spot
      INTEGER(ppm_kind_int64) :: fspot
      INTEGER(ppm_kind_int64) :: sspot

      info = 0

      IF (PRESENT(existed)) THEN
         IF (existed) THEN
            jump=0_ppm_kind_int64

            ! Get the other key that results in same
            ! hash key as for the inputkey.
            fspot=HashFunc(key,seed1,table%nrow)
            spot=fspot+1_ppm_kind_int64

            ! Keep on searching withing bounds of hash table.
            DO WHILE (jump.LT.table%nrow)
               ! If an empty slot found ...
               IF (table%keys(spot).EQ.key) THEN
                  !Remove the key and the corresponding value and RETURN.
                  table%borders_pos(spot)=MCMCParticle(-1,zero)
                  table%keys(spot)=htable_null_li
                  RETURN
               ENDIF
               ! If the current slot is occupied, jump to next key that results
               ! in same hash function.
               jump=jump+1_ppm_kind_int64
               IF (jump.EQ.1_ppm_kind_int64) THEN
                  sspot=HashFunc(key,seed2,table%nrow)
               ENDIF

               spot=table%h_key(fspot,sspot,jump)
            ENDDO

            ! If NOT returned within the while-loop, that means the key was NOT found
            info = ppm_error_fatal
         ENDIF
         RETURN
      ELSE
         IF (ANY(key.EQ.table%keys)) THEN
            jump=0_ppm_kind_int64

            ! Get the other key that results in same
            ! hash key as for the inputkey.
            fspot=HashFunc(key,seed1,table%nrow)
            spot=fspot+1_ppm_kind_int64

            ! Keep on searching withing bounds of hash table.
            DO WHILE (jump.LT.table%nrow)
               ! If an empty slot found ...
               IF (table%keys(spot).EQ.key) THEN
                  !Remove the key and the corresponding value and RETURN.
                  table%borders_pos(spot)=MCMCParticle(-1,zero)
                  table%keys(spot)=htable_null_li
                  RETURN
               ENDIF
               ! If the current slot is occupied, jump to next key that results
               ! in same hash function.
               jump=jump+1_ppm_kind_int64
               IF (jump.EQ.1_ppm_kind_int64) THEN
                  sspot=HashFunc(key,seed2,table%nrow)
               ENDIF

               spot=table%h_key(fspot,sspot,jump)
            ENDDO
            ! If NOT returned within the while-loop, that means the key was NOT found
            info = ppm_error_fatal
         ENDIF
      ENDIF
      RETURN
      END SUBROUTINE MCMCParticle_hash_remove

      SUBROUTINE MCMCParticle_hash_remove_2D3D(table,key_2D3D,candlabel,info,existed)

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_MCMCParticlehtable)     :: table
      !!! The hashtable
      INTEGER, DIMENSION(:), INTENT(IN   ) :: key_2D3D
      !!! Key to be removed
      INTEGER,               INTENT(IN   ) :: candlabel
      !!! Input label
      INTEGER,               INTENT(  OUT) :: info
      !!! Info for whether removal was successful or not. 0 if SUCCESSFUL
      !!! and -1 otherwise.
      LOGICAL, OPTIONAL,     INTENT(IN   ) :: existed

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: jump
      INTEGER(ppm_kind_int64) :: spot
      INTEGER(ppm_kind_int64) :: fspot
      INTEGER(ppm_kind_int64) :: sspot
      INTEGER(ppm_kind_int64) :: key
      INTEGER                 :: ssize

      info = 0

      IF (PRESENT(existed)) THEN
         IF (existed) THEN
            jump=0_ppm_kind_int64

            ssize=SIZE(key_2D3D)
            SELECT CASE (ssize)
            CASE (2)
               SELECT CASE (candlabel)
               CASE (0)
                  ! Get the address corresponding to given key
                  fspot=HashFunc_XY(key_2D3D(1),key_2D3D(2),seed1,table%nrow,key)
               CASE DEFAULT
                  ! Get the address corresponding to given key
                  fspot=HashFunc_XYLabel(key_2D3D(1),key_2D3D(2),candlabel,seed1,table%nrow,key)
               END SELECT
            CASE (3)
               SELECT CASE (candlabel)
               CASE (0)
                  ! Get the address corresponding to given key
                  fspot=HashFunc_XYZ(key_2D3D(1),key_2D3D(2),key_2D3D(3),seed1,table%nrow,key)
               CASE DEFAULT
                  ! Get the address corresponding to given key
                  fspot=HashFunc_XYZLabel(key_2D3D(1),key_2D3D(2),key_2D3D(3),candlabel,seed1,table%nrow,key)
               END SELECT
            CASE DEFAULT
               info=ppm_error_fatal
               RETURN
            END SELECT

            spot=fspot+1_ppm_kind_int64

            ! Keep on searching withing bounds of hash table.
            DO WHILE (jump.LT.table%nrow)
               ! If an empty slot found ...
               IF (table%keys(spot).EQ.key) THEN
                  !Remove the key and the corresponding value and RETURN.
                  table%borders_pos(spot)=MCMCParticle(-1,zero)
                  table%keys(spot)=htable_null_li
                  RETURN
               ENDIF
               ! If the current slot is occupied, jump to next key that results
               ! in same hash function.
               jump=jump+1_ppm_kind_int64
               IF (jump.EQ.1_ppm_kind_int64) THEN
                  sspot=HashFunc(key,seed2,table%nrow)
               ENDIF

               spot=table%h_key(fspot,sspot,jump)
            ENDDO

            ! If NOT returned within the while-loop, that means the key was NOT found
            info = ppm_error_fatal
         ENDIF
         RETURN
      ELSE
         ssize=SIZE(key_2D3D)
         SELECT CASE (ssize)
         CASE (2)
            SELECT CASE (candlabel)
            CASE (0)
               key=IndexHashFunctor_label_2d(key_2D3D(1),key_2D3D(2))
            CASE DEFAULT
               key=IndexHashFunctor_label_2d(key_2D3D(1),key_2D3D(2),candlabel)
            END SELECT
         CASE (3)
            SELECT CASE (candlabel)
            CASE (0)
               key=IndexHashFunctor_label_3d(key_2D3D(1),key_2D3D(2),key_2D3D(3))
            CASE DEFAULT
               key=IndexHashFunctor_label_3d(key_2D3D(1),key_2D3D(2),key_2D3D(3),candlabel)
            END SELECT
         CASE DEFAULT
            info=ppm_error_fatal
            RETURN
         END SELECT

         IF (ANY(key.EQ.table%keys)) THEN
            jump=0_ppm_kind_int64

            ! Get the other key that results in same
            ! hash key as for the inputkey.
            fspot=HashFunc(key,seed1,table%nrow)
            spot=fspot+1_ppm_kind_int64

            ! Keep on searching withing bounds of hash table.
            DO WHILE (jump.LT.table%nrow)
               ! If an empty slot found ...
               IF (table%keys(spot).EQ.key) THEN
                  !Remove the key and the corresponding value and RETURN.
                  table%borders_pos(spot)=MCMCParticle(-1,zero)
                  table%keys(spot)=htable_null_li
                  RETURN
               ENDIF
               ! If the current slot is occupied, jump to next key that results
               ! in same hash function.
               jump=jump+1_ppm_kind_int64
               IF (jump.EQ.1_ppm_kind_int64) THEN
                  sspot=HashFunc(key,seed2,table%nrow)
               ENDIF

               spot=table%h_key(fspot,sspot,jump)
            ENDDO
            ! If NOT returned within the while-loop, that means the key was NOT found
            info = ppm_error_fatal
         ENDIF
      ENDIF
      RETURN
      END SUBROUTINE MCMCParticle_hash_remove_2D3D

      SUBROUTINE MCMCParticle_hash_remove_2d(table,key_1,key_2,candlabel,info,existed)

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_MCMCParticlehtable) :: table
      !!! The hashtable
      INTEGER,           INTENT(IN   ) :: key_1
      INTEGER,           INTENT(IN   ) :: key_2
      !!! Key to be removed
      INTEGER,           INTENT(IN   ) :: candlabel
      !!! Input label
      INTEGER,           INTENT(  OUT) :: info
      !!! Info for whether removal was successful or not. 0 if SUCCESSFUL
      !!! and -1 otherwise.
      LOGICAL, OPTIONAL, INTENT(IN   ) :: existed

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: jump
      INTEGER(ppm_kind_int64) :: spot
      INTEGER(ppm_kind_int64) :: fspot
      INTEGER(ppm_kind_int64) :: sspot
      INTEGER(ppm_kind_int64) :: key

      info = 0

      IF (PRESENT(existed)) THEN
         IF (existed) THEN
            jump=0_ppm_kind_int64

            SELECT CASE (candlabel)
            CASE (0)
               ! Get the address corresponding to given key
               fspot=HashFunc_XY(key_1,key_2,seed1,table%nrow,key)
            CASE DEFAULT
               ! Get the address corresponding to given key
               fspot=HashFunc_XYLabel(key_1,key_2,candlabel,seed1,table%nrow,key)
            END SELECT

            spot=fspot+1_ppm_kind_int64

            ! Keep on searching withing bounds of hash table.
            DO WHILE (jump.LT.table%nrow)
               ! If an empty slot found ...
               IF (table%keys(spot).EQ.key) THEN
                  !Remove the key and the corresponding value and RETURN.
                  table%borders_pos(spot)=MCMCParticle(-1,zero)
                  table%keys(spot)=htable_null_li
                  RETURN
               ENDIF
               ! If the current slot is occupied, jump to next key that results
               ! in same hash function.
               jump=jump+1_ppm_kind_int64
               IF (jump.EQ.1_ppm_kind_int64) THEN
                  sspot=HashFunc(key,seed2,table%nrow)
               ENDIF

               spot=table%h_key(fspot,sspot,jump)
            ENDDO

            ! If NOT returned within the while-loop, that means the key was NOT found
            info = ppm_error_fatal
         ENDIF
         RETURN
      ELSE
         SELECT CASE (candlabel)
         CASE (0)
            key=IndexHashFunctor_label_2d(key_1,key_2)
         CASE DEFAULT
            key=IndexHashFunctor_label_2d(key_1,key_2,candlabel)
         END SELECT

         IF (ANY(key.EQ.table%keys)) THEN
            jump=0_ppm_kind_int64

            ! Get the other key that results in same
            ! hash key as for the inputkey.
            fspot=HashFunc(key,seed1,table%nrow)
            spot=fspot+1_ppm_kind_int64

            ! Keep on searching withing bounds of hash table.
            DO WHILE (jump.LT.table%nrow)
               ! If an empty slot found ...
               IF (table%keys(spot).EQ.key) THEN
                  !Remove the key and the corresponding value and RETURN.
                  table%borders_pos(spot)=MCMCParticle(-1,zero)
                  table%keys(spot)=htable_null_li
                  RETURN
               ENDIF
               ! If the current slot is occupied, jump to next key that results
               ! in same hash function.
               jump=jump+1_ppm_kind_int64
               IF (jump.EQ.1_ppm_kind_int64) THEN
                  sspot=HashFunc(key,seed2,table%nrow)
               ENDIF

               spot=table%h_key(fspot,sspot,jump)
            ENDDO
            ! If NOT returned within the while-loop, that means the key was NOT found
            info = ppm_error_fatal
         ENDIF
      ENDIF
      RETURN
      END SUBROUTINE MCMCParticle_hash_remove_2d

      SUBROUTINE MCMCParticle_hash_remove_3d(table,key_1,key_2,key_3,candlabel,info,existed)

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
      INTEGER,           INTENT(IN   ) :: candlabel
      !!! Input label
      INTEGER,           INTENT(  OUT) :: info
      !!! Info for whether removal was successful or not. 0 if SUCCESSFUL
      !!! and -1 otherwise.
      LOGICAL, OPTIONAL, INTENT(IN   ) :: existed

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: jump
      INTEGER(ppm_kind_int64) :: spot
      INTEGER(ppm_kind_int64) :: fspot
      INTEGER(ppm_kind_int64) :: sspot
      INTEGER(ppm_kind_int64) :: key

      info = 0

      IF (PRESENT(existed)) THEN
         IF (existed) THEN
            jump=0_ppm_kind_int64

            SELECT CASE (candlabel)
            CASE (0)
               ! Get the address corresponding to given key
               fspot=HashFunc_XYZ(key_1,key_2,key_3,seed1,table%nrow,key)
            CASE DEFAULT
               ! Get the address corresponding to given key
               fspot=HashFunc_XYZLabel(key_1,key_2,key_3,candlabel,seed1,table%nrow,key)
            END SELECT

            spot=fspot+1_ppm_kind_int64

            ! Keep on searching withing bounds of hash table.
            DO WHILE (jump.LT.table%nrow)
               ! If an empty slot found ...
               IF (table%keys(spot).EQ.key) THEN
                  !Remove the key and the corresponding value and RETURN.
                  table%borders_pos(spot)=MCMCParticle(-1,zero)
                  table%keys(spot)=htable_null_li
                  RETURN
               ENDIF
               ! If the current slot is occupied, jump to next key that results
               ! in same hash function.
               jump=jump+1_ppm_kind_int64
               IF (jump.EQ.1_ppm_kind_int64) THEN
                  sspot=HashFunc(key,seed2,table%nrow)
               ENDIF

               spot=table%h_key(fspot,sspot,jump)
            ENDDO

            ! If NOT returned within the while-loop, that means the key was NOT found
            info = ppm_error_fatal
         ENDIF
         RETURN
      ELSE
         SELECT CASE (candlabel)
         CASE (0)
            key=IndexHashFunctor_label_3d(key_1,key_2,key_3)
         CASE DEFAULT
            key=IndexHashFunctor_label_3d(key_1,key_2,key_3,candlabel)
         END SELECT

         IF (ANY(key.EQ.table%keys)) THEN
            jump=0_ppm_kind_int64

            ! Get the other key that results in same
            ! hash key as for the inputkey.
            fspot=HashFunc(key,seed1,table%nrow)
            spot=fspot+1_ppm_kind_int64

            ! Keep on searching withing bounds of hash table.
            DO WHILE (jump.LT.table%nrow)
               ! If an empty slot found ...
               IF (table%keys(spot).EQ.key) THEN
                  !Remove the key and the corresponding value and RETURN.
                  table%borders_pos(spot)=MCMCParticle(-1,zero)
                  table%keys(spot)=htable_null_li
                  RETURN
               ENDIF
               ! If the current slot is occupied, jump to next key that results
               ! in same hash function.
               jump=jump+1_ppm_kind_int64
               IF (jump.EQ.1_ppm_kind_int64) THEN
                  sspot=HashFunc(key,seed2,table%nrow)
               ENDIF

               spot=table%h_key(fspot,sspot,jump)
            ENDDO
            ! If NOT returned within the while-loop, that means the key was NOT found
            info = ppm_error_fatal
         ENDIF
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

      IF (table%nrow*2_ppm_kind_int64.GE.INT(ppm_big_i-1,ppm_kind_int64)) THEN
         !TOCHCECK
         fail("hashtable with more than 2^31-1 elements will fail",ppm_error=ppm_error_fatal)
      ENDIF

      nsize=INT(table%nrow)

      SELECT CASE (nsize)
      CASE (0)
         CALL table%create(1,info)
         or_fail("table%create")
      CASE DEFAULT
         ALLOCATE(keys_tmp(nsize),SOURCE=table%keys,STAT=info)
         or_fail_alloc("keys_tmp")

         ALLOCATE(borders_pos_tmp(nsize),SOURCE=table%borders_pos,STAT=info)
         or_fail_alloc("keys_tmp & borders_pos_tmp")

         NULLIFY(ranklist)
         CALL ppm_util_qsort(keys_tmp,ranklist,info,nsize)
         or_fail("ppm_util_qsort")

         CALL table%destroy(info)
         or_fail("table%destroy")

         CALL table%create(nsize*2,info)
         or_fail("table%create")

         DO i=1,nsize
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
      END SELECT

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

      nsize=INT(table%nrow)
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
      !!! The hashtable.
      REAL(MK), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: proposal

      INTEGER,                             INTENT(IN   ) :: proposalSize
      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      REAL(ppm_kind_double) :: t0

      INTEGER(ppm_kind_int64) :: i
      INTEGER                 :: nsize
      INTEGER                 :: info

      CHARACTER(LEN=*), PARAMETER :: caller="MCMCParticle_hash_packproposal"

      !-------------------------------------------------------------------------
      !  Initialize
      !-------------------------------------------------------------------------
      CALL substart(caller,t0,info)

      IF (table%nrow.GT.0_ppm_kind_int64) THEN
         IF (ALLOCATED(proposal)) THEN
            IF (proposalSize.NE.SIZE(proposal)) THEN
               DEALLOCATE(proposal,STAT=info)
               or_fail_dealloc("proposal",ppm_error=ppm_error_fatal)
               ALLOCATE(proposal(proposalSize),STAT=info)
               or_fail_alloc("proposal",ppm_error=ppm_error_fatal)
            ENDIF
         ELSE
            ALLOCATE(proposal(proposalSize),STAT=info)
            or_fail_alloc("proposal",ppm_error=ppm_error_fatal)
         ENDIF
         nsize=0
         DO i=1_ppm_kind_int64,table%nrow
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

      FUNCTION MCMCParticle_hash_packproposal_(table,proposal,startindex,endindex) RESULT(info)
      !---------------------------------------------------------------------
      !  Note :
      !  This subroutine does not check the size of the proposal !
      !  It is user responsibility to make sure that the allocated proposal
      !  Has been allocated from startindex:endindex
      !---------------------------------------------------------------------
      USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_fatal
      USE ppm_module_substart, ONLY : substart
      USE ppm_module_substop, ONLY : substop
      USE ppm_module_error, ONLY : ppm_error,ppm_err_alloc,ppm_err_dealloc
      IMPLICIT NONE
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_MCMCParticlehtable)                   :: table
      !!! The hashtable.
      REAL(MK), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: proposal

      INTEGER,                             INTENT(IN   ) :: startindex
      INTEGER,                             INTENT(IN   ) :: endindex
      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      REAL(ppm_kind_double) :: t0

      INTEGER(ppm_kind_int64) :: i
      INTEGER                 :: nsize
      INTEGER                 :: info

      CHARACTER(LEN=*), PARAMETER :: caller="MCMCParticle_hash_packproposal_"

      !-------------------------------------------------------------------------
      !  Initialize
      !-------------------------------------------------------------------------
      CALL substart(caller,t0,info)

      IF (table%nrow.GT.0_ppm_kind_int64) THEN
         nsize=startindex-1
         DO i=1_ppm_kind_int64,table%nrow
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
      END FUNCTION MCMCParticle_hash_packproposal_

      FUNCTION MCMCParticle_hash_elementAt_2D3D(table,elementIndex,elementKey) RESULT(value)
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
      INTEGER :: nsize
      !-------------------------------------------------------------------------
      !  Initialize
      !-------------------------------------------------------------------------
      nsize=INT(table%nrow)
      IF (nsize.GT.0) THEN
         IF (elementIndex.LE.nsize) THEN
            elementIndex_=0
            DO i=1,nsize
               IF (table%keys(i).NE.htable_null_li) THEN
                  elementIndex_=elementIndex_+1
                  IF (elementIndex_.EQ.elementIndex) THEN
                     value=table%borders_pos(i)
                     IF (PRESENT(elementKey)) THEN
                        SELECT CASE (SIZE(elementKey))
                        CASE (2)
                           CALL HashIndexFunctor_label_2d(table%keys(i),elementKey(1),elementKey(2))
                        CASE (3)
                           CALL HashIndexFunctor_label_3d(table%keys(i),elementKey(1),elementKey(2),elementKey(3))
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
      END FUNCTION MCMCParticle_hash_elementAt_2D3D

      FUNCTION MCMCParticle_hash_elementAt_2D(table,elementIndex,coord1,coord2) RESULT(value)
      !!! This function returns the elementIndex element in the hash table
      !!! If the elementKey is present it gets the coordinates from Hash key
      IMPLICIT NONE
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_MCMCParticlehtable) :: table
      !!! The hashtable. The pointer must not be NULL
      INTEGER,          INTENT(IN   ) :: elementIndex
      INTEGER,          INTENT(  OUT) :: coord1
      INTEGER,          INTENT(  OUT) :: coord2
      !!! coordinates which has been hashed as a hash key

      TYPE(MCMCParticle)              :: Value
      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER :: i,elementIndex_
      INTEGER :: nsize
      !-------------------------------------------------------------------------
      !  Initialize
      !-------------------------------------------------------------------------
      nsize=INT(table%nrow)
      IF (nsize.GT.0) THEN
         IF (elementIndex.LE.nsize) THEN
            elementIndex_=0
            DO i=1,nsize
               IF (table%keys(i).NE.htable_null_li) THEN
                  elementIndex_=elementIndex_+1
                  IF (elementIndex_.EQ.elementIndex) THEN
                     value=table%borders_pos(i)
                     CALL HashIndexFunctor_label_2d(table%keys(i),coord1,coord2)
                     RETURN
                  ENDIF !elementIndex_.EQ.elementIndex
               ENDIF !table%keys(i).NE.htable_null_li
            ENDDO !i=1,table%nrow
         ENDIF !elementIndex.LE.table%nrow
      ENDIF !table%nrow.GT.0
      Value=MCMCParticle(-1,zero)
      coord1=-1
      coord2=-1
      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      RETURN
      END FUNCTION MCMCParticle_hash_elementAt_2D

      FUNCTION MCMCParticle_hash_elementAt_3D(table,elementIndex,coord1,coord2,coord3) RESULT(value)
      !!! This function returns the elementIndex element in the hash table
      !!! If the elementKey is present it gets the coordinates from Hash key
      IMPLICIT NONE
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_MCMCParticlehtable) :: table
      !!! The hashtable. The pointer must not be NULL
      INTEGER,          INTENT(IN   ) :: elementIndex
      INTEGER,          INTENT(  OUT) :: coord1
      INTEGER,          INTENT(  OUT) :: coord2
      INTEGER,          INTENT(  OUT) :: coord3
      !!! coordinates which has been hashed as a hash key

      TYPE(MCMCParticle)              :: Value
      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER :: i,elementIndex_
      INTEGER :: nsize
      !-------------------------------------------------------------------------
      !  Initialize
      !-------------------------------------------------------------------------
      nsize=INT(table%nrow)
      IF (nsize.GT.0) THEN
         IF (elementIndex.LE.nsize) THEN
            elementIndex_=0
            DO i=1,nsize
               IF (table%keys(i).NE.htable_null_li) THEN
                  elementIndex_=elementIndex_+1
                  IF (elementIndex_.EQ.elementIndex) THEN
                     value=table%borders_pos(i)
                     CALL HashIndexFunctor_label_3d(table%keys(i),coord1,coord2,coord3)
                     RETURN
                  ENDIF !elementIndex_.EQ.elementIndex
               ENDIF !table%keys(i).NE.htable_null_li
            ENDDO !i=1,table%nrow
         ENDIF !elementIndex.LE.table%nrow
      ENDIF !table%nrow.GT.0
      Value=MCMCParticle(-1,zero)
      coord1=-1
      coord2=-1
      coord3=-1
      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      RETURN
      END FUNCTION MCMCParticle_hash_elementAt_3D

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

      CHARACTER(LEN=ppm_char) :: caller='HashIndex_create_htable'

      CALL substart(caller,t0,info)

      IF (nelement.LE.0) THEN
         table%nrow = 1_ppm_kind_int64
      ELSE
         table%nrow = 2_ppm_kind_int64**(CEILING(LOG(REAL(nelement))/LOG(2.0)))
      ENDIF

      !---------------------------------------------------------------------
      !  Allocate array for hash table keys.
      !  Set everything to NULL.
      !---------------------------------------------------------------------
      ALLOCATE(table%keys(table%nrow),SOURCE=htable_null_li,STAT=Info)
      or_fail_alloc('Failed to alocate htable_keys!',ppm_error=ppm_error_fatal)

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
      IF (ALLOCATED(table%keys)) THEN
         DEALLOCATE(table%keys,STAT=info)
         or_fail_dealloc('Failed to deallocate htable_keys!',ppm_error=ppm_error_fatal)
      ENDIF

      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      END SUBROUTINE HashIndex_destroy_htable

      FUNCTION HashIndex_h_key(table,spot1,spot2,jump) RESULT(address)
      !!! Given the spots and jump value, returns corresponding address on
      !!! "borders" array.

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_HashIndextable), INTENT(IN   ) :: table
      !!! The hashtable to create. The pointer must not be NULL

      INTEGER(ppm_kind_int64),      INTENT(IN   ) :: spot1
      !!! First spot value to avoid double computing
      INTEGER(ppm_kind_int64),      INTENT(IN   ) :: spot2
      !!! Second spot value to avoid double computing
      INTEGER(ppm_kind_int64),      INTENT(IN   ) :: jump
      !!! Jump value for double hashing
      INTEGER(ppm_kind_int64)                     :: address
      !!! Address that corresponds to given key
      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      ! we need to work with 64 bit integers because
      ! jump*h_func might overflow, and as there is no unsigned
      ! integer the resulting value might be negative and not
      ! a valid array index
      address=MOD(spot1+spot2*jump,table%nrow)+1_ppm_kind_int64
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
      INTEGER(ppm_kind_int64) :: jump
      INTEGER(ppm_kind_int64) :: spot
      INTEGER(ppm_kind_int64) :: fspot
      INTEGER(ppm_kind_int64) :: sspot

      info = 0
      jump = 0_ppm_kind_int64

      ! Get the address corresponding to given key
      fspot=HashFunc(key,seed1,table%nrow)
      spot=fspot+1_ppm_kind_int64

      ! Keep on searching withing bounds of hash table.
      DO WHILE (jump.LT.table%nrow)
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
         jump=jump + 1_ppm_kind_int64
         IF (jump.EQ.1_ppm_kind_int64) THEN
            sspot=HashFunc(key,seed2,table%nrow)
         ENDIF

         spot=table%h_key(fspot,sspot,jump)
      ENDDO

      ! If NOT returned within the while-loop, that means our hash table is
      ! not sufficiently large. Hence, regrowth will be needed!
      CALL table%grow(info,key)
      RETURN
      END SUBROUTINE HashIndex_hash_insert

      SUBROUTINE HashIndex_hash_insert_2D3D(table,key_2D3D,info)

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_HashIndextable)         :: table
      !!! The hashtable

      INTEGER, DIMENSION(:), INTENT(IN   ) :: key_2D3D
      !!! Key to be stored
      INTEGER,               INTENT(  OUT) :: info
      !!! Info for whether insertion was successful or not. 0 if SUCCESSFUL
      !!! and -1 otherwise.

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: key
      INTEGER                 :: ssize

      ssize=SIZE(key_2D3D)
      SELECT CASE (ssize)
      CASE (2)
         key=IndexHashFunctor64_2d(key_2D3D)
      CASE (3)
         key=IndexHashFunctor64_3d(key_2D3D)
      CASE DEFAULT
         info=ppm_error_fatal
         RETURN
      END SELECT
      CALL table%HashIndex_hash_insert(key,info)
      RETURN
      END SUBROUTINE HashIndex_hash_insert_2D3D

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
      INTEGER(ppm_kind_int64) :: jump
      INTEGER(ppm_kind_int64) :: spot
      INTEGER(ppm_kind_int64) :: fspot
      INTEGER(ppm_kind_int64) :: sspot

      jump=0_ppm_kind_int64

      ! Get the other key that results in same
      ! hash key as for the inputkey.
      fspot=HashFunc(key,seed1,table%nrow)
      spot=fspot+1_ppm_kind_int64

      ! Keep on searching while we don't come across a NULL value or we don't
      ! exceed bounds of hash table.
      loop: DO WHILE(jump.LT.table%nrow)
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
         jump=jump+1_ppm_kind_int64
         IF (jump.EQ.1_ppm_kind_int64) THEN
            sspot=HashFunc(key,seed2,table%nrow)
         ENDIF

         spot=table%h_key(fspot,sspot,jump)
      ENDDO loop
      HashIndex_hash_search = .FALSE.
      RETURN
      END FUNCTION HashIndex_hash_search

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
      INTEGER(ppm_kind_int64) :: jump
      INTEGER(ppm_kind_int64) :: spot
      INTEGER(ppm_kind_int64) :: fspot
      INTEGER(ppm_kind_int64) :: sspot

      info = 0
      IF (PRESENT(existed)) THEN
         IF (existed) THEN
            jump=0_ppm_kind_int64

            ! Get the address corresponding to given key
            fspot=HashFunc(key,seed1,table%nrow)
            spot=fspot+1_ppm_kind_int64

            ! Keep on searching withing bounds of hash table.
            DO WHILE (jump.LT.table%nrow)
               ! If an empty slot found ...
               IF (table%keys(spot).EQ.key) THEN
                  !Remove the key and RETURN.
                  table%keys(spot)=htable_null_li
                  RETURN
               ENDIF
               ! If the current slot is occupied, jump to next key that results
               ! in same hash function.
               jump=jump+1_ppm_kind_int64
               IF (jump.EQ.1_ppm_kind_int64) THEN
                  sspot=HashFunc(key,seed2,table%nrow)
               ENDIF

               spot=table%h_key(fspot,sspot,jump)
            ENDDO
            ! If NOT returned within the while-loop, that means the key was NOT found
            info = ppm_error_fatal
         ENDIf
         RETURN
      ELSE
         IF (ANY(key.EQ.table%keys)) THEN
            jump=0_ppm_kind_int64

            ! Get the address corresponding to given key
            fspot=HashFunc(key,seed1,table%nrow)
            spot=fspot+1_ppm_kind_int64

            ! Keep on searching withing bounds of hash table.
            DO WHILE (jump.LT.table%nrow)
               ! If an empty slot found ...
               IF (table%keys(spot).EQ.key) THEN
                  !Remove the key and RETURN.
                  table%keys(spot)=htable_null_li
                  RETURN
               ENDIF
               ! If the current slot is occupied, jump to next key that results
               ! in same hash function.
               jump=jump+1_ppm_kind_int64
               IF (jump.EQ.1_ppm_kind_int64) THEN
                  sspot=HashFunc(key,seed2,table%nrow)
               ENDIF

               spot=table%h_key(fspot,sspot,jump)
            ENDDO
            ! If NOT returned within the while-loop, that means the key was NOT found
            info = ppm_error_fatal
         ENDIF
      ENDIF
      RETURN
      END SUBROUTINE HashIndex_hash_remove

      SUBROUTINE HashIndex_hash_remove_2D3D(table,key_2D3D, info,existed)

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_HashIndextable)         :: table
      !!! The hashtable

      INTEGER, DIMENSION(:), INTENT(IN   ) :: key_2D3D
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

      ssize=SIZE(key_2D3D)

      SELECT CASE (ssize)
      CASE (2)
         key=IndexHashFunctor64_2d(key_2D3D)
      CASE (3)
         key=IndexHashFunctor64_3d(key_2D3D)
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
      END SUBROUTINE HashIndex_hash_remove_2D3D

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

      IF (table%nrow*2_ppm_kind_int64.GE.INT(ppm_big_i-1,ppm_kind_int64)) THEN
         !TOCHCECK
         fail("hashtable with more than 2^31-1 elements will fail",ppm_error=ppm_error_fatal)
      ENDIF

      nsize=INT(table%nrow)

      SELECT CASE (nsize)
      CASE (0)
         CALL table%create(1,info)
         or_fail("table%create")
      CASE DEFAULT
         ALLOCATE(keys_tmp(nsize),SOURCE=table%keys,STAT=info)
         or_fail_alloc("keys_tmp")

         NULLIFY(ranklist)
         CALL ppm_util_qsort(keys_tmp,ranklist,info,nsize)
         or_fail("ppm_util_qsort")

         CALL table%destroy(info)
         or_fail("table%destroy")

         CALL table%create(nsize*2,info)
         or_fail("table%create")

         DO i=1,nsize
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
      END SELECT

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

      nsize=INT(table%nrow)
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

         IF (PRESENT(nsize)) THEN
            IF (nsize.LT.ssize) THEN
               l=0
               DO i=1,ssize
                  IF (Particlehtable(i)%nrow.GT.0_ppm_kind_int64) CYCLE
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

         l=INT(Particlehtable(0)%nrow)
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
            l=INT(Particlehtable(i)%nrow)
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

         l=INT(Particlehtabletmp(0)%nrow)
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
               l=INT(Particlehtabletmp(i)%nrow)
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
               l=INT(Particlehtabletmp(i)%nrow)
               IF (l.LE.0) CYCLE

               CALL Particlehtable(i)%create(l,info)
               or_fail("Particlehtable(i)%create")

               FORALL (j=1:l)
                  Particlehtable(i)%keys(j)          =Particlehtabletmp(i)%keys(j)
                  Particlehtable(i)%borders_pos(j)=Particlehtabletmp(i)%borders_pos(j)
               END FORALL

               CALL Particlehtabletmp(i)%destroy(info)
               or_fail("Particlehtabletmp(i)%destroy")
            ENDDO
         ENDIF

         DEALLOCATE(Particlehtabletmp,STAT=info)
         or_fail_dealloc("Particlehtabletmp")
      ELSE
         nsize_=15

         ALLOCATE(Particlehtable(0:nsize_),STAT=info)
         or_fail_alloc("Particlehtable")
      ENDIF !(ALLOCATED(Particlehtable))

      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      END SUBROUTINE MCMCParticle_realloc

      SUBROUTINE MCMCResults_create_htable(table,nelement,info)
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
      CLASS(ppm_rc_MCMCResultshtable) :: table
      !!! The hashtable to create.
      INTEGER,         INTENT(IN   ) :: nelement
      !!! Number of desired elements.
      INTEGER,         INTENT(  OUT) :: info

      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double) :: t0

      CHARACTER(LEN=ppm_char) :: caller='MCMCResults_create_htable'

      CALL substart(caller,t0,info)

      IF (nelement.LE.0) THEN
         table%nrow = 1_ppm_kind_int64
      ELSE
         table%nrow = 2_ppm_kind_int64**(CEILING(LOG(REAL(nelement))/LOG(2.0)))
      ENDIF

      !---------------------------------------------------------------------
      !  Allocate array for hash table keys and array for positions on "borders" array.
      ! &
      !  Set everything to NULL.
      !---------------------------------------------------------------------
      ALLOCATE(table%keys(table%nrow),SOURCE=htable_null_li,STAT=Info)
      or_fail_alloc('Failed to alocate htable_keys!',ppm_error=ppm_error_fatal)

      ALLOCATE(table%borders_pos(table%nrow),STAT=Info)
      or_fail_alloc('Failed to alocate htable_borders_pos!',ppm_error=ppm_error_fatal)

      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      END SUBROUTINE MCMCResults_create_htable

      SUBROUTINE MCMCResults_destroy_htable(table,info,lcomplete)

      USE ppm_module_data
      USE ppm_module_alloc
      USE ppm_module_error
      USE ppm_module_substart
      USE ppm_module_substop
      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      CLASS(ppm_rc_MCMCResultshtable)  :: table
      !!! The hashtable to create. The pointer must not be NULL
      INTEGER,           INTENT(  OUT) :: info
      LOGICAL, OPTIONAL, INTENT(IN   ) :: lcomplete
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double) :: t0

      INTEGER(ppm_kind_int64) :: i

      CHARACTER(LEN=ppm_char) :: caller='MCMCResults_destroy_htable'

      CALL substart(caller,t0,info)

      !---------------------------------------------------------------------
      !  Deallocate array for hash table keys and array for positions on "borders" array.
      !---------------------------------------------------------------------
      IF (ALLOCATED(table%keys)) THEN
         IF (PRESENT(lcomplete)) THEN
            IF (lcomplete) THEN
               DO i=1_ppm_kind_int64,table%nrow
                  CALL table%borders_pos(i)%destroy()
               ENDDO
            ENDIF
         ENDIF

         DEALLOCATE(table%keys,table%borders_pos,STAT=info)
         or_fail_dealloc('Failed to deallocate htable_keys & htable_borders_pos!',ppm_error=ppm_error_fatal)
      ENDIF

      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      END SUBROUTINE MCMCResults_destroy_htable

      FUNCTION MCMCResults_h_key(table,spot1,spot2,jump) RESULT(address)
      !!! Given the spots and jump value, returns corresponding address on
      !!! "borders" array.

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_MCMCResultshtable), INTENT(IN   ) :: table
      !!! The hashtable to create. The pointer must not be NULL

      INTEGER(ppm_kind_int64),          INTENT(IN   ) :: spot1
      !!! First spot value to avoid double computing
      INTEGER(ppm_kind_int64),          INTENT(IN   ) :: spot2
      !!! Second spot value to avoid double computing
      INTEGER(ppm_kind_int64),          INTENT(IN   ) :: jump
      !!! Jump value for double hashing
      INTEGER(ppm_kind_int64)                         :: address
      !!! Address that corresponds to given key
      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      ! we need to work with 64 bit integers because
      ! jump*h_func might overflow, and as there is no unsigned
      ! integer the resulting value might be negative and not
      ! a valid array index
      address=MOD(spot1+spot2*jump,table%nrow)+1_ppm_kind_int64
      RETURN
      END FUNCTION MCMCResults_h_key

      SUBROUTINE MCMCResults_hash_insert(table,key,value1,value2,info)
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
      CLASS(ppm_rc_MCMCResultshtable)       :: table
      !!! The hashtable
      INTEGER(ppm_kind_int64), INTENT(IN   ) :: key
      !!! Key to be stored
      INTEGER(ppm_kind_int64), INTENT(IN   ) :: value1
      INTEGER,                 INTENT(IN   ) :: value2
      !!! Value that corresponds to given key
      INTEGER,                 INTENT(  OUT) :: info
      !!! Info for whether insertion was successful or not. 0 if SUCCESSFUL
      !!! and -1 otherwise.

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: jump
      INTEGER(ppm_kind_int64) :: spot
      INTEGER(ppm_kind_int64) :: fspot
      INTEGER(ppm_kind_int64) :: sspot

      info = 0
      jump = 0_ppm_kind_int64

      ! Get the address corresponding to given key
      fspot=HashFunc(key,seed1,table%nrow)
      spot=fspot+1_ppm_kind_int64

      ! Keep on searching withing bounds of hash table.
      DO WHILE (jump.LT.table%nrow)
         ! If an empty slot found ...
         IF (table%keys(spot).EQ.htable_null_li) THEN
            ! Store the key and the corresponding value and RETURN.
            table%keys(spot) = key
            CALL table%borders_pos(spot)%add(value1,value2)
            RETURN
         !If the key is the same the values should be added at the end
         ELSE IF (table%keys(spot).EQ.key) THEN
            CALL table%borders_pos(spot)%add(value1,value2)
            RETURN
         ENDIF
         ! If the current slot is occupied, jump to next key that results
         ! in same hash function.
         jump=jump + 1_ppm_kind_int64
         IF (jump.EQ.1_ppm_kind_int64) THEN
            sspot=HashFunc(key,seed2,table%nrow)
         ENDIF

         spot=table%h_key(fspot,sspot,jump)
      ENDDO

      ! If NOT returned within the while-loop, that means our hash table is
      ! not sufficiently large. Hence, regrowth will be needed!
      CALL table%grow(info,key,value1,value2)
      RETURN
      END SUBROUTINE MCMCResults_hash_insert

      SUBROUTINE MCMCResults_hash_insert_2D3D(table,key_2D3D,value1,value2,info)

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_MCMCResultshtable)     :: table
      !!! The hashtable
      INTEGER,   DIMENSION(:), INTENT(IN   ) :: key_2D3D
      !!! Key to be stored
      INTEGER(ppm_kind_int64), INTENT(IN   ) :: value1
      INTEGER,                 INTENT(IN   ) :: value2
      !!! Value that corresponds to given key
      INTEGER,                 INTENT(  OUT) :: info
      !!! Info for whether insertion was successful or not. 0 if SUCCESSFUL
      !!! and -1 otherwise.

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: jump
      INTEGER(ppm_kind_int64) :: spot
      INTEGER(ppm_kind_int64) :: fspot
      INTEGER(ppm_kind_int64) :: sspot
      INTEGER(ppm_kind_int64) :: key
      INTEGER                 :: ssize

      ssize=SIZE(key_2D3D)
      SELECT CASE (ssize)
      CASE (2)
         info=0

         ! Get the address corresponding to given key
         fspot=HashFunc_XY(key_2D3D(1),key_2D3D(2),seed1,table%nrow,key)
      CASE (3)
         info=0

         ! Get the address corresponding to given key
         fspot=HashFunc_XYZ(key_2D3D(1),key_2D3D(2),key_2D3D(3),seed1,table%nrow,key)
      CASE DEFAULT
         info=ppm_error_fatal
         RETURN
      END SELECT

      jump = 0_ppm_kind_int64

      spot=fspot+1_ppm_kind_int64

      ! Keep on searching withing bounds of hash table.
      DO WHILE (jump.LT.table%nrow)
         ! If an empty slot found ...
         IF (table%keys(spot).EQ.htable_null_li) THEN
            ! Store the key and the corresponding value and RETURN.
            table%keys(spot) = key
            CALL table%borders_pos(spot)%add(value1,value2)
            RETURN
         !If the key is the same the value should be updated
         ELSE IF (table%keys(spot).EQ.key) THEN
            CALL table%borders_pos(spot)%add(value1,value2)
            RETURN
         ENDIF
         ! If the current slot is occupied, jump to next key that results
         ! in same hash function.
         jump=jump + 1_ppm_kind_int64
         IF (jump.EQ.1_ppm_kind_int64) THEN
            sspot=HashFunc(key,seed2,table%nrow)
         ENDIF

         spot=table%h_key(fspot,sspot,jump)
      ENDDO

      ! If NOT returned within the while-loop, that means our hash table is
      ! not sufficiently large. Hence, regrowth will be needed!
      CALL table%grow(info,key,value1,value2)
      RETURN
      END SUBROUTINE MCMCResults_hash_insert_2D3D

      SUBROUTINE MCMCResults_hash_insert_2d(table,key_1,key_2,value1,value2,info)

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_MCMCResultshtable)        :: table
      !!! The hashtable
      INTEGER,                 INTENT(IN   ) :: key_1
      INTEGER,                 INTENT(IN   ) :: key_2
      !!! Key to be stored
      INTEGER(ppm_kind_int64), INTENT(IN   ) :: value1
      INTEGER,                 INTENT(IN   ) :: value2
      !!! Value that corresponds to given key
      INTEGER,                 INTENT(  OUT) :: info
      !!! Info for whether insertion was successful or not. 0 if SUCCESSFUL
      !!! and -1 otherwise.

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: jump
      INTEGER(ppm_kind_int64) :: spot
      INTEGER(ppm_kind_int64) :: fspot
      INTEGER(ppm_kind_int64) :: sspot
      INTEGER(ppm_kind_int64) :: key

      info = 0
      jump = 0_ppm_kind_int64

      ! Get the address corresponding to given key
      fspot=HashFunc_XY(key_1,key_2,seed1,table%nrow,key)

      spot=fspot+1_ppm_kind_int64

      ! Keep on searching withing bounds of hash table.
      DO WHILE (jump.LT.table%nrow)
         ! If an empty slot found ...
         IF (table%keys(spot).EQ.htable_null_li) THEN
            ! Store the key and the corresponding value and RETURN.
            table%keys(spot) = key
            CALL table%borders_pos(spot)%add(value1,value2)
            RETURN
         !If the key is the same the value should be updated
         ELSE IF (table%keys(spot).EQ.key) THEN
            CALL table%borders_pos(spot)%add(value1,value2)
            RETURN
         ENDIF
         ! If the current slot is occupied, jump to next key that results
         ! in same hash function.
         jump=jump + 1_ppm_kind_int64
         IF (jump.EQ.1_ppm_kind_int64) THEN
            sspot=HashFunc(key,seed2,table%nrow)
         ENDIF

         spot=table%h_key(fspot,sspot,jump)
      ENDDO

      ! If NOT returned within the while-loop, that means our hash table is
      ! not sufficiently large. Hence, regrowth will be needed!
      CALL table%grow(info,key,value1,value2)
      RETURN
      END SUBROUTINE MCMCResults_hash_insert_2d

      SUBROUTINE MCMCResults_hash_insert_3d(table,key_1,key_2,key_3,value1,value2,info)

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_MCMCResultshtable)  :: table
      !!! The hashtable
      INTEGER,                 INTENT(IN   ) :: key_1
      INTEGER,                 INTENT(IN   ) :: key_2
      INTEGER,                 INTENT(IN   ) :: key_3
      !!! Key to be stored
      INTEGER(ppm_kind_int64), INTENT(IN   ) :: value1
      INTEGER,                 INTENT(IN   ) :: value2
      !!! Value that corresponds to given key
      INTEGER,                 INTENT(  OUT) :: info
      !!! Info for whether insertion was successful or not. 0 if SUCCESSFUL
      !!! and -1 otherwise.

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: jump
      INTEGER(ppm_kind_int64) :: spot
      INTEGER(ppm_kind_int64) :: fspot
      INTEGER(ppm_kind_int64) :: sspot
      INTEGER(ppm_kind_int64) :: key

      info = 0
      jump = 0_ppm_kind_int64

      ! Get the address corresponding to given key
      fspot=HashFunc_XYZ(key_1,key_2,key_3,seed1,table%nrow,key)

      spot=fspot+1_ppm_kind_int64

      ! Keep on searching withing bounds of hash table.
      DO WHILE (jump.LT.table%nrow)
         ! If an empty slot found ...
         IF (table%keys(spot).EQ.htable_null_li) THEN
            ! Store the key and the corresponding value and RETURN.
            table%keys(spot) = key
            CALL table%borders_pos(spot)%add(value1,value2)
            RETURN
         !If the key is the same the value should be updated
         ELSE IF (table%keys(spot).EQ.key) THEN
            CALL table%borders_pos(spot)%add(value1,value2)
            RETURN
         ENDIF
         ! If the current slot is occupied, jump to next key that results
         ! in same hash function.
         jump=jump + 1_ppm_kind_int64
         IF (jump.EQ.1_ppm_kind_int64) THEN
            sspot=HashFunc(key,seed2,table%nrow)
         ENDIF

         spot=table%h_key(fspot,sspot,jump)
      ENDDO

      ! If NOT returned within the while-loop, that means our hash table is
      ! not sufficiently large. Hence, regrowth will be needed!
      CALL table%grow(info,key,value1,value2)
      RETURN
      END SUBROUTINE MCMCResults_hash_insert_3d

      FUNCTION MCMCResults_hash_search(table,key) RESULT(value)
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
      CLASS(ppm_rc_MCMCResultshtable), INTENT(IN   ) :: table
      !!! The hashtable to create. The pointer must not be NULL
      INTEGER(ppm_kind_int64),          INTENT(IN   ) :: key
      !!! Input key, which the corresponding value is asked for
      TYPE(ppm_rc_list_2)                             :: value
      !!! Value corresponding to the input key

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: jump
      INTEGER(ppm_kind_int64) :: spot
      INTEGER(ppm_kind_int64) :: fspot
      INTEGER(ppm_kind_int64) :: sspot

      LOGICAL :: KeyExist

      jump=0_ppm_kind_int64

      KeyExist=.TRUE.

      ! Get the other key that results in same
      ! hash key as for the inputkey.
      fspot=HashFunc(key,seed1,table%nrow)
      spot=fspot+1_ppm_kind_int64

      ! Keep on searching while we don't come across a NULL value or we don't
      ! exceed bounds of hash table.
      loop: DO WHILE(jump.LT.table%nrow)
          ! If key matches ...
          IF (table%keys(spot).EQ.key) THEN
             ! Set the return value and return
             value%first => table%borders_pos(spot)%first
             value%last  => table%borders_pos(spot)%last
             RETURN
          ELSE IF (table%keys(spot).EQ.htable_null_li) THEN
             IF (KeyExist) THEN
                KeyExist=.FALSE.
                IF (.NOT.ANY(key.EQ.table%keys)) EXIT loop
             ENDIF
          ENDIF
          ! Otherwise, keep on incrementing jump distance
          jump=jump+1_ppm_kind_int64
          IF (jump.EQ.1_ppm_kind_int64) THEN
             sspot=HashFunc(key,seed2,table%nrow)
          ENDIF

          spot=table%h_key(fspot,sspot,jump)
      ENDDO loop
      value%first => NULL()
      value%last  => NULL()
      RETURN
      END FUNCTION MCMCResults_hash_search

      FUNCTION MCMCResults_hash_search_2d(table,key_1,key_2) RESULT(value)

      IMPLICIT NONE
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_MCMCResultshtable), INTENT(IN   ) :: table
      !!! The hashtable

      INTEGER,                          INTENT(IN   ) :: key_1
      INTEGER,                          INTENT(IN   ) :: key_2
      !!! Input key, which the corresponding value is asked for

      TYPE(ppm_rc_list_2)                             :: value
      !!! Value corresponding to the input key
      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: jump
      INTEGER(ppm_kind_int64) :: spot
      INTEGER(ppm_kind_int64) :: fspot
      INTEGER(ppm_kind_int64) :: sspot
      INTEGER(ppm_kind_int64) :: key

      LOGICAL :: KeyExist

      jump=0_ppm_kind_int64

      KeyExist=.TRUE.

      ! Get the address corresponding to given key
      fspot=HashFunc_XY(key_1,key_2,seed1,table%nrow,key)

      spot=fspot+1_ppm_kind_int64

      ! Keep on searching while we don't come across a NULL value or we don't
      ! exceed bounds of hash table.
      loop: DO WHILE(jump.LT.table%nrow)
          ! If key matches ...
          IF (table%keys(spot).EQ.key) THEN
             ! Set the return value and return
             value%first => table%borders_pos(spot)%first
             value%last  => table%borders_pos(spot)%last
             RETURN
          ELSE IF (table%keys(spot).EQ.htable_null_li) THEN
             IF (KeyExist) THEN
                KeyExist=.FALSE.
                IF (.NOT.ANY(key.EQ.table%keys)) EXIT loop
             ENDIF
          ENDIF
          ! Otherwise, keep on incrementing jump distance
          jump=jump+1_ppm_kind_int64
          IF (jump.EQ.1_ppm_kind_int64) THEN
             sspot=HashFunc(key,seed2,table%nrow)
          ENDIF

          spot=table%h_key(fspot,sspot,jump)
      ENDDO loop
      value%first => NULL()
      value%last  => NULL()
      RETURN
      END FUNCTION MCMCResults_hash_search_2d

      FUNCTION MCMCResults_hash_search_3d(table,key_1,key_2,key_3) RESULT(value)

      IMPLICIT NONE
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_MCMCResultshtable), INTENT(IN   ) :: table
      !!! The hashtable

      INTEGER,                          INTENT(IN   ) :: key_1
      INTEGER,                          INTENT(IN   ) :: key_2
      INTEGER,                          INTENT(IN   ) :: key_3
      !!! Input key, which the corresponding value is asked for

      TYPE(ppm_rc_list_2)                              :: value
      !!! Value corresponding to the input key
      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: jump
      INTEGER(ppm_kind_int64) :: spot
      INTEGER(ppm_kind_int64) :: fspot
      INTEGER(ppm_kind_int64) :: sspot
      INTEGER(ppm_kind_int64) :: key

      LOGICAL :: KeyExist

      jump=0_ppm_kind_int64

      KeyExist=.TRUE.

      ! Get the address corresponding to given key
      fspot=HashFunc_XYZ(key_1,key_2,key_3,seed1,table%nrow,key)

      spot=fspot+1_ppm_kind_int64

      ! Keep on searching while we don't come across a NULL value or we don't
      ! exceed bounds of hash table.
      loop: DO WHILE(jump.LT.table%nrow)
          ! If key matches ...
          IF (table%keys(spot).EQ.key) THEN
             ! Set the return value and return
             value%first => table%borders_pos(spot)%first
             value%last  => table%borders_pos(spot)%last
             RETURN
          ELSE IF (table%keys(spot).EQ.htable_null_li) THEN
             IF (KeyExist) THEN
                KeyExist=.FALSE.
                IF (.NOT.ANY(key.EQ.table%keys)) EXIT loop
             ENDIF
          ENDIF
          ! Otherwise, keep on incrementing jump distance
          jump=jump+1_ppm_kind_int64
          IF (jump.EQ.1_ppm_kind_int64) THEN
             sspot=HashFunc(key,seed2,table%nrow)
          ENDIF

          spot=table%h_key(fspot,sspot,jump)
      ENDDO loop
      value%first => NULL()
      value%last  => NULL()
      RETURN
      END FUNCTION MCMCResults_hash_search_3d

      !TODO check this sub
      !Yaser
      !This function PROBABLY suffers from a bug!
      SUBROUTINE MCMCResults_hash_remove(table,key,info,existed)
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
      CLASS(ppm_rc_MCMCResultshtable)        :: table
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
      INTEGER(ppm_kind_int64) :: jump
      INTEGER(ppm_kind_int64) :: spot
      INTEGER(ppm_kind_int64) :: fspot
      INTEGER(ppm_kind_int64) :: sspot

      info = 0

      IF (PRESENT(existed)) THEN
         IF (existed) THEN
            jump=0_ppm_kind_int64

            ! Get the other key that results in same
            ! hash key as for the inputkey.
            fspot=HashFunc(key,seed1,table%nrow)
            spot=fspot+1_ppm_kind_int64

            ! Keep on searching withing bounds of hash table.
            DO WHILE (jump.LT.table%nrow)
               ! If an empty slot found ...
               IF (table%keys(spot).EQ.key) THEN
                  !Remove the key and the corresponding value and RETURN.
                  CALL table%borders_pos(spot)%destroy()
                  table%keys(spot)=htable_null_li
                  RETURN
               ENDIF
               ! If the current slot is occupied, jump to next key that results
               ! in same hash function.
               jump=jump+1_ppm_kind_int64
               IF (jump.EQ.1_ppm_kind_int64) THEN
                  sspot=HashFunc(key,seed2,table%nrow)
               ENDIF

               spot=table%h_key(fspot,sspot,jump)
            ENDDO

            ! If NOT returned within the while-loop, that means the key was NOT found
            info = ppm_error_fatal
         ENDIF
         RETURN
      ELSE
         IF (ANY(key.EQ.table%keys)) THEN
            jump=0_ppm_kind_int64

            ! Get the other key that results in same
            ! hash key as for the inputkey.
            fspot=HashFunc(key,seed1,table%nrow)
            spot=fspot+1_ppm_kind_int64

            ! Keep on searching withing bounds of hash table.
            DO WHILE (jump.LT.table%nrow)
               ! If an empty slot found ...
               IF (table%keys(spot).EQ.key) THEN
                  !Remove the key and the corresponding value and RETURN.
                  CALL table%borders_pos(spot)%destroy()
                  table%keys(spot)=htable_null_li
                  RETURN
               ENDIF
               ! If the current slot is occupied, jump to next key that results
               ! in same hash function.
               jump=jump+1_ppm_kind_int64
               IF (jump.EQ.1_ppm_kind_int64) THEN
                  sspot=HashFunc(key,seed2,table%nrow)
               ENDIF

               spot=table%h_key(fspot,sspot,jump)
            ENDDO
            ! If NOT returned within the while-loop, that means the key was NOT found
            info = ppm_error_fatal
         ENDIF
      ENDIF
      RETURN
      END SUBROUTINE MCMCResults_hash_remove

      SUBROUTINE MCMCResults_hash_remove_2D3D(table,key_2D3D,info,existed)

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_MCMCResultshtable)     :: table
      !!! The hashtable
      INTEGER, DIMENSION(:), INTENT(IN   ) :: key_2D3D
      !!! Key to be removed
      INTEGER,               INTENT(  OUT) :: info
      !!! Info for whether removal was successful or not. 0 if SUCCESSFUL
      !!! and -1 otherwise.
      LOGICAL, OPTIONAL,     INTENT(IN   ) :: existed

      !---------------------------------------------------------------------
      !  Local variables
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64) :: jump
      INTEGER(ppm_kind_int64) :: spot
      INTEGER(ppm_kind_int64) :: fspot
      INTEGER(ppm_kind_int64) :: sspot
      INTEGER(ppm_kind_int64) :: key
      INTEGER                 :: ssize

      info = 0

      IF (PRESENT(existed)) THEN
         IF (existed) THEN
            jump=0_ppm_kind_int64

            ssize=SIZE(key_2D3D)
            SELECT CASE (ssize)
            CASE (2)
               ! Get the address corresponding to given key
               fspot=HashFunc_XY(key_2D3D(1),key_2D3D(2),seed1,table%nrow,key)
            CASE (3)
               ! Get the address corresponding to given key
               fspot=HashFunc_XYZ(key_2D3D(1),key_2D3D(2),key_2D3D(3),seed1,table%nrow,key)
            CASE DEFAULT
               info=ppm_error_fatal
               RETURN
            END SELECT

            spot=fspot+1_ppm_kind_int64

            ! Keep on searching withing bounds of hash table.
            DO WHILE (jump.LT.table%nrow)
               ! If an empty slot found ...
               IF (table%keys(spot).EQ.key) THEN
                  !Remove the key and the corresponding value and RETURN.
                  CALL table%borders_pos(spot)%destroy()
                  table%keys(spot)=htable_null_li
                  RETURN
               ENDIF
               ! If the current slot is occupied, jump to next key that results
               ! in same hash function.
               jump=jump+1_ppm_kind_int64
               IF (jump.EQ.1_ppm_kind_int64) THEN
                  sspot=HashFunc(key,seed2,table%nrow)
               ENDIF

               spot=table%h_key(fspot,sspot,jump)
            ENDDO

            ! If NOT returned within the while-loop, that means the key was NOT found
            info = ppm_error_fatal
         ENDIF
         RETURN
      ELSE
         ssize=SIZE(key_2D3D)
         SELECT CASE (ssize)
         CASE (2)
            key=IndexHashFunctor_label_2d(key_2D3D(1),key_2D3D(2))
         CASE (3)
            key=IndexHashFunctor_label_3d(key_2D3D(1),key_2D3D(2),key_2D3D(3))
         CASE DEFAULT
            info=ppm_error_fatal
            RETURN
         END SELECT

         IF (ANY(key.EQ.table%keys)) THEN
            jump=0_ppm_kind_int64

            ! Get the other key that results in same
            ! hash key as for the inputkey.
            fspot=HashFunc(key,seed1,table%nrow)
            spot=fspot+1_ppm_kind_int64

            ! Keep on searching withing bounds of hash table.
            DO WHILE (jump.LT.table%nrow)
               ! If an empty slot found ...
               IF (table%keys(spot).EQ.key) THEN
                  !Remove the key and the corresponding value and RETURN.
                  CALL table%borders_pos(spot)%destroy()
                  table%keys(spot)=htable_null_li
                  RETURN
               ENDIF
               ! If the current slot is occupied, jump to next key that results
               ! in same hash function.
               jump=jump+1_ppm_kind_int64
               IF (jump.EQ.1_ppm_kind_int64) THEN
                  sspot=HashFunc(key,seed2,table%nrow)
               ENDIF

               spot=table%h_key(fspot,sspot,jump)
            ENDDO
            ! If NOT returned within the while-loop, that means the key was NOT found
            info = ppm_error_fatal
         ENDIF
      ENDIF
      RETURN
      END SUBROUTINE MCMCResults_hash_remove_2D3D

      SUBROUTINE MCMCResults_hash_remove_2d(table,key_1,key_2,info,existed)

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_MCMCResultshtable) :: table
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
      INTEGER(ppm_kind_int64) :: jump
      INTEGER(ppm_kind_int64) :: spot
      INTEGER(ppm_kind_int64) :: fspot
      INTEGER(ppm_kind_int64) :: sspot
      INTEGER(ppm_kind_int64) :: key

      info = 0

      IF (PRESENT(existed)) THEN
         IF (existed) THEN
            jump=0_ppm_kind_int64

            ! Get the address corresponding to given key
            fspot=HashFunc_XY(key_1,key_2,seed1,table%nrow,key)

            spot=fspot+1_ppm_kind_int64

            ! Keep on searching withing bounds of hash table.
            DO WHILE (jump.LT.table%nrow)
               ! If an empty slot found ...
               IF (table%keys(spot).EQ.key) THEN
                  !Remove the key and the corresponding value and RETURN.
                  CALL table%borders_pos(spot)%destroy()
                  table%keys(spot)=htable_null_li
                  RETURN
               ENDIF
               ! If the current slot is occupied, jump to next key that results
               ! in same hash function.
               jump=jump+1_ppm_kind_int64
               IF (jump.EQ.1_ppm_kind_int64) THEN
                  sspot=HashFunc(key,seed2,table%nrow)
               ENDIF

               spot=table%h_key(fspot,sspot,jump)
            ENDDO

            ! If NOT returned within the while-loop, that means the key was NOT found
            info = ppm_error_fatal
         ENDIF
         RETURN
      ELSE
         key=IndexHashFunctor_label_2d(key_1,key_2)

         IF (ANY(key.EQ.table%keys)) THEN
            jump=0_ppm_kind_int64

            ! Get the other key that results in same
            ! hash key as for the inputkey.
            fspot=HashFunc(key,seed1,table%nrow)
            spot=fspot+1_ppm_kind_int64

            ! Keep on searching withing bounds of hash table.
            DO WHILE (jump.LT.table%nrow)
               ! If an empty slot found ...
               IF (table%keys(spot).EQ.key) THEN
                  !Remove the key and the corresponding value and RETURN.
                  CALL table%borders_pos(spot)%destroy()
                  table%keys(spot)=htable_null_li
                  RETURN
               ENDIF
               ! If the current slot is occupied, jump to next key that results
               ! in same hash function.
               jump=jump+1_ppm_kind_int64
               IF (jump.EQ.1_ppm_kind_int64) THEN
                  sspot=HashFunc(key,seed2,table%nrow)
               ENDIF

               spot=table%h_key(fspot,sspot,jump)
            ENDDO
            ! If NOT returned within the while-loop, that means the key was NOT found
            info = ppm_error_fatal
         ENDIF
      ENDIF
      RETURN
      END SUBROUTINE MCMCResults_hash_remove_2d

      SUBROUTINE MCMCResults_hash_remove_3d(table,key_1,key_2,key_3,info,existed)

      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_MCMCResultshtable) :: table
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
      INTEGER(ppm_kind_int64) :: jump
      INTEGER(ppm_kind_int64) :: spot
      INTEGER(ppm_kind_int64) :: fspot
      INTEGER(ppm_kind_int64) :: sspot
      INTEGER(ppm_kind_int64) :: key

      info = 0

      IF (PRESENT(existed)) THEN
         IF (existed) THEN
            jump=0_ppm_kind_int64

            ! Get the address corresponding to given key
            fspot=HashFunc_XYZ(key_1,key_2,key_3,seed1,table%nrow,key)

            spot=fspot+1_ppm_kind_int64

            ! Keep on searching withing bounds of hash table.
            DO WHILE (jump.LT.table%nrow)
               ! If an empty slot found ...
               IF (table%keys(spot).EQ.key) THEN
                  !Remove the key and the corresponding value and RETURN.
                  CALL table%borders_pos(spot)%destroy()
                  table%keys(spot)=htable_null_li
                  RETURN
               ENDIF
               ! If the current slot is occupied, jump to next key that results
               ! in same hash function.
               jump=jump+1_ppm_kind_int64
               IF (jump.EQ.1_ppm_kind_int64) THEN
                  sspot=HashFunc(key,seed2,table%nrow)
               ENDIF

               spot=table%h_key(fspot,sspot,jump)
            ENDDO

            ! If NOT returned within the while-loop, that means the key was NOT found
            info = ppm_error_fatal
         ENDIF
         RETURN
      ELSE
         key=IndexHashFunctor_label_3d(key_1,key_2,key_3)

         IF (ANY(key.EQ.table%keys)) THEN
            jump=0_ppm_kind_int64

            ! Get the other key that results in same
            ! hash key as for the inputkey.
            fspot=HashFunc(key,seed1,table%nrow)
            spot=fspot+1_ppm_kind_int64

            ! Keep on searching withing bounds of hash table.
            DO WHILE (jump.LT.table%nrow)
               ! If an empty slot found ...
               IF (table%keys(spot).EQ.key) THEN
                  !Remove the key and the corresponding value and RETURN.
                  CALL table%borders_pos(spot)%destroy()
                  table%keys(spot)=htable_null_li
                  RETURN
               ENDIF
               ! If the current slot is occupied, jump to next key that results
               ! in same hash function.
               jump=jump+1_ppm_kind_int64
               IF (jump.EQ.1_ppm_kind_int64) THEN
                  sspot=HashFunc(key,seed2,table%nrow)
               ENDIF

               spot=table%h_key(fspot,sspot,jump)
            ENDDO
            ! If NOT returned within the while-loop, that means the key was NOT found
            info = ppm_error_fatal
         ENDIF
      ENDIF
      RETURN
      END SUBROUTINE MCMCResults_hash_remove_3d

      SUBROUTINE MCMCResults_grow_htable(table,info,key_,value1_,value2_)
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
      CLASS(ppm_rc_MCMCResultshtable)                 :: table
      !!! The hashtable to grow.

      INTEGER,                           INTENT(  OUT) :: info

      INTEGER(ppm_kind_int64), OPTIONAL, INTENT(IN   ) :: key_
      !!! Key to be stored
      INTEGER(ppm_kind_int64), OPTIONAL, INTENT(IN   ) :: value1_
      INTEGER,                 OPTIONAL, INTENT(IN   ) :: value2_
      !!! Value that corresponds to given key
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      TYPE(ppm_rc_list_2), DIMENSION(:), ALLOCATABLE :: borders_pos_tmp

      REAL(ppm_kind_double) :: t0

      INTEGER(ppm_kind_int64), DIMENSION(:), ALLOCATABLE :: keys_tmp
      INTEGER,                 DIMENSION(:), POINTER     :: ranklist
      INTEGER                                            :: nsize
      INTEGER                                            :: nsize2
      INTEGER                                            :: i,j
      INTEGER(ppm_kind_int64)                            :: jump
      INTEGER(ppm_kind_int64)                            :: spot
      INTEGER(ppm_kind_int64)                            :: fspot
      INTEGER(ppm_kind_int64)                            :: sspot

      CHARACTER(LEN=ppm_char) :: caller='MCMCResults_grow_htable'

      CALL substart(caller,t0,info)

      IF (table%nrow*2_ppm_kind_int64.GE.INT(ppm_big_i-1,ppm_kind_int64)) THEN
         !TOCHCECK
         fail("hashtable with more than 2^31-1 elements will fail",ppm_error=ppm_error_fatal)
      ENDIF

      nsize=INT(table%nrow)

      SELECT CASE (nsize)
      CASE (0)
         CALL table%create(1,info)
         or_fail("table%create")
      CASE DEFAULT
         ALLOCATE(keys_tmp(nsize),SOURCE=table%keys,STAT=info)
         or_fail_alloc("keys_tmp")

         ALLOCATE(borders_pos_tmp(nsize),STAT=info)
         or_fail_alloc("borders_pos_tmp")

         DO i=1,nsize
            borders_pos_tmp(i)%first => table%borders_pos(i)%first
            borders_pos_tmp(i)%last  => table%borders_pos(i)%last
         ENDDO

         NULLIFY(ranklist)
         CALL ppm_util_qsort(keys_tmp,ranklist,info,nsize)
         or_fail("ppm_util_qsort")

         CALL table%destroy(info)
         or_fail("table%destroy")

         CALL table%create(nsize*2,info)
         or_fail("table%create")

      5555 CONTINUE

         DO i=1,nsize
            j=ranklist(i)
            IF (keys_tmp(j).EQ.htable_null_li) CYCLE

            jump = 0_ppm_kind_int64

            ! Get the address corresponding to given key
            fspot=HashFunc(keys_tmp(j),seed1,table%nrow)
            spot=fspot+1_ppm_kind_int64

            ! Keep on searching withing bounds of hash table.
            DO WHILE (jump.LT.table%nrow)
               ! If an empty slot found ...
               IF (table%keys(spot).EQ.htable_null_li) THEN
                  ! Store the key and the corresponding value and RETURN.
                  table%keys(spot) = keys_tmp(j)
                  table%borders_pos(spot)%first => borders_pos_tmp(j)%first
                  table%borders_pos(spot)%last  => borders_pos_tmp(j)%last
                  EXIT
               ELSE IF (table%keys(spot).EQ.keys_tmp(j)) THEN
                  IF (table%nrow*2_ppm_kind_int64.GE.INT(ppm_big_i-1,ppm_kind_int64)) THEN
                     fail("hashtable with more than 2^31-1 elements will fail",ppm_error=ppm_error_fatal)
                  ENDIF

                  nsize2=INT(table%nrow)*2

                  CALL table%destroy(info)
                  or_fail("table%destroy")

                  CALL table%create(nsize2,info)
                  or_fail("table%create")

                  GOTO 5555
               ENDIF
               ! If the current slot is occupied, jump to next key that results
               ! in same hash function.
               jump=jump + 1_ppm_kind_int64
               IF (jump.EQ.1_ppm_kind_int64) THEN
                  sspot=HashFunc(keys_tmp(j),seed2,table%nrow)
               ENDIF
               spot=table%h_key(fspot,sspot,jump)
            ENDDO ! WHILE (jump.LT.table%nrow)

            IF (jump.GE.table%nrow) THEN
               IF (table%nrow*2_ppm_kind_int64.GE.INT(ppm_big_i-1,ppm_kind_int64)) THEN
                  fail("hashtable with more than 2^31-1 elements will fail",ppm_error=ppm_error_fatal)
               ENDIF

               nsize2=INT(table%nrow)*2

               CALL table%destroy(info)
               or_fail("table%destroy")

               CALL table%create(nsize2,info)
               or_fail("table%create")

               GOTO 5555
            ENDIF
         ENDDO !i=1,nsize

         DEALLOCATE(keys_tmp,borders_pos_tmp,STAT=info)
         or_fail_dealloc("keys_tmp & borders_pos_tmp")

         !---------------------------------------------------------------------
         !  Deallocate ranklist array.
         !---------------------------------------------------------------------
         CALL ppm_alloc(ranklist,[0],ppm_param_dealloc,info)
         or_fail_dealloc('htable_keys',ppm_error=ppm_error_fatal)
      END SELECT

      IF (PRESENT(key_)) THEN
         IF (PRESENT(value1_).AND.PRESENT(value2_)) THEN
            CALL table%MCMCResults_hash_insert(key_,value1_,value2_,info)
         ENDIF
      ENDIF

      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      END SUBROUTINE MCMCResults_grow_htable

      FUNCTION MCMCResults_hash_size(table) RESULT(value)
      IMPLICIT NONE
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      CLASS(ppm_rc_MCMCResultshtable) :: table
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
      END FUNCTION MCMCResults_hash_size

#ifdef __F2003
#else
      ! Note - This code makes an assumption about how your machine behaves -

      ! 1. sizeof(INTEGER(ppm_kind_int64)) == 8

      ! Limitation
      ! It will not produce the same results on little-endian and big-endian machines.
      FUNCTION HashFunc(key,seed,tablenrow) RESULT(hash_val)
      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64), INTENT(IN   ) :: key
      !!! Input key
      INTEGER(ppm_kind_int64), INTENT(IN   ) :: seed
      !!! Seed to be used for mixing
      INTEGER(ppm_kind_int64), INTENT(IN   ) :: tablenrow
      !!! Number of table rows
      INTEGER(ppm_kind_int64)                :: hash_val
      !!! Result of the hash function
      !---------------------------------------------------------------------
      !  Local variables and parameters
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64), PARAMETER :: two32=4294967296_ppm_kind_int64
      INTEGER(ppm_kind_int64), PARAMETER :: m = 1540483477_ppm_kind_int64 !0x5bd1e995
      INTEGER(ppm_kind_int64)            :: h,k

      k=IBITS(key,0,32)

      ! This multiplication is safe.
      ! For the biggest 32 bit integer times m, it does not overflow.
      k=MOD(k*m,two32)
      k=IEOR(k,ISHFT(k,-24))
      k=MOD(k*m,two32)

      ! Initialize the hash to a 'random' value
      ! 8 is for 8 Bytes key len
      h=IEOR(seed,8_ppm_kind_int64)
      h=MOD(h*m,two32)
      h=IEOR(h,k)

      k=IBITS(key,32,32)

      k=MOD(k*m,two32)
      k=IEOR(k,ISHFT(k,-24))
      k=MOD(k*m,two32)

      h=MOD(h*m,two32)
      h=IEOR(h,k)

      ! Do a few final mixes of the hash to ensure the last few
      ! bytes are well-incorporated.
      h=IEOR(h,ISHFT(h,-13))
      h=MOD(h*m,two32)
      h=IEOR(h,ISHFT(h,-15))

      hash_val=IAND(h,tablenrow-1_ppm_kind_int64)
      RETURN
      END FUNCTION HashFunc

      FUNCTION HashFunc_XY(X32,Y32,seed,tablenrow,XY) RESULT(hash_val)
      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      INTEGER,                 INTENT(IN   ) :: X32
      !!! Input key
      INTEGER,                 INTENT(IN   ) :: Y32
      !!! Input key
      INTEGER(ppm_kind_int64), INTENT(IN   ) :: seed
      !!! Seed to be used for mixing
      INTEGER(ppm_kind_int64), INTENT(IN   ) :: tablenrow
      !!! Number of table rows
      INTEGER(ppm_kind_int64), INTENT(  OUT) :: XY
      !!! Input key
      INTEGER(ppm_kind_int64)                :: hash_val
      !!! Result of the hash function
      !---------------------------------------------------------------------
      !  Local variables and parameters
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64), PARAMETER :: two32=4294967296_ppm_kind_int64
      INTEGER(ppm_kind_int64), PARAMETER :: m = 1540483477_ppm_kind_int64 !0x5bd1e995
      INTEGER(ppm_kind_int64)            :: h,k
      INTEGER(ppm_kind_int64)            :: X64,Y64

      X64=INT(X32,ppm_kind_int64)
      Y64=INT(Y32,ppm_kind_int64)
      ! this is fine as *X32 and *Y32 are positive numbers and less than 2^16 (<=65535)
      XY=IOR(ISHFT(Y64,16),X64)

      k=IBITS(XY,0,32)

      ! This multiplication is safe.
      ! For the biggest 32 bit integer times m, it does not overflow.
      k=MOD(k*m,two32)
      k=IEOR(k,ISHFT(k,-24))
      k=MOD(k*m,two32)

      ! Initialize the hash to a 'random' value
      ! 8 is for 8 Bytes key len
      h=IEOR(seed,8_ppm_kind_int64)
      h=MOD(h*m,two32)
      h=IEOR(h,k)

      k=IBITS(XY,32,32)

      k=MOD(k*m,two32)
      k=IEOR(k,ISHFT(k,-24))
      k=MOD(k*m,two32)

      h=MOD(h*m,two32)
      h=IEOR(h,k)

      ! Do a few final mixes of the hash to ensure the last few
      ! bytes are well-incorporated.
      h=IEOR(h,ISHFT(h,-13))
      h=MOD(h*m,two32)
      h=IEOR(h,ISHFT(h,-15))

      hash_val=IAND(h,tablenrow-1_ppm_kind_int64)
      RETURN
      END FUNCTION HashFunc_XY

      FUNCTION HashFunc_XYLabel(X32,Y32,clabel,seed,tablenrow,XYLabel) RESULT(hash_val)
      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      INTEGER,                 INTENT(IN   ) :: X32
      !!! Input key
      INTEGER,                 INTENT(IN   ) :: Y32
      !!! Input key
      INTEGER,                 INTENT(IN   ) :: clabel
      !!! Input key
      INTEGER(ppm_kind_int64), INTENT(IN   ) :: seed
      !!! Seed to be used for mixing
      INTEGER(ppm_kind_int64), INTENT(IN   ) :: tablenrow
      !!! Number of table rows
      INTEGER(ppm_kind_int64), INTENT(  OUT) :: XYLabel
      !!! Input key
      INTEGER(ppm_kind_int64)                :: hash_val
      !!! Result of the hash function
      !---------------------------------------------------------------------
      !  Local variables and parameters
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64), PARAMETER :: two32=4294967296_ppm_kind_int64
      INTEGER(ppm_kind_int64), PARAMETER :: m = 1540483477_ppm_kind_int64 !0x5bd1e995
      INTEGER(ppm_kind_int64)            :: h,k
      INTEGER(ppm_kind_int64)            :: X64,Y64

      X64=INT(X32,ppm_kind_int64)
      Y64=INT(Y32,ppm_kind_int64)
      ! this is fine as *X32 and *Y32 are positive numbers and less than 2^16 (<=65535)
      XYLabel=INT(clabel,ppm_kind_int64)
      XYLabel=IOR(ISHFT(XYLabel,32),IOR(ISHFT(Y64,16),X64))

      k=IBITS(XYLabel,0,32)

      ! This multiplication is safe.
      ! For the biggest 32 bit integer times m, it does not overflow.
      k=MOD(k*m,two32)
      k=IEOR(k,ISHFT(k,-24))
      k=MOD(k*m,two32)

      ! Initialize the hash to a 'random' value
      ! 8 is for 8 Bytes key len
      h=IEOR(seed,8_ppm_kind_int64)
      h=MOD(h*m,two32)
      h=IEOR(h,k)

      k=IBITS(XYLabel,32,32)

      k=MOD(k*m,two32)
      k=IEOR(k,ISHFT(k,-24))
      k=MOD(k*m,two32)

      h=MOD(h*m,two32)
      h=IEOR(h,k)

      ! Do a few final mixes of the hash to ensure the last few
      ! bytes are well-incorporated.
      h=IEOR(h,ISHFT(h,-13))
      h=MOD(h*m,two32)
      h=IEOR(h,ISHFT(h,-15))

      hash_val=IAND(h,tablenrow-1_ppm_kind_int64)
      RETURN
      END FUNCTION HashFunc_XYLabel

      FUNCTION HashFunc_XYZ(X32,Y32,Z32,seed,tablenrow,XYZ) RESULT(hash_val)
      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      INTEGER,                 INTENT(IN   ) :: X32
      !!! Input key
      INTEGER,                 INTENT(IN   ) :: Y32
      !!! Input key
      INTEGER,                 INTENT(IN   ) :: Z32
      !!! Input key
      INTEGER(ppm_kind_int64), INTENT(IN   ) :: seed
      !!! Seed to be used for mixing
      INTEGER(ppm_kind_int64), INTENT(IN   ) :: tablenrow
      !!! Number of table rows
      INTEGER(ppm_kind_int64), INTENT(  OUT) :: XYZ
      !!! Input key
      INTEGER(ppm_kind_int64)                :: hash_val
      !!! Result of the hash function
      !---------------------------------------------------------------------
      !  Local variables and parameters
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64), PARAMETER :: two32=4294967296_ppm_kind_int64
      INTEGER(ppm_kind_int64), PARAMETER :: m = 1540483477_ppm_kind_int64 !0x5bd1e995
      INTEGER(ppm_kind_int64)            :: h,k
      INTEGER(ppm_kind_int64)            :: X64,Y64

      X64=INT(X32,ppm_kind_int64)
      Y64=INT(Y32,ppm_kind_int64)

      ! this is fine as X32,Y32 & Z32 are positive numbers and less than 2^11 (<=2047)
      XYZ=INT(Z32,ppm_kind_int64)
      XYZ=IOR(ISHFT(XYZ,22),IOR(ISHFT(Y64,11),X64))

      k=IBITS(XYZ,0,32)

      ! This multiplication is safe.
      ! For the biggest 32 bit integer times m, it does not overflow.
      k=MOD(k*m,two32)
      k=IEOR(k,ISHFT(k,-24))
      k=MOD(k*m,two32)

      ! Initialize the hash to a 'random' value
      ! 8 is for 8 Bytes key len
      h=IEOR(seed,8_ppm_kind_int64)
      h=MOD(h*m,two32)
      h=IEOR(h,k)

      k=IBITS(XYZ,32,32)

      k=MOD(k*m,two32)
      k=IEOR(k,ISHFT(k,-24))
      k=MOD(k*m,two32)

      h=MOD(h*m,two32)
      h=IEOR(h,k)

      ! Do a few final mixes of the hash to ensure the last few
      ! bytes are well-incorporated.
      h=IEOR(h,ISHFT(h,-13))
      h=MOD(h*m,two32)
      h=IEOR(h,ISHFT(h,-15))

      hash_val=IAND(h,tablenrow-1_ppm_kind_int64)
      RETURN
      END FUNCTION HashFunc_XYZ

      FUNCTION HashFunc_XYZLabel(X32,Y32,Z32,clabel,seed,tablenrow,XYZLabel) RESULT(hash_val)
      IMPLICIT NONE

      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
      INTEGER,                 INTENT(IN   ) :: X32
      !!! Input key
      INTEGER,                 INTENT(IN   ) :: Y32
      !!! Input key
      INTEGER,                 INTENT(IN   ) :: Z32
      !!! Input key
      INTEGER,                 INTENT(IN   ) :: clabel
      !!! Input key
      INTEGER(ppm_kind_int64), INTENT(IN   ) :: seed
      !!! Seed to be used for mixing
      INTEGER(ppm_kind_int64), INTENT(IN   ) :: tablenrow
      !!! Number of table rows
      INTEGER(ppm_kind_int64), INTENT(  OUT) :: XYZLabel
      !!! Input key
      INTEGER(ppm_kind_int64)                :: hash_val
      !!! Result of the hash function
      !---------------------------------------------------------------------
      !  Local variables and parameters
      !---------------------------------------------------------------------
      INTEGER(ppm_kind_int64), PARAMETER :: two32=4294967296_ppm_kind_int64
      INTEGER(ppm_kind_int64), PARAMETER :: m = 1540483477_ppm_kind_int64 !0x5bd1e995
      INTEGER(ppm_kind_int64)            :: h,k
      INTEGER(ppm_kind_int64)            :: X64,Y64,Z64

      X64=INT(X32,ppm_kind_int64)
      Y64=INT(Y32,ppm_kind_int64)
      Z64=INT(Z32,ppm_kind_int64)

      ! this is fine as X32,Y32 & Z are positive numbers and less than 2^11 (<=2047)
      XYZLabel=INT(clabel,ppm_kind_int64)
      XYZLabel=IOR(ISHFT(XYZLabel,11),Z64)
      XYZLabel=IOR(ISHFT(XYZLabel,22),IOR(ISHFT(Y64,11),X64))

      k=IBITS(XYZLabel,0,32)

      ! This multiplication is safe.
      ! For the biggest 32 bit integer times m, it does not overflow.
      k=MOD(k*m,two32)
      k=IEOR(k,ISHFT(k,-24))
      k=MOD(k*m,two32)

      ! Initialize the hash to a 'random' value
      ! 8 is for 8 Bytes key len
      h=IEOR(seed,8_ppm_kind_int64)
      h=MOD(h*m,two32)
      h=IEOR(h,k)

      k=IBITS(XYZLabel,32,32)

      k=MOD(k*m,two32)
      k=IEOR(k,ISHFT(k,-24))
      k=MOD(k*m,two32)

      h=MOD(h*m,two32)
      h=IEOR(h,k)

      ! Do a few final mixes of the hash to ensure the last few
      ! bytes are well-incorporated.
      h=IEOR(h,ISHFT(h,-13))
      h=MOD(h*m,two32)
      h=IEOR(h,ISHFT(h,-15))

      hash_val=IAND(h,tablenrow-1_ppm_kind_int64)
      RETURN
      END FUNCTION HashFunc_XYZLabel
#endif
#endif

#if   __DIME == __2D
      ELEMENTAL FUNCTION DTYPE(IndexHashFunctor64_label0)(index_1,index_2) RESULT(key)
#elif __DIME == __3D
      ELEMENTAL FUNCTION DTYPE(IndexHashFunctor64_label0)(index_1,index_2,index_3) RESULT(key)
#endif
        !!!
        IMPLICIT NONE
        !-------------------------------------------------------------------------
        !  Includes
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        INTEGER,  INTENT(IN   ) :: index_1
        INTEGER,  INTENT(IN   ) :: index_2
#if   __DIME == __3D
        INTEGER,  INTENT(IN   ) :: index_3
#endif
        INTEGER(ppm_kind_int64) :: key
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
        INTEGER(ppm_kind_int64) :: index1,index2
#if   __DIME == __3D
        INTEGER(ppm_kind_int64) :: index3
#endif

        index1=INT(index_1,ppm_kind_int64)
        index2=INT(index_2,ppm_kind_int64)
#if   __DIME == __2D
        key=IOR(ISHFT(index2,16),index1)
#elif __DIME == __3D
        index3=INT(index_3,ppm_kind_int64)

        key=IOR(ISHFT(index3,22),IOR(ISHFT(index2,11),index1))
#endif
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        RETURN
      END FUNCTION DTYPE(IndexHashFunctor64_label0)

#if   __DIME == __2D
      ELEMENTAL FUNCTION DTYPE(IndexHashFunctor64_label)(index_1,index_2,clabel) RESULT(key)
#elif __DIME == __3D
      ELEMENTAL FUNCTION DTYPE(IndexHashFunctor64_label)(index_1,index_2,index_3,clabel) RESULT(key)
#endif
        !!!
        IMPLICIT NONE
        !-------------------------------------------------------------------------
        !  Includes
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        INTEGER,  INTENT(IN   ) :: index_1
        INTEGER,  INTENT(IN   ) :: index_2
#if   __DIME == __3D
        INTEGER,  INTENT(IN   ) :: index_3
#endif
        INTEGER,  INTENT(IN   ) :: clabel
        INTEGER(ppm_kind_int64) :: key

        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
        INTEGER(ppm_kind_int64) :: index1,index2
#if   __DIME == __3D
        INTEGER(ppm_kind_int64) :: index3
#endif

        index1=INT(index_1,ppm_kind_int64)
        index2=INT(index_2,ppm_kind_int64)
#if   __DIME == __2D

        key=INT(clabel,ppm_kind_int64)
        key=IOR(ISHFT(key,32),IOR(ISHFT(index2,16),index1))
#elif __DIME == __3D
        index3=INT(index_3,ppm_kind_int64)

        ! Label is less than 2^31 and this shifting is safe here
        ! this is fine as X,Y & Z are positive numbers and less than 2^11 (<=2047)
        key=INT(clabel,ppm_kind_int64)
        key=IOR(ISHFT(key,11),index3)
        key=IOR(ISHFT(key,22),IOR(ISHFT(index2,11),index1))
#endif
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        RETURN
      END FUNCTION DTYPE(IndexHashFunctor64_label)

#if   __DIME == __2D
      SUBROUTINE DTYPE(HashIndexFunctor64_label0)(key,index_1,index_2)
#elif __DIME == __3D
      SUBROUTINE DTYPE(HashIndexFunctor64_label0)(key,index_1,index_2,index_3)
#endif
        !!!
        IMPLICIT NONE
        !-------------------------------------------------------------------------
        !  Includes
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        INTEGER(ppm_kind_int64), INTENT(IN   ) :: key
        INTEGER,                 INTENT(  OUT) :: index_1
        INTEGER,                 INTENT(  OUT) :: index_2
#if   __DIME == __3D
        INTEGER,                 INTENT(  OUT) :: index_3
#endif
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
#if   __DIME == __2D
        index_1=INT(IBITS(key, 0,16))
        index_2=INT(IBITS(key,16,16))
#elif __DIME == __3D
        index_1=INT(IBITS(key, 0,11))
        index_2=INT(IBITS(key,11,11))
        index_3=INT(IBITS(key,22,11))
#endif
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        RETURN
      END SUBROUTINE DTYPE(HashIndexFunctor64_label0)


#if   __DIME == __2D
      SUBROUTINE DTYPE(HashIndexFunctor64_label)(key,index_1,index_2,clabel)
#elif __DIME == __3D
      SUBROUTINE DTYPE(HashIndexFunctor64_label)(key,index_1,index_2,index_3,clabel)
#endif

        !!!
        IMPLICIT NONE
        !-------------------------------------------------------------------------
        !  Includes
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        INTEGER(ppm_kind_int64), INTENT(IN   ) :: key
        INTEGER,                 INTENT(  OUT) :: index_1
        INTEGER,                 INTENT(  OUT) :: index_2
#if   __DIME == __3D
        INTEGER,                 INTENT(  OUT) :: index_3
#endif
        INTEGER,                 INTENT(  OUT) :: clabel

        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
#if   __DIME == __2D
        index_1=INT(IBITS(key, 0,16))
        index_2=INT(IBITS(key,16,16))
        clabel =INT(IBITS(key,32,32))
#elif __DIME == __3D
        index_1=INT(IBITS(key, 0,11))
        index_2=INT(IBITS(key,11,11))
        index_3=INT(IBITS(key,22,11))
        clabel =INT(IBITS(key,33,31))
#endif
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        RETURN
      END SUBROUTINE DTYPE(HashIndexFunctor64_label)

      PURE FUNCTION DTYPE(IndexHashFunctor64)(index_) RESULT(key)
        !!!
        IMPLICIT NONE
        !-------------------------------------------------------------------------
        !  Includes
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        INTEGER, DIMENSION(:), INTENT(IN   ) :: index_
        INTEGER(ppm_kind_int64)              :: key
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
        INTEGER(ppm_kind_int64), DIMENSION(__DIME) :: index_64
        index_64=INT(index_,ppm_kind_int64)
#if   __DIME == __2D
        key=IOR(ISHFT(index_64(2),16),index_64(1))
#elif __DIME == __3D
        key=IOR(ISHFT(index_64(3),32),IOR(ISHFT(index_64(2),16),index_64(1)))
#endif
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        RETURN
      END FUNCTION DTYPE(IndexHashFunctor64)

#if   __DIME == __2D
      ELEMENTAL FUNCTION DTYPE(IndexHashFunctor64_)(index_1,index_2) RESULT(key)
#elif __DIME == __3D
      ELEMENTAL FUNCTION DTYPE(IndexHashFunctor64_)(index_1,index_2,index_3) RESULT(key)
#endif
        !!!
        IMPLICIT NONE
        !-------------------------------------------------------------------------
        !  Includes
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        INTEGER,  INTENT(IN   ) :: index_1
        INTEGER,  INTENT(IN   ) :: index_2
#if   __DIME == __3D
        INTEGER,  INTENT(IN   ) :: index_3
#endif
        INTEGER(ppm_kind_int64) :: key
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
        INTEGER(ppm_kind_int64) :: index1,index2
#if   __DIME == __2D
        index1=INT(index_1,ppm_kind_int64)
        index2=INT(index_2,ppm_kind_int64)
        key=IOR(ISHFT(index2,16),index1)
#elif __DIME == __3D
        INTEGER(ppm_kind_int64) :: index3
        index1=INT(index_1,ppm_kind_int64)
        index2=INT(index_2,ppm_kind_int64)
        index3=INT(index_3,ppm_kind_int64)
        key=IOR(ISHFT(index3,32),IOR(ISHFT(index2,16),index1))
#endif
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        RETURN
      END FUNCTION DTYPE(IndexHashFunctor64_)

      PURE FUNCTION DTYPE(IndexHashFunctor32)(index_) RESULT(key)
        !!!
        IMPLICIT NONE
        !-------------------------------------------------------------------------
        !  Includes
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        INTEGER, DIMENSION(:), INTENT(IN   ) :: index_
        INTEGER                              :: key
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
#if   __DIME == __2D
        key=IOR(ISHFT(index_(2),16),index_(1))
#elif __DIME == __3D
        key=IOR(ISHFT(IAND(index_(3),1023),22),IOR(ISHFT(IAND(index_(2),2047),11),IAND(index_(1),2047)))
#endif
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        RETURN
      END FUNCTION DTYPE(IndexHashFunctor32)

#if   __DIME == __2D
      ELEMENTAL FUNCTION DTYPE(IndexHashFunctor32_)(index_1,index_2) RESULT(key)
#elif __DIME == __3D
      ELEMENTAL FUNCTION DTYPE(IndexHashFunctor32_)(index_1,index_2,index_3) RESULT(key)
#endif
        !!!
        IMPLICIT NONE
        !-------------------------------------------------------------------------
        !  Includes
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        INTEGER, INTENT(IN   ) :: index_1
        INTEGER, INTENT(IN   ) :: index_2
#if   __DIME == __3D
        INTEGER, INTENT(IN   ) :: index_3
#endif
        INTEGER                :: key
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
#if   __DIME == __2D
        key=IOR(ISHFT(index_2,16),IAND(index_1,65535))
#elif __DIME == __3D
        key=IOR(ISHFT(IAND(index_3,1023),22),IOR(ISHFT(IAND(index_2,2047),11),IAND(index_1,2047)))
#endif
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        RETURN
      END FUNCTION DTYPE(IndexHashFunctor32_)

      SUBROUTINE DTYPE(HashIndexFunctor32)(key,index_)
        !!!
        IMPLICIT NONE
        !-------------------------------------------------------------------------
        !  Includes
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        INTEGER,               INTENT(IN   ) :: key
        INTEGER, DIMENSION(:), INTENT(  OUT) :: index_
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
#if   __DIME == __2D
        index_(1)=IBITS(key,0,16)
        index_(2)=IBITS(key,16,16)
#elif __DIME == __3D
        index_(1)=IBITS(key,0,11)
        index_(2)=IBITS(key,11,11)
        index_(3)=IBITS(key,22,10)
#endif
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        RETURN
      END SUBROUTINE DTYPE(HashIndexFunctor32)

#if   __DIME == __2D
      SUBROUTINE DTYPE(HashIndexFunctor32_)(key,index_1,index_2)
#elif __DIME == __3D
      SUBROUTINE DTYPE(HashIndexFunctor32_)(key,index_1,index_2,index_3)
#endif
        !!!
        IMPLICIT NONE
        !-------------------------------------------------------------------------
        !  Includes
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        INTEGER, INTENT(IN   ) :: key
        INTEGER, INTENT(  OUT) :: index_1
        INTEGER, INTENT(  OUT) :: index_2
#if   __DIME == __3D
        INTEGER, INTENT(  OUT) :: index_3
#endif
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
#if   __DIME == __2D
        index_1=IBITS(key,0,16)
        index_2=IBITS(key,16,16)
#elif __DIME == __3D
        index_1=IBITS(key,0,11)
        index_2=IBITS(key,11,11)
        index_3=IBITS(key,22,10)
#endif
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        RETURN
      END SUBROUTINE DTYPE(HashIndexFunctor32_)

      SUBROUTINE DTYPE(HashIndexFunctor64)(key,index_)
        !!!
        IMPLICIT NONE
        !-------------------------------------------------------------------------
        !  Includes
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        INTEGER(ppm_kind_int64), INTENT(IN   ) :: key
        INTEGER, DIMENSION(:),   INTENT(  OUT) :: index_
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
        index_(1)=INT(IBITS(key, 0,16))
        index_(2)=INT(IBITS(key,16,16))
#if   __DIME == __3D
        index_(3)=INT(IBITS(key,32,16))
#endif
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        RETURN
      END SUBROUTINE DTYPE(HashIndexFunctor64)

#if   __DIME == __2D
      SUBROUTINE DTYPE(HashIndexFunctor64_)(key,index_1,index_2)
#elif __DIME == __3D
      SUBROUTINE DTYPE(HashIndexFunctor64_)(key,index_1,index_2,index_3)
#endif
        !!!
        IMPLICIT NONE
        !-------------------------------------------------------------------------
        !  Includes
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        INTEGER(ppm_kind_int64), INTENT(IN   ) :: key
        INTEGER,                 INTENT(  OUT) :: index_1
        INTEGER,                 INTENT(  OUT) :: index_2
#if   __DIME == __3D
        INTEGER,                 INTENT(  OUT) :: index_3
#endif
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
        index_1=INT(IBITS(key, 0,16))
        index_2=INT(IBITS(key,16,16))
#if   __DIME == __3D
        index_3=INT(IBITS(key,32,16))
#endif
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        RETURN
      END SUBROUTINE DTYPE(HashIndexFunctor64_)



