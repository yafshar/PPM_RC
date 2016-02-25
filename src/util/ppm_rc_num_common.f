      !-------------------------------------------------------------------------
      !  Subroutine   :                    compute_num_common
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

      SUBROUTINE compute_num_common(a,b,c,n,info)

        IMPLICIT NONE

        INTEGER, DIMENSION(:), POINTER       :: a
        !first input
        INTEGER, DIMENSION(:), POINTER       :: b
        !second input
        INTEGER, DIMENSION(:), POINTER       :: c
        !common output
        INTEGER,               INTENT(  OUT) :: n
        !number of common
        INTEGER,               INTENT(  OUT) :: info

        REAL(ppm_kind_double) :: t0

        INTEGER :: i,iopt

        CHARACTER(LEN=ppm_char) :: caller='compute_num_common'

        LOGICAL, DIMENSION(SIZE(a)) :: MASK

        CALL substart(caller,t0,info)

        MASK=[(ANY(a(i).EQ.b),i=1,SIZE(a))]

        n = COUNT(MASK)

        iopt=ppm_param_alloc_grow
        CALL ppm_alloc(c,(/n/),iopt,info)
        or_fail_alloc("compute_num_common c")

        c=PACK(a,MASK)
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
      9999 CONTINUE
        CALL substop(caller,t0,info)
      END SUBROUTINE compute_num_common

