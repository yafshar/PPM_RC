      !-------------------------------------------------------------------------
      !  Subroutine   :                    ppm_rc_mem_usage
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

      SUBROUTINE ppm_rc_mem_usage(valueRSS,info)

!         USE ifport !if on intel compiler
        IMPLICIT NONE

        INTEGER, INTENT(INOUT) :: valueRSS
        !input
        INTEGER, INTENT(  OUT) :: info

        REAL(ppm_kind_double) :: t0

        INTEGER :: gpid

        CHARACTER(LEN=30)       :: dummy
        CHARACTER(LEN=ppm_char) :: rss_use_name,command
        CHARACTER(LEN=ppm_char) :: caller ='ppm_rc_mem_usage'

        CALL substart(caller,t0,info)

        gpid=GETPID()

        WRITE(dummy,'(I0)') gpid
        WRITE(rss_use_name,'(A,I0)')"/tmp/rss_use_",rank

        WRITE(command,'(4A)') "cat /proc/",TRIM(ADJUSTL(dummy)),"/status | grep RSS | grep -Eo '[0-9]{0,9}' > ",TRIM(ADJUSTL(rss_use_name))

        info=SYSTEM(command)

        OPEN(UNIT=777,FILE=TRIM(ADJUSTL(rss_use_name)))
        READ(777,FMT=*,END=7777,ERR=8888) valueRSS
        CLOSE(777,STATUS='DELETE')

        valueRSS=valueRSS/1024-valueRSS0
        !Convert memory use in KB to MB

        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        GOTO 9999
      7777 CLOSE(777,STATUS='DELETE')
        valueRSS=-1
        GOTO 9999
      8888 CLOSE(777,STATUS='DELETE')
        valueRSS=-2
      9999 CONTINUE
        CALL substop(caller,t0,info)
      END SUBROUTINE ppm_rc_mem_usage
