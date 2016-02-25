      !-------------------------------------------------------------------------
      !  Subroutine   :                  ppm_rc_FilterCandidates
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
      !  Author           - y.afshar           June   2014
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Subroutine   :                  ppm_rc_FilterCandidates
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Loops through all the particle Candidates, and remove
      !                 particles which are illegal
      !
      !  Input        :
      !
      !  Input/output :
      !
      !  Output       : info       (I) return status. 0 on success.
      !
      !  Routines     :
      !
      !
      !  References   :
      !-------------------------------------------------------------------------
      SUBROUTINE DTYPE(ppm_rc_FilterCandidates)(info)

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
        INTEGER, INTENT(   OUT)  :: info
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
        CLASS(ppm_t_subpatch_), POINTER :: sbpitr

        TYPE(ppm_rc_list), POINTER :: seed

        REAL(ppm_kind_double) :: t0

#if   __DIME == __2D
        INTEGER, CONTIGUOUS, DIMENSION(:,:),   POINTER :: DTYPE(wpp)
#elif __DIME == __3D
        INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER :: DTYPE(wpp)
#endif
        INTEGER,             DIMENSION(:),     POINTER :: seedn
        INTEGER                                        :: ipatch,ppart,nsize,i

        CHARACTER(LEN=ppm_char) :: caller="ppm_rc_FilterCandidates"
        !-------------------------------------------------------------------------
        !  Initialize
        !-------------------------------------------------------------------------
        CALL substart(caller,t0,info)

        !According to the previous routines Candidates should not contain
        !any ghost particles
        ipatch=mesh%subpatch%nb
        nsize=SUM(Candidates(1:ipatch)%nb)
        ALLOCATE(bufi(nsize),STAT=info)
        or_fail_alloc("bufi")

        NULLIFY(DTYPE(wpp))

        sbpitr => mesh%subpatch%begin()
        ipatch=1
        i=0
        DO WHILE (ASSOCIATED(sbpitr))
           CALL sbpitr%get_field(pind,DTYPE(wpp),info)
           or_fail("Failed to get field i_wp data.")

           seed => Candidates(ipatch)%last()
           DO WHILE (ASSOCIATED(seed))
              seedn => seed%first%getValue()
#if   __DIME == __2D
              ppart=DTYPE(wpp)(seedn(1),seedn(2))
#elif __DIME == __3D
              ppart=DTYPE(wpp)(seedn(1),seedn(2),seedn(3))
#endif
              !Filter all candidates with the illegal indices
              !Filter candidates according to their energy
              SELECT CASE (accepteda(ppart))
              CASE (0)
                 !If it is a new particle, its index should be removed
                 IF (ppart.GT.Mpart) THEN
#if   __DIME == __2D
                    DTYPE(wpp)(seedn(1),seedn(2))=0
#elif __DIME == __3D
                    DTYPE(wpp)(seedn(1),seedn(2),seedn(3))=0
#endif
                 ENDIF

                 CALL Candidates(ipatch)%remove(info)
                 or_fail("Candidates(ipatch)%remove")

                 DEALLOCATE(seed,STAT=info)
                 or_fail_dealloc("Failed to deallocate seed")

                 seed => Candidates(ipatch)%at(Candidates(ipatch)%iter_id)
                 CYCLE

              CASE (1)
                 i=i+1
                 bufi(i)=ppart

              END SELECT

              seed => Candidates(ipatch)%prev()
           ENDDO !ASSOCIATED(seed)

           sbpitr => mesh%subpatch%next()
           ipatch=ipatch+1
        ENDDO !WHILE (ASSOCIATED(sbpitr))

        ipatch=mesh%subpatch%nb
        nsize=SUM(Candidates(1:ipatch)%nb)
        ALLOCATE(Candidates_list(nsize),SOURCE=bufi(1:nsize),STAT=info)
        or_fail_alloc("Candidates_list")

        DEALLOCATE(bufi,accepteda,STAT=info)
        or_fail_dealloc("bufi & accepteda")

        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
      9999 CONTINUE
        CALL substop(caller,t0,info)
        RETURN
      END SUBROUTINE DTYPE(ppm_rc_FilterCandidates)
