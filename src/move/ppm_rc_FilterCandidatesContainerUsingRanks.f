      !-------------------------------------------------------------------------
      !  Subroutine   :  ppm_rc_FilterCandidatesContainerUsingRanks
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
      !  Subroutine   :  ppm_rc_FilterCandidatesContainerUsingRanks
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
      SUBROUTINE DTYPE(ppm_rc_FilterCandidatesContainerUsingRanks)(info)

        !-------------------------------------------------------------------------
        !  Modules
        !-------------------------------------------------------------------------
        USE ppm_module_mpi
        USE ppm_module_util_qsort, ONLY : ppm_util_qsort

        USE ppm_rc_module_util, ONLY : label_exist
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

        REAL(ppm_kind_double)               :: t0
#ifdef __MPI
        REAL(MK)                            :: countsmin
        REAL(MK), DIMENSION(2)              :: counts2
        REAL(MK), DIMENSION(:), ALLOCATABLE :: counts
#endif

#if   __DIME == __2D
        INTEGER, CONTIGUOUS, DIMENSION(:,:),   POINTER :: DTYPE(wpp)
#elif __DIME == __3D
        INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER :: DTYPE(wpp)
#endif
        INTEGER, DIMENSION(:),     POINTER :: seedn
        INTEGER                            :: ipatch,ppart
        INTEGER                            :: nsize,nsize_,i
        INTEGER                            :: iter_id
#ifdef __MPI
        INTEGER                            :: request
#endif

        CHARACTER(LEN=ppm_char) :: caller="ppm_rc_FilterCandidatesContainerUsingRanks"
        !-------------------------------------------------------------------------
        !  Initialize
        !-------------------------------------------------------------------------
        CALL substart(caller,t0,info)

        ! Accept all the move, so nothing to be done here
        IF (AcceptedPointsFactor.GE.one) THEN
#ifdef __MPI
           IF (ppm_nproc.GT.1.AND..NOT.ConvergenceMASK(2)) THEN
              ipatch=mesh%subpatch%nb
              nsize_=SUM(Candidates(1:ipatch)%nb)

              SELECT CASE (nsize_)
              CASE (0)
                 counts2(1)=zero
                 counts2(2)=bigs+one

              CASE (1)
                 counts2(1)=one
                 counts2(2)=energya(Candidates_list(1))

              CASE DEFAULT
                 counts2=REAL(nsize_,MK)

              END SELECT

              ALLOCATE(counts(ppm_nproc*2),STAT=info)
              or_fail_alloc("counts")

              CALL MPI_Iallgather(counts2,2,ppm_mpi_kind,counts,2,ppm_mpi_kind,comm,request,info)
              or_fail_MPI("MPI_Iallgather")
           ELSE
              nsize_=2
           ENDIF
#endif
           GOTO 9999
        ENDIF

        ipatch=mesh%subpatch%nb
        nsize=SUM(Candidates(1:ipatch)%nb)
        nsize_=CEILING(nsize*AcceptedPointsFactor)

        !Sort the candidates according to their energy gradients, in
        !an ascending order
        CALL ppm_util_qsort(energya(Candidates_list),energyrank,info)
        or_fail("ppm_util_qsort")

#ifdef __MPI
        IF (ppm_nproc.GT.1.AND..NOT.ConvergenceMASK(2)) THEN
           ALLOCATE(counts(ppm_nproc*2),STAT=info)
           or_fail_alloc("counts")

           counts2=REAL(nsize_,MK)
           IF (nsize_.EQ.1) THEN
              counts2(2)=energya(Candidates_list(energyrank(1)))
           ENDIF

           CALL MPI_Iallgather(counts2,2,ppm_mpi_kind,counts,2,ppm_mpi_kind,comm,request,info)
           or_fail_MPI("MPI_Iallgather")
        ENDIF
#endif

        ! There might be several particles with the same energy level
        ! [NOTE]
        ! This is not available in the original algorithm
        DO i=nsize_+1,nsize
           IF (i.EQ.1) CYCLE
           IF (energya(Candidates_list(energyrank(i))).GT. &
           &   energya(Candidates_list(energyrank(i-1)))) EXIT
        ENDDO

        nsize=MERGE(i-1,0,i.GT.2)

        ALLOCATE(vSortedList(nsize),bufi(nsize),STAT=info)
        or_fail_alloc("vSortedList & bufi")

        !Arrange the sorted candidates from high to lowest energy
        !This has been done to help removing candidate from Candidates list
        !the point is in the removing step as it will change the pointer array
        !order by swaping the last element and the deleted one
        FORALL (i=1:nsize) vSortedList(i)=Candidates_list(energyrank(nsize-i+1))

        bufi=vSortedList
        !Sort the candidates according to their labels
        CALL ppm_util_qsort(bufi,info,nsize)
        or_fail("ppm_util_qsort")

        NULLIFY(DTYPE(wpp))

        sbpitr => mesh%subpatch%begin()
        ipatch=1
        i=0
        DO WHILE (ASSOCIATED(sbpitr))
           CALL sbpitr%get_field(pind,DTYPE(wpp),info)
           or_fail("Failed to get field i_wp data.")

           seed => Candidates(ipatch)%last()
           DO WHILE (ASSOCIATED(seed))
              i=i+1
              ppart=Candidates_list(i)

              IF (label_exist(ppart,bufi,nsize)) THEN
                 seed => Candidates(ipatch)%prev()
                 CYCLE
              ELSE IF (ppart.GT.Mpart) THEN
                 seedn => seed%first%getValue()

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
           ENDDO !ASSOCIATED(seed)
           !Removing extra particles from the list

           seed => Candidates(ipatch)%begin()
           DO WHILE (ASSOCIATED(seed))
              seedn => seed%first%getValue()

#if   __DIME == __2D
              ppart=DTYPE(wpp)(seedn(1),seedn(2))
#elif __DIME == __3D
              ppart=DTYPE(wpp)(seedn(1),seedn(2),seedn(3))
#endif

              iter_id=Candidates(ipatch)%iter_id

              !This way, I would also sort the Candidates
              !from the highest energy to the lowest
              DO i=iter_id,nsize
                 IF (ppart.EQ.vSortedList(i)) EXIT
              ENDDO

              IF (iter_id.EQ.i) THEN
                 !Means that the candidate is in the right place
                 seed => Candidates(ipatch)%next()
                 CYCLE
              ENDIF

              !just a swap between the current link and
              tmp_Candidates => Candidates(ipatch)%vec(i)%t
              Candidates(ipatch)%vec(i)%t => seed
              Candidates(ipatch)%vec(iter_id)%t => tmp_Candidates

              seed => Candidates(ipatch)%vec(iter_id)%t
           ENDDO !ASSOCIATED(seed)

           sbpitr => mesh%subpatch%next()
           ipatch=ipatch+1
        ENDDO !ASSOCIATED(sbpitr)
        !Candidates also has been sorted in descending order of energy

        DEALLOCATE(bufi,STAT=info)
        or_fail_dealloc("bufi")

        CALL MOVE_ALLOC(vSortedList,Candidates_list)
        !Destroy the vSortedList array and allocate it to Candidates_list
        !Now Candidates_list has been sorted in descending order of energy

        CALL ppm_alloc(energyrank,(/0/),ppm_param_dealloc,info)
        or_fail_dealloc("energyrank")

        NULLIFY(tmp_Candidates)

        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
      9999 CONTINUE

#ifdef __MPI
        IF (ppm_nproc.GT.1.AND..NOT.ConvergenceMASK(2)) THEN
           !wait for counts to be avilable
           CALL MPI_Wait(request,MPI_STATUS_IGNORE,info)
           or_fail_MPI("MPI_Wait")

           IF (ALL(counts(1:ppm_nproc*2:2).LE.oneplus)) THEN
              countsmin=bigs+one
              DO i=2,ppm_nproc*2,2
                 IF (counts(i).LE.countsmin.AND.counts(i).LT.zero) countsmin=counts(i)
              ENDDO

              IF (ABS(countsmin-counts2(2)).GT.lmyeps) THEN
                 nsize=0

                 DEALLOCATE(Candidates_list,STAT=info)
                 or_fail_dealloc("Candidates_list")

                 ALLOCATE(Candidates_list(nsize),STAT=info)
                 or_fail_alloc("Candidates_list")

                 sbpitr => mesh%subpatch%begin()
                 ipatch=1
                 DO WHILE (ASSOCIATED(sbpitr))
                    CALL Candidates(ipatch)%destroy(info)
                    or_fail("Candidates(ipatch)%destroy")
                    sbpitr => mesh%subpatch%next()
                    ipatch=ipatch+1
                 ENDDO !WHILE (ASSOCIATED(sbpitr))
              ENDIF
           ENDIF

           DEALLOCATE(counts,STAT=info)
           or_fail_dealloc("counts")
        ENDIF

        ConvergenceMASK(2)=nsize_.GT.1
#endif

        CALL substop(caller,t0,info)
        RETURN
      END SUBROUTINE DTYPE(ppm_rc_FilterCandidatesContainerUsingRanks)
