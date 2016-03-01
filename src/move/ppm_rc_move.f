      !-------------------------------------------------------------------------
      !  Subroutine   :                    ppm_rc_move
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
      !  Subroutine   :                    ppm_rc_move
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Samples the moves to accept and executes them.
      !                 Handles fusion and fission of regions.
      !
      !               Move all the points that are simple. Non simple points
      !               remain in the candidates list.
      !
      !               We first move all the FG-simple points. This we do because
      !               it happens  that points that are not simple at the first
      !               place get simple after the change of other points.
      !               The non-simple points will be treated in a separate loop
      !               afterward
      !
      !  Input        :
      !
      !  Input/output :
      !
      !  Output       : detot      (F) total energy gain from all moves
      !                 nchange    (I) total number of accepted moves
      !                 info       (I) return status. 0 on success.
      !
      !  Routines     :
      !
      !  Remarks      :
      !
      !  References   :
      !-------------------------------------------------------------------------
      LOGICAL FUNCTION DTYPE(ppm_rc_move)(info)

        !-------------------------------------------------------------------------
        !  Modules
        !-------------------------------------------------------------------------
        USE ppm_module_mpi
        USE ppm_module_util_qsort, ONLY : ppm_util_qsort
        USE ppm_module_util_time, ONLY : ppm_util_time

        USE ppm_rc_module_util, ONLY : label_exist
        IMPLICIT NONE

        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        INTEGER, INTENT(  OUT) :: info

        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
        CLASS(ppm_t_subpatch_),    POINTER :: sbpitr

        CLASS(ppm_t_part_prop_s_), POINTER :: prop

        TYPE(ppm_rc_list), POINTER :: seed
        TYPE(ppm_rc_list), POINTER :: vseed
        TYPE(ppm_rc_list), POINTER :: tseed
        TYPE(ppm_rc_link), POINTER :: seedlnk
        TYPE(ppm_rc_list_)         :: vLabelsToCheck

#if   __DIME == __2D
        REAL(MK), CONTIGUOUS, DIMENSION(:,:),   POINTER :: DTYPE(wpi)
#elif __DIME == __3D
        REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER :: DTYPE(wpi)
#endif
        REAL(MK),             DIMENSION(:,:),   POINTER :: xp
        REAL(ppm_kind_double)                           :: t0,tm0,tm1,tc0,tc1

        INTEGER,             DIMENSION(:), ALLOCATABLE :: old_ghost
        INTEGER,             DIMENSION(:),     POINTER :: seedn
        INTEGER,             DIMENSION(:),     POINTER :: vseedn
        INTEGER,             DIMENSION(:,:),   POINTER :: wpl
#if   __DIME == __2D
        INTEGER, CONTIGUOUS, DIMENSION(:,:),   POINTER :: DTYPE(wpl)
        INTEGER, CONTIGUOUS, DIMENSION(:,:),   POINTER :: DTYPE(wpp)
#elif __DIME == __3D
        INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER :: DTYPE(wpl)
        INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER :: DTYPE(wpp)
#endif
        INTEGER                                        :: ipatch
        INTEGER,             DIMENSION(__DIME)              :: ll,dd
        INTEGER,             DIMENSION(2)              :: ld
        INTEGER,             DIMENSION(:),     POINTER :: Nm
        INTEGER,             DIMENSION(:,:),   POINTER :: vFGTNvector
        INTEGER                                        :: vTopoNbItr
        INTEGER                                        :: nsize,iopt
        INTEGER                                        :: ppart,qpart
        INTEGER                                        :: ipart
        INTEGER,             DIMENSION(:), ALLOCATABLE :: vCheckedLabels
        INTEGER,             DIMENSION(:), POINTER     :: vCheckedLabelsranking
        INTEGER,             DIMENSION(:), POINTER     :: list_del_parts_ranking
        INTEGER                                        :: vCurrentLabel,vCandidateLabel
        INTEGER                                        :: vLabel,vLabel1,vLabel2,val
        INTEGER                                        :: vNSplits
        INTEGER                                        :: nneigh
        INTEGER                                        :: i,j,l,ii
#if   __DIME == __3D
        INTEGER                                        :: k
#endif
        INTEGER                                        :: iter_id,nb,il1
        INTEGER                                        :: Candidates_list_iter
#ifdef __MPI
        INTEGER                                        :: request
#endif

        CHARACTER(LEN=ppm_char) :: caller='ppm_rc_move'

        LOGICAL :: vChange
        LOGICAL :: vConvergence
        LOGICAL :: vSimple
        LOGICAL :: vValidPoint
        LOGICAL :: vSplit
        LOGICAL :: vMerge
        LOGICAL :: IsEnclosedByLabelFGConnectivity
        LOGICAL :: MASK
#ifdef __MPI
        LOGICAL :: MASKG
#endif

        !-------------------------------------------------------------------------
        !  Initialize
        !-------------------------------------------------------------------------
        CALL substart(caller,t0,info)

        IF (debug.GT.0) THEN
           CALL ppm_util_time(t0)
        ENDIF

        !-------------------------------------------------------------------------
        !  Initialize
        !-------------------------------------------------------------------------
        DTYPE(ppm_rc_move)=.FALSE.

        vConvergence=.TRUE.
        vChange=.TRUE.

        !-------------------------------------------------------------------------
        !  Create the hash table for global to local index
        !-------------------------------------------------------------------------
        CALL htableg%create(NpartNew+32,info)
        or_fail("create htableg")

        NpartNewtmp=0

        nsize=Mpart+NpartNew
        ALLOCATE(processeda(nsize),STAT=info)
        or_fail_alloc("processeda")

        processeda=1

        Candidates_list_iter=SIZE(Candidates_list)

        NULLIFY(DTYPE(wpi),DTYPE(wpl),DTYPE(wpp))

        !-------------------------------------------------------------------------
        ! We first move all the FG-simple points. This happens
        ! as points that are not simple at the first place get simple after
        ! the change of the other points. The non-simple points will be treated
        ! in a separate loop afterwards.
        !
        ! Firstly, we are trying to move all inner FG-simple points 2<= and <= Nm-1
        ! to avoid unnecassary ghost communication as much as possible
        !-------------------------------------------------------------------------
        sbpitr => mesh%subpatch%begin()
        ipatch=1
        DO WHILE (ASSOCIATED(sbpitr))
           MASK=vChange.AND.Candidates(ipatch)%nb.GT.0

           Nm => sbpitr%nnodes

           CALL sbpitr%get_field(image,DTYPE(wpi),info)
           or_fail("Failed to get field r_wp data.")

           CALL sbpitr%get_field(labels,DTYPE(wpl),info)
           or_fail("Failed to get field i_wp data.")

           CALL sbpitr%get_field(pind,DTYPE(wpp),info)
           or_fail("Failed to get field i_wp data.")

           inner_mask_loop: DO WHILE (MASK)
              vChange=.FALSE.

              !-------------------------------------------------------------------------
              !particle with the lowest energy
              !-------------------------------------------------------------------------
              seed => Candidates(ipatch)%last()
              inner_seed_loop: DO WHILE (ASSOCIATED(seed))
                 seedn => seed%first%getValue()

                 IF (ANY(seedn.LE.2.OR.seedn.GE.Nm-1)) THEN
                    seed => Candidates(ipatch)%prev()
                    CYCLE inner_seed_loop
                 ENDIF

#if   __DIME == __2D
                 ppart=DTYPE(wpp)(seedn(1),seedn(2))
#elif __DIME == __3D
                 ppart=DTYPE(wpp)(seedn(1),seedn(2),seedn(3))
#endif

                 IF (ppart.GT.nsize.OR.ppart.LE.0) THEN
                    CALL Candidates(ipatch)%remove(info)
                    or_fail("Candidates(ipatch)%remove")

                    DEALLOCATE(seed,STAT=info)
                    or_fail_dealloc("Failed to deallocate seed")

                    DO i=Candidates_list_iter,1,-1
                       IF (Candidates_list(i).EQ.ppart) THEN
                          Candidates_list(i)=-1
                          EXIT
                       ENDIF
                    ENDDO

                    seed => Candidates(ipatch)%at(Candidates(ipatch)%iter_id)
                    CYCLE inner_seed_loop
                 ELSE
                    !-------------------------------------------------------------------------
                    ! we will reuse the processed flag to indicate if a particle is
                    ! a seed
                    !-------------------------------------------------------------------------
                    processeda(ppart)=0
                 ENDIF

                 vFGTNvector => TopologicalNumberFunction%EvaluateAdjacentRegionsFGTNAtIndex(seedn,DTYPE(wpl))

                 vSimple=.TRUE.
                 DO vTopoNbItr=1,SIZE(vFGTNvector,DIM=2)
                    IF (vFGTNvector(2,vTopoNbItr).NE.1.OR.vFGTNvector(3,vTopoNbItr).NE.1) THEN
                       vSimple=.FALSE.
                       EXIT
                    ENDIF
                 ENDDO

                 DEALLOCATE(vFGTNvector,STAT=info)
                 or_fail_dealloc("vFGTNvector")
                 NULLIFY(vFGTNvector)

                 IF (vSimple) THEN
                    vChange     =.TRUE.
                    vConvergence=.FALSE.

                    CALL DTYPE(ppm_rc_ChangeContourParticleLabelToCandidateLabel) &
                    &   (ipatch,DTYPE(wpi),DTYPE(wpl),DTYPE(wpp),        &
                    &    seedn(1:__DIME),Nm,ppart,info)
                    or_fail("ppm_rc_ChangeContourParticleLabelToCandidateLabel")

                    CALL Candidates(ipatch)%remove(info)
                    or_fail("Candidates(ipatch)%remove")

                    DEALLOCATE(seed,STAT=info)
                    or_fail_dealloc("Failed to deallocate seed")

                    DO i=Candidates_list_iter,1,-1
                       IF (Candidates_list(i).EQ.ppart) THEN
                          Candidates_list(i)=-1
                          EXIT
                       ENDIF
                    ENDDO

                    seed => Candidates(ipatch)%at(Candidates(ipatch)%iter_id)
                    CYCLE inner_seed_loop
                 ENDIF !vSimple

                 seed => Candidates(ipatch)%prev()
              ENDDO inner_seed_loop

              MASK=vChange.AND.Candidates(ipatch)%nb.GT.0

           ENDDO inner_mask_loop !(MASK)

           sbpitr => mesh%subpatch%next()
           ipatch=ipatch+1
        ENDDO !(ASSOCIATED(sbpitr))

#ifdef __MPI
        ipatch=mesh%subpatch%nb
        MASK=ANY(Candidates(1:ipatch)%nb.GT.0)

        IF (debug.GT.0) THEN
           CALL ppm_util_time(tc0)
        ENDIF

        MASKG=.FALSE.

        CALL MPI_Iallreduce(MASK,MASKG,1,MPI_LOGICAL,MPI_LOR,comm,request,info)
        or_fail_MPI("MPI_Iallreduce")
#endif

        NULLIFY(DTYPE(wpi),DTYPE(wpl),DTYPE(wpp))

        !-------------------------------------------------------------------------
        !  destroy the memory which will not be used
        !-------------------------------------------------------------------------
        IF (ALLOCATED(energya)) THEN
           DEALLOCATE(energya,STAT=info)
           or_fail_dealloc("energya")
        ENDIF

        IF (ppm_nproc.GT.1) THEN
           ALLOCATE(old_ghost(ghost_size),STAT=info)
           or_fail_alloc("old_ghost")

           !-------------------------------------------------------------------------
           ! Ghost get
           !-------------------------------------------------------------------------
           ghostsize=1

           CALL mesh%map_ghost_get(info,ghostsize=ghostsize)
           or_fail("mesh%map_ghost_get")
        ENDIF

        vChange =.TRUE.

#ifdef __MPI
        !-------------------------------------------------------------------------
        ! wait to have MASKG
        !-------------------------------------------------------------------------
        CALL MPI_Wait(request,MPI_STATUS_IGNORE,info)
        or_fail_MPI("MPI_Wait")

        IF (debug.GT.0) THEN
           CALL ppm_util_time(tc1)
           tmove_SimpleComm=tmove_SimpleComm+tc1-tc0
        ENDIF
#endif

        !-------------------------------------------------------------------------
        ! We first move all the FG-simple points. This happens
        ! as points that are not simple at the first place get simple after
        ! the change of the other points. The non-simple points will be treated
        ! in a separate loop afterwards.
        !-------------------------------------------------------------------------
        sbpitr => mesh%subpatch%begin()
        ipatch=1
        DO WHILE (ASSOCIATED(sbpitr))
           Nm => sbpitr%nnodes

           CALL sbpitr%get_field(image,DTYPE(wpi),info)
           or_fail("Failed to get field r_wp data.")

           CALL sbpitr%get_field(labels,DTYPE(wpl),info)
           or_fail("Failed to get field i_wp data.")

           CALL sbpitr%get_field(pind,DTYPE(wpp),info)
           or_fail("Failed to get field i_wp data.")

           MASK=vChange.AND.Candidates(ipatch)%nb.GT.0

           mask_loop: DO WHILE (MASK)
              vChange=.FALSE.

              seed => Candidates(ipatch)%last()
              seed_loop: DO WHILE (ASSOCIATED(seed))
                 seedn => seed%first%getValue()

#if   __DIME == __2D
                 ppart=DTYPE(wpp)(seedn(1),seedn(2))
#elif __DIME == __3D
                 ppart=DTYPE(wpp)(seedn(1),seedn(2),seedn(3))
#endif

                 IF (ppart.GT.nsize.OR.ppart.LE.0) THEN
                    CALL Candidates(ipatch)%remove(info)
                    or_fail("Candidates(ipatch)%remove")

                    DEALLOCATE(seed,STAT=info)
                    or_fail_dealloc("Failed to deallocate seed")

                    DO i=Candidates_list_iter,1,-1
                       IF (Candidates_list(i).EQ.ppart) THEN
                          Candidates_list(i)=-1
                          EXIT
                       ENDIF
                    ENDDO

                    seed => Candidates(ipatch)%at(Candidates(ipatch)%iter_id)
                    CYCLE seed_loop
                 ELSE
                    !-------------------------------------------------------------------------
                    ! we will reuse the processed flag to indicate if a particle is
                    ! a seed
                    !-------------------------------------------------------------------------
                    processeda(ppart)=0
                 ENDIF

                 vFGTNvector => TopologicalNumberFunction%EvaluateAdjacentRegionsFGTNAtIndex(seedn,DTYPE(wpl))

                 vSimple=.TRUE.
                 DO vTopoNbItr=1,SIZE(vFGTNvector,DIM=2)
                    IF (vFGTNvector(2,vTopoNbItr).NE.1.OR.vFGTNvector(3,vTopoNbItr).NE.1) THEN
                       vSimple=.FALSE.
                       EXIT
                    ENDIF
                 ENDDO

                 DEALLOCATE(vFGTNvector,STAT=info)
                 or_fail_dealloc("vFGTNvector")
                 NULLIFY(vFGTNvector)

                 IF (vSimple) THEN
                    vChange     =.TRUE.
                    vConvergence=.FALSE.

                    CALL DTYPE(ppm_rc_ChangeContourParticleLabelToCandidateLabel) &
                    &   (ipatch,DTYPE(wpi),DTYPE(wpl),DTYPE(wpp),seedn(1:__DIME), &
                    &   Nm,ppart,info)
                    or_fail("ppm_rc_ChangeContourParticleLabelToCandidateLabel")

                    CALL Candidates(ipatch)%remove(info)
                    or_fail("Candidates(ipatch)%remove")

                    DEALLOCATE(seed,STAT=info)
                    or_fail_dealloc("Failed to deallocate seed")

                    DO i=Candidates_list_iter,1,-1
                       IF (Candidates_list(i).EQ.ppart) THEN
                          Candidates_list(i)=-1
                          EXIT
                       ENDIF
                    ENDDO

                    seed => Candidates(ipatch)%at(Candidates(ipatch)%iter_id)
                    CYCLE seed_loop
                 ENDIF !vSimple

                 seed => Candidates(ipatch)%prev()
              ENDDO seed_loop

              MASK=vChange.AND.Candidates(ipatch)%nb.GT.0

#ifdef __MPI
              IF (debug.GT.0) THEN
                 CALL ppm_util_time(tc0)
              ENDIF

              MASKG=.FALSE.

              CALL MPI_Iallreduce(MASK,MASKG,1,MPI_LOGICAL,MPI_LOR,comm,request,info)
              or_fail_MPI("MPI_Iallreduce")

              IF (ppm_nproc.GT.1) THEN
                 !-------------------------------------------------------------------------
                 !  ghost update for labels
                 ! Any label update at this part should be seen at the other processor
                 !-------------------------------------------------------------------------
                 CALL labels%map_ghost_push(mesh,info)
                 or_fail("labels%map_ghost_push")
                 CALL mesh%map_isend(info,sendrecv=.TRUE.)
                 or_fail("mesh%map_isend")

                 !-------------------------------------------------------------------------
                 ! Copy the ghost for comparison with the new values and impose the changes
                 ! to inside domain
                 !-------------------------------------------------------------------------
                 CALL DTYPE(ppm_rc_ghost_copy)(sbpitr,old_ghost,info)
                 or_fail("ppm_rc_ghost_copy")

                 CALL mesh%map_isend(info,sendrecv=.FALSE.)
                 or_fail("mesh%map_isend")
                 CALL labels%map_ghost_pop(mesh,info)
                 or_fail("labels%map_ghost_pop")

#include "./ChangeGhostContourParticleLabelToCandidateLabel.inc"
              ENDIF

              !-------------------------------------------------------------------------
              !Wait to have MASKG
              !-------------------------------------------------------------------------
              CALL MPI_Wait(request,MPI_STATUS_IGNORE,info)
              or_fail_MPI("MPI_Wait")

              IF (debug.GT.0) THEN
                 CALL ppm_util_time(tc1)
                 tmove_SimpleComm=tmove_SimpleComm+tc1-tc0
              ENDIF
#endif
!MPI
           ENDDO mask_loop !(MASK)

#ifdef __MPI
           IF (debug.GT.0) THEN
              CALL ppm_util_time(tc0)
           ENDIF
           !-------------------------------------------------------------------------
           ! Note: here MASK is false
           !
           ! processors whose MASK is .FALSE. should also participate in
           ! ghost push/send/pop
           !-------------------------------------------------------------------------
           DO WHILE (MASKG)
              MASKG=.FALSE.

              CALL MPI_Iallreduce(MASK,MASKG,1,MPI_LOGICAL,MPI_LOR,comm,request,info)
              or_fail_MPI("MPI_Iallreduce")

              IF (ppm_nproc.GT.1) THEN
                 !-------------------------------------------------------------------------
                 !  ghost update for labels
                 ! Any change at this part should be seen at the other processor
                 !-------------------------------------------------------------------------
                 CALL labels%map_ghost_push(mesh,info)
                 or_fail("labels%map_ghost_push")
                 CALL mesh%map_isend(info,sendrecv=.TRUE.)
                 or_fail("mesh%map_isend")

                 !-------------------------------------------------------------------------
                 ! Copy the ghost for comparison with the new values and impose the changes
                 ! to inside domain
                 !-------------------------------------------------------------------------
                 CALL DTYPE(ppm_rc_ghost_copy)(sbpitr,old_ghost(1:ghost_size),info)
                 or_fail("ppm_rc_ghost_copy")

                 CALL mesh%map_isend(info,sendrecv=.FALSE.)
                 or_fail("mesh%map_isend")
                 CALL labels%map_ghost_pop(mesh,info)
                 or_fail("labels%map_ghost_pop")

#include "./ChangeGhostContourParticleLabelToCandidateLabel.inc"
              ENDIF

              CALL MPI_Wait(request,MPI_STATUS_IGNORE,info)
              or_fail_MPI("MPI_Wait")
           ENDDO !MASKG

           IF (debug.GT.0) THEN
              CALL ppm_util_time(tc1)
              tmove_SimpleComm=tmove_SimpleComm+tc1-tc0
           ENDIF
#endif

           sbpitr => mesh%subpatch%next()
           ipatch=ipatch+1
        ENDDO !(ASSOCIATED(sbpitr))

        NULLIFY(DTYPE(wpi),DTYPE(wpl),DTYPE(wpp))

        IF (debug.GT.0) THEN
           CALL ppm_util_time(tm1)
           tmove_Simple=tmove_Simple+tm1-t0
           CALL ppm_util_time(tm0)
        ENDIF

        !-------------------------------------------------------------------------
        ! Now we know that all the points in the list are 'currently' not simple.
        ! We move them anyway (if topological constraints allow) but record
        ! (for every particle) where to relabel (using the seed set). Placing
        ! the seed is necessary for every particle to ensure relabeling even
        ! if a bunch of neighboring particles change. The seed will be ignored
        ! later on if the corresponding FG region is not present in the
        ! neighborhood anymore.
        ! TODO:
        ! The following code is dependent on the iteration order if splits/handles
        ! are not allowed. A solution would be to sort the candidates beforehand.
        ! This should be computationally not too expensive since we assume there
        ! are not many non-simple points.
        !-------------------------------------------------------------------------
        sbpitr => mesh%subpatch%begin()
        ipatch=1
        DO WHILE (ASSOCIATED(sbpitr))
           Nm => sbpitr%nnodes

           CALL sbpitr%get_field(image,DTYPE(wpi),info)
           or_fail("Failed to get field r_wp data.")

           CALL sbpitr%get_field(labels,DTYPE(wpl),info)
           or_fail("Failed to get field i_wp data.")

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

              IF (ppart.LE.0) THEN
                 CALL Candidates(ipatch)%remove(info)
                 or_fail("Candidates(ipatch)%remove")

                 DEALLOCATE(seed,STAT=info)
                 or_fail_dealloc("Failed to deallocate seed")

                 seed => Candidates(ipatch)%at(Candidates(ipatch)%iter_id)
                 CYCLE
              ENDIF

              vCurrentLabel  =labela(ppart)
              vCandidateLabel=candlabela(ppart)

              vValidPoint=.TRUE.

              vFGTNvector => TopologicalNumberFunction%EvaluateAdjacentRegionsFGTNAtIndex(seedn,DTYPE(wpl))

              !-------------------------------------------------------------------------
              ! Check for handles:
              !
              ! if the point was not disqualified already and we disallow
              ! introducing handles (not only self fusion!), we check if
              ! there is an introduction of a handle.
              !-------------------------------------------------------------------------
              IF (.NOT.AllowHandles) THEN
                 vValidPoint_loop: DO vTopoNbItr=1,SIZE(vFGTNvector,DIM=2)
                    IF (vFGTNvector(1,vTopoNbItr).EQ.vCandidateLabel) THEN
                       IF (vFGTNvector(2,vTopoNbItr).GT.1) THEN
                          vValidPoint=.FALSE.
                          EXIT vValidPoint_loop
                       ENDIF
                    ENDIF
                    !-------------------------------------------------------------------------
                    ! criterion to detect surface points or surface junctions.
                    ! or surface-curve junctions. Changing a surface-curve
                    ! junction will also cause a split.
                    !-------------------------------------------------------------------------
                    IF (vFGTNvector(3,vTopoNbItr).GT.1) THEN
                       vValidPoint=.FALSE.
                       EXIT vValidPoint_loop
                    ENDIF
                 ENDDO vValidPoint_loop
              ENDIF !(vValidPoint.AND..NOT.AllowHandles)


              !-------------------------------------------------------------------------
              ! Check for splits:
              !
              ! This we have to do either to forbid
              ! the change in topology or to register the seed point for
              ! relabelling.
              ! if the point was not disqualified already and we disallow
              ! splits, then we check if the 'old' label undergoes a split.
              !-------------------------------------------------------------------------
              IF (vValidPoint) THEN
                 ! - "allow introducing holes": T_FG(x, L = l') > 1
                 ! - "allow splits": T_FG > 2 && T_BG == 1
                 vSplit=.FALSE.
                 vTopoNbItr_loop: DO vTopoNbItr=1,SIZE(vFGTNvector,DIM=2)
                    IF (vFGTNvector(1,vTopoNbItr).EQ.vCurrentLabel) THEN
                       IF (vFGTNvector(2,vTopoNbItr).GT.1) THEN
                          vSplit=.TRUE.
                          EXIT vTopoNbItr_loop
                       ENDIF
                    ENDIF
                 ENDDO vTopoNbItr_loop

                 IF (vSplit) THEN
                    IF (AllowFission) THEN
                       DO i=1,FG_ConnectivityType%NumberOfNeighbors
                          ll=seedn(1:__DIME)+FG_ConnectivityType%NeighborsPoints(:,i)

                          IF (ANY(ll.LT.1.OR.ll.GT.Nm)) CYCLE
#if   __DIME == __2D
                          vLabel=ABS(DTYPE(wpl)(ll(1),ll(2)))
#elif __DIME == __3D
                          vLabel=ABS(DTYPE(wpl)(ll(1),ll(2),ll(3)))
#endif
                          IF (vLabel.EQ.vCurrentLabel) THEN
                             ALLOCATE(vseed,STAT=info)
                             or_fail_alloc("vseed")
#if   __DIME == __2D
                             CALL vseed%add(ll(1),ll(2),vLabel)
#elif __DIME == __3D
                             CALL vseed%add(ll(1),ll(2),ll(3),vLabel)
#endif
                             CALL m_Seeds(ipatch)%push(vseed,info)
                             or_fail("could not add new seed to the collection")

                             !-------------------------------------------------------------------------
                             ! At the position where we put the seed, inform the particle
                             ! that it has to inform its neighbor in case it moves (if there
                             ! is a particle at all at this spot; else we don't have a problem
                             ! because the label will not move at the spot and therefore the
                             ! seed will be effective).
                             !-------------------------------------------------------------------------
#if   __DIME == __2D
                             qpart=DTYPE(wpp)(ll(1),ll(2))
#elif __DIME == __3D
                             qpart=DTYPE(wpp)(ll(1),ll(2),ll(3))
#endif

                             IF (ANY(qpart.EQ.Candidates_list)) THEN
                                processeda(qpart)=1
                             ENDIF
                          ENDIF !(vLabel.EQ.vCurrentLabel)
                       ENDDO !i=1,FG_ConnectivityType%NumberOfNeighbors
                    ELSE !(.NOT.AllowFission)
                       !disallow the move.
                       vValidPoint=.FALSE.
                    ENDIF !(AllowFission)
                 ENDIF !(vSplit)
              ENDIF !(vValidPoint)

              DEALLOCATE(vFGTNvector,STAT=info)
              or_fail_dealloc("vFGTNvector")
              NULLIFY(vFGTNvector)

              !-------------------------------------------------------------------------
              ! If the move doesn't change topology or is allowed (and registered
              ! as seed) to change the topology, perform the move (in the
              ! second iteration; in the first iteration seed points need to
              ! be collected):
              !-------------------------------------------------------------------------
              IF (vValidPoint) THEN
                 vConvergence = .FALSE.

                 CALL DTYPE(ppm_rc_ChangeContourParticleLabelToCandidateLabel) &
                 &   (ipatch,DTYPE(wpi),DTYPE(wpl),DTYPE(wpp),seedn(1:__DIME), &
                 &   Nm,ppart,info)
                 or_fail("ppm_rc_ChangeContourParticleLabelToCandidateLabel")

                 !-------------------------------------------------------------------------
                 ! check if the point was a seed and if so, hand the scapegoat
                 ! to the other particle
                 !-------------------------------------------------------------------------
                 IF (processeda(ppart).EQ.1) THEN
                    DO i=1,FG_ConnectivityType%NumberOfNeighbors
                       ll=seedn(1:__DIME)+FG_ConnectivityType%NeighborsPoints(:,i)

                       IF (ANY(ll.LT.1.OR.ll.GT.Nm)) CYCLE

#if   __DIME == __2D
                       vLabel=ABS(DTYPE(wpl)(ll(1),ll(2)))
#elif __DIME == __3D
                       vLabel=ABS(DTYPE(wpl)(ll(1),ll(2),ll(3)))
#endif
                       IF (vLabel.EQ.vCurrentLabel) THEN
                          ALLOCATE(vseed,STAT=info)
                          or_fail_alloc("vseed")
#if   __DIME == __2D
                          CALL vseed%add(ll(1),ll(2),vLabel)
#elif __DIME == __3D
                          CALL vseed%add(ll(1),ll(2),ll(3),vLabel)
#endif

                          !-------------------------------------------------------------------------
                          ! At the position where we put the seed, inform the particle
                          ! that it has to inform its neighbor in case it moves (if there
                          ! is a particle at all at this spot; else we don't have a problem
                          ! because the label will not move at the spot and therefore the
                          ! seed will be effective).
                          !-------------------------------------------------------------------------
                          CALL m_Seeds(ipatch)%push(vseed,info)
                          or_fail("could not add new seed to the collection")

#if   __DIME == __2D
                          qpart=DTYPE(wpp)(ll(1),ll(2))
#elif __DIME == __3D
                          qpart=DTYPE(wpp)(ll(1),ll(2),ll(3))
#endif

                          IF (ANY(qpart.EQ.Candidates_list)) THEN
                             processeda(qpart)=1
                          ENDIF
                       ENDIF !(vLabel.EQ.vCurrentLabel)
                    ENDDO !i=1,FG_ConnectivityType%NumberOfNeighbors

                    vseed => m_Seeds(ipatch)%begin()
                    DO WHILE (ASSOCIATED(vseed))
                       vseedn => vseed%first%getValue()
                       IF (ALL(seedn(1:__DIME).EQ.vseedn(1:__DIME))) THEN
                          CALL m_Seeds(ipatch)%remove(info)
                          or_fail("m_Seeds(ipatch)%remove")

                          DEALLOCATE(vseed,STAT=info)
                          or_fail_dealloc("Failed to deallocate vseed")
                       ENDIF
                       vseed => m_Seeds(ipatch)%next()
                    ENDDO !ASSOCIATED(vseed)
                 ENDIF !(processeda(ppart).EQ.1)
              ENDIF !(vValidPoint)

              !-------------------------------------------------------------------------
              !TOCHECK
              !safely remove the last element
              !-------------------------------------------------------------------------
              CALL Candidates(ipatch)%remove(info)
              or_fail("Candidates(ipatch)%remove")

              DEALLOCATE(seed,STAT=info)
              or_fail_dealloc("Failed to deallocate seed")

              DO i=Candidates_list_iter,1,-1
                 IF (Candidates_list(i).EQ.ppart) THEN
                    Candidates_list(i)=-1
                    EXIT
                 ENDIF
              ENDDO

              seed => Candidates(ipatch)%at(Candidates(ipatch)%iter_id)
           ENDDO !(ASSOCIATED(seed))

           IF (ppm_nproc.GT.1) THEN
              !-------------------------------------------------------------------------
              ! Ghost update for labels
              ! Any label update at this part should be seen at the other processor
              !-------------------------------------------------------------------------
              CALL labels%map_ghost_push(mesh,info)
              or_fail("labels%map_ghost_push")
              CALL mesh%map_isend(info,sendrecv=.TRUE.)
              or_fail("mesh%map_send")

              !-------------------------------------------------------------------------
              ! Copy the ghost for comparison with the new values and impose the changes
              ! to inside domain
              !-------------------------------------------------------------------------
              CALL DTYPE(ppm_rc_ghost_copy)(sbpitr,old_ghost(1:ghost_size),info)
              or_fail("ppm_rc_ghost_copy")

              CALL mesh%map_isend(info,sendrecv=.FALSE.)
              or_fail("mesh%map_send")
              CALL labels%map_ghost_pop(mesh,info)
              or_fail("labels%map_ghost_pop")

#include "./ChangeGhostContourParticleLabelToCandidateLabel.inc"
           ENDIF

           sbpitr => mesh%subpatch%next()
           ipatch=ipatch+1
        ENDDO !(ASSOCIATED(sbpitr))

        NULLIFY(DTYPE(wpi),DTYPE(wpl),DTYPE(wpp))

        IF (debug.GT.0) THEN
           CALL ppm_util_time(tm1)
           tmove_NotSimple=tmove_NotSimple+tm1-tm0
        ENDIF

#ifdef __MPI
        ConvergenceMASK(1)=vConvergence
        CALL MPI_Iallreduce(MPI_IN_PLACE,ConvergenceMASK,2, &
        &    MPI_LOGICAL,MPI_LAND,comm,request,info)
        or_fail_MPI("MPI_Iallreduce")
#endif

        IF (ppm_nproc.GT.1) THEN
           DEALLOCATE(old_ghost,STAT=info)
           or_fail_dealloc("old_ghost")
        ENDIF

        !-------------------------------------------------------------------------
        ! Free memory which is not used anymore
        !-------------------------------------------------------------------------
        DEALLOCATE(Candidates_list,processeda,STAT=info)
        or_fail_dealloc("Candidates_list & processeda")

        !-------------------------------------------------------------------------
        ! Perform relabeling of the regions that did a split:
        !-------------------------------------------------------------------------
        vSplit=.FALSE.
        vMerge=.FALSE.

        sbpitr => mesh%subpatch%begin()
        ipatch=1
        DO WHILE (ASSOCIATED(sbpitr))
           CALL sbpitr%get_field(image,DTYPE(wpi),info)
           or_fail("Failed to get field r_wp data.")

           CALL sbpitr%get_field(labels,DTYPE(wpl),info)
           or_fail("Failed to get field i_wp data.")

           seed => m_Seeds(ipatch)%begin()
           DO WHILE (ASSOCIATED(seed))
              seedn => seed%first%getValue()

              !-------------------------------------------------------------------------
              ! Maybe it has been relabelled already.
              ! It is not necessary but should speed up
              !-------------------------------------------------------------------------
#if   __DIME == __2D
              vLabel1=ABS(DTYPE(wpl)(seedn(1),seedn(2)))
#elif __DIME == __3D
              vLabel1=ABS(DTYPE(wpl)(seedn(1),seedn(2),seedn(3)))
#endif

              IF (vLabel1.EQ.seedn(__DIME+1)) THEN
                 nneigh=0
                 DO i=1,FG_ConnectivityType%NumberOfNeighbors
                    ll=seedn(1:__DIME)+FG_ConnectivityType%NeighborsPoints(:,i)
#if   __DIME == __2D
                    IF (ABS(DTYPE(wpl)(ll(1),ll(2))).EQ.seedn(3)) THEN
#elif __DIME == __3D
                    IF (ABS(DTYPE(wpl)(ll(1),ll(2),ll(3))).EQ.seedn(4)) THEN
#endif
                       nneigh=nneigh+1
                       EXIT
                    ENDIF
                 ENDDO

                 IF (nneigh.EQ.0) THEN
                    !-------------------------------------------------------------------------
                    ! If this is a single point region, we need to remove it,
                    ! it is non signifacant, we still keep it in m_Seeds to remove the
                    ! particle
                    !-------------------------------------------------------------------------
#if   __DIME == __2D
                    DTYPE(wpl)(seedn(1),seedn(2))=0
#elif __DIME == __3D
                    DTYPE(wpl)(seedn(1),seedn(2),seedn(3))=0
#endif

                    CALL e_data%UpdateStatisticsWhenJump(DTYPE(wpi),seedn,vLabel1,0,info)
                    or_fail("e_data%UpdateStatisticsWhenJump")
                 ELSE
                    ALLOCATE(vseed,STAT=info)
                    or_fail_alloc("vseed")
#if   __DIME == __2D
                    CALL vseed%add(seedn(1),seedn(2))
#elif __DIME == __3D
                    CALL vseed%add(seedn(1),seedn(2),seedn(3))
#endif
                    CALL ppm_rc_seeds(ipatch)%push(vseed,info)
                    or_fail("could not add new seed to the collection")

                    !-------------------------------------------------------------------------
                    ! If it has no neighbor, it will cause single point region
                    ! creation, which is not desired, so we should know about this
                    !-------------------------------------------------------------------------
                    CALL m_Seeds(ipatch)%remove(info)
                    or_fail("m_Seeds(ipatch)%remove")

                    DEALLOCATE(seed,STAT=info)
                    or_fail_dealloc("Failed to deallocate seed")
                 ENDIF !(nneigh.EQ.0)
              ELSE
                 !-------------------------------------------------------------------------
                 ! If the label has already changed we need to remove the seed
                 !-------------------------------------------------------------------------
                 CALL m_Seeds(ipatch)%remove(info)
                 or_fail("m_Seeds(ipatch)%remove")

                 DEALLOCATE(seed,STAT=info)
                 or_fail_dealloc("Failed to deallocate seed")
              ENDIF !vLabel1.EQ.seedn(__DIME+1)
              seed => m_Seeds(ipatch)%next()
           ENDDO
           sbpitr => mesh%subpatch%next()
           ipatch=ipatch+1
        ENDDO

        NULLIFY(DTYPE(wpi),DTYPE(wpl))

        ipatch=mesh%subpatch%nb

        vNSplits = SUM(ppm_rc_seeds(1:ipatch)%nb)
        IF (vNSplits.GT.0) THEN
           vSplit=.TRUE.

           ALLOCATE(nlabels(vNSplits),STAT=info)
           or_fail_alloc("nlabels")

           DO i=1,vNSplits
              nlabels(i)=loc_label
              loc_label=loc_label-1
           ENDDO
        ENDIF

        !-------------------------------------------------------------------------
        ! Merge competing regions if they meet the merging criterion.
        !-------------------------------------------------------------------------
        IF (AllowFusion) THEN
           !-------------------------------------------------------------------------
           ! vLabelsToCheck contains all labels whose neighboring
           ! regions we still have to check.
           !-------------------------------------------------------------------------
           nsize=16

           ALLOCATE(vCheckedLabels(nsize),STAT=info)
           or_fail_alloc("vCheckedLabels")

           NULLIFY(vCheckedLabelsranking)

           nb=vNSplits

           sbpitr => mesh%subpatch%begin()
           ipatch=1
           DO WHILE (ASSOCIATED(sbpitr))
              CALL sbpitr%get_field(labels,DTYPE(wpl),info)
              or_fail("Failed to get field i_wp data.")

              !-------------------------------------------------------------------------
              ! If the candidate has already changed as a result of split and
              ! went to zero then no fusion will happen
              !-------------------------------------------------------------------------
              seed => CompetingRegions(ipatch)%last()
              DO WHILE (ASSOCIATED(seed))
                 seedn => seed%first%getValue()

#if   __DIME == __2D
                 IF (DTYPE(wpl)(seedn(1),seedn(2)).EQ.0.OR. &
                 &   DTYPE(wpl)(seedn(__DIME+2),seedn(__DIME+3)).EQ.0) THEN
#elif __DIME == __3D
                 IF (DTYPE(wpl)(seedn(1),seedn(2),seedn(3)).EQ.0.OR. &
                 &   DTYPE(wpl)(seedn(__DIME+2),seedn(__DIME+3),seedn(__DIME+4)).EQ.0) THEN
#endif
                    CALL CompetingRegions(ipatch)%remove(info)
                    or_fail("CompetingRegions(ipatch)%remove")

                    DEALLOCATE(seed,STAT=info)
                    or_fail_dealloc("Failed to deallocate seed")

                    seed => CompetingRegions(ipatch)%at(CompetingRegions(ipatch)%iter_id)
                    CYCLE
                 ENDIF
                 seed => CompetingRegions(ipatch)%prev()
              ENDDO !ASSOCIATED(seed)

              Nm => sbpitr%nnodes

              il1=0
              vCheckedLabels=0

              seed => CompetingRegions(ipatch)%begin()
              DO WHILE (ASSOCIATED(seed))
                 seedn => seed%first%getValue()

                 vLabel1 = seedn(__DIME+1)
                 vLabel2 = seedn(2*(__DIME+1))

                 IF (vLabel1.EQ.0.OR.vLabel2.EQ.0) THEN
                    stdout("Warning there is a zero label in CompetingRegions!!!")
                    seed => CompetingRegions(ipatch)%next()
                    CYCLE
                 ENDIF

                 IF (label_exist(vLabel1,vCheckedLabels,il1).OR. &
                 &   label_exist(vLabel2,vCheckedLabels,il1)) THEN
                    !-------------------------------------------------------------------------
                    ! do nothing, since this label was already in a fusion-chain
                    !-------------------------------------------------------------------------
                 ELSE
                    iter_id=CompetingRegions(ipatch)%iter_id

                    !-------------------------------------------------------------------------
                    ! Insert the first element (then iterate)
                    !-------------------------------------------------------------------------
                    CALL vLabelsToCheck%add(vLabel1)

                    il1=il1+1
                    IF (il1.GT.nsize) THEN
                       ALLOCATE(bufi(nsize+1),STAT=info)
                       or_fail_alloc("Tmp array allocation failed!")

                       FORALL (i=1:nsize) bufi(i)=vCheckedLabels(i)
                       nsize=nsize*2

                       CALL MOVE_ALLOC(bufi,vCheckedLabels)
                    ENDIF
                    !-------------------------------------------------------------------------
                    ! CheckedLabels contains all labels that are already in.
                    !-------------------------------------------------------------------------
                    vCheckedLabels(il1)=vLabel1

                    CALL ppm_util_qsort(vCheckedLabels,vCheckedLabelsranking,info,il1)
                    or_fail("ppm_util_qsort")

                    ALLOCATE(bufi(nsize),STAT=info)
                    or_fail_alloc("bufi")

                    DO i=1,il1
                       bufi(i)=vCheckedLabels(vCheckedLabelsranking(i))
                    ENDDO

                    CALL MOVE_ALLOC(bufi,vCheckedLabels)

                    ALLOCATE(tseed,STAT=info)
                    or_fail_alloc("tseed")
#if   __DIME == __2D
                    CALL tseed%add(seedn(1),seedn(2))
#elif __DIME == __3D
                    CALL tseed%add(seedn(1),seedn(2),seedn(3))
#endif
                    CALL ppm_rc_seeds(ipatch)%push(tseed,info)
                    or_fail("could not add new seed to the collection")

                    nb=nb+1
                    IF (nb.EQ.1) THEN
                       ALLOCATE(nlabels(1),STAT=info)
                       or_fail_alloc("nlabels")
                    ELSE
                       ld(1)=SIZE(nlabels,1)
                       IF (nb.GT.ld(1)) THEN
                          ALLOCATE(bufi(1:ld(1)*2),STAT=info)
                          or_fail_alloc("bufi")

                          FORALL (i=1:ld(1)) bufi(i)=nlabels(i)

                          CALL MOVE_ALLOC(bufi,nlabels)
                       ENDIF
                    ENDIF
                    nlabels(nb)=loc_label

                    DO WHILE (ASSOCIATED(vLabelsToCheck%first))
                       vLabel=vLabelsToCheck%first%getValue()
                       CALL vLabelsToCheck%remove()

                       !-------------------------------------------------------------------------
                       ! if one of the labels is involved in another fusion, we'll also relabel
                       ! this third/fourth... region
                       !-------------------------------------------------------------------------
                       vseed => CompetingRegions(ipatch)%begin()
                       DO WHILE (ASSOCIATED(vseed))
                          vseedn => vseed%first%getValue()

                          vLabel1 = vseedn(__DIME+1)
                          vLabel2 = vseedn(2*(__DIME+1))

                          IF (vLabel1.EQ.vLabel) THEN
!                              IF (.NOT.ANY(vLabel2.EQ.vCheckedLabels)) THEN
                             IF (.NOT.label_exist(vLabel2,vCheckedLabels,il1)) THEN
                                ALLOCATE(tseed,STAT=info)
                                or_fail_alloc("tseed")
#if   __DIME == __2D
                                CALL tseed%add(vseedn(4),vseedn(5))
#elif __DIME == __3D
                                CALL tseed%add(vseedn(5),vseedn(6),vseedn(7))
#endif
                                CALL ppm_rc_seeds(ipatch)%push(tseed,info)
                                or_fail("could not add new seed to the collection")

                                nb=nb+1
                                ld(1)=SIZE(nlabels,1)
                                IF (nb.GT.ld(1)) THEN
                                   ALLOCATE(bufi(ld(1)*2),STAT=info)
                                   or_fail_alloc("bufi")

                                   FORALL (i=1:ld(1)) bufi(i)=nlabels(i)

                                   CALL MOVE_ALLOC(bufi,nlabels)
                                ENDIF
                                nlabels(nb)=loc_label

                                il1=il1+1
                                IF (il1.GT.nsize) THEN
                                   ALLOCATE(bufi(nsize+1),STAT=info)
                                   or_fail_alloc("Tmp array allocation failed!")

                                   FORALL (i=1:nsize) bufi(i)=vCheckedLabels(i)
                                   nsize=nsize*2

                                   CALL MOVE_ALLOC(bufi,vCheckedLabels)
                                ENDIF
                                !-------------------------------------------------------------------------
                                ! CheckedLabels contains all labels that are already in.
                                !-------------------------------------------------------------------------
                                vCheckedLabels(il1)=vLabel2

                                CALL ppm_util_qsort(vCheckedLabels,vCheckedLabelsranking,info,il1)
                                or_fail("ppm_util_qsort")

                                ALLOCATE(bufi(nsize),STAT=info)
                                or_fail_alloc("bufi")

                                DO i=1,il1
                                   bufi(i)=vCheckedLabels(vCheckedLabelsranking(i))
                                ENDDO

                                CALL MOVE_ALLOC(bufi,vCheckedLabels)

                                CALL vLabelsToCheck%add(vLabel2)
                             ELSE IF (ANY(vseedn(__DIME+2:2*__DIME+1).LT.1.OR.vseedn(__DIME+2:2*__DIME+1).GT.Nm)) THEN
                                !-------------------------------------------------------------------------
                                !
                                !-------------------------------------------------------------------------
                                ALLOCATE(tseed,STAT=info)
                                or_fail_alloc("tseed")
#if   __DIME == __2D
                                CALL tseed%add(vseedn(4),vseedn(5))
#elif __DIME == __3D
                                CALL tseed%add(vseedn(5),vseedn(6),vseedn(7))
#endif
                                CALL ppm_rc_seeds(ipatch)%push(tseed,info)
                                or_fail("could not add new seed to the collection")

                                nb=nb+1
                                ld(1)=SIZE(nlabels,1)
                                IF (nb.GT.ld(1)) THEN
                                   ALLOCATE(bufi(1:ld(1)*2),STAT=info)
                                   or_fail_alloc("bufi")

                                   FORALL (i=1:ld(1)) bufi(i)=nlabels(i)

                                   CALL MOVE_ALLOC(bufi,nlabels)
                                ENDIF
                                nlabels(nb)=loc_label
                             ENDIF
                          ENDIF
                          IF (vLabel2.EQ.vLabel) THEN
!                              IF (.NOT.ANY(vLabel1.EQ.vCheckedLabels)) THEN
                             IF (.NOT.label_exist(vLabel1,vCheckedLabels,il1)) THEN
                                ALLOCATE(tseed,STAT=info)
                                or_fail_alloc("tseed")
#if   __DIME == __2D
                                CALL tseed%add(vseedn(1),vseedn(2))
#elif __DIME == __3D
                                CALL tseed%add(vseedn(1),vseedn(2),vseedn(3))
#endif
                                CALL ppm_rc_seeds(ipatch)%push(tseed,info)
                                or_fail("could not add new seed to the collection")

                                nb=nb+1
                                ld(1)=SIZE(nlabels,1)
                                IF (nb.GT.ld(1)) THEN
                                   ALLOCATE(bufi(1:ld(1)*2),STAT=info)
                                   or_fail_alloc("bufi")

                                   FORALL (i=1:ld(1)) bufi(i)=nlabels(i)

                                   CALL MOVE_ALLOC(bufi,nlabels)
                                ENDIF
                                nlabels(nb)=loc_label

                                il1=il1+1
                                IF (il1.GT.nsize) THEN
                                   ALLOCATE(bufi(nsize+1),STAT=info)
                                   or_fail_alloc("Tmp array allocation failed!")

                                   FORALL (i=1:nsize) bufi(i)=vCheckedLabels(i)
                                   nsize=nsize*2

                                   CALL MOVE_ALLOC(bufi,vCheckedLabels)
                                ENDIF
                                !-------------------------------------------------------------------------
                                ! CheckedLabels contains all labels that are already in.
                                !-------------------------------------------------------------------------
                                vCheckedLabels(il1)=vLabel1

                                CALL ppm_util_qsort(vCheckedLabels,vCheckedLabelsranking,info,il1)
                                or_fail("ppm_util_qsort")

                                ALLOCATE(bufi(nsize),STAT=info)
                                or_fail_alloc("bufi")

                                DO i=1,il1
                                   bufi(i)=vCheckedLabels(vCheckedLabelsranking(i))
                                ENDDO

                                CALL MOVE_ALLOC(bufi,vCheckedLabels)

                                CALL vLabelsToCheck%add(vLabel1)
                             ELSE IF (ANY(vseedn(1:__DIME).LT.1.OR.vseedn(1:__DIME).GT.Nm)) THEN
                                !-------------------------------------------------------------------------
                                !
                                !-------------------------------------------------------------------------
                                ALLOCATE(tseed,STAT=info)
                                or_fail_alloc("tseed")
#if   __DIME == __2D
                                CALL tseed%add(vseedn(1),vseedn(2))
#elif __DIME == __3D
                                CALL tseed%add(vseedn(1),vseedn(2),vseedn(3))
#endif
                                CALL ppm_rc_seeds(ipatch)%push(tseed,info)
                                or_fail("could not add new seed to the collection")

                                nb=nb+1
                                ld(1)=SIZE(nlabels,1)
                                IF (nb.GT.ld(1)) THEN
                                   ALLOCATE(bufi(1:ld(1)*2),STAT=info)
                                   or_fail_alloc("bufi")

                                   FORALL (i=1:ld(1)) bufi(i)=nlabels(i)

                                   CALL MOVE_ALLOC(bufi,nlabels)
                                ENDIF
                                nlabels(nb)=loc_label
                             ENDIF
                          ENDIF !vLabel2.EQ.vLabel
                          vseed => CompetingRegions(ipatch)%next()
                       ENDDO !ASSOCIATED(vseed)
                    ENDDO !ASSOCIATED(vLabelsToCheck%first)

                    CompetingRegions(ipatch)%iter_id=iter_id
                    loc_label=loc_label-1
                 ENDIF !ANY(vLabel1.EQ.vCheckedLabels).OR.ANY(vLabel2.EQ.vCheckedLabels)

                 seed => CompetingRegions(ipatch)%next()
              ENDDO !ASSOCIATED(seed)
              vMerge = .TRUE.
              sbpitr => mesh%subpatch%next()
              ipatch=ipatch+1
           ENDDO !ASSOCIATED(sbpitr)
           NULLIFY(DTYPE(wpl))

           DEALLOCATE(vCheckedLabels,STAT=info)
           or_fail_dealloc("vCheckedLabels")

           CALL ppm_alloc(vCheckedLabelsranking,ld,ppm_param_dealloc,info)
           or_fail_dealloc("vCheckedLabelsranking")

           CALL vLabelsToCheck%destroy()

           IF (ALLOCATED(nlabels)) THEN
              IF (SIZE(nlabels).GE.nb+1) THEN
                 FORALL (i=nb+1:SIZE(nlabels)) nlabels(i)=-1
              ENDIF
           ENDIF
        ENDIF !AllowFusion

        !-------------------------------------------------------------------------
        ! After, we checked the fusion and split, we know about the
        ! seed points which will fire and change the region label
        !
        ! Fire from the known seed points
        ! The fire will make some of the boundary particles a seed point
        ! to fire to the other processor through the ghost update
        !-------------------------------------------------------------------------
        CALL DTYPE(forestfire)(Part,mesh,.TRUE.,info)
        or_fail("forestfire")


        nsize=SIZE(list_del_parts)
        !-------------------------------------------------------------------------
        ! After firing from seed poionts, we got rid off the single point region,
        ! which was newly created, now we need to remove the particles of those
        ! regoins from the contourPoints
        !-------------------------------------------------------------------------
        sbpitr => mesh%subpatch%begin()
        ipatch=1
        DO WHILE (ASSOCIATED(sbpitr))
           CALL sbpitr%get_field(pind,DTYPE(wpp),info)
           or_fail("Failed to get field i_wp data.")

           seed => m_Seeds(ipatch)%begin()
           DO WHILE (ASSOCIATED(seed))
              seedn => seed%first%getValue()
#if   __DIME == __2D
              qpart=DTYPE(wpp)(seedn(1),seedn(2))
              DTYPE(wpp)(seedn(1),seedn(2))=0
#elif __DIME == __3D
              qpart=DTYPE(wpp)(seedn(1),seedn(2),seedn(3))
              DTYPE(wpp)(seedn(1),seedn(2),seedn(3))=0
#endif
              IF (qpart.GT.0.AND.qpart.LE.Npart) THEN
                 del_parts=del_parts+1
                 IF (del_parts.GT.nsize) THEN
                    ALLOCATE(bufi(nsize*2),STAT=info)
                    or_fail_alloc("Tmp array allocation failed!")

                    FORALL (i=1:nsize) bufi(i)=list_del_parts(i)

                    nsize=nsize*2

                    CALL MOVE_ALLOC(bufi,list_del_parts)
                 ENDIF
                 list_del_parts(del_parts)=qpart
              ELSE IF (qpart.GT.Mpart) THEN
                 val=htableg%search(qpart)
                 IF (val.NE.htable_null) THEN
                    del_parts=del_parts+1
                    IF (del_parts.GT.nsize) THEN
                       ALLOCATE(bufi(nsize*2),STAT=info)
                       or_fail_alloc("Tmp array allocation failed!")

                       FORALL (i=1:nsize) bufi(i)=list_del_parts(i)

                       nsize=nsize*2

                       CALL MOVE_ALLOC(bufi,list_del_parts)
                    ENDIF
                    list_del_parts(del_parts)=val
                    !-------------------------------------------------------------------------
                    ! val is the vec index in particle container which should
                    ! be removed
                    !-------------------------------------------------------------------------
                 ENDIF !val.NE.htable_null
              ENDIF !qpart.GT.0.AND.qpart.LE.Npart
              seed => m_Seeds(ipatch)%next()
           ENDDO !ASSOCIATED(seed)

           CALL m_Seeds(ipatch)%destroy(info)
           or_fail("m_Seeds(ipatch)%destroy")

           sbpitr => mesh%subpatch%next()
           ipatch=ipatch+1
        ENDDO !ASSOCIATED(sbpitr)

        NULLIFY(DTYPE(wpp))

        !-------------------------------------------------------------------------
        ! We need to remove inner contourPoints particles and
        ! the deleted ones from the particle container list
        !-------------------------------------------------------------------------
        nsize=SIZE(list_del_parts)


        sbpitr => mesh%subpatch%begin()
        ipatch=1
        DO WHILE (ASSOCIATED(sbpitr))
           CALL sbpitr%get_field(labels,DTYPE(wpl),info)
           or_fail("Failed to get field i_wp data.")

           CALL sbpitr%get_field(pind,DTYPE(wpp),info)
           or_fail("Failed to get field i_wp data.")

           seed => ppm_rc_seeds(ipatch)%begin()
           DO WHILE (ASSOCIATED(seed))
              seedlnk => seed%first
              DO WHILE (ASSOCIATED(seedlnk))
                 seedn => seedlnk%getValue()

#if   __DIME == __2D
                 vLabel=ABS(DTYPE(wpl)(seedn(1),seedn(2)))
#elif __DIME == __3D
                 vLabel=ABS(DTYPE(wpl)(seedn(1),seedn(2),seedn(3)))
#endif

                 IsEnclosedByLabelFGConnectivity=.TRUE.
                 DO i=1,FG_ConnectivityType%NumberOfNeighbors
                    ll=seedn+FG_ConnectivityType%NeighborsPoints(:,i)
#if   __DIME == __2D
                    IF (ABS(DTYPE(wpl)(ll(1),ll(2))).NE.vLabel) THEN
#elif __DIME == __3D
                    IF (ABS(DTYPE(wpl)(ll(1),ll(2),ll(3))).NE.vLabel) THEN
#endif
                       IsEnclosedByLabelFGConnectivity=.FALSE.
                       EXIT
                    ENDIF
                 ENDDO !i=1,FG_ConnectivityType%NumberOfNeighbors
                 IF (IsEnclosedByLabelFGConnectivity) THEN
#if   __DIME == __2D
                    DTYPE(wpl)(seedn(1),seedn(2))=vLabel
                    qpart=DTYPE(wpp)(seedn(1),seedn(2))
                    DTYPE(wpp)(seedn(1),seedn(2))=0
#elif __DIME == __3D
                    DTYPE(wpl)(seedn(1),seedn(2),seedn(3))=vLabel
                    qpart=DTYPE(wpp)(seedn(1),seedn(2),seedn(3))
                    DTYPE(wpp)(seedn(1),seedn(2),seedn(3))=0
#endif
                    IF (qpart.GT.0.AND.qpart.LE.Npart) THEN
                       del_parts=del_parts+1
                       IF (del_parts.GT.nsize) THEN
                          ALLOCATE(bufi(nsize*2),STAT=info)
                          or_fail_alloc("Tmp array allocation failed!")

                          FORALL (i=1:nsize) bufi(i)=list_del_parts(i)

                          nsize=nsize*2

                          CALL MOVE_ALLOC(bufi,list_del_parts)
                       ENDIF
                       list_del_parts(del_parts)=qpart
                    ELSE IF (qpart.GT.Mpart) THEN
                       val=htableg%search(qpart)
                       IF (val.NE.htable_null) THEN
                          del_parts=del_parts+1
                          IF (del_parts.GT.nsize) THEN
                             ALLOCATE(bufi(nsize*2),STAT=info)
                             or_fail_alloc("Tmp array allocation failed!")

                             FORALL (i=1:nsize) bufi(i)=list_del_parts(i)

                             nsize=nsize*2

                             CALL MOVE_ALLOC(bufi,list_del_parts)
                          ENDIF
                          list_del_parts(del_parts)=val
                          !-------------------------------------------------------------------------
                          ! val is the vec index in particle container which should
                          ! be removed
                          !-------------------------------------------------------------------------
                       ENDIF !val.NE.htable_null
                    ENDIF !qpart.GT.0.AND.qpart.LE.Npart
                 ENDIF !(IsEnclosedByLabelFGConnectivity)

                 seedlnk => seedlnk%nextLink()
              ENDDO !WHILE (ASSOCIATED(seedlnk))
              seed => ppm_rc_seeds(ipatch)%next()
           ENDDO
           sbpitr => mesh%subpatch%next()
           ipatch=ipatch+1
        ENDDO

        NULLIFY(DTYPE(wpl),DTYPE(wpp))

        CALL htableg%destroy(info)
        or_fail("destroying the hash table has Failed!")

        !-------------------------------------------------------------------------
        ! Sort the candidates according to their number
        !-------------------------------------------------------------------------
        NULLIFY(list_del_parts_ranking)
        CALL ppm_util_qsort(list_del_parts,list_del_parts_ranking,info,del_parts)
        or_fail("ppm_util_qsort")

        !TOCHECK for ipatch
        !TODO
        sbpitr => mesh%subpatch%begin()
        ipatch=1
        DO WHILE (ASSOCIATED(sbpitr))
           DO ipart=del_parts,1,-1
              iter_id=list_del_parts(list_del_parts_ranking(ipart))

              InnerContourContainer(ipatch)%iter_id=iter_id

              seed => InnerContourContainer(ipatch)%vec(iter_id)%t

              CALL InnerContourContainer(ipatch)%remove(info)
              or_fail("InnerContourContainer(ipatch)%remove")

              DEALLOCATE(seed,STAT=info)
              or_fail_dealloc("Failed to deallocate seed")
           ENDDO
           EXIT
           sbpitr => mesh%subpatch%next()
           ipatch=ipatch+1
        ENDDO

        !-------------------------------------------------------------------------
        ! We are removing the particles using the list_del_parts
        !-------------------------------------------------------------------------
        DEALLOCATE(list_del_parts,STAT=info)
        or_fail_dealloc("list_del_parts")

        CALL ppm_alloc(list_del_parts_ranking,ld,ppm_param_dealloc,info)
        or_fail_dealloc("list_del_parts_ranking")

        IF (debug.GT.0) THEN
           CALL ppm_util_time(tm0)
        ENDIF

        sbpitr => mesh%subpatch%begin()
        ipatch=1
        DO WHILE (ASSOCIATED(sbpitr))
           CALL ppm_rc_seeds(ipatch)%destroy(info)
           or_fail("ppm_rc_seeds(ipatch)%destroy")

           sbpitr => mesh%subpatch%next()
           ipatch=ipatch+1
        ENDDO
        !Freeing the memory which we do not need any more

        nsize=SUM(InnerContourContainer(:)%nb)

        !-------------------------------------------------------------------------
        ! Now we need to update the particles
        !-------------------------------------------------------------------------

        Part%Npart=nsize
        Part%Mpart=nsize

        Part%flags(ppm_part_ghosts)=.FALSE.

        IF (Mpart.NE.nsize) THEN
           iopt=ppm_param_alloc_fit
           ld(1)=__DIME
           ld(2)=nsize
           CALL ppm_alloc(Part%xp,ld,iopt,info)
           or_fail_alloc("Part%xp")
        ENDIF

        !-------------------------------------------------------------------------
        ! create special particle labels, it is important to send and recieve
        ! hot ghosts information and the old labels update
        !-------------------------------------------------------------------------
        CALL plabels%create(2,info,dtype=ppm_type_int,name="particle_labels")
        or_fail("Create plabels field failed!" )

        CALL plabels%discretize_on(Part,info)
        or_fail("plabels discretization on Particles failed!")

        NULLIFY(xp,wpl)
        CALL Part%get(plabels,wpl,info)
        or_fail("Part%get labels is failed!")

        CALL Part%get_xp(xp,info)
        or_fail("Part%get_xp")

        sbpitr => mesh%subpatch%begin()
        ipart=0
        ipatch=1
        DO WHILE (ASSOCIATED(sbpitr))
           Nm => sbpitr%nnodes

           CALL sbpitr%get_field(labels,DTYPE(wpl),info)
           or_fail("Failed to get field wpl data.")

           CALL sbpitr%get_field(pind,DTYPE(wpp),info)
           or_fail("Failed to get field wpp data.")

           !-------------------------------------------------------------------------
           ! clean the particle index list
           !-------------------------------------------------------------------------
           DTYPE(wpp)=0

           seed => InnerContourContainer(ipatch)%begin()
           DO WHILE (ASSOCIATED(seed))
              seedn => seed%first%getValue()

              ipart=ipart+1

              !TODO change to get_pos2-or 3d
              xp(1,ipart)=REAL(seedn(1)+sbpitr%istart(1)-2,MK)
              xp(2,ipart)=REAL(seedn(2)+sbpitr%istart(2)-2,MK)
#if   __DIME == __2D

              DTYPE(wpp)(seedn(1),seedn(2))=ipart
              wpl(1,ipart)=DTYPE(wpl)(seedn(1),seedn(2))
#elif __DIME == __3D
              xp(3,ipart)=REAL(seedn(3)+sbpitr%istart(3)-2,MK)
              !particle position

              DTYPE(wpp)(seedn(1),seedn(2),seedn(3))=ipart
              wpl(1,ipart)=DTYPE(wpl)(seedn(1),seedn(2),seedn(3))
#endif

              !-------------------------------------------------------------------------
              ! If there is any new seed
              !-------------------------------------------------------------------------
              IF (ALLOCATED(nlabels)) THEN
                 IF (ANY(seedn.EQ.1.OR.seedn.EQ.Nm)) THEN
                    !if the label of the ghost particle is equal to the newlabel
                    !then it is a hot particle which should fire in another
                    !processor
                    IF (label_exist(ABS(wpl(1,ipart)),nlabels,nb,.TRUE.)) THEN
                       wpl(1,ipart)=ABS(wpl(1,ipart))
                       ! The hot particles are positive in label
                       wpl(2,ipart)=0

                       neigh_loop: DO ii = 1,FG_ConnectivityType%NumberOfNeighbors
                          ll=seedn+FG_ConnectivityType%NeighborsPoints(:,ii)
                          IF (ANY(ll.EQ.1.OR.ll.EQ.Nm)) THEN
#if   __DIME == __2D
                             IF (DTYPE(wpl)(ll(1),ll(2)).EQ.wpl(1,ipart)) THEN
#elif __DIME == __3D
                             IF (DTYPE(wpl)(ll(1),ll(2),ll(3)).EQ.wpl(1,ipart)) THEN
#endif
                                IF (ANY(ll.LT.1.OR.ll.GT.Nm)) CYCLE neigh_loop
                                wpl(2,ipart)=ii
                                !-------------------------------------------------------------------------
                                ! the index of a label which would carry the
                                ! old particle label
                                !-------------------------------------------------------------------------
                                EXIT neigh_loop
                             ENDIF
                          ENDIF !ANY(ll.EQ.1.OR.ll.EQ.Nm)
                       ENDDO neigh_loop

                       IF (wpl(2,ipart).EQ.0) THEN
                          neigh_loop2: DO ii = 1,FG_ConnectivityType%NumberOfNeighbors
                             ll=seedn+FG_ConnectivityType%NeighborsPoints(:,ii)
                             IF (ANY(ll.LT.1.OR.ll.GT.Nm)) THEN
#if   __DIME == __2D
                                IF (ABS(DTYPE(wpl)(ll(1),ll(2))).EQ.wpl(1,ipart)) THEN
#elif __DIME == __3D
                                IF (ABS(DTYPE(wpl)(ll(1),ll(2),ll(3))).EQ.wpl(1,ipart)) THEN
#endif
                                   wpl(2,ipart)=ii
                                   !-------------------------------------------------------------------------
                                   ! the index of a label which would carry the
                                   ! old particle label
                                   !-------------------------------------------------------------------------
                                   EXIT neigh_loop2
                                ENDIF
                             ENDIF !ANY(ll.EQ.1.OR.ll.EQ.Nm)
                          ENDDO neigh_loop2
                       ENDIF

                       IF (wpl(2,ipart).EQ.0) THEN
                          !-------------------------------------------------------------------------
                          ! This seed should not participate in firing, as it
                          ! does not have a topologically connected region neighbor
                          !-------------------------------------------------------------------------
                          wpl(1,ipart)=-wpl(1,ipart)
                       ENDIF
                    ENDIF !ANY(ABS(wpl(1,ipart)).EQ.nlabels)
                 ENDIF !ANY(seedn.EQ.1.OR.seedn.EQ.Nm)
              ENDIF !(ALLOCATED(nlabels))

              seed => InnerContourContainer(ipatch)%next()
           ENDDO !ASSOCIATED(seed)

           partnm(ipatch)=InnerContourContainer(ipatch)%nb

           sbpitr => mesh%subpatch%next()
           ipatch=ipatch+1
        ENDDO

        NULLIFY(DTYPE(wpl),DTYPE(wpp))

        CALL Part%set_xp(xp,info)
        or_fail("Part%set_xp")

        !-------------------------------------------------------------------------
        !do not change any flag
        !-------------------------------------------------------------------------
        CALL Part%set(plabels,wpl,info)
        or_fail("Part%set for <plabels> field is failed")

        !-------------------------------------------------------------------------
        !We are using seed later
        !-------------------------------------------------------------------------
        NULLIFY(seed)

        !-------------------------------------------------------------------------
        !Number of particles
        !-------------------------------------------------------------------------
        Npart=Part%Npart

        IF (Mpart.NE.nsize) THEN
           CALL Part%realloc_prop(Part%gi,info)
           or_fail("Part%realloc_prop")
        ENDIF

        !-------------------------------------------------------------------------
        ! This is a hack to make it work
        !-------------------------------------------------------------------------
        Part%gi%flags(ppm_ppt_map_ghosts)=.TRUE.
        Part%gi%flags(ppm_ppt_map_parts) =.TRUE.
        Part%flags(ppm_part_areinside)   =.TRUE.
        Part%flags(ppm_part_partial)     =.TRUE.

        prop => part%props%begin()
        DO WHILE (ASSOCIATED(prop))
           prop%flags(ppm_ppt_partial)=.TRUE.
           prop => part%props%next()
        ENDDO

        CALL Part%comp_global_index(info)
        or_fail("Part%comp_global_index")

        !-------------------------------------------------------------------------
        ! Particle ghost mapping
        !-------------------------------------------------------------------------
        CALL Part%map_ghosts(info)
        or_fail("Part%map_ghosts")


        CALL Part%get(plabels,wpl,info,with_ghosts=.TRUE.,read_only=.TRUE.)
        or_fail("Part%get for <labels> field is failed!")

        !-------------------------------------------------------------------------
        ! get_xp does not change any flag
        !-------------------------------------------------------------------------
        CALL Part%get_xp(xp,info,with_ghosts=.TRUE.)
        or_fail("Part%get_xp")

        sbpitr => mesh%subpatch%begin()
        ipatch=1
        DO WHILE (ASSOCIATED(sbpitr))
           Nm => sbpitr%nnodes

           CALL sbpitr%get_field(labels,DTYPE(wpl),info)
           or_fail("Failed to get field wpl data.")

           CALL sbpitr%get_field(pind,DTYPE(wpp),info)
           or_fail("Failed to get field wpp data.")

           part_loop: DO ipart=Part%Npart+1,Part%Mpart
              i=NINT(xp(1,ipart))+2-sbpitr%istart(1)
              j=NINT(xp(2,ipart))+2-sbpitr%istart(2)
              ll(1)=i
              ll(2)=j

#if   __DIME == __2D
              IF (ANY(ll.LT.0.OR.ll.GT.Nm+1)) THEN
                 DTYPE(wpp)(i,j)=ipart
                 !-------------------------------------------------------------------------
                 ! We are only interested in the ghost particles which are
                 ! adjacent to the domain, the rest of the ghost particles are
                 ! not part of the InnerContourContainer list.
                 !-------------------------------------------------------------------------
                 CYCLE part_loop
              ENDIF
              !-------------------------------------------------------------------------
              ! TODO
              ! TOCHECK
              ! The implementation is not general
              !-------------------------------------------------------------------------
              IF (wpl(1,ipart).GT.0) THEN
                 dd=ll+FG_ConnectivityType%NeighborsPoints(:,wpl(2,ipart))
                 IF (DTYPE(wpl)(dd(1),dd(2)).NE.0) THEN
                    !I do not fire the background
                    !-------------------------------------------------------------------------
                    ! if the ghost particle label is positive it is a hot particle
                    ! and it should fire in the other processor domain
                    !-------------------------------------------------------------------------
                    ALLOCATE(tseed,STAT=info)
                    or_fail_alloc("tseed")

                    !-------------------------------------------------------------------------
                    ! coordinates of the hot particle and
                    ! neighbor index of the old label
                    !-------------------------------------------------------------------------
                    CALL tseed%add(i,j,wpl(2,ipart))

                    CALL ppm_rc_seeds(ipatch)%push(tseed,info)
                    or_fail("could not add new seed to the collection")

                    DTYPE(wpl)(i,j)=-wpl(1,ipart)
!               ELSE
!                  !If there is any new seed
!                  IF (ALLOCATED(nlabels)) THEN
!                     IF (ANY(ABS(DTYPE(wpl)(i,j)).EQ.nlabels)) THEN
!                     ELSE
!                        DTYPE(wpl)(i,j)=wpl(1,ipart)
!                     ENDIF
!                  ELSE
!                     DTYPE(wpl)(i,j)=wpl(1,ipart)
!                  ENDIF
                 ENDIF
              ENDIF

              DTYPE(wpp)(i,j)=ipart
#elif __DIME == __3D
              k=NINT(xp(3,ipart))+2-sbpitr%istart(3)
              ll(3)=k

              IF (ANY(ll.LT.0.OR.ll.GT.Nm+1)) THEN
                 DTYPE(wpp)(i,j,k)=ipart
                 !-------------------------------------------------------------------------
                 ! We are only interested in the ghost particles which are
                 ! adjacent to the domain, the rest of the ghost particles are
                 ! not part of the InnerContourContainer list.
                 !-------------------------------------------------------------------------
                 CYCLE part_loop
              ENDIF

              IF (wpl(1,ipart).GT.0) THEN
                 dd=ll+FG_ConnectivityType%NeighborsPoints(:,wpl(2,ipart))
                 IF (DTYPE(wpl)(dd(1),dd(2),dd(3)).NE.0) THEN
                    !I do not fire the background

                    ALLOCATE(tseed,STAT=info)
                    or_fail_alloc("tseed")

                    !-------------------------------------------------------------------------
                    ! coordinates of the hot particle and
                    ! neighbor index of the old particle
                    !-------------------------------------------------------------------------
                    CALL tseed%add(i,j,k,wpl(2,ipart))

                    CALL ppm_rc_seeds(ipatch)%push(tseed,info)
                    or_fail("could not add new seed to the collection")

                    DTYPE(wpl)(i,j,k)=-wpl(1,ipart)
!               ELSE
!                  !If there is any new seed
!                  IF (ALLOCATED(nlabels)) THEN
!                     IF (ANY(ABS(DTYPE(wpl)(i,j,k)).EQ.nlabels)) THEN
!                     ELSE
!                        DTYPE(wpl)(i,j,k)=wpl(1,ipart)
!                     ENDIF
!                  ELSE
!                     DTYPE(wpl)(i,j,k)=wpl(1,ipart)
!                  ENDIF
                 ENDIF
              ENDIF

              DTYPE(wpp)(i,j,k)=ipart
#endif

              ALLOCATE(seed,STAT=info)
              or_fail_alloc("seed")
              CALL seed%add(ll)
              CALL InnerContourContainer(ipatch)%push(seed,info)
              or_fail("could not add new seed to the collection")
           ENDDO part_loop

           sbpitr => mesh%subpatch%next()
           ipatch=ipatch+1
        ENDDO

        NULLIFY(DTYPE(wpl),DTYPE(wpp))

        CALL Part%set_xp(xp,info,read_only=.TRUE.)
        or_fail("Part%set_xp")

        CALL Part%set(plabels,wpl,info,read_only=.TRUE.)
        or_fail("Part%set for <wpl> pointer is failed!")

        IF (ALLOCATED(nlabels)) THEN
           DEALLOCATE(nlabels,STAT=info)
           or_fail_dealloc("nlabels")
        ENDIF

        !-------------------------------------------------------------------------
        ! Total number of particles (Npart+ghost particles)
        !-------------------------------------------------------------------------
        Mpart=Part%Mpart

        ALLOCATE(list_del_parts(Npart),STAT=info)
        or_fail_alloc("list_del_parts")

        del_parts=0

        IF (debug.GT.0) THEN
           CALL ppm_util_time(tm1)
           tmove_Part=tmove_Part+tm1-tm0
           CALL ppm_util_time(tm0)
        ENDIF

        !-------------------------------------------------------------------------
        ! After the forestfire the label is uptodate
        !-------------------------------------------------------------------------
        CALL DTYPE(forestfire)(Part,mesh,.FALSE.,info)
        or_fail("forestfire")

        IF (debug.GT.0) THEN
           CALL ppm_util_time(tm1)
           tmove_ghostfire=tmove_ghostfire+tm1-tm0
        ENDIF

        !-------------------------------------------------------------------------
        ! We need to free memory which will not be used later
        !-------------------------------------------------------------------------
        CALL plabels%destroy(info)
        or_fail("plabels%destroy")

#ifdef __MPI
        !-------------------------------------------------------------------------
        ! Wait till every processor has the vConvergence!
        !-------------------------------------------------------------------------
        CALL MPI_Wait(request,MPI_STATUS_IGNORE,info)
        or_fail_MPI("MPI_Wait")
#endif

        DTYPE(ppm_rc_move)=ConvergenceMASK(1)

        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
      9999 CONTINUE
        CALL substop(caller,t0,info)
        RETURN
      END FUNCTION DTYPE(ppm_rc_move)