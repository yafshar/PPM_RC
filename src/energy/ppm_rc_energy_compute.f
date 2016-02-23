      !-------------------------------------------------------------------------
      !  Subroutine   :                  ppm_rc_energy_compute
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Loops through all the particles, finds the mother-
      !                 daughter relationships and stores them, and computes
      !                 the energy differences of label changes.
      !
      !  Input        :
      !
      !  Input/output :
      !
      !  Output       : info       (I) return status. 0 on success.
      !
      !  Routines     :
      !
      !  Remarks      : A daughter of a particle is any particle in the
      !                 6-neighborhood that belongs to a different region.
      !
      !  DIMENSION(3+17) :: ld
      !  for 3D case
      !  1     :: i
      !  2     :: j
      !  3     :: k
      !  4     :: candlabel
      !  5     :: referncecount
      !  6     :: nmothers
      !  7     :: ndaughters
      !  8-13  :: daughters
      ! 14     :: accepted
      !
      !
      !  DIMENSION(2+14) :: ld
      !  for 2D case
      !  1     :: i
      !  2     :: j
      !  3     :: candlabel
      !  4     :: referneccount
      !  5     :: nmothers
      !  6     :: ndaughters
      !  7-10  :: mothers
      ! 11-14  :: daughters
      ! 15     :: accepted
      !
      !  References   :
      !-------------------------------------------------------------------------
      !  MOSAIC Group
      !  Max Planck Institute of Molecular Cell Biology and Genetics
      !  Pfotenhauerstr. 108, 01307 Dresden, Germany
      !
      !  Author           - y.afshar           June   2014
      !-------------------------------------------------------------------------

      SUBROUTINE DTYPE(ppm_rc_energy_compute)(info)

        !-------------------------------------------------------------------------
        !  Modules
        !-------------------------------------------------------------------------
        USE ppm_module_mpi

        USE ppm_rc_module_linkedlist, ONLY : ppm_rc_list,InnerContourContainer, &
        &   Candidates,CompetingRegions
        USE ppm_rc_module_topologicalnumber, ONLY : FG_ConnectivityType
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
        TYPE(ppm_rc_list), POINTER :: seedi

#if   __DIME == __2D
        REAL(MK), CONTIGUOUS, DIMENSION(:,:),   POINTER :: DTYPE(wpi)
        !!!pointer to the intensity field data
#elif __DIME == __3D
        REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER :: DTYPE(wpi)
#endif
        REAL(ppm_kind_double)                           :: t0
        REAL(MK)                                        :: etemp,etemp1,etemp2
        REAL(MK)                                        :: TotalEnergy_,TotalEnergy

#if   __DIME == __2D
        INTEGER, CONTIGUOUS, DIMENSION(:,:),   POINTER :: DTYPE(wpl)
        !!!pointer to the label data
        INTEGER, CONTIGUOUS, DIMENSION(:,:),   POINTER :: DTYPE(wpp)
        !!!pointer to the particle index data
#elif __DIME == __3D
        INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER :: DTYPE(wpl)
        INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER :: DTYPE(wpp)
#endif
        INTEGER,             DIMENSION(:),     POINTER :: seedn
        INTEGER,             DIMENSION(:),     POINTER :: Nm
        INTEGER,             DIMENSION(__DIME)              :: ll
        INTEGER,             DIMENSION(9)              :: ierr
        INTEGER                                        :: i,ii,jj
        INTEGER                                        :: ipatch,nsize,nsizet
        INTEGER                                        :: ipart,ppart,qpart
        INTEGER                                        :: mylabel,newlabel

        CHARACTER(LEN=ppm_char) :: caller="ppm_rc_energy_compute"

        LOGICAL :: e_merge,l_e_merge
        LOGICAL :: IsEnclosedByLabelFGConnectivity
        !-------------------------------------------------------------------------
        !  Initialize
        !-------------------------------------------------------------------------
        CALL substart(caller,t0,info)

        IF (debug.GT.3) THEN
           SELECT TYPE (e_data)
           TYPE IS (E_PC)
              TotalEnergy_=zero
              NULLIFY(DTYPE(wpi),DTYPE(wpl))
              sbpitr => mesh%subpatch%begin()
              ipatch=1
              DO WHILE (ASSOCIATED(sbpitr))
                 Nm => sbpitr%nnodes

                 CALL sbpitr%get_field(image,DTYPE(wpi),info)
                 or_fail("Failed to get field r_wp data.")

                 CALL sbpitr%get_field(labels,DTYPE(wpl),info)
                 or_fail("Failed to get field i_wp data.")

                 TotalEnergy_=TotalEnergy_+ &
                 & e_data%CalculateTotalEnergy(DTYPE(wpi),DTYPE(wpl),Nm,info)

                 sbpitr => mesh%subpatch%next()
                 ipatch=ipatch+1
              ENDDO !WHILE (ASSOCIATED(sbpitr))

              CALL MPI_Reduce(TotalEnergy_,TotalEnergy, &
              &    1,ppm_mpi_kind,MPI_SUM,0,comm,info)
              or_fail_MPI("MPI_Reduce")

              IF (rank.EQ.0) THEN
                 stdout(istep,TotalEnergy)
                 WRITE(UNIT=3333,FMT='(I16,f20.10)')istep,TotalEnergy
                 CALL ppm_log(caller,cbuf,info)
              ENDIF

           TYPE IS (E_PS)
              TotalEnergy_=zero
              NULLIFY(DTYPE(wpi),DTYPE(wpl))
              sbpitr => mesh%subpatch%begin()
              ipatch=1
              DO WHILE (ASSOCIATED(sbpitr))
                 Nm => sbpitr%nnodes

                 CALL sbpitr%get_field(image,DTYPE(wpi),info)
                 or_fail("Failed to get field r_wp data.")

                 CALL sbpitr%get_field(labels,DTYPE(wpl),info)
                 or_fail("Failed to get field i_wp data.")

                 TotalEnergy_=TotalEnergy_+ &
                 & e_data%CalculateTotalEnergy(DTYPE(wpi),DTYPE(wpl), &
                 & (/Nm(1:__DIME),sbpitr%istart(1:__DIME),sbpitr%iend(1:__DIME)/),info)

                 sbpitr => mesh%subpatch%next()
                 ipatch=ipatch+1
              ENDDO !WHILE (ASSOCIATED(sbpitr))

              CALL MPI_Reduce(TotalEnergy_,TotalEnergy, &
              &    1,ppm_mpi_kind,MPI_SUM,0,comm,info)
              or_fail_MPI("MPI_Reduce")

              IF (rank.EQ.0) THEN
                 stdout(istep,TotalEnergy)
                 WRITE(UNIT=3333,FMT='(I16,f20.10)')istep,TotalEnergy
                 CALL ppm_log(caller,cbuf,info)
              ENDIF

           END SELECT
        ENDIF

        !!! WE are using l_e_merge to differentiate between energy functionals
        SELECT CASE (e_data%m_EnergyFunctional)
        CASE (1)
           l_e_merge=.FALSE.
        CASE DEFAULT
           l_e_merge=.TRUE.
        END SELECT

        ierr=0
        IF (ALLOCATED(labela))      DEALLOCATE(labela,STAT=ierr(1))
        IF (ALLOCATED(candlabela))  DEALLOCATE(candlabela,STAT=ierr(2))
        IF (ANY(ierr(1:2).NE.0)) info=ppm_error_fatal
        or_fail_dealloc("Array deallocation is failed", ppm_error=ppm_error_fatal)

        IF (Mpart.GT.0) THEN
           nsize=2**(CEILING(LOG(REAL(Mpart))/LOG(2.0)))
        ELSE
           nsize=0

           sbpitr => mesh%subpatch%begin()
           ipatch=1
           DO WHILE (ASSOCIATED(sbpitr))
              check_true(<#InnerContourContainer(ipatch)%nb.LE.0#>,"Error there is a particle while Maprt=0!")
              sbpitr => mesh%subpatch%next()
              ipatch=ipatch+1
           ENDDO !WHILE (ASSOCIATED(sbpitr))
        ENDIF
        ALLOCATE(energya(nsize),       STAT=ierr(1))
        ALLOCATE(labela(nsize),        STAT=ierr(2))
        ALLOCATE(candlabela(nsize),    STAT=ierr(3))
        ALLOCATE(ccandlabela(nsize),   STAT=ierr(4))
        ALLOCATE(ndaughtersa(nsize),   STAT=ierr(5))
        ALLOCATE(daughtersa(2*__DIME,nsize),STAT=ierr(6))
        ALLOCATE(nmothersa(nsize),     STAT=ierr(7))
        ALLOCATE(accepteda(nsize),     STAT=ierr(8))
        ALLOCATE(mothersa(2*__DIME,nsize),  STAT=ierr(9))
        IF (ANY(ierr.NE.0)) info=ppm_error_fatal
        or_fail_alloc("Array allocation to the new size is failed", ppm_error=ppm_error_fatal)

        ccandlabela=0
        ndaughtersa=0
        nmothersa=0
        accepteda=-1

        ALLOCATE(missedparticles(1:Npart),STAT=info)
        or_fail_alloc("missedparticles")
        missedparticles=0

        !By DEFAULT no fusion between regions is happening
        e_merge=.FALSE.

        NULLIFY(DTYPE(wpi),DTYPE(wpl),DTYPE(wpp))
        !TODO
        !TOCHECK
        !there is an optimization here to not compute energy for ghost particles
        !and update them during next step computaion
        sbpitr => mesh%subpatch%begin()
        ipatch=1
        DO WHILE (ASSOCIATED(sbpitr))
           CALL sbpitr%get_field(image,DTYPE(wpi),info)
           or_fail("Failed to get field r_wp data.")

           CALL sbpitr%get_field(labels,DTYPE(wpl),info)
           or_fail("Failed to get field i_wp data.")

           CALL sbpitr%get_field(pind,DTYPE(wpp),info)
           or_fail("Failed to get field i_wp data.")

           seed => InnerContourContainer(ipatch)%begin()
           icontour_loop: DO WHILE (ASSOCIATED(seed))
              seedn => seed%first%getValue()

#if   __DIME == __2D
              ppart=DTYPE(wpp)(seedn(1),seedn(2))
#elif __DIME == __3D
              ppart=DTYPE(wpp)(seedn(1),seedn(2),seedn(3))
#endif

              IF (ppart.LE.0) THEN
                 seed => InnerContourContainer(ipatch)%next()
                 CYCLE icontour_loop
              ENDIF

#if   __DIME == __2D
              mylabel = ABS(DTYPE(wpl)(seedn(1),seedn(2)))
#elif __DIME == __3D
              mylabel = ABS(DTYPE(wpl)(seedn(1),seedn(2),seedn(3)))
#endif

              !Assign the particle label and candidate label
              labela(ppart)=mylabel
              candlabela(ppart)=0

              ![NOTE]
              !Important
              !e_length energy for ghost particle is right with the current
              !labels ghost update method
              !within the current implementation at least two layers of ghost
              !would be updated which is OK for this compuation, but for future
              !optimization one should be careful about this!
              etemp=e_data%EvaluateEnergyDifference                 &
              &     (DTYPE(wpi),DTYPE(wpl),seedn,mylabel,0,e_merge) &
              &    +e_length%EvaluateEnergyDifference               &
              &     (DTYPE(wpi),DTYPE(wpl),seedn,mylabel,0,e_merge)

              energya(ppart)=etemp

              !TODO
              !TOCHECK
              IF (etemp.LT.zero.AND.ppart.LE.Npart) THEN
                 accepteda(ppart)=1

                 ALLOCATE(seedi,STAT=info)
                 or_fail_alloc("seedi")
                 CALL seedi%add(seedn)
                 CALL Candidates(ipatch)%push(seedi,info)
                 or_fail("Candidates(ipatch)%push")
              ENDIF
              seed => InnerContourContainer(ipatch)%next()
           ENDDO icontour_loop

           sbpitr => mesh%subpatch%next()
           ipatch=ipatch+1
        ENDDO !WHILE (ASSOCIATED(sbpitr))

        NULLIFY(DTYPE(wpi),DTYPE(wpl),DTYPE(wpp))

        sbpitr => mesh%subpatch%begin()
        qpart=Mpart
        ipatch=1
        DO WHILE (ASSOCIATED(sbpitr))
           Nm => sbpitr%nnodes

           CALL sbpitr%get_field(image,DTYPE(wpi),info)
           or_fail("Failed to get field r_wp data.")

           CALL sbpitr%get_field(labels,DTYPE(wpl),info)
           or_fail("Failed to get field i_wp data.")

           CALL sbpitr%get_field(pind,DTYPE(wpp),info)
           or_fail("Failed to get field i_wp data.")

           seed => InnerContourContainer(ipatch)%begin()
           contour_loop: DO WHILE (ASSOCIATED(seed))
              seedn => seed%first%getValue()

#if   __DIME == __2D
              ppart=DTYPE(wpp)(seedn(1),seedn(2))
#elif __DIME == __3D
              ppart=DTYPE(wpp)(seedn(1),seedn(2),seedn(3))
#endif

              IF (ppart.LE.0) THEN
                 seed => InnerContourContainer(ipatch)%next()
                 CYCLE contour_loop
              ELSE IF (ppart.GT.Npart) THEN
                 IsEnclosedByLabelFGConnectivity=.FALSE.
              ELSE
                 IsEnclosedByLabelFGConnectivity=.TRUE.
              ENDIF

              ! store the label at this point
#if   __DIME == __2D
              mylabel = ABS(DTYPE(wpl)(seedn(1),seedn(2)))
#elif __DIME == __3D
              mylabel = ABS(DTYPE(wpl)(seedn(1),seedn(2),seedn(3)))
#endif

              neigh_loop: DO i=1,FG_ConnectivityType%NumberOfNeighbors
                 ll=seedn+FG_ConnectivityType%NeighborsPoints(:,i)
                 !TOCHECK
                 !TODO
                 !check the pointer which could point to the local memory
                 !to be sure that the cash misses is the least
                 IF (ALL(ll.GE.1.AND.ll.LE.Nm)) THEN
#if   __DIME == __2D
                    newlabel = ABS(DTYPE(wpl)(ll(1),ll(2)))
#elif __DIME == __3D
                    newlabel = ABS(DTYPE(wpl)(ll(1),ll(2),ll(3)))
#endif

                    IF (newlabel.EQ.mylabel) THEN
                       CYCLE neigh_loop
                    ELSE IF (newlabel.EQ.FORBIDDEN) THEN
                       IsEnclosedByLabelFGConnectivity=.FALSE.
                       CYCLE neigh_loop
                    ENDIF

                    IsEnclosedByLabelFGConnectivity=.FALSE.

#if   __DIME == __2D
                    ipart = DTYPE(wpp)(ll(1),ll(2))
#elif __DIME == __3D
                    ipart = DTYPE(wpp)(ll(1),ll(2),ll(3))
#endif

                    ! if there is not particle yet, create one
                    IF (ipart.EQ.0) THEN
                       !TOCHECK:
                       !whether the label is zero or not!
                       !set label of q = 0 in case ppart is the real particle
                       !There is a very special case, when I have part of two region
                       !between two processors with no boundaries particles
                       !In one procesor fusion will happen, while the other
                       !have no idea about
! ! !                        IF (ppart.LE.Npart) THEN
! ! ! #if   __DIME == __2D
! ! !                           IF (DTYPE(wpl)(ll(1),ll(2)).NE.0) THEN
! ! !                              stdout("************** St is wrong here******************")
! ! !                              stdout(Npart,Mpart,ipart,ppart,seedn,istep)
! ! !                              stdout('wpl_2d(ll(1),ll(2))',ll,mylabel,newlabel)
! ! !                           ENDIF
! ! !                           DTYPE(wpl)(ll(1),ll(2))=0
! ! !                           newlabel=0
! ! ! #elif __DIME == __3D
! ! !                           DTYPE(wpl)(ll(1),ll(2),LL(3))=0
! ! !                           newlabel=0
! ! ! #endif
! ! !                        ENDIF

                       IF (newlabel.NE.0) CYCLE neigh_loop

                       e_merge=.FALSE.

                       etemp=e_data%EvaluateEnergyDifference                     &
                       &     (DTYPE(wpi),DTYPE(wpl),ll,newlabel,mylabel,e_merge) &
                       &    +e_length%EvaluateEnergyDifference                   &
                       &     (DTYPE(wpi),DTYPE(wpl),ll,newlabel,mylabel,e_merge)

                       !TOCHECK
                       !this is what I think is right you should be sure about this
                       IF (etemp.GE.zero) CYCLE neigh_loop

                       qpart=qpart+1
#if   __DIME == __2D
                       DTYPE(wpp)(ll(1),ll(2))=qpart
#elif __DIME == __3D
                       DTYPE(wpp)(ll(1),ll(2),ll(3))=qpart
#endif

                       ALLOCATE(seedi,STAT=info)
                       or_fail_alloc("seedi")
                       CALL seedi%add(ll)
                       CALL Candidates(ipatch)%push(seedi,info)
                       or_fail("Candidates(ipatch)%push")

                       IF (nsize.LT.qpart) THEN
#include "./ppm_rc_energy_compute_resize.inc"
                       ENDIF

                       ! register q in p's daughter list
                       ndaughtersa(ppart)=ndaughtersa(ppart)+1
                       daughtersa(ndaughtersa(ppart),ppart)=qpart

                       candlabela(qpart)=mylabel
                       ccandlabela(qpart)=1
                       nmothersa(qpart)=1
                       !register p in q's parent list
                       mothersa(nmothersa(qpart),qpart)=ppart
!                        ndaughtersa(qpart)=0
                       accepteda(qpart)=1
                       energya(qpart)=etemp

                       labela(qpart)=0
                    ! there is a particle at this point
                    ELSE
                       IF (mylabel.EQ.candlabela(ipart)) THEN
                          !register p in q's parent list
                          nmothersa(ipart)=nmothersa(ipart)+1
                          mothersa(nmothersa(ipart),ipart)=ppart

                          !register q in p's daughter list
                          ndaughtersa(ppart)=ndaughtersa(ppart)+1
                          daughtersa(ndaughtersa(ppart),ppart)=ipart

                          ccandlabela(ipart)=ccandlabela(ipart)+1

                          IF (ppart.GT.Npart) THEN
                             IF (ipart.LE.Npart) missedparticles(ipart)=0
                          ENDIF
                       ELSE
                          !I did this crazy thing which seems really nice
                          !why I should assign a mother daughter relation for
                          !a particle which will not be one,
                          !so I just removed this part and put an if condition

                          e_merge=l_e_merge
                          ! the particle exists and the energy for moving to
                          ! candlabel has previously been computed
                          etemp1=energya(ipart)

                          etemp2=e_data%EvaluateEnergyDifference                     &
                          &      (DTYPE(wpi),DTYPE(wpl),ll,newlabel,mylabel,e_merge) &
                          &     +e_length%EvaluateEnergyDifference                   &
                          &      (DTYPE(wpi),DTYPE(wpl),ll,newlabel,mylabel,e_merge)

                          !I have seen a very strange behavior,
                          !I can have a single particle among 4 other particles
                          !as it is single it has an energy of -bigs
                          !the energy for changing to any of its neighbors is also -bigs
                          !and it will cause topology violation so it is not accepted to move
                          !by changing the GT to GE it seems that it resolves the problem
                          IF (etemp1.GE.etemp2) THEN
                             !register p in q's parent list
                             nmothersa(ipart)=nmothersa(ipart)+1
                             mothersa(nmothersa(ipart),ipart)=ppart

                             !register q in p's daughter list
                             ndaughtersa(ppart)=ndaughtersa(ppart)+1
                             daughtersa(ndaughtersa(ppart),ppart)=ipart

                             IF (ppart.GT.Npart) THEN
                                IF (ipart.LE.Npart) missedparticles(ipart)=0
                             ENDIF

                             candlabela(ipart)=mylabel
                             !TODO
                             !TOCHECK
                             ccandlabela(ipart)=COUNT(labela(mothersa(1:nmothersa(ipart),ipart)).EQ.mylabel)
                             energya(ipart)=etemp2
                             IF (etemp2.LT.zero) THEN
                                !If this move has been rejected before, we should not change the decision
                                !as this point could be chosen for some merging and should not move at all
                                IF (accepteda(ipart).EQ.-1) THEN
                                   accepteda(ipart)=1
                                   IF (etemp1.GE.zero) THEN
                                      !If etemp1 greater or equal than zero, it means that
                                      !this particle is not among the candidate list, and now
                                      !its energy value has been changed and it should be in
                                      !the list
                                      ALLOCATE(seedi,STAT=info)
                                      or_fail_alloc("seedi")
                                      CALL seedi%add(ll)
                                      CALL Candidates(ipatch)%push(seedi,info)
                                      or_fail("Candidates(ipatch)%push")
                                   ENDIF
                                ENDIF
                             ENDIF
                             !TOCHECK
                             !TODO
                             !find the better way
                             IF (AllowFusion.AND.labela(ipart).NE.0.AND.candlabela(ipart).NE.0) THEN
                                IF (.NOT.l_e_merge) THEN
                                   e_merge=.TRUE.
                                   etemp=e_data%EvaluateEnergyDifference(DTYPE(wpi), &
                                   &     DTYPE(wpl),ll,newlabel,mylabel,e_merge)
                                ENDIF
                                IF (e_merge) THEN
                                   IF (AllowFusionZ) THEN
#if   __DIME == __3D
                                      IF (ABS(ll(3)-seedn(3)).EQ.1) THEN
                                         ALLOCATE(seedi,STAT=info)
                                         or_fail_alloc("seedi")

                                         CALL seedi%add(ll(1),ll(2),ll(3),newlabel, &
                                         &    seedn(1),seedn(2),seedn(3),mylabel)

                                         CALL CompetingRegions(ipatch)%push(seedi,info)
                                         or_fail("CompetingRegions%push")

                                         !The seed should be exempted from any move
                                         accepteda(ipart)=0
                                         accepteda(ppart)=0
                                      ENDIF
#endif
                                   ELSE
                                      ALLOCATE(seedi,STAT=info)
                                      or_fail_alloc("seedi")
#if   __DIME == __2D
                                      CALL seedi%add(ll(1),ll(2),newlabel, &
                                      &    seedn(1),seedn(2),mylabel)
#elif __DIME == __3D
                                      CALL seedi%add(ll(1),ll(2),ll(3),newlabel, &
                                      &    seedn(1),seedn(2),seedn(3),mylabel)
#endif

                                      CALL CompetingRegions(ipatch)%push(seedi,info)
                                      or_fail("CompetingRegions%push")

                                      !The seed should be exempted from any move
                                      accepteda(ipart)=0
                                      accepteda(ppart)=0
                                   ENDIF
                                ENDIF
                             ENDIF
                          ENDIF !(etemp1.GT.etemp2)
                       ENDIF !(mylabel.EQ.candlabela(ipart))
                    ENDIF !(ipart.LT.0)
                 ELSE IF (ppart.LE.Npart) THEN
#if   __DIME == __2D
                    ipart = DTYPE(wpp)(ll(1),ll(2))
#elif __DIME == __3D
                    ipart = DTYPE(wpp)(ll(1),ll(2),ll(3))
#endif

                    ! if there is not particle yet, this one could be a missedparticles
                    IF (ipart.GT.0) THEN !CYCLE neigh_loop
                       !!!It is possible that a real particle has a candidate in
                       !!!background labelled ghost and will be missed from
                       !!!the list for violation of topological movement
                       missedparticles(ppart)=1

! #if   __DIME == __2D
!                        newlabel = ABS(DTYPE(wpl)(ll(1),ll(2)))
! #elif __DIME == __3D
!                        newlabel = ABS(DTYPE(wpl)(ll(1),ll(2),ll(3)))
! #endif

                       IF (IsEnclosedByLabelFGConnectivity) THEN
#if   __DIME == __2D
                          newlabel = -DTYPE(wpl)(ll(1),ll(2))
#elif __DIME == __3D
                          newlabel = -DTYPE(wpl)(ll(1),ll(2),ll(3))
#endif
                          IF (newlabel.NE.mylabel) THEN
                             IsEnclosedByLabelFGConnectivity=.FALSE.
                          ENDIF
                       ENDIF


!                        e_merge=l_e_merge
!                        ! the particle exists and the energy for moving to
!                        ! candlabel has previously been computed
!                        etemp1=energya(ipart)
!
!                        etemp2=e_data%EvaluateEnergyDifference                     &
!                        &      (DTYPE(wpi),DTYPE(wpl),ll,newlabel,mylabel,e_merge) &
!                        &     +e_length%EvaluateEnergyDifference                   &
!                        &      (DTYPE(wpi),DTYPE(wpl),ll,newlabel,mylabel,e_merge)
!
!                        IF (etemp1.GE.etemp2) THEN
!                           IF (AllowFusion) THEN
!                              IF (.NOT.l_e_merge) THEN
!                                 e_merge=.TRUE.
!                                 etemp=e_data%EvaluateEnergyDifference(DTYPE(wpi), &
!                                 &     DTYPE(wpl),ll,newlabel,mylabel,e_merge)
!                              ENDIF
!                              IF (e_merge) THEN
!                                 IF (AllowFusionZ) THEN
! #if   __DIME == __3D
!                                    IF (ABS(ll(3)-seedn(3)).EQ.1) THEN
!                                       ALLOCATE(seedi,STAT=info)
!                                       or_fail_alloc("seedi")
!
!                                       CALL seedi%add(ll(1),ll(2),ll(3),newlabel, &
!                                       &    seedn(1),seedn(2),seedn(3),mylabel)
!
!                                       CALL CompetingRegions(ipatch)%push(seedi,info)
!                                       or_fail("CompetingRegions%push")
!
!                                       !The seed should be exempted from any move
!                                       accepteda(ipart)=0
!                                       accepteda(ppart)=0
!                                    ENDIF
! #endif
!                                 ELSE
!                                    ALLOCATE(seedi,STAT=info)
!                                    or_fail_alloc("seedi")
! #if   __DIME == __2D
!                                    CALL seedi%add(ll(1),ll(2),newlabel, &
!                                    &    seedn(1),seedn(2),mylabel)
! #elif __DIME == __3D
!                                    CALL seedi%add(ll(1),ll(2),ll(3),newlabel, &
!                                    &    seedn(1),seedn(2),seedn(3),mylabel)
! #endif
!
!                                    CALL CompetingRegions(ipatch)%push(seedi,info)
!                                    or_fail("CompetingRegions%push")
!
!                                    !The seed should be exempted from any move
!                                    accepteda(ipart)=0
!                                    accepteda(ppart)=0
!                                 ENDIF
!                              ENDIF
!                           ENDIF
!                        ENDIF !(etemp1.GT.etemp2)

                       CYCLE neigh_loop
                    ENDIF

#if   __DIME == __2D
                    newlabel = DTYPE(wpl)(ll(1),ll(2))
#elif __DIME == __3D
                    newlabel = DTYPE(wpl)(ll(1),ll(2),ll(3))
#endif

                    ! If thers is a label, then it is inside
                    !one region
                    IF (newlabel.NE.0) THEN
                       IF (newlabel.NE.mylabel) THEN
                          IsEnclosedByLabelFGConnectivity=.FALSE.
                       ENDIF
                       CYCLE neigh_loop
                    ENDIF

                    IsEnclosedByLabelFGConnectivity=.FALSE.

                    e_merge=.FALSE.

                    etemp=e_data%EvaluateEnergyDifference              &
                    &     (DTYPE(wpi),DTYPE(wpl),ll,0,mylabel,e_merge) &
                    &    +e_length%EvaluateEnergyDifference            &
                    &     (DTYPE(wpi),DTYPE(wpl),ll,0,mylabel,e_merge)

                    !TOCHECK
                    !TODO
                    !this is what I think is right, you should be sure about this
                    IF (etemp.GE.zero) CYCLE neigh_loop

                    !!!It is possible that a real particle has a candidate in
                    !!!background labelled ghost and will be missed from
                    !!!the list for violation of topological movement
                    missedparticles(ppart)=1
                 ENDIF !(ALL(ll.GE.1.AND.ll.LE.Nm)) THEN
              ENDDO neigh_loop

              IF (IsEnclosedByLabelFGConnectivity) THEN
#if   __DIME == __2D
                 DTYPE(wpl)(seedn(1),seedn(2))=mylabel
                 DTYPE(wpp)(seedn(1),seedn(2))=0
#elif __DIME == __3D
                 DTYPE(wpl)(seedn(1),seedn(2),seedn(3))=mylabel
                 DTYPE(wpp)(seedn(1),seedn(2),seedn(3))=0
#endif
                 del_parts=del_parts+1
                 list_del_parts(del_parts)=ppart
              ENDIF

              seed => InnerContourContainer(ipatch)%next()
           ENDDO contour_loop

           sbpitr => mesh%subpatch%next()
           ipatch=ipatch+1
        ENDDO !WHILE (ASSOCIATED(sbpitr))

        NpartNew=qpart-Mpart

        DEALLOCATE(mothersa,STAT=info)
        or_fail_dealloc("mothersa")

        IF (nsize.GT.qpart) THEN
           nsize=qpart
#include "./ppm_rc_energy_compute_correctsize.inc"
        ENDIF

        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
      9999 CONTINUE
        CALL substop(caller,t0,info)
        RETURN
      END SUBROUTINE DTYPE(ppm_rc_energy_compute)

