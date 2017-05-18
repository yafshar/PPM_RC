        LOGICAL FUNCTION MCMCOffBoundarySample_2d(Growth,LabelIn)
          !!!  For the insertion and deletion of floating particles, we use a
          !!!  constant energy, therefore the MH ratio be equal to the FBR.
          !!!
          !!!  Note that we also could use the energy as if we would do the job as
          !!!  proposal distribution, but we would have to calculate the backward
          !!!  energy as well.
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_fatal
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_error,ppm_err_sub_failed
          USE ppm_module_mesh_typedef, ONLY : ppm_t_subpatch_

          USE ppm_rc_module_global, ONLY : half,mesh,labels
          USE ppm_rc_module_rnd, ONLY : ppm_rc_Saru_RPRNG
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          LOGICAL, INTENT(IN   ) :: Growth
          INTEGER, INTENT(IN   ) :: LabelIn
          !!! Absolute value label from regionsLabel
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(MCMCParticle) :: tmpParticle

          CLASS(ppm_t_subpatch_), POINTER :: sbpitr

          REAL(ppm_kind_double) :: t0

          INTEGER, CONTIGUOUS, DIMENSION(:,:), POINTER :: labels_2d
          INTEGER,             DIMENSION(:),   POINTER :: Nm
          INTEGER,             DIMENSION(2)            :: index_
          INTEGER                                      :: info

          LOGICAL :: Child

          CHARACTER(LEN=ppm_char) :: caller='MCMCOffBoundarySample'
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          NULLIFY(labels_2d)

          sbpitr => mesh%subpatch%begin()
          DO WHILE (ASSOCIATED(sbpitr))
             NM => sbpitr%nnodes

             CALL sbpitr%get_field(labels,labels_2d,info)
             or_fail("Failed to get field labels_2d data.",ppm_error=ppm_error_fatal)

             ! get index from Edge Density
             index_=MCMCGetIndexFromEdgeDensity(labels_2d,Nm)

             ! search for the particle at this index
             tmpParticle=MCMCFloatingParticles%search(index_(1),index_(2))

             ! Particle does not exist
             IF (tmpParticle%candlabel.EQ.-1) THEN
                ! insert particle
                IF (Growth) THEN
                   ! the target-distribution probability ratio
                   Child=ppm_rc_Saru_RPRNG().GT.half

                   tmpParticle%candlabel=MERGE(LabelIn,0,Child)
                   ! TO CHECK
                   ! TO DO
                   ! This one does not exist in the original code (to compute the proposal)
                   tmpParticle%proposal=MCMCproposal(labels_2d,index_)

                   MCMCOffBoundarySample_2d= &
                   &  MCMCInsertFloatingParticle(index_(1),index_(2),tmpParticle,.FALSE.,.FALSE.)
                   GOTO 9999
                ! reject
                ELSE
                   MCMCOffBoundarySample_2d=.FALSE.
                   GOTO 9999
                ENDIF
             ELSE
                ! Particle exists which means reject
                IF (Growth) THEN
                   MCMCOffBoundarySample_2d=.FALSE.
                   GOTO 9999
                ! .NOT.Growth.AND.Particle exists
                ELSE
                   MCMCOffBoundarySample_2d= &
                   & MCMCEraseFloatingParticle(index_(1),index_(2),tmpParticle,tmpParticle,.FALSE.)
                   GOTO 9999
                ENDIF
             ENDIF

             sbpitr => mesh%subpatch%next()
          ENDDO

        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        END FUNCTION MCMCOffBoundarySample_2d

        LOGICAL FUNCTION MCMCOffBoundarySample_3d(Growth,LabelIn)
          !!!  For the insertion and deletion of floating particles, we use a
          !!!  constant energy, therefore the MH ratio be equal to the FBR.
          !!!
          !!!  Note that we also could use the energy as if we would do the job as
          !!!  proposal distribution, but we would have to calculate the backward
          !!!  energy as well.
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_fatal
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_error,ppm_err_sub_failed
          USE ppm_module_mesh_typedef, ONLY : ppm_t_subpatch_

          USE ppm_rc_module_global, ONLY : half,mesh,labels
          USE ppm_rc_module_rnd, ONLY : ppm_rc_Saru_RPRNG
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          LOGICAL, INTENT(IN   ) :: Growth
          INTEGER, INTENT(IN   ) :: LabelIn
          !!! Absolute value label from regionsLabel
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(MCMCParticle) :: tmpParticle

          CLASS(ppm_t_subpatch_), POINTER :: sbpitr

          REAL(ppm_kind_double) :: t0

          INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER :: labels_3d
          INTEGER,             DIMENSION(:),     POINTER :: Nm
          INTEGER,             DIMENSION(3)              :: index_
          INTEGER                                        :: info

          LOGICAL :: Child

          CHARACTER(LEN=ppm_char) :: caller='MCMCOffBoundarySample'
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          NULLIFY(labels_3d)

          sbpitr => mesh%subpatch%begin()
          DO WHILE (ASSOCIATED(sbpitr))
             NM => sbpitr%nnodes

             CALL sbpitr%get_field(labels,labels_3d,info)
             or_fail("Failed to get field labels_3d data.",ppm_error=ppm_error_fatal)

             ! get index from Edge Density
             index_=MCMCGetIndexFromEdgeDensity(labels_3d,Nm)

             ! search for the particle at this index
             tmpParticle=MCMCFloatingParticles%search(index_(1),index_(2),index_(3))

             ! Particle does not exist
             IF (tmpParticle%candlabel.EQ.-1) THEN
                ! insert particle
                IF (Growth) THEN
                   ! the target-distribution probability ratio
                   Child=ppm_rc_Saru_RPRNG().GT.half

                   tmpParticle%candlabel=MERGE(LabelIn,0,Child)
                   ! TO CHECK
                   ! TO DO
                   ! This one does not exist in the original code (to compute the proposal)
                   tmpParticle%proposal=MCMCproposal(labels_3d,index_)

                   MCMCOffBoundarySample_3d= &
                   &  MCMCInsertFloatingParticle(index_(1),index_(2),index_(3),tmpParticle,.FALSE.,.FALSE.)
                   GOTO 9999
                ! reject
                ELSE
                   MCMCOffBoundarySample_3d=.FALSE.
                   GOTO 9999
                ENDIF
             ! Particle exists
             ELSE
                ! reject
                IF (Growth) THEN
                   MCMCOffBoundarySample_3d=.FALSE.
                   GOTO 9999
                ! .NOT.Growth.AND.Particle exists
                ELSE
                   MCMCOffBoundarySample_3d= &
                   & MCMCEraseFloatingParticle(index_(1),index_(2),index_(3),tmpParticle,tmpParticle,.FALSE.)
                   GOTO 9999
                ENDIF
             ENDIF

             sbpitr => mesh%subpatch%next()
          ENDDO

        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        END FUNCTION MCMCOffBoundarySample_3d

        !----------------------------------------------------------------------
        ! Find particle B
        ! In case of pair proposals, we find a partner for each proposed particle.
        ! Now, we know A and Q(A), and we need to build another (2nd step) discrete proposal distribution.
        ! We sample from it to determine the partner particle B.
        ! Furthermore we calculate the conditional proposal probability Q(B|A).
        ! In a second step we calculate the conditional Q(A|B).
        ! The same then needs to be done for the backward probabilities Qb(A), Qb(B), Qb(A|B) and Qb(B|A).
        !
        ! Notation:
        ! Q is the forward and Qb the backward probability.
        ! A is a forward praticle and B' the backward particle.
        ! Qb always assumes backward particles as its arguments!
        ! Hence, Qb_A_B is the probabily Qb(A'|B').
        !----------------------------------------------------------------------
        FUNCTION MCMCFindParticleB(ppmrcdim) RESULT(info)
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_fatal,ppm_error_error
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_error,ppm_err_sub_failed, &
          &   ppm_err_alloc,ppm_err_dealloc
          USE ppm_module_mesh_typedef, ONLY : ppm_t_subpatch_

          USE ppm_rc_module_global, ONLY : oned,mesh,image,labels,MCMCuseBiasedProposal
          USE ppm_rc_module_rnd, ONLY : ppm_rc_GetPartDistrIndex,ppm_rc_Saru_IPRNG,  &
          &   ppm_rc_DestroyParticlesDiscrDistr,ppm_rc_GenerateParticlesFwdProposalsDiscrDistr
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER, INTENT(IN   ) :: ppmrcdim
          !!! Problem dimension

          INTEGER                :: info
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(MCMCParticle), DIMENSION(:), ALLOCATABLE :: Parts_Q_BgivenA
          TYPE(MCMCParticle), DIMENSION(:), ALLOCATABLE :: Parts_Qb_AgivenB
          TYPE(MCMCParticle)                            :: Part_A
          TYPE(MCMCParticle)                            :: RPart_A
          TYPE(MCMCParticle)                            :: Part_B

          CLASS(ppm_t_subpatch_), POINTER :: sbpitr

          REAL(MK), CONTIGUOUS, DIMENSION(:,:),   POINTER     :: image_2d
          REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER     :: image_3d
          REAL(MK),             DIMENSION(:),     ALLOCATABLE :: ProposalsVector
          REAL(MK)                                            :: Normalizer_Q_B_A
          REAL(MK)                                            :: Normalizer_Qb_A_B
          REAL(MK)                                            :: qb_A_B_unnorm
          REAL(ppm_kind_double)                               :: t0

          INTEGER, CONTIGUOUS, DIMENSION(:,:),   POINTER     :: labels_2d
          INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER     :: labels_3d
          INTEGER,             DIMENSION(:),     POINTER     :: Nm
          INTEGER,             DIMENSION(:,:),   ALLOCATABLE :: Parts_Q_BgivenA_coords
          INTEGER,             DIMENSION(:,:),   ALLOCATABLE :: Parts_Qb_AgivenB_coords
          INTEGER,             DIMENSION(ppmrcdim)           :: index_A
          INTEGER,             DIMENSION(ppmrcdim)           :: index_B
          INTEGER                                            :: nsize
          INTEGER                                            :: i
          INTEGER                                            :: cbox
          INTEGER                                            :: nx,ny,nz
          INTEGER                                            :: ParticleIndex
          INTEGER                                            :: CandidateLabel

          CHARACTER(LEN=ppm_char) :: caller='MCMCFindParticleB'
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          SELECT CASE (ppmrcdim)
          CASE (2)
             NULLIFY(labels_2d,image_2d)
             sbpitr => mesh%subpatch%begin()
             DO WHILE (ASSOCIATED(sbpitr))
                Nm => sbpitr%nnodes

                CALL sbpitr%get_field(image,image_2d,info)
                or_fail("Failed to get field image_2d data.",ppm_error=ppm_error_fatal)

                CALL sbpitr%get_field(labels,labels_2d,info)
                or_fail("Failed to get field labels_2d data.",ppm_error=ppm_error_fatal)

                !-------------------------------------------------------------------------
                ! Using pair potential with step size > 1 leads to slightly wrong proposals
                ! and hence the detailed balance is not guaranteed.
                !-------------------------------------------------------------------------
                ! Find particle A:
                !-------------------------------------------------------------------------
                index_A=MCMC_CandidateMove_Index(:,1)

                info=MCMCApplyParticle_2d(image_2d,labels_2d,index_A,Nm,MCMC_CandidateMove(1),.TRUE.)
                or_fail("MCMCApplyParticle")

                CandidateLabel=MCMC_CandidateMove(1)%candlabel

                !-------------------------------------------------------------------------
                ! BTW we have to remember if A' is a floating particle in the state x -> A.
                !-------------------------------------------------------------------------
                MCMC_Particle_Ab_IsFloating(1)=MCMCParticleHasFloatingProperty(labels_2d,index_A,CandidateLabel)

                !-------------------------------------------------------------------------
                ! Get the particles involved in the second step of the
                ! proposal (with updated proposals)
                ! nsize is at least 1 to count for particle A
                !-------------------------------------------------------------------------
                nsize=MCMCgetParticlesInBGNeighborhood(labels_2d,index_A,Nm,Parts_Q_BgivenA,Parts_Q_BgivenA_coords)

                !-------------------------------------------------------------------------
                ! Insert particle A in the set for B such that it is possible
                ! that we have a single particle move.
                !-------------------------------------------------------------------------
                Parts_Q_BgivenA(nsize)%candlabel=CandidateLabel

                Part_A=Parts_Q_BgivenA(nsize)

                !-------------------------------------------------------------------------
                ! Find B
                !-------------------------------------------------------------------------

                !-------------------------------------------------------------------------
                ! Choose B from Q(B|A) and calculate Q(B|A).
                !-------------------------------------------------------------------------
                IF (MCMCuseBiasedProposal) THEN
                   ALLOCATE(ProposalsVector(nsize),STAT=info)
                   or_fail_alloc("ProposalsVector")

                   FORALL (i=1:nsize) ProposalsVector(i)=Parts_Q_BgivenA(i)%proposal

                   Normalizer_Q_B_A=SUM(ProposalsVector)

                   !-------------------------------------------------------------------------
                   ! Create a discrete distribution over particles
                   !-------------------------------------------------------------------------
                   info=ppm_rc_GenerateParticlesFwdProposalsDiscrDistr(ProposalsVector,nsize)
                   or_fail("ppm_rc_GenerateParticlesFwdProposalsDiscrDistr",ppm_error=ppm_error_fatal)

                   ParticleIndex=ppm_rc_GetPartDistrIndex()

                   Part_B=Parts_Q_BgivenA(ParticleIndex)

                   !-------------------------------------------------------------------------
                   ! The value m_Proposal of B is currently equal to Q_B_A.
                   !-------------------------------------------------------------------------
                   MCMC_q_B_A(1)=REAL(Part_B%proposal,ppm_kind_double)/REAL(Normalizer_Q_B_A,ppm_kind_double)
                ELSE
                   ParticleIndex=ppm_rc_Saru_IPRNG(nsize)

                   Part_B=Parts_Q_BgivenA(ParticleIndex)

                   !-------------------------------------------------------------------------
                   ! The value m_Proposal of B is currently equal to Q_B_A.
                   !-------------------------------------------------------------------------
                   MCMC_q_B_A(1)=oned/REAL(nsize,ppm_kind_double)
                ENDIF !MCMCuseBiasedProposal

                !-------------------------------------------------------------------------
                ! store B (and its original label).
                !-------------------------------------------------------------------------
                MCMC_PartnerMove(1)=Part_B

                index_B=Parts_Q_BgivenA_coords(:,ParticleIndex)

                !-------------------------------------------------------------------------
                ! TOCHECK
                !-------------------------------------------------------------------------
                MCMC_PartnerMove_Index(:,1)=index_B

                DEALLOCATE(Parts_Q_BgivenA,Parts_Q_BgivenA_coords,STAT=info)
                or_fail_dealloc("Parts_Q_BgivenA,Parts_Q_BgivenA_coords")

                IF (ALL(index_A.EQ.index_B).AND.Part_A%candlabel.EQ.Part_B%candlabel) THEN
                   MCMC_SingleParticleMoveForPairProposals=.TRUE.
                   MCMC_LabelsBeforeJump_B(1)=MCMC_LabelsBeforeJump_A(1)
                ELSE
                   MCMC_LabelsBeforeJump_B(1)=ABS(labels_2d(index_B(1),index_B(2)))
                ENDIF

                !-------------------------------------------------------------------------
                ! Get the reverse particle of A (without proposal update as it is not necessary):
                !-------------------------------------------------------------------------
                RPart_A=MCMCParticle(MCMC_LabelsBeforeJump_A(1),MCMC_CandidateMove(1)%proposal)

                !-------------------------------------------------------------------------
                ! In case that Part_A == Part_B, we must already undo the simulated
                ! move in order to calculate Q'(B'|A') (== Q'(A'|B'))
                !-------------------------------------------------------------------------
                IF (MCMC_SingleParticleMoveForPairProposals) THEN
                   info=MCMCApplyParticle_2d(image_2d,labels_2d,index_A,Nm,RPart_A,.TRUE.)
                   or_fail("MCMCApplyParticle")
                ENDIF

                !-------------------------------------------------------------------------
                ! Get the reverse particle of B:
                ! in the current state of the label image and the
                ! containers we can calculate qb_A'_B' as well. We assume now
                ! that B' was applied and we calculate the probability for A'.
                !-------------------------------------------------------------------------
                nsize=MCMCgetParticlesInFGNeighborhood(labels_2d,index_B,Nm,Parts_Qb_AgivenB,Parts_Qb_AgivenB_coords)

                IF (MCMCuseBiasedProposal) THEN
                   Normalizer_Qb_A_B=SUM(Parts_Qb_AgivenB(:)%proposal)
                   qb_A_B_unnorm=MCMCproposal(labels_2d,index_A)

                   MCMC_qb_A_B(1)=REAL(qb_A_B_unnorm,ppm_kind_double)/REAL(Normalizer_Qb_A_B,ppm_kind_double)
                ELSE
                   MCMC_qb_A_B(1)=oned/REAL(nsize,ppm_kind_double)
                ENDIF

                DEALLOCATE(Parts_Qb_AgivenB,Parts_Qb_AgivenB_coords,STAT=info)
                or_fail_dealloc("Parts_Qb_AgivenB,Parts_Qb_AgivenB_coords")

                IF (.NOT.MCMC_SingleParticleMoveForPairProposals) THEN
                   !-------------------------------------------------------------------------
                   ! undo the simulated move.
                   !-------------------------------------------------------------------------
                   info=MCMCApplyParticle_2d(image_2d,labels_2d,index_A,Nm,RPart_A,.TRUE.)
                   or_fail("MCMCApplyParticle")


                   nx=INT(REAL(index_B(1),MK)/MCMCcellsize(1))
                   IF (nx.EQ.MCMCcellNm(1)) nx=nx-1

                   ny=INT(REAL(index_B(2),MK)/MCMCcellsize(2))
                   IF (ny.EQ.MCMCcellNm(2)) ny=ny-1

                   cbox=1+nx+MCMCcellNm(1)*ny

                   !-------------------------------------------------------------------------
                   ! Now we can calculate Q_B (as we now know B and the original
                   ! state has been recovered).
                   !-------------------------------------------------------------------------
                   IF (MCMCIsRegularParticle(index_B(1),index_B(2),cbox,Part_B)) THEN
                      IF (MCMCuseBiasedProposal) THEN
                         !-------------------------------------------------------------------------
                         ! B Proposal needs to be recalculated here:
                         !-------------------------------------------------------------------------
                         MCMC_q_B(1)=REAL(MCMCproposal(labels_2d,index_B),ppm_kind_double)/(MCMCTotalNormalizer+MCMCTotalNormalizerlocal)
                      ELSE
                         MCMC_q_B(1)=oned/(MCMCTotalNormalizer+MCMCTotalNormalizerlocal)
                      ENDIF
                   ELSE
                      MCMC_q_B(1)=zerod
                   ENDIF
                ENDIF

                IF (MCMCuseBiasedProposal) THEN
                   DEALLOCATE(ProposalsVector,STAT=info)
                   or_fail_dealloc("ProposalsVector")

                   !-------------------------------------------------------------------------
                   ! make sure that there is no discrete distribution available in memory
                   !-------------------------------------------------------------------------
                   info=ppm_rc_DestroyParticlesDiscrDistr()
                   or_fail("ppm_rc_DestroyParticlesDiscrDistr")
                ENDIF

                sbpitr => mesh%subpatch%next()
             ENDDO !ASSOCIATED(sbpitr)
             NULLIFY(labels_2d,image_2d)
          CASE (3)
             NULLIFY(labels_3d,image_3d)
             sbpitr => mesh%subpatch%begin()
             DO WHILE (ASSOCIATED(sbpitr))
                CALL sbpitr%get_field(labels,labels_3d,info)
                or_fail("Failed to get field labels_3d data.",ppm_error=ppm_error_fatal)

                CALL sbpitr%get_field(image,image_3d,info)
                or_fail("Failed to get field labels_3d data.",ppm_error=ppm_error_fatal)

                !-------------------------------------------------------------------------
                ! Find particle A:
                !-------------------------------------------------------------------------
                index_A=MCMC_CandidateMove_Index(:,1)

                info=MCMCApplyParticle_3d(image_3d,labels_3d,index_A,Nm,MCMC_CandidateMove(1),.TRUE.)
                or_fail("MCMCApplyParticle")

                CandidateLabel=MCMC_CandidateMove(1)%candlabel

                !-------------------------------------------------------------------------
                ! BTW we have to remember if A' is a floating particle in the state x -> A.
                !-------------------------------------------------------------------------
                MCMC_Particle_Ab_IsFloating(1)=MCMCParticleHasFloatingProperty(labels_3d,index_A,CandidateLabel)

                !-------------------------------------------------------------------------
                ! Get the particles involved in the second step of the
                ! proposal (with updated proposals)
                ! nsize is at least 1 to count for particle A
                !-------------------------------------------------------------------------
                nsize=MCMCgetParticlesInBGNeighborhood(labels_3d,index_A,Nm,Parts_Q_BgivenA,Parts_Q_BgivenA_coords)

                !-------------------------------------------------------------------------
                ! Insert particle A in the set for B such that it is possible
                ! that we have a single particle move.
                !-------------------------------------------------------------------------
                Parts_Q_BgivenA(nsize)%candlabel=CandidateLabel

                Part_A=Parts_Q_BgivenA(nsize)

                !-------------------------------------------------------------------------
                ! Find B
                !-------------------------------------------------------------------------

                !-------------------------------------------------------------------------
                ! Choose B from Q(B|A) and calculate Q(B|A).
                !-------------------------------------------------------------------------
                IF (MCMCuseBiasedProposal) THEN
                   ALLOCATE(ProposalsVector(nsize),STAT=info)
                   or_fail_alloc("ProposalsVector")

                   FORALL (i=1:nsize) ProposalsVector(i)=Parts_Q_BgivenA(i)%proposal

                   Normalizer_Q_B_A=SUM(ProposalsVector)

                   !-------------------------------------------------------------------------
                   ! Create a discrete distribution over particles
                   !-------------------------------------------------------------------------
                   info=ppm_rc_GenerateParticlesFwdProposalsDiscrDistr(ProposalsVector,nsize)
                   or_fail("ppm_rc_GenerateParticlesFwdProposalsDiscrDistr", &
                   & ppm_error=ppm_error_fatal)

                   ParticleIndex=ppm_rc_GetPartDistrIndex()

                   Part_B=Parts_Q_BgivenA(ParticleIndex)

                   !-------------------------------------------------------------------------
                   ! The value m_Proposal of B is currently equal to Q_B_A.
                   !-------------------------------------------------------------------------
                   MCMC_q_B_A(1)=REAL(Part_B%proposal,ppm_kind_double)/REAL(Normalizer_Q_B_A,ppm_kind_double)
                ELSE
                   ParticleIndex=ppm_rc_Saru_IPRNG(nsize)

                   Part_B=Parts_Q_BgivenA(ParticleIndex)

                   !-------------------------------------------------------------------------
                   ! The value m_Proposal of B is currently equal to Q_B_A.
                   !-------------------------------------------------------------------------
                   MCMC_q_B_A(1)=oned/REAL(nsize,ppm_kind_double)
                ENDIF !MCMCuseBiasedProposal

                !-------------------------------------------------------------------------
                ! store B (and its original label).
                !-------------------------------------------------------------------------
                MCMC_PartnerMove(1)=Part_B

                index_B=Parts_Q_BgivenA_coords(:,ParticleIndex)

                !-------------------------------------------------------------------------
                ! TOCHECK
                !-------------------------------------------------------------------------
                MCMC_PartnerMove_Index(:,1)=index_B

                DEALLOCATE(Parts_Q_BgivenA,Parts_Q_BgivenA_coords,STAT=info)
                or_fail_dealloc("Parts_Q_BgivenA,Parts_Q_BgivenA_coords")

                IF (ALL(index_A.EQ.index_B).AND.Part_A%candlabel.EQ.Part_B%candlabel) THEN
                   MCMC_SingleParticleMoveForPairProposals=.TRUE.
                   MCMC_LabelsBeforeJump_B(1)=MCMC_LabelsBeforeJump_A(1)
                ELSE
                   MCMC_LabelsBeforeJump_B(1)=ABS(labels_3d(index_B(1),index_B(2),index_B(3)))
                ENDIF

                !-------------------------------------------------------------------------
                ! Get the reverse particle of A (without proposal update as it
                ! is not necessary):
                !-------------------------------------------------------------------------
                RPart_A=MCMCParticle(MCMC_LabelsBeforeJump_A(1),MCMC_CandidateMove(1)%proposal)

                !-------------------------------------------------------------------------
                ! In case that Part_A == Part_B, we must already undo the simulated
                ! move in order to calculate Q'(B'|A') (== Q'(A'|B'))
                !-------------------------------------------------------------------------
                IF (MCMC_SingleParticleMoveForPairProposals) THEN
                   info=MCMCApplyParticle_3d(image_3d,labels_3d,index_A,Nm,RPart_A,.TRUE.)
                   or_fail("MCMCApplyParticle")
                ENDIF

                !-------------------------------------------------------------------------
                ! Get the reverse particle of B:
                ! in the current state of the label image and the
                ! containers we can calculate qb_A'_B' as well. We assume now
                ! that B' was applied and we calculate the probability for A'.
                !-------------------------------------------------------------------------
                nsize=MCMCgetParticlesInFGNeighborhood(labels_3d,index_B,Nm,Parts_Qb_AgivenB,Parts_Qb_AgivenB_coords)

                IF (MCMCuseBiasedProposal) THEN
                   Normalizer_Qb_A_B=SUM(Parts_Qb_AgivenB(:)%proposal)
                   qb_A_B_unnorm=MCMCproposal(labels_3d,index_A)

                   MCMC_qb_A_B(1)=REAL(qb_A_B_unnorm,ppm_kind_double)/REAL(Normalizer_Qb_A_B,ppm_kind_double)
                ELSE
                   MCMC_qb_A_B(1)=oned/REAL(nsize,ppm_kind_double)
                ENDIF

                DEALLOCATE(Parts_Qb_AgivenB,Parts_Qb_AgivenB_coords,STAT=info)
                or_fail_dealloc("Parts_Qb_AgivenB,Parts_Qb_AgivenB_coords")

                IF (.NOT.MCMC_SingleParticleMoveForPairProposals) THEN
                   !-------------------------------------------------------------------------
                   ! undo the simulated move.
                   !-------------------------------------------------------------------------
                   info=MCMCApplyParticle_3d(image_3d,labels_3d,index_A,Nm,RPart_A,.TRUE.)
                   or_fail("MCMCApplyParticle")

                   nx=INT(REAL(index_B(1),MK)/MCMCcellsize(1))
                   IF (nx.EQ.MCMCcellNm(1)) nx=nx-1

                   ny=INT(REAL(index_B(2),MK)/MCMCcellsize(2))
                   IF (ny.EQ.MCMCcellNm(2)) ny=ny-1

                   nz=INT(REAL(index_B(3),MK)/MCMCcellsize(3))
                   IF (nz.EQ.MCMCcellNm(3)) nz=nz-1

                   cbox=1+nx+MCMCcellNm(1)*(ny+nz*MCMCcellNm(2))

                   !-------------------------------------------------------------------------
                   ! Now we can calculate Q_B (as we now know B and the original
                   ! state has been recovered).
                   !-------------------------------------------------------------------------
                   IF (MCMCIsRegularParticle(index_B(1),index_B(2),index_B(3),cbox,Part_B)) THEN
                      IF (MCMCuseBiasedProposal) THEN
                         !-------------------------------------------------------------------------
                         ! B Proposal needs to be recalculated here:
                         !-------------------------------------------------------------------------
                         MCMC_q_B(1)=REAL(MCMCproposal(labels_3d,index_B),ppm_kind_double)/(MCMCTotalNormalizer+MCMCTotalNormalizerlocal)
                      ELSE
                         MCMC_q_B(1)=oned/(MCMCTotalNormalizer+MCMCTotalNormalizerlocal)
                      ENDIF
                   ELSE
                      MCMC_q_B(1)=zerod
                   ENDIF
                ENDIF

                IF (MCMCuseBiasedProposal) THEN
                   DEALLOCATE(ProposalsVector,STAT=info)
                   or_fail_dealloc("ProposalsVector")

                   !-------------------------------------------------------------------------
                   ! make sure that there is no discrete distribution available in memory
                   !-------------------------------------------------------------------------
                   info=ppm_rc_DestroyParticlesDiscrDistr()
                   or_fail("ppm_rc_DestroyParticlesDiscrDistr")
                ENDIF

                sbpitr => mesh%subpatch%next()
             ENDDO
             NULLIFY(labels_3d,image_3d)
          END SELECT
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        END FUNCTION MCMCFindParticleB












        LOGICAL FUNCTION MCMCMove(ppmrcdim,info)
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_fatal,ppm_error_error
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_error,ppm_err_sub_failed, &
          &   ppm_err_dealloc
          USE ppm_module_mesh_typedef, ONLY : ppm_t_subpatch_

          USE ppm_rc_module_global, ONLY : one,mesh,image,labels,MCMCstepsize, &
          &   MCMCuseBiasedProposal,MCMCusePairProposal
          USE ppm_rc_module_rnd, ONLY : ppm_rc_GetPartDistrIndex,ppm_rc_Saru_IPRNG,  &
          &   ppm_rc_DestroyParticlesDiscrDistr,ppm_rc_GenerateParticlesFwdProposalsDiscrDistr
          USE ppm_rc_module_linkedlist, ONLY : MCMCAppliedParticleOrigLabels, &
          &   ppm_rc_link
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER, INTENT(IN   ) :: ppmrcdim
          !!! Problem dimension

          INTEGER, INTENT(  OUT) :: info
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(MCMCParticle), DIMENSION(:), ALLOCATABLE :: Parts_Q_AgivenB
          TYPE(MCMCParticle), DIMENSION(:), ALLOCATABLE :: Parts_Qb_BgivenA
          TYPE(MCMCParticle)                            :: Part_A
          TYPE(MCMCParticle)                            :: Part_B
          TYPE(MCMCParticle)                            :: ReverseFloatingP
!
          CLASS(ppm_t_subpatch_), POINTER :: sbpitr
!
! !           TYPE(ppm_rc_link), POINTER :: seedlink
!
          REAL(MK), CONTIGUOUS, DIMENSION(:,:),   POINTER     :: image_2d
          REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER     :: image_3d
          REAL(ppm_kind_double)                               :: t0
          REAL(MK)                                            :: Normalizer_Q_A_B
          REAL(MK)                                            :: Normalizer_Qb_B_A
          REAL(MK)                                            :: Normalizer_qb_A
          REAL(MK)                                            :: Normalizer_qb_B
!           REAL(MK)                                            :: ForwardBackwardRatio
          REAL(MK)                                            :: etemp
          REAL(MK)                                            :: proposalA
          REAL(MK)                                            :: proposalBb
!           REAL(ppm_kind_double)                               :: ProposeAFloatInXb
!           REAL(MK)                                            :: HastingsRatio

          INTEGER,             DIMENSION(:,:),   ALLOCATABLE :: Parts_Q_AgivenB_coords
          INTEGER,             DIMENSION(:,:),   ALLOCATABLE :: Parts_Qb_BgivenA_coords
          INTEGER, CONTIGUOUS, DIMENSION(:,:),   POINTER     :: labels_2d
          INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER     :: labels_3d
          INTEGER,             DIMENSION(:),     POINTER     :: Nm
!           INTEGER,             DIMENSION(:),     POINTER     :: seedn
          INTEGER,             DIMENSION(ppmrcdim)           :: index_A
          INTEGER,             DIMENSION(ppmrcdim)           :: index_B
          INTEGER,             DIMENSION(ppmrcdim)           :: index_
          INTEGER                                            :: nsize
          INTEGER                                            :: PC
          INTEGER                                            :: OriginalLabel

          LOGICAL :: e_merge
!           LOGICAL :: Accept

          CHARACTER(LEN=ppm_char) :: caller='MCMCMoves'
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          e_merge=.FALSE.

          SELECT CASE (ppmrcdim)
          CASE (2)
             NULLIFY(labels_2d,image_2d)
             sbpitr => mesh%subpatch%begin()
             DO WHILE (ASSOCIATED(sbpitr))
                Nm => sbpitr%nnodes

                CALL sbpitr%get_field(image,image_2d,info)
                or_fail("Failed to get field image_2d data.",ppm_error=ppm_error_fatal)

                CALL sbpitr%get_field(labels,labels_2d,info)
                or_fail("Failed to get field labels_2d data.",ppm_error=ppm_error_fatal)

                IF (MCMCusePairProposal.AND..NOT.MCMC_SingleParticleMoveForPairProposals) THEN
                   !-------------------------------------------------------------------------
                   ! Iterate the candidates, calculate the energy and perform the moves.
                   !-------------------------------------------------------------------------
                   DO PC=1,MCMCstepsize
                      !-------------------------------------------------------------------------
                      ! apply particle A and B, start with B:
                      ! it is necessary that we start with particle B as we have
                      ! to calculate Q(A|B) and Qb(B|A).
                      !-------------------------------------------------------------------------
                      Part_B =MCMC_PartnerMove(PC)
                      OriginalLabel=MCMC_LabelsBeforeJump_B(PC)
                      index_B=MCMC_PartnerMove_Index(:,PC)
                      !-------------------------------------------------------------------------
                      ! We calculate the energy and apply them move if
                      ! - the move has not been performed beforehand (a particle
                      !   was sampled twice)
                      ! - THE FOLLOWING IS CURRENTLY A NON ISSUE AS B CANNOT BE A':
                      !   particle B is not the reverse particle of particle A (in
                      !   case of m_usePairProposals, i.e. vN > 1). This is important
                      !   because the energy update gets corrupted as particle B is
                      !   not a valid particle before A was applied. To resolve this
                      !   we just perform a one particle move (we don't apply B).
                      !-------------------------------------------------------------------------
                      IF (.NOT.MCMCAppliedParticles%containselement(index_B(1),index_B(2),Part_B)) THEN
                         !-------------------------------------------------------------------------
                         ! Calculate the energy difference when changing this candidate:
                         !-------------------------------------------------------------------------
                         etemp=e_data%EvaluateEnergyDifference(image_2d,labels_2d,   &
                         &     index_B,OriginalLabel,Part_B%candlabel,e_merge)       &
                         &    +e_length%EvaluateEnergyDifference(image_2d,labels_2d, &
                         &     index_B,OriginalLabel,Part_B%candlabel,e_merge)

                         MCMCTotEnergyDiff=MCMCTotEnergyDiff+etemp

                         !-------------------------------------------------------------------------
                         ! Finally, perform the (particle-)move
                         !-------------------------------------------------------------------------
                         info=MCMCApplyParticle_2d(image_2d,index_B,Nm,labels_2d,Part_B,.FALSE.)
                         or_fail("MCMCApplyParticle!")

                         CALL MCMCAppliedParticles%insert(index_B(1),index_B(2),Part_B,info)
                         or_fail("MCMCAppliedParticles%insert")

                         CALL MCMCAppliedParticleOrigLabels%add(index_B(1),index_B(2),OriginalLabel)
                      ENDIF
                      !-------------------------------------------------------------------------
                      ! Calculate Q(A|B) and Qb(B|A) in case we moved B only; this is
                      ! when vN == 2.
                      ! Get the neighbors (conditional particles) and sum up
                      ! their proposal values; this is the normalizer for the
                      ! discrete probability Q(A|B)
                      !-------------------------------------------------------------------------
                      nsize=MCMCgetParticlesInFGNeighborhood(labels_2d,index_B,Nm, &
                      &     Parts_Q_AgivenB,Parts_Q_AgivenB_coords)
                      !-------------------------------------------------------------------------
                      ! add particle B as this is always a candidate as well
                      !-------------------------------------------------------------------------

                      index_A=MCMC_CandidateMove_Index(:,PC)
                      !-------------------------------------------------------------------------
                      ! Part_A%Proposal is not valid anymore. Part_A
                      ! got a new proposal when applying particle B.
                      !-------------------------------------------------------------------------
                      proposalA=MCMCproposal(labels_2d,index_A)

                      IF (MCMCuseBiasedProposal) THEN
                         Normalizer_Q_A_B=SUM(Parts_Q_AgivenB(:)%proposal)
                         MCMC_q_A_B(PC)=proposalA/Normalizer_Q_A_B
                      ELSE
                         MCMC_q_A_B(PC)=one/REAL(nsize,MK)
                      ENDIF !MCMCuseBiasedProposal

                      !-------------------------------------------------------------------------
                      ! create A'
                      ! Calculate Qb(B'|A')
                      !-------------------------------------------------------------------------
                      nsize=MCMCgetParticlesInFGNeighborhood(labels_2d,index_A,Nm, &
                      &     Parts_Qb_BgivenA,Parts_Qb_BgivenA_coords)

                      IF (MCMCuseBiasedProposal) THEN
                         Normalizer_Qb_B_A=SUM(Parts_Qb_BgivenA(:)%proposal)

                         proposalBb=MCMCproposal(labels_2d,index_B)

                         MCMC_qb_B_A(PC)=proposalBb/Normalizer_Qb_B_A
                      ELSE
                         MCMC_qb_B_A(PC)=one/REAL(nsize,MK)
                      ENDIF !MCMCuseBiasedProposal

                      ! apply particle A
                      Part_A =MCMC_CandidateMove(PC)
                      OriginalLabel=MCMC_LabelsBeforeJump_A(PC)
                      index_A=MCMC_CandidateMove_Index(:,PC)

                      !-------------------------------------------------------------------------
                      ! We calculate the energy and apply them move if
                      ! - the move has not been performed beforehand (a particle
                      !   was sampled twice)
                      ! - THE FOLLOWING IS CURRENTLY A NON ISSUE AS B CANNOT BE A':
                      !   particle B is not the reverse particle of particle A (in
                      !   case of m_usePairProposals, i.e. vN > 1). This is important
                      !   because the energy update gets corrupted as particle B is
                      !   not a valid particle before A was applied. To resolve this
                      !   we just perform a one particle move (we don't apply B).
                      !-------------------------------------------------------------------------
                      IF (.NOT.MCMCAppliedParticles%containselement(index_A(1),index_A(2),Part_A)) THEN
                         ! Calculate the energy difference when changing this candidate:
                         etemp=e_data%EvaluateEnergyDifference(image_2d,labels_2d,   &
                         &     index_A,OriginalLabel,Part_A%candlabel,e_merge)       &
                         &    +e_length%EvaluateEnergyDifference(image_2d,labels_2d, &
                         &     index_A,OriginalLabel,Part_A%candlabel,e_merge)

                         MCMCTotEnergyDiff=MCMCTotEnergyDiff+etemp
                         !-------------------------------------------------------------------------
                         ! Finally, perform the (particle-)move
                         !-------------------------------------------------------------------------
                         info=MCMCApplyParticle_2d(image_2d,index_A,Nm,labels_2d,Part_A,.FALSE.)
                         or_fail("MCMCApplyParticle!")

                         CALL MCMCAppliedParticles%insert(index_A(1),index_A(2),Part_A,info)
                         or_fail("MCMCAppliedParticles%insert")

                         CALL MCMCAppliedParticleOrigLabels%add(index_A(1),index_A(2),OriginalLabel)
                      ENDIF
                   ENDDO !PC=1,MCMCstepsize
                ELSE
                   !-------------------------------------------------------------------------
                   ! Iterate the candidates, calculate the energy and perform the moves.
                   !-------------------------------------------------------------------------
                   DO PC=1,MCMCstepsize
                      Part_A = MCMC_CandidateMove(PC)
                      !-------------------------------------------------------------------------
                      ! apply particle A and B, start with B:
                      ! it is necessary that we start with particle B as we have
                      ! to calculate Q(A|B) and Qb(B|A).
                      !-------------------------------------------------------------------------
                      OriginalLabel=MCMC_LabelsBeforeJump_A(PC)
                      index_A=MCMC_CandidateMove_Index(:,PC)

                      !-------------------------------------------------------------------------
                      ! We calculate the energy and apply them move if
                      ! - the move has not been performed beforehand (a particle
                      !   was sampled twice)
                      ! - THE FOLLOWING IS CURRENTLY A NON ISSUE AS B CANNOT BE A':
                      !   particle B is not the reverse particle of particle A (in
                      !   case of m_usePairProposals, i.e. vN > 1). This is important
                      !   because the energy update gets corrupted as particle B is
                      !   not a valid particle before A was applied. To resolve this
                      !   we just perform a one particle move (we don't apply B).
                      !-------------------------------------------------------------------------
                      IF (.NOT.MCMCAppliedParticles%containselement(index_A(1),index_A(2),Part_A)) THEN
                         ! Calculate the energy difference when changing this candidate:
                         etemp=e_data%EvaluateEnergyDifference(image_2d,labels_2d,   &
                         &     index_A,OriginalLabel,Part_A%candlabel,e_merge)       &
                         &    +e_length%EvaluateEnergyDifference(image_2d,labels_2d, &
                         &     index_A,OriginalLabel,Part_A%candlabel,e_merge)

                         MCMCTotEnergyDiff=MCMCTotEnergyDiff+etemp
                         !-------------------------------------------------------------------------
                         ! Finally, perform the (particle-)move
                         !-------------------------------------------------------------------------
                         info=MCMCApplyParticle_2d(image_2d,index_A,Nm,labels_2d,Part_A,.FALSE.)
                         or_fail("MCMCApplyParticle!")

                         CALL MCMCAppliedParticles%insert(index_A(1),index_A(2),Part_A,info)
                         or_fail("MCMCAppliedParticles%insert")

                         CALL MCMCAppliedParticleOrigLabels%add(index_A(1),index_A(2),OriginalLabel)
                      ENDIF
                   ENDDO !PC=1,MCMCstepsize
                ENDIF !MCMCusePairProposal.AND..NOT.MCMC_SingleParticleMoveForPairProposals

                !-------------------------------------------------------------------------
                ! Free memory
                !-------------------------------------------------------------------------
                IF (ALLOCATED(Parts_Q_AgivenB)) THEN
                   DEALLOCATE(Parts_Q_AgivenB,Parts_Q_AgivenB_coords,STAT=info)
                   or_fail_dealloc("Parts_Q_AgivenB & Parts_Q_AgivenB_coords")
                ENDIF
                IF (ALLOCATED(Parts_Qb_BgivenA)) THEN
                   DEALLOCATE(Parts_Qb_BgivenA,Parts_Qb_BgivenA_coords,STAT=info)
                   or_fail_dealloc("Parts_Qb_BgivenA & Parts_Qb_BgivenA_coords")
                ENDIF

                DO PC=1,MCMCstepsize
                   !-------------------------------------------------------------------------
                   ! Correct the containers whenever floating particles were involved:
                   ! The method moveParticles, for simplicity, only works on the regular
                   ! particle set.
                   !-------------------------------------------------------------------------

                   !-------------------------------------------------------------------------
                   ! First, figure out if A' or B' is floating:
                   ! Figure out if the backward particles are floating:
                   !-------------------------------------------------------------------------
                   IF (MCMCusePairProposal) THEN
                      index_=MCMC_PartnerMove_Index(:,PC)

                      MCMC_Particle_Bb_IsFloating(PC)= &
                      & MCMCParticleHasFloatingProperty(labels_2d,index_,MCMC_LabelsBeforeJump_B(PC))
                   ELSE
                      index_=MCMC_CandidateMove_Index(:,PC)
                      !-------------------------------------------------------------------------
                      ! if we're not in pair proposal mode we did not yet check if
                      ! A's reverse particle is floating (else we did already):
                      !-------------------------------------------------------------------------
                      MCMC_Particle_Ab_IsFloating(PC)= &
                      & MCMCParticleHasFloatingProperty(labels_2d,index_,MCMC_LabelsBeforeJump_A(PC))
                   ENDIF

                   !-------------------------------------------------------------------------
                   ! the first condition is needed when not using pair proposal mode
                   !-------------------------------------------------------------------------
                   IF (MCMC_Particle_Ab_IsFloating(PC)) THEN
                      ReverseFloatingP=MCMCParticle(MCMC_LabelsBeforeJump_A(PC),MCMCproposal(labels_2d,index_))
                   ENDIF

                   !-------------------------------------------------------------------------
                   ! in pair proposal, if A' is floating, B' is as well (they are the same particle):
                   !-------------------------------------------------------------------------
                   IF (MCMCusePairProposal.AND.MCMC_Particle_Bb_IsFloating(PC)) THEN
                      !-------------------------------------------------------------------------
                      ! only possible in pair proposal mode
                      !-------------------------------------------------------------------------
                      ReverseFloatingP=MCMCParticle(MCMC_LabelsBeforeJump_B(PC),MCMCproposal(labels_2d,index_))
                   ENDIF

                   !-------------------------------------------------------------------------
                   ! finally convert the regular particle into a floating particle,
                   ! i.e. insert it in the floating DS and remove it from the regular:
                   !-------------------------------------------------------------------------
                   IF (MCMC_Particle_Ab_IsFloating(PC).OR.MCMC_Particle_Bb_IsFloating(PC)) THEN
                      !-------------------------------------------------------------------------
                      ! insert the reverse particle in the appropriate container. If
                      ! there is no space, we reject the move.
                      !-------------------------------------------------------------------------
                      IF (.NOT.MCMCInsertFloatingParticle(index_(1),index_(2),ReverseFloatingP,.TRUE.)) THEN
                         !-------------------------------------------------------------------------
                         ! TODO: calling MCMCReject from here invalidates the result.
                         ! MCMCReject(&vAppliedParticles,&vAppliedParticleOrigLabels);
                         !-------------------------------------------------------------------------
                         HardReject=.TRUE.
                      ENDIF
                   ENDIF
                ENDDO !PC=1,MCMCstepsize

                !-------------------------------------------------------------------------
                ! We are now in the state x'.
                ! Calculate Q'(A) and maybe Q'(B). Note that this has to be done after
                ! all particles were applied.
                !-------------------------------------------------------------------------
                DO PC=1,MCMCstepsize
                   !-------------------------------------------------------------------------
                   ! Calculate MCMC_qb_A and MCMC_qb_B
                   !-------------------------------------------------------------------------
                   IF (.NOT.MCMCuseBiasedProposal) THEN
                      MCMC_qb_A(PC)=one
                      MCMC_qb_B(PC)=one
                   ELSE
                      index_A=MCMC_CandidateMove_Index(:,PC)
                      MCMC_qb_A(PC)=MCMCproposal(labels_2d,index_A)
                      IF (MCMCusePairProposal.AND..NOT.MCMC_SingleParticleMoveForPairProposals) THEN
                         index_B=MCMC_PartnerMove_Index(:,PC)
                         MCMC_qb_B(PC)=MCMCproposal(labels_2d,index_B)
                      ENDIF
                   ENDIF

                   !-------------------------------------------------------------------------
                   ! Normalize MCMC_qb_A and MCMC_qb_B
                   !-------------------------------------------------------------------------
                   IF (MCMC_Particle_Ab_IsFloating(PC)) THEN
                      Normalizer_qb_A=REAL(MCMCFloatingParticlesProposalNormalizer,MK)
                   ELSE
                      Normalizer_qb_A=MCMCGetProposalNormalizer(MCMC_CandidateMove(PC)%candlabel,MCMC_LabelsBeforeJump_A(PC))
                   ENDIF
                   MCMC_qb_A(PC)=MCMC_qb_A(PC)/Normalizer_qb_A

                   IF (MCMCusePairProposal.AND..NOT.MCMC_SingleParticleMoveForPairProposals) THEN
                      IF (MCMC_Particle_Bb_IsFloating(PC)) THEN
                         Normalizer_qb_B=REAL(MCMCFloatingParticlesProposalNormalizer,MK)
                      ELSE
                         Normalizer_qb_B=MCMCGetProposalNormalizer(MCMC_PartnerMove(PC)%candlabel,MCMC_LabelsBeforeJump_B(PC))
                      ENDIF
                      MCMC_qb_B(PC)=MCMC_qb_B(PC)/Normalizer_qb_B
                   ENDIF
                   !-------------------------------------------------------------------------
                   ! Finally, we omit half of the calculations if particle A == B
                   !-------------------------------------------------------------------------
                   IF (MCMC_SingleParticleMoveForPairProposals) THEN
                      MCMC_q_B(PC)   =MCMC_q_A(PC)
                      MCMC_q_A_B(PC) =MCMC_q_B_A(PC);
                      MCMC_qb_B_A(PC)=MCMC_qb_A_B(PC);
                      MCMC_qb_B(PC)  =MCMC_qb_A(PC);
                   ENDIF
                ENDDO !PC=1,MCMCstepsize

!                 ! Calculate the forward-backward ratio:
!                 ForwardBackwardRatio=one
!                 DO PC=1,MCMCstepsize
!                    IF (MCMCusePairProposal) THEN
!                       IF (MCMC_Particle_Ab_IsFloating(PC).OR.MCMC_Particle_Bb_IsFloating(PC).OR.MCMC_Particle_A_IsFloating(PC)) THEN
!                          ForwardBackwardRatio=ForwardBackwardRatio* &
!                          & (MCMC_qb_B(PC)*MCMC_qb_A_B(PC))/(MCMC_q_A(PC)*MCMC_q_B_A(PC))
!                       ELSE
!                          ForwardBackwardRatio=ForwardBackwardRatio* &
!                          & (MCMC_qb_A(PC)*MCMC_qb_B_A(PC)+MCMC_qb_B(PC)*MCMC_qb_A_B(PC))/(MCMC_q_A(PC)*MCMC_q_B_A(PC)+MCMC_q_B(PC)*MCMC_q_A_B(PC))
!                       ENDIF
!                    ELSE
!                       IF (MCMC_Particle_A_IsFloating(PC)) THEN
!                          ! we distroy a floating particle, in the next iteration there
!                          ! will be one floating particle less, hence the probability
!                          ! in the x' to sample a floating particle is (note that
!                          ! both normalizers are in state x'):
!                          ProposeAFloatInXb=MCMCTotalNormalizer/(MCMCFloatingParticlesProposalNormalizer+MCMCTotalNormalizer)
!                          ProposeAFloatInXb=halfd*ProposeAFloatInXb/MCMCProbabilityToProposeAFloatingParticle
!
!                          ForwardBackwardRatio=ForwardBackwardRatio*REAL(ProposeAFloatInXb,MK)*MCMC_qb_A(PC)/MCMC_q_A(PC)
!                       ELSE IF (MCMC_Particle_Ab_IsFloating(PC)) THEN
!                          ! we create a floating particle, in the next iteration there
!                          ! will be one floating particle more, hence the probability
!                          ! in the x' to sample a floating particle is (note that
!                          ! m_MCMCTotalNormalizer is updated to x'):
!                          ProposeAFloatInXb=MCMCFloatingParticlesProposalNormalizer/(MCMCFloatingParticlesProposalNormalizer+MCMCTotalNormalizer)
!                          ProposeAFloatInXb=ProposeAFloatInXb/(halfd*(oned-MCMCProbabilityToProposeAFloatingParticle))
!
!                          ForwardBackwardRatio=ForwardBackwardRatio*REAL(ProposeAFloatInXb,MK)*MCMC_qb_A(PC)/MCMC_q_A(PC)
!                       ELSE
!                          ! Shrinkage and growth events have the same probability.
!                          ! We hence only need to compare the individual particle
!                          ! ratios.
!                          ForwardBackwardRatio=ForwardBackwardRatio*MCMC_qb_A(PC)/MCMC_q_A(PC)
!                       ENDIF
!                    ENDIF
!                 ENDDO !PC=1,MCMCstepsize
!
!                 ! Compute the Hastingsratio:
!                 HastingsRatio = EXP(-MCMCTotEnergyDiff/MCMCtemperature)*ForwardBackwardRatio
!                 ! debug: check if the average of the FBR is equal to one.
!                 ! HastingsRatio = ForwardBackwardRatio
!                 ! Should I stay or shoud I go; the Metropolis-Hastings algorithm:
!                 IF (HastingsRatio.GE.one) THEN
!                    Accept=.TRUE.
!                 ELSE
!                    IF (HastingsRatio.GT.ppm_rc_Saru_RPRNG()) THEN
!                       Accept=.TRUE.
!                    ELSE
!                       Accept=.FALSE.
!                    ENDIF
!                 ENDIF
!
!                 ! Register the result (if we accept) or rollback to the previous state.
!                 IF (Accept.AND..NOT.HardReject) THEN
!                    seedlink => MCMCAppliedParticleOrigLabels%first
!                    DO WHILE (ASSOCIATED(seedlink))
!                       seedn => seedlnk%getValue()
!                       ! store the results and finish (next iteration).
!                       Part_A=MCMCAppliedParticles%search(seedn(1),seedn(2))
!                       !MCMCStore index,originallabel,step
!                       ! something changed, particles were inserted and deleted. We need
!                       ! to update the edge-map constant.
!                       !MCMCUpdateRegularParticleMapInNeighborhood(vAppliedParticles[vM].m_Index);
!                       seedlnk => seedlnk%nextLink()
!                    ENDDO
!                 ELSE
!                    !MCMCReject(&vAppliedParticles,&vAppliedParticleOrigLabels);
!                 ENDIF
!
!                 MCMCMove=Accept

                sbpitr => mesh%subpatch%next()
             ENDDO !ASSOCIATED(sbpitr)
             NULLIFY(labels_2d,image_2d)
          CASE (3)
             NULLIFY(labels_3d,image_3d)
             sbpitr => mesh%subpatch%begin()
             DO WHILE (ASSOCIATED(sbpitr))
                CALL sbpitr%get_field(labels,labels_3d,info)
                or_fail("Failed to get field labels_3d data.",ppm_error=ppm_error_fatal)

                CALL sbpitr%get_field(image,image_3d,info)
                or_fail("Failed to get field labels_3d data.",ppm_error=ppm_error_fatal)

                IF (MCMCusePairProposal.AND..NOT.MCMC_SingleParticleMoveForPairProposals) THEN
                   !-------------------------------------------------------------------------
                   ! Iterate the candidates, calculate the energy and perform the moves.
                   !-------------------------------------------------------------------------
                   DO PC=1,MCMCstepsize
                      !-------------------------------------------------------------------------
                      ! apply particle A and B, start with B:
                      !-------------------------------------------------------------------------
                      Part_B=MCMC_PartnerMove(PC)
                      !-------------------------------------------------------------------------
                      ! it is necessary that we start with particle B as we have
                      ! to calculate Q(A|B) and Qb(B|A).
                      !-------------------------------------------------------------------------
                      OriginalLabel=MCMC_LabelsBeforeJump_B(PC)
                      index_B=MCMC_PartnerMove_Index(:,PC)
                      !-------------------------------------------------------------------------
                      ! We calculate the energy and apply the move if
                      ! - the move has not been performed beforehand (a particle
                      !   was sampled twice)
                      ! - THE FOLLOWING IS CURRENTLY A NON ISSUE AS B CANNOT BE A':
                      !   particle B is not the reverse particle of particle A (in
                      !   case of m_usePairProposals, i.e. vN > 1). This is important
                      !   because the energy update gets corrupted as particle B is
                      !   not a valid particle before A was applied. To resolve this
                      !   we just perform a one particle move (we don't apply B).
                      !-------------------------------------------------------------------------
                      IF (.NOT.MCMCAppliedParticles%containselement(index_B(1),index_B(2),index_B(3),Part_B)) THEN
                         !-------------------------------------------------------------------------
                         ! Calculate the energy difference when changing this candidate:
                         !-------------------------------------------------------------------------
                         etemp=e_data%EvaluateEnergyDifference(image_3d,labels_3d,   &
                         &     index_B,OriginalLabel,Part_B%candlabel,e_merge)       &
                         &    +e_length%EvaluateEnergyDifference(image_3d,labels_3d, &
                         &     index_B,OriginalLabel,Part_B%candlabel,e_merge)

                         MCMCTotEnergyDiff=MCMCTotEnergyDiff+etemp
                         !-------------------------------------------------------------------------
                         ! Finally, perform the (particle-)move
                         !-------------------------------------------------------------------------
                         info=MCMCApplyParticle_3d(image_3d,index_B,Nm,labels_3d,Part_B,.FALSE.)
                         or_fail("MCMCApplyParticle")

                         CALL MCMCAppliedParticles%insert(index_B(1),index_B(2),index_B(3),Part_B,info)
                         or_fail("MCMCAppliedParticles%insert")

                         CALL MCMCAppliedParticleOrigLabels%add(index_B(1),index_B(2),index_B(3),OriginalLabel)
                      ENDIF
                      !-------------------------------------------------------------------------
                      ! Calculate Q(A|B) and Qb(B|A) in case we moved B only;
                      ! Get the neighbors (conditional particles) and sum up
                      ! their proposal values; this is the normalizer for the
                      ! discrete probability Q(A|B)
                      !-------------------------------------------------------------------------
                      nsize=MCMCgetParticlesInFGNeighborhood(labels_3d,index_B,Nm, &
                      &     Parts_Q_AgivenB,Parts_Q_AgivenB_coords)
                      Part_B%proposal=Parts_Q_AgivenB(nsize)%proposal

                      index_A=MCMC_CandidateMove_Index(:,PC)
                      !-------------------------------------------------------------------------
                      ! Part_A%proposal is not valid anymore. Particle A
                      ! got a new proposal when applying particle B.
                      !-------------------------------------------------------------------------
                      proposalA=MCMCproposal(labels_3d,index_A)

                      IF (MCMCuseBiasedProposal) THEN
                         Normalizer_Q_A_B=SUM(Parts_Q_AgivenB(:)%proposal)
                         MCMC_q_A_B(PC)=proposalA/Normalizer_Q_A_B
                      ELSE
                         MCMC_q_A_B(PC)=one/REAL(nsize,MK)
                      ENDIF

                      !-------------------------------------------------------------------------
                      ! create A'
                      ! Calculate Qb(B'|A')
                      !-------------------------------------------------------------------------
                      nsize=MCMCgetParticlesInFGNeighborhood(labels_3d,index_A,Nm, &
                      &     Parts_Qb_BgivenA,Parts_Qb_BgivenA_coords)

                      IF (MCMCuseBiasedProposal) THEN
                         Normalizer_Qb_B_A=SUM(Parts_Qb_BgivenA(:)%proposal)
                         !-------------------------------------------------------------------------
                         ! the proposal of the backward particle (given A) is:
                         !-------------------------------------------------------------------------
                         proposalBb=MCMCproposal(labels_3d,index_B)

                         MCMC_qb_B_A(PC)=proposalBb/Normalizer_Qb_B_A
                      ELSE
                         MCMC_qb_B_A(PC)=one/REAL(nsize,MK)
                      ENDIF

                      !-------------------------------------------------------------------------
                      ! apply particle A
                      !-------------------------------------------------------------------------
                      Part_A=MCMC_CandidateMove(PC)
                      OriginalLabel=MCMC_LabelsBeforeJump_A(PC)
                      index_A=MCMC_CandidateMove_Index(:,PC)
                      !-------------------------------------------------------------------------
                      ! We calculate the energy and apply them move iff
                      ! - the move has not been performed beforehand (a particle
                      !   was sampled twice)
                      ! - THE FOLLOWING IS CURRENTLY A NON ISSUE AS B CANNOT BE A':
                      !   particle B is not the reverse particle of particle A (in
                      !   case of m_usePairProposals, i.e. vN > 1). This is important
                      !   because the energy update gets corrupted as particle B is
                      !   not a valid particle before A was applied. To resolve this
                      !   we just perform a one particle move (we don't apply B).
                      !-------------------------------------------------------------------------
                      IF (.NOT.MCMCAppliedParticles%containselement(index_A(1),index_A(2),index_A(3),Part_A)) THEN
                         !-------------------------------------------------------------------------
                         ! Calculate the energy difference when changing this candidate:
                         !-------------------------------------------------------------------------
                         etemp=e_data%EvaluateEnergyDifference(image_3d,labels_3d,   &
                         &     index_A,OriginalLabel,Part_A%candlabel,e_merge)       &
                         &    +e_length%EvaluateEnergyDifference(image_3d,labels_3d, &
                         &     index_A,OriginalLabel,Part_A%candlabel,e_merge)

                         MCMCTotEnergyDiff=MCMCTotEnergyDiff+etemp
                         !-------------------------------------------------------------------------
                         ! Finally, perform the (particle-)move
                         !-------------------------------------------------------------------------
                         info=MCMCApplyParticle_3d(image_3d,index_A,Nm,labels_3d,Part_A,.FALSE.)
                         or_fail("MCMCApplyParticle")

                         CALL MCMCAppliedParticles%insert(index_A(1),index_A(2),index_A(3),Part_A,info)
                         or_fail("MCMCAppliedParticles%insert")

                         CALL MCMCAppliedParticleOrigLabels%add(index_A(1),index_A(2),index_A(3),OriginalLabel)
                      ENDIF
                   ENDDO !PC=1,MCMCstepsize
                ELSE
                   !-------------------------------------------------------------------------
                   ! Iterate the candidates, calculate the energy and perform the moves.
                   !-------------------------------------------------------------------------
                   DO PC=1,MCMCstepsize
                      !-------------------------------------------------------------------------
                      ! apply particle A:
                      !-------------------------------------------------------------------------
                      Part_A=MCMC_CandidateMove(PC)
                      OriginalLabel=MCMC_LabelsBeforeJump_A(PC)
                      index_A=MCMC_CandidateMove_Index(:,PC)
                      !-------------------------------------------------------------------------
                      ! We calculate the energy and apply the move if
                      ! - the move has not been performed beforehand (a particle
                      !   was sampled twice)
                      ! - THE FOLLOWING IS CURRENTLY A NON ISSUE AS B CANNOT BE A':
                      !   particle B is not the reverse particle of particle A (in
                      !   case of m_usePairProposals, i.e. vN > 1). This is important
                      !   because the energy update gets corrupted as particle B is
                      !   not a valid particle before A was applied. To resolve this
                      !   we just perform a one particle move (we don't apply B).
                      !-------------------------------------------------------------------------
                      IF (.NOT.MCMCAppliedParticles%containselement(index_A(1),index_A(2),index_A(3),Part_A)) THEN
                         !-------------------------------------------------------------------------
                         ! Calculate the energy difference when changing this candidate:
                         !-------------------------------------------------------------------------
                         etemp=e_data%EvaluateEnergyDifference(image_3d,labels_3d,   &
                         &     index_A,OriginalLabel,Part_A%candlabel,e_merge)       &
                         &    +e_length%EvaluateEnergyDifference(image_3d,labels_3d, &
                         &     index_A,OriginalLabel,Part_A%candlabel,e_merge)

                         MCMCTotEnergyDiff=MCMCTotEnergyDiff+etemp
                         !-------------------------------------------------------------------------
                         ! Finally, perform the (particle-)move
                         !-------------------------------------------------------------------------
                         info=MCMCApplyParticle_3d(image_3d,index_A,Nm,labels_3d,Part_A,.FALSE.)
                         or_fail("MCMCApplyParticle")

                         CALL MCMCAppliedParticles%insert(index_A(1),index_A(2),index_A(3),Part_A,info)
                         or_fail("MCMCAppliedParticles%insert")

                         CALL MCMCAppliedParticleOrigLabels%add(index_A(1),index_A(2),index_A(3),OriginalLabel)
                      ENDIF
                   ENDDO !PC=1,MCMCstepsize
                ENDIF !MCMCusePairProposal.AND..NOT.MCMC_SingleParticleMoveForPairProposals

                ! Free memory
                IF (ALLOCATED(Parts_Q_AgivenB)) THEN
                   DEALLOCATE(Parts_Q_AgivenB,Parts_Q_AgivenB_coords,STAT=info)
                   or_fail_dealloc("Parts_Q_AgivenB & Parts_Q_AgivenB_coords")
                ENDIF
                IF (ALLOCATED(Parts_Qb_BgivenA)) THEN
                   DEALLOCATE(Parts_Qb_BgivenA,Parts_Qb_BgivenA_coords,STAT=info)
                   or_fail_dealloc("Parts_Qb_BgivenA & Parts_Qb_BgivenA_coords")
                ENDIF

                DO PC=1,MCMCstepsize
                   !-------------------------------------------------------------------------
                   ! Correct the containers whenever floating particles were involved:
                   ! The method moveParticles, for simplicity, only works on the regular
                   ! particle set.
                   !-------------------------------------------------------------------------

                   !-------------------------------------------------------------------------
                   ! First, figure out if A' or B' is floating:
                   ! Figure out if the backward particles are floating:
                   !-------------------------------------------------------------------------
                   IF (MCMCusePairProposal) THEN
                      index_=MCMC_PartnerMove_Index(:,PC)

                      MCMC_Particle_Bb_IsFloating(PC)= &
                      & MCMCParticleHasFloatingProperty(labels_3d,index_,MCMC_LabelsBeforeJump_B(PC))
                   ELSE
                      index_=MCMC_CandidateMove_Index(:,PC)
                      !-------------------------------------------------------------------------
                      ! if we're not in pair proposal mode we did not yet check if
                      ! A's reverse particle is floating (else we did already):
                      !-------------------------------------------------------------------------
                      MCMC_Particle_Ab_IsFloating(PC)= &
                      & MCMCParticleHasFloatingProperty(labels_3d,index_,MCMC_LabelsBeforeJump_A(PC))
                   ENDIF

                   !-------------------------------------------------------------------------
                   ! the first condition is needed when not using pair proposal mode
                   !-------------------------------------------------------------------------
                   IF (MCMC_Particle_Ab_IsFloating(PC)) THEN
                      ReverseFloatingP=MCMCParticle(MCMC_LabelsBeforeJump_A(PC),MCMCproposal(labels_3d,index_))
                   ENDIF

                   !-------------------------------------------------------------------------
                   ! in pair proposal, if A' is floating, B' is as well (they are the same particle):
                   !-------------------------------------------------------------------------
                   IF (MCMCusePairProposal.AND.MCMC_Particle_Bb_IsFloating(PC)) THEN
                      !-------------------------------------------------------------------------
                      ! only possible in pair proposal mode
                      !-------------------------------------------------------------------------
                      ReverseFloatingP=MCMCParticle(MCMC_LabelsBeforeJump_B(PC),MCMCproposal(labels_3d,index_))
                   ENDIF

                   !-------------------------------------------------------------------------
                   ! finally convert the regular particle into a floating particle,
                   ! i.e. insert it in the floating DS and remove it from the regular:
                   !-------------------------------------------------------------------------
                   IF (MCMC_Particle_Ab_IsFloating(PC).OR.MCMC_Particle_Bb_IsFloating(PC)) THEN
                      !-------------------------------------------------------------------------
                      ! insert the reverse particle in the appropriate container. If
                      ! there is no space, we reject the move.
                      !-------------------------------------------------------------------------
                      IF (.NOT.MCMCInsertFloatingParticle(index_(1),index_(2),index_(3),ReverseFloatingP,.TRUE.)) THEN
                         !-------------------------------------------------------------------------
                         ! TODO: calling MCMCReject from here invalidates the result.
                         ! MCMCReject(&vAppliedParticles,&vAppliedParticleOrigLabels);
                         !-------------------------------------------------------------------------
                         HardReject=.TRUE.
                      ENDIF
                   ENDIF
                ENDDO !PC=1,MCMCstepsize

                !-------------------------------------------------------------------------
                ! We are now in the state x'.
                ! Calculate Q'(A) and maybe Q'(B). Note that this has to be done after
                ! all particles were applied.
                !-------------------------------------------------------------------------
                DO PC=1,MCMCstepsize
                   !-------------------------------------------------------------------------
                   ! Calculate MCMC_qb_A and MCMC_qb_B
                   !-------------------------------------------------------------------------
                   IF (.NOT.MCMCuseBiasedProposal) THEN
                      MCMC_qb_A(PC)=one
                      MCMC_qb_B(PC)=one
                   ELSE
                      index_A=MCMC_CandidateMove_Index(:,PC)
                      MCMC_qb_A(PC)=MCMCproposal(labels_3d,index_A)
                      IF (MCMCusePairProposal.AND..NOT.MCMC_SingleParticleMoveForPairProposals) THEN
                         index_B=MCMC_PartnerMove_Index(:,PC)
                         MCMC_qb_B(PC)=MCMCproposal(labels_3d,index_B)
                      ENDIF
                   ENDIF

                   !-------------------------------------------------------------------------
                   ! Normalize MCMC_qb_A and MCMC_qb_B
                   !-------------------------------------------------------------------------
                   IF (MCMC_Particle_Ab_IsFloating(PC)) THEN
                      Normalizer_qb_A=REAL(MCMCFloatingParticlesProposalNormalizer,MK)
                   ELSE
                      Normalizer_qb_A=MCMCGetProposalNormalizer(MCMC_CandidateMove(PC)%candlabel,MCMC_LabelsBeforeJump_A(PC))
                   ENDIF
                   MCMC_qb_A(PC)=MCMC_qb_A(PC)/Normalizer_qb_A

                   IF (MCMCusePairProposal.AND..NOT.MCMC_SingleParticleMoveForPairProposals) THEN
                      IF (MCMC_Particle_Bb_IsFloating(PC)) THEN
                         Normalizer_qb_B=REAL(MCMCFloatingParticlesProposalNormalizer,MK)
                      ELSE
                         Normalizer_qb_B=MCMCGetProposalNormalizer(MCMC_PartnerMove(PC)%candlabel,MCMC_LabelsBeforeJump_B(PC))
                      ENDIF
                      MCMC_qb_B(PC)=MCMC_qb_B(PC)/Normalizer_qb_B
                   ENDIF
                   !-------------------------------------------------------------------------
                   ! Finally, we omit half of the calculations if particle A == B
                   !-------------------------------------------------------------------------
                   IF (MCMC_SingleParticleMoveForPairProposals) THEN
                      MCMC_q_B(PC)   =MCMC_q_A(PC)
                      MCMC_q_A_B(PC) =MCMC_q_B_A(PC);
                      MCMC_qb_B_A(PC)=MCMC_qb_A_B(PC);
                      MCMC_qb_B(PC)  =MCMC_qb_A(PC);
                   ENDIF
                ENDDO !PC=1,MCMCstepsize

                sbpitr => mesh%subpatch%next()
             ENDDO
             NULLIFY(labels_3d,image_3d)
          END SELECT

          MCMCMove=.TRUE.
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        END FUNCTION MCMCMove









