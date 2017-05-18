        SUBROUTINE MCMCCheckParameters(info)
          !!! Check the input parameters
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_fatal
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_error,ppm_err_argument

          USE ppm_rc_module_global, ONLY : MCMCstepsize,MCMCusePairProposal, &
          &   MCMCcontinue,nsteql
          IMPLICIT NONE

          INTEGER, INTENT(INOUT) :: info
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          REAL(ppm_kind_double) :: t0

          CHARACTER(LEN=ppm_char) :: caller='MCMCCheckParameters'
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          IF (MCMCstepsize.GT.1.AND.MCMCusePairProposal) THEN
             fail("Using pair potential with step size > 1 leads to slightly wrong proposals and hence the detailed balance is not guaranteed.", &
             & ppm_error=ppm_error_fatal)
          ENDIF

          IF (.NOT.MCMCcontinue.AND.nsteql.GT.0) THEN
             fail("Energy equilibrium state and changing the energy functional does not work with MCMC!", &
             & ppm_error=ppm_error_fatal)
          ENDIF
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        END SUBROUTINE MCMCCheckParameters

        SUBROUTINE MCMCFreeMemory(info)
          !!! Free memroy which is not used anymore
          !!! This subroutine is useful especially if we are continuing from
          !!! RC client for uncertainty quantification
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_error
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_error,ppm_err_argument,ppm_err_dealloc, &
          &   ppm_err_sub_failed
          USE ppm_module_mesh_typedef, ONLY : ppm_t_subpatch_

          USE ppm_rc_module_global, ONLY : pind,plabels,energya,labela,candlabela, &
          &   ccandlabela,ndaughtersa,nmothersa,mothersa,daughtersa,accepteda,     &
          &   processeda,mesh,Part,ineighproc,procflag
          USE ppm_rc_module_linkedlist, ONLY : Candidates,m_Seeds, &
          &   CompetingRegions,ppm_rc_seeds,ppm_rc_seeds_to_remove, &
          &   InnerContourContainer
          USE ppm_rc_module_energy, ONLY : e_data
          IMPLICIT NONE

          INTEGER, INTENT(INOUT) :: info
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          CLASS(ppm_t_subpatch_), POINTER :: sbpitr

          REAL(ppm_kind_double) :: t0

          INTEGER :: ipatch

          CHARACTER(LEN=ppm_char) :: caller='MCMCFreeMemory'
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          !-----------------------------------------------------------------------
          !  Particle index term
          !-----------------------------------------------------------------------
          IF (ASSOCIATED(pind)) THEN
             CALL pind%destroy(info)
             or_fail("pind%destroy")
             DEALLOCATE(pind,STAT=info)
             or_fail_dealloc("pind")
             NULLIFY(pind)
          ENDIF
          !-----------------------------------------------------------------------
          !  Particle label
          !-----------------------------------------------------------------------
          IF (ASSOCIATED(plabels)) THEN
             CALL plabels%destroy(info)
             or_fail("plabels%destroy")
             DEALLOCATE(plabels,STAT=info)
             or_fail_dealloc("plabels")
             NULLIFY(plabels)
          ENDIF
          !---------------------------------------------------------------------
          ! Destroy the particle Part
          !---------------------------------------------------------------------
          IF (ASSOCIATED(Part)) THEN
             CALL Part%destroy(info)
             or_fail("Failed to destroy Part.")

             dealloc_pointer("Part")
          ENDIF
          !-----------------------------------------------------------------------
          !  Energy term
          !-----------------------------------------------------------------------
          IF (ALLOCATED(energya)) THEN
             DEALLOCATE(energya,STAT=info)
             or_fail_dealloc("energya")
          ENDIF
          IF (ALLOCATED(labela)) THEN
             DEALLOCATE(labela,STAT=info)
             or_fail_dealloc("labela")
          ENDIF
          IF (ALLOCATED(candlabela)) THEN
             DEALLOCATE(candlabela,STAT=info)
             or_fail_dealloc("candlabela")
          ENDIF
          IF (ALLOCATED(ccandlabela)) THEN
             DEALLOCATE(ccandlabela,STAT=info)
             or_fail_dealloc("ccandlabela")
          ENDIF
          IF (ALLOCATED(ndaughtersa)) THEN
             DEALLOCATE(ndaughtersa,STAT=info)
             or_fail_dealloc("ndaughtersa")
          ENDIF
          IF (ALLOCATED(nmothersa)) THEN
             DEALLOCATE(nmothersa,STAT=info)
             or_fail_dealloc("nmothersa")
          ENDIF
          IF (ALLOCATED(mothersa)) THEN
             DEALLOCATE(mothersa,STAT=info)
             or_fail_dealloc("mothersa")
          ENDIF
          IF (ALLOCATED(daughtersa)) THEN
             DEALLOCATE(daughtersa,STAT=info)
             or_fail_dealloc("daughtersa")
          ENDIF
          IF (ALLOCATED(accepteda)) THEN
             DEALLOCATE(accepteda,STAT=info)
             or_fail_dealloc("accepteda")
          ENDIF
          IF (ALLOCATED(processeda)) THEN
             DEALLOCATE(processeda,STAT=info)
             or_fail_dealloc("processeda")
          ENDIF

          !----------------------------------------------------------------------
          !  Clear the old candidate set
          !----------------------------------------------------------------------
          sbpitr => mesh%subpatch%begin()
          ipatch=1
          DO WHILE (ASSOCIATED(sbpitr))
             IF (ALLOCATED(ppm_rc_seeds)) THEN
                CALL ppm_rc_seeds(ipatch)%destroy(info)
                or_fail("ppm_rc_seeds(ipatch)%destroy")
             ENDIF
             IF (ALLOCATED(ppm_rc_seeds_to_remove)) THEN
                CALL ppm_rc_seeds_to_remove(ipatch)%destroy(info)
                or_fail("ppm_rc_seeds_to_remove(ipatch)%destroy")
             ENDIF
             IF (ALLOCATED(InnerContourContainer)) THEN
                CALL InnerContourContainer(ipatch)%destroy(info)
                or_fail("InnerContourContainer(ipatch)%destroy")
             ENDIF
             IF (ALLOCATED(Candidates)) THEN
                CALL Candidates(ipatch)%destroy(info)
                or_fail("Candidates(ipatch)%destroy")
             ENDIF
             IF (ALLOCATED(CompetingRegions)) THEN
                CALL CompetingRegions(ipatch)%destroy(info)
                or_fail("CompetingRegions(ipatch)%destroy")
             ENDIF
             IF (ALLOCATED(m_Seeds)) THEN
                CALL m_Seeds(ipatch)%destroy(info)
                or_fail("m_Seeds(ipatch)%destroy")
             ENDIF
             sbpitr => mesh%subpatch%next()
             ipatch=ipatch+1
          ENDDO

          IF (ALLOCATED(ppm_rc_seeds)) THEN
             DEALLOCATE(ppm_rc_seeds,STAT=info)
             or_fail_dealloc("ppm_rc_seeds")
          ENDIF
          IF (ALLOCATED(ppm_rc_seeds_to_remove)) THEN
             DEALLOCATE(ppm_rc_seeds_to_remove,STAT=info)
             or_fail_dealloc("ppm_rc_seeds_to_remove")
          ENDIF
          IF (ALLOCATED(InnerContourContainer)) THEN
             DEALLOCATE(InnerContourContainer,STAT=info)
             or_fail_dealloc("InnerContourContainer")
          ENDIF
          IF (ALLOCATED(Candidates)) THEN
             DEALLOCATE(Candidates,STAT=info)
             or_fail_dealloc("Candidates")
          ENDIF
          IF (ALLOCATED(CompetingRegions)) THEN
             DEALLOCATE(CompetingRegions,STAT=info)
             or_fail_dealloc("CompetingRegions")
          ENDIF
          IF (ALLOCATED(m_Seeds)) THEN
             DEALLOCATE(m_Seeds,STAT=info)
             or_fail_dealloc("m_Seeds")
          ENDIF
          IF (ALLOCATED(ineighproc)) THEN
             DEALLOCATE(ineighproc,STAT=info)
             or_fail_dealloc("ineighproc")
          ENDIF
          IF (ALLOCATED(procflag)) THEN
             DEALLOCATE(procflag,STAT=info)
             or_fail_dealloc("procflag")
          ENDIF

          CALL e_data%shrink(info)
          or_fail("Failed to shrink the Extra memory array")

        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        END SUBROUTINE MCMCFreeMemory

        SUBROUTINE CreateMCMClengthProposalMask(info)
          !!! Create lengthProposalMask for fast proposal computation in biased
          !!! MCMC
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_error
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_error,ppm_err_alloc,ppm_err_dealloc, &
          &   ppm_err_argument

          USE ppm_rc_module_global, ONLy : one,ppm_rc_dim,MCMCuseBiasedProposal
          USE ppm_rc_module_topologicalnumber, ONLY : FG_ConnectivityType,BG_ConnectivityType
          IMPLICIT NONE

          INTEGER, INTENT(INOUT) :: info

          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          REAL(ppm_kind_double) :: t0

          INTEGER :: vLength,NumberOfNeighbors
          INTEGER :: i

          CHARACTER(LEN=ppm_char) :: caller='CreateMCMClengthProposalMask'
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          IF (.NOT.MCMCuseBiasedProposal) GOTO 9999

          CALL check

          IF (ALLOCATED(MCMClengthProposalMask)) THEN
             DEALLOCATE(MCMClengthProposalMask,STAT=info)
             or_fail_dealloc("MCMClengthProposalMask")
          ENDIF

          NumberOfNeighbors=BG_ConnectivityType%NumberOfNeighbors

          ALLOCATE(MCMClengthProposalMask(NumberOfNeighbors),STAT=info)
          or_fail_alloc("Failed to allocate MCMClengthProposalMask")

          !-------------------------------------------------------------------------
          ! Incorporate a smooth contour bias into a discrete proposal
          ! For smooth boundary proposal; weight is related to the contour length
          ! For weights we use the length approximation for contours and surfaces
          ! by Boykov and Kolmogorov
          ! 0.7=sqrt(2)/2
          !
          !     0.7  1  0.7
          !        \ | /
          !     1  - . - 1
          !        / | \
          !     0.7  1  0.7
          !
          ! Prepare a fast proposal computation:
          !-------------------------------------------------------------------------
          SELECT CASE (ppm_rc_dim)
          CASE (2)
             DO i=1,NumberOfNeighbors
                vLength=BG_ConnectivityType%NeighborsPoints(1,i)* &
                &       BG_ConnectivityType%NeighborsPoints(1,i)+ &
                &       BG_ConnectivityType%NeighborsPoints(2,i)* &
                &       BG_ConnectivityType%NeighborsPoints(2,i)

                MCMClengthProposalMask(i)=one/SQRT(REAL(vLength,MK))
             ENDDO
          CASE (3)
             DO i=1,NumberOfNeighbors
                vLength=BG_ConnectivityType%NeighborsPoints(1,i)* &
                &       BG_ConnectivityType%NeighborsPoints(1,i)+ &
                &       BG_ConnectivityType%NeighborsPoints(2,i)* &
                &       BG_ConnectivityType%NeighborsPoints(2,i)+ &
                &       BG_ConnectivityType%NeighborsPoints(3,i)* &
                &       BG_ConnectivityType%NeighborsPoints(3,i)

                MCMClengthProposalMask(i)=one/SQRT(REAL(vLength,MK))
             ENDDO
          END SELECT
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        CONTAINS
        SUBROUTINE check
        IMPLICIT NONE
          check_true(<#ASSOCIATED(FG_ConnectivityType)#>, &
          & "Topology Connectivity is not defined yet!",exit_point=8888)
        8888 CONTINUE
        END SUBROUTINE check
        END SUBROUTINE CreateMCMClengthProposalMask

        SUBROUTINE DestroyMCMClengthProposalMask(info)
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_error
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_error,ppm_err_dealloc
          IMPLICIT NONE

          INTEGER, INTENT(INOUT) :: info

          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          REAL(ppm_kind_double) :: t0

          CHARACTER(LEN=ppm_char) :: caller='DestroyMCMClengthProposalMask'
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)
          IF (ALLOCATED(MCMClengthProposalMask)) THEN
             DEALLOCATE(MCMClengthProposalMask,STAT=info)
             or_fail_dealloc("MCMClengthProposalMask")
          ENDIF
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        END SUBROUTINE DestroyMCMClengthProposalMask

        FUNCTION MCMCInsertBGLabelInRegionLabel() RESULT(info)
          !!! Craete MCMCRegionLabel array and insert the Background label in it
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_fatal
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_error,ppm_err_alloc,ppm_err_argument
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
          INTEGER :: info
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          REAL(ppm_kind_double) :: t0

          CHARACTER(LEN=ppm_char) :: caller='MCMCInsertLabelInRegionLabel'
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          IF (MCMCRegionLabelSize.EQ.0) THEN
             MCMCRegionLabelSize=1
             ALLOCATE(MCMCRegionLabel(MCMCRegionLabelSize),STAT=info)
             or_fail_alloc("MCMCRegionLabel",ppm_error=ppm_error_fatal)
             MCMCRegionLabel(1)=0
          ELSE
             fail("MCMCRegionLabel is not allocated yet and can not have MCMCRegionLabelSize bigger than 0!", &
             & ppm_error=ppm_error_fatal)
          ENDIF
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        END FUNCTION MCMCInsertBGLabelInRegionLabel

        FUNCTION MCMCInsertLabelInRegionLabel(LabelIn) RESULT(info)
          !!! Insert the new label in the MCMCRegionLabel array if
          !!! it does not already exist there. (local on the processor)
          !!! The new label will be in the right place in an ascending sorted
          !!! array order
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_fatal
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_error,ppm_err_alloc,ppm_err_argument

          USE ppm_rc_module_util, ONLY : ppm_rc_label_exist,ppm_rc_label_index
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER, INTENT(IN   ) :: LabelIn
          !!! A positive label (we put in regionsLabel, if it does not already exist)
          INTEGER                :: info
          !!! 0 status on success
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          REAL(ppm_kind_double) :: t0

          INTEGER :: labelindex
          INTEGER :: nsize
          INTEGER :: i

          CHARACTER(LEN=ppm_char) :: caller='MCMCInsertLabelInRegionLabel'
          CHARACTER(LEN=ppm_char) :: cbuf
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          IF (.NOT.ppm_rc_label_exist(LabelIn,MCMCRegionLabel,MCMCRegionLabelSize)) THEN
             ! Find the index, where the new label will be inside the array of labels
             labelindex=ppm_rc_label_index(LabelIn,MCMCRegionLabel,MCMCRegionLabelSize)

             ! Here nsize is at least 1 which is the background label
             nsize=SIZE(MCMCRegionLabel,DIM=1)

             MCMCRegionLabelSize=MCMCRegionLabelSize+1

             IF (MCMCRegionLabelSize.GT.nsize) THEN
                ALLOCATE(tmp1_i(nsize*2),STAT=info)
                or_fail_alloc("tmp1_i",ppm_error=ppm_error_fatal)

                FORALL (i=1:labelindex-1) tmp1_i(i)=MCMCRegionLabel(i)
                tmp1_i(labelindex)=LabelIn
                FORALL (i=labelindex:nsize) tmp1_i(i+1)=MCMCRegionLabel(i)

                CALL MOVE_ALLOC(tmp1_i,MCMCRegionLabel)
             ELSE
                IF      (labelindex.EQ.MCMCRegionLabelSize) THEN
                   MCMCRegionLabel(MCMCRegionLabelSize)=LabelIn
                ELSE IF (labelindex.GT.1) THEN
                   DO i=MCMCRegionLabelSize,labelindex+1,-1
                      MCMCRegionLabel(i)=MCMCRegionLabel(i-1)
                   ENDDO
                   MCMCRegionLabel(labelindex)=LabelIn
                ELSE IF (labelindex.EQ.1) THEN
                   DO i=MCMCRegionLabelSize,2,-1
                      MCMCRegionLabel(i)=MCMCRegionLabel(i-1)
                   ENDDO
                   MCMCRegionLabel(1)=LabelIn
                ! Here labelindex.EQ.0 only means that LabelIn is negative
                ! which is not acceptable
                ELSE IF (labelindex.EQ.0) THEN
                   WRITE(cbuf,'(A,I10,A)')"Inserting negative label ",LabelIn," is not possible!"
                   fail(cbuf,ppm_error=ppm_error_fatal)
                ENDIF
             ENDIF
          ENDIF
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        END FUNCTION MCMCInsertLabelInRegionLabel

        FUNCTION MCMCUpdateRegionLabel() RESULT(info)
          !!! Set up a vector MCMCRegionLabel that maps natural numbers, index of
          !!! the vector to the region labels (which are local on this processor).
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_fatal
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_error,ppm_err_dealloc,ppm_err_sub_failed
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER :: info
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          REAL(ppm_kind_double) :: t0

          INTEGER :: nsize
          INTEGER :: i,RLabel

          CHARACTER(LEN=ppm_char) :: caller='MCMCUpdateRegionLabel'
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          DEALLOCATE(MCMCRegionLabel,STAT=info)
          or_fail_dealloc("MCMCRegionLabel",ppm_error=ppm_error_fatal)

          MCMCRegionLabelSize=0

          ! The first element is the background label
          info=MCMCInsertLabelInRegionLabel()
          or_fail("MCMCInsertLabelInRegionLabel",ppm_error=ppm_error_fatal)

          ! Now we put the other labels
          nsize=SIZE(MCMCRegularParticles)

          DO i=1,nsize
             ! If the region exists and contains Particles
             IF (MCMCRegularParticles(i)%size().GT.0) THEN

                ! Find the region Label
                RLabel=e_data%Rlabel(i)

                ! Insert the label inside the array
                info=MCMCInsertLabelInRegionLabel(RLabel)
                or_fail("MCMCInsertLabelInRegionLabel",ppm_error=ppm_error_fatal)
             ENDIF
          ENDDO
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        END FUNCTION MCMCUpdateRegionLabel

        FUNCTION MCMCgetRegularParticlesAtIndex_2d(labels_,coord,Nm,MCMCParticles) RESULT(nsize)
          !!! This function returns nsize regular particles at label index
          !!! nsize particles are :
          !!! children of neighboring regions plus
          !!! if label at label index is non zero (one parent)
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          USE ppm_rc_module_global, ONLY : FORBIDDEN
          USE ppm_rc_module_topologicalnumber, ONLY : FG_ConnectivityType, &
          &   BG_ConnectivityType
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER, CONTIGUOUS, DIMENSION(:,:), POINTER       :: labels_
          !!! image labels
          INTEGER,             DIMENSION(:),   INTENT(IN   ) :: coord
          !!! Coordinates of the particle position
          INTEGER,             DIMENSION(:),   INTENT(IN   ) :: Nm
          !!! Domain size
          TYPE(MCMCParticle),  DIMENSION(:),   INTENT(INOUT) :: MCMCParticles
          !!! Regular particles at coordinate and its neighbors
          INTEGER                                            :: nsize
          !!! Number of Regular particle (size of MCMCParticles)
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          REAL(MK) :: Proposal

          INTEGER, DIMENSION(:,:), POINTER :: tmplabels
          INTEGER, DIMENSION(8)            :: LabelNghtmp
          INTEGER, DIMENSION(2)            :: ll
          INTEGER                          :: cLabel
          INTEGER                          :: LabelNgh
          INTEGER                          :: i

          LOGICAL :: NotParentInserted

          nsize=0

          ! We are outside real domain, we should not put particle here
          IF (ANY(coord.LT.1.OR.coord.GT.Nm)) RETURN

          IF (labels_(coord(1),coord(2)).EQ.FORBIDDEN) RETURN

          tmplabels => labels_(coord(1)-1:coord(1)+1,coord(2)-1:coord(2)+1)

          cLabel=ABS(tmplabels(2,2))

          Proposal=MCMCproposal(tmplabels)

          SELECT CASE (cLabel)
          ! It can not be a parent
          CASE (0)
             DO i=1,FG_ConnectivityType%NumberOfNeighbors
                ll=FG_ConnectivityType%NeighborsPoints(:,i)+2
                LabelNgh=ABS(tmplabels(ll(1),ll(2)))
                ! here there should be a particle since two neighboring pixels
                ! have different labels.
                IF (LabelNgh.NE.0.AND.LabelNgh.NE.FORBIDDEN) THEN
                   IF (ANY(LabelNgh.EQ.LabelNghtmp(1:nsize))) CYCLE
                   ! FG labels have a daughter placed at this spot.
                   nsize=nsize+1
                   MCMCParticles(nsize)=MCMCParticle(LabelNgh,Proposal)
                   LabelNghtmp(nsize)=LabelNgh
                ENDIF
             ENDDO
          CASE DEFAULT
             NotParentInserted=.TRUE.

             DO i=1,FG_ConnectivityType%NumberOfNeighbors
                ll=FG_ConnectivityType%NeighborsPoints(:,i)+2
                LabelNgh=ABS(tmplabels(ll(1),ll(2)))
                ! here there should be a particle since two neighboring pixels
                ! have different labels.
                IF (LabelNgh.NE.cLabel.AND.LabelNgh.NE.FORBIDDEN) THEN
                   ! FG labels have a daughter placed at this spot.
                   IF (LabelNgh.NE.0) THEN
                      IF (ANY(LabelNgh.EQ.LabelNghtmp(1:nsize))) CYCLE
                      nsize=nsize+1
                      MCMCParticles(nsize)=MCMCParticle(LabelNgh,Proposal)
                      LabelNghtmp(nsize)=LabelNgh
                   ENDIF
                   ! this is a non-zero pixel with different neighbors,
                   ! hence there must be a mother in the list:
                   IF (NotParentInserted) THEN
                      NotParentInserted=.FALSE.
                      nsize=nsize+1
                      MCMCParticles(nsize)=MCMCParticle(0,Proposal)
                      LabelNghtmp(nsize)=0
                   ENDIF
                ENDIF
             ENDDO !i=1,FG_ConnectivityType%NumberOfNeighbors

             ! Check the BG neighborhood now if we need to insert a parent.
             IF (NotParentInserted) THEN
                DO i=1,BG_ConnectivityType%NumberOfNeighbors
                   ll=BG_ConnectivityType%NeighborsPoints(:,i)+2
                   LabelNgh=ABS(tmplabels(ll(1),ll(2)))
                   ! This is a FG pixel with a neighbor of a different label.
                   ! finally, insert a parent particle
                   IF (cLabel.NE.LabelNgh) THEN
                      nsize=nsize+1
                      MCMCParticles(nsize)=MCMCParticle(0,Proposal)
                      EXIT
                   ENDIF
                ENDDO !i=1,BG_ConnectivityType%NumberOfNeighbors
             ENDIF !NotParentInserted
          END SELECT !cLabel
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCgetRegularParticlesAtIndex_2d

        FUNCTION MCMCgetRegularParticlesAtIndex__2d(tmplabels,coord1,coord2,Nm,MCMCParticles) RESULT(nsize)
          !!! This function returns nsize regular particles at label index
          !!! nsize particles are children of neighboring regions and
          !!! if label at label index is non zero (one parent)
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          USE ppm_rc_module_global, ONLY : FORBIDDEN
          USE ppm_rc_module_topologicalnumber, ONLY : FG_ConnectivityType, &
          &   BG_ConnectivityType
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER,            DIMENSION(:,:), POINTER       :: tmplabels
          !!! image labels
          INTEGER,                            INTENT(IN   ) :: coord1
          !!! X Coordinates of the particle position
          INTEGER,                            INTENT(IN   ) :: coord2
          !!! Y Coordinates of the particle position
          INTEGER,            DIMENSION(:),   INTENT(IN   ) :: Nm
          !!! Domain Size
          TYPE(MCMCParticle), DIMENSION(:),   INTENT(INOUT) :: MCMCParticles
          !!! Regular particles at coordinate and its neighbors
          INTEGER                                           :: nsize
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          REAL(MK) :: Proposal

          INTEGER, DIMENSION(2) :: ll
          INTEGER, DIMENSION(8) :: LabelNghtmp
          INTEGER               :: cLabel
          INTEGER               :: LabelNgh
          INTEGER               :: i

          LOGICAL :: NotParentInserted

          nsize=0

          ! We are outside real domain, we should not put particle here
          IF (coord1.LT.1.OR.coord2.LT.1.OR.coord1.GT.Nm(1).OR.coord2.GT.Nm(2)) RETURN

          IF (tmplabels(2,2).EQ.FORBIDDEN) RETURN

          cLabel=ABS(tmplabels(2,2))

          proposal=MCMCproposal(tmplabels)

          SELECT CASE (cLabel)
          ! we can put a child particle here
          CASE (0)
             DO i=1,FG_ConnectivityType%NumberOfNeighbors
                ll=FG_ConnectivityType%NeighborsPoints(:,i)+2
                LabelNgh=ABS(tmplabels(ll(1),ll(2)))
                IF (LabelNgh.NE.0.AND.LabelNgh.NE.FORBIDDEN) THEN
                   IF (ANY(LabelNgh.EQ.LabelNghtmp(1:nsize))) CYCLE
                   ! here there should be a particle since two neighboring pixels
                   ! have different labels.
                   ! FG labels have a daughter placed at this spot.
                   nsize=nsize+1
                   MCMCParticles(nsize)=MCMCParticle(LabelNgh,proposal)
                   LabelNghtmp(nsize)=LabelNgh
                ENDIF
             ENDDO
          ! Either a parent or a child
          CASE DEFAULT
             NotParentInserted=.TRUE.
             DO i=1,FG_ConnectivityType%NumberOfNeighbors
                ll=FG_ConnectivityType%NeighborsPoints(:,i)+2
                LabelNgh=ABS(tmplabels(ll(1),ll(2)))
                IF (LabelNgh.NE.cLabel.AND.LabelNgh.NE.FORBIDDEN) THEN
                   ! here there should be a particle since two neighboring pixels
                   ! have different labels.
                   IF (LabelNgh.NE.0) THEN
                      IF (ANY(LabelNgh.EQ.LabelNghtmp(1:nsize))) CYCLE
                      ! FG labels have a daughter placed at this spot.
                      nsize=nsize+1
                      MCMCParticles(nsize)=MCMCParticle(LabelNgh,proposal)
                      LabelNghtmp(nsize)=LabelNgh
                   ENDIF
                   IF (NotParentInserted) THEN
                      ! this is a non-zero pixel with different neighbors,
                      ! hence there must be a mother in the list:
                      NotParentInserted=.FALSE.
                      nsize=nsize+1
                      MCMCParticles(nsize)=MCMCParticle(0,proposal)
                      LabelNghtmp(nsize)=0
                   ENDIF
                ENDIF
             ENDDO !i=1,FG_ConnectivityType%NumberOfNeighbors

             ! Check the BG neighborhood now if we need to insert a parent.
             IF (NotParentInserted) THEN
                DO i=1,BG_ConnectivityType%NumberOfNeighbors
                   ll=BG_ConnectivityType%NeighborsPoints(:,i)+2
                   LabelNgh=ABS(tmplabels(ll(1),ll(2)))
                   IF (cLabel.NE.LabelNgh) THEN
                      ! This is a FG pixel with a neighbor of a different label.
                      ! finally, insert a parent particle
                      nsize=nsize+1
                      MCMCParticles(nsize)=MCMCParticle(0,proposal)
                      EXIT
                   ENDIF
                ENDDO !i=1,BG_ConnectivityType%NumberOfNeighbors
             ENDIF !NotParentInserted
          END SELECT !cLabel
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCgetRegularParticlesAtIndex__2d

        FUNCTION MCMCgetRegularParticlesAtIndex_3d(labels_,coord,Nm,MCMCParticles) RESULT(nsize)
          !!! This function returns nsize regular particles at label index
          !!! nsize particles are children of neighboring regions and
          !!! if label at label index is non zero (one parent)
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          USE ppm_rc_module_global, ONLY : FORBIDDEN
          USE ppm_rc_module_topologicalnumber, ONLY : FG_ConnectivityType, &
          &   BG_ConnectivityType
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: labels_
          !!! image labels
          INTEGER,             DIMENSION(:),     INTENT(IN   ) :: coord
          !!! Coordinates of the particle position
          INTEGER,             DIMENSION(:),     INTENT(IN   ) :: Nm
          !!! Domain Size
          TYPE(MCMCParticle),  DIMENSION(:),     INTENT(INOUT) :: MCMCParticles
          !!! Regular particles at coordinate and its neighbors
          INTEGER                                              :: nsize
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          REAL(MK) :: proposal

          INTEGER, DIMENSION(:,:,:), POINTER :: tmplabels
          INTEGER, DIMENSION(27)             :: LabelNghtmp
          INTEGER, DIMENSION(3)              :: ll
          INTEGER                            :: cLabel
          INTEGER                            :: LabelNgh
          INTEGER                            :: i

          LOGICAL :: NotParentInserted

          nsize=0

          ! We are outside real domain, we should not put a particle here
          IF (ANY(coord.LT.1.OR.coord.GT.Nm)) RETURN

          IF (labels_(coord(1),coord(2),coord(3)).EQ.FORBIDDEN) RETURN

          tmplabels => labels_(coord(1)-1:coord(1)+1,coord(2)-1:coord(2)+1,coord(3)-1:coord(3)+1)

          cLabel=ABS(tmplabels(2,2,2))

          proposal=MCMCproposal(tmplabels)

          SELECT CASE (cLabel)
          CASE (0)
             DO i=1,FG_ConnectivityType%NumberOfNeighbors
                ll=FG_ConnectivityType%NeighborsPoints(:,i)+2
                LabelNgh=ABS(tmplabels(ll(1),ll(2),ll(3)))
                IF (LabelNgh.NE.0.AND.LabelNgh.NE.FORBIDDEN) THEN
                   IF (ANY(LabelNgh.EQ.LabelNghtmp(1:nsize))) CYCLE
                   ! here there should be a particle since two neighboring pixels
                   ! have different labels.
                   ! FG labels have a daughter placed at this spot.
                   nsize=nsize+1
                   MCMCParticles(nsize)=MCMCParticle(LabelNgh,proposal)
                   LabelNghtmp(nsize)=LabelNgh
                ENDIF
             ENDDO
          CASE DEFAULT
             NotParentInserted=.TRUE.
             ! If coord is outside the domain, we should not put a parent
             ! there
             DO i=1,FG_ConnectivityType%NumberOfNeighbors
                ll=FG_ConnectivityType%NeighborsPoints(:,i)+2
                LabelNgh=ABS(tmplabels(ll(1),ll(2),ll(3)))
                IF (LabelNgh.NE.cLabel.AND.LabelNgh.NE.FORBIDDEN) THEN
                ! here there should be a particle since two neighboring pixels
                ! have different labels.
                   IF (LabelNgh.NE.0) THEN
                      IF (ANY(LabelNgh.EQ.LabelNghtmp(1:nsize))) CYCLE
                      ! FG labels have a daughter placed at this spot.
                      nsize=nsize+1
                      MCMCParticles(nsize)=MCMCParticle(LabelNgh,proposal)
                      LabelNghtmp(nsize)=LabelNgh
                   ENDIF
                   IF (NotParentInserted) THEN
                   ! this is a non-zero pixel with different neighbors,
                   ! hence there must be a mother in the list:
                      NotParentInserted=.FALSE.
                      nsize=nsize+1
                      MCMCParticles(nsize)=MCMCParticle(0,proposal)
                      LabelNghtmp(nsize)=0
                   ENDIF
                ENDIF
             ENDDO !i=1,FG_ConnectivityType%NumberOfNeighbors

             ! Check the BG neighborhood now if we need to insert a parent.
             IF (NotParentInserted) THEN
                DO i=1,BG_ConnectivityType%NumberOfNeighbors
                   ll=BG_ConnectivityType%NeighborsPoints(:,i)+2
                   LabelNgh=ABS(tmplabels(ll(1),ll(2),ll(3)))
                   IF (cLabel.NE.LabelNgh) THEN
                      ! This is a FG pixel with a neighbor of a different label.
                      ! finally, insert a parent particle
                      nsize=nsize+1
                      MCMCParticles(nsize)=MCMCParticle(0,proposal)
                      EXIT
                   ENDIF
                ENDDO !i=1,BG_ConnectivityType%NumberOfNeighbors
             ENDIF !NotParentInserted
          END SELECT !cLabel
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCgetRegularParticlesAtIndex_3d

        FUNCTION MCMCgetRegularParticlesAtIndex__3d(tmplabels,coord1,coord2,coord3,Nm,MCMCParticles) RESULT(nsize)
          !!! This function returns nsize regular particles at label index
          !!! nsize particles are children of neighboring regions and
          !!! if label at label index is non zero (one parent)
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          USE ppm_rc_module_global, ONLY : FORBIDDEN
          USE ppm_rc_module_topologicalnumber, ONLY : FG_ConnectivityType, &
          &   BG_ConnectivityType
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER,            DIMENSION(:,:,:), POINTER       :: tmplabels
          !!! image labels
          INTEGER,                              INTENT(IN   ) :: coord1
          !!! X Coordinates of the particle position
          INTEGER,                              INTENT(IN   ) :: coord2
          !!! Y Coordinates of the particle position
          INTEGER,                              INTENT(IN   ) :: coord3
          !!! Z Coordinates of the particle position
          INTEGER,            DIMENSION(:),     INTENT(IN   ) :: Nm
          !!! Domain size
          TYPE(MCMCParticle), DIMENSION(:),     INTENT(INOUT) :: MCMCParticles
          !!! Regular particles at coordinate and its neighbors
          INTEGER                                             :: nsize
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          REAL(MK) :: proposal

          INTEGER, DIMENSION(3)  :: ll
          INTEGER, DIMENSION(27) :: LabelNghtmp
          INTEGER                :: cLabel
          INTEGER                :: LabelNgh
          INTEGER                :: i

          LOGICAL :: NotParentInserted

          nsize=0

          ! We are outside real domain, we should not put a particle here
          IF (coord1.LT.1    .OR.coord2.LT.1    .OR.coord3.LT.1.OR. &
          &   coord1.GT.Nm(1).OR.coord2.GT.Nm(2).OR.coord3.GT.Nm(3)) RETURN

          IF (tmplabels(2,2,2).EQ.FORBIDDEN) RETURN

          cLabel=ABS(tmplabels(2,2,2))

          proposal=MCMCproposal(tmplabels)

          SELECT CASE (cLabel)
          CASE (0)
             DO i=1,FG_ConnectivityType%NumberOfNeighbors
                ll=FG_ConnectivityType%NeighborsPoints(:,i)+2
                LabelNgh=ABS(tmplabels(ll(1),ll(2),ll(3)))
                IF (LabelNgh.NE.0.AND.LabelNgh.NE.FORBIDDEN) THEN
                   IF (ANY(LabelNgh.EQ.LabelNghtmp(1:nsize))) CYCLE
                   ! here there should be a particle since two neighboring pixels
                   ! have different labels.
                   ! FG labels have a daughter placed at this spot.
                   nsize=nsize+1
                   MCMCParticles(nsize)=MCMCParticle(LabelNgh,proposal)
                   LabelNghtmp(nsize)=LabelNgh
                ENDIF
             ENDDO
          CASE DEFAULT
             NotParentInserted=.NOT.(coord1.LT.1.OR.coord2.LT.1.OR.coord3.LT.1.OR. &
             &                       coord1.GT.Nm(1).OR.coord2.GT.Nm(2).OR.coord3.GT.Nm(3))
             ! If coord is outside the domain, we should not put a parent
             ! there
             DO i=1,FG_ConnectivityType%NumberOfNeighbors
                ll=FG_ConnectivityType%NeighborsPoints(:,i)+2
                LabelNgh=ABS(tmplabels(ll(1),ll(2),ll(3)))
                IF (LabelNgh.NE.cLabel.AND.LabelNgh.NE.FORBIDDEN) THEN
                ! here there should be a particle since two neighboring pixels
                ! have different labels.
                   IF (LabelNgh.NE.0) THEN
                      IF (ANY(LabelNgh.EQ.LabelNghtmp(1:nsize))) CYCLE
                      ! FG labels have a daughter placed at this spot.
                      nsize=nsize+1
                      MCMCParticles(nsize)=MCMCParticle(LabelNgh,proposal)
                      LabelNghtmp(nsize)=LabelNgh
                   ENDIF
                   IF (NotParentInserted) THEN
                   ! this is a non-zero pixel with different neighbors,
                   ! hence there must be a mother in the list:
                      NotParentInserted=.FALSE.
                      nsize=nsize+1
                      MCMCParticles(nsize)=MCMCParticle(0,proposal)
                      LabelNghtmp(nsize)=0
                   ENDIF
                ENDIF
             ENDDO !i=1,FG_ConnectivityType%NumberOfNeighbors

             ! Check the BG neighborhood now if we need to insert a parent.
             IF (NotParentInserted) THEN
                DO i=1,BG_ConnectivityType%NumberOfNeighbors
                   ll=BG_ConnectivityType%NeighborsPoints(:,i)+2
                   LabelNgh=ABS(tmplabels(ll(1),ll(2),ll(3)))
                   IF (cLabel.NE.LabelNgh) THEN
                      ! This is a FG pixel with a neighbor of a different label.
                      ! finally, insert a parent particle
                      nsize=nsize+1
                      MCMCParticles(nsize)=MCMCParticle(0,proposal)
                      EXIT
                   ENDIF
                ENDDO !i=1,BG_ConnectivityType%NumberOfNeighbors
             ENDIF !NotParentInserted
          END SELECT !cLabel
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCgetRegularParticlesAtIndex__3d

        FUNCTION MCMCproposal_2d(labels_,coord)
          !!! Compute the Length-prior driven proposal at coord
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          USE ppm_rc_module_global, ONLy : zero,half,one,MCMCuseBiasedProposal
          USE ppm_rc_module_topologicalnumber, ONLY : BG_ConnectivityType
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER, CONTIGUOUS, DIMENSION(:,:), POINTER       :: labels_
          INTEGER,             DIMENSION(:),   INTENT(IN   ) :: coord

          REAL(MK)                                           :: MCMCproposal_2d
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          INTEGER, DIMENSION(:,:), POINTER :: tmplabels
          INTEGER, DIMENSION(2)            :: ll
          INTEGER                          :: cLabel
          INTEGER                          :: LabelNgh
          INTEGER                          :: i

          LOGICAL :: vFloating
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          IF (.NOT.MCMCuseBiasedProposal) THEN
             MCMCproposal_2d=one
             RETURN
          ENDIF
          ! The choice of a data-driven proposal distribution depends a lot on
          ! the contrast in the image (?) and the current segmentation.

          ! L1: 1 / ( 1 + abs(d) )
          ! MCMCproposal_2d=1./(1.+abs(Means[Particle.m_CandidateLabel]-Image(Particle.Index)))
          ! L2: 1 / d^2
          ! MCMCproposal_2d=1./(1.+abs(Means[Particle.m_CandidateLabel]-Image(Particle.Index)))
          ! exp(-d^2)
          ! MCMCproposal_2d=exp((Means[Particle.CandidateLabel] - Image(Particle.Index)) * &
          ! &                   (Means[Particle.CandidateLabel] - Image(Particle.Index)))

          ! Length-prior driven proposal:
          MCMCproposal_2d=zero
          vFloating=.TRUE.

          tmplabels => labels_(coord(1)-1:coord(1)+1,coord(2)-1:coord(2)+1)
          cLabel=ABS(tmplabels(2,2))

          DO i=1,BG_ConnectivityType%NumberOfNeighbors
             ll=BG_ConnectivityType%NeighborsPoints(:,i)+2
             LabelNgh=ABS(tmplabels(ll(1),ll(2)))
             IF (cLabel.NE.LabelNgh) THEN
                MCMCproposal_2d=MCMCproposal_2d+MCMClengthProposalMask(i)
                vFloating=.FALSE.
             ENDIF
          ENDDO

          IF (vFloating) THEN
             ! floating particles need a proposal > 0: We take half of the smallest
             ! element in the mask (hack:we assume the smallest element to be at
             ! position 1).
             MCMCproposal_2d = MCMClengthProposalMask(1)*half
          ENDIF
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCproposal_2d

        FUNCTION MCMCproposal__2d(tmplabels)
          !!! Compute the Length-prior driven proposal at coord
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          USE ppm_rc_module_global, ONLy : zero,half,one,MCMCuseBiasedProposal
          USE ppm_rc_module_topologicalnumber, ONLY : BG_ConnectivityType
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER, DIMENSION(:,:), POINTER :: tmplabels
          !!! Temporary label index around coordinate

          REAL(MK)                         :: MCMCproposal__2d
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          INTEGER, DIMENSION(2) :: ll
          INTEGER               :: cLabel
          INTEGER               :: LabelNgh
          INTEGER               :: i

          LOGICAL :: vFloating
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          IF (.NOT.MCMCuseBiasedProposal) THEN
             MCMCproposal__2d=one
             RETURN
          ENDIF
          ! The choice of a data-driven proposal distribution depends a lot on
          ! the contrast in the image (?) and the current segmentation.

          ! L1: 1 / ( 1 + abs(d) )
          ! MCMCproposal__2d=1./(1.+abs(Means[Particle.m_CandidateLabel]-Image(Particle.Index)))
          ! L2: 1 / d^2
          ! MCMCproposal__2d=1./(1.+abs(Means[Particle.m_CandidateLabel]-Image(Particle.Index)))
          ! exp(-d^2)
          ! MCMCproposal__2d=exp((Means[Particle.CandidateLabel] - Image(Particle.Index)) * &
          ! &                   (Means[Particle.CandidateLabel] - Image(Particle.Index)))

          ! Length-prior driven proposal:
          MCMCproposal__2d=zero
          vFloating=.TRUE.

          cLabel=ABS(tmplabels(2,2))

          DO i=1,BG_ConnectivityType%NumberOfNeighbors
             ll=BG_ConnectivityType%NeighborsPoints(:,i)+2
             LabelNgh=ABS(tmplabels(ll(1),ll(2)))
             IF (cLabel.NE.LabelNgh) THEN
                MCMCproposal__2d=MCMCproposal__2d+MCMClengthProposalMask(i)
                vFloating=.FALSE.
             ENDIF
          ENDDO

          IF (vFloating) THEN
             ! floating particles need a proposal > 0: We take half of the smallest
             ! element in the mask (hack:we assume the smallest element to be at
             ! position 0).
             MCMCproposal__2d = MCMClengthProposalMask(1)*half
          ENDIF
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCproposal__2d

        FUNCTION MCMCproposal_3d(labels_,coord)
          !!! Compute the Length-prior driven proposal at coord
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          USE ppm_rc_module_global, ONLy : zero,half,one,MCMCuseBiasedProposal
          USE ppm_rc_module_topologicalnumber, ONLY : BG_ConnectivityType
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: labels_
          INTEGER,             DIMENSION(:),     INTENT(IN   ) :: coord

          REAL(MK)                                             :: MCMCproposal_3d
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          INTEGER, DIMENSION(:,:,:), POINTER :: tmplabels
          INTEGER, DIMENSION(3)              :: ll
          INTEGER                            :: cLabel
          INTEGER                            :: LabelNgh
          INTEGER                            :: i

          LOGICAL :: vFloating
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          IF (.NOT.MCMCuseBiasedProposal) THEN
             MCMCproposal_3d=one
             RETURN
          ENDIF
          ! The choice of a data-driven proposal distribution depends a lot on
          ! the contrast in the image (?) and the current segmentation.

          ! L1: 1 / ( 1 + abs(d) )
          ! MCMCproposal_3d=1./(1.+abs(Means[Particle.m_CandidateLabel]-Image(Particle.Index)))
          ! L2: 1 / d^2
          ! MCMCproposal_3d=1./(1.+abs(Means[Particle.m_CandidateLabel]-Image(Particle.Index)))
          ! exp(-d^2)
          ! MCMCproposal_3d=exp((Means[Particle.CandidateLabel] - Image(Particle.Index)) * &
          ! &                   (Means[Particle.CandidateLabel] - Image(Particle.Index)))

          ! Length-prior driven proposal:
          MCMCproposal_3d=zero
          vFloating=.TRUE.

          tmplabels => labels_(coord(1)-1:coord(1)+1,coord(2)-1:coord(2)+1,coord(3)-1:coord(3)+1)

          cLabel=ABS(tmplabels(2,2,2))

          DO i=1,BG_ConnectivityType%NumberOfNeighbors
             ll=BG_ConnectivityType%NeighborsPoints(:,i)+2
             LabelNgh=ABS(tmplabels(ll(1),ll(2),ll(3)))
             IF (cLabel.NE.LabelNgh) THEN
                MCMCproposal_3d=MCMCproposal_3d+MCMClengthProposalMask(i)
                vFloating=.FALSE.
             ENDIF
          ENDDO

          IF (vFloating) THEN
             ! floating particles need a proposal > 0: We take half of the smallest
             ! element in the mask (hack:we assume the smallest element to be at
             ! position 0).
             MCMCproposal_3d = MCMClengthProposalMask(1)*half
          ENDIF
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCproposal_3d

        FUNCTION MCMCproposal__3d(tmplabels)
          !!! Compute the Length-prior driven proposal at coord
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          USE ppm_rc_module_global, ONLy : zero,half,one,MCMCuseBiasedProposal
          USE ppm_rc_module_topologicalnumber, ONLY : BG_ConnectivityType
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER, DIMENSION(:,:,:), POINTER :: tmplabels
          !!! Temporary label index around coordinate

          REAL(MK)                           :: MCMCproposal__3d
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          INTEGER, DIMENSION(3) :: ll
          INTEGER               :: cLabel
          INTEGER               :: LabelNgh
          INTEGER               :: i

          LOGICAL :: vFloating
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          IF (.NOT.MCMCuseBiasedProposal) THEN
             MCMCproposal__3d=one
             RETURN
          ENDIF
          ! The choice of a data-driven proposal distribution depends a lot on
          ! the contrast in the image (?) and the current segmentation.

          ! L1: 1 / ( 1 + abs(d) )
          ! MCMCproposal__3d=1./(1.+abs(Means[Particle.m_CandidateLabel]-Image(Particle.Index)))
          ! L2: 1 / d^2
          ! MCMCproposal__3d=1./(1.+abs(Means[Particle.m_CandidateLabel]-Image(Particle.Index)))
          ! exp(-d^2)
          ! MCMCproposal__3d=exp((Means[Particle.CandidateLabel] - Image(Particle.Index)) * &
          ! &                   (Means[Particle.CandidateLabel] - Image(Particle.Index)))

          ! Length-prior driven proposal:
          MCMCproposal__3d=zero
          vFloating=.TRUE.

          cLabel=ABS(tmplabels(2,2,2))

          DO i=1,BG_ConnectivityType%NumberOfNeighbors
             ll=BG_ConnectivityType%NeighborsPoints(:,i)+2
             LabelNgh=ABS(tmplabels(ll(1),ll(2),ll(3)))
             IF (cLabel.NE.LabelNgh) THEN
                MCMCproposal__3d=MCMCproposal__3d+MCMClengthProposalMask(i)
                vFloating=.FALSE.
             ENDIF
          ENDDO

          IF (vFloating) THEN
             ! floating particles need a proposal > 0: We take half of the smallest
             ! element in the mask (hack:we assume the smallest element to be at
             ! position 0).
             MCMCproposal__3d = MCMClengthProposalMask(1)*half
          ENDIF
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCproposal__3d

        ! TOCHECK : The original implementation is wrong!!!!
        ! Returns false if topology is changed when applying this particle. This
        ! is based on the current state of the label image.
        ! TODO: this method should be changed to achieve full topological control.
        !      now the particle is rejected if it changes somehow the topology.
        LOGICAL FUNCTION MCMCIsParticleTopoValid_2d(labels_,coord,CandidateLabel)
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          USE ppm_rc_module_global, ONLy : AllowFission,AllowHandles
          USE ppm_rc_module_topologicalnumber, ONLY : FGandBGTopoNbPairType, &
          &   TopologicalNumberFunction
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER, CONTIGUOUS, DIMENSION(:,:), POINTER       :: labels_
          INTEGER,             DIMENSION(:),   INTENT(IN   ) :: coord
          INTEGER,                             INTENT(IN   ) :: CandidateLabel
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(FGandBGTopoNbPairType) :: FGBGTNP

          INTEGER :: ContainerLabel
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------

          IF (.NOT.AllowFission.OR..NOT.AllowHandles) THEN
             ! If the particle is a parent
             IF (CandidateLabel.EQ.0) THEN
                ! Get the correct label to access the container of the particle
                ContainerLabel=ABS(labels_(coord(1),coord(2)))

                FGBGTNP = TopologicalNumberFunction%EvaluateFGTNOfLabelAtIndex(coord,labels_,ContainerLabel)

                IF (FGBGTNP%FGNumber.NE.1.OR.FGBGTNP%BGNumber.NE.1) THEN
                   MCMCIsParticleTopoValid_2d=.FALSE.
                   RETURN
                ENDIF
             ELSE
                ! if the both labels are not 0, we must calculate the
                ! topo numbers for the current and the candidate label.
                FGBGTNP = TopologicalNumberFunction%EvaluateFGTNOfLabelAtIndex(coord,labels_,CandidateLabel)

                IF (FGBGTNP%FGNumber.NE.1.OR.FGBGTNP%BGNumber.NE.1) THEN
                   MCMCIsParticleTopoValid_2d=.FALSE.
                   RETURN
                ENDIF

                ContainerLabel=ABS(labels_(coord(1),coord(2)))

                FGBGTNP = TopologicalNumberFunction%EvaluateFGTNOfLabelAtIndex(ContainerLabel)

                IF (FGBGTNP%FGNumber.NE.1.OR.FGBGTNP%BGNumber.NE.1) THEN
                   MCMCIsParticleTopoValid_2d=.FALSE.
                   RETURN
                ENDIF
             ENDIF
          ENDIF !.NOT.AllowFission.OR..NOT.AllowHandles

          MCMCIsParticleTopoValid_2d=.TRUE.
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCIsParticleTopoValid_2d

        LOGICAL FUNCTION MCMCIsParticleTopoValid__2d(tmplabels,CandidateLabel)
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          USE ppm_rc_module_global, ONLy : AllowFission,AllowHandles
          USE ppm_rc_module_topologicalnumber, ONLY : FGandBGTopoNbPairType, &
          &   TopologicalNumberFunction
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER, DIMENSION(:,:), POINTER       :: tmplabels
          INTEGER,                 INTENT(IN   ) :: CandidateLabel
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(FGandBGTopoNbPairType) :: FGBGTNP

          INTEGER :: ContainerLabel
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------

          IF (.NOT.AllowFission.OR..NOT.AllowHandles) THEN
             IF (CandidateLabel.EQ.0) THEN
                ! Get the correct label to access the container of the particle
                ContainerLabel=ABS(tmplabels(2,2))

                FGBGTNP = TopologicalNumberFunction%EvaluateFGTNOfLabelAtIndex(tmplabels,ContainerLabel)

                IF (FGBGTNP%FGNumber.NE.1.OR.FGBGTNP%BGNumber.NE.1) THEN
                   MCMCIsParticleTopoValid__2d=.FALSE.
                   RETURN
                ENDIF
             ELSE
                ! if the both labels are not 0, we must calculate the
                ! topo numbers for the current and the candidate label.
                FGBGTNP = TopologicalNumberFunction%EvaluateFGTNOfLabelAtIndex(tmplabels,CandidateLabel)

                IF (FGBGTNP%FGNumber.NE.1.OR.FGBGTNP%BGNumber.NE.1) THEN
                   MCMCIsParticleTopoValid__2d=.FALSE.
                   RETURN
                ENDIF

                ContainerLabel=ABS(tmplabels(2,2))

                FGBGTNP = TopologicalNumberFunction%EvaluateFGTNOfLabelAtIndex(ContainerLabel)

                IF (FGBGTNP%FGNumber.NE.1.OR.FGBGTNP%BGNumber.NE.1) THEN
                   MCMCIsParticleTopoValid__2d=.FALSE.
                   RETURN
                ENDIF
             ENDIF
          ENDIF !.NOT.AllowFission.OR..NOT.AllowHandles

          MCMCIsParticleTopoValid__2d=.TRUE.
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCIsParticleTopoValid__2d

        LOGICAL FUNCTION MCMCIsParticleTopoValid_3d(labels_,coord,CandidateLabel)
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          USE ppm_rc_module_global, ONLy : AllowFission,AllowHandles
          USE ppm_rc_module_topologicalnumber, ONLY : FGandBGTopoNbPairType, &
          &   TopologicalNumberFunction
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: labels_
          INTEGER,             DIMENSION(:),     INTENT(IN   ) :: coord
          INTEGER,                               INTENT(IN   ) :: CandidateLabel
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(FGandBGTopoNbPairType) :: FGBGTNP

          INTEGER :: ContainerLabel

          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------

          IF (.NOT.AllowFission.OR..NOT.AllowHandles) THEN
             IF (CandidateLabel.EQ.0) THEN
                ! Get the correct label to access the container of the particle
                ContainerLabel=ABS(labels_(coord(1),coord(2),coord(3)))

                FGBGTNP = TopologicalNumberFunction%EvaluateFGTNOfLabelAtIndex(coord,labels_,ContainerLabel)

                IF (FGBGTNP%FGNumber.NE.1.OR.FGBGTNP%BGNumber.NE.1) THEN
                   MCMCIsParticleTopoValid_3d=.FALSE.
                   RETURN
                ENDIF
             ELSE
                FGBGTNP = TopologicalNumberFunction%EvaluateFGTNOfLabelAtIndex(coord,labels_,CandidateLabel)

                IF (FGBGTNP%FGNumber.NE.1.OR.FGBGTNP%BGNumber.NE.1) THEN
                   MCMCIsParticleTopoValid_3d=.FALSE.
                   RETURN
                ENDIF

                ! if the both labels are not 0, we must calculate the
                ! topo numbers for the current and the candidate label.
                ContainerLabel=ABS(labels_(coord(1),coord(2),coord(3)))

                FGBGTNP = TopologicalNumberFunction%EvaluateFGTNOfLabelAtIndex(ContainerLabel)

                IF (FGBGTNP%FGNumber.NE.1.OR.FGBGTNP%BGNumber.NE.1) THEN
                   MCMCIsParticleTopoValid_3d=.FALSE.
                   RETURN
                ENDIF
             ENDIF
          ENDIF !.NOT.AllowFission.OR..NOT.AllowHandles

          MCMCIsParticleTopoValid_3d=.TRUE.
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCIsParticleTopoValid_3d

        LOGICAL FUNCTION MCMCIsParticleTopoValid__3d(tmplabels,CandidateLabel)
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          USE ppm_rc_module_global, ONLy : AllowFission,AllowHandles
          USE ppm_rc_module_topologicalnumber, ONLY : FGandBGTopoNbPairType, &
          &   TopologicalNumberFunction
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER, DIMENSION(:,:,:), POINTER       :: tmplabels
          INTEGER,                   INTENT(IN   ) :: CandidateLabel
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(FGandBGTopoNbPairType) :: FGBGTNP

          INTEGER :: ContainerLabel

          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------

          IF (.NOT.AllowFission.OR..NOT.AllowHandles) THEN
             IF (CandidateLabel.EQ.0) THEN
                ! Get the correct label to access the container of the particle
                ContainerLabel=ABS(tmplabels(2,2,2))

                FGBGTNP = TopologicalNumberFunction%EvaluateFGTNOfLabelAtIndex(tmplabels,ContainerLabel)

                IF (FGBGTNP%FGNumber.NE.1.OR.FGBGTNP%BGNumber.NE.1) THEN
                   MCMCIsParticleTopoValid__3d=.FALSE.
                   RETURN
                ENDIF
             ELSE
                FGBGTNP = TopologicalNumberFunction%EvaluateFGTNOfLabelAtIndex(tmplabels,CandidateLabel)

                IF (FGBGTNP%FGNumber.NE.1.OR.FGBGTNP%BGNumber.NE.1) THEN
                   MCMCIsParticleTopoValid__3d=.FALSE.
                   RETURN
                ENDIF

                ! if the both labels are not 0, we must calculate the
                ! topo numbers for the current and the candidate label.
                ContainerLabel=ABS(tmplabels(2,2,2))

                FGBGTNP = TopologicalNumberFunction%EvaluateFGTNOfLabelAtIndex(ContainerLabel)

                IF (FGBGTNP%FGNumber.NE.1.OR.FGBGTNP%BGNumber.NE.1) THEN
                   MCMCIsParticleTopoValid__3d=.FALSE.
                   RETURN
                ENDIF
             ENDIF
          ENDIF !.NOT.AllowFission.OR..NOT.AllowHandles

          MCMCIsParticleTopoValid__3d=.TRUE.
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCIsParticleTopoValid__3d

        FUNCTION MCMCgetParticlesInFGNeighborhood_2d(labels_,coord,Nm,MCMCParticles,MCMCParticlesCoords) RESULT(nsize)
          !!! The method returns a list of all particles in the neighborhood INCLUSIVE
          !!! the proposal for the particle at coord, which comes at the end!
          !!! If coord is outside domain, it is still Okay
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_error_fatal
          USE ppm_module_error, ONLY : ppm_err_alloc,ppm_error

          USE ppm_rc_module_topologicalnumber, ONLY : FG_ConnectivityType
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER, CONTIGUOUS, DIMENSION(:,:),              POINTER       :: labels_
          INTEGER,             DIMENSION(:),                INTENT(IN   ) :: coord
          !!! input coord
          INTEGER,             DIMENSION(:),                INTENT(IN   ) :: Nm
          TYPE(MCMCParticle),  DIMENSION(:),   ALLOCATABLE, INTENT(INOUT) :: MCMCParticles
          INTEGER,             DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: MCMCParticlesCoords
          INTEGER                                                         :: nsize
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(MCMCParticle), DIMENSION(:), ALLOCATABLE :: tmpParticles
          TYPE(MCMCParticle), DIMENSION(9)              :: MCMCParticles_

          REAL(MK) :: proposal

          INTEGER, DIMENSION(:,:), ALLOCATABLE :: tmpcoords
          INTEGER, DIMENSION(:,:), POINTER     :: tmplabels
          INTEGER, DIMENSION(:,:), POINTER     :: labelstmp
          INTEGER, DIMENSION(2)                :: ll,ld
          INTEGER                              :: nsize_
          INTEGER                              :: i,j
          INTEGER                              :: info

          CHARACTER(LEN=*), PARAMETER :: caller="MCMCgetParticlesInFGNeighborhood"
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          nsize=0

          labelstmp => labels_(coord(1)-2:coord(1)+2,coord(2)-2:coord(2)+2)

          DO i=1,FG_ConnectivityType%NumberOfNeighbors
             ld=FG_ConnectivityType%NeighborsPoints(:,i)+coord
             IF (ANY(ld.LT.1.OR.ld.GT.Nm)) CYCLE

             ll=FG_ConnectivityType%NeighborsPoints(:,i)+3

             tmplabels => labelstmp(ll(1)-1:ll(1)+1,ll(2)-1:ll(2)+1)

             nsize_=MCMCgetRegularParticlesAtIndex(tmplabels,ld(1),ld(2),Nm,MCMCParticles_)

             IF (nsize_.GT.0) THEN
                ALLOCATE(tmpParticles(nsize+nsize_),tmpcoords(2,nsize+nsize_),STAT=info)
                or_fail_alloc("Failed to allocate tmpParticles & tmpcoords!", &
                & ppm_error=ppm_error_fatal,exit_point=no)

                FORALL (j=1:nsize)
                   tmpParticles(j)=MCMCParticles(j)
                   tmpcoords(:,j)=MCMCParticlesCoords(:,j)
                END FORALL
                FORALL (j=1:nsize_)
                   tmpParticles(nsize+j)=MCMCParticles_(j)
                   tmpcoords(:,nsize+j)=ld
                END FORALL

                CALL MOVE_ALLOC(tmpParticles,MCMCParticles)
                CALL MOVE_ALLOC(tmpcoords,MCMCParticlesCoords)

                nsize=nsize+nsize_
             ENDIF
          ENDDO !i=1,FG_ConnectivityType%NumberOfNeighbors

          ! Insert particle such that it is possible
          ! that we have a single particle move.
          tmplabels => labelstmp(2:4,2:4)
          proposal=MCMCproposal(tmplabels)

          ALLOCATE(tmpParticles(nsize+1),tmpcoords(2,nsize+1),STAT=info)
          or_fail_alloc("Failed to allocate tmpParticles & tmpcoords!", &
          & ppm_error=ppm_error_fatal,exit_point=no)

          FORALL (j=1:nsize)
             tmpParticles(j)=MCMCParticles(j)
             tmpcoords(:,j)=MCMCParticlesCoords(:,j)
          END FORALL

          tmpParticles(nsize+1)=MCMCParticle(-1,proposal)
          tmpcoords(:,nsize+1)=coord

          nsize=nsize+1

          CALL MOVE_ALLOC(tmpParticles,MCMCParticles)
          CALL MOVE_ALLOC(tmpcoords,MCMCParticlesCoords)
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCgetParticlesInFGNeighborhood_2d

        FUNCTION MCMCgetParticlesInFGNeighborhood_3d(labels_,coord,Nm,MCMCParticles,MCMCParticlesCoords) RESULT(nsize)
          !!! The method returns a list of all particles in the neighborhood INCLUSIVE
          !!! the particles at coord, which comes at the end!
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_error_fatal
          USE ppm_module_error, ONLY : ppm_err_alloc,ppm_error

          USE ppm_rc_module_topologicalnumber, ONLY : FG_ConnectivityType
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER, CONTIGUOUS, DIMENSION(:,:,:),            POINTER       :: labels_
          INTEGER,             DIMENSION(:),                INTENT(IN   ) :: coord
          INTEGER,             DIMENSION(:),                INTENT(IN   ) :: Nm
          TYPE(MCMCParticle),  DIMENSION(:),   ALLOCATABLE, INTENT(INOUT) :: MCMCParticles
          INTEGER,             DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: MCMCParticlesCoords
          INTEGER                                                         :: nsize
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(MCMCParticle), DIMENSION(:), ALLOCATABLE :: tmpParticles
          TYPE(MCMCParticle), DIMENSION(27)             :: MCMCParticles_

          REAL(MK) :: proposal

          INTEGER, DIMENSION(:,:),   ALLOCATABLE :: tmpcoords
          INTEGER, DIMENSION(:,:,:), POINTER     :: tmplabels
          INTEGER, DIMENSION(:,:,:), POINTER     :: labelstmp
          INTEGER, DIMENSION(3)                  :: ll,ld
          INTEGER                                :: nsize_
          INTEGER                                :: i,j
          INTEGER                                :: info

          CHARACTER(LEN=*), PARAMETER :: caller="MCMCgetParticlesInFGNeighborhood"
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          nsize=0

          labelstmp => labels_(coord(1)-2:coord(1)+2,coord(2)-2:coord(2)+2,coord(3)-2:coord(3)+2)

          DO i=1,FG_ConnectivityType%NumberOfNeighbors
             ld=FG_ConnectivityType%NeighborsPoints(:,i)+coord

             IF (ANY(ld.LT.1.OR.ld.GT.Nm)) CYCLE

             ll=FG_ConnectivityType%NeighborsPoints(:,i)+3

             tmplabels => labelstmp(ll(1)-1:ll(1)+1,ll(2)-1:ll(2)+1,ll(3)-1:ll(3)+1)

             nsize_=MCMCgetRegularParticlesAtIndex(tmplabels,ld(1),ld(2),ld(3),Nm,MCMCParticles_)

             IF (nsize_.GT.0) THEN
                ALLOCATE(tmpParticles(nsize+nsize_),tmpcoords(3,nsize+nsize_),STAT=info)
                or_fail_alloc("Failed to allocate tmpParticles & tmpcoords!", &
                & ppm_error=ppm_error_fatal,exit_point=no)

                FORALL (j=1:nsize)
                   tmpParticles(j)=MCMCParticles(j)
                   tmpcoords(:,j)=MCMCParticlesCoords(:,j)
                END FORALL
                FORALL (j=1:nsize_)
                   tmpParticles(nsize+j)=MCMCParticles_(j)
                   tmpcoords(:,nsize+j)=ld
                END FORALL

                CALL MOVE_ALLOC(tmpParticles,MCMCParticles)
                CALL MOVE_ALLOC(tmpcoords,MCMCParticlesCoords)

                nsize=nsize+nsize_
             ENDIF
          ENDDO !i=1,FG_ConnectivityType%NumberOfNeighbors

          ! Insert particle such that it is possible
          ! that we have a single particle move.
          tmplabels => labelstmp(2:4,2:4,2:4)
          proposal=MCMCproposal(tmplabels)

          ALLOCATE(tmpParticles(nsize+1),tmpcoords(3,nsize+1),STAT=info)
          or_fail_alloc("Failed to allocate tmpParticles & tmpcoords!", &
          & ppm_error=ppm_error_fatal,exit_point=no)

          FORALL (j=1:nsize)
             tmpParticles(j)=MCMCParticles(j)
             tmpcoords(:,j)=MCMCParticlesCoords(:,j)
          END FORALL

          tmpParticles(nsize+1)=MCMCParticle(-1,proposal)
          tmpcoords(:,nsize+1)=coord

          nsize=nsize+1

          CALL MOVE_ALLOC(tmpParticles,MCMCParticles)
          CALL MOVE_ALLOC(tmpcoords,MCMCParticlesCoords)
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCgetParticlesInFGNeighborhood_3d

        FUNCTION MCMCgetParticlesInBGNeighborhood_2d(labels_,coord,Nm,MCMCParticles,MCMCParticlesCoords) RESULT(nsize)
          !!! The method returns a list of all particles in the neighborhood INCLUSIVE
          !!! the particles at coord, which comes at the end!
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_error_fatal
          USE ppm_module_error, ONLY : ppm_err_alloc,ppm_error

          USE ppm_rc_module_global, ONLY : AllowFission,AllowHandles
          USE ppm_rc_module_topologicalnumber, ONLY : BG_ConnectivityType
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER, CONTIGUOUS, DIMENSION(:,:),              POINTER       :: labels_
          INTEGER,             DIMENSION(:),                INTENT(IN   ) :: coord
          INTEGER,             DIMENSION(:),                INTENT(IN   ) :: Nm
          TYPE(MCMCParticle),  DIMENSION(:),   ALLOCATABLE, INTENT(INOUT) :: MCMCParticles
          INTEGER,             DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: MCMCParticlesCoords
          INTEGER                                                         :: nsize
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(MCMCParticle), DIMENSION(:), ALLOCATABLE :: tmpParticles
          TYPE(MCMCParticle), DIMENSION(9)              :: MCMCParticles_

          REAL(MK) :: proposal

          INTEGER, DIMENSION(:,:), ALLOCATABLE :: tmpcoords
          INTEGER, DIMENSION(:,:), POINTER     :: tmplabels
          INTEGER, DIMENSION(:,:), POINTER     :: labelstmp
          INTEGER, DIMENSION(2)                :: ll,ld
          INTEGER                              :: nsize_
          INTEGER                              :: i,j
          INTEGER                              :: info

          CHARACTER(LEN=*), PARAMETER :: caller="MCMCgetParticlesInBGNeighborhood"
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          nsize=0

          labelstmp => labels_(coord(1)-2:coord(1)+2,coord(2)-2:coord(2)+2)

          DO i=1,BG_ConnectivityType%NumberOfNeighbors
             ld=BG_ConnectivityType%NeighborsPoints(:,i)+coord

             IF (ANY(ld.LT.1.OR.ld.GT.Nm)) CYCLE

             ll=BG_ConnectivityType%NeighborsPoints(:,i)+3

             tmplabels => labelstmp(ll(1)-1:ll(1)+1,ll(2)-1:ll(2)+1)

             nsize_=MCMCgetRegularParticlesAtIndex(tmplabels,ld(1),ld(2),Nm,MCMCParticles_)

             IF (nsize_.GT.0) THEN
                IF (.NOT.AllowFission.OR..NOT.AllowHandles) THEN
                   !-------------------------------------------------------------------------
                   ! filter the points to meet topological constraints:
                   !-------------------------------------------------------------------------
                   DO j=nsize_,1,-1
                      !-------------------------------------------------------------------------
                      ! Removes all non-simple points from the particle set
                      !-------------------------------------------------------------------------
                      IF (.NOT.MCMCIsParticleTopoValid(tmplabels,MCMCParticles_(j)%candlabel)) THEN
                         MCMCParticles_(j)=MCMCParticles_(nsize_)
                         nsize_=nsize_-1
                      ENDIF
                   ENDDO
                ENDIF
             ENDIF !nsize_.GT.0

             IF (nsize_.GT.0) THEN
                ALLOCATE(tmpParticles(nsize+nsize_),tmpcoords(2,nsize+nsize_),STAT=info)
                or_fail_alloc("Failed to allocate tmpParticles & tmpcoords!", &
                & ppm_error=ppm_error_fatal,exit_point=no)

                FORALL (j=1:nsize)
                   tmpParticles(j)=MCMCParticles(j)
                   tmpcoords(:,j)=MCMCParticlesCoords(:,j)
                END FORALL
                FORALL (j=1:nsize_)
                   tmpParticles(nsize+j)=MCMCParticles_(j)
                   tmpcoords(:,nsize+j)=ld
                END FORALL

                CALL MOVE_ALLOC(tmpParticles,MCMCParticles)
                CALL MOVE_ALLOC(tmpcoords,MCMCParticlesCoords)

                nsize=nsize+nsize_
             ENDIF
          ENDDO !i=1,BG_ConnectivityType%NumberOfNeighbors

          !-------------------------------------------------------------------------
          ! Insert particle A in the set for B such that it is possible
          ! that we have a single particle move.
          !-------------------------------------------------------------------------
          tmplabels => labelstmp(2:4,2:4)
          proposal=MCMCproposal(tmplabels)

          ALLOCATE(tmpParticles(nsize+1),tmpcoords(2,nsize+1),STAT=info)
          or_fail_alloc("Failed to allocate tmpParticles & tmpcoords!", &
          & ppm_error=ppm_error_fatal,exit_point=no)

          FORALL (j=1:nsize)
             tmpParticles(j)=MCMCParticles(j)
             tmpcoords(:,j)=MCMCParticlesCoords(:,j)
          END FORALL

          tmpParticles(nsize+1)=MCMCParticle(-1,proposal)
          tmpcoords(:,nsize+1)=coord

          nsize=nsize+1

          CALL MOVE_ALLOC(tmpParticles,MCMCParticles)
          CALL MOVE_ALLOC(tmpcoords,MCMCParticlesCoords)

          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCgetParticlesInBGNeighborhood_2d

        FUNCTION MCMCgetParticlesInBGNeighborhood_3d(labels_,coord,Nm,MCMCParticles,MCMCParticlesCoords) RESULT(nsize)
          !!! The method returns a list of all particles in the neighborhood INCLUSIVE
          !!! the particles at coord, which comes at the end!
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_error_fatal
          USE ppm_module_error, ONLY : ppm_err_alloc,ppm_error

          USE ppm_rc_module_global, ONLY : AllowFission,AllowHandles
          USE ppm_rc_module_topologicalnumber, ONLY : BG_ConnectivityType
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER, CONTIGUOUS, DIMENSION(:,:,:),            POINTER       :: labels_
          INTEGER,             DIMENSION(:),                INTENT(IN   ) :: coord
          INTEGER,             DIMENSION(:),                INTENT(IN   ) :: Nm
          TYPE(MCMCParticle),  DIMENSION(:),   ALLOCATABLE, INTENT(INOUT) :: MCMCParticles
          INTEGER,             DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: MCMCParticlesCoords
          INTEGER                                                         :: nsize
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(MCMCParticle), DIMENSION(:), ALLOCATABLE :: tmpParticles
          TYPE(MCMCParticle), DIMENSION(27)             :: MCMCParticles_

          REAL(MK) :: proposal

          INTEGER, DIMENSION(:,:),   ALLOCATABLE :: tmpcoords
          INTEGER, DIMENSION(:,:,:), POINTER     :: tmplabels
          INTEGER, DIMENSION(:,:,:), POINTER     :: labelstmp
          INTEGER, DIMENSION(3)                  :: ll,ld
          INTEGER                                :: nsize_
          INTEGER                                :: i,j
          INTEGER                                :: info

          CHARACTER(LEN=*), PARAMETER :: caller="MCMCgetParticlesInBGNeighborhood"
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          nsize=0

          labelstmp => labels_(coord(1)-2:coord(1)+2,coord(2)-2:coord(2)+2,coord(3)-2:coord(3)+2)

          DO i=1,BG_ConnectivityType%NumberOfNeighbors
             ld=BG_ConnectivityType%NeighborsPoints(:,i)+coord

             IF (ANY(ld.LT.1.OR.ld.GT.Nm)) CYCLE

             ll=BG_ConnectivityType%NeighborsPoints(:,i)+3

             tmplabels => labelstmp(ll(1)-1:ll(1)+1,ll(2)-1:ll(2)+1,ll(3)-1:ll(3)+1)

             nsize_=MCMCgetRegularParticlesAtIndex(tmplabels,ld(1),ld(2),ld(3),Nm,MCMCParticles_)

             IF (nsize_.GT.0) THEN
                IF (.NOT.AllowFission.OR..NOT.AllowHandles) THEN
                ! filter the points to meet topological constraints:
                   DO j=nsize_,1,-1
                      ! Removes all non-simple points from the particle set
                      IF (.NOT.MCMCIsParticleTopoValid(tmplabels,MCMCParticles_(j)%candlabel)) THEN
                         MCMCParticles_(j)=MCMCParticles_(nsize_)
                         nsize_=nsize_-1
                      ENDIF
                   ENDDO
                ENDIF
             ENDIF !nsize_.GT.0

             IF (nsize_.GT.0) THEN
                ALLOCATE(tmpParticles(nsize+nsize_),tmpcoords(3,nsize+nsize_),STAT=info)
                or_fail_alloc("Failed to allocate tmpParticles & tmpcoords!", &
                & ppm_error=ppm_error_fatal,exit_point=no)

                FORALL (j=1:nsize)
                   tmpParticles(j)=MCMCParticles(j)
                   tmpcoords(:,j)=MCMCParticlesCoords(:,j)
                END FORALL
                FORALL (j=1:nsize_)
                   tmpParticles(nsize+j)=MCMCParticles_(j)
                   tmpcoords(:,nsize+j)=ld
                END FORALL

                CALL MOVE_ALLOC(tmpParticles,MCMCParticles)
                CALL MOVE_ALLOC(tmpcoords,MCMCParticlesCoords)

                nsize=nsize+nsize_
             ENDIF
          ENDDO !i=1,BG_ConnectivityType%NumberOfNeighbors

          ! Insert particle A in the set for B such that it is possible
          ! that we have a single particle move.
          tmplabels => labelstmp(2:4,2:4,2:4)
          proposal=MCMCproposal(tmplabels)

          ALLOCATE(tmpParticles(nsize+1),tmpcoords(3,nsize+1),STAT=info)
          or_fail_alloc("Failed to allocate tmpParticles & tmpcoords!", &
          & ppm_error=ppm_error_fatal,exit_point=no)

          FORALL (j=1:nsize)
             tmpParticles(j)=MCMCParticles(j)
             tmpcoords(:,j)=MCMCParticlesCoords(:,j)
          END FORALL

          tmpParticles(nsize+1)=MCMCParticle(-1,proposal)
          tmpcoords(:,nsize+1)=coord

          nsize=nsize+1

          CALL MOVE_ALLOC(tmpParticles,MCMCParticles)
          CALL MOVE_ALLOC(tmpcoords,MCMCParticlesCoords)

          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCgetParticlesInBGNeighborhood_3d

        FUNCTION MCMCGetIndexFromEdgeDensity_2d(labels_,Nm) RESULT(index_)
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          !!!
          USE ppm_module_data, ONLY : ppm_kind_int64

          USE ppm_rc_module_global, ONLY : FORBIDDEN
          USE ppm_rc_module_rnd, ONLY : ppm_rc_GetImageDistrIndex
          USE ppm_rc_module_util, ONLY : id_gtl_2d
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER, CONTIGUOUS, DIMENSION(:,:), POINTER       :: labels_
          INTEGER,             DIMENSION(:),   INTENT(IN   ) :: Nm
          INTEGER,             DIMENSION(2)                  :: index_
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          INTEGER(ppm_kind_int64) :: lindex

          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          DO
             lindex=ppm_rc_GetImageDistrIndex()
             CALL id_gtl_2d(Nm,lindex,index_)
             IF (labels_(index_(1),index_(2)).NE.FORBIDDEN) EXIT
          ENDDO
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCGetIndexFromEdgeDensity_2d

        FUNCTION MCMCGetIndexFromEdgeDensity_3d(labels_,Nm) RESULT(index_)
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          !!!
          USE ppm_module_data, ONLY : ppm_kind_int64

          USE ppm_rc_module_global, ONLY : FORBIDDEN
          USE ppm_rc_module_rnd, ONLY : ppm_rc_GetImageDistrIndex
          USE ppm_rc_module_util, ONLY : id_gtl_3d
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: labels_
          INTEGER,             DIMENSION(:),     INTENT(IN   ) :: Nm
          INTEGER,             DIMENSION(3)                    :: index_
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          INTEGER(ppm_kind_int64) :: lindex

          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          DO
             lindex=ppm_rc_GetImageDistrIndex()
             CALL id_gtl_3d(Nm,lindex,index_)
             IF (labels_(index_(1),index_(2),index_(3)).NE.FORBIDDEN) EXIT
          ENDDO
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCGetIndexFromEdgeDensity_3d

        FUNCTION MCMCInsertCandidatesToContainers_2d(coord1,coord2,cbox, &
        &        MCMCParticle_,CurrentLabel,DoRecord,ParticleReplaced) RESULT (info)
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_fatal,ppm_error_error
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_err_alloc,ppm_error,ppm_err_argument, &
          &   ppm_err_sub_failed

          USE ppm_rc_module_linkedlist, ONLY : ppm_rc_list,MCMCParticleInContainerHistory
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER,            INTENT(IN   ) :: coord1
          INTEGER,            INTENT(IN   ) :: coord2
          INTEGER,            INTENT(IN   ) :: cbox
          TYPE(MCMCParticle), INTENT(IN   ) :: MCMCParticle_
          INTEGER,            INTENT(IN   ) :: CurrentLabel
          LOGICAL,            INTENT(IN   ) :: DoRecord
          LOGICAL,            INTENT(  OUT) :: ParticleReplaced
          INTEGER                           :: info
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(MCMCParticle) :: tmpParticle
          TYPE(MCMCHistoryParticle) :: tmpHistoryParticle

          TYPE(ppm_rc_list), POINTER ::  seed

          REAL(ppm_kind_double) :: t0
          REAL(MK)              :: ProposalDiff

          CHARACTER(LEN=ppm_char) :: caller='MCMCInsertCandidatesToContainers'
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          ParticleReplaced=.FALSE.

          tmpParticle=MCMCRegularParticlesInCell(cbox)%search(coord1,coord2,MCMCParticle_%candlabel)

          CALL MCMCRegularParticlesInCell(cbox)%insert(coord1,coord2,MCMCParticle_,info)
          or_fail("MCMCRegularParticlesInCell(cbox)%insert")

          IF (tmpParticle%candlabel.EQ.-1) THEN
             !-------------------------------------------------------------------------
             ! The particle does not exist
             !-------------------------------------------------------------------------
             MCMCTotalNormalizerlocal=MCMCTotalNormalizerlocal+REAL(MCMCParticle_%proposal,ppm_kind_double)
          ELSE
             !-------------------------------------------------------------------------
             ! this is a replacement:
             !-------------------------------------------------------------------------
             ProposalDiff=MCMCParticle_%proposal-tmpParticle%proposal

             MCMCTotalNormalizerlocal=MCMCTotalNormalizerlocal+REAL(ProposalDiff,ppm_kind_double)

             ParticleReplaced=.TRUE.
          ENDIF

          IF (DoRecord) THEN
             !-------------------------------------------------------------------------
             ! Note that the order here is important: first insert the particle
             ! that gets replaced. When restoring the state afterwards, the
             ! particle history is iterated in reverse order.
             !-------------------------------------------------------------------------
             IF (ParticleReplaced) THEN
                !-------------------------------------------------------------------------
                ! This is a replaced particle (tmpParticle%candlabel.eq.MCMCParticle_%candlabel)
                !-------------------------------------------------------------------------
                tmpHistoryParticle=MCMCHistoryParticle(tmpParticle%candlabel, &
                & CurrentLabel,tmpParticle%proposal,.FALSE.)
                ALLOCATE(seed,STAT=info)
                or_fail_alloc("seed")
                CALL seed%add(coord1,coord2,tmpHistoryParticle)
                CALL MCMCParticleInContainerHistory%push(seed,info)
                or_fail("MCMCParticleInContainerHistory%push")
             ENDIF

             tmpHistoryParticle=MCMCHistoryParticle(MCMCParticle_%candlabel, &
             & CurrentLabel,MCMCParticle_%proposal,.TRUE.)
             ALLOCATE(seed,STAT=info)
             or_fail_alloc("seed")
             CALL seed%add(coord1,coord2,tmpHistoryParticle)
             CALL MCMCParticleInContainerHistory%push(seed,info)
             or_fail("MCMCParticleInContainerHistory%push")
          ENDIF
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        END FUNCTION MCMCInsertCandidatesToContainers_2d

        FUNCTION MCMCInsertCandidatesToContainers_3d(coord1,coord2,coord3,cbox, &
        &        MCMCParticle_,CurrentLabel,DoRecord,ParticleReplaced) RESULT(info)
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_fatal,ppm_error_error
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_err_alloc,ppm_error,ppm_err_argument, &
          &   ppm_err_sub_failed

          USE ppm_rc_module_linkedlist, ONLY : ppm_rc_list,MCMCParticleInContainerHistory
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER,            INTENT(IN   ) :: coord1
          INTEGER,            INTENT(IN   ) :: coord2
          INTEGER,            INTENT(IN   ) :: coord3
          INTEGER,            INTENT(IN   ) :: cbox
          TYPE(MCMCParticle), INTENT(IN   ) :: MCMCParticle_
          INTEGER,            INTENT(IN   ) :: CurrentLabel
          LOGICAL,            INTENT(IN   ) :: DoRecord
          LOGICAL,            INTENT(  OUT) :: ParticleReplaced
          INTEGER                           :: info
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(MCMCParticle) :: tmpParticle

          TYPE(MCMCHistoryParticle) :: tmpHistoryParticle

          TYPE(ppm_rc_list), POINTER ::  seed

          REAL(ppm_kind_double) :: t0
          REAL(MK)              :: ProposalDiff

          CHARACTER(LEN=ppm_char) :: caller='MCMCInsertCandidatesToContainers'
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          ParticleReplaced=.FALSE.

          tmpParticle=MCMCRegularParticlesInCell(cbox)%search(coord1,coord2,coord3,MCMCParticle_%candlabel)

          CALL MCMCRegularParticlesInCell(cbox)%insert(coord1,coord2,coord3,MCMCParticle_,info)
          or_fail("MCMCRegularParticlesInCell(cbox)%insert")

          IF (tmpParticle%candlabel.EQ.-1) THEN
             !-------------------------------------------------------------------------
             ! The particle does not exist
             !-------------------------------------------------------------------------
             MCMCTotalNormalizerlocal=MCMCTotalNormalizerlocal+REAL(MCMCParticle_%proposal,ppm_kind_double)
          ELSE
             !-------------------------------------------------------------------------
             ! this is a replacement:
             !-------------------------------------------------------------------------
             ProposalDiff=MCMCParticle_%proposal-tmpParticle%proposal

             MCMCTotalNormalizerlocal=MCMCTotalNormalizerlocal+REAL(ProposalDiff,ppm_kind_double)

             ParticleReplaced=.TRUE.
          ENDIF

          IF (DoRecord) THEN
             !-------------------------------------------------------------------------
             ! Note that the order here is important: first insert the particle
             ! that gets replaced. When restoring the state afterwards, the
             ! particle history is iterated in reverse order.
             !-------------------------------------------------------------------------
             IF (ParticleReplaced) THEN
                !-------------------------------------------------------------------------
                ! This is a replacement (tmpParticle%candlabel.eq.MCMCParticle_%candlabel)
                !-------------------------------------------------------------------------
                tmpHistoryParticle=MCMCHistoryParticle(tmpParticle%candlabel, &
                & CurrentLabel,tmpParticle%proposal,.FALSE.)
                ALLOCATE(seed,STAT=info)
                or_fail_alloc("seed")
                CALL seed%add(coord1,coord2,coord3,tmpHistoryParticle)
                CALL MCMCParticleInContainerHistory%push(seed,info)
                or_fail("MCMCParticleInContainerHistory%push")
             ENDIF

             tmpHistoryParticle=MCMCHistoryParticle(MCMCParticle_%candlabel, &
             & CurrentLabel,MCMCParticle_%proposal,.TRUE.)
             ALLOCATE(seed,STAT=info)
             or_fail_alloc("seed")
             CALL seed%add(coord1,coord2,coord3,tmpHistoryParticle)
             CALL MCMCParticleInContainerHistory%push(seed,info)
             or_fail("MCMCParticleInContainerHistory%push")
          ENDIF

        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        END FUNCTION MCMCInsertCandidatesToContainers_3d

        FUNCTION MCMCEraseCandidatesFromContainers_2d(coord1,coord2,cbox, &
        &        MCMCParticle_,CurrentLabel,DoRecord) RESULT (info)
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_fatal,ppm_error_error
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_err_alloc,ppm_error,ppm_err_argument, &
          &   ppm_err_sub_failed

          USE ppm_rc_module_linkedlist, ONLY : ppm_rc_list,MCMCParticleInContainerHistory
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER,            INTENT(IN   ) :: coord1
          INTEGER,            INTENT(IN   ) :: coord2
          INTEGER,            INTENT(IN   ) :: cbox
          TYPE(MCMCParticle), INTENT(IN   ) :: MCMCParticle_
          INTEGER,            INTENT(IN   ) :: CurrentLabel
          LOGICAL,            INTENT(IN   ) :: DoRecord
          INTEGER                           :: info
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(MCMCParticle) :: tmpParticle
          TYPE(MCMCHistoryParticle) :: tmpHistoryParticle

          TYPE(ppm_rc_list), POINTER ::  seed

          REAL(ppm_kind_double) :: t0

          CHARACTER(LEN=ppm_char) :: caller='MCMCEraseCandidatesFromContainers'
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          tmpParticle=MCMCRegularParticlesInCell(cbox)%search(coord1,coord2,MCMCParticle_%candlabel)

          IF (tmpParticle%candlabel.GT.-1) THEN
             !-------------------------------------------------------------------------
             ! Particle exists
             !-------------------------------------------------------------------------
             CALL MCMCRegularParticlesInCell(cbox)%remove(coord1,coord2,MCMCParticle_%candlabel,info,.TRUE.)
             or_fail("MCMCRegularParticlesInCell(cbox)%remove")

             MCMCTotalNormalizerlocal=MCMCTotalNormalizerlocal-REAL(tmpParticle%proposal,ppm_kind_double)

             IF (DoRecord) THEN
                tmpHistoryParticle=MCMCHistoryParticle(tmpParticle%candlabel, &
                & CurrentLabel,tmpParticle%proposal,.FALSE.)
                ALLOCATE(seed,STAT=info)
                or_fail_alloc("seed")
                CALL seed%add(coord1,coord2,tmpHistoryParticle)
                CALL MCMCParticleInContainerHistory%push(seed,info)
                or_fail("MCMCParticleInContainerHistory%push")
             ENDIF
          ENDIF !tmpParticle%candlabel.GT.-1
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        END FUNCTION MCMCEraseCandidatesFromContainers_2d

        FUNCTION MCMCEraseCandidatesFromContainers_3d(coord1,coord2,coord3,cbox, &
        &        MCMCParticle_,CurrentLabel,DoRecord) RESULT (info)
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_fatal,ppm_error_error
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_err_alloc,ppm_error,ppm_err_argument, &
          &   ppm_err_sub_failed

          USE ppm_rc_module_linkedlist, ONLY : ppm_rc_list,MCMCParticleInContainerHistory
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER,            INTENT(IN   ) :: coord1
          INTEGER,            INTENT(IN   ) :: coord2
          INTEGER,            INTENT(IN   ) :: coord3
          INTEGER,            INTENT(IN   ) :: cbox
          TYPE(MCMCParticle), INTENT(IN   ) :: MCMCParticle_
          INTEGER,            INTENT(IN   ) :: CurrentLabel
          LOGICAL,            INTENT(IN   ) :: DoRecord
          INTEGER                           :: info
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(MCMCParticle) :: tmpParticle
          TYPE(MCMCHistoryParticle) :: tmpHistoryParticle

          TYPE(ppm_rc_list), POINTER ::  seed

          REAL(ppm_kind_double) :: t0

          CHARACTER(LEN=ppm_char) :: caller='MCMCEraseCandidatesFromContainers'
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          tmpParticle=MCMCRegularParticlesInCell(cbox)%search(coord1,coord2,coord3,MCMCParticle_%candlabel)

          IF (tmpParticle%candlabel.GT.-1) THEN
             ! Particle exists
             CALL MCMCRegularParticlesInCell(cbox)%remove(coord1,coord2,coord3,MCMCParticle_%candlabel,info,.TRUE.)
             or_fail("MCMCRegularParticlesInCell(cbox)%remove")

             MCMCTotalNormalizerlocal=MCMCTotalNormalizerlocal-REAL(tmpParticle%proposal,ppm_kind_double)

             IF (DoRecord) THEN
                tmpHistoryParticle=MCMCHistoryParticle(tmpParticle%candlabel, &
                & CurrentLabel,tmpParticle%proposal,.FALSE.)
                ALLOCATE(seed,STAT=info)
                or_fail_alloc("seed")
                CALL seed%add(coord1,coord2,coord3,tmpHistoryParticle)
                CALL MCMCParticleInContainerHistory%push(seed,info)
                or_fail("MCMCParticleInContainerHistory%push")
             ENDIF
          ENDIF !tmpParticle%candlabel.GT.-1
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        END FUNCTION MCMCEraseCandidatesFromContainers_3d

        LOGICAL FUNCTION MCMCInsertFloatingParticle_2d(coord1,coord2,cbox, &
        &                MCMCParticle_,DoRecord,MCMCParticleExist)
          !!! Insert a floating particle only if it doesn't exist.
          !!! If it is existed, the return value will be false.
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_fatal
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_error,ppm_err_alloc,ppm_err_sub_failed

          USE ppm_rc_module_linkedlist, ONLY : ppm_rc_list, &
          &   MCMCFloatingParticleInContainerHistory
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER,            INTENT(IN   ) :: coord1
          INTEGER,            INTENT(IN   ) :: coord2
          INTEGER,            INTENT(IN   ) :: cbox
          TYPE(MCMCParticle), INTENT(IN   ) :: MCMCParticle_
          LOGICAL,            INTENT(IN   ) :: DoRecord
          LOGICAL, OPTIONAL,  INTENT(IN   ) :: MCMCParticleExist
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(MCMCParticle) :: tmpParticle
          TYPE(MCMCHistoryParticle) :: tmpHistoryParticle

          TYPE(ppm_rc_list), POINTER ::  seed

          REAL(ppm_kind_double) :: t0

          INTEGER :: info

          CHARACTER(LEN=ppm_char) :: caller="MCMCInsertFloatingParticle"
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          IF (PRESENT(MCMCParticleExist)) THEN
             ! the particle already exists
             IF (MCMCParticleExist) THEN
                MCMCInsertFloatingParticle_2d=.FALSE.
                GOTO 9999
             ! there is no floating particle at this position
             ELSE
                CALL MCMCFloatingParticlesInCell(cbox)%insert(coord1,coord2,MCMCParticle_,info)
                or_fail("MCMCFloatingParticlesInCell(cbox)%insert",ppm_error=ppm_error_fatal)

                IF (DoRecord) THEN
                   tmpHistoryParticle=MCMCHistoryParticle(MCMCParticle_%candlabel, &
                   & 0,MCMCParticle_%proposal,.TRUE.)
                   ALLOCATE(seed,STAT=info)
                   or_fail_alloc("seed",ppm_error=ppm_error_fatal)
                   CALL seed%add(coord1,coord2,tmpHistoryParticle)
                   CALL MCMCFloatingParticleInContainerHistory%push(seed,info)
                   or_fail("MCMCFloatingParticleInContainerHistory%push",ppm_error=ppm_error_fatal)
                ENDIF

                MCMCTotalNormalizerlocal=MCMCTotalNormalizerlocal+REAL(MCMCParticle_%proposal,ppm_kind_double)

                MCMCInsertFloatingParticle_2d=.TRUE.
                GOTO 9999
             ENDIF
          ENDIF

          tmpParticle=MCMCFloatingParticlesInCell(cbox)%search(coord1,coord2,MCMCParticle_%candlabel)
          ! there is no floating particle at this position
          IF (tmpParticle%candlabel.EQ.-1) THEN
             CALL MCMCFloatingParticlesInCell(cbox)%insert(coord1,coord2,MCMCParticle_,info)
             or_fail("MCMCFloatingParticlesInCell(cbox)%insert",ppm_error=ppm_error_fatal)

             IF (DoRecord) THEN
                tmpHistoryParticle=MCMCHistoryParticle(MCMCParticle_%candlabel, &
                & 0,MCMCParticle_%proposal,.TRUE.)
                ALLOCATE(seed,STAT=info)
                or_fail_alloc("seed",ppm_error=ppm_error_fatal)
                CALL seed%add(coord1,coord2,tmpHistoryParticle)
                CALL MCMCFloatingParticleInContainerHistory%push(seed,info)
                or_fail("MCMCFloatingParticleInContainerHistory%push",ppm_error=ppm_error_fatal)
             ENDIF

             MCMCTotalNormalizerlocal=MCMCTotalNormalizerlocal+REAL(MCMCParticle_%proposal,ppm_kind_double)

             MCMCInsertFloatingParticle_2d=.TRUE.
          ! the particle already exists
          ELSE
             MCMCInsertFloatingParticle_2d=.FALSE.
          ENDIF
        9999 CONTINUE
          CALL substop(caller,t0,info)
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCInsertFloatingParticle_2d

        LOGICAL FUNCTION MCMCInsertFloatingParticle_3d(coord1,coord2,coord3,cbox, &
        &                MCMCParticle_,DoRecord,MCMCParticleExist)
          !!! Insert a floating particle only if it doesn't exist.
          !!! If it is existed, the return value will be false.
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_fatal
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_error,ppm_err_alloc,ppm_err_sub_failed

          USE ppm_rc_module_linkedlist, ONLY : ppm_rc_list, &
          &   MCMCFloatingParticleInContainerHistory
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER,            INTENT(IN   ) :: coord1
          INTEGER,            INTENT(IN   ) :: coord2
          INTEGER,            INTENT(IN   ) :: coord3
          INTEGER,            INTENT(IN   ) :: cbox
          TYPE(MCMCParticle), INTENT(IN   ) :: MCMCParticle_
          LOGICAL,            INTENT(IN   ) :: DoRecord
          LOGICAL, OPTIONAL,  INTENT(IN   ) :: MCMCParticleExist
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(MCMCParticle) :: tmpParticle
          TYPE(MCMCHistoryParticle) :: tmpHistoryParticle

          TYPE(ppm_rc_list), POINTER ::  seed

          REAL(ppm_kind_double) :: t0

          INTEGER :: info

          CHARACTER(LEN=ppm_char) :: caller="MCMCInsertFloatingParticle"
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          IF (PRESENT(MCMCParticleExist)) THEN
             ! the particle already exists
             IF (MCMCParticleExist) THEN
                MCMCInsertFloatingParticle_3d=.FALSE.
                GOTO 9999
             ! there is no floating particle at this position
             ELSE
                CALL MCMCFloatingParticlesInCell(cbox)%insert(coord1,coord2,coord3,MCMCParticle_,info)
                or_fail("MCMCFloatingParticlesInCell(cbox)%insert",ppm_error=ppm_error_fatal)

                IF (DoRecord) THEN
                   tmpHistoryParticle=MCMCHistoryParticle(MCMCParticle_%candlabel, &
                   & 0,MCMCParticle_%proposal,.TRUE.)
                   ALLOCATE(seed,STAT=info)
                   or_fail_alloc("seed",ppm_error=ppm_error_fatal)
                   CALL seed%add(coord1,coord2,coord3,tmpHistoryParticle)
                   CALL MCMCFloatingParticleInContainerHistory%push(seed,info)
                   or_fail("MCMCFloatingParticleInContainerHistory%push",ppm_error=ppm_error_fatal)
                ENDIF

                MCMCTotalNormalizerlocal=MCMCTotalNormalizerlocal+REAL(MCMCParticle_%proposal,ppm_kind_double)

                MCMCInsertFloatingParticle_3d=.TRUE.
                GOTO 9999
             ENDIF
          ENDIF

          tmpParticle=MCMCFloatingParticlesInCell(cbox)%search(coord1,coord2,coord3,MCMCParticle_%candlabel)
          ! there is no floating particle at this position
          IF (tmpParticle%candlabel.EQ.-1) THEN
             CALL MCMCFloatingParticlesInCell(cbox)%insert(coord1,coord2,coord3,MCMCParticle_,info)
             or_fail("MCMCFloatingParticlesInCell(cbox)%insert",ppm_error=ppm_error_fatal)

             IF (DoRecord) THEN
                tmpHistoryParticle=MCMCHistoryParticle(MCMCParticle_%candlabel, &
                & 0,MCMCParticle_%proposal,.TRUE.)
                ALLOCATE(seed,STAT=info)
                or_fail_alloc("seed",ppm_error=ppm_error_fatal)
                CALL seed%add(coord1,coord2,coord3,tmpHistoryParticle)
                CALL MCMCFloatingParticleInContainerHistory%push(seed,info)
                or_fail("MCMCFloatingParticleInContainerHistory%push",ppm_error=ppm_error_fatal)
             ENDIF

             MCMCTotalNormalizerlocal=MCMCTotalNormalizerlocal+REAL(MCMCParticle_%proposal,ppm_kind_double)

             MCMCInsertFloatingParticle_3d=.TRUE.
          ! the particle already exists
          ELSE
             MCMCInsertFloatingParticle_3d=.FALSE.
          ENDIF
        9999 CONTINUE
          CALL substop(caller,t0,info)
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCInsertFloatingParticle_3d

        FUNCTION MCMCInsertFloatingParticleCumulative_2d(coord1,coord2,cbox,MCMCParticle_,DoRecord) RESULT(info)
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_fatal
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_error,ppm_err_alloc,ppm_err_sub_failed

          USE ppm_rc_module_linkedlist, ONLY : ppm_rc_list, &
          &   MCMCFloatingParticleInContainerHistory
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER,            INTENT(IN   ) :: coord1
          INTEGER,            INTENT(IN   ) :: coord2
          INTEGER,            INTENT(IN   ) :: cbox
          TYPE(MCMCParticle), INTENT(IN   ) :: MCMCParticle_
          LOGICAL,            INTENT(IN   ) :: DoRecord
          INTEGER                           :: info
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(MCMCParticle) :: tmpParticle
          TYPE(MCMCHistoryParticle) :: tmpHistoryParticle

          TYPE(ppm_rc_list), POINTER ::  seed

          REAL(ppm_kind_double) :: t0

          CHARACTER(LEN=ppm_char) :: caller="MCMCInsertFloatingParticleCumulative"
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          tmpParticle=MCMCFloatingParticlesInCell(cbox)%search(coord1,coord2,MCMCParticle_%candlabel)

          !-----------------------------------------------------------------------
          ! there is no floating particle at this position
          !-----------------------------------------------------------------------
          IF (tmpParticle%candlabel.EQ.-1) THEN
             CALL MCMCFloatingParticlesInCell(cbox)%insert(coord1,coord2,MCMCParticle_,info)
          !-----------------------------------------------------------------------
          ! the particle already exists
          !-----------------------------------------------------------------------
          ELSE
             !-----------------------------------------------------------------------
             ! The element did already exist. We add up the proposal and insert
             ! the element again (in order to overwrite).
             !-----------------------------------------------------------------------
             tmpParticle%proposal=tmpParticle%proposal+MCMCParticle_%proposal

             CALL MCMCFloatingParticlesInCell(cbox)%insert(coord1,coord2,tmpParticle,info)
          ENDIF
          or_fail("MCMCFloatingParticlesInCell(cbox)%insert",ppm_error=ppm_error_fatal)

          IF (DoRecord) THEN
             tmpHistoryParticle=MCMCHistoryParticle(MCMCParticle_%candlabel, &
             & 0,MCMCParticle_%proposal,.TRUE.)
             ALLOCATE(seed,STAT=info)
             or_fail_alloc("seed",ppm_error=ppm_error_fatal)
             CALL seed%add(coord1,coord2,tmpHistoryParticle)
             CALL MCMCFloatingParticleInContainerHistory%push(seed,info)
             or_fail("MCMCFloatingParticleInContainerHistory%push",ppm_error=ppm_error_fatal)
          ENDIF

          MCMCTotalNormalizerlocal=MCMCTotalNormalizerlocal+REAL(MCMCParticle_%proposal,ppm_kind_double)

        9999 CONTINUE
          CALL substop(caller,t0,info)
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCInsertFloatingParticleCumulative_2d

        FUNCTION MCMCInsertFloatingParticleCumulative_3d(coord1,coord2,coord3,cbox,MCMCParticle_,DoRecord) RESULT(info)
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_fatal
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_error,ppm_err_alloc,ppm_err_sub_failed

          USE ppm_rc_module_linkedlist, ONLY : ppm_rc_list, &
          &   MCMCFloatingParticleInContainerHistory
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER,            INTENT(IN   ) :: coord1
          INTEGER,            INTENT(IN   ) :: coord2
          INTEGER,            INTENT(IN   ) :: coord3
          INTEGER,            INTENT(IN   ) :: cbox
          TYPE(MCMCParticle), INTENT(IN   ) :: MCMCParticle_
          LOGICAL,            INTENT(IN   ) :: DoRecord
          INTEGER                           :: info
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(MCMCParticle) :: tmpParticle
          TYPE(MCMCHistoryParticle) :: tmpHistoryParticle

          TYPE(ppm_rc_list), POINTER ::  seed

          REAL(ppm_kind_double) :: t0

          CHARACTER(LEN=ppm_char) :: caller="MCMCInsertFloatingParticleCumulative"
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          tmpParticle=MCMCFloatingParticlesInCell(cbox)%search(coord1,coord2,coord3,MCMCParticle_%candlabel)

          !-----------------------------------------------------------------------
          ! there is no floating particle at this position
          !-----------------------------------------------------------------------
          IF (tmpParticle%candlabel.EQ.-1) THEN
             CALL MCMCFloatingParticlesInCell(cbox)%insert(coord1,coord2,coord3,MCMCParticle_,info)
          !-----------------------------------------------------------------------
          ! the particle already exists
          !-----------------------------------------------------------------------
          ELSE
             !-----------------------------------------------------------------------
             ! The element did already exist. We add up the proposal and insert
             ! the element again (in order to overwrite).
             !-----------------------------------------------------------------------
             tmpParticle%proposal=tmpParticle%proposal+MCMCParticle_%proposal

             CALL MCMCFloatingParticlesInCell(cbox)%insert(coord1,coord2,coord3,tmpParticle,info)
          ENDIF
          or_fail("MCMCFloatingParticlesInCell(cbox)%insert",ppm_error=ppm_error_fatal)

          IF (DoRecord) THEN
             tmpHistoryParticle=MCMCHistoryParticle(MCMCParticle_%candlabel, &
             & 0,MCMCParticle_%proposal,.TRUE.)
             ALLOCATE(seed,STAT=info)
             or_fail_alloc("seed",ppm_error=ppm_error_fatal)
             CALL seed%add(coord1,coord2,coord3,tmpHistoryParticle)
             CALL MCMCFloatingParticleInContainerHistory%push(seed,info)
             or_fail("MCMCFloatingParticleInContainerHistory%push",ppm_error=ppm_error_fatal)
          ENDIF

          MCMCTotalNormalizerlocal=MCMCTotalNormalizerlocal+REAL(MCMCParticle_%proposal,ppm_kind_double)

        9999 CONTINUE
          CALL substop(caller,t0,info)
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCInsertFloatingParticleCumulative_3d

        LOGICAL FUNCTION MCMCEraseFloatingParticle_2d(coord1,coord2,cbox,MCMCParticle_,DoRecord)
          !!! Erase a floating particle only if it exists.
          !!! If it does not exist, the return value will be true.
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_fatal
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_error,ppm_err_alloc,ppm_err_sub_failed

          USE ppm_rc_module_global, ONLY : small,ten
          USE ppm_rc_module_linkedlist, ONLY : ppm_rc_list, &
          &   MCMCFloatingParticleInContainerHistory
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER,            INTENT(IN   ) :: coord1
          INTEGER,            INTENT(IN   ) :: coord2
          INTEGER,            INTENT(IN   ) :: cbox
          TYPE(MCMCParticle), INTENT(IN   ) :: MCMCParticle_
          LOGICAL,            INTENT(IN   ) :: DoRecord
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(MCMCParticle) :: tmpParticle
          TYPE(MCMCHistoryParticle) :: tmpHistoryParticle

          TYPE(ppm_rc_list), POINTER ::  seed

          REAL(ppm_kind_double) :: t0

          INTEGER :: info

          CHARACTER(LEN=ppm_char) :: caller="MCMCEraseFloatingParticle"
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          tmpParticle=MCMCFloatingParticlesInCell(cbox)%search(coord1,coord2,MCMCParticle_%candlabel)

          ! there is no floating particle at this position
          IF (tmpParticle%candlabel.EQ.-1) THEN
             MCMCEraseFloatingParticle_2d=.TRUE.
          ! there is a floating particle at this position and we erase that
          ELSE
             MCMCEraseFloatingParticle_2d=.FALSE.

             IF (DoRecord) THEN
                tmpHistoryParticle=MCMCHistoryParticle(MCMCParticle_%candlabel, &
                & 0,MCMCParticle_%proposal,.FALSE.)
                ALLOCATE(seed,STAT=info)
                or_fail_alloc("seed",ppm_error=ppm_error_fatal)
                CALL seed%add(coord1,coord2,tmpHistoryParticle)
                CALL MCMCFloatingParticleInContainerHistory%push(seed,info)
                or_fail("MCMCFloatingParticleInContainerHistory%push",ppm_error=ppm_error_fatal)
             ENDIF

             ! the particle has still some proposal left
             IF ((tmpParticle%proposal-MCMCParticle_%proposal).GT.small*ten) THEN
                MCMCTotalNormalizerlocal=MCMCTotalNormalizerlocal-REAL(MCMCParticle_%proposal,ppm_kind_double)

                tmpParticle%proposal=tmpParticle%proposal-MCMCParticle_%proposal

                CALL MCMCFloatingParticlesInCell(cbox)%insert(coord1,coord2,tmpParticle,info)
                or_fail("MCMCFloatingParticlesInCell(cbox)%insert",ppm_error=ppm_error_fatal)
             ELSE
                ! the particle gets deleted. We only remove the amount from the
                ! normalizer that has been stored in the set.
                MCMCTotalNormalizerlocal=MCMCTotalNormalizerlocal-REAL(tmpParticle%proposal,ppm_kind_double)

                CALL MCMCFloatingParticlesInCell(cbox)%remove(coord1,coord2,MCMCParticle_%candlabel,info,.TRUE.)
                or_fail("MCMCFloatingParticlesInCell(cbox)%remove",ppm_error=ppm_error_fatal)
             ENDIF

             MCMCEraseFloatingParticle_2d=.TRUE.
          ENDIF
        9999 CONTINUE
          CALL substop(caller,t0,info)
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCEraseFloatingParticle_2d

        LOGICAL FUNCTION MCMCEraseFloatingParticle__2d(coord1,coord2,cbox,searchedParticle,MCMCParticle_,DoRecord)
          !!! Erase a floating particle only if it exists. The particle at coords
          !!! has been searched as searchedParticle.
          !!! If it does not exist, the return value will be true.
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_fatal
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_error,ppm_err_alloc,ppm_err_sub_failed

          USE ppm_rc_module_global, ONLY : small,ten
          USE ppm_rc_module_linkedlist, ONLY : ppm_rc_list, &
          &   MCMCFloatingParticleInContainerHistory
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER,            INTENT(IN   ) :: coord1
          INTEGER,            INTENT(IN   ) :: coord2
          INTEGER,            INTENT(IN   ) :: cbox
          TYPE(MCMCParticle), INTENT(INOUT) :: searchedParticle
          TYPE(MCMCParticle), INTENT(IN   ) :: MCMCParticle_
          LOGICAL,            INTENT(IN   ) :: DoRecord
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(MCMCHistoryParticle) :: tmpHistoryParticle

          TYPE(ppm_rc_list), POINTER ::  seed

          REAL(ppm_kind_double) :: t0

          INTEGER :: info

          CHARACTER(LEN=ppm_char) :: caller="MCMCEraseFloatingParticle"
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          ! there is no floating particle at this position
          IF (searchedParticle%candlabel.EQ.-1) THEN
             MCMCEraseFloatingParticle__2d=.TRUE.
          ELSE
             IF (DoRecord) THEN
                tmpHistoryParticle=MCMCHistoryParticle(MCMCParticle_%candlabel, &
                & 0,MCMCParticle_%proposal,.FALSE.)
                ALLOCATE(seed,STAT=info)
                or_fail_alloc("seed",ppm_error=ppm_error_fatal)
                CALL seed%add(coord1,coord2,tmpHistoryParticle)
                CALL MCMCFloatingParticleInContainerHistory%push(seed,info)
                or_fail("MCMCFloatingParticleInContainerHistory%push",ppm_error=ppm_error_fatal)
             ENDIF

             ! the particle has still some proposal left
             IF ((searchedParticle%proposal-MCMCParticle_%proposal).GT.small*ten) THEN
                MCMCTotalNormalizerlocal=MCMCTotalNormalizerlocal-REAL(MCMCParticle_%proposal,ppm_kind_double)

                searchedParticle%proposal=searchedParticle%proposal-MCMCParticle_%proposal

                CALL MCMCFloatingParticlesInCell(cbox)%insert(coord1,coord2,searchedParticle,info)
                or_fail("MCMCFloatingParticlesInCell(cbox)%insert",ppm_error=ppm_error_fatal)
             ELSE
                ! the particle gets deleted. We only remove the amount from the
                ! normalizer that has been stored in the set.
                MCMCTotalNormalizerlocal=MCMCTotalNormalizerlocal-REAL(searchedParticle%proposal,ppm_kind_double)

                CALL MCMCFloatingParticlesInCell(cbox)%remove(coord1,coord2,MCMCParticle_%candlabel,info,.TRUE.)
                or_fail("MCMCFloatingParticlesInCell(cbox)%remove",ppm_error=ppm_error_fatal)
             ENDIF
             MCMCEraseFloatingParticle__2d=.TRUE.
          ENDIF
        9999 CONTINUE
          CALL substop(caller,t0,info)
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCEraseFloatingParticle__2d

        LOGICAL FUNCTION MCMCEraseFloatingParticle_3d(coord1,coord2,coord3,cbox,MCMCParticle_,DoRecord)
          !!! Erase a floating particle only if it exists.
          !!! If it does not exist, the return value will be true.
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_fatal
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_error,ppm_err_alloc,ppm_err_sub_failed

          USE ppm_rc_module_global, ONLY : small,ten
          USE ppm_rc_module_linkedlist, ONLY : ppm_rc_list, &
          &   MCMCFloatingParticleInContainerHistory
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER,            INTENT(IN   ) :: coord1
          INTEGER,            INTENT(IN   ) :: coord2
          INTEGER,            INTENT(IN   ) :: coord3
          INTEGER,            INTENT(IN   ) :: cbox
          TYPE(MCMCParticle), INTENT(IN   ) :: MCMCParticle_
          LOGICAL,            INTENT(IN   ) :: DoRecord
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(MCMCParticle) :: tmpParticle
          TYPE(MCMCHistoryParticle) :: tmpHistoryParticle

          TYPE(ppm_rc_list), POINTER ::  seed

          REAL(ppm_kind_double) :: t0

          INTEGER :: info

          CHARACTER(LEN=ppm_char) :: caller="MCMCEraseFloatingParticle"
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          tmpParticle=MCMCFloatingParticlesInCell(cbox)%search(coord1,coord2,coord3,MCMCParticle_%candlabel)

          ! there is no floating particle at this position
          IF (tmpParticle%candlabel.EQ.-1) THEN
             MCMCEraseFloatingParticle_3d=.TRUE.
          ELSE
             IF (DoRecord) THEN
                tmpHistoryParticle=MCMCHistoryParticle(MCMCParticle_%candlabel, &
                & 0,MCMCParticle_%proposal,.FALSE.)
                ALLOCATE(seed,STAT=info)
                or_fail_alloc("seed",ppm_error=ppm_error_fatal)
                CALL seed%add(coord1,coord2,coord3,tmpHistoryParticle)
                CALL MCMCFloatingParticleInContainerHistory%push(seed,info)
                or_fail("MCMCFloatingParticleInContainerHistory%push",ppm_error=ppm_error_fatal)
             ENDIF

             ! the particle has still some proposal left
             IF ((tmpParticle%proposal-MCMCParticle_%proposal).GT.small*ten) THEN
                MCMCTotalNormalizerlocal=MCMCTotalNormalizerlocal-REAL(MCMCParticle_%proposal,ppm_kind_double)

                tmpParticle%proposal=tmpParticle%proposal-MCMCParticle_%proposal

                CALL MCMCFloatingParticlesInCell(cbox)%insert(coord1,coord2,coord3,tmpParticle,info)
                or_fail("MCMCFloatingParticlesInCell(cbox)%insert",ppm_error=ppm_error_fatal)
             ELSE
                ! the particle gets deleted. We only remove the amount from the
                ! normalizer that has been stored in the set.
                MCMCTotalNormalizerlocal=MCMCTotalNormalizerlocal-REAL(tmpParticle%proposal,ppm_kind_double)

                CALL MCMCFloatingParticlesInCell(cbox)%remove(coord1,coord2,coord3,MCMCParticle_%candlabel,info,.TRUE.)
                or_fail("MCMCFloatingParticlesInCell(cbox)%remove",ppm_error=ppm_error_fatal)
             ENDIF
             MCMCEraseFloatingParticle_3d=.TRUE.
          ENDIF
        9999 CONTINUE
          CALL substop(caller,t0,info)
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCEraseFloatingParticle_3d

        LOGICAL FUNCTION MCMCEraseFloatingParticle__3d(coord1,coord2,coord3,cbox,searchedParticle,MCMCParticle_,DoRecord)
          !!! Erase a floating particle only if it exists. The particle at coords
          !!! has been searched as searchedParticle.
          !!! If it does not exist, the return value will be true.
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_fatal
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_error,ppm_err_alloc,ppm_err_sub_failed

          USE ppm_rc_module_global, ONLY : small,ten
          USE ppm_rc_module_linkedlist, ONLY : ppm_rc_list, &
          &   MCMCFloatingParticleInContainerHistory
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER,            INTENT(IN   ) :: coord1
          INTEGER,            INTENT(IN   ) :: coord2
          INTEGER,            INTENT(IN   ) :: coord3
          INTEGER,            INTENT(IN   ) :: cbox
          TYPE(MCMCParticle), INTENT(INOUT) :: searchedParticle
          TYPE(MCMCParticle), INTENT(IN   ) :: MCMCParticle_
          LOGICAL,            INTENT(IN   ) :: DoRecord
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(MCMCHistoryParticle) :: tmpHistoryParticle

          TYPE(ppm_rc_list), POINTER ::  seed

          REAL(ppm_kind_double) :: t0

          INTEGER :: info

          CHARACTER(LEN=ppm_char) :: caller="MCMCEraseFloatingParticle"
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          ! there is no floating particle at this position
          IF (searchedParticle%candlabel.EQ.-1) THEN
             MCMCEraseFloatingParticle__3d=.FALSE.
          ELSE
             IF (DoRecord) THEN
                tmpHistoryParticle=MCMCHistoryParticle(MCMCParticle_%candlabel, &
                & 0,MCMCParticle_%proposal,.FALSE.)
                ALLOCATE(seed,STAT=info)
                or_fail_alloc("seed",ppm_error=ppm_error_fatal)
                CALL seed%add(coord1,coord2,coord3,tmpHistoryParticle)
                CALL MCMCFloatingParticleInContainerHistory%push(seed,info)
                or_fail("MCMCFloatingParticleInContainerHistory%push",ppm_error=ppm_error_fatal)
             ENDIF

             ! the particle has still some proposal left
             IF ((searchedParticle%proposal-MCMCParticle_%proposal).GT.small*ten) THEN
                MCMCTotalNormalizerlocal=MCMCTotalNormalizerlocal-REAL(MCMCParticle_%proposal,ppm_kind_double)

                searchedParticle%proposal=searchedParticle%proposal-MCMCParticle_%proposal

                CALL MCMCFloatingParticlesInCell(cbox)%insert(coord1,coord2,coord3,searchedParticle,info)
                or_fail("MCMCFloatingParticlesInCell(cbox)%insert",ppm_error=ppm_error_fatal)
             ELSE
                ! the particle gets deleted. We only remove the amount from the
                ! normalizer that has been stored in the set.
                MCMCTotalNormalizerlocal=MCMCTotalNormalizerlocal-REAL(searchedParticle%proposal,ppm_kind_double)

                CALL MCMCFloatingParticlesInCell(cbox)%remove(coord1,coord2,coord3,MCMCParticle_%candlabel,info,.TRUE.)
                or_fail("MCMCFloatingParticlesInCell(cbox)%remove",ppm_error=ppm_error_fatal)
             ENDIF
             MCMCEraseFloatingParticle__3d=.TRUE.
          ENDIF
        9999 CONTINUE
          CALL substop(caller,t0,info)
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCEraseFloatingParticle__3d

        FUNCTION ppm_rc_MCMCinitmove() RESULT(info)
          USE ppm_module_data, ONLY : ppm_char,ppm_kind_double,ppm_kind_int64, &
          &   ppm_nproc,ppm_error_error
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_error,ppm_err_sub_failed, &
          &   ppm_err_alloc

          USE ppm_rc_module_global, ONLY : mesh,ghostsize_run,ghost_size
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER :: info
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          REAL(ppm_kind_double) :: t0

          CHARACTER(LEN=ppm_char) :: caller='ppm_rc_MCMCinitmove'

          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          IF (ppm_nproc.GT.1) THEN
             ALLOCATE(MCMColdghost(ghost_size),STAT=info)
             or_fail_alloc("MCMColdghost")

             !-------------------------------------------------------------------------
             ! Ghost get
             !-------------------------------------------------------------------------
             CALL mesh%map_ghost_get(info,ghostsize=ghostsize_run)
             or_fail("mesh%map_ghost_get")
          ENDIF
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        END FUNCTION ppm_rc_MCMCinitmove

        FUNCTION MCMCAddAndRemoveParticlesWhenMove_2d(labels_,coord,Nm,MCMCParticle_) RESULT(info)
          !!! Updates the DS (label image and the regular particle containers) when
          !!! applying a particle. Note that floating particles will not be updated!
          !!! The method only ensures that L and the regular  particles are correct.
          !!! The method expects the label image NOT to be updated already. Else the
          !!! operations performed will be wrong.
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_error
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_err_alloc,ppm_error,ppm_err_sub_failed

          USE ppm_rc_module_global, ONLY : FORBIDDEN
          USE ppm_rc_module_topologicalnumber, ONLY : FG_ConnectivityType, &
          &   BG_ConnectivityType,IsEnclosedByLabel_BGConnectivity
          USE ppm_rc_module_linkedlist, ONLY : ppm_rc_list,MCMCLabelImageHistory
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER, CONTIGUOUS, DIMENSION(:,:), POINTER       :: labels_
          INTEGER,             DIMENSION(:),   INTENT(IN   ) :: coord
          INTEGER,             DIMENSION(:),   INTENT(IN   ) :: Nm
          TYPE(MCMCParticle),                  INTENT(IN   ) :: MCMCParticle_
          INTEGER                                            :: info
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(MCMCParticle) :: tmpParticle
          TYPE(MCMCParticle) :: ReverseParticle

          TYPE(ppm_rc_list), POINTER ::  seed

          REAL(ppm_kind_double) :: t0
          REAL(MK)              :: proposal

          INTEGER, DIMENSION(:,:), POINTER :: labelstmp
          INTEGER, DIMENSION(:,:), POINTER :: tmplabels

          INTEGER, DIMENSION(2)            :: ll,ld,dd
          INTEGER                          :: LabelFrom
          INTEGER                          :: LabelTo
          INTEGER                          :: StoreLabel
          INTEGER                          :: LabelNgh
          INTEGER                          :: i,j
          INTEGER                          :: nx,ny
          INTEGER                          :: cbox

          LOGICAL :: ReverseParticleIsFloating
          LOGICAL :: Replaced
          LOGICAL :: HasOtherMother

          CHARACTER(LEN=ppm_char) :: caller='MCMCAddAndRemoveParticlesWhenMove'
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          !-------------------------------------------------------------------------
          ! Initilize a neighborhood iterator for fast access to the label image
          ! in the vicinity of the particle of interest.
          !-------------------------------------------------------------------------
          labelstmp => labels_(coord(1)-2:coord(1)+2,coord(2)-2:coord(2)+2)

          StoreLabel=labelstmp(3,3)

          !-------------------------------------------------------------------------
          ! Initialize locals from the particle
          !-------------------------------------------------------------------------
          LabelFrom=ABS(StoreLabel)
          LabelTo=MCMCParticle_%candlabel

          !-------------------------------------------------------------------------
          ! We remove the particle and insert the reverse particle:
          ! Here we cannot decide if the reverse particle is indeed
          ! a floating particle (this is because we apply B first, then A).
          ! The solution is that this method only works on the regular particle
          ! set. Floating particles will be detected and treated outside of
          ! this method. Here we omit the potential insertion of floating
          ! particles.
          ! In order to (maybe) replace the particle with its reverse particle we:
          ! Simulate the move, calculate the proposal for the backward particle,
          ! create the backward particle, check if the backward particle is
          ! floating, and finally restore the label image:
          !-------------------------------------------------------------------------

          labelstmp(3,3)=-LabelTo

          tmplabels => labelstmp(2:4,2:4)

          proposal=MCMCproposal(tmplabels)

          ReverseParticle=MCMCParticle(LabelFrom,proposal)
          !                                                        (tmplabels,LabelIn,CandidateLabel)
          !                                                         tmplabels,       ,ReverseParticle%candlabel
          ReverseParticleIsFloating=MCMCParticleHasFloatingProperty(tmplabels,LabelTo,LabelFrom)

          labelstmp(3,3)=StoreLabel


          nx=INT(REAL(coord(1),MK)/MCMCcellsize(1))
          IF (nx.EQ.MCMCcellNm(1)) nx=nx-1

          ny=INT(REAL(coord(2),MK)/MCMCcellsize(2))
          IF (ny.EQ.MCMCcellNm(2)) ny=ny-1

          cbox=1+nx+ny*MCMCcellNm(1)

          !-------------------------------------------------------------------------
          ! Insert the reverse particle (if its not floating, see above):
          !-------------------------------------------------------------------------
          IF (.NOT.ReverseParticleIsFloating) THEN
             info=MCMCInsertCandidatesToContainers(coord(1),coord(2),cbox, &
             &    ReverseParticle,LabelTo,.TRUE.,Replaced)
             or_fail("MCMCInsertCandidatesToContainers")
          ENDIF

          !-------------------------------------------------------------------------
          ! erase the currently applied particle (if its floating this will not hurt
          ! to call erase here).
          !-------------------------------------------------------------------------
          info=MCMCEraseCandidatesFromContainers(coord(1),coord(2),cbox, &
          &    MCMCParticle_,LabelFrom,.TRUE.)
          or_fail("MCMCEraseCandidatesFromContainers")

          !-------------------------------------------------------------------------
          ! What particle would be added or removed to the contour lists
          ! of the currently shrinking region:
          ! (compare to also AddNeighborsAtRemove(..))
          !-------------------------------------------------------------------------
          IF (LabelFrom.NE.0) THEN
             !-------------------------------------------------------------------------
             ! A FG region is shrinking (Moving particle is either a parent or child)
             !-------------------------------------------------------------------------
             DO i=1,BG_ConnectivityType%NumberOfNeighbors
                ld=coord+BG_ConnectivityType%NeighborsPoints(:,i)
                IF (ANY(ld.LT.1.OR.ld.GT.Nm)) CYCLE
                ll=BG_ConnectivityType%NeighborsPoints(:,i)+3
                LabelNgh=labelstmp(ll(1),ll(2))
                !-------------------------------------------------------------------------
                ! Check if new points enter the contour (internal point becomes a parent):
                ! Internal points of this region are positive:
                !-------------------------------------------------------------------------
                IF (LabelNgh.EQ.LabelFrom) THEN
                   tmplabels => labelstmp(ll(1)-1:ll(1)+1,ll(2)-1:ll(2)+1)
                   proposal=MCMCproposal(tmplabels)

                   !-------------------------------------------------------------------------
                   ! Internal point becomes a parent
                   !-------------------------------------------------------------------------
                   tmpParticle=MCMCParticle(0,proposal)

                   nx=INT(REAL(ld(1),MK)/MCMCcellsize(1))
                   IF (nx.EQ.MCMCcellNm(1)) nx=nx-1

                   ny=INT(REAL(ld(2),MK)/MCMCcellsize(2))
                   IF (ny.EQ.MCMCcellNm(2)) ny=ny-1

                   cbox=1+nx+ny*MCMCcellNm(1)

                   !-------------------------------------------------------------------------
                   ! Insert new particle in container
                   !-------------------------------------------------------------------------
                   info=MCMCInsertCandidatesToContainers(ld(1),ld(2),cbox, &
                   &    tmpParticle,LabelFrom,.TRUE.,Replaced)
                   or_fail("MCMCInsertCandidatesToContainers")

                   IF (.NOT.Replaced) THEN
                      labelstmp(ll(1),ll(2))=-LabelFrom

                      ALLOCATE(seed,STAT=info)
                      or_fail_alloc("seed")
                      CALL seed%add(ld(1),ld(2),LabelFrom)
                      CALL MCMCLabelImageHistory%push(seed,info)
                      or_fail("MCMCLabelImageHistory%push")
                   ENDIF
                ENDIF !LabelNgh.GT.0.AND.LabelNgh.EQ.LabelFrom
             ENDDO !i=1,BG_ConnectivityType%NumberOfNeighbors

             DO i=1,FG_ConnectivityType%NumberOfNeighbors
                ld=coord+FG_ConnectivityType%NeighborsPoints(:,i)
                IF (ANY(ld.LT.1.OR.ld.GT.Nm)) CYCLE
                ll=FG_ConnectivityType%NeighborsPoints(:,i)+3
                LabelNgh=ABS(labelstmp(ll(1),ll(2)))

                !-------------------------------------------------------------------------
                ! check if there are FG-neighbors with no other mother of the
                ! same label --> orphan.
                ! we first 'simulate' the move:
                !-------------------------------------------------------------------------
                IF (LabelNgh.NE.LabelFrom) THEN
                   StoreLabel=labelstmp(3,3)

                   labelstmp(3,3)=-LabelTo

                   !-------------------------------------------------------------------------
                   ! Check if this neighbor has other mothers from this label:
                   !-------------------------------------------------------------------------
                   HasOtherMother=.FALSE.
                   DO j=1,FG_ConnectivityType%NumberOfNeighbors
                      dd=ll+FG_ConnectivityType%NeighborsPoints(:,i)
                      IF (ABS(labelstmp(dd(1),dd(2))).EQ.LabelFrom) THEN
                         HasOtherMother=.TRUE.
                         EXIT
                      ENDIF
                   ENDDO
                   IF (.NOT.HasOtherMother) THEN
                      !-------------------------------------------------------------------------
                      ! The orphan has label equal to what we read from the
                      ! label image and has a candidate label of the
                      ! currently shrinking region.
                      !
                      ! Note:
                      ! When erasing a particle we do not need to compute the proposal again
                      ! Only coords, and particle candidate label are important for the key in
                      ! the hash tabel
                      !-------------------------------------------------------------------------
                      tmpParticle=MCMCParticle(LabelFrom,zero)

                      nx=INT(REAL(ld(1),MK)/MCMCcellsize(1))
                      IF (nx.EQ.MCMCcellNm(1)) nx=nx-1

                      ny=INT(REAL(ld(2),MK)/MCMCcellsize(2))
                      IF (ny.EQ.MCMCcellNm(2)) ny=ny-1

                      cbox=1+nx+ny*MCMCcellNm(1)

                      info=MCMCEraseCandidatesFromContainers(ld(1),ld(2),cbox, &
                      &    tmpParticle,LabelNgh,.TRUE.)
                      or_fail("MCMCEraseCandidatesFromContainers")
                   ENDIF !.NOT.HasOtherMother

                   labelstmp(3,3)=StoreLabel
                ENDIF !LabelNgh.NE.LabelFrom
             ENDDO !i=1,FG_ConnectivityType%NumberOfNeighbors
          ENDIF ! LabelFrom.NE.0

          !-------------------------------------------------------------------------
          ! Growing: figure out the changes of candidates for the expanding region
          !-------------------------------------------------------------------------
          IF (LabelTo.NE.0) THEN
             !-------------------------------------------------------------------------
             ! we are growing
             !
             ! Neighbors: Figure out what (neighboring)mother points are going
             ! to be interior points:
             !
             ! simulate the move
             !-------------------------------------------------------------------------
             StoreLabel=labelstmp(3,3)

             labelstmp(3,3)=-LabelTo

             DO i=1,BG_ConnectivityType%NumberOfNeighbors
                ld=coord+BG_ConnectivityType%NeighborsPoints(:,i)
                IF (ANY(ld.LT.1.OR.ld.GT.Nm)) CYCLE
                ll=BG_ConnectivityType%NeighborsPoints(:,i)+3
                LabelNgh=labelstmp(ll(1),ll(2))

                tmplabels => labelstmp(ll(1)-1:ll(1)+1,ll(2)-1:ll(2)+1)

                !-------------------------------------------------------------------------
                ! TODO: the isEnclosedByLabel method could use the above iterator
                !-------------------------------------------------------------------------
                IF (LabelNgh.EQ.-LabelTo.AND.IsEnclosedByLabel_BGConnectivity(tmplabels,LabelTo)) THEN
                   !-------------------------------------------------------------------------
                   ! Remove the parent that got enclosed; it had a the label
                   ! of the currently expanding region and a candidate label of 0.
                   !-------------------------------------------------------------------------
                   tmpParticle=MCMCParticle(0,zero)

                   nx=INT(REAL(ld(1),MK)/MCMCcellsize(1))
                   IF (nx.EQ.MCMCcellNm(1)) nx=nx-1

                   ny=INT(REAL(ld(2),MK)/MCMCcellsize(2))
                   IF (ny.EQ.MCMCcellNm(2)) ny=ny-1

                   cbox=1+nx+ny*MCMCcellNm(1)

                   info=MCMCEraseCandidatesFromContainers(ld(1),ld(2),cbox, &
                   &    tmpParticle,LabelTo,.TRUE.)
                   or_fail("MCMCEraseCandidatesFromContainers")

                   !-------------------------------------------------------------------------
                   ! update the label image
                   !-------------------------------------------------------------------------
                   labelstmp(ll(1),ll(2))=LabelTo

                   ALLOCATE(seed,STAT=info)
                   or_fail_alloc("seed")
                   CALL seed%add(ld(1),ld(2),LabelNgh)
                   CALL MCMCLabelImageHistory%push(seed,info)
                   or_fail("MCMCLabelImageHistory%push")
                ENDIF
             ENDDO !i=1,BG_ConnectivityType%NumberOfNeighbors

             labelstmp(3,3)=StoreLabel

             !-------------------------------------------------------------------------
             ! Figure out if a point renders to a candidate. These are
             ! all the FG-neighbors with a different label that are not yet
             ! candidates of the currently expanding region.
             !-------------------------------------------------------------------------
             DO i=1,FG_ConnectivityType%NumberOfNeighbors
                ld=coord+FG_ConnectivityType%NeighborsPoints(:,i)
                IF (ANY(ld.LT.1.OR.ld.GT.Nm)) CYCLE
                ll=FG_ConnectivityType%NeighborsPoints(:,i)+3
                LabelNgh=labelstmp(ll(1),ll(2))

                IF (ABS(LabelNgh).NE.LabelTo.AND.LabelNgh.NE.FORBIDDEN) THEN
                   !-------------------------------------------------------------------------
                   ! check if there is no other mother (hence the particle is
                   ! not in the container yet). This we could do by checking the
                   ! neighborhood of the label image or by checking the
                   ! cointainers.
                   ! Here: we check the (not yet updated!) label image.
                   !-------------------------------------------------------------------------
                   HasOtherMother=.FALSE.
                   DO j=1,FG_ConnectivityType%NumberOfNeighbors
                      dd=ll+FG_ConnectivityType%NeighborsPoints(:,i)
                      IF (ABS(labelstmp(dd(1),dd(2))).EQ.LabelTo) THEN
                         HasOtherMother=.TRUE.
                         EXIT
                      ENDIF
                   ENDDO
                   IF (.NOT.HasOtherMother) THEN
                      !-------------------------------------------------------------------------
                      ! This is a new child. It's current label we have to read
                      ! from the label image, the candidate label is the label of
                      ! the currently expanding region.
                      !-------------------------------------------------------------------------
                      tmplabels => labelstmp(ll(1)-1:ll(1)+1,ll(2)-1:ll(2)+1)
                      proposal=MCMCproposal(tmplabels)

                      tmpParticle=MCMCParticle(LabelTo,proposal)

                      nx=INT(REAL(ld(1),MK)/MCMCcellsize(1))
                      IF (nx.EQ.MCMCcellNm(1)) nx=nx-1

                      ny=INT(REAL(ld(2),MK)/MCMCcellsize(2))
                      IF (ny.EQ.MCMCcellNm(2)) ny=ny-1

                      cbox=1+nx+ny*MCMCcellNm(1)

                      info=MCMCInsertCandidatesToContainers(ld(1),ld(2),cbox, &
                      &    tmpParticle,ABS(LabelNgh),.TRUE.,Replaced)
                      or_fail("MCMCInsertCandidatesToContainers")

                      IF (.NOT.Replaced) THEN
                         labelstmp(ll(1),ll(2))=SIGN(LabelNgh,-1)

                         ALLOCATE(seed,STAT=info)
                         or_fail_alloc("seed")
                         CALL seed%add(ld(1),ld(2),LabelNgh)
                         CALL MCMCLabelImageHistory%push(seed,info)
                         or_fail("MCMCLabelImageHistory%push")
                      ENDIF
                   ENDIF !.NOT.HasOtherMother
                ENDIF !LabelNgh.NE.LabelTo.AND.LabelNgh.NE.FORBIDDEN
             ENDDO !i=1,FG_ConnectivityType%NumberOfNeighbors
          ENDIF ! LabelTo.NE.0
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        END FUNCTION MCMCAddAndRemoveParticlesWhenMove_2d

        FUNCTION MCMCAddAndRemoveParticlesWhenMove_3d(labels_,coord,Nm,MCMCParticle_) RESULT(info)
          ! Updates the DS (label image and the regular particle containers) when
          ! applying a particle. Note that floating particles will not be updated!
          ! The method only ensures that L and the regular  particles are correct.
          ! The method expects the label image NOT to be updated already. Else the
          ! operations performed will be wrong.
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_error
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_err_alloc,ppm_error,ppm_err_sub_failed

          USE ppm_rc_module_global, ONLY : FORBIDDEN
          USE ppm_rc_module_topologicalnumber, ONLY : FG_ConnectivityType, &
          &   BG_ConnectivityType,IsEnclosedByLabel_BGConnectivity
          USE ppm_rc_module_linkedlist, ONLY : ppm_rc_list,MCMCLabelImageHistory
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: labels_
          INTEGER,             DIMENSION(:),     INTENT(IN   ) :: coord
          INTEGER,             DIMENSION(:),     INTENT(IN   ) :: Nm
          TYPE(MCMCParticle),                    INTENT(IN   ) :: MCMCParticle_
          INTEGER                                              :: info
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(MCMCParticle) :: tmpParticle
          TYPE(MCMCParticle) :: ReverseParticle

          TYPE(ppm_rc_list), POINTER ::  seed

          REAL(ppm_kind_double) :: t0
          REAL(MK)              :: proposal

          INTEGER, DIMENSION(:,:,:), POINTER :: labelstmp
          INTEGER, DIMENSION(:,:,:), POINTER :: tmplabels

          INTEGER, DIMENSION(3)              :: ll,ld,dd
          INTEGER                            :: LabelFrom
          INTEGER                            :: LabelTo
          INTEGER                            :: StoreLabel
          INTEGER                            :: LabelNgh
          INTEGER                            :: i,j
          INTEGER                            :: nx,ny,nz
          INTEGER                            :: cbox


          LOGICAL :: ReverseParticleIsFloating
          LOGICAL :: Replaced
          LOGICAL :: HasOtherMother

          CHARACTER(LEN=ppm_char) :: caller='MCMCAddAndRemoveParticlesWhenMove'
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          !-----------------------------------------------------------------------
          ! Initilize an neighborhooditerator for fast access to the label image
          ! in the vicinity of the particle of interest.
          !-----------------------------------------------------------------------
          labelstmp => labels_(coord(1)-2:coord(1)+2,coord(2)-2:coord(2)+2,coord(3)-2:coord(3)+2)

          StoreLabel=labelstmp(3,3,3)

          !-----------------------------------------------------------------------
          ! Initialize locals from the particle
          !-----------------------------------------------------------------------
          LabelFrom=ABS(StoreLabel)
          LabelTo=MCMCParticle_%candlabel

          !-----------------------------------------------------------------------
          ! We remove the particle and insert the reverse particle:
          ! Here we cannot decide if the reverse particle is indeed
          ! a floating particle (this is because we apply B first, then A).
          ! The solution is that this method only works on the regular particle
          ! set. Floating particles will be detected and treated outside of
          ! this method. Here we omit the potential insertion of floating
          ! particles.
          ! In order to (maybe) replace the particle with its reverse particle we:
          ! Simulate the move, calculate the proposal for the backward particle,
          ! create the backward particle, check if the backward particle is
          ! floating, and finally restore the label image:
          !-----------------------------------------------------------------------

          labelstmp(3,3,3)=-LabelTo

          tmplabels => labelstmp(2:4,2:4,2:4)
          proposal=MCMCproposal(tmplabels)

          ReverseParticle=MCMCParticle(LabelFrom,proposal)
          !                                                        (tmplabels,LabelIn,CandidateLabel)
          !                                                         tmplabels,       ,ReverseParticle%candlabel
          ReverseParticleIsFloating=MCMCParticleHasFloatingProperty(tmplabels,LabelTo,LabelFrom)

          labelstmp(3,3,3)=StoreLabel


          nx=INT(REAL(coord(1),MK)/MCMCcellsize(1))
          IF (nx.EQ.MCMCcellNm(1)) nx=nx-1

          ny=INT(REAL(coord(2),MK)/MCMCcellsize(2))
          IF (ny.EQ.MCMCcellNm(2)) ny=ny-1

          nz=INT(REAL(coord(3),MK)/MCMCcellsize(3))
          IF (nz.EQ.MCMCcellNm(3)) nz=nz-1

          cbox=1+nx+MCMCcellNm(1)*(ny+nz*MCMCcellNm(2))

          !-----------------------------------------------------------------------
          ! Insert the reverse particle (if its not floating, see above):
          !-----------------------------------------------------------------------
          IF (.NOT.ReverseParticleIsFloating) THEN
             info=MCMCInsertCandidatesToContainers(coord(1),coord(2),coord(3),cbox, &
             &    ReverseParticle,LabelTo,.TRUE.,Replaced)
             or_fail("MCMCInsertCandidatesToContainers")
          ENDIF

          !-----------------------------------------------------------------------
          ! erase the currently applied particle (if its floating this will not hurt
          ! to call erase here).
          !-----------------------------------------------------------------------
          info=MCMCEraseCandidatesFromContainers(coord(1),coord(2),coord(3),cbox, &
          &    MCMCParticle_,LabelFrom,.TRUE.)
          or_fail("MCMCEraseCandidatesFromContainers")

          !-----------------------------------------------------------------------
          ! What particle would be added or removed to the contour lists
          ! of the currently shrinking region:
          ! (compare to also AddNeighborsAtRemove(..))
          !-----------------------------------------------------------------------
          IF (LabelFrom.NE.0) THEN
             !-----------------------------------------------------------------------
             ! A FG region is shrinking
             !-----------------------------------------------------------------------
             DO i=1,BG_ConnectivityType%NumberOfNeighbors
                ld=coord+BG_ConnectivityType%NeighborsPoints(:,i)
                IF (ANY(ld.LT.1.OR.ld.GT.Nm)) CYCLE
                ll=BG_ConnectivityType%NeighborsPoints(:,i)+3
                LabelNgh=labelstmp(ll(1),ll(2),ll(3))
                !-----------------------------------------------------------------------
                ! Check if new points enter the contour (internal point becomes a parent):
                ! Internal points of this region are positive:
                !-----------------------------------------------------------------------
                IF (LabelNgh.EQ.LabelFrom) THEN
                   tmplabels => labelstmp(ll(1)-1:ll(1)+1,ll(2)-1:ll(2)+1,ll(3)-1:ll(3)+1)
                   proposal=MCMCproposal(tmplabels)

                   tmpParticle=MCMCParticle(0,proposal)

                   nx=INT(REAL(ld(1),MK)/MCMCcellsize(1))
                   IF (nx.EQ.MCMCcellNm(1)) nx=nx-1

                   ny=INT(REAL(ld(2),MK)/MCMCcellsize(2))
                   IF (ny.EQ.MCMCcellNm(2)) ny=ny-1

                   nz=INT(REAL(ld(3),MK)/MCMCcellsize(3))
                   IF (nz.EQ.MCMCcellNm(3)) nz=nz-1

                   cbox=1+nx+MCMCcellNm(1)*(ny+nz*MCMCcellNm(2))

                   info=MCMCInsertCandidatesToContainers(ld(1),ld(2),ld(3),cbox, &
                   &    tmpParticle,LabelFrom,.TRUE.,Replaced)
                   or_fail("MCMCInsertCandidatesToContainers")

                   IF (.NOT.Replaced) THEN
                      labelstmp(ll(1),ll(2),ll(3))=-LabelFrom

                      ALLOCATE(seed,STAT=info)
                      or_fail_alloc("seed")
                      CALL seed%add(ld(1),ld(2),ld(3),LabelFrom)
                      CALL MCMCLabelImageHistory%push(seed,info)
                      or_fail("MCMCLabelImageHistory%push")
                   ENDIF
                ENDIF !LabelNgh.GT.0.AND.LabelNgh.EQ.LabelFrom
             ENDDO !i=1,BG_ConnectivityType%NumberOfNeighbors

             DO i=1,FG_ConnectivityType%NumberOfNeighbors
                ld=coord+FG_ConnectivityType%NeighborsPoints(:,i)
                IF (ANY(ld.LT.1.OR.ld.GT.Nm)) CYCLE
                ll=FG_ConnectivityType%NeighborsPoints(:,i)+3
                LabelNgh=ABS(labelstmp(ll(1),ll(2),ll(3)))

                !-----------------------------------------------------------------------
                ! check if there are FG-neighbors with no other mother of the
                ! same label --> orphan.
                ! we first 'simulate' the move:
                !-----------------------------------------------------------------------
                IF (LabelNgh.NE.LabelFrom) THEN
                   StoreLabel=labelstmp(3,3,3)

                   labelstmp(3,3,3)=-LabelTo

                   !-----------------------------------------------------------------------
                   ! check if this neighbor has other mothers from this label:
                   !-----------------------------------------------------------------------
                   HasOtherMother=.FALSE.
                   DO j=1,FG_ConnectivityType%NumberOfNeighbors
                      dd=ll+FG_ConnectivityType%NeighborsPoints(:,i)
                      IF (ABS(labelstmp(dd(1),dd(2),dd(3))).EQ.LabelFrom) THEN
                         HasOtherMother=.TRUE.
                         EXIT
                      ENDIF
                   ENDDO
                   IF (.NOT.HasOtherMother) THEN
                      !-----------------------------------------------------------------------
                      ! The orphan has label equal to what we read from the
                      ! label image and has a candidate label of the
                      ! currently shrinking region.
                      !-----------------------------------------------------------------------
                      tmpParticle=MCMCParticle(LabelFrom,zero)

                      nx=INT(REAL(ld(1),MK)/MCMCcellsize(1))
                      IF (nx.EQ.MCMCcellNm(1)) nx=nx-1

                      ny=INT(REAL(ld(2),MK)/MCMCcellsize(2))
                      IF (ny.EQ.MCMCcellNm(2)) ny=ny-1

                      nz=INT(REAL(ld(3),MK)/MCMCcellsize(3))
                      IF (nz.EQ.MCMCcellNm(3)) nz=nz-1

                      cbox=1+nx+MCMCcellNm(1)*(ny+nz*MCMCcellNm(2))

                      info=MCMCEraseCandidatesFromContainers(ld(1),ld(2),ld(3),cbox, &
                      &    tmpParticle,LabelNgh,.TRUE.)
                      or_fail("MCMCEraseCandidatesFromContainers")
                   ENDIF !.NOT.HasOtherMother

                   labelstmp(3,3,3)=StoreLabel
                ENDIF !LabelNgh.NE.LabelFrom
             ENDDO !i=1,FG_ConnectivityType%NumberOfNeighbors
          ENDIF ! LabelFrom.NE.0

          !-----------------------------------------------------------------------
          ! Growing: figure out the changes of candidates for the expanding region
          !-----------------------------------------------------------------------
          IF (LabelTo.NE.0) THEN
             !-----------------------------------------------------------------------
             ! we are growing
             !
             ! Neighbors: Figure out what (neighboring)mother points are going
             ! to be interior points:
             !
             ! simulate the move
             !-----------------------------------------------------------------------
             StoreLabel=labelstmp(3,3,3)

             labelstmp(3,3,3)=-LabelTo

             DO i=1,BG_ConnectivityType%NumberOfNeighbors
                ld=coord+BG_ConnectivityType%NeighborsPoints(:,i)
                IF (ANY(ld.LT.1.OR.ld.GT.Nm)) CYCLE
                ll=BG_ConnectivityType%NeighborsPoints(:,i)+3
                LabelNgh=labelstmp(ll(1),ll(2),ll(3))

                tmplabels => labelstmp(ll(1)-1:ll(1)+1,ll(2)-1:ll(2)+1,ll(3)-1:ll(3)+1)

                !-----------------------------------------------------------------------
                ! TODO: the isEnclosedByLabel method could use the above iterator
                !-----------------------------------------------------------------------
                IF (LabelNgh.EQ.-LabelTo.AND.IsEnclosedByLabel_BGConnectivity(tmplabels,LabelTo)) THEN
                   !-----------------------------------------------------------------------
                   ! Remove the parent that got enclosed; it had a the label
                   ! of the currently expanding region and a candidate label of 0.
                   !-----------------------------------------------------------------------
                   tmpParticle=MCMCParticle(0,zero)

                   nx=INT(REAL(ld(1),MK)/MCMCcellsize(1))
                   IF (nx.EQ.MCMCcellNm(1)) nx=nx-1

                   ny=INT(REAL(ld(2),MK)/MCMCcellsize(2))
                   IF (ny.EQ.MCMCcellNm(2)) ny=ny-1

                   nz=INT(REAL(ld(3),MK)/MCMCcellsize(3))
                   IF (nz.EQ.MCMCcellNm(3)) nz=nz-1

                   cbox=1+nx+MCMCcellNm(1)*(ny+nz*MCMCcellNm(2))

                   info=MCMCEraseCandidatesFromContainers(ld(1),ld(2),ld(3),cbox, &
                   &    tmpParticle,LabelTo,.TRUE.)
                   or_fail("MCMCEraseCandidatesFromContainers")

                   !-----------------------------------------------------------------------
                   ! update the label image (we're not using the update
                   ! mechanism of the optimizer anymore):
                   !-----------------------------------------------------------------------
                   labelstmp(ll(1),ll(2),ll(3))=LabelTo

                   ALLOCATE(seed,STAT=info)
                   or_fail_alloc("seed")
                   CALL seed%add(ld(1),ld(2),ld(3),LabelNgh)
                   CALL MCMCLabelImageHistory%push(seed,info)
                   or_fail("MCMCLabelImageHistory%push")
                ENDIF
             ENDDO !i=1,BG_ConnectivityType%NumberOfNeighbors

             labelstmp(3,3,3)=StoreLabel

             !-----------------------------------------------------------------------
             ! Figure out if a point renders to a candidate. These are
             ! all the FG-neighbors with a different label that are not yet
             ! candidates of the currently expanding region.
             !-----------------------------------------------------------------------
             DO i=1,FG_ConnectivityType%NumberOfNeighbors
                ld=coord+FG_ConnectivityType%NeighborsPoints(:,i)
                IF (ANY(ld.LT.1.OR.ld.GT.Nm)) CYCLE
                ll=FG_ConnectivityType%NeighborsPoints(:,i)+3
                LabelNgh=labelstmp(ll(1),ll(2),ll(3))

                IF (ABS(LabelNgh).NE.LabelTo.AND.LabelNgh.NE.FORBIDDEN) THEN
                   !-----------------------------------------------------------------------
                   ! check if there is no other mother (hence the particle is
                   ! not in the container yet). This we could do by checking the
                   ! neighborhood of the label image or by checking the
                   ! cointainers.
                   ! Here: we check the (not yet updated!) label image.
                   !-----------------------------------------------------------------------
                   HasOtherMother=.FALSE.
                   DO j=1,FG_ConnectivityType%NumberOfNeighbors
                      dd=ll+FG_ConnectivityType%NeighborsPoints(:,i)
                      IF (ABS(labelstmp(dd(1),dd(2),dd(3))).EQ.LabelTo) THEN
                         HasOtherMother=.TRUE.
                         EXIT
                      ENDIF
                   ENDDO
                   IF (.NOT.HasOtherMother) THEN
                      !-----------------------------------------------------------------------
                      ! This is a new child. It's current label we have to read
                      ! from the label image, the candidate label is the label of
                      ! the currently expanding region.
                      !-----------------------------------------------------------------------
                      tmplabels => labelstmp(ll(1)-1:ll(1)+1,ll(2)-1:ll(2)+1,ll(3)-1:ll(3)+1)
                      proposal=MCMCproposal(tmplabels)

                      tmpParticle=MCMCParticle(LabelTo,proposal)

                      nx=INT(REAL(ld(1),MK)/MCMCcellsize(1))
                      IF (nx.EQ.MCMCcellNm(1)) nx=nx-1

                      ny=INT(REAL(ld(2),MK)/MCMCcellsize(2))
                      IF (ny.EQ.MCMCcellNm(2)) ny=ny-1

                      nz=INT(REAL(ld(3),MK)/MCMCcellsize(3))
                      IF (nz.EQ.MCMCcellNm(3)) nz=nz-1

                      cbox=1+nx+MCMCcellNm(1)*(ny+nz*MCMCcellNm(2))

                      info=MCMCInsertCandidatesToContainers(ld(1),ld(2),ld(3),cbox, &
                      &    tmpParticle,ABS(LabelNgh),.TRUE.,Replaced)
                      or_fail("MCMCInsertCandidatesToContainers")

                      IF (.NOT.Replaced) THEN
                         labelstmp(ll(1),ll(2),ll(3))=SIGN(LabelNgh,-1)

                         ALLOCATE(seed,STAT=info)
                         or_fail_alloc("seed")
                         CALL seed%add(ld(1),ld(2),ld(3),LabelNgh)
                         CALL MCMCLabelImageHistory%push(seed,info)
                         or_fail("MCMCLabelImageHistory%push")
                      ENDIF
                   ENDIF !.NOT.HasOtherMother
                ENDIF !LabelNgh.NE.LabelTo.AND.LabelNgh.NE.FORBIDDEN
             ENDDO !i=1,FG_ConnectivityType%NumberOfNeighbors
          ENDIF ! LabelTo.NE.0
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        END FUNCTION MCMCAddAndRemoveParticlesWhenMove_3d

        FUNCTION MCMCupdateProposalsAndFilterTopologyInNeighborhood_2d(labels_,coord,Nm) RESULT(info)
          !!! Proposal update and topology control are combined in this method for
          !!! efficiency reason.
          !!! All particles in the 3x3x3 unit cube to be updated. We use the
          !!! label image to find the particles. Their proposal is updated. If the
          !!! particles do not fulfill the topological constraints, they will be
          !!! removed from the containers.
          !!! Note: floating particles will not be updated!
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_error
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_error,ppm_err_sub_failed

          USE ppm_rc_module_topologicalnumber, ONLY : BG_ConnectivityType
          USE ppm_rc_module_global, ONLY : MCMCParticle
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER, CONTIGUOUS, DIMENSION(:,:), POINTER       :: labels_
          INTEGER,             DIMENSION(:),   INTENT(IN   ) :: coord
          INTEGER,             DIMENSION(:),   INTENT(IN   ) :: Nm
          INTEGER                                            :: info
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(MCMCParticle), DIMENSION(9) :: MCMCParticles

          REAL(ppm_kind_double) :: t0

          INTEGER, DIMENSION(:,:), POINTER :: labelstmp
          INTEGER, DIMENSION(:,:), POINTER :: tmplabels
          INTEGER, DIMENSION(2)            :: ll,ld
          INTEGER                          :: ContainerLabel
          INTEGER                          :: cLabel
          INTEGER                          :: nsize
          INTEGER                          :: i,j
          INTEGER                          :: nx,ny
          INTEGER                          :: cbox

          LOGICAL :: Replaced

          CHARACTER(LEN=ppm_char) :: caller='MCMCupdateProposalsAndFilterTopologyInNeighborhood'
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          labelstmp => labels_(coord(1)-2:coord(1)+2,coord(2)-2:coord(2)+2)
          tmplabels  => labelstmp(2:4,2:4)

          cLabel=ABS(tmplabels(2,2))

          ! nsize will be zero if coords are outside domain
          nsize=MCMCgetRegularParticlesAtIndex(tmplabels,coord(1),coord(2),Nm,MCMCParticles)

          IF (nsize.GT.0) THEN
             nx=INT(REAL(coord(1),MK)/MCMCcellsize(1))
             IF (nx.EQ.MCMCcellNm(1)) nx=nx-1

             ny=INT(REAL(coord(2),MK)/MCMCcellsize(2))
             IF (ny.EQ.MCMCcellNm(2)) ny=ny-1

             cbox=1+nx+ny*MCMCcellNm(1)

             DO i=1,nsize
                ContainerLabel=MERGE(cLabel,MCMCParticles(i)%candlabel,MCMCParticles(i)%candlabel.EQ.0)
                IF (MCMCIsParticleTopoValid(tmplabels,MCMCParticles(i)%candlabel)) THEN
                   ! replace the particle with the new propsal
                   info=MCMCInsertCandidatesToContainers(coord(1),coord(2),cbox, &
                   &    MCMCParticles(i),ContainerLabel,.TRUE.,Replaced)
                   or_fail("MCMCInsertCandidatesToContainers")
                ELSE
                   ! remove the particle
                   info=MCMCEraseCandidatesFromContainers(coord(1),coord(2),cbox, &
                   &    MCMCParticles(i),ContainerLabel,.TRUE.)
                   or_fail("MCMCEraseCandidatesFromContainers")
                ENDIF
             ENDDO !i=1,nsize
          ENDIF !(nsize.GT.0)

          DO i=1,BG_ConnectivityType%NumberOfNeighbors
             ld=coord+BG_ConnectivityType%NeighborsPoints(:,i)
             IF (ANY(ld.LT.1.OR.ld.GT.Nm)) CYCLE
             ll=BG_ConnectivityType%NeighborsPoints(:,i)+3

             tmplabels => labelstmp(ll(1)-1:ll(1)+1,ll(2)-1:ll(2)+1)

             cLabel=ABS(tmplabels(2,2))

             nsize=MCMCgetRegularParticlesAtIndex(tmplabels,ld(1),ld(2),Nm,MCMCParticles)

             IF (nsize.GT.0) THEN
                nx=INT(REAL(ld(1),MK)/MCMCcellsize(1))
                IF (nx.EQ.MCMCcellNm(1)) nx=nx-1

                ny=INT(REAL(ld(2),MK)/MCMCcellsize(2))
                IF (ny.EQ.MCMCcellNm(2)) ny=ny-1

                cbox=1+nx+ny*MCMCcellNm(1)

                DO j=1,nsize
                   ContainerLabel=MERGE(cLabel,MCMCParticles(j)%candlabel,MCMCParticles(j)%candlabel.EQ.0)
                   IF (MCMCIsParticleTopoValid(tmplabels,MCMCParticles(j)%candlabel)) THEN
                      ! replace the particle with the new propsal
                      info=MCMCInsertCandidatesToContainers(ld(1),ld(2),cbox, &
                      &    MCMCParticles(j),ContainerLabel,.TRUE.,Replaced)
                      or_fail("MCMCInsertCandidatesToContainers")
                   ELSE
                      ! remove the particle
                      info=MCMCEraseCandidatesFromContainers(ld(1),ld(2),cbox, &
                      &    MCMCParticles(j),ContainerLabel,.TRUE.)
                      or_fail("MCMCEraseCandidatesFromContainers")
                   ENDIF
                ENDDO !j=1,nsize
             ENDIF !(nsize.GT.0)
          ENDDO !i=1,BG_ConnectivityType%NumberOfNeighbors
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        END FUNCTION MCMCupdateProposalsAndFilterTopologyInNeighborhood_2d

        FUNCTION MCMCupdateProposalsAndFilterTopologyInNeighborhood_3d(labels_,coord,Nm) RESULT(info)
          !!! Proposal update and topology control are combined in this method for
          !!! efficiency reason.
          !!! All particles in the 3x3x3 unit cube to be updated. We use the
          !!! label image to find the particles. Their proposal is updated. If the
          !!! particles do not fulfill the topological constraints, they will be
          !!! removed from the containers.
          !!! Note: floating particles will not be updated!
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_error
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_error,ppm_err_sub_failed

          USE ppm_rc_module_topologicalnumber, ONLY : BG_ConnectivityType
          USE ppm_rc_module_global, ONLY : MCMCParticle
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: labels_
          INTEGER,             DIMENSION(:),     INTENT(IN   ) :: coord
          INTEGER,             DIMENSION(:),     INTENT(IN   ) :: Nm
          INTEGER                                              :: info
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(MCMCParticle), DIMENSION(27) :: MCMCParticles

          REAL(ppm_kind_double) :: t0

          INTEGER, DIMENSION(:,:,:), POINTER :: labelstmp
          INTEGER, DIMENSION(:,:,:), POINTER :: tmplabels
          INTEGER, DIMENSION(3)              :: ll,ld
          INTEGER                            :: ContainerLabel
          INTEGER                            :: cLabel
          INTEGER                            :: nsize
          INTEGER                            :: i,j
          INTEGER                            :: nx,ny,nz
          INTEGER                            :: cbox

          LOGICAL :: Replaced

          CHARACTER(LEN=ppm_char) :: caller='MCMCupdateProposalsAndFilterTopologyInNeighborhood'
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          labelstmp => labels_(coord(1)-2:coord(1)+2,coord(2)-2:coord(2)+2,coord(3)-2:coord(3)+2)
          tmplabels  => labelstmp(2:4,2:4,2:4)

          cLabel=ABS(tmplabels(2,2,2))

          nsize=MCMCgetRegularParticlesAtIndex(tmplabels,coord(1),coord(2),coord(3),Nm,MCMCParticles)

          IF (nsize.GT.0) THEN
             nx=INT(REAL(coord(1),MK)/MCMCcellsize(1))
             IF (nx.EQ.MCMCcellNm(1)) nx=nx-1

             ny=INT(REAL(coord(2),MK)/MCMCcellsize(2))
             IF (ny.EQ.MCMCcellNm(2)) ny=ny-1

             nz=INT(REAL(coord(3),MK)/MCMCcellsize(3))
             IF (nz.EQ.MCMCcellNm(3)) nz=nz-1

             cbox=1+nx+MCMCcellNm(1)*(ny+nz*MCMCcellNm(2))

             DO i=1,nsize
                ContainerLabel=MERGE(cLabel,MCMCParticles(i)%candlabel,MCMCParticles(i)%candlabel.EQ.0)
                IF (MCMCIsParticleTopoValid(tmplabels,MCMCParticles(i)%candlabel)) THEN
                   ! replace the particle with the new propsal
                   info=MCMCInsertCandidatesToContainers(coord(1),coord(2),coord(3),cbox, &
                   &    MCMCParticles(i),ContainerLabel,.TRUE.,Replaced)
                   or_fail("MCMCInsertCandidatesToContainers")
                ELSE
                   ! remove the particle
                   info=MCMCEraseCandidatesFromContainers(coord(1),coord(2),coord(3),cbox, &
                   &    MCMCParticles(i),ContainerLabel,.TRUE.)
                   or_fail("MCMCEraseCandidatesFromContainers")
                ENDIF
             ENDDO !i=1,nsize
          ENDIF !(nsize.GE.1)

          DO i=1,BG_ConnectivityType%NumberOfNeighbors
             ld=coord+BG_ConnectivityType%NeighborsPoints(:,i)
             ! TODO
             ! TOCHECK: I think this one should be changed in parallel version
             IF (ANY(ld.LT.1.OR.ld.GT.Nm)) CYCLE
             ll=BG_ConnectivityType%NeighborsPoints(:,i)+3

             tmplabels => labelstmp(ll(1)-1:ll(1)+1,ll(2)-1:ll(2)+1,ll(3)-1:ll(3)+1)

             cLabel=ABS(tmplabels(2,2,2))

             nsize=MCMCgetRegularParticlesAtIndex(tmplabels,ld(1),ld(2),ld(3),Nm,MCMCParticles)

             IF (nsize.GT.0) THEN
                nx=INT(REAL(ld(1),MK)/MCMCcellsize(1))
                IF (nx.EQ.MCMCcellNm(1)) nx=nx-1

                ny=INT(REAL(ld(2),MK)/MCMCcellsize(2))
                IF (ny.EQ.MCMCcellNm(2)) ny=ny-1

                nz=INT(REAL(ld(3),MK)/MCMCcellsize(3))
                IF (nz.EQ.MCMCcellNm(3)) nz=nz-1

                cbox=1+nx+MCMCcellNm(1)*(ny+nz*MCMCcellNm(2))

                DO j=1,nsize
                   ContainerLabel=MERGE(cLabel,MCMCParticles(j)%candlabel,MCMCParticles(j)%candlabel.EQ.0)
                   IF (MCMCIsParticleTopoValid(tmplabels,MCMCParticles(j)%candlabel)) THEN
                      ! replace the particle with the new propsal
                      info=MCMCInsertCandidatesToContainers(ld(1),ld(2),ld(3),cbox, &
                      &    MCMCParticles(j),ContainerLabel,.TRUE.,Replaced)
                      or_fail("MCMCInsertCandidatesToContainers")
                   ELSE
                      ! remove the particle
                      info=MCMCEraseCandidatesFromContainers(ld(1),ld(2),ld(3),cbox, &
                      &    MCMCParticles(j),ContainerLabel,.TRUE.)
                      or_fail("MCMCEraseCandidatesFromContainers")
                   ENDIF
                ENDDO !j=1,nsize
             ENDIF !(nsize.GT.0)
          ENDDO !i=1,BG_ConnectivityType%NumberOfNeighbors
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        END FUNCTION MCMCupdateProposalsAndFilterTopologyInNeighborhood_3d

        LOGICAL FUNCTION MCMCIsRegularParticle_2d (coord1,coord2,cbox,MCMCParticle_)
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER,            INTENT(IN   ) :: coord1
          INTEGER,            INTENT(IN   ) :: coord2
          INTEGER,            INTENT(IN   ) :: cbox
          TYPE(MCMCParticle), INTENT(IN   ) :: MCMCParticle_
          !!! Absolute value label from regionsLabel at Coordinates (coord1,coord2)
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !
          !-------------------------------------------------------------------------
          MCMCIsRegularParticle_2d=MCMCRegularParticlesInCell(cbox)%containselement(coord1,coord2,MCMCParticle_)
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCIsRegularParticle_2d

        LOGICAL FUNCTION MCMCIsRegularParticle_3d (coord1,coord2,coord3,cbox,MCMCParticle_)
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER,            INTENT(IN   ) :: coord1
          INTEGER,            INTENT(IN   ) :: coord2
          INTEGER,            INTENT(IN   ) :: coord3
          INTEGER,            INTENT(IN   ) :: cbox
          TYPE(MCMCParticle), INTENT(IN   ) :: MCMCParticle_
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !
          !-------------------------------------------------------------------------
          MCMCIsRegularParticle_3d=MCMCRegularParticlesInCell(cbox)%containselement(coord1,coord2,coord3,MCMCParticle_)
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCIsRegularParticle_3d

        LOGICAL FUNCTION MCMCParticleHasFloatingProperty_2d(labels_,coords,CandidateLabel)
          USE ppm_rc_module_topologicalnumber, ONLY : IsSingleFGPoint, &
          &   IsEnclosedByLabel_BGConnectivity
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER, CONTIGUOUS, DIMENSION(:,:), POINTER       :: labels_
          INTEGER,             DIMENSION(:),   INTENT(IN   ) :: coords
          INTEGER,                             INTENT(IN   ) :: CandidateLabel
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          INTEGER :: cLabel

          IF (CandidateLabel.GT.0) THEN
             ! This is a daughter. It is floating if there is no mother, i.e. if
             ! there is no corresponding region in the FG neighborhood.
             MCMCParticleHasFloatingProperty_2d=IsSingleFGPoint(coords,labels_,CandidateLabel)
             RETURN
          ELSE
             ! else: this is a potential mother (candidate label is 0). According
             ! to the BG connectivity it fulfills the floating property
             ! only if there is no other region in the BG neighborhood. Otherwise
             ! this pixel might well go to the BG label without changing the topo.
             cLabel=ABS(labels_(coords(1),coords(2)))
             MCMCParticleHasFloatingProperty_2d=IsEnclosedByLabel_BGConnectivity(coords,labels_,cLabel)
             RETURN
          ENDIF
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCParticleHasFloatingProperty_2d

        LOGICAL FUNCTION MCMCParticleHasFloatingProperty___2d(labels_,coords,LabelIn,CandidateLabel)
          USE ppm_rc_module_topologicalnumber, ONLY : IsSingleFGPoint, &
          &   IsEnclosedByLabel_BGConnectivity
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER, CONTIGUOUS, DIMENSION(:,:), POINTER       :: labels_
          INTEGER,             DIMENSION(:),   INTENT(IN   ) :: coords
          INTEGER,                             INTENT(IN   ) :: LabelIn
          INTEGER,                             INTENT(IN   ) :: CandidateLabel
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          IF (CandidateLabel.GT.0) THEN
             ! This is a daughter. It is floating if there is no mother, i.e. if
             ! there is no corresponding region in the FG neighborhood.
             MCMCParticleHasFloatingProperty___2d=IsSingleFGPoint(coords,labels_,CandidateLabel)
             RETURN
          ELSE
             ! else: this is a potential mother (candidate label is 0). According
             ! to the BG connectivity it fulfills the floating property
             ! only if there is no other region in the BG neighborhood. Otherwise
             ! this pixel might well go to the BG label without changing the topo.
             MCMCParticleHasFloatingProperty___2d=IsEnclosedByLabel_BGConnectivity(coords,labels_,LabelIn)
             RETURN
          ENDIF
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCParticleHasFloatingProperty___2d

        LOGICAL FUNCTION MCMCParticleHasFloatingProperty__2d(tmplabels,CandidateLabel)
          USE ppm_rc_module_topologicalnumber, ONLY : IsSingleFGPoint, &
          &   IsEnclosedByLabel_BGConnectivity
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER, DIMENSION(:,:), POINTER       :: tmplabels
          INTEGER,                 INTENT(IN   ) :: CandidateLabel
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          INTEGER :: cLabel

          IF (CandidateLabel.GT.0) THEN
             ! This is a daughter. It is floating if there is no mother, i.e. if
             ! there is no corresponding region in the FG neighborhood.
             MCMCParticleHasFloatingProperty__2d=IsSingleFGPoint(tmplabels,CandidateLabel)
             RETURN
          ELSE
             ! else: this is a potential mother (candidate label is 0). According
             ! to the BG connectivity it fulfills the floating property
             ! only if there is no other region in the BG neighborhood. Otherwise
             ! this pixel might well go to the BG label without changing the topo.
             cLabel=ABS(tmplabels(2,2))
             MCMCParticleHasFloatingProperty__2d=IsEnclosedByLabel_BGConnectivity(tmplabels,cLabel)
             RETURN
          ENDIF
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCParticleHasFloatingProperty__2d

        LOGICAL FUNCTION MCMCParticleHasFloatingProperty____2d(tmplabels,LabelIn,CandidateLabel)
          USE ppm_rc_module_topologicalnumber, ONLY : IsSingleFGPoint, &
          &   IsEnclosedByLabel_BGConnectivity
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER, DIMENSION(:,:), POINTER       :: tmplabels
          INTEGER,                 INTENT(IN   ) :: LabelIn
          INTEGER,                 INTENT(IN   ) :: CandidateLabel
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          IF (CandidateLabel.GT.0) THEN
             !-----------------------------------------------------------------------
             ! This is a daughter. It is floating if there is no mother, i.e. if
             ! there is no corresponding region in the FG neighborhood.
             !-----------------------------------------------------------------------
             MCMCParticleHasFloatingProperty____2d=IsSingleFGPoint(tmplabels,CandidateLabel)
             RETURN
          ELSE
             !-----------------------------------------------------------------------
             ! else: this is a potential mother (candidate label is 0). According
             ! to the BG connectivity it fulfills the floating property
             ! only if there is no other region in the BG neighborhood. Otherwise
             ! this pixel might well go to the BG label without changing the topo.
             !-----------------------------------------------------------------------
             MCMCParticleHasFloatingProperty____2d=IsEnclosedByLabel_BGConnectivity(tmplabels,LabelIn)
             RETURN
          ENDIF
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCParticleHasFloatingProperty____2d

        LOGICAL FUNCTION MCMCParticleHasFloatingProperty_3d(labels_,coords,CandidateLabel)
          USE ppm_rc_module_topologicalnumber, ONLY : IsSingleFGPoint, &
          &   IsEnclosedByLabel_BGConnectivity
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: labels_
          INTEGER,             DIMENSION(:),     INTENT(IN   ) :: coords
          INTEGER,                               INTENT(IN   ) :: CandidateLabel
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          INTEGER :: cLabel

          IF (CandidateLabel.GT.0) THEN
             ! This is a daughter. It is floating if there is no mother, i.e. if
             ! there is no corresponding region in the FG neighborhood.
             MCMCParticleHasFloatingProperty_3d=IsSingleFGPoint(coords,labels_,CandidateLabel)
             RETURN
          ELSE
             ! else: this is a potential mother (candidate label is 0). According
             ! to the BG connectivity it fulfills the floating property
             ! only if there is no other region in the BG neighborhood. Otherwise
             ! this pixel might well go to the BG label without changing the topo.
             cLabel=ABS(labels_(coords(1),coords(2),coords(3)))
             MCMCParticleHasFloatingProperty_3d=IsEnclosedByLabel_BGConnectivity(coords,labels_,cLabel)
             RETURN
          ENDIF
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCParticleHasFloatingProperty_3d

        LOGICAL FUNCTION MCMCParticleHasFloatingProperty___3d(labels_,coords,LabelIn,CandidateLabel)
          USE ppm_rc_module_topologicalnumber, ONLY : IsSingleFGPoint, &
          &   IsEnclosedByLabel_BGConnectivity
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: labels_
          INTEGER,             DIMENSION(:),     INTENT(IN   ) :: coords
          INTEGER,                               INTENT(IN   ) :: LabelIn
          INTEGER,                               INTENT(IN   ) :: CandidateLabel
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          IF (CandidateLabel.GT.0) THEN
             ! This is a daughter. It is floating if there is no mother, i.e. if
             ! there is no corresponding region in the FG neighborhood.
             MCMCParticleHasFloatingProperty___3d=IsSingleFGPoint(coords,labels_,CandidateLabel)
             RETURN
          ELSE
             ! else: this is a potential mother (candidate label is 0). According
             ! to the BG connectivity it fulfills the floating property
             ! only if there is no other region in the BG neighborhood. Otherwise
             ! this pixel might well go to the BG label without changing the topo.
             MCMCParticleHasFloatingProperty___3d=IsEnclosedByLabel_BGConnectivity(coords,labels_,LabelIn)
             RETURN
          ENDIF
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCParticleHasFloatingProperty___3d

        LOGICAL FUNCTION MCMCParticleHasFloatingProperty__3d(tmplabels,CandidateLabel)
          USE ppm_rc_module_topologicalnumber, ONLY : IsSingleFGPoint, &
          &   IsEnclosedByLabel_BGConnectivity
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER, DIMENSION(:,:,:), POINTER       :: tmplabels
          INTEGER,                   INTENT(IN   ) :: CandidateLabel
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          INTEGER :: cLabel

          IF (CandidateLabel.GT.0) THEN
             ! This is a daughter. It is floating if there is no mother, i.e. if
             ! there is no corresponding region in the FG neighborhood.
             MCMCParticleHasFloatingProperty__3d=IsSingleFGPoint(tmplabels,CandidateLabel)
             RETURN
          ELSE
             ! else: this is a potential mother (candidate label is 0). According
             ! to the BG connectivity it fulfills the floating property
             ! only if there is no other region in the BG neighborhood. Otherwise
             ! this pixel might well go to the BG label without changing the topo.
             cLabel=ABS(tmplabels(2,2,2))
             MCMCParticleHasFloatingProperty__3d=IsEnclosedByLabel_BGConnectivity(tmplabels,cLabel)
             RETURN
          ENDIF
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCParticleHasFloatingProperty__3d

        LOGICAL FUNCTION MCMCParticleHasFloatingProperty____3d(tmplabels,LabelIn,CandidateLabel)
          USE ppm_rc_module_topologicalnumber, ONLY : IsSingleFGPoint, &
          &   IsEnclosedByLabel_BGConnectivity
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER, DIMENSION(:,:,:), POINTER       :: tmplabels
          INTEGER,                   INTENT(IN   ) :: LabelIn
          INTEGER,                   INTENT(IN   ) :: CandidateLabel
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          IF (CandidateLabel.GT.0) THEN
             ! This is a daughter. It is floating if there is no mother, i.e. if
             ! there is no corresponding region in the FG neighborhood.
             MCMCParticleHasFloatingProperty____3d=IsSingleFGPoint(tmplabels,CandidateLabel)
             RETURN
          ELSE
             ! else: this is a potential mother (candidate label is 0). According
             ! to the BG connectivity it fulfills the floating property
             ! only if there is no other region in the BG neighborhood. Otherwise
             ! this pixel might well go to the BG label without changing the topo.
             MCMCParticleHasFloatingProperty____3d=IsEnclosedByLabel_BGConnectivity(tmplabels,LabelIn)
             RETURN
          ENDIF
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
          RETURN
        END FUNCTION MCMCParticleHasFloatingProperty____3d

        FUNCTION MCMCApplyParticle_2d(image_,labels_,coord,Nm,MCMCParticle_,DoSimulate) RESULT(info)
          !!! Updates the DS (label image and the regular particle containers) when
          !!! applying a particle. Note that floating particles will not be updated!
          !!! The method only ensures that L and the regular  particles are correct.
          !!! The method expects the label image NOT to be updated already. Else the
          !!! operations performed will be wrong.
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_error
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_err_alloc,ppm_error,ppm_err_sub_failed, &
          &   ppm_err_argument

          USE ppm_rc_module_global, ONLY : MCMCuseBiasedProposal,AllowFission, &
          &   AllowHandles
          USE ppm_rc_module_linkedlist, ONLY : ppm_rc_list,MCMCLabelImageHistory
          USE ppm_rc_module_topologicalnumber, ONLY : IsEnclosedByLabel_BGConnectivity
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          REAL(MK), CONTIGUOUS, DIMENSION(:,:), POINTER       :: image_
          INTEGER,  CONTIGUOUS, DIMENSION(:,:), POINTER       :: labels_
          INTEGER,              DIMENSION(:),   INTENT(IN   ) :: coord
          INTEGER,              DIMENSION(:),   INTENT(IN   ) :: Nm
          TYPE(MCMCParticle),                   INTENT(IN   ) :: MCMCParticle_
          LOGICAL,                              INTENT(IN   ) :: DoSimulate
          INTEGER                                             :: info
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(ppm_rc_list), POINTER ::  seed

          REAL(ppm_kind_double) :: t0

          INTEGER :: FromLabel
          INTEGER :: ToLabel
          INTEGER :: cLabel

          CHARACTER(LEN=ppm_char) :: caller='MCMCApplyParticle'
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          !-------------------------------------------------------------------------
          ! Maintain the regular particle container and the label image
          !-------------------------------------------------------------------------
          info=MCMCAddAndRemoveParticlesWhenMove(labels_,coord,Nm,MCMCParticle_)
          or_fail("MCMCAddAndRemoveParticlesWhenMove")

          !-------------------------------------------------------------------------
          ! Update the label image. The new point is either a contour particle or 0,
          ! therefor the negative label value is set.
          !-------------------------------------------------------------------------
          cLabel=labels_(coord(1),coord(2))

          FromLabel=ABS(cLabel)
          ToLabel=MCMCParticle_%candlabel

          check_true(<#FromLabel.NE.ToLabel#>,"particle and its candidate have the same label!")

          ALLOCATE(seed,STAT=info)
          or_fail_alloc("seed")
          CALL seed%add(coord(1),coord(2),cLabel)
          CALL MCMCLabelImageHistory%push(seed,info)
          or_fail("MCMCLabelImageHistory%push")

          labels_(coord(1),coord(2))=ToLabel

          !-------------------------------------------------------------------------
          ! Update the statistics of the propagating and the loser region.
          !-------------------------------------------------------------------------
          IF (.NOT.DoSimulate) THEN
             CALL e_data%UpdateStatisticsWhenJump(image_,coord,FromLabel,ToLabel,info)
             or_fail("e_data%UpdateStatisticsWhenJump")
          ENDIF

          !-------------------------------------------------------------------------
          ! Update the proposals for all particles in the neighborhood (as they
          ! might have changed).
          !-------------------------------------------------------------------------
          IF (MCMCuseBiasedProposal.OR..NOT.AllowFission.OR..NOT.AllowHandles) THEN
             info=MCMCupdateProposalsAndFilterTopologyInNeighborhood(labels_,coord,Nm)
             or_fail("MCMCupdateProposalsAndFilterTopologyInNeighborhood")
          ENDIF
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        END FUNCTION MCMCApplyParticle_2d

        FUNCTION MCMCApplyParticle_3d(image_,labels_,coord,Nm,MCMCParticle_,DoSimulate) RESULT(info)
          !!! Updates the DS (label image and the regular particle containers) when
          !!! applying a particle. Note that floating particles will not be updated!
          !!! The method only ensures that L and the regular  particles are correct.
          !!! The method expects the label image NOT to be updated already. Else the
          !!! operations performed will be wrong.
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_error
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_err_alloc,ppm_error,ppm_err_sub_failed, &
          &   ppm_err_argument

          USE ppm_rc_module_global, ONLY : MCMCuseBiasedProposal,AllowFission, &
          &   AllowHandles
          USE ppm_rc_module_linkedlist, ONLY : ppm_rc_list,MCMCLabelImageHistory
          USE ppm_rc_module_topologicalnumber, ONLY : IsEnclosedByLabel_BGConnectivity
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: image_
          INTEGER,  CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: labels_
          INTEGER,              DIMENSION(:),     INTENT(IN   ) :: coord
          INTEGER,              DIMENSION(:),     INTENT(IN   ) :: Nm
          TYPE(MCMCParticle),                     INTENT(IN   ) :: MCMCParticle_
          LOGICAL,                                INTENT(IN   ) :: DoSimulate
          INTEGER                                               :: info
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(ppm_rc_list), POINTER ::  seed

          REAL(ppm_kind_double) :: t0

          INTEGER :: FromLabel
          INTEGER :: ToLabel
          INTEGER :: cLabel

          CHARACTER(LEN=ppm_char) :: caller='MCMCApplyParticle'
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          !-------------------------------------------------------------------------
          ! Maintain the regular particle container and the label image
          !-------------------------------------------------------------------------
          info=MCMCAddAndRemoveParticlesWhenMove(labels_,coord,Nm,MCMCParticle_)
          or_fail("MCMCAddAndRemoveParticlesWhenMove")

          !-------------------------------------------------------------------------
          ! Update the label image. The new point is either a contour particle or 0,
          ! therefor the negative label value is set.
          !-------------------------------------------------------------------------
          cLabel=labels_(coord(1),coord(2),coord(3))

          FromLabel=ABS(cLabel)
          ToLabel=MCMCParticle_%candlabel

          check_true(<#FromLabel.NE.ToLabel#>,"particle and its candidate have the same label!")

          ALLOCATE(seed,STAT=info)
          or_fail_alloc("seed")
          CALL seed%add(coord(1),coord(2),coord(3),cLabel)
          CALL MCMCLabelImageHistory%push(seed,info)
          or_fail("MCMCLabelImageHistory%push")

          labels_(coord(1),coord(2),coord(3))=ToLabel

          !-------------------------------------------------------------------------
          ! Update the statistics of the propagating and the loser region.
          !-------------------------------------------------------------------------
          IF (.NOT.DoSimulate) THEN
             CALL e_data%UpdateStatisticsWhenJump(image_,coord,FromLabel,ToLabel,info)
             or_fail("e_data%UpdateStatisticsWhenJump")
          ENDIF

          !-------------------------------------------------------------------------
          ! Update the proposals for all particles in the neighborhood (as they
          ! might have changed).
          !-------------------------------------------------------------------------
          IF (MCMCuseBiasedProposal.OR..NOT.AllowFission.OR..NOT.AllowHandles) THEN
             info=MCMCupdateProposalsAndFilterTopologyInNeighborhood(labels_,coord,Nm)
             or_fail("MCMCupdateProposalsAndFilterTopologyInNeighborhood")
          ENDIF
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        END FUNCTION MCMCApplyParticle_3d

        FUNCTION MCMCAddAndRemoveParticlesFromGhostMove_2d(labels_,coord,Nm,oldghost) RESULT(info)
          !!! Updates the DS (label image and the regular particle containers) when
          !!! applying a particle. Note that floating particles will not be updated!
          !!! The method only ensures that L and the regular  particles are correct.
          !!! The method expects the label image NOT to be updated already. Else the
          !!! operations performed will be wrong.
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_error
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_err_alloc,ppm_error,ppm_err_sub_failed

          USE ppm_rc_module_global, ONLY : FORBIDDEN
          USE ppm_rc_module_topologicalnumber, ONLY : FG_ConnectivityType, &
          &   BG_ConnectivityType,IsEnclosedByLabel_BGConnectivity
          USE ppm_rc_module_linkedlist, ONLY : ppm_rc_list,MCMCLabelImageHistory
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER, CONTIGUOUS, DIMENSION(:,:), POINTER       :: labels_
          INTEGER,             DIMENSION(:),   INTENT(IN   ) :: coord
          INTEGER,             DIMENSION(:),   INTENT(IN   ) :: Nm
          INTEGER,                             INTENT(IN   ) :: oldghost
          INTEGER                                            :: info
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(MCMCParticle) :: tmpParticle

          TYPE(ppm_rc_list), POINTER ::  seed

          REAL(ppm_kind_double) :: t0
          REAL(MK)              :: proposal

          INTEGER, DIMENSION(:,:), POINTER :: labelstmp
          INTEGER, DIMENSION(:,:), POINTER :: tmplabels

          INTEGER, DIMENSION(2)            :: ll,ld,dd
          INTEGER                          :: LabelFrom
          INTEGER                          :: LabelTo
          INTEGER                          :: LabelNgh
          INTEGER                          :: i,j
          INTEGER                          :: nx,ny
          INTEGER                          :: cbox

          LOGICAL :: Replaced
          LOGICAL :: HasOtherMother

          CHARACTER(LEN=ppm_char) :: caller='MCMCAddAndRemoveParticlesFromGhostMove'
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          ! Initilize an neighborhooditerator for fast access to the label image
          ! in the vicinity of the particle of interest.
          labelstmp => labels_(coord(1)-2:coord(1)+2,coord(2)-2:coord(2)+2)

          ! Initialize locals from the particle
          LabelFrom=ABS(oldghost)
          LabelTo=ABS(labelstmp(3,3))

          ! What particle would be added or removed to the contour lists
          ! of the currently shrinking region:
          ! (compare to also AddNeighborsAtRemove(..))
          IF (LabelFrom.NE.0) THEN
             ! a FG region is shrinking (Moving particle is either a parent or child)
             DO i=1,BG_ConnectivityType%NumberOfNeighbors
                ld=coord+BG_ConnectivityType%NeighborsPoints(:,i)
                IF (ANY(ld.LT.1.OR.ld.GT.Nm)) CYCLE
                ll=BG_ConnectivityType%NeighborsPoints(:,i)+3
                ! We do not use ABS here, to find internal points
                LabelNgh=labelstmp(ll(1),ll(2))
                ! Check if new points enter the contour (internal point becomes a parent):
                ! Internal points of this region are positive:
                IF (LabelNgh.EQ.LabelFrom) THEN
                   tmplabels => labelstmp(ll(1)-1:ll(1)+1,ll(2)-1:ll(2)+1)
                   proposal=MCMCproposal(tmplabels)

                   ! Internal point becomes a parent
                   tmpParticle=MCMCParticle(0,proposal)

                   nx=INT(REAL(ld(1),MK)/MCMCcellsize(1))
                   IF (nx.EQ.MCMCcellNm(1)) nx=nx-1

                   ny=INT(REAL(ld(2),MK)/MCMCcellsize(2))
                   IF (ny.EQ.MCMCcellNm(2)) ny=ny-1

                   cbox=1+nx+ny*MCMCcellNm(1)

                   ! Insert new particle in container
                   info=MCMCInsertCandidatesToContainers(ld(1),ld(2),cbox, &
                   &    tmpParticle,LabelFrom,.TRUE.,Replaced)
                   or_fail("MCMCInsertCandidatesToContainers")

                   IF (.NOT.Replaced) THEN
                      labelstmp(ll(1),ll(2))=-LabelFrom

                      ALLOCATE(seed,STAT=info)
                      or_fail_alloc("seed")
                      CALL seed%add(ld(1),ld(2),LabelFrom)
                      CALL MCMCLabelImageHistory%push(seed,info)
                      or_fail("MCMCLabelImageHistory%push")
                   ENDIF
                ENDIF !LabelNgh.GT.0.AND.LabelNgh.EQ.LabelFrom
             ENDDO !i=1,BG_ConnectivityType%NumberOfNeighbors

             DO i=1,FG_ConnectivityType%NumberOfNeighbors
                ld=coord+FG_ConnectivityType%NeighborsPoints(:,i)
                IF (ANY(ld.LT.1.OR.ld.GT.Nm)) CYCLE
                ll=FG_ConnectivityType%NeighborsPoints(:,i)+3
                LabelNgh=ABS(labelstmp(ll(1),ll(2)))

                ! check if there are FG-neighbors with no other mother of the
                ! same label --> orphan.
                ! we first 'simulate' the move:
                IF (LabelNgh.NE.LabelFrom) THEN
                   ! check if this neighbor has other mothers from this label:
                   HasOtherMother=.FALSE.
                   DO j=1,FG_ConnectivityType%NumberOfNeighbors
                      dd=ll+FG_ConnectivityType%NeighborsPoints(:,i)
                      IF (ABS(labelstmp(dd(1),dd(2))).EQ.LabelFrom) THEN
                         HasOtherMother=.TRUE.
                         EXIT
                      ENDIF
                   ENDDO
                   IF (.NOT.HasOtherMother) THEN
                      ! The orphan has label equal to what we read from the
                      ! label image and has a candidate label of the
                      ! currently shrinking region.
                      tmplabels => labelstmp(ll(1)-1:ll(1)+1,ll(2)-1:ll(2)+1)
                      proposal=MCMCproposal(tmplabels)

                      tmpParticle=MCMCParticle(LabelFrom,proposal)

                      nx=INT(REAL(ld(1),MK)/MCMCcellsize(1))
                      IF (nx.EQ.MCMCcellNm(1)) nx=nx-1

                      ny=INT(REAL(ld(2),MK)/MCMCcellsize(2))
                      IF (ny.EQ.MCMCcellNm(2)) ny=ny-1

                      cbox=1+nx+ny*MCMCcellNm(1)

                      info=MCMCEraseCandidatesFromContainers(ld(1),ld(2),cbox, &
                      &    tmpParticle,LabelNgh,.TRUE.)
                      or_fail("MCMCEraseCandidatesFromContainers")
                   ENDIF !.NOT.HasOtherMother
                ENDIF !LabelNgh.NE.LabelFrom
             ENDDO !i=1,FG_ConnectivityType%NumberOfNeighbors
          ENDIF ! LabelFrom.NE.0

          ! Growing: figure out the changes of candidates for the expanding region
          IF (LabelTo.NE.0) THEN
             ! we are growing

             ! Neighbors: Figure out what (neighboring)mother points are going
             ! to be interior points:

             DO i=1,BG_ConnectivityType%NumberOfNeighbors
                ld=coord+BG_ConnectivityType%NeighborsPoints(:,i)
                IF (ANY(ld.LT.1.OR.ld.GT.Nm)) CYCLE
                ll=BG_ConnectivityType%NeighborsPoints(:,i)+3
                LabelNgh=labelstmp(ll(1),ll(2))

                tmplabels => labelstmp(ll(1)-1:ll(1)+1,ll(2)-1:ll(2)+1)

                ! TODO: the isEnclosedByLabel method could use the above iterator
                IF (LabelNgh.EQ.-LabelTo.AND.IsEnclosedByLabel_BGConnectivity(tmplabels,LabelTo)) THEN
                   ! Remove the parent that got enclosed; it had a the label
                   ! of the currently expanding region and a candidate label
                   ! of 0.
                   tmpParticle=MCMCParticle(0,zero)

                   nx=INT(REAL(ld(1),MK)/MCMCcellsize(1))
                   IF (nx.EQ.MCMCcellNm(1)) nx=nx-1

                   ny=INT(REAL(ld(2),MK)/MCMCcellsize(2))
                   IF (ny.EQ.MCMCcellNm(2)) ny=ny-1

                   cbox=1+nx+ny*MCMCcellNm(1)

                   info=MCMCEraseCandidatesFromContainers(ld(1),ld(2),cbox, &
                   &    tmpParticle,LabelTo,.TRUE.)
                   or_fail("MCMCEraseCandidatesFromContainers")

                   ! update the label image (we're not using the update
                   ! mechanism of the optimizer anymore):
                   labelstmp(ll(1),ll(2))=LabelTo

                   ALLOCATE(seed,STAT=info)
                   or_fail_alloc("seed")
                   CALL seed%add(ld(1),ld(2),LabelNgh)
                   CALL MCMCLabelImageHistory%push(seed,info)
                   or_fail("MCMCLabelImageHistory%push")
                ENDIF
             ENDDO !i=1,BG_ConnectivityType%NumberOfNeighbors

             ! Figure out if a point renders to a candidate. These are
             ! all the FG-neighbors with a different label that are not yet
             ! candidates of the currently expanding region.
             DO i=1,FG_ConnectivityType%NumberOfNeighbors
                ld=coord+FG_ConnectivityType%NeighborsPoints(:,i)
                IF (ANY(ld.LT.1.OR.ld.GT.Nm)) CYCLE
                ll=FG_ConnectivityType%NeighborsPoints(:,i)+3
                LabelNgh=labelstmp(ll(1),ll(2))

                IF (ABS(LabelNgh).NE.LabelTo.AND.LabelNgh.NE.FORBIDDEN) THEN
                   ! check if there is no other mother (hence the particle is
                   ! not in the container yet). This we could do by checking the
                   ! neighborhood of the label image or by checking the
                   ! cointainers.
                   ! Here: we check the (not yet updated!) label image.
                   HasOtherMother=.FALSE.
                   DO j=1,FG_ConnectivityType%NumberOfNeighbors
                      dd=ll+FG_ConnectivityType%NeighborsPoints(:,i)
                      IF (ABS(labelstmp(dd(1),dd(2))).EQ.LabelTo) THEN
                         HasOtherMother=.TRUE.
                         EXIT
                      ENDIF
                   ENDDO
                   IF (.NOT.HasOtherMother) THEN
                      ! This is a new child. It's current label we have to read
                      ! from the label image, the candidate label is the label of
                      ! the currently expanding region.
                      tmplabels => labelstmp(ll(1)-1:ll(1)+1,ll(2)-1:ll(2)+1)
                      proposal=MCMCproposal(tmplabels)

                      tmpParticle=MCMCParticle(LabelTo,proposal)

                      nx=INT(REAL(ld(1),MK)/MCMCcellsize(1))
                      IF (nx.EQ.MCMCcellNm(1)) nx=nx-1

                      ny=INT(REAL(ld(2),MK)/MCMCcellsize(2))
                      IF (ny.EQ.MCMCcellNm(2)) ny=ny-1

                      cbox=1+nx+ny*MCMCcellNm(1)

                      info=MCMCInsertCandidatesToContainers(ld(1),ld(2),cbox, &
                      &    tmpParticle,ABS(LabelNgh),.TRUE.,Replaced)
                      or_fail("MCMCInsertCandidatesToContainers")

                      IF (.NOT.Replaced) THEN
                         labelstmp(ll(1),ll(2))=SIGN(LabelNgh,-1)

                         ALLOCATE(seed,STAT=info)
                         or_fail_alloc("seed")
                         CALL seed%add(ld(1),ld(2),LabelNgh)
                         CALL MCMCLabelImageHistory%push(seed,info)
                         or_fail("MCMCLabelImageHistory%push")
                      ENDIF
                   ENDIF !.NOT.HasOtherMother
                ENDIF !LabelNgh.NE.LabelTo.AND.LabelNgh.NE.FORBIDDEN
             ENDDO !i=1,FG_ConnectivityType%NumberOfNeighbors
          ENDIF ! LabelTo.NE.0
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        END FUNCTION MCMCAddAndRemoveParticlesFromGhostMove_2d

        FUNCTION MCMCAddAndRemoveParticlesFromGhostMove_3d(labels_,coord,Nm,oldghost) RESULT(info)
          ! Updates the DS (label image and the regular particle containers) when
          ! applying a particle. Note that floating particles will not be updated!
          ! The method only ensures that L and the regular  particles are correct.
          ! The method expects the label image NOT to be updated already. Else the
          ! operations performed will be wrong.
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_error
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_err_alloc,ppm_error,ppm_err_sub_failed

          USE ppm_rc_module_global, ONLY : FORBIDDEN
          USE ppm_rc_module_topologicalnumber, ONLY : FG_ConnectivityType, &
          &   BG_ConnectivityType,IsEnclosedByLabel_BGConnectivity
          USE ppm_rc_module_linkedlist, ONLY : ppm_rc_list,MCMCLabelImageHistory
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: labels_
          INTEGER,             DIMENSION(:),     INTENT(IN   ) :: coord
          INTEGER,             DIMENSION(:),     INTENT(IN   ) :: Nm
          INTEGER,                               INTENT(IN   ) :: oldghost
          INTEGER                                              :: info
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(MCMCParticle) :: tmpParticle

          TYPE(ppm_rc_list), POINTER ::  seed

          REAL(ppm_kind_double) :: t0
          REAL(MK)              :: proposal

          INTEGER, DIMENSION(:,:,:), POINTER :: labelstmp
          INTEGER, DIMENSION(:,:,:), POINTER :: tmplabels

          INTEGER, DIMENSION(3)              :: ll,ld,dd
          INTEGER                            :: LabelFrom
          INTEGER                            :: LabelTo
          INTEGER                            :: LabelNgh
          INTEGER                            :: i,j
          INTEGER                            :: nx,ny,nz
          INTEGER                            :: cbox


          LOGICAL :: Replaced
          LOGICAL :: HasOtherMother

          CHARACTER(LEN=ppm_char) :: caller='MCMCAddAndRemoveParticlesFromGhostMove'
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          ! Initilize an neighborhooditerator for fast access to the label image
          ! in the vicinity of the particle of interest.
          labelstmp => labels_(coord(1)-2:coord(1)+2,coord(2)-2:coord(2)+2,coord(3)-2:coord(3)+2)

          ! Initialize locals from the particle
          LabelFrom=ABS(oldghost)
          LabelTo=ABS(labelstmp(3,3,3))

          ! What particle would be added or removed to the contour lists
          ! of the currently shrinking region:
          ! (compare to also AddNeighborsAtRemove(..))
          IF (LabelFrom.NE.0) THEN
             ! a FG region is shrinking
             DO i=1,BG_ConnectivityType%NumberOfNeighbors
                ld=coord+BG_ConnectivityType%NeighborsPoints(:,i)
                IF (ANY(ld.LT.1.OR.ld.GT.Nm)) CYCLE
                ll=BG_ConnectivityType%NeighborsPoints(:,i)+3
                LabelNgh=labelstmp(ll(1),ll(2),ll(3))
                ! Check if new points enter the contour (internal point becomes a parent):
                ! Internal points of this region are positive:
                IF (LabelNgh.EQ.LabelFrom) THEN
                   tmplabels => labelstmp(ll(1)-1:ll(1)+1,ll(2)-1:ll(2)+1,ll(3)-1:ll(3)+1)
                   proposal=MCMCproposal(tmplabels)

                   tmpParticle=MCMCParticle(0,proposal)

                   nx=INT(REAL(ld(1),MK)/MCMCcellsize(1))
                   IF (nx.EQ.MCMCcellNm(1)) nx=nx-1

                   ny=INT(REAL(ld(2),MK)/MCMCcellsize(2))
                   IF (ny.EQ.MCMCcellNm(2)) ny=ny-1

                   nz=INT(REAL(ld(3),MK)/MCMCcellsize(3))
                   IF (nz.EQ.MCMCcellNm(3)) nz=nz-1

                   cbox=1+nx+MCMCcellNm(1)*(ny+nz*MCMCcellNm(2))

                   info=MCMCInsertCandidatesToContainers(ld(1),ld(2),ld(3),cbox, &
                   &    tmpParticle,LabelFrom,.TRUE.,Replaced)
                   or_fail("MCMCInsertCandidatesToContainers")

                   IF (.NOT.Replaced) THEN
                      labelstmp(ll(1),ll(2),ll(3))=-LabelFrom

                      ALLOCATE(seed,STAT=info)
                      or_fail_alloc("seed")
                      CALL seed%add(ld(1),ld(2),ld(3),LabelFrom)
                      CALL MCMCLabelImageHistory%push(seed,info)
                      or_fail("MCMCLabelImageHistory%push")
                   ENDIF
                ENDIF !LabelNgh.GT.0.AND.LabelNgh.EQ.LabelFrom
             ENDDO !i=1,BG_ConnectivityType%NumberOfNeighbors

             DO i=1,FG_ConnectivityType%NumberOfNeighbors
                ld=coord+FG_ConnectivityType%NeighborsPoints(:,i)
                IF (ANY(ld.LT.1.OR.ld.GT.Nm)) CYCLE
                ll=FG_ConnectivityType%NeighborsPoints(:,i)+3
                LabelNgh=ABS(labelstmp(ll(1),ll(2),ll(3)))

                ! check if there are FG-neighbors with no other mother of the
                ! same label --> orphan.
                ! we first 'simulate' the move:
                IF (LabelNgh.NE.LabelFrom) THEN
                   ! check if this neighbor has other mothers from this label:
                   HasOtherMother=.FALSE.
                   DO j=1,FG_ConnectivityType%NumberOfNeighbors
                      dd=ll+FG_ConnectivityType%NeighborsPoints(:,i)
                      IF (ABS(labelstmp(dd(1),dd(2),dd(3))).EQ.LabelFrom) THEN
                         HasOtherMother=.TRUE.
                         EXIT
                      ENDIF
                   ENDDO
                   IF (.NOT.HasOtherMother) THEN
                      ! The orphan has label equal to what we read from the
                      ! label image and has a candidate label of the
                      ! currently shrinking region.
                      tmplabels => labelstmp(ll(1)-1:ll(1)+1,ll(2)-1:ll(2)+1,ll(3)-1:ll(3)+1)
                      proposal=MCMCproposal(tmplabels)

                      tmpParticle=MCMCParticle(LabelFrom,proposal)

                      nx=INT(REAL(ld(1),MK)/MCMCcellsize(1))
                      IF (nx.EQ.MCMCcellNm(1)) nx=nx-1

                      ny=INT(REAL(ld(2),MK)/MCMCcellsize(2))
                      IF (ny.EQ.MCMCcellNm(2)) ny=ny-1

                      nz=INT(REAL(ld(3),MK)/MCMCcellsize(3))
                      IF (nz.EQ.MCMCcellNm(3)) nz=nz-1

                      cbox=1+nx+MCMCcellNm(1)*(ny+nz*MCMCcellNm(2))

                      info=MCMCEraseCandidatesFromContainers(ld(1),ld(2),ld(3),cbox, &
                      &    tmpParticle,LabelNgh,.TRUE.)
                      or_fail("MCMCEraseCandidatesFromContainers")
                   ENDIF !.NOT.HasOtherMother
                ENDIF !LabelNgh.NE.LabelFrom
             ENDDO !i=1,FG_ConnectivityType%NumberOfNeighbors
          ENDIF ! LabelFrom.NE.0

          ! Growing: figure out the changes of candidates for the expanding region
          IF (LabelTo.NE.0) THEN
             ! we are growing

             ! Neighbors: Figure out what (neighboring)mother points are going
             ! to be interior points:
             DO i=1,BG_ConnectivityType%NumberOfNeighbors
                ld=coord+BG_ConnectivityType%NeighborsPoints(:,i)
                IF (ANY(ld.LT.1.OR.ld.GT.Nm)) CYCLE
                ll=BG_ConnectivityType%NeighborsPoints(:,i)+3
                LabelNgh=labelstmp(ll(1),ll(2),ll(3))

                tmplabels => labelstmp(ll(1)-1:ll(1)+1,ll(2)-1:ll(2)+1,ll(3)-1:ll(3)+1)

                ! TODO: the isEnclosedByLabel method could use the above iterator
                IF (LabelNgh.EQ.-LabelTo.AND.IsEnclosedByLabel_BGConnectivity(tmplabels,LabelTo)) THEN
                   ! Remove the parent that got enclosed; it had a the label
                   ! of the currently expanding region and a candidate label
                   ! of 0.
                   tmpParticle=MCMCParticle(0,zero)

                   nx=INT(REAL(ld(1),MK)/MCMCcellsize(1))
                   IF (nx.EQ.MCMCcellNm(1)) nx=nx-1

                   ny=INT(REAL(ld(2),MK)/MCMCcellsize(2))
                   IF (ny.EQ.MCMCcellNm(2)) ny=ny-1

                   nz=INT(REAL(ld(3),MK)/MCMCcellsize(3))
                   IF (nz.EQ.MCMCcellNm(3)) nz=nz-1

                   cbox=1+nx+MCMCcellNm(1)*(ny+nz*MCMCcellNm(2))

                   info=MCMCEraseCandidatesFromContainers(ld(1),ld(2),ld(3),cbox, &
                   &    tmpParticle,LabelTo,.TRUE.)
                   or_fail("MCMCEraseCandidatesFromContainers")

                   ! update the label image (we're not using the update
                   ! mechanism of the optimizer anymore):
                   labelstmp(ll(1),ll(2),ll(3))=LabelTo

                   ALLOCATE(seed,STAT=info)
                   or_fail_alloc("seed")
                   CALL seed%add(ld(1),ld(2),ld(3),LabelNgh)
                   CALL MCMCLabelImageHistory%push(seed,info)
                   or_fail("MCMCLabelImageHistory%push")
                ENDIF
             ENDDO !i=1,BG_ConnectivityType%NumberOfNeighbors

             ! Figure out if a point renders to a candidate. These are
             ! all the FG-neighbors with a different label that are not yet
             ! candidates of the currently expanding region.
             DO i=1,FG_ConnectivityType%NumberOfNeighbors
                ld=coord+FG_ConnectivityType%NeighborsPoints(:,i)
                IF (ANY(ld.LT.1.OR.ld.GT.Nm)) CYCLE
                ll=FG_ConnectivityType%NeighborsPoints(:,i)+3
                LabelNgh=labelstmp(ll(1),ll(2),ll(3))

                IF (ABS(LabelNgh).NE.LabelTo.AND.LabelNgh.NE.FORBIDDEN) THEN
                   ! check if there is no other mother (hence the particle is
                   ! not in the container yet). This we could do by checking the
                   ! neighborhood of the label image or by checking the
                   ! cointainers.
                   ! Here: we check the (not yet updated!) label image.
                   HasOtherMother=.FALSE.
                   DO j=1,FG_ConnectivityType%NumberOfNeighbors
                      dd=ll+FG_ConnectivityType%NeighborsPoints(:,i)
                      IF (ABS(labelstmp(dd(1),dd(2),dd(3))).EQ.LabelTo) THEN
                         HasOtherMother=.TRUE.
                         EXIT
                      ENDIF
                   ENDDO
                   IF (.NOT.HasOtherMother) THEN
                      ! This is a new child. It's current label we have to read
                      ! from the label image, the candidate label is the label of
                      ! the currently expanding region.
                      tmplabels => labelstmp(ll(1)-1:ll(1)+1,ll(2)-1:ll(2)+1,ll(3)-1:ll(3)+1)
                      proposal=MCMCproposal(tmplabels)

                      tmpParticle=MCMCParticle(LabelTo,proposal)

                      nx=INT(REAL(ld(1),MK)/MCMCcellsize(1))
                      IF (nx.EQ.MCMCcellNm(1)) nx=nx-1

                      ny=INT(REAL(ld(2),MK)/MCMCcellsize(2))
                      IF (ny.EQ.MCMCcellNm(2)) ny=ny-1

                      nz=INT(REAL(ld(3),MK)/MCMCcellsize(3))
                      IF (nz.EQ.MCMCcellNm(3)) nz=nz-1

                      cbox=1+nx+MCMCcellNm(1)*(ny+nz*MCMCcellNm(2))

                      info=MCMCInsertCandidatesToContainers(ld(1),ld(2),ld(3),cbox, &
                      &    tmpParticle,ABS(LabelNgh),.TRUE.,Replaced)
                      or_fail("MCMCInsertCandidatesToContainers")

                      IF (.NOT.Replaced) THEN
                         labelstmp(ll(1),ll(2),ll(3))=SIGN(LabelNgh,-1)

                         ALLOCATE(seed,STAT=info)
                         or_fail_alloc("seed")
                         CALL seed%add(ld(1),ld(2),ld(3),LabelNgh)
                         CALL MCMCLabelImageHistory%push(seed,info)
                         or_fail("MCMCLabelImageHistory%push")
                      ENDIF
                   ENDIF !.NOT.HasOtherMother
                ENDIF !LabelNgh.NE.LabelTo.AND.LabelNgh.NE.FORBIDDEN
             ENDDO !i=1,FG_ConnectivityType%NumberOfNeighbors
          ENDIF ! LabelTo.NE.0
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        END FUNCTION MCMCAddAndRemoveParticlesFromGhostMove_3d

        FUNCTION MCMCApplyParticleFromGhost_2d(labels_,coord,Nm,oldghost) RESULT(info)
          !!! Updates the DS (label image and the regular particle containers) when
          !!! applying a particle. Note that floating particles will not be updated!
          !!! The method only ensures that L and the regular  particles are correct.
          !!! The method expects the label image NOT to be updated already. Else the
          !!! operations performed will be wrong.
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_error
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_error,ppm_err_sub_failed,ppm_err_argument

          USE ppm_rc_module_global, ONLY : MCMCuseBiasedProposal,AllowFission, &
          &   AllowHandles
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER,  CONTIGUOUS, DIMENSION(:,:), POINTER       :: labels_
          INTEGER,              DIMENSION(:),   INTENT(IN   ) :: coord
          INTEGER,              DIMENSION(:),   INTENT(IN   ) :: Nm
          INTEGER,                              INTENT(IN   ) :: oldghost
          INTEGER                                             :: info
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          REAL(ppm_kind_double) :: t0

          CHARACTER(LEN=ppm_char) :: caller='MCMCApplyParticleFromGhost'
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          ! Maintain the regular particle container and the label image
          info=MCMCAddAndRemoveParticlesFromGhostMove(labels_,coord,Nm,oldghost)
          or_fail("MCMCAddAndRemoveParticlesFromGhostMove")

          ! Update the proposals for all particles in the neighborhood (as they
          ! might have changed).
          IF (MCMCuseBiasedProposal.OR..NOT.AllowFission.OR..NOT.AllowHandles) THEN
             info=MCMCupdateProposalsAndFilterTopologyInNeighborhood(labels_,coord,Nm)
             or_fail("MCMCupdateProposalsAndFilterTopologyInNeighborhood")
          ENDIF
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        END FUNCTION MCMCApplyParticleFromGhost_2d

        FUNCTION MCMCApplyParticleFromGhost_3d(labels_,coord,Nm,oldghost) RESULT(info)
          !!! Updates the DS (label image and the regular particle containers) when
          !!! applying a particle. Note that floating particles will not be updated!
          !!! The method only ensures that L and the regular  particles are correct.
          !!! The method expects the label image NOT to be updated already. Else the
          !!! operations performed will be wrong.
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_error
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_error,ppm_err_sub_failed,ppm_err_argument

          USE ppm_rc_module_global, ONLY : MCMCuseBiasedProposal,AllowFission, &
          &   AllowHandles
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER,  CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: labels_
          INTEGER,              DIMENSION(:),     INTENT(IN   ) :: coord
          INTEGER,              DIMENSION(:),     INTENT(IN   ) :: Nm
          INTEGER,                                INTENT(IN   ) :: oldghost
          INTEGER                                               :: info
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          REAL(ppm_kind_double) :: t0

          CHARACTER(LEN=ppm_char) :: caller='MCMCApplyParticleFromGhost'
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          ! Maintain the regular particle container and the label image
          info=MCMCAddAndRemoveParticlesFromGhostMove(labels_,coord,Nm,oldghost)
          or_fail("MCMCAddAndRemoveParticlesFromGhostMove")

          ! Update the proposals for all particles in the neighborhood (as they
          ! might have changed).
          IF (MCMCuseBiasedProposal.OR..NOT.AllowFission.OR..NOT.AllowHandles) THEN
             info=MCMCupdateProposalsAndFilterTopologyInNeighborhood(labels_,coord,Nm)
             or_fail("MCMCupdateProposalsAndFilterTopologyInNeighborhood")
          ENDIF
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        END FUNCTION MCMCApplyParticleFromGhost_3d

        LOGICAL FUNCTION MCMCFindParticleA(ppmrcdim,CCIndx,RPICSize,FPICSize,info)
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          ! TRUE result of the function means the move is rejected or not
          !-------------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_fatal,ppm_error_error
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_error,ppm_err_sub_failed, &
          &   ppm_err_argument,ppm_err_dealloc
          USE ppm_module_mesh_typedef, ONLY : ppm_t_subpatch_

          USE ppm_rc_module_global, ONLY : oned,mesh,labels,MCMCstepsize, &
          &   MCMCuseBiasedProposal
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
          INTEGER, INTENT(IN   ) :: CCIndx
          !!! Current Cell Box index
          INTEGER, INTENT(IN   ) :: RPICSize
          !!! Number of Regular particles in the current cell cbox
          INTEGER, INTENT(IN   ) :: FPICSize
          !!! Number of Floating particles in the current cell cbox
          INTEGER, INTENT(  OUT) :: info
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          CLASS(ppm_t_subpatch_), POINTER :: sbpitr

          REAL(ppm_kind_double)               :: t0
          REAL(MK), DIMENSION(:), ALLOCATABLE :: AllParticlesFwdProposals
          !!! For each particle within the region, calculate the proposal

          INTEGER, CONTIGUOUS, DIMENSION(:,:),   POINTER :: labels_2d
          INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER :: labels_3d
          INTEGER,             DIMENSION(ppmrcdim)       :: index_
          INTEGER                                        :: PC
          INTEGER                                        :: ParticleIndex
          INTEGER                                        :: cLabel
          INTEGER                                        :: ActiveCandidatesSize
          !!! Size of the AllParticlesFwdProposals

          LOGICAL :: ParticleAIsFloating
          LOGICAL :: IsDelParticleA

          CHARACTER(LEN=ppm_char) :: caller='MCMCFindParticleA'
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)
          !-------------------------------------------------------------------------
          ! For the insertion and deletion of floating particles, we use a
          ! constant energy, therefore the MH ratio be equal to the FBR.
          ! Note that we also could use the energy as if we would do the job as
          ! proposal distribution, but we would have to calculate the backward
          ! energy as well.
          !-------------------------------------------------------------------------
          MCMCFindParticleA=.TRUE.

          !-------------------------------------------------------------------------
          ! Number of particles in the current cell
          !-------------------------------------------------------------------------
          ActiveCandidatesSize=RPICSize+FPICSize

          IF (MCMCuseBiasedProposal) THEN
             info=MCMCRegularParticlesInCell(CCIndx)%packproposal(AllParticlesFwdProposals,ActiveCandidatesSize)
             or_fail("MCMCRegularParticlesInCell(CCIndx)%packproposal",ppm_error=ppm_error_fatal)

             info=MCMCFloatingParticlesInCell(CCIndx)%packproposal(AllParticlesFwdProposals,RPICSize+1,ActiveCandidatesSize)
             or_fail("MCMCFloatingParticlesInCell(CCIndx)%packproposal",ppm_error=ppm_error_fatal)

             !-------------------------------------------------------------------------
             ! Create a discrete distribution over particles
             !-------------------------------------------------------------------------
             info=ppm_rc_GenerateParticlesFwdProposalsDiscrDistr(AllParticlesFwdProposals,ActiveCandidatesSize)
             or_fail("ppm_rc_GenerateParticlesFwdProposalsDiscrDistr",ppm_error=ppm_error_fatal)
          ENDIF

          SELECT CASE (ppmrcdim)
          CASE (2)
             NULLIFY(labels_2d)
             sbpitr => mesh%subpatch%begin()
             DO WHILE (ASSOCIATED(sbpitr))
                CALL sbpitr%get_field(labels,labels_2d,info)
                or_fail("Failed to get field labels_2d data.",ppm_error=ppm_error_fatal)

                !-------------------------------------------------------------------------
                ! Find particle A:
                !-------------------------------------------------------------------------
                DO PC=1,MCMCstepsize
                   !-------------------------------------------------------------------------
                   ! Using biased proposal
                   !-------------------------------------------------------------------------
                   IF (MCMCuseBiasedProposal) THEN
                      ParticleIndex=ppm_rc_GetPartDistrIndex()
                      MCMC_q_A(PC)=REAL(AllParticlesFwdProposals(ParticleIndex),ppm_kind_double)
                   !-------------------------------------------------------------------------
                   ! Using unbiased proposal
                   !-------------------------------------------------------------------------
                   ELSE
                      ParticleIndex=ppm_rc_Saru_IPRNG(ActiveCandidatesSize)
                      MCMC_q_A(PC)=oned
                   ENDIF

                   !-------------------------------------------------------------------------
                   ! If particle is a regular particle
                   !-------------------------------------------------------------------------
                   IF (ParticleIndex.LE.RPICSize) THEN
                      MCMCActiveCandidates => MCMCRegularParticlesInCell(CCIndx)

                      ParticleAIsFloating=.FALSE.
                   !-------------------------------------------------------------------------
                   ! If particle is a floating particle
                   !-------------------------------------------------------------------------
                   ELSE
                      MCMCActiveCandidates => MCMCFloatingParticlesInCell(CCIndx)

                      !-------------------------------------------------------------------------
                      ! correction on the index
                      !-------------------------------------------------------------------------
                      ParticleIndex=ParticleIndex-RPICSize

                      ParticleAIsFloating=.TRUE.
                   ENDIF

                   !-------------------------------------------------------------------------
                   ! Get the particle
                   !-------------------------------------------------------------------------
                   MCMC_CandidateMove(PC)=MCMCActiveCandidates%elementAt(ParticleIndex,index_)
                   IF (MCMC_CandidateMove(PC)%candlabel.EQ.-1) THEN
                      fail("There is no particle in this position!!!",ppm_error=ppm_error_fatal)
                   ENDIF

                   !-------------------------------------------------------------------------
                   ! particle position
                   !-------------------------------------------------------------------------
                   MCMC_CandidateMove_Index(:,PC)=index_
                   !-------------------------------------------------------------------------
                   ! label at this index
                   !-------------------------------------------------------------------------
                   cLabel=ABS(labels_2d(index_(1),index_(2)))

                   IF (ParticleAIsFloating) THEN
                      !-------------------------------------------------------------------------
                      ! A floating particle was selected.
                      !-------------------------------------------------------------------------
                      MCMC_Particle_A_IsFloating(PC)=.TRUE.

                      !-------------------------------------------------------------------------
                      ! Immediately reject the move (self transition) if the label
                      ! at the particular position is the same. Since this is a
                      ! self transition we keep the particle in the proposals.
                      !-------------------------------------------------------------------------
                      IF (cLabel.EQ.MCMC_CandidateMove(PC)%candlabel) THEN
                         !!! reject
                         MCMCFindParticleA=.FALSE.
                         GOTO 8888
                      ENDIF
                      !-------------------------------------------------------------------------
                      ! METHOD 5: do it (the regular particle will be deleted) --> 49.06 % / avg. 430 fp
                      !-------------------------------------------------------------------------

                      !-------------------------------------------------------------------------
                      ! METHOD 4: reject (because the backward prob is 0) --> 49.2 %
                      !-------------------------------------------------------------------------
                      IF (MCMCIsRegularParticle(index_(1),index_(2),CCIndx,MCMC_CandidateMove(PC))) THEN
                         !-------------------------------------------------------------------------
                         ! reject
                         !-------------------------------------------------------------------------
                         MCMCFindParticleA=.FALSE.
                         GOTO 8888
                      ENDIF

                      !-------------------------------------------------------------------------
                      ! if we are here, do the job
                      !-------------------------------------------------------------------------
                      IsDelParticleA=MCMCEraseFloatingParticle(index_(1),index_(2),CCIndx, &
                      & MCMC_CandidateMove(PC),MCMC_CandidateMove(PC),.TRUE.)

                      !-------------------------------------------------------------------------
                      ! METHOD 3 = METHOD 4
                      !-------------------------------------------------------------------------
                      ! METHOD 2: apply the floating particle if there is no regular
                      ! one. If there is a regular one, we apply this and keep the
                      ! float. (convert the float to regular), i.e. we just don't delete
                      ! the floating p. --> gets bi-modal with mean 48.75%, avg f.p. 8250
                      !-------------------------------------------------------------------------

                      !-------------------------------------------------------------------------
                      ! METHOD 1: just do it (remove both particles)
                      ! convert it to a real particle (floating particles always store
                      ! the region label they belong to). We need to check if we the
                      ! need to set the candidate label to 0.
                      !-------------------------------------------------------------------------
                   ENDIF !ParticleAIsFloating

                   !-------------------------------------------------------------------------
                   ! Fill in some properties for this particles
                   !-------------------------------------------------------------------------
                   MCMC_LabelsBeforeJump_A(PC)=cLabel

                   MCMC_q_A(PC)=MCMC_q_A(PC)/(MCMCTotalNormalizer+MCMCTotalNormalizerlocal)
                ENDDO !PC=1,MCMCstepsize

                sbpitr => mesh%subpatch%next()
             ENDDO !ASSOCIATED(sbpitr)
             NULLIFY(labels_2d)
          CASE (3)
             NULLIFY(labels_3d)
             sbpitr => mesh%subpatch%begin()
             DO WHILE (ASSOCIATED(sbpitr))
                CALL sbpitr%get_field(labels,labels_3d,info)
                or_fail("Failed to get field labels_3d data.",ppm_error=ppm_error_fatal)

                !-------------------------------------------------------------------------
                ! Find particle A:
                !-------------------------------------------------------------------------
                DO PC=1,MCMCstepsize
                   IF (MCMCuseBiasedProposal) THEN
                      ParticleIndex=ppm_rc_GetPartDistrIndex()
                      MCMC_q_A(PC)=REAL(AllParticlesFwdProposals(ParticleIndex),ppm_kind_double)
                   ELSE
                      ParticleIndex=ppm_rc_Saru_IPRNG(ActiveCandidatesSize)
                      MCMC_q_A(PC)=oned
                   ENDIF

                   !-------------------------------------------------------------------------
                   ! If particle is a regular particle
                   !-------------------------------------------------------------------------
                   IF (ParticleIndex.LE.RPICSize) THEN
                      MCMCActiveCandidates => MCMCRegularParticlesInCell(CCIndx)

                      ParticleAIsFloating=.FALSE.
                   !-------------------------------------------------------------------------
                   ! If particle is a floating particle
                   !-------------------------------------------------------------------------
                   ELSE
                      MCMCActiveCandidates => MCMCFloatingParticlesInCell(CCIndx)

                      !-------------------------------------------------------------------------
                      ! correction on the index
                      !-------------------------------------------------------------------------
                      ParticleIndex=ParticleIndex-RPICSize

                      ParticleAIsFloating=.TRUE.
                   ENDIF

                   MCMC_CandidateMove(PC)=MCMCActiveCandidates%elementAt(ParticleIndex,index_)
                   IF (MCMC_CandidateMove(PC)%candlabel.EQ.-1) THEN
                      fail("There is no particle in this position!!!",ppm_error=ppm_error_fatal)
                   ENDIF

                   !-------------------------------------------------------------------------
                   ! particle position
                   !-------------------------------------------------------------------------
                   MCMC_CandidateMove_Index(:,PC)=index_
                   !-------------------------------------------------------------------------
                   ! labels at this index
                   !-------------------------------------------------------------------------
                   cLabel=ABS(labels_3d(index_(1),index_(2),index_(3)))

                   IF (ParticleAIsFloating) THEN
                      !-------------------------------------------------------------------------
                      ! A floating particle was selected.
                      !-------------------------------------------------------------------------
                      MCMC_Particle_A_IsFloating(PC)=.TRUE.

                      !-------------------------------------------------------------------------
                      ! Immediately accept the move (self transition) if the label
                      ! at the particular position is the same. Since this is a
                      ! self transition we keep the particle in the proposals.
                      !-------------------------------------------------------------------------
                      IF (cLabel.EQ.MCMC_CandidateMove(PC)%candlabel) THEN
                         ! reject
                         MCMCFindParticleA=.FALSE.
                         GOTO 8888
                      ENDIF
                      !-------------------------------------------------------------------------
                      ! METHOD 5: do it (the regular particle will be deleted) --> 49.06 % / avg. 430 fp
                      !-------------------------------------------------------------------------

                      !-------------------------------------------------------------------------
                      ! METHOD 4: reject (because the backward prob is 0) --> 49.2 %
                      !-------------------------------------------------------------------------
                      IF (MCMCIsRegularParticle(index_(1),index_(2),index_(3),CCIndx,MCMC_CandidateMove(PC))) THEN
                         !-------------------------------------------------------------------------
                         ! reject
                         !-------------------------------------------------------------------------
                         MCMCFindParticleA=.FALSE.
                         GOTO 8888
                      ENDIF

                      !-------------------------------------------------------------------------
                      ! if we are here, do the job
                      !-------------------------------------------------------------------------
                      IsDelParticleA=MCMCEraseFloatingParticle(index_(1),index_(2),index_(3),CCIndx, &
                      & MCMC_CandidateMove(PC),MCMC_CandidateMove(PC),.TRUE.)

                      !-------------------------------------------------------------------------
                      ! METHOD 3 = METHOD 4
                      !-------------------------------------------------------------------------
                      !-------------------------------------------------------------------------
                      ! METHOD 2: apply the floating particle if there is no regular
                      ! one. If there is a regular one, we apply this and keep the
                      ! float. (convert the float to regular), i.e. we just don't delete
                      ! the floating p. --> gets bi-modal with mean 48.75%, avg f.p. 8250
                      !-------------------------------------------------------------------------
                      !-------------------------------------------------------------------------
                      ! METHOD 1: just do it (remove both particles)
                      ! convert it to a real particle (floating particles always store
                      ! the region label they belong to). We need to check if we the
                      ! need to set the candidate label to 0.
                      !-------------------------------------------------------------------------
                   ENDIF !ParticleAIsFloating

                   !-------------------------------------------------------------------------
                   ! Fill in some properties for this particles
                   !-------------------------------------------------------------------------
                   MCMC_LabelsBeforeJump_A(PC)=cLabel

                   MCMC_q_A(PC)=MCMC_q_A(PC)/(MCMCTotalNormalizer+MCMCTotalNormalizerlocal)
                ENDDO !PC=1,MCMCstepsize

                sbpitr => mesh%subpatch%next()
             ENDDO
             NULLIFY(labels_3d)
          END SELECT

        8888 CONTINUE
          IF (MCMCuseBiasedProposal) THEN
             !-------------------------------------------------------------------------
             ! Make sure that there is no discrete distribution available in memory
             !-------------------------------------------------------------------------
             info=ppm_rc_DestroyParticlesDiscrDistr()
             or_fail("ppm_rc_DestroyParticlesDiscrDistr")
             !-------------------------------------------------------------------------
             ! Free memory
             !-------------------------------------------------------------------------
             DEALLOCATE(AllParticlesFwdProposals,STAT=info)
             or_fail_dealloc("AllParticlesFwdProposals")
          ENDIF
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        END FUNCTION MCMCFindParticleA

        FUNCTION MCMCReject_2d(image_,labels_,AppliedParticleOrigLabels) RESULT(info)
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_fatal
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_error,ppm_err_dealloc,ppm_err_sub_failed

          USE ppm_rc_module_linkedlist, ONLY : ppm_rc_link,ppm_rc_list, &
          &   MCMCParticleInContainerHistory,MCMCLabelImageHistory,     &
          &   MCMCFloatingParticleInContainerHistory
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          REAL(MK), CONTIGUOUS, DIMENSION(:,:), POINTER       :: image_

          INTEGER,  CONTIGUOUS, DIMENSION(:,:), POINTER       :: labels_

          TYPE(ppm_rc_list),                    INTENT(INOUT) :: AppliedParticleOrigLabels

          INTEGER                                             :: info
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(ppm_rc_list),  POINTER :: seed
          TYPE(ppm_rc_link),  POINTER :: seedlink

          TYPE(MCMCHistoryParticle) :: tmpHistoryParticle

          TYPE(MCMCParticle) :: tmpParticle

          REAL(ppm_kind_double) :: t0

          INTEGER, DIMENSION(:), POINTER :: seedn
          INTEGER                        :: coord1
          INTEGER                        :: coord2
          INTEGER                        :: cbox
          INTEGER                        :: nx,ny

          CHARACTER(LEN=ppm_char) :: caller='MCMCReject'

          LOGICAL :: Replaced
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)
          !-------------------------------------------------------------------------
          ! First, recover the theoretical state:
          ! - recover the label image
          ! - recover the particle set
          ! As we have some (redundant) speedup statistics we also need to:
          ! - recover the statistics (by simulation of the backward particle)
          ! - recover the proposal normalizers (taken care of within the
          !   insertion/deletion methods).
          !-------------------------------------------------------------------------
          seed => MCMCParticleInContainerHistory%last()
          DO WHILE (ASSOCIATED(seed))
             seedn => seed%first%getValue(coord1,coord2,tmpHistoryParticle)

             nx=INT(REAL(coord1,MK)/MCMCcellsize(1))
             IF (nx.EQ.MCMCcellNm(1)) nx=nx-1

             ny=INT(REAL(coord2,MK)/MCMCcellsize(2))
             IF (ny.EQ.MCMCcellNm(2)) ny=ny-1

             cbox=1+nx+MCMCcellNm(1)*ny

             tmpParticle=MCMCParticle(tmpHistoryParticle%candlabel,tmpHistoryParticle%proposal)

             IF (tmpHistoryParticle%wasadded) THEN
                info=MCMCEraseCandidatesFromContainers(coord1,coord2,cbox, &
                &    tmpParticle,tmpHistoryParticle%orglabel,.FALSE.)
                or_fail("MCMCEraseCandidatesFromContainers",ppm_error=ppm_error_fatal)
             ELSE
                info=MCMCInsertCandidatesToContainers(coord1,coord2,cbox, &
                &    tmpParticle,tmpHistoryParticle%orglabel,.FALSE.,Replaced)
                or_fail("MCMCInsertCandidatesToContainers",ppm_error=ppm_error_fatal)
             ENDIF

             seed => MCMCParticleInContainerHistory%prev()
          ENDDO

          seed => MCMCFloatingParticleInContainerHistory%last()
          DO WHILE (ASSOCIATED(seed))
             seedn => seed%first%getValue(coord1,coord2,tmpHistoryParticle)

             nx=INT(REAL(coord1,MK)/MCMCcellsize(1))
             IF (nx.EQ.MCMCcellNm(1)) nx=nx-1

             ny=INT(REAL(coord2,MK)/MCMCcellsize(2))
             IF (ny.EQ.MCMCcellNm(2)) ny=ny-1

             cbox=1+nx+MCMCcellNm(1)*ny

             tmpParticle=MCMCParticle(tmpHistoryParticle%candlabel,tmpHistoryParticle%proposal)

             IF (tmpHistoryParticle%wasadded) THEN
                Replaced=MCMCEraseFloatingParticle(coord1,coord2,cbox,tmpParticle,.FALSE.)
             ELSE
                info=MCMCInsertFloatingParticleCumulative(coord1,coord2,cbox,tmpParticle,.FALSE.)
                or_fail("MCMCInsertFloatingParticleCumulative",ppm_error=ppm_error_fatal)
             ENDIF

             seed => MCMCFloatingParticleInContainerHistory%prev()
          ENDDO

          !-------------------------------------------------------------------------
          ! recover the label image:
          !-------------------------------------------------------------------------
          seed => MCMCLabelImageHistory%last()
          DO WHILE (ASSOCIATED(seed))
             seedn => seed%first%getValue()

             labels_(seedn(1),seedn(2))=seedn(3)

             seed => MCMCLabelImageHistory%prev()
          ENDDO

          !-------------------------------------------------------------------------
          ! recover the statistics:
          !-------------------------------------------------------------------------
          seedlink => AppliedParticleOrigLabels%first
          DO WHILE (ASSOCIATED(seedlink))
             seedn => seedlink%getValue()

             !-------------------------------------------------------------------------
             ! coord1        =seedn(1)
             ! coord2        =seedn(2)
             ! part%candlabel=seedn(3)  -- aFromLabel
             ! OriginalLabel =seedn(4)  -- aToLabel
             !-------------------------------------------------------------------------
             CALL e_data%UpdateStatisticsWhenJump(image_,seedn,seedn(3),seedn(4),info)
             or_fail("e_data%UpdateStatisticsWhenJump",ppm_error=ppm_error_fatal)

             seedlink => seedlink%nextLink()
          ENDDO
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        END FUNCTION MCMCReject_2d

        FUNCTION MCMCReject_3d(image_,labels_,AppliedParticleOrigLabels) RESULT(info)
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_error_fatal
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_error,ppm_err_dealloc,ppm_err_sub_failed

          USE ppm_rc_module_linkedlist, ONLY : ppm_rc_link,ppm_rc_list, &
          &   MCMCParticleInContainerHistory,MCMCLabelImageHistory,     &
          &   MCMCFloatingParticleInContainerHistory
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: image_

          INTEGER,  CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: labels_

          TYPE(ppm_rc_list),                      INTENT(IN   ) :: AppliedParticleOrigLabels

          INTEGER                                               :: info
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(ppm_rc_list),  POINTER :: seed
          TYPE(ppm_rc_link),  POINTER :: seedlink

          TYPE(MCMCHistoryParticle) :: tmpHistoryParticle

          TYPE(MCMCParticle) :: tmpParticle

          REAL(ppm_kind_double) :: t0

          INTEGER, DIMENSION(:), POINTER :: seedn
          INTEGER                        :: coord1
          INTEGER                        :: coord2
          INTEGER                        :: coord3
          INTEGER                        :: cbox
          INTEGER                        :: nx,ny,nz

          CHARACTER(LEN=ppm_char) :: caller='MCMCReject'

          LOGICAL :: Replaced
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)
          !-------------------------------------------------------------------------
          ! First, recover the theoretical state:
          ! - recover the label image
          ! - recover the particle set
          ! As we have some (redundant) speedup statistics we also need to:
          ! - recover the statistics (by simulation of the backward particle)
          ! - recover the proposal normalizers (taken care of within the
          !   insertion/deletion methods).
          !-------------------------------------------------------------------------
          seed => MCMCParticleInContainerHistory%last()
          DO WHILE (ASSOCIATED(seed))
             seedn => seed%first%getValue(coord1,coord2,coord3,tmpHistoryParticle)

             nx=INT(REAL(coord1,MK)/MCMCcellsize(1))
             IF (nx.EQ.MCMCcellNm(1)) nx=nx-1

             ny=INT(REAL(coord2,MK)/MCMCcellsize(2))
             IF (ny.EQ.MCMCcellNm(2)) ny=ny-1

             nz=INT(REAL(coord3,MK)/MCMCcellsize(3))
             IF (nz.EQ.MCMCcellNm(3)) nz=nz-1

             cbox=1+nx+MCMCcellNm(1)*(ny+nz*MCMCcellNm(2))

             tmpParticle=MCMCParticle(tmpHistoryParticle%candlabel,tmpHistoryParticle%proposal)

             IF (tmpHistoryParticle%wasadded) THEN
                info=MCMCEraseCandidatesFromContainers(coord1,coord2,coord3,cbox, &
                &    tmpParticle,tmpHistoryParticle%orglabel,.FALSE.)
                or_fail("MCMCEraseCandidatesFromContainers",ppm_error=ppm_error_fatal)
             ELSE
                info=MCMCInsertCandidatesToContainers(coord1,coord2,coord3,cbox, &
                &    tmpParticle,tmpHistoryParticle%orglabel,.FALSE.,Replaced)
                or_fail("MCMCInsertCandidatesToContainers",ppm_error=ppm_error_fatal)
             ENDIF

             seed => MCMCParticleInContainerHistory%prev()
          ENDDO

          seed => MCMCFloatingParticleInContainerHistory%last()
          DO WHILE (ASSOCIATED(seed))
             seedn => seed%first%getValue(coord1,coord2,coord3,tmpHistoryParticle)

             nx=INT(REAL(coord1,MK)/MCMCcellsize(1))
             IF (nx.EQ.MCMCcellNm(1)) nx=nx-1

             ny=INT(REAL(coord2,MK)/MCMCcellsize(2))
             IF (ny.EQ.MCMCcellNm(2)) ny=ny-1

             nz=INT(REAL(coord3,MK)/MCMCcellsize(3))
             IF (nz.EQ.MCMCcellNm(3)) nz=nz-1

             cbox=1+nx+MCMCcellNm(1)*(ny+nz*MCMCcellNm(2))

             tmpParticle=MCMCParticle(tmpHistoryParticle%candlabel,tmpHistoryParticle%proposal)

             IF (tmpHistoryParticle%wasadded) THEN
                Replaced=MCMCEraseFloatingParticle(coord1,coord2,coord3,cbox,tmpParticle,.FALSE.)
             ELSE
                info=MCMCInsertFloatingParticleCumulative(coord1,coord2,coord3,cbox,tmpParticle,.FALSE.)
                or_fail("MCMCInsertFloatingParticleCumulative",ppm_error=ppm_error_fatal)
             ENDIF

             seed => MCMCFloatingParticleInContainerHistory%prev()
          ENDDO

          !-------------------------------------------------------------------------
          ! recover the label image:
          !-------------------------------------------------------------------------
          seed => MCMCLabelImageHistory%last()
          DO WHILE (ASSOCIATED(seed))
             seedn => seed%first%getValue()

             labels_(seedn(1),seedn(2),seedn(3))=seedn(4)

             seed => MCMCLabelImageHistory%prev()
          ENDDO

          !-------------------------------------------------------------------------
          ! recover the statistics:
          !-------------------------------------------------------------------------
          seedlink => AppliedParticleOrigLabels%first
          DO WHILE (ASSOCIATED(seedlink))
             seedn => seedlink%getValue()

             !-------------------------------------------------------------------------
             ! coord1        =seedn(1)
             ! coord2        =seedn(2)
             ! coord3        =seedn(3)
             ! part%candlabel=seedn(4)
             ! OriginalLabel =seedn(5)
             !-------------------------------------------------------------------------
             CALL e_data%UpdateStatisticsWhenJump(image_,seedn,seedn(4),seedn(5),info)
             or_fail("e_data%UpdateStatisticsWhenJump",ppm_error=ppm_error_fatal)

             seedlink => seedlink%nextLink()
          ENDDO
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        END FUNCTION MCMCReject_3d

        LOGICAL FUNCTION MCMCMove(ppmrcdim,LocalMoves,info)
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_kind_double,ppm_kind_int64,ppm_error_fatal,ppm_error_error
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_error,ppm_err_sub_failed, &
          &   ppm_err_dealloc
          USE ppm_module_mesh_typedef, ONLY : ppm_t_subpatch_

          USE ppm_rc_module_global, ONLY : zerod,oned,mesh,image,labels,MCMCstepsize, &
          &   MCMCuseBiasedProposal,TotalMoves,MCMCtemperature
          USE ppm_rc_module_rnd, ONLY : ppm_rc_Saru_RPRNGD
          USE ppm_rc_module_linkedlist, ONLY : ppm_rc_link,ppm_rc_list
          USE ppm_rc_module_hash, ONLY : ppm_rc_MCMCParticlehtable
          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER,                 INTENT(IN   ) :: ppmrcdim
          !!! Problem dimension
          INTEGER(ppm_kind_int64), INTENT(IN   ) :: LocalMoves
          INTEGER,                 INTENT(  OUT) :: info
          !-----------------------------------------------------------------------
          !  Local variables
          !-----------------------------------------------------------------------
          TYPE(ppm_rc_MCMCParticlehtable) :: MCMCAppliedParticles

          TYPE(ppm_rc_list) :: MCMCAppliedParticleOrigLabels

          TYPE(ppm_rc_link), POINTER :: seedlnk

          TYPE(MCMCParticle) :: Part_A
          TYPE(MCMCParticle) :: ReverseFloatingP

          CLASS(ppm_t_subpatch_), POINTER :: sbpitr

          REAL(MK), CONTIGUOUS, DIMENSION(:,:),   POINTER     :: image_2d
          REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER     :: image_3d
          REAL(ppm_kind_double)                               :: t0
          REAL(ppm_kind_double)                               :: ForwardBackwardRatio
          REAL(ppm_kind_double)                               :: HastingsRatio
          REAL(ppm_kind_double)                               :: TotEnergyDiff
          REAL(MK)                                            :: etemp
          !!! MCMC Energy difference

          INTEGER(ppm_kind_int64)                            :: Iteration
          INTEGER, CONTIGUOUS, DIMENSION(:,:),   POINTER     :: labels_2d
          INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER     :: labels_3d
          INTEGER,             DIMENSION(:),     POINTER     :: Nm
          INTEGER,             DIMENSION(:),     POINTER     :: seedn
          INTEGER,             DIMENSION(ppmrcdim)           :: index_A
          INTEGER                                            :: PC
          INTEGER                                            :: cbox
          INTEGER                                            :: nx,ny,nz
          INTEGER                                            :: OriginalLabel

          LOGICAL :: e_merge
          LOGICAL :: Accept
          LOGICAL :: HardReject

          CHARACTER(LEN=ppm_char) :: caller='MCMCMoves'
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          e_merge      =.FALSE.
          HardReject   =.FALSE.
          TotEnergyDiff=zerod

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
                ! Iterate the candidates, calculate the energy and perform the moves.
                !-------------------------------------------------------------------------
                DO PC=1,MCMCstepsize
                   Part_A = MCMC_CandidateMove(PC)
                   !-------------------------------------------------------------------------
                   ! Apply particle A and B, start with B:
                   ! it is necessary that we start with particle B as we have
                   ! to calculate Q(A|B) and Qb(B|A).
                   !-------------------------------------------------------------------------
                   OriginalLabel=MCMC_LabelsBeforeJump_A(PC)

                   index_A=MCMC_CandidateMove_Index(:,PC)
                   !-------------------------------------------------------------------------
                   ! We calculate the energy and apply the move if
                   ! - the move has not been performed beforehand (a particle   was sampled twice)
                   ! - THE FOLLOWING IS CURRENTLY A NON ISSUE AS B CANNOT BE A':
                   !   particle B is not the reverse particle of particle A (in
                   !   case of m_usePairProposals, i.e. vN > 1). This is important
                   !   because the energy update gets corrupted as particle B is
                   !   not a valid particle before A was applied. To resolve this
                   !   we just perform a one particle move (we don't apply B).
                   !-------------------------------------------------------------------------
                   IF (.NOT.MCMCAppliedParticles%containselement(index_A(1),index_A(2),Part_A)) THEN
                      !-------------------------------------------------------------------------
                      ! Calculate the energy difference when changing this candidate:
                      !-------------------------------------------------------------------------
                      etemp=e_data%EvaluateEnergyDifference(image_2d,labels_2d,   &
                      &     index_A,OriginalLabel,Part_A%candlabel,e_merge)       &
                      &    +e_length%EvaluateEnergyDifference(image_2d,labels_2d, &
                      &     index_A,OriginalLabel,Part_A%candlabel,e_merge)

                      TotEnergyDiff=TotEnergyDiff+REAL(etemp,ppm_kind_double)

                      !-------------------------------------------------------------------------
                      ! Finally, perform the (particle-)move
                      !-------------------------------------------------------------------------
                      info=MCMCApplyParticle_2d(image_2d,labels_2d,index_A,Nm,Part_A,.FALSE.)
                      or_fail("MCMCApplyParticle!")

                      CALL MCMCAppliedParticles%insert(index_A(1),index_A(2),Part_A,info)
                      or_fail("MCMCAppliedParticles%insert")

                      CALL MCMCAppliedParticleOrigLabels%add(index_A(1),index_A(2),Part_A%candlabel,OriginalLabel)
                   ENDIF
                ENDDO !PC=1,MCMCstepsize

                CALL MCMCAppliedParticles%destroy(info)
                or_fail("MCMCAppliedParticles%destroy",ppm_error=ppm_error_fatal)

                DO PC=1,MCMCstepsize
                   !-------------------------------------------------------------------------
                   ! Correct the containers whenever floating particles were involved:
                   ! The method moveParticles, for simplicity, only works on the regular
                   ! particle set.
                   !-------------------------------------------------------------------------

                   !-------------------------------------------------------------------------
                   ! First, figure out if A' is floating:
                   ! Figure out if the backward particles are floating:
                   !-------------------------------------------------------------------------
                   OriginalLabel=MCMC_LabelsBeforeJump_A(PC)

                   index_A=MCMC_CandidateMove_Index(:,PC)
                   !-------------------------------------------------------------------------
                   ! if we're not in pair proposal mode we did not yet check if
                   ! A's reverse particle is floating (else we did already):
                   !-------------------------------------------------------------------------
                   MCMC_Particle_Ab_IsFloating(PC)=MCMCParticleHasFloatingProperty(labels_2d,index_A,OriginalLabel)

                   !-------------------------------------------------------------------------
                   ! the first condition is needed when not using pair proposal mode
                   !-------------------------------------------------------------------------
                   IF (MCMC_Particle_Ab_IsFloating(PC)) THEN
                      ReverseFloatingP=MCMCParticle(OriginalLabel,MCMCproposal(labels_2d,index_A))

                      !-------------------------------------------------------------------------
                      ! finally convert the regular particle into a floating particle,
                      ! i.e. insert it in the floating DS and remove it from the regular:
                      !-------------------------------------------------------------------------
                      !-------------------------------------------------------------------------
                      ! insert the reverse particle in the appropriate container. If
                      ! there is no space, we reject the move.
                      !-------------------------------------------------------------------------
                      nx=INT(REAL(index_A(1),MK)/MCMCcellsize(1))
                      IF (nx.EQ.MCMCcellNm(1)) nx=nx-1

                      ny=INT(REAL(index_A(2),MK)/MCMCcellsize(2))
                      IF (ny.EQ.MCMCcellNm(2)) ny=ny-1

                      cbox=1+nx+MCMCcellNm(1)*ny

                      IF (.NOT.MCMCInsertFloatingParticle(index_A(1),index_A(2),cbox,ReverseFloatingP,.TRUE.)) THEN
                         !-------------------------------------------------------------------------
                         ! TODO: calling MCMCReject from here invalidates the result.
                         !-------------------------------------------------------------------------
                         HardReject=.TRUE.
                      ENDIF
                   ENDIF !(MCMC_Particle_Ab_IsFloating(PC))
                ENDDO !PC=1,MCMCstepsize

                !-------------------------------------------------------------------------
                ! We are now in the state x'.
                ! Calculate Q'(A) and maybe Q'(B). Note that this has to be done after
                ! all particles were applied.
                !-------------------------------------------------------------------------
                DO PC=1,MCMCstepsize
                   !-------------------------------------------------------------------------
                   ! Calculate MCMC_qb_A
                   !-------------------------------------------------------------------------
                   IF (MCMCuseBiasedProposal) THEN
                      index_A=MCMC_CandidateMove_Index(:,PC)
                      MCMC_qb_A(PC)=REAL(MCMCproposal(labels_2d,index_A),ppm_kind_double)
                   ELSE
                      MCMC_qb_A(PC)=oned
                   ENDIF

                   !-------------------------------------------------------------------------
                   ! Normalize MCMC_qb_A
                   !-------------------------------------------------------------------------
                   MCMC_qb_A(PC)=MCMC_qb_A(PC)/(MCMCTotalNormalizer+MCMCTotalNormalizerlocal)
                ENDDO !PC=1,MCMCstepsize

                !-------------------------------------------------------------------------
                ! Calculate the forward-backward ratio:
                !-------------------------------------------------------------------------
                ForwardBackwardRatio=oned

                DO PC=1,MCMCstepsize
                   !-------------------------------------------------------------------------
                   ! Shrinkage and growth events have the same probability.
                   ! We hence only need to compare the individual particle ratios.
                   !-------------------------------------------------------------------------
                   ForwardBackwardRatio=ForwardBackwardRatio*MCMC_qb_A(PC)/MCMC_q_A(PC)
                ENDDO !PC=1,MCMCstepsize

                !-------------------------------------------------------------------------
                ! Compute the Hastingsratio:
                !-------------------------------------------------------------------------
                HastingsRatio = EXP(-TotEnergyDiff/MCMCtemperature)*ForwardBackwardRatio
                !-------------------------------------------------------------------------
                ! debug: check if the average of the FBR is equal to one.
                ! HastingsRatio = ForwardBackwardRatio
                ! Should I stay or shoud I go; the Metropolis-Hastings algorithm:
                !-------------------------------------------------------------------------
                IF (HastingsRatio.GE.oned) THEN
                   Accept=.TRUE.
                ELSE
                   IF (ppm_rc_Saru_RPRNGD().LT.HastingsRatio) THEN
                      Accept=.TRUE.
                   ELSE
                      Accept=.FALSE.
                   ENDIF
                ENDIF

                !-------------------------------------------------------------------------
                ! Register the result (if we accept) or rollback to the previous state.
                !-------------------------------------------------------------------------
                IF (Accept.AND..NOT.HardReject) THEN
                   Iteration=TotalMoves+LocalMoves

                   seedlnk => MCMCAppliedParticleOrigLabels%first
                   DO WHILE (ASSOCIATED(seedlnk))
                      seedn => seedlnk%getValue()

                      !-------------------------------------------------------------------------
                      ! store the results and finish (next iteration).
                      ! MCMCStore index,step,originallabel
                      !-------------------------------------------------------------------------
                      CALL MCMCResults%insert(seedn(1),seedn(2),Iteration,seedn(4),info)
                      or_fail("MCMCResults%insert",ppm_error=ppm_error_fatal)

                      !-------------------------------------------------------------------------
                      ! something changed, particles were inserted and deleted. We need
                      ! to update the edge-map constant.
                      !-------------------------------------------------------------------------
                      ! MCMCUpdateRegularParticleMapInNeighborhood(seedn(1),seedn(2));
                      !-------------------------------------------------------------------------
                      seedlnk => seedlnk%nextLink()
                   ENDDO
                ELSE
                   info=MCMCReject(image_2d,labels_2d,MCMCAppliedParticleOrigLabels)
                   or_fail("MCMCReject")
                ENDIF

                MCMCMove=Accept

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
                ! Iterate the candidates, calculate the energy and perform the moves.
                !-------------------------------------------------------------------------
                DO PC=1,MCMCstepsize
                   Part_A = MCMC_CandidateMove(PC)
                   !-------------------------------------------------------------------------
                   ! Apply particle A and B, start with B:
                   ! it is necessary that we start with particle B as we have
                   ! to calculate Q(A|B) and Qb(B|A).
                   !-------------------------------------------------------------------------
                   OriginalLabel=MCMC_LabelsBeforeJump_A(PC)

                   index_A=MCMC_CandidateMove_Index(:,PC)
                   !-------------------------------------------------------------------------
                   ! We calculate the energy and apply the move if
                   ! - the move has not been performed beforehand (a particle   was sampled twice)
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

                      TotEnergyDiff=TotEnergyDiff+REAL(etemp,ppm_kind_double)

                      !-------------------------------------------------------------------------
                      ! Finally, perform the (particle-)move
                      !-------------------------------------------------------------------------
                      info=MCMCApplyParticle_3d(image_3d,labels_3d,index_A,Nm,Part_A,.FALSE.)
                      or_fail("MCMCApplyParticle!")

                      CALL MCMCAppliedParticles%insert(index_A(1),index_A(2),index_A(3),Part_A,info)
                      or_fail("MCMCAppliedParticles%insert")

                      CALL MCMCAppliedParticleOrigLabels%add(index_A(1),index_A(2),index_A(3),Part_A%candlabel,OriginalLabel)
                   ENDIF
                ENDDO !PC=1,MCMCstepsize

                CALL MCMCAppliedParticles%destroy(info)
                or_fail("MCMCAppliedParticles%destroy",ppm_error=ppm_error_fatal)

                DO PC=1,MCMCstepsize
                   !-------------------------------------------------------------------------
                   ! Correct the containers whenever floating particles were involved:
                   ! The method moveParticles, for simplicity, only works on the regular
                   ! particle set.
                   !-------------------------------------------------------------------------

                   !-------------------------------------------------------------------------
                   ! First, figure out if A' is floating:
                   ! Figure out if the backward particles are floating:
                   !-------------------------------------------------------------------------
                   OriginalLabel=MCMC_LabelsBeforeJump_A(PC)

                   index_A=MCMC_CandidateMove_Index(:,PC)
                   !-------------------------------------------------------------------------
                   ! if we're not in pair proposal mode we did not yet check if
                   ! A's reverse particle is floating (else we did already):
                   !-------------------------------------------------------------------------
                   MCMC_Particle_Ab_IsFloating(PC)=MCMCParticleHasFloatingProperty(labels_3d,index_A,OriginalLabel)

                   !-------------------------------------------------------------------------
                   ! the first condition is needed when not using pair proposal mode
                   !-------------------------------------------------------------------------
                   IF (MCMC_Particle_Ab_IsFloating(PC)) THEN
                      ReverseFloatingP=MCMCParticle(OriginalLabel,MCMCproposal(labels_3d,index_A))

                      !-------------------------------------------------------------------------
                      ! finally convert the regular particle into a floating particle,
                      ! i.e. insert it in the floating DS and remove it from the regular:
                      !-------------------------------------------------------------------------
                      !-------------------------------------------------------------------------
                      ! insert the reverse particle in the appropriate container. If
                      ! there is no space, we reject the move.
                      !-------------------------------------------------------------------------
                      nx=INT(REAL(index_A(1),MK)/MCMCcellsize(1))
                      IF (nx.EQ.MCMCcellNm(1)) nx=nx-1

                      ny=INT(REAL(index_A(2),MK)/MCMCcellsize(2))
                      IF (ny.EQ.MCMCcellNm(2)) ny=ny-1

                      nz=INT(REAL(index_A(3),MK)/MCMCcellsize(3))
                      IF (nz.EQ.MCMCcellNm(3)) nz=nz-1

                      cbox=1+nx+MCMCcellNm(1)*(ny+nz*MCMCcellNm(2))

                      IF (.NOT.MCMCInsertFloatingParticle(index_A(1),index_A(2),index_A(3),cbox,ReverseFloatingP,.TRUE.)) THEN
                         !-------------------------------------------------------------------------
                         ! TODO: calling MCMCReject from here invalidates the result.
                         !-------------------------------------------------------------------------
                         HardReject=.TRUE.
                      ENDIF
                   ENDIF !(MCMC_Particle_Ab_IsFloating(PC))
                ENDDO !PC=1,MCMCstepsize

                !-------------------------------------------------------------------------
                ! We are now in the state x'.
                ! Calculate Q'(A) and maybe Q'(B). Note that this has to be done after
                ! all particles were applied.
                !-------------------------------------------------------------------------
                DO PC=1,MCMCstepsize
                   !-------------------------------------------------------------------------
                   ! Calculate MCMC_qb_A
                   !-------------------------------------------------------------------------
                   IF (MCMCuseBiasedProposal) THEN
                      index_A=MCMC_CandidateMove_Index(:,PC)
                      MCMC_qb_A(PC)=REAL(MCMCproposal(labels_3d,index_A),ppm_kind_double)
                   ELSE
                      MCMC_qb_A(PC)=oned
                   ENDIF

                   !-------------------------------------------------------------------------
                   ! Normalize MCMC_qb_A
                   !-------------------------------------------------------------------------
                   MCMC_qb_A(PC)=MCMC_qb_A(PC)/(MCMCTotalNormalizer+MCMCTotalNormalizerlocal)
                ENDDO !PC=1,MCMCstepsize

                !-------------------------------------------------------------------------
                ! Calculate the forward-backward ratio:
                !-------------------------------------------------------------------------
                ForwardBackwardRatio=oned

                DO PC=1,MCMCstepsize
                   !-------------------------------------------------------------------------
                   ! Shrinkage and growth events have the same probability.
                   ! We hence only need to compare the individual particle ratios.
                   !-------------------------------------------------------------------------
                   ForwardBackwardRatio=ForwardBackwardRatio*MCMC_qb_A(PC)/MCMC_q_A(PC)
                ENDDO !PC=1,MCMCstepsize

                !-------------------------------------------------------------------------
                ! Compute the Hastingsratio:
                !-------------------------------------------------------------------------
                HastingsRatio = EXP(-TotEnergyDiff/MCMCtemperature)*ForwardBackwardRatio
                !-------------------------------------------------------------------------
                ! debug: check if the average of the FBR is equal to one.
                ! HastingsRatio = ForwardBackwardRatio
                ! Should I stay or shoud I go; the Metropolis-Hastings algorithm:
                !-------------------------------------------------------------------------
                IF (HastingsRatio.GE.oned) THEN
                   Accept=.TRUE.
                ELSE
                   IF (ppm_rc_Saru_RPRNGD().LT.HastingsRatio) THEN
                      Accept=.TRUE.
                   ELSE
                      Accept=.FALSE.
                   ENDIF
                ENDIF

                !-------------------------------------------------------------------------
                ! Register the result (if we accept) or rollback to the previous state.
                !-------------------------------------------------------------------------
                IF (Accept.AND..NOT.HardReject) THEN
                   Iteration=TotalMoves+LocalMoves

                   seedlnk => MCMCAppliedParticleOrigLabels%first
                   DO WHILE (ASSOCIATED(seedlnk))
                      seedn => seedlnk%getValue()

                      !-------------------------------------------------------------------------
                      ! store the results and finish (next iteration).
                      ! MCMCStore index,step,originallabel
                      !-------------------------------------------------------------------------
                      CALL MCMCResults%insert(seedn(1),seedn(2),seedn(3),Iteration,seedn(5),info)
                      or_fail("MCMCResults%insert",ppm_error=ppm_error_fatal)

                      !-------------------------------------------------------------------------
                      ! something changed, particles were inserted and deleted. We need
                      ! to update the edge-map constant.
                      !-------------------------------------------------------------------------
                      ! MCMCUpdateRegularParticleMapInNeighborhood(seedn(1),seedn(2));
                      !-------------------------------------------------------------------------
                      seedlnk => seedlnk%nextLink()
                   ENDDO
                ELSE
                   info=MCMCReject(image_3d,labels_3d,MCMCAppliedParticleOrigLabels)
                   or_fail("MCMCReject")
                ENDIF

                MCMCMove=Accept

                sbpitr => mesh%subpatch%next()
             ENDDO
             NULLIFY(labels_3d,image_3d)
          END SELECT

          CALL MCMCAppliedParticleOrigLabels%destroy()
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        END FUNCTION MCMCMove











