      !-------------------------------------------------------------------------
      !  Subroutine   :                    ppm_rc_finalize
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
      !  Subroutine   :                    ppm_rc_finalize
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Clean up global memory and close open units. This
      !                 also finalizes ppm and mpi.
      !
      !  Input        : errstat       (I) error status. 0 if this is a
      !                                   normal termination. .NE. 0 if
      !                                   this is a program abort.
      !
      !  Input/output :
      !
      !  Output       : info          (I) return status. 0 on success.
      !
      !  Routines     : substart
      !                 substop
      !                 ppm_finalize
      !                 MPI_Finalize
      !
      !  Remarks      :
      !
      !  References   :
      !-------------------------------------------------------------------------
      SUBROUTINE ppm_rc_finalize(errstat,info)

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_finalize, ONLY : ppm_finalize
      USE ppm_module_mpi

      USE ppm_rc_module_linkedlist
      USE ppm_rc_module_topologicalnumber, ONLY : TopologicalNumberFunction, &
      &   FG_ConnectivityType
      USE ppm_rc_module_energy, ONLY : e_data,e_length,E_ContourLengthApprox
      USE ppm_rc_module_rnd, ONLY : ppm_rc_DestroyImageDiscrDistr, &
      &   ppm_rc_DestroyParticlesDiscrDistr
      USE ppm_rc_module_mcmc, ONLY : MCMClengthProposalMask,            &
      &   MCMCAllParticlesFwdProposals,MCMCchildren,MCMCparents,        &
      &   MCMC_CandidateMove,MCMC_CandidateMove_Index,MCMC_PartnerMove, &
      &   MCMC_qb_B,MCMC_q_A,MCMC_q_B,MCMC_q_A_B,MCMC_q_B_A,            &
      &   MCMC_qb_A,MCMC_qb_A_B,MCMC_qb_B_A,MCMC_LabelsBeforeJump_A,    &
      &   MCMC_LabelsBeforeJump_B,MCMC_Particle_Ab_IsFloating,          &
      &   MCMC_Particle_Bb_IsFloating,MCMC_Particle_A_IsFloating,       &
      &   MCMC_PartnerMove_Index
      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      INTEGER, INTENT(IN   ) :: errstat
      INTEGER, INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      CLASS(ppm_t_subpatch_), POINTER :: sbpitr

      REAL(ppm_kind_double) :: t0

      INTEGER :: i,ipatch

      CHARACTER(LEN=ppm_char) :: caller = 'ppm_rc_finalize'
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialize
      !-------------------------------------------------------------------------
      CALL substart(caller,t0,info)

      !-------------------------------------------------------------------------
      !  Destroy the hash table for the region labels
      !-------------------------------------------------------------------------
      IF (ASSOCIATED(htable)) THEN
         CALL htable%destroy(info)
         or_fail("htable%destroy")
         DEALLOCATE(htable,STAT=info)
         or_fail_dealloc("htable")
         NULLIFY(htable)
      ENDIF

      !-------------------------------------------------------------------------
      !  Deallocate Fields
      !-------------------------------------------------------------------------
      IF (ASSOCIATED(image)) THEN
         CALL image%destroy(info)
         or_fail("image%destroy")
         DEALLOCATE(image,STAT=info)
         or_fail_dealloc("image")
         NULLIFY(image)
      ENDIF
      IF (ASSOCIATED(labels)) THEN
         CALL labels%destroy(info)
         or_fail("labels%destroy")
         DEALLOCATE(labels,STAT=info)
         or_fail_dealloc("labels")
         NULLIFY(labels)
      ENDIF
      IF (ASSOCIATED(pind)) THEN
         CALL pind%destroy(info)
         or_fail("pind%destroy")
         DEALLOCATE(pind,STAT=info)
         or_fail_dealloc("pind")
         NULLIFY(pind)
      ENDIF
      IF (ASSOCIATED(plabels)) THEN
         DEALLOCATE(plabels,STAT=info)
         or_fail_dealloc("plabels")
         NULLIFY(plabels)
      ENDIF
      IF (ASSOCIATED(Backuplabels)) THEN
         CALL Backuplabels%destroy(info)
         or_fail("Backuplabels%destroy")
         DEALLOCATE(Backuplabels,STAT=info)
         or_fail_dealloc("Backuplabels")
         NULLIFY(Backuplabels)
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
      IF (ALLOCATED(energya)) THEN
         DEALLOCATE(energya,STAT=info)
         or_fail_dealloc("energya")
      ENDIF

      IF (ASSOCIATED(e_data)) THEN
         CALL e_data%destroy(info)
         or_fail_dealloc("e_data%destroy")
         DEALLOCATE(e_data,STAT=info)
         or_fail_dealloc("e_data")
         NULLIFY(e_data)
      ENDIF

      IF (ASSOCIATED(e_length)) THEN
         SELECT TYPE (e_length)
         TYPE IS (E_ContourLengthApprox)
            CALL e_length%destroy(info)
            or_fail("e_length%destroy")
         END SELECT
         DEALLOCATE(e_length,STAT=info)
         or_fail_dealloc("e_length")
         NULLIFY(e_length)
      ENDIF

      !-------------------------------------------------------------------------
      !  Topology related members
      !-------------------------------------------------------------------------
      CALL TopologicalNumberFunction%destroy(info)
      or_fail("TopologicalNumberFunction%destroy")

      IF (ASSOCIATED(FG_ConnectivityType)) THEN
         CALL FG_ConnectivityType%destroy(info)
         or_fail("FG_ConnectivityType%destroy")
         DEALLOCATE(FG_ConnectivityType,STAT=info)
         or_fail_dealloc("FG_ConnectivityType")
         NULLIFY(FG_ConnectivityType)
      ENDIF
      !-------------------------------------------------------------------------
      !  Deallocate Domain decomposition variables
      !-------------------------------------------------------------------------
      IF (ASSOCIATED(proc_speed)) THEN
         CALL ppm_alloc(proc_speed,(/1/),ppm_param_dealloc,info)
         or_fail_dealloc("proc_speed")
         NULLIFY(proc_speed)
      ENDIF
      IF (ASSOCIATED(cost)) THEN
         CALL ppm_alloc(cost,(/1/),ppm_param_dealloc,info)
         or_fail_dealloc("cost")
         NULLIFY(cost)
      ENDIF
      IF (ASSOCIATED(min_phys)) THEN
         DEALLOCATE(min_phys,STAT=info)
         or_fail_dealloc("min_phys")
         NULLIFY(min_phys)
      ENDIF
      IF (ASSOCIATED(max_phys)) THEN
         DEALLOCATE(max_phys,STAT=info)
         or_fail_dealloc("max_phys")
         NULLIFY(max_phys)
      ENDIF
      IF (ALLOCATED(ghostsize)) THEN
         DEALLOCATE(ghostsize,STAT=info)
         or_fail_dealloc("ghostsize")
      ENDIF
      IF (ALLOCATED(ioghostsize)) THEN
         DEALLOCATE(ioghostsize,STAT=info)
         or_fail_dealloc("ioghostsize")
      ENDIF
      IF (ALLOCATED(inighostsize)) THEN
         DEALLOCATE(inighostsize,STAT=info)
         or_fail_dealloc("inighostsize")
      ENDIF
      IF (ALLOCATED(ghostsize_equil)) THEN
         DEALLOCATE(ghostsize_equil,STAT=info)
         or_fail_dealloc("ghostsize_equil")
      ENDIF
      IF (ALLOCATED(ghostsize_run)) THEN
         DEALLOCATE(ghostsize_run,STAT=info)
         or_fail_dealloc("ghostsize_run")
      ENDIF
      IF (ALLOCATED(bcdef)) THEN
         DEALLOCATE(bcdef,STAT=info)
         or_fail_dealloc("bcdef")
      ENDIF
      IF (ALLOCATED(ineighproc)) THEN
         DEALLOCATE(ineighproc,STAT=info)
         or_fail_dealloc("ineighproc")
      ENDIF
      IF (ALLOCATED(procflag)) THEN
         DEALLOCATE(procflag,STAT=info)
         or_fail_dealloc("procflag")
      ENDIF
      IF (ALLOCATED(ineighproc)) THEN
         DEALLOCATE(ineighproc,STAT=info)
         or_fail_dealloc("ineighproc")
      ENDIF

      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      IF (ALLOCATED(init_rd)) THEN
         DEALLOCATE(init_rd,STAT=info)
         or_fail_dealloc("init_rd")
      ENDIF
      IF (ALLOCATED(init_sp)) THEN
         DEALLOCATE(init_sp,STAT=info)
         or_fail_dealloc("init_sp")
      ENDIF

      !----------------------------------------------------------------------
      !  Clear the old candidate set
      !----------------------------------------------------------------------
      IF (ASSOCIATED(mesh)) THEN
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
            IF (ALLOCATED(MCMCParticleInContainerHistory)) THEN
               CALL MCMCParticleInContainerHistory(ipatch)%destroy(info)
               or_fail("MCMCParticleInContainerHistory(ipatch)%destroy")
            ENDIF
            IF (ALLOCATED(MCMCFloatingParticleInContainerHistory)) THEN
               CALL MCMCFloatingParticleInContainerHistory(ipatch)%destroy(info)
               or_fail("MCMCFloatingParticleInContainerHistory(ipatch)%destroy")
            ENDIF
            IF (ALLOCATED(MCMCLabelImageHistory)) THEN
               CALL MCMCLabelImageHistory(ipatch)%destroy(info)
               or_fail("MCMCLabelImageHistory(ipatch)%destroy")
            ENDIF
            IF (ALLOCATED(MCMCAppliedParticleOrigLabels)) THEN
               !!! This is a list structure not a collection
               !!! ppm_rc_list
               CALL MCMCAppliedParticleOrigLabels(ipatch)%destroy()
            ENDIF
            sbpitr => mesh%subpatch%next()
            ipatch=ipatch+1
         ENDDO
      ENDIF

      IF (ALLOCATED(ppm_rc_seeds)) THEN
         DEALLOCATE(ppm_rc_seeds,STAT=info)
         or_fail_dealloc("ppm_rc_seeds")
      ENDIF
      IF (ALLOCATED(Candidates)) THEN
         DEALLOCATE(Candidates,STAT=info)
         or_fail_dealloc("Candidates")
      ENDIF
      IF (ALLOCATED(InnerContourContainer)) THEN
         DEALLOCATE(InnerContourContainer,STAT=info)
         or_fail_dealloc("InnerContourContainer")
      ENDIF
      IF (ALLOCATED(CompetingRegions)) THEN
         DEALLOCATE(CompetingRegions,STAT=info)
         or_fail_dealloc("CompetingRegions")
      ENDIF
      IF (ALLOCATED(m_Seeds)) THEN
         DEALLOCATE(m_Seeds,STAT=info)
         or_fail_dealloc("m_Seeds")
      ENDIF
      IF (ALLOCATED(Candidates_list)) THEN
         DEALLOCATE(Candidates_list,STAT=info)
         or_fail_dealloc("Candidates_list")
      ENDIF
      IF (ALLOCATED(ppm_rc_seeds_to_remove)) THEN
         DEALLOCATE(ppm_rc_seeds_to_remove,STAT=info)
         or_fail_dealloc("ppm_rc_seeds_to_remove")
      ENDIF
      IF (ALLOCATED(MCMCParticleInContainerHistory)) THEN
         DEALLOCATE(MCMCParticleInContainerHistory,STAT=info)
         or_fail_dealloc("MCMCParticleInContainerHistory")
      ENDIF
      IF (ALLOCATED(MCMCFloatingParticleInContainerHistory)) THEN
         DEALLOCATE(MCMCFloatingParticleInContainerHistory,STAT=info)
         or_fail_dealloc("MCMCFloatingParticleInContainerHistory")
      ENDIF
      IF (ALLOCATED(MCMCLabelImageHistory)) THEN
         DEALLOCATE(MCMCLabelImageHistory,STAT=info)
         or_fail_dealloc("MCMCLabelImageHistory")
      ENDIF
      IF (ALLOCATED(MCMCAppliedParticleOrigLabels)) THEN
         DEALLOCATE(MCMCAppliedParticleOrigLabels,STAT=info)
         or_fail_dealloc("MCMCAppliedParticleOrigLabels")
      ENDIF
      !---------------------------------------------------------------------
      ! Destroy MCMC memory
      !---------------------------------------------------------------------
      SELECT CASE (UseMCMC)
      CASE (.TRUE.)
         IF (ALLOCATED(MCMClengthProposalMask)) THEN
            DEALLOCATE(MCMClengthProposalMask,STAT=info)
            or_fail_dealloc("MCMClengthProposalMask")
         ENDIF
         IF (ALLOCATED(MCMCAllParticlesFwdProposals)) THEN
            DEALLOCATE(MCMCAllParticlesFwdProposals,STAT=info)
            or_fail_dealloc("MCMCAllParticlesFwdProposals")
         ENDIF

         info=ppm_rc_DestroyImageDiscrDistr()
         or_fail("ppm_rc_DestroyImageDiscrDistr")

         info=ppm_rc_DestroyParticlesDiscrDistr()
         or_fail("ppm_rc_DestroyParticlesDiscrDistr")

         IF (ALLOCATED(MCMCparents)) THEN
            DO i=0,SIZE(MCMCparents)-1
               CALL MCMCparents(i)%destroy(info)
               or_fail("MCMCparents(i)%destroy")
            ENDDO
            DEALLOCATE(MCMCparents,STAT=info)
            or_fail_dealloc("MCMCparents(i)%destroy")
         ENDIF
         IF (ALLOCATED(MCMCchildren)) THEN
            DO i=0,SIZE(MCMCchildren)-1
               CALL MCMCchildren(i)%destroy(info)
               or_fail("MCMCchildren(i)%destroy")
            ENDDO
            DEALLOCATE(MCMCchildren,STAT=info)
            or_fail_dealloc("MCMCchildren(i)%destroy")
         ENDIF
         IF (ALLOCATED(MCMC_CandidateMove)) THEN
            DEALLOCATE(MCMC_CandidateMove,STAT=info)
            or_fail_dealloc("MCMCCandidateMove")
         ENDIF
         IF (ALLOCATED(MCMC_CandidateMove_Index)) THEN
            DEALLOCATE(MCMC_CandidateMove_Index,STAT=info)
            or_fail_dealloc("MCMC_CandidateMove_Index")
         ENDIF
         IF (ALLOCATED(MCMC_PartnerMove_Index)) THEN
            DEALLOCATE(MCMC_PartnerMove_Index,STAT=info)
            or_fail_dealloc("MCMC_PartnerMove_Index")
         ENDIF
         IF (ALLOCATED(MCMC_PartnerMove)) THEN
            DEALLOCATE(MCMC_PartnerMove,STAT=info)
            or_fail_dealloc("MCMCPartnerMove")
         ENDIF
         IF (ALLOCATED(MCMC_q_A)) THEN
            DEALLOCATE(MCMC_q_A,STAT=info)
            or_fail_dealloc("MCMC_q_A")
         ENDIF
         IF (ALLOCATED(MCMC_q_B)) THEN
            DEALLOCATE(MCMC_q_B,STAT=info)
            or_fail_dealloc("MCMC_q_B")
         ENDIF
         IF (ALLOCATED(MCMC_q_A_B)) THEN
            DEALLOCATE(MCMC_q_A_B,STAT=info)
            or_fail_dealloc("MCMC_q_A_B")
         ENDIF
         IF (ALLOCATED(MCMC_q_B_A)) THEN
            DEALLOCATE(MCMC_q_B_A,STAT=info)
            or_fail_dealloc("MCMC_q_B_A")
         ENDIF
         IF (ALLOCATED(MCMC_qb_A)) THEN
            DEALLOCATE(MCMC_qb_A,STAT=info)
            or_fail_dealloc("MCMC_qb_A")
         ENDIF
         IF (ALLOCATED(MCMC_qb_B)) THEN
            DEALLOCATE(MCMC_qb_B,STAT=info)
            or_fail_dealloc("MCMC_qb_B")
         ENDIF
         IF (ALLOCATED(MCMC_qb_A_B)) THEN
            DEALLOCATE(MCMC_qb_A_B,STAT=info)
            or_fail_dealloc("MCMC_qb_A_B")
         ENDIF
         IF (ALLOCATED(MCMC_qb_B_A)) THEN
            DEALLOCATE(MCMC_qb_B_A,STAT=info)
            or_fail_dealloc("MCMC_qb_B_A")
         ENDIF
         IF (ALLOCATED(MCMC_LabelsBeforeJump_A)) THEN
            DEALLOCATE(MCMC_LabelsBeforeJump_A,STAT=info)
            or_fail_dealloc("MCMC_LabelsBeforeJump_A")
         ENDIF
         IF (ALLOCATED(MCMC_LabelsBeforeJump_B)) THEN
            DEALLOCATE(MCMC_LabelsBeforeJump_B,STAT=info)
            or_fail_dealloc("MCMC_LabelsBeforeJump_B")
         ENDIF
         IF (ALLOCATED(MCMC_Particle_Ab_IsFloating)) THEN
            DEALLOCATE(MCMC_Particle_Ab_IsFloating,STAT=info)
            or_fail_dealloc("MCMC_Particle_Ab_IsFloating")
         ENDIF
         IF (ALLOCATED(MCMC_Particle_Bb_IsFloating)) THEN
            DEALLOCATE(MCMC_Particle_Bb_IsFloating,STAT=info)
            or_fail_dealloc("MCMC_Particle_Bb_IsFloating")
         ENDIF
         IF (ALLOCATED(MCMC_Particle_A_IsFloating)) THEN
            DEALLOCATE(MCMC_Particle_A_IsFloating,STAT=info)
            or_fail_dealloc("MCMC_Particle_A_IsFloating")
         ENDIF
      END SELECT

      !---------------------------------------------------------------------
      ! Destroy the part
      !---------------------------------------------------------------------
      IF (ASSOCIATED(Part)) THEN
         CALL Part%destroy(info)
         or_fail("Failed to destroy Part.")

         dealloc_pointer("Part")
      ENDIF

      !---------------------------------------------------------------------
      ! Destroy the mesh
      !---------------------------------------------------------------------
      IF (ASSOCIATED(mesh)) THEN
         CALL mesh%destroy(info)
         or_fail('Failed to destroy mesh.')

         dealloc_pointer("mesh")
      ENDIF

      !-------------------------------------------------------------------------
      !  Finalize PPM
      !-------------------------------------------------------------------------
      CALL ppm_finalize(info)

      SELECT CASE (errstat)
      CASE (0)
#ifdef __MPI
         CALL MPI_Barrier(comm,info)
         or_fail_MPI("MPI_Barrier")

         !---------------------------------------------------------------------
         !  Finalize MPI
         !---------------------------------------------------------------------
         CALL MPI_Finalize(info)
         or_fail_MPI('FAILED TO FINALIZE MPI. BAILING OUT!')

      CASE DEFAULT
         !---------------------------------------------------------------------
         !  Abort MPI
         !---------------------------------------------------------------------
         CALL MPI_Abort(comm,errstat,info)
         or_fail_MPI('FAILED TO ABORT MPI!')
#endif
      END SELECT

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      END SUBROUTINE ppm_rc_finalize
