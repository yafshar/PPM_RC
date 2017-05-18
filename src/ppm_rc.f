      !-------------------------------------------------------------------------
      !  Program      :                    ppm_rc
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
      !  Author           - y.afshar           Feb    2016
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  Please do cite:
      !
      !  Y. Afshar, and I. F. Sbalzarini. A Parallel Distributed-Memory Particle
      !  Method Enables Acquisition-Rate Segmentation of Large Fluorescence
      !  Microscopy Images. PLoS ONE 11(4):e0152528, (2016).
      !
      !  when publishing research data obtained using PPM_RC
      !-------------------------------------------------------------------------

      PROGRAM ppm_rc
        !----------------------------------------------------------------------
        !  Modules
        !----------------------------------------------------------------------
        USE ppm_module_mpi

        USE ppm_rc_module_global
        USE ppm_rc_module_linkedlist, ONLY : Candidates,m_Seeds, &
        &   CompetingRegions
        USE ppm_rc_module_init, ONLY : ppm_rc_init_arg,ppm_rc_init_seg_2d, &
        &   ppm_rc_init_seg_3d,ppm_rc_init_mcmc_2d,ppm_rc_init_mcmc_3d
        USE ppm_rc_module_energy, ONLY : ppm_rc_energy_compute_2d, &
        &   ppm_rc_energy_compute_3d,e_data,ppm_rc_energy_parameter_redefine, &
        &   ppm_rc_energy_parameter_redefine_mcmc
        USE ppm_rc_module_move
        USE ppm_rc_module_finalize, ONLY : ppm_rc_finalize
        USE ppm_rc_module_write
#ifdef __Linux
        USE ppm_rc_module_util, ONLY : ppm_rc_mem_usage
#endif
        USE ppm_rc_module_hash, ONLY : ppm_rc_MCMCParticlehtable
        USE ppm_rc_module_rnd, ONLY : ppm_rc_Saru_IPRNG,ppm_rc_Saru_RPRNG, &
        &   ppm_rc_Saru_RPRNGD,ppm_rc_GetPartDistrIndex
        USE ppm_rc_module_mcmc
        IMPLICIT NONE

        !-----------------------------------------------------------------------
        !  Includes
        !-----------------------------------------------------------------------
        !-----------------------------------------------------------------------
        !  Local variables
        !-----------------------------------------------------------------------
        CLASS(ppm_t_subpatch_), POINTER :: sbpitr

        TYPE(MCMCParticle) :: tmpParticle

        ! variables for timing
        REAL(ppm_kind_double) :: t0,tinit,ttime(11)
        REAL(ppm_kind_double) :: tenergy, tcontour, tfilter1,tfilter2
        REAL(ppm_kind_double) :: tenergys,tcontours,tfilters
        REAL(ppm_kind_double) :: tosc, tmove, timg
        REAL(ppm_kind_double) :: toscs,tmoves,timgs
        REAL(ppm_kind_double) :: rnd

        INTEGER(ppm_kind_int64)                            :: LocalMoves,LocalAcceptedMoves
#ifdef __MPI
        INTEGER(ppm_kind_int64), DIMENSION(2)              :: TmpMoves
        INTEGER                                            :: request1,request2
#endif
        INTEGER                                            :: info
#ifdef __Linux
        INTEGER                                            :: memory
#endif
        INTEGER                                            :: i,j,ipatch
        INTEGER                                            :: RLabel,label_region,PC,ParticleIndex
        INTEGER,                 DIMENSION(:), ALLOCATABLE :: index_


        CHARACTER(LEN=ppm_char) :: caller='ppm_rc'

        LOGICAL :: OffBoundarySampling
        LOGICAL :: Growth
        LOGICAL :: done
        LOGICAL :: Convergence
        LOGICAL :: Accepted
        LOGICAL :: ParticleAIsFloating
        LOGICAL :: rejected

        !-------------------------------------------------------------------------
        !  Initializing
        !-------------------------------------------------------------------------
        tenergys  = zerod
        tcontours = zerod
        tfilters  = zerod
        toscs     = zerod
        tmoves    = zerod
        timgs     = zerod

        !-------------------------------------------------------------------------
        ! Initialize the ppm library and ppm_rc arguments & variables
        !-------------------------------------------------------------------------
        CALL ppm_rc_init_arg(info)
        SELECT CASE (info)
        CASE (0)
        CASE (exit_gracefully)
           GOTO 7000
        CASE DEFAULT
           or_fail("ppm_rc_init_arg",ppm_error=ppm_error_fatal)
        END SELECT

        CALL ppm_util_time(t0)
        !-------------------------------------------------------------------------
        ! Initialization starts
        !-------------------------------------------------------------------------
        SELECT CASE (ppm_rc_dim)
        CASE (2)
           SELECT CASE (UseMCMC)
           CASE (.TRUE.)
              SELECT CASE (MCMCcontinue)
              CASE (.TRUE.)
                 !-------------------------------------------------------------------------
                 ! We continue MCMC after Segmentation is done by RC
                 !-------------------------------------------------------------------------
                 CALL ppm_rc_init_seg_2d(info)
              END SELECT
           CASE (.FALSE.)
              !-------------------------------------------------------------------------
              ! The segmentation algorithm
              !-------------------------------------------------------------------------
              CALL ppm_rc_init_seg_2d(info)
           END SELECT
        CASE (3)
           SELECT CASE (UseMCMC)
           CASE (.TRUE.)
              SELECT CASE (MCMCcontinue)
              CASE (.TRUE.)
                 !-------------------------------------------------------------------------
                 ! We continue MCMC after Segmentation is done by RC
                 !-------------------------------------------------------------------------
                 CALL ppm_rc_init_seg_3d(info)
              END SELECT
           CASE (.FALSE.)
              !-------------------------------------------------------------------------
              ! The segmentation algorithm
              !-------------------------------------------------------------------------
              CALL ppm_rc_init_seg_3d(info)
           END SELECT
        END SELECT
        SELECT CASE (info)
        CASE (0)
        CASE (exit_gracefully)
           GOTO 7000
        CASE DEFAULT
           or_fail('Failed to initialize!',exit_point=8000)
        END SELECT

        ! This is the DRC segmentation, which can be followed by MCMC sampling
        IF (.NOT.UseMCMC.OR.(UseMCMC.AND.MCMCcontinue)) THEN
           ! TODO remove this from here
           ! This is a temporary hack for producing some output
           debug=2
           CALL ppm_util_time(tinit)

           ! Make sure initialization is done on every processor
           IF (ppm_rc_debug.GT.0) THEN
              ! If we are on the debug mode every prcessor should print the Initialization time
              stdout("Initialization took: ",'tinit-t0'," secs.")
           ELSE
              IF (rank.EQ.0) THEN
                 stdout("Initialization took: ",'tinit-t0'," secs.")
              ENDIF
           ENDIF

           !-------------------------------------------------------------------------
           !!! RC Segmentation
           !-------------------------------------------------------------------------
           done=istep.GE.maxiter
           Convergence=.FALSE.

           !-------------------------------------------------------------------------
           !  Do iterations
           !-------------------------------------------------------------------------
           RC_itrloop: DO WHILE (.NOT.done)
              istep=istep+1

              IF (rank.EQ.0) THEN
                 stdout("istep=",istep)
              ENDIF

              IF (istep.EQ.nsteql) THEN
                 CALL ppm_rc_energy_parameter_redefine(info)
                 or_fail("ppm_rc_energy_parameter_redefine")
              ENDIF

              !----------------------------------------------------------------------
              !  Clear the old candidate set
              !----------------------------------------------------------------------
              sbpitr => mesh%subpatch%begin()
              ipatch=1
              DO WHILE (ASSOCIATED(sbpitr))
                 CALL Candidates(ipatch)%destroy(info)
                 or_fail("Candidates(ipatch)%destroy")

                 CALL m_Seeds(ipatch)%destroy(info)
                 or_fail("m_Seeds(ipatch)%destroy")

                 CALL CompetingRegions(ipatch)%destroy(info)
                 or_fail("CompetingRegions(ipatch)%destroy")

                 sbpitr => mesh%subpatch%next()
                 ipatch=ipatch+1
              ENDDO

              CALL ppm_util_time(t0)

              SELECT CASE (ppm_rc_dim)
              CASE (2)
                 !-----------------------------------------------------------------
                 !  Compute the energies and build the mother/daughter lists
                 !-----------------------------------------------------------------
                 CALL ppm_rc_energy_compute_2d(info)
                 or_fail("ppm_rc_energy_compute")
                 CALL ppm_util_time(tenergy)
#ifdef __Linux
                 IF (debug.GT.1) THEN
                    CALL ppm_rc_mem_usage(memory,info)
                    stdout("mem_usage at ppm_rc_energy_compute=",memory)
                 ENDIF
#endif
                 !-----------------------------------------------------------------
                 !  Find topologically compatible candidates and store their
                 !  indices as an accepted or not
                 !-----------------------------------------------------------------
                 CALL ppm_rc_contour_propagation_2d(info)
                 or_fail("ppm_rc_contour_propagation")
                 CALL ppm_util_time(tcontour)
#ifdef __Linux
                 IF (debug.GT.1) THEN
                    CALL ppm_rc_mem_usage(memory,info)
                    stdout("mem_usage at ppm_rc_contour_propagation=",memory)
                 ENDIF
#endif
                 !-----------------------------------------------------------------
                 !  Filtering: Filter unaccepted particles
                 !-----------------------------------------------------------------
                 CALL ppm_rc_FilterCandidates_2d(info)
                 or_fail("ppm_rc_FilterCandidates")
                 CALL ppm_util_time(tfilter1)
#ifdef __Linux
                 IF (debug.GT.1) THEN
                    CALL ppm_rc_mem_usage(memory,info)
                    stdout("mem_usage at ppm_rc_FilterCandidates=",memory)
                 ENDIF
#endif

                 !-----------------------------------------------------------------
                 ! Oscillations: Detecting Oscillations and reduce the moves.
                 !-----------------------------------------------------------------
                 CALL ppm_rc_DetectOscillations2(info)
                 or_fail("ppm_rc_DetectOscillations")
                 CALL ppm_util_time(tosc)
#ifdef __Linux
                 IF (debug.GT.1) THEN
                    CALL ppm_rc_mem_usage(memory,info)
                    stdout("mem_usage at ppm_rc_DetectOscillations=",memory)
                 ENDIF
#endif
                 !-----------------------------------------------------------------
                 ! Filtering: Intermediate step: filter the candidates according
                 ! to their rank and spacial position.
                 !-----------------------------------------------------------------
                 CALL ppm_rc_FilterCandidatesContainerUsingRanks_2d(info)
                 or_fail("ppm_rc_FilterCandidatesContainerUsingRanks")
                 CALL ppm_util_time(tfilter2)
#ifdef __Linux
                 IF (debug.GT.1) THEN
                    CALL ppm_rc_mem_usage(memory,info)
                    stdout("mem_usage at ppm_rc_FilterCandidatesContainerUsingRanks=",memory)
                 ENDIF
#endif
                 !-----------------------------------------------------------------
                 ! Move: Move all the points that are simple. Non simple points
                 ! remain in the candidates list.
                 !-----------------------------------------------------------------
                 Convergence=ppm_rc_move_2d(info)
                 or_fail("ppm_rc_move")
                 CALL ppm_util_time(tmove)
#ifdef __Linux
                 IF (debug.GT.1) THEN
                    CALL ppm_rc_mem_usage(memory,info)
                    stdout("mem_usage at ppm_rc_move=",memory)
                 ENDIF
#endif

                 !-----------------------------------------------------------------
                 ! Write the output.
                 !-----------------------------------------------------------------
                 IF (MOD(istep,freqoutput).EQ.0) THEN
                    CALL ppm_rc_write_image_label_2d(labels,mesh,outputfile,info) !bitsPerSampleW=16
                    or_fail("ppm_rc_write_image")
!                      call mesh%print_vtk(outputname,info,Field=labels)
                    CALL ppm_util_time(timg)
                    timgs=timgs+timg-tmove
                 ENDIF

              CASE (3)
                 !-----------------------------------------------------------------
                 !  Compute the energies and build the mother/daughter lists
                 !-----------------------------------------------------------------
                 CALL ppm_rc_energy_compute_3d(info)
                 or_fail("ppm_rc_energy_compute")
                 CALL ppm_util_time(tenergy)
#ifdef __Linux
                 IF (debug.GT.1) THEN
                    CALL ppm_rc_mem_usage(memory,info)
                    stdout("mem_usage at ppm_rc_energy_compute=",memory)
                 ENDIF
#endif
                 !-----------------------------------------------------------------
                 !  Find topologically compatible candidates and store their
                 !  indices as an accepted or not
                 !-----------------------------------------------------------------
                 CALL ppm_rc_contour_propagation_3d(info)
                 or_fail("ppm_rc_contour_propagation")
                 CALL ppm_util_time(tcontour)
#ifdef __Linux
                 IF (debug.GT.1) THEN
                    CALL ppm_rc_mem_usage(memory,info)
                    stdout("mem_usage at ppm_rc_contour_propagation=",memory)
                 ENDIF
#endif
                 !-----------------------------------------------------------------
                 !  Filtering: Filter unaccepted particles
                 !-----------------------------------------------------------------
                 CALL ppm_rc_FilterCandidates_3d(info)
                 or_fail("ppm_rc_FilterCandidates")
                 CALL ppm_util_time(tfilter1)
#ifdef __Linux
                 IF (debug.GT.1) THEN
                    CALL ppm_rc_mem_usage(memory,info)
                    stdout("mem_usage at ppm_rc_FilterCandidates=",memory)
                 ENDIF
#endif

                 !-----------------------------------------------------------------
                 ! Oscillations: Detecting Oscillations and reduce the moves.
                 !-----------------------------------------------------------------
                 CALL ppm_rc_DetectOscillations2(info)
                 or_fail("ppm_rc_DetectOscillations")
                 CALL ppm_util_time(tosc)
#ifdef __Linux
                 IF (debug.GT.1) THEN
                    CALL ppm_rc_mem_usage(memory,info)
                    stdout("mem_usage at ppm_rc_DetectOscillations=",memory)
                 ENDIF
#endif
                 !-----------------------------------------------------------------
                 ! Filtering: Intermediate step: filter the candidates according
                 ! to their rank and spacial position.
                 !-----------------------------------------------------------------
                 CALL ppm_rc_FilterCandidatesContainerUsingRanks_3d(info)
                 or_fail("ppm_rc_FilterCandidatesContainerUsingRanks")
                 CALL ppm_util_time(tfilter2)
#ifdef __Linux
                 IF (debug.GT.1) THEN
                    CALL ppm_rc_mem_usage(memory,info)
                    stdout("mem_usage at ppm_rc_FilterCandidatesContainerUsingRanks=",memory)
                 ENDIF
#endif
                 !-----------------------------------------------------------------
                 ! Move: Move all the points that are simple. Non simple points
                 ! remain in the candidates list.
                 !-----------------------------------------------------------------
                 Convergence=ppm_rc_move_3d(info)
                 or_fail("ppm_rc_move")
                 CALL ppm_util_time(tmove)
#ifdef __Linux
                 IF (debug.GT.1) THEN
                    CALL ppm_rc_mem_usage(memory,info)
                    stdout("mem_usage at ppm_rc_move=",memory)
                 ENDIF
#endif

                 !-----------------------------------------------------------------
                 ! Write the output.
                 !-----------------------------------------------------------------
                 IF (MOD(istep,freqoutput).EQ.0) THEN
                    CALL ppm_rc_write_image_label_3d(labels,mesh,outputfile,info,idn=istep) !bitsPerSampleW=16
                    or_fail("ppm_rc_write_image")
!                      call mesh%print_vtk(outputname,info,Field=labels)
!                      call ppm_vtk_particles("yaser",Part,info)
                    CALL ppm_util_time(timg)
                    timgs=timgs+timg-tmove
                 ENDIF

              END SELECT

              IF (rank.EQ.0) THEN
                 stdout("The algorithm found ",'e_data%size()'," connected FG region.")
              ENDIF

              tenergys =tenergys+tenergy-t0
              tcontours=tcontours+tcontour-tenergy
              tfilters =tfilters+tfilter1-tcontour+tfilter2-tosc
              toscs    =toscs+tosc-tfilter1
              tmoves   =tmoves+tmove-tfilter2

              done=(istep.GE.maxiter.OR.Convergence)
           ENDDO RC_itrloop !(.NOT.done)

#ifdef  __MPI
           CALL MPI_Ibarrier(comm,request1,info)
           or_fail_MPI("MPI_Ibarrier")
#endif

           CALL e_data%RemoveNotSignificantRegions(info)
           or_fail("RemoveNotSignificantRegions")

           IF (Convergence) THEN
              CALL ppm_util_time(t0)
              SELECT CASE (ppm_rc_dim)
              CASE (2)
                 !-------------------------------------------------------------------------
                 ! check if all went well on image
                 !-------------------------------------------------------------------------
                 CALL ppm_rc_write_image_label_2d(labels,mesh,"converged",info)
                 IF (debug.GT.1) THEN
                    CALL ppm_rc_write_image_label_2d(labels,mesh,"converged",info,idn=FORBIDDEN)
                 ENDIF

              CASE (3)
                 !-------------------------------------------------------------------------
                 ! check if all went well on image
                 !-------------------------------------------------------------------------
                 CALL ppm_rc_write_image_label_3d(labels,mesh,"converged",info)
                 IF (debug.GT.1) THEN
                    CALL ppm_rc_write_image_label_3d(labels,mesh,"converged",info,idn=FORBIDDEN)
                 ENDIF

              END SELECT
              or_fail("ppm_rc_write_image")

              CALL ppm_util_time(timg)
              timgs=timgs+timg-t0
           ELSE
              CALL ppm_util_time(t0)
              SELECT CASE (ppm_rc_dim)
              CASE (2)
                 IF (debug.GT.1) THEN
                    CALL ppm_rc_write_image_label_2d(labels,mesh,outputfile,info,idn=FORBIDDEN)
                 ENDIF

              CASE (3)
                 IF (debug.GT.1) THEN
                    CALL ppm_rc_write_image_label_3d(labels,mesh,outputfile,info,idn=FORBIDDEN)
                 ENDIF

              END SELECT
              or_fail("ppm_rc_write_image")

              CALL ppm_util_time(timg)
              timgs=timgs+timg-t0
           ENDIF

           !-------------------------------------------------------------------------
           ! Writing the labels using IO topology for further analysis
           !-------------------------------------------------------------------------
           IF (ppm_nproc.GT.1) THEN
              CALL ppm_util_time(t0)
              SELECT CASE (ppm_rc_dim)
              CASE (2)
                 CALL ppm_rc_write_image_label_2d(labels,mesh,"iolabel",info,liotopo=.TRUE.)
                 !CALL ppm_rc_write_image_label_2d(labels,mesh,"iolabel",info,idn=FORBIDDEN,liotopo=.TRUE.)
              CASE (3)
                 CALL ppm_rc_write_image_label_3d(labels,mesh,"iolabel",info,liotopo=.TRUE.)
                 !CALL ppm_rc_write_image_label_3d(labels,mesh,"iolabel",info,idn=FORBIDDEN,liotopo=.TRUE.)
              END SELECT
              or_fail("ppm_rc_write_image")

              CALL ppm_util_time(timg)
              timgs=timgs+timg-t0
           ENDIF

#ifdef  __MPI
           ttime( 1)=tenergys
           ttime( 2)=tcontours
           ttime( 3)=tfilters
           ttime( 4)=toscs
           ttime( 5)=timgs
           ttime( 6)=tmoves
           IF (debug.GT.0) THEN
              ttime( 7)=tmove_Simple
              ttime( 8)=tmove_SimpleComm
              ttime( 9)=tmove_NotSimple
              ttime(10)=tmove_Part
              ttime(11)=tmove_ghostfire
           ELSE
              ttime( 7)=zerod
              ttime( 8)=zerod
              ttime( 9)=zerod
              ttime(10)=zerod
              ttime(11)=zerod
           ENDIF

           CALL MPI_Iallreduce(MPI_IN_PLACE,ttime,11,MPI_DOUBLE_PRECISION,MPI_SUM,comm,request2,info)
           or_fail_MPI("MPI_Iallreduce")

           CALL MPI_Wait(request2,MPI_STATUS_IGNORE,info)
           or_fail_MPI("MPI_Wait")

           ttime=ttime/REAL(ppm_nproc,ppm_kind_double)

           tenergys=ttime(1)
           tcontours=ttime(2)
           tfilters=ttime(3)
           toscs=ttime(4)
           timgs=ttime(5)
           tmoves=ttime(6)
           IF (debug.GT.0) THEN
              tmove_Simple=ttime(7)
              tmove_SimpleComm=ttime(8)
              tmove_NotSimple=ttime(9)
              tmove_Part=ttime(10)
              tmove_ghostfire=ttime(11)
           ENDIF
#endif

           !-------------------------------------------------------------------------
           !  Finalize
           !-------------------------------------------------------------------------
           IF (rank.EQ.0) THEN
              stdout("")
              CALL ppm_log(caller,cbuf,info)
              stdout("Segmentation took:           ",'tenergys+tcontours+tfilters+tmoves'," secs.")
              CALL ppm_log(caller,cbuf,info)
              stdout("")
              CALL ppm_log(caller,cbuf,info)
              stdout("Energy Compute took:         ",tenergys, " secs.")
              CALL ppm_log(caller,cbuf,info)
              stdout("Contour propagation took:    ",tcontours," secs.")
              CALL ppm_log(caller,cbuf,info)
              stdout("Filtering took:              ",tfilters, " secs.")
              CALL ppm_log(caller,cbuf,info)
              stdout("Detecting oscillations took: ",toscs,    " secs.")
              CALL ppm_log(caller,cbuf,info)
              stdout("Writing output images took:  ",timgs,    " secs.")
              CALL ppm_log(caller,cbuf,info)
              stdout("Data structure update took:  ",tmoves,   " secs.")
              CALL ppm_log(caller,cbuf,info)
              stdout("")
              !-------------------------------------------------------------------------
              ! Number of connected FG regions without Background
              !-------------------------------------------------------------------------
              stdout("The algorithm found ",'e_data%size()'," connected FG regions.")
              CALL ppm_log(caller,cbuf,info)
              stdout("")
              CALL ppm_log(caller,cbuf,info)
              j=0
              DO i=0,SIZE(e_data%gCount)-1
                 IF (e_data%gCount(i).GT.oned) THEN
                    j=j+1
                    stdout(j,'e_data%Rlabel(i)','e_data%gCount(i)','e_data%gSums(i)')
                    CALL ppm_log(caller,cbuf,info)
                 ENDIF
              ENDDO
              IF (debug.GT.0) THEN
                 stdout("")
                 CALL ppm_log(caller,cbuf,info)
                 stdout("Data structure update took:  ",tmoves,          " secs.")
                 CALL ppm_log(caller,cbuf,info)
                 stdout("tmove_Simple              :  ",tmove_Simple,    " secs.")
                 CALL ppm_log(caller,cbuf,info)
                 stdout("tmove_SimpleComm          :  ",tmove_SimpleComm," secs.")
                 CALL ppm_log(caller,cbuf,info)
                 stdout("tmove_NotSimple           :  ",tmove_NotSimple, " secs.")
                 CALL ppm_log(caller,cbuf,info)
                 stdout("tmove_Part                :  ",tmove_Part,      " secs.")
                 CALL ppm_log(caller,cbuf,info)
                 stdout("tmove_ghostfire           :  ",tmove_ghostfire, " secs.")
                 CALL ppm_log(caller,cbuf,info)
              ENDIF
              stdout("")
              CALL ppm_log(caller,cbuf,info)
              IF (UseMCMC) THEN
                 stdout("Segmentation is Done.")
              ELSE
                 stdout("Segmentation is Done. Goodbye...")
              ENDIF
              CALL ppm_log(caller,cbuf,info)
              stdout("")
              CALL ppm_log(caller,cbuf,info)

              IF (debug.GT.3) THEN
                 CLOSE(3333)
              ENDIF
           ENDIF !(rank.EQ.0)
#ifdef  __MPI
           !-------------------------------------------------------------------------
           ! Wait for everyone to pass the barrier
           !-------------------------------------------------------------------------
           CALL MPI_Wait(request1,MPI_STATUS_IGNORE,info)
           or_fail_MPI("MPI_Wait")
#endif
           CALL ppm_util_time(t0)
        ENDIF !.NOT.UseMCMC.OR.(UseMCMC.AND.MCMCcontinue)

        IF (UseMCMC) THEN
           !-------------------------------------------------------------------------
           ! MCMC sampling
           !-------------------------------------------------------------------------

           CALL ppm_rc_energy_parameter_redefine_mcmc(info)
           or_fail("ppm_rc_energy_parameter_redefine_mcmc")

           !-------------------------------------------------------------------------
           ! MCMC initialization starts
           !-------------------------------------------------------------------------
           SELECT CASE (ppm_rc_dim)
           CASE (2)
              CALL ppm_rc_init_mcmc_2d(info)
           CASE (3)
              CALL ppm_rc_init_mcmc_3d(info)
           END SELECT
           SELECT CASE (info)
           CASE (0)
           CASE (exit_gracefully)
              stdout("exit_gracefully")
              GOTO 7000
           CASE DEFAULT
              or_fail('Failed to initialize!',exit_point=8000)
           END SELECT

           CALL ppm_util_time(tinit)
           ! Make sure initialization is done on every processor
           IF (ppm_rc_debug.GT.0) THEN
              ! If we are on the debug mode every prcessor should print the Initialization time
              stdout("MCMC initialization took: ",'tinit-t0'," secs.")
           ELSE
              IF (rank.EQ.0) THEN
                 stdout("MCMC initialization took: ",'tinit-t0'," secs.")
              ENDIF
           ENDIF

           ! Allocate the MCMColdghost memory for changes in the ghost
           ! Also ask for the mesh%ghost_get
           info=ppm_rc_MCMCinitmove()
           or_fail("ppm_rc_MCMCinitmove")



           Accepted=.FALSE.

           ALLOCATE(index_(ppm_rc_dim),STAT=info)
           or_fail_alloc("index_")

           !-------------------------------------------------------------------------
           !  Do iterations
           !-------------------------------------------------------------------------
           istep=0

           TotalMoves=0_ppm_kind_int64
           AcceptedMoves=0_ppm_kind_int64

           done=.FALSE.

           MCMC_itrloop: DO WHILE (.NOT.done)
              ! Iteration counter
              istep=istep+1

              IF (rank.EQ.0) THEN
                 stdout("Number of steps =",istep," Total number of samples =",TotalMoves)
              ENDIF

              ! The following check is theoretically not needed as there should
              ! always be a floating particle somewhere (given the algorithm
              ! got initialized with more than one region). But in the burn-in phase
              ! we delete floating particles.

              CALL ppm_util_time(t0)

              info=ppm_rc_mcmc_move_2d(LocalMoves,LocalAcceptedMoves)






!
!               IF (MCMCRegionLabelSize.GT.1) THEN
!                  ! Sample a region number (a FG region number; without BG)
!                  i = ppm_rc_Saru_IPRNG(MCMCRegionLabelSize-1)+1
!
!                  ! Find the corresponding label
!                  RLabel = MCMCRegionLabel(i)
!
!                  IF (AllowFission.AND.AllowFusion) THEN
!                     ! using uniform random number
!                     OffBoundarySampling=ppm_rc_Saru_RPRNG().LT.MCMCsampleOffBoundaryPercentage
!
!                     IF (OffBoundarySampling) THEN
!                        ! using uniform random number
!                        Growth=ppm_rc_Saru_RPRNG().LT.half
!
!                        SELECT CASE (ppm_rc_dim)
!                        CASE (2)
!                           Accepted=MCMCOffBoundarySample_2d(Growth,RLabel)
!                        CASE (3)
!                           Accepted=MCMCOffBoundarySample_3d(Growth,RLabel)
!                        END SELECT !ppm_rc_dim
!
!                        AcceptedMoves=AcceptedMoves+MERGE(1,0,Accepted)
!                        done=(istep.GE.maxiter)
!                        CYCLE MCMC_itrloop
!                     ENDIF !OffBoundarySampling
!                  ENDIF !AllowFission.AND.AllowFusion
!
!                  ! finding the label_region index in its array
!                  label_region=htable%search(RLabel)
!                  IF (label_region.EQ.htable_null) THEN
!                     fail("The region does not exist!!!",ppm_error=ppm_error_fatal)
!                  ENDIF
!
!                  IF (MCMCparents(label_region)%size().EQ.0) THEN
!                     ! This is an empty region. Maybe there exists a floating particle
!                     ! with no future for this region. But if Count == 0, it will not
!                     ! be accepted according to the definition of the energy. We hence
!                     ! cleanup the statistics (kill the region).
!                     info=MCMCUpdateRegionLabel()
!                     or_fail("MCMCUpdateRegionLabel")
!
!                     done=(istep.GE.maxiter)
!                     CYCLE MCMC_itrloop
!                  ENDIF
!
!                  ! Figure out if Particle A will cause growth, shrinkage or
!                  ! will be a floating particle
!                  MCMCProbabilityToProposeAFloatingParticle=MCMCFloatingParticlesProposalNormalizer/ &
!                  & (MCMCFloatingParticlesProposalNormalizer + MCMCTotalNormalizer)
!
!                  ! create double precision uniform random number
!                  rnd=ppm_rc_Saru_RPRNGD()
!
!                  ! We will choose one out of the floating particles
!                  ! rnd < F/(F+T)
!                  IF (rnd.LT.MCMCProbabilityToProposeAFloatingParticle) THEN
!                     ParticleAIsFloating=.TRUE.
!
!                     MCMCActiveCandidates => MCMCFloatingParticles
!                  ! rnd < (F+T/2)/(F+T)
!                  ELSE IF (rnd.LT.half*(MCMCProbabilityToProposeAFloatingParticle+oned)) THEN
!                     ParticleAIsFloating=.FALSE.
!
!                     MCMCActiveCandidates => MCMCchildren(label_region)
!                  ! rnd < 1
!                  ELSE
!                     ParticleAIsFloating=.FALSE.
!
!                     MCMCActiveCandidates => MCMCparents(label_region)
!                  ENDIF
!
!                  ! Draw n particles from the discrete distribution
!                  MCMC_Particle_Ab_IsFloating=.FALSE.
!                  MCMC_Particle_Bb_IsFloating=.FALSE.
!                  MCMC_Particle_A_IsFloating=.FALSE.
!                  MCMC_SingleParticleMoveForPairProposals=.FALSE.
!
!                  ! Find particle A: (q_A,CandidateMove,LabelsBeforeJump_A)
!                  info=MCMCFindParticleA(ParticleAIsFloating,ppm_rc_dim,rejected)
!                  or_fail("MCMCFindParticleA")
!
!                  IF (rejected) THEN
!                     done=(istep.GE.MCMCmaxsamples)
!                     CYCLE MCMC_itrloop
!                  ENDIF
!
!                  !!! In case of pair proposals, we find a partner for each proposed particle.
!                  !!! We now know A and Q(A). Now it needs to build
!                  !!! another (2nd step) discrete proposal distribution. We sample from it
!                  !!! to determine the partner particle B. Furthermore we calculate the
!                  !!! conditional proposal probability Q(B|A). In a second step we calculate
!                  !!! the conditional Q(A|B). The same then needs to be done for the backward
!                  !!! probabilities Qb(A), Qb(B), Qb(A|B) and Qb(B|A).
!                  !!! Notation:
!                  !!! - Q is the forward and Qb the backward probability. A is
!                  !!!   a forward praticle and B' the backward particle.
!                  !!! - Qb always assumes backward particles as its arguments! Hence,
!                  !!!   Qb_A_B is the probabily Qb(A'|B').
!                  IF (MCMCusePairProposal) THEN
!                     info=MCMCFindParticleB(ppm_rc_dim)
!                     or_fail("MCMCFindParticleB")
!                  ENDIF !MCMCusePairProposal
!
!
!
!               ENDIF !MCMCRegionLabelSize.GT.1
!
!               AcceptedMoves=AcceptedMoves+MERGE(1,0,Accepted)
!
#ifdef  __MPI
              TmpMoves(1)=LocalMoves
              TmpMoves(2)=LocalAcceptedMoves

              CALL MPI_Iallreduce(MPI_IN_PLACE,TmpMoves,2, &
              &    MPI_INTEGER8,MPI_SUM,comm,request1,info)
              or_fail_MPI("MPI_Iallreduce")

              CALL MPI_Wait(request1,MPI_STATUS_IGNORE,info)
              or_fail_MPI("MPI_Wait")

              LocalMoves=TmpMoves(1)
              LocalAcceptedMoves=TmpMoves(2)
#endif

              TotalMoves=TotalMoves+LocalMoves
              AcceptedMoves=AcceptedMoves+LocalAcceptedMoves

              done=(istep.GE.MCMCmaxiter.OR.TotalMoves.GE.MCMCmaxsamples)
           ENDDO MCMC_itrloop

           DEALLOCATE(index_,STAT=info)
           or_fail_dealloc("index_")
        ENDIF !(UseMCMC)

      7000 CONTINUE
        CALL ppm_rc_finalize(0,info)

        !-------------------------------------------------------------------------
        !  Skip error handler
        !-------------------------------------------------------------------------
        GOTO 9999

        !-------------------------------------------------------------------------
        !  Error handler
        !-------------------------------------------------------------------------
      8000 CONTINUE
        CALL ppm_rc_finalize(1,info)

        !-------------------------------------------------------------------------
        !  End
        !-------------------------------------------------------------------------
      9999 CONTINUE
        STOP
      END PROGRAM ppm_rc
