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
#ifdef __MPI
        USE ppm_module_mpi
#endif

        USE ppm_rc_module_global
        USE ppm_rc_module_linkedlist, ONLY : Candidates,m_Seeds, &
        &   CompetingRegions
        USE ppm_rc_module_init, ONLY : ppm_rc_init_arg,ppm_rc_init_seg_2d, &
        &   ppm_rc_init_seg_3d
        USE ppm_rc_module_energy, ONLY : ppm_rc_energy_compute_2d, &
        &   ppm_rc_energy_compute_3d,e_data
        USE ppm_rc_module_move
        USE ppm_rc_module_finalize, ONLY : ppm_rc_finalize
        USE ppm_rc_module_write
#ifdef __Linux
        USE ppm_rc_module_util, ONLY : ppm_rc_mem_usage
#endif
        IMPLICIT NONE

        !-----------------------------------------------------------------------
        !  Includes
        !-----------------------------------------------------------------------
        !-----------------------------------------------------------------------
        !  Local variables
        !-----------------------------------------------------------------------
        CLASS(ppm_t_subpatch_), POINTER :: sbpitr

        ! variables for timing
        REAL(ppm_kind_double) :: t0,tinit,ttime(11)
        REAL(ppm_kind_double) :: tenergy, tcontour, tfilter1,tfilter2
        REAL(ppm_kind_double) :: tenergys,tcontours,tfilters
        REAL(ppm_kind_double) :: tosc, tmove, timg
        REAL(ppm_kind_double) :: toscs,tmoves,timgs

        INTEGER :: info
#ifdef __Linux
        INTEGER :: memory
#endif
        INTEGER :: i,j,ipatch
#ifdef __MPI
        INTEGER :: request1,request2
#endif
        CHARACTER(LEN=ppm_char) :: caller='ppm_rc'

        LOGICAL :: done
        LOGICAL :: Convergence

        !-------------------------------------------------------------------------
        !  initializing
        !-------------------------------------------------------------------------
        tenergys  = zerod
        tcontours = zerod
        tfilters  = zerod
        toscs     = zerod
        tmoves    = zerod
        timgs     = zerod

        istep = 0

        CALL ppm_rc_init_arg(info)
        SELECT CASE (info)
        CASE (0)
        CASE (exit_gracefully)
           GOTO 7000
        CASE DEFAULT
           or_fail("ppm_rc_init_arg",ppm_error=ppm_error_fatal)
        END SELECT

        !Initialization starts
        CALL ppm_util_time(t0)

        SELECT CASE (ppm_rc_dim)
        CASE (2)
           CALL ppm_rc_init_seg_2d(info)
        CASE (3)
           CALL ppm_rc_init_seg_3d(info)
        END SELECT
        SELECT CASE (info)
        CASE (0)
        CASE (exit_gracefully)
           GOTO 7000
        CASE DEFAULT
           or_fail('Failed to initialize!',exit_point=8000)
        END SELECT

        !-------------------------------------------------------------------------
        !  Write initial diagnostics
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
        !  Write initial condition to output file
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
        !  debugging
        !  debug = 0, no debugging
        !  debug = 1, activate debugging time variables
        !  debug = 2, on LINUX trying to compute the memory use by every processor
        !  debug = 3, writing the labels using IO topology for further analysis
        !  debug = 4, computing and writing the total energy for PC & PS energy functional
        !             the outputfile will be created with number of processors as a suffix
        !             and at every step the total energy will be written there
        !-------------------------------------------------------------------------
        debug=2
        IF (debug.GT.3) THEN
           IF (rank.EQ.0) THEN
              WRITE(cbuf,'(A,A1,I0,A1)') TRIM(inputimage),"_",ppm_nproc,CHAR(0)
              OPEN(3333,FILE=cbuf,STATUS='REPLACE')
           ENDIF
        ENDIF

        CALL ppm_util_time(tinit)

        IF (rank.EQ.0) THEN
           stdout("Initialization took: ",'tinit-t0'," secs.")
        ENDIF

        IF (UseMCMC) THEN
        ELSE !UseMCMC
        !!! RC Segmentation
           ghostsize=1

           done=.FALSE.
           Convergence=.FALSE.

           !-------------------------------------------------------------------------
           !  Do iterations
           !-------------------------------------------------------------------------
           RC_itrloop: DO WHILE (.NOT.done)
              istep=istep+1

              IF (rank.EQ.0) THEN
                 stdout("istep=",istep)
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

                 IF (MOD(istep,freqoutput).EQ.0) THEN
                    CALL ppm_rc_write_image_label_2d(labels, &
                    &    mesh,outputfile,info) !,idn=istep)
                    or_fail("ppm_rc_write_image")
!                      call mesh%print_vtk(outputname,info,Field=labels)
                    CALL ppm_util_time(timg)
                    timgs=timgs+timg-tmove
                 ENDIF
                 IF (rank.EQ.0) THEN
                    stdout("The algorithm found ",'COUNT(e_data%gCount.GT.one)-1'," connected FG region.")
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
                 IF (MOD(istep,freqoutput).EQ.0) THEN
                    CALL ppm_rc_write_image_label_3d(labels, &
                    &    mesh,outputfile,info,idn=istep)
                    or_fail("ppm_rc_write_image")
!                      call mesh%print_vtk(outputname,info,Field=labels)
!                      call ppm_vtk_particles("yaser",Part,info)
                    CALL ppm_util_time(timg)
                    timgs=timgs+timg-tmove
                 ENDIF
                 IF (rank.EQ.0) THEN
                    stdout("The algorithm found ",'COUNT(e_data%gCount.GT.one)-1'," connected FG region.")
                 ENDIF

              END SELECT

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
                 ! check if all went well on image
                 CALL ppm_rc_write_image_label_2d(labels,mesh,"converged",info)
                 IF (debug.GT.1) THEN
                    CALL ppm_rc_write_image_label_2d(labels,mesh,"converged",info,idn=FORBIDDEN)
                 ENDIF

              CASE (3)
                 ! check if all went well on image
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

           IF (debug.GT.2) THEN
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
              stdout("Segmentation took:           ",'tenergys+tcontours+tfilters+tmoves'," secs.")
              stdout("")
              stdout("Energy Compute took:         ",tenergys, " secs.")
              stdout("Contour propagation took:    ",tcontours," secs.")
              stdout("Filtering took:              ",tfilters, " secs.")
              stdout("Detecting oscillations took: ",toscs,    " secs.")
              stdout("Writing output images took:  ",timgs,    " secs.")
              stdout("Data structure update took:  ",tmoves,   " secs.")
              stdout("")
              stdout("The algorithm found ",'COUNT(e_data%gCount.GT.one)-1'," connected FG region.")
              stdout("")
              j=0
              DO i=0,SIZE(e_data%gCount)-1
                 IF (e_data%gCount(i).GT.oned) THEN
                    j=j+1
                    stdout(j,'e_data%Rlabel(i)','e_data%gCount(i)','e_data%gSums(i)')
                 ENDIF
              ENDDO
              IF (debug.GT.0) THEN
                 stdout("")
                 stdout("Data structure update took:  ",tmoves,          " secs.")
                 stdout("tmove_Simple              :  ",tmove_Simple,    " secs.")
                 stdout("tmove_SimpleComm          :  ",tmove_SimpleComm," secs.")
                 stdout("tmove_NotSimple           :  ",tmove_NotSimple, " secs.")
                 stdout("tmove_Part                :  ",tmove_Part,      " secs.")
                 stdout("tmove_ghostfire           :  ",tmove_ghostfire, " secs.")
              ENDIF
              stdout("")
              stdout("Done. Goodbye...")
              stdout("")
              stdout_f('(A)',"****************************************************************")
              CALL ppm_log(caller,cbuf,info)
              stdout_f('(A)',"**        Please do cite:                                     **")
              CALL ppm_log(caller,cbuf,info)
              stdout_f('(A)',"**        `PLoS ONE 11(4):e0152528, (2016)`                   **")
              CALL ppm_log(caller,cbuf,info)
              stdout_f('(A)',"**        when publishing research data obtained using PPM_RC **")
              CALL ppm_log(caller,cbuf,info)
              stdout_f('(A)',"****************************************************************")
              CALL ppm_log(caller,cbuf,info)
           ENDIF
        ENDIF !(UseMCMC)

        IF (debug.GT.3) THEN
           IF (rank.EQ.0) THEN
              CLOSE(3333)
           ENDIF
        ENDIF

#ifdef  __MPI
        !wait for everyone to pass the barrier
        CALL MPI_Wait(request1,MPI_STATUS_IGNORE,info)
        or_fail_MPI("MPI_Wait")
#endif

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
