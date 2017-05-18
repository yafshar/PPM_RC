        !-------------------------------------------------------------------------
        !  Subroutine   :                    ppm_rc_mcmc_move
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
        !  Author           - y.afshar           May   2016
        !-------------------------------------------------------------------------

        !-------------------------------------------------------------------------
        !  Subroutine   :                    ppm_rc_mcmc_move
        !-------------------------------------------------------------------------
        !
        !  Purpose      : Samples the moves to accept and executes them.
        !                 Handles fusion and fission of regions.
        !
        !
        !  Input        :
        !
        !  Input/output :
        !
        !  Output       : LocalMoves         (LI) total number of samples
        !                 LocalAcceptedMoves (LI) total number of accepted moves
        !                 info                (I) return status. 0 on success.
        !
        !  Routines     :
        !
        !  Remarks      :
        !
        !  References   :
        !-------------------------------------------------------------------------
        FUNCTION DTYPE(ppm_rc_mcmc_move)(LocalMoves,LocalAcceptedMoves) RESULT(info)

          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_data, ONLY : ppm_char,ppm_kind_double,ppm_kind_int64, &
          &   ppm_nproc,ppm_error_error,ppm_rank
          USE ppm_module_substart, ONLY : substart
          USE ppm_module_substop, ONLY : substop
          USE ppm_module_error, ONLY : ppm_error,ppm_err_sub_failed, &
          &   ppm_err_alloc,ppm_err_dealloc,ppm_err_mpi_fail
          USE ppm_module_write
          USE ppm_module_topo_typedef, ONLY : ppm_t_topo,ppm_topo
          USE ppm_module_mesh_typedef, ONLY : ppm_t_subpatch_

          USE ppm_module_mpi
          USE ppm_module_util_qsort, ONLY : ppm_util_qsort
          USE ppm_module_util_time, ONLY : ppm_util_time

          USE ppm_rc_module_global, ONLY : MK,procflag,proccolor,mesh,labels, &
          &   ppm_rc_dim,istep,topoid,comm,usrseed
          USE ppm_rc_module_linkedlist, ONLY : MCMCParticleInContainerHistory, &
          &   MCMCFloatingParticleInContainerHistory,MCMCLabelImageHistory
          USE ppm_rc_module_util, ONLY : ppm_rc_label_exist,ppm_rc_label_index, &
          &   ppm_rc_shuffle
          USE ppm_rc_module_rnd, ONLY : ppm_rc_Saru_SEED,ppm_rc_Saru_RPRNG
          USE ppm_rc_module_fire, ONLY : DTYPE(ppm_rc_ghost_copy)
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER(ppm_kind_int64), INTENT(  OUT) :: LocalMoves
          INTEGER(ppm_kind_int64), INTENT(  OUT) :: LocalAcceptedMoves
          INTEGER                                :: info

          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          CLASS(ppm_t_subpatch_), POINTER :: sbpitr

          TYPE(ppm_t_topo), POINTER :: topo

          REAL(ppm_kind_double) :: t0

          ! Local variables for ghost communications
          INTEGER, DIMENSION(:), ALLOCATABLE :: sendlist
          INTEGER, DIMENSION(:), ALLOCATABLE :: recvlist
          INTEGER                            :: nsendlist
          INTEGER                            :: nrecvlist
          ! Dummy pointer for number of mesh points in each direction
          INTEGER, DIMENSION(:), POINTER     :: Nm
          ! Active processor
          INTEGER                            :: activecolor
          INTEGER                            :: cellcolor
          ! Iterator
          INTEGER                            :: nsweeps
          ! Subsweep iterator
          INTEGER                            :: subsweep
          ! Maximum number of particles in any cell
          INTEGER                            :: NpMax
          ! Number of particles in currebt cbox cell
          INTEGER                            :: Np
          ! Current box index
          INTEGER                            :: cbox
          ! Start box number
          INTEGER                            :: sbox
          ! End box number
          INTEGER                            :: ebox
          INTEGER                            :: i
          INTEGER                            :: j
          INTEGER                            :: k
#ifdef  __MPI
          INTEGER                            :: request
#endif

          CHARACTER(LEN=ppm_char) :: caller='ppm_rc_mcmc_move'

          LOGICAL :: Accept
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          !-------------------------------------------------------------------------
          !  Find the maximum number of particles in any cell
          !-------------------------------------------------------------------------
          NpMax=0
          DO i=1,PRODUCT(MCMCcellNm)
             j=MCMCRegularParticlesInCell(i)%size()
             k=MCMCFloatingParticlesInCell(i)%size()
             NpMax=MAX(NpMax,j+k)
          ENDDO

#ifdef  __MPI
          IF (ppm_nproc.GT.1) THEN
             CALL MPI_Iallreduce(MPI_IN_PLACE,NpMax,1,MPI_INTEGER,MPI_MAX,comm,request,info)
             or_fail_MPI("MPI_Iallreduce")
          ENDIF
#endif

          LocalMoves=0_ppm_kind_int64
          LocalAcceptedMoves=0_ppm_kind_int64

          !-------------------------------------------------------------------------
          ! Initialize RANDOM_SEED with the istep so every processor produce
          ! the same sequence
          !-------------------------------------------------------------------------
          info=ppm_rc_Saru_SEED(istep,usrseed)

          !-------------------------------------------------------------------------
          ! To avoid bias in our sampling, we shuffle the independent cell sets.
          ! Every processor produces the same sequence
          !-------------------------------------------------------------------------
          CALL ppm_rc_shuffle(procflag)
          !-------------------------------------------------------------------------
          ! Boundary cell Indices do not need to be the same on every processor
          !-------------------------------------------------------------------------
          CALL ppm_rc_shuffle(MCMCboundarycellIndex,MCMCboundarycellSize)

          !-------------------------------------------------------------------------
          ! Pointer to the current topo
          !-------------------------------------------------------------------------
          topo => ppm_topo(topoid)%t

          !-------------------------------------------------------------------------
          ! Dummy iterator
          !-------------------------------------------------------------------------
          subsweep=0

#ifdef  __MPI
          IF (ppm_nproc.GT.1) THEN
             CALL MPI_Wait(request,MPI_STATUS_IGNORE,info)
             or_fail_MPI("MPI_Wait")
          ENDIF
#endif

          cell_loop: DO nsweeps=1,ISHFT(1,ppm_rc_dim) ! 4 or 8 in 2D & 3D respectively
             ! Find out the active color from the shuffled independent cell sets
             activecolor=procflag(nsweeps)

             nsendlist=0
             nrecvlist=0

             !-------------------------------------------------------------------------
             ! This is an active processor
             ! The active processor is sampling in its boundary cells
             !-------------------------------------------------------------------------
             IF (proccolor.EQ.activecolor) THEN
                !-------------------------------------------------------------------------
                ! This processor now does the border for ghost communications
                ! At first samples the border area
                ! Then non-blocking send/recive
                ! Wait
                ! Go to the inside
                !-------------------------------------------------------------------------
                IF (ppm_nproc.GT.1) THEN
                   DO i=1,topo%nneighproc
                      IF (topo%ineighcolor(i).NE.activecolor) nsendlist=nsendlist+1
                   ENDDO

                   ALLOCATE(sendlist(nsendlist),recvlist(nrecvlist),STAT=info)
                   or_fail_alloc("sendlist & recvlist")

                   j=0
                   DO i=1,topo%nneighproc
                      IF (topo%ineighcolor(i).NE.activecolor) THEN
                         j=j+1
                         sendlist(j)=topo%ineighcolor(i)
                      ENDIF
                   ENDDO
                ENDIF !(ppm_nproc.GT.1)

                ! Sub-sweep in boundary cells
                DO i=1,MCMCboundarycellSize
                   ! Current box index
                   cbox=MCMCboundarycellIndex(i)

                   j=MCMCRegularParticlesInCell(cbox)%size()
                   k=MCMCFloatingParticlesInCell(cbox)%size()
                   ! Number of particles in the current cell cbox
                   Np=j+k
                   ! TODO
                   ! I am not sure how should I do this
                   ! For now I consider regular particles (both of regular and floating particles)

                   IF (Np.EQ.0) CYCLE

                   LocalMoves=LocalMoves+1_ppm_kind_int64

                   IF (ppm_rc_Saru_RPRNG().LT.REAL(Np,MK)/REAL(NpMax,MK)) THEN
                      !----------------------------------------------------------------------
                      ! Clear the old sets
                      ! (These sets will help to revert the move in case it gets rejected).
                      !----------------------------------------------------------------------
                      CALL MCMCParticleInContainerHistory%destroy(info)
                      or_fail("MCMCParticleInContainerHistory%destroy")

                      CALL MCMCFloatingParticleInContainerHistory%destroy(info)
                      or_fail("MCMCFloatingParticleInContainerHistory%destroy")

                      CALL MCMCLabelImageHistory%destroy(info)
                      or_fail("MCMCLabelImageHistory%destroy")

                      !----------------------------------------------------------------------
                      ! Currently it is possible that the same candidate is in the move set.
                      ! Hence we store the applied moves to avoid duplicates.
                      !----------------------------------------------------------------------

                      !----------------------------------------------------------------------
                      ! Draw n particles from the discrete distribution
                      !----------------------------------------------------------------------
                      MCMC_Particle_Ab_IsFloating =.FALSE.
                      MCMC_Particle_A_IsFloating  =.FALSE.

                      !----------------------------------------------------------------------
                      ! Find particle A
                      !----------------------------------------------------------------------
                      IF (.NOT.MCMCFindParticleA(__DIME,cbox,j,k,info)) CYCLE
                      or_fail("MCMCFindParticleA")

                      Accept=MCMCMove(__DIME,LocalMoves,info)
                      or_fail("MCMCMove")

                      LocalAcceptedMoves=LocalAcceptedMoves+MERGE(1_ppm_kind_int64,0_ppm_kind_int64,Accept)

                      !----------------------------------------------------------------------
                      ! Iterate the candidates, calculate the energy and perform the moves
                      !----------------------------------------------------------------------
                   ENDIF !(ppm_rc_Saru_RPRNG().LT.REAL(Np,MK)/REAL(NpMax,MK))
                ENDDO !i=1,MCMCboundarycellSize

                !-------------------------------------------------------------------------
                ! Ghost update for labels
                ! Any label update at this part should be seen at the other processor
                !-------------------------------------------------------------------------
                IF (ppm_nproc.GT.1) THEN
                   CALL mesh%map_ghost_get(info,sendlist=sendlist,recvlist=recvlist)
                   or_fail("labels%map_ghost_push")

                   CALL labels%map_ghost_push(mesh,info)
                   or_fail("labels%map_ghost_push")
                   CALL mesh%map_isend(info,sendrecv=.TRUE.)
                   or_fail("mesh%map_isend")

                   DEALLOCATE(sendlist,recvlist,STAT=info)
                   or_fail_dealloc("sendlist & recvlist")

                   CALL mesh%map_isend(info,sendrecv=.FALSE.)
                   or_fail("mesh%map_isend")
                   CALL labels%map_ghost_pop(mesh,info)
                   or_fail("labels%map_ghost_pop")
                ENDIF
             !-------------------------------------------------------------------------
             ! This is an inactive processor sampling in interior cells
             ! the proccolor is diferent from activecolor
             !------------------------------------------------------------------------
             ELSE
                !------------------------------------------------------------------------
                ! This processor now works on the inside cells
                ! First non-blocking recieve
                ! Work inside on cellcolor cells
                ! Wait
                !------------------------------------------------------------------------
                IF (ppm_nproc.GT.1) THEN
                   DO i=1,topo%nneighproc
                      IF (topo%ineighcolor(i).EQ.activecolor) nrecvlist=nrecvlist+1
                   ENDDO

                   ALLOCATE(sendlist(nsendlist),recvlist(nrecvlist),STAT=info)
                   or_fail_alloc("sendlist & recvlist")

                   j=0
                   DO i=1,topo%nneighproc
                      IF (topo%ineighcolor(i).EQ.activecolor) THEN
                         j=j+1
                         recvlist(j)=topo%ineighcolor(i)
                      ENDIF
                   ENDDO

                   !-------------------------------------------------------------------------
                   ! Ghost update for labels
                   ! Any label update at this part should be seen at the other processor
                   !-------------------------------------------------------------------------
                   CALL mesh%map_ghost_get(info,sendlist=sendlist,recvlist=recvlist)
                   or_fail("labels%map_ghost_push")

                   CALL labels%map_ghost_push(mesh,info)
                   or_fail("labels%map_ghost_push")
                   CALL mesh%map_isend(info,sendrecv=.TRUE.)
                   or_fail("mesh%map_isend")

                   DEALLOCATE(sendlist,recvlist,STAT=info)
                   or_fail_dealloc("sendlist & recvlist")
                ENDIF !ppm_nproc.GT.1

                subsweep=subsweep+1
                cellcolor=procflag(subsweep)

                sbox=MCMCinteriorcellDisp(cellcolor-1)+1
                ebox=MCMCinteriorcellDisp(cellcolor)
                ! Sub-sweep in interior cells
                DO i=sbox,ebox
                   ! Index of the current box
                   cbox=MCMCinteriorcellIndex(i)
                   ! Count number of particles in cell cbox
                   j=MCMCRegularParticlesInCell(i)%size()
                   k=MCMCFloatingParticlesInCell(i)%size()
                   Np=j+k
                   ! TODO
                   ! I am not sure how should I do this
                   ! Should I only consider regular particles or both of regular and floating particles
                   IF (Np.EQ.0) CYCLE
                   IF (ppm_rc_Saru_RPRNG().LT.REAL(Np,MK)/REAL(NpMax,MK)) THEN
                      !----------------------------------------------------------------------
                      !  Clear the old sets
                      !  These list will help to revert the move in case it gets rejected.
                      !----------------------------------------------------------------------
                      CALL MCMCParticleInContainerHistory%destroy(info)
                      or_fail("MCMCParticleInContainerHistory%destroy")

                      CALL MCMCFloatingParticleInContainerHistory%destroy(info)
                      or_fail("MCMCFloatingParticleInContainerHistory%destroy")

                      CALL MCMCLabelImageHistory%destroy(info)
                      or_fail("MCMCLabelImageHistory%destroy")

                      ! Currently it is possible that the same candidate is in the move set.
                      ! Hence we store the applied moves to avoid duplicates.
                   ELSE
                      ! TODO
                      ! Now I consider this as a rejection
                      LocalMoves=LocalMoves+1_ppm_kind_int64
                   ENDIF
                ENDDO

                !-------------------------------------------------------------------------
                ! TODO
                ! Here I can do an optimization, by comparison of the MCMCghostcellSize
                ! With the active MCMCcoloringcellSize and if there is ahuge ratio,
                ! I can stop the loop in the middle get the ghost and then continue the
                ! loop Then All the processors are kind of synchronization
                !-------------------------------------------------------------------------

#if   __DIME == __2D
#elif __DIME == __3D
#endif
                IF (ppm_nproc.GT.1) THEN
                   sbpitr => mesh%subpatch%begin()
                   DO WHILE (ASSOCIATED(sbpitr))
                      !-------------------------------------------------------------------------
                      ! Copy the ghost for comparison with the new values and impose the changes
                      ! to inside domain
                      !-------------------------------------------------------------------------
                      CALL DTYPE(ppm_rc_ghost_copy)(sbpitr,MCMColdghost,1,info)
                      or_fail("ppm_rc_ghost_copy")

                      sbpitr => mesh%subpatch%next()
                   ENDDO !WHILE(ASSOCIATED(sbpitr))

                   CALL mesh%map_isend(info,sendrecv=.FALSE.)
                   or_fail("mesh%map_isend")
                   CALL labels%map_ghost_pop(mesh,info)
                   or_fail("labels%map_ghost_pop")

                   sbpitr => mesh%subpatch%begin()
                   DO WHILE (ASSOCIATED(sbpitr))
                      Nm => sbpitr%nnodes
! #include "./ChangeGhostContourParticleLabelToCandidateLabel.inc"
                      sbpitr => mesh%subpatch%next()
                   ENDDO !WHILE(ASSOCIATED(sbpitr))
                ENDIF !ppm_nproc.GT.1
             ENDIF !(proccolor.EQ.activecolor)
          ENDDO cell_loop

          !-------------------------------------------------------------------------
          ! The last sub-sweep, where all processors sample in active interior cells.
          !-------------------------------------------------------------------------
          subsweep=ISHFT(1,ppm_rc_dim)
          cellcolor=procflag(subsweep)
          sbox=MCMCinteriorcellDisp(cellcolor-1)+1
          ebox=MCMCinteriorcellDisp(cellcolor)
          DO i=sbox,ebox
             ! Index of the current box
             cbox=MCMCinteriorcellIndex(i)
             ! Count number of particles in cell cbox
             j=MCMCRegularParticlesInCell(i)%size()
             k=MCMCFloatingParticlesInCell(i)%size()
             Np=j+k
             ! TODO
             ! I am not sure how should I do this
             ! Should I only consider regular particles or both of regular and floating particles
             IF (Np.EQ.0) CYCLE
             IF (ppm_rc_Saru_RPRNG().LT.REAL(Np,MK)/REAL(NpMax,MK)) THEN
                !----------------------------------------------------------------------
                !  Clear the old sets
                !  These list will help to revert the move in case it gets rejected.
                !----------------------------------------------------------------------
                CALL MCMCParticleInContainerHistory%destroy(info)
                or_fail("MCMCParticleInContainerHistory%destroy")

                CALL MCMCFloatingParticleInContainerHistory%destroy(info)
                or_fail("MCMCFloatingParticleInContainerHistory%destroy")

                CALL MCMCLabelImageHistory%destroy(info)
                or_fail("MCMCLabelImageHistory%destroy")

                ! Currently it is possible that the same candidate is in the move set.
                ! Hence we store the applied moves to avoid duplicates.
             ELSE
                ! TODO
                ! Now I consider this as a rejection
                LocalMoves=LocalMoves+1_ppm_kind_int64
             ENDIF
          ENDDO
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        END FUNCTION DTYPE(ppm_rc_mcmc_move)



