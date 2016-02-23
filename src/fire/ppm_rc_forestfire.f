
        SUBROUTINE DTYPE(initforestfire)(MeshIn,initghost_size,info)

          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_mpi
          USE ppm_module_util_qsort, ONLY : ppm_util_qsort
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !------------------------------------------------------------------------
          !  Arguments
          !------------------------------------------------------------------------
          CLASS(ppm_t_equi_mesh_), POINTER       :: MeshIn

          INTEGER,                 INTENT(IN   ) :: initghost_size
          ! the size of the whole 1 layer ghost around all subpatches
          INTEGER,                 INTENT(  OUT) :: info

          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          CLASS(ppm_t_subpatch_), POINTER :: sbpitr

          TYPE(ppm_rc_list), POINTER :: seed
          TYPE(ppm_rc_stat), POINTER :: trstat

#if   __DIME == __2D
          REAL(MK), CONTIGUOUS, DIMENSION(:,:),   POINTER :: DTYPE(wpi)
#elif __DIME == __3D
          REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER :: DTYPE(wpi)
#endif
          REAL(ppm_kind_double),DIMENSION(:),     POINTER :: value
          REAL(ppm_kind_double),DIMENSION(3)              :: val
          REAL(ppm_kind_double)                           :: t0,dummy

#if   __DIME == __2D
          INTEGER, CONTIGUOUS, DIMENSION(:,:),   POINTER :: DTYPE(wpl)
#elif __DIME == __3D
          INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER :: DTYPE(wpl)
#endif
          INTEGER, DIMENSION(:), ALLOCATABLE :: old_ghost
          INTEGER, DIMENSION(:),     POINTER :: Nm
          INTEGER, DIMENSION(:),     POINTER :: unique_new_regions
          INTEGER, DIMENSION(__DIME)              :: ld,ldu
          INTEGER                            :: iopt,isize,nsize
          INTEGER                            :: i,j,l,ipatch
#if   __DIME == __3D
          INTEGER                            :: k
#endif
          INTEGER                            :: nrgs,sum_nrgs
          INTEGER                            :: htable_size
          INTEGER                            :: Nregions
#ifdef __MPI
          INTEGER, DIMENSION(:), ALLOCATABLE :: displ,counts
          INTEGER                            :: request
#endif

          INTEGER(ppm_kind_int64)            :: ll

          CHARACTER(ppm_char) :: caller="initforestfire"

          LOGICAL                            :: MASK
          LOGICAL, DIMENSION(:), ALLOCATABLE :: MASKL

          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          ALLOCATE(tmp_region_stat,STAT=info)
          or_fail_alloc("tmp_region_stat")

          ldu(1)=MeshIn%subpatch%nb
          ALLOCATE(seednm(1:ldu(1)),STAT=info)
          or_fail_alloc("seednm")

          ! first step in recognizig connected regions
          CALL DTYPE(initfire)(MeshIn,info)
          or_fail("fire")

          !-------------------------------------------------------------------------
          !  ghost update.
          !-------------------------------------------------------------------------
          IF (ppm_nproc.GT.1) THEN
             ghostsize=1
             !We only need one layer of ghosts to complete clustering of
             !different regions
             CALL MeshIn%map_ghost_get(info,ghostsize=ghostsize)
             or_fail("MeshIn%map_ghost_get")

             CALL labels%map_ghost_push(MeshIn,info)
             or_fail("labels%map_ghost_push")
             CALL MeshIn%map_isend(info,sendrecv=.TRUE.)
             or_fail("MeshIn%map_send")
          ENDIF

          ALLOCATE(old_ghost(initghost_size),STAT=info)
          or_fail_alloc("old_ghost")

          IF (initghost_size.GT.0) THEN
             sbpitr => MeshIn%subpatch%begin()
             isize=0
             nsize=0
             DO WHILE (ASSOCIATED(sbpitr))
                Nm => sbpitr%nnodes
#if   __DIME == __2D
                IF (sbpitr%istart(1).NE.1) THEN
                   nsize=nsize+Nm(2)
                ENDIF
                IF (sbpitr%iend(1).NE.Ngrid(1)) THEN
                   nsize=nsize+Nm(2)
                ENDIF
                IF (sbpitr%istart(2).NE.1) THEN
                   nsize=nsize+Nm(1)+2
                ENDIF
                IF (sbpitr%iend(2).NE.Ngrid(2)) THEN
                   nsize=nsize+Nm(1)+2
                ENDIF
#elif __DIME == __3D
                IF (sbpitr%istart(1).NE.1) THEN
                   nsize=nsize+Nm(2)*Nm(3)
                ENDIF
                IF (sbpitr%iend(1).NE.Ngrid(1)) THEN
                   nsize=nsize+Nm(2)*Nm(3)
                ENDIF
                IF (sbpitr%istart(2).NE.1) THEN
                   nsize=nsize+Nm(3)*(Nm(1)+2)
                ENDIF
                IF (sbpitr%iend(2).NE.Ngrid(2)) THEN
                   nsize=nsize+Nm(3)*(Nm(1)+2)
                ENDIF
                IF (sbpitr%istart(3).NE.1) THEN
                   nsize=nsize+(Nm(1)+2)*(Nm(2)+2)
                ENDIF
                IF (sbpitr%iend(3).NE.Ngrid(3)) THEN
                   nsize=nsize+(Nm(1)+2)*(Nm(2)+2)
                ENDIF
#endif
                CALL DTYPE(ppm_rc_ghost_copy)(sbpitr,old_ghost(isize+1:nsize),info)
                or_fail("ppm_rc_ghost_copy")

                isize=nsize
                sbpitr => MeshIn%subpatch%next()
             ENDDO
          ENDIF !initghost_size.GT.0

          IF (ppm_nproc.GT.1) THEN
             CALL MeshIn%map_isend(info,sendrecv=.FALSE.)
             or_fail("MeshIn%map_send")
             CALL labels%map_ghost_pop(MeshIn,info)
             or_fail("labels%map_ghost_pop")
          ENDIF

          NULLIFY(DTYPE(wpl))

          sbpitr => MeshIn%subpatch%begin()
          nsize=0
          ipatch=1
          DO WHILE (ASSOCIATED(sbpitr))
             Nm => sbpitr%nnodes

             CALL sbpitr%get_field(labels,DTYPE(wpl),info)
             or_fail("Failed to get field i_wp data.")

#if   __DIME == __2D
             !-x west
             i=0
             l=nsize
             IF (sbpitr%istart(1).NE.1) THEN
                DO j=1,Nm(2)
                   l=l+1
                   IF ((ABS(old_ghost(l)).GT.ABS(DTYPE(wpl)(i,j)) &
                   &   .AND.ABS(DTYPE(wpl)(i,j)).GT.1)) THEN
                      ALLOCATE(seed,STAT=info)
                      or_fail_alloc("seed")
                      CALL seed%add(i,j)
                      CALL ppm_rc_seeds(ipatch)%push(seed,info)
                      or_fail("could not add new seed to the collection")
                   ELSE
                      DTYPE(wpl)(i,j)=old_ghost(l)
                   ENDIF
                ENDDO !j

                nsize=nsize+Nm(2)
             ENDIF
             !+x east
             i=Nm(1)+1
             l=nsize
             IF (sbpitr%iend(1).NE.Ngrid(1)) THEN
                DO j=1,Nm(2)
                   l=l+1
                   IF ((ABS(old_ghost(l)).GT.ABS(DTYPE(wpl)(i,j)) &
                   &   .AND.ABS(DTYPE(wpl)(i,j)).GT.1)) THEN
                      ALLOCATE(seed,STAT=info)
                      or_fail_alloc("seed")
                      CALL seed%add(i,j)
                      CALL ppm_rc_seeds(ipatch)%push(seed,info)
                      or_fail("could not add new seed to the collection")
                   ELSE
                      DTYPE(wpl)(i,j)=old_ghost(l)
                   ENDIF
                ENDDO !j
                nsize=nsize+Nm(2)
             ENDIF
             !-y south
             j=0
             l=nsize
             IF (sbpitr%istart(2).NE.1) THEN
                i=0
                l=l+1
                DTYPE(wpl)(i,j)=old_ghost(l)
                DO i=1,Nm(1)
                   l=l+1
                   IF ((ABS(old_ghost(l)).GT.ABS(DTYPE(wpl)(i,j)) &
                   &   .AND.ABS(DTYPE(wpl)(i,j)).GT.1)) THEN
                      ALLOCATE(seed,STAT=info)
                      or_fail_alloc("seed")
                      CALL seed%add(i,j)
                      CALL ppm_rc_seeds(ipatch)%push(seed,info)
                      or_fail("could not add new seed to the collection")
                   ELSE
                      DTYPE(wpl)(i,j)=old_ghost(l)
                   ENDIF
                ENDDO !i
                i=Nm(1)+1
                l=l+1
                DTYPE(wpl)(i,j)=old_ghost(l)
                nsize=nsize+Nm(1)+2
             ENDIF
             !+y north
             j=Nm(2)+1
             l=nsize
             IF (sbpitr%iend(2).NE.Ngrid(2)) THEN
                i=0
                l=l+1
                DTYPE(wpl)(i,j)=old_ghost(l)
                DO i=1,Nm(1)
                   l=l+1
                   IF ((ABS(old_ghost(l)).GT.ABS(DTYPE(wpl)(i,j)) &
                   &   .AND.ABS(DTYPE(wpl)(i,j)).GT.1)) THEN
                      ALLOCATE(seed,STAT=info)
                      or_fail_alloc("seed")
                      CALL seed%add(i,j)
                      CALL ppm_rc_seeds(ipatch)%push(seed,info)
                      or_fail("could not add new seed to the collection")
                   ELSE
                      DTYPE(wpl)(i,j)=old_ghost(l)
                   ENDIF
                ENDDO !j
                i=Nm(1)+1
                l=l+1
                DTYPE(wpl)(i,j)=old_ghost(l)
                nsize=nsize+Nm(1)+2
             ENDIF
#elif __DIME == __3D
             !-x west
             i=0
             l=nsize
             IF (sbpitr%istart(1).NE.1) THEN
                DO k=1,Nm(3)
                   DO j=1,Nm(2)
                      l=l+1
                      IF ((ABS(old_ghost(l)).GT.ABS(DTYPE(wpl)(i,j,k)) &
                      &   .AND.ABS(DTYPE(wpl)(i,j,k)).GT.1)) THEN
                         ALLOCATE(seed,STAT=info)
                         or_fail_alloc("seed")
                         CALL seed%add(i,j,k)
                         CALL ppm_rc_seeds(ipatch)%push(seed,info)
                         or_fail("could not add new seed to the collection")
                      ELSE
                         DTYPE(wpl)(i,j,k)=old_ghost(l)
                      ENDIF
                   ENDDO !j
                ENDDO !k
                nsize=nsize+Nm(2)*Nm(3)
             ENDIF
             !+x east
             i=Nm(1)+1
             l=nsize
             IF (sbpitr%iend(1).NE.Ngrid(1)) THEN
                DO k=1,Nm(3)
                   DO j=1,Nm(2)
                      l=l+1
                      IF ((ABS(old_ghost(l)).GT.ABS(DTYPE(wpl)(i,j,k)) &
                      &   .AND.ABS(DTYPE(wpl)(i,j,k)).GT.1)) THEN
                         ALLOCATE(seed,STAT=info)
                         or_fail_alloc("seed")
                         CALL seed%add(i,j,k)
                         CALL ppm_rc_seeds(ipatch)%push(seed,info)
                         or_fail("could not add new seed to the collection")
                      ELSE
                         DTYPE(wpl)(i,j,k)=old_ghost(l)
                      ENDIF
                   ENDDO !j
                ENDDO !k
                nsize=nsize+Nm(2)*Nm(3)
             ENDIF
             !-y south
             j=0
             l=nsize
             IF (sbpitr%istart(2).NE.1) THEN
                DO k=1,Nm(3)
                   i=0
                   l=l+1
                   DTYPE(wpl)(i,j,k)=old_ghost(l)
                   DO i=1,Nm(1)
                      l=l+1
                      IF ((ABS(old_ghost(l)).GT.ABS(DTYPE(wpl)(i,j,k)) &
                      &   .AND.ABS(DTYPE(wpl)(i,j,k)).GT.1)) THEN
                         ALLOCATE(seed,STAT=info)
                         or_fail_alloc("seed")
                         CALL seed%add(i,j,k)
                         CALL ppm_rc_seeds(ipatch)%push(seed,info)
                         or_fail("could not add new seed to the collection")
                      ELSE
                         DTYPE(wpl)(i,j,k)=old_ghost(l)
                      ENDIF
                   ENDDO !i
                   i=Nm(1)+1
                   l=l+1
                   DTYPE(wpl)(i,j,k)=old_ghost(l)
                ENDDO !k
                nsize=nsize+(Nm(1)+2)*Nm(3)
             ENDIF
             !+y north
             j=Nm(2)+1
             l=nsize
             IF (sbpitr%iend(2).NE.Ngrid(2)) THEN
                DO k=1,Nm(3)
                   i=0
                   l=l+1
                   DTYPE(wpl)(i,j,k)=old_ghost(l)
                   DO i=1,Nm(1)
                      l=l+1
                      IF ((ABS(old_ghost(l)).GT.ABS(DTYPE(wpl)(i,j,k)) &
                      &   .AND.ABS(DTYPE(wpl)(i,j,k)).GT.1)) THEN
                         ALLOCATE(seed,STAT=info)
                         or_fail_alloc("seed")
                         CALL seed%add(i,j,k)
                         CALL ppm_rc_seeds(ipatch)%push(seed,info)
                         or_fail("could not add new seed to the collection")
                      ELSE
                         DTYPE(wpl)(i,j,k)=old_ghost(l)
                      ENDIF
                   ENDDO !i
                   i=Nm(1)+1
                   l=l+1
                   DTYPE(wpl)(i,j,k)=old_ghost(l)
                ENDDO !k
                nsize=nsize+(Nm(1)+2)*Nm(3)
             ENDIF
             !-z bottom
             k=0
             l=nsize
             IF (sbpitr%istart(3).NE.1) THEN
                j=0
                DO i=0,Nm(1)+1
                   l=l+1
                   DTYPE(wpl)(i,j,k)=old_ghost(l)
                ENDDO
                DO j=1,Nm(2)
                   i=0
                   l=l+1
                   DTYPE(wpl)(i,j,k)=old_ghost(l)
                   DO i=1,Nm(1)
                      l=l+1
                      IF ((ABS(old_ghost(l)).GT.ABS(DTYPE(wpl)(i,j,k)) &
                      &   .AND.ABS(DTYPE(wpl)(i,j,k)).GT.1)) THEN
                         ALLOCATE(seed,STAT=info)
                         or_fail_alloc("seed")
                         CALL seed%add(i,j,k)
                         CALL ppm_rc_seeds(ipatch)%push(seed,info)
                         or_fail("could not add new seed to the collection")
                      ELSE
                         DTYPE(wpl)(i,j,k)=old_ghost(l)
                      ENDIF
                   ENDDO !i
                   i=Nm(1)+1
                   l=l+1
                   DTYPE(wpl)(i,j,k)=old_ghost(l)
                ENDDO !j
                j=Nm(2)+1
                DO i=0,Nm(1)+1
                   l=l+1
                   DTYPE(wpl)(i,j,k)=old_ghost(l)
                ENDDO
                nsize=nsize+(Nm(1)+2)*(Nm(2)+2)
             ENDIF
             !+z top
             k=Nm(3)+1
             l=nsize
             IF (sbpitr%iend(3).NE.Ngrid(3)) THEN
                j=0
                DO i=0,Nm(1)+1
                   l=l+1
                   DTYPE(wpl)(i,j,k)=old_ghost(l)
                ENDDO
                DO j=1,Nm(2)
                   i=0
                   l=l+1
                   DTYPE(wpl)(i,j,k)=old_ghost(l)
                   DO i=1,Nm(1)
                      l=l+1
                      IF ((ABS(old_ghost(l)).GT.ABS(DTYPE(wpl)(i,j,k)) &
                      &   .AND.ABS(DTYPE(wpl)(i,j,k)).GT.1)) THEN
                         ALLOCATE(seed,STAT=info)
                         or_fail_alloc("seed")
                         CALL seed%add(i,j,k)
                         CALL ppm_rc_seeds(ipatch)%push(seed,info)
                         or_fail("could not add new seed to the collection")
                      ELSE
                         DTYPE(wpl)(i,j,k)=old_ghost(l)
                      ENDIF
                   ENDDO !i
                   i=Nm(1)+1
                   l=l+1
                   DTYPE(wpl)(i,j,k)=old_ghost(l)
                ENDDO !j
                j=Nm(2)+1
                DO i=0,Nm(1)+1
                   l=l+1
                   DTYPE(wpl)(i,j,k)=old_ghost(l)
                ENDDO
                nsize=nsize+(Nm(1)+2)*(Nm(2)+2)
             ENDIF
#endif

             sbpitr => MeshIn%subpatch%next()
             ipatch=ipatch+1
          ENDDO

          NULLIFY(DTYPE(wpl))

          l=MeshIn%subpatch%nb
          ldu(1)=SUM(ppm_rc_seeds(1:l)%nb-seednm(1:l))

          MASK=(ldu(1).GT.0)
          ghost_on_fire=MERGE(1,0,MASK)

#ifdef __MPI
          CALL MPI_Iallreduce(MPI_IN_PLACE,ghost_on_fire,1,&
          &    MPI_INTEGER,MPI_SUM,comm,request,info)
          or_fail_MPI("MPI_Iallreduce")
#endif

          ALLOCATE(ghost_nlabels(ldu(1)),STAT=info)
          or_fail_alloc("ghost_nlabels")

          j=SIZE(nlabels,DIM=1)
          DO i=j,1,-1
             IF (nlabels(i).GT.0) EXIT
          ENDDO
          IF (i.EQ.0) THEN
             IF (j.GT.0) THEN
                ALLOCATE(tmp1_i(i),STAT=info)
                or_fail_alloc("tmp1_i")
                CALL MOVE_ALLOC(tmp1_i,nlabels)
             ENDIF
          ELSE
             IF (i.LT.j) THEN
                ALLOCATE(tmp1_i(i),STAT=info)
                or_fail_alloc("tmp1_i")
                FORALL (l=1:i) tmp1_i(l)=nlabels(l)
                CALL MOVE_ALLOC(tmp1_i,nlabels)
             ENDIF
          ENDIF

          CALL ppm_util_qsort(nlabels,info,i)
          or_fail("ppm_util_qsort")

#ifdef __MPI
          !wait for ghost_on_fire
          CALL MPI_Wait(request,MPI_STATUS_IGNORE,info)
          or_fail_MPI("MPI_Wait")
#endif

          DO WHILE (ghost_on_fire.GT.0)
             ! first step in recognizig connected regions
             CALL DTYPE(initghostfire)(MeshIn,info)
             or_fail("ghostfire")

             IF (ppm_nproc.GT.1) THEN
                !-------------------------------------------------------------------------
                !  ghost update.
                !-------------------------------------------------------------------------
                CALL labels%map_ghost_push(MeshIn,info)
                or_fail("labels%map_ghost_push")
                CALL MeshIn%map_isend(info,sendrecv=.TRUE.)
                or_fail("MeshIn%map_send")
             ENDIF

             sbpitr => MeshIn%subpatch%begin()
             isize=0
             nsize=0
             DO WHILE (ASSOCIATED(sbpitr))
                Nm => sbpitr%nnodes

#if   __DIME == __2D
                IF (sbpitr%istart(1).NE.1) THEN
                   nsize=nsize+Nm(2)
                ENDIF
                IF (sbpitr%iend(1).NE.Ngrid(1)) THEN
                   nsize=nsize+Nm(2)
                ENDIF
                IF (sbpitr%istart(2).NE.1) THEN
                   nsize=nsize+Nm(1)+2
                ENDIF
                IF (sbpitr%iend(2).NE.Ngrid(2)) THEN
                   nsize=nsize+Nm(1)+2
                ENDIF
#elif __DIME == __3D
                IF (sbpitr%istart(1).NE.1) THEN
                   nsize=nsize+Nm(2)*Nm(3)
                ENDIF
                IF (sbpitr%iend(1).NE.Ngrid(1)) THEN
                   nsize=nsize+Nm(2)*Nm(3)
                ENDIF
                IF (sbpitr%istart(2).NE.1) THEN
                   nsize=nsize+Nm(3)*(Nm(1)+2)
                ENDIF
                IF (sbpitr%iend(2).NE.Ngrid(2)) THEN
                   nsize=nsize+Nm(3)*(Nm(1)+2)
                ENDIF
                IF (sbpitr%istart(3).NE.1) THEN
                   nsize=nsize+(Nm(1)+2)*(Nm(2)+2)
                ENDIF
                IF (sbpitr%iend(3).NE.Ngrid(3)) THEN
                   nsize=nsize+(Nm(1)+2)*(Nm(2)+2)
                ENDIF
#endif

                CALL DTYPE(ppm_rc_ghost_copy)(sbpitr,old_ghost(isize+1:nsize),info)
                or_fail("ppm_rc_ghost_copy")

                isize=nsize
                sbpitr => MeshIn%subpatch%next()
             ENDDO

             IF (ppm_nproc.GT.1) THEN
                CALL MeshIn%map_isend(info,sendrecv=.FALSE.)
                or_fail("MeshIn%map_send")
                CALL labels%map_ghost_pop(MeshIn,info)
                or_fail("labels%map_ghost_pop")
             ENDIF

             sbpitr => MeshIn%subpatch%begin()
             nsize=0
             ipatch=1
             DO WHILE (ASSOCIATED(sbpitr))
                Nm => sbpitr%nnodes

                CALL sbpitr%get_field(labels,DTYPE(wpl),info)
                or_fail("Failed to get field i_wp data.")

#if   __DIME == __2D
                !-x west
                i=0
                l=nsize
                IF (sbpitr%istart(1).NE.1) THEN
                   DO j=1,Nm(2)
                      l=l+1
                      IF ((ABS(old_ghost(l)).GT.ABS(DTYPE(wpl)(i,j)) &
                      &   .AND.ABS(DTYPE(wpl)(i,j)).GT.1)) THEN
                         ALLOCATE(seed,STAT=info)
                         or_fail_alloc("seed")
                         CALL seed%add(i,j)
                         CALL ppm_rc_seeds(ipatch)%push(seed,info)
                         or_fail("could not add new seed to the collection")
                      ELSE
                         DTYPE(wpl)(i,j)=old_ghost(l)
                      ENDIF
                   ENDDO !j
                   nsize=nsize+Nm(2)
                ENDIF
                !+x east
                i=Nm(1)+1
                l=nsize
                IF (sbpitr%iend(1).NE.Ngrid(1)) THEN
                   DO j=1,Nm(2)
                      l=l+1
                      IF ((ABS(old_ghost(l)).GT.ABS(DTYPE(wpl)(i,j)) &
                      &   .AND.ABS(DTYPE(wpl)(i,j)).GT.1)) THEN
                         ALLOCATE(seed,STAT=info)
                         or_fail_alloc("seed")
                         CALL seed%add(i,j)
                         CALL ppm_rc_seeds(ipatch)%push(seed,info)
                         or_fail("could not add new seed to the collection")
                      ELSE
                         DTYPE(wpl)(i,j)=old_ghost(l)
                      ENDIF
                   ENDDO !i
                   nsize=nsize+Nm(2)
                ENDIF
                !-y south
                j=0
                l=nsize
                IF (sbpitr%istart(2).NE.1) THEN
                   i=0
                   l=l+1
                   DTYPE(wpl)(i,j)=old_ghost(l)
                   DO i=1,Nm(1)
                      l=l+1
                      IF ((ABS(old_ghost(l)).GT.ABS(DTYPE(wpl)(i,j)) &
                      &   .AND.ABS(DTYPE(wpl)(i,j)).GT.1)) THEN
                         ALLOCATE(seed,STAT=info)
                         or_fail_alloc("seed")
                         CALL seed%add(i,j)
                         CALL ppm_rc_seeds(ipatch)%push(seed,info)
                         or_fail("could not add new seed to the collection")
                      ELSE
                         DTYPE(wpl)(i,j)=old_ghost(l)
                      ENDIF
                   ENDDO !i
                   i=Nm(1)+1
                   l=l+1
                   DTYPE(wpl)(i,j)=old_ghost(l)
                   nsize=nsize+Nm(1)+2
                ENDIF
                !+y north
                j=Nm(2)+1
                l=nsize
                IF (sbpitr%iend(2).NE.Ngrid(2)) THEN
                   i=0
                   l=l+1
                   DTYPE(wpl)(i,j)=old_ghost(l)
                   DO i=1,Nm(1)
                      l=l+1
                      IF ((ABS(old_ghost(l)).GT.ABS(DTYPE(wpl)(i,j)) &
                      &   .AND.ABS(DTYPE(wpl)(i,j)).GT.1)) THEN
                         ALLOCATE(seed,STAT=info)
                         or_fail_alloc("seed")
                         CALL seed%add(i,j)
                         CALL ppm_rc_seeds(ipatch)%push(seed,info)
                         or_fail("could not add new seed to the collection")
                      ELSE
                         DTYPE(wpl)(i,j)=old_ghost(l)
                      ENDIF
                   ENDDO !i
                   i=Nm(1)+1
                   l=l+1
                   DTYPE(wpl)(i,j)=old_ghost(l)
                   nsize=nsize+Nm(1)+2
                ENDIF
#elif __DIME == __3D
                !-x wset
                i=0
                l=nsize
                IF (sbpitr%istart(1).NE.1) THEN
                   DO k=1,Nm(3)
                      DO j=1,Nm(2)
                         l=l+1
                         IF ((ABS(old_ghost(l)).GT.ABS(DTYPE(wpl)(i,j,k)) &
                         &   .AND.ABS(DTYPE(wpl)(i,j,k)).GT.1)) THEN
                            ALLOCATE(seed,STAT=info)
                            or_fail_alloc("seed")
                            CALL seed%add(i,j,k)
                            CALL ppm_rc_seeds(ipatch)%push(seed,info)
                            or_fail("could not add new seed to the collection")
                         ELSE
                            DTYPE(wpl)(i,j,k)=old_ghost(l)
                         ENDIF
                      ENDDO !j
                   ENDDO !k
                   nsize=nsize+Nm(2)*Nm(3)
                ENDIF
                !+x east
                i=Nm(1)+1
                l=nsize
                IF (sbpitr%iend(1).NE.Ngrid(1)) THEN
                   DO k=1,Nm(3)
                      DO j=1,Nm(2)
                         l=l+1
                         IF ((ABS(old_ghost(l)).GT.ABS(DTYPE(wpl)(i,j,k)) &
                         &   .AND.ABS(DTYPE(wpl)(i,j,k)).GT.1)) THEN
                            ALLOCATE(seed,STAT=info)
                            or_fail_alloc("seed")
                            CALL seed%add(i,j,k)
                            CALL ppm_rc_seeds(ipatch)%push(seed,info)
                            or_fail("could not add new seed to the collection")
                         ELSE
                            DTYPE(wpl)(i,j,k)=old_ghost(l)
                         ENDIF
                      ENDDO !j
                   ENDDO !k
                   nsize=nsize+Nm(2)*Nm(3)
                ENDIF
                !-y south
                j=0
                l=nsize
                IF (sbpitr%istart(2).NE.1) THEN
                   DO k=1,Nm(3)
                      i=0
                      l=l+1
                      DTYPE(wpl)(i,j,k)=old_ghost(l)
                      DO i=1,Nm(1)
                         l=l+1
                         IF ((ABS(old_ghost(l)).GT.ABS(DTYPE(wpl)(i,j,k)) &
                         &   .AND.ABS(DTYPE(wpl)(i,j,k)).GT.1)) THEN
                            ALLOCATE(seed,STAT=info)
                            or_fail_alloc("seed")
                            CALL seed%add(i,j,k)
                            CALL ppm_rc_seeds(ipatch)%push(seed,info)
                            or_fail("could not add new seed to the collection")
                         ELSE
                            DTYPE(wpl)(i,j,k)=old_ghost(l)
                         ENDIF
                      ENDDO !i
                      i=Nm(1)+1
                      l=l+1
                      DTYPE(wpl)(i,j,k)=old_ghost(l)
                   ENDDO !k
                   nsize=nsize+(Nm(1)+2)*Nm(3)
                ENDIF
                !+y north
                j=Nm(2)+1
                l=nsize
                IF (sbpitr%iend(2).NE.Ngrid(2)) THEN
                   DO k=1,Nm(3)
                      i=0
                      l=l+1
                      DTYPE(wpl)(i,j,k)=old_ghost(l)
                      DO i=1,Nm(1)
                         l=l+1
                         IF ((ABS(old_ghost(l)).GT.ABS(DTYPE(wpl)(i,j,k)) &
                         &   .AND.ABS(DTYPE(wpl)(i,j,k)).GT.1)) THEN
                            ALLOCATE(seed,STAT=info)
                            or_fail_alloc("seed")
                            CALL seed%add(i,j,k)
                            CALL ppm_rc_seeds(ipatch)%push(seed,info)
                            or_fail("could not add new seed to the collection")
                         ELSE
                            DTYPE(wpl)(i,j,k)=old_ghost(l)
                         ENDIF
                      ENDDO !i
                      i=Nm(1)+1
                      l=l+1
                      DTYPE(wpl)(i,j,k)=old_ghost(l)
                   ENDDO !k
                   nsize=nsize+(Nm(1)+2)*Nm(3)
                ENDIF
                !-z bottom
                k=0
                l=nsize
                IF (sbpitr%istart(3).NE.1) THEN
                   j=0
                   DO i=0,Nm(1)+1
                      l=l+1
                      DTYPE(wpl)(i,j,k)=old_ghost(l)
                   ENDDO
                   DO j=1,Nm(2)
                      i=0
                      l=l+1
                      DTYPE(wpl)(i,j,k)=old_ghost(l)
                      DO i=1,Nm(1)
                         l=l+1
                         IF ((ABS(old_ghost(l)).GT.ABS(DTYPE(wpl)(i,j,k)) &
                         &   .AND.ABS(DTYPE(wpl)(i,j,k)).GT.1)) THEN
                            ALLOCATE(seed,STAT=info)
                            or_fail_alloc("seed")
                            CALL seed%add(i,j,k)
                            CALL ppm_rc_seeds(ipatch)%push(seed,info)
                            or_fail("could not add new seed to the collection")
                         ELSE
                            DTYPE(wpl)(i,j,k)=old_ghost(l)
                         ENDIF
                      ENDDO !i
                      i=Nm(1)+1
                      l=l+1
                      DTYPE(wpl)(i,j,k)=old_ghost(l)
                   ENDDO !j
                   j=Nm(2)+1
                   DO i=0,Nm(1)+1
                      l=l+1
                      DTYPE(wpl)(i,j,k)=old_ghost(l)
                   ENDDO
                   nsize=nsize+(Nm(1)+2)*(Nm(2)+2)
                ENDIF
                !+z top
                k=Nm(3)+1
                l=nsize
                IF (sbpitr%iend(3).NE.Ngrid(3)) THEN
                   j=0
                   DO i=0,Nm(1)+1
                      l=l+1
                      DTYPE(wpl)(i,j,k)=old_ghost(l)
                   ENDDO
                   DO j=1,Nm(2)
                      i=0
                      l=l+1
                      DTYPE(wpl)(i,j,k)=old_ghost(l)
                      DO i=1,Nm(1)
                         l=l+1
                         IF ((ABS(old_ghost(l)).GT.ABS(DTYPE(wpl)(i,j,k)) &
                         &   .AND.ABS(DTYPE(wpl)(i,j,k)).GT.1)) THEN
                            ALLOCATE(seed,STAT=info)
                            or_fail_alloc("seed")
                            CALL seed%add(i,j,k)
                            CALL ppm_rc_seeds(ipatch)%push(seed,info)
                            or_fail("could not add new seed to the collection")
                         ELSE
                            DTYPE(wpl)(i,j,k)=old_ghost(l)
                         ENDIF
                      ENDDO !i
                      i=Nm(1)+1
                      l=l+1
                      DTYPE(wpl)(i,j,k)=old_ghost(l)
                   ENDDO !j
                   j=Nm(2)+1
                   DO i=0,Nm(1)+1
                      l=l+1
                      DTYPE(wpl)(i,j,k)=old_ghost(l)
                   ENDDO
                   nsize=nsize+(Nm(1)+2)*(Nm(2)+2)
                ENDIF
#endif

                sbpitr => MeshIn%subpatch%next()
                ipatch=ipatch+1
             ENDDO !WHILE (ASSOCIATED(sbpitr))

             NULLIFY(DTYPE(wpl))

             l=MeshIn%subpatch%nb
             ldu(1)=SUM(ppm_rc_seeds(1:l)%nb-seednm(1:l))

             MASK=(ldu(1).GT.0)
             ghost_on_fire=MERGE(1,0,MASK)

#ifdef __MPI
             CALL MPI_Iallreduce(MPI_IN_PLACE,ghost_on_fire,1,&
             &    MPI_INTEGER,MPI_SUM,comm,request,info)
             or_fail_MPI("MPI_Iallreduce")
#endif

             IF (ldu(1).GT.SIZE(ghost_nlabels)) THEN
                DEALLOCATE(ghost_nlabels,STAT=info)
                or_fail_dealloc("ghost_nlabels array deallocation failed!")

                ALLOCATE(ghost_nlabels(ldu(1)),STAT=info)
                or_fail_alloc("Temp array allocation failed!")
             ENDIF

#ifdef __MPI
             !wait for ghost_on_fire
             CALL MPI_Wait(request,MPI_STATUS_IGNORE,info)
             or_fail_MPI("MPI_Wait")
#endif
          ENDDO !WHILE (ghost_on_fire > 0)


          ALLOCATE(trstat,STAT=info)
          or_fail_alloc("trstat")

          CALL trstat%add(zerod,zerod,zerod,zerod)
          !add the background as a region

          CALL tmp_region_stat%push(trstat,info)
          or_fail("tmp_region_stat%push")

          nrgs=tmp_region_stat%nb
          !number of regions on each processor
          !plus background region

#ifdef __MPI
          ALLOCATE(counts(ppm_nproc),STAT=info)
          or_fail_alloc("counts")

          counts=0

          CALL MPI_Iallgather(nrgs,1,MPI_INTEGER,counts,1,MPI_INTEGER,comm,request,info)
          or_fail_MPI("MPI_Iallgather")
#endif

          !this will point to the background region which is added to the
          !stat collection
          trstat => tmp_region_stat%last()

          NULLIFY(DTYPE(wpi))

          sbpitr => MeshIn%subpatch%begin()
          DO WHILE (ASSOCIATED(sbpitr))
             Nm => sbpitr%nnodes

             CALL sbpitr%get_field(image,DTYPE(wpi),info)
             or_fail("Failed to get field i_wp data.")

             CALL sbpitr%get_field(labels,DTYPE(wpl),info)
             or_fail("Failed to get field i_wp data.")

             val=zerod
             ll=0_ppm_kind_int64

#if   __DIME == __2D
             DO j=1,Nm(2)
                DO i=1,Nm(1)
                   IF (DTYPE(wpl)(i,j).EQ.0.OR.DTYPE(wpl)(i,j).EQ.FORBIDDEN) THEN
                      ll=ll+1_ppm_kind_int64
                      dummy=REAL(DTYPE(wpi)(i,j),ppm_kind_double)
                      val(2)=val(2)+dummy
                      val(3)=val(3)+dummy*dummy
                   ENDIF
                ENDDO
             ENDDO
#elif __DIME == __3D
             DO k=1,Nm(3)
                DO j=1,Nm(2)
                   DO i=1,Nm(1)
                      IF (DTYPE(wpl)(i,j,k).EQ.0.OR.DTYPE(wpl)(i,j,k).EQ.FORBIDDEN) THEN
                         ll=ll+1_ppm_kind_int64
                         dummy=REAL(DTYPE(wpi)(i,j,k),ppm_kind_double)
                         val(2)=val(2)+dummy
                         val(3)=val(3)+dummy*dummy
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
#endif

             val(1)=REAL(ll,ppm_kind_double)
             CALL trstat%add(val)

             sbpitr => MeshIn%subpatch%next()
          ENDDO !WHILE (ASSOCIATED(sbpitr))

          NULLIFY(DTYPE(wpi),DTYPE(wpl))

          ALLOCATE(tmp_region_stat_compact(4,nrgs),STAT=info)
          or_fail_alloc("tmp_region_stat_compact")

          trstat => tmp_region_stat%begin()
          i=1
          DO WHILE (ASSOCIATED(trstat))
             value => trstat%getValue()

             tmp_region_stat_compact(1,i)=value(1)
             tmp_region_stat_compact(2,i)=value(2)
             tmp_region_stat_compact(3,i)=value(3)
             tmp_region_stat_compact(4,i)=value(4)

             trstat => tmp_region_stat%next()
             i=i+1
          ENDDO

          CALL tmp_region_stat%destroy(info)
          or_fail("tmp_region_stat%destroy")

          DEALLOCATE(tmp_region_stat,STAT=info)
          or_fail_dealloc("tmp_region_stat")

#ifdef __MPI
          !wait for ghost_on_fire
          CALL MPI_Wait(request,MPI_STATUS_IGNORE,info)
          or_fail_MPI("MPI_Wait")

          sum_nrgs=SUM(counts)

          counts=4*counts
#else
          sum_nrgs=nrgs
#endif
          ALLOCATE(tmp_region_stat_aggregate(4,sum_nrgs),STAT=info)
          or_fail_alloc("tmp_region_stat_aggregate")

#ifdef __MPI
          ALLOCATE(displ(ppm_nproc),STAT=info)
          or_fail_alloc("displ")

          displ(1)=0
          DO i=2,ppm_nproc
             displ(i)=displ(i-1)+counts(i-1)
          ENDDO

          CALL MPI_Iallgatherv(tmp_region_stat_compact,4*nrgs,MPI_DOUBLE_PRECISION, &
          &    tmp_region_stat_aggregate,counts,displ,MPI_DOUBLE_PRECISION,comm,    &
          &    request,info)
          or_fail_MPI("MPI_Iallgatherv")
#endif

          DEALLOCATE(ghost_nlabels,STAT=info)
          or_fail_dealloc("ghost_nlabels")

          DEALLOCATE(seednm,STAT=info)
          or_fail_dealloc("seednm")

          DEALLOCATE(old_ghost,STAT=info)
          or_fail_dealloc("old_ghost")

          sbpitr => MeshIn%subpatch%begin()
          ipatch=1
          DO WHILE (ASSOCIATED(sbpitr))
             CALL ppm_rc_seeds(ipatch)%destroy(info)
             or_fail("ppm_rc_seeds(ipatch)%destroy")

             sbpitr => MeshIn%subpatch%next()
             ipatch=ipatch+1
          ENDDO

#ifdef __MPI
          CALL MPI_Wait(request,MPI_STATUS_IGNORE,info)
          or_fail_MPI("MPI_Wait")

          DEALLOCATE(displ,counts,STAT=info)
          or_fail_dealloc("displ & counts")
#else
          tmp_region_stat_aggregate=tmp_region_stat_compact
#endif

          DEALLOCATE(tmp_region_stat_compact,STAT=info)
          or_fail_dealloc("tmp_region_stat_compact")

          NULLIFY(unique_new_regions)

          CALL unique(tmp_region_stat_aggregate(1,:),unique_new_regions,info)
          or_fail("unique")

          IF (ASSOCIATED(unique_new_regions)) THEN
             Nregions = SIZE(unique_new_regions)
          ELSE
             Nregions = 0
          ENDIF

          CALL e_data%galloc(Nregions-1,info)
          or_fail_alloc("e_data")
          !Minus one is for the background

          !-------------------------------------------------------------------------
          !  Create the hash table for the region labels
          !-------------------------------------------------------------------------
          ALLOCATE(htable,STAT=info)
          or_fail_alloc("Failed to allocate htable")

          ALLOCATE(MASKL(sum_nrgs),STAT=info)
          or_fail_alloc("MASKL")

          ! TODO
          htable_size = (Nregions+32) !Size of hash table

          CALL htable%create(htable_size,info)
          or_fail("create htable")

          nrgs=0

          DO i=1,Nregions
             SELECT CASE (unique_new_regions(i))
             CASE (0)
                CALL htable%insert(0_ppm_kind_int64,0,info)
                or_fail("hash insert")

                e_data%Rlabel(0)=0

                FORALL (l=1:sum_nrgs)
                   MASKL(l)=INT(tmp_region_stat_aggregate(1,l)).EQ.unique_new_regions(i)
                END FORALL

                e_data%gCount(0)=SUM(tmp_region_stat_aggregate(2,:),MASKL)
                e_data%gSums(0) =SUM(tmp_region_stat_aggregate(3,:),MASKL)
                e_data%gSumsq(0)=SUM(tmp_region_stat_aggregate(4,:),MASKL)

                nrgs=i

             CASE DEFAULT
                j=MERGE(i,i-1,nrgs.EQ.0)
                !Correcting the index, in case BG region is in the middle
                !of the unique_new_regions list

                CALL htable%insert(unique_new_regions(i),j,info)
                or_fail("hash insert")

                e_data%Rlabel(j)=unique_new_regions(i)

                FORALL (l=1:sum_nrgs)
                   MASKL(l)=INT(tmp_region_stat_aggregate(1,l)).EQ.unique_new_regions(i)
                END FORALL

                e_data%gCount(j)=SUM(tmp_region_stat_aggregate(2,:),MASKL)
                e_data%gSums(j) =SUM(tmp_region_stat_aggregate(3,:),MASKL)
                e_data%gSumsq(j)=SUM(tmp_region_stat_aggregate(4,:),MASKL)

             END SELECT
          ENDDO

#if   __DIME == __2D
          check_true(<#ABS(e_data%gCount(0)-REAL(Ngrid(1),ppm_kind_double)*REAL(Ngrid(2),ppm_kind_double)+SUM(e_data%gCount(1:Nregions-1))).LT.smalld#>, &
          & "Background pixels populations!")
#elif __DIME == __3D
!           stdout('e_data%gCount(0)')
!           stdout('REAL(Ngrid(1),ppm_kind_double)*REAL(Ngrid(2),ppm_kind_double)*REAL(Ngrid(3),ppm_kind_double)')
!           stdout('SUM(e_data%gCount(1:Nregions-1))')
!           stdout('ABS(e_data%gCount(0)-REAL(Ngrid(1),ppm_kind_double)*REAL(Ngrid(2),ppm_kind_double)*REAL(Ngrid(3),ppm_kind_double)+SUM(e_data%gCount(1:Nregions-1)))')
          check_true(<#ABS(e_data%gCount(0)-REAL(Ngrid(1),ppm_kind_double)*REAL(Ngrid(2),ppm_kind_double)*REAL(Ngrid(3),ppm_kind_double)+SUM(e_data%gCount(1:Nregions-1))).LT.smalld#>, &
          & "Background pixels populations!")
#endif

          iopt=ppm_param_dealloc
          CALL ppm_alloc(unique_new_regions,ld,iopt,info)
          or_fail_alloc("unique_new_regions")

          DEALLOCATE(tmp_region_stat_aggregate,MASKL,STAT=info)
          or_fail_dealloc("tmp_region_stat_aggregate & MASKL")

          CALL e_data%lalloc(Nregions-1,info)
          or_fail_alloc("e_data")

          !---------------------------------------------------------------------
          !  Return
          !---------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
        END SUBROUTINE DTYPE(initforestfire)

        SUBROUTINE DTYPE(forestfire)(PartIn,MeshIn,lfire,info)
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_mpi
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !------------------------------------------------------------------------
          !  Arguments
          !------------------------------------------------------------------------
          CLASS(ppm_t_particles_s_), POINTER       :: PartIn
          CLASS(ppm_t_equi_mesh_),   POINTER       :: MeshIn

          LOGICAL,                   INTENT(IN   ) :: lfire

          INTEGER,                   INTENT(  OUT) :: info

          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          CLASS(ppm_t_subpatch_), POINTER :: sbpitr

          TYPE(ppm_rc_list), POINTER :: seed
          TYPE(ppm_rc_link), POINTER :: seedlnk
          TYPE(ppm_rc_stat), POINTER :: trstat

          REAL(ppm_kind_double), DIMENSION(:), POINTER :: value
          REAL(ppm_kind_double)                        :: t0

#if   __DIME == __2D
          INTEGER, CONTIGUOUS, DIMENSION(:,:),   POINTER :: DTYPE(wpl)
          INTEGER, CONTIGUOUS, DIMENSION(:,:),   POINTER :: DTYPE(wpp)
#elif __DIME == __3D
          INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER :: DTYPE(wpl)
          INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER :: DTYPE(wpp)
#endif
          INTEGER, DIMENSION(:), ALLOCATABLE :: old_ghost
          INTEGER, DIMENSION(:),     POINTER :: seedn
          INTEGER, DIMENSION(:),     POINTER :: unique_new_regions
          INTEGER, CONTIGUOUS, DIMENSION(:,:),   POINTER :: wpl
          INTEGER, DIMENSION(__DIME+1)              :: ld
          INTEGER, DIMENSION(__DIME)              :: ll
          INTEGER                            :: iopt,isize,nsize
          INTEGER                            :: i,j,k,l,ipatch
          INTEGER                            :: iseed
          INTEGER                            :: nrgs,sum_nrgs
          INTEGER                            :: ipart,opart
          INTEGER                            :: vLabel
          INTEGER                            :: Nregions
#ifdef __MPI
          INTEGER, DIMENSION(:), ALLOCATABLE :: displ,counts
          INTEGER                            :: request
#endif

          CHARACTER(ppm_char) :: caller="forestfire"

          LOGICAL                            :: MASK
          LOGICAL, DIMENSION(:), ALLOCATABLE :: MASKL
          LOGICAL                            :: IsEnclosedByLabelFGConnectivity

          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          IF (lfire) THEN


! if (rank.eq.0) then
!              sbpitr => MeshIn%subpatch%begin()
!              DO WHILE (ASSOCIATED(sbpitr))
!                 CALL sbpitr%get_field(labels,DTYPE(wpl),info)
!                 or_fail("Failed to get field wpl data.")
! #if   __DIME == __2D
!                  stdout("---------f0--------")
!                  stdout('wpl_2d(35,35)','wpl_2d(35,55)','wpl_2d(55,35)')
!
!            DO i=0,SIZE(e_data%gCount)-1
!               IF (e_data%Rlabel(i).eq.wpl_2d(35,35).or.e_data%Rlabel(i).eq.wpl_2d(35,55).or.e_data%Rlabel(i).eq.wpl_2d(55,35)) THEN
!                  stdout('e_data%Rlabel(i)','e_data%gCount(i)','e_data%gSums(i)','e_data%gSumsq(i)')
!               ENDIF
!            ENDDO
!
!
! #endif
!              sbpitr => MeshIn%subpatch%next()
!           ENDDO !ASSOCIATED(sbpitr)
!
! else
! do i=1,1000000
! t0=log(real(i,8))
! enddo
! endif
!
! call mpi_barrier(comm,info)
! call mpi_barrier(comm,info)
! call mpi_barrier(comm,info)
! call mpi_barrier(comm,info)
! call mpi_barrier(comm,info)
! call mpi_barrier(comm,info)





             !We want to gather statistical information about regions
             !with new color
             ALLOCATE(tmp_region_stat,STAT=info)
             or_fail_alloc("tmp_region_stat")

             !an array to keep the record of the seeds and particles in main domain
             ld(1)=MeshIn%subpatch%nb
             ALLOCATE(partnm(ld(1)),STAT=info)
             or_fail_alloc("partnm")

             !if there is any firing seed, fire from the seed
             IF (ANY(ppm_rc_seeds(:)%nb.GT.0)) THEN
                ! first step in recognizig connected regions
                CALL DTYPE(fire)(MeshIn,info)
                or_fail("fire")
             ENDIF

! if (rank.eq.0) then
!              sbpitr => MeshIn%subpatch%begin()
!              DO WHILE (ASSOCIATED(sbpitr))
!                 CALL sbpitr%get_field(labels,DTYPE(wpl),info)
!                 or_fail("Failed to get field wpl data.")
! #if   __DIME == __2D
!                  stdout("---------ff0--------")
!                  stdout('wpl_2d(35,35)','wpl_2d(35,55)','wpl_2d(55,35)')
!
!
!            DO i=0,SIZE(e_data%gCount)-1
!               IF (e_data%Rlabel(i).eq.wpl_2d(35,35).or.e_data%Rlabel(i).eq.wpl_2d(35,55).or.e_data%Rlabel(i).eq.wpl_2d(55,35)) THEN
!                  stdout('e_data%Rlabel(i)','e_data%gCount(i)','e_data%gSums(i)','e_data%gSumsq(i)')
!               ENDIF
!            ENDDO
!
!
!
!
! #endif
!              sbpitr => MeshIn%subpatch%next()
!           ENDDO !ASSOCIATED(sbpitr)
!
! else
! do i=1,1000000
! t0=log(real(i,8))
! enddo
! endif
!
! call mpi_barrier(comm,info)
! call mpi_barrier(comm,info)
! call mpi_barrier(comm,info)
! call mpi_barrier(comm,info)
! call mpi_barrier(comm,info)
! call mpi_barrier(comm,info)



             GOTO 9999
          ENDIF

          ll(1)=SUM(ppm_rc_seeds(:)%nb)

          MASK=(ll(1).GT.0)
          ghost_on_fire=MERGE(1,0,MASK)

#ifdef __MPI
          CALL MPI_Iallreduce(MPI_IN_PLACE,ghost_on_fire,1, &
          &    MPI_INTEGER,MPI_SUM,comm,request,info)
          or_fail_MPI("MPI_Iallreduce")
#endif

          ll(1)=SUM(InnerContourContainer(:)%nb)-Npart
          !Not all of the ghost particles are in this list
          ALLOCATE(old_ghost(ll(1)),STAT=info)
          or_fail_alloc("old_ghost")

          NULLIFY(wpl,DTYPE(wpl),DTYPE(wpp))

          CALL PartIn%get(plabels,wpl,info,with_ghosts=.TRUE.)
          or_fail("PartIn%get for <plabels> field is failed!")

          sbpitr => MeshIn%subpatch%begin()
          ipatch=1
          opart=0
          DO WHILE (ASSOCIATED(sbpitr))
             CALL sbpitr%get_field(pind,DTYPE(wpp),info)
             or_fail("Failed to get field wpp data.")

             iseed=partnm(ipatch)+1
             !the number of real particles at this patch
             seed => InnerContourContainer(ipatch)%at(iseed)
             DO WHILE (ASSOCIATED(seed))
                seedn => seed%first%getValue()
#if   __DIME == __2D
                ipart=DTYPE(wpp)(seedn(1),seedn(2))
#elif __DIME == __3D
                ipart=DTYPE(wpp)(seedn(1),seedn(2),seedn(3))
#endif

                opart=opart+1
                old_ghost(opart)=SIGN(wpl(1,ipart),-1)
                !copy the ghost particles which have the negative labels
                !SIGN operation is to correct the hot particles labels

                seed => InnerContourContainer(ipatch)%next()
             ENDDO
             sbpitr => MeshIn%subpatch%next()
             ipatch=ipatch+1
          ENDDO !ASSOCIATED(sbpitr)
          !Copy the ghost particles in old_ghost array

          NULLIFY(DTYPE(wpp))

          CALL PartIn%set(plabels,wpl,info)
          or_fail("PartIn%set for <plabels> field is failed!")

          ghostfirenewjob=.TRUE.

#ifdef __MPI
          !wait for ghost_on_fire
          CALL MPI_Wait(request,MPI_STATUS_IGNORE,info)
          or_fail_MPI("MPI_Wait")
#endif

          DO WHILE (ghost_on_fire.GT.0)
             ! first step in recognizig connected regions
             CALL DTYPE(ghostfire)(PartIn,MeshIn,info)
             or_fail("ghostfire")

             ghostfirenewjob=.FALSE.

             CALL PartIn%get(plabels,wpl,info)
             or_fail("PartIn%get for <plabels> field is failed!")

             !TOCHECK
             !TODO
             !update particle labels with new labels from ghostfire
             sbpitr => MeshIn%subpatch%begin()
             ipatch=1
             DO WHILE (ASSOCIATED(sbpitr))
                CALL sbpitr%get_field(labels,DTYPE(wpl),info)
                or_fail("Failed to get field wpl data.")

                CALL sbpitr%get_field(pind,DTYPE(wpp),info)
                or_fail("Failed to get field wpp data.")

                seed => InnerContourContainer(ipatch)%begin()
                seed_loop: DO WHILE (ASSOCIATED(seed))
                   IF (InnerContourContainer(ipatch)%iter_id.GT.partnm(ipatch)) EXIT seed_loop

                   seedn => seed%first%getValue()

#if   __DIME == __2D
                   ipart=DTYPE(wpp)(seedn(1),seedn(2))
#elif __DIME == __3D
                   ipart=DTYPE(wpp)(seedn(1),seedn(2),seedn(3))
#endif

                   IF (ipart.LE.0) THEN
                      seed => InnerContourContainer(ipatch)%next()
                      CYCLE seed_loop
                   ENDIF

#if   __DIME == __2D
                   wpl(1,ipart)=DTYPE(wpl)(seedn(1),seedn(2))
#elif __DIME == __3D
                   wpl(1,ipart)=DTYPE(wpl)(seedn(1),seedn(2),seedn(3))
#endif

                   seed => InnerContourContainer(ipatch)%next()
                ENDDO seed_loop
                sbpitr => MeshIn%subpatch%next()
                ipatch=ipatch+1
             ENDDO

             NULLIFY(DTYPE(wpl),DTYPE(wpp))

             CALL PartIn%set(plabels,wpl,info)
             or_fail("PartIn%set for <plabels> field is failed!")

             IF (ppm_nproc.GT.1) THEN
                !-------------------------------------------------------------------------
                !  ghost update.
                !-------------------------------------------------------------------------
                CALL PartIn%map_ghost_push(info,plabels)
                or_fail("PartIn%map_ghost_push")
                CALL PartIn%map_ghost_isend(info)
                or_fail("PartIn%map_ghost_send")
                CALL PartIn%map_ghost_pop(info,plabels)
                or_fail("PartIn%map_ghost_pop")
             ENDIF
             !-------------------------------------------------------------------------
             !
             !-------------------------------------------------------------------------

             CALL PartIn%get(plabels,wpl,info,with_ghosts=.TRUE.)
             or_fail("PartIn%get for <plabels> field is failed!")

             ghost_on_fire=0

             sbpitr => MeshIn%subpatch%begin()
             ipatch=1
             opart=0
             DO WHILE (ASSOCIATED(sbpitr))
                CALL sbpitr%get_field(pind,DTYPE(wpp),info)
                or_fail("Failed to get field wpp data.")

                !Loop only through the ghost particles
                iseed=partnm(ipatch)+1
                seed => InnerContourContainer(ipatch)%at(iseed)
                DO WHILE (ASSOCIATED(seed))
                   seedn => seed%first%getValue()
#if   __DIME == __2D
                   ipart=DTYPE(wpp)(seedn(1),seedn(2))
#elif __DIME == __3D
                   ipart=DTYPE(wpp)(seedn(1),seedn(2),seedn(3))
#endif
                   opart=opart+1
                   !labels are negative!
                   IF (old_ghost(opart).NE.wpl(1,ipart)) THEN
                      ALLOCATE(seed,STAT=info)
                      or_fail_alloc("seed")

#if   __DIME == __2D
                      CALL seed%add(seedn(1),seedn(2),-wpl(1,ipart))
#elif __DIME == __3D
                      CALL seed%add(seedn(1),seedn(2),seedn(3),-wpl(1,ipart))
#endif

                      CALL ppm_rc_seeds(ipatch)%push(seed,info)
                      or_fail("could not add new seed to the collection")

                      old_ghost(opart)=wpl(1,ipart)

                      IF (ghost_on_fire.EQ.0) THEN
                         ghost_on_fire=1
#ifdef __MPI
                         CALL MPI_Iallreduce(MPI_IN_PLACE,ghost_on_fire,1,&
                         &    MPI_INTEGER,MPI_SUM,comm,request,info)
                         or_fail_MPI("MPI_Iallreduce")
#endif
                      ENDIF
                   ENDIF
                   !check the ghost labels, if there is any change which could
                   !be effective, meaning the new label has the lower value
                   !and could affect the current label, the lower value is to
                   !avoid oscillation
                   seed => InnerContourContainer(ipatch)%next()
                ENDDO
                sbpitr => MeshIn%subpatch%next()
                ipatch=ipatch+1
             ENDDO

             IF (ghost_on_fire.EQ.0) THEN
#ifdef __MPI
                CALL MPI_Iallreduce(MPI_IN_PLACE,ghost_on_fire,1,&
                &    MPI_INTEGER,MPI_SUM,comm,request,info)
                or_fail_MPI("MPI_Iallreduce")
#endif
             ENDIF

             NULLIFY(DTYPE(wpp))

             CALL PartIn%set(plabels,wpl,info)
             or_fail("PartIn%set for <plabels> field is failed!")

#ifdef __MPI
             !wait for ghost_on_fire
             CALL MPI_Wait(request,MPI_STATUS_IGNORE,info)
             or_fail_MPI("MPI_Wait")
#endif
          ENDDO !WHILE (ghost_on_fire > 0)

          !Number of regions
          nrgs=tmp_region_stat%nb

#ifdef __MPI
          ALLOCATE(counts(ppm_nproc),STAT=info)
          or_fail_alloc("counts")

          CALL MPI_Iallgather(nrgs,1,MPI_INTEGER,counts,1,MPI_INTEGER,comm,request,info)
          or_fail_MPI("MPI_Iallgather")
#endif

          !I have moved UpdateStatistics here so take advantage of
          !updating the latest changes, in case we have split or merge
          CALL e_data%UpdateStatistics(info,sendrecv=.TRUE.)
          or_fail("e_data%UpdateStatistics")

          ALLOCATE(tmp_region_stat_compact(4,nrgs),STAT=info)
          or_fail_alloc("tmp_region_stat_compact")

          trstat => tmp_region_stat%begin()
          i=1
          DO WHILE (ASSOCIATED(trstat))
             value => trstat%getValue()

             tmp_region_stat_compact(1,i)=value(1)
             tmp_region_stat_compact(2,i)=value(2)
             tmp_region_stat_compact(3,i)=value(3)
             tmp_region_stat_compact(4,i)=value(4)

             trstat => tmp_region_stat%next()
             i=i+1
          ENDDO

          DEALLOCATE(partnm,STAT=info)
          or_fail_dealloc("partnm")
          !Freeing the memory which we do not need any more

          DEALLOCATE(old_ghost,STAT=info)
          or_fail_dealloc("old_ghost")
          !Freeing the memory which we do not need any more

          sbpitr => MeshIn%subpatch%begin()
          ipatch=1
          DO WHILE (ASSOCIATED(sbpitr))
             CALL ppm_rc_seeds(ipatch)%destroy(info)
             or_fail("ppm_rc_seeds(ipatch)%destroy")

             sbpitr => MeshIn%subpatch%next()
             ipatch=ipatch+1
          ENDDO
          !Freeing the memory which we do not need any more

          CALL tmp_region_stat%destroy(info)
          or_fail("tmp_region_stat%destroy")
          !Freeing the memory which we do not need any more

          DEALLOCATE(tmp_region_stat,STAT=info)
          or_fail_dealloc("tmp_region_stat")
          !Freeing the memory which we do not need any more
          NULLIFY(tmp_region_stat)

#ifdef __MPI
          ALLOCATE(displ(ppm_nproc),STAT=info)
          or_fail_alloc("displ")

          !wait for tmp_region_stat_compact
          CALL MPI_Wait(request,MPI_STATUS_IGNORE,info)
          or_fail_MPI("MPI_Wait")

          sum_nrgs=SUM(counts)

          counts=4*counts
#else
          sum_nrgs=nrgs
#endif
          ALLOCATE(tmp_region_stat_aggregate(4,sum_nrgs),STAT=info)
          or_fail_alloc("tmp_region_stat_aggregate")

#ifdef __MPI
          displ(1)=0
          DO i=2,ppm_nproc
             displ(i)=displ(i-1)+counts(i-1)
          ENDDO

          CALL MPI_Iallgatherv(tmp_region_stat_compact,4*nrgs,MPI_DOUBLE_PRECISION, &
          &    tmp_region_stat_aggregate,counts,displ,MPI_DOUBLE_PRECISION,comm,    &
          &    request,info)
          or_fail_MPI("MPI_Iallgatherv")
#else
          tmp_region_stat_aggregate=tmp_region_stat_compact
#endif

          NULLIFY(unique_new_regions)

          IF (ppm_nproc.GT.1) THEN
             !-------------------------------------------------------------------------
             !  ghost update for labels
             ! Remark :
             !          taking advantage of computation over communication
             !          by starting the nonblocking send and recv and waiting
             !          for completion at the end
             !-------------------------------------------------------------------------
             CALL MeshIn%map_ghost_get(info)
             or_fail("MeshIn%map_ghost_get")

             CALL labels%map_ghost_push(MeshIn,info)
             or_fail("labels%map_ghost_push")
             CALL MeshIn%map_isend(info,sendrecv=.TRUE.)
             or_fail("MeshIn%map_isend")
          ENDIF

#ifdef __MPI
          !wait for tmp_region_stat_compact
          CALL MPI_Wait(request,MPI_STATUS_IGNORE,info)
          or_fail_MPI("MPI_Wait")

          DEALLOCATE(displ,counts,STAT=info)
          or_fail_dealloc("displ & counts")
#endif

          DEALLOCATE(tmp_region_stat_compact,STAT=info)
          or_fail_dealloc("tmp_region_stat_compact")

          CALL unique(tmp_region_stat_aggregate(1,:),unique_new_regions,info)
          or_fail("unique")

          IF (ASSOCIATED(unique_new_regions)) THEN
             Nregions = SIZE(unique_new_regions)
          ELSE
             Nregions = 0
          ENDIF

          !-------------------------------------------------------------------------
          !  grow the size of hash table if necessary
          !-------------------------------------------------------------------------
          nsize=SIZE(e_data%gCount)

          IF (Nregions+nsize.GT.htable%nrow) THEN
             CALL htable%grow(info)
             or_fail_alloc("Failed to grow the htable")
          ENDIF

          !NOW every processor should have regions global value available
          CALL e_data%UpdateStatistics(info,sendrecv=.FALSE.)
          or_fail("e_data%UpdateStatistics")

          isize=COUNT(e_data%gCount.LT.oned)
          IF (Nregions.GT.isize) THEN
             nsize=nsize+Nregions-isize-1
             CALL e_data%grow(nsize,info)
             or_fail("e_data%grow")
          ELSE
             !Bounds OF the array are 0:nsize-1
             nsize=nsize-1
          ENDIF

          ALLOCATE(MASKL(sum_nrgs),STAT=info)
          or_fail_alloc("MASKL")

          k=1
          DO i=1,Nregions
             DO j=k,nsize
                IF (e_data%gCount(j).LT.oned) EXIT
             ENDDO
             k=j

             CALL htable%insert(unique_new_regions(i),k,info)
             or_fail("hash insert")

             e_data%Rlabel(k)=unique_new_regions(i)

             DO l=1,sum_nrgs
                MASKL(l)=INT(tmp_region_stat_aggregate(1,l)).EQ.unique_new_regions(i)
             ENDDO

             e_data%gCount(k)=SUM(tmp_region_stat_aggregate(2,:),MASKL)
             e_data%gSums(k) =SUM(tmp_region_stat_aggregate(3,:),MASKL)
             e_data%gSumsq(k)=SUM(tmp_region_stat_aggregate(4,:),MASKL)
          ENDDO

          iopt=ppm_param_dealloc
          CALL ppm_alloc(unique_new_regions,ld,iopt,info)
          or_fail_alloc("unique_new_regions")

          DEALLOCATE(tmp_region_stat_aggregate,MASKL,STAT=info)
          or_fail_dealloc("tmp_region_stat_aggregate & MASKL")

          CALL e_data%lalloc(nsize,info)
          or_fail("e_data%lalloc")

          IF (ppm_nproc.GT.1) THEN
             CALL MeshIn%map_isend(info,sendrecv=.FALSE.)
             or_fail("MeshIn%map_isend")
             CALL labels%map_ghost_pop(MeshIn,info)
             or_fail("labels%map_ghost_pop")
          ENDIF

          sbpitr => MeshIn%subpatch%begin()
          ipatch=1
          DO WHILE (ASSOCIATED(sbpitr))
             CALL sbpitr%get_field(labels,DTYPE(wpl),info)
             or_fail("Failed to get field wpl data.")

             CALL sbpitr%get_field(pind,DTYPE(wpp),info)
             or_fail("Failed to get field wpp data.")

             seed => ppm_rc_seeds_to_remove(ipatch)%begin()

             IF (ASSOCIATED(seed)) THEN
                seedlnk => seed%first
             ELSE
                NULLIFY(seedlnk)
             ENDIF

             DO WHILE (ASSOCIATED(seedlnk))
                seedn => seedlnk%getValue()

#if   __DIME == __2D
                IF (DTYPE(wpl)(seedn(1),seedn(2)).GE.0) THEN
                   seedlnk => seedlnk%nextLink()
                   CYCLE
                ENDIF

                vLabel=ABS(DTYPE(wpl)(seedn(1),seedn(2)))
#elif __DIME == __3D
                IF (DTYPE(wpl)(seedn(1),seedn(2),seedn(3)).GE.0) THEN
                   seedlnk => seedlnk%nextLink()
                   CYCLE
                ENDIF

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
                   ipart=DTYPE(wpp)(seedn(1),seedn(2))
                   DTYPE(wpp)(seedn(1),seedn(2))=0
#elif __DIME == __3D
                   DTYPE(wpl)(seedn(1),seedn(2),seedn(3))=vLabel
                   ipart=DTYPE(wpp)(seedn(1),seedn(2),seedn(3))
                   DTYPE(wpp)(seedn(1),seedn(2),seedn(3))=0
#endif
                   IF (ipart.GT.0.AND.ipart.LE.Npart) THEN
                      del_parts=del_parts+1
                      list_del_parts(del_parts)=ipart
                   ENDIF
                ENDIF !(IsEnclosedByLabelFGConnectivity)

                seedlnk => seedlnk%nextLink()
             ENDDO !WHILE (ASSOCIATED(seedlnk))

             CALL ppm_rc_seeds_to_remove(ipatch)%destroy(info)
             or_fail("ppm_rc_seeds(ipatch)%destroy")

             sbpitr => MeshIn%subpatch%next()
             ipatch=ipatch+1
          ENDDO
          NULLIFY(DTYPE(wpl),DTYPE(wpp))

          !---------------------------------------------------------------------
          !  Return
          !---------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
        END SUBROUTINE DTYPE(forestfire)

