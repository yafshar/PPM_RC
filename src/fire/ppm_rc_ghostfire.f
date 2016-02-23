
      SUBROUTINE DTYPE(initghostfire)(MeshIn,info)

        !-------------------------------------------------------------------------
        !  Modules
        !-------------------------------------------------------------------------
        USE ppm_module_util_qsort, ONLY : ppm_util_qsort

        USE ppm_rc_module_util, ONLY : label_exist
        IMPLICIT NONE
        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        CLASS(ppm_t_equi_mesh_), POINTER       :: MeshIn

        INTEGER,                 INTENT(  OUT) :: info
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
        CLASS(ppm_t_subpatch_), POINTER :: sbpitr

        TYPE(ppm_rc_stat), POINTER :: trstat
        TYPE(ppm_rc_list), POINTER :: seed
        TYPE(ppm_rc_list), POINTER :: seedlst
        TYPE(ppm_rc_list), POINTER :: seedlsti
        TYPE(ppm_rc_link), POINTER :: seedlnk

#if   __DIME == __2D
        REAL(MK), CONTIGUOUS, DIMENSION(:,:),   POINTER :: DTYPE(wpi)
#elif __DIME == __3D
        REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER :: DTYPE(wpi)
#endif
        REAL(ppm_kind_double), DIMENSION(:),    POINTER :: value
        REAL(ppm_kind_double), DIMENSION(3)             :: val
        REAL(ppm_kind_double)                           :: t0,dummy

#if   __DIME == __2D
        INTEGER, CONTIGUOUS, DIMENSION(:,:),   POINTER :: DTYPE(wpl)
#elif __DIME == __3D
        INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER :: DTYPE(wpl)
#endif
        INTEGER, DIMENSION(:,:),   ALLOCATABLE :: labelled
        INTEGER, DIMENSION(:),     POINTER :: seedn
        INTEGER, DIMENSION(:),     POINTER :: Nm
        INTEGER, DIMENSION(__DIME)              :: ld,ld_
        INTEGER                            :: i,j,kk,nn
        INTEGER                            :: idgs,idge,lngn,ghost_nlabels_
        INTEGER                            :: nsize,ipatch,iseed

        CHARACTER(LEN=ppm_char) :: caller="initghostfire"

        LOGICAL, SAVE :: newjob=.TRUE.
        LOGICAL       :: lghost

        !-------------------------------------------------------------------------
        !  Externals
        !-------------------------------------------------------------------------

        !-------------------------------------------------------------------------
        !  Initialize
        !-------------------------------------------------------------------------
         CALL substart(caller,t0,info)

         NULLIFY(DTYPE(wpi),DTYPE(wpl))
         NULLIFY(seedlst,seedlsti)

         nsize=SIZE(nlabels,DIM=1)

         sbpitr => MeshIn%subpatch%begin()
         idgs=0
         idge=0
         ipatch=1
         patch_loop: DO WHILE (ASSOCIATED(sbpitr))
            Nm => sbpitr%nnodes

            CALL sbpitr%get_field(image,DTYPE(wpi),info)
            or_fail("Failed to get field r_wp data.")

            CALL sbpitr%get_field(labels,DTYPE(wpl),info)
            or_fail("Failed to get field i1_wp data.")

            iseed=seednm(ipatch)+1
            seed => ppm_rc_seeds(ipatch)%at(iseed)

            idgs=idge+1

            seed_loop: DO WHILE (ASSOCIATED(seed))
               seedlnk => seed%first
               seedn   => seedlnk%getValue()

               lghost=.TRUE.

               idge=idge+1

#if   __DIME == __2D
               ghost_nlabels_=ABS(DTYPE(wpl)(seedn(1),seedn(2)))
#elif __DIME == __3D
               ghost_nlabels_=ABS(DTYPE(wpl)(seedn(1),seedn(2),seedn(3)))
#endif
               ghost_nlabels(idge)=ghost_nlabels_
               lngn=idge-idgs+1

               safe_loop: DO

                  IF (.NOT.ASSOCIATED(seedlst)) THEN
                     ALLOCATE(seedlst,STAT=info)
                     or_fail_alloc("seedlst")
                  ENDIF
                  IF (.NOT.ASSOCIATED(seedlsti)) THEN
                     ALLOCATE(seedlsti,STAT=info)
                     or_fail_alloc("seedlsti")
                  ENDIF

                  DO WHILE (ASSOCIATED(seedlnk))
                     seedn => seedlnk%getValue()
#if   __DIME == __2D
                     IF (DTYPE(wpl)(seedn(1),seedn(2)).GT.0.OR.lghost) THEN
#elif __DIME == __3D
                     IF (DTYPE(wpl)(seedn(1),seedn(2),seedn(3)).GT.0.OR.lghost) THEN
#endif
                        lghost=.FALSE.

                        neigh_loop: DO kk = 1,FG_ConnectivityType%NumberOfNeighbors
                           ld=seedn(1:__DIME)+FG_ConnectivityType%NeighborsPoints(:,kk)

                           !check to see whether we are inside the domain-+1 ghost
                           main_if: IF (ALL(ld.GE.1.AND.ld.LE.Nm)) THEN
#if   __DIME == __2D
                              IF (ANY(ABS(DTYPE(wpl)(ld(1),ld(2))).EQ.ghost_nlabels(idgs:idge))) THEN
                                 IF (ABS(DTYPE(wpl)(ld(1),ld(2))).LE.ghost_nlabels(idge)) CYCLE neigh_loop
                              ENDIF

                              main_if2: IF (label_exist(ABS(DTYPE(wpl)(ld(1),ld(2))),nlabels,nsize,i)) THEN
                                 CALL DTYPE(ppm_rc_floodFill)(labelled,DTYPE(wpl),Nm,ld(1:__DIME), &
                                 &    nlabels(i),ghost_nlabels(idge),0,info)

                                 nn=0
                                 DO j=1,COUNT(labelled(1,:).GE.0)
                                    ld_=labelled(:,j)
                                    IF (ANY(ld_.LT.1.OR.ld_.GT.Nm)) CYCLE
                                    nn=nn+1
                                 ENDDO

                                 DO j=1,COUNT(labelled(1,:).GE.0)
                                    ld_=labelled(:,j)
                                    IF (ALL(ld_.GE.1.AND.ld_.LE.Nm)) EXIT
                                 ENDDO

                                 IF (ALL(ld_.GE.1.AND.ld_.LE.Nm)) THEN
                                    CALL seed%add(ld_)
                                 ENDIF

                                 DEALLOCATE(labelled,STAT=info)
                                 or_fail_dealloc("labelled")

                                 !update the statistic
                                 trstat => tmp_region_stat%begin()
                                 stat_loop: DO WHILE (ASSOCIATED(trstat))
                                    value => trstat%getValue()
                                    IF (INT(value(1)).EQ.nlabels(i).AND.INT(value(2)).EQ.nn) THEN
                                       value(1)=REAL(ghost_nlabels(idge),ppm_kind_double)
                                       nlabels(i)=ghost_nlabels(idge)

                                       CALL ppm_util_qsort(nlabels,info,nsize)
                                       or_fail("ppm_util_qsort")

                                       EXIT stat_loop
                                    ENDIF
                                    trstat => tmp_region_stat%next()
                                 ENDDO stat_loop
                                 NULLIFY(trstat)

                                 CYCLE neigh_loop
#elif __DIME == __3D
                              IF (ANY(ABS(DTYPE(wpl)(ld(1),ld(2),ld(3))).EQ.ghost_nlabels(idgs:idge))) THEN
                                 IF (ABS(DTYPE(wpl)(ld(1),ld(2),ld(3))).LE.ghost_nlabels(idge)) CYCLE neigh_loop
                              ENDIF

                              main_if2: IF (label_exist(ABS(DTYPE(wpl)(ld(1),ld(2),ld(3))),nlabels,nsize,i)) THEN
                                 CALL DTYPE(ppm_rc_floodFill)(labelled,DTYPE(wpl),Nm,ld(1:__DIME), &
                                 &    nlabels(i),ghost_nlabels(idge),0,info)

                                 nn=0
                                 DO j=1,COUNT(labelled(1,:).GE.0)
                                    ld_=labelled(:,j)
                                    IF (ANY(ld_.LT.1.OR.ld_.GT.Nm)) CYCLE
                                    nn=nn+1
                                 ENDDO

                                 DO j=1,COUNT(labelled(1,:).GE.0)
                                    ld_=labelled(:,j)
                                    IF (ALL(ld_.GE.1.AND.ld_.LE.Nm)) EXIT
                                 ENDDO

                                 IF (ALL(ld_.GE.1.AND.ld_.LE.Nm)) THEN
                                    CALL seed%add(ld_)
                                 ENDIF

                                 DEALLOCATE(labelled,STAT=info)
                                 or_fail_dealloc("labelled")

                                 !update the statistic
                                 trstat => tmp_region_stat%begin()
                                 stat_loop: DO WHILE (ASSOCIATED(trstat))
                                    value => trstat%getValue()
                                    IF (INT(value(1)).EQ.nlabels(i).AND.INT(value(2)).EQ.nn) THEN
                                       value(1)=REAL(ghost_nlabels(idge),ppm_kind_double)
                                       nlabels(i)=ghost_nlabels(idge)

                                       CALL ppm_util_qsort(nlabels,info,nsize)
                                       or_fail("ppm_util_qsort")

                                       EXIT stat_loop
                                    ENDIF
                                    trstat => tmp_region_stat%next()
                                 ENDDO stat_loop
                                 NULLIFY(trstat)

                                 CYCLE neigh_loop
#endif

                              ELSE main_if2
#if   __DIME == __2D
                                 SELECT CASE (DTYPE(wpl)(ld(1),ld(2)))
#elif __DIME == __3D
                                 SELECT CASE (DTYPE(wpl)(ld(1),ld(2),ld(3)))
#endif
                                 CASE (1)
#if   __DIME == __2D
                                    DTYPE(wpl)(ld(1),ld(2))=ghost_nlabels(idge)
#elif __DIME == __3D
                                    DTYPE(wpl)(ld(1),ld(2),ld(3))=ghost_nlabels(idge)
#endif

                                    CALL seedlst%add(ld)
                                    CALL seedlsti%add(ld)

                                    IF (ASSOCIATED(seedlsti%first,seedlsti%last)) THEN
                                       NULLIFY(trstat)
                                       ALLOCATE(trstat,STAT=info)
                                       or_fail_alloc("trstat")

                                       dummy=REAL(ghost_nlabels(idge),ppm_kind_double)

                                       val(1)=oned
#if   __DIME == __2D
                                       val(2)=REAL(DTYPE(wpi)(ld(1),ld(2)),ppm_kind_double)
#elif __DIME == __3D
                                       val(2)=REAL(DTYPE(wpi)(ld(1),ld(2),ld(3)),ppm_kind_double)
#endif
                                       val(3)=val(2)*val(2)

                                       CALL trstat%add(dummy,val(1),val(2),val(3))
                                       CALL tmp_region_stat%push(trstat,info)
                                       or_fail("tmp_region_stat%push")

                                       trstat => tmp_region_stat%last()
                                    ELSE
                                       val(1)=oned
#if   __DIME == __2D
                                       val(2)=REAL(DTYPE(wpi)(ld(1),ld(2)),ppm_kind_double)
#elif __DIME == __3D
                                       val(2)=REAL(DTYPE(wpi)(ld(1),ld(2),ld(3)),ppm_kind_double)
#endif
                                       val(3)=val(2)*val(2)

                                       CALL trstat%add(val)
                                    ENDIF
                                 CASE DEFAULT
                                    IF (ASSOCIATED(seedlsti%first)) THEN
#if   __DIME == __2D
                                       IF (DTYPE(wpl)(seedn(1),seedn(2)).NE.FORBIDDEN) THEN
                                          DTYPE(wpl)(seedn(1),seedn(2))=-ghost_nlabels(idge)
#elif __DIME == __3D
                                       IF (DTYPE(wpl)(seedn(1),seedn(2),seedn(3)).NE.FORBIDDEN) THEN
                                          DTYPE(wpl)(seedn(1),seedn(2),seedn(3))=-ghost_nlabels(idge)
#endif
                                       ENDIF
                                    ENDIF

                                 END SELECT
                              ENDIF main_if2

                           ENDIF main_if

                        ENDDO neigh_loop
                     ENDIF !(i1_wp.GT.0)

                     seedlnk => seedlnk%nextLink()
                  ENDDO !WHILE (ASSOCIATED(seedlnk))

                  IF (ASSOCIATED(seedlst%first)) THEN
                     seedlnk => seedlst%first
                     CALL seed%merge(seedlst)

                     IF (ASSOCIATED(seedlsti%first)) THEN
                        CALL seedlsti%destroy()
                        DEALLOCATE(seedlsti,STAT=info)
                        or_fail_dealloc("seedlsti")
                        NULLIFY(seedlsti)
                     ENDIF
                  ELSE
                     IF (ASSOCIATED(seedlst)) THEN
                        DEALLOCATE(seedlst,STAT=info)
                        or_fail_dealloc("seedlst")
                     ENDIF
                     IF (ASSOCIATED(seedlsti)) THEN
                        DEALLOCATE(seedlsti,STAT=info)
                        or_fail_dealloc("seedlsti")
                     ENDIF
                     NULLIFY(seedlst,seedlsti)
                     EXIT safe_loop
                  ENDIF

               ENDDO safe_loop

               seed => ppm_rc_seeds(ipatch)%next()
            ENDDO seed_loop

            iseed=seednm(ipatch)+1
            seed => ppm_rc_seeds(ipatch)%at(iseed)

            IF (newjob) THEN
               newjob=.FALSE.
               newjob_loop: DO WHILE (ASSOCIATED(seed))
                  seedlnk => seed%first
                  DO WHILE (ASSOCIATED(seedlnk))
                     seedn => seedlnk%getValue()
                     IF (ALL(seedn.GE.1.AND.seedn.LE.Nm)) THEN
                        CALL seed%swap(seedlnk)

                        seednm(ipatch)=seednm(ipatch)+1
                        iseed=SUM(seednm)

                        IF (nsize.LT.iseed) THEN
                           ld(1)=nsize+16

                           ALLOCATE (tmp1_i(ld(1)),STAT=info)
                           or_fail_alloc("Temp array allocation failed!")

                           FORALL (i=1:nsize) tmp1_i(i)=nlabels(i)

                           CALL MOVE_ALLOC(tmp1_i,nlabels)

                           nsize=ld(1)

                           FORALL (i=iseed:nsize) nlabels(i)=-bigi
                        ENDIF
#if   __DIME == __2D
                        nlabels(iseed)=ABS(DTYPE(wpl)(seedn(1),seedn(2)))
#elif __DIME == __3D
                        nlabels(iseed)=ABS(DTYPE(wpl)(seedn(1),seedn(2),seedn(3)))
#endif

                        seedlnk => seed%first
                        CALL seed%destroy(seedlnk)

                        seed => ppm_rc_seeds(ipatch)%next()
                        CYCLE newjob_loop
                     ENDIF
                     seedlnk => seedlnk%nextLink()
                  ENDDO !(ASSOCIATED(seedlnk))

                  CALL ppm_rc_seeds(ipatch)%remove(info)
                  or_fail("ppm_rc_seeds(ipatch)%remove")

                  DEALLOCATE(seed,STAT=info)
                  or_fail_dealloc("Failed to deallocate seed")

                  seed => ppm_rc_seeds(ipatch)%next()
               ENDDO newjob_loop

               DO i=nsize,1,-1
                  IF (nlabels(i).GT.0) EXIT
               ENDDO
               IF (i.EQ.0) THEN
                  IF (nsize.GT.0) THEN
                     ALLOCATE(tmp1_i(i),STAT=info)
                     or_fail_alloc("tmp1_i")
                     CALL MOVE_ALLOC(tmp1_i,nlabels)
                  ENDIF
               ELSE
                  IF (i.LT.nsize) THEN
                     ALLOCATE(tmp1_i(i),STAT=info)
                     or_fail_alloc("tmp1_i")
                     FORALL (j=1:i) tmp1_i(j)=nlabels(j)
                     CALL MOVE_ALLOC(tmp1_i,nlabels)
                  ENDIF
               ENDIF

               CALL ppm_util_qsort(nlabels,info,i)
               or_fail("ppm_util_qsort")
            ELSE
               DO WHILE (ASSOCIATED(seed))
                  CALL ppm_rc_seeds(ipatch)%remove(info)
                  or_fail("ppm_rc_seeds(ipatch)%remove")

                  DEALLOCATE(seed,STAT=info)
                  or_fail_dealloc("Failed to deallocate seed")

                  seed => ppm_rc_seeds(ipatch)%next()
               ENDDO
            ENDIF !(newjob)

            sbpitr => MeshIn%subpatch%next()
            ipatch=ipatch+1
         ENDDO patch_loop

        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
      9999 CONTINUE
        CALL substop(caller,t0,info)
      END SUBROUTINE DTYPE(initghostfire)


      SUBROUTINE DTYPE(ghostfire)(PartIn,MeshIn,info)

        USE ppm_rc_module_energy, ONLY : e_data
        IMPLICIT NONE
        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        CLASS(ppm_t_particles_s_), POINTER       :: PartIn

        CLASS(ppm_t_equi_mesh_),   POINTER       :: MeshIn

        INTEGER,                   INTENT(  OUT) :: info
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
        CLASS(ppm_t_subpatch_), POINTER :: sbpitr

        TYPE(ppm_rc_stat),  POINTER :: trstat
        TYPE(ppm_rc_stat),  POINTER :: trstat_
        TYPE(ppm_rc_list),  POINTER :: seed
        TYPE(ppm_rc_list_), POINTER :: removeseed
        TYPE(ppm_rc_list),  POINTER :: seedlst
        TYPE(ppm_rc_link_), POINTER :: seedlnk

#if   __DIME == __2D
        REAL(MK), CONTIGUOUS, DIMENSION(:,:),   POINTER :: DTYPE(wpi)
#elif __DIME == __3D
        REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER :: DTYPE(wpi)
#endif
        REAL(ppm_kind_double), DIMENSION(:),    POINTER :: value
        REAL(ppm_kind_double), DIMENSION(:),    POINTER :: value_
        REAL(ppm_kind_double), DIMENSION(3)             :: val
        REAL(ppm_kind_double)                           :: t0,dummy

#if   __DIME == __2D
        INTEGER, CONTIGUOUS, DIMENSION(:,:),   POINTER :: DTYPE(wpl)
        INTEGER, CONTIGUOUS, DIMENSION(:,:),   POINTER :: DTYPE(wpp)
#elif __DIME == __3D
        INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER :: DTYPE(wpl)
        INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER :: DTYPE(wpp)
#endif
        INTEGER, CONTIGUOUS, DIMENSION(:,:),   POINTER :: wpl
        INTEGER, DIMENSION(:),     POINTER     :: seedn
        INTEGER, DIMENSION(:),     POINTER     :: Nm
        INTEGER, DIMENSION(:,:),   ALLOCATABLE :: labelled
        INTEGER, DIMENSION(__DIME)                  :: ld,dd
        INTEGER                                :: llsum
        INTEGER                                :: i,ii
        INTEGER                                :: nsize,ipatch
        INTEGER                                :: oldlabel,newlabel
        INTEGER                                :: label_region
        INTEGER                                :: j,l,vLabel,ppart
#if   __DIME == __3D
        INTEGER                                :: k
#endif
        CHARACTER(LEN=ppm_char) :: caller='ghostfire'

        LOGICAL :: IsEnclosedByLabelFGConnectivity,lremove

        !-------------------------------------------------------------------------
        !  Externals
        !-------------------------------------------------------------------------

        !-------------------------------------------------------------------------
        !  Initialize
        !-------------------------------------------------------------------------
         CALL substart(caller,t0,info)

         NULLIFY(wpl,DTYPE(wpi),DTYPE(wpl),DTYPE(wpp))

         NULLIFY(removeseed)

         CALL PartIn%get(plabels,wpl,info)
         or_fail("PartIn%get for <plabels> field is failed!")

         sbpitr => MeshIn%subpatch%begin()
         ipatch=1
         patch_loop: DO WHILE (ASSOCIATED(sbpitr))
            CALL sbpitr%get_field(image,DTYPE(wpi),info)
            or_fail("Failed to get field wpi data.")

            CALL sbpitr%get_field(labels,DTYPE(wpl),info)
            or_fail("Failed to get field wpl data.")

            CALL sbpitr%get_field(pind,DTYPE(wpp),info)
            or_fail("Failed to get field wpp data.")

            Nm => sbpitr%nnodes

            ALLOCATE(seedlst,STAT=info)
            or_fail_alloc("seedlst")

            CALL ppm_rc_seeds_to_remove(ipatch)%push(seedlst,info)
            or_fail("ppm_rc_seeds_to_remove(ipatch)%push")

            seedlst => ppm_rc_seeds_to_remove(ipatch)%begin()

            seed => ppm_rc_seeds(ipatch)%begin()
            seed_loop: DO WHILE (ASSOCIATED(seed))
               seedn => seed%first%getValue()

               IF (ghostfirenewjob) THEN
#if   __DIME == __2D
                  newlabel=ABS(DTYPE(wpl)(seedn(1),seedn(2)))
                  ii = seedn(3)
                  !old label neighbor index number according to particle index
#elif __DIME == __3D
                  newlabel=ABS(DTYPE(wpl)(seedn(1),seedn(2),seedn(3)))
                  ii = seedn(4)
                  !old label neighbor index number according to particle index
#endif

                  ld=seedn(1:__DIME)+FG_ConnectivityType%NeighborsPoints(:,ii)
                  !coordinate of the neighbor cell which carries the old label

#if   __DIME == __2D
                  oldlabel=ABS(DTYPE(wpl)(ld(1),ld(2)))
#elif __DIME == __3D
                  oldlabel=ABS(DTYPE(wpl)(ld(1),ld(2),ld(3)))
#endif
               ELSE
#if   __DIME == __2D
                  oldlabel=ABS(DTYPE(wpl)(seedn(1),seedn(2)))
                  newlabel=seedn(3)
#elif __DIME == __3D
                  oldlabel=ABS(DTYPE(wpl)(seedn(1),seedn(2),seedn(3)))
                  newlabel=seedn(4)
#endif
               ENDIF

               IF (oldlabel.EQ.0) THEN
                  !If the oldlabel is BACKGROUND, we need to remove this seed
                  CALL ppm_rc_seeds(ipatch)%remove(info)
                  or_fail("ppm_rc_seeds(ipatch)%remove")

                  DEALLOCATE(seed,STAT=info)
                  or_fail_dealloc("Failed to deallocate seed")

                  seed => ppm_rc_seeds(ipatch)%next()
                  CYCLE seed_loop
               ENDIF

               label_region=htable%search(oldlabel)

               IF (label_region.EQ.htable_null) THEN
                  ! this region itself is a new region, it seems that two merge
                  ! or merge and split at two processors are happeneing at the
                  ! same time, so I should choose between the two of them
                  IF (oldlabel.LE.newlabel) THEN
                  !If the oldlabel has the smaller value it should be
                  !color of the region, so we need to remove this seed
                     CALL ppm_rc_seeds(ipatch)%remove(info)
                     or_fail("ppm_rc_seeds(ipatch)%remove")

                     DEALLOCATE(seed,STAT=info)
                     or_fail_dealloc("Failed to deallocate seed")

                     seed => ppm_rc_seeds(ipatch)%next()
                     CYCLE seed_loop
                  ENDIF
               ENDIF

               IF (oldlabel.EQ.newlabel) THEN
                  CALL ppm_rc_seeds(ipatch)%remove(info)
                  or_fail("ppm_rc_seeds(ipatch)%remove")

                  DEALLOCATE(seed,STAT=info)
                  or_fail_dealloc("Failed to deallocate seed")

                  seed => ppm_rc_seeds(ipatch)%next()
                  CYCLE seed_loop
               ENDIF

               CALL DTYPE(ppm_rc_floodFill)(labelled,DTYPE(wpl), &
               &    Nm,seedn(1:__DIME),oldlabel,newlabel,0,info)
               or_fail("ppm_rc_floodFill")

               ALLOCATE(trstat,STAT=info)
               or_fail_alloc("trstat")

               dummy=REAL(newlabel,ppm_kind_double)

               CALL trstat%add(dummy,zerod,zerod,zerod)

               CALL tmp_region_stat%push(trstat,info)
               or_fail("tmp_region_stat%push")

               val=zerod

               trstat => tmp_region_stat%last()

               nsize=COUNT(labelled(1,:).GE.0)

               DO i=1,nsize
                  ld=labelled(:,i)
                  IF (ANY(ld.LT.1.OR.ld.GT.Nm)) CYCLE

                  val(1)=val(1)+oned

#if   __DIME == __2D
                  dummy=REAL(DTYPE(wpi)(ld(1),ld(2)),ppm_kind_double)
#elif __DIME == __3D
                  dummy=REAL(DTYPE(wpi)(ld(1),ld(2),ld(3)),ppm_kind_double)
#endif

                  val(2)=val(2)+dummy
                  val(3)=val(3)+dummy*dummy
               ENDDO

               CALL trstat%add(val)

               value => trstat%getValue()

               IF (label_region.EQ.htable_null) THEN
                  IF (nsize.GT.0) THEN
                     lremove=.TRUE.
                  ELSE
                     lremove=.FALSE.
                  ENDIF

                  !the old label does not exist, it seems that two merge
                  !or merge and split at two processors are happeneing at the
                  !same time, so I should choose between the two of them.
                  !The oldlabel has the bigger value it should be removed from
                  !statistics and should be replaced by the newlabel color
                  llsum=MERGE(-bigi,0,INT(value(2)).EQ.0)

                  trstat_ => tmp_region_stat%begin()
                  DO WHILE (ASSOCIATED(trstat_))
                     value_ => trstat_%getValue()
                     IF (INT(value_(1)).EQ.oldlabel) THEN
                        IF (INT(value_(2)).EQ.INT(value(2))) THEN
                           lremove=.FALSE.

                           CALL tmp_region_stat%remove(info)
                           or_fail("tmp_region_stat%remove")

                           DEALLOCATE(trstat_,STAT=info)
                           or_fail_dealloc("Failed to deallocate trstat_")

                           EXIT
                        ELSE
                           llsum=llsum+INT(value_(2))
                        ENDIF
                     ENDIF
                     trstat_ => tmp_region_stat%next()
                  ENDDO !WHILE (ASSOCIATED(trstat_))

                  IF (llsum.EQ.INT(value(2))) THEN
                     trstat_ => tmp_region_stat%begin()
                     DO WHILE (ASSOCIATED(trstat_))
                        value_ => trstat_%getValue()
                        IF (INT(value_(1)).EQ.oldlabel) THEN
                           lremove=.FALSE.

                           CALL tmp_region_stat%remove(info)
                           or_fail("tmp_region_stat%remove")

                           DEALLOCATE(trstat_,STAT=info)
                           or_fail_dealloc("Failed to deallocate trstat_")
                        ENDIF
                        trstat_ => tmp_region_stat%next()
                     ENDDO !WHILE (ASSOCIATED(trstat_))
                  ENDIF !(llsum.EQ.INT(value(2)))

                  IF (lremove) THEN
                     IF (ASSOCIATED(removeseed)) THEN
                        CALL removeseed%add(oldlabel)
                     ELSE
                        ALLOCATE(removeseed,STAT=info)
                        or_fail_alloc("removeseed")

                        CALL removeseed%add(oldlabel)
                     ENDIF
                  ENDIF
               ELSE
                  e_data%lCount(label_region)=e_data%lCount(label_region)-value(2)
                  i=label_region*2
                  e_data%lSumslSumsq(i)  =e_data%lSumslSumsq(i)  -value(3)
                  e_data%lSumslSumsq(i+1)=e_data%lSumslSumsq(i+1)-value(4)
                  !update the old region statistics
               ENDIF !label_region.EQ.htable_null

               NULLIFY(trstat)

               ! Remove the particles which has changed to internal points
               ! because of merging
               DO i=1,nsize
                  ld=labelled(:,i)
                  IF (ANY(ld.LT.1.OR.ld.GT.Nm)) CYCLE
#if   __DIME == __2D
                  IF (DTYPE(wpl)(ld(1),ld(2)).GE.0) CYCLE
#elif __DIME == __3D
                  IF (DTYPE(wpl)(ld(1),ld(2),ld(3)).GE.0) CYCLE
#endif

#if   __DIME == __2D
                  vLabel=-DTYPE(wpl)(ld(1),ld(2))
#elif __DIME == __3D
                  vLabel=-DTYPE(wpl)(ld(1),ld(2),ld(3))
#endif

                  IsEnclosedByLabelFGConnectivity=.TRUE.
                  DO j=1,FG_ConnectivityType%NumberOfNeighbors
                     dd=ld+FG_ConnectivityType%NeighborsPoints(:,j)
#if   __DIME == __2D
                     IF (ABS(DTYPE(wpl)(dd(1),dd(2))).NE.vLabel) THEN
#elif __DIME == __3D
                     IF (ABS(DTYPE(wpl)(dd(1),dd(2),dd(3))).NE.vLabel) THEN
#endif
                        IsEnclosedByLabelFGConnectivity=.FALSE.
                        EXIT
                     ENDIF
                  ENDDO !j=1,FG_ConnectivityType%NumberOfNeighbors
                  IF (IsEnclosedByLabelFGConnectivity) THEN
#if   __DIME == __2D
                     DTYPE(wpl)(ld(1),ld(2))=vLabel
                     ppart=DTYPE(wpp)(ld(1),ld(2))
                     wpl(1,ppart)=-bigi
                     DTYPE(wpp)(ld(1),ld(2))=0
#elif __DIME == __3D
                     DTYPE(wpl)(ld(1),ld(2),ld(3))=vLabel
                     ppart=DTYPE(wpp)(ld(1),ld(2),ld(3))
                     wpl(1,ppart)=-bigi
                     DTYPE(wpp)(ld(1),ld(2),ld(3))=0
#endif
                     del_parts=del_parts+1
                     !ppart is definitely less than Npart as no new
                     !particle is created yet
                     list_del_parts(del_parts)=ppart
                  ELSE
                     CALL seedlst%add(ld)
                  ENDIF !(IsEnclosedByLabelFGConnectivity)
               ENDDO !i=1,nsize

               CALL ppm_rc_seeds(ipatch)%remove(info)
               or_fail("ppm_rc_seeds(ipatch)%remove")

               DEALLOCATE(seed,STAT=info)
               or_fail_dealloc("Failed to deallocate seed")

               seed => ppm_rc_seeds(ipatch)%next()
            ENDDO seed_loop

            NULLIFY(seedlst)

            sbpitr => MeshIn%subpatch%next()
            ipatch=ipatch+1
         ENDDO patch_loop

         IF (ALLOCATED(labelled)) THEN
            DEALLOCATE(labelled,STAT=info)
            or_fail_dealloc("labelled")
         ENDIF

         !TODO
         !This is very expensive, I should get rid of it
         IF (ASSOCIATED(removeseed)) THEN
            sbpitr => MeshIn%subpatch%begin()
            DO WHILE (ASSOCIATED(sbpitr))
               CALL sbpitr%get_field(image,DTYPE(wpi),info)
               or_fail("Failed to get field wpi data.")

               CALL sbpitr%get_field(labels,DTYPE(wpl),info)
               or_fail("Failed to get field wpl data.")

               Nm => sbpitr%nnodes

               seedlnk => removeseed%first
               DO WHILE (ASSOCIATED(seedlnk))
                  oldlabel=seedlnk%getValue()

                  trstat_ => tmp_region_stat%begin()
                  DO WHILE (ASSOCIATED(trstat_))
                     value_ => trstat_%getValue()
                     IF (INT(value_(1)).EQ.oldlabel) THEN
                        CALL tmp_region_stat%remove(info)
                        or_fail("tmp_region_stat%remove")

                        DEALLOCATE(trstat_,STAT=info)
                        or_fail_dealloc("Failed to deallocate trstat_")
                     ENDIF
                     trstat_ => tmp_region_stat%next()
                  ENDDO !WHILE (ASSOCIATED(trstat_))

#if   __DIME == __2D
                  IF (ANY(ABS(DTYPE(wpl)(1:Nm(1),1:Nm(2))).EQ.oldlabel)) THEN
                     NULLIFY(trstat)

                     ALLOCATE(trstat,STAT=info)
                     or_fail_alloc("trstat")

                     dummy=REAL(oldlabel,ppm_kind_double)

                     CALL trstat%add(dummy,zerod,zerod,zerod)

                     CALL tmp_region_stat%push(trstat,info)
                     or_fail("tmp_region_stat%push")

                     val=zerod

                     trstat => tmp_region_stat%last()

                     l=0

                     DO j=1,Nm(2)
                        DO i=1,Nm(1)
                           IF (ABS(DTYPE(wpl)(i,j)).EQ.oldlabel) THEN
                              l=l+1
                              dummy=REAL(DTYPE(wpi)(i,j),ppm_kind_double)
                              val(2)=val(2)+dummy
                              val(3)=val(3)+dummy*dummy
                           ENDIF
                        ENDDO
                     ENDDO

                     val(1)=REAL(l,ppm_kind_double)

                     CALL trstat%add(val)
                  ENDIF
#elif __DIME == __3D
                  IF (ANY(ABS(DTYPE(wpl)(1:Nm(1),1:Nm(2),1:Nm(3))).EQ.oldlabel)) THEN
                     NULLIFY(trstat)

                     ALLOCATE(trstat,STAT=info)
                     or_fail_alloc("trstat")

                     dummy=REAL(oldlabel,ppm_kind_double)

                     CALL trstat%add(dummy,zerod,zerod,zerod)

                     CALL tmp_region_stat%push(trstat,info)
                     or_fail("tmp_region_stat%push")

                     val=zerod

                     trstat => tmp_region_stat%last()

                     l=0

                     DO k=1,Nm(3)
                        DO j=1,Nm(2)
                           DO i=1,Nm(1)
                              IF (ABS(DTYPE(wpl)(i,j,k)).EQ.oldlabel) THEN
                                 l=l+1
                                 dummy=REAL(DTYPE(wpi)(i,j,k),ppm_kind_double)
                                 val(2)=val(2)+dummy
                                 val(3)=val(3)+dummy*dummy
                              ENDIF
                           ENDDO
                        ENDDO
                     ENDDO

                     val(1)=REAL(l,ppm_kind_double)

                     CALL trstat%add(val)
                  ENDIF
#endif
                  seedlnk => seedlnk%nextLink()
               ENDDO !WHILE (ASSOCIATED(seedlnk))

               sbpitr => MeshIn%subpatch%next()
            ENDDO

            CALL removeseed%destroy()

            DEALLOCATE(removeseed,STAT=info)
            or_fail_dealloc("removeseed")

            NULLIFY(removeseed)
         ENDIF !ASSOCIATED(removeseed)

        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
      9999 CONTINUE
        CALL substop(caller,t0,info)
      END SUBROUTINE DTYPE(ghostfire)

