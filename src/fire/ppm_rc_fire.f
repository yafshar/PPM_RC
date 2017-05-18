
      SUBROUTINE DTYPE(initfire)(MeshIn,info)

        USE ppm_rc_module_util, ONLY : ppm_rc_label_exist
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
        TYPE(ppm_rc_link), POINTER :: seedlnk

#if   __DIME == __2D
        REAL(MK), CONTIGUOUS, DIMENSION(:,:),   POINTER :: wpi
#elif __DIME == __3D
        REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER :: wpi
#endif
        REAL(ppm_kind_double), DIMENSION(:),     POINTER :: value
        REAL(ppm_kind_double), DIMENSION(3)              :: val
        REAL(ppm_kind_double)                            :: t0,dummy

#if   __DIME == __2D
        INTEGER, CONTIGUOUS, DIMENSION(:,:), POINTER :: wpl
#elif __DIME == __3D
        INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER :: wpl
#endif
        INTEGER, DIMENSION(:),     POINTER :: seedn
        INTEGER, DIMENSION(:),     POINTER :: Nm
        INTEGER, DIMENSION(__DIME)              :: ld,ld_
        INTEGER                            :: kk,ipatch,iseed
        INTEGER                            :: nsize,nrg,arrayindex
        INTEGER                            :: clabel
        INTEGER                            :: i,j,z
#if   __DIME == __3D
        INTEGER                            :: k
#endif

        CHARACTER(LEN=ppm_char) :: caller='initfire'

        LOGICAL :: IsSeedForgotten

        !-------------------------------------------------------------------------
        !  Externals
        !-------------------------------------------------------------------------

        !-------------------------------------------------------------------------
        !  Initialize
        !-------------------------------------------------------------------------
         CALL substart(caller,t0,info)

         NULLIFY(seedlst,trstat)
         NULLIFY(wpi,wpl)

         sbpitr => MeshIn%subpatch%begin()
         ipatch=1
         iseed=1
         patch_loop: DO WHILE (ASSOCIATED(sbpitr))
            Nm => sbpitr%nnodes

            CALL sbpitr%get_field(image,wpi,info)
            or_fail("Failed to get field r_wp data.")

            CALL sbpitr%get_field(labels,wpl,info)
            or_fail("Failed to get field i1_wp data.")

            seed => ppm_rc_seeds(ipatch)%begin()
            seed_loop: DO WHILE (ASSOCIATED(seed))
               seedlnk => seed%first
               seedn   => seedlnk%getValue()

#if   __DIME == __2D
               IF (ppm_rc_label_exist(ABS(wpl(seedn(1),seedn(2))),nlabels,iseed-1,.TRUE.)) THEN
#elif __DIME == __3D
               IF (ppm_rc_label_exist(ABS(wpl(seedn(1),seedn(2),seedn(3))),nlabels,iseed-1,.TRUE.)) THEN
#endif
                   CALL ppm_rc_seeds(ipatch)%remove(info)
                   or_fail("ppm_rc_seeds(ipatch)%remove")

                   DEALLOCATE(seed,STAT=info)
                   or_fail_dealloc("Failed to deallocate seed")

                   seed => ppm_rc_seeds(ipatch)%next()
                   CYCLE seed_loop
               ENDIF

#if   __DIME == __2D
               IF (wpl(seedn(1),seedn(2)).NE.FORBIDDEN) THEN
                  wpl(seedn(1),seedn(2))=nlabels(iseed)
#elif __DIME == __3D
               IF (wpl(seedn(1),seedn(2),seedn(3)).NE.FORBIDDEN) THEN
                  wpl(seedn(1),seedn(2),seedn(3))=nlabels(iseed)
#endif
                  val(1)=REAL(nlabels(iseed),ppm_kind_double)
                  val(2)=oned
#if   __DIME == __2D
                  val(3)=REAL(wpi(seedn(1),seedn(2)),ppm_kind_double)
#elif __DIME == __3D
                  val(3)=REAL(wpi(seedn(1),seedn(2),seedn(3)),ppm_kind_double)
#endif
               ELSE
                  !It is possible that the seed is in the padded layer at image
                  !border and should not be forgotten
                  IsSeedForgotten=.FALSE.
                  DO kk = 1,FG_ConnectivityType%NumberOfNeighbors
                     ld=seedn(1:__DIME)+FG_ConnectivityType%NeighborsPoints(:,kk)
                     IF (ALL(ld.GE.1.AND.ld.LE.Nm)) THEN
#if   __DIME == __2D
                        IF (wpl(ld(1),ld(2)).EQ.1) THEN
#elif __DIME == __3D
                        IF (wpl(ld(1),ld(2),ld(3)).EQ.1) THEN
#endif
                           IsSeedForgotten=.TRUE.
                           EXIT
                        ENDIF
                     ENDIF
                  ENDDO !j=1,FG_ConnectivityType%NumberOfNeighbors

                  IF (IsSeedForgotten) THEN
                     val(1)=REAL(nlabels(iseed),ppm_kind_double)
                     val(2)=zerod
                     val(3)=zerod
                  ELSE
                     CALL ppm_rc_seeds(ipatch)%remove(info)
                     or_fail("ppm_rc_seeds(ipatch)%remove")

                     DEALLOCATE(seed,STAT=info)
                     or_fail_dealloc("Failed to deallocate seed")

                     seed => ppm_rc_seeds(ipatch)%next()
                     CYCLE seed_loop
                  ENDIF
               ENDIF

               IF (ASSOCIATED(trstat)) NULLIFY(trstat)
               ALLOCATE(trstat,STAT=info)
               or_fail_alloc("trstat")

               dummy=val(3)*val(3)

               CALL trstat%add(val(1),val(2),val(3),dummy)
               CALL tmp_region_stat%push(trstat,info)
               or_fail("tmp_region_stat%push")

               trstat => tmp_region_stat%last()

               safe_loop: DO

                  IF (.NOT.ASSOCIATED(seedlst)) THEN
                     ALLOCATE(seedlst,STAT=info)
                     or_fail_alloc("seedlst")
                  ENDIF

                  DO WHILE (ASSOCIATED(seedlnk))
                     seedn => seedlnk%getValue()
                     IF (ALL(seedn.GE.1.AND.seedn.LE.Nm)) THEN
#if   __DIME == __2D
                     IF (wpl(seedn(1),seedn(2)).NE.0) THEN
#elif __DIME == __3D
                     IF (wpl(seedn(1),seedn(2),seedn(3)).NE.0) THEN
#endif
                        DO kk = 1,FG_ConnectivityType%NumberOfNeighbors
                           ld=seedn(1:__DIME)+FG_ConnectivityType%NeighborsPoints(:,kk)

                           !check to see whether we are inside the domain-+1 ghost
                           IF (ALL(ld.GE.0.AND.ld.LE.Nm+1)) THEN
#if   __DIME == __2D
                              IF (ABS(wpl(ld(1),ld(2))).EQ.nlabels(iseed)) CYCLE
                              IF (ppm_rc_label_exist(ABS(wpl(seedn(1),seedn(2))),nlabels,iseed-1,.TRUE.)) THEN
                                 IF (wpl(seedn(1),seedn(2)).NE.FORBIDDEN) THEN
                                    wpl(seedn(1),seedn(2))=-nlabels(iseed)
                                 ENDIF
#elif __DIME == __3D
                              IF (ABS(wpl(ld(1),ld(2),ld(3))).EQ.nlabels(iseed)) CYCLE
                              IF (ppm_rc_label_exist(ABS(wpl(seedn(1),seedn(2),seedn(3))),nlabels,iseed-1,.TRUE.)) THEN
                                 IF (wpl(seedn(1),seedn(2),seedn(3)).NE.FORBIDDEN) THEN
                                    wpl(seedn(1),seedn(2),seedn(3))=-nlabels(iseed)
                                 ENDIF
#endif
                                 CYCLE
                              ENDIF
#if   __DIME == __2D
                              SELECT CASE (wpl(ld(1),ld(2)))
#elif __DIME == __3D
                              SELECT CASE (wpl(ld(1),ld(2),ld(3)))
#endif
                              CASE (1)
#if   __DIME == __2D
                                 wpl(ld(1),ld(2))=nlabels(iseed)
#elif __DIME == __3D
                                 wpl(ld(1),ld(2),ld(3))=nlabels(iseed)
#endif
                                 CALL seedlst%add(ld)

                                 IF (ALL(ld.GE.1.AND.ld.LE.Nm)) THEN
                                    val(1)=oned
#if   __DIME == __2D
                                    val(2)=REAL(wpi(ld(1),ld(2)),ppm_kind_double)
#elif __DIME == __3D
                                    val(2)=REAL(wpi(ld(1),ld(2),ld(3)),ppm_kind_double)
#endif
                                    val(3)=val(2)*val(2)

                                    CALL trstat%add(val)
                                 ENDIF

                              !What ever else the original seed would be on the
                              !region border
                              CASE DEFAULT
#if   __DIME == __2D
                                 IF (wpl(seedn(1),seedn(2)).NE.FORBIDDEN) THEN
                                     wpl(seedn(1),seedn(2))=-nlabels(iseed)
#elif __DIME == __3D
                                 IF (wpl(seedn(1),seedn(2),seedn(3)).NE.FORBIDDEN) THEN
                                     wpl(seedn(1),seedn(2),seedn(3))=-nlabels(iseed)
#endif
                                 ENDIF

                              END SELECT

                           ENDIF !(ALL(ld.GE.0).AND.ALL(ld.LE.Nm+1))

                        ENDDO !kk = 1,
                     ENDIF !(i1_wp_DTYPE(__DIME).GT.0)
                     ENDIF !(ALL(seedn.GE.1.AND.seedn.LE.Nm)) THEN
                     seedlnk => seedlnk%nextLink()
                  ENDDO !WHILE (ASSOCIATED(seedlnk))

                  IF (ASSOCIATED(seedlst%first)) THEN
                     seedlnk => seedlst%first
                     CALL seed%merge(seedlst)
                  ELSE
                     IF (ASSOCIATED(seedlst)) THEN
                        DEALLOCATE(seedlst,STAT=info)
                        or_fail_dealloc("seedlst")
                     ENDIF
                     NULLIFY(seedlst)
                     EXIT safe_loop
                  ENDIF

               ENDDO safe_loop

               seedlnk => seed%first
               CALL seed%destroy(seedlnk)

               seed => ppm_rc_seeds(ipatch)%next()
               iseed=iseed+1
            ENDDO seed_loop

            seednm(ipatch)=ppm_rc_seeds(ipatch)%nb

            sbpitr => MeshIn%subpatch%next()
            ipatch=ipatch+1
         ENDDO patch_loop

         NULLIFY(seed)
         NULLIFY(wpi,wpl)

         SELECT CASE (vInitKind)
         CASE (e_fromfile)
            nsize=SIZE(nlabels,DIM=1)
            check_true(<#nsize.EQ.0#>,"Error SIZE(nlabels)/=0!")

            nrg=0

            sbpitr => MeshIn%subpatch%begin()
            ipatch=1
            DO WHILE (ASSOCIATED(sbpitr))
               Nm => sbpitr%nnodes

               seednm(ipatch)=0

               CALL sbpitr%get_field(image,wpi,info)
               or_fail("Failed to get field r_wp data.")

               CALL sbpitr%get_field(labels,wpl,info)
               or_fail("Failed to get field i1_wp data.")

#if   __DIME == __2D
               DO j=1,Nm(2)
                  DO i=1,Nm(1)
                     ! Current label
                     clabel=wpl(i,j)

                     IF (clabel.GT.0.AND.clabel.NE.FORBIDDEN) THEN
                        IF (ppm_rc_label_exist(clabel,nlabels,nrg)) THEN
                           val(1)=oned
                           val(2)=REAL(wpi(i,j),ppm_kind_double)
                           val(3)=val(2)*val(2)

                           IF (clabel.EQ.INT(value(1))) THEN
                              CALL trstat%add(val)
                           ELSE
                              trstat => tmp_region_stat%last()
                              trstat_loop: DO WHILE (ASSOCIATED(trstat))
                                 value => trstat%getValue()
                                 IF (clabel.EQ.INT(value(1))) THEN
                                    CALL trstat%add(val)
                                    EXIT trstat_loop
                                 ENDIF
                                 trstat => tmp_region_stat%prev()
                              ENDDO trstat_loop
                           ENDIF
                        ELSE
                           arrayindex=ppm_rc_label_index(clabel,nlabels,nrg)

                           ! New label
                           nrg=nrg+1
                           IF (nrg.GT.nsize) THEN
                              ALLOCATE(tmp1_i(nrg*2),STAT=info)
                              or_fail_alloc("Tmp array allocation failed!")

                              IF (arrayindex.EQ.0) THEN
                                 tmp1_i(nrg)=clabel
                              ELSE
                                 FORALL (z=1:arrayindex-1) tmp1_i(z)=nlabels(z)
                                 tmp1_i(arrayindex)=clabel
                                 FORALL (z=arrayindex:nsize) tmp1_i(z+1)=nlabels(z)
                              ENDIF

                              nsize=nrg*2

                              CALL MOVE_ALLOC(tmp1_i,nlabels)
                           ELSE
                              IF      (arrayindex.EQ.nrg) THEN
                                 nlabels(nrg)=clabel
                              ELSE IF (arrayindex.GT.1) THEN
                                 DO z=nrg,arrayindex+1,-1
                                    nlabels(z)=nlabels(z-1)
                                 ENDDO
                                 nlabels(arrayindex)=clabel
                              ELSE IF (arrayindex.EQ.1) THEN
                                 DO z=nrg,2,-1
                                    nlabels(z)=nlabels(z-1)
                                 ENDDO
                                 nlabels(1)=clabel
                              ELSE IF (arrayindex.EQ.0) THEN
                                 nlabels(nrg)=clabel
                              ENDIF
                           ENDIF


                           val(1)=REAL(clabel,ppm_kind_double)
                           val(2)=oned
                           val(3)=REAL(wpi(i,j),ppm_kind_double)
                           dummy=val(3)*val(3)

                           NULLIFY(trstat)
                           ALLOCATE(trstat,STAT=info)
                           or_fail_alloc("trstat")
                           CALL trstat%add(val(1),val(2),val(3),dummy)
                           CALL tmp_region_stat%push(trstat,info)
                           or_fail("tmp_region_stat%push")

                           trstat => tmp_region_stat%last()
                           value => trstat%getValue()
                        ENDIF

                        DO kk =1,FG_ConnectivityType%NumberOfNeighbors
                           ld(1)=i+FG_ConnectivityType%NeighborsPoints(1,kk)
                           ld(2)=j+FG_ConnectivityType%NeighborsPoints(2,kk)
                           IF (ABS(wpl(ld(1),ld(2))).NE.clabel) THEN
                               wpl(i,j)=-clabel
                               EXIT
                           ENDIF
                        ENDDO !kk =1,FG_ConnectivityType%NumberOfNeighbors
                     ENDIF !clabel.GT.0.AND.clabel.NE.FORBIDDEN
                  ENDDO !i=1,Nm(1)
               ENDDO !j=1,Nm(2)
#elif __DIME == __3D
               DO k=1,Nm(3)
                  DO j=1,Nm(2)
                     DO i=1,Nm(1)
                        ! Current label
                        clabel=wpl(i,j,k)

                        IF (clabel.GT.0.AND.clabel.NE.FORBIDDEN) THEN
                           IF (ppm_rc_label_exist(clabel,nlabels,nrg)) THEN
                              val(1)=oned
                              val(2)=REAL(wpi(i,j,k),ppm_kind_double)
                              val(3)=val(2)*val(2)

                              IF (clabel.EQ.INT(value(1))) THEN
                                 CALL trstat%add(val)
                              ELSE
                                 trstat => tmp_region_stat%last()
                                 trstat_loop: DO WHILE (ASSOCIATED(trstat))
                                    value => trstat%getValue()
                                    IF (clabel.EQ.INT(value(1))) THEN
                                       CALL trstat%add(val)
                                       EXIT trstat_loop
                                    ENDIF
                                    trstat => tmp_region_stat%prev()
                                 ENDDO trstat_loop
                              ENDIF
                           ELSE
                              arrayindex=ppm_rc_label_index(clabel,nlabels,nrg)

                              ! New label
                              nrg=nrg+1
                              IF (nrg.GT.nsize) THEN
                                 ALLOCATE(tmp1_i(nrg*2),STAT=info)
                                 or_fail_alloc("Tmp array allocation failed!")

                                 IF (arrayindex.EQ.0) THEN
                                    tmp1_i(nrg)=clabel
                                 ELSE
                                    FORALL (z=1:arrayindex-1) tmp1_i(z)=nlabels(z)
                                    tmp1_i(arrayindex)=clabel
                                    FORALL (z=arrayindex:nsize) tmp1_i(z+1)=nlabels(z)
                                 ENDIF

                                 nsize=nrg*2

                                 CALL MOVE_ALLOC(tmp1_i,nlabels)
                              ELSE
                                 IF      (arrayindex.EQ.nrg) THEN
                                    nlabels(nrg)=clabel
                                 ELSE IF (arrayindex.GT.1) THEN
                                    DO z=nrg,arrayindex+1,-1
                                       nlabels(z)=nlabels(z-1)
                                    ENDDO
                                    nlabels(arrayindex)=clabel
                                 ELSE IF (arrayindex.EQ.1) THEN
                                    DO z=nrg,2,-1
                                       nlabels(z)=nlabels(z-1)
                                    ENDDO
                                    nlabels(1)=clabel
                                 ELSE IF (arrayindex.EQ.0) THEN
                                    nlabels(nrg)=clabel
                                 ENDIF
                              ENDIF


                              val(1)=REAL(clabel,ppm_kind_double)
                              val(2)=oned
                              val(3)=REAL(wpi(i,j,k),ppm_kind_double)
                              dummy=val(3)*val(3)

                              NULLIFY(trstat)
                              ALLOCATE(trstat,STAT=info)
                              or_fail_alloc("trstat")
                              CALL trstat%add(val(1),val(2),val(3),dummy)
                              CALL tmp_region_stat%push(trstat,info)
                              or_fail("tmp_region_stat%push")

                              trstat => tmp_region_stat%last()
                              value => trstat%getValue()
                           ENDIF

                           DO kk =1,FG_ConnectivityType%NumberOfNeighbors
                              ld(1)=i+FG_ConnectivityType%NeighborsPoints(1,kk)
                              ld(2)=j+FG_ConnectivityType%NeighborsPoints(2,kk)
                              ld(3)=k+FG_ConnectivityType%NeighborsPoints(3,kk)
                              IF (ABS(wpl(ld(1),ld(2),ld(3))).NE.clabel) THEN
                                 wpl(i,j,k)=-clabel
                                 EXIT
                              ENDIF
                           ENDDO !kk =1,FG_ConnectivityType%NumberOfNeighbors
                        ENDIF !clabel.GT.0.AND.clabel.NE.FORBIDDEN
                     ENDDO !i=1,Nm(1)
                  ENDDO !j=1,Nm(2)
               ENDDO !k=1,Nm(3)
#endif

               sbpitr => MeshIn%subpatch%next()
               ipatch=ipatch+1
            ENDDO !WHILE (ASSOCIATED(sbpitr))

            FORALL (i=nrg+1:nsize) nlabels(i)=-bigi

            IF (ppm_rc_label_exist(1,nlabels,nrg)) THEN
               stdout("In initializiation from a file there is a region with label intensity of 1.")
               fail("Region with label of 1! Label 1 is not allowed! (Please at first use createuniquelabel)", &
               & ppm_error=ppm_error_fatal)
            ENDIF

         CASE DEFAULT
            sbpitr => MeshIn%subpatch%begin()
            ipatch=1
            DO WHILE (ASSOCIATED(sbpitr))
               Nm => sbpitr%nnodes

               CALL sbpitr%get_field(image,wpi,info)
               or_fail("Failed to get field r_wp data.")

               CALL sbpitr%get_field(labels,wpl,info)
               or_fail("Failed to get field i1_wp data.")

               ld_=1
#if   __DIME == __2D
               DO WHILE (ANY(wpl(1:Nm(1),ld_(2):Nm(2)).EQ.1))
                  j_loop: DO j=ld_(2),Nm(2)
                     DO i=1,Nm(1)
                        IF (wpl(i,j).EQ.1) EXIT j_loop
                     ENDDO
                  ENDDO j_loop
                  ld_(1)=i
                  ld_(2)=j
#elif __DIME == __3D
               DO WHILE (ANY(wpl(1:Nm(1),1:Nm(2),ld_(3):Nm(3)).EQ.1))
                  k_loop: DO k=ld_(3),Nm(3)
                     DO j=1,Nm(2)
                        DO i=1,Nm(1)
                           IF (wpl(i,j,k).EQ.1) EXIT k_loop
                        ENDDO
                     ENDDO
                  ENDDO k_loop
                  ld_(1)=i
                  ld_(2)=j
                  ld_(3)=k
#endif

                  ALLOCATE(seed,STAT=info)
                  or_fail_alloc("seed")
                  CALL seed%add(ld_)
                  CALL ppm_rc_seeds(ipatch)%push(seed,info)
                  or_fail("could not add new seed to the collection")

                  nsize=SIZE(nlabels,DIM=1)
                  seednm(ipatch)=seednm(ipatch)+1
                  !one seed has been added
                  iseed=SUM(seednm)

                  IF (nsize.LT.iseed) THEN
                     ld(1)=nsize+16

                     ALLOCATE(tmp1_i(ld(1)),STAT=info)
                     or_fail_alloc("Temp array allocation failed!")

                     FORALL (i=1:nsize) tmp1_i(i)=nlabels(i)

                     CALL MOVE_ALLOC(tmp1_i,nlabels)

                     nsize=ld(1)

                     FORALL (i=iseed:nsize) nlabels(i)=-bigi
                  ENDIF

                  IF (nlabels(iseed).EQ.-bigi) THEN
                     nlabels(iseed)=loc_label
                     loc_label=loc_label-1
                  ENDIF

#if   __DIME == __2D
                  wpl(ld_(1),ld_(2))=nlabels(iseed)
#elif __DIME == __3D
                  wpl(ld_(1),ld_(2),ld_(3))=nlabels(iseed)
#endif

                  ALLOCATE(trstat,STAT=info)
                  or_fail_alloc("trstat")

                  val(1)=REAL(nlabels(iseed),ppm_kind_double)
                  val(2)=oned
#if   __DIME == __2D
                  val(3)=REAL(wpi(ld_(1),ld_(2)),ppm_kind_double)
#elif __DIME == __3D
                  val(3)=REAL(wpi(ld_(1),ld_(2),ld_(3)),ppm_kind_double)
#endif
                  dummy=val(3)*val(3)

                  CALL trstat%add(val(1),val(2),val(3),dummy)
                  CALL tmp_region_stat%push(trstat,info)
                  or_fail("tmp_region_stat%push")

                  trstat => tmp_region_stat%last()

                  seed    => ppm_rc_seeds(ipatch)%last()
                  seedlnk => seed%first

                  safe_loop_: DO

                     IF (.NOT.ASSOCIATED(seedlst)) THEN
                        ALLOCATE(seedlst,STAT=info)
                        or_fail_alloc("seedlst")
                     ENDIF

                     DO WHILE (ASSOCIATED(seedlnk))
                        seedn => seedlnk%getValue()
                        IF (ALL(seedn.GE.1.AND.seedn.LE.Nm)) THEN
#if   __DIME == __2D
                        IF (wpl(seedn(1),seedn(2)).NE.0) THEN
#elif __DIME == __3D
                        IF (wpl(seedn(1),seedn(2),seedn(3)).NE.0) THEN
#endif
                           DO kk = 1,FG_ConnectivityType%NumberOfNeighbors
                              ld=seedn(1:__DIME)+FG_ConnectivityType%NeighborsPoints(:,kk)

                              !check to see whether we are inside the domain-+1 ghost
                              IF (ALL(ld.GE.0.AND.ld.LE.Nm+1)) THEN
#if   __DIME == __2D
                                 IF (ABS(wpl(ld(1),ld(2))).EQ.nlabels(iseed)) CYCLE
                                 IF (ppm_rc_label_exist(ABS(wpl(ld(1),ld(2))),nlabels,iseed-1,.TRUE.)) THEN
                                    wpl(seedn(1),seedn(2))=-nlabels(iseed)
#elif __DIME == __3D
                                 IF (ABS(wpl(ld(1),ld(2),ld(3))).EQ.nlabels(iseed)) CYCLE
                                 IF (ppm_rc_label_exist(ABS(wpl(ld(1),ld(2),ld(3))),nlabels,iseed-1,.TRUE.)) THEN
                                    wpl(seedn(1),seedn(2),seedn(3))=-nlabels(iseed)
#endif
                                    CYCLE
                                 ENDIF

#if   __DIME == __2D
                                 SELECT CASE (wpl(ld(1),ld(2)))
#elif __DIME == __3D
                                 SELECT CASE (wpl(ld(1),ld(2),ld(3)))
#endif
                                 CASE (1)
#if   __DIME == __2D
                                    wpl(ld(1),ld(2))=nlabels(iseed)
#elif __DIME == __3D
                                    wpl(ld(1),ld(2),ld(3))=nlabels(iseed)
#endif
                                    CALL seedlst%add(ld)

                                    IF (ALL(ld.GE.1.AND.ld.LE.Nm)) THEN
                                       val(1)=oned
#if   __DIME == __2D
                                       val(2)=REAL(wpi(ld(1),ld(2)),ppm_kind_double)
#elif __DIME == __3D
                                       val(2)=REAL(wpi(ld(1),ld(2),ld(3)),ppm_kind_double)
#endif
                                       val(3)=val(2)*val(2)

                                       CALL trstat%add(val)
                                    ENDIF

                                 CASE DEFAULT
#if   __DIME == __2D
                                    wpl(seedn(1),seedn(2))=-nlabels(iseed)
#elif __DIME == __3D
                                    wpl(seedn(1),seedn(2),seedn(3))=-nlabels(iseed)
#endif

                                 END SELECT

                              ENDIF !(ALL(ld.GE.0.AND.ld.LE.Nm+1))

                           ENDDO !kk = 1,
                        ENDIF !wpl(seedn(1),seedn(2)).GT.0
                        ENDIF !(ALL(seedn.GE.1.AND.seedn.LE.Nm)) THEN
                        seedlnk => seedlnk%nextLink()
                     ENDDO !WHILE (ASSOCIATED(seedlnk))

                     IF (ASSOCIATED(seedlst%first)) THEN
                        seedlnk => seedlst%first
                        CALL seed%merge(seedlst)
                     ELSE
                        IF (ASSOCIATED(seedlst)) THEN
                           DEALLOCATE(seedlst,STAT=info)
                           or_fail_dealloc("seedlst")
                        ENDIF
                        NULLIFY(seedlst)
                        EXIT safe_loop_
                     ENDIF

                  ENDDO safe_loop_

                  seedlnk => seed%first
                  CALL seed%destroy(seedlnk)

               ENDDO ! WHILE (ANY(wpl(1:Nm(1),ld_(2):Nm(2)).EQ.1))

               sbpitr => MeshIn%subpatch%next()
               ipatch=ipatch+1
            ENDDO !WHILE (ASSOCIATED(sbpitr))

         END SELECT

         NULLIFY(wpi,wpl)

         !Some thresholding which you should check later!
         !Here I get rid of all the regions which have intensity
         !below 1
         sbpitr => MeshIn%subpatch%begin()
         DO WHILE (ASSOCIATED(sbpitr))
            CALL sbpitr%get_field(labels,wpl,info)
            or_fail("Failed to get field i1_wp data.")

            trstat => tmp_region_stat%begin()
            DO WHILE (ASSOCIATED(trstat))
               value => trstat%getValue()

               !TODO
               !TOCHECK
               !check this criteria later
               IF (value(3).LT.smalld) THEN
                  wpl=MERGE(0,wpl,ABS(wpl).EQ.INT(value(1)))

                  CALL tmp_region_stat%remove(info)
                  or_fail("tmp_region_stat%remove")

                  DEALLOCATE(trstat,STAT=info)
                  or_fail_dealloc("Failed to deallocate trstat")
               ENDIF

               trstat => tmp_region_stat%next()
            ENDDO
            sbpitr => MeshIn%subpatch%next()
         ENDDO

         !-------------------------------------------------------------------------
         !  Return
         !-------------------------------------------------------------------------
      9999 CONTINUE
         CALL substop(caller,t0,info)
      END SUBROUTINE DTYPE(initfire)


      SUBROUTINE DTYPE(fire)(MeshIn,info)

        USE ppm_rc_module_energy, ONLY : e_data
        USE ppm_rc_module_util, ONLY : ppm_rc_label_exist
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
        TYPE(ppm_rc_list), POINTER :: vseed

#if   __DIME == __2D
        REAL(MK), CONTIGUOUS, DIMENSION(:,:),   POINTER :: wpi
#elif __DIME == __3D
        REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER :: wpi
#endif
        REAL(ppm_kind_double), DIMENSION(:),    POINTER :: value
        REAL(ppm_kind_double), DIMENSION(3)             :: val
        REAL(ppm_kind_double)                           :: t0,dummy

#if   __DIME == __2D
        INTEGER, CONTIGUOUS, DIMENSION(:,:),   POINTER     :: wpl
#elif __DIME == __3D
        INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER     :: wpl
#endif
        INTEGER, DIMENSION(:),     POINTER     :: seedn
        INTEGER, DIMENSION(:),     POINTER     :: Nm
        INTEGER, DIMENSION(:,:),   ALLOCATABLE :: labelled
        INTEGER, DIMENSION(__DIME)                  :: ld,ld_

        INTEGER                                :: i,ipatch,nsize
        INTEGER                                :: iseede,iseed
        INTEGER                                :: oldlabel,newlabel,vLabel
        INTEGER                                :: label_region
        INTEGER                                :: old_label_region,new_label_region

        CHARACTER(LEN=ppm_char) :: caller="fire"

        LOGICAL :: IsEnclosedByLabelFGConnectivity

        !-------------------------------------------------------------------------
        !  Externals
        !-------------------------------------------------------------------------

        !-------------------------------------------------------------------------
        !  Initialize
        !-------------------------------------------------------------------------
         CALL substart(caller,t0,info)

         nsize=64
         ALLOCATE(labelled(__DIME,nsize),STAT=info)
         or_fail_alloc("labelled")

         NULLIFY(wpi,wpl)

         sbpitr => MeshIn%subpatch%begin()
         ipatch=1
         i=MeshIn%subpatch%nb
         iseed=SUM(ppm_rc_seeds(1:i)%nb)
         iseede=iseed
         patch_loop: DO WHILE (ASSOCIATED(sbpitr))
            Nm => sbpitr%nnodes

            CALL sbpitr%get_field(image,wpi,info)
            or_fail("Failed to get field r_wp data.")

            CALL sbpitr%get_field(labels,wpl,info)
            or_fail("Failed to get field i_wp data.")

            seed => ppm_rc_seeds(ipatch)%last()
            seed_loop: DO WHILE (ASSOCIATED(seed))
               seedn => seed%first%getValue()

#if   __DIME == __2D
               oldlabel=ABS(wpl(seedn(1),seedn(2)))
#elif __DIME == __3D
               oldlabel=ABS(wpl(seedn(1),seedn(2),seedn(3)))
#endif
               IF (oldlabel.EQ.0) THEN
                  CALL ppm_rc_seeds(ipatch)%remove(info)
                  or_fail("ppm_rc_seeds(ipatch)%remove")

                  DEALLOCATE(seed,STAT=info)
                  or_fail_dealloc("Failed to deallocate seed")

                  seed => ppm_rc_seeds(ipatch)%at(ppm_rc_seeds(ipatch)%iter_id)
                  iseed=iseed-1
                  CYCLE seed_loop
               ENDIF
               IF (ppm_rc_label_exist(oldlabel,nlabels(iseed:iseede),iseede-iseed+1,.TRUE.)) THEN
                  CALL ppm_rc_seeds(ipatch)%remove(info)
                  or_fail("ppm_rc_seeds(ipatch)%remove")

                  DEALLOCATE(seed,STAT=info)
                  or_fail_dealloc("Failed to deallocate seed")

                  seed => ppm_rc_seeds(ipatch)%at(ppm_rc_seeds(ipatch)%iter_id)
                  iseed=iseed-1
                  CYCLE seed_loop
               ENDIF

               newlabel=nlabels(iseed)

               CALL ppm_rc_floodFill(labelled,wpl,Nm,seedn(1:__DIME),oldlabel,newlabel,0,info)
               or_fail("ppm_rc_floodFill")

               nsize=COUNT(labelled(1,:).GE.0)

               SELECT CASE (nsize)
               CASE (1)
                  !If this is a single point region, we need to remove it,
                  !it is non signifacant
                  !
                  !We turn the enclosed point to the enclosing region
                  !In case of being neighbor to
                  ld=labelled(:,1)

                  IF (ANY(ld.LT.1.OR.ld.GT.Nm)) THEN
                     !We are in the ghost layer
                  ELSE
                     ld_=ld+FG_ConnectivityType%NeighborsPoints(:,1)
#if   __DIME == __2D
                     vLabel=ABS(wpl(ld_(1),ld_(2)))
#elif __DIME == __3D
                     vLabel=ABS(wpl(ld_(1),ld_(2),ld_(3)))
#endif

                     newlabel=0

                     IF (vLabel.EQ.0) THEN
                        IsEnclosedByLabelFGConnectivity=.FALSE.
                     ELSE
                        IsEnclosedByLabelFGConnectivity=.TRUE.

                        IF (vLabel.EQ.FORBIDDEN) THEN
                           ld_=ld+FG_ConnectivityType%NeighborsPoints(:,2)
#if   __DIME == __2D
                           vLabel=ABS(wpl(ld_(1),ld_(2)))
#elif __DIME == __3D
                           vLabel=ABS(wpl(ld_(1),ld_(2),ld_(3)))
#endif
                        ENDIF

                        DO i=2,FG_ConnectivityType%NumberOfNeighbors
                           ld_=ld+FG_ConnectivityType%NeighborsPoints(:,i)
#if   __DIME == __2D
                           IF (ABS(wpl(ld_(1),ld_(2))).NE.vLabel) THEN
                              IF (ABS(wpl(ld_(1),ld_(2))).EQ.FORBIDDEN) CYCLE
#elif __DIME == __3D
                           IF (ABS(wpl(ld_(1),ld_(2),ld_(3))).NE.vLabel) THEN
                              IF (ABS(wpl(ld_(1),ld_(2),ld_(3))).EQ.FORBIDDEN) CYCLE
#endif
                              IsEnclosedByLabelFGConnectivity=.FALSE.
                              EXIT
                           ENDIF
                        ENDDO
                     ENDIF

                     IF (IsEnclosedByLabelFGConnectivity) THEN
                        newlabel=vLabel
#if   __DIME == __2D
                        wpl(ld(1),ld(2))=newlabel
#elif __DIME == __3D
                        wpl(ld(1),ld(2),ld(3))=newlabel
#endif

                        old_label_region=htable%search(oldlabel)
                        new_label_region=htable%search(newlabel)
                        IF (old_label_region.NE.htable_null.AND.new_label_region.NE.htable_null) THEN
                           CALL e_data%UpdateStatisticsWhenJump(wpi,ld,oldlabel,newlabel,info)
                           or_fail("e_data%UpdateStatisticsWhenJump")
                        ELSE IF (old_label_region.NE.htable_null.AND.new_label_region.EQ.htable_null) THEN
                           e_data%lCount(old_label_region)=e_data%lCount(old_label_region)-oned
                           i=old_label_region*2
#if   __DIME == __2D
                           dummy=REAL(wpi(ld(1),ld(2)),ppm_kind_double)
#elif __DIME == __3D
                           dummy=REAL(wpi(ld(1),ld(2),ld(3)),ppm_kind_double)
#endif
                           e_data%lSumslSumsq(i)  =e_data%lSumslSumsq(i)  -dummy
                           e_data%lSumslSumsq(i+1)=e_data%lSumslSumsq(i+1)-dummy*dummy

                           trstat => tmp_region_stat%begin()
                           DO WHILE (ASSOCIATED(trstat))
                              value => trstat%getValue()
                              IF (INT(value(1)).EQ.newlabel) THEN
                                 value(2)=value(2)+oned
                                 value(3)=value(3)+dummy
                                 value(4)=value(4)+dummy*dummy
                                 EXIT
                              ENDIF
                              trstat => tmp_region_stat%next()
                           ENDDO
                        ELSE IF (old_label_region.EQ.htable_null.AND.new_label_region.NE.htable_null) THEN
                           e_data%lCount(new_label_region)=e_data%lCount(new_label_region)+oned
                           i=new_label_region*2
#if   __DIME == __2D
                           dummy=REAL(wpi(ld(1),ld(2)),ppm_kind_double)
#elif __DIME == __3D
                           dummy=REAL(wpi(ld(1),ld(2),ld(3)),ppm_kind_double)
#endif
                           e_data%lSumslSumsq(i)  =e_data%lSumslSumsq(i)  +dummy
                           e_data%lSumslSumsq(i+1)=e_data%lSumslSumsq(i+1)+dummy*dummy

                           trstat => tmp_region_stat%begin()
                           DO WHILE (ASSOCIATED(trstat))
                              value => trstat%getValue()
                              IF (INT(value(1)).EQ.oldlabel) THEN
                                 value(2)=value(2)-oned
                                 value(3)=value(3)-dummy
                                 value(4)=value(4)-dummy*dummy
                                 EXIT
                              ENDIF
                              trstat => tmp_region_stat%next()
                           ENDDO
                        ELSE
#if   __DIME == __2D
                           dummy=REAL(wpi(ld(1),ld(2)),ppm_kind_double)
#elif __DIME == __3D
                           dummy=REAL(wpi(ld(1),ld(2),ld(3)),ppm_kind_double)
#endif

                           trstat => tmp_region_stat%begin()
                           DO WHILE (ASSOCIATED(trstat))
                              value => trstat%getValue()
                              IF (INT(value(1)).EQ.oldlabel) THEN
                                 value(2)=value(2)-oned
                                 value(3)=value(3)-dummy
                                 value(4)=value(4)-dummy*dummy
                                 EXIT
                              ENDIF
                              trstat => tmp_region_stat%next()
                           ENDDO
                           trstat => tmp_region_stat%begin()
                           DO WHILE (ASSOCIATED(trstat))
                              value => trstat%getValue()
                              IF (INT(value(1)).EQ.newlabel) THEN
                                 value(2)=value(2)+oned
                                 value(3)=value(3)+dummy
                                 value(4)=value(4)+dummy*dummy
                                 EXIT
                              ENDIF
                              trstat => tmp_region_stat%next()
                           ENDDO
                        ENDIF !old_label_region.NE.htable_null.AND.new_label_region.NE.htable_null

                        ALLOCATE(vseed,STAT=info)
                        or_fail_alloc("vseed")
                        CALL vseed%add(ld)
                        !Adding the coordinate for removing the particle from
                        !contourPoints list
                        CALL m_Seeds(ipatch)%push(vseed,info)
                        or_fail("could not add new seed to the collection")
                     ELSE
#if   __DIME == __2D
                        wpl(ld(1),ld(2))=0
#elif __DIME == __3D
                        wpl(ld(1),ld(2),ld(3))=0
#endif

                        old_label_region=htable%search(oldlabel)
                        IF (old_label_region.NE.htable_null) THEN
                           CALL e_data%UpdateStatisticsWhenJump(wpi,ld,oldlabel,0,info)
                           or_fail("e_data%UpdateStatisticsWhenJump")
                        ELSE
                           e_data%lCount(0)=e_data%lCount(0)+oned
#if   __DIME == __2D
                           dummy=REAL(wpi(ld(1),ld(2)),ppm_kind_double)
#elif __DIME == __3D
                           dummy=REAL(wpi(ld(1),ld(2),ld(3)),ppm_kind_double)
#endif
                           e_data%lSumslSumsq(0)=e_data%lSumslSumsq(0)+dummy
                           e_data%lSumslSumsq(1)=e_data%lSumslSumsq(1)+dummy*dummy

                           trstat => tmp_region_stat%begin()
                           DO WHILE (ASSOCIATED(trstat))
                              value => trstat%getValue()
                              IF (INT(value(1)).EQ.oldlabel) THEN
                                 value(2)=value(2)-oned
                                 value(3)=value(3)-dummy
                                 value(4)=value(4)-dummy*dummy
                                 EXIT
                              ENDIF
                              trstat => tmp_region_stat%next()
                           ENDDO
                        ENDIF
                     ENDIF !IsEnclosedByLabelFGConnectivity
                  ENDIF !ANY(ld.LT.1.OR.ld.GT.Nm)

               CASE DEFAULT
                  ALLOCATE(trstat,STAT=info)
                  or_fail_alloc("trstat")
                  !at this point there is no stat, so any point that is start
                  !firing, will be part of the new region and we do not
                  !need to kill the old stat

                  dummy=REAL(newlabel,ppm_kind_double)

                  CALL trstat%add(dummy,zerod,zerod,zerod)
                  CALL tmp_region_stat%push(trstat,info)
                  or_fail("tmp_region_stat%push")

                  val=zerod

                  trstat => tmp_region_stat%last()

                  SELECT CASE (AllowFusion)
                  CASE (.TRUE.)
                     DO i=1,nsize
                        ld=labelled(:,i)
                        IF (ANY(ld.LT.1.OR.ld.GT.Nm)) CYCLE
#if   __DIME == __2D
                        IF (wpl(ld(1),ld(2)).LT.0) THEN
#elif __DIME == __3D
                        IF (wpl(ld(1),ld(2),ld(3)).LT.0) THEN
#endif
                           CALL seed%add(ld)
                           !IF Fusion is allowed, later we wanna check the
                           !border particles to see whether they are enclosed
                           !by other label particles or not, so we will keep
                           !their indexes
                        ENDIF

                        val(1)=val(1)+oned
#if   __DIME == __2D
                        dummy=REAL(wpi(ld(1),ld(2)),ppm_kind_double)
#elif __DIME == __3D
                        dummy=REAL(wpi(ld(1),ld(2),ld(3)),ppm_kind_double)
#endif
                        val(2)=val(2)+dummy
                        val(3)=val(3)+dummy*dummy
                     ENDDO !i=1,nsize

                     CALL trstat%add(val)

                  CASE DEFAULT
                     DO i=1,nsize
                        ld=labelled(:,i)
                        IF (ANY(ld.LT.1.OR.ld.GT.Nm)) CYCLE

                        val(1)=val(1)+oned
#if   __DIME == __2D
                        dummy=REAL(wpi(ld(1),ld(2)),ppm_kind_double)
#elif __DIME == __3D
                        dummy=REAL(wpi(ld(1),ld(2),ld(3)),ppm_kind_double)
#endif
                        val(2)=val(2)+dummy
                        val(3)=val(3)+dummy*dummy
                     ENDDO !i=1,nsize

                     CALL trstat%add(val)

                  END SELECT !AllowFusion

                  label_region=htable%search(oldlabel)

                  check_true(<#label_region.NE.htable_null#>, &
                  & "Fail!!!, There should be an oldlabel available.")

                  value => trstat%getValue()

                  e_data%lCount(label_region)=e_data%lCount(label_region)-value(2)
                  i=label_region*2
                  e_data%lSumslSumsq(i)  =e_data%lSumslSumsq(i)  -value(3)
                  e_data%lSumslSumsq(i+1)=e_data%lSumslSumsq(i+1)-value(4)
                  !update the old region statistics

                  NULLIFY(trstat)

               END SELECT !nsize

               IF (2*nsize.LE.SIZE(labelled,DIM=2)) THEN
                  DEALLOCATE(labelled,STAT=info)
                  or_fail_dealloc("labelled")
                  ALLOCATE(labelled(__DIME,nsize),STAT=info)
                  or_fail_alloc("labelled")
               ENDIF

               seed => ppm_rc_seeds(ipatch)%prev()
               iseed=iseed-1
            ENDDO seed_loop

            sbpitr => MeshIn%subpatch%next()
            ipatch=ipatch+1
         ENDDO patch_loop

         IF (ALLOCATED(labelled)) THEN
            DEALLOCATE(labelled,STAT=info)
            or_fail_dealloc("labelled")
         ENDIF

         !-------------------------------------------------------------------------
         !  Return
         !-------------------------------------------------------------------------
      9999 CONTINUE
         CALL substop(caller,t0,info)
      END SUBROUTINE DTYPE(fire)

