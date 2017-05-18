        SUBROUTINE DTYPE(RCInitClass_GetOutput)(this,FieldIn,MeshIn,FieldOut,info)

          USE ppm_module_interfaces, ONLY : ppm_t_equi_mesh_,ppm_t_field_
          USE ppm_module_mesh_typedef, ONLY : ppm_t_equi_mesh
          USE ppm_module_field_typedef, ONLY : ppm_t_field

#ifdef __Linux
          USE ppm_rc_module_util, ONLY : ppm_rc_mem_usage
#endif
          USE ppm_rc_module_fire, ONLY : ppm_rc_floodFill
          USE ppm_rc_module_linkedlist, ONLY : ppm_rc_list,ppm_rc_link
          USE ppm_rc_module_write, ONLY : DTYPE(ppm_rc_write_image), &
          &   DTYPE(ppm_rc_write_image_label)
          USE ppm_rc_module_filter, ONLY : OtsuThresholdImageFilter, &
          &   DTYPE(ppm_rc_GaussianImageFilter),DTYPE(ppm_rc_SobelImageFilter)
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(RCInitClass)                     :: this

          CLASS(ppm_t_field_),     POINTER       :: FieldIn
          !!! the source filed

          CLASS(ppm_t_equi_mesh_), POINTER       :: MeshIn

          CLASS(ppm_t_field_),     POINTER       :: FieldOut
          !!! the source filed

          INTEGER,                 INTENT(  OUT) :: info

          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          TYPE(OtsuThresholdImageFilter) :: vOtsuThsFilter

          CLASS(ImageSource),     POINTER :: elem

          CLASS(ppm_t_subpatch_), POINTER :: sbpitr

          CLASS(ppm_t_field_), POINTER :: image_gf
          !blured pixel intensities of the image

          TYPE(ppm_rc_list), POINTER ::  seed
          TYPE(ppm_rc_list), POINTER :: rseed
          TYPE(ppm_rc_link), POINTER :: seedlnk

#if   __DIME == __2D
          REAL(MK),              CONTIGUOUS, DIMENSION(:,:),   POINTER :: wpr
          REAL(ppm_kind_single), CONTIGUOUS, DIMENSION(:,:),   POINTER :: wprs
#elif __DIME == __3D
          REAL(MK),              CONTIGUOUS, DIMENSION(:,:,:), POINTER :: wpr
          REAL(ppm_kind_single), CONTIGUOUS, DIMENSION(:,:,:), POINTER :: wprs
          REAL(MK)                                                     :: sz
#endif
          REAL(MK)                                                     :: s,sm,sp,sx,sy
          REAL(MK),                          DIMENSION(__DIME)         :: elem_size

#if   __DIME == __2D
          INTEGER, CONTIGUOUS, DIMENSION(:,:),   POINTER :: wpi
#elif __DIME == __3D
          INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER :: wpi
#endif
          INTEGER,             DIMENSION(:),     POINTER :: seedn
          INTEGER,             DIMENSION(:),     POINTER :: Nm
          INTEGER,             DIMENSION(:),     POINTER :: Nmm
          INTEGER,             DIMENSION(:),     POINTER :: hi_a
          INTEGER,             DIMENSION(:),     POINTER :: lo_a
          INTEGER,             DIMENSION(:),     POINTER :: istart
          INTEGER,             DIMENSION(:),     POINTER :: iend
          INTEGER,             DIMENSION(__DIME)         :: bx,nx
          INTEGER,             DIMENSION(__DIME)         :: ld,ll,lt
          INTEGER,             DIMENSION(__DIME)         :: ls,le
          INTEGER                                        :: i,j,ii,jj
#if   __DIME == __3D
          INTEGER                                        :: k,kk
#endif
          INTEGER                                        :: ipatch
#ifdef __Linux
          INTEGER                                        :: memory
#endif
          INTEGER                                        :: BGValue
          INTEGER                                        :: FGValue

          LOGICAL :: lbox

          start_subroutine("GetOutput")

          ipatch=MeshIn%subpatch%nb
          IF (ipatch.GT.0) THEN
             ALLOCATE(ppm_rc_seeds(ipatch),STAT=info)
             or_fail_alloc("ppm_rc_seeds")
          ELSE
             fail("There is no patch defined for the initial mesh.")
          ENDIF

          elem => this%elem
          BGValue=elem%m_BackgroundValue
          FGValue=elem%m_ForegroundValue

          SELECT CASE (vInitKind)
          CASE (e_fromfile) ! From a file
          CASE (e_rect) ! rect
          !2 - rect     initializes with one or several rectangular regions.
             ALLOCATE(rseed,STAT=info)
             or_fail_alloc("rseed")

             ll=elem%m_Size-1
#if   __DIME == __2D
             j=0
             DO i=0,ll(1)
                CALL rseed%add(i,j)
             ENDDO !i
             j=ll(2)
             DO i=0,ll(1)
                CALL rseed%add(i,j)
             ENDDO !i
             i=0
             DO j=1,ll(2)-1
                CALL rseed%add(i,j)
             ENDDO !j
             i=ll(1)
             DO j=1,ll(2)-1
                CALL rseed%add(i,j)
             ENDDO !j
#elif __DIME == __3D
             k=0
             DO j=0,ll(2)
                DO i=0,ll(1)
                   CALL rseed%add(i,j,k)
                ENDDO !i
             ENDDO ! j
             k=ll(3)
             DO j=0,ll(2)
                DO i=0,ll(1)
                   CALL rseed%add(i,j,k)
                ENDDO !i
             ENDDO ! j

             j=0
             DO k=1,ll(3)-1
                DO i=0,ll(1)
                   CALL rseed%add(i,j,k)
                ENDDO !i
             ENDDO ! k
             j=ll(2)
             DO k=1,ll(3)-1
                DO i=0,ll(1)
                   CALL rseed%add(i,j,k)
                ENDDO !i
             ENDDO ! k

             i=0
             DO k=1,ll(3)-1
                DO j=1,ll(2)-1
                   CALL rseed%add(i,j,k)
                ENDDO !j
             ENDDO ! k
             i=ll(1)
             DO k=1,ll(3)-1
                DO j=1,ll(2)-1
                   CALL rseed%add(i,j,k)
                ENDDO !j
             ENDDO ! k
#endif


             Nm => MeshIn%Nm
             !number of global mesh nodes

             nx(1:__DIME)=Nm(1:__DIME)/elem%m_Spacing(1:__DIME)
             !Max number of cells in each direction

             NULLIFY(wpi)
             ipatch=0
             sbpitr => MeshIn%subpatch%begin()
             DO WHILE (ASSOCIATED(sbpitr))
                Nmm    => sbpitr%nnodes
                istart => sbpitr%istart
                iend   => sbpitr%iend

                ls=istart(1:__DIME)-elem%m_Size(1:__DIME)
                le=iend(1:__DIME)  +elem%m_Size(1:__DIME)

                ipatch=ipatch+1

                CALL sbpitr%get_field(FieldOut,wpi,info)
                or_fail("Failed to get field wpl data.")

                wpi=BGValue
                !initial FieldOut value to the Background Value

#if   __DIME == __3D
                !distance from the border
                bx(3)=-elem%m_Size(3)
                loop_k: DO k=1,nx(3)
                   bx(3)=bx(3)+elem%m_Spacing(3)
                   IF (bx(3).LT.ls(3).OR.bx(3).GT.le(3)) CYCLE loop_k
#endif
                !distance from the border
                bx(2)=-elem%m_Size(2)
                loop_j: DO j=1,nx(2)
                   bx(2)=bx(2)+elem%m_Spacing(2)
                   IF (bx(2).LT.ls(2).OR.bx(2).GT.le(2)) CYCLE loop_j
                   !distance from the border
                   bx(1)=-elem%m_Size(1)
                   loop_i: DO i=1,nx(1)
                      bx(1)=bx(1)+elem%m_Spacing(1)
                      IF (bx(1).LT.ls(1).OR.bx(1).GT.le(1)) CYCLE loop_i

                      ld=bx-istart+1
                      !change to local index

                      lbox=.TRUE.

                      IF (ALL(ld.GE.1.AND.ld.LE.Nmm)) THEN
                         IF (ALL(ld.LT.Nmm)) THEN
                            ALLOCATE(seed,STAT=info)
                            or_fail_alloc("seed")
#if   __DIME == __2D
                            CALL seed%add(ld(1)+1,ld(2)+1)
#elif __DIME == __3D
                            CALL seed%add(ld(1)+1,ld(2)+1,ld(3)+1)
#endif
                            CALL ppm_rc_seeds(ipatch)%push(seed,info)
                            or_fail("could not add new seed to the collection")
                            lbox=.FALSE.
                         ENDIF
                      ELSE IF (ALL(ld.LE.1).AND.ALL(ld+elem%m_Size-1.GE.Nmm)) THEN
                         ALLOCATE(seed,STAT=info)
                         or_fail_alloc("seed")
                         ll=1
                         CALL seed%add(ll)
                         CALL ppm_rc_seeds(ipatch)%push(seed,info)
                         or_fail("could not add new seed to the collection")

                         seed => ppm_rc_seeds(ipatch)%last()
                         IF (ASSOCIATED(seed)) THEN
                            seedn => seed%first%getValue()
                            CALL ppm_rc_floodFill(wpi,Nmm,seedn(1:__DIME),BGValue,FGValue,1,info)
                         ENDIF
#if   __DIME == __2D
                         EXIT loop_j
#elif __DIME == __3D
                         EXIT loop_k
#endif
                      ENDIF

                      seedlnk => rseed%first
                      DO WHILE (ASSOCIATED(seedlnk))
                         seedn => seedlnk%getValue()

                         ll=ld+seedn
                         !Origin + the border points

                         IF (ALL(ll.GE.1.AND.ll.LE.Nmm)) THEN
                            IF (lbox) THEN
                               IF (ALL(ll.GT.1.AND.ll.LT.Nmm)) THEN
#if   __DIME == __2D
                                  IF (.NOT. &
                                  &  ((seedn(1).EQ.0.AND.seedn(2).EQ.elem%m_Size(2)-1).OR. &
                                  &   (seedn(2).EQ.0.AND.seedn(1).EQ.elem%m_Size(1)-1)))  THEN
#elif __DIME == __3D
                                  IF (.NOT. &
                                  &  ((seedn(1).EQ.0.AND.seedn(2).EQ.0.AND.seedn(3).EQ.elem%m_Size(3)-1).OR. &
                                  &   (seedn(1).EQ.0.AND.seedn(3).EQ.0.AND.seedn(2).EQ.elem%m_Size(2)-1).OR. &
                                  &   (seedn(2).EQ.0.AND.seedn(3).EQ.0.AND.seedn(1).EQ.elem%m_Size(1)-1)))  THEN
#endif
                                  ALLOCATE(seed,STAT=info)
                                  or_fail_alloc("seed")

                                  IF (ALL(ll.GT.ld)) THEN
                                     IF (ALL(ll-1.GT.ld)) THEN
#if   __DIME == __2D
                                        CALL seed%add(ll(1)-1,ll(2)-1)
                                     ELSE IF (ll(1)-1.EQ.ld(1)) THEN
                                        CALL seed%add(ll(1),ll(2)-1)
                                     ELSE
                                        CALL seed%add(ll(1)-1,ll(2))
                                     ENDIF
#elif __DIME == __3D
                                        CALL seed%add(ll(1)-1,ll(2)-1,ll(3)-1)
                                     ELSE IF (ll(1)-1.EQ.ld(1)) THEN
                                        IF (ll(3)-1.EQ.ld(3)) THEN
                                           CALL seed%add(ll(1),ll(2)-1,ll(3))
                                        ELSE
                                           CALL seed%add(ll(1),ll(2)-1,ll(3)-1)
                                        ENDIF
                                     ELSE IF (ll(2)-1.EQ.ld(2)) THEN
                                        IF (ll(3)-1.EQ.ld(3)) THEN
                                           CALL seed%add(ll(1)-1,ll(2),ll(3))
                                        ELSE
                                           CALL seed%add(ll(1)-1,ll(2),ll(3)-1)
                                        ENDIF
                                     ELSE IF (ll(3)-1.EQ.ld(3)) THEN
                                        CALL seed%add(ll(1)-1,ll(2)-1,ll(3))
                                     ENDIF
#endif
                                  ELSE IF (ll(1).GT.ld(1)) THEN
#if   __DIME == __2D
                                     CALL seed%add(ll(1),ll(2)+1)
                                  ELSE IF (ll(2).GT.ld(2)) THEN
                                     CALL seed%add(ll(1)+1,ll(2))
#elif __DIME == __3D
                                     IF (ll(2).GT.ld(2)) THEN
                                        CALL seed%add(ll(1),ll(2),ll(3)+1)
                                     ELSE IF (ll(3).GT.ld(3)) THEN
                                        CALL seed%add(ll(1),ll(2)+1,ll(3))
                                     ELSE
                                        CALL seed%add(ll(1),ll(2)+1,ll(3)+1)
                                     ENDIF
                                  ELSE IF (ll(2).GT.ld(2)) THEN
                                  !When you are here it means that ll(1).eq.ld(1)
                                     IF (ll(3).GT.ld(3)) THEN
                                        CALL seed%add(ll(1)+1,ll(2),ll(3))
                                     ELSE
                                        CALL seed%add(ll(1)+1,ll(2),ll(3)+1)
                                     ENDIF
                                  ELSE IF (ll(3).GT.ld(3)) THEN
                                     CALL seed%add(ll(1)+1,ll(2)+1,ll(3))
#endif
                                  ENDIF

                                  CALL ppm_rc_seeds(ipatch)%push(seed,info)
                                  or_fail("could not add new seed to the collection")

                                  lbox=.FALSE.
                                  ENDIF
                               ENDIF !ALL(ll.GT.1.AND.ll.LT.Nmm)
                            ENDIF !lbox

#if   __DIME == __2D
                            wpi(ll(1),ll(2))=FGValue
#elif __DIME == __3D
                            wpi(ll(1),ll(2),ll(3))=FGValue
#endif
                         ENDIF !ALL(ll.GE.1.AND.ll.LE.Nmm)

                         seedlnk => seedlnk%nextLink()
                      ENDDO !ASSOCIATED(seedlnk)

                      !If we have a start point
                      IF (.NOT.lbox) THEN
                         seed => ppm_rc_seeds(ipatch)%last()
                         IF (ASSOCIATED(seed)) THEN
                            seedn => seed%first%getValue()
                            CALL ppm_rc_floodFill(wpi,Nmm,seedn(1:__DIME),BGValue,FGValue,1,info)
                         ENDIF
                      ENDIF
                   ENDDO loop_i
                ENDDO loop_j
#if   __DIME == __3D
                ENDDO loop_k
#endif

                sbpitr => MeshIn%subpatch%next()
             ENDDO !ASSOCIATED(sbpitr)

             NULLIFY(wpi)

          CASE (e_sphere) ! spheres
          !3 - sphere   initializes with one or several spherical regions.
             ALLOCATE(rseed,STAT=info)
             or_fail_alloc("rseed")

             SELECT CASE (MAXVAL(elem%m_Size))
             CASE (1)
                sm=oneminus
                sp=oneplus

             CASE (2:15)
                sm=0.5_MK
                sp=oneplus

             CASE DEFAULT
                sm=0.8_MK
                sp=oneplus

             END SELECT

             sx=REAL(elem%m_Size(1)*elem%m_Size(1))
             sy=REAL(elem%m_Size(2)*elem%m_Size(2))
#if   __DIME == __2D
             DO j=-elem%m_Size(2),elem%m_Size(2)
                DO i=-elem%m_Size(1),elem%m_Size(1)
                   s=REAL(i*i)/sx+REAL(j*j)/sy
                   IF (s.GE.sm.AND.s.LE.sp) THEN
                      CALL rseed%add(i,j)
                   ENDIF
                ENDDO !i
             ENDDO !j
#elif __DIME == __3D
             sz=REAL(elem%m_Size(3)*elem%m_Size(3))

             DO k=-elem%m_Size(3),elem%m_Size(3)
                DO j=-elem%m_Size(2),elem%m_Size(2)
                   DO i=-elem%m_Size(1),elem%m_Size(1)
                      s=REAL(i*i)/sx+REAL(j*j)/sy+REAL(k*k)/sz
                      IF (s.GE.sm.AND.s.LE.sp) THEN
                         CALL rseed%add(i,j,k)
                      ENDIF
                   ENDDO !i
                ENDDO !j
             ENDDO ! k
#endif


             Nm => MeshIn%Nm

    !         stdout('elem%m_Size')
    !         stdout(Nm)

             nx(1:__DIME)=Nm(1:__DIME)/elem%m_Spacing(1:__DIME)
             !Max number of cells in each direction

             NULLIFY(wpi)
             ipatch=0
             sbpitr => MeshIn%subpatch%begin()
             DO WHILE (ASSOCIATED(sbpitr))
                Nmm    => sbpitr%nnodes
                istart => sbpitr%istart
                iend   => sbpitr%iend

                ls=istart(1:__DIME)-elem%m_Size(1:__DIME)
                le=iend(1:__DIME)  +elem%m_Size(1:__DIME)

                ipatch=ipatch+1

                CALL sbpitr%get_field(FieldOut,wpi,info)
                or_fail("Failed to get field wpl data.")

                wpi=BGValue

#if   __DIME == __3D
                !distance from the border
                bx(3)=-elem%m_Size(3)
                loop_z: DO k=1,nx(3)
                   bx(3)=bx(3)+elem%m_Spacing(3)
                   IF (bx(3).LT.ls(3).OR.bx(3).GT.le(3)) CYCLE loop_z
#endif
                !distance from the border
                bx(2)=-elem%m_Size(2)
                loop_y: DO j=1,nx(2)
                   bx(2)=bx(2)+elem%m_Spacing(2)
                   IF (bx(2).LT.ls(2).OR.bx(2).GT.le(2)) CYCLE loop_y
                   !distance from the border
                   bx(1)=-elem%m_Size(1)
                   loop_x: DO i=1,nx(1)
                      bx(1)=bx(1)+elem%m_Spacing(1)
                      IF (bx(1).LT.ls(1).OR.bx(1).GT.le(1)) CYCLE loop_x

                      !IF ((i.LT.2.OR.i.GT.7).OR.(j.LT.2.OR.j.GT.7)) CYCLE

                      !change to local index
                      ld=bx-istart+1

                      lbox=.TRUE.

                      IF (ALL(ld.GE.1.AND.ld.LE.Nmm)) THEN
                         ALLOCATE(seed,STAT=info)
                         or_fail_alloc("seed")

                         CALL seed%add(ld)
                         CALL ppm_rc_seeds(ipatch)%push(seed,info)
                         or_fail("could not add new seed to the collection")
                         lbox=.FALSE.
                      ENDIF

                      IF (ALL(ld-elem%m_Size+1.LE.1).AND.ALL(ld+elem%m_Size-1.GE.Nmm)) THEN
                         IF (lbox) THEN
                            ALLOCATE(seed,STAT=info)
                            or_fail_alloc("seed")
                            ll=1
                            CALL seed%add(ll)
                            CALL ppm_rc_seeds(ipatch)%push(seed,info)
                            or_fail("could not add new seed to the collection")
                         ENDIF

                         seed => ppm_rc_seeds(ipatch)%last()
                         IF (ASSOCIATED(seed)) THEN
                            seedn => seed%first%getValue()
                            CALL ppm_rc_floodFill(wpi,Nmm,seedn(1:__DIME),BGValue,FGValue,1,info)
                         ENDIF
#if   __DIME == __2D
                         EXIT loop_y
#elif __DIME == __3D
                         EXIT loop_z
#endif
                      ENDIF

                      seedlnk => rseed%first
                      DO WHILE (ASSOCIATED(seedlnk))
                         seedn => seedlnk%getValue()
                         !Origin + the border points
                         ll=ld+seedn

                         IF (ALL(ll.GE.1.AND.ll.LE.Nmm)) THEN
                            IF (lbox) THEN
                               IF (ALL(ll.GT.1.AND.ll.LT.Nmm)) THEN
                                  ALLOCATE(seed,STAT=info)
                                  or_fail_alloc("seed")

                                  IF (ALL(ll.GT.ld)) THEN
                                     lt=ll-1
                                  ELSE IF (ALL(ll.LT.ld)) THEN
                                     lt=ll+1
                                  ELSE
                                     lt=ll
                                  ENDIF

                                  DO WHILE (lt(1).GT.ld(1).AND.lt(1).GT.1.AND.lt(1).LT.Nmm(1))
                                     lt(1)=lt(1)-1
                                  ENDDO
                                  DO WHILE (lt(2).GT.ld(2).AND.lt(2).GT.1.AND.lt(2).LT.Nmm(2))
                                     lt(2)=lt(2)-1
                                  ENDDO
#if   __DIME == __3D
                                  DO WHILE (lt(3).GT.ld(3).AND.lt(3).GT.1.AND.lt(3).LT.Nmm(3))
                                     lt(3)=lt(3)-1
                                  ENDDO
#endif
                                  DO WHILE (lt(1).LT.ld(1).AND.lt(1).GT.1.AND.lt(1).LT.Nmm(1))
                                     lt(1)=lt(1)+1
                                  ENDDO
                                  DO WHILE (lt(2).LT.ld(2).AND.lt(2).GT.1.AND.lt(2).LT.Nmm(2))
                                     lt(2)=lt(2)+1
                                  ENDDO
#if   __DIME == __3D
                                  DO WHILE (lt(3).LT.ld(3).AND.lt(3).GT.1.AND.lt(3).LT.Nmm(3))
                                     lt(3)=lt(3)+1
                                  ENDDO
#endif

                                  CALL seed%add(lt)
                                  CALL ppm_rc_seeds(ipatch)%push(seed,info)
                                  or_fail("could not add new seed to the collection")

                                  lbox=.FALSE.
                               ENDIF !ALL(ll.GT.1.AND.ll.LT.Nmm)
                            ENDIF !lbox

#if   __DIME == __2D
                            wpi(ll(1),ll(2))=FGValue
#elif __DIME == __3D
                            wpi(ll(1),ll(2),ll(3))=FGValue
#endif
                         ENDIF !ALL(ll.GE.1.AND.ll.LE.Nmm)

                         seedlnk => seedlnk%nextLink()
                      ENDDO !ASSOCIATED(seedlnk)

                      seed => ppm_rc_seeds(ipatch)%last()
                      IF (ASSOCIATED(seed)) THEN
                         seedn => seed%first%getValue()
                         CALL ppm_rc_floodFill(wpi,Nmm,seedn(1:__DIME),BGValue,FGValue,1,info)
                      ENDIF

                   ENDDO loop_x
                ENDDO loop_y
#if   __DIME == __3D
                ENDDO loop_z
#endif

                sbpitr => MeshIn%subpatch%next()
             ENDDO !ASSOCIATED(sbpitr)

             NULLIFY(wpi)

          CASE (e_otsu)
             !4 - otsu  performs an otsu thresholding on the read image label
             !which is FieldOut
             IF (LGT(TRIM(initimage),"")) THEN
                CALL vOtsuThsFilter%DTYPE(GenerateData)(FieldOut,       &
                &    MeshIn,histSize,info,lowerBound=INT(m_lowerBound), &
                &    upperBound=INT(m_upperBound),MaskValue=FGValue)
             ELSE
                fail("!!!This initialization mode needs a label input image for otsu thresholding!!!", &
                & ppm_error=ppm_error_fatal)
             ENDIF

          CASE (e_localmax)
          !6 - localmax initializes a region at local intensity maxima
             IF (ALL(elem%m_Size.EQ.0)) THEN
                !We are getting the localmaxima as an input, so we do not blur
                IF (ppm_nproc.GT.1) THEN
                   ghostsize=1
                   !-------------------------------------------------------------------------
                   !  Get field ghosts for image
                   !-------------------------------------------------------------------------
                   CALL MeshIn%map_ghost_get(info,ghostsize=ghostsize)
                   or_fail("MeshIn%map_ghost_get")

                   CALL FieldIn%map_ghost_push(MeshIn,info)
                   or_fail("FieldIn%map_ghost_push")
                   CALL MeshIn%map_isend(info)
                   or_fail("MeshIn%map_isend")
                   CALL FieldIn%map_ghost_pop(MeshIn,info)
                   or_fail("FieldIn%map_ghost_pop")
                ENDIF !ppm_nproc.GT.1


                NULLIFY(wprs)

                sbpitr => MeshIn%subpatch%begin()
                ipatch=1
                DO WHILE (ASSOCIATED(sbpitr))
                   Nm => sbpitr%nnodes
                   !-------------------------------------------------------------------------
                   !  Find all local maxima in the intensity field which
                   !  are regions with intensity > intensity_ths.
                   !  Put spheres at the position of the seeds (label image).
                   !-------------------------------------------------------------------------
                   SELECT CASE (FieldIn%data_type)
                   CASE (ppm_type_real_single)
                      CALL sbpitr%get_field(FieldIn,wprs,info)
                      or_fail("Failed to get field data tmp.")
                   CASE DEFAULT
                      fail("!!!This type is not implemented yet for this initialization!!!", &
                      & ppm_error=ppm_error_fatal)
                   END SELECT

                   ls=MERGE(2,1,sbpitr%istart.EQ.1)
                   le=MERGE(Nm-1,Nm,sbpitr%iend.EQ.Ngrid)

#if   __DIME == __2D
                   DO j=ls(2),le(2)
                      DO i=ls(1),le(1)
                         IF (wprs(i,j).GT.intensity_ths) THEN
                            IF (wprs(i,j).GE.wprs(i-1,j).AND. &
                            &   wprs(i,j).GE.wprs(i+1,j).AND. &
                            &   wprs(i,j).GE.wprs(i,j-1).AND. &
                            &   wprs(i,j).GE.wprs(i,j+1)) THEN
                               ALLOCATE(seed,STAT=info)
                               or_fail_alloc("seed")
                               CALL seed%add(i,j)
                               CALL ppm_rc_seeds(ipatch)%push(seed,info)
                               or_fail("could not add new seed to the collection")
                            ENDIF
                         ENDIF
                      ENDDO !i
                   ENDDO !j
#elif __DIME == __3D
                   DO k=ls(3),le(3)
                      DO j=ls(2),le(2)
                         DO i=ls(1),le(1)
                            IF (wprs(i,j,k).GT.intensity_ths) THEN
                               IF (wprs(i,j,k).GE.wprs(i-1,j,k).AND. &
                               &   wprs(i,j,k).GE.wprs(i+1,j,k).AND. &
                               &   wprs(i,j,k).GE.wprs(i,j-1,k).AND. &
                               &   wprs(i,j,k).GE.wprs(i,j+1,k).AND. &
                               &   wprs(i,j,k).GE.wprs(i,j,k-1).AND. &
                               &   wprs(i,j,k).GE.wprs(i,j,k+1)) THEN
                                  ALLOCATE(seed,STAT=info)
                                  or_fail_alloc("seed")
                                  CALL seed%add(i,j,k)
                                  CALL ppm_rc_seeds(ipatch)%push(seed,info)
                                  or_fail("could not add new seed to the collection")
                               ENDIF
                            ENDIF
                         ENDDO !i
                      ENDDO !j
                   ENDDO !k
#endif

                   sbpitr => MeshIn%subpatch%next()
                   ipatch = ipatch+1
                ENDDO

#ifdef __Linux
                IF (debug.GT.1) THEN
                   CALL ppm_rc_mem_usage(memory,info)
                   stdout("mem_usage with local creation=",memory)
                ENDIF
#endif
             ELSE
                !after blurring with a Gaussian filter

                !---------------------------------------------------------------------
                !  blur the image using Gaussian filter
                !---------------------------------------------------------------------
                SELECT CASE (lNormalize)
                CASE (.FALSE.)
                   fail("!!!This initialization mode will use the normalized image!!!", &
                   & ppm_error=ppm_error_fatal)

                END SELECT

                !---------------------------------------------------------------------
                !  blur the image using Gaussian filter
                !---------------------------------------------------------------------
                ! copy the image:
                ALLOCATE(ppm_t_field::image_gf,STAT=info)
                or_fail_alloc('Failed to allocate image_gf.')

                CALL image_gf%create(1,info,dtype=ppm_type_real_single,name="blured_pixel_intensity")
                or_fail("Create field failed!" )

                elem_size=REAL(elem%m_Size,MK)
                CALL DTYPE(ppm_rc_GaussianImageFilter)(elem_size,FieldIn,MeshIn,info, &
                &    image_gf,KernelFactor=1.95_MK)

!                 CALL DTYPE(ppm_rc_write_image)(image_gf,MeshIn,"Gaussian",info)

                IF (ppm_nproc.GT.1) THEN
                   ghostsize=1
                   !-------------------------------------------------------------------------
                   !  Get field ghosts for image
                   !-------------------------------------------------------------------------
                   CALL MeshIn%map_ghost_get(info,ghostsize=ghostsize)
                   or_fail("MeshIn%map_ghost_get")

                   CALL image_gf%map_ghost_push(MeshIn,info)
                   or_fail("image_gf%map_ghost_push")
                   CALL MeshIn%map_isend(info)
                   or_fail("MeshIn%map_isend")
                   CALL image_gf%map_ghost_pop(MeshIn,info)
                   or_fail("image_gf%map_ghost_pop")
                ENDIF !ppm_nproc.GT.1


                NULLIFY(wprs)

                sbpitr => MeshIn%subpatch%begin()
                ipatch=1
                DO WHILE (ASSOCIATED(sbpitr))
                   Nm => sbpitr%nnodes

                   !-------------------------------------------------------------------------
                   !  Find all local maxima in the intensity field which
                   !  are regions with intensity > intensity_ths.
                   !  Put spheres at the position of the seeds (label image).
                   !-------------------------------------------------------------------------
                   CALL sbpitr%get_field(image_gf,wprs,info)
                   or_fail("Failed to get field data tmp.")

#if   __DIME == __3D
                   !padding the left xy image face
                   IF (sbpitr%istart(3).EQ.1) THEN
                      IF (inighostsize(3).GT.0) THEN
                         DO j=0,Nm(2)+1
                            DO i=0,Nm(1)+1
                               wprs(i,j,0)=wprs(i,j,1)
                            ENDDO !i=0,Nm(1)+1
                         ENDDO !j=0,Nm(2)+1
                      ENDIF !(inighostsize(3).GT.0)
                   ENDIF !(sbpitr%istart(3).EQ.1)
                   !padding the right xy image face
                   IF (sbpitr%iend(3).EQ.Ngrid(3)) THEN
                      IF (inighostsize(3).GT.0) THEN
                         k=Nm(3)+1
                         DO j=0,Nm(2)+1
                            DO i=0,Nm(1)+1
                               wprs(i,j,k)=wprs(i,j,Nm(3))
                            ENDDO !i=0,Nm(1)+1
                         ENDDO !j=0,Nm(2)+1
                      ENDIF !(inighostsize(3).GT.0)
                   ENDIF !(sbpitr%iend(3).EQ.Ngrid(3))
#endif

                   !padding the left xz image face
                   IF (sbpitr%istart(2).EQ.1) THEN
                      IF (inighostsize(2).GT.0) THEN
#if   __DIME == __2D
                         DO i=0,Nm(1)+1
                            wprs(i,0)=wprs(i,1)
                         ENDDO !i=0,Nm(1)+1
#elif __DIME == __3D
                         DO k=0,Nm(3)+1
                            DO i=0,Nm(1)+1
                               wprs(i,0,k)=wprs(i,1,k)
                            ENDDO !i=0,Nm(1)+1
                         ENDDO !k=0,Nm(3)+1
#endif
                      ENDIF !(inighostsize(2).GT.0)
                   ENDIF !(sbpitr%istart(2).EQ.1)
                   !padding the right xz image face
                   IF (sbpitr%iend(2).EQ.Ngrid(2)) THEN
                      IF (inighostsize(2).GT.0) THEN
                         j=Nm(2)+1
#if   __DIME == __2D
                         DO i=0,Nm(1)+1
                            wprs(i,j)= wprs(i,Nm(2))
                         ENDDO !i=0,Nm(1)+1
#elif __DIME == __3D
                         DO k=0,Nm(3)+1
                            DO i=0,Nm(1)+1
                               wprs(i,j,k)=wprs(i,Nm(2),k)
                            ENDDO !i=0,Nm(1)+1
                         ENDDO !k=0,Nm(3)+1
#endif
                      ENDIF !(inighostsize(2).GT.0)
                   ENDIF !(sbpitr%iend(2).EQ.Ngrid(2))

                   !padding the left yz image face
                   IF (sbpitr%istart(1).EQ.1) THEN
                      IF (inighostsize(1).GT.0) THEN
#if   __DIME == __2D
                         DO j=0,Nm(2)+1
                            wprs(0,j)=wprs(1,j)
                         ENDDO !j=0,Nm(2)+1
#elif __DIME == __3D
                         DO k=0,Nm(3)+1
                            DO j=0,Nm(2)+1
                               wprs(0,j,k)=wprs(1,j,k)
                            ENDDO !j=0,Nm(2)+1
                         ENDDO !k=0,Nm(3)+1
#endif
                      ENDIF !(inighostsize(1).GT.0)
                   ENDIF !(sbpitr%istart(1).EQ.1)
                   !padding the right yz image face
                   IF (sbpitr%iend(1).EQ.Ngrid(1)) THEN
                      IF (inighostsize(1).GT.0) THEN
                         i=Nm(1)+1
#if   __DIME == __2D
                         DO j=0,Nm(2)+1
                            wprs(i,j)=wprs(Nm(1),j)
                         ENDDO !j=0,Nm(2)+1
#elif __DIME == __3D
                         DO k=0,Nm(3)+1
                            DO j=0,Nm(2)+1
                               wprs(i,j,k)=wprs(Nm(1),j,k)
                            ENDDO !j=0,Nm(2)+1
                         ENDDO !k=0,Nm(3)+1
#endif
                      ENDIF !(inighostsize(1).GT.0)
                   ENDIF !(sbpitr%iend(1).EQ.Ngrid(1))

#if   __DIME == __2D
                   DO j=1,Nm(2)
                      DO i=1,Nm(1)
                         IF (wprs(i,j).GT.intensity_ths) THEN
                            IF (wprs(i,j).GE.wprs(i-1,j).AND. &
                            &   wprs(i,j).GE.wprs(i+1,j).AND. &
                            &   wprs(i,j).GE.wprs(i,j-1).AND. &
                            &   wprs(i,j).GE.wprs(i,j+1)) THEN
                               ALLOCATE(seed,STAT=info)
                               or_fail_alloc("seed")
                               CALL seed%add(i,j)
                               CALL ppm_rc_seeds(ipatch)%push(seed,info)
                               or_fail("could not add new seed to the collection")
                            ENDIF
                         ENDIF
                      ENDDO !i
                   ENDDO !j
#elif __DIME == __3D
                   DO k=1,Nm(3)
                      DO j=1,Nm(2)
                         DO i=1,Nm(1)
                            IF (wprs(i,j,k).GT.intensity_ths) THEN
                               IF (wprs(i,j,k).GE.wprs(i-1,j,k).AND. &
                               &   wprs(i,j,k).GE.wprs(i+1,j,k).AND. &
                               &   wprs(i,j,k).GE.wprs(i,j-1,k).AND. &
                               &   wprs(i,j,k).GE.wprs(i,j+1,k).AND. &
                               &   wprs(i,j,k).GE.wprs(i,j,k-1).AND. &
                               &   wprs(i,j,k).GE.wprs(i,j,k+1)) THEN
                                  ALLOCATE(seed,STAT=info)
                                  or_fail_alloc("seed")
                                  CALL seed%add(i,j,k)
                                  CALL ppm_rc_seeds(ipatch)%push(seed,info)
                                  or_fail("could not add new seed to the collection")
                               ENDIF
                            ENDIF
                         ENDDO !i
                      ENDDO !j
                   ENDDO !k
#endif

                   sbpitr => MeshIn%subpatch%next()
                   ipatch = ipatch+1
                ENDDO

                NULLIFY(wprs)

#ifdef __Linux
                IF (debug.GT.1) THEN
                   CALL ppm_rc_mem_usage(memory,info)
                   stdout("mem_usage with local creation=",memory)
                ENDIF
#endif

                CALL image_gf%destroy(info)
                or_fail("Destroy image_gf field failed!")

                ! deallocate the pointer:
                dealloc_pointer("image_gf")

#ifdef __Linux
                IF (debug.GT.1) THEN
                   CALL ppm_rc_mem_usage(memory,info)
                   stdout("mem_usage at image_gf destroy=",memory)
                ENDIF
#endif
             ENDIF !(ALL(elem%m_Size.EQ.0))

             !TODO
             !TOCHECK
             !This is temporary solution
             ALLOCATE(rseed,STAT=info)
             or_fail_alloc("rseed")

             sm=0.5_MK
             sp=oneplus

             ii=MIN(elem%m_Size(1),5)
             jj=MIN(elem%m_Size(2),5)
             IF (ii.LE.1) ii=2
             IF (jj.LE.1) jj=2

             sx=REAL(ii*ii)
             sy=REAL(jj*jj)
#if   __DIME == __2D
             DO j=-jj,jj
                DO i=-ii,ii
                   s=REAL(i*i)/sx+REAL(j*j)/sy
                   IF (s.LE.sp) THEN               !s.GE.sm.AND.s.LE.sp) THEN
                      CALL rseed%add(i,j)
                   ENDIF
                ENDDO !i
             ENDDO !j
#elif __DIME == __3D
             kk=MIN(elem%m_Size(3),5)
             IF (kk.LE.1) kk=2

             sz=REAL(kk*kk)
             DO k=-kk,kk
                DO j=-jj,jj
                   DO i=-ii,ii
                      s=REAL(i*i)/sx+REAL(j*j)/sy+REAL(k*k)/sz
                      IF (s.LE.sp) THEN               !s.GE.sm.AND.s.LE.sp) THEN
                         CALL rseed%add(i,j,k)
                      ENDIF
                   ENDDO !i
                ENDDO !j
             ENDDO ! k
#endif

             NULLIFY(wpi,wpr)

             sbpitr => MeshIn%subpatch%begin()
             ipatch=1
             DO WHILE (ASSOCIATED(sbpitr))
                Nm => sbpitr%nnodes
                hi_a => sbpitr%hi_a
                lo_a => sbpitr%lo_a

!                 CALL sbpitr%get_field(image,wpr,info)
!                 or_fail("Failed to get field wpi data.")

                CALL sbpitr%get_field(FieldOut,wpi,info)
                or_fail("Failed to get field wpl data.")

                wpi=BGValue

                seed => ppm_rc_seeds(ipatch)%begin()
                DO WHILE (ASSOCIATED(seed))
                   seedn => seed%first%getValue()

! #if   __DIME == __2D
!                    IF (wpi(seedn(1),seedn(2)).EQ.BGValue) THEN
! #elif __DIME == __3D
!                    IF (wpi(seedn(1),seedn(2),seedn(3)).EQ.BGValue) THEN
! #endif
                      seedlnk => rseed%first
                      DO WHILE (ASSOCIATED(seedlnk))
                         Nmm => seedlnk%getValue()
                         !Origin + the border points
                         ll=seedn+Nmm

                         IF (ALL(ll.GE.lo_a.AND.ll.LE.hi_a)) THEN
#if   __DIME == __2D
                            wpi(ll(1),ll(2))=FGValue
#elif __DIME == __3D
                            wpi(ll(1),ll(2),ll(3))=FGValue
#endif
                         ENDIF !ALL(ll.GE.1.AND.ll.LE.Nmm)

                         seedlnk => seedlnk%nextLink()
                      ENDDO !ASSOCIATED(seedlnk)

!                       CALL ppm_rc_floodFill(wpi,Nm,seedn(1:__DIME), &
!                       &    BGValue,FGValue,0,info)
!                       CALL ppm_rc_floodFill(wpi,wpr,Nm,   &
!                       &    seedn(1:__DIME),BGValue,FGValue, &
!                       &    Tolerance_,Tolerance_,0,info)
!                       or_fail("ppm_rc_floodFill")
!                    ELSE
!                       CALL ppm_rc_seeds(ipatch)%remove(info)
!                       or_fail("ppm_rc_seeds(ipatch)%remove")
!                    ENDIF
                   seed => ppm_rc_seeds(ipatch)%next()
                ENDDO

                sbpitr => MeshIn%subpatch%next()
                ipatch = ipatch+1
             ENDDO

             NULLIFY(wpi,wpr)

          END SELECT !vInitKind

          SELECT CASE (vInitKind)
          CASE (e_fromfile,e_otsu)
          ! 1- From File and 4 - otsu
          ! If we read the init label from a file, we have everything.
          !
          ! Also for an otsu thresholding we do not need to put the ghost back
          ! on the real components

             sbpitr => MeshIn%subpatch%begin()
             DO WHILE (ASSOCIATED(sbpitr))
                Nm   => sbpitr%nnodes
                hi_a => sbpitr%hi_a
                lo_a => sbpitr%lo_a

                CALL sbpitr%get_field(FieldOut,wpi,info)
                or_fail("Failed to get field i_tmp data.")

#if   __DIME == __2D
                IF (sbpitr%istart(1).EQ.1) THEN
                   DO j=lo_a(2),hi_a(2)
                      wpi(1,j)=FORBIDDEN
                   ENDDO
                ENDIF
                IF (sbpitr%istart(2).EQ.1) THEN
                   DO i=lo_a(1),hi_a(1)
                      wpi(i,1)=FORBIDDEN
                   ENDDO
                ENDIF
                IF (sbpitr%iend(1).EQ.Ngrid(1)) THEN
                   DO j=lo_a(2),hi_a(2)
                      wpi(Nm(1),j)=FORBIDDEN
                   ENDDO
                ENDIF
                IF (sbpitr%iend(2).EQ.Ngrid(2)) THEN
                   DO i=lo_a(1),hi_a(1)
                      wpi(i,Nm(2))=FORBIDDEN
                   ENDDO
                ENDIF
#elif __DIME == __3D
                IF (sbpitr%istart(1).EQ.1) THEN
                   DO k=lo_a(3),hi_a(3)
                      DO j=lo_a(2),hi_a(2)
                         wpi(1,j,k)=FORBIDDEN
                      ENDDO
                   ENDDO
                ENDIF
                IF (sbpitr%istart(2).EQ.1) THEN
                   DO k=lo_a(3),hi_a(3)
                      DO i=lo_a(1),hi_a(1)
                         wpi(i,1,k)=FORBIDDEN
                      ENDDO
                   ENDDO
                ENDIF
                IF (sbpitr%istart(3).EQ.1) THEN
                   DO j=lo_a(2),hi_a(2)
                      DO i=lo_a(1),hi_a(1)
                         wpi(i,j,1)=FORBIDDEN
                      ENDDO
                   ENDDO
                ENDIF
                IF (sbpitr%iend(1).EQ.Ngrid(1)) THEN
                   DO k=lo_a(3),hi_a(3)
                      DO j=lo_a(2),hi_a(2)
                         wpi(Nm(1),j,k)=FORBIDDEN
                      ENDDO
                   ENDDO
                ENDIF
                IF (sbpitr%iend(2).EQ.Ngrid(2)) THEN
                   DO k=lo_a(3),hi_a(3)
                      DO i=lo_a(1),hi_a(1)
                         wpi(i,Nm(2),k)=FORBIDDEN
                      ENDDO
                   ENDDO
                ENDIF
                IF (sbpitr%iend(3).EQ.Ngrid(3)) THEN
                   DO j=lo_a(2),hi_a(2)
                      DO i=lo_a(1),hi_a(1)
                         wpi(i,j,Nm(3))=FORBIDDEN
                      ENDDO
                   ENDDO
                ENDIF
#endif

                sbpitr => MeshIn%subpatch%next()
             ENDDO !WHILE (ASSOCIATED(sbpitr))

          CASE DEFAULT
             IF (ppm_nproc.GT.1) THEN
                !-------------------------------------------------------------------------
                !  Ghost Put for FieldOut field
                !-------------------------------------------------------------------------
                CALL MeshIn%map_ghost_put(info)
                or_fail("MeshIn%map_ghost_put")

                CALL FieldOut%map_ghost_push(MeshIn,info)
                or_fail("FieldOut%map_ghost_push")
                CALL MeshIn%map_isend(info,sendrecv=.TRUE.)
                or_fail("MeshIn%map_isend")
             ENDIF

             CALL rseed%destroy()

             dealloc_pointer("rseed")

             IF (ppm_nproc.GT.1) THEN
                CALL MeshIn%map_isend(info,sendrecv=.FALSE.)
                or_fail("MeshIn%map_isend")
                CALL FieldOut%map_ghost_pop(MeshIn,info)
                or_fail("FieldOut%map_ghost_pop")
             ENDIF

             sbpitr => MeshIn%subpatch%begin()
             DO WHILE (ASSOCIATED(sbpitr))
                Nm   => sbpitr%nnodes
                hi_a => sbpitr%hi_a
                lo_a => sbpitr%lo_a

                CALL sbpitr%get_field(FieldOut,wpi,info)
                or_fail("Failed to get field i_tmp data.")

#if   __DIME == __2D
                DO j=1,Nm(2)
                   DO i=1,Nm(1)
                      IF (wpi(i,j).GT.FGValue) THEN
                         wpi(i,j)=FGValue
                      ENDIF
                   ENDDO !i=1,Nm(1)
                ENDDO !j=1,Nm(2)

                IF (sbpitr%istart(1).EQ.1) THEN
                   DO j=lo_a(2),hi_a(2)
                      wpi(1,j)=FORBIDDEN
                   ENDDO
                ENDIF
                IF (sbpitr%istart(2).EQ.1) THEN
                   DO i=lo_a(1),hi_a(1)
                      wpi(i,1)=FORBIDDEN
                   ENDDO
                ENDIF
                IF (sbpitr%iend(1).EQ.Ngrid(1)) THEN
                   DO j=lo_a(2),hi_a(2)
                      wpi(Nm(1),j)=FORBIDDEN
                   ENDDO
                ENDIF
                IF (sbpitr%iend(2).EQ.Ngrid(2)) THEN
                   DO i=lo_a(1),hi_a(1)
                      wpi(i,Nm(2))=FORBIDDEN
                   ENDDO
                ENDIF
#elif __DIME == __3D
                DO k=1,Nm(3)
                   DO j=1,Nm(2)
                      DO i=1,Nm(1)
                         IF (wpi(i,j,k).GT.FGValue) THEN
                            wpi(i,j,k)=FGValue
                         ENDIF
                      ENDDO !i=1,Nm(1)
                   ENDDO !j=1,Nm(2)
                ENDDO !k=1,Nm(3)

                IF (sbpitr%istart(1).EQ.1) THEN
                   DO k=lo_a(3),hi_a(3)
                      DO j=lo_a(2),hi_a(2)
                         wpi(1,j,k)=FORBIDDEN
                      ENDDO
                   ENDDO
                ENDIF
                IF (sbpitr%istart(2).EQ.1) THEN
                   DO k=lo_a(3),hi_a(3)
                      DO i=lo_a(1),hi_a(1)
                         wpi(i,1,k)=FORBIDDEN
                      ENDDO
                   ENDDO
                ENDIF
                IF (sbpitr%istart(3).EQ.1) THEN
                   DO j=lo_a(2),hi_a(2)
                      DO i=lo_a(1),hi_a(1)
                         wpi(i,j,1)=FORBIDDEN
                      ENDDO
                   ENDDO
                ENDIF
                IF (sbpitr%iend(1).EQ.Ngrid(1)) THEN
                   DO k=lo_a(3),hi_a(3)
                      DO j=lo_a(2),hi_a(2)
                         wpi(Nm(1),j,k)=FORBIDDEN
                      ENDDO
                   ENDDO
                ENDIF
                IF (sbpitr%iend(2).EQ.Ngrid(2)) THEN
                   DO k=lo_a(3),hi_a(3)
                      DO i=lo_a(1),hi_a(1)
                         wpi(i,Nm(2),k)=FORBIDDEN
                      ENDDO
                   ENDDO
                ENDIF
                IF (sbpitr%iend(3).EQ.Ngrid(3)) THEN
                   DO j=lo_a(2),hi_a(2)
                      DO i=lo_a(1),hi_a(1)
                         wpi(i,j,Nm(3))=FORBIDDEN
                      ENDDO
                   ENDDO
                ENDIF
#endif

                sbpitr => MeshIn%subpatch%next()
             ENDDO

          END SELECT !vInitKind

          end_subroutine()

        END SUBROUTINE DTYPE(RCInitClass_GetOutput)

