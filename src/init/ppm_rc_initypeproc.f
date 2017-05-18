        ! destroy
        SUBROUTINE ImageSource_destroy(this,info)

          IMPLICIT NONE

          CLASS(ImageSource)     :: this

          INTEGER, INTENT(  OUT) :: info

          start_subroutine("ImageSource_destroy")

          IF (ALLOCATED(this%m_Size)) THEN
             DEALLOCATE(this%m_Size,STAT=info)
             or_fail_dealloc("this%m_Size")
          ENDIF
          IF (ALLOCATED(this%m_Origin)) THEN
             DEALLOCATE(this%m_Origin,STAT=info)
             or_fail_dealloc("this%m_Origin")
          ENDIF
          IF (ALLOCATED(this%m_Spacing)) THEN
             DEALLOCATE(this%m_Spacing,STAT=info)
             or_fail_dealloc("this%m_Spacing")
          ENDIF

          end_subroutine()

        END SUBROUTINE ImageSource_destroy

        ! Constructor
        SUBROUTINE FFileImageSource_create(this,info,     &
        &          ForegroundValue,BackgroundValue,size_, &
        &          origin_,spacing_)

          IMPLICIT NONE

          CLASS(FFileImageSource)                         :: this

          INTEGER,                         INTENT(  OUT) :: info
          INTEGER,               OPTIONAL, INTENT(IN   ) :: ForegroundValue
          INTEGER,               OPTIONAL, INTENT(IN   ) :: BackgroundValue
          INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN   ) :: size_
          INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN   ) :: origin_
          INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN   ) :: spacing_

          start_subroutine("FFileImageSource_create")

          IF (PRESENT(ForegroundValue)) THEN
             IF (ForegroundValue.LT.0) THEN
                this%m_ForegroundValue=1
             ELSE
                this%m_ForegroundValue=ForegroundValue
             ENDIF
          ELSE
             this%m_ForegroundValue=1
          ENDIF
          IF (PRESENT(BackgroundValue)) THEN
             IF (BackgroundValue.LT.0) THEN
                this%m_BackgroundValue=0
             ELSE
                this%m_BackgroundValue=BackgroundValue
             ENDIF
          ELSE
             this%m_BackgroundValue=0
          ENDIF

          end_subroutine()

        END SUBROUTINE FFileImageSource_create

        ! Constructor
        SUBROUTINE RectangularImageSource_create(this,info, &
        &          ForegroundValue,BackgroundValue,size_,   &
        &          origin_,spacing_)

          IMPLICIT NONE

          CLASS(RectangularImageSource)                  :: this

          INTEGER,                         INTENT(  OUT) :: info
          INTEGER,               OPTIONAL, INTENT(IN   ) :: ForegroundValue
          INTEGER,               OPTIONAL, INTENT(IN   ) :: BackgroundValue
          INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN   ) :: size_
          INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN   ) :: origin_
          INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN   ) :: spacing_

          start_subroutine("RectangularImageSource_create")

          IF (PRESENT(ForegroundValue)) THEN
             IF (ForegroundValue.LT.0) THEN
                this%m_ForegroundValue=1
             ELSE
                this%m_ForegroundValue=ForegroundValue
             ENDIF
          ELSE
             this%m_ForegroundValue=1
          ENDIF

          IF (PRESENT(BackgroundValue)) THEN
             IF (BackgroundValue.LT.0) THEN
                this%m_BackgroundValue=0
             ELSE
                this%m_BackgroundValue=BackgroundValue
             ENDIF
          ELSE
             this%m_BackgroundValue=0
          ENDIF

          ALLOCATE(this%m_Size(ppm_rc_dim),this%m_Origin(ppm_rc_dim), &
          &        this%m_Spacing(ppm_rc_dim),STAT=info)
          or_fail_alloc("this%m_Size, this%m_Spacing & this%m_Origin")

          IF (PRESENT(size_)) THEN
             IF (ANY(size_.LT.0)) THEN
                this%m_Size=32
             ELSE
                this%m_Size=size_
             ENDIF
          ELSE
             this%m_Size=32
          ENDIF
          ! Initial image is 32 wide pixels in each direction.

          check_true(<#ALL(this%m_Size.LT.(max_phys(1:ppm_rc_dim)-min_phys(1:ppm_rc_dim)))#>, &
          & "Initialization size is bigger than the domain size")

          IF (PRESENT(spacing_)) THEN
             IF (ANY(spacing_.LT.0)) THEN
                this%m_Spacing=SUM(this%m_Size)/ppm_rc_dim+10
             ELSE
                this%m_Spacing=spacing_
             ENDIF
          ELSE
             this%m_Spacing=SUM(this%m_Size)/ppm_rc_dim+10
          ENDIF

          check_true(<#ALL(this%m_Spacing.GT.this%m_Size)#>, &
          & "displacement between two origins is not enough")

          this%m_Origin=this%m_Size/2

          end_subroutine()

        END SUBROUTINE RectangularImageSource_create

        ! Constructor
        SUBROUTINE SphereImageSource_create(this,info,    &
        &          ForegroundValue,BackgroundValue,size_, &
        &          origin_,spacing_)

          IMPLICIT NONE

          CLASS(SphereImageSource)                       :: this

          INTEGER,                         INTENT(  OUT) :: info
          INTEGER,               OPTIONAL, INTENT(IN   ) :: ForegroundValue
          INTEGER,               OPTIONAL, INTENT(IN   ) :: BackgroundValue
          INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN   ) :: size_
          INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN   ) :: origin_
          INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN   ) :: spacing_

          start_subroutine("SphereImageSource_create")

          IF (PRESENT(ForegroundValue)) THEN
             IF (ForegroundValue.LT.0) THEN
                this%m_ForegroundValue=1
             ELSE
                this%m_ForegroundValue=ForegroundValue
             ENDIF
          ELSE
             this%m_ForegroundValue=1
          ENDIF
          IF (PRESENT(BackgroundValue)) THEN
             IF (BackgroundValue.LT.0) THEN
                this%m_BackgroundValue=0
             ELSE
                this%m_BackgroundValue=BackgroundValue
             ENDIF
          ELSE
             this%m_BackgroundValue=0
          ENDIF

          ALLOCATE(this%m_Size(ppm_rc_dim),this%m_Origin(ppm_rc_dim), &
          &        this%m_Spacing(ppm_rc_dim),STAT=info)
          or_fail_alloc("this%m_Size, this%m_Spacing & this%m_Origin")

          IF (PRESENT(size_)) THEN
             IF (ANY(size_.LT.0)) THEN
                this%m_Size=16
             ELSE
                this%m_Size=size_
             ENDIF
          ELSE
             this%m_Size=16
          ENDIF
          ! Initial image has a 16 wide pixels radius.

          check_true(<#ALL(2*this%m_Size(1:ppm_rc_dim).LT.(max_phys(1:ppm_rc_dim)-min_phys(1:ppm_rc_dim)))#>, &
          & "Initial size is bigger than the domain size")

          IF (PRESENT(spacing_)) THEN
             IF (ANY(spacing_.LT.0)) THEN
                this%m_Spacing=2*SUM(this%m_Size)/ppm_rc_dim+10
             ELSE
                this%m_Spacing=spacing_
             ENDIF
          ELSE
             this%m_Spacing=2*SUM(this%m_Size)/ppm_rc_dim+10
          ENDIF

          check_true(<#ALL(this%m_Spacing.GT.2*this%m_Size)#>, &
          & "displacement between two origins is not enough")

          this%m_Origin=0

          end_subroutine()

        END SUBROUTINE SphereImageSource_create


        ! Constructor
        SUBROUTINE OtsuImageSource_create(this,info,      &
        &          ForegroundValue,BackgroundValue,size_, &
        &          origin_,spacing_)

          IMPLICIT NONE

          CLASS(OtsuImageSource)                         :: this

          INTEGER,                         INTENT(  OUT) :: info
          INTEGER,               OPTIONAL, INTENT(IN   ) :: ForegroundValue
          INTEGER,               OPTIONAL, INTENT(IN   ) :: BackgroundValue
          INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN   ) :: size_
          INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN   ) :: origin_
          INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN   ) :: spacing_

          start_subroutine("OtsuImageSource_create")

          IF (PRESENT(ForegroundValue)) THEN
             IF (ForegroundValue.LT.0) THEN
                this%m_ForegroundValue=1
             ELSE
                this%m_ForegroundValue=ForegroundValue
             ENDIF
          ELSE
             this%m_ForegroundValue=1
          ENDIF
          IF (PRESENT(BackgroundValue)) THEN
             IF (BackgroundValue.LT.0) THEN
                this%m_BackgroundValue=0
             ELSE
                this%m_BackgroundValue=BackgroundValue
             ENDIF
          ELSE
             this%m_BackgroundValue=0
          ENDIF

          end_subroutine()

        END SUBROUTINE OtsuImageSource_create



        SUBROUTINE RCInitClass_CreateStructuringElement(this, &
        &          dim_,info,ForegroundValue,BackgroundValue, &
        &          size_,origin_,spacing_)

          IMPLICIT NONE

          CLASS(RCInitClass)                             :: this

          INTEGER,                         INTENT(IN   ) :: dim_
          INTEGER,                         INTENT(  OUT) :: info
          INTEGER,               OPTIONAL, INTENT(IN   ) :: ForegroundValue
          INTEGER,               OPTIONAL, INTENT(IN   ) :: BackgroundValue
          INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN   ) :: size_
          INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN   ) :: origin_
          INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN   ) :: spacing_

          CLASS(ImageSource), POINTER :: elem

          INTEGER                  :: ForegroundValue_
          INTEGER                  :: BackgroundValue_
          INTEGER, DIMENSION(dim_) :: size__
          INTEGER, DIMENSION(dim_) :: origin__
          INTEGER, DIMENSION(dim_) :: spacing__

          start_subroutine("CreateStructuringElement")

          SELECT CASE (vInitKind)
          CASE (e_fromfile)
             check_true(<#LGT(TRIM(initimage),"")#>, &
             & "The initimage filename is not given in the Ctrl file!")

             ALLOCATE(FFileImageSource::this%elem,STAT=info)
          CASE (e_rect)
             ALLOCATE(RectangularImageSource::this%elem,STAT=info)
          CASE (e_sphere)
             ALLOCATE(SphereImageSource::this%elem,STAT=info)
          CASE (e_otsu)
             ALLOCATE(OtsuImageSource::this%elem,STAT=info)
          CASE (e_localmax)
             ALLOCATE(SphereImageSource::this%elem,STAT=info)
          CASE DEFAULT
             fail("This init mode has not been implemented yet!", &
             & ppm_error=ppm_error_fatal)
          END SELECT
          or_fail_alloc("elem")

          elem => this%elem

          ForegroundValue_=MERGE(ForegroundValue,-1,PRESENT(ForegroundValue))
          BackgroundValue_=MERGE(BackgroundValue,-1,PRESENT(BackgroundValue))
          IF (PRESENT(size_)) THEN
             size__=size_(1:dim_)
          ELSE
             size__=-1
          ENDIF
          IF (PRESENT(origin_)) THEN
             origin__=origin_(1:dim_)
          ELSE
             origin__=-1
          ENDIF
          IF (PRESENT(spacing_)) THEN
             spacing__=spacing_(1:dim_)
          ELSE
             spacing__=-1
          ENDIF

          CALL elem%create(info,ForegroundValue_,BackgroundValue_, &
          &    size__,origin__,spacing__)
          or_fail("elem%create")

          !DEALLOCATE read values from Ctrl file
          DEALLOCATE(init_rd,init_sp,STAT=info)
          or_fail_dealloc("init_rd & init_sp")

          end_subroutine()

        END SUBROUTINE RCInitClass_CreateStructuringElement


       SUBROUTINE RCInitClass_DestroyStructuringElement(this,info)

          IMPLICIT NONE

          CLASS(RCInitClass)     :: this

          INTEGER, INTENT(  OUT) :: info

          start_subroutine("DestroyStructuringElement")

          ASSOCIATE(elem => this%elem)
             CALL elem%destroy(info)
             or_fail("elem%destroy")
          END ASSOCIATE

          DEALLOCATE(this%elem,STAT=info)
          or_fail_dealloc("this%elem")

          NULLIFY(this%elem)

          end_subroutine()

        END SUBROUTINE RCInitClass_DestroyStructuringElement


