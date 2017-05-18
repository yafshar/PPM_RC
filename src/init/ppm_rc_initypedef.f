
        TYPE,ABSTRACT :: ImageSource
          INTEGER                            :: m_ForegroundValue
          INTEGER                            :: m_BackgroundValue
          INTEGER, DIMENSION(:), ALLOCATABLE :: m_Size
          !!! size of the output image for initialization task
          INTEGER, DIMENSION(:), ALLOCATABLE :: m_Origin
          !!! origin
          !!! Parameters for the FG & BG
          INTEGER, DIMENSION(:), ALLOCATABLE :: m_Spacing
          !!! spacing between two inital objects origin
          !!! the spacing is measured from the start point
          !!! to start point ~ m_Spacing = m_Size+d
        CONTAINS
          PROCEDURE(ImageSource_create), DEFERRED :: create
          PROCEDURE :: destroy => ImageSource_destroy
        END TYPE ImageSource

        TYPE,EXTENDS(ImageSource) :: FFileImageSource
        !m_Size is the rectangle length
        CONTAINS
          PROCEDURE :: create => FFileImageSource_create
        END TYPE FFileImageSource

        TYPE,EXTENDS(ImageSource) :: RectangularImageSource
        !m_Size is the rectangle length
        CONTAINS
          PROCEDURE :: create => RectangularImageSource_create
        END TYPE RectangularImageSource

        TYPE,EXTENDS(ImageSource) :: SphereImageSource
        !m_Size is the radius
        CONTAINS
          PROCEDURE :: create => SphereImageSource_create
        END TYPE SphereImageSource

        TYPE,EXTENDS(ImageSource) :: OtsuImageSource
        CONTAINS
          PROCEDURE :: create => OtsuImageSource_create
        END TYPE OtsuImageSource

        ! Base class for initialization functions for
        ! the Region Competition optimizer and sampler
        TYPE :: RCInitClass
          CLASS(ImageSource), POINTER :: elem => NULL()

        CONTAINS
          PROCEDURE :: CreateElement  => RCInitClass_CreateStructuringElement
          PROCEDURE :: DestroyElement => RCInitClass_DestroyStructuringElement
          PROCEDURE :: GetOutput_2d   => RCInitClass_GetOutput_2d
          PROCEDURE :: GetOutput_3d   => RCInitClass_GetOutput_3d
        END TYPE RCInitClass

        !----------------------------------------------------------------------
        !  INTERFACES
        !----------------------------------------------------------------------
        INTERFACE
          ! Constructor
          SUBROUTINE ImageSource_create(this,info,ForegroundValue, &
          &          BackgroundValue,size_,origin_,spacing_)
            IMPORT :: ImageSource
            IMPLICIT NONE
            CLASS(ImageSource)                             :: this
            INTEGER,                         INTENT(  OUT) :: info
            INTEGER,               OPTIONAL, INTENT(IN   ) :: ForegroundValue
            INTEGER,               OPTIONAL, INTENT(IN   ) :: BackgroundValue
            INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN   ) :: size_
            INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN   ) :: origin_
            INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN   ) :: spacing_
          END SUBROUTINE ImageSource_create
        END INTERFACE
