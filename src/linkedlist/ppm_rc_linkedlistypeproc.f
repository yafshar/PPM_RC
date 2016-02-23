        FUNCTION nextLink(this)
          CLASS(ppm_rc_link)         :: this

          TYPE(ppm_rc_link), POINTER :: nextLink

          nextLink => this%next
        END FUNCTION nextLink

        FUNCTION nextLink_(this)
          CLASS(ppm_rc_link_)         :: this

          TYPE(ppm_rc_link_), POINTER :: nextLink_

          nextLink_ => this%next
        END FUNCTION nextLink_

        FUNCTION nextLink_2(this)
          CLASS(ppm_rc_link_2)         :: this
          TYPE(ppm_rc_link_2), POINTER :: nextLink_2
          nextLink_2 => this%next
        END FUNCTION nextLink_2

        SUBROUTINE setNextLink(this,next)
          CLASS(ppm_rc_link)         :: this

          TYPE(ppm_rc_link), POINTER :: next

          this%next => next
        END SUBROUTINE setNextLink

        SUBROUTINE setNextLink_(this,next)
          CLASS(ppm_rc_link_)         :: this

          TYPE(ppm_rc_link_), POINTER :: next

          this%next => next
        END SUBROUTINE setNextLink_

        SUBROUTINE setNextLink_2(this,next)
          CLASS(ppm_rc_link_2)         :: this
          TYPE(ppm_rc_link_2), POINTER :: next
          this%next => next
        END SUBROUTINE setNextLink_2

        FUNCTION getValue_(this)
          CLASS(ppm_rc_link)             :: this

          INTEGER, DIMENSION(:), POINTER :: getValue_

          getValue_ => this%value
        END FUNCTION getValue_

        FUNCTION getValueMCMCParticle_2d(this,value1,value2,MCMCParticle_)
          USE ppm_rc_module_util, ONLY : MCMCParticle,HashIndexFunctor_2d
          IMPLICIT NONE
          CLASS(ppm_rc_link)                   :: this

          INTEGER,               INTENT(  OUT) :: value1
          INTEGER,               INTENT(  OUT) :: value2
          TYPE(MCMCParticle),    INTENT(  OUT) :: MCMCParticle_

          INTEGER, DIMENSION(:), POINTER       :: getValueMCMCParticle_2d

          getValueMCMCParticle_2d => this%value
          CALL HashIndexFunctor_2d(getValueMCMCParticle_2d(1),value1,value2)
          MCMCParticle_%candlabel=getValueMCMCParticle_2d(2)
          MCMCParticle_%proposal=TRANSFER(getValueMCMCParticle_2d(3),MCMCParticle_%proposal)
        END FUNCTION getValueMCMCParticle_2d

        FUNCTION getValueMCMCParticle_3d(this,value1,value2,value3,MCMCParticle_)
          USE ppm_rc_module_util, ONLY : MCMCParticle,HashIndexFunctor_3d
          IMPLICIT NONE
          CLASS(ppm_rc_link)                   :: this

          INTEGER,               INTENT(  OUT) :: value1
          INTEGER,               INTENT(  OUT) :: value2
          INTEGER,               INTENT(  OUT) :: value3
          TYPE(MCMCParticle),    INTENT(  OUT) :: MCMCParticle_

          INTEGER, DIMENSION(:), POINTER       :: getValueMCMCParticle_3d

          getValueMCMCParticle_3d => this%value
          CALL HashIndexFunctor_3d(getValueMCMCParticle_3d(1),value1,value2,value3)
          MCMCParticle_%candlabel=getValueMCMCParticle_3d(2)
          MCMCParticle_%proposal=TRANSFER(getValueMCMCParticle_3d(3),MCMCParticle_%proposal)
        END FUNCTION getValueMCMCParticle_3d

        FUNCTION getValueMCMCHistoryParticle_2d(this,value1,value2,MCMCHistoryParticle_)
          USE ppm_rc_module_util, ONLY : MCMCHistoryParticle,HashIndexFunctor_2d
          IMPLICIT NONE
          CLASS(ppm_rc_link)                       :: this

          INTEGER,                   INTENT(  OUT) :: value1
          INTEGER,                   INTENT(  OUT) :: value2
          TYPE(MCMCHistoryParticle), INTENT(  OUT) :: MCMCHistoryParticle_

          INTEGER, DIMENSION(:),     POINTER       :: getValueMCMCHistoryParticle_2d

          getValueMCMCHistoryParticle_2d => this%value
          CALL HashIndexFunctor_2d(getValueMCMCHistoryParticle_2d(1),value1,value2)
          MCMCHistoryParticle_%candlabel=getValueMCMCHistoryParticle_2d(2)
          MCMCHistoryParticle_%orglabel=getValueMCMCHistoryParticle_2d(3)
          MCMCHistoryParticle_%proposal=TRANSFER(getValueMCMCHistoryParticle_2d(4),MCMCHistoryParticle_%proposal)
          MCMCHistoryParticle_%wasadded=TRANSFER(getValueMCMCHistoryParticle_2d(5),.TRUE.)
        END FUNCTION getValueMCMCHistoryParticle_2d

        FUNCTION getValueMCMCHistoryParticle_3d(this,value1,value2,value3,MCMCHistoryParticle_)
          USE ppm_rc_module_util, ONLY : MCMCHistoryParticle,HashIndexFunctor_3d
          IMPLICIT NONE
          CLASS(ppm_rc_link)                       :: this

          INTEGER,                   INTENT(  OUT) :: value1
          INTEGER,                   INTENT(  OUT) :: value2
          INTEGER,                   INTENT(  OUT) :: value3
          TYPE(MCMCHistoryParticle), INTENT(  OUT) :: MCMCHistoryParticle_

          INTEGER, DIMENSION(:),     POINTER       :: getValueMCMCHistoryParticle_3d

          getValueMCMCHistoryParticle_3d => this%value
          CALL HashIndexFunctor_3d(getValueMCMCHistoryParticle_3d(1),value1,value2,value3)
          MCMCHistoryParticle_%candlabel=getValueMCMCHistoryParticle_3d(2)
          MCMCHistoryParticle_%orglabel=getValueMCMCHistoryParticle_3d(3)
          MCMCHistoryParticle_%proposal=TRANSFER(getValueMCMCHistoryParticle_3d(4),MCMCHistoryParticle_%proposal)
          MCMCHistoryParticle_%wasadded=TRANSFER(getValueMCMCHistoryParticle_3d(5),.TRUE.)
        END FUNCTION getValueMCMCHistoryParticle_3d

        FUNCTION getValue__(this)
          CLASS(ppm_rc_link_) :: this

          INTEGER             :: getValue__

          getValue__= this%value
        END FUNCTION getValue__

        FUNCTION getXYValue_(this)
          USE ppm_rc_module_util, ONLY : HashIndexFunctor_2d
          IMPLICIT NONE
          CLASS(ppm_rc_link_)   :: this
          INTEGER, DIMENSION(2) :: getXYValue_
          CALL HashIndexFunctor_2d(this%value,getXYValue_)
        END FUNCTION getXYValue_

        FUNCTION getXYZValue_(this)
          USE ppm_rc_module_util, ONLY : HashIndexFunctor_3d
          IMPLICIT NONE
          CLASS(ppm_rc_link_)   :: this
          INTEGER, DIMENSION(3) :: getXYZValue_
          CALL HashIndexFunctor_3d(this%value,getXYZValue_)
        END FUNCTION getXYZValue_

        FUNCTION getValue_2_2d(this,key1,key2) RESULT(value)
          USE ppm_rc_module_util, ONLY : HashIndexFunctor_2d
          IMPLICIT NONE
          CLASS(ppm_rc_link_2)   :: this
          INTEGER, INTENT(  OUT) :: key1
          INTEGER, INTENT(  OUT) :: key2
          INTEGER                :: value
          CALL HashIndexFunctor_2d(this%key,key1,key2)
          value=this%value
        END FUNCTION getValue_2_2d

        FUNCTION getValue_2_3d(this,key1,key2,key3) RESULT(value)
          USE ppm_rc_module_util, ONLY : HashIndexFunctor_3d
          IMPLICIT NONE
          CLASS(ppm_rc_link_2)   :: this
          INTEGER, INTENT(  OUT) :: key1
          INTEGER, INTENT(  OUT) :: key2
          INTEGER, INTENT(  OUT) :: key3
          INTEGER                :: value
          CALL HashIndexFunctor_3d(this%key,key1,key2,key3)
          value=this%value
        END FUNCTION getValue_2_3d

        FUNCTION getValue_2(this,key) RESULT(value)
          USE ppm_rc_module_util, ONLY : HashIndexFunctor_2d,HashIndexFunctor_3d
          IMPLICIT NONE
          CLASS(ppm_rc_link_2)                 :: this
          INTEGER, DIMENSION(:), INTENT(  OUT) :: key
          INTEGER                              :: value
          IF (SIZE(key).EQ.2) THEN
             CALL HashIndexFunctor_2d(this%key,key(1),key(2))
          ELSE
             CALL HashIndexFunctor_3d(this%key,key(1),key(2),key(3))
          ENDIF
          value=this%value
        END FUNCTION getValue_2


        FUNCTION constructor(value, next)
          INTEGER, DIMENSION(:), INTENT(IN   ) :: value

          TYPE(ppm_rc_link),     POINTER       :: next

          TYPE(ppm_rc_link),     POINTER       :: constructor

          INTEGER :: info,vsize

          info=0
          vsize=SIZE(value,DIM=1)

          ALLOCATE(constructor)
          ALLOCATE(constructor%value(vsize),SOURCE=value,STAT=info)
          IF (info.NE.0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc, "constructor", &
             &    'Failed to allocate value.',  __LINE__, info)
          ENDIF

          constructor%next => next

        END FUNCTION constructor

        FUNCTION constructor1(value1,next)
          INTEGER,           INTENT(IN   ) :: value1

          TYPE(ppm_rc_link), POINTER       :: next

          TYPE(ppm_rc_link), POINTER       :: constructor1

          ALLOCATE(constructor1)
          ALLOCATE(constructor1%value(1))
          constructor1%value(1)=value1
          constructor1%next => next
        END FUNCTION constructor1
        FUNCTION constructor2(value1,value2,next)
          INTEGER,           INTENT(IN   ) :: value1
          INTEGER,           INTENT(IN   ) :: value2

          TYPE(ppm_rc_link), POINTER       :: next

          TYPE(ppm_rc_link), POINTER       :: constructor2

          ALLOCATE(constructor2)
          ALLOCATE(constructor2%value(2))
          constructor2%value(1)=value1
          constructor2%value(2)=value2
          constructor2%next => next
        END FUNCTION constructor2
        FUNCTION constructor3(value1,value2,value3,next)
          INTEGER,           INTENT(IN   ) :: value1
          INTEGER,           INTENT(IN   ) :: value2
          INTEGER,           INTENT(IN   ) :: value3

          TYPE(ppm_rc_link), POINTER       :: next

          TYPE(ppm_rc_link), POINTER       :: constructor3

          ALLOCATE(constructor3)
          ALLOCATE(constructor3%value(3))
          constructor3%value(1)=value1
          constructor3%value(2)=value2
          constructor3%value(3)=value3
          constructor3%next => next
        END FUNCTION constructor3
        FUNCTION constructor4(value1,value2,value3,value4,next)
          INTEGER,           INTENT(IN   ) :: value1
          INTEGER,           INTENT(IN   ) :: value2
          INTEGER,           INTENT(IN   ) :: value3
          INTEGER,           INTENT(IN   ) :: value4

          TYPE(ppm_rc_link), POINTER       :: next

          TYPE(ppm_rc_link), POINTER       :: constructor4

          ALLOCATE(constructor4)
          ALLOCATE(constructor4%value(4))
          constructor4%value(1)=value1
          constructor4%value(2)=value2
          constructor4%value(3)=value3
          constructor4%value(4)=value4
          constructor4%next => next
        END FUNCTION constructor4
        FUNCTION constructor5(value1,value2,value3,value4,value5,next)
          INTEGER,           INTENT(IN   ) :: value1
          INTEGER,           INTENT(IN   ) :: value2
          INTEGER,           INTENT(IN   ) :: value3
          INTEGER,           INTENT(IN   ) :: value4
          INTEGER,           INTENT(IN   ) :: value5

          TYPE(ppm_rc_link), POINTER       :: next

          TYPE(ppm_rc_link), POINTER       :: constructor5

          ALLOCATE(constructor5)
          ALLOCATE(constructor5%value(5))
          constructor5%value(1)=value1
          constructor5%value(2)=value2
          constructor5%value(3)=value3
          constructor5%value(4)=value4
          constructor5%value(5)=value5
          constructor5%next => next
        END FUNCTION constructor5
        FUNCTION constructor6(value1,value2,value3,value4,value5,value6,next)
          INTEGER,           INTENT(IN   ) :: value1
          INTEGER,           INTENT(IN   ) :: value2
          INTEGER,           INTENT(IN   ) :: value3
          INTEGER,           INTENT(IN   ) :: value4
          INTEGER,           INTENT(IN   ) :: value5
          INTEGER,           INTENT(IN   ) :: value6

          TYPE(ppm_rc_link), POINTER       :: next

          TYPE(ppm_rc_link), POINTER       :: constructor6

          ALLOCATE(constructor6)
          ALLOCATE(constructor6%value(6))
          constructor6%value(1)=value1
          constructor6%value(2)=value2
          constructor6%value(3)=value3
          constructor6%value(4)=value4
          constructor6%value(5)=value5
          constructor6%value(6)=value6
          constructor6%next => next
        END FUNCTION constructor6
        FUNCTION constructor7(value1,value2,value3,value4,value5,value6,value7,next)
          INTEGER,           INTENT(IN   ) :: value1
          INTEGER,           INTENT(IN   ) :: value2
          INTEGER,           INTENT(IN   ) :: value3
          INTEGER,           INTENT(IN   ) :: value4
          INTEGER,           INTENT(IN   ) :: value5
          INTEGER,           INTENT(IN   ) :: value6
          INTEGER,           INTENT(IN   ) :: value7

          TYPE(ppm_rc_link), POINTER       :: next

          TYPE(ppm_rc_link), POINTER       :: constructor7

          ALLOCATE(constructor7)
          ALLOCATE(constructor7%value(7))
          constructor7%value(1)=value1
          constructor7%value(2)=value2
          constructor7%value(3)=value3
          constructor7%value(4)=value4
          constructor7%value(5)=value5
          constructor7%value(6)=value6
          constructor7%value(7)=value7
          constructor7%next => next
        END FUNCTION constructor7
        FUNCTION constructor8(value1,value2,value3,value4,value5,value6,value7,value8,next)
          INTEGER,           INTENT(IN   ) :: value1
          INTEGER,           INTENT(IN   ) :: value2
          INTEGER,           INTENT(IN   ) :: value3
          INTEGER,           INTENT(IN   ) :: value4
          INTEGER,           INTENT(IN   ) :: value5
          INTEGER,           INTENT(IN   ) :: value6
          INTEGER,           INTENT(IN   ) :: value7
          INTEGER,           INTENT(IN   ) :: value8

          TYPE(ppm_rc_link), POINTER       :: next

          TYPE(ppm_rc_link), POINTER       :: constructor8

          ALLOCATE(constructor8)
          ALLOCATE(constructor8%value(8))
          constructor8%value(1)=value1
          constructor8%value(2)=value2
          constructor8%value(3)=value3
          constructor8%value(4)=value4
          constructor8%value(5)=value5
          constructor8%value(6)=value6
          constructor8%value(7)=value7
          constructor8%value(8)=value8
          constructor8%next => next
        END FUNCTION constructor8

        FUNCTION constructorMCMCParticle_2d(value1,value2,MCMCParticle_,next)
          USE ppm_rc_module_util, ONLY : MCMCParticle,IndexHashFunctor_2d
          IMPLICIT NONE
          INTEGER,            INTENT(IN   ) :: value1
          INTEGER,            INTENT(IN   ) :: value2

          TYPE(MCMCParticle), INTENT(IN   ) :: MCMCParticle_

          TYPE(ppm_rc_link),  POINTER       :: next

          TYPE(ppm_rc_link),  POINTER       :: constructorMCMCParticle_2d

          ALLOCATE(constructorMCMCParticle_2d)
          ALLOCATE(constructorMCMCParticle_2d%value(3))
          constructorMCMCParticle_2d%value(1)=IndexHashFunctor_2d(value1,value2)
          constructorMCMCParticle_2d%value(2)=MCMCParticle_%candlabel
          constructorMCMCParticle_2d%value(3)=TRANSFER(MCMCParticle_%proposal,1)
          constructorMCMCParticle_2d%next => next
        END FUNCTION constructorMCMCParticle_2d

        FUNCTION constructorMCMCParticle_3d(value1,value2,value3,MCMCParticle_,next)
          USE ppm_rc_module_util, ONLY : MCMCParticle,IndexHashFunctor_3d
          IMPLICIT NONE
          INTEGER,            INTENT(IN   ) :: value1
          INTEGER,            INTENT(IN   ) :: value2
          INTEGER,            INTENT(IN   ) :: value3

          TYPE(MCMCParticle), INTENT(IN   ) :: MCMCParticle_

          TYPE(ppm_rc_link),  POINTER       :: next

          TYPE(ppm_rc_link),  POINTER       :: constructorMCMCParticle_3d

          ALLOCATE(constructorMCMCParticle_3d)
          ALLOCATE(constructorMCMCParticle_3d%value(3))
          constructorMCMCParticle_3d%value(1)=IndexHashFunctor_3d(value1,value2,value3)
          constructorMCMCParticle_3d%value(2)=MCMCParticle_%candlabel
          constructorMCMCParticle_3d%value(3)=TRANSFER(MCMCParticle_%proposal,1)
          constructorMCMCParticle_3d%next => next
        END FUNCTION constructorMCMCParticle_3d

        FUNCTION constructorMCMCHistoryParticle_2d(value1,value2,MCMCHistoryParticle_,next)
          USE ppm_rc_module_util, ONLY : MCMCHistoryParticle,IndexHashFunctor_2d
          IMPLICIT NONE
          INTEGER,                   INTENT(IN   ) :: value1
          INTEGER,                   INTENT(IN   ) :: value2

          TYPE(MCMCHistoryParticle), INTENT(IN   ) :: MCMCHistoryParticle_

          TYPE(ppm_rc_link),         POINTER       :: next

          TYPE(ppm_rc_link),         POINTER       :: constructorMCMCHistoryParticle_2d

          ALLOCATE(constructorMCMCHistoryParticle_2d)
          ALLOCATE(constructorMCMCHistoryParticle_2d%value(5))
          constructorMCMCHistoryParticle_2d%value(1)=IndexHashFunctor_2d(value1,value2)
          constructorMCMCHistoryParticle_2d%value(2)=MCMCHistoryParticle_%candlabel
          constructorMCMCHistoryParticle_2d%value(3)=MCMCHistoryParticle_%orglabel
          constructorMCMCHistoryParticle_2d%value(4)=TRANSFER(MCMCHistoryParticle_%proposal,1)
          constructorMCMCHistoryParticle_2d%value(5)=TRANSFER(MCMCHistoryParticle_%wasadded,1)
          constructorMCMCHistoryParticle_2d%next => next
        END FUNCTION constructorMCMCHistoryParticle_2d

        FUNCTION constructorMCMCHistoryParticle_3d(value1,value2,value3,MCMCHistoryParticle_,next)
          USE ppm_rc_module_util, ONLY : MCMCHistoryParticle,IndexHashFunctor_3d
          IMPLICIT NONE
          INTEGER,                   INTENT(IN   ) :: value1
          INTEGER,                   INTENT(IN   ) :: value2
          INTEGER,                   INTENT(IN   ) :: value3

          TYPE(MCMCHistoryParticle), INTENT(IN   ) :: MCMCHistoryParticle_

          TYPE(ppm_rc_link),         POINTER       :: next

          TYPE(ppm_rc_link),         POINTER       :: constructorMCMCHistoryParticle_3d

          ALLOCATE(constructorMCMCHistoryParticle_3d)
          ALLOCATE(constructorMCMCHistoryParticle_3d%value(5))
          constructorMCMCHistoryParticle_3d%value(1)=IndexHashFunctor_3d(value1,value2,value3)
          constructorMCMCHistoryParticle_3d%value(2)=MCMCHistoryParticle_%candlabel
          constructorMCMCHistoryParticle_3d%value(3)=MCMCHistoryParticle_%orglabel
          constructorMCMCHistoryParticle_3d%value(4)=TRANSFER(MCMCHistoryParticle_%proposal,1)
          constructorMCMCHistoryParticle_3d%value(5)=TRANSFER(MCMCHistoryParticle_%wasadded,1)
          constructorMCMCHistoryParticle_3d%next => next
        END FUNCTION constructorMCMCHistoryParticle_3d

        FUNCTION constructor_(value,next)
          IMPLICIT NONE
          INTEGER,            INTENT(IN   ) :: value

          TYPE(ppm_rc_link_), POINTER       :: next

          TYPE(ppm_rc_link_), POINTER       :: constructor_

          ALLOCATE(constructor_)
          constructor_%value=value
          constructor_%next => next

        END FUNCTION constructor_
        FUNCTION constructor_2d(value1,value2,next)
          USE ppm_rc_module_util, ONLY : IndexHashFunctor_2d
          IMPLICIT NONE
          INTEGER,            INTENT(IN   ) :: value1
          INTEGER,            INTENT(IN   ) :: value2

          TYPE(ppm_rc_link_), POINTER       :: next

          TYPE(ppm_rc_link_), POINTER       :: constructor_2d

          ALLOCATE(constructor_2d)
          constructor_2d%value=IndexHashFunctor_2d(value1,value2)
          constructor_2d%next => next
        END FUNCTION constructor_2d
        FUNCTION constructor_3d(value1,value2,value3,next)
          USE ppm_rc_module_util, ONLY : IndexHashFunctor_3d
          IMPLICIT NONE
          INTEGER,            INTENT(IN   ) :: value1
          INTEGER,            INTENT(IN   ) :: value2
          INTEGER,            INTENT(IN   ) :: value3

          TYPE(ppm_rc_link_), POINTER       :: next

          TYPE(ppm_rc_link_), POINTER       :: constructor_3d

          ALLOCATE(constructor_3d)
          constructor_3d%value=IndexHashFunctor_3d(value1,value2,value3)
          constructor_3d%next => next
        END FUNCTION constructor_3d

        FUNCTION constructor_2_2d(key1,key2,value,next)
          USE ppm_rc_module_util, ONLY : IndexHashFunctor64_2d
          IMPLICIT NONE
          INTEGER,             INTENT(IN   ) :: key1
          INTEGER,             INTENT(IN   ) :: key2
          INTEGER,             INTENT(IN   ) :: value
          TYPE(ppm_rc_link_2), POINTER       :: next
          TYPE(ppm_rc_link_2), POINTER       :: constructor_2_2d
          ALLOCATE(constructor_2_2d)
          constructor_2_2d%key=IndexHashFunctor64_2d(key1,key2)
          constructor_2_2d%value=value
          constructor_2_2d%next => next
        END FUNCTION constructor_2_2d

        FUNCTION constructor_2_3d(key1,key2,key3,value,next)
          USE ppm_rc_module_util, ONLY : IndexHashFunctor64_3d
          IMPLICIT NONE
          INTEGER,             INTENT(IN   ) :: key1
          INTEGER,             INTENT(IN   ) :: key2
          INTEGER,             INTENT(IN   ) :: key3
          INTEGER,             INTENT(IN   ) :: value
          TYPE(ppm_rc_link_2), POINTER       :: next
          TYPE(ppm_rc_link_2), POINTER       :: constructor_2_3d
          ALLOCATE(constructor_2_3d)
          constructor_2_3d%key=IndexHashFunctor64_3d(key1,key2,key3)
          constructor_2_3d%value=value
          constructor_2_3d%next => next
        END FUNCTION constructor_2_3d

        SUBROUTINE addValue(this,value)
          !!! Add a value to a link (by creating a link) and the link to the list
          CLASS(ppm_rc_list)                   :: this

          INTEGER, DIMENSION(:), INTENT(IN   ) :: value

          TYPE(ppm_rc_link), POINTER :: lastLink
          TYPE(ppm_rc_link), POINTER :: newLink

          IF (.NOT. ASSOCIATED(this%first)) then
             this%first => ppm_rc_link(value, this%first)
             this%last  => this%first
          ELSE
             lastLink => this%last%nextLink()
             newLink  => ppm_rc_link(value,lastLink)
             CALL this%last%setNextLink(newLink)
             this%last => newLink
          ENDIF
        END SUBROUTINE addValue
        SUBROUTINE addValue1(this,value1)
          !!! Add a value to a link (by creating a link) and the link to the list
          CLASS(ppm_rc_list)     :: this

          INTEGER, INTENT(IN   ) :: value1

          TYPE(ppm_rc_link), POINTER :: lastLink
          TYPE(ppm_rc_link), POINTER :: newLink

          IF (.NOT.ASSOCIATED(this%first)) then
             this%first => ppm_rc_link(value1,this%first)
             this%last  => this%first
          ELSE
             lastLink => this%last%nextLink()
             newLink  => ppm_rc_link(value1,lastLink)
             CALL this%last%setNextLink(newLink)
             this%last => newLink
          ENDIF
        END SUBROUTINE addValue1
        SUBROUTINE addValue2(this,value1,value2)
          !!! Add a value to a link (by creating a link) and the link to the list
          CLASS(ppm_rc_list)     :: this

          INTEGER, INTENT(IN   ) :: value1
          INTEGER, INTENT(IN   ) :: value2

          TYPE(ppm_rc_link), POINTER :: lastLink
          TYPE(ppm_rc_link), POINTER :: newLink

          IF (.NOT.ASSOCIATED(this%first)) then
             this%first => ppm_rc_link(value1,value2,this%first)
             this%last  => this%first
          ELSE
             lastLink => this%last%nextLink()
             newLink  => ppm_rc_link(value1,value2,lastLink)
             CALL this%last%setNextLink(newLink)
             this%last => newLink
          ENDIF
        END SUBROUTINE addValue2
        SUBROUTINE addValue3(this,value1,value2,value3)
          !!! Add a value to a link (by creating a link) and the link to the list
          CLASS(ppm_rc_list)     :: this

          INTEGER, INTENT(IN   ) :: value1
          INTEGER, INTENT(IN   ) :: value2
          INTEGER, INTENT(IN   ) :: value3

          TYPE(ppm_rc_link), POINTER :: lastLink
          TYPE(ppm_rc_link), POINTER :: newLink

          IF (.NOT. ASSOCIATED(this%first)) then
             this%first => ppm_rc_link(value1,value2,value3,this%first)
             this%last  => this%first
          ELSE
             lastLink => this%last%nextLink()
             newLink  => ppm_rc_link(value1,value2,value3,lastLink)
             CALL this%last%setNextLink(newLink)
             this%last => newLink
          ENDIF
        END SUBROUTINE addValue3
        SUBROUTINE addValue4(this,value1,value2,value3,value4)
          !!! Add a value to a link (by creating a link) and the link to the list
          CLASS(ppm_rc_list)     :: this

          INTEGER, INTENT(IN   ) :: value1
          INTEGER, INTENT(IN   ) :: value2
          INTEGER, INTENT(IN   ) :: value3
          INTEGER, INTENT(IN   ) :: value4

          TYPE(ppm_rc_link), POINTER :: lastLink
          TYPE(ppm_rc_link), POINTER :: newLink

          IF (.NOT. ASSOCIATED(this%first)) then
             this%first => ppm_rc_link(value1,value2,value3,value4,this%first)
             this%last  => this%first
          ELSE
             lastLink => this%last%nextLink()
             newLink  => ppm_rc_link(value1,value2,value3,value4,lastLink)
             CALL this%last%setNextLink(newLink)
             this%last => newLink
          ENDIF
        END SUBROUTINE addValue4
        SUBROUTINE addValue5(this,value1,value2,value3,value4,value5)
          !!! Add a value to a link (by creating a link) and the link to the list
          CLASS(ppm_rc_list)     :: this

          INTEGER, INTENT(IN   ) :: value1
          INTEGER, INTENT(IN   ) :: value2
          INTEGER, INTENT(IN   ) :: value3
          INTEGER, INTENT(IN   ) :: value4
          INTEGER, INTENT(IN   ) :: value5

          TYPE(ppm_rc_link), POINTER :: lastLink
          TYPE(ppm_rc_link), POINTER :: newLink

          IF (.NOT. ASSOCIATED(this%first)) then
             this%first => ppm_rc_link(value1,value2,value3,value4,value5,this%first)
             this%last  => this%first
          ELSE
             lastLink => this%last%nextLink()
             newLink  => ppm_rc_link(value1,value2,value3,value4,value5,lastLink)
             CALL this%last%setNextLink(newLink)
             this%last => newLink
          ENDIF
        END SUBROUTINE addValue5
        SUBROUTINE addValue6(this,value1,value2,value3,value4,value5,value6)
          !!! Add a value to a link (by creating a link) and the link to the list
          CLASS(ppm_rc_list)     :: this

          INTEGER, INTENT(IN   ) :: value1
          INTEGER, INTENT(IN   ) :: value2
          INTEGER, INTENT(IN   ) :: value3
          INTEGER, INTENT(IN   ) :: value4
          INTEGER, INTENT(IN   ) :: value5
          INTEGER, INTENT(IN   ) :: value6

          TYPE(ppm_rc_link), POINTER :: lastLink
          TYPE(ppm_rc_link), POINTER :: newLink

          IF (.NOT. ASSOCIATED(this%first)) then
             this%first => ppm_rc_link(value1,value2,value3,value4,value5,value6,this%first)
             this%last  => this%first
          ELSE
             lastLink => this%last%nextLink()
             newLink  => ppm_rc_link(value1,value2,value3,value4,value5,value6,lastLink)
             CALL this%last%setNextLink(newLink)
             this%last => newLink
          ENDIF
        END SUBROUTINE addValue6
        SUBROUTINE addValue7(this,value1,value2,value3,value4,value5,value6,value7)
          !!! Add a value to a link (by creating a link) and the link to the list
          CLASS(ppm_rc_list)     :: this

          INTEGER, INTENT(IN   ) :: value1
          INTEGER, INTENT(IN   ) :: value2
          INTEGER, INTENT(IN   ) :: value3
          INTEGER, INTENT(IN   ) :: value4
          INTEGER, INTENT(IN   ) :: value5
          INTEGER, INTENT(IN   ) :: value6
          INTEGER, INTENT(IN   ) :: value7

          TYPE(ppm_rc_link), POINTER :: lastLink
          TYPE(ppm_rc_link), POINTER :: newLink

          IF (.NOT. ASSOCIATED(this%first)) then
             this%first => ppm_rc_link(value1,value2,value3,value4,value5,value6,value7,this%first)
             this%last  => this%first
          ELSE
             lastLink => this%last%nextLink()
             newLink  => ppm_rc_link(value1,value2,value3,value4,value5,value6,value7,lastLink)
             CALL this%last%setNextLink(newLink)
             this%last => newLink
          ENDIF
        END SUBROUTINE addValue7
        SUBROUTINE addValue8(this,value1,value2,value3,value4,value5,value6,value7,value8)
          !!! Add a value to a link (by creating a link) and the link to the list
          CLASS(ppm_rc_list)     :: this

          INTEGER, INTENT(IN   ) :: value1
          INTEGER, INTENT(IN   ) :: value2
          INTEGER, INTENT(IN   ) :: value3
          INTEGER, INTENT(IN   ) :: value4
          INTEGER, INTENT(IN   ) :: value5
          INTEGER, INTENT(IN   ) :: value6
          INTEGER, INTENT(IN   ) :: value7
          INTEGER, INTENT(IN   ) :: value8

          TYPE(ppm_rc_link), POINTER :: lastLink
          TYPE(ppm_rc_link), POINTER :: newLink

          IF (.NOT. ASSOCIATED(this%first)) then
             this%first => ppm_rc_link(value1,value2,value3,value4,value5,value6,value7,value8,this%first)
             this%last  => this%first
          ELSE
             lastLink => this%last%nextLink()
             newLink  => ppm_rc_link(value1,value2,value3,value4,value5,value6,value7,value8,lastLink)
             CALL this%last%setNextLink(newLink)
             this%last => newLink
          ENDIF
        END SUBROUTINE addValue8

        SUBROUTINE addValueMCMCParticle_2d(this,value1,value2,MCMCParticle_)
          !!! Add a value to a link (by creating a link) and the link to the list
          USE ppm_rc_module_util, ONLY : MCMCParticle
          IMPLICIT NONE
          CLASS(ppm_rc_list)                :: this

          INTEGER,            INTENT(IN   ) :: value1
          INTEGER,            INTENT(IN   ) :: value2
          TYPE(MCMCParticle), INTENT(IN   ) :: MCMCParticle_

          TYPE(ppm_rc_link), POINTER :: lastLink
          TYPE(ppm_rc_link), POINTER :: newLink

          IF (.NOT. ASSOCIATED(this%first)) then
             this%first => ppm_rc_link(value1,value2,MCMCParticle_,this%first)
             this%last  => this%first
          ELSE
             lastLink => this%last%nextLink()
             newLink  => ppm_rc_link(value1,value2,MCMCParticle_,lastLink)
             CALL this%last%setNextLink(newLink)
             this%last => newLink
          ENDIF
        END SUBROUTINE addValueMCMCParticle_2d

        SUBROUTINE addValueMCMCParticle_3d(this,value1,value2,value3,MCMCParticle_)
          !!! Add a value to a link (by creating a link) and the link to the list
          USE ppm_rc_module_util, ONLY : MCMCParticle
          IMPLICIT NONE
          CLASS(ppm_rc_list)                :: this

          INTEGER,            INTENT(IN   ) :: value1
          INTEGER,            INTENT(IN   ) :: value2
          INTEGER,            INTENT(IN   ) :: value3
          TYPE(MCMCParticle), INTENT(IN   ) :: MCMCParticle_

          TYPE(ppm_rc_link), POINTER :: lastLink
          TYPE(ppm_rc_link), POINTER :: newLink

          IF (.NOT. ASSOCIATED(this%first)) then
             this%first => ppm_rc_link(value1,value2,value3,MCMCParticle_,this%first)
             this%last  => this%first
          ELSE
             lastLink => this%last%nextLink()
             newLink  => ppm_rc_link(value1,value2,value3,MCMCParticle_,lastLink)
             CALL this%last%setNextLink(newLink)
             this%last => newLink
          ENDIF
        END SUBROUTINE addValueMCMCParticle_3d

        SUBROUTINE addValueMCMCHistoryParticle_2d(this,value1,value2,MCMCHistoryParticle_)
          !!! Add a value to a link (by creating a link) and the link to the list
          USE ppm_rc_module_util, ONLY : MCMCHistoryParticle
          IMPLICIT NONE
          CLASS(ppm_rc_list)                       :: this

          INTEGER,                   INTENT(IN   ) :: value1
          INTEGER,                   INTENT(IN   ) :: value2
          TYPE(MCMCHistoryParticle), INTENT(IN   ) :: MCMCHistoryParticle_

          TYPE(ppm_rc_link), POINTER :: lastLink
          TYPE(ppm_rc_link), POINTER :: newLink

          IF (.NOT. ASSOCIATED(this%first)) then
             this%first => ppm_rc_link(value1,value2,MCMCHistoryParticle_,this%first)
             this%last  => this%first
          ELSE
             lastLink => this%last%nextLink()
             newLink  => ppm_rc_link(value1,value2,MCMCHistoryParticle_,lastLink)
             CALL this%last%setNextLink(newLink)
             this%last => newLink
          ENDIF
        END SUBROUTINE addValueMCMCHistoryParticle_2d

        SUBROUTINE addValueMCMCHistoryParticle_3d(this,value1,value2,value3,MCMCHistoryParticle_)
          !!! Add a value to a link (by creating a link) and the link to the list
          USE ppm_rc_module_util, ONLY : MCMCHistoryParticle
          IMPLICIT NONE
          CLASS(ppm_rc_list)                       :: this

          INTEGER,                   INTENT(IN   ) :: value1
          INTEGER,                   INTENT(IN   ) :: value2
          INTEGER,                   INTENT(IN   ) :: value3
          TYPE(MCMCHistoryParticle), INTENT(IN   ) :: MCMCHistoryParticle_

          TYPE(ppm_rc_link), POINTER :: lastLink
          TYPE(ppm_rc_link), POINTER :: newLink

          IF (.NOT. ASSOCIATED(this%first)) then
             this%first => ppm_rc_link(value1,value2,value3,MCMCHistoryParticle_,this%first)
             this%last  => this%first
          ELSE
             lastLink => this%last%nextLink()
             newLink  => ppm_rc_link(value1,value2,value3,MCMCHistoryParticle_,lastLink)
             CALL this%last%setNextLink(newLink)
             this%last => newLink
          ENDIF
        END SUBROUTINE addValueMCMCHistoryParticle_3d

        SUBROUTINE addValue_(this, value)
          !!! Add a value to a link (by creating a link) and the link to the list
          CLASS(ppm_rc_list_)    :: this

          INTEGER, INTENT(IN   ) :: value

          TYPE(ppm_rc_link_), POINTER :: lastLink
          TYPE(ppm_rc_link_), POINTER :: newLink

          IF (.NOT. ASSOCIATED(this%first)) THEN
             this%first => ppm_rc_link_(value, this%first)
             this%last  => this%first
          ELSE
             lastLink => this%last%nextLink()
             newLink  => ppm_rc_link_(value,lastLink)
             CALL this%last%setNextLink(newLink)
             this%last => newLink
          ENDIF
        END SUBROUTINE addValue_

        SUBROUTINE addValue__(this, value)
          !!! Add a value to a link (by creating a link) and the link to the list
          CLASS(ppm_rc_list_)                  :: this

          INTEGER, DIMENSION(:), INTENT(IN   ) :: value

          TYPE(ppm_rc_link_), POINTER :: lastLink
          TYPE(ppm_rc_link_), POINTER :: newLink

          IF (.NOT. ASSOCIATED(this%first)) THEN
             SELECT CASE (SIZE(value))
             CASE (2)
                this%first => ppm_rc_link_(value(1),value(2),this%first)
             CASE (3)
                this%first => ppm_rc_link_(value(1),value(2),value(3),this%first)
             END SELECT
             this%last => this%first
          ELSE
             lastLink => this%last%nextLink()
             SELECT CASE (SIZE(value))
             CASE (2)
                newLink => ppm_rc_link_(value(1),value(2),lastLink)
             CASE (3)
                newLink => ppm_rc_link_(value(1),value(2),value(3),lastLink)
             END SELECT
             CALL this%last%setNextLink(newLink)
             this%last => newLink
          ENDIF
        END SUBROUTINE addValue__

        SUBROUTINE addValue_2d(this,value1,value2)
          !!! Add a value to a link (by creating a link) and the link to the list
          CLASS(ppm_rc_list_)    :: this

          INTEGER, INTENT(IN   ) :: value1
          INTEGER, INTENT(IN   ) :: value2

          TYPE(ppm_rc_link_), POINTER :: lastLink
          TYPE(ppm_rc_link_), POINTER :: newLink

          IF (.NOT. ASSOCIATED(this%first)) THEN
             this%first => ppm_rc_link_(value1,value2,this%first)
             this%last  => this%first
          ELSE
             lastLink  => this%last%nextLink()
             newLink   => ppm_rc_link_(value1,value2,lastLink)
             CALL this%last%setNextLink(newLink)
             this%last => newLink
          ENDIF
        END SUBROUTINE addValue_2d

        SUBROUTINE addValue_3d(this,value1,value2,value3)
          !!! Add a value to a link (by creating a link) and the link to the list
          CLASS(ppm_rc_list_)    :: this

          INTEGER, INTENT(IN   ) :: value1
          INTEGER, INTENT(IN   ) :: value2
          INTEGER, INTENT(IN   ) :: value3

          TYPE(ppm_rc_link_), POINTER :: lastLink
          TYPE(ppm_rc_link_), POINTER :: newLink

          IF (.NOT. ASSOCIATED(this%first)) THEN
             this%first => ppm_rc_link_(value1,value2,value3,this%first)
             this%last  => this%first
          ELSE
             lastLink  => this%last%nextLink()
             newLink   => ppm_rc_link_(value1,value2,value3,lastLink)
             CALL this%last%setNextLink(newLink)
             this%last => newLink
          ENDIF
        END SUBROUTINE addValue_3d

        SUBROUTINE addValue_2_2d(this,key1,key2,value)
          !!! Add a value to a link (by creating a link) and the link to the list
          CLASS(ppm_rc_list_2)   :: this
          INTEGER, INTENT(IN   ) :: key1
          INTEGER, INTENT(IN   ) :: key2
          INTEGER, INTENT(IN   ) :: value
          TYPE(ppm_rc_link_2), POINTER :: lastLink
          TYPE(ppm_rc_link_2), POINTER :: newLink
          IF (.NOT.ASSOCIATED(this%first)) THEN
             this%first => ppm_rc_link_2(key1,key2,value,this%first)
             this%last  => this%first
          ELSE
             lastLink  => this%last%nextLink()
             newLink   => ppm_rc_link_2(key1,key2,value,lastLink)
             CALL this%last%setNextLink(newLink)
             this%last => newLink
          ENDIF
        END SUBROUTINE addValue_2_2d

        SUBROUTINE addValue_2_3d(this,key1,key2,key3,value)
          !!! Add a value to a link (by creating a link) and the link to the list
          CLASS(ppm_rc_list_2)   :: this
          INTEGER, INTENT(IN   ) :: key1
          INTEGER, INTENT(IN   ) :: key2
          INTEGER, INTENT(IN   ) :: key3
          INTEGER, INTENT(IN   ) :: value
          TYPE(ppm_rc_link_2), POINTER :: lastLink
          TYPE(ppm_rc_link_2), POINTER :: newLink
          IF (.NOT.ASSOCIATED(this%first)) THEN
             this%first => ppm_rc_link_2(key1,key2,key3,value,this%first)
             this%last  => this%first
          ELSE
             lastLink  => this%last%nextLink()
             newLink   => ppm_rc_link_2(key1,key2,key3,value,lastLink)
             CALL this%last%setNextLink(newLink)
             this%last => newLink
          ENDIF
        END SUBROUTINE addValue_2_3d

        SUBROUTINE addValue_2_(this,key,value)
          !!! Add a value to a link (by creating a link) and the link to the list
          CLASS(ppm_rc_list_2)                 :: this
          INTEGER, DIMENSION(:), INTENT(IN   ) :: key
          INTEGER,               INTENT(IN   ) :: value
          TYPE(ppm_rc_link_2), POINTER :: lastLink
          TYPE(ppm_rc_link_2), POINTER :: newLink
          IF (.NOT.ASSOCIATED(this%first)) THEN
             IF (SIZE(key).EQ.2) THEN
                this%first => ppm_rc_link_2(key(1),key(2),value,this%first)
             ELSE
                this%first => ppm_rc_link_2(key(1),key(2),key(3),value,this%first)
             ENDIF
             this%last  => this%first
          ELSE
             lastLink  => this%last%nextLink()
             IF (SIZE(key).EQ.2) THEN
                newLink => ppm_rc_link_2(key(1),key(2),value,lastLink)
             ELSE
                newLink => ppm_rc_link_2(key(1),key(2),key(3),value,lastLink)
             ENDIF
             CALL this%last%setNextLink(newLink)
             this%last => newLink
          ENDIF
        END SUBROUTINE addValue_2_

        SUBROUTINE removeLink(this,link)
          !!! Remove a link from the list (slow...)
          CLASS(ppm_rc_list)                   :: this

          TYPE(ppm_rc_link), OPTIONAL, POINTER :: link

          TYPE(ppm_rc_link), POINTER :: prevLink
          TYPE(ppm_rc_link), POINTER :: currLink

          currLink => this%first
          SELECT CASE (PRESENT(link))
          CASE (.FALSE.)
             this%first => currLink%next
             IF (ASSOCIATED(currLink%value)) DEALLOCATE(currLink%value)
             NULLIFY(currLink%value)
             DEALLOCATE(currLink)
             NULLIFY(currLink)
             RETURN

          CASE (.TRUE.)
             IF (ASSOCIATED(currLink,link)) THEN
                this%first => currLink%next
                IF (ASSOCIATED(currLink%value)) DEALLOCATE(currLink%value)
                NULLIFY(currLink%value)
                DEALLOCATE(currLink)
                NULLIFY(currLink)
                RETURN
             ENDIF
             prevLink => this%first
             DO WHILE(ASSOCIATED(currLink))
                IF (ASSOCIATED(currLink,link)) THEN
                   prevLink%next => currLink%next
                   IF (ASSOCIATED(currLink%value)) DEALLOCATE(currLink%value)
                   NULLIFY(currLink%value)
                   DEALLOCATE(currLink)
                   NULLIFY(currLink)
                   RETURN
                ENDIF
                prevLink => currLink
                currLink => currLink%next
             ENDDO

          END SELECT
        END SUBROUTINE removeLink

        SUBROUTINE removeLink_(this,link)
          !!! Remove a link from the list (slow...)
          CLASS(ppm_rc_list_)                   :: this

          TYPE(ppm_rc_link_), OPTIONAL, POINTER :: link

          TYPE(ppm_rc_link_), POINTER :: prevLink
          TYPE(ppm_rc_link_), POINTER :: currLink

          currLink => this%first
          SELECT CASE (PRESENT(link))
          CASE (.FALSE.)
             this%first => currLink%next
             DEALLOCATE(currLink)
             NULLIFY(currLink)
             RETURN

          CASE (.TRUE.)
             IF (ASSOCIATED(currLink,link)) THEN
                this%first => currLink%next
                DEALLOCATE(currLink)
                NULLIFY(currLink)
                RETURN
             ENDIF
             prevLink => this%first
             DO WHILE(ASSOCIATED(currLink))
                IF (ASSOCIATED(currLink,link)) THEN
                   prevLink%next => currLink%next
                   DEALLOCATE(currLink)
                   NULLIFY(currLink)
                   RETURN
                ENDIF
                prevLink => currLink
                currLink => currLink%next
             ENDDO

          END SELECT
        END SUBROUTINE removeLink_

        SUBROUTINE removeLink_2(this,link)
          !!! Remove a link from the list (slow...)
          CLASS(ppm_rc_list_2)                   :: this
          TYPE(ppm_rc_link_2), OPTIONAL, POINTER :: link
          TYPE(ppm_rc_link_2), POINTER :: prevLink
          TYPE(ppm_rc_link_2), POINTER :: currLink
          currLink => this%first
          SELECT CASE (PRESENT(link))
          CASE (.FALSE.)
             this%first => currLink%next
             DEALLOCATE(currLink)
             NULLIFY(currLink)
             RETURN
          CASE (.TRUE.)
             IF (ASSOCIATED(currLink,link)) THEN
                this%first => currLink%next
                DEALLOCATE(currLink)
                NULLIFY(currLink)
                RETURN
             ENDIF
             prevLink => this%first
             DO WHILE(ASSOCIATED(currLink))
                IF (ASSOCIATED(currLink,link)) THEN
                   prevLink%next => currLink%next
                   DEALLOCATE(currLink)
                   NULLIFY(currLink)
                   RETURN
                ENDIF
                prevLink => currLink
                currLink => currLink%next
             ENDDO
          END SELECT
        END SUBROUTINE removeLink_2

        SUBROUTINE destroy(this,link_)
          !!! Destroy a list
          !!! If there is no link_, it will destroy the whole list
          !!! if there is a link it will destroy the next link afterwards
          CLASS(ppm_rc_list)                   :: this
          TYPE(ppm_rc_link), OPTIONAL, POINTER :: link_

          TYPE(ppm_rc_link), POINTER :: curr
          TYPE(ppm_rc_link), POINTER :: next

          IF (PRESENT(link_)) THEN
             IF (ASSOCIATED(link_)) THEN
                curr => link_%next
             ELSE
                NULLIFY(curr)
             ENDIF
          ELSE
             curr => this%first
          ENDIF

          DO WHILE(ASSOCIATED(curr))
             next => curr%next
             IF (ASSOCIATED(curr%value)) DEALLOCATE(curr%value)
             NULLIFY(curr%value)
             DEALLOCATE(curr)
             NULLIFY(curr)
             curr => next
          ENDDO

          IF (PRESENT(link_)) THEN
             this%last => link_
             NULLIFY(link_%next)
          ELSE
             IF (ASSOCIATED(this%last))  NULLIFY(this%last)
             IF (ASSOCIATED(this%first)) NULLIFY(this%first)
          ENDIF

        END SUBROUTINE destroy

        SUBROUTINE destroy_(this,link_)
          !!! Destroy a list
          !!! If there is no link_, it will destroy the whole list
          !!! if there is a link it will destroy the next link afterwards
          CLASS(ppm_rc_list_)                   :: this
          TYPE(ppm_rc_link_), OPTIONAL, POINTER :: link_

          TYPE(ppm_rc_link_), POINTER :: curr
          TYPE(ppm_rc_link_), POINTER :: next

          IF (PRESENT(link_)) THEN
             IF (ASSOCIATED(link_)) THEN
                curr => link_%next
             ELSE
                NULLIFY(curr)
             ENDIF
          ELSE
             curr => this%first
          ENDIF

          DO WHILE(ASSOCIATED(curr))
             next => curr%next
             DEALLOCATE(curr)
             NULLIFY(curr)
             curr => next
          ENDDO

          IF (PRESENT(link_)) THEN
             this%last => link_
             NULLIFY(link_%next)
          ELSE
             IF (ASSOCIATED(this%last))  NULLIFY(this%last)
             IF (ASSOCIATED(this%first)) NULLIFY(this%first)
          ENDIF

        END SUBROUTINE destroy_

        SUBROUTINE destroy_2(this,link_)
          !!! Destroy a list
          !!! If there is no link_, it will destroy the whole list
          !!! if there is a link it will destroy the next link afterwards
          CLASS(ppm_rc_list_2)                   :: this
          TYPE(ppm_rc_link_2), OPTIONAL, POINTER :: link_
          TYPE(ppm_rc_link_2), POINTER :: curr
          TYPE(ppm_rc_link_2), POINTER :: next
          IF (PRESENT(link_)) THEN
             IF (ASSOCIATED(link_)) THEN
                curr => link_%next
             ELSE
                NULLIFY(curr)
             ENDIF
          ELSE
             curr => this%first
          ENDIF
          DO WHILE(ASSOCIATED(curr))
             next => curr%next
             DEALLOCATE(curr)
             NULLIFY(curr)
             curr => next
          ENDDO
          IF (PRESENT(link_)) THEN
             this%last => link_
             NULLIFY(link_%next)
          ELSE
             IF (ASSOCIATED(this%last))  NULLIFY(this%last)
             IF (ASSOCIATED(this%first)) NULLIFY(this%first)
          ENDIF
        END SUBROUTINE destroy_2

        SUBROUTINE mergeLinks(this,links)
          !!! Merge links from the list to the available list
          CLASS(ppm_rc_list)         :: this

          TYPE(ppm_rc_list), POINTER :: links

          TYPE(ppm_rc_link), POINTER :: currLink

          IF (ASSOCIATED(links%first)) THEN
             currLink => links%first
             CALL this%last%setNextLink(currLink)
             this%last => links%last
             NULLIFY(links%first,links%last)
             DEALLOCATE(links)
             NULLIFY(links)
          ENDIF

        END SUBROUTINE mergeLinks

        SUBROUTINE mergeLinks_(this,links)
          !!! Merge links from the list to the available list
          CLASS(ppm_rc_list_)         :: this

          TYPE(ppm_rc_list_), POINTER :: links

          TYPE(ppm_rc_link_), POINTER :: currLink

          IF (ASSOCIATED(links%first)) THEN
             currLink => links%first
             CALL this%last%setNextLink(currLink)
             this%last => links%last
             NULLIFY(links%first,links%last)
             DEALLOCATE(links)
             NULLIFY(links)
          ENDIF

        END SUBROUTINE mergeLinks_

        SUBROUTINE mergeLinks_2(this,links)
          !!! Merge links from the list to the available list
          CLASS(ppm_rc_list_2)         :: this
          TYPE(ppm_rc_list_2), POINTER :: links
          TYPE(ppm_rc_link_2), POINTER :: currLink
          IF (ASSOCIATED(links%first)) THEN
             currLink => links%first
             CALL this%last%setNextLink(currLink)
             this%last => links%last
             NULLIFY(links%first,links%last)
             DEALLOCATE(links)
             NULLIFY(links)
          ENDIF
        END SUBROUTINE mergeLinks_2

        SUBROUTINE swap1(this,curr)
          !!! swap a link with the root
          CLASS(ppm_rc_list)         :: this

          TYPE(ppm_rc_link), POINTER :: curr

!           INTEGER, DIMENSION(SIZE(curr%value)) :: tmp
          INTEGER, DIMENSION(:), POINTER :: tmp

          LOGICAL, SAVE :: new=.TRUE.

          IF (new) THEN
             new=.FALSE.
             tmp              => this%first%value
             this%first%value => curr%value
             curr%value       => tmp
          ENDIF

        END SUBROUTINE swap1

        SUBROUTINE swap1_(this,curr)
          !!! swap a link with the root
          CLASS(ppm_rc_list_)         :: this

          TYPE(ppm_rc_link_), POINTER :: curr

          INTEGER :: tmp

          LOGICAL, SAVE :: new=.TRUE.

          IF (new) THEN
             new=.FALSE.
             tmp              = this%first%value
             this%first%value = curr%value
             curr%value       = tmp
          ENDIF

        END SUBROUTINE swap1_

        SUBROUTINE swap1_2(this,curr)
          !!! swap a link with the root
          CLASS(ppm_rc_list_2)         :: this
          TYPE(ppm_rc_link_2), POINTER :: curr
          INTEGER(ppm_kind_int64) :: tmpkey
          INTEGER                 :: tmpvalue
          LOGICAL, SAVE :: new=.TRUE.
          IF (new) THEN
             new=.FALSE.
             tmpkey           = this%first%key
             tmpvalue         = this%first%value
             this%first%key   = curr%key
             this%first%value = curr%value
             curr%key         = tmpkey
             curr%value       = tmpvalue
          ENDIF
        END SUBROUTINE swap1_2

        SUBROUTINE swap2(this,curr1,curr2)
          !!! swap a link with the root
          CLASS(ppm_rc_list)         :: this

          TYPE(ppm_rc_link), POINTER :: curr1
          TYPE(ppm_rc_link), POINTER :: curr2

!           INTEGER, DIMENSION(SIZE(curr1%value)) :: tmp
          INTEGER, DIMENSION(:), POINTER :: tmp

          tmp => curr1%value
          curr1%value => curr2%value
          curr2%value => tmp

        END SUBROUTINE swap2

        SUBROUTINE swap2_(this,curr1,curr2)
          !!! swap a link with the root
          CLASS(ppm_rc_list_)         :: this

          TYPE(ppm_rc_link_), POINTER :: curr1
          TYPE(ppm_rc_link_), POINTER :: curr2

          INTEGER :: tmp

          tmp = curr1%value
          curr1%value = curr2%value
          curr2%value = tmp

        END SUBROUTINE swap2_

        SUBROUTINE swap2_2(this,curr1,curr2)
          !!! swap a link with the root
          CLASS(ppm_rc_list_2)         :: this
          TYPE(ppm_rc_link_2), POINTER :: curr1
          TYPE(ppm_rc_link_2), POINTER :: curr2
          INTEGER(ppm_kind_int64) :: tmpkey
          INTEGER                 :: tmpvalue
          tmpkey      =curr1%key
          tmpvalue    = curr1%value
          curr1%key   = curr2%key
          curr1%value = curr2%value
          curr2%key   = tmpkey
          curr2%value = tmpvalue
        END SUBROUTINE swap2_2

        !BEGIN
        FUNCTION ppm_rc_c_list_begin(this) RESULT (iterator)

            IMPLICIT NONE
            !-------------------------------------------------------------------------
            !  Arguments
            !-------------------------------------------------------------------------
            CLASS(ppm_rc_c_list), INTENT(INOUT) :: this

            TYPE(ppm_rc_list),    POINTER       :: iterator

            IF (this%nb.GT.0) THEN
               this%iter_id = this%min_id
               IF (this%iter_id.LE.this%max_id) THEN
                  iterator => this%vec(this%min_id)%t
                  RETURN
               ENDIF
            ENDIF
            iterator => NULL()
        END FUNCTION ppm_rc_c_list_begin
        FUNCTION ppm_rc_c_list_begin_(this) RESULT (iterator)

            IMPLICIT NONE
            !-------------------------------------------------------------------------
            !  Arguments
            !-------------------------------------------------------------------------
            CLASS(ppm_rc_c_list_), INTENT(INOUT) :: this

            TYPE(ppm_rc_list_),    POINTER       :: iterator

            IF (this%nb.GT.0) THEN
               this%iter_id = this%min_id
               IF (this%iter_id.LE.this%max_id) THEN
                  iterator => this%vec(this%min_id)%t
                  RETURN
               ENDIF
            ENDIF
            iterator => NULL()
        END FUNCTION ppm_rc_c_list_begin_
        FUNCTION ppm_rc_c_list_2_begin(this) RESULT (iterator)
            IMPLICIT NONE
            !-------------------------------------------------------------------------
            !  Arguments
            !-------------------------------------------------------------------------
            CLASS(ppm_rc_c_list_2), INTENT(INOUT) :: this
            TYPE(ppm_rc_list_2),    POINTER       :: iterator
            IF (this%nb.GT.0) THEN
               this%iter_id = this%min_id
               IF (this%iter_id.LE.this%max_id) THEN
                  iterator => this%vec(this%min_id)%t
                  RETURN
               ENDIF
            ENDIF
            iterator => NULL()
        END FUNCTION ppm_rc_c_list_2_begin

        !NEXT
        FUNCTION ppm_rc_c_list_next(this) RESULT (iterator)

            IMPLICIT NONE
            !-------------------------------------------------------------------------
            !  Arguments
            !-------------------------------------------------------------------------
            CLASS(ppm_rc_c_list), INTENT(INOUT) :: this

            TYPE(ppm_rc_list),    POINTER       :: iterator

            IF (this%nb.GT.0) THEN
               this%iter_id = this%iter_id + 1
               IF (this%iter_id.GE.this%min_id .AND. this%iter_id.LE.this%max_id) THEN
                  iterator => this%vec(this%iter_id)%t
                  RETURN
               ENDIF
            ENDIF
            iterator => NULL()
        END FUNCTION ppm_rc_c_list_next
        FUNCTION ppm_rc_c_list_next_(this) RESULT (iterator)

            IMPLICIT NONE
            !-------------------------------------------------------------------------
            !  Arguments
            !-------------------------------------------------------------------------
            CLASS(ppm_rc_c_list_), INTENT(INOUT) :: this

            TYPE(ppm_rc_list_),    POINTER       :: iterator

            IF (this%nb.GT.0) THEN
               this%iter_id = this%iter_id + 1
               IF (this%iter_id.GE.this%min_id .AND. this%iter_id.LE.this%max_id) THEN
                  iterator => this%vec(this%iter_id)%t
                  RETURN
               ENDIF
            ENDIF
            iterator => NULL()
        END FUNCTION ppm_rc_c_list_next_
        FUNCTION ppm_rc_c_list_2_next(this) RESULT (iterator)
            IMPLICIT NONE
            !-------------------------------------------------------------------------
            !  Arguments
            !-------------------------------------------------------------------------
            CLASS(ppm_rc_c_list_2), INTENT(INOUT) :: this
            TYPE(ppm_rc_list_2),    POINTER       :: iterator
            IF (this%nb.GT.0) THEN
               this%iter_id = this%iter_id + 1
               IF (this%iter_id.GE.this%min_id .AND. this%iter_id.LE.this%max_id) THEN
                  iterator => this%vec(this%iter_id)%t
                  RETURN
               ENDIF
            ENDIF
            iterator => NULL()
        END FUNCTION ppm_rc_c_list_2_next

        !PREVIOUS
        FUNCTION ppm_rc_c_list_prev(this) RESULT (iterator)

            IMPLICIT NONE
            !-------------------------------------------------------------------------
            !  Arguments
            !-------------------------------------------------------------------------
            CLASS(ppm_rc_c_list), INTENT(INOUT) :: this

            TYPE(ppm_rc_list),    POINTER       :: iterator

            IF (this%nb.GT.0) THEN
               this%iter_id = this%iter_id - 1
               IF (this%iter_id.GE.this%min_id .AND. this%iter_id.LE.this%max_id) THEN
                  iterator => this%vec(this%iter_id)%t
                  RETURN
               ENDIF
            ENDIF
            iterator => NULL()
        END FUNCTION ppm_rc_c_list_prev
        FUNCTION ppm_rc_c_list_prev_(this) RESULT (iterator)

            IMPLICIT NONE
            !-------------------------------------------------------------------------
            !  Arguments
            !-------------------------------------------------------------------------
            CLASS(ppm_rc_c_list_), INTENT(INOUT) :: this

            TYPE(ppm_rc_list_),    POINTER       :: iterator

            IF (this%nb.GT.0) THEN
               this%iter_id = this%iter_id - 1
               IF (this%iter_id.GE.this%min_id .AND. this%iter_id.LE.this%max_id) THEN
                  iterator => this%vec(this%iter_id)%t
                  RETURN
               ENDIF
            ENDIF
            iterator => NULL()
        END FUNCTION ppm_rc_c_list_prev_
        FUNCTION ppm_rc_c_list_2_prev(this) RESULT (iterator)
            IMPLICIT NONE
            !-------------------------------------------------------------------------
            !  Arguments
            !-------------------------------------------------------------------------
            CLASS(ppm_rc_c_list_2), INTENT(INOUT) :: this
            TYPE(ppm_rc_list_2),    POINTER       :: iterator
            IF (this%nb.GT.0) THEN
               this%iter_id = this%iter_id - 1
               IF (this%iter_id.GE.this%min_id .AND. this%iter_id.LE.this%max_id) THEN
                  iterator => this%vec(this%iter_id)%t
                  RETURN
               ENDIF
            ENDIF
            iterator => NULL()
        END FUNCTION ppm_rc_c_list_2_prev
        !LAST
        FUNCTION ppm_rc_c_list_last(this) RESULT (iterator)

            IMPLICIT NONE
            !-------------------------------------------------------------------------
            !  Arguments
            !-------------------------------------------------------------------------
            CLASS(ppm_rc_c_list), INTENT(INOUT) :: this

            TYPE(ppm_rc_list),    POINTER       :: iterator

            IF (this%nb.GT.0) THEN
               this%iter_id = this%max_id
               IF (this%iter_id.GE.this%min_id) THEN
                  iterator => this%vec(this%max_id)%t
                  RETURN
               ENDIF
            ENDIF
            iterator => NULL()
        END FUNCTION ppm_rc_c_list_last
        FUNCTION ppm_rc_c_list_last_(this) RESULT (iterator)

            IMPLICIT NONE
            !-------------------------------------------------------------------------
            !  Arguments
            !-------------------------------------------------------------------------
            CLASS(ppm_rc_c_list_), INTENT(INOUT) :: this

            TYPE(ppm_rc_list_),    POINTER       :: iterator

            IF (this%nb.GT.0) THEN
               this%iter_id = this%max_id
               IF (this%iter_id.GE.this%min_id) THEN
                  iterator => this%vec(this%max_id)%t
                  RETURN
               ENDIF
            ENDIF
            iterator => NULL()
        END FUNCTION ppm_rc_c_list_last_
        FUNCTION ppm_rc_c_list_2_last(this) RESULT (iterator)
            IMPLICIT NONE
            !-------------------------------------------------------------------------
            !  Arguments
            !-------------------------------------------------------------------------
            CLASS(ppm_rc_c_list_2), INTENT(INOUT) :: this
            TYPE(ppm_rc_list_2),    POINTER       :: iterator
            IF (this%nb.GT.0) THEN
               this%iter_id = this%max_id
               IF (this%iter_id.GE.this%min_id) THEN
                  iterator => this%vec(this%max_id)%t
                  RETURN
               ENDIF
            ENDIF
            iterator => NULL()
        END FUNCTION ppm_rc_c_list_2_last
        !DESTROY CONTAINER
        SUBROUTINE ppm_rc_c_list_destroy(this,info)

            IMPLICIT NONE
            !-------------------------------------------------------------------------
            !  Arguments
            !-------------------------------------------------------------------------
            CLASS(ppm_rc_c_list), INTENT(INOUT) :: this

            INTEGER,              INTENT(  OUT) :: info

            !-------------------------------------------------------------------------
            ! local variables
            !-------------------------------------------------------------------------
            TYPE(ppm_rc_list), POINTER :: p

            start_subroutine("ppm_rc_c_list_destroy")

            p => this%begin()
            DO WHILE(ASSOCIATED(p))
               !yaser I think this is necessary for complete destroyment
               CALL p%destroy()

               DEALLOCATE(P,STAT=info)
               or_fail_dealloc("Could not deallocate collection element")

               p => this%next()
            ENDDO

            IF (ASSOCIATED(this%vec)) THEN
               DEALLOCATE(this%vec,STAT=info)
               or_fail_dealloc("Could not deallocate collection array")
            ENDIF

            this%iter_id = 0
            this%min_id = 0
            this%max_id = 0
            this%nb = 0
            this%vec_size=0

            end_subroutine()
        END SUBROUTINE ppm_rc_c_list_destroy

        SUBROUTINE ppm_rc_c_list_destroy_(this,info)

            IMPLICIT NONE
            !-------------------------------------------------------------------------
            !  Arguments
            !-------------------------------------------------------------------------
            CLASS(ppm_rc_c_list_), INTENT(INOUT) :: this

            INTEGER,              INTENT(  OUT) :: info

            !-------------------------------------------------------------------------
            ! local variables
            !-------------------------------------------------------------------------
            TYPE(ppm_rc_list_), POINTER :: p

            start_subroutine("ppm_rc_c_list_destroy")

            p => this%begin()
            DO WHILE(ASSOCIATED(p))
               !yaser I think this is necessary for complete destroyment
               CALL p%destroy()

               DEALLOCATE(P,STAT=info)
               or_fail_dealloc("Could not deallocate collection element")

               p => this%next()
            ENDDO

            IF (ASSOCIATED(this%vec)) THEN
               DEALLOCATE(this%vec,STAT=info)
               or_fail_dealloc("Could not deallocate collection array")
            ENDIF

            this%iter_id = 0
            this%min_id = 0
            this%max_id = 0
            this%nb = 0
            this%vec_size=0

            end_subroutine()
        END SUBROUTINE ppm_rc_c_list_destroy_
        SUBROUTINE ppm_rc_c_list_2_destroy(this,info)
            IMPLICIT NONE
            !-------------------------------------------------------------------------
            !  Arguments
            !-------------------------------------------------------------------------
            CLASS(ppm_rc_c_list_2), INTENT(INOUT) :: this
            INTEGER,                INTENT(  OUT) :: info
            !-------------------------------------------------------------------------
            ! local variables
            !-------------------------------------------------------------------------
            TYPE(ppm_rc_list_2), POINTER :: p
            start_subroutine("ppm_rc_c_list_destroy")
            p => this%begin()
            DO WHILE(ASSOCIATED(p))
               !yaser I think this is necessary for complete destroyment
               CALL p%destroy()
               DEALLOCATE(P,STAT=info)
               or_fail_dealloc("Could not deallocate collection element")
               p => this%next()
            ENDDO
            IF (ASSOCIATED(this%vec)) THEN
               DEALLOCATE(this%vec,STAT=info)
               or_fail_dealloc("Could not deallocate collection array")
            ENDIF
            this%iter_id = 0
            this%min_id = 0
            this%max_id = 0
            this%nb = 0
            this%vec_size=0
            end_subroutine()
        END SUBROUTINE ppm_rc_c_list_2_destroy

        !EXISTS
        FUNCTION ppm_rc_c_list_exists(this,id) RESULT(exists)
            !!! Check whether an element exists and can be accessed at this id

            IMPLICIT NONE
            !-------------------------------------------------------------------------
            !  Arguments
            !-------------------------------------------------------------------------
            CLASS(ppm_rc_c_list), INTENT(IN   ) :: this
            !!! Data structure containing the particles
            INTEGER,              INTENT(IN   ) :: id
            !!! id where the data is stored
            LOGICAL                             :: exists
            !!! Return status, on success 0.

            !-------------------------------------------------------------------------
            ! local variables
            !-------------------------------------------------------------------------

            !-------------------------------------------------------------------------
            ! Check arguments
            !-------------------------------------------------------------------------
            IF (id.LE.0 .OR. id.LT.this%min_id .OR. id.GT.this%max_id) THEN
               exists = .FALSE.
               RETURN
            ENDIF
            IF (.NOT.ASSOCIATED(this%vec)) THEN
               exists = .FALSE.
               RETURN
            ENDIF
            IF (.NOT.ASSOCIATED(this%vec(id)%t)) THEN
               exists = .FALSE.
               RETURN
            ENDIF
            exists = .TRUE.

            RETURN
        END FUNCTION ppm_rc_c_list_exists

        FUNCTION ppm_rc_c_list_exists_(this,id) RESULT(exists)
            !!! Check whether an element exists and can be accessed at this id

            IMPLICIT NONE
            !-------------------------------------------------------------------------
            !  Arguments
            !-------------------------------------------------------------------------
            CLASS(ppm_rc_c_list_), INTENT(IN   ) :: this
            !!! Data structure containing the particles
            INTEGER,               INTENT(IN   ) :: id
            !!! id where the data is stored
            LOGICAL                              :: exists
            !!! Return status, on success 0.

            !-------------------------------------------------------------------------
            ! local variables
            !-------------------------------------------------------------------------

            !-------------------------------------------------------------------------
            ! Check arguments
            !-------------------------------------------------------------------------
            IF (id.LE.0 .OR. id.LT.this%min_id .OR. id.GT.this%max_id) THEN
               exists = .FALSE.
               RETURN
            ENDIF
            IF (.NOT.ASSOCIATED(this%vec)) THEN
               exists = .FALSE.
               RETURN
            ENDIF
            IF (.NOT.ASSOCIATED(this%vec(id)%t)) THEN
               exists = .FALSE.
               RETURN
            ENDIF
            exists = .TRUE.

            RETURN
        END FUNCTION ppm_rc_c_list_exists_

        FUNCTION ppm_rc_c_list_2_exists(this,id) RESULT(exists)
            !!! Check whether an element exists and can be accessed at this id
            IMPLICIT NONE
            !-------------------------------------------------------------------------
            !  Arguments
            !-------------------------------------------------------------------------
            CLASS(ppm_rc_c_list_2), INTENT(IN   ) :: this
            !!! Data structure containing the particles
            INTEGER,                INTENT(IN   ) :: id
            !!! id where the data is stored
            LOGICAL                               :: exists
            !!! Return status, on success 0.
            !-------------------------------------------------------------------------
            ! local variables
            !-------------------------------------------------------------------------
            !-------------------------------------------------------------------------
            ! Check arguments
            !-------------------------------------------------------------------------
            IF (id.LE.0 .OR. id.LT.this%min_id .OR. id.GT.this%max_id) THEN
               exists = .FALSE.
               RETURN
            ENDIF
            IF (.NOT.ASSOCIATED(this%vec)) THEN
               exists = .FALSE.
               RETURN
            ENDIF
            IF (.NOT.ASSOCIATED(this%vec(id)%t)) THEN
               exists = .FALSE.
               RETURN
            ENDIF
            exists = .TRUE.
            RETURN
        END FUNCTION ppm_rc_c_list_2_exists
        !HAS
        FUNCTION ppm_rc_c_list_has(this,element) RESULT(has)
            !!! Check whether an element is present in the collection (slow...)

            IMPLICIT NONE
            !-------------------------------------------------------------------------
            !  Arguments
            !-------------------------------------------------------------------------
            CLASS(ppm_rc_c_list),       INTENT(INOUT) :: this
            !!! Data structure containing the particles
            TYPE(ppm_rc_list),  TARGET, INTENT(IN   ) :: element
            !!! element which is being searched for
            LOGICAL                                   :: has
            !!! true if this element belongs to the collection

            !-------------------------------------------------------------------------
            ! local variables
            !-------------------------------------------------------------------------
            TYPE(ppm_rc_list), POINTER :: p

            has = .TRUE.

            p => this%begin()
            DO WHILE(ASSOCIATED(p))
               IF (ASSOCIATED(p,element)) RETURN
               p => this%next()
            ENDDO

            has = .FALSE.

            RETURN
        END FUNCTION ppm_rc_c_list_has

        FUNCTION ppm_rc_c_list_has_(this,element) RESULT(has)
            !!! Check whether an element is present in the collection (slow...)

            IMPLICIT NONE
            !-------------------------------------------------------------------------
            !  Arguments
            !-------------------------------------------------------------------------
            CLASS(ppm_rc_c_list_),       INTENT(INOUT) :: this
            !!! Data structure containing the particles
            TYPE(ppm_rc_list_),  TARGET, INTENT(IN   ) :: element
            !!! element which is being searched for
            LOGICAL                                    :: has
            !!! true if this element belongs to the collection

            !-------------------------------------------------------------------------
            ! local variables
            !-------------------------------------------------------------------------
            TYPE(ppm_rc_list_), POINTER :: p

            has = .TRUE.

            p => this%begin()
            DO WHILE(ASSOCIATED(p))
               IF (ASSOCIATED(p,element)) RETURN
               p => this%next()
            ENDDO

            has = .FALSE.

            RETURN
        END FUNCTION ppm_rc_c_list_has_

        FUNCTION ppm_rc_c_list_2_has(this,element) RESULT(has)
            !!! Check whether an element is present in the collection (slow...)
            IMPLICIT NONE
            !-------------------------------------------------------------------------
            !  Arguments
            !-------------------------------------------------------------------------
            CLASS(ppm_rc_c_list_2),       INTENT(INOUT) :: this
            !!! Data structure containing the particles
            TYPE(ppm_rc_list_2),  TARGET, INTENT(IN   ) :: element
            !!! element which is being searched for
            LOGICAL                                     :: has
            !!! true if this element belongs to the collection
            !-------------------------------------------------------------------------
            ! local variables
            !-------------------------------------------------------------------------
            TYPE(ppm_rc_list_2), POINTER :: p
            has = .TRUE.
            p => this%begin()
            DO WHILE(ASSOCIATED(p))
               IF (ASSOCIATED(p,element)) RETURN
               p => this%next()
            ENDDO
            has = .FALSE.
            RETURN
        END FUNCTION ppm_rc_c_list_2_has

        !GET_ID
        FUNCTION ppm_rc_c_list_get_id(this,element) RESULT(id)
            !!! Returns the id of an element in the collection (slow...)
            !!! Returs -1 if the element is not found.

            IMPLICIT NONE
            !-------------------------------------------------------------------------
            !  Arguments
            !-------------------------------------------------------------------------
            CLASS(ppm_rc_c_list),      INTENT(INOUT) :: this
            !!! Data structure containing the particles
            TYPE(ppm_rc_list), TARGET, INTENT(IN   ) :: element
            !!! Element which is being searched for
            INTEGER                                  :: id
            !!! id where the data is stored

            !-------------------------------------------------------------------------
            ! local variables
            !-------------------------------------------------------------------------
            TYPE(ppm_rc_list), POINTER :: p

            p => this%begin()
            DO WHILE(ASSOCIATED(p))
               IF (ASSOCIATED(p,element)) THEN
                  id = this%iter_id
                  RETURN
               ENDIF
               p => this%next()
            ENDDO

            id = -1

            RETURN
        END FUNCTION ppm_rc_c_list_get_id

        FUNCTION ppm_rc_c_list_get_id_(this,element) RESULT(id)
            !!! Returns the id of an element in the collection (slow...)
            !!! Returs -1 if the element is not found.

            IMPLICIT NONE
            !-------------------------------------------------------------------------
            !  Arguments
            !-------------------------------------------------------------------------
            CLASS(ppm_rc_c_list_),      INTENT(INOUT) :: this
            !!! Data structure containing the particles
            TYPE(ppm_rc_list_), TARGET, INTENT(IN   ) :: element
            !!! Element which is being searched for
            INTEGER                                   :: id
            !!! id where the data is stored

            !-------------------------------------------------------------------------
            ! local variables
            !-------------------------------------------------------------------------
            TYPE(ppm_rc_list_), POINTER :: p

            p => this%begin()
            DO WHILE(ASSOCIATED(p))
               IF (ASSOCIATED(p,element)) THEN
                  id = this%iter_id
                  RETURN
               ENDIF
               p => this%next()
            ENDDO

            id = -1

            RETURN
        END FUNCTION ppm_rc_c_list_get_id_

        FUNCTION ppm_rc_c_list_2_get_id(this,element) RESULT(id)
            !!! Returns the id of an element in the collection (slow...)
            !!! Returs -1 if the element is not found.
            IMPLICIT NONE
            !-------------------------------------------------------------------------
            !  Arguments
            !-------------------------------------------------------------------------
            CLASS(ppm_rc_c_list_2),      INTENT(INOUT) :: this
            !!! Data structure containing the particles
            TYPE(ppm_rc_list_2), TARGET, INTENT(IN   ) :: element
            !!! Element which is being searched for
            INTEGER                                    :: id
            !!! id where the data is stored
            !-------------------------------------------------------------------------
            ! local variables
            !-------------------------------------------------------------------------
            TYPE(ppm_rc_list_2), POINTER :: p
            p => this%begin()
            DO WHILE(ASSOCIATED(p))
               IF (ASSOCIATED(p,element)) THEN
                  id = this%iter_id
                  RETURN
               ENDIF
               p => this%next()
            ENDDO
            id = -1
            RETURN
        END FUNCTION ppm_rc_c_list_2_get_id

        !PUSH
        SUBROUTINE ppm_rc_c_list_push(this,element,info,id)
            !!! add an element into the collection

            IMPLICIT NONE
            !-------------------------------------------------------------------------
            !  Arguments
            !-------------------------------------------------------------------------
            CLASS(ppm_rc_c_list),           INTENT(INOUT) :: this
            TYPE(ppm_rc_list),    POINTER                 :: element
            INTEGER,                        INTENT(  OUT) :: info
            INTEGER,              OPTIONAL, INTENT(  OUT) :: id
            !!! index of the element in the collection

            !-------------------------------------------------------------------------
            ! local variables
            !-------------------------------------------------------------------------
            start_subroutine("ppm_rc_c_list_push")

            !add the element at the end of the array
            this%min_id = 1
            this%max_id = this%max_id + 1
            this%nb = this%nb + 1
            IF (PRESENT(id)) id = this%max_id

            IF (this%max_id.GT.this%vec_size) THEN
               CALL this%grow_size(info)
               or_fail("could not grow ppm_rc_c_list to a larger size")
            ENDIF

            IF (ASSOCIATED(this%vec(this%max_id)%t)) THEN
               fail("Pointer at position of new element is already associated. Something wrong in the Collection data structure")
            ENDIF

            this%vec(this%max_id)%t => element

            check_associated_noscope(<#this%vec(this%max_id)%t#>,"Pushing element into collection failed unexpectedly")

            element => NULL()

            end_subroutine()
        END SUBROUTINE ppm_rc_c_list_push
        SUBROUTINE ppm_rc_c_list_push_(this,element,info,id)
            !!! add an element into the collection

            IMPLICIT NONE
            !-------------------------------------------------------------------------
            !  Arguments
            !-------------------------------------------------------------------------
            CLASS(ppm_rc_c_list_),           INTENT(INOUT) :: this
            TYPE(ppm_rc_list_),    POINTER                 :: element
            INTEGER,                         INTENT(  OUT) :: info
            INTEGER,               OPTIONAL, INTENT(  OUT) :: id
            !!! index of the element in the collection

            !-------------------------------------------------------------------------
            ! local variables
            !-------------------------------------------------------------------------
            start_subroutine("ppm_rc_c_list_push")

            !add the element at the end of the array
            this%min_id = 1
            this%max_id = this%max_id + 1
            this%nb = this%nb + 1
            IF (PRESENT(id)) id = this%max_id

            IF (this%max_id.GT.this%vec_size) THEN
               CALL this%grow_size(info)
               or_fail("could not grow ppm_rc_c_list to a larger size")
            ENDIF

            IF (ASSOCIATED(this%vec(this%max_id)%t)) THEN
               fail("Pointer at position of new element is already associated. Something wrong in the Collection data structure")
            ENDIF

            this%vec(this%max_id)%t => element

            check_associated_noscope(<#this%vec(this%max_id)%t#>,"Pushing element into collection failed unexpectedly")

            element => NULL()

            end_subroutine()
        END SUBROUTINE ppm_rc_c_list_push_

        SUBROUTINE ppm_rc_c_list_2_push(this,element,info,id)
            !!! add an element into the collection
            IMPLICIT NONE
            !-------------------------------------------------------------------------
            !  Arguments
            !-------------------------------------------------------------------------
            CLASS(ppm_rc_c_list_2),           INTENT(INOUT) :: this
            TYPE(ppm_rc_list_2),    POINTER                 :: element
            INTEGER,                          INTENT(  OUT) :: info
            INTEGER,                OPTIONAL, INTENT(  OUT) :: id
            !!! index of the element in the collection
            !-------------------------------------------------------------------------
            ! local variables
            !-------------------------------------------------------------------------
            start_subroutine("ppm_rc_c_list_push")
            !add the element at the end of the array
            this%min_id = 1
            this%max_id = this%max_id + 1
            this%nb = this%nb + 1
            IF (PRESENT(id)) id = this%max_id
            IF (this%max_id.GT.this%vec_size) THEN
               CALL this%grow_size(info)
               or_fail("could not grow ppm_rc_c_list to a larger size")
            ENDIF
            IF (ASSOCIATED(this%vec(this%max_id)%t)) THEN
               fail("Pointer at position of new element is already associated. Something wrong in the Collection data structure")
            ENDIF
            this%vec(this%max_id)%t => element
            check_associated_noscope(<#this%vec(this%max_id)%t#>,"Pushing element into collection failed unexpectedly")
            element => NULL()
            end_subroutine()
        END SUBROUTINE ppm_rc_c_list_2_push
        !REMOVE
        SUBROUTINE ppm_rc_c_list_remove(this,info,element,id)
            !!! If element is present, remove it from the collection
            !!! else, remove the current element (as defined by the iterator pointer)

            IMPLICIT NONE
            !-------------------------------------------------------------------------
            !  Arguments
            !-------------------------------------------------------------------------
            CLASS(ppm_rc_c_list),         INTENT(INOUT) :: this

            INTEGER,                      INTENT(  OUT) :: info

            TYPE(ppm_rc_list), OPTIONAL, INTENT(INOUT) :: element
            INTEGER,           OPTIONAL, INTENT(IN   ) :: id

            !-------------------------------------------------------------------------
            ! local variables
            !-------------------------------------------------------------------------
            INTEGER :: del_id
            INTEGER :: iter_id_save

            iter_id_save = this%iter_id

            IF (PRESENT(element)) THEN
               del_id = MERGE(id,this%get_id(element),PRESENT(id))
               IF (del_id.LT.0) RETURN
            ELSE
               del_id = MERGE(id,this%iter_id,PRESENT(id))
            ENDIF

            !deallocate the element
            CALL this%vec(del_id)%t%destroy()

            IF (PRESENT(id)) THEN
               IF (del_id+1.LE.this%max_id) THEN
                  this%vec(del_id)%t => this%vec(del_id+1)%t
               ELSE
                  this%vec(del_id)%t => NULL()
                  this%nb = this%nb - 1
                  this%max_id = this%max_id - 1
                  this%iter_id = iter_id_save - 1
                  IF (this%nb.EQ.0 .OR. this%max_id.EQ.0) this%min_id = 0
               ENDIF
            ELSE
               !swap with the last non-empty element of the collection
               IF (this%max_id.GT.this%min_id) THEN
                  this%vec(del_id)%t => this%vec(this%max_id)%t
                  this%vec(this%max_id)%t => NULL()
               ELSE
                  this%vec(del_id)%t => NULL()
               ENDIF

               this%nb = this%nb - 1
               this%max_id = this%max_id - 1
               this%iter_id = iter_id_save - 1
               IF (this%nb.EQ.0 .OR. this%max_id.EQ.0) this%min_id = 0
            ENDIF

            info=0
            RETURN
        END SUBROUTINE ppm_rc_c_list_remove
        SUBROUTINE ppm_rc_c_list_remove_(this,info,element,id)
            !!! If element is present, remove it from the collection
            !!! else, remove the current element (as defined by the iterator pointer)

            IMPLICIT NONE
            !-------------------------------------------------------------------------
            !  Arguments
            !-------------------------------------------------------------------------
            CLASS(ppm_rc_c_list_),         INTENT(INOUT) :: this

            INTEGER,                       INTENT(  OUT) :: info

            TYPE(ppm_rc_list_), OPTIONAL, INTENT(INOUT) :: element
            INTEGER,            OPTIONAL, INTENT(IN   ) :: id

            !-------------------------------------------------------------------------
            ! local variables
            !-------------------------------------------------------------------------
            INTEGER :: del_id
            INTEGER :: iter_id_save

            iter_id_save = this%iter_id

            IF (PRESENT(element)) THEN
               del_id = MERGE(id,this%get_id(element),PRESENT(id))
               IF (del_id.LT.0) RETURN
            ELSE
               del_id = MERGE(id,this%iter_id,PRESENT(id))
            ENDIF

            !deallocate the element
            CALL this%vec(del_id)%t%destroy()

            IF (PRESENT(id)) THEN
               IF (del_id+1.LE.this%max_id) THEN
                  this%vec(del_id)%t => this%vec(del_id+1)%t
               ELSE
                  this%vec(del_id)%t => NULL()
                  this%nb = this%nb - 1
                  this%max_id = this%max_id - 1
                  this%iter_id = iter_id_save - 1
                  IF (this%nb.EQ.0 .OR. this%max_id.EQ.0) this%min_id = 0
               ENDIF
            ELSE
               !swap with the last non-empty element of the collection
               IF (this%max_id.GT.this%min_id) THEN
                  this%vec(del_id)%t => this%vec(this%max_id)%t
                  this%vec(this%max_id)%t => NULL()
               ELSE
                  this%vec(del_id)%t => NULL()
               ENDIF

               this%nb = this%nb - 1
               this%max_id = this%max_id - 1
               this%iter_id = iter_id_save - 1
               IF (this%nb.EQ.0 .OR. this%max_id.EQ.0) this%min_id = 0
            ENDIF

            info=0

            RETURN
        END SUBROUTINE ppm_rc_c_list_remove_
        SUBROUTINE ppm_rc_c_list_2_remove(this,info,element,id)
            !!! If element is present, remove it from the collection
            !!! else, remove the current element (as defined by the iterator pointer)
            IMPLICIT NONE
            !-------------------------------------------------------------------------
            !  Arguments
            !-------------------------------------------------------------------------
            CLASS(ppm_rc_c_list_2),        INTENT(INOUT) :: this
            INTEGER,                       INTENT(  OUT) :: info
            TYPE(ppm_rc_list_2), OPTIONAL, INTENT(INOUT) :: element
            INTEGER,             OPTIONAL, INTENT(IN   ) :: id
            !-------------------------------------------------------------------------
            ! local variables
            !-------------------------------------------------------------------------
            INTEGER :: del_id
            INTEGER :: iter_id_save
            iter_id_save = this%iter_id
            IF (PRESENT(element)) THEN
               del_id = MERGE(id,this%get_id(element),PRESENT(id))
               IF (del_id.LT.0) RETURN
            ELSE
               del_id = MERGE(id,this%iter_id,PRESENT(id))
            ENDIF
            !deallocate the element
            CALL this%vec(del_id)%t%destroy()
            IF (PRESENT(id)) THEN
               IF (del_id+1.LE.this%max_id) THEN
                  this%vec(del_id)%t => this%vec(del_id+1)%t
               ELSE
                  this%vec(del_id)%t => NULL()
                  this%nb = this%nb - 1
                  this%max_id = this%max_id - 1
                  this%iter_id = iter_id_save - 1
                  IF (this%nb.EQ.0 .OR. this%max_id.EQ.0) this%min_id = 0
               ENDIF
            ELSE
               !swap with the last non-empty element of the collection
               IF (this%max_id.GT.this%min_id) THEN
                  this%vec(del_id)%t => this%vec(this%max_id)%t
                  this%vec(this%max_id)%t => NULL()
               ELSE
                  this%vec(del_id)%t => NULL()
               ENDIF
               this%nb = this%nb - 1
               this%max_id = this%max_id - 1
               this%iter_id = iter_id_save - 1
               IF (this%nb.EQ.0 .OR. this%max_id.EQ.0) this%min_id = 0
            ENDIF
            info=0
            RETURN
        END SUBROUTINE ppm_rc_c_list_2_remove

        !GROW COLLECTION SIZE
        SUBROUTINE ppm_rc_c_list_grow_size(this,info)
            !!! Reallocate collection to a larger size
            !!! (twice the size seems good)

            IMPLICIT NONE
            !-------------------------------------------------------------------------
            !  Arguments
            !-------------------------------------------------------------------------
            CLASS(ppm_rc_c_list), INTENT(INOUT) :: this

            INTEGER,              INTENT(  OUT) :: info

            !-------------------------------------------------------------------------
            !  Local variables
            !-------------------------------------------------------------------------
            TYPE(ppm_rc_t_ptr_list), DIMENSION(:), POINTER :: vec_temp

            INTEGER :: i

            start_subroutine("ppm_rc_c_list_grow_size")

            IF (this%vec_size.LE.0 .OR. .NOT.ASSOCIATED(this%vec)) THEN
                !if the array is empty, allocate with a reasonable size
                this%vec_size = 10
                ALLOCATE(this%vec(this%vec_size),STAT=info)
                or_fail_alloc("could not allocate collection array")
            ELSE
                !if the array is full, double its size
                ! allocate a temporary array to store the current element
                ALLOCATE(vec_temp(this%vec_size),STAT=info)
                or_fail_alloc("could not allocate temporary array")

                ! copy the elements to the temporary array
                DO i=1,this%vec_size
                   vec_temp(i)%t => this%vec(i)%t
                ENDDO

                !reallocate the collection to a larger size
                DEALLOCATE(this%vec,STAT=info)
                or_fail_dealloc("could not deallocate collection array")

                ALLOCATE(this%vec(2*this%vec_size),STAT=info)
                or_fail_alloc("could not allocate collection array")

                !copy the elements back from the temporary array to the collection
                DO i=1,this%vec_size
                   this%vec(i)%t => vec_temp(i)%t
                ENDDO
                this%vec_size = 2 * this%vec_size

                DEALLOCATE(vec_temp,STAT=info)
                or_fail_dealloc("could not deallocate temporary array")
                NULLIFY(vec_temp)
            ENDIF

            end_subroutine()
        END SUBROUTINE ppm_rc_c_list_grow_size
        SUBROUTINE ppm_rc_c_list_grow_size_(this,info)
            !!! Reallocate collection to a larger size
            !!! (twice the size seems good)

            IMPLICIT NONE
            !-------------------------------------------------------------------------
            !  Arguments
            !-------------------------------------------------------------------------
            CLASS(ppm_rc_c_list_), INTENT(INOUT) :: this

            INTEGER,               INTENT(  OUT) :: info

            !-------------------------------------------------------------------------
            !  Local variables
            !-------------------------------------------------------------------------
            TYPE(ppm_rc_t_ptr_list_), DIMENSION(:), POINTER :: vec_temp

            INTEGER :: i

            start_subroutine("ppm_rc_c_list_grow_size")

            IF (this%vec_size.LE.0 .OR. .NOT.ASSOCIATED(this%vec)) THEN
                !if the array is empty, allocate with a reasonable size
                this%vec_size = 10
                ALLOCATE(this%vec(this%vec_size),STAT=info)
                or_fail_alloc("could not allocate collection array")
            ELSE
                !if the array is full, double its size
                ! allocate a temporary array to store the current element
                ALLOCATE(vec_temp(this%vec_size),STAT=info)
                or_fail_alloc("could not allocate temporary array")

                ! copy the elements to the temporary array
                DO i=1,this%vec_size
                   vec_temp(i)%t => this%vec(i)%t
                ENDDO

                !reallocate the collection to a larger size
                DEALLOCATE(this%vec,STAT=info)
                or_fail_dealloc("could not deallocate collection array")

                ALLOCATE(this%vec(2*this%vec_size),STAT=info)
                or_fail_alloc("could not allocate collection array")

                !copy the elements back from the temporary array to the collection
                DO i=1,this%vec_size
                   this%vec(i)%t => vec_temp(i)%t
                ENDDO
                this%vec_size = 2 * this%vec_size

                DEALLOCATE(vec_temp,STAT=info)
                or_fail_dealloc("could not deallocate temporary array")
                NULLIFY(vec_temp)
            ENDIF

            end_subroutine()
        END SUBROUTINE ppm_rc_c_list_grow_size_
        SUBROUTINE ppm_rc_c_list_2_grow_size(this,info)
            !!! Reallocate collection to a larger size
            !!! (twice the size seems good)
            IMPLICIT NONE
            !-------------------------------------------------------------------------
            !  Arguments
            !-------------------------------------------------------------------------
            CLASS(ppm_rc_c_list_2), INTENT(INOUT) :: this
            INTEGER,                INTENT(  OUT) :: info
            !-------------------------------------------------------------------------
            !  Local variables
            !-------------------------------------------------------------------------
            TYPE(ppm_rc_t_ptr_list_2), DIMENSION(:), POINTER :: vec_temp
            INTEGER :: i
            start_subroutine("ppm_rc_c_list_grow_size")
            IF (this%vec_size.LE.0.OR..NOT.ASSOCIATED(this%vec)) THEN
                !if the array is empty, allocate with a reasonable size
                this%vec_size = 10
                ALLOCATE(this%vec(this%vec_size),STAT=info)
                or_fail_alloc("could not allocate collection array")
            ELSE
                !if the array is full, double its size
                ! allocate a temporary array to store the current element
                ALLOCATE(vec_temp(this%vec_size),STAT=info)
                or_fail_alloc("could not allocate temporary array")
                ! copy the elements to the temporary array
                DO i=1,this%vec_size
                   vec_temp(i)%t => this%vec(i)%t
                ENDDO
                !reallocate the collection to a larger size
                DEALLOCATE(this%vec,STAT=info)
                or_fail_dealloc("could not deallocate collection array")
                ALLOCATE(this%vec(2*this%vec_size),STAT=info)
                or_fail_alloc("could not allocate collection array")
                !copy the elements back from the temporary array to the collection
                DO i=1,this%vec_size
                   this%vec(i)%t => vec_temp(i)%t
                ENDDO
                this%vec_size = 2 * this%vec_size
                DEALLOCATE(vec_temp,STAT=info)
                or_fail_dealloc("could not deallocate temporary array")
                NULLIFY(vec_temp)
            ENDIF
            end_subroutine()
        END SUBROUTINE ppm_rc_c_list_2_grow_size

        !AT
        FUNCTION ppm_rc_c_list_at(this,i) RESULT (element)
            !!! Access the i-th element of the collection
            !!! (the ordering is the same as that given by the iterators next()
            !!! begin(), prev() and last() )
            !!! If there is less than i elements in the collection, the function
            !!! returns a NULL pointer.

            IMPLICIT NONE
            !-------------------------------------------------------------------------
            !  Arguments
            !-------------------------------------------------------------------------
            CLASS(ppm_rc_c_list), INTENT(INOUT) :: this

            INTEGER,              INTENT(IN   ) :: i

            TYPE(ppm_rc_list),    POINTER       :: element

            IF (i.GE.this%min_id.AND.i.LE.this%max_id.AND.i.NE.0) THEN
               element => this%vec(i)%t
               this%iter_id=i
               DO WHILE (.NOT.ASSOCIATED(element))
                  IF (this%iter_id.GT.this%max_id) THEN
                     element => NULL()
                     RETURN
                  ENDIF
                  element => this%next()
               ENDDO
            ELSE
               element => NULL()
            ENDIF
            RETURN
        END FUNCTION ppm_rc_c_list_at
        FUNCTION ppm_rc_c_list_at_(this,i) RESULT (element)
            !!! Access the i-th element of the collection
            !!! (the ordering is the same as that given by the iterators next()
            !!! begin(), prev() and last() )
            !!! If there is less than i elements in the collection, the function
            !!! returns a NULL pointer.

            IMPLICIT NONE
            !-------------------------------------------------------------------------
            !  Arguments
            !-------------------------------------------------------------------------
            CLASS(ppm_rc_c_list_), INTENT(INOUT) :: this

            INTEGER,               INTENT(IN   ) :: i

            TYPE(ppm_rc_list_),    POINTER       :: element

            IF (i.GE.this%min_id.AND.i.LE.this%max_id.AND.i.NE.0) THEN
               element => this%vec(i)%t
               this%iter_id=i
               DO WHILE (.NOT.ASSOCIATED(element))
                  IF (this%iter_id.GT.this%max_id) THEN
                     element => NULL()
                     RETURN
                  ENDIF
                  element => this%next()
               ENDDO
            ELSE
               element => NULL()
            ENDIF
            RETURN
        END FUNCTION ppm_rc_c_list_at_
        FUNCTION ppm_rc_c_list_2_at(this,i) RESULT (element)
            !!! Access the i-th element of the collection
            !!! (the ordering is the same as that given by the iterators next()
            !!! begin(), prev() and last() )
            !!! If there is less than i elements in the collection, the function
            !!! returns a NULL pointer.
            IMPLICIT NONE
            !-------------------------------------------------------------------------
            !  Arguments
            !-------------------------------------------------------------------------
            CLASS(ppm_rc_c_list_2), INTENT(INOUT) :: this
            INTEGER,                INTENT(IN   ) :: i
            TYPE(ppm_rc_list_2),    POINTER       :: element
            IF (i.GE.this%min_id.AND.i.LE.this%max_id.AND.i.NE.0) THEN
               element => this%vec(i)%t
               this%iter_id=i
               DO WHILE (.NOT.ASSOCIATED(element))
                  IF (this%iter_id.GT.this%max_id) THEN
                     element => NULL()
                     RETURN
                  ENDIF
                  element => this%next()
               ENDDO
            ELSE
               element => NULL()
            ENDIF
            RETURN
        END FUNCTION ppm_rc_c_list_2_at
