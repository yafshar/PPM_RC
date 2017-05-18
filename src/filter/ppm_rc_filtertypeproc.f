#if   __DIME == __3D
        !
        SUBROUTINE ThresholdImageFilter_ThresholdAbove(this,thresh,info,OutsideValue)

          IMPLICIT NONE

          CLASS(ThresholdImageFilter)       :: this

          REAL(MK),           INTENT(IN   ) :: thresh

          INTEGER,            INTENT(  OUT) :: info

          REAL(MK), OPTIONAL, INTENT(IN   ) :: OutsideValue

          start_subroutine("ThresholdImageFilter_ThresholdAbove")
          this%m_Upper=REAL(thresh,ppm_kind_double)
          this%m_Lower=zerod
          this%m_OutsideValue=MERGE(REAL(OutsideValue,ppm_kind_double),zerod,PRESENT(OutsideValue))
          end_subroutine()

        END SUBROUTINE ThresholdImageFilter_ThresholdAbove

        !
        SUBROUTINE ThresholdImageFilter_ThresholdBelow(this,thresh,info,OutsideValue)

          IMPLICIT NONE

          CLASS(ThresholdImageFilter)       :: this

          REAL(MK),           INTENT(IN   ) :: thresh

          INTEGER,            INTENT(  OUT) :: info

          REAL(MK), OPTIONAL, INTENT(IN   ) :: OutsideValue

          start_subroutine("ThresholdImageFilter_ThresholdBelow")
          this%m_Upper=REAL(big,ppm_kind_double)
          this%m_Lower=REAL(thresh,ppm_kind_double)
          this%m_OutsideValue=MERGE(REAL(OutsideValue,ppm_kind_double),zerod,PRESENT(OutsideValue))
          end_subroutine()

        END SUBROUTINE ThresholdImageFilter_ThresholdBelow

        !
        SUBROUTINE ThresholdImageFilter_ThresholdOutside(this,lower,upper,info,OutsideValue)

          IMPLICIT NONE

          CLASS(ThresholdImageFilter)       :: this

          REAL(MK),           INTENT(IN   ) :: lower
          REAL(MK),           INTENT(IN   ) :: upper

          INTEGER,            INTENT(  OUT) :: info

          REAL(MK), OPTIONAL, INTENT(IN   ) :: OutsideValue

          start_subroutine("ThresholdImageFilter_ThresholdOutside")
          check_true(<#lower.LE.upper#>,'Lower threshold cannot be greater than upper threshold.')
          this%m_Upper=REAL(upper,ppm_kind_double)
          this%m_Lower=REAL(lower,ppm_kind_double)
          this%m_OutsideValue=MERGE(REAL(OutsideValue,ppm_kind_double),zerod,PRESENT(OutsideValue))
          end_subroutine()

        END SUBROUTINE ThresholdImageFilter_ThresholdOutside
#endif

        !
        SUBROUTINE DTYPE(ThresholdImageFilter_GenerateData)(this, &
        &          FieldIn,MeshIn,info,FieldOut,MaskValue)

          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(ThresholdImageFilter)                  :: this

          CLASS(ppm_t_field_),           POINTER       :: FieldIn
          !!! the source filed which is supposed to be

          CLASS(ppm_t_equi_mesh_),       POINTER       :: MeshIn
          !!! source mesh on which fields are descritized

          INTEGER,                       INTENT(  OUT) :: info

          CLASS(ppm_t_field_), OPTIONAL, POINTER       :: FieldOut
          !!! the source filed which is supposed to be

          INTEGER,             OPTIONAL, INTENT(IN   ) :: MaskValue

          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          CLASS(ppm_t_subpatch_), POINTER :: sbpitr

#if   __DIME == __2D
          REAL(MK), CONTIGUOUS, DIMENSION(:,:),   POINTER :: DTYPE(wpin_r)
          REAL(MK), CONTIGUOUS, DIMENSION(:,:),   POINTER :: DTYPE(wpout_r)
#elif __DIME == __3D
          REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER :: DTYPE(wpin_r)
          REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER :: DTYPE(wpout_r)
#endif

#if   __DIME == __2D
          INTEGER, CONTIGUOUS, DIMENSION(:,:),   POINTER :: DTYPE(wpin_i)
          INTEGER, CONTIGUOUS, DIMENSION(:,:),   POINTER :: DTYPE(wpout_i)
#elif __DIME == __3D
          INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER :: DTYPE(wpin_i)
          INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER :: DTYPE(wpout_i)
#endif
          INTEGER,             DIMENSION(:),     POINTER :: Nn

          LOGICAL                                        :: m_Threshold
          !To decide whether it is just a simple thresholding or
          !it is a mask

          start_subroutine("ThresholdImageFilter_GenerateData")

          CALL check()

          m_Threshold=MERGE(.FALSE.,.TRUE.,PRESENT(MaskValue))

          sbpitr => MeshIn%subpatch%begin()
          DO WHILE (ASSOCIATED(sbpitr))
             Nn => sbpitr%nnodes

             SELECT CASE (FieldIn%data_type)
             CASE (ppm_type_int)
                CALL sbpitr%get_field(FieldIn,DTYPE(wpin_i),info)
                or_fail("Failed to get field wpin_i data.")

                IF (PRESENT(FieldOut)) THEN
                   CALL sbpitr%get_field(FieldOut,DTYPE(wpout_i),info)
                   or_fail("Failed to get field wpout data.")
                ELSE
                   DTYPE(wpout_i) => DTYPE(wpin_i)
                ENDIF

                SELECT CASE (m_Threshold)
                CASE (.TRUE.)
                   CALL DTYPE(Update_i)(DTYPE(wpin_i),DTYPE(wpout_i),Nn,info)
                CASE DEFAULT
                   CALL DTYPE(Update_i)(DTYPE(wpin_i),DTYPE(wpout_i),Nn,info,MaskValue)
                END SELECT
                or_fail("Update")

             CASE (ppm_type_real,ppm_type_real_single)
                CALL sbpitr%get_field(FieldIn,DTYPE(wpin_r),info)
                or_fail("Failed to get field wpin_i data.")

                IF (PRESENT(FieldOut)) THEN
                   CALL sbpitr%get_field(FieldOut,DTYPE(wpout_r),info)
                   or_fail("Failed to get field wpout data.")
                ELSE
                   DTYPE(wpout_r) => DTYPE(wpin_r)
                ENDIF

                SELECT CASE (m_Threshold)
                CASE (.TRUE.)
                   CALL DTYPE(Update_r)(DTYPE(wpin_r),DTYPE(wpout_r),Nn,info)
                CASE DEFAULT
                   CALL DTYPE(Update_r)(DTYPE(wpin_r),DTYPE(wpout_r),Nn,info,REAL(MaskValue,MK))
                END SELECT
                or_fail("Update")

             END SELECT

             sbpitr => MeshIn%subpatch%next()
          ENDDO !(ASSOCIATED(sbpitr))

          end_subroutine()

          CONTAINS
          SUBROUTINE check
          IMPLICIT NONE
            check_true(<#FieldIn%is_discretized_on(MeshIn)#>, &
            & "The input Field has not been descritized on the mesh!",exit_point=8888)

            IF (PRESENT(FieldOut)) THEN
               check_true(<#FieldOut%is_discretized_on(MeshIn)#>, &
               & "The output Field is present, but it has not been descritized on the mesh!",exit_point=8888)

               check_true(<#FieldOut%data_type.EQ.FieldIn%data_type#>, &
               & "The Input and output fields have different data types!",exit_point=8888)
            ENDIF
          8888 CONTINUE
          END SUBROUTINE check

          SUBROUTINE DTYPE(Update_i)(wpin,wpout,nnodes,info_,m_MaskValue)
            IMPLICIT NONE
#if   __DIME == __2D
            INTEGER, CONTIGUOUS, DIMENSION(:,:),   POINTER       :: wpin
            INTEGER, CONTIGUOUS, DIMENSION(:,:),   POINTER       :: wpout
#elif __DIME == __3D
            INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: wpin
            INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: wpout
#endif

            INTEGER,             DIMENSION(:),     INTENT(IN   ) :: nnodes
            INTEGER,                               INTENT(  OUT) :: info_
            INTEGER,             OPTIONAL,         INTENT(IN   ) :: m_MaskValue

            INTEGER :: lower,upper,OutsideValue_
            INTEGER :: i,j
#if   __DIME == __3D
            INTEGER :: k
#endif

            info_=0

            lower=INT(this%m_Lower)
            IF (this%m_Upper.GE.REAL(bigi-1,MK)) THEN
               upper=bigi-1
            ELSE
               upper=INT(this%m_Upper)
            ENDIF
            OutsideValue_=INT(this%m_OutsideValue)

            IF (ASSOCIATED(wpout,wpin)) THEN
               SELECT CASE (lower)
               CASE (0)
#if   __DIME == __2D
                  DO j=1,nnodes(2)
                     DO i=1,nnodes(1)
                        IF (wpin(i,j).GE.upper) THEN
                           wpout(i,j)=OutsideValue_
                        ENDIF
#elif __DIME == __3D
                  DO k=1,nnodes(3)
                     DO j=1,nnodes(2)
                        DO i=1,nnodes(1)
                           IF (wpin(i,j,k).GE.upper) THEN
                              wpout(i,j,k)=OutsideValue_
                           ENDIF
                        ENDDO
#endif
                     ENDDO
                  ENDDO

               CASE DEFAULT
#if   __DIME == __2D
                  DO j=1,nnodes(2)
                     DO i=1,nnodes(1)
                        IF (wpin(i,j).LE.lower.OR.wpin(i,j).GE.upper) THEN
                           wpout(i,j)=OutsideValue_
                        ENDIF
#elif __DIME == __3D
                  DO k=1,nnodes(3)
                     DO j=1,nnodes(2)
                        DO i=1,nnodes(1)
                           IF (wpin(i,j,k).LE.lower.OR.wpin(i,j,k).GE.upper) THEN
                              wpout(i,j,k)=OutsideValue_
                           ENDIF
                        ENDDO
#endif
                     ENDDO
                  ENDDO

               END SELECT
            ELSE
               SELECT CASE (lower)
               CASE (0)
#if   __DIME == __2D
                  DO CONCURRENT (i=1:nnodes(1),j=1:nnodes(2))
                     IF (wpin(i,j).GE.upper) THEN
                        wpout(i,j)=OutsideValue_
                     ELSE
                        wpout(i,j)=wpin(i,j)
                     ENDIF
#elif __DIME == __3D
                  DO CONCURRENT (i=1:nnodes(1),j=1:nnodes(2),k=1:nnodes(3))
                     IF (wpin(i,j,k).GE.upper) THEN
                        wpout(i,j,k)=OutsideValue_
                     ELSE
                        wpout(i,j,k)=wpin(i,j,k)
                     ENDIF
#endif
                  ENDDO
               CASE DEFAULT
#if   __DIME == __2D
                  DO CONCURRENT (i=1:nnodes(1),j=1:nnodes(2))
                     IF (wpin(i,j).LE.lower.OR.wpin(i,j).GE.upper) THEN
                        wpout(i,j)=OutsideValue_
                     ELSE
                        wpout(i,j)=wpin(i,j)
                     ENDIF
#elif __DIME == __3D
                  DO CONCURRENT (i=1:nnodes(1),j=1:nnodes(2),k=1:nnodes(3))
                     IF (wpin(i,j,k).LE.lower.OR.wpin(i,j,k).GE.upper) THEN
                        wpout(i,j,k)=OutsideValue_
                     ELSE
                        wpout(i,j,k)=wpin(i,j,k)
                     ENDIF
#endif
                  ENDDO

               END SELECT
            ENDIF !(ASSOCIATED(wpout,wpin))

            IF (PRESENT(m_MaskValue)) THEN
#if   __DIME == __2D
               DO j=1,nnodes(2)
                  DO i=1,nnodes(1)
                     IF (wpout(i,j).NE.OutsideValue_) THEN
                        wpout(i,j)=m_MaskValue
                     ENDIF
#elif __DIME == __3D
               DO k=1,nnodes(3)
                  DO j=1,nnodes(2)
                     DO i=1,nnodes(1)
                        IF (wpout(i,j,k).NE.OutsideValue_) THEN
                           wpout(i,j,k)=m_MaskValue
                        ENDIF
                     ENDDO
#endif
                  ENDDO
               ENDDO
            ENDIF

          CONTINUE
          END SUBROUTINE DTYPE(Update_i)

          SUBROUTINE DTYPE(Update_r)(wpin,wpout,nnodes,info_,m_MaskValue)
            IMPLICIT NONE
#if   __DIME == __2D
            REAL(MK), CONTIGUOUS, DIMENSION(:,:),   POINTER       :: wpin
            REAL(MK), CONTIGUOUS, DIMENSION(:,:),   POINTER       :: wpout
#elif __DIME == __3D
            REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: wpin
            REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: wpout
#endif
            INTEGER,              DIMENSION(:),     INTENT(IN   ) :: nnodes
            INTEGER,                                INTENT(  OUT) :: info_

            REAL(MK),             OPTIONAL,         INTENT(IN   ) :: m_MaskValue

            REAL(MK) :: lower,upper,OutsideValue_

            INTEGER :: i,j
#if   __DIME == __3D
            INTEGER :: k
#endif

            info_=0

            lower=REAL(this%m_Lower,MK)
            upper=REAL(this%m_Upper,MK)
            OutsideValue_=REAL(this%m_OutsideValue,MK)

            IF (ASSOCIATED(wpout,wpin)) THEN
               SELECT CASE (INT(lower))
               CASE (0)
#if   __DIME == __2D
                  DO j=1,nnodes(2)
                     DO i=1,nnodes(1)
                        IF (wpin(i,j).GE.upper) THEN
                           wpout(i,j)=OutsideValue_
                        ENDIF
#elif __DIME == __3D
                  DO k=1,nnodes(3)
                     DO j=1,nnodes(2)
                        DO i=1,nnodes(1)
                           IF (wpin(i,j,k).GE.upper) THEN
                              wpout(i,j,k)=OutsideValue_
                           ENDIF
                        ENDDO
#endif
!                   END FORALL
                     ENDDO
                  ENDDO
               CASE DEFAULT
#if   __DIME == __2D
                  DO j=1,nnodes(2)
                     DO i=1,nnodes(1)
                        IF (wpin(i,j).LE.lower.OR.wpin(i,j).GE.upper) THEN
                           wpout(i,j)=OutsideValue_
                        ENDIF
#elif __DIME == __3D
                  DO k=1,nnodes(3)
                     DO j=1,nnodes(2)
                        DO i=1,nnodes(1)
                           IF (wpin(i,j,k).LE.lower.OR.wpin(i,j,k).GE.upper) THEN
                              wpout(i,j,k)=OutsideValue_
                           ENDIF
                        ENDDO
#endif
                     ENDDO
                  ENDDO
               END SELECT
            ELSE
               SELECT CASE (INT(lower))
               CASE (0)
#if   __DIME == __2D
                  DO CONCURRENT (i=1:nnodes(1),j=1:nnodes(2))
                     IF (wpin(i,j).GE.upper) THEN
                        wpout(i,j)=OutsideValue_
                     ELSE
                        wpout(i,j)=wpin(i,j)
                     ENDIF
#elif __DIME == __3D
                  DO CONCURRENT (i=1:nnodes(1),j=1:nnodes(2),k=1:nnodes(3))
                     IF (wpin(i,j,k).GE.upper) THEN
                        wpout(i,j,k)=OutsideValue_
                     ELSE
                        wpout(i,j,k)=wpin(i,j,k)
                     ENDIF
#endif
                  ENDDO
               CASE DEFAULT
#if   __DIME == __2D
                  DO CONCURRENT (i=1:nnodes(1),j=1:nnodes(2))
                     IF (wpin(i,j).LE.lower.OR.wpin(i,j).GE.upper) THEN
                        wpout(i,j)=OutsideValue_
                     ELSE
                        wpout(i,j)=wpin(i,j)
                     ENDIF
#elif __DIME == __3D
                  DO CONCURRENT (i=1:nnodes(1),j=1:nnodes(2),k=1:nnodes(3))
                     IF (wpin(i,j,k).LE.lower.OR.wpin(i,j,k).GE.upper) THEN
                        wpout(i,j,k)=OutsideValue_
                     ELSE
                        wpout(i,j,k)=wpin(i,j,k)
                     ENDIF
#endif
                  ENDDO
               END SELECT
            ENDIF !(ASSOCIATED(wpout,wpin))

            IF (PRESENT(m_MaskValue)) THEN
#if   __DIME == __2D
               DO j=1,nnodes(2)
                  DO i=1,nnodes(1)
                     IF (ABS(wpout(i,j)-OutsideValue_).GT.lmyeps) THEN
                        wpout(i,j)=m_MaskValue
                     ENDIF
#elif __DIME == __3D
               DO k=1,nnodes(3)
                  DO j=1,nnodes(2)
                     DO i=1,nnodes(1)
                        IF (ABS(wpout(i,j,k)-OutsideValue_).GT.lmyeps) THEN
                           wpout(i,j,k)=m_MaskValue
                        ENDIF
                     ENDDO
#endif
                  ENDDO
               ENDDO
            ENDIF !

          CONTINUE
          END SUBROUTINE DTYPE(Update_r)
        END SUBROUTINE DTYPE(ThresholdImageFilter_GenerateData)

#if   __DIME == __3D
        !
        SUBROUTINE Histogram_Initialize(this,hSize,info,lowerBound,upperBound)

          IMPLICIT NONE

          CLASS(Histogram)                  :: this

          INTEGER,            INTENT(IN   ) :: hSize
          INTEGER,            INTENT(  OUT) :: info

          REAL(MK), OPTIONAL, INTENT(IN   ) :: lowerBound
          REAL(MK), OPTIONAL, INTENT(IN   ) :: upperBound

          REAL(MK) :: lowerBound_
          REAL(MK) :: upperBound_

          INTEGER :: i

          start_subroutine("Histogram_Initialize")

          this%histSize=hSize+1
          !!!The upperBound will be in the new cell
          !!!that is why we need one more

          ALLOCATE(this%FrequencyContainer(hSize+1), &
          &        this%OffsetTable(hSize+1),STAT=info)
          or_fail_alloc("FrequencyContainer & this%OffsetTable")

          this%FrequencyContainer=zerod

          IF (PRESENT(lowerBound)) THEN
             lowerBound_=MERGE(lowerBound,zero,lowerBound.GT.zero)
          ELSE
             lowerBound_=zero
          ENDIF

          IF (PRESENT(upperBound)) THEN
             upperBound_=MERGE(upperBound,255._MK,upperBound.GT.zero)
          ELSE
             upperBound_=255._MK
          ENDIF

          this%histOffset=REAL(upperBound_-lowerBound_,ppm_kind_double)/REAL(hSize,ppm_kind_double)

          this%OffsetTable(1)=REAL(lowerBound_,ppm_kind_double)+this%histOffset*halfd

          DO i=2,hSize+1
             this%OffsetTable(i)=this%OffsetTable(i-1)+this%histOffset
          ENDDO

          end_subroutine()

        END SUBROUTINE Histogram_Initialize

        SUBROUTINE Histogram_destroy(this,info)

          IMPLICIT NONE

          CLASS(Histogram)       :: this

          INTEGER, INTENT(  OUT) :: info

          start_subroutine("Histogram_destroy")

          IF (ALLOCATED(this%FrequencyContainer)) THEN
             DEALLOCATE(this%FrequencyContainer,STAT=info)
             or_fail_dealloc("FrequencyContainer")
          ENDIF
          IF (ALLOCATED(this%OffsetTable)) THEN
             DEALLOCATE(this%OffsetTable,STAT=info)
             or_fail_dealloc("OffsetTable")
          ENDIF

          end_subroutine()

        END SUBROUTINE Histogram_destroy
#endif

        SUBROUTINE DTYPE(Histogram_create)(this,FieldIn,MeshIn,info)

          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(Histogram)                       :: this

          CLASS(ppm_t_field_),     POINTER       :: FieldIn
          !!! the source filed which is supposed to be

          CLASS(ppm_t_equi_mesh_), POINTER       :: MeshIn
          !!! source mesh on which fields are descritized

          INTEGER,                 INTENT(  OUT) :: info

          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          CLASS(ppm_t_subpatch_), POINTER :: sbpitr

#if   __DIME == __2D
          REAL(MK), CONTIGUOUS, DIMENSION(:,:),   POINTER :: DTYPE(wp_r)
#elif __DIME == __3D
          REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER :: DTYPE(wp_r)
#endif

#if   __DIME == __2D
          INTEGER, CONTIGUOUS, DIMENSION(:,:),   POINTER :: DTYPE(wp_i)
#elif __DIME == __3D
          INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER :: DTYPE(wp_i)
#endif
          INTEGER,             DIMENSION(:),     POINTER :: Nn

          start_subroutine("Histogram_create")

          CALL check()

          sbpitr => MeshIn%subpatch%begin()
          DO WHILE (ASSOCIATED(sbpitr))
             Nn => sbpitr%nnodes

             SELECT CASE (FieldIn%data_type)
             CASE (ppm_type_int)
                CALL sbpitr%get_field(FieldIn,DTYPE(wp_i),info)
                or_fail("Failed to get field wpin_i data.")

                CALL DTYPE(Update_i)(DTYPE(wp_i),Nn,info)
                or_fail("Update")

             CASE (ppm_type_real,ppm_type_real_single)
                CALL sbpitr%get_field(FieldIn,DTYPE(wp_r),info)
                or_fail("Failed to get field wpin_i data.")

                CALL DTYPE(Update_r)(DTYPE(wp_r),Nn,info)
                or_fail("Update")

             END SELECT

             sbpitr => MeshIn%subpatch%next()
          ENDDO !(ASSOCIATED(sbpitr))

          end_subroutine()

          CONTAINS

          SUBROUTINE check
          IMPLICIT NONE
            check_true(<#FieldIn%is_discretized_on(MeshIn)#>, &
            & "The input Field has not been descritized on the mesh!",exit_point=8888)

            check_true(<#ALLOCATED(this%FrequencyContainer).AND.ALLOCATED(this%OffsetTable)#>, &
            & "The histogram is not initialized yet!",exit_point=8888)

          8888 CONTINUE
          END SUBROUTINE check

          SUBROUTINE DTYPE(Update_i)(wpin,nnodes,info)
            IMPLICIT NONE
#if   __DIME == __2D
            INTEGER, CONTIGUOUS, DIMENSION(:,:),   POINTER       :: wpin
#elif __DIME == __3D
            INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: wpin
#endif
            INTEGER,             DIMENSION(:),     INTENT(IN   ) :: nnodes
            INTEGER,                               INTENT(  OUT) :: info

            INTEGER :: i,j,l
#if   __DIME == __3D
            INTEGER :: k
#endif

            info=0

            ALLOCATE(tmp1_i(1:nnodes(1)),STAT=info)
            or_fail_alloc("tmp1_i",exit_point=7777)

#if   __DIME == __2D
            DO j=1,nnodes(2)
               FORALL (i=1:nnodes(1))
                  tmp1_i(i)=INT(REAL(wpin(i,j),ppm_kind_double)/this%histOffset)+1
               END FORALL

               DO i=1,nnodes(1)
                  l=tmp1_i(i)
                  this%FrequencyContainer(l)= this%FrequencyContainer(l)+oned
               ENDDO
            ENDDO
#elif __DIME == __3D
            DO k=1,nnodes(3)
               DO j=1,nnodes(2)
                  FORALL (i=1:nnodes(1))
                     tmp1_i(i)=INT(REAL(wpin(i,j,k),ppm_kind_double)/this%histOffset)+1
                  END FORALL

                  DO i=1,nnodes(1)
                     l=tmp1_i(i)
                     this%FrequencyContainer(l)=this%FrequencyContainer(l)+oned
                  ENDDO
               ENDDO
            ENDDO
#endif

            DEALLOCATE(tmp1_i,STAT=info)
            or_fail_dealloc("tmp1_i",exit_point=7777)

          7777 CONTINUE
          END SUBROUTINE DTYPE(Update_i)

          SUBROUTINE DTYPE(Update_r)(wpin,nnodes,info)
            IMPLICIT NONE
#if   __DIME == __2D
            REAL(MK), CONTIGUOUS, DIMENSION(:,:),   POINTER       :: wpin
#elif __DIME == __3D
            REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: wpin
#endif
            INTEGER,              DIMENSION(:),     INTENT(IN   ) :: nnodes
            INTEGER,                                INTENT(  OUT) :: info

            INTEGER :: i,j,l
#if   __DIME == __3D
            INTEGER :: k
#endif

            info=0

            ALLOCATE(tmp1_i(1:nnodes(1)),STAT=info)
            or_fail_alloc("tmp1_i",exit_point=7777)

#if   __DIME == __2D
            DO j=1,nnodes(2)
               FORALL (i=1:nnodes(1))
                  tmp1_i(i)=INT(REAL(wpin(i,j),ppm_kind_double)/this%histOffset)+1
               END FORALL

               DO i=1,nnodes(1)
                  l=tmp1_i(i)
                  this%FrequencyContainer(l)= this%FrequencyContainer(l)+oned
               ENDDO
            ENDDO
#elif __DIME == __3D
            DO k=1,nnodes(3)
               DO j=1,nnodes(2)
                  FORALL (i=1:nnodes(1))
                     tmp1_i(i)=INT(REAL(wpin(i,j,k),ppm_kind_double)/this%histOffset)+1
                  END FORALL

                  DO i=1,nnodes(1)
                     l=tmp1_i(i)
                     this%FrequencyContainer(l)=this%FrequencyContainer(l)+oned
                  ENDDO
               ENDDO
            ENDDO
#endif

            DEALLOCATE(tmp1_i,STAT=info)
            or_fail_dealloc("tmp1_i",exit_point=7777)

          7777 CONTINUE
          END SUBROUTINE DTYPE(Update_r)
        END SUBROUTINE DTYPE(Histogram_create)


#if   __DIME == __3D
        !
        SUBROUTINE Histogram_Mean(this,mean,info,totalFrequency)

          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          USE ppm_module_mpi
          IMPLICIT NONE

          CLASS(Histogram)                               :: this

          REAL(ppm_kind_double),           INTENT(  OUT) :: mean

          INTEGER,                         INTENT(  OUT) :: info

          REAL(ppm_kind_double), OPTIONAL, INTENT(  OUT) :: totalFrequency

          REAL(ppm_kind_double), DIMENSION(2) :: sum_

          start_subroutine("Histogram_Mean")

#ifdef __MPI
          CALL MPI_Allreduce(MPI_IN_PLACE,this%FrequencyContainer, &
          &    this%histSize,MPI_DOUBLE_PRECISION,MPI_SUM,comm,info)
          or_fail_MPI("MPI_Allreduce")
#endif

          sum_(1)=SUM(this%FrequencyContainer)
          sum_(2)=DOT_PRODUCT(this%OffsetTable,this%FrequencyContainer)

          mean=sum_(2)/sum_(1)

          IF (PRESENT(totalFrequency)) totalFrequency=sum_(1)

          end_subroutine()

        END SUBROUTINE Histogram_Mean

        !
        SUBROUTINE Histogram_Quantile(this,p,Quantile,info)

          IMPLICIT NONE

          CLASS(Histogram)                     :: this

          REAL(MK),              INTENT(IN   ) :: p
          REAL(ppm_kind_double), INTENT(  OUT) :: Quantile

          INTEGER,               INTENT(  OUT) :: info

          REAL(ppm_kind_double) :: p_n_prev,p_n,f_n,p_
          REAL(ppm_kind_double) :: cumulated
          REAL(ppm_kind_double) :: totalFrequency
          REAL(ppm_kind_double) :: mean
          REAL(ppm_kind_double) :: binProportion

          INTEGER :: m,n

          start_subroutine("Histogram_Quantile")

          CALL this%Mean(mean,info,totalFrequency)
          or_fail("this%Mean")

          cumulated=zerod

          p_=REAL(p,ppm_kind_double)

          IF (p.LT.half) THEN
             n  =1
             p_n=zerod
             DO WHILE (n.LE.this%histSize.AND.p_n.LT.p_)
                f_n = this%FrequencyContainer(n)
                cumulated=cumulated+f_n
                p_n_prev=p_n
                p_n = cumulated/totalFrequency
                n=n+1
             ENDDO

             binProportion = f_n / totalFrequency

             Quantile=this%OffsetTable(n-1)-this%histOffset*halfd+((p_-p_n_prev)/binProportion)*this%histOffset
          ELSE
             n  =this%histSize
             m  =1
             p_n=oned
             DO WHILE (m.LT.this%histSize.AND.p_n.GT.p_)
                f_n = this%FrequencyContainer(n)
                cumulated=cumulated+f_n
                p_n_prev=p_n
                p_n = oned - cumulated/totalFrequency
                n=n-1
                m=m+1
             ENDDO

             binProportion = f_n / totalFrequency

             Quantile=this%OffsetTable(n+1)-this%histOffset*halfd+((p_n_prev-p_)/binProportion)*this%histOffset
          ENDIF

          end_subroutine()

        END SUBROUTINE Histogram_Quantile
#endif

        !
        SUBROUTINE DTYPE(OtsuThresholdCalculator_SetInput)(this, &
        &          FieldIn,MeshIn,hSize,lowerBound,upperBound,info)

          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(OtsuThresholdCalculator)         :: this

          CLASS(ppm_t_field_),     POINTER       :: FieldIn
          !!! the source filed which is supposed to be

          CLASS(ppm_t_equi_mesh_), POINTER       :: MeshIn
          !!! source mesh on which fields are descritized

          INTEGER,                 INTENT(IN   ) :: hSize
          !!! the histogram size
          INTEGER,                 INTENT(IN   ) :: lowerBound
          INTEGER,                 INTENT(IN   ) :: upperBound
          INTEGER,                 INTENT(  OUT) :: info

          start_subroutine("OtsuThresholdCalculator_SetInput")

          CALL this%Initialize(hSize,info,REAL(lowerBound,MK),REAL(upperBound,MK))
          or_fail("Initialize")

          CALL this%DTYPE(create)(FieldIn,MeshIn,info)
          or_fail("create")

          end_subroutine()

        END SUBROUTINE DTYPE(OtsuThresholdCalculator_SetInput)


        SUBROUTINE DTYPE(OtsuThresholdImageFilter_GenerateData)(this, &
        &          FieldIn,MeshIn,hSize,info,FieldOut,lowerBound,     &
        &          upperBound,MaskValue)

          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(OtsuThresholdImageFilter)              :: this

          CLASS(ppm_t_field_),           POINTER       :: FieldIn
          !!! the source filed which is supposed to be

          CLASS(ppm_t_equi_mesh_),       POINTER       :: MeshIn
          !!! source mesh on which fields are descritized

          INTEGER,                       INTENT(IN   ) :: hSize
          !!!Histogram size(number of Bins)

          INTEGER,                       INTENT(  OUT) :: info

          CLASS(ppm_t_field_), OPTIONAL, POINTER       :: FieldOut
          !!! the source filed which is supposed to be

          INTEGER,             OPTIONAL, INTENT(IN   ) :: lowerBound
          !!!Histogram lowerBound
          INTEGER,             OPTIONAL, INTENT(IN   ) :: upperBound
          !!!Histogram upperBound
          INTEGER,             OPTIONAL, INTENT(IN   ) :: MaskValue

          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          TYPE(OtsuThresholdCalculator) :: otsuhistogram

          REAL(ppm_kind_double)            :: mean,totalFrequency
          REAL(ppm_kind_double)            :: totalMean
          REAL(ppm_kind_double), PARAMETER :: tol = 0.00001_ppm_kind_double
          REAL(ppm_kind_double)            :: freqLeft
          REAL(ppm_kind_double)            :: meanLeft
          REAL(ppm_kind_double)            :: meanRight
          REAL(ppm_kind_double)            :: maxVarBetween
          REAL(ppm_kind_double)            :: varBetween
          REAL(ppm_kind_double)            :: freqLeftOld
          REAL(ppm_kind_double)            :: meanLeftOld

          INTEGER :: i
          INTEGER :: lowerBound_,upperBound_
          INTEGER :: maxBinNumber

          LOGICAL :: m_Threshold
          !logical value to decide whether only doing thresholding
          !or create a Mask based on threshold

          start_subroutine("OtsuThresholdImageFilter_GenerateData")

          CALL check()

          IF (PRESENT(lowerBound)) THEN
             lowerBound_=MERGE(lowerBound,0,lowerBound.GT.0)
          ELSE
             lowerBound_=0
          ENDIF

          IF (PRESENT(upperBound)) THEN
             upperBound_=MERGE(upperBound,255,upperBound.GT.0)
          ELSE
             upperBound_=255
          ENDIF

          m_Threshold=MERGE(.FALSE.,.TRUE.,PRESENT(MaskValue))

          CALL otsuhistogram%DTYPE(SetInput)(FieldIn,MeshIn,hSize,lowerBound_,upperBound_,info)
          or_fail("otsuhistogram%SetInput")

          CALL otsuhistogram%Mean(mean,info,totalFrequency)
          or_fail("otsuhistogram%Mean")

          check_true(<#totalFrequency.GT.zerod#>,"Histogram is empty")

          ALLOCATE(tmp1_rd(otsuhistogram%histSize),STAT=info)
          or_fail_alloc("tmp1_rd")

          tmp1_rd=otsuhistogram%FrequencyContainer/totalFrequency
          !compute the probability distribution

          totalMean=zerod
          !totalMean level of the original picture

          DO i=1,otsuhistogram%histSize
             totalMean=totalMean+REAL(i,ppm_kind_double)*tmp1_rd(i)
             !muT=sum_1_L(i*wi)
          ENDDO

          !compute Otsu's threshold by maximizing the between-class
          !variance
          freqLeft=tmp1_rd(1)
          !w0
          meanLeft=oned
          !mu0

          IF (freqLeft.LT.oned) THEN
             meanRight=(totalMean-meanLeft*freqLeft)/(oned-freqLeft)
             !mu1=(muT-mu0*w0)/(1-w0)
          ELSE
             meanRight=zerod
          ENDIF

          maxVarBetween=freqLeft*(oned-freqLeft)*(meanLeft-meanRight)**2
          !maxVarBetween=w0(1-w0)(mu0-mu1)^2

          maxBinNumber=1

          freqLeftOld=freqLeft
          meanLeftOld=meanLeft

          DO i=2,otsuhistogram%histSize
             freqLeft=freqLeft+tmp1_rd(i)
             !w0
             IF (freqLeft.GT.zerod) THEN
                meanLeft=(meanLeftOld*freqLeftOld+REAL(i,ppm_kind_double)*tmp1_rd(i))/freqLeft
                !mu0=(sum_1_i(i*w0i))/(sum_1_i(w0))
             ENDIF

             IF (freqLeft.GE.oned) THEN
                meanRight=zerod
                !mu1=0
             ELSE
                meanRight=(totalMean-meanLeft*freqLeft)/(oned-freqLeft)
                !mu1=(muT-mu0*w0)/(1-w0)
             ENDIF

             varBetween=freqLeft*(oned-freqLeft)*(meanLeft-meanRight)**2
             !varBetween=w0*(1-w0)(mu0-mu1)^2

             IF ((varBetween-tol).GT.maxVarBetween) THEN
                maxVarBetween=varBetween
                maxBinNumber =i
             ENDIF

             freqLeftOld=freqLeft
             meanLeftOld=meanLeft
             !cache old values
          ENDDO

          DEALLOCATE(tmp1_rd,STAT=info)
          or_fail_dealloc("tmp1_rd")

          ! should be this for backward compatibility
          ! this->GetOutput()->Set( static_cast<OutputType>(
          ! histogram->GetBinMin( 0, maxBinNumber + 1 ) ) );

          CALL this%Calculator%ThresholdBelow(REAL(otsuhistogram%OffsetTable(maxBinNumber),MK),info)
          or_fail("this%Calculator%ThresholdBelow")

          CALL otsuhistogram%destroy(info)
          or_fail("otsuhistogram%destroy")

          SELECT CASE (m_Threshold)
          CASE (.TRUE.)
          !We only will do threshold
             SELECT CASE (PRESENT(FieldOut))
             CASE (.TRUE.)
                CALL this%Calculator%DTYPE(GenerateData)(FieldIn,MeshIn,info,FieldOut)
             CASE (.FALSE.)
                CALL this%Calculator%DTYPE(GenerateData)(FieldIn,MeshIn,info)
             END SELECT
          CASE DEFAULT
          !create a MASK based on thresholding
             SELECT CASE (PRESENT(FieldOut))
             CASE (.TRUE.)
                CALL this%Calculator%DTYPE(GenerateData)(FieldIn, &
                &    MeshIn,info,FieldOut,MaskValue=MaskValue)
             CASE (.FALSE.)
                CALL this%Calculator%DTYPE(GenerateData)(FieldIn, &
                &    MeshIn,info,MaskValue=MaskValue)
             END SELECT
          END SELECT
          or_fail("this%Calculator%GenerateData")

          end_subroutine()
        CONTAINS
          SUBROUTINE check
          IMPLICIT NONE
            check_true(<#hSize.GT.0#>,"Error otsu thresholding requires Histogram size (number of Bins) greater than zero!",exit_point=8888)

            check_true(<#((PRESENT(lowerBound).AND.PRESENT(upperBound)).OR.(.NOT.(PRESENT(lowerBound).AND.PRESENT(upperBound))))#>, &
            & "Either both lowerBound & upperBound should be available or non!",exit_point=8888)

            check_true(<#FieldIn%is_discretized_on(MeshIn)#>, &
            & "The input Field has not been descritized on the mesh!",exit_point=8888)

            IF (PRESENT(FieldOut)) THEN
               check_true(<#FieldOut%is_discretized_on(MeshIn)#>, &
               & "The output Field is present, but it has not been descritized on the mesh!",exit_point=8888)

               check_true(<#FieldOut%data_type.EQ.FieldIn%data_type#>, &
               & "The Input and output fields have different data types!",exit_point=8888)
            ENDIF
          8888 CONTINUE
          END SUBROUTINE check
        END SUBROUTINE DTYPE(OtsuThresholdImageFilter_GenerateData)


        SUBROUTINE DTYPE(OtsuMultipleThresholdsImageFilter_GenerateData)(this, &
        &          FieldIn,MeshIn,hSize,info,FieldOut,lowerBound,upperBound,   &
        &          MaskValue,NumberOfThresholds)

          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(OtsuMultipleThresholdsImageFilter)     :: this

          CLASS(ppm_t_field_),           POINTER       :: FieldIn
          !!! the source filed which is supposed to be

          CLASS(ppm_t_equi_mesh_),       POINTER       :: MeshIn
          !!! source mesh on which fields are descritized

          INTEGER,                       INTENT(IN   ) :: hSize
          !!!Histogram size(number of Bins)

          INTEGER,                       INTENT(  OUT) :: info

          CLASS(ppm_t_field_), OPTIONAL, POINTER       :: FieldOut
          !!! the source filed which is supposed to be

          INTEGER,             OPTIONAL, INTENT(IN   ) :: lowerBound
          !!!Histogram lowerBound
          INTEGER,             OPTIONAL, INTENT(IN   ) :: upperBound
          !!!Histogram upperBound
          INTEGER,             OPTIONAL, INTENT(IN   ) :: MaskValue
          INTEGER,             OPTIONAL, INTENT(IN   ) :: NumberOfThresholds

          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          TYPE(OtsuThresholdCalculator) :: otsuhistogram

          REAL(ppm_kind_double)                            :: globalMean,globalFrequencySum
          REAL(ppm_kind_double)                            :: freqSum
          REAL(ppm_kind_double), DIMENSION(:), ALLOCATABLE :: classFrequency
          REAL(ppm_kind_double), DIMENSION(:), ALLOCATABLE :: classMean
          REAL(ppm_kind_double)                            :: meanSum
          REAL(ppm_kind_double)                            :: maxVarBetween
          REAL(ppm_kind_double)                            :: varBetween

          INTEGER                            :: i
          INTEGER                            :: numberOfClasses
          INTEGER, DIMENSION(:), ALLOCATABLE :: thresholdIndexes
          INTEGER, DIMENSION(:), ALLOCATABLE :: maxVarThresholdIndexes
          INTEGER                            :: lowerBound_,upperBound_,NumberOfThresholds_
          !INTEGER :: maxBinNumber

          LOGICAL :: m_Threshold
          !logical value to decide whether only doing thresholding
          !or create a Mask based on threshold

          start_subroutine("OtsuMultipleThresholdsImageFilter_GenerateData")

          CALL check()

          IF (PRESENT(lowerBound)) THEN
             lowerBound_=MERGE(lowerBound,0,lowerBound.GT.0)
          ELSE
             lowerBound_=0
          ENDIF

          IF (PRESENT(upperBound)) THEN
             upperBound_=MERGE(upperBound,255,upperBound.GT.0)
          ELSE
             upperBound_=255
          ENDIF

          m_Threshold=MERGE(.FALSE.,.TRUE.,PRESENT(MaskValue))

          IF (PRESENT(NumberOfThresholds)) THEN
             NumberOfThresholds_=MERGE(NumberOfThresholds,1,NumberOfThresholds.GT.0)
          ELSE
             NumberOfThresholds_=1
          ENDIF

          CALL otsuhistogram%DTYPE(SetInput)(FieldIn,MeshIn,hSize, &
          &    lowerBound_,upperBound_,info)
          or_fail("otsuhistogram%SetInput")

          CALL otsuhistogram%Mean(globalMean,info,globalFrequencySum)
          or_fail("otsuhistogram%Mean")

          check_true(<#globalFrequencySum.GT.zerod#>,"Histogram is empty")

          ALLOCATE(thresholdIndexes(NumberOfThresholds_),STAT=info)
          or_fail_alloc("thresholdIndexes")

          !initialize thresholds
          FORALL (i=1:NumberOfThresholds_) thresholdIndexes(i)=i

          ALLOCATE(maxVarThresholdIndexes(NumberOfThresholds_),SOURCE=thresholdIndexes,STAT=info)
          or_fail_alloc("maxVarThresholdIndexes")

          !compute frequency and mean of initial classes
          numberOfClasses = NumberOfThresholds_ + 1

          ALLOCATE(classFrequency(numberOfClasses),STAT=info)
          or_fail_alloc("classFrequency")

          freqSum = zerod

          DO i=1,numberOfClasses-1
             classFrequency(i)=otsuhistogram%FrequencyContainer(thresholdIndexes(i))
             freqSum = freqSum + classFrequency(i)
          ENDDO

          classFrequency(numberOfClasses)=globalFrequencySum - freqSum

          ALLOCATE(classMean(numberOfClasses),STAT=info)
          or_fail_alloc("classMean")

          meanSum=zerod

          DO i=1,numberOfClasses-1
             IF (classFrequency(i).GT.zerod) THEN
                classMean(i) = otsuhistogram%OffsetTable(i)
             ELSE
                classMean(i) = zerod
             ENDIF
             meanSum = meanSum + classMean(i) * classFrequency(i)
          ENDDO

          IF (classFrequency(numberOfClasses).GT.zerod) THEN
             classMean(numberOfClasses) = (globalMean*globalFrequencySum-meanSum)/classFrequency(numberOfClasses)
          ELSE
             classMean(numberOfClasses) = zerod
          ENDIF

          maxVarBetween=zerod
          DO i=1,numberOfClasses
             maxVarBetween=maxVarBetween+classFrequency(i)*((globalMean-classMean(i))*(globalMean-classMean(i)))
          ENDDO

          !explore all possible threshold configurations and choose the one that
          !yields maximum between-class variance
          DO WHILE (IncrementThresholds(thresholdIndexes,classMean,classFrequency))
             varBetween=zerod
             DO i=1,numberOfClasses
                varBetween = varBetween+classFrequency(i)*((globalMean-classMean(i))*(globalMean-classMean(i)))
             ENDDO

             IF (varBetween.GT.maxVarBetween) THEN
                maxVarBetween = varBetween
                maxVarThresholdIndexes = thresholdIndexes
             ENDIF
          ENDDO

          ! should be this for backward compatibility
          ! this->GetOutput()->Set( static_cast<OutputType>(
          ! histogram->GetBinMin( 0, maxBinNumber + 1 ) ) );
          CALL this%Calculator%ThresholdBelow(REAL(otsuhistogram%OffsetTable(maxVarThresholdIndexes(1)),MK),info)
          or_fail("this%Calculator%ThresholdBelow")

          DEALLOCATE(thresholdIndexes,maxVarThresholdIndexes,classFrequency,classMean,STAT=info)
          or_fail_dealloc("thresholdIndexes,maxVarThresholdIndexes,classFrequency & classMean")

          CALL otsuhistogram%destroy(info)
          or_fail("otsuhistogram%destroy")

          SELECT CASE (m_Threshold)
          CASE (.TRUE.)
          !We only will do threshold
             SELECT CASE (PRESENT(FieldOut))
             CASE (.TRUE.)
                CALL this%Calculator%DTYPE(GenerateData)(FieldIn,MeshIn,info,FieldOut)
             CASE (.FALSE.)
                CALL this%Calculator%DTYPE(GenerateData)(FieldIn,MeshIn,info)
             END SELECT
          CASE DEFAULT
          !create a MASK based on thresholding
             SELECT CASE (PRESENT(FieldOut))
             CASE (.TRUE.)
                CALL this%Calculator%DTYPE(GenerateData)(FieldIn, &
                &    MeshIn,info,FieldOut,MaskValue=MaskValue)
             CASE (.FALSE.)
                CALL this%Calculator%DTYPE(GenerateData)(FieldIn, &
                &    MeshIn,info,MaskValue=MaskValue)
             END SELECT
          END SELECT
          or_fail("this%Calculator%GenerateData")

          end_subroutine()
        CONTAINS
          SUBROUTINE check
          IMPLICIT NONE
            check_true(<#hSize.GT.0#>,"Error otsu thresholding requires Histogram size (number of Bins) greater than zero!",exit_point=8888)

            check_true(<#((PRESENT(lowerBound).AND.PRESENT(upperBound)).OR.(.NOT.(PRESENT(lowerBound).AND.PRESENT(upperBound))))#>, &
            & "Either both lowerBound & upperBound should be available or non!",exit_point=8888)

            check_true(<#FieldIn%is_discretized_on(MeshIn)#>, &
            & "The input Field has not been descritized on the mesh!",exit_point=8888)

            IF (PRESENT(FieldOut)) THEN
               check_true(<#FieldOut%is_discretized_on(MeshIn)#>, &
               & "The output Field is present, but it has not been descritized on the mesh!",exit_point=8888)

               check_true(<#FieldOut%data_type.EQ.FieldIn%data_type#>, &
               & "The Input and output fields have different data types!",exit_point=8888)
            ENDIF
          8888 CONTINUE
          END SUBROUTINE check

          FUNCTION IncrementThresholds(thresholdIndexes_,classMean_,classFrequency_)
          IMPLICIT NONE

          INTEGER,               DIMENSION(:), INTENT(INOUT) :: thresholdIndexes_

          REAL(ppm_kind_double), DIMENSION(:), INTENT(INOUT) :: classMean_
          REAL(ppm_kind_double), DIMENSION(:), INTENT(INOUT) :: classFrequency_

          LOGICAL                                            :: IncrementThresholds

          REAL(ppm_kind_double) :: meanOld,freqOld

          INTEGER :: j
          INTEGER :: numberOfHistogramBins,numberOfClasses

          numberOfHistogramBins = otsuhistogram%histSize
          numberOfClasses = SIZE(classMean_)

          !from the upper threshold down
          DO i=NumberOfThresholds_,1,-1
             !if this threshold can be incremented (i.e.
             !we're not at the end of the histogram)
             IF (thresholdIndexes_(i).LT.numberOfHistogramBins-1-(NumberOfThresholds_-i)) THEN
                !increment it and update mean and frequency of the class bounded by the
                !threshold
                thresholdIndexes_(i)=thresholdIndexes_(i)+1

                meanOld = classMean_(i)
                freqOld = classFrequency_(i)

                classFrequency_(i)=classFrequency_(i)+otsuhistogram%FrequencyContainer(thresholdIndexes_(i))

                IF (classFrequency_(i).GT.zerod) THEN
                   classMean_(i)=(meanOld*freqOld+otsuhistogram%OffsetTable(thresholdIndexes_(i))*otsuhistogram%FrequencyContainer(thresholdIndexes_(i)))/classFrequency_(i)
                ELSE
                   classMean_(i)=zerod
                ENDIF

                !set higher thresholds adjacent to their previous ones, and update mean
                !and frequency of the respective classes

                DO j=i+1,NumberOfThresholds_
                   thresholdIndexes_(j)=thresholdIndexes_(j-1)+1
                   classFrequency_(j)=otsuhistogram%FrequencyContainer(thresholdIndexes_(j))

                   IF (classFrequency_(j).GT.zerod) THEN
                      classMean_(j)=otsuhistogram%OffsetTable(thresholdIndexes_(j))
                   ELSE
                      classMean_(j)=zerod
                   ENDIF
                ENDDO

                !update mean and frequency of the highest class
                classFrequency_(numberOfClasses)=globalFrequencySum
                classMean_(numberOfClasses)=globalMean*globalFrequencySum

                DO j=1,numberOfClasses
                   classFrequency_(numberOfClasses)=classFrequency_(numberOfClasses)-classFrequency_(j)
                   classMean_(numberOfClasses)=classMean_(numberOfClasses)-classMean_(j)*classFrequency_(j)
                ENDDO

                IF (classFrequency_(numberOfClasses).GT.zerod) THEN
                   classMean_(numberOfClasses)=classMean_(numberOfClasses)/classFrequency_(numberOfClasses)
                ELSE
                   classMean_(numberOfClasses)=zerod
                ENDIF

                !exit the loop if a threshold has been incremented
                EXIT
             !if this threshold can't be incremented
             ELSE
                !if it's the lowest threshold
                IF (i.EQ.0) THEN
                   !we couldn't increment because we're done
                   IncrementThresholds=.FALSE.
                   RETURN
                ENDIF
             ENDIF
          ENDDO !i=NumberOfThresholds_,1,-1

          !we incremented
          IncrementThresholds=.TRUE.
          !7777 CONTINUE
          RETURN
          END FUNCTION IncrementThresholds

        END SUBROUTINE DTYPE(OtsuMultipleThresholdsImageFilter_GenerateData)












        !
        SUBROUTINE DTYPE(FilterDataArray)(this,wpin,wpout,nnodes,direction,info)
        !!! Apply the Recursive Filter to an array of data.
        !!! This routine does not work for in-place data buffer
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(RecursiveSeparableImageFilter)                  :: this
#if   __DIME == __2D
          REAL(MK), CONTIGUOUS, DIMENSION(:,:),   POINTER       :: wpin
          !!! the input source filed
          REAL(MK), CONTIGUOUS, DIMENSION(:,:),   POINTER       :: wpout
          !!! the output source filed
#elif __DIME == __3D
          REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: wpin
          !!! the input source filed
          REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: wpout
          !!! the output source filed
#endif
          INTEGER,              DIMENSION(:),     INTENT(IN   ) :: nnodes
          !!! Number of mesh points in every direction
          INTEGER,                                INTENT(IN   ) :: direction
          !!! direction in which the filter is to be applied X,Y,Z ~ 1,2,3
          INTEGER,                                INTENT(  OUT) :: info
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          REAL(MK)              :: outV
          REAL(ppm_kind_double) ::  t0

          INTEGER               :: i,j,l,m,n,o,p,Nm1,Nm2
#if   __DIME == __3D
          INTEGER               :: k,Nm3
#endif
          INTEGER, DIMENSION(5) :: ll

          CHARACTER(LEN=ppm_char) :: caller = 'FilterDataArray'

          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          IF (ASSOCIATED(wpin,wpout)) THEN
             fail("Recursive Filter failed to work with inplace data buffer!", &
             & ppm_error=ppm_error_fatal)
          ENDIF

          Nm1=nnodes(1)
          Nm2=nnodes(2)

#if   __DIME == __2D
          SELECT CASE (direction)
          CASE (1)
          !!! X direction
             ALLOCATE(tmp1_r(Nm1),STAT=info)
             or_fail_alloc("Failed to allocate tmp1_r!")

             DO j=1,Nm2
                ! this value is assumed to exist from the border to infinity.
                outV = wpin(1,j)

                ! Initialize borders
                tmp1_r(1) =      outV*this%N0 +      outV*this%N1 +      outV*this%N2 + outV*this%N3
                tmp1_r(2) = wpin(2,j)*this%N0 +      outV*this%N1 +      outV*this%N2 + outV*this%N3
                tmp1_r(3) = wpin(3,j)*this%N0 + wpin(2,j)*this%N1 +      outV*this%N2 + outV*this%N3
                tmp1_r(4) = wpin(4,j)*this%N0 + wpin(3,j)*this%N1 + wpin(2,j)*this%N2 + outV*this%N3

                ! note that the outV value is multiplied by the Boundary coefficients BNi
                tmp1_r(1) = tmp1_r(1) - (    outV*this%BN1 +     outV*this%BN2 +     outV*this%BN3 + outV*this%BN4)
                tmp1_r(2) = tmp1_r(2) - (tmp1_r(1)*this%D1 +     outV*this%BN2 +     outV*this%BN3 + outV*this%BN4)
                tmp1_r(3) = tmp1_r(3) - (tmp1_r(2)*this%D1 + tmp1_r(1)*this%D2 +     outV*this%BN3 + outV*this%BN4)
                tmp1_r(4) = tmp1_r(4) - (tmp1_r(3)*this%D1 + tmp1_r(2)*this%D2 + tmp1_r(1)*this%D3 + outV*this%BN4)

                !Recursively filter the rest
                DO i=5,Nm1
                   tmp1_r(i) = wpin(i,j)*this%N0 + wpin(i-1,j)*this%N1 + wpin(i-2,j)*this%N2 + wpin(i-3,j)*this%N3
                   tmp1_r(i) = tmp1_r(i) - (tmp1_r(i-1)*this%D1 + tmp1_r(i-2)*this%D2 + tmp1_r(i-3)*this%D3 + tmp1_r(i-4)*this%D4)
                ENDDO !i=1,Nm1

                FORALL (i=1:Nm1)
                   wpout(i,j)=tmp1_r(i)
                END FORALL

                ! AntiCausal direction pass
                ! this value is assumed to exist from the border to infinity.
                outV = wpin(Nm1,j)

                ! Initialize borders
                tmp1_r(Nm1  ) =          outV*this%M1 +          outV*this%M2 +        outV*this%M3 + outV*this%M4
                tmp1_r(Nm1-1) = wpin(Nm1  ,j)*this%M1 +          outV*this%M2 +        outV*this%M3 + outV*this%M4
                tmp1_r(Nm1-2) = wpin(Nm1-1,j)*this%M1 + wpin(Nm1  ,j)*this%M2 +        outV*this%M3 + outV*this%M4
                tmp1_r(Nm1-3) = wpin(Nm1-2,j)*this%M1 + wpin(Nm1-1,j)*this%M2 + wpin(Nm1,j)*this%M3 + outV*this%M4

                ! note that the outV value is multiplied by the Boundary coefficients BMi
                tmp1_r(Nm1  ) = tmp1_r(Nm1  ) - (        outV*this%BM1 +         outV*this%BM2 +       outV*this%BM3 + outV*this%BM4)
                tmp1_r(Nm1-1) = tmp1_r(Nm1-1) - (tmp1_r(Nm1  )*this%D1 +         outV*this%BM2 +       outV*this%BM3 + outV*this%BM4)
                tmp1_r(Nm1-2) = tmp1_r(Nm1-2) - (tmp1_r(Nm1-1)*this%D1 + tmp1_r(Nm1  )*this%D2 +       outV*this%BM3 + outV*this%BM4)
                tmp1_r(Nm1-3) = tmp1_r(Nm1-3) - (tmp1_r(Nm1-2)*this%D1 + tmp1_r(Nm1-1)*this%D2 + tmp1_r(Nm1)*this%D3 + outV*this%BM4)

                ! Recursively filter the rest
                DO i=Nm1-3,2,-1
                   tmp1_r(i-1) = wpin(i,j)*this%M1 + wpin(i+1,j)*this%M2 + wpin(i+2,j)*this%M3 + wpin(i+3,j)*this%M4
                   tmp1_r(i-1) = tmp1_r(i-1) - (tmp1_r(i)*this%D1 + tmp1_r(i+1)*this%D2 + tmp1_r(i+2)*this%D3 + tmp1_r(i+3)*this%D4)
                ENDDO

                ! Roll the antiCausal part into the output
                FORALL (i=1:Nm1)
                   wpout(i,j)=wpout(i,j)+tmp1_r(i)
                END FORALL
             ENDDO !j=ldl(2),ldu(2)

             DEALLOCATE(tmp1_r,STAT=info)
             or_fail_dealloc("Failed to deallocate tmp1_r!")
          CASE (2)
          !!! Y direction

             !At each pixel, I should keep four values of the first, second
             !third and the fourth row to avoid cash missses
             ALLOCATE(tmp2_r(Nm1,5),STAT=info)
             or_fail_alloc("tmp2_r")

             ! This optimisation is done to avoid recomputing what has been computed once
             FORALL (i=1:Nm1)
                ! Initialize borders
                tmp2_r(i,1) = wpin(i,1)*this%N0 + wpin(i,1)*this%N1 + wpin(i,1)*this%N2 + wpin(i,1)*this%N3
                tmp2_r(i,2) = wpin(i,2)*this%N0 + wpin(i,1)*this%N1 + wpin(i,1)*this%N2 + wpin(i,1)*this%N3
                tmp2_r(i,3) = wpin(i,3)*this%N0 + wpin(i,2)*this%N1 + wpin(i,1)*this%N2 + wpin(i,1)*this%N3
                tmp2_r(i,4) = wpin(i,4)*this%N0 + wpin(i,3)*this%N1 + wpin(i,2)*this%N2 + wpin(i,1)*this%N3

                ! note that the wpin(i,1) value is multiplied by the Boundary coefficients BNi
                tmp2_r(i,1) = tmp2_r(i,1) - ( wpin(i,1)*this%BN1 +  wpin(i,1)*this%BN2 +  wpin(i,1)*this%BN3 + wpin(i,1)*this%BN4)
                tmp2_r(i,2) = tmp2_r(i,2) - (tmp2_r(i,1)*this%D1 +  wpin(i,1)*this%BN2 +  wpin(i,1)*this%BN3 + wpin(i,1)*this%BN4)
                tmp2_r(i,3) = tmp2_r(i,3) - (tmp2_r(i,2)*this%D1 + tmp2_r(i,1)*this%D2 +  wpin(i,1)*this%BN3 + wpin(i,1)*this%BN4)
                tmp2_r(i,4) = tmp2_r(i,4) - (tmp2_r(i,3)*this%D1 + tmp2_r(i,2)*this%D2 + tmp2_r(i,1)*this%D3 + wpin(i,1)*this%BN4)
             END FORALL

             FORALL(i=1:5) ll(i)=i

             !Recursively filter the rest
             DO j=5,Nm2
                l=ll(5)
                m=ll(4)
                n=ll(3)
                o=ll(2)
                p=ll(1)
                FORALL (i=1:Nm1)
                   tmp2_r(i,l)=wpin(i,j)*this%N0 + wpin(i,j-1)*this%N1 + wpin(i,j-2)*this%N2 + wpin(i,j-3)*this%N3
                   tmp2_r(i,l)=tmp2_r(i,l) - (tmp2_r(i,m)*this%D1 + tmp2_r(i,n)*this%D2 + tmp2_r(i,o)*this%D3 + tmp2_r(i,p)*this%D4)
                END FORALL

                !Store the causal result
                l=j-4
                FORALL (i=1:Nm1)
                   wpout(i,l)=tmp2_r(i,p)
                END FORALL
                ll=CSHIFT(ll,SHIFT=1,DIM=1)
             ENDDO !j=5,Nm2

             !Store the causal result for the last 4 lines
             p=ll(1)
             j=Nm2-3
             FORALL (i=1:Nm1)
                wpout(i,j)=tmp2_r(i,p)
             END FORALL
             o=ll(2)
             j=Nm2-2
             FORALL (i=1:Nm1)
                wpout(i,j)=tmp2_r(i,o)
             END FORALL
             n=ll(3)
             j=Nm2-1
             FORALL (i=1:Nm1)
                wpout(i,j)=tmp2_r(i,n)
             END FORALL
             m=ll(4)
             j=Nm2
             FORALL (i=1:Nm1)
                wpout(i,j)=tmp2_r(i,m)
             END FORALL


             !AntiCausal direction pass
             FORALL (i=1:Nm1)
                ! Initialize borders
                tmp2_r(i,1) = wpin(i,Nm2  )*this%M1 + wpin(i,Nm2  )*this%M2 + wpin(i,Nm2)*this%M3 + wpin(i,Nm2)*this%M4
                tmp2_r(i,2) = wpin(i,Nm2  )*this%M1 + wpin(i,Nm2  )*this%M2 + wpin(i,Nm2)*this%M3 + wpin(i,Nm2)*this%M4
                tmp2_r(i,3) = wpin(i,Nm2-1)*this%M1 + wpin(i,Nm2  )*this%M2 + wpin(i,Nm2)*this%M3 + wpin(i,Nm2)*this%M4
                tmp2_r(i,4) = wpin(i,Nm2-2)*this%M1 + wpin(i,Nm2-1)*this%M2 + wpin(i,Nm2)*this%M3 + wpin(i,Nm2)*this%M4

                ! note that the outV value is multiplied by the Boundary coefficients BMi
                tmp2_r(i,1) = tmp2_r(i,1) - (wpin(i,Nm2)*this%BM1 + wpin(i,Nm2)*this%BM2 + wpin(i,Nm2)*this%BM3 + wpin(i,Nm2)*this%BM4)
                tmp2_r(i,2) = tmp2_r(i,2) - (tmp2_r(i, 1)*this%D1 + wpin(i,Nm2)*this%BM2 + wpin(i,Nm2)*this%BM3 + wpin(i,Nm2)*this%BM4)
                tmp2_r(i,3) = tmp2_r(i,3) - (tmp2_r(i, 2)*this%D1 + tmp2_r(i, 1)*this%D2 + wpin(i,Nm2)*this%BM3 + wpin(i,Nm2)*this%BM4)
                tmp2_r(i,4) = tmp2_r(i,4) - (tmp2_r(i, 3)*this%D1 + tmp2_r(i, 2)*this%D2 + tmp2_r(i, 1)*this%D3 + wpin(i,Nm2)*this%BM4)
             END FORALL

             FORALL(i=1:5) ll(i)=i

             !Recursively filter the rest
             DO j=Nm2-3,2,-1
                l=ll(5)
                m=ll(4)
                n=ll(3)
                o=ll(2)
                p=ll(1)
                FORALL (i=1:Nm1)
                   tmp2_r(i,l)=wpin(i,j)*this%M1 + wpin(i,j+1)*this%M2 + wpin(i,j+2)*this%M3 + wpin(i,j+3)*this%M4
                   tmp2_r(i,l)=tmp2_r(i,l) - (tmp2_r(i,m)*this%D1 + tmp2_r(i,n)*this%D2 + tmp2_r(i,o)*this%D3 + tmp2_r(i,p)*this%D4)
                END FORALL

                l=j+3
                !Store the causal result
                FORALL (i=1:Nm1)
                   wpout(i,l)=wpout(i,l)+tmp2_r(i,p)
                END FORALL
                ll=CSHIFT(ll,SHIFT=1,DIM=1)
             ENDDO !j=Nm2-3,2,-1

             !Store the causal result for the last 4 lines
             p=ll(1)
             j=4
             FORALL (i=1:Nm1)
                wpout(i,j)=tmp2_r(i,p)
             END FORALL
             o=ll(2)
             j=3
             FORALL (i=1:Nm1)
                wpout(i,j)=tmp2_r(i,o)
             END FORALL
             n=ll(3)
             j=2
             FORALL (i=1:Nm1)
                wpout(i,j)=tmp2_r(i,n)
             END FORALL
             m=ll(4)
             j=1
             FORALL (i=1:Nm1)
                wpout(i,j)=tmp2_r(i,m)
             END FORALL

             DEALLOCATE(tmp2_r,STAT=info)
             or_fail_dealloc("Failed to deallocate tmp2_r!")
          CASE DEFAULT
             fail("Wrong direction for 2D problem!!!")
          END SELECT
#elif __DIME == __3D
          Nm3=nnodes(3)

          SELECT CASE (direction)
          CASE (1)
          !!! X direction
             IF (Nm1.LT.4) THEN
                fail("The number of pixels along dimension is less than 4!!! This filter requires a minimum of four pixels along the dimension to be processed!", &
                & ppm_error=ppm_error_fatal)
             ENDIF

             ALLOCATE(tmp1_r(Nm1),STAT=info)
             or_fail_alloc("Failed to allocate tmp1_r!")

             DO k=1,Nm3
                DO j=1,Nm2
                   ! this value is assumed to exist from the border to infinity.
                   outV = wpin(1,j,k)

                   ! Initialize borders
                   tmp1_r(1) =        outV*this%N0 +        outV*this%N1 +        outV*this%N2 + outV*this%N3
                   tmp1_r(2) = wpin(2,j,k)*this%N0 +        outV*this%N1 +        outV*this%N2 + outV*this%N3
                   tmp1_r(3) = wpin(3,j,k)*this%N0 + wpin(2,j,k)*this%N1 +        outV*this%N2 + outV*this%N3
                   tmp1_r(4) = wpin(4,j,k)*this%N0 + wpin(3,j,k)*this%N1 + wpin(2,j,k)*this%N2 + outV*this%N3

                   ! note that the outV value is multiplied by the Boundary coefficients BNi
                   tmp1_r(1) = tmp1_r(1) - (    outV*this%BN1 +     outV*this%BN2 +     outV*this%BN3 + outV*this%BN4)
                   tmp1_r(2) = tmp1_r(2) - (tmp1_r(1)*this%D1 +     outV*this%BN2 +     outV*this%BN3 + outV*this%BN4)
                   tmp1_r(3) = tmp1_r(3) - (tmp1_r(2)*this%D1 + tmp1_r(1)*this%D2 +     outV*this%BN3 + outV*this%BN4)
                   tmp1_r(4) = tmp1_r(4) - (tmp1_r(3)*this%D1 + tmp1_r(2)*this%D2 + tmp1_r(1)*this%D3 + outV*this%BN4)

                   !Recursively filter the rest
                   DO i=5,Nm1
                      tmp1_r(i) = wpin(i,j,k)*this%N0 + wpin(i-1,j,k)*this%N1 + wpin(i-2,j,k)*this%N2 + wpin(i-3,j,k)*this%N3
                      tmp1_r(i) = tmp1_r(i) - (tmp1_r(i-1)*this%D1 + tmp1_r(i-2)*this%D2 + tmp1_r(i-3)*this%D3 + tmp1_r(i-4)*this%D4)
                   ENDDO !i=1,Nm1

                   FORALL (i=1:Nm1)
                      wpout(i,j,k)=tmp1_r(i)
                   END FORALL

                   ! AntiCausal direction pass
                   ! this value is assumed to exist from the border to infinity.
                   outV = wpin(Nm1,j,k)

                   ! Initialize borders
                   tmp1_r(Nm1  ) =            outV*this%M1 +            outV*this%M2 +          outV*this%M3 + outV*this%M4
                   tmp1_r(Nm1-1) = wpin(Nm1  ,j,k)*this%M1 +            outV*this%M2 +          outV*this%M3 + outV*this%M4
                   tmp1_r(Nm1-2) = wpin(Nm1-1,j,k)*this%M1 + wpin(Nm1  ,j,k)*this%M2 +          outV*this%M3 + outV*this%M4
                   tmp1_r(Nm1-3) = wpin(Nm1-2,j,k)*this%M1 + wpin(Nm1-1,j,k)*this%M2 + wpin(Nm1,j,k)*this%M3 + outV*this%M4

                   ! note that the outV value is multiplied by the Boundary coefficients BMi
                   tmp1_r(Nm1  ) = tmp1_r(Nm1  ) - (        outV*this%BM1 +         outV*this%BM2 +       outV*this%BM3 + outV*this%BM4)
                   tmp1_r(Nm1-1) = tmp1_r(Nm1-1) - (tmp1_r(Nm1  )*this%D1 +         outV*this%BM2 +       outV*this%BM3 + outV*this%BM4)
                   tmp1_r(Nm1-2) = tmp1_r(Nm1-2) - (tmp1_r(Nm1-1)*this%D1 + tmp1_r(Nm1  )*this%D2 +       outV*this%BM3 + outV*this%BM4)
                   tmp1_r(Nm1-3) = tmp1_r(Nm1-3) - (tmp1_r(Nm1-2)*this%D1 + tmp1_r(Nm1-1)*this%D2 + tmp1_r(Nm1)*this%D3 + outV*this%BM4)

                   ! Recursively filter the rest
                   DO i=Nm1-3,2,-1
                      tmp1_r(i-1) = wpin(i,j,k)*this%M1 + wpin(i+1,j,k)*this%M2 + wpin(i+2,j,k)*this%M3 + wpin(i+3,j,k)*this%M4
                      tmp1_r(i-1) = tmp1_r(i-1) - (tmp1_r(i)*this%D1 + tmp1_r(i+1)*this%D2 + tmp1_r(i+2)*this%D3 + tmp1_r(i+3)*this%D4)
                   ENDDO

                   ! Roll the antiCausal part into the output
                   FORALL (i=1:Nm1)
                      wpout(i,j,k)=wpout(i,j,k)+tmp1_r(i)
                   END FORALL
                ENDDO !j=1,Nm2
             ENDDO !k=1,Nm3

             DEALLOCATE(tmp1_r,STAT=info)
             or_fail_dealloc("Failed to deallocate tmp1_r!")
          CASE (2)
          !!! Y direction
             IF (Nm2.LT.4) THEN
                fail("The number of pixels along dimension is less than 4!!! This filter requires a minimum of four pixels along the dimension to be processed!", &
                & ppm_error=ppm_error_fatal)
             ENDIF

             !At each pixel, I should keep three values of the first, second
             !and third row
             ALLOCATE(tmp2_r(Nm1,5),STAT=info)
             or_fail_alloc("tmp2_r")

             DO K=1,NM3
                ! This optimisation is done to avoid recomputing what has been computed once
                FORALL (i=1:Nm1)
                   ! Initialize borders
                   tmp2_r(i,1) = wpin(i,1,k)*this%N0 + wpin(i,1,k)*this%N1 + wpin(i,1,k)*this%N2 + wpin(i,1,k)*this%N3
                   tmp2_r(i,2) = wpin(i,2,k)*this%N0 + wpin(i,1,k)*this%N1 + wpin(i,1,k)*this%N2 + wpin(i,1,k)*this%N3
                   tmp2_r(i,3) = wpin(i,3,k)*this%N0 + wpin(i,2,k)*this%N1 + wpin(i,1,k)*this%N2 + wpin(i,1,k)*this%N3
                   tmp2_r(i,4) = wpin(i,4,k)*this%N0 + wpin(i,3,k)*this%N1 + wpin(i,2,k)*this%N2 + wpin(i,1,k)*this%N3

                   ! note that the wpin(i,1) value is multiplied by the Boundary coefficients BNi
                   tmp2_r(i,1) = tmp2_r(i,1) - (wpin(i,1,k)*this%BN1 + wpin(i,1,k)*this%BN2 +  wpin(i,1,k)*this%BN3 + wpin(i,1,k)*this%BN4)
                   tmp2_r(i,2) = tmp2_r(i,2) - (tmp2_r(i,1)*this%D1  + wpin(i,1,k)*this%BN2 +  wpin(i,1,k)*this%BN3 + wpin(i,1,k)*this%BN4)
                   tmp2_r(i,3) = tmp2_r(i,3) - (tmp2_r(i,2)*this%D1  + tmp2_r(i,1)*this%D2  +  wpin(i,1,k)*this%BN3 + wpin(i,1,k)*this%BN4)
                   tmp2_r(i,4) = tmp2_r(i,4) - (tmp2_r(i,3)*this%D1  + tmp2_r(i,2)*this%D2  +  tmp2_r(i,1)*this%D3  + wpin(i,1,k)*this%BN4)
                END FORALL

                FORALL(i=1:5) ll(i)=i

                !Recursively filter the rest
                DO j=5,Nm2
                   l=ll(5)
                   m=ll(4)
                   n=ll(3)
                   o=ll(2)
                   p=ll(1)
                   FORALL (i=1:Nm1)
                      tmp2_r(i,l)=wpin(i,j,k)*this%N0 + wpin(i,j-1,k)*this%N1 + wpin(i,j-2,k)*this%N2 + wpin(i,j-3,k)*this%N3
                      tmp2_r(i,l)=tmp2_r(i,l) - (tmp2_r(i,m)*this%D1 + tmp2_r(i,n)*this%D2 + tmp2_r(i,o)*this%D3 + tmp2_r(i,p)*this%D4)
                   END FORALL

                   !Store the causal result
                   l=j-4
                   FORALL (i=1:Nm1)
                      wpout(i,l,k)=tmp2_r(i,p)
                   END FORALL
                   ll=CSHIFT(ll,SHIFT=1,DIM=1)
                ENDDO !j=5,Nm2

                !Store the causal result for the last 4 lines
                p=ll(1)
                j=Nm2-3
                FORALL (i=1:Nm1)
                   wpout(i,j,k)=tmp2_r(i,p)
                END FORALL
                o=ll(2)
                j=Nm2-2
                FORALL (i=1:Nm1)
                   wpout(i,j,k)=tmp2_r(i,o)
                END FORALL
                n=ll(3)
                j=Nm2-1
                FORALL (i=1:Nm1)
                   wpout(i,j,k)=tmp2_r(i,n)
                END FORALL
                m=ll(4)
                j=Nm2
                FORALL (i=1:Nm1)
                   wpout(i,j,k)=tmp2_r(i,m)
                END FORALL

                !AntiCausal direction pass
                FORALL (i=1:Nm1)
                   ! Initialize borders
                   tmp2_r(i,1) = wpin(i,Nm2  ,k)*this%M1 + wpin(i,Nm2  ,k)*this%M2 + wpin(i,Nm2,k)*this%M3 + wpin(i,Nm2,k)*this%M4
                   tmp2_r(i,2) = wpin(i,Nm2  ,k)*this%M1 + wpin(i,Nm2  ,k)*this%M2 + wpin(i,Nm2,k)*this%M3 + wpin(i,Nm2,k)*this%M4
                   tmp2_r(i,3) = wpin(i,Nm2-1,k)*this%M1 + wpin(i,Nm2  ,k)*this%M2 + wpin(i,Nm2,k)*this%M3 + wpin(i,Nm2,k)*this%M4
                   tmp2_r(i,4) = wpin(i,Nm2-2,k)*this%M1 + wpin(i,Nm2-1,k)*this%M2 + wpin(i,Nm2,k)*this%M3 + wpin(i,Nm2,k)*this%M4

                   ! note that the outV value is multiplied by the Boundary coefficients BMi
                   tmp2_r(i,1) = tmp2_r(i,1) - (wpin(i,Nm2,k)*this%BM1 + wpin(i,Nm2,k)*this%BM2 + wpin(i,Nm2,k)*this%BM3 + wpin(i,Nm2,k)*this%BM4)
                   tmp2_r(i,2) = tmp2_r(i,2) - (tmp2_r(i,  1)*this%D1  + wpin(i,Nm2,k)*this%BM2 + wpin(i,Nm2,k)*this%BM3 + wpin(i,Nm2,k)*this%BM4)
                   tmp2_r(i,3) = tmp2_r(i,3) - (tmp2_r(i,  2)*this%D1  + tmp2_r(i,  1)*this%D2  + wpin(i,Nm2,k)*this%BM3 + wpin(i,Nm2,k)*this%BM4)
                   tmp2_r(i,4) = tmp2_r(i,4) - (tmp2_r(i,  3)*this%D1  + tmp2_r(i,  2)*this%D2  + tmp2_r(i,  1)*this%D3  + wpin(i,Nm2,k)*this%BM4)
                END FORALL

                FORALL(i=1:5) ll(i)=i

                !Recursively filter the rest
                DO j=Nm2-3,2,-1
                   l=ll(5)
                   m=ll(4)
                   n=ll(3)
                   o=ll(2)
                   p=ll(1)
                   FORALL (i=1:Nm1)
                      tmp2_r(i,l)=wpin(i,j,k)*this%M1 + wpin(i,j+1,k)*this%M2 + wpin(i,j+2,k)*this%M3 + wpin(i,j+3,k)*this%M4
                      tmp2_r(i,l)=tmp2_r(i,l) - (tmp2_r(i,m)*this%D1 + tmp2_r(i,n)*this%D2 + tmp2_r(i,o)*this%D3 + tmp2_r(i,p)*this%D4)
                   END FORALL

                   l=j+3
                   !Store the causal result
                   FORALL (i=1:Nm1)
                      wpout(i,l,k)=wpout(i,l,k)+tmp2_r(i,p)
                   END FORALL
                   ll=CSHIFT(ll,SHIFT=1,DIM=1)
                ENDDO !j=Nm2-3,2,-1

                !Store the causal result for the last 4 lines
                m=ll(4)
                j=1
                FORALL (i=1:Nm1)
                   wpout(i,j,k)=tmp2_r(i,m)
                END FORALL
                n=ll(3)
                j=2
                FORALL (i=1:Nm1)
                   wpout(i,j,k)=tmp2_r(i,n)
                END FORALL
                o=ll(2)
                j=3
                FORALL (i=1:Nm1)
                   wpout(i,j,k)=tmp2_r(i,o)
                END FORALL
                p=ll(1)
                j=4
                FORALL (i=1:Nm1)
                   wpout(i,j,k)=tmp2_r(i,p)
                END FORALL
             ENDDO !K=1,NM3

             DEALLOCATE(tmp2_r,STAT=info)
             or_fail_dealloc("Failed to deallocate tmp2_r!")
          CASE (3)
          !!! Z direction
             IF (Nm3.LT.4) THEN
                fail("The number of pixels along dimension is less than 4!!! This filter requires a minimum of four pixels along the dimension to be processed!", &
                & ppm_error=ppm_error_fatal)
             ENDIF

             !At each pixel, I should keep three values of the first, second
             !and third row
             ALLOCATE(tmp2_r(Nm1,5),STAT=info)
             or_fail_alloc("tmp2_r")

             DO j=1,NM2
                ! This optimisation is done to avoid recomputing what has been computed once
                FORALL (i=1:Nm1)
                   ! Initialize borders
                   tmp2_r(i,2) = wpin(i,j,2)*this%N0
                   tmp2_r(i,3) = wpin(i,j,2)*this%N1
                   tmp2_r(i,4) = wpin(i,j,2)*this%N2
                END FORALL

                FORALL (i=1:Nm1)
                   ! Initialize borders
                   tmp2_r(i,3) = tmp2_r(i,3)+wpin(i,j,3)*this%N0
                   tmp2_r(i,4) = tmp2_r(i,4)+wpin(i,j,3)*this%N1
                END FORALL

                FORALL (i=1:Nm1)
                   ! Initialize borders
                   tmp2_r(i,4) = tmp2_r(i,4)+wpin(i,j,4)*this%N0
                END FORALL

                FORALL (i=1:Nm1)
                   ! Initialize borders
                   tmp2_r(i,1) = wpin(i,j,1)*this%N0 + wpin(i,j,1)*this%N1 + wpin(i,j,1)*this%N2 + wpin(i,j,1)*this%N3
                   tmp2_r(i,2) = tmp2_r(i,2)         + wpin(i,j,1)*this%N1 + wpin(i,j,1)*this%N2 + wpin(i,j,1)*this%N3
                   tmp2_r(i,3) = tmp2_r(i,3)                               + wpin(i,j,1)*this%N2 + wpin(i,j,1)*this%N3
                   tmp2_r(i,4) = tmp2_r(i,4)                                                     + wpin(i,j,1)*this%N3

                   ! note that the wpin(i,1) value is multiplied by the Boundary coefficients BNi
                   tmp2_r(i,1) = tmp2_r(i,1) - (wpin(i,j,1)*this%BN1 + wpin(i,j,1)*this%BN2 + wpin(i,j,1)*this%BN3 + wpin(i,j,1)*this%BN4)
                   tmp2_r(i,2) = tmp2_r(i,2) - (tmp2_r(i,1)*this%D1  + wpin(i,j,1)*this%BN2 + wpin(i,j,1)*this%BN3 + wpin(i,j,1)*this%BN4)
                   tmp2_r(i,3) = tmp2_r(i,3) - (tmp2_r(i,2)*this%D1  + tmp2_r(i,1)*this%D2  + wpin(i,j,1)*this%BN3 + wpin(i,j,1)*this%BN4)
                   tmp2_r(i,4) = tmp2_r(i,4) - (tmp2_r(i,3)*this%D1  + tmp2_r(i,2)*this%D2  + tmp2_r(i,1)*this%D3  + wpin(i,j,1)*this%BN4)
                END FORALL

                FORALL(i=1:5) ll(i)=i

                !Recursively filter the rest
                DO k=5,Nm3
                   l=ll(5)
                   FORALL (i=1:Nm1)
                      tmp2_r(i,l)=wpin(i,j,k)*this%N0
                   END FORALL
                   m=k-1
                   FORALL (i=1:Nm1)
                      tmp2_r(i,l)=tmp2_r(i,l)+wpin(i,j,m)*this%N1
                   END FORALL
                   m=k-2
                   FORALL (i=1:Nm1)
                      tmp2_r(i,l)=tmp2_r(i,l)+wpin(i,j,m)*this%N2
                   END FORALL
                   m=k-3
                   FORALL (i=1:Nm1)
                      tmp2_r(i,l)=tmp2_r(i,l)+wpin(i,j,m)*this%N3
                   END FORALL

                   m=ll(4)
                   n=ll(3)
                   o=ll(2)
                   p=ll(1)
                   FORALL (i=1:Nm1)
                      tmp2_r(i,l)=tmp2_r(i,l) - (tmp2_r(i,m)*this%D1 + tmp2_r(i,n)*this%D2 + tmp2_r(i,o)*this%D3 + tmp2_r(i,p)*this%D4)
                   END FORALL

                   !Store the causal result
                   l=k-4
                   FORALL (i=1:Nm1)
                      wpout(i,j,l)=tmp2_r(i,p)
                   END FORALL
                   ll=CSHIFT(ll,SHIFT=1,DIM=1)
                ENDDO !k=5,Nm3

                !Store the causal result for the last 4 lines
                p=ll(1)
                k=Nm3-3
                FORALL (i=1:Nm1)
                   wpout(i,j,k)=tmp2_r(i,p)
                END FORALL
                o=ll(2)
                k=Nm3-2
                FORALL (i=1:Nm1)
                   wpout(i,j,k)=tmp2_r(i,o)
                END FORALL
                n=ll(3)
                k=Nm3-1
                FORALL (i=1:Nm1)
                   wpout(i,j,k)=tmp2_r(i,n)
                END FORALL
                m=ll(4)
                k=Nm3
                FORALL (i=1:Nm1)
                   wpout(i,j,k)=tmp2_r(i,m)
                END FORALL

                !------------------------
                !AntiCausal direction pass
                !------------------------
                m=Nm3-1
                FORALL (i=1:Nm1)
                   ! Initialize borders
                   tmp2_r(i,3) = wpin(i,j,m)*this%M1
                   tmp2_r(i,4) = wpin(i,j,m)*this%M2
                END FORALL

                m=Nm3-2
                FORALL (i=1:Nm1)
                   ! Initialize borders
                   tmp2_r(i,4) = tmp2_r(i,4)+wpin(i,j,m)*this%M1
                END FORALL

                FORALL (i=1:Nm1)
                   ! Initialize borders
                   tmp2_r(i,1) = wpin(i,j,Nm3  )*this%M1 + wpin(i,j,Nm3  )*this%M2 + wpin(i,j,Nm3)*this%M3 + wpin(i,j,Nm3)*this%M4
                   tmp2_r(i,2) = wpin(i,j,Nm3  )*this%M1 + wpin(i,j,Nm3  )*this%M2 + wpin(i,j,Nm3)*this%M3 + wpin(i,j,Nm3)*this%M4
                   tmp2_r(i,3) = tmp2_r(i,3)             + wpin(i,j,Nm3  )*this%M2 + wpin(i,j,Nm3)*this%M3 + wpin(i,j,Nm3)*this%M4
                   tmp2_r(i,4) = tmp2_r(i,4)                                       + wpin(i,j,Nm3)*this%M3 + wpin(i,j,Nm3)*this%M4

                   ! note that the outV value is multiplied by the Boundary coefficients BMi
                   tmp2_r(i,1) = tmp2_r(i,1) - (wpin(i,j,Nm3)*this%BM1 + wpin(i,j,Nm3)*this%BM2 + wpin(i,j,Nm3)*this%BM3 + wpin(i,j,Nm3)*this%BM4)
                   tmp2_r(i,2) = tmp2_r(i,2) - (tmp2_r(i,  1)*this%D1  + wpin(i,j,Nm3)*this%BM2 + wpin(i,j,Nm3)*this%BM3 + wpin(i,j,Nm3)*this%BM4)
                   tmp2_r(i,3) = tmp2_r(i,3) - (tmp2_r(i,  2)*this%D1  + tmp2_r(i,  1)*this%D2  + wpin(i,j,Nm3)*this%BM3 + wpin(i,j,Nm3)*this%BM4)
                   tmp2_r(i,4) = tmp2_r(i,4) - (tmp2_r(i,  3)*this%D1  + tmp2_r(i,  2)*this%D2  + tmp2_r(i,  1)*this%D3  + wpin(i,j,Nm3)*this%BM4)
                END FORALL

                FORALL(i=1:5) ll(i)=i

                !Recursively filter the rest
                DO k=Nm3-3,2,-1
                   l=ll(5)
                   FORALL (i=1:Nm1)
                      tmp2_r(i,l)=wpin(i,j,k)*this%M1
                   END FORALL
                   m=k+1
                   FORALL (i=1:Nm1)
                      tmp2_r(i,l)=tmp2_r(i,l)+wpin(i,j,m)*this%M2
                   END FORALL
                   m=k+2
                   FORALL (i=1:Nm1)
                      tmp2_r(i,l)=tmp2_r(i,l)+wpin(i,j,m)*this%M3
                   END FORALL
                   m=k+3
                   FORALL (i=1:Nm1)
                      tmp2_r(i,l)=tmp2_r(i,l)+wpin(i,j,m)*this%M4
                   END FORALL

                   m=ll(4)
                   n=ll(3)
                   o=ll(2)
                   p=ll(1)
                   FORALL (i=1:Nm1)
                      tmp2_r(i,l)=tmp2_r(i,l) - (tmp2_r(i,m)*this%D1 + tmp2_r(i,n)*this%D2 + tmp2_r(i,o)*this%D3 + tmp2_r(i,p)*this%D4)
                   END FORALL

                   l=k+3
                   !Store the causal result
                   FORALL (i=1:Nm1)
                      wpout(i,j,l)=wpout(i,j,l)+tmp2_r(i,p)
                   END FORALL
                   ll=CSHIFT(ll,SHIFT=1,DIM=1)
                ENDDO !k=Nm3-3,2,-1

                !Store the causal result for the last 4 lines
                m=ll(4)
                FORALL (i=1:Nm1)
                   wpout(i,j,1)=tmp2_r(i,m)
                END FORALL
                n=ll(3)
                FORALL (i=1:Nm1)
                   wpout(i,j,2)=tmp2_r(i,n)
                END FORALL
                o=ll(2)
                FORALL (i=1:Nm1)
                   wpout(i,j,3)=tmp2_r(i,o)
                END FORALL
                p=ll(1)
                FORALL (i=1:Nm1)
                   wpout(i,j,4)=tmp2_r(i,p)
                END FORALL
             ENDDO !j=1,NM2

             DEALLOCATE(tmp2_r,STAT=info)
             or_fail_dealloc("Failed to deallocate tmp2_r!")
          CASE DEFAULT
             fail("Wrong direction for 3D problem!!!")
          END SELECT
#endif
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        END SUBROUTINE DTYPE(FilterDataArray)

#if   __DIME == __3D
        !
        SUBROUTINE ComputeNCoefficients(this,A1,B1,W1,L1,A2,B2,W2,L2, &
        &          N0,N1,N2,N3,SN,DN,EN)
        !!!Compute the N coefficients in the recursive filter.
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(RecursiveGaussianImageFilter) :: this

          REAL(MK),             INTENT(IN   ) :: A1
          REAL(MK),             INTENT(IN   ) :: B1
          REAL(MK),             INTENT(IN   ) :: W1
          REAL(MK),             INTENT(IN   ) :: L1
          REAL(MK),             INTENT(IN   ) :: A2
          REAL(MK),             INTENT(IN   ) :: B2
          REAL(MK),             INTENT(IN   ) :: W2
          REAL(MK),             INTENT(IN   ) :: L2
          REAL(MK),             INTENT(INOUT) :: N0
          REAL(MK),             INTENT(INOUT) :: N1
          REAL(MK),             INTENT(INOUT) :: N2
          REAL(MK),             INTENT(INOUT) :: N3
          REAL(MK),             INTENT(INOUT) :: SN
          REAL(MK),             INTENT(INOUT) :: DN
          REAL(MK),             INTENT(INOUT) :: EN

          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          REAL(MK) :: Sin1,Sin2,Cos1,Cos2,Exp1,Exp2

          Sin1 = SIN(W1/this%Sigma)
          Sin2 = SIN(W2/this%Sigma)
          Cos1 = COS(W1/this%Sigma)
          Cos2 = COS(W2/this%Sigma)
          Exp1 = EXP(L1/this%Sigma)
          Exp2 = EXP(L2/this%Sigma)

          N0 = A1 + A2
          N1 = Exp2 * ( B2 * Sin2 - ( A2 + two * A1 ) * Cos2 )
          N1 = N1 + Exp1 * ( B1 * Sin1 - ( A1 + two * A2 ) * Cos1 )
          N2 = ( A1 + A2 ) * Cos2 * Cos1
          N2 = N2 - (B1 * Cos2 * Sin1 + B2 * Cos1 * Sin2)
          N2 = N2 * two * Exp1 * Exp2
          N2 = N2 + (A2 * Exp1 * Exp1 + A1 * Exp2 * Exp2)
          N3 = Exp2 * Exp1 * Exp1 * ( B2 * Sin2 - A2 * Cos2 )
          N3 = N3 + Exp1 * Exp2 * Exp2 * ( B1 * Sin1 - A1 * Cos1 )

          SN = N0 + N1 + N2 + N3
          DN = N1 + two * N2 + three * N3
          EN = N1 + 4.0_MK * N2 + 9._MK * N3

        END SUBROUTINE ComputeNCoefficients

        !
        SUBROUTINE ComputeDCoefficients(this,W1,L1,W2,L2,SD,DD,ED)
        !!!Compute the D coefficients in the recursive filter.
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(RecursiveGaussianImageFilter) :: this

          REAL(MK),             INTENT(IN   ) :: W1
          REAL(MK),             INTENT(IN   ) :: L1
          REAL(MK),             INTENT(IN   ) :: W2
          REAL(MK),             INTENT(IN   ) :: L2
          REAL(MK),             INTENT(INOUT) :: SD
          REAL(MK),             INTENT(INOUT) :: DD
          REAL(MK),             INTENT(INOUT) :: ED

          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          REAL(MK) :: Cos1,Cos2,Exp1,Exp2

          Cos1 = COS(W1/this%Sigma)
          Cos2 = COS(W2/this%Sigma)
          Exp1 = EXP(L1/this%Sigma)
          Exp2 = EXP(L2/this%Sigma)

          this%D4 = Exp1 * Exp1 * Exp2 * Exp2
          this%D3 = -two * Cos1 * Exp1 * Exp2 * Exp2
          this%D3 = this%D3 - two * Cos2 * Exp2 * Exp1 * Exp1
          this%D2 = 4.0_MK * Cos2 * Cos1 * Exp1 * Exp2
          this%D2 = this%D2 + (Exp1 * Exp1 + Exp2 * Exp2)
          this%D1 = -two * ( Exp2 * Cos2 + Exp1 * Cos1 )

          SD = one + this%D1 + this%D2 + this%D3 + this%D4
          DD = this%D1 + two * this%D2 + three * this%D3 + 4.0_MK * this%D4
          ED = this%D1 + 4.0_MK * this%D2 + 9._MK * this%D3 + 16._MK * this%D4

        END SUBROUTINE ComputeDCoefficients

        !
        SUBROUTINE ComputeRemainingCoefficients(this,symmetric)
        !!! Compute the M coefficients and the boundary coefficients in the
        !!! recursive filter.
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(RecursiveGaussianImageFilter) :: this

          LOGICAL,              INTENT(IN   ) :: symmetric

          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          REAL(MK) :: SN,SM,SD

          SELECT CASE (symmetric)
          CASE (.TRUE.)
             this%M1 = this%N1 - this%D1 * this%N0
             this%M2 = this%N2 - this%D2 * this%N0
             this%M3 = this%N3 - this%D3 * this%N0
             this%M4 =         - this%D4 * this%N0

          CASE DEFAULT
             this%M1 = -(this%N1 - this%D1 * this%N0)
             this%M2 = -(this%N2 - this%D2 * this%N0)
             this%M3 = -(this%N3 - this%D3 * this%N0)
             this%M4 =             this%D4 * this%N0

          END SELECT

          ! Compute coefficients to be used at the boundaries
          ! in order to simulate edge extension boundary conditions.
          SN = this%N0 + this%N1 + this%N2 + this%N3
          SM = this%M1 + this%M2 + this%M3 + this%M4
          SD = one + this%D1 + this%D2 + this%D3 + this%D4

          this%BN1 = this%D1 * SN / SD
          this%BN2 = this%D2 * SN / SD
          this%BN3 = this%D3 * SN / SD
          this%BN4 = this%D4 * SN / SD

          this%BM1 = this%D1 * SM / SD
          this%BM2 = this%D2 * SM / SD
          this%BM3 = this%D3 * SM / SD
          this%BM4 = this%D4 * SM / SD

        END SUBROUTINE ComputeRemainingCoefficients

        !
        SUBROUTINE RecursiveGaussianImageFilter_SetUp(this,info)
        !!!Compute filter for Gaussian kernel
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(RecursiveGaussianImageFilter) :: this

          INTEGER,              INTENT(  OUT) :: info
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          !Parameters of exponential series
          REAL(MK), PARAMETER :: A1(3)=(/ 1.3530_MK,-0.6724_MK,-1.3563_MK/)
          REAL(MK), PARAMETER :: B1(3)=(/ 1.8151_MK,-3.4327_MK, 5.2318_MK/)
          REAL(MK), PARAMETER :: W1   =   0.6681_MK
          REAL(MK), PARAMETER :: L1   =  -1.3932_MK
          REAL(MK), PARAMETER :: A2(3)=(/-0.3531_MK, 0.6724_MK, 0.3446_MK/)
          REAL(MK), PARAMETER :: B2(3)=(/ 0.0902_MK, 0.6100_MK,-2.2355_MK/)
          REAL(MK), PARAMETER :: W2   =   2.0787_MK
          REAL(MK), PARAMETER :: L2   =  -1.3732_MK

          REAL(MK)            :: SD, DD, ED
          REAL(MK)            :: SN, DN, EN
          REAL(MK)            :: alpha,beta

          REAL(MK)            :: N0_0, N1_0, N2_0, N3_0
          REAL(MK)            :: N0_2, N1_2, N2_2, N3_2
          REAL(MK)            :: SN0, DN0, EN0
          REAL(MK)            :: SN2, DN2, EN2

          LOGICAL :: symmetric

          start_subroutine("RecursiveGaussianImageFilter_SetUp")

          CALL check()

          CALL this%ComputeDCoefficients(W1,L1,W2,L2,SD,DD,ED)

          SELECT CASE (this%Order)
          CASE (0)
             ! Approximation of convolution with a gaussian.
             CALL this%ComputeNCoefficients(A1(1),B1(1),W1,L1,A2(1),B2(1), &
             &    W2,L2,this%N0,this%N1,this%N2,this%N3,SN,DN,EN)

             alpha = 2 * SN / SD - this%N0
             this%N0=this%N0 / alpha
             this%N1=this%N1 / alpha
             this%N2=this%N2 / alpha
             this%N3=this%N3 / alpha

             symmetric=.TRUE.

             CALL this%ComputeRemainingCoefficients(symmetric)
          CASE (1)
             !Approximation of convolution with the first derivative of a Gaussian
             CALL this%ComputeNCoefficients(A1(2),B1(2),W1,L1,A2(2),B2(2), &
             &    W2,L2,this%N0,this%N1,this%N2,this%N3,SN,DN,EN)

             alpha = 2 * ( SN * DD - DN * SD ) / ( SD * SD )

             IF (this%NormalizeAcrossScale) THEN
                this%N0=this%N0 * this%Sigma / alpha
                this%N1=this%N1 * this%Sigma / alpha
                this%N2=this%N2 * this%Sigma / alpha
                this%N3=this%N3 * this%Sigma / alpha
             ELSE
                this%N0=this%N0 / alpha
                this%N1=this%N1 / alpha
                this%N2=this%N2 / alpha
                this%N3=this%N3 / alpha
             ENDIF

             symmetric=.FALSE.

             CALL this%ComputeRemainingCoefficients(symmetric)
          CASE (2)
             !Approximation of convolution with the second derivative of a
             !Gaussian.
             CALL this%ComputeNCoefficients(A1(1),B1(1),W1,L1,A2(1),B2(1), &
             &    W2,L2,N0_0,N1_0,N2_0,N3_0,SN0,DN0,EN0)

             CALL this%ComputeNCoefficients(A1(3),B1(3),W1,L1,A2(3),B2(3), &
             &    W2,L2,N0_2,N1_2,N2_2,N3_2,SN2,DN2,EN2)

             beta = -( two * SN2 - SD * N0_2 ) / ( two * SN0 - SD * N0_0 )

             this%N0 = N0_2 + beta * N0_0
             this%N1 = N1_2 + beta * N1_0
             this%N2 = N2_2 + beta * N2_0
             this%N3 = N3_2 + beta * N3_0

             SN = SN2 + beta * SN0
             DN = DN2 + beta * DN0
             EN = EN2 + beta * EN0

             alpha = EN * SD * SD - ED * SN * SD - 2 * DN * DD * SD + 2 * DD * DD * SN
             alpha = alpha / (SD * SD * SD)

             IF (this%NormalizeAcrossScale) THEN
                this%N0=this%N0 * this%Sigma * this%Sigma / alpha
                this%N1=this%N1 * this%Sigma * this%Sigma / alpha
                this%N2=this%N2 * this%Sigma * this%Sigma / alpha
                this%N3=this%N3 * this%Sigma * this%Sigma / alpha
             ELSE
                this%N0=this%N0 / alpha
                this%N1=this%N1 / alpha
                this%N2=this%N2 / alpha
                this%N3=this%N3 / alpha
             ENDIF

             symmetric=.TRUE.

             CALL this%ComputeRemainingCoefficients(symmetric)
          CASE DEFAULT
             fail("Unknown Order!")
          END SELECT

          end_subroutine()
        CONTAINS
          SUBROUTINE check
            IF (this%Sigma.LE.zero) THEN
               fail("Sigma must be greater than zero.",exit_point=7777)
            ENDIF
            7777 CONTINUE
            RETURN
          END SUBROUTINE check
        END SUBROUTINE RecursiveGaussianImageFilter_SetUp


        SUBROUTINE GenerateDiscreteGaussianCoefficients(direction,Mask,Sigma,MaximumKernelWidth,KernelWidth,info,MaximumError)
          !!! Compute filter for Gaussian kernel
          !!!
          !!!
          !!!
          !!! References:
          !!! The function is adapted and amended from ITK software implementation
          !!!
          !!! The Gaussian kernel contained in this operator was described
          !!! by Tony Lindeberg (Discrete Scale-Space Theory and the Scale-Space
          !!! Primal Sketch.  Dissertation. Royal Institute of Technology, Stockholm,
          !!! Sweden. May 1991.).
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          INTEGER,                             INTENT(IN   ) :: direction

          REAL(MK), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: Mask
          REAL(MK),                            INTENT(IN   ) :: Sigma

          INTEGER,                             INTENT(IN   ) :: MaximumKernelWidth
          INTEGER,                             INTENT(  OUT) :: KernelWidth
          !!! The final size of the Kernel which is equal to MaximumKernelWidth
          !!! This is chnaging when there is MaximumError
          INTEGER,                             INTENT(  OUT) :: info

          REAL(MK), OPTIONAL,                  INTENT(IN   ) :: MaximumError
          !!! The "maximum error" allowed in the discrete Gaussian function.
          !!! "Maximum errror" is defined as the difference between the area
          !!! under the discrete Gaussian curve and the area under the continuous
          !!! Gaussian. Maximum error affects the Gaussian operator size. Care should
          !!! be taken not to make this value too small relative to the standard
          !!! deviation, lest the operator size become unreasonably large.
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------

          REAL(ppm_kind_double), DIMENSION(:), ALLOCATABLE :: Maskd
          REAL(ppm_kind_double)                            :: et,cap
          REAL(ppm_kind_double)                            :: Sigmad
          REAL(ppm_kind_double)                            :: msumd
          REAL(MK)                                         :: dd
          REAL(MK)                                         :: coef,coef1
          REAL(MK)                                         :: msum

          INTEGER :: i

          start_subroutine("GenerateDiscreteGaussianCoefficients")

          IF (MaximumKernelWidth.LT.1) THEN
             fail("MaximumKernelWidth can not be less than 1",ppm_error=ppm_error_fatal)
          ENDIF

          IF (PRESENT(MaximumError)) THEN
             IF (MaximumKernelWidth.LT.2) THEN
                fail("MaximumKernelWidth can not be less than 2 when using MaximumError to adjustthe Kernel size", &
                & ppm_error=ppm_error_fatal)
             ENDIF

             Sigmad=REAL(Sigma,ppm_kind_double)

             et=EXP(-Sigmad*Sigmad)
             cap=oned-REAL(MaximumError,ppm_kind_double)

             ALLOCATE(Maskd(0:MaximumKernelWidth),STAT=info)
             or_fail_alloc("Failed to allocate Maskd")

             Maskd(0)=et*ModifiedBesselI0(Sigmad)
             Maskd(1)=et*ModifiedBesselI1(Sigmad)
             msumd=Maskd(0)+Maskd(1)
             i=1
             DO WHILE (msumd.LT.cap)
                i=i+1
                Maskd(i)=et*ModifiedBesselI(i,Sigmad)
                msumd=msumd+Maskd(i)*twod
                IF (Maskd(i).LE.zerod) EXIT
                IF (i.EQ.MaximumKernelWidth) THEN
                   IF (rank.EQ.0) THEN
                      stdout("Kernel size has exceeded the specified maximum width of ",MaximumKernelWidth, &
                      & " and has been truncated to ",MaximumKernelWidth," elements.")
                   ENDIF
                   EXIT
                ENDIF
             ENDDO
             ! Normalize the coefficients so that their sum is one.
             Maskd=Maskd/msumd
             KernelWidth=i

             IF (ALLOCATED(Mask)) THEN
                IF (SIZE(Mask).NE.KernelWidth*2+1) THEN
                   DEALLOCATE(Mask,STAT=info)
                   or_fail_dealloc("Mask")

                   ALLOCATE(Mask(-KernelWidth:KernelWidth),STAT=info)
                   or_fail_alloc("Failed to allocate Mask")
                ENDIF
             ELSE
                ALLOCATE(Mask(-KernelWidth:KernelWidth),STAT=info)
                or_fail_alloc("Failed to allocate Mask")
             ENDIF

             Mask(0)=REAL(Maskd(0),MK)
             DO i=1,KernelWidth
                Mask( i)=REAL(Maskd(i),MK)
                Mask(-i)=Mask(i)
             ENDDO

             DEALLOCATE(Maskd,STAT=info)
             or_fail_dealloc("Failed to Deallocate Mask")
          ELSE
             KernelWidth=MaximumKernelWidth

             IF (ALLOCATED(Mask)) THEN
                IF (SIZE(Mask).NE.KernelWidth*2+1) THEN
                   DEALLOCATE(Mask,STAT=info)
                   or_fail_dealloc("Mask")

                   ALLOCATE(Mask(-KernelWidth:KernelWidth),STAT=info)
                   or_fail_alloc("Failed to allocate Mask")
                ENDIF
             ELSE
                ALLOCATE(Mask(-KernelWidth:KernelWidth),STAT=info)
                or_fail_alloc("Failed to allocate Mask")
             ENDIF

             dd=pixel(1)/pixel(direction)*pixel(1)/pixel(direction)

             !Gaussians = 1/(sqrt(2*pi)*sigma).*exp(-x.^2/(2*sigma.^2))
             coef =one/(SQRT(two*pi)*Sigma)
             coef1=one/(two*Sigma*Sigma)/dd
             FORALL (i=-KernelWidth:KernelWidth) Mask(-i)= coef * EXP(-REAL(i*i,MK)*coef1)
             msum=SUM(Mask)
             Mask=Mask/msum
          ENDIF

          end_subroutine()

        CONTAINS
          FUNCTION ModifiedBesselI0(Sigma_) RESULT(accumulator)
            IMPLICIT NONE
            REAL(ppm_kind_double), INTENT(IN   ) :: Sigma_
            REAL(ppm_kind_double)                :: accumulator
            !-------------------------------------------------------------------------
            !  Local variables
            !-------------------------------------------------------------------------
            REAL(ppm_kind_double), PARAMETER :: a1= 3.5156229_ppm_kind_double
            REAL(ppm_kind_double), PARAMETER :: a2= 3.0899424_ppm_kind_double
            REAL(ppm_kind_double), PARAMETER :: a3= 1.2067492_ppm_kind_double
            REAL(ppm_kind_double), PARAMETER :: a4= 1.2659732_ppm_kind_double
            REAL(ppm_kind_double), PARAMETER :: a5= 0.0360768_ppm_kind_double
            REAL(ppm_kind_double), PARAMETER :: a6= 0.0045813_ppm_kind_double
            REAL(ppm_kind_double), PARAMETER :: b1= 0.39894228_ppm_kind_double
            REAL(ppm_kind_double), PARAMETER :: b2= 0.01328592_ppm_kind_double
            REAL(ppm_kind_double), PARAMETER :: b3= 0.00225319_ppm_kind_double
            REAL(ppm_kind_double), PARAMETER :: b4=-0.00157565_ppm_kind_double
            REAL(ppm_kind_double), PARAMETER :: b5= 0.00916281_ppm_kind_double
            REAL(ppm_kind_double), PARAMETER :: b6=-0.02057706_ppm_kind_double
            REAL(ppm_kind_double), PARAMETER :: b7=-0.02635537_ppm_kind_double
            REAL(ppm_kind_double), PARAMETER :: b8=-0.01647633_ppm_kind_double
            REAL(ppm_kind_double), PARAMETER :: b9= 0.00392377_ppm_kind_double
            REAL(ppm_kind_double)            :: m

            IF (Sigma_*Sigma_.LT.3.75_ppm_kind_double) THEN
               m=Sigma_*Sigma_*Sigma_*Sigma_/14.0625_ppm_kind_double
               accumulator=oned+m*(a1+m*(a2+m*(a3+m*(a4+m*(a5+m*a6)))))
            ELSE
               m=3.75_ppm_kind_double/Sigma_/Sigma_
               accumulator=(EXP(Sigma_*Sigma_)/Sigma_)*(b1+m*(b2+m*(b3+m*(b4+m*(b5+m*(b6+m*(b7+m*(b8+m*b9))))))))
            ENDIF
            RETURN
          END FUNCTION ModifiedBesselI0

          FUNCTION ModifiedBesselI1(Sigma_) RESULT(accumulator)
            IMPLICIT NONE
            REAL(ppm_kind_double), INTENT(IN   ) :: Sigma_
            REAL(ppm_kind_double)                :: accumulator
            !-------------------------------------------------------------------------
            !  Local variables
            !-------------------------------------------------------------------------
            REAL(ppm_kind_double), PARAMETER :: a1= 0.87890594_ppm_kind_double
            REAL(ppm_kind_double), PARAMETER :: a2= 0.51498869_ppm_kind_double
            REAL(ppm_kind_double), PARAMETER :: a3= 0.15084934_ppm_kind_double
            REAL(ppm_kind_double), PARAMETER :: a4= 0.02658733_ppm_kind_double
            REAL(ppm_kind_double), PARAMETER :: a5= 0.00301532_ppm_kind_double
            REAL(ppm_kind_double), PARAMETER :: a6= 0.00032411_ppm_kind_double
            REAL(ppm_kind_double), PARAMETER :: b1= 0.39894228_ppm_kind_double
            REAL(ppm_kind_double), PARAMETER :: b2=-0.03988024_ppm_kind_double
            REAL(ppm_kind_double), PARAMETER :: b3=-0.00362018_ppm_kind_double
            REAL(ppm_kind_double), PARAMETER :: b4= 0.00163801_ppm_kind_double
            REAL(ppm_kind_double), PARAMETER :: b5=-0.01031555_ppm_kind_double
            REAL(ppm_kind_double), PARAMETER :: b6= 0.02282967_ppm_kind_double
            REAL(ppm_kind_double), PARAMETER :: b7=-0.02895312_ppm_kind_double
            REAL(ppm_kind_double), PARAMETER :: b8= 0.01787654_ppm_kind_double
            REAL(ppm_kind_double), PARAMETER :: b9= 0.00420059_ppm_kind_double
            REAL(ppm_kind_double)            :: m

            IF (Sigma_*Sigma_.LT.3.75_ppm_kind_double) THEN
               m=Sigma_*Sigma_*Sigma_*Sigma_/14.0625_ppm_kind_double
               accumulator=Sigma_*Sigma_*(halfd+m*(a1+m*(a2+m*(a3+m*(a4+m*(a5+m*a6))))))
            ELSE
               m=3.75_ppm_kind_double/Sigma_/Sigma_
               accumulator=(EXP(Sigma_*Sigma_)/Sigma_)*(b1+m*(b2+m*(b3+m*(b4+m*(b5+m*b6+m*(b7+m*(b8-m*b9)))))))
            ENDIF
            RETURN
          END FUNCTION ModifiedBesselI1

          FUNCTION ModifiedBesselI(n,Sigma_) RESULT(accumulator)
            IMPLICIT NONE

            INTEGER,               INTENT(IN   ) :: n

            REAL(ppm_kind_double), INTENT(IN   ) :: Sigma_
            REAL(ppm_kind_double)                :: accumulator
            !-------------------------------------------------------------------------
            !  Local variables
            !-------------------------------------------------------------------------
            REAL(ppm_kind_double), PARAMETER :: a=1.0E10_ppm_kind_double
            REAL(ppm_kind_double), PARAMETER :: b=1.0E-10_ppm_kind_double
            REAL(ppm_kind_double)            :: qim,qi,qip,toy

            INTEGER :: l

            toy=twod/Sigma_/Sigma_
            qip=zerod
            qi=oned
            accumulator=zerod

            DO l=2*(n+INT(SQRT(40._MK*REAL(n,MK)))),2,-1
               qim=qip+REAL(l,ppm_kind_double)*toy*qi
               qip=qi
               qi=qim
               IF (ABS(qi).GT.a) THEN
                  accumulator=accumulator*b
                  qi=qi*b
                  qip=qip*b
               ENDIF
               IF (l.EQ.n) accumulator=qip
            ENDDO
            accumulator=accumulator*ModifiedBesselI0(Sigma_)/qi
            RETURN
          END FUNCTION ModifiedBesselI
        END SUBROUTINE GenerateDiscreteGaussianCoefficients

#endif




