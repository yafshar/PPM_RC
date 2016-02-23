
        SUBROUTINE ppm_rc_DetectOscillations(info)

        !TODO
        !TOCHECK
        !!!!It should be checked later

          IMPLICIT NONE
          INTEGER, INTENT(  OUT) :: info
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          REAL(MK), DIMENSION(1:OscillationHistoryLength), SAVE :: m_OscillationsEnergyHist=0.0_MK
          REAL(MK), PARAMETER                                     :: c=0.00001_MK
          REAL(MK)                                                :: vSum,vSumOld

          INTEGER, DIMENSION(1:OscillationHistoryLength), SAVE :: m_OscillationsNumberHist=0
          INTEGER                                                :: nsize,i

          start_subroutine("ppm_rc_DetectOscillations")

          vSum=SUM(energya(Candidates_list))
          nsize=SIZE(Candidates_list)

          DO i=1,OscillationHistoryLength
             vSumOld=m_OscillationsEnergyHist(i)
             IF (nsize.EQ.m_OscillationsNumberHist(i).AND. &
             &   ABS(vSum-vSumOld).LE.c*ABS(vSum)) THEN
                !here we assume that we're oscillating,so we decrease the acceptance factor:
                AcceptedPointsFactor=AcceptedPointsFactor*AcceptedPointsReductionFactor
             ENDIF
          ENDDO
          m_OscillationsEnergyHist=EOSHIFT(m_OscillationsEnergyHist,SHIFT=1)
          m_OscillationsNumberHist=EOSHIFT(m_OscillationsNumberHist,SHIFT=1)

          !Fill the new elements:
          m_OscillationsEnergyHist(OscillationHistoryLength)=vSum
          m_OscillationsNumberHist(OscillationHistoryLength)=nsize

          end_subroutine()

        END SUBROUTINE ppm_rc_DetectOscillations

        SUBROUTINE ppm_rc_DetectOscillations2(info)

          IMPLICIT NONE
          INTEGER, INTENT(  OUT) :: info
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          REAL(ppm_kind_double), DIMENSION(:), ALLOCATABLE, SAVE :: vSums
          REAL(ppm_kind_double), PARAMETER                       :: alpha = 0.1_ppm_kind_double
          REAL(ppm_kind_double)                                  :: vSum,vSumOld,vSumNew
          REAL(ppm_kind_double)                                  :: winvar,totvar,fac

          INTEGER       :: i,s
          INTEGER, SAVE :: iter_id=0

          LOGICAL, SAVE :: new=.TRUE.

          start_subroutine("ppm_rc_DetectOscillations2")

          s=SIZE(Candidates_list)

          vSum=0.0_ppm_kind_double
          DO i=1,s
             vSum=vSum+REAL(energya(Candidates_list(i)),ppm_kind_double)
          ENDDO

          IF (new) THEN
             ALLOCATE(vSums(maxiter),STAT=info)
             or_fail_alloc("vSums",ppm_error=ppm_error_fatal)
             vSums  =zerod
             vSumOld=vSum
          ELSE
             vSumOld=vSums(iter_id)
          ENDIF

          iter_id=iter_id+1
          vSumNew=alpha*vSum+(oned-alpha)*vSumOld
          vSums(iter_id)=vSumNew

          IF (new) THEN
             new=.FALSE.
          ELSE
             totvar=var(vSums,1,iter_id)
             IF (totvar.LE.smalld) RETURN

             i=MAX(1,iter_id-OscillationHistoryLength)
             winvar=var(vSums,i,iter_id)
             fac=SQRT(winvar/totvar)
             IF (fac.LT.OscillationThreshold) THEN
                IF (AcceptedPointsFactor.LT.lmyeps) THEN
                   AcceptedPointsFactor=zero
                ELSE
                   AcceptedPointsFactor=AcceptedPointsFactor*AcceptedPointsReductionFactor
                ENDIF
             ENDIF
          ENDIF

          end_subroutine()

        CONTAINS

          !Computes the variance of a(st:en) variables
          FUNCTION var(a,st,en)
            IMPLICIT NONE
            REAL(ppm_kind_double), DIMENSION(:), INTENT(IN   ) :: a
            INTEGER,                             INTENT(IN   ) :: st
            INTEGER,                             INTENT(IN   ) :: en
            REAL(ppm_kind_double)                              :: var

            REAL(ppm_kind_double), DIMENSION(en-st+1) :: b
            REAL(ppm_kind_double)                     :: nr,mean
            nr=REAL(en-st+1,ppm_kind_double)
            mean=SUM(a(st:en))/nr
            b=mean-a(st:en)
            var=DOT_PRODUCT(b,b)/nr
          END FUNCTION var

        END SUBROUTINE ppm_rc_DetectOscillations2

