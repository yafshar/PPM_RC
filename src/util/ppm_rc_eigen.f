      !-------------------------------------------------------------------------
      !  Subroutine   :                    ppm_rc_id
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Calculates the eigenvalues of a symmetric 3x3 matrix A
      !               : for symmetric 3x3 matrices.
      !
      !  Input        :
      !
      !  Input/output :
      !
      !  Output       : info (I) return status. 0 on success.
      !
      !  Routines     :
      !
      !
      !  Remarks      :
      !
      !  References   : Joachim Kopp
      !                 Numerical diagonalization of hermitian 3x3 matrices
      !                 arXiv.org preprint: physics/0610206
      !                 Int. J. Mod. Phys. C19 (2008) 523-548
      !
      !  Revisions    : Y.Afshar
      !-------------------------------------------------------------------------
      !  MOSAIC Group
      !  Max Planck Institute of Molecular Cell Biology and Genetics
      !  Pfotenhauerstr. 108, 01307 Dresden, Germany
      !
      !  Author           - y.afshar           June   2014
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_rc_CalcEigenvalueAnalytical(Matrix,EigenVals)
        !!! Calculates the eigenvalues of a symmetric 3x3 matrix M using Cardano's
        !!! analytical algorithm.
        !!! Only the diagonal and upper triangular parts of M are accessed.
        !!! The access is read-only.
        IMPLICIT NONE
        !-------------------------------------------------------------------------
        !  Includes
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        REAL(MK), DIMENSION(:), INTENT(IN   ) :: Matrix
        !!! M: The symmetric input matrix
        !!! M Elements are M11,M12,M21,M13,M23,M33
        REAL(MK), DIMENSION(:), INTENT(  OUT) :: EigenVals
        !!! W: Storage buffer for eigenvalues
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
        REAL(ppm_kind_double), PARAMETER    :: rsqrt3=0.57735026918962576451_ppm_kind_double
        REAL(ppm_kind_double), PARAMETER    :: onethird=0.33333333333333333333_ppm_kind_double
        REAL(ppm_kind_double), DIMENSION(6) :: A
        REAL(ppm_kind_double), DIMENSION(3) :: W
        REAL(ppm_kind_double)               :: A12A23
        REAL(ppm_kind_double)               :: A12A12
        REAL(ppm_kind_double)               :: A23A23
        REAL(ppm_kind_double)               :: A13A13
        REAL(ppm_kind_double)               :: TR,C1,C0
        REAL(ppm_kind_double)               :: P,sqrtP
        REAL(ppm_kind_double)               :: Q,C
        REAL(ppm_kind_double)               :: S,PHI
        REAL(ppm_kind_double)               :: tmp
        !!! Determine coefficients of characteristic poynomial. We write
        !!!           | A11   A12   A13 | | A(1) A(2) A(4) |
        !!!      A =  |       A22   A23 | |      A(3) A(5) |
        !!!           |             A33 | |           A(6) |
        A=REAL(Matrix,ppm_kind_double)

        A11A22=A(1)*A(3)
        A12A12=A(2)*A(2)
        A23A23=A(5)*A(5)
        A13A13=A(4)*A(4)

        TR=A(1)+A(3)+A(6)

        C1=(A11A22+A(1)*A(6)+A(3)*A(6))-(A12A12+A23A23+A13A13)
        C0=A(6)*A12A12+A(1)*A23A23+A(3)*A13A13-A11A22*A(6)-twod*A(4)*A(2)*A(5)

        P=TR*TR-threed*C1
        Q=TR*(P-1.5_ppm_kind_double*C1)-13.5_ppm_kind_double*C0

        sqrtP=SQRT(ABS(P))

        PHI=27.0_ppm_kind_double*(0.25_ppm_kind_double*C1*C1*(P-C1)+C0*(Q+6.75_ppm_kind_double*C0))
        PHI=onethird*ATAN2(SQRT(ABS(PHI)),Q)
        C=sqrtP*COS(PHI)
        S=rsqrt3*sqrtP*SIN(PHI)

        W(2)=onethird*(TR-C)
        W(3)=W(2)+S
        W(1)=W(2)+C
        W(2)=W(2)-S

        !!! Order eigen values
        !!! lambda_1 < lambda_2 < lambda_3
        tmp=W(2)
        IF (W(1).GT.W(2)) THEN
           W(2)=W(1)
           W(1)=tmp
           tmp=W(2)
        ENDIF
        IF (tmp.GT.W(3)) THEN
           IF (W(1).GT.W(3)) THEN
              W(2)=W(1)
              W(1)=W(3)
           ELSE
              W(2)=W(3)
           ENDIF
           W(3)=tmp
        ENDIF

        EigenVals=REAL(W,MK)
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        RETURN
      END SUBROUTINE ppm_rc_CalcEigenvalueAnalytical

