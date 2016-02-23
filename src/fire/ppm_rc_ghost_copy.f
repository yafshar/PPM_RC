      !-------------------------------------------------------------------------
      !  Subroutine   :                    ppm_rc_id
      !-------------------------------------------------------------------------
      !
      !  Purpose      :
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
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  MOSAIC Group
      !  Max Planck Institute of Molecular Cell Biology and Genetics
      !  Pfotenhauerstr. 108, 01307 Dresden, Germany
      !
      !  Author           - y.afshar           June   2014
      !-------------------------------------------------------------------------

      SUBROUTINE DTYPE(ppm_rc_ghost_copy)(sbpitr,reshapedghost,info)

        IMPLICIT NONE
        !-------------------------------------------------------------------------
        !  Includes
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        CLASS(ppm_t_subpatch_), POINTER       :: sbpitr

        INTEGER,  DIMENSION(:), INTENT(INOUT) :: reshapedghost
        INTEGER,                INTENT(  OUT) :: info
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
        REAL(ppm_kind_double) :: t0

        INTEGER,             DIMENSION(:),     POINTER :: Nm
#if   __DIME == __2D
        INTEGER, CONTIGUOUS, DIMENSION(:,:),   POINTER :: DTYPE(wp)
#elif __DIME == __3D
        INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER :: DTYPE(wp)
        INTEGER                                        :: k
#endif
        INTEGER                                        :: i,j,l,m

        CHARACTER(ppm_char) :: caller="ppm_rc_ghost_copy"

        !-------------------------------------------------------------------------
        !  Initialise
        !-------------------------------------------------------------------------
        CALL substart(caller,t0,info)

        Nm => sbpitr%nnodes

        NULLIFY(DTYPE(wp))

        CALL sbpitr%get_field(labels,DTYPE(wp),info)
        or_fail("Failed to get field wp data.")

        m=1

#if   __DIME == __2D
        !-x west
        IF (sbpitr%istart(1).NE.1) THEN
           l=m
           i=0
           DO j=1,Nm(2)
              reshapedghost(l)=wp_2d(i,j)
              l=l+1
           ENDDO !j=1,Nm(2)
           m=m+Nm(2)
        ENDIF
        !+x east
        IF (sbpitr%iend(1).NE.Ngrid(1)) THEN
           l=m
           i=Nm(1)+1
           DO j=1,Nm(2)
              reshapedghost(l)=wp_2d(i,j)
              l=l+1
           ENDDO !j=1,Nm(2)
           m=m+Nm(2)
        ENDIF
        !-y south
        IF (sbpitr%istart(2).NE.1) THEN
           l=m
           j=0
           DO i=0,Nm(1)+1
              reshapedghost(l)=wp_2d(i,j)
              l=l+1
           ENDDO !j=1,Nm(2)
           m=m+Nm(1)+2
        ENDIF
        !+y north
        IF (sbpitr%iend(2).NE.Ngrid(2)) THEN
           l=m
           j=Nm(2)+1
           DO i=0,Nm(1)+1
              reshapedghost(l)=wp_2d(i,j)
              l=l+1
           ENDDO !j=1,Nm(2)
        ENDIF
#elif __DIME == __3D
        !-x west
        IF (sbpitr%istart(1).NE.1) THEN
           l=m
           i=0
           DO k=1,Nm(3)
              DO j=1,Nm(2)
                 reshapedghost(l)=wp_3d(i,j,k)
                 l=l+1
              ENDDO !i=1,Nm(2)
           ENDDO !j=1,Nm(3)
           m=m+Nm(2)*Nm(3)
        ENDIF
        !+x east
        IF (sbpitr%iend(1).NE.Ngrid(1)) THEN
           l=m
           i=Nm(1)+1
           DO k=1,Nm(3)
              DO j=1,Nm(2)
                 reshapedghost(l)=wp_3d(i,j,k)
                 l=l+1
              ENDDO !i=1,Nm(2)
           ENDDO !j=1,Nm(3)
           m=m+Nm(2)*Nm(3)
        ENDIF
        !-y south
        IF (sbpitr%istart(2).NE.1) THEN
           l=m
           j=0
           DO k=1,Nm(3)
              DO i=0,Nm(1)+1
                 reshapedghost(l)=wp_3d(i,j,k)
                 l=l+1
              ENDDO !i=1,Nm(2)
           ENDDO !j=1,Nm(3)
           m=m+(Nm(1)+2)*Nm(3)
        ENDIF
        !+y north
        IF (sbpitr%iend(2).NE.Ngrid(2)) THEN
           l=m
           j=Nm(2)+1
           DO k=1,Nm(3)
              DO i=0,Nm(1)+1
                 reshapedghost(l)=wp_3d(i,j,k)
                 l=l+1
              ENDDO !i=1,Nm(2)
           ENDDO !j=1,Nm(3)
           m=m+(Nm(1)+2)*Nm(3)
        ENDIF
        !-z bottom
        IF (sbpitr%istart(3).NE.1) THEN
           l=m
           k=0
           DO j=0,Nm(2)+1
              DO i=0,Nm(1)+1
                 reshapedghost(l)=wp_3d(i,j,k)
                 l=l+1
              ENDDO !i=1,Nm(2)
           ENDDO !j=1,Nm(3)
           m=m+(Nm(1)+2)*(Nm(2)+2)
        ENDIF
        !+z top
        IF (sbpitr%iend(3).NE.Ngrid(3)) THEN
           l=m
           k=Nm(3)+1
           DO j=0,Nm(2)+1
              DO i=0,Nm(1)+1
                 reshapedghost(l)=wp_3d(i,j,k)
                 l=l+1
              ENDDO !i=1,Nm(2)
           ENDDO !j=1,Nm(3)
        ENDIF
#endif
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
      9999 CONTINUE
        CALL substop(caller,t0,info)
        RETURN
      END SUBROUTINE DTYPE(ppm_rc_ghost_copy)