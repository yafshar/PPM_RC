#if   __DIME == __2D
             !-x west
             l=0
             IF (sbpitr%istart(1).NE.1) THEN
                i=0
                DO j=1,Nm(2)
                   l=l+1
                   IF (old_ghost(l).NE.DTYPE(wpl)(i,j)) THEN
                      CALL DTYPE(ppm_rc_ChangeContourParticleLabelToCandidateLabel)   &
                      &   (ipatch,DTYPE(wpl),DTYPE(wpp),(/i,j/),Nm,ABS(old_ghost(l)), &
                      &   ABS(DTYPE(wpl)(i,j)),info)
                      or_fail("ppm_rc_ChangeContourParticleLabelToCandidateLabel")
                   ENDIF
                ENDDO !j
             ENDIF
             !+x east
             IF (sbpitr%iend(1).NE.Ngrid(1)) THEN
                i=Nm(1)+1
                DO j=1,Nm(2)
                   l=l+1
                   IF (old_ghost(l).NE.DTYPE(wpl)(i,j)) THEN
                      CALL DTYPE(ppm_rc_ChangeContourParticleLabelToCandidateLabel)   &
                      &   (ipatch,DTYPE(wpl),DTYPE(wpp),(/i,j/),Nm,ABS(old_ghost(l)), &
                      &   ABS(DTYPE(wpl)(i,j)),info)
                      or_fail("ppm_rc_ChangeContourParticleLabelToCandidateLabel")
                   ENDIF
                ENDDO !j
             ENDIF
             !-y south
             IF (sbpitr%istart(2).NE.1) THEN
                j=0
                DO i=0,Nm(1)+1
                   l=l+1
                   IF (old_ghost(l).NE.DTYPE(wpl)(i,j)) THEN
                      CALL DTYPE(ppm_rc_ChangeContourParticleLabelToCandidateLabel)   &
                      &   (ipatch,DTYPE(wpl),DTYPE(wpp),(/i,j/),Nm,ABS(old_ghost(l)), &
                      &   ABS(DTYPE(wpl)(i,j)),info)
                      or_fail("ppm_rc_ChangeContourParticleLabelToCandidateLabel")
                   ENDIF
                ENDDO !i
             ENDIF
             !+y north
             IF (sbpitr%iend(2).NE.Ngrid(2)) THEN
                j=Nm(2)+1
                DO i=0,Nm(1)+1
                   l=l+1
                   IF (old_ghost(l).NE.DTYPE(wpl)(i,j)) THEN
                      CALL DTYPE(ppm_rc_ChangeContourParticleLabelToCandidateLabel)   &
                      &   (ipatch,DTYPE(wpl),DTYPE(wpp),(/i,j/),Nm,ABS(old_ghost(l)), &
                      &   ABS(DTYPE(wpl)(i,j)),info)
                      or_fail("ppm_rc_ChangeContourParticleLabelToCandidateLabel")
                   ENDIF
                ENDDO !j
             ENDIF
#elif __DIME == __3D
             !-x west
             l=0
             IF (sbpitr%istart(1).NE.1) THEN
                i=0
                DO k=1,Nm(3)
                   DO j=1,Nm(2)
                      l=l+1
                      IF (old_ghost(l).NE.DTYPE(wpl)(i,j,k)) THEN
                         CALL DTYPE(ppm_rc_ChangeContourParticleLabelToCandidateLabel)     &
                         &   (ipatch,DTYPE(wpl),DTYPE(wpp),(/i,j,k/),Nm,ABS(old_ghost(l)), &
                         &   ABS(DTYPE(wpl)(i,j,k)),info)
                         or_fail("ppm_rc_ChangeContourParticleLabelToCandidateLabel")
                      ENDIF
                   ENDDO !j
                ENDDO !k
             ENDIF
             !+x east
             IF (sbpitr%iend(1).NE.Ngrid(1)) THEN
                i=Nm(1)+1
                DO k=1,Nm(3)
                   DO j=1,Nm(2)
                      l=l+1
                      IF (old_ghost(l).NE.DTYPE(wpl)(i,j,k)) THEN
                         CALL DTYPE(ppm_rc_ChangeContourParticleLabelToCandidateLabel)     &
                         &   (ipatch,DTYPE(wpl),DTYPE(wpp),(/i,j,k/),Nm,ABS(old_ghost(l)), &
                         &   ABS(DTYPE(wpl)(i,j,k)),info)
                         or_fail("ppm_rc_ChangeContourParticleLabelToCandidateLabel")
                      ENDIF
                   ENDDO !j
                ENDDO !k
             ENDIF
             !-y south
             IF (sbpitr%istart(2).NE.1) THEN
                j=0
                DO k=1,Nm(3)
                   DO i=0,Nm(1)+1
                      l=l+1
                      IF (old_ghost(l).NE.DTYPE(wpl)(i,j,k)) THEN
                         CALL DTYPE(ppm_rc_ChangeContourParticleLabelToCandidateLabel)     &
                         &   (ipatch,DTYPE(wpl),DTYPE(wpp),(/i,j,k/),Nm,ABS(old_ghost(l)), &
                         &   ABS(DTYPE(wpl)(i,j,k)),info)
                         or_fail("ppm_rc_ChangeContourParticleLabelToCandidateLabel")
                      ENDIF
                   ENDDO !i
                ENDDO !k
             ENDIF
             !+y north
             IF (sbpitr%iend(2).NE.Ngrid(2)) THEN
                j=Nm(2)+1
                DO k=1,Nm(3)
                   DO i=0,Nm(1)+1
                      l=l+1
                      IF (old_ghost(l).NE.DTYPE(wpl)(i,j,k)) THEN
                         CALL DTYPE(ppm_rc_ChangeContourParticleLabelToCandidateLabel)     &
                         &   (ipatch,DTYPE(wpl),DTYPE(wpp),(/i,j,k/),Nm,ABS(old_ghost(l)), &
                         &   ABS(DTYPE(wpl)(i,j,k)),info)
                         or_fail("ppm_rc_ChangeContourParticleLabelToCandidateLabel")
                      ENDIF
                   ENDDO !i
                ENDDO !k
             ENDIF
             !-z bottom
             IF (sbpitr%istart(3).NE.1) THEN
                k=0
                DO j=0,Nm(2)+1
                   DO i=0,Nm(1)+1
                      l=l+1
                      IF (old_ghost(l).NE.DTYPE(wpl)(i,j,k)) THEN
                         CALL DTYPE(ppm_rc_ChangeContourParticleLabelToCandidateLabel)     &
                         &   (ipatch,DTYPE(wpl),DTYPE(wpp),(/i,j,k/),Nm,ABS(old_ghost(l)), &
                         &   ABS(DTYPE(wpl)(i,j,k)),info)
                         or_fail("ppm_rc_ChangeContourParticleLabelToCandidateLabel")
                      ENDIF
                   ENDDO !i
                ENDDO !j
             ENDIF
             !+z top
             IF (sbpitr%iend(3).NE.Ngrid(3)) THEN
                k=Nm(3)+1
                DO j=0,Nm(2)+1
                   DO i=0,Nm(1)+1
                      l=l+1
                      IF (old_ghost(l).NE.DTYPE(wpl)(i,j,k)) THEN
                         CALL DTYPE(ppm_rc_ChangeContourParticleLabelToCandidateLabel)     &
                         &   (ipatch,DTYPE(wpl),DTYPE(wpp),(/i,j,k/),Nm,ABS(old_ghost(l)), &
                         &   ABS(DTYPE(wpl)(i,j,k)),info)
                         or_fail("ppm_rc_ChangeContourParticleLabelToCandidateLabel")
                      ENDIF
                   ENDDO !i
                ENDDO !j
             ENDIF
#endif
