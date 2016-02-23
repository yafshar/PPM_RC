                          ALLOCATE(energyt(nsize),SOURCE=energya(1:nsize),STAT=ierr(1))
                          CALL MOVE_ALLOC(energyt,energya)

                          ALLOCATE(labelt(nsize),SOURCE=labela(1:nsize),STAT=ierr(2))
                          CALL MOVE_ALLOC(labelt,labela)

                          ALLOCATE(candlabelt(nsize),SOURCE=candlabela(1:nsize),STAT=ierr(3))
                          CALL MOVE_ALLOC(candlabelt,candlabela)

                          ALLOCATE(ccandlabelt(nsize),SOURCE=ccandlabela(1:nsize),STAT=ierr(4))
                          CALL MOVE_ALLOC(ccandlabelt,ccandlabela)

                          ALLOCATE(ndaughterst(nsize),SOURCE=ndaughtersa(1:nsize),STAT=ierr(5))
                          CALL MOVE_ALLOC(ndaughterst,ndaughtersa)

                          ALLOCATE(daughterst(2*__DIME,nsize),STAT=ierr(6))
                          FORALL (ii=1:2*__DIME,jj=1:nsize) daughterst(ii,jj)=daughtersa(ii,jj)
                          CALL MOVE_ALLOC(daughterst,daughtersa)

                          ALLOCATE(nmotherst(nsize),SOURCE=nmothersa(1:nsize),STAT=ierr(7))
                          CALL MOVE_ALLOC(nmotherst,nmothersa)

                          ALLOCATE(acceptedt(nsize),STAT=ierr(8))
                          DO i=1,nsize
                             IF (accepteda(i).EQ.-1) THEN
                                acceptedt(i)=0
                             ELSE
                                acceptedt(i)=accepteda(i)
                             ENDIF
                          ENDDO
                          CALL MOVE_ALLOC(acceptedt,accepteda)

                          IF (ANY(ierr(1:8).NE.0)) info=ppm_error_fatal
                          or_fail_alloc("Array allocation to the correct size is failed", &
                          & ppm_error=ppm_error_fatal)
