                          nsizet=2*nsize
                          ALLOCATE(energyt(nsizet),STAT=ierr(1))
                          energyt(1:nsize)=energya(1:nsize)
                          CALL MOVE_ALLOC(energyt,energya)

                          ALLOCATE(labelt(nsizet),STAT=ierr(2))
                          labelt(1:nsize)=labela(1:nsize)
                          CALL MOVE_ALLOC(labelt,labela)

                          ALLOCATE(candlabelt(nsizet),STAT=ierr(3))
                          candlabelt(1:nsize)=candlabela(1:nsize)
                          CALL MOVE_ALLOC(candlabelt,candlabela)

                          ALLOCATE(ccandlabelt(nsizet),STAT=ierr(4))
                          ccandlabelt(1:nsize)=ccandlabela(1:nsize)
                          CALL MOVE_ALLOC(ccandlabelt,ccandlabela)
                          ccandlabela(nsize+1:nsizet)=0

                          ALLOCATE(ndaughterst(nsizet),STAT=ierr(5))
                          ndaughterst(1:nsize)=ndaughtersa(1:nsize)
                          CALL MOVE_ALLOC(ndaughterst,ndaughtersa)
                          ndaughtersa(nsize+1:nsizet)=0

                          ALLOCATE(daughterst(2*__DIME,nsizet),STAT=ierr(6))
                          FORALL (ii=1:2*__DIME,jj=1:nsize) daughterst(ii,jj)=daughtersa(ii,jj)
                          CALL MOVE_ALLOC(daughterst,daughtersa)

                          ALLOCATE(nmotherst(nsizet),STAT=ierr(7))
                          nmotherst(1:nsize)=nmothersa(1:nsize)
                          CALL MOVE_ALLOC(nmotherst,nmothersa)
                          nmothersa(nsize+1:nsizet)=0

                          ALLOCATE(acceptedt(nsizet),STAT=ierr(8))
                          acceptedt(1:nsize)=accepteda(1:nsize)
                          CALL MOVE_ALLOC(acceptedt,accepteda)
                          accepteda(nsize+1:nsizet)=-1

                          ALLOCATE(motherst(2*__DIME,nsizet),STAT=ierr(9))
                          FORALL (ii=1:2*__DIME,jj=1:nsize) motherst(ii,jj)=mothersa(ii,jj)
                          CALL MOVE_ALLOC(motherst,mothersa)

                          IF (ANY(ierr.NE.0)) info=ppm_error_fatal
                          or_fail_alloc("Array allocation to the new size is failed", &
                          & ppm_error=ppm_error_fatal)

                          nsize=nsizet

