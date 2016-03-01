      !-------------------------------------------------------------------------
      !  Subroutine   :                  ppm_rc_contour_propagation
      !-------------------------------------------------------------------------
      ! Copyright (c) 2016 MOSAIC Group (MPI-CBG Dresden)
      !
      !
      ! This file is part of the PPM_RC program.
      !
      ! PPM_RC is free software: you can redistribute it and/or modify
      ! it under the terms of the GNU Lesser General Public License
      ! as published by the Free Software Foundation, either
      ! version 3 of the License, or (at your option) any later
      ! version.
      !
      ! PPM_RC is distributed in the hope that it will be useful,
      ! but WITHOUT ANY WARRANTY; without even the implied warranty of
      ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
      ! GNU General Public License for more details.
      !
      ! You should have received a copy of the GNU General Public License
      ! and the GNU Lesser General Public License along with PPM_RC. If not,
      ! see <http://www.gnu.org/licenses/>.
      !-------------------------------------------------------------------------
      !  MOSAIC Group
      !  Max Planck Institute of Molecular Cell Biology and Genetics
      !  Pfotenhauerstr. 108, 01307 Dresden, Germany
      !
      !  Author           - y.afshar           June   2014
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Subroutine   :                  ppm_rc_contour_propagation
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Loops through all the particles,
      !
      !  Input        :
      !
      !  Input/output :
      !
      !  Output       : info       (I) return status. 0 on success.
      !
      !  Routines     :
      !
      !  Remarks      : A daughter of a particle is any particle in the
      !                 4-6-neighborhood that belongs to a different region.
      !
      !
      !  References   :
      !-------------------------------------------------------------------------
      SUBROUTINE DTYPE(ppm_rc_contour_propagation)(info)

        !-------------------------------------------------------------------------
        !  Modules
        !-------------------------------------------------------------------------
        USE ppm_module_mpi
        USE ppm_module_util_qsort, ONLY : ppm_util_qsort
        IMPLICIT NONE

        !-------------------------------------------------------------------------
        !  Includes
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        INTEGER, INTENT(   OUT)  :: info

        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
        TYPE(ppm_rc_list), POINTER :: seed

        CLASS(ppm_t_subpatch_), POINTER :: sbpitr

        REAL(ppm_kind_double) :: t0

        INTEGER, DIMENSION(:),     POINTER ::  subgraph
        INTEGER, DIMENSION(:),     POINTER ::  subgraph1
        INTEGER, DIMENSION(:),     POINTER :: ssubgraph
        INTEGER, DIMENSION(:),     POINTER :: gi
        INTEGER, DIMENSION(1)              :: ldu
        INTEGER, DIMENSION(9)              :: ierr
        INTEGER                            :: numpart
        INTEGER                            :: ipart,i,j,k,l,m
        INTEGER                            :: nsize,nsizea
        INTEGER                            :: nsizebuf,nsizesub
        INTEGER                            :: missed
        INTEGER                            :: ip,iopt,nd,d,ipatch
#ifdef __MPI
        INTEGER                            :: partOffset
        INTEGER, DIMENSION(:), ALLOCATABLE :: git
        INTEGER, DIMENSION(:), ALLOCATABLE :: displ,counts
        INTEGER                            :: tag,ind,dsp
        INTEGER                            :: ns,nr
        INTEGER                            :: ns1,nr1
        INTEGER                            :: ns2,nr2
        INTEGER                            :: ns3,nr3
        INTEGER                            :: root
        INTEGER, DIMENSION(2)              :: request
        INTEGER, DIMENSION(:), ALLOCATABLE :: requestsr
#endif
        INTEGER                            :: iter_id

        CHARACTER(LEN=ppm_char) :: caller="ppm_rc_contour_propagation"

        LOGICAL :: clg
        LOGICAL :: missedpart
        LOGICAL :: lgraph
        !-------------------------------------------------------------------------
        !  Initialize
        !-------------------------------------------------------------------------
        CALL substart(caller,t0,info)

        !-------------------------------------------------------------------------
        !  assigning the global index to the newly created particles
        !-------------------------------------------------------------------------
#ifdef __MPI
        CALL MPI_Iallreduce(Npart,totNpart,1,MPI_INTEGER,MPI_SUM,comm,request(1),info)
        or_fail_MPI("MPI_Iallreduce")

        partOffset=0
        CALL MPI_Iexscan(NpartNew,partOffset,1,MPI_INTEGER,MPI_SUM,comm,request(2),info)
        or_fail_MPI("MPI_Exscan")
#endif

        !-------------------------------------------------------------------------
        ! Loop through all the patches and remove the ghost particles
        ! from the InnerContourContainer
        !-------------------------------------------------------------------------
        sbpitr => mesh%subpatch%begin()
        ipatch=1
        DO WHILE (ASSOCIATED(sbpitr))
           IF (InnerContourContainer(ipatch)%nb.GT.0) THEN
              iter_id=InnerContourContainer(ipatch)%max_id
              InnerContourContainer(ipatch)%iter_id=iter_id
              DO WHILE (iter_id.GT.Npart)
                 seed => InnerContourContainer(ipatch)%vec(iter_id)%t

                 !get rid of all the ghost particles in InnerContourContainer
                 CALL InnerContourContainer(ipatch)%remove(info)
                 or_fail("InnerContourContainer(ipatch)%remove")

                 DEALLOCATE(seed,STAT=info)
                 or_fail_dealloc("Failed to deallocate seed")

                 iter_id=InnerContourContainer(ipatch)%iter_id
              ENDDO
           ENDIF

           sbpitr => mesh%subpatch%next()
           ipatch=ipatch+1
        ENDDO

        !-------------------------------------------------------------------------
        !  Initialize subgraph pointer to NULL
        !-------------------------------------------------------------------------
        NULLIFY(subgraph)

        !Total number of particles on this processor
        nsize=Mpart+NpartNew

        !If there is any graph
        lgraph=.FALSE.
        DO i=1,nsize
           IF (ndaughtersa(i).GT.0) THEN
              lgraph=.TRUE.
              EXIT
           ENDIF
        ENDDO

        IF (lgraph) THEN
#ifdef __MPI
           !-------------------------------------------------------------------------
           ! ssigning global index
           !-------------------------------------------------------------------------
           NULLIFY(gi)

           CALL Part%get(Part%gi,gi,info,read_only=.TRUE.,with_ghosts=.TRUE.)
           or_fail("get particles global index pointer")

           ALLOCATE(git(nsize),STAT=info)
           or_fail_alloc("git")

           FORALL (i=1:Mpart) git(i)=gi(i)

           !does not change any flag
           CALL Part%set(Part%gi,gi,info)
           or_fail("set particles global index pointer")
#endif

           !-------------------------------------------------------------------------
           ! Construct undirectde graph from candidate list of
           ! all particles we consider moving to another region
           !-------------------------------------------------------------------------
           CALL ppm_rc_create_graph(nsize,daughtersa,ndaughtersa,subgraph,info)
           or_fail("ppm_rc_create_graph")

#ifdef __MPI
           ALLOCATE(bufa(32),STAT=info)
           or_fail_alloc("bufa")

           nsizebuf=32

           nsize=SIZE(subgraph,DIM=1)

           l=0
           j=2
           DO i=2,nsize
              IF (subgraph(i).NE.-1) CYCLE

              !--------------------------------------------------------------------
              ! maximal connected subgraph
              !--------------------------------------------------------------------
              !Upperbound of the maximal connected subgraph
              k=i-1
              !Maximal connected subgraph
              ssubgraph => subgraph(j:k)
              !
              nsizesub=SIZE(ssubgraph,DIM=1)
              !Lowerbound of the maximal connected subgraph
              j=i+1

              !there is no graph
              IF (nsizesub.EQ.1) THEN
                 IF (ssubgraph(1).LE.Npart) THEN
                    IF (missedparticles(ssubgraph(1)).EQ.1) THEN
                       ! grow the memory if it is needed
                       IF (l+1.GT.nsizebuf) THEN
                          DO WHILE(l+nsizesub.GT.nsizebuf)
                             nsizebuf=nsizebuf*2
                          ENDDO
                          ALLOCATE(buft(nsizebuf),STAT=info)
                          or_fail_alloc("buft")

                          FORALL (m=1:l) buft(m)=bufa(m)

                          CALL MOVE_ALLOC(buft,bufa)
                       ENDIF

                       l=l+1
                       bufa(l)=ssubgraph(1)
                    ENDIF !(missedpart)
                 ENDIF !(ssubgraph(1).LE.Npart)
                 CYCLE
              ENDIF !(nsizesub.EQ.1)

              !Check to see whether the subgraph will cross the border or not
              IF (ANY(ssubgraph.GT.Npart.AND.ssubgraph.LE.Mpart)) THEN
                 ! grow the memory if it is needed
                 IF (l+nsizesub.GT.nsizebuf) THEN
                    DO WHILE(l+nsizesub.GT.nsizebuf)
                       nsizebuf=nsizebuf*2
                    ENDDO
                    ALLOCATE(buft(nsizebuf),STAT=info)
                    or_fail_alloc("buft")

                    FORALL (m=1:l) buft(m)=bufa(m)

                    CALL MOVE_ALLOC(buft,bufa)
                 ENDIF

                 DO ip=1,nsizesub
                    IF (ssubgraph(ip).LE.Npart) THEN
                       l=l+1
                       bufa(l)=ssubgraph(ip)
                    ELSE IF (ssubgraph(ip).GT.Mpart) THEN
                       l=l+1
                       bufa(l)=ssubgraph(ip)
                    ELSE
                       !This is only for 4,6-connectivity
                       IF (ndaughtersa(ssubgraph(ip)).EQ.1) THEN
                          d=daughtersa(1,ssubgraph(ip))
                          !Who is the daughter
                          ndaughtersa(d)=ndaughtersa(d)+1
                          daughtersa(ndaughtersa(d),d)=ssubgraph(ip)
                          !definitely this daughter component is empty
                          nmothersa(d)=nmothersa(d)+1000
                          !make it obvious which one is carrying extra information
                       ENDIF
                    ENDIF
                 ENDDO !ip=1,nsizesub

                 ssubgraph(1)=-2
              !local on this processor
              ELSE
                 missedpart=.FALSE.

                 DO ip=1,nsizesub
                    IF (ssubgraph(ip).GT.Npart) CYCLE
                    IF (missedparticles(ssubgraph(ip)).EQ.1) THEN
                       missedpart=.TRUE.
                       EXIT
                    ENDIF
                 ENDDO

                 IF (missedpart) THEN
                    ! grow the memory if it is needed
                    IF (l+nsizesub.GT.nsizebuf) THEN
                       DO WHILE(l+nsizesub.GT.nsizebuf)
                          nsizebuf=nsizebuf*2
                       ENDDO
                       ALLOCATE(buft(nsizebuf),STAT=info)
                       or_fail_alloc("buft")

                       FORALL (m=1:l) buft(m)=bufa(m)

                       CALL MOVE_ALLOC(buft,bufa)
                    ENDIF

                    DO ip=1,nsizesub
                       l=l+1
                       bufa(l)=ssubgraph(ip)
                    ENDDO

                    ssubgraph(1)=-2
                 ENDIF !(missedpart)
              ENDIF !(ANY(ssubgraph.GT.Npart.AND.ssubgraph.LE.Mpart))
           ENDDO !i=2,nsize

           !Upperbound of the maximal connected subgraph
           k=nsize
           !Maximal connected subgraph
           ssubgraph => subgraph(j:k)
           !
           nsizesub=SIZE(ssubgraph,DIM=1)

           !Is there is a graph or a single missed particle
           IF (nsizesub.EQ.1) THEN
              IF (ssubgraph(1).LE.Npart) THEN
                 IF (missedparticles(ssubgraph(1)).EQ.1) THEN
                    ! grow the memory if it is needed
                    IF (l+1.GT.nsizebuf) THEN
                       ALLOCATE(buft(l+1),STAT=info)
                       or_fail_alloc("buft")

                       FORALL (m=1:l) buft(m)=bufa(m)

                       CALL MOVE_ALLOC(buft,bufa)
                    ENDIF

                    l=l+1
                    bufa(l)=ssubgraph(1)
                 ENDIF !(missedpart)
              ENDIF !(ssubgraph(1).LE.Npart)
           !Is there is a graph or a single missed particle
           ELSE IF (nsizesub.GT.1) THEN
              !Check to see whether the subgraph will cross the border or not
              IF (ANY(ssubgraph.GT.Npart.AND.ssubgraph.LE.Mpart)) THEN
                 ! grow the memory if it is needed
                 IF (l+nsizesub.GT.nsizebuf) THEN
                    ALLOCATE(buft(l+nsizesub),STAT=info)
                    or_fail_alloc("buft")

                    FORALL (m=1:l) buft(m)=bufa(m)

                    CALL MOVE_ALLOC(buft,bufa)
                 ENDIF

                 DO ip=1,nsizesub
                    IF (ssubgraph(ip).LE.Npart) THEN
                       l=l+1
                       bufa(l)=ssubgraph(ip)
                    ELSE IF (ssubgraph(ip).GT.Mpart) THEN
                       l=l+1
                       bufa(l)=ssubgraph(ip)
                    ELSE
                       IF (ndaughtersa(ssubgraph(ip)).EQ.1) THEN
                          d=daughtersa(1,ssubgraph(ip))
                          !Who is the daughter
                          ndaughtersa(d)=ndaughtersa(d)+1
                          daughtersa(ndaughtersa(d),d)=ssubgraph(ip)
                          !definitely this daughter component is empty
                          nmothersa(d)=nmothersa(d)+1000
                          !make it obvious which one is carrying extra information
                       ENDIF
                    ENDIF
                 ENDDO

                 ssubgraph(1)=-2
              !local on this processor
              ELSE
                 missedpart=.FALSE.

                 DO ip=1,nsizesub
                    IF (ssubgraph(ip).GT.Npart) CYCLE
                    IF (missedparticles(ssubgraph(ip)).EQ.1) THEN
                       missedpart=.TRUE.
                       EXIT
                    ENDIF
                 ENDDO

                 IF (missedpart) THEN
                    ! grow the memory if it is needed
                    IF (l+nsizesub.GT.nsizebuf) THEN
                       ALLOCATE(buft(l+nsizesub),STAT=info)
                       or_fail_alloc("buft")

                       FORALL (m=1:l) buft(m)=bufa(m)

                       CALL MOVE_ALLOC(buft,bufa)
                    ENDIF

                    DO ip=1,nsizesub
                       l=l+1
                       bufa(l)=ssubgraph(ip)
                    ENDDO

                    ssubgraph(1)=-2
                 ENDIF !(missedpart)
              ENDIF !(ANY(ssubgraph.GT.Npart.AND.ssubgraph.LE.Mpart))
           ENDIF !(nsizesub.GT.1)
        ELSE !lgraph=.FALSE. there is no graph
           nsize=0
           l=0

           IF (ANY(missedparticles(1:Npart).EQ.1)) THEN
              missed=COUNT(missedparticles.EQ.1)

              ALLOCATE(bufa(1:missed),git(1:Npart),STAT=info)
              or_fail_alloc("bufa & git")

              DO i=1,Npart
                 IF (missedparticles(i).EQ.1) THEN
                    l=l+1
                    bufa(l)=i
                 ENDIF
              ENDDO !i=1,Npart

              NULLIFY(gi)
              CALL Part%get(Part%gi,gi,info,read_only=.TRUE.)
              or_fail("get particles global index pointer")

              FORALL (i=1:Npart) git(i)=gi(i)

              !does not change any flag
              CALL Part%set(Part%gi,gi,info)
              or_fail("set particles global index pointer")
           ELSE
              ALLOCATE(bufa(1:nsize),git(1:nsize),STAT=info)
              or_fail_alloc("bufa & git")
           ENDIF
#endif
        ENDIF !(lgraph)

#ifdef __MPI
        ALLOCATE(counts(nneighproc),STAT=info)
        or_fail_alloc("counts")

        counts=0

        ! size of the graph that has subgraphs crossed the borders
        nsizea=l

        nsize=nneighproc*2
        ALLOCATE(requestsr(nsize),STAT=info)
        or_fail_alloc("requestsr")

        nr=0
        DO l=1,nneighproc
           ind=ineighproc(l)
           !convert index to rank
           root=ind-1
           IF (rank.LT.root) THEN
              tag=root*(root-1)/2+rank
           ELSE IF (rank.GT.root) THEN
              tag=rank*(rank-1)/2+root
           ENDIF

           nr=nr+1
           CALL MPI_Irecv(counts(l),1,MPI_INTEGER,root,tag,comm,requestsr(nr),info)
           or_fail_MPI("MPI_Irecv")

           nr=nr+1
           CALL MPI_Isend(nsizea,1,MPI_INTEGER,root,tag,comm,requestsr(nr),info)
           or_fail_MPI("MPI_Isend")
        ENDDO !l=1,nneighproc

        !Make sure that we have partOffset & totNpart on this processor
        CALL MPI_Waitall(2,request,MPI_STATUSES_IGNORE,info)
        or_fail_MPI("MPI_Waitall")

        IF (lgraph) THEN
           !Total number of particles on this processor
           nsize=Mpart+NpartNew

           !To make the numbering correct
           partOffset=partOffset+totNpart-Mpart

           FORALL (i=Mpart+1:nsize) git(i)=partOffset+i
        ENDIF

        ALLOCATE(displ(nneighproc),STAT=info)
        or_fail_alloc("displ")

        DEALLOCATE(missedparticles,STAT=info)
        or_fail_dealloc("missedparticles")

        !wait for counts to be avilable on every rank
        CALL MPI_Waitall(nr,requestsr,MPI_STATUSES_IGNORE,info)
        or_fail_MPI("MPI_Waitall")

        DEALLOCATE(requestsr,STAT=info)
        or_fail_dealloc("requestsr")

        !Total number of particles which their information will be gathered
        !on this processor
        numpart=nsizea
        DO i=1,nneighproc
           IF (procflag(i+1).NE.procflag(1)) THEN
              numpart=numpart+counts(i)
           ENDIF
        ENDDO

        IF (nneighproc.GT.0) THEN
           displ(1)=nsizea
           DO i=1,nneighproc-1
              IF (procflag(i+1).NE.procflag(1)) THEN
                 displ(i+1)=displ(i)+counts(i)
              ELSE
                 displ(i+1)=displ(i)
              ENDIF
           ENDDO
        ENDIF !(nneighproc.GT.0)

        !If we are root, we would allocate enough memory to gather the
        !required information for global graph, and we gather information
        !from other processors.
        nsize=COUNT(procflag.NE.procflag(1))
        ns=MERGE(nsize,0,nsizea.GT.0)
        nr=0
        DO i=1,nneighproc
           IF (counts(i).GT.0.AND.procflag(i+1).NE.procflag(1)) THEN
              nr=nr+1
           ENDIF
        ENDDO

        nsize=(ns+nr)
        ALLOCATE(requestsr(nsize),STAT=info)
        or_fail_alloc("requestsr")

        !number of received messages
        nr=0

        ALLOCATE(sndrcvbuf(numpart*(8+2*__DIME)),STAT=info)
        or_fail_alloc("Failed to allocate sndrcvbuf!")

        ALLOCATE(ndaughterst(numpart),STAT=info)
        or_fail_alloc("Failed to allocate ndaughterst!")

        SELECT CASE (procflag(1))
        CASE (1)
           !
           DO l=1,nneighproc
              IF (counts(l).LE.0.OR.procflag(l+1).EQ.1) CYCLE

              !convert index to rank
              root=ineighproc(l)-1
              IF (rank.LT.root) THEN
                 tag=root*(root-1)/2+rank
              ELSE IF (rank.GT.root) THEN
                 tag=rank*(rank-1)/2+root
              ENDIF
              dsp=displ(l)*(8+2*__DIME)+1

              nr=nr+1
              CALL MPI_Irecv(sndrcvbuf(dsp),counts(l)*(8+2*__DIME),MPI_INTEGER,root,tag+1,comm,requestsr(nr),info)
              or_fail_MPI("MPI_Irecv(ltg)")
           ENDDO !l=1,nneighproc

           !Pack the data in sendrecv buffer
           DO i=1,nsizea
              sndrcvbuf(i)=git(bufa(i))
           ENDDO
           l=nsizea
           DO i=1,nsizea
              sndrcvbuf(i+l)=nmothersa(bufa(i))
           ENDDO
           !
           DO i=1,nsizea
              ndaughterst(i)=ndaughtersa(bufa(i))
           ENDDO
           l=2*nsizea
           DO i=1,nsizea
              sndrcvbuf(i+l)=ndaughterst(i)
           ENDDO
           l=3*nsizea-2*__DIME
           DO i=1,nsizea
#if    __DIME == __2D
              k=i*4+l
#elif  __DIME == __3D
              k=i*6+l
#endif
              DO j=1,ndaughterst(i)
                 sndrcvbuf(k+j)=git(daughtersa(j,bufa(i)))
              ENDDO
           ENDDO

           !TODO CHECK
           DEALLOCATE(git,STAT=info)
           or_fail_dealloc("git")

           l=(3+2*__DIME)*nsizea
           DO i=1,nsizea
              sndrcvbuf(i+l)=TRANSFER(energya(bufa(i)),1)
           ENDDO
           l=(4+2*__DIME)*nsizea
           DO i=1,nsizea
              sndrcvbuf(i+l)=labela(bufa(i))
           ENDDO
           l=(5+2*__DIME)*nsizea
           DO i=1,nsizea
              sndrcvbuf(i+l)=candlabela(bufa(i))
           ENDDO

        CASE(2)
           !
           DO l=1,nneighproc
              IF (counts(l).LE.0.OR.procflag(l+1).EQ.2) CYCLE
              !convert index to rank
              root=ineighproc(l)-1
              IF (rank.LT.root) THEN
                 tag=root*(root-1)/2+rank
              ELSE IF (rank.GT.root) THEN
                 tag=rank*(rank-1)/2+root
              ENDIF
              dsp=displ(l)*(8+2*__DIME)+1

              nr=nr+1
              CALL MPI_Irecv(sndrcvbuf(dsp),counts(l)*(8+2*__DIME),MPI_INTEGER,root,tag+1,comm,requestsr(nr),info)
              or_fail_MPI("MPI_Irecv(sndrcvbuf)")
           ENDDO !l=1,nneighproc

           !Pack the data in sendrecv buffer
           DO i=1,nsizea
              sndrcvbuf(i)=git(bufa(i))
           ENDDO
           l=nsizea
           DO i=1,nsizea
              sndrcvbuf(i+l)=nmothersa(bufa(i))
           ENDDO
           !
           DO i=1,nsizea
              ndaughterst(i)=ndaughtersa(bufa(i))
           ENDDO
           l=2*nsizea
           DO i=1,nsizea
              sndrcvbuf(i+l)=ndaughterst(i)
           ENDDO
           l=3*nsizea-2*__DIME
           DO i=1,nsizea
#if    __DIME == __2D
              k=i*4+l
#elif  __DIME == __3D
              k=i*6+l
#endif
              DO j=1,ndaughterst(i)
                 sndrcvbuf(k+j)=git(daughtersa(j,bufa(i)))
              ENDDO
           ENDDO

           !TODO CHECK
           DEALLOCATE(git,STAT=info)
           or_fail_dealloc("git")

           l=(3+2*__DIME)*nsizea
           DO i=1,nsizea
              sndrcvbuf(i+l)=TRANSFER(energya(bufa(i)),1)
           ENDDO
           l=(4+2*__DIME)*nsizea
           DO i=1,nsizea
              sndrcvbuf(i+l)=labela(bufa(i))
           ENDDO
           l=(5+2*__DIME)*nsizea
           DO i=1,nsizea
              sndrcvbuf(i+l)=candlabela(bufa(i))
           ENDDO
           l=(6+2*__DIME)*nsizea
           DO i=1,nsizea
              sndrcvbuf(i+l)=ccandlabela(bufa(i))
           ENDDO
           l=(7+2*__DIME)*nsizea
           DO i=1,nsizea
              sndrcvbuf(i+l)=accepteda(bufa(i))
           ENDDO

           ns=nr
           IF (nsizea.GT.0) THEN
              DO l=1,nneighproc
                 IF (procflag(l+1).EQ.2) CYCLE
                 !convert index to rank
                 root=ineighproc(l)-1
                 IF (rank.LT.root) THEN
                    tag=root*(root-1)/2+rank
                 ELSE IF (rank.GT.root) THEN
                    tag=rank*(rank-1)/2+root
                 ENDIF

                 ns=ns+1
                 CALL MPI_Isend(sndrcvbuf,nsizea*(8+2*__DIME),MPI_INTEGER,root,tag+1,comm,requestsr(ns),info)
                 or_fail_MPI("MPI_Isend(sndrcvbuf)")
              ENDDO !l=1,nneighproc
           ENDIF !(nsizea.GT.0)

        END SELECT







        !Fill the local array on this processor
        ALLOCATE(        ltg(         numpart),STAT=ierr(1))
        ALLOCATE(  nmotherst(         numpart),STAT=ierr(2))
        ALLOCATE( daughterst(2*__DIME,numpart),STAT=ierr(3))
        ALLOCATE(    energyt(         numpart),STAT=ierr(4))
        ALLOCATE(     labelt(         numpart),STAT=ierr(5))
        ALLOCATE( candlabelt(         numpart),STAT=ierr(6))
        ALLOCATE(ccandlabelt(         numpart),STAT=ierr(7))
        ALLOCATE(  acceptedt(         numpart),STAT=ierr(8))
        IF (ANY(ierr(1:8).NE.0)) info=ppm_error_fatal
        or_fail_alloc("energyt,labelt,candlabelt,ccandlabelt,nmotherst,acceptedt,ltg & daughterst")

        !Local on this processor
        DO i=1,nsizea
           ltg(i)=sndrcvbuf(i)
        ENDDO
        l=nsizea
        DO i=1,nsizea
           nmotherst(i)=sndrcvbuf(i+l)
        ENDDO
        l=3*nsizea-2*__DIME
        DO i=1,nsizea
#if    __DIME == __2D
           k=i*4+l
#elif  __DIME == __3D
           k=i*6+l
#endif
           DO j=1,ndaughterst(i)
              k=k+1
              daughterst(j,i)=sndrcvbuf(k)
           ENDDO
        ENDDO
        DO i=1,nsizea
           energyt(i)=energya(bufa(i))
        ENDDO
        l=(4+2*__DIME)*nsizea
        DO i=1,nsizea
           labelt(i)=sndrcvbuf(i+l)
        ENDDO
        l=(5+2*__DIME)*nsizea
        DO i=1,nsizea
           candlabelt(i)=sndrcvbuf(i+l)
        ENDDO
        SELECT CASE (procflag(1))
        CASE (1)
           DO i=1,nsizea
              ccandlabelt(i)=ccandlabela(bufa(i))
           ENDDO
           DO i=1,nsizea
              acceptedt(i)=accepteda(bufa(i))
           ENDDO
        CASE (2)
           l=(6+2*__DIME)*nsizea
           DO i=1,nsizea
              ccandlabelt(i)=sndrcvbuf(i+l)
           ENDDO
           l=(7+2*__DIME)*nsizea
           DO i=1,nsizea
              acceptedt(i)=sndrcvbuf(i+l)
           ENDDO
        END SELECT










        SELECT CASE (procflag(1))
        CASE (1)
           !-------------------------------------------------------------------------
           !  Create the hash table for global to local index
           !-------------------------------------------------------------------------
           CALL htableg%create(numpart,info)
           or_fail("create htableg")

           !Now we should have ltg, nmotherst, ndaughterst, daughterst
           IF (nr.GT.0) THEN
              CALL MPI_Waitall(nr,requestsr(1:nr),MPI_STATUSES_IGNORE,info)
              or_fail_MPI("MPI_Waitall")
           ENDIF
           !
           DO l=1,nneighproc
              IF (counts(l).LE.0.OR.procflag(l+1).EQ.1) CYCLE

              !convert index to rank
              root=ineighproc(l)-1
              IF (rank.LT.root) THEN
                 tag=root*(root-1)/2+rank
              ELSE IF (rank.GT.root) THEN
                 tag=rank*(rank-1)/2+root
              ENDIF

              dsp=displ(l)*(8+2*__DIME)

              k=displ(l)

              DO i=1,counts(l)
                 dsp=dsp+1
                 ltg(k+i)=sndrcvbuf(dsp)
              ENDDO
              DO i=1,counts(l)
                 dsp=dsp+1
                 nmotherst(k+i)=sndrcvbuf(dsp)
              ENDDO
              DO i=1,counts(l)
                 dsp=dsp+1
                 ndaughterst(k+i)=sndrcvbuf(dsp)
              ENDDO
              DO i=1,counts(l)
                 DO j=1,ndaughterst(k+i)
                    daughterst(j,k+i)=sndrcvbuf(dsp+j)
                 ENDDO
                 dsp=dsp+2*__DIME
              ENDDO
              DO i=1,counts(l)
                 dsp=dsp+1
                 energyt(k+i)=TRANSFER(sndrcvbuf(dsp),1._MK)
              ENDDO
              DO i=1,counts(l)
                 dsp=dsp+1
                 labelt(k+i)=sndrcvbuf(dsp)
              ENDDO
              DO i=1,counts(l)
                 dsp=dsp+1
                 candlabelt(k+i)=sndrcvbuf(dsp)
              ENDDO
              DO i=1,counts(l)
                 dsp=dsp+1
                 ccandlabelt(k+i)=sndrcvbuf(dsp)
              ENDDO
              DO i=1,counts(l)
                 dsp=dsp+1
                 acceptedt(k+i)=sndrcvbuf(dsp)
              ENDDO
           ENDDO !l=1,nneighproc

           DO i=1,numpart
              !All the numpart particles are unique real particles
              CALL htableg%insert(ltg(i),i,info)
              or_fail("hash insert")
           ENDDO

           !Converting global index in daughters list to local index
           DO ip=1,numpart
              IF (nmotherst(ip).GT.100) THEN
                 !!!It means some information have been encrypted here, and this
                 !!!particle is the daughter of another particle which was on the
                 !!!ghost layer
                 nmotherst(ip)=nmotherst(ip)-1000
                 !!!correct the number of mothers for this particle
                 missed=ndaughterst(ip)
                 !!!index of the mother which has been encrypted here
                 ndaughterst(ip)=ndaughterst(ip)-1
                 !!!correct number of daughters for this particle

                 d=htableg%search(daughterst(missed,ip))

                 IF (d.NE.htable_null) THEN
                    IF (ndaughterst(d).LT.2*__DIME) THEN
                       !!!find the real mother which was on the ghost
                       ndaughterst(d)=ndaughterst(d)+1

                       IF (d.LT.ip) THEN
                          daughterst(ndaughterst(d),d)=ip
                       ELSE
                          daughterst(ndaughterst(d),d)=ltg(ip)
                       ENDIF
                    ENDIF !(ndaughterst(d).LT.2*__DIME)
                 ENDIF
              ENDIF !(nmotherst(ip).GT.100)
              ind=0
              DO nd=1,ndaughterst(ip)
                 d=htableg%search(daughterst(nd,ip))
                 IF (d.EQ.htable_null) THEN
                    ndaughterst(ip)=ndaughterst(ip)-1
                    CYCLE
                 ENDIF
                 ind=ind+1
                 daughterst(ind,ip)=d
              ENDDO !nd=1,ndaughterst(ip)
           ENDDO !ip=1,numpart

           !-------------------------------------------------------------------------
           !  destroy the hash table
           !-------------------------------------------------------------------------
           CALL htableg%destroy(info)
           or_fail("destroy_htable")

           !-------------------------------------------------------------------------
           !  destroy not necessary memory
           !-------------------------------------------------------------------------
           DEALLOCATE(ltg,STAT=info)
           or_fail_dealloc("ltg")

           !-------------------------------------------------------------------------
           !  Create graph on every node
           !-------------------------------------------------------------------------
           NULLIFY(subgraph1)

           CALL ppm_rc_create_graph(numpart,daughterst,ndaughterst,subgraph1,info)
           or_fail("ppm_rc_create_graph")

           !TO CHECK
           !This part will be done on every rank
           nsize=MERGE(SIZE(subgraph1,DIM=1),0,ASSOCIATED(subgraph1))
           l=0
           j=2
           subgraph1_loop_b : DO i=2,nsize
              IF (subgraph1(i).NE.-1) CYCLE subgraph1_loop_b
              !Upperbound of the maximal connected subgraph1
              k=i-1
              !Maximal connected subgraph1
              ssubgraph => subgraph1(j:k)
              !
              nsizesub=SIZE(ssubgraph,DIM=1)
              !Lowerbound of the maximal connected subgraph1
              j=i+1

              !There is no graph
              IF (nsizesub.EQ.1) CYCLE subgraph1_loop_b

              CALL ppm_util_qsort(energyt(ssubgraph),energyrank,info)
              or_fail("ppm_util_qsort")

              maximal_connected_b : DO ip=1,nsizesub
                 ipart=ssubgraph(energyrank(ip))

                 IF (energyt(ipart).GE.zero) EXIT maximal_connected_b
                 !It's already been rejected, no need to chceck further
                 IF (acceptedt(ipart).EQ.0) CYCLE maximal_connected_b

                 !C1: if p is a child its referncecount is greater equal than one
                 IF (nmotherst(ipart).GT.0) THEN
                    IF (ccandlabelt(ipart).LT.1) THEN
                       ! C1 is not TRUE
                       acceptedt(ipart)=0
                       CYCLE maximal_connected_b
                    ENDIF
                 ENDIF
                 !C1 is TRUE

                 !C2: if p is a parent all of its children which are accepted
                 ! should should have referncecount greater than 1
                 IF (ndaughterst(ipart).GT.0) THEN
                    c2_b: DO nd=1,ndaughterst(ipart)
                       d=daughterst(nd,ipart)

                       IF (acceptedt(d).EQ.1) THEN
                          IF (ccandlabelt(d).LE.1) THEN
                             !C2 is not TRUE
                             acceptedt(ipart)=0
                             CYCLE maximal_connected_b
                          ENDIF
                       ENDIF
                    ENDDO c2_b
                    !C2 is TRUE

                    clg=.FALSE.
                    !C3: if p is a parent at least one of its children
                    ! is not yet accepted or has a candlabel equal to p
                    c3_b: DO nd=1,ndaughterst(ipart)
                       d=daughterst(nd,ipart)
                       IF (acceptedt(d).EQ.0.OR.candlabelt(d).NE.labelt(ipart)) THEN
                          clg=.TRUE.
                          EXIT c3_b
                       ENDIF
                    ENDDO c3_b

                    IF (.NOT.clg) THEN
                       !C3 is not TRUE
                       acceptedt(ipart)=0
                       CYCLE maximal_connected_b
                    ENDIF

                    !C1, C2, and C3 are TRUE
                    DO nd=1,ndaughterst(ipart)
                       d=daughterst(nd,ipart)
                       IF (candlabelt(d).EQ.labelt(ipart)) THEN
                          ccandlabelt(d)=ccandlabelt(d)-1
                       ENDIF
                    ENDDO
                 ENDIF !(ndaughterst(ipart).GT.0)
              ENDDO maximal_connected_b
           ENDDO subgraph1_loop_b

           IF (nsize.GT.0) THEN
              !Upperbound of the maximal connected subgraph1
              k=nsize
              !Maximal connected subgraph1
              ssubgraph => subgraph1(j:k)
              !
              nsizesub=SIZE(ssubgraph,DIM=1)

              !Is there any graph?
              IF (nsizesub.GT.1) THEN
                 CALL ppm_util_qsort(energyt(ssubgraph),energyrank,info)
                 or_fail("ppm_util_qsort")

                 maximal_connected_o_b : DO ip=1,nsizesub
                    ipart=ssubgraph(energyrank(ip))

                    IF (energyt(ipart).GE.zero) EXIT maximal_connected_o_b
                    !It's already been rejected, no need to chceck further
                    IF (acceptedt(ipart).EQ.0) CYCLE maximal_connected_o_b

                    !C1: if p is a child its referncecount is greater or equal than one
                    IF (nmotherst(ipart).GT.0) THEN
                       IF (ccandlabelt(ipart).LT.1) THEN
                          ! C1 is not TRUE
                          acceptedt(ipart)=0
                          CYCLE maximal_connected_o_b
                       ENDIF
                    ENDIF
                    !C1 is TRUE

                    !C2: if p is a parent all of its children which are accepted
                    ! should should have referncecount greater than 1
                    IF (ndaughterst(ipart).GT.0) THEN
                       c2_o_b: DO nd=1,ndaughterst(ipart)
                          d=daughterst(nd,ipart)

                          IF (acceptedt(d).EQ.1) THEN
                             IF (ccandlabelt(d).LE.1) THEN
                                !C2 is not TRUE
                                acceptedt(ipart)=0
                                CYCLE maximal_connected_o_b
                             ENDIF
                          ENDIF
                       ENDDO c2_o_b
                       !C2 is TRUE

                       clg=.FALSE.
                       !C3: if p is a parent at least one of its children
                       ! is not yet accepted or has a candlabel equal to p
                       c3_o_b: DO nd=1,ndaughterst(ipart)
                          d=daughterst(nd,ipart)
                          IF (acceptedt(d).EQ.0.OR.candlabelt(d).NE.labelt(ipart)) THEN
                             clg=.TRUE.
                             EXIT c3_o_b
                          ENDIF
                       ENDDO c3_o_b

                       IF (.NOT.clg) THEN
                          !C3 is not TRUE
                          acceptedt(ipart)=0
                          CYCLE maximal_connected_o_b
                       ENDIF

                       !C1, C2, and C3 are TRUE
                       DO nd=1,ndaughterst(ipart)
                          d=daughterst(nd,ipart)
                          IF (candlabelt(d).EQ.labelt(ipart)) THEN
                             ccandlabelt(d)=ccandlabelt(d)-1
                          ENDIF
                       ENDDO
                    ENDIF !(ndaughterst(ipart).GT.0)
                 ENDDO maximal_connected_o_b
              ENDIF !(nsizesub.GT.1)
           ENDIF !nsize.GT.0

           l=(6+2*__DIME)*nsizea
           DO i=1,nsizea
              sndrcvbuf(i+l)=ccandlabelt(i)
           ENDDO
           l=(7+2*__DIME)*nsizea
           DO i=1,nsizea
              sndrcvbuf(i+l)=acceptedt(i)
           ENDDO

           ns=nr
           IF (nsizea.GT.0) THEN
              DO l=1,nneighproc
                 IF (procflag(l+1).EQ.1) CYCLE

                 !convert index to rank
                 root=ineighproc(l)-1
                 IF (rank.LT.root) THEN
                    tag=root*(root-1)/2+rank
                 ELSE IF (rank.GT.root) THEN
                    tag=rank*(rank-1)/2+root
                 ENDIF

                 ns=ns+1
                 CALL MPI_Isend(sndrcvbuf,nsizea*(8+2*__DIME),MPI_INTEGER,root,tag+1,comm,requestsr(ns),info)
                 or_fail_MPI("MPI_Isend(sndrcvbuf)")
              ENDDO !l=1,nneighproc
           ENDIF !(nsizea.GT.0)

           DO i=1,nsizea
              accepteda(bufa(i))=acceptedt(i)
           ENDDO

           CALL ppm_rc_destroy_graph(subgraph1,info)
           or_fail("ppm_rc_destroy_graph")

           NULLIFY(subgraph1)

           iopt=ppm_param_dealloc
           CALL ppm_alloc(energyrank,ldu,iopt,info)
           or_fail_dealloc("energyrank")

           NULLIFY(energyrank)

           DEALLOCATE(nmotherst,STAT=ierr(1))
           DEALLOCATE(ndaughterst,STAT=ierr(2))
           DEALLOCATE(daughterst,STAT=ierr(3))
           DEALLOCATE(energyt,STAT=ierr(4))
           DEALLOCATE(labelt,STAT=ierr(5))
           DEALLOCATE(candlabelt,STAT=ierr(6))
           DEALLOCATE(ccandlabelt,STAT=ierr(7))
           DEALLOCATE(acceptedt,STAT=ierr(8))
           IF (ANY(ierr(1:8).NE.0)) info=ppm_error_fatal
           or_fail_dealloc("nmotherst,ndaughterst,daughterst,energyt,labelt,candlabelt,ccandlabelt & acceptedt")

        ENd SELECT
#endif

        !Now we do the local part on every processor
        !to take advantage of computation over communication
        !If there is any graph
        IF (lgraph) THEN
           nsize=SIZE(subgraph,DIM=1)

           l=0
           j=2
           subgraph_loop : DO i=2,nsize
              IF (subgraph(i).NE.-1) CYCLE subgraph_loop

              !--------------------------------------------------------------------
              ! maximal connected subgraph
              !--------------------------------------------------------------------
              !Upperbound of the maximal connected subgraph
              k=i-1
              !Maximal connected subgraph
              ssubgraph => subgraph(j:k)
              !
              nsizesub=SIZE(ssubgraph,DIM=1)
              !Lowerbound of the maximal connected subgraph
              j=i+1

              !there is no graph
              IF (nsizesub.EQ.1) CYCLE subgraph_loop

              !If it crosses the border do nothing, as it has been mapped to
              !the root processor
              IF (ssubgraph(1).EQ.-2) CYCLE subgraph_loop

              CALL ppm_util_qsort(energya(ssubgraph),energyrank,info)
              or_fail("ppm_util_qsort")

              maximal_connected : DO ip=1,nsizesub
                 ipart=ssubgraph(energyrank(ip))

                 ! this is the sorted array based on energy so
                 ! if any array is getting bigger than zero the rest is so
                 IF (energya(ipart).GE.zero) EXIT maximal_connected
                 !It's already been rejected, no need to chceck further
                 IF (accepteda(ipart).EQ.0) CYCLE maximal_connected

                 !C1: if p is a child its referncecount is greater equal than one
                 IF (nmothersa(ipart).GT.0) THEN
                    IF (ccandlabela(ipart).LT.1) THEN
                       ! C1 is not TRUE
                       accepteda(ipart)=0
                       CYCLE maximal_connected
                    ENDIF
                 ENDIF
                 !C1 is TRUE

                 !C2: if p is a parent all of its children which are accepted
                 ! should have referncecount greater than 1
                 IF (ndaughtersa(ipart).GT.0) THEN
                    c2: DO nd=1,ndaughtersa(ipart)
                       d=daughtersa(nd,ipart)
                       IF (accepteda(d).EQ.1) THEN
                          IF (ccandlabela(d).LE.1) THEN
                          !C2 is not TRUE
                             accepteda(ipart)=0
                             CYCLE maximal_connected
                          ENDIF
                       ENDIF
                    ENDDO c2
                    !C2 is TRUE

                    clg=.FALSE.
                    !C3: if p is a parent at least one of its children is not
                    !yet accepted or has a candlabel which is not equal to p
                    c3: DO nd=1,ndaughtersa(ipart)
                       d=daughtersa(nd,ipart)
                       IF (accepteda(d).EQ.0.OR.candlabela(d).NE.labela(ipart)) THEN
                          clg=.TRUE.
                          EXIT c3
                       ENDIF
                    ENDDO c3

                    IF (.NOT.clg) THEN
                       !C3 is not TRUE
                       accepteda(ipart)=0
                       CYCLE maximal_connected
                    ENDIF

                    !C1, C2, and C3 are TRUE
                    DO nd=1,ndaughtersa(ipart)
                       d=daughtersa(nd,ipart)
                       IF (candlabela(d).EQ.labela(ipart)) THEN
                          ccandlabela(d)=ccandlabela(d)-1
                       ENDIF
                    ENDDO
                 ENDIF ! (ndaughtersa(ipart).GT.0)
              ENDDO maximal_connected
           ENDDO subgraph_loop

           !Upperbound of the maximal connected subgraph
           k=nsize
           !Maximal connected subgraph
           ssubgraph => subgraph(j:k)
           !
           nsizesub=SIZE(ssubgraph,DIM=1)

           !Is there is a graph
           IF (nsizesub.GT.1) THEN
              !If it crosses the border do nothing, it has been mapped to the master processor
              !local on this processor
              IF (ssubgraph(1).NE.-2) THEN
                 CALL ppm_util_qsort(energya(ssubgraph),energyrank,info)
                 or_fail("ppm_util_qsort")

                 maximal_connected_o : DO ip=1,nsizesub
                    ipart=ssubgraph(energyrank(ip))

                    IF (energya(ipart).GE.zero) EXIT maximal_connected_o
                    !It's already been rejected, no need to chceck further
                    IF (accepteda(ipart).EQ.0) CYCLE maximal_connected_o

                    !C1: if p is a child its referncecount is greater equal than one
                    IF (nmothersa(ipart).GT.0) THEN
                       IF (ccandlabela(ipart).LT.1) THEN
                          ! C1 is not TRUE
                          accepteda(ipart)=0
                          CYCLE maximal_connected_o
                       ENDIF
                    ENDIF
                    !C1 is TRUE

                    !C2: if p is a parent all of its children which are accepted
                    ! should should have referncecount greater than 1
                    IF (ndaughtersa(ipart).GT.0) THEN
                       c2_o: DO nd=1,ndaughtersa(ipart)
                          d=daughtersa(nd,ipart)
                          IF (accepteda(d).EQ.1) THEN
                             IF (ccandlabela(d).LE.1) THEN
                                !C2 is not TRUE
                                accepteda(ipart)=0
                                CYCLE maximal_connected_o
                             ENDIF
                          ENDIF
                       ENDDO c2_o
                       !C2 is TRUE

                       clg=.FALSE.
                       !C3: if p is a parent at least one of its children
                       ! is not yet accepted or has a candlabel equal to p
                       c3_o: DO nd=1,ndaughtersa(ipart)
                          d=daughtersa(nd,ipart)
                          IF (accepteda(d).EQ.0.OR.candlabela(d).NE.labela(ipart)) THEN
                             clg=.TRUE.
                             EXIT c3_o
                          ENDIF
                       ENDDO c3_o

                       IF (.NOT.clg) THEN
                          !C3 is not TRUE
                          accepteda(ipart)=0
                          CYCLE maximal_connected_o
                       ENDIF

                       !C1, C2, and C3 are TRUE
                       DO nd=1,ndaughtersa(ipart)
                          d=daughtersa(nd,ipart)
                          IF (candlabela(d).EQ.labela(ipart)) THEN
                             ccandlabela(d)=ccandlabela(d)-1
                          ENDIF
                       ENDDO
                    ENDIF !(ndaughtersa(ipart).GT.0)
                 ENDDO maximal_connected_o
              ENDIF !(ssubgraph(1).NE.-2))
           ENDIF !(nsizesub.GT.1)
        ENDIF !(lgraph)

        DEALLOCATE(nmothersa,STAT=ierr(1))
        DEALLOCATE(ndaughtersa,STAT=ierr(2))
        DEALLOCATE(daughtersa,STAT=ierr(3))
        DEALLOCATE(ccandlabela,STAT=ierr(4))
        IF (ANY(ierr(1:4).NE.0)) info=ppm_error_fatal
        or_fail_dealloc("nmothersa,ndaughtersa,daughtersa & ccandlabela")

#ifdef __MPI
        SELECT CASE (procflag(1))
        CASE (2)
           !-------------------------------------------------------------------------
           !  Create the hash table for global to local index
           !-------------------------------------------------------------------------
           CALL htableg%create(numpart,info)
           or_fail("create htableg")

           !Now we should have finished receiving data
           IF (nr.GT.0) THEN
              CALL MPI_Waitall(nr,requestsr(1:nr),MPI_STATUSES_IGNORE,info)
              or_fail_MPI("MPI_Waitall")
           ENDIF

           !
           DO l=1,nneighproc
              IF (counts(l).LE.0.OR.procflag(l+1).EQ.2) CYCLE
              !convert index to rank
              root=ineighproc(l)-1
              IF (rank.LT.root) THEN
                 tag=root*(root-1)/2+rank
              ELSE IF (rank.GT.root) THEN
                 tag=rank*(rank-1)/2+root
              ENDIF

              dsp=displ(l)*(8+2*__DIME)

              k=displ(l)

              DO i=1,counts(l)
                 dsp=dsp+1
                 ltg(k+i)=sndrcvbuf(dsp)
              ENDDO
              DO i=1,counts(l)
                 dsp=dsp+1
                 nmotherst(k+i)=sndrcvbuf(dsp)
              ENDDO
              DO i=1,counts(l)
                 dsp=dsp+1
                 ndaughterst(k+i)=sndrcvbuf(dsp)
              ENDDO
              DO i=1,counts(l)
                 DO j=1,ndaughterst(k+i)
                    daughterst(j,k+i)=sndrcvbuf(dsp+j)
                 ENDDO
                 dsp=dsp+2*__DIME
              ENDDO
              DO i=1,counts(l)
                 dsp=dsp+1
                 energyt(k+i)=TRANSFER(sndrcvbuf(dsp),1._MK)
              ENDDO
              DO i=1,counts(l)
                 dsp=dsp+1
                 labelt(k+i)=sndrcvbuf(dsp)
              ENDDO
              DO i=1,counts(l)
                 dsp=dsp+1
                 candlabelt(k+i)=sndrcvbuf(dsp)
              ENDDO
              DO i=1,counts(l)
                 dsp=dsp+1
                 ccandlabelt(k+i)=sndrcvbuf(dsp)
              ENDDO
              DO i=1,counts(l)
                 dsp=dsp+1
                 acceptedt(k+i)=sndrcvbuf(dsp)
              ENDDO
           ENDDO !l=1,nneighproc

           DO i=1,numpart
              !All the numpart particles are unique real particles
              CALL htableg%insert(ltg(i),i,info)
              or_fail("hash insert")
           ENDDO

           !Converting global index in daughters list to local index
           DO ip=1,numpart
              IF (nmotherst(ip).GT.100) THEN
                 !!!It means some information have been encrypted here, and this
                 !!!particle is the daughter of another particle which was on the
                 !!!ghost layer
                 nmotherst(ip)=nmotherst(ip)-1000
                 !!!correct the number of mothers for this particle
                 missed=ndaughterst(ip)
                 !!!index of the mother which has been encrypted here
                 ndaughterst(ip)=ndaughterst(ip)-1
                 !!!correct number of daughters for this particle

                 d=htableg%search(daughterst(missed,ip))

                 IF (d.NE.htable_null) THEN
                    IF (ndaughterst(d).LT.2*__DIME) THEN
                       !!!find the real mother which was on the ghost
                       ndaughterst(d)=ndaughterst(d)+1

                       IF (d.LT.ip) THEN
                          daughterst(ndaughterst(d),d)=ip
                       ELSE
                          daughterst(ndaughterst(d),d)=ltg(ip)
                       ENDIF
                    ENDIF !(ndaughterst(d).LT.2*__DIME)
                 ENDIF
              ENDIF !(nmotherst(ip).GT.100)
              ind=0
              DO nd=1,ndaughterst(ip)
                 d=htableg%search(daughterst(nd,ip))
                 IF (d.EQ.htable_null) THEN
                    ndaughterst(ip)=ndaughterst(ip)-1
                    CYCLE
                 ENDIF
                 ind=ind+1
                 daughterst(ind,ip)=d
              ENDDO !nd=1,ndaughterst(ip)
           ENDDO !ip=1,numpart

           !-------------------------------------------------------------------------
           ! free memory which I do not need
           !-------------------------------------------------------------------------
           DEALLOCATE(ltg,STAT=info)
           or_fail_dealloc("ltg")

           !-------------------------------------------------------------------------
           ! destroy the hash table
           !-------------------------------------------------------------------------
           CALL htableg%destroy(info)
           or_fail("destroy_htable")

           !-------------------------------------------------------------------------
           !  Create graph on every node
           !-------------------------------------------------------------------------
           CALL ppm_rc_destroy_graph(subgraph,info)
           or_fail("ppm_rc_destroy_graph")

           NULLIFY(subgraph)

           CALL ppm_rc_create_graph(numpart,daughterst,ndaughterst,subgraph,info)
           or_fail("ppm_rc_create_graph")

           !TO CHECK
           !This part will be done on every rank
           nsize=MERGE(SIZE(subgraph,DIM=1),0,ASSOCIATED(subgraph))
           l=0
           j=2
           subgraph_loop_b_2 : DO i=2,nsize
              IF (subgraph(i).NE.-1) CYCLE subgraph_loop_b_2
              !Upperbound of the maximal connected subgraph
              k=i-1
              !Maximal connected subgraph
              ssubgraph => subgraph(j:k)
              !
              nsizesub=SIZE(ssubgraph,DIM=1)
              !Lowerbound of the maximal connected subgraph
              j=i+1

              !There is no graph
              IF (nsizesub.EQ.1) CYCLE subgraph_loop_b_2

              CALL ppm_util_qsort(energyt(ssubgraph),energyrank,info)
              or_fail("ppm_util_qsort")

              maximal_connected_b_2 : DO ip=1,nsizesub
                 ipart=ssubgraph(energyrank(ip))

                 IF (energyt(ipart).GE.zero) EXIT maximal_connected_b_2
                 !It's already been rejected, no need to chceck further
                 IF (acceptedt(ipart).EQ.0) CYCLE maximal_connected_b_2

                 !C1: if p is a child its referncecount is greater equal than one
                 IF (nmotherst(ipart).GT.0) THEN
                    IF (ccandlabelt(ipart).LT.1) THEN
                       ! C1 is not TRUE
                       acceptedt(ipart)=0
                       CYCLE maximal_connected_b_2
                    ENDIF
                 ENDIF
                 !C1 is TRUE

                 !C2: if p is a parent all of its children which are accepted
                 ! should should have referncecount greater than 1
                 IF (ndaughterst(ipart).GT.0) THEN
                    c2_b_2: DO nd=1,ndaughterst(ipart)
                       d=daughterst(nd,ipart)

                       IF (acceptedt(d).EQ.1) THEN
                          IF (ccandlabelt(d).LE.1) THEN
                             !C2 is not TRUE
                             acceptedt(ipart)=0
                             CYCLE maximal_connected_b_2
                          ENDIF
                       ENDIF
                    ENDDO c2_b_2
                    !C2 is TRUE

                    clg=.FALSE.
                    !C3: if p is a parent at least one of its children
                    ! is not yet accepted or has a candlabel equal to p
                    c3_b_2: DO nd=1,ndaughterst(ipart)
                       d=daughterst(nd,ipart)
                       IF (acceptedt(d).EQ.0.OR.candlabelt(d).NE.labelt(ipart)) THEN
                          clg=.TRUE.
                          EXIT c3_b_2
                       ENDIF
                    ENDDO c3_b_2

                    IF (.NOT.clg) THEN
                       !C3 is not TRUE
                       acceptedt(ipart)=0
                       CYCLE maximal_connected_b_2
                    ENDIF

                    !C1, C2, and C3 are TRUE
                    DO nd=1,ndaughterst(ipart)
                       d=daughterst(nd,ipart)
                       IF (candlabelt(d).EQ.labelt(ipart)) THEN
                          ccandlabelt(d)=ccandlabelt(d)-1
                       ENDIF
                    ENDDO
                 ENDIF !(ndaughterst(ipart).GT.0)
              ENDDO maximal_connected_b_2
           ENDDO subgraph_loop_b_2

           IF (nsize.GT.0) THEN
              !Upperbound of the maximal connected subgraph
              k=nsize
              !Maximal connected subgraph
              ssubgraph => subgraph(j:k)
              !
              nsizesub=SIZE(ssubgraph,DIM=1)

              !Is there any graph?
              IF (nsizesub.GT.1) THEN
                 CALL ppm_util_qsort(energyt(ssubgraph),energyrank,info)
                 or_fail("ppm_util_qsort")

                 maximal_connected_o_b_2 : DO ip=1,nsizesub
                    ipart=ssubgraph(energyrank(ip))

                    IF (energyt(ipart).GE.zero) EXIT maximal_connected_o_b_2
                    !It's already been rejected, no need to chceck further
                    IF (acceptedt(ipart).EQ.0) CYCLE maximal_connected_o_b_2

                    !C1: if p is a child its referncecount is greater or equal than one
                    IF (nmotherst(ipart).GT.0) THEN
                       IF (ccandlabelt(ipart).LT.1) THEN
                          ! C1 is not TRUE
                          acceptedt(ipart)=0
                          CYCLE maximal_connected_o_b_2
                       ENDIF
                    ENDIF
                    !C1 is TRUE

                    !C2: if p is a parent all of its children which are accepted
                    ! should should have referncecount greater than 1
                    IF (ndaughterst(ipart).GT.0) THEN
                       c2_o_b_2: DO nd=1,ndaughterst(ipart)
                          d=daughterst(nd,ipart)

                          IF (acceptedt(d).EQ.1) THEN
                             IF (ccandlabelt(d).LE.1) THEN
                                !C2 is not TRUE
                                acceptedt(ipart)=0
                                CYCLE maximal_connected_o_b_2
                             ENDIF
                          ENDIF
                       ENDDO c2_o_b_2
                       !C2 is TRUE

                       clg=.FALSE.
                       !C3: if p is a parent at least one of its children
                       ! is not yet accepted or has a candlabel equal to p
                       c3_o_b_2: DO nd=1,ndaughterst(ipart)
                          d=daughterst(nd,ipart)
                          IF (acceptedt(d).EQ.0.OR.candlabelt(d).NE.labelt(ipart)) THEN
                             clg=.TRUE.
                             EXIT c3_o_b_2
                          ENDIF
                       ENDDO c3_o_b_2

                       IF (.NOT.clg) THEN
                          !C3 is not TRUE
                          acceptedt(ipart)=0
                          CYCLE maximal_connected_o_b_2
                       ENDIF

                       !C1, C2, and C3 are TRUE
                       DO nd=1,ndaughterst(ipart)
                          d=daughterst(nd,ipart)
                          IF (candlabelt(d).EQ.labelt(ipart)) THEN
                             ccandlabelt(d)=ccandlabelt(d)-1
                          ENDIF
                       ENDDO
                    ENDIF !(ndaughterst(ipart).GT.0)
                 ENDDO maximal_connected_o_b_2
              ENDIF !(nsizesub.GT.1)
           ENDIF !nsize.GT.0

           DO i=1,nsizea
              accepteda(bufa(i))=acceptedt(i)
           ENDDO

           DEALLOCATE(nmotherst,STAT=ierr(1))
           DEALLOCATE(ndaughterst,STAT=ierr(2))
           DEALLOCATE(daughterst,STAT=ierr(3))
           DEALLOCATE(energyt,STAT=ierr(4))
           DEALLOCATE(labelt,STAT=ierr(5))
           DEALLOCATE(candlabelt,STAT=ierr(6))
           DEALLOCATE(ccandlabelt,STAT=ierr(7))
           DEALLOCATE(acceptedt,STAT=ierr(8))
           IF (ANY(ierr(1:8).NE.0)) info=ppm_error_fatal
           or_fail_dealloc("nmotherst,ndaughterst,daughterst,energyt,labelt,candlabelt,ccandlabelt & acceptedt")

        END SELECT

        CALL ppm_rc_destroy_graph(subgraph,info)
        or_fail("ppm_rc_destroy_graph")

        NULLIFY(subgraph)

        iopt=ppm_param_dealloc
        CALL ppm_alloc(energyrank,ldu,iopt,info)
        or_fail_dealloc("energyrank")

        NULLIFY(energyrank)

        DEALLOCATE(bufa,counts,displ,STAT=info)
        or_fail_dealloc("bufa,counts & displ",ppm_error=ppm_error_fatal)

        !Wait till the unfinished send are complete
        IF (ns-nr.GT.0) THEN
           CALL MPI_Waitall(ns-nr,requestsr(nr+1:ns),MPI_STATUSES_IGNORE,info)
           or_fail_MPI("MPI_Waitall")
        ENDIF

        DEALLOCATE(sndrcvbuf,STAT=info)
        or_fail_dealloc("Failed to deallocate sndrcvbuf!")

        DEALLOCATE(requestsr,STAT=info)
        or_fail_dealloc("requestsr")
#else
        CALL ppm_rc_destroy_graph(subgraph,info)
        or_fail("ppm_rc_destroy_graph")

        NULLIFY(subgraph)

        iopt=ppm_param_dealloc
        CALL ppm_alloc(energyrank,ldu,iopt,info)
        or_fail_dealloc("energyrank")

        NULLIFY(energyrank)
#endif

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      9999 CONTINUE
        CALL substop(caller,t0,info)
        RETURN
      END SUBROUTINE DTYPE(ppm_rc_contour_propagation)
