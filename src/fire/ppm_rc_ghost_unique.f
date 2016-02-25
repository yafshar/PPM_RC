      !-------------------------------------------------------------------------
      !  Subroutine   :                    ppm_rc_ghost_unique
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
      !  Subroutine   :                    ppm_rc_ghost_unique
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
      SUBROUTINE DTYPE(ppm_rc_ghost_unique)(MeshIn,FieldIn,info)

        !-------------------------------------------------------------------------
        !  Modules
        !-------------------------------------------------------------------------
        USE ppm_module_mpi
        IMPLICIT NONE

        !-------------------------------------------------------------------------
        !  Includes
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        CLASS(ppm_t_equi_mesh_), POINTER      :: MeshIn
        CLASS(ppm_t_field_),     POINTER      :: FieldIn

        INTEGER,                INTENT(  OUT) :: info
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
        TYPE(ppm_t_topo), POINTER :: topo

        CLASS(ppm_t_subpatch_), POINTER :: sbpitr

        TYPE(ppm_rc_c_list), DIMENSION(:), ALLOCATABLE :: ppm_rc_border_seeds

        TYPE(ppm_rc_list), POINTER :: seed
        TYPE(ppm_rc_list), POINTER :: seedlst
        TYPE(ppm_rc_link), POINTER :: seedlnk

        REAL(ppm_kind_double) :: t0

        INTEGER, DIMENSION(:),         POINTER :: Nm
        INTEGER, DIMENSION(:),         POINTER :: seedn
        INTEGER, DIMENSION(:),         POINTER :: bc
        INTEGER, DIMENSION(:),         POINTER :: iistart,jistart
        !!! Lower-left coordinates on the global mesh
        INTEGER, DIMENSION(:),         POINTER :: iiend,jend



#if   __DIME == __2D
        INTEGER, CONTIGUOUS, DIMENSION(:,:),       POINTER :: DTYPE(i_wp) => NULL()
        INTEGER, DIMENSION(:,:),   ALLOCATABLE :: border
        INTEGER, DIMENSION(:,:),   ALLOCATABLE :: unique_border
        INTEGER, DIMENSION(:),     ALLOCATABLE :: buff
#elif __DIME == __3D
        INTEGER, CONTIGUOUS, DIMENSION(:,:,:),     POINTER :: DTYPE(i_wp) => NULL()
        INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: border
        INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: unique_border
        INTEGER, DIMENSION(:,:),   ALLOCATABLE :: buff
#endif
        INTEGER, DIMENSION(:,:),   ALLOCATABLE :: ineighsubs,content
        INTEGER, DIMENSION(:),     ALLOCATABLE :: displ,counts


        INTEGER                                :: i,j,k,l,m,n,c,s,iopt
        INTEGER                                :: isub,jsub,iisub
        INTEGER                                :: max_id,iter_id
        INTEGER                                :: lmin,lmax
        INTEGER, DIMENSION(__DIME)                  :: ld
        INTEGER                                :: maxnneighsubs
        INTEGER                                :: nsubs,nsize,isize
        INTEGER                                :: is2p,js2p,iis2p

        CHARACTER(ppm_char) :: caller="ppm_rc_ghost_copy"

        !-------------------------------------------------------------------------
        !  Initialize
        !-------------------------------------------------------------------------
        CALL substart(caller,t0,info)

        topo => ppm_topo(MeshIn%topoid)%t
        nsubs=topo%nsubs

        maxnneighsubs=0

        CALL MPI_Reduce(topo%nneighsubs,maxnneighsubs,1,MPI_INTEGER,MPI_MAX,0,comm,info)
        or_fail_MPI("MPI_Reduce")

        ld(1)=maxnneighsubs
        ld(2)=nsubs

        ALLOCATE(ineighsubs(ld(1),ld(2)),STAT=info)
        or_fail_alloc("ineighsubs")

        ineighsubs=0

        ALLOCATE(displ(nsubs),counts(nsubs),STAT=info)
        or_fail_alloc("displ, counts")

        DO isub=1,nsubs
           displ(isub)=(isub-1)*maxnneighsubs
        ENDDO
        counts=maxnneighsubs

        CALL MPI_Gatherv(topo%ineighsubs,topo%nneighsubs,MPI_INTEGER, &
        &    ineighsubs,counts,displ,MPI_INTEGER,0,comm,info)
        or_fail_MPI("MPI_Gatherv")

        nsize=0
        DO isub=1,nsubs
           SELECT TYPE(p => MeshIn%subpatch_by_sub(isub)%vec(1)%t)
           TYPE IS (ppm_t_subpatch)
              nsize=MAX(nsize,MAXVAL(p%nnodes))
           END SELECT
        ENDDO !isub=1,nsubs

        ld(1:__DIME-1)=MERGE(nsize,0,rank.EQ.0)
        ld(__DIME) =MERGE(__DIME*2*nsubs,0,rank.EQ.0)
#if   __DIME == __2D
        ALLOCATE(border(ld(1),ld(2)),STAT=info)
        or_fail_alloc("border")

        ALLOCATE(buff(nsize*4),STAT=info)
        or_fail_alloc("buff")

        DO isub=1,nsubs
           displ(isub)=(isub-1)*4*nsize
        ENDDO
        counts=4*nsize
#elif __DIME == __3D
        ALLOCATE(border(ld(1),ld(2),ld(3)),STAT=info)
        or_fail_alloc("border")

        ALLOCATE(buff(nsize,nsize*6),STAT=info)
        or_fail_alloc("buff")

        DO isub=1,topo%nsubs
           displ(isub)=(isub-1)*6*nsize*nsize
        ENDDO
        counts=6*nsize*nsize
#endif
        buff=0

        sbpitr => MeshIn%subpatch%begin()
        IF (ASSOCIATED(sbpitr)) THEN
           Nm => sbpitr%nnodes
           iistart => sbpitr%istart
           iiend   => sbpitr%iend

           CALL sbpitr%get_field(FieldIn,DTYPE(i_wp),info)
           or_fail("Failed to get field i_wp data.")

#if   __DIME == __2D
           !west
           IF (iistart(1).NE.1) THEN
              DO j=2,Nm(2)-1
                 IF (i_wp_2d(1,j).LE.0) CYCLE
                 IF (i_wp_2d(0,j  ).GT.0.OR. &
                 &   i_wp_2d(0,j-1).GT.0.OR. &
                 &   i_wp_2d(0,j+1).GT.0) buff(j)=i_wp_2d(1,j)
              ENDDO
           ENDIF
           !east
           IF (iiend(1).NE.Ngrid(1)) THEN
              DO j=2,Nm(2)-1
                 IF (i_wp_2d(Nm(1),j).LE.0) CYCLE
                 IF (i_wp_2d(Nm(1)+1,j  ).GT.0.OR. &
                 &   i_wp_2d(Nm(1)+1,j-1).GT.0.OR. &
                 &   i_wp_2d(Nm(1)+1,j+1).GT.0) buff(nsize+j)=i_wp_2d(Nm(1),j)
              ENDDO
!               FORALL (j=2:Nm(2)-1) buff(nsize+j)=i_wp_2d(Nm(1),j)
           ENDIF
           !south
           IF (iistart(2).NE.1) THEN
              l=2*nsize
              DO i=1,Nm(1)
                 IF (i_wp_2d(i,1).LE.0) CYCLE
                 IF (i_wp_2d(i-1,0).GT.0.OR. &
                 &   i_wp_2d(i,  0).GT.0.OR. &
                 &   i_wp_2d(i+1,0).GT.0) buff(l+i)=i_wp_2d(i,1)
              ENDDO
!               FORALL (i=1:Nm(1)) buff(l+i)=i_wp_2d(i,1)
           ENDIF
           !north
           IF (iiend(2).NE.Ngrid(2)) THEN
              l=3*nsize
              DO i=1,Nm(1)
                 IF (i_wp_2d(i,Nm(2)).LE.0) CYCLE
                 IF (i_wp_2d(i-1,Nm(2)+1).GT.0.OR. &
                 &   i_wp_2d(i,  Nm(2)+1).GT.0.OR. &
                 &   i_wp_2d(i+1,Nm(2)+1).GT.0) buff(l+i)=i_wp_2d(i,Nm(2))
              ENDDO
!               FORALL (i=1:Nm(1)) buff(l+i)=i_wp_2d(i,Nm(2))
           ENDIF
#elif __DIME == __3D
           !west
           IF (iistart(1).NE.1) THEN
              DO k=2,Nm(3)-1
                 DO j=2,Nm(2)-1
                    IF (i_wp_3d(1,j,k).LE.0) CYCLE
                    IF (i_wp_3d(0,j-1,k-1).GT.0.OR. &
                    &   i_wp_3d(0,j  ,k-1).GT.0.OR. &
                    &   i_wp_3d(0,j+1,k-1).GT.0.OR. &
                    &   i_wp_3d(0,j-1,k  ).GT.0.OR. &
                    &   i_wp_3d(0,j  ,k  ).GT.0.OR. &
                    &   i_wp_3d(0,j+1,k  ).GT.0.OR. &
                    &   i_wp_3d(0,j-1,k+1).GT.0.OR. &
                    &   i_wp_3d(0,j  ,k+1).GT.0.OR. &
                    &   i_wp_3d(0,j+1,k+1).GT.0) buff(j,k)=i_wp_3d(1,j,k)
                 ENDDO
              ENDDO
!               FORALL (j=2:Nm(2)-1,k=2:Nm(3)-1) buff(j,k)=i_wp_3d(1,j,k)
           ENDIF
           !east
           IF (iiend(1).NE.Ngrid(1)) THEN
              DO k=2,Nm(3)-1
                 DO j=2,Nm(2)-1
                    IF (i_wp_3d(Nm(1),j,k).LE.0) CYCLE
                    IF (i_wp_3d(Nm(1)+1,j-1,k-1).GT.0.OR. &
                    &   i_wp_3d(Nm(1)+1,j  ,k-1).GT.0.OR. &
                    &   i_wp_3d(Nm(1)+1,j+1,k-1).GT.0.OR. &
                    &   i_wp_3d(Nm(1)+1,j-1,k  ).GT.0.OR. &
                    &   i_wp_3d(Nm(1)+1,j  ,k  ).GT.0.OR. &
                    &   i_wp_3d(Nm(1)+1,j+1,k  ).GT.0.OR. &
                    &   i_wp_3d(Nm(1)+1,j-1,k+1).GT.0.OR. &
                    &   i_wp_3d(Nm(1)+1,j  ,k+1).GT.0.OR. &
                    &   i_wp_3d(Nm(1)+1,j+1,k+1).GT.0) buff(j,nsize+k)=i_wp_3d(Nm(1),j,k)
                 ENDDO
              ENDDO
!               FORALL (j=2:Nm(2)-1,k=2:Nm(3)-1) buff(j,nsize+k)=i_wp_3d(Nm(1),j,k)
           ENDIF
           !south
           IF (iistart(2).NE.1) THEN
              l=2*nsize
              DO k=2,Nm(3)-1
                 DO i=1,Nm(1)
                    IF (i_wp_3d(i,1,k).LE.0) CYCLE
                    IF (i_wp_3d(i-1,0,k-1).GT.0.OR. &
                    &   i_wp_3d(i  ,0,k-1).GT.0.OR. &
                    &   i_wp_3d(i+1,0,k-1).GT.0.OR. &
                    &   i_wp_3d(i-1,0,k  ).GT.0.OR. &
                    &   i_wp_3d(i  ,0,k  ).GT.0.OR. &
                    &   i_wp_3d(i+1,0,k  ).GT.0.OR. &
                    &   i_wp_3d(i-1,0,k+1).GT.0.OR. &
                    &   i_wp_3d(i  ,0,k+1).GT.0.OR. &
                    &   i_wp_3d(i+1,0,k+1).GT.0) buff(i,l+k)=i_wp_3d(i,1,k)
                 ENDDO
              ENDDO
!               FORALL (i=1:Nm(1),k=2:Nm(3)-1) buff(i,l+k)=i_wp_3d(i,1,k)
           ENDIF
           !north
           IF (iiend(2).NE.Ngrid(2)) THEN
              l=3*nsize
              DO k=2,Nm(3)-1
                 DO i=1,Nm(1)
                    IF (i_wp_3d(i,Nm(2),k).LE.0) CYCLE
                    IF (i_wp_3d(i-1,Nm(2)+1,k-1).GT.0.OR. &
                    &   i_wp_3d(i  ,Nm(2)+1,k-1).GT.0.OR. &
                    &   i_wp_3d(i+1,Nm(2)+1,k-1).GT.0.OR. &
                    &   i_wp_3d(i-1,Nm(2)+1,k  ).GT.0.OR. &
                    &   i_wp_3d(i  ,Nm(2)+1,k  ).GT.0.OR. &
                    &   i_wp_3d(i+1,Nm(2)+1,k  ).GT.0.OR. &
                    &   i_wp_3d(i-1,Nm(2)+1,k+1).GT.0.OR. &
                    &   i_wp_3d(i  ,Nm(2)+1,k+1).GT.0.OR. &
                    &   i_wp_3d(i+1,Nm(2)+1,k+1).GT.0) buff(i,l+k)=i_wp_3d(i,Nm(2),k)
                 ENDDO
              ENDDO
!               FORALL (i=1:Nm(1),k=2:Nm(3)-1) buff(i,l+k)=i_wp_3d(i,Nm(2),k)
           ENDIF
           !bottom
           IF (iistart(3).NE.1) THEN
              l=4*nsize
              DO j=1,Nm(2)
                 DO i=1,Nm(1)
                    IF (i_wp_3d(i,j,1).LE.0) CYCLE
                    IF (i_wp_3d(i-1,j-1,0).GT.0.OR. &
                    &   i_wp_3d(i  ,j-1,0).GT.0.OR. &
                    &   i_wp_3d(i+1,j-1,0).GT.0.OR. &
                    &   i_wp_3d(i-1,j  ,0).GT.0.OR. &
                    &   i_wp_3d(i  ,j  ,0).GT.0.OR. &
                    &   i_wp_3d(i+1,j  ,0).GT.0.OR. &
                    &   i_wp_3d(i-1,j+1,0).GT.0.OR. &
                    &   i_wp_3d(i  ,j+1,0).GT.0.OR. &
                    &   i_wp_3d(i+1,j+1,0).GT.0) buff(i,l+j)=i_wp_3d(i,j,1)
                 ENDDO
              ENDDO
!               FORALL (i=1:Nm(1),j=1:Nm(2)) buff(i,l+j)=i_wp_3d(i,j,1)
           ENDIF
           !top
           IF (iiend(3).NE.Ngrid(3)) THEN
              l=5*nsize
              DO j=1,Nm(2)
                 DO i=1,Nm(1)
                    IF (i_wp_3d(i,j,Nm(3)).LE.0) CYCLE
                    IF (i_wp_3d(i-1,j-1,Nm(3)+1).GT.0.OR. &
                    &   i_wp_3d(i  ,j-1,Nm(3)+1).GT.0.OR. &
                    &   i_wp_3d(i+1,j-1,Nm(3)+1).GT.0.OR. &
                    &   i_wp_3d(i-1,j  ,Nm(3)+1).GT.0.OR. &
                    &   i_wp_3d(i  ,j  ,Nm(3)+1).GT.0.OR. &
                    &   i_wp_3d(i+1,j  ,Nm(3)+1).GT.0.OR. &
                    &   i_wp_3d(i-1,j+1,Nm(3)+1).GT.0.OR. &
                    &   i_wp_3d(i  ,j+1,Nm(3)+1).GT.0.OR. &
                    &   i_wp_3d(i+1,j+1,Nm(3)+1).GT.0) buff(i,l+j)=i_wp_3d(i,j,Nm(3))
                 ENDDO
              ENDDO
!               FORALL (i=1:Nm(1),j=1:Nm(2)) buff(i,l+j)=i_wp_3d(i,j,Nm(3))
           ENDIF
#endif

           CALL MPI_Gatherv(buff,counts(1),MPI_INTEGER,&
           & border,counts,displ,MPI_INTEGER,0,comm,info)
           or_fail_MPI("MPI_Gatherv")

        ENDIF

        IF (rank.EQ.0) THEN
           ALLOCATE(ppm_rc_border_seeds(ppm_nproc),STAT=info)
           or_fail_alloc("ppm_rc_border_seeds")

           m=0
           DO isub=1,nsubs
#if   __DIME == __2D
              c=(isub-1)*4+1
              l=COUNT(border(:,c:c+3).GT.0))
#elif __DIME == __3D
              !last column in the border data
              c=(isub-1)*6+1
              l=COUNT(border(:,:,c:c+5).GT.0)
#endif
              m=MAX(m,l)
           ENDDO

           ALLOCATE(content(m,nsubs),STAT=info)
           or_fail_alloc("content")

           content=0

           DO isub=1,nsubs
              !processor id
              is2p=topo%sub2proc(isub)+1

#if   __DIME == __2D
              c=displ(is2p)/nsize+1

              l=0
              DO j=c,c+3
                 DO i=1,nsize
                    !I do not care about negative labels which are the border
                    IF (border(i,j).LE.0) CYCLE
                    l=l+1
!                     DO m=1,l-1
!                        IF (content(m,is2p).EQ.border(i,j)) EXIT
!                     ENDDO
                    m=MAXLOC(content(1:l-1,is2p),DIM=1,MASK=content(:,is2p).EQ.content(l,is2p))

                    ld=(/i,j/)

                    SELECT CASE (m)
                    CASE (0)
                       content(l,is2p)=border(i,j)

                       ALLOCATE(seed,STAT=info)
                       or_fail_alloc("seed")
                       CALL seed%add(ld)
                       CALL ppm_rc_border_seeds(is2p)%push(seed,info)
                       or_fail("could not add new seed to the collection")

                    CASE DEFAULT
                       seedlst => ppm_rc_border_seeds(is2p)%vec(m)%t
                       CALL seedlst%add(ld)

                    END SELECT

                 ENDDO !i
              ENDDO !j
#elif __DIME == __3D
              !last column in the border data
              c=displ(is2p)/(nsize*nsize)+1

              l=0
              DO k=c,c+5
                 DO j=1,nsize
                    DO i=1,nsize
                       IF (border(i,j,k).LE.0) CYCLE
                       l=l+1

!                        DO m=1,l-1
!                           IF (content(m,is2p).EQ.border(i,j,k)) EXIT
!                        ENDDO

                       m=MAXLOC(content(1:l-1,is2p),DIM=1,MASK=content(:,is2p).EQ.content(l,is2p))

                       ld=(/i,j,k/)

                       SELECT CASE (m)
                       CASE (0)
                          content(l,is2p)=border(i,j,k)

                          ALLOCATE(seed,STAT=info)
                          or_fail_alloc("seed")
                          CALL seed%add(ld)
                          CALL ppm_rc_seeds(is2p)%push(seed,info)
                          or_fail("could not add new seed to the collection")

                       CASE DEFAULT
                          seedlst => ppm_rc_seeds(is2p)%vec(m)%t
                          CALL seedlst%add(ld)

                       END SELECT

                    ENDDO !i
                 ENDDO !j
              ENDDO !k
#endif
           ENDDO !isub=1,nsubs

           DO isub=1,nsubs
              is2p=topo%sub2proc(isub)+1

              IF (ppm_rc_seeds(is2p)%nb.LE.0) CYCLE
              seed => ppm_rc_seeds(is2p)%vec(1)%t
              content(1,is2p)=0
              DO WHILE (ASSOCIATED(seed))
                 seedlnk => seed%first

                 DO WHILE (ASSOCIATED(seedlnk))
                    seedn => seedlnk%getValue()

                    iis2p=(seedn(__DIME)-1)/(2*__DIME)
                    DO iisub=1,nsubs
                       IF (iis2p.EQ.topo%sub2proc(iisub)) EXIT
                    ENDDO

                    SELECT TYPE(p => MeshIn%subpatch_by_sub(iisub)%vec(1)%t)
                    TYPE IS (ppm_t_subpatch)
                       iistart => p%istart
                       iiend   => p%iend
                       bc      => p%bc
                    END SELECT

                    SELECT CASE (MOD(seedn(__DIME),2*__DIME))
                    CASE (1)
                       !west
                       IF (bc(1).EQ.1) THEN
                          seedlnk => seed%nextLink()
                          CYCLE
                       ENDIF
                       west_neigh: DO n=1,maxnneighsubs
                          IF (ineighsubs(n,iis2p).EQ.0) EXIT west_neigh
                          jsub=ineighsubs(n,iis2p)
                          SELECT TYPE(p => MeshIn%subpatch_by_sub(jsub)%vec(1)%t)
                          TYPE IS (ppm_t_subpatch)
                             jiend => p%iend
                          END SELECT
                          IF (iistart(1).NE.jiend(1)+1) CYCLE west_neigh
                          js2p=topo%sub2proc(jsub)+1
#if   __DIME == __2D
                          !east_neigh
                          c=displ(js2p)/nsize+2
!                           IF (border(seedn(1),c).LE.0) EXIT west_neigh
                          l=MAX(border(seedn(1)-1,c), &
                          &     border(seedn(1)  ,c), &
                          &     border(seedn(1)+1,c))
#elif __DIME == __3D
                          !east_neigh
                          c=displ(js2p)/(nsize*nsize)+2
!                           IF (border(seedn(1),seedn(2),c).LE.0) EXIT west_neigh
                          l=MAX(border(seedn(1)-1,seedn(2)-1,c), &
                          &     border(seedn(1)  ,seedn(2)-1,c), &
                          &     border(seedn(1)+1,seedn(2)-1,c), &
                          &     border(seedn(1)-1,seedn(2)  ,c), &
                          &     border(seedn(1)  ,seedn(2)  ,c), &
                          &     border(seedn(1)+1,seedn(2)  ,c), &
                          &     border(seedn(1)-1,seedn(2)+1,c), &
                          &     border(seedn(1)  ,seedn(2)+1,c), &
                          &     border(seedn(1)+1,seedn(2)+1,c))
#endif

!                           DO m=1,SIZE(content,DIM=1)
!                              IF (content(m,js2p).EQ.l) EXIT
!                           ENDDO
                          m=MAXLOC(content(:,js2p),DIM=1,MASK=content(:,js2p).EQ.l)
!                           IF (m.GT.SIZE(content,DIM=1)) EXIT west_neigh
                          IF (m.EQ.0) EXIT west_neigh

                          content(m,js2p)=0
                          seedlst => ppm_rc_seeds(js2p)%vec(m)%t
                          CALL seed%merge(seedlst)

                          max_id=ppm_rc_seeds(js2p)%max_id
                          ppm_rc_seeds(js2p)%vec(m)%t => ppm_rc_seeds(js2p)%vec(max_id)%t
                          ppm_rc_seeds(js2p)%vec(max_id)%t => NULL()
                          ppm_rc_seeds(js2p)%nb=ppm_rc_seeds(js2p)%nb-1
                          ppm_rc_seeds(js2p)%max_id=ppm_rc_seeds(js2p)%max_id-1

                          EXIT west_neigh
                       ENDDO west_neigh

                    CASE (2)
                       !east
                       IF (bc(2).EQ.1) THEN
                          seedlnk => seed%nextLink()
                          CYCLE
                       ENDIF
                       east_neigh: DO n=1,maxnneighsubs
                          IF (ineighsubs(n,iis2p).EQ.0) EXIT east_neigh
                          jsub=ineighsubs(n,iis2p)
                          SELECT TYPE(p => MeshIn%subpatch_by_sub(jsub)%vec(1)%t)
                          TYPE IS (ppm_t_subpatch)
                             jistart => p%istart
                          END SELECT
                          IF (iiend(1).NE.jistart(1)-1) CYCLE east_neigh
                          js2p=topo%sub2proc(jsub)+1
#if   __DIME == __2D
                          !west_neigh
                          c=displ(js2p)/nsize+1
!                           IF (border(seedn(1),c).LE.0) EXIT east_neigh
!                           l=border(seedn(1),c)
                          l=MAX(border(seedn(1)-1,c), &
                          &     border(seedn(1)  ,c), &
                          &     border(seedn(1)+1,c))
#elif __DIME == __3D
                          !west_neigh
                          c=displ(js2p)/(nsize*nsize)+1
!                           IF (border(seedn(1),seedn(2),c).LE.0) EXIT east_neigh
!                           l=border(seedn(1),seedn(2),c)
                          l=MAX(border(seedn(1)-1,seedn(2)-1,c), &
                          &     border(seedn(1)  ,seedn(2)-1,c), &
                          &     border(seedn(1)+1,seedn(2)-1,c), &
                          &     border(seedn(1)-1,seedn(2)  ,c), &
                          &     border(seedn(1)  ,seedn(2)  ,c), &
                          &     border(seedn(1)+1,seedn(2)  ,c), &
                          &     border(seedn(1)-1,seedn(2)+1,c), &
                          &     border(seedn(1)  ,seedn(2)+1,c), &
                          &     border(seedn(1)+1,seedn(2)+1,c))
#endif

!                           DO m=1,SIZE(content,DIM=1)
!                              IF (content(m,js2p).EQ.l) EXIT
!                           ENDDO
                          m=MAXLOC(content(:,js2p),DIM=1,MASK=content(:,js2p).EQ.l)
!                           IF (m.GT.SIZE(content,DIM=1)) EXIT east_neigh
                          IF (m.EQ.0) EXIT east_neigh

                          content(m,js2p)=0
                          seedlst => ppm_rc_seeds(js2p)%vec(m)%t
                          CALL seed%merge(seedlst)

                          max_id=ppm_rc_seeds(js2p)%max_id
                          ppm_rc_seeds(js2p)%vec(m)%t => ppm_rc_seeds(js2p)%vec(max_id)%t
                          ppm_rc_seeds(js2p)%vec(max_id)%t => NULL()
                          ppm_rc_seeds(js2p)%nb=ppm_rc_seeds(js2p)%nb-1
                          ppm_rc_seeds(js2p)%max_id=ppm_rc_seeds(js2p)%max_id-1

                          EXIT east_neigh
                       ENDDO east_neigh

                    CASE (3)
                       !south
                       IF (bc(3).EQ.1) THEN
                          seedlnk => seed%nextLink()
                          CYCLE
                       ENDIF
                       south_neigh: DO n=1,maxnneighsubs
                          IF (ineighsubs(n,iis2p).EQ.0) EXIT south_neigh
                          jsub=ineighsubs(n,iis2p)
                          SELECT TYPE(p => MeshIn%subpatch_by_sub(jsub)%vec(1)%t)
                          TYPE IS (ppm_t_subpatch)
                             jiend => p%iend
                          END SELECT
                          IF (iistart(2).NE.jiend(2)+1) CYCLE south_neigh
                          js2p=topo%sub2proc(jsub)+1
#if   __DIME == __2D
                          !north_neigh
                          c=displ(js2p)/nsize+4
!                           IF (border(seedn(1),c).LE.0) EXIT south_neigh
!                           l=border(seedn(1),c)
                          l=MAX(border(seedn(1)-1,c), &
                          &     border(seedn(1)  ,c), &
                          &     border(seedn(1)+1,c))
#elif __DIME == __3D
                          !north_neigh
                          c=displ(js2p)/(nsize*nsize)+4
!                           IF (border(seedn(1),seedn(2),c).LE.0) EXIT south_neigh
!                           l=border(seedn(1),seedn(2),c)
                          l=MAX(border(seedn(1)-1,seedn(2)-1,c), &
                          &     border(seedn(1)  ,seedn(2)-1,c), &
                          &     border(seedn(1)+1,seedn(2)-1,c), &
                          &     border(seedn(1)-1,seedn(2)  ,c), &
                          &     border(seedn(1)  ,seedn(2)  ,c), &
                          &     border(seedn(1)+1,seedn(2)  ,c), &
                          &     border(seedn(1)-1,seedn(2)+1,c), &
                          &     border(seedn(1)  ,seedn(2)+1,c), &
                          &     border(seedn(1)+1,seedn(2)+1,c))
#endif
!                           DO m=1,SIZE(content,DIM=1)
!                              IF (content(m,js2p).EQ.l) EXIT
!                           ENDDO
                          m=MAXLOC(content(:,js2p),DIM=1,MASK=content(:,js2p).EQ.l)
!                           IF (m.GT.SIZE(content,DIM=1)) EXIT south_neigh
                          IF (m.EQ.0) EXIT south_neigh

                          content(m,js2p)=0
                          seedlst => ppm_rc_seeds(js2p)%vec(m)%t
                          CALL seed%merge(seedlst)

                          max_id=ppm_rc_seeds(js2p)%max_id
                          ppm_rc_seeds(js2p)%vec(m)%t => ppm_rc_seeds(js2p)%vec(max_id)%t
                          ppm_rc_seeds(js2p)%vec(max_id)%t => NULL()
                          ppm_rc_seeds(js2p)%nb=ppm_rc_seeds(js2p)%nb-1
                          ppm_rc_seeds(js2p)%max_id=ppm_rc_seeds(js2p)%max_id-1

                          EXIT south_neigh
                       ENDDO south_neigh

                    CASE (4)
                       !north
                       IF (bc(4).EQ.1) THEN
                          seedlnk => seed%nextLink()
                          CYCLE
                       ENDIF
                       north_neigh: DO n=1,maxnneighsubs
                          IF (ineighsubs(n,iis2p).EQ.0) EXIT north_neigh
                          jsub=ineighsubs(n,iis2p)
                          SELECT TYPE(p => MeshIn%subpatch_by_sub(jsub)%vec(1)%t)
                          TYPE IS (ppm_t_subpatch)
                             jistart => p%istart
                          END SELECT
                          IF (iiend(2).NE.jistart(2)-1) CYCLE north_neigh
                          js2p=topo%sub2proc(jsub)+1
#if   __DIME == __2D
                          !north_neigh
                          c=displ(js2p)/nsize+3
!                           IF (border(seedn(1),c).LE.0) EXIT north_neigh
!                           l=border(seedn(1),c)
                          l=MAX(border(seedn(1)-1,c), &
                          &     border(seedn(1)  ,c), &
                          &     border(seedn(1)+1,c))
#elif __DIME == __3D
                          !north_neigh
                          c=displ(js2p)/(nsize*nsize)+3
!                           IF (border(seedn(1),seedn(2),c).LE.0) EXIT north_neigh
!                           l=border(seedn(1),seedn(2),c)
                          l=MAX(border(seedn(1)-1,seedn(2)-1,c), &
                          &     border(seedn(1)  ,seedn(2)-1,c), &
                          &     border(seedn(1)+1,seedn(2)-1,c), &
                          &     border(seedn(1)-1,seedn(2)  ,c), &
                          &     border(seedn(1)  ,seedn(2)  ,c), &
                          &     border(seedn(1)+1,seedn(2)  ,c), &
                          &     border(seedn(1)-1,seedn(2)+1,c), &
                          &     border(seedn(1)  ,seedn(2)+1,c), &
                          &     border(seedn(1)+1,seedn(2)+1,c))
#endif
!                           DO m=1,SIZE(content,DIM=1)
!                              IF (content(m,js2p).EQ.l) EXIT
!                           ENDDO
                          m=MAXLOC(content(:,js2p),DIM=1,MASK=content(:,js2p).EQ.l)

!                           IF (m.GT.SIZE(content,DIM=1)) EXIT north_neigh
                          IF (m.EQ.0) EXIT north_neigh

                          content(m,js2p)=0
                          seedlst => ppm_rc_seeds(js2p)%vec(m)%t
                          CALL seed%merge(seedlst)

                          max_id=ppm_rc_seeds(js2p)%max_id
                          ppm_rc_seeds(js2p)%vec(m)%t => ppm_rc_seeds(js2p)%vec(max_id)%t
                          ppm_rc_seeds(js2p)%vec(max_id)%t => NULL()
                          ppm_rc_seeds(js2p)%nb=ppm_rc_seeds(js2p)%nb-1
                          ppm_rc_seeds(js2p)%max_id=ppm_rc_seeds(js2p)%max_id-1

                          EXIT south_neigh
                       ENDDO north_neigh

#if   __DIME == __3D
                    CASE (5)
                       !bottom
                       IF (bc(5).EQ.1) THEN
                          seedlnk => seed%nextLink()
                          CYCLE
                       ENDIF
                       bottom_neigh: DO n=1,maxnneighsubs
                          IF (ineighsubs(n,iis2p).EQ.0) EXIT bottom_neigh
                          jsub=ineighsubs(n,iis2p)
                          SELECT TYPE(p => MeshIn%subpatch_by_sub(jsub)%vec(1)%t)
                          TYPE IS (ppm_t_subpatch)
                             jiend => p%iend
                          END SELECT
                          IF (iistart(3).NE.jiend(3)+1) CYCLE bottom_neigh
                          js2p=topo%sub2proc(jsub)+1
                          !top_neigh
                          c=displ(js2p)/(nsize*nsize)+6
!                           IF (border(seedn(1),seedn(2),c).LE.0) EXIT bottom_neigh
!                           l=border(seedn(1),seedn(2),c)
                          l=MAX(border(seedn(1)-1,seedn(2)-1,c), &
                          &     border(seedn(1)  ,seedn(2)-1,c), &
                          &     border(seedn(1)+1,seedn(2)-1,c), &
                          &     border(seedn(1)-1,seedn(2)  ,c), &
                          &     border(seedn(1)  ,seedn(2)  ,c), &
                          &     border(seedn(1)+1,seedn(2)  ,c), &
                          &     border(seedn(1)-1,seedn(2)+1,c), &
                          &     border(seedn(1)  ,seedn(2)+1,c), &
                          &     border(seedn(1)+1,seedn(2)+1,c))

!                           DO m=1,SIZE(content,DIM=1)
!                              IF (content(m,js2p).EQ.l) EXIT
!                           ENDDO
                          m=MAXLOC(content(:,js2p),DIM=1,MASK=content(:,js2p).EQ.l)
!                           IF (m.GT.SIZE(content,DIM=1)) EXIT bottom_neigh
                          IF (m.EQ.0) EXIT bottom_neigh

                          content(m,js2p)=0

                          seedlst => ppm_rc_seeds(js2p)%vec(m)%t
                          CALL seed%merge(seedlst)

                          max_id=ppm_rc_seeds(js2p)%max_id
                          ppm_rc_seeds(js2p)%vec(m)%t => ppm_rc_seeds(js2p)%vec(max_id)%t
                          ppm_rc_seeds(js2p)%vec(max_id)%t => NULL()
                          ppm_rc_seeds(js2p)%nb=ppm_rc_seeds(js2p)%nb-1
                          ppm_rc_seeds(js2p)%max_id=ppm_rc_seeds(js2p)%max_id-1

                          EXIT bottom_neigh
                       ENDDO bottom_neigh

                    CASE (6)
                       !top
                       IF (bc(6).EQ.1) THEN
                          seedlnk => seed%nextLink()
                          CYCLE
                       ENDIF
                       top_neigh: DO n=1,maxnneighsubs
                          IF (ineighsubs(n,iis2p).EQ.0) EXIT top_neigh
                          jsub=ineighsubs(n,iis2p)
                          SELECT TYPE(p => MeshIn%subpatch_by_sub(jsub)%vec(1)%t)
                          TYPE IS (ppm_t_subpatch)
                             jistart => p%istart
                          END SELECT
                          IF (iiend(3).NE.jistart(3)-1) CYCLE top_neigh
                          js2p=topo%sub2proc(jsub)+1
                          !bottom_neigh
                          c=displ(js2p)/(nsize*nsize)+5
!                           IF (border(seedn(1),seedn(2),c).LE.0) EXIT top_neigh
!                           l=border(seedn(1),seedn(2),c)
                          l=MAX(border(seedn(1)-1,seedn(2)-1,c), &
                          &     border(seedn(1)  ,seedn(2)-1,c), &
                          &     border(seedn(1)+1,seedn(2)-1,c), &
                          &     border(seedn(1)-1,seedn(2)  ,c), &
                          &     border(seedn(1)  ,seedn(2)  ,c), &
                          &     border(seedn(1)+1,seedn(2)  ,c), &
                          &     border(seedn(1)-1,seedn(2)+1,c), &
                          &     border(seedn(1)  ,seedn(2)+1,c), &
                          &     border(seedn(1)+1,seedn(2)+1,c))

!                           DO m=1,SIZE(content,DIM=1)
!                              IF (content(m,js2p).EQ.l) EXIT
!                           ENDDO
                          m=MAXLOC(content(:,js2p),DIM=1,MASK=content(:,js2p).EQ.l)
!                           IF (m.GT.SIZE(content,DIM=1)) EXIT top_neigh
                          IF (m.EQ.0) EXIT top_neigh

                          content(m,js2p)=0

                          seedlst => ppm_rc_seeds(js2p)%vec(m)%t
                          CALL seed%merge(seedlst)

                          max_id=ppm_rc_seeds(js2p)%max_id
                          ppm_rc_seeds(js2p)%vec(m)%t => ppm_rc_seeds(js2p)%vec(max_id)%t
                          ppm_rc_seeds(js2p)%vec(max_id)%t => NULL()
                          ppm_rc_seeds(js2p)%nb=ppm_rc_seeds(js2p)%nb-1
                          ppm_rc_seeds(js2p)%max_id=ppm_rc_seeds(js2p)%max_id-1

                          EXIT top_neigh
                       ENDDO top_neigh
#endif

                    END SELECT

                    seedlnk => seed%nextLink()
                 ENDDO !WHILE(ASSOCIATED(seedlnk))

                 DO m=1,SIZE(content,DIM=1)
                    IF (content(m,is2p).NE.0) EXIT
                 ENDDO
                 IF (m.GT.SIZE(content,DIM=1)) EXIT
                 seed => ppm_rc_seeds(is2p)%vec(m)%t
              ENDDO !WHILE(ASSOCIATED(seed))

           ENDDO !isub=1,nsubs

           DO isub=1,nsubs
              seed => ppm_rc_seeds(isub)%begin()
              DO WHILE (ASSOCIATED(seed))
                 seedlnk => seed%first
                 DO WHILE (ASSOCIATED(seedlnk))
                    seedn => seedlnk%getValue()

#if   __DIME == __2D
                    lmin=MIN(lmin,border(seedn(1),seedn(2)))
#elif __DIME == __3D
                    lmin=MIN(lmin,border(seedn(1),seedn(2),seedn(3)))
#endif
                    seedlnk => seed%nextLink()
                 ENDDO !WHILE(ASSOCIATED(seedlnk))

                 seedlnk => seed%first
                 DO WHILE (ASSOCIATED(seedlnk))
                    seedn => seedlnk%getValue()

#if   __DIME == __2D
                    border(seedn(1),seedn(2))=lmin
#elif __DIME == __3D
                    border(seedn(1),seedn(2),seedn(3))=lmin
#endif
                    seedlnk => seed%nextLink()
                 ENDDO !WHILE(ASSOCIATED(seedlnk))

                 lmin=0
                 seed => ppm_rc_seeds(isub)%next()
              ENDDO !WHILE(ASSOCIATED(seed))
           ENDDO !isub=1,nsubs

        ELSE
           buf=0
        ENDIF !(rank.EQ.0)

        CALL MPI_Scatterv(border,counts,displ,MPI_INTEGER,&
        &    buff,counts(1),MPI_INTEGER,0,comm,info)
        or_fail_MPI("MPI_Scatterv")

        sbpitr => MeshIn%subpatch%begin()
        IF (ASSOCIATED(sbpitr)) THEN
           Nm => sbpitr%nnodes
           iistart => sbpitr%istart
           iiend   => sbpitr%iend

           CALL sbpitr%get_field(FieldIn,DTYPE(i_wp),info)
           or_fail("Failed to get field i_wp data.")
#if   __DIME == __2D
           !west
           IF (iistart(1).NE.1) THEN
              FORALL (j=2:Nm(2)-1) i_wp_2d(1,j)=buff(j)
           ENDIF
           !east
           IF (iiend(1).NE.Ngrid(1)) THEN
              FORALL (j=2:Nm(2)-1) i_wp_2d(Nm(1),j)=buff(nsize+j)
           ENDIF
           !south
           IF (iistart(2).NE.1) THEN
              l=2*nsize
              FORALL (i=1:Nm(1)) i_wp_2d(i,1)=buff(l+i)
           ENDIF
           !north
           IF (iiend(2).NE.Ngrid(2)) THEN
              l=3*nsize
              FORALL (i=1:Nm(1)) i_wp_2d(i,Nm(2))=buff(l+i)
           ENDIF
#elif __DIME == __3D
           !west
           IF (iistart(1).NE.1) THEN
              FORALL (j=2:Nm(2)-1,k=2:Nm(3)-1) i_wp_3d(1,j,k)=buff(j,k)
           ENDIF
           !east
           IF (iiend(1).NE.Ngrid(1)) THEN
              FORALL (j=2:Nm(2)-1,k=2:Nm(3)-1) i_wp_3d(Nm(1),j,k)=buff(j,nsize+k)
           ENDIF
           !south
           IF (iistart(2).NE.1) THEN
              l=2*nsize
              FORALL (i=1:Nm(1),k=2:Nm(3)-1) i_wp_3d(i,1,k)=buff(i,l+k)
           ENDIF
           !north
           IF (iiend(2).NE.Ngrid(2)) THEN
              l=3*nsize
              FORALL (i=1:Nm(1),k=2:Nm(3)-1) i_wp_3d(i,Nm(2),k)=buff(i,l+k)
           ENDIF
           !bottom
           IF (iistart(3).NE.1) THEN
              l=4*nsize
              FORALL (i=1:Nm(1),j=1:Nm(2)) i_wp_3d(i,j,1)=buff(i,l+j)
           ENDIF
           !top
           IF (iiend(3).NE.Ngrid(3)) THEN
              l=5*nsize
              FORALL (i=1:Nm(1),j=1:Nm(2)) i_wp_3d(i,j,Nm(3))=buff(i,l+j)
           ENDIF
#endif
        ENDIF

        IF (rank.EQ.0) THEN
           DO isub=1,nsubs
              CALL ppm_rc_border_seeds(isub)%destroy(info)
              or_fail("ppm_rc_border_seeds(isub)%destroy")
           ENDDO
           DEALLOCATE(ppm_rc_border_seeds,STAT=info)
           or_fail_dealloc("ppm_rc_border_seeds")
        ENDIF

        DEALLOCATE(ineighsubs,STAT=info)
        or_fail_dealloc("ineighsubs")
        DEALLOCATE(displ,counts,STAT=info)
        or_fail_dealloc("displ,counts")
        DEALLOCATE(border,STAT=info)
        or_fail_dealloc("border")
        DEALLOCATE(buff,STAT=info)
        or_fail_dealloc("buff")

        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
      9999 CONTINUE
        CALL substop(caller,t0,info)
        RETURN
      END SUBROUTINE DTYPE(ppm_rc_ghost_unique)


