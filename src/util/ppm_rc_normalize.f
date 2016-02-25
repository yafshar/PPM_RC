      !-------------------------------------------------------------------------
      !  Subroutine   :                    ppm_rc_normalize
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
      !  Subroutine   :                    ppm_rc_normalize
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Shifts and scales the intensity values in the image
      !                 to [0,1].
      !
      !
      !  Input        : FieldIn1, FieldIn2, MeshIn
      !
      !  Input/output :
      !
      !  Output       : info  (I) return status. 0 on success.
      !
      !  Routines     :
      !
      !
      !  Remarks      : in this process ghost layers are not considered
      !               : !TODO in case lmax=lmin it will produce seg fault
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      SUBROUTINE DTYPE(ppm_rc_normalize1)(FieldIn1,FieldIn2,MeshIn,info)

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
        CLASS(ppm_t_field_),     POINTER       :: FieldIn1
        !!! the source filed which is supposed to be normalized
        CLASS(ppm_t_field_),     POINTER       :: FieldIn2
        !!! the output normalized field

        CLASS(ppm_t_equi_mesh_), POINTER       :: MeshIn
        !!! source mesh on which fields are descritized

        INTEGER,                 INTENT(  OUT) :: info
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
        CLASS(ppm_t_subpatch_), POINTER :: sbpitr

#if   __DIME == __2D
        REAL(MK),              CONTIGUOUS, DIMENSION(:,:),   POINTER :: DTYPE(wp1)
        REAL(MK),              CONTIGUOUS, DIMENSION(:,:),   POINTER :: DTYPE(wp2)
        REAL(ppm_kind_single), CONTIGUOUS, DIMENSION(:,:),   POINTER :: DTYPE(wp2s)
#elif __DIME == __3D
        REAL(MK),              CONTIGUOUS, DIMENSION(:,:,:), POINTER :: DTYPE(wp1)
        REAL(MK),              CONTIGUOUS, DIMENSION(:,:,:), POINTER :: DTYPE(wp2)
        REAL(ppm_kind_single), CONTIGUOUS, DIMENSION(:,:,:), POINTER :: DTYPE(wp2s)
#endif
        REAL(ppm_kind_double)                                        :: t0
        REAL(MK)                                                     :: coef
        REAL(MK)                                                     :: lmin,lmax

        INTEGER, DIMENSION(:), POINTER :: Nm
        INTEGER, DIMENSION(:), POINTER :: hi_a
        INTEGER, DIMENSION(:), POINTER :: lo_a
        INTEGER                        :: i,j
#if   __DIME == __3D
        INTEGER                        :: k
#endif
#ifdef __MPI
        INTEGER, DIMENSION(2)          :: request
#endif

        CHARACTER(LEN=ppm_char) :: caller = 'ppm_rc_normalize'

        !-------------------------------------------------------------------------
        !  Initialize
        !-------------------------------------------------------------------------
        CALL substart(caller,t0,info)

        NULLIFY(DTYPE(wp1),DTYPE(wp2))

        SELECT CASE (lNormalize)
        CASE (.FALSE.)
        !It means thet FieldIn1 has not been previously normalized
           lmin=big
           lmax=-big

           sbpitr => MeshIn%subpatch%begin()
           DO WHILE (ASSOCIATED(sbpitr))
              Nm => sbpitr%nnodes

              CALL sbpitr%get_field(FieldIn1,DTYPE(wp1),info)
              or_fail("Failed to get field wp data.")

#if   __DIME == __2D
              DO j=1,Nm(2)
                 DO i=1,Nm(1)
                    IF (DTYPE(wp1)(i,j).LT.lmin) lmin=DTYPE(wp1)(i,j)
                    IF (DTYPE(wp1)(i,j).GT.lmax) lmax=DTYPE(wp1)(i,j)
                 ENDDO
              ENDDO
#elif __DIME == __3D
              DO k=1,Nm(3)
                 DO j=1,Nm(2)
                    DO i=1,Nm(1)
                       IF (DTYPE(wp1)(i,j,k).LT.lmin) lmin=DTYPE(wp1)(i,j,k)
                       IF (DTYPE(wp1)(i,j,k).GT.lmax) lmax=DTYPE(wp1)(i,j,k)
                    ENDDO
                 ENDDO
              ENDDO
#endif

              sbpitr => MeshIn%subpatch%next()
           ENDDO

#ifdef __MPI
           CALL MPI_Iallreduce(MPI_IN_PLACE,lmin,1, &
           &    ppm_mpi_kind,MPI_MIN,comm,request(1),info)
           or_fail_MPI("MPI_Iallreduce failed")

           CALL MPI_Iallreduce(MPI_IN_PLACE,lmax,1, &
           &    ppm_mpi_kind,MPI_MAX,comm,request(2),info)
           or_fail_MPI("MPI_Iallreduce failed")

           CALL MPI_Waitall(2,request,MPI_STATUSES_IGNORE,info)
           or_fail_MPI("MPI_Waitall")
#endif

           !-------------------------------------------------------------------------
           !  Scale and shift the FieldIn2
           !-------------------------------------------------------------------------
           Normalfac=lmax-lmin
           coef=one/Normalfac

           Tolerance_=Tolerance*coef
           !correct the tolerance, when we are normalizing the image

           sbpitr => MeshIn%subpatch%begin()
           DO WHILE (ASSOCIATED(sbpitr))
              Nm => sbpitr%nnodes

              CALL sbpitr%get_field(FieldIn1,DTYPE(wp1),info)
              or_fail("Failed to get field wp data.")

              SELECT CASE (FieldIn2%data_type)
              CASE (ppm_type_real)
                 CALL sbpitr%get_field(FieldIn2,DTYPE(wp2),info)
                 or_fail("Failed to get field wp data.")

#if   __DIME == __2D
                 FORALL (i=1:Nm(1),j=1:Nm(2))
                    DTYPE(wp2)(i,j)=(DTYPE(wp1)(i,j)-lmin)*coef
#elif __DIME == __3D
                 FORALL (i=1:Nm(1),j=1:Nm(2),k=1:Nm(3))
                    DTYPE(wp2)(i,j,k)=(DTYPE(wp1)(i,j,k)-lmin)*coef
#endif
                 END FORALL

              CASE (ppm_type_real_single)
                 IF (ppm_kind.EQ.ppm_kind_double) THEN
                    CALL sbpitr%get_field(FieldIn2,DTYPE(wp2s),info)
                    or_fail("Failed to get field wp data.")

#if   __DIME == __2D
                    FORALL (i=1:Nm(1),j=1:Nm(2))
                       DTYPE(wp2s)(i,j)=REAL((DTYPE(wp1)(i,j)-lmin)*coef,ppm_kind_single)
#elif __DIME == __3D
                    FORALL (i=1:Nm(1),j=1:Nm(2),k=1:Nm(3))
                       DTYPE(wp2s)(i,j,k)=REAL((DTYPE(wp1)(i,j,k)-lmin)*coef,ppm_kind_single)
#endif
                    END FORALL

                 ELSE IF (ppm_kind.EQ.ppm_kind_single) THEN
                    CALL sbpitr%get_field(FieldIn2,DTYPE(wp2),info)
                    or_fail("Failed to get field wp data.")

#if   __DIME == __2D
                    FORALL (i=1:Nm(1),j=1:Nm(2))
                       DTYPE(wp2)(i,j)=(DTYPE(wp1)(i,j)-lmin)*coef
#elif __DIME == __3D
                    FORALL (i=1:Nm(1),j=1:Nm(2),k=1:Nm(3))
                       DTYPE(wp2)(i,j,k)=(DTYPE(wp1)(i,j,k)-lmin)*coef
#endif
                    END FORALL

                 ENDIF

              END SELECT

              sbpitr => MeshIn%subpatch%next()
           ENDDO

        CASE (.TRUE.)
           sbpitr => MeshIn%subpatch%begin()
           DO WHILE (ASSOCIATED(sbpitr))
              Nm   => sbpitr%nnodes
              hi_a => sbpitr%hi_a
              lo_a => sbpitr%lo_a

              CALL sbpitr%get_field(FieldIn1,DTYPE(wp1),info)
              or_fail("Failed to get field i_wp data.")

              SELECT CASE (FieldIn2%data_type)
              CASE (ppm_type_real)
                 CALL sbpitr%get_field(FieldIn2,DTYPE(wp2),info)
                 or_fail("Failed to get field wp data.")

#if   __DIME == __2D
                 DO j=lo_a(2),hi_a(2)
                    DO i=lo_a(1),hi_a(1)
!                  FORALL (i=lo_a(1):hi_a(1),j=lo_a(2):hi_a(2))
                       DTYPE(wp2)(i,j)=DTYPE(wp1)(i,j)
                    ENDDO !i=lo_a(1),hi_a(1)
                 ENDDO !j=lo_a(2),hi_a(2)
#elif __DIME == __3D
                 DO k=lo_a(3),hi_a(3)
                    DO j=lo_a(2),hi_a(2)
                       DO i=lo_a(1),hi_a(1)
!                  FORALL (i=lo_a(1):hi_a(1),j=lo_a(2):hi_a(2),k=lo_a(3):hi_a(3))
                          DTYPE(wp2)(i,j,k)=DTYPE(wp1)(i,j,k)
                       ENDDO !i=lo_a(1),hi_a(1)
                    ENDDO !j=lo_a(2),hi_a(2)
                 ENDDO !k=lo_a(3),hi_a(3)
#endif
!                  END FORALL

              CASE (ppm_type_real_single)
                 IF (ppm_kind.EQ.ppm_kind_double) THEN
                    CALL sbpitr%get_field(FieldIn2,DTYPE(wp2s),info)
                    or_fail("Failed to get field wp data.")

#if   __DIME == __2D
                    DO j=lo_a(2),hi_a(2)
                       DO i=lo_a(1),hi_a(1)
                          DTYPE(wp2s)(i,j)=REAL(DTYPE(wp1)(i,j),ppm_kind_single)
                       ENDDO !i=lo_a(1),hi_a(1)
                    ENDDO !j=lo_a(2),hi_a(2)
#elif __DIME == __3D
                    DO k=lo_a(3),hi_a(3)
                       DO j=lo_a(2),hi_a(2)
                          DO i=lo_a(1),hi_a(1)
                             DTYPE(wp2s)(i,j,k)=REAL(DTYPE(wp1)(i,j,k),ppm_kind_single)
                          ENDDO !i=lo_a(1),hi_a(1)
                       ENDDO !j=lo_a(2),hi_a(2)
                    ENDDO !k=lo_a(3),hi_a(3)
#endif
                 ELSE IF (ppm_kind.EQ.ppm_kind_single) THEN
                    CALL sbpitr%get_field(FieldIn2,DTYPE(wp2),info)
                    or_fail("Failed to get field wp data.")

#if   __DIME == __2D
                    DO j=lo_a(2),hi_a(2)
                       DO i=lo_a(1),hi_a(1)
                          DTYPE(wp2)(i,j)=DTYPE(wp1)(i,j)
                       ENDDO !i=lo_a(1),hi_a(1)
                    ENDDO !j=lo_a(2),hi_a(2)
#elif __DIME == __3D
                    DO k=lo_a(3),hi_a(3)
                       DO j=lo_a(2),hi_a(2)
                          DO i=lo_a(1),hi_a(1)
                             DTYPE(wp2)(i,j,k)=DTYPE(wp1)(i,j,k)
                          ENDDO !i=lo_a(1),hi_a(1)
                       ENDDO !j=lo_a(2),hi_a(2)
                    ENDDO !k=lo_a(3),hi_a(3)
#endif
                 ENDIF

              END SELECT

              sbpitr => MeshIn%subpatch%next()
           ENDDO

        END SELECT

        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
      9999 CONTINUE
        CALL substop(caller,t0,info)
        RETURN
      CONTAINS
      SUBROUTINE check
      IMPLICIT NONE
        IF (FieldIn1%data_type.NE.ppm_type_real) THEN
           fail("FieldIn1 data type is not correct!",exit_point=8888)
        ENDIF
        IF (FieldIn2%data_type.NE.ppm_type_real) THEN
           fail("FieldIn2 data type is not correct!",exit_point=8888)
        ENDIF
      8888 CONTINUE
      END SUBROUTINE check
      END SUBROUTINE DTYPE(ppm_rc_normalize1)


      SUBROUTINE DTYPE(ppm_rc_normalize2)(FieldIn,MeshIn,info,MinThr,MaxThr)

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
        CLASS(ppm_t_field_),     POINTER       :: FieldIn
        !!! the source filed which is supposed to be normalized

        CLASS(ppm_t_equi_mesh_), POINTER       :: MeshIn
        !!! source mesh on which fields are descritized

        INTEGER,                 INTENT(  OUT) :: info

        REAL(MK), OPTIONAL,      INTENT(IN   ) :: MinThr
        !!! Lower threshhold value in normalizing the image
        REAL(MK), OPTIONAL,      INTENT(IN   ) :: MaxThr
        !!! Upper threshhold value in normalizeing the image
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
        CLASS(ppm_t_subpatch_), POINTER :: sbpitr

#if   __DIME == __2D
        REAL(MK), CONTIGUOUS, DIMENSION(:,:),   POINTER :: DTYPE(wp)
#elif __DIME == __3D
        REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER :: DTYPE(wp)
#endif
        REAL(ppm_kind_double)                           :: t0
        REAL(MK)                                        :: coef
        REAL(MK)                                        :: lmin,lmax

        INTEGER, DIMENSION(:), POINTER :: Nm
        INTEGER                        :: i,j
#if   __DIME == __3D
        INTEGER                        :: k
#endif
#ifdef __MPI
        INTEGER, DIMENSION(2)          :: request
#endif

        CHARACTER(LEN=ppm_char) :: caller = 'ppm_rc_normalize'

        LOGICAL :: lthr

        !-------------------------------------------------------------------------
        !  Initialize
        !-------------------------------------------------------------------------
        CALL substart(caller,t0,info)

        IF (PRESENT(MinThr).AND.PRESENT(MaxThr)) THEN
           lmin=MinThr
           lmax=MaxThr
           lthr=.TRUE.
        ELSE IF (PRESENT(MinThr).OR.PRESENT(MaxThr)) THEN
           lmin=zero
           lmax=MERGE(MinThr,MaxThr,PRESENT(MinThr))
           lthr=.TRUE.
        ELSE
           lmin=big
           lmax=-big
           lthr=.FALSE.
        ENDIF

        NULLIFY(DTYPE(wp))

        sbpitr => MeshIn%subpatch%begin()
        DO WHILE (ASSOCIATED(sbpitr))
           Nm => sbpitr%nnodes

           CALL sbpitr%get_field(FieldIn,DTYPE(wp),info)
           or_fail("Failed to get field wp data.")

           SELECT CASE (lthr)
           CASE (.TRUE.)
#if   __DIME == __2D
              DO j=1,Nm(2)
                 DO i=1,Nm(1)
                    IF (DTYPE(wp)(i,j).GT.lmax) THEN
                       DTYPE(wp)(i,j)=lmax
                    ELSE IF (DTYPE(wp)(i,j).LT.lmin) THEN
                       DTYPE(wp)(i,j)=lmin
                    ENDIF
                 ENDDO
              ENDDO
#elif __DIME == __3D
              DO k=1,Nm(3)
                 DO j=1,Nm(2)
                    DO i=1,Nm(1)
                       IF (DTYPE(wp)(i,j,k).GT.lmax) THEN
                          DTYPE(wp)(i,j,k)=lmax
                       ELSE IF (DTYPE(wp)(i,j,k).LT.lmin) THEN
                          DTYPE(wp)(i,j,k)=lmin
                       ENDIF
                    ENDDO
                 ENDDO
              ENDDO
#endif

           CASE (.FALSE.)
#if   __DIME == __2D
              DO j=1,Nm(2)
                 DO i=1,Nm(1)
                    IF (DTYPE(wp)(i,j).LT.lmin) lmin=DTYPE(wp)(i,j)
                    IF (DTYPE(wp)(i,j).GT.lmax) lmax=DTYPE(wp)(i,j)
                 ENDDO
              ENDDO
#elif __DIME == __3D
              DO k=1,Nm(3)
                 DO j=1,Nm(2)
                    DO i=1,Nm(1)
                       IF (DTYPE(wp)(i,j,k).LT.lmin) lmin=DTYPE(wp)(i,j,k)
                       IF (DTYPE(wp)(i,j,k).GT.lmax) lmax=DTYPE(wp)(i,j,k)
                    ENDDO
                 ENDDO
              ENDDO
#endif

           END SELECT

           sbpitr => MeshIn%subpatch%next()
        ENDDO

        SELECT CASE (lthr)
        CASE (.TRUE.)
        CASE (.FALSE.)
#ifdef __MPI
           CALL MPI_Iallreduce(MPI_IN_PLACE,lmin,1, &
           &    ppm_mpi_kind,MPI_MIN,comm,request(1),info)
           or_fail_MPI("MPI_Iallreduce failed")

           CALL MPI_Iallreduce(MPI_IN_PLACE,lmax,1, &
           &    ppm_mpi_kind,MPI_MAX,comm,request(2),info)
           or_fail_MPI("MPI_Iallreduce failed")

           CALL MPI_Waitall(2,request,MPI_STATUSES_IGNORE,info)
           or_fail_MPI("MPI_Waitall")
#endif

        END SELECT

!         subimageintensity=zerod
!         imageintensity=zerod

        !-------------------------------------------------------------------------
        !  Scale and shift the FieldIn2
        !-------------------------------------------------------------------------
        Normalfac=lmax-lmin
        coef=one/Normalfac

        Tolerance_=Tolerance*coef
        !correct the tolerance, when we are normalizing the image

        sbpitr => MeshIn%subpatch%begin()
        DO WHILE (ASSOCIATED(sbpitr))
           Nm => sbpitr%nnodes

           CALL sbpitr%get_field(FieldIn,DTYPE(wp),info)
           or_fail("Failed to get field wp data.")

#if   __DIME == __2D
           DO j=1,Nm(2)
              DO i=1,Nm(1)
                 DTYPE(wp)(i,j)=(DTYPE(wp)(i,j)-lmin)*coef
!                  subimageintensity=subimageintensity+DTYPE(wp)(i,j)
              ENDDO
           ENDDO
#elif __DIME == __3D
           DO k=1,Nm(3)
              DO j=1,Nm(2)
                 DO i=1,Nm(1)
                    DTYPE(wp)(i,j,k)=(DTYPE(wp)(i,j,k)-lmin)*coef
!                     subimageintensity=subimageintensity+DTYPE(wp)(i,j,k)
                 ENDDO
              ENDDO
           ENDDO
#endif

! #ifdef __MPI
!            CALL MPI_Iallreduce(subimageintensity,imageintensity,1, &
!            &    MPI_DOUBLE_PRECISION,MPI_SUM,comm,request(1),info)
!            or_fail_MPI("MPI_Iallreduce failed")
!
!            CALL MPI_Wait(request(1),MPI_STATUS_IGNORE,info)
!            or_fail_MPI("MPI_Wait")
! #endif

           sbpitr => MeshIn%subpatch%next()
        ENDDO

        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
      9999 CONTINUE
        CALL substop(caller,t0,info)
        RETURN
      CONTAINS
      SUBROUTINE check
      IMPLICIT NONE
        IF (FieldIn%data_type.NE.ppm_type_real) THEN
           fail("Field data type is not correct!",exit_point=8888)
        ENDIF
      8888 CONTINUE
      END SUBROUTINE check
      END SUBROUTINE DTYPE(ppm_rc_normalize2)
