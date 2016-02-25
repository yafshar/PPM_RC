      !-------------------------------------------------------------------------
      !  Subroutine   :                    ppm_rc_CopyImageAndNormalize
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
      !  Author           - y.afshar           Dec   2015
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Subroutine   :                    ppm_rc_CopyImageAndNormalize
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Shifts and scales the intensity values in the image
      !                 to [0,1] and store the result in the output
      !
      !
      !  Input        : FieldIn1, FieldIn2, MeshIn
      !
      !  Input/output :
      !
      !  Output       : info  (I) return status. 0 on success.
      !
      !  Note         : This routine per default only copy the input image to the
      !                 output one and do not use ghosts, do not normalize!
      !                 This is the user responsibility to know whether the ghosts
      !                 are updated before calling this routine or they will get
      !                 updated after!
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      SUBROUTINE DTYPE(ppm_rc_CopyImageAndNormalize)(FieldIn1,FieldIn2,MeshIn,info, &
      &          withGhost,Normalize,FieldMinVal,FieldMaxVal)

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

        LOGICAL,  OPTIONAL,      INTENT(IN   ) :: withGhost
        !!! Do operation with or without ghost
        !!! (If FieldIn1 ghosts are uptodate, it is safe to do operation with ghosts)
        LOGICAL,  OPTIONAL,      INTENT(IN   ) :: Normalize
        !!! Whether normalize the image or just copy the content

        REAL(MK), OPTIONAL,      INTENT(  OUT) :: FieldMinVal
        REAL(MK), OPTIONAL,      INTENT(  OUT) :: FieldMaxVal
        !!! Return the global min and max value
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
        REAL(MK)                                                     :: FieldMinVal_,FieldMaxVal_
#ifdef __MPI
        REAL(MK), DIMENSION(2)                                       :: FieldMinMax
#endif

#if   __DIME == __2D
        INTEGER, CONTIGUOUS, DIMENSION(:,:),   POINTER :: DTYPE(wpi1)
        INTEGER, CONTIGUOUS, DIMENSION(:,:),   POINTER :: DTYPE(wpi2)
#elif __DIME == __3D
        INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER :: DTYPE(wpi1)
        INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER :: DTYPE(wpi2)
#endif
        INTEGER,             DIMENSION(:),     POINTER :: Nm
        INTEGER,             DIMENSION(:),     POINTER :: hi_a
        INTEGER,             DIMENSION(:),     POINTER :: lo_a
        INTEGER                                        :: i,j
#if   __DIME == __3D
        INTEGER                                        :: k
#endif
#ifdef __MPI
        INTEGER                                        :: request
#endif

        CHARACTER(LEN=ppm_char) :: caller = 'ppm_rc_CopyImageAndNormalize'

        LOGICAL ::Normalize_,withGhost_

        !-------------------------------------------------------------------------
        !  Initialize
        !-------------------------------------------------------------------------
        CALL substart(caller,t0,info)

        withGhost_=MERGE(withGhost,.FALSE.,PRESENT(withGhost))
        Normalize_=MERGE(Normalize,.FALSE.,PRESENT(Normalize))

        NULLIFY(DTYPE(wp1),DTYPE(wp2))

        SELECT CASE (Normalize_)
        CASE (.TRUE.)
           !It means thet FieldIn1 has not been previously normalized
           FieldMinVal_=big
           FieldMaxVal_=-big

           sbpitr => MeshIn%subpatch%begin()
           DO WHILE (ASSOCIATED(sbpitr))
              Nm => sbpitr%nnodes

              CALL sbpitr%get_field(FieldIn1,DTYPE(wp1),info)
              or_fail("Failed to get field wp data.")

#if   __DIME == __2D
              DO j=1,Nm(2)
                 DO i=1,Nm(1)
                    IF (DTYPE(wp1)(i,j).LT.FieldMinVal_) FieldMinVal_=DTYPE(wp1)(i,j)
                    IF (DTYPE(wp1)(i,j).GT.FieldMaxVal_) FieldMaxVal_=DTYPE(wp1)(i,j)
                 ENDDO
              ENDDO
#elif __DIME == __3D
              DO k=1,Nm(3)
                 DO j=1,Nm(2)
                    DO i=1,Nm(1)
                       IF (DTYPE(wp1)(i,j,k).LT.FieldMinVal_) FieldMinVal_=DTYPE(wp1)(i,j,k)
                       IF (DTYPE(wp1)(i,j,k).GT.FieldMaxVal_) FieldMaxVal_=DTYPE(wp1)(i,j,k)
                    ENDDO
                 ENDDO
              ENDDO
#endif

              sbpitr => MeshIn%subpatch%next()
           ENDDO

#ifdef __MPI

           FieldMinMax(1)= FieldMinVal_
           FieldMinMax(2)=-FieldMaxVal_

           CALL MPI_Iallreduce(MPI_IN_PLACE,FieldMinMax,2, &
           &    ppm_mpi_kind,MPI_MIN,comm,request,info)
           or_fail_MPI("MPI_Iallreduce failed")

           CALL MPI_Wait(request,MPI_STATUS_IGNORE,info)
           or_fail_MPI("MPI_Wait")

           FieldMinVal_= FieldMinMax(1)
           FieldMaxVal_=-FieldMinMax(2)
#endif

           !-------------------------------------------------------------------------
           !  Scale and shift the FieldIn2
           !-------------------------------------------------------------------------
           coef=one/(FieldMaxVal_-FieldMinVal_)

           sbpitr => MeshIn%subpatch%begin()
           DO WHILE (ASSOCIATED(sbpitr))
              Nm   => sbpitr%nnodes
              hi_a => sbpitr%hi_a
              lo_a => sbpitr%lo_a

              CALL sbpitr%get_field(FieldIn1,DTYPE(wp1),info)
              or_fail("Failed to get field wp data.")

              SELECT CASE (FieldIn2%data_type)
              CASE (ppm_type_real)
                 CALL sbpitr%get_field(FieldIn2,DTYPE(wp2),info)
                 or_fail("Failed to get field wp data.")

                 SELECT CASE (withGhost_)
                 CASE (.FALSE.)
#if   __DIME == __2D
                    FORALL (i=1:Nm(1),j=1:Nm(2))
                       DTYPE(wp2)(i,j)=(DTYPE(wp1)(i,j)-FieldMinVal_)*coef
#elif __DIME == __3D
                    FORALL (i=1:Nm(1),j=1:Nm(2),k=1:Nm(3))
                       DTYPE(wp2)(i,j,k)=(DTYPE(wp1)(i,j,k)-FieldMinVal_)*coef
#endif
                    END FORALL
                 CASE (.TRUE.)
#if   __DIME == __2D
                    DO j=lo_a(2),hi_a(2)
                       DO i=lo_a(1),hi_a(1)
                          DTYPE(wp2)(i,j)=(DTYPE(wp1)(i,j)-FieldMinVal_)*coef
                       ENDDO !i=lo_a(1),hi_a(1)
                    ENDDO !j=lo_a(2),hi_a(2)
#elif __DIME == __3D
                    DO k=lo_a(3),hi_a(3)
                       DO j=lo_a(2),hi_a(2)
                          DO i=lo_a(1),hi_a(1)
                             DTYPE(wp2)(i,j,k)=(DTYPE(wp1)(i,j,k)-FieldMinVal_)*coef
                          ENDDO !i=lo_a(1),hi_a(1)
                       ENDDO !j=lo_a(2),hi_a(2)
                    ENDDO !k=lo_a(3),hi_a(3)
#endif

                 END SELECT

              CASE (ppm_type_real_single)
                 IF (ppm_kind.EQ.ppm_kind_double) THEN
                    CALL sbpitr%get_field(FieldIn2,DTYPE(wp2s),info)
                    or_fail("Failed to get field wp data.")

                    SELECT CASE (withGhost_)
                    CASE (.FALSE.)
#if   __DIME == __2D
                       FORALL (i=1:Nm(1),j=1:Nm(2))
                          DTYPE(wp2s)(i,j)=REAL((DTYPE(wp1)(i,j)-FieldMinVal_)*coef,ppm_kind_single)
#elif __DIME == __3D
                       FORALL (i=1:Nm(1),j=1:Nm(2),k=1:Nm(3))
                          DTYPE(wp2s)(i,j,k)=REAL((DTYPE(wp1)(i,j,k)-FieldMinVal_)*coef,ppm_kind_single)
#endif
                       END FORALL
                    CASE (.TRUE.)
#if   __DIME == __2D
                       DO j=lo_a(2),hi_a(2)
                          DO i=lo_a(1),hi_a(1)
                             DTYPE(wp2s)(i,j)=REAL((DTYPE(wp1)(i,j)-FieldMinVal_)*coef,ppm_kind_single)
                          ENDDO !i=lo_a(1),hi_a(1)
                       ENDDO !j=lo_a(2),hi_a(2)
#elif __DIME == __3D
                       DO k=lo_a(3),hi_a(3)
                          DO j=lo_a(2),hi_a(2)
                             DO i=lo_a(1),hi_a(1)
                                DTYPE(wp2s)(i,j,k)=REAL((DTYPE(wp1)(i,j,k)-FieldMinVal_)*coef,ppm_kind_single)
                             ENDDO !i=lo_a(1),hi_a(1)
                          ENDDO !j=lo_a(2),hi_a(2)
                       ENDDO !k=lo_a(3),hi_a(3)
#endif
                    END SELECT

                 ELSE IF (ppm_kind.EQ.ppm_kind_single) THEN
                    CALL sbpitr%get_field(FieldIn2,DTYPE(wp2),info)
                    or_fail("Failed to get field wp data.")

                    SELECT CASE (withGhost_)
                    CASE (.FALSE.)
#if   __DIME == __2D
                       FORALL (i=1:Nm(1),j=1:Nm(2))
                          DTYPE(wp2)(i,j)=(DTYPE(wp1)(i,j)-FieldMinVal_)*coef
#elif __DIME == __3D
                       FORALL (i=1:Nm(1),j=1:Nm(2),k=1:Nm(3))
                          DTYPE(wp2)(i,j,k)=(DTYPE(wp1)(i,j,k)-FieldMinVal_)*coef
#endif
                       END FORALL
                    CASE (.TRUE.)
#if   __DIME == __2D
                       DO j=lo_a(2),hi_a(2)
                          DO i=lo_a(1),hi_a(1)
                             DTYPE(wp2)(i,j)=(DTYPE(wp1)(i,j)-FieldMinVal_)*coef
                          ENDDO !i=lo_a(1),hi_a(1)
                       ENDDO !j=lo_a(2),hi_a(2)
#elif __DIME == __3D
                       DO k=lo_a(3),hi_a(3)
                          DO j=lo_a(2),hi_a(2)
                             DO i=lo_a(1),hi_a(1)
                                DTYPE(wp2)(i,j,k)=(DTYPE(wp1)(i,j,k)-FieldMinVal_)*coef
                             ENDDO !i=lo_a(1),hi_a(1)
                          ENDDO !j=lo_a(2),hi_a(2)
                       ENDDO !k=lo_a(3),hi_a(3)
#endif
                    END SELECT
                 ENDIF

              END SELECT

              sbpitr => MeshIn%subpatch%next()
           ENDDO

           IF (PRESENT(FieldMinVal)) FieldMinVal=FieldMinVal_
           IF (PRESENT(FieldMaxVal)) FieldMaxVal=FieldMaxVal_

        CASE (.FALSE.) !(Normalize_)
           sbpitr => MeshIn%subpatch%begin()
           DO WHILE (ASSOCIATED(sbpitr))
              Nm   => sbpitr%nnodes
              hi_a => sbpitr%hi_a
              lo_a => sbpitr%lo_a

              SELECT CASE (FieldIn2%data_type)
              CASE (ppm_type_real)
                 CALL sbpitr%get_field(FieldIn1,DTYPE(wp1),info)
                 or_fail("Failed to get field i_wp data.")

                 CALL sbpitr%get_field(FieldIn2,DTYPE(wp2),info)
                 or_fail("Failed to get field wp data.")

                 SELECT CASE (withGhost_)
                 CASE (.FALSE.)
#if   __DIME == __2D
                    FORALL (i=1:Nm(1),j=1:Nm(2))
                       DTYPE(wp2)(i,j)=DTYPE(wp1)(i,j)
#elif __DIME == __3D
                    FORALL (i=1:Nm(1),j=1:Nm(2),k=1:Nm(3))
                       DTYPE(wp2)(i,j,k)=DTYPE(wp1)(i,j,k)
#endif
                    END FORALL
                 CASE (.TRUE.)
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
                 END SELECT

              CASE (ppm_type_real_single)
                 CALL sbpitr%get_field(FieldIn1,DTYPE(wp1),info)
                 or_fail("Failed to get field i_wp data.")

                 IF (ppm_kind.EQ.ppm_kind_double) THEN
                    CALL sbpitr%get_field(FieldIn2,DTYPE(wp2s),info)
                    or_fail("Failed to get field wp data.")

                    SELECT CASE (withGhost_)
                    CASE (.FALSE.)
#if   __DIME == __2D
                       FORALL (i=1:Nm(1),j=1:Nm(2))
                          DTYPE(wp2s)(i,j)=REAL(DTYPE(wp1)(i,j),ppm_kind_single)
#elif __DIME == __3D
                       FORALL (i=1:Nm(1),j=1:Nm(2),k=1:Nm(3))
                          DTYPE(wp2s)(i,j,k)=REAL(DTYPE(wp1)(i,j,k),ppm_kind_single)
#endif
                       END FORALL
                    CASE (.TRUE.)
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
                    END SELECT

                 ELSE IF (ppm_kind.EQ.ppm_kind_single) THEN
                    CALL sbpitr%get_field(FieldIn2,DTYPE(wp2),info)
                    or_fail("Failed to get field wp data.")

                    SELECT CASE (withGhost_)
                    CASE (.FALSE.)
#if   __DIME == __2D
                       FORALL (i=1:Nm(1),j=1:Nm(2))
                          DTYPE(wp2)(i,j)=DTYPE(wp1)(i,j)
#elif __DIME == __3D
                       FORALL (i=1:Nm(1),j=1:Nm(2),k=1:Nm(3))
                          DTYPE(wp2)(i,j,k)=DTYPE(wp1)(i,j,k)
#endif
                       END FORALL
                    CASE (.TRUE.)
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
                    END SELECT
                 ENDIF

              CASE (ppm_type_int)
                 CALL sbpitr%get_field(FieldIn1,DTYPE(wpi1),info)
                 or_fail("Failed to get field i_wp data.")

                 CALL sbpitr%get_field(FieldIn2,DTYPE(wpi2),info)
                 or_fail("Failed to get field wp data.")

                 SELECT CASE (withGhost_)
                 CASE (.FALSE.)
#if   __DIME == __2D
                    FORALL (i=1:Nm(1),j=1:Nm(2))
                       DTYPE(wpi2)(i,j)=DTYPE(wpi1)(i,j)
#elif __DIME == __3D
                    FORALL (i=1:Nm(1),j=1:Nm(2),k=1:Nm(3))
                       DTYPE(wpi2)(i,j,k)=DTYPE(wpi1)(i,j,k)
#endif
                    END FORALL
                 CASE (.TRUE.)
#if   __DIME == __2D
                    DO j=lo_a(2),hi_a(2)
                       DO i=lo_a(1),hi_a(1)
                          DTYPE(wpi2)(i,j)=DTYPE(wpi1)(i,j)
                       ENDDO !i=lo_a(1),hi_a(1)
                    ENDDO !j=lo_a(2),hi_a(2)
#elif __DIME == __3D
                    DO k=lo_a(3),hi_a(3)
                       DO j=lo_a(2),hi_a(2)
                          DO i=lo_a(1),hi_a(1)
                             DTYPE(wpi2)(i,j,k)=DTYPE(wpi1)(i,j,k)
                          ENDDO !i=lo_a(1),hi_a(1)
                       ENDDO !j=lo_a(2),hi_a(2)
                    ENDDO !k=lo_a(3),hi_a(3)
#endif
                 END SELECT !withGhost_

              END SELECT !FieldIn2%data_type

              sbpitr => MeshIn%subpatch%next()
           ENDDO

        END SELECT !Normalize_

        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
      9999 CONTINUE
        CALL substop(caller,t0,info)
        RETURN
      CONTAINS
      SUBROUTINE check
      IMPLICIT NONE
      SELECT CASE (Normalize_)
      CASE (.TRUE.)
         IF (FieldIn2%data_type.EQ.ppm_type_int) THEN
            fail("FieldIn2 data type is not correct for the normalized data!",exit_point=8888)
         ENDIF
         IF (FieldIn1%data_type.NE.ppm_type_real) THEN
            fail("FieldIn1 data type is not correct!",exit_point=8888)
         ENDIF
         IF (FieldIn2%data_type.NE.ppm_type_real) THEN
            fail("FieldIn2 data type is not correct!",exit_point=8888)
         ENDIF
      CASE DEFAULT
         IF (FieldIn1%data_type.NE.FieldIn2%data_type) THEN
            fail("FieldIn1 data type is not the same as FieldIn2 data type!",exit_point=8888)
         ENDIF
         IF (PRESENT(FieldMinVal).OR.PRESENT(FieldMaxVal)) THEN
            stdout("Warning !!! The Min and Max values for FieldIn1 are not computed !!! The returned values are wrong!!!")
         ENDIF
      END SELECT
      8888 CONTINUE
      END SUBROUTINE check
      END SUBROUTINE DTYPE(ppm_rc_CopyImageAndNormalize)