      !-------------------------------------------------------------------------
      !  Subroutine   :                    ppm_rc_median
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
      !  Author           - y.afshar           October   2014
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Subroutine   :                    ppm_rc_median
      !-------------------------------------------------------------------------
      !
      !  Purpose      :
      !
      !
      !  Input        : MaskSize,FieldIn,MeshIn,info,FieldOut,
      !                 boundary_condition,is_normalized
      !
      !  Input/output :
      !
      !  Output       : info  (I) return status. 0 on success.
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

      !TODDO
      !Optimization can be done here to avoid multipl copy when we are passing
      !image indexes

      SUBROUTINE DTYPE(ppm_rc_median)(MaskSize,FieldIn,MeshIn,info, &
      &          FieldOut,boundary_condition,is_normalized)

        IMPLICIT NONE
        !-------------------------------------------------------------------------
        !  Includes
        !-------------------------------------------------------------------------

        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        INTEGER,                       INTENT(IN   ) :: MaskSize

        CLASS(ppm_t_field_),           POINTER       :: FieldIn
        !!! the source filed which is supposed to be

        CLASS(ppm_t_equi_mesh_),       POINTER       :: MeshIn
        !!! source mesh on which fields are descritized

        INTEGER,                       INTENT(  OUT) :: info

        CLASS(ppm_t_field_), OPTIONAL, POINTER       :: FieldOut
        !!! the output filed which is supposed to be

        INTEGER,             OPTIONAL, INTENT(IN   ) :: boundary_condition

        LOGICAL,             OPTIONAL, INTENT(IN   ) :: is_normalized

        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
        CLASS(ppm_t_subpatch_), POINTER :: sbpitr

        REAL(ppm_kind_double) :: t0

#if   __DIME == __2D
        REAL(MK), CONTIGUOUS, DIMENSION(:,:),   POINTER :: DTYPE(wpin_r)
        REAL(MK), CONTIGUOUS, DIMENSION(:,:),   POINTER :: DTYPE(wpout_r)
#elif __DIME == __3D
        REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER :: DTYPE(wpin_r)
        REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER :: DTYPE(wpout_r)
#endif

#if   __DIME == __2D
        INTEGER, CONTIGUOUS, DIMENSION(:,:),   POINTER :: DTYPE(wpin_i)
        INTEGER, CONTIGUOUS, DIMENSION(:,:),   POINTER :: DTYPE(wpout_i)
#elif __DIME == __3D
        INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER :: DTYPE(wpin_i)
        INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER :: DTYPE(wpout_i)
#endif

        INTEGER, DIMENSION(:), POINTER :: nnodes

        CHARACTER(LEN=ppm_char) :: caller = 'ppm_rc_median'

        !-------------------------------------------------------------------------
        !  Initialize
        !-------------------------------------------------------------------------
        CALL substart(caller,t0,info)

        CALL check()

        !currently we suppose that the field has been padded at the borders
        sbpitr => MeshIn%subpatch%begin()
        DO WHILE (ASSOCIATED(sbpitr))
           nnodes => sbpitr%nnodes

           SELECT CASE (FieldIn%data_type)
           CASE (ppm_type_int)
              CALL sbpitr%get_field(FieldIn,DTYPE(wpin_i),info)
              or_fail("Failed to get field wpin_i data.")

              IF (PRESENT(FieldOut)) THEN
                 CALL sbpitr%get_field(FieldOut,DTYPE(wpout_i),info)
                 or_fail("Failed to get field wpout data.")
              ELSE
                 DTYPE(wpout_i) => DTYPE(wpin_i)
              ENDIF

              CALL Update_i(DTYPE(wpin_i),DTYPE(wpout_i),nnodes,MaskSize,info)
              or_fail("Update")

           CASE (ppm_type_real)
              CALL sbpitr%get_field(FieldIn,DTYPE(wpin_r),info)
              or_fail("Failed to get field wpin_i data.")

              IF (PRESENT(FieldOut)) THEN
                 CALL sbpitr%get_field(FieldOut,DTYPE(wpout_r),info)
                 or_fail("Failed to get field wpout data.")
              ELSE
                 DTYPE(wpout_r) => DTYPE(wpin_r)
              ENDIF

              CALL Update_r(DTYPE(wpin_r),DTYPE(wpout_r),nnodes,MaskSize,info)
              or_fail("Update")

           END SELECT

           sbpitr => MeshIn%subpatch%next()
        ENDDO !(ASSOCIATED(sbpitr))

        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
      9999 CONTINUE
        CALL substop(caller,t0,info)
        RETURN
      CONTAINS
      SUBROUTINE check
      IMPLICIT NONE
        check_true(<#ALL(MeshIn%ghostsize(1:__DIME).GE.MaskSize)#>, &
        & "The ghost size is smaller than the mask size, so there would be problem at the borders!",exit_point=8888)

        IF (PRESENT(FieldOut)) THEN
           check_true(<#FieldOut%is_discretized_on(MeshIn)#>, &
           & "The output Field is present, but it has not been descritized on the mesh!",exit_point=8888)

           check_true(<#FieldOut%data_type.EQ.FieldIn%data_type#>, &
           & "The Input and output fields have different data types!",exit_point=8888)
        ENDIF

      8888 CONTINUE
      END SUBROUTINE check
#if   __DIME == __2D
      SUBROUTINE Update_i(wpin,wpout,nnodes_,msize_,info)
        USE ppm_module_util_qsort
        IMPLICIT NONE
        INTEGER, CONTIGUOUS, DIMENSION(:,:), POINTER       :: wpin
        INTEGER, CONTIGUOUS, DIMENSION(:,:), POINTER       :: wpout
        INTEGER,             DIMENSION(:),   INTENT(IN   ) :: nnodes_
        INTEGER,                             INTENT(IN   ) :: msize_
        INTEGER,                             INTENT(  OUT) :: info

        INTEGER, CONTIGUOUS, DIMENSION(:), POINTER :: Mask_
        INTEGER                                    :: i,j,l,m
        INTEGER                                    :: ii,jj,ll

        IF (msize_.EQ.0) THEN
        ELSE
           m=2*msize_+1
           m=m*m

           ALLOCATE(tmp1_i(m),tmp2_i(1:nnodes_(1),1:nnodes_(2)),STAT=info)
           or_fail_alloc("tmp1_i,tmp2_i",exit_point=7777)

           l=2*msize_+1+msize_+1
           !index of the center cell

           NULLIFY(Mask_)

           DO j=1,nnodes_(2)
              DO i=1,nnodes_(1)
                 ll=1
                 DO jj=j-msize_,j+msize_
                    DO ii=i-msize_,i+msize_
                       tmp1_i(ll)=wpin(ii,jj)
                       ll=ll+1
                    ENDDO !ii=i-msize_:i+msize_
                 ENDDO !jj=j-msize_:j+msize_

                 CALL ppm_util_qsort(tmp1_i,Mask_,info,m)

                 tmp2_i(i,j)=tmp1_i(Mask_(l))
              ENDDO !i=1,nnodes_(1)
           ENDDO !j=1,nnodes_(2)

           FORALL (i=1:nnodes_(1),j=1:nnodes_(2))
              wpout(i,j)=tmp2_i(i,j)
           END FORALL

           DEALLOCATE(tmp1_i,tmp2_i,STAT=info)
           or_fail_dealloc("tmp2_i",exit_point=7777)

           CALL ppm_alloc(Mask_,ldc,ppm_param_dealloc,info)
           or_fail_dealloc("Mask_",exit_point=7777)
        ENDIF
      7777 CONTINUE
      END SUBROUTINE Update_i

      SUBROUTINE Update_r(wpin,wpout,nnodes_,msize_,info)
        USE ppm_module_util_qsort
        IMPLICIT NONE
        REAL(MK), CONTIGUOUS, DIMENSION(:,:), POINTER       :: wpin
        REAL(MK), CONTIGUOUS, DIMENSION(:,:), POINTER       :: wpout
        INTEGER,              DIMENSION(:),   INTENT(IN   ) :: nnodes_
        INTEGER,                              INTENT(IN   ) :: msize_
        INTEGER,                              INTENT(  OUT) :: info

        INTEGER, CONTIGUOUS, DIMENSION(:), POINTER :: Mask_
        INTEGER                                    :: i,j,l,m
        INTEGER                                    :: ii,jj,ll

        IF (msize_.EQ.0) THEN
        ELSE
           m=2*msize_+1
           m=m*m

           ALLOCATE(tmp1_r(m),tmp2_r(1:nnodes_(1),1:nnodes_(2)),STAT=info)
           or_fail_alloc("tmp1_i & tmp2_r",exit_point=7777)

           l=2*msize_+1+msize_+1
           !index of the center cell

           NULLIFY(Mask_)

           DO j=1,nnodes_(2)
              DO i=1,nnodes_(1)
                 ll=1
                 DO jj=j-msize_,j+msize_
                    DO ii=i-msize_,i+msize_
                       tmp1_r(ll)=wpin(ii,jj)
                       ll=ll+1
                    ENDDO !ii=i-msize_,i+msize_
                 ENDDO !jj=j-msize_,j+msize_

                 CALL ppm_util_qsort(tmp1_r,Mask_,info,m)

                 tmp2_r(i,j)=tmp1_r(Mask_(l))
              ENDDO !i=1,nnodes_(1)
           ENDDO !j=1,nnodes_(2)

           FORALL (i=1:nnodes_(1),j=1:nnodes_(2))
              wpout(i,j)=tmp2_r(i,j)
           END FORALL

           DEALLOCATE(tmp1_r,tmp2_r,STAT=info)
           or_fail_dealloc("tmp1_r,tmp2_r",exit_point=7777)

           CALL ppm_alloc(Mask_,ldc,ppm_param_dealloc,info)
           or_fail_dealloc("Mask_",exit_point=7777)
        ENDIF
      7777 CONTINUE
      END SUBROUTINE Update_r
#elif __DIME == __3D
      SUBROUTINE Update_i(wpin,wpout,nnodes_,msize_,info)
        USE ppm_module_util_qsort
        IMPLICIT NONE
        INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: wpin
        INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: wpout
        INTEGER,             DIMENSION(:),     INTENT(IN   ) :: nnodes_
        INTEGER,                               INTENT(IN   ) :: msize_
        INTEGER,                               INTENT(  OUT) :: info

        INTEGER, CONTIGUOUS, DIMENSION(:), POINTER :: Mask_
        INTEGER                                    :: i,j,k,l,m
        INTEGER                                    :: ii,jj,kk,ll

        IF (msize_.EQ.0) THEN
        ELSE
           m=2*msize_+1
           m=m*m*m

           ALLOCATE(tmp1_i(m),tmp3_i(1:nnodes_(1),1:nnodes_(2),1:nnodes_(3)),STAT=info)
           or_fail_alloc("tmp1_i & tmp3_i",exit_point=7777)

           l=(2*msize_+1)**2+2*msize_+1+msize_+1
           !middle cell

           NULLIFY(Mask_)

           DO k=1,nnodes_(3)
              DO j=1,nnodes_(2)
                 DO i=1,nnodes_(1)
                    ll=1
                    DO kk=k-msize_,k+msize_
                       DO jj=j-msize_,j+msize_
                          DO ii=i-msize_,i+msize_
                             tmp1_i(ll)=wpin(ii,jj,kk)
                             ll=ll+1
                          ENDDO !ii=i-msize_,i+msize_
                       ENDDO !jj=j-msize_,j+msize_
                    ENDDO !kk=k-msize_,k+msize_

                    CALL ppm_util_qsort(tmp1_i,Mask_,info,m)

                    tmp3_i(i,j,k)=tmp1_i(Mask_(l))
                 ENDDO
              ENDDO
           ENDDO

           FORALL (i=1:nnodes_(1),j=1:nnodes_(2),k=1:nnodes_(3))
              wpout(i,j,k)=tmp3_i(i,j,k)
           END FORALL

           DEALLOCATE(tmp1_i,tmp3_i,STAT=info)
           or_fail_dealloc("tmp1_i & tmp3_i",exit_point=7777)

           CALL ppm_alloc(Mask_,ldc,ppm_param_dealloc,info)
           or_fail_dealloc("Mask_",exit_point=7777)
        ENDIF
      7777 CONTINUE
      END SUBROUTINE Update_i

      SUBROUTINE Update_r(wpin,wpout,nnodes_,msize_,info)
        USE ppm_module_util_qsort
        IMPLICIT NONE
        REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: wpin
        REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: wpout
        INTEGER,              DIMENSION(:),     INTENT(IN   ) :: nnodes_
        INTEGER,                                INTENT(IN   ) :: msize_
        INTEGER,                                INTENT(  OUT) :: info

        INTEGER, CONTIGUOUS, DIMENSION(:), POINTER :: Mask_
        INTEGER                                    :: i,j,k,l,m
        INTEGER                                    :: ii,jj,kk,ll

        IF (msize_.EQ.0) THEN
        ELSE
           m=2*msize_+1
           m=m*m*m

           ALLOCATE(tmp1_r(m),tmp3_r(1:nnodes_(1),1:nnodes_(2),1:nnodes_(3)),STAT=info)
           or_fail_alloc("tmp1_r & tmp3_r",exit_point=7777)

           l=(2*msize_+1)**2+2*msize_+1+msize_+1
           !middle cell

           NULLIFY(Mask_)

           DO k=1,nnodes_(3)
              DO j=1,nnodes_(2)
                 DO i=1,nnodes_(1)
                    ll=1
                    DO kk=k-msize_,k+msize_
                       DO jj=j-msize_,j+msize_
                          DO ii=i-msize_,i+msize_
                             tmp1_r(ll)=wpin(ii,jj,kk)
                             ll=ll+1
                          ENDDO !ii=i-msize_,i+msize_
                       ENDDO !jj=j-msize_,j+msize_
                    ENDDO !kk=k-msize_,k+msize_

                    CALL ppm_util_qsort(tmp1_r,Mask_,info,m)

                    tmp3_r(i,j,k)=tmp1_r(Mask_(l))
                 ENDDO
              ENDDO
           ENDDO

           FORALL (i=1:nnodes_(1),j=1:nnodes_(2),k=1:nnodes_(3))
              wpout(i,j,k)=tmp3_r(i,j,k)
           END FORALL

           DEALLOCATE(tmp1_r,tmp3_r,STAT=info)
           or_fail_dealloc("tmp1_r & tmp3_r",exit_point=7777)

           CALL ppm_alloc(Mask_,ldc,ppm_param_dealloc,info)
           or_fail_dealloc("Mask_",exit_point=7777)
        ENDIF
      7777 CONTINUE
      END SUBROUTINE Update_r
#endif
      END SUBROUTINE DTYPE(ppm_rc_median)

