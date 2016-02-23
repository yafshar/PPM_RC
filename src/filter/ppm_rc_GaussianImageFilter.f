        !-------------------------------------------------------------------------
        !  Subroutine   :                    ppm_rc_GaussianImageFilter
        !-------------------------------------------------------------------------
        !
        !  Purpose      : Convolve image by a Gaussian filter.
        !
        !
        !  Input        : MaskIn, FieldIn, MeshIn, FieldOut
        !
        !  Input/output :
        !
        !  Output       : info  (I) return status. 0 on success.
        !
        !  Routines     :
        !
        !
        !  Remarks      :Lots of room for improvement
        !                e.g. providing blocks for loop through the data
        !
        !                This Routine is working with the same input and output!
        !                It needs an input requested region with the ghost sizes
        !                (larger or equal by the size of the Gaussian kernel)
        !
        !  References   :
        !
        !  Revisions    :
        !-------------------------------------------------------------------------
        !  MOSAIC Group
        !  Max Planck Institute of Molecular Cell Biology and Genetics
        !  Pfotenhauerstr. 108, 01307 Dresden, Germany
        !
        !  Author           - y.afshar           Dec       2015
        !-------------------------------------------------------------------------

        SUBROUTINE DTYPE(ppm_rc_GaussianImageFilter)(G_Sigma, &
        &          FieldIn,MeshIn,info,FieldOut,KernelFactor)

          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------

          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          REAL(MK), DIMENSION(:),        INTENT(IN   ) :: G_Sigma
          !!! Standard deviation Sigma for Gaussian convolution

          CLASS(ppm_t_field_),           POINTER       :: FieldIn
          !!! the source filed which is supposed to be

          CLASS(ppm_t_equi_mesh_),       POINTER       :: MeshIn
          !!! source mesh on which fields are descritized

          INTEGER,                       INTENT(  OUT) :: info

          CLASS(ppm_t_field_), OPTIONAL, POINTER       :: FieldOut
          !!! the output filed which is supposed to be

          REAL(MK),            OPTIONAL, INTENT(IN   ) :: KernelFactor
          !!! factor for the Gaussian kernel
          !!! The convolution mask will have the size of : FLOOR(KernelFactor*G_Sigma)
          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          CLASS(ppm_t_subpatch_), POINTER :: sbpitr

          REAL(ppm_kind_double)                           :: t0
#if   __DIME == __2D
          REAL(MK), CONTIGUOUS, DIMENSION(:,:),   POINTER :: DTYPE(wpin)
          REAL(MK), CONTIGUOUS, DIMENSION(:,:),   POINTER :: DTYPE(wpout)
#elif __DIME == __3D
          REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER :: DTYPE(wpin)
          REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER :: DTYPE(wpout)
#endif
          REAL(MK),             DIMENSION(:), ALLOCATABLE :: Mask
          REAL(MK)                                        :: coef,coef1,sum1

          INTEGER, DIMENSION(:), POINTER :: Nm
          INTEGER, DIMENSION(:), POINTER :: lo_a
          INTEGER, DIMENSION(:), POINTER :: hi_a
          INTEGER, DIMENSION(__DIME)          :: ghostsize_
          INTEGER                        :: i,j,s1,s2
#if   __DIME == __3D
          INTEGER                        :: k,s3
#endif

          CHARACTER(LEN=ppm_char) :: caller = 'ppm_rc_GaussianImageFilter'

          LOGICAL :: lsymmetric
          LOGICAL :: is_discretized_on
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          !Get the G_Sigma vector dimension. It should be 1 for
          !symmetric Gaussian and 2 or 3 for Anisotropic Gaussian
          s1=SIZE(G_Sigma,1)

          SELECT CASE (s1)
          CASE (1)
          !!! isotropic (i.e. circularly symmetric) Gaussian
             lsymmetric=.TRUE.
          CASE (2)
             IF (INT(G_Sigma(1)).EQ.INT(G_Sigma(2))) THEN
                lsymmetric=.TRUE.
             ELSE
                lsymmetric=.FALSE.
             ENDIF
          !!! Anisotropic 2D Gaussian
          CASE (3)
             IF (INT(G_Sigma(1)).EQ.INT(G_Sigma(2)).AND.INT(G_Sigma(1)).EQ.INT(G_Sigma(3))) THEN
                lsymmetric=.TRUE.
             ELSE
                lsymmetric=.FALSE.
             ENDIF
             !!! Anisotropic 3D Gaussian
          CASE DEFAULT
             fail("Wrong dimension for GaussianImageFilter G_Sigma!")
          END SELECT

          !!! We are considering two times G_Sigma as a default
          !!! kernel size for convolution
          coef= MERGE(KernelFactor,two,PRESENT(KernelFactor))

          CALL check()

          !!! TODO
          !!! I should change this routine find the index if ghostsize>FLOOR(coef*G_Sigma)
          !!! we are safe other wise compute the new sigmas that the whole image should
          !!! be convolved with to produce the original sigma effect
          !!! nsigma=ghostsize/coef ..... X*nsigma1^2+nsigma2^2=G_Sigma^2 then convolve
          !!! the image X times with sigma1 and one time with sigma2

          IF (ppm_nproc.GT.1) THEN
             !!! calculate the required ghost size for X direction
             ghostsize_(1)=FLOOR(coef*G_Sigma(1))
             ghostsize_(2)=0
             ghostsize_(__DIME)=0

             !-------------------------------------------------------------------------
             !  Get field ghosts for image
             !-------------------------------------------------------------------------
             CALL MeshIn%map_ghost_get(info,ghostsize=ghostsize_)
             or_fail("MeshIn%map_ghost_get")

             CALL FieldIn%map_ghost_push(MeshIn,info)
             or_fail("FieldIn%map_ghost_push")
             CALL MeshIn%map_isend(info,sendrecv=.TRUE.)
             or_fail("MeshIn%map_isend")

             s1=ghostsize_(1)
          ELSE
             s1=FLOOR(coef*G_Sigma(1))
          ENDIF

          IF (.NOT.is_discretized_on) THEN
             IF (.NOT.ASSOCIATED(FieldOut%discr_info)) THEN
                CALL FieldOut%create(1,info,dtype=FieldIn%data_type,name="blured_image")
                or_fail("Failed to create field!" )
             ENDIF

             CALL FieldOut%discretize_on(MeshIn,info)
             or_fail("Failed to discretize Field on Mesh!")
          ENDIF

          IF (.NOT.lsymmetric) THEN
             s2=FLOOR(coef*G_Sigma(2))
#if   __DIME == __3D
             s3=FLOOR(coef*G_Sigma(3))
#endif
          ENDIF

          !-------------------------------------------------------------------------
          !  Pads the input at the boarder of the global domain
          !  pad_size: how many times the most outer plane is copied
          !  (~ ghost layer size)
          !-------------------------------------------------------------------------
          NULLIFY(DTYPE(wpin))
          sbpitr => MeshIn%subpatch%begin()
          DO WHILE (ASSOCIATED(sbpitr))
             Nm   => sbpitr%nnodes
             hi_a => sbpitr%hi_a
             lo_a => sbpitr%lo_a

             CALL sbpitr%get_field(FieldIn,DTYPE(wpin),info)
             or_fail("Failed to get field wpin data.")

#if __DIME == __2D
             !padding the left yz image face
             IF (sbpitr%istart(1).EQ.1) THEN
                DO j=1,Nm(2)
                   DO i=lo_a(1),0
                      DTYPE(wpin)(i,j)=DTYPE(wpin)(1,j)
                   ENDDO !i=lo_a(1),0
                ENDDO !j=lo_a(2),hi_a(2)
             ENDIF !(sbpitr%istart(1).EQ.1)
             !padding the right yz image face
             IF (sbpitr%iend(1).EQ.Ngrid(1)) THEN
                DO j=1,Nm(2)
                   DO i=Nm(1)+1,sbpitr%hi_a(1)
                      DTYPE(wpin)(i,j)=DTYPE(wpin)(Nm(1),j)
                   ENDDO !i=Nm(1)+1,sbpitr%hi_a(1)
                ENDDO !j=lo_a(2),hi_a(2)
             ENDIF !(sbpitr%iend(1).EQ.Ngrid(1))
#elif __DIME == __3D
             !padding the left yz image face
             IF (sbpitr%istart(1).EQ.1) THEN
                DO k=1,Nm(3)
                   DO j=1,Nm(2)
                      DO i=lo_a(1),0
                         DTYPE(wpin)(i,j,k)= DTYPE(wpin)(1,j,k)
                      ENDDO !i=lo_a(1),0
                   ENDDO !j=1,Nm(2)
                ENDDO !k=1,Nm(3)
             ENDIF !(sbpitr%istart(1).EQ.1)
             !padding the right yz image face
             IF (sbpitr%iend(1).EQ.Ngrid(1)) THEN
                DO k=1,Nm(3)
                   DO j=1,Nm(2)
                      DO i=Nm(1)+1,sbpitr%hi_a(1)
                         DTYPE(wpin)(i,j,k)= DTYPE(wpin)(Nm(1),j,k)
                      ENDDO !i=Nm(1)+1,sbpitr%hi_a(1)
                   ENDDO !j=1,Nm(2)
                ENDDO !k=1,Nm(3)
             ENDIF !(sbpitr%iend(1).EQ.Ngrid(1))
#endif
             sbpitr => MeshIn%subpatch%next()
          ENDDO
          NULLIFY(DTYPE(wpin))

          ALLOCATE(Mask(-s1:s1),STAT=info)
          or_fail_alloc("Failed to allocate Mask")

          !Gaussians = 1/(sqrt(2*pi)*sigma).*exp(-x.^2/(2*sigma.^2))
          coef =one/(SQRT(two*pi)*G_Sigma(1))
          coef1=one/(two*G_Sigma(1)*G_Sigma(1))
          FORALL (i=-s1:s1) Mask(-i)= coef * EXP(-REAL(i*i,MK)*coef1)
          sum1=SUM(Mask)
          Mask=Mask/sum1

          IF (ppm_nproc.GT.1) THEN
             CALL MeshIn%map_isend(info,sendrecv=.FALSE.)
             or_fail("MeshIn%map_isend")
             CALL FieldIn%map_ghost_pop(MeshIn,info)
             or_fail("FieldIn%map_ghost_pop")
          ENDIF

          !The field padded at the borders in the X direction
          !and the ghost is received
          sbpitr => MeshIn%subpatch%begin()
          DO WHILE (ASSOCIATED(sbpitr))
             Nm => sbpitr%nnodes

             CALL sbpitr%get_field(FieldIn,DTYPE(wpin),info)
             or_fail("Failed to get field wpin data.")

             IF (PRESENT(FieldOut)) THEN
                NULLIFY(DTYPE(wpout))
                CALL sbpitr%get_field(FieldOut,DTYPE(wpout),info)
                or_fail("Failed to get field wpout data.")
             ELSE
                DTYPE(wpout) => DTYPE(wpin)
             ENDIF

             !!! Gaussian convolution in X direction
             CALL Update(Mask,s1,DTYPE(wpin),DTYPE(wpout),Nm,1,info)
             or_fail("Update")

             sbpitr => MeshIn%subpatch%next()
          ENDDO !(ASSOCIATED(sbpitr))
          NULLIFY(DTYPE(wpin))

          IF (ppm_nproc.GT.1) THEN
             !!! calculate the required ghost size for Y direction
             ghostsize_(1)=0
             ghostsize_(2)=MERGE(s1,s2,lsymmetric)
#if   __DIME == __3D
             ghostsize_(3)=0
#endif

             !-------------------------------------------------------------------------
             !  Get field ghosts for image
             !-------------------------------------------------------------------------
             CALL MeshIn%map_ghost_get(info,ghostsize=ghostsize_)
             or_fail("MeshIn%map_ghost_get")

             IF (PRESENT(FieldOut)) THEN
                CALL FieldOut%map_ghost_push(MeshIn,info)
                or_fail("FieldOut%map_ghost_push")
             ELSE
                CALL FieldIn%map_ghost_push(MeshIn,info)
                or_fail("FieldIn%map_ghost_push")
             ENDIF
             CALL MeshIn%map_isend(info,sendrecv=.TRUE.)
             or_fail("MeshIn%map_isend")
          ENDIF

          !-------------------------------------------------------------------------
          !  Pads the input at the boarder of the global domain
          !  pad_size: how many times the most outer plane is copied
          !  (~ ghost layer size)
          !-------------------------------------------------------------------------
          sbpitr => MeshIn%subpatch%begin()
          DO WHILE (ASSOCIATED(sbpitr))
             Nm   => sbpitr%nnodes
             hi_a => sbpitr%hi_a
             lo_a => sbpitr%lo_a

             IF (PRESENT(FieldOut)) THEN
                CALL sbpitr%get_field(FieldOut,DTYPE(wpin),info)
             ELSE
                CALL sbpitr%get_field(FieldIn,DTYPE(wpin),info)
             ENDIF
             or_fail("Failed to get field wpin data.")

#if __DIME == __2D
             !padding the left xz image face
             IF (sbpitr%istart(2).EQ.1) THEN
                DO j=lo_a(2),0
                   DO i=1,Nm(1)
                      DTYPE(wpin)(i,j)=DTYPE(wpin)(i,1)
                   ENDDO !i=1,Nm(1)
                ENDDO !j=lo_a(2),0
             ENDIF !(sbpitr%istart(2).EQ.1)
             !padding the right xz image face
             IF (sbpitr%iend(2).EQ.Ngrid(2)) THEN
                DO j=Nm(2)+1,hi_a(2)
                   DO i=1,Nm(1)
                      DTYPE(wpin)(i,j)= DTYPE(wpin)(i,Nm(2))
                   ENDDO !i=1,Nm(1)
                ENDDO !Nm(2)+1,hi_a(2)
             ENDIF !(sbpitr%iend(2).EQ.Ngrid(2))
#elif __DIME == __3D
             !padding the left xz image face
             IF (sbpitr%istart(2).EQ.1) THEN
                DO k=1,Nm(3)
                   DO j=lo_a(2),0
                      DO i=1,Nm(1)
                         DTYPE(wpin)(i,j,k)=DTYPE(wpin)(i,1,k)
                      ENDDO !i=1,Nm(1)
                   ENDDO !j=lo_a(2),0
                ENDDO !k=1,Nm(3)
             ENDIF !(sbpitr%istart(2).EQ.1)
             !padding the right xz image face
             IF (sbpitr%iend(2).EQ.Ngrid(2)) THEN
                DO k=1,Nm(3)
                   DO j=Nm(2)+1,hi_a(2)
                      DO i=1,Nm(1)
                         DTYPE(wpin)(i,j,k)= DTYPE(wpin)(i,Nm(2),k)
                      ENDDO !i=1,Nm(1)
                   ENDDO !j=Nm(2)+1,hi_a(2)
                ENDDO !k=1,Nm(3)
             ENDIF !(sbpitr%iend(2).EQ.Ngrid(2))
#endif

             sbpitr => MeshIn%subpatch%next()
          ENDDO
          NULLIFY(DTYPE(wpin))

          IF (lsymmetric) THEN
             s2=s1
          ELSE
             DEALLOCATE(Mask,STAT=info)
             or_fail_dealloc("Failed to deallocate Mask")

             ALLOCATE(Mask(-s2:s2),STAT=info)
             or_fail_alloc("Failed to allocate Mask")

             !Gaussians = 1/(sqrt(2*pi)*sigma).*exp(-y.^2/(2*sigma.^2))
             coef =one/(SQRT(two*pi)*G_Sigma(2))
             coef1=one/(two*G_Sigma(2)*G_Sigma(2))
             FORALL (i=-s2:s2) Mask(-i)= coef * EXP(-REAL(i*i,MK)*coef1)
             sum1=SUM(Mask)
             Mask=Mask/sum1
          ENDIF

          IF (ppm_nproc.GT.1) THEN
             CALL MeshIn%map_isend(info,sendrecv=.FALSE.)
             or_fail("MeshIn%map_isend")

             IF (PRESENT(FieldOut)) THEN
                CALL FieldOut%map_ghost_pop(MeshIn,info)
                or_fail("FieldOut%map_ghost_pop")
             ELSE
                CALL FieldIn%map_ghost_pop(MeshIn,info)
                or_fail("FieldIn%map_ghost_pop")
             ENDIF
          ENDIF

          !The field padded at the borders in the X direction
          sbpitr => MeshIn%subpatch%begin()
          DO WHILE (ASSOCIATED(sbpitr))
             Nm => sbpitr%nnodes

             IF (PRESENT(FieldOut)) THEN
                CALL sbpitr%get_field(FieldOut,DTYPE(wpin),info)
             ELSE
                CALL sbpitr%get_field(FieldIn,DTYPE(wpin),info)
             ENDIF
             or_fail("Failed to get field wpin data.")

             DTYPE(wpout) => DTYPE(wpin)

             CALL Update(Mask,s2,DTYPE(wpin),DTYPE(wpout),Nm,2,info)
             or_fail("Update")

             sbpitr => MeshIn%subpatch%next()
          ENDDO !(ASSOCIATED(sbpitr))
          NULLIFY(DTYPE(wpin))

#if __DIME == __3D
          IF (ppm_nproc.GT.1) THEN
             !!! calculate the required ghost size for Z direction
             ghostsize_(1)=0
             ghostsize_(2)=0
             ghostsize_(3)=MERGE(s1,s3,lsymmetric)

             !-------------------------------------------------------------------------
             !  Get field ghosts for image
             !-------------------------------------------------------------------------
             CALL MeshIn%map_ghost_get(info,ghostsize=ghostsize_)
             or_fail("MeshIn%map_ghost_get")

             IF (PRESENT(FieldOut)) THEN
                CALL FieldOut%map_ghost_push(MeshIn,info)
                or_fail("FieldOut%map_ghost_push")
             ELSE
                CALL FieldIn%map_ghost_push(MeshIn,info)
                or_fail("FieldIn%map_ghost_push")
             ENDIF
             CALL MeshIn%map_isend(info,sendrecv=.TRUE.)
             or_fail("MeshIn%map_isend")
          ENDIF

          !-------------------------------------------------------------------------
          !  Pads the input at the boarder of the global domain
          !  pad_size: how many times the most outer plane is copied
          !  (~ ghost layer size)
          !-------------------------------------------------------------------------
          sbpitr => MeshIn%subpatch%begin()
          DO WHILE (ASSOCIATED(sbpitr))
             Nm   => sbpitr%nnodes
             hi_a => sbpitr%hi_a
             lo_a => sbpitr%lo_a

             IF (PRESENT(FieldOut)) THEN
                CALL sbpitr%get_field(FieldOut,DTYPE(wpin),info)
             ELSE
                CALL sbpitr%get_field(FieldIn,DTYPE(wpin),info)
             ENDIF
             or_fail("Failed to get field wpin data.")

             !padding the bottom xy image face
             IF (sbpitr%istart(3).EQ.1) THEN
                DO k=lo_a(3),0
                   DO j=1,Nm(2)
                      DO i=1,Nm(1)
                         DTYPE(wpin)(i,j,k)=DTYPE(wpin)(i,j,1)
                      ENDDO !i=1,Nm(1)
                   ENDDO !j=1,Nm(2)
                ENDDO !k=lo_a(3),0
             ENDIF !(sbpitr%istart(3).EQ.1)
             !padding the top xy image face
             IF (sbpitr%iend(3).EQ.Ngrid(3)) THEN
                DO k=Nm(3)+1,hi_a(3)
                   DO j=1,Nm(2)
                      DO i=1,Nm(1)
                         DTYPE(wpin)(i,j,k)= DTYPE(wpin)(i,j,Nm(3))
                      ENDDO !i=1,Nm(1)
                   ENDDO !j=1,Nm(2)
                ENDDO !k=Nm(3)+1,hi_a(3)
             ENDIF !(sbpitr%iend(3).EQ.Ngrid(3))

             sbpitr => MeshIn%subpatch%next()
          ENDDO
          NULLIFY(DTYPE(wpin))

          IF (lsymmetric) THEN
             s3=s1
          ELSE
             DEALLOCATE(Mask,STAT=info)
             or_fail_dealloc("Failed to deallocate Mask")

             ALLOCATE(Mask(-s3:s3),STAT=info)
             or_fail_alloc("Failed to allocate Mask")

             !Gaussians = 1/(sqrt(2*pi)*sigma).*exp(-z.^2/(2*sigma.^2))
             coef =one/(SQRT(two*pi)*G_Sigma(3))
             coef1=one/(two*G_Sigma(3)*G_Sigma(3))
             FORALL (i=-s3:s3) Mask(-i)= coef * EXP(-REAL(i*i,MK)*coef1)
             sum1=SUM(Mask)
             Mask=Mask/sum1
          ENDIF

          IF (ppm_nproc.GT.1) THEN
             CALL MeshIn%map_isend(info,sendrecv=.FALSE.)
             or_fail("MeshIn%map_isend")

             IF (PRESENT(FieldOut)) THEN
                CALL FieldOut%map_ghost_pop(MeshIn,info)
                or_fail("FieldOut%map_ghost_pop")
             ELSE
                CALL FieldIn%map_ghost_pop(MeshIn,info)
                or_fail("FieldIn%map_ghost_pop")
             ENDIF
          ENDIF

          !The field padded at the borders in the X direction
          sbpitr => MeshIn%subpatch%begin()
          DO WHILE (ASSOCIATED(sbpitr))
             Nm => sbpitr%nnodes

             IF (PRESENT(FieldOut)) THEN
                CALL sbpitr%get_field(FieldOut,DTYPE(wpin),info)
             ELSE
                CALL sbpitr%get_field(FieldIn,DTYPE(wpin),info)
             ENDIF
             or_fail("Failed to get field wpin data.")

             DTYPE(wpout) => DTYPE(wpin)

             CALL Update(Mask,s3,DTYPE(wpin),DTYPE(wpout),Nm,3,info)
             or_fail("Update")

             sbpitr => MeshIn%subpatch%next()
          ENDDO !(ASSOCIATED(sbpitr))
          NULLIFY(DTYPE(wpin))
#endif

          DEALLOCATE(Mask,STAT=info)
          or_fail_dealloc("Failed to deallocate Mask")

          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        CONTAINS
        SUBROUTINE check
        IMPLICIT NONE
          ghostsize_=FLOOR(coef*G_Sigma)

          check_true(<#ALL(MeshIn%ghostsize.GE.ghostsize_)#>, &
          & "The ghost size is smaller than the kernel size, so there would be problem at the borders!",exit_point=8888)

          IF (PRESENT(FieldOut)) THEN
             IF (FieldOut%is_discretized_on(MeshIn)) THEN
                is_discretized_on=.TRUE.
                check_true(<#FieldOut%data_type.EQ.FieldIn%data_type#>, &
                & "The Input and output fields have different data types!",exit_point=8888)
             ELSE
                is_discretized_on=.FALSE.
             ENDIF
          ENDIF
        8888 CONTINUE
        END SUBROUTINE check
#if   __DIME == __2D
        SUBROUTINE Update(Mask,MaskSize,wpin,wpout,nnodes,direction,info)
          IMPLICIT NONE
          REAL(MK),             DIMENSION(:),   INTENT(IN   ) :: Mask
          INTEGER,                              INTENT(IN   ) :: MaskSize
          REAL(MK), CONTIGUOUS, DIMENSION(:,:), POINTER       :: wpin
          REAL(MK), CONTIGUOUS, DIMENSION(:,:), POINTER       :: wpout
          INTEGER,              DIMENSION(:),   INTENT(IN   ) :: nnodes
          INTEGER,                              INTENT(IN   ) :: direction
          INTEGER,                              INTENT(  OUT) :: info

          INTEGER                            :: l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11
          INTEGER                            :: l12,l13,l14,l15,l16,l17,l18,l19,l20,l21
          INTEGER                            :: i,j,l,m,n
          INTEGER                            :: jp
          INTEGER, DIMENSION(:), ALLOCATABLE :: ll

          SELECT CASE (direction)
          CASE (1)
          !!! X direction
             ALLOCATE(tmp1_r(nnodes(1)),STAT=info)
             or_fail_alloc("Failed to allocate tmp1_r!",exit_point=7777)

             SELECT CASE (MaskSize)
             CASE (1)
                DO j=1,nnodes(2)
                   DO i=1,nnodes(1)
                      tmp1_r(i)=wpin(i-1,j)*Mask(1)+wpin(i,j)*Mask(2)+wpin(i+1,j)*Mask(3)
                   ENDDO !i=1,Nm(1)
                   FORALL (i=1:nnodes(1))
                      wpout(i,j)=tmp1_r(i)
                   END FORALL
                ENDDO !j=ldl(2),ldu(2)
             CASE (2)
                DO j=1,nnodes(2)
                   DO i=1,nnodes(1)
                      tmp1_r(i)=wpin(i-2,j)*Mask(1)+wpin(i-1,j)*Mask(2)+wpin(i,j)*Mask(3)+wpin(i+1,j)*Mask(4)+ &
                      &         wpin(i+2,j)*Mask(5)
                   ENDDO !i=1,Nm(1)
                   FORALL (i=1:nnodes(1))
                      wpout(i,j)=tmp1_r(i)
                   END FORALL
                ENDDO !j=ldl(2),ldu(2)
             CASE (3)
                DO j=1,nnodes(2)
                   DO i=1,nnodes(1)
                      tmp1_r(i)=wpin(i-3,j)*Mask(1)+wpin(i-2,j)*Mask(2)+wpin(i-1,j)*Mask(3)+wpin(i,j)*Mask(4)+  &
                      &         wpin(i+1,j)*Mask(5)+wpin(i+2,j)*Mask(6)+wpin(i+3,j)*Mask(7)
                   ENDDO !i=1,Nm(1)
                   FORALL (i=1:nnodes(1))
                      wpout(i,j)=tmp1_r(i)
                   END FORALL
                ENDDO !j=ldl(2),ldu(2)
             CASE (4)
                DO j=1,nnodes(2)
                   DO i=1,nnodes(1)
                      tmp1_r(i)=wpin(i-4,j)*Mask(1)+wpin(i-3,j)*Mask(2)+wpin(i-2,j)*Mask(3)+wpin(i-1,j)*Mask(4)+ &
                      &         wpin(i  ,j)*Mask(5)+wpin(i+1,j)*Mask(6)+wpin(i+2,j)*Mask(7)+wpin(i+3,j)*Mask(8)+ &
                      &         wpin(i+4,j)*Mask(9)
                   ENDDO !i=1,Nm(1)
                   FORALL (i=1:nnodes(1))
                      wpout(i,j)=tmp1_r(i)
                   END FORALL
                ENDDO !j=ldl(2),ldu(2)
             CASE (5)
                DO j=1,nnodes(2)
                   DO i=1,nnodes(1)
                      tmp1_r(i)=wpin(i-5,j)*Mask(1)+wpin(i-4,j)*Mask( 2)+wpin(i-3,j)*Mask(3)+wpin(i-2,j)*Mask(4)+ &
                      &         wpin(i-1,j)*Mask(5)+wpin(i  ,j)*Mask( 6)+wpin(i+1,j)*Mask(7)+wpin(i+2,j)*Mask(8)+ &
                      &         wpin(i+3,j)*Mask(9)+wpin(i+4,j)*Mask(10)+wpin(i+5,j)*Mask(11)
                   ENDDO !i=1,Nm(1)
                   FORALL (i=1:nnodes(1))
                      wpout(i,j)=tmp1_r(i)
                   END FORALL
                ENDDO !j=ldl(2),ldu(2)
             CASE (6)
                DO j=1,nnodes(2)
                   DO i=1,nnodes(1)
                      tmp1_r(i)=wpin(i-6,j)*Mask( 1)+wpin(i-5,j)*Mask( 2)+wpin(i-4,j)*Mask( 3)+wpin(i-3,j)*Mask( 4)+ &
                      &         wpin(i-2,j)*Mask( 5)+wpin(i-1,j)*Mask( 6)+wpin(i  ,j)*Mask( 7)+wpin(i+1,j)*Mask( 8)+ &
                      &         wpin(i+2,j)*Mask( 9)+wpin(i+3,j)*Mask(10)+wpin(i+4,j)*Mask(11)+wpin(i+5,j)*Mask(12)+ &
                      &         wpin(i+6,j)*Mask(13)
                   ENDDO !i=1,Nm(1)
                   FORALL (i=1:nnodes(1))
                      wpout(i,j)=tmp1_r(i)
                   END FORALL
                ENDDO !j=ldl(2),ldu(2)
             CASE (7)
                DO j=1,nnodes(2)
                   DO i=1,nnodes(1)
                      tmp1_r(i)=wpin(i-7,j)*Mask( 1)+wpin(i-6,j)*Mask( 2)+wpin(i-5,j)*Mask( 3)+wpin(i-4,j)*Mask( 4)+ &
                      &         wpin(i-3,j)*Mask( 5)+wpin(i-2,j)*Mask( 6)+wpin(i-1,j)*Mask( 7)+wpin(i  ,j)*Mask( 8)+ &
                      &         wpin(i+1,j)*Mask( 9)+wpin(i+2,j)*Mask(10)+wpin(i+3,j)*Mask(11)+wpin(i+4,j)*Mask(12)+ &
                      &         wpin(i+5,j)*Mask(13)+wpin(i+6,j)*Mask(14)+wpin(i+7,j)*Mask(15)
                   ENDDO !i=1,Nm(1)
                   FORALL (i=1:nnodes(1))
                      wpout(i,j)=tmp1_r(i)
                   END FORALL
                ENDDO !j=ldl(2),ldu(2)
             CASE (8)
                DO j=1,nnodes(2)
                   DO i=1,nnodes(1)
                      tmp1_r(i)=wpin(i-8,j)*Mask( 1)+wpin(i-7,j)*Mask( 2)+wpin(i-6,j)*Mask( 3)+wpin(i-5,j)*Mask( 4)+ &
                      &         wpin(i-4,j)*Mask( 5)+wpin(i-3,j)*Mask( 6)+wpin(i-2,j)*Mask( 7)+wpin(i-1,j)*Mask( 8)+ &
                      &         wpin(i  ,j)*Mask( 9)+wpin(i+1,j)*Mask(10)+wpin(i+2,j)*Mask(11)+wpin(i+3,j)*Mask(12)+ &
                      &         wpin(i+4,j)*Mask(13)+wpin(i+5,j)*Mask(14)+wpin(i+6,j)*Mask(15)+wpin(i+7,j)*Mask(16)+ &
                      &         wpin(i+8,j)*Mask(17)
                   ENDDO !i=1,Nm(1)
                   FORALL (i=1:nnodes(1))
                      wpout(i,j)=tmp1_r(i)
                   END FORALL
                ENDDO !j=ldl(2),ldu(2)
             CASE (9)
                DO j=1,nnodes(2)
                   DO i=1,nnodes(1)
                      tmp1_r(i)=wpin(i-9,j)*Mask( 1)+wpin(i-8,j)*Mask( 2)+wpin(i-7,j)*Mask( 3)+wpin(i-6,j)*Mask( 4)+ &
                      &         wpin(i-5,j)*Mask( 5)+wpin(i-4,j)*Mask( 6)+wpin(i-3,j)*Mask( 7)+wpin(i-2,j)*Mask( 8)+ &
                      &         wpin(i-1,j)*Mask( 9)+wpin(i  ,j)*Mask(10)+wpin(i+1,j)*Mask(11)+wpin(i+2,j)*Mask(12)+ &
                      &         wpin(i+3,j)*Mask(13)+wpin(i+4,j)*Mask(14)+wpin(i+5,j)*Mask(15)+wpin(i+6,j)*Mask(16)+ &
                      &         wpin(i+7,j)*Mask(17)+wpin(i+8,j)*Mask(18)+wpin(i+9,j)*Mask(19)
                   ENDDO !i=1,Nm(1)
                   FORALL (i=1:nnodes(1))
                      wpout(i,j)=tmp1_r(i)
                   END FORALL
                ENDDO !j=ldl(2),ldu(2)
             CASE (10)
                DO j=1,nnodes(2)
                   DO i=1,nnodes(1)
                      tmp1_r(i)=wpin(i-10,j)*Mask( 1)+wpin(i-9,j)*Mask( 2)+wpin(i-8,j)*Mask( 3)+wpin(i-7,j)*Mask( 4)+ &
                      &         wpin(i- 6,j)*Mask( 5)+wpin(i-5,j)*Mask( 6)+wpin(i-4,j)*Mask( 7)+wpin(i-3,j)*Mask( 8)+ &
                      &         wpin(i- 2,j)*Mask( 9)+wpin(i-1,j)*Mask(10)+wpin(i  ,j)*Mask(11)+wpin(i+1,j)*Mask(12)+ &
                      &         wpin(i+ 2,j)*Mask(13)+wpin(i+3,j)*Mask(14)+wpin(i+4,j)*Mask(15)+wpin(i+5,j)*Mask(16)+ &
                      &         wpin(i+ 6,j)*Mask(17)+wpin(i+7,j)*Mask(18)+wpin(i+8,j)*Mask(19)+wpin(i+9,j)*Mask(20)+ &
                      &         wpin(i+10,j)*Mask(21)
                   ENDDO !i=1,Nm(1)
                   FORALL (i=1:nnodes(1))
                      wpout(i,j)=tmp1_r(i)
                   END FORALL
                ENDDO !j=ldl(2),ldu(2)
             CASE (11)
                DO j=1,nnodes(2)
                   DO i=1,nnodes(1)
                      tmp1_r(i)=wpin(i-11,j)*Mask( 1)+wpin(i-10,j)*Mask( 2)+wpin(i- 9,j)*Mask( 3)+wpin(i-8,j)*Mask( 4)+ &
                      &         wpin(i- 7,j)*Mask( 5)+wpin(i- 6,j)*Mask( 6)+wpin(i- 5,j)*Mask( 7)+wpin(i-4,j)*Mask( 8)+ &
                      &         wpin(i- 3,j)*Mask( 9)+wpin(i- 2,j)*Mask(10)+wpin(i- 1,j)*Mask(11)+wpin(i  ,j)*Mask(12)+ &
                      &         wpin(i+ 1,j)*Mask(13)+wpin(i+ 2,j)*Mask(14)+wpin(i+ 3,j)*Mask(15)+wpin(i+4,j)*Mask(16)+ &
                      &         wpin(i+ 5,j)*Mask(17)+wpin(i+ 6,j)*Mask(18)+wpin(i+ 7,j)*Mask(19)+wpin(i+8,j)*Mask(20)+ &
                      &         wpin(i+ 9,j)*Mask(21)+wpin(i+10,j)*Mask(22)+wpin(i+11,j)*Mask(23)
                   ENDDO !i=1,Nm(1)
                   FORALL (i=1:nnodes(1))
                      wpout(i,j)=tmp1_r(i)
                   END FORALL
                ENDDO !j=ldl(2),ldu(2)
             CASE (12)
                DO j=1,nnodes(2)
                   DO i=1,nnodes(1)
                      tmp1_r(i)=wpin(i-12,j)*Mask( 1)+wpin(i-11,j)*Mask( 2)+wpin(i-10,j)*Mask( 3)+wpin(i- 9,j)*Mask( 4)+ &
                      &         wpin(i- 8,j)*Mask( 5)+wpin(i- 7,j)*Mask( 6)+wpin(i- 6,j)*Mask( 7)+wpin(i- 5,j)*Mask( 8)+ &
                      &         wpin(i- 4,j)*Mask( 9)+wpin(i- 3,j)*Mask(10)+wpin(i- 2,j)*Mask(11)+wpin(i- 1,j)*Mask(12)+ &
                      &         wpin(i   ,j)*Mask(13)+wpin(i+ 1,j)*Mask(14)+wpin(i+ 2,j)*Mask(15)+wpin(i+ 3,j)*Mask(16)+ &
                      &         wpin(i+ 4,j)*Mask(17)+wpin(i+ 5,j)*Mask(18)+wpin(i+ 6,j)*Mask(19)+wpin(i+ 7,j)*Mask(20)+ &
                      &         wpin(i+ 8,j)*Mask(21)+wpin(i+ 9,j)*Mask(22)+wpin(i+10,j)*Mask(23)+wpin(i+11,j)*Mask(24)+ &
                      &         wpin(i+12,j)*Mask(25)
                   ENDDO !i=1,Nm(1)
                   FORALL (i=1:nnodes(1))
                      wpout(i,j)=tmp1_r(i)
                   END FORALL
                ENDDO !j=ldl(2),ldu(2)
             CASE (13)
                DO j=1,nnodes(2)
                   DO i=1,nnodes(1)
                      tmp1_r(i)=wpin(i-13,j)*Mask( 1)+wpin(i-12,j)*Mask( 2)+wpin(i-11,j)*Mask( 3)+wpin(i-10,j)*Mask( 4)+ &
                      &         wpin(i- 9,j)*Mask( 5)+wpin(i- 8,j)*Mask( 6)+wpin(i- 7,j)*Mask( 7)+wpin(i- 6,j)*Mask( 8)+ &
                      &         wpin(i- 5,j)*Mask( 9)+wpin(i- 4,j)*Mask(10)+wpin(i- 3,j)*Mask(11)+wpin(i- 2,j)*Mask(12)+ &
                      &         wpin(i- 1,j)*Mask(13)+wpin(i   ,j)*Mask(14)+wpin(i+ 1,j)*Mask(15)+wpin(i+ 2,j)*Mask(16)+ &
                      &         wpin(i+ 3,j)*Mask(17)+wpin(i+ 4,j)*Mask(18)+wpin(i+ 5,j)*Mask(19)+wpin(i+ 6,j)*Mask(20)+ &
                      &         wpin(i+ 7,j)*Mask(21)+wpin(i+ 8,j)*Mask(22)+wpin(i+ 9,j)*Mask(23)+wpin(i+10,j)*Mask(24)+ &
                      &         wpin(i+11,j)*Mask(25)+wpin(i+12,j)*Mask(26)+wpin(i+13,j)*Mask(27)
                   ENDDO !i=1,Nm(1)
                   FORALL (i=1:nnodes(1))
                      wpout(i,j)=tmp1_r(i)
                   END FORALL
                ENDDO !j=ldl(2),ldu(2)
             CASE (14)
                DO j=1,nnodes(2)
                   DO i=1,nnodes(1)
                      tmp1_r(i)=wpin(i-14,j)*Mask( 1)+wpin(i-13,j)*Mask( 2)+wpin(i-12,j)*Mask( 3)+wpin(i-11,j)*Mask( 4)+ &
                      &         wpin(i-10,j)*Mask( 5)+wpin(i- 9,j)*Mask( 6)+wpin(i- 8,j)*Mask( 7)+wpin(i- 7,j)*Mask( 8)+ &
                      &         wpin(i- 6,j)*Mask( 9)+wpin(i- 5,j)*Mask(10)+wpin(i- 4,j)*Mask(11)+wpin(i- 3,j)*Mask(12)+ &
                      &         wpin(i- 2,j)*Mask(13)+wpin(i- 1,j)*Mask(14)+wpin(i   ,j)*Mask(15)+wpin(i+ 1,j)*Mask(16)+ &
                      &         wpin(i+ 2,j)*Mask(17)+wpin(i+ 3,j)*Mask(18)+wpin(i+ 4,j)*Mask(19)+wpin(i+ 5,j)*Mask(20)+ &
                      &         wpin(i+ 6,j)*Mask(21)+wpin(i+ 7,j)*Mask(22)+wpin(i+ 8,j)*Mask(23)+wpin(i+ 9,j)*Mask(24)+ &
                      &         wpin(i+10,j)*Mask(25)+wpin(i+11,j)*Mask(26)+wpin(i+12,j)*Mask(27)+wpin(i+13,j)*Mask(28)+ &
                      &         wpin(i+14,j)*Mask(29)
                   ENDDO !i=1,Nm(1)
                   FORALL (i=1:nnodes(1))
                      wpout(i,j)=tmp1_r(i)
                   END FORALL
                ENDDO !j=ldl(2),ldu(2)
             CASE (15)
                DO j=1,nnodes(2)
                   DO i=1,nnodes(1)
                      tmp1_r(i)=wpin(i-15,j)*Mask( 1)+wpin(i-14,j)*Mask( 2)+wpin(i-13,j)*Mask( 3)+wpin(i-12,j)*Mask( 4)+ &
                      &         wpin(i-11,j)*Mask( 5)+wpin(i-10,j)*Mask( 6)+wpin(i- 9,j)*Mask( 7)+wpin(i- 8,j)*Mask( 8)+ &
                      &         wpin(i- 7,j)*Mask( 9)+wpin(i- 6,j)*Mask(10)+wpin(i- 5,j)*Mask(11)+wpin(i- 4,j)*Mask(12)+ &
                      &         wpin(i- 3,j)*Mask(13)+wpin(i- 2,j)*Mask(14)+wpin(i- 1,j)*Mask(15)+wpin(i   ,j)*Mask(16)+ &
                      &         wpin(i+ 1,j)*Mask(17)+wpin(i+ 2,j)*Mask(18)+wpin(i+ 3,j)*Mask(19)+wpin(i+ 4,j)*Mask(20)+ &
                      &         wpin(i+ 5,j)*Mask(21)+wpin(i+ 6,j)*Mask(22)+wpin(i+ 7,j)*Mask(23)+wpin(i+ 8,j)*Mask(24)+ &
                      &         wpin(i+ 9,j)*Mask(25)+wpin(i+10,j)*Mask(26)+wpin(i+11,j)*Mask(27)+wpin(i+12,j)*Mask(28)+ &
                      &         wpin(i+13,j)*Mask(29)+wpin(i+14,j)*Mask(30)+wpin(i+15,j)*Mask(31)
                   ENDDO !i=1,Nm(1)
                   FORALL (i=1:nnodes(1))
                      wpout(i,j)=tmp1_r(i)
                   END FORALL
                ENDDO !j=ldl(2),ldu(2)
             CASE (16)
                DO j=1,nnodes(2)
                   DO i=1,nnodes(1)
                      tmp1_r(i)=wpin(i-16,j)*Mask( 1)+wpin(i-15,j)*Mask( 2)+wpin(i-14,j)*Mask( 3)+wpin(i-13,j)*Mask( 4)+ &
                      &         wpin(i-12,j)*Mask( 5)+wpin(i-11,j)*Mask( 6)+wpin(i-10,j)*Mask( 7)+wpin(i- 9,j)*Mask( 8)+ &
                      &         wpin(i- 8,j)*Mask( 9)+wpin(i- 7,j)*Mask(10)+wpin(i- 6,j)*Mask(11)+wpin(i- 5,j)*Mask(12)+ &
                      &         wpin(i- 4,j)*Mask(13)+wpin(i- 3,j)*Mask(14)+wpin(i- 2,j)*Mask(15)+wpin(i- 1,j)*Mask(16)+ &
                      &         wpin(i   ,j)*Mask(17)+wpin(i+ 1,j)*Mask(18)+wpin(i+ 2,j)*Mask(19)+wpin(i+ 3,j)*Mask(20)+ &
                      &         wpin(i+ 4,j)*Mask(21)+wpin(i+ 5,j)*Mask(22)+wpin(i+ 6,j)*Mask(23)+wpin(i+ 7,j)*Mask(24)+ &
                      &         wpin(i+ 8,j)*Mask(25)+wpin(i+ 9,j)*Mask(26)+wpin(i+10,j)*Mask(27)+wpin(i+11,j)*Mask(28)+ &
                      &         wpin(i+12,j)*Mask(29)+wpin(i+13,j)*Mask(30)+wpin(i+14,j)*Mask(31)+wpin(i+15,j)*Mask(32)+ &
                      &         wpin(i+16,j)*Mask(33)
                   ENDDO !i=1,Nm(1)
                   FORALL (i=1:nnodes(1))
                      wpout(i,j)=tmp1_r(i)
                   END FORALL
                ENDDO !j=ldl(2),ldu(2)
             CASE (17)
                DO j=1,nnodes(2)
                   DO i=1,nnodes(1)
                      tmp1_r(i)=wpin(i-17,j)*Mask( 1)+wpin(i-16,j)*Mask( 2)+wpin(i-15,j)*Mask( 3)+wpin(i-14,j)*Mask( 4)+ &
                      &         wpin(i-13,j)*Mask( 5)+wpin(i-12,j)*Mask( 6)+wpin(i-11,j)*Mask( 7)+wpin(i-10,j)*Mask( 8)+ &
                      &         wpin(i- 9,j)*Mask( 9)+wpin(i- 8,j)*Mask(10)+wpin(i- 7,j)*Mask(11)+wpin(i- 6,j)*Mask(12)+ &
                      &         wpin(i- 5,j)*Mask(13)+wpin(i- 4,j)*Mask(14)+wpin(i- 3,j)*Mask(15)+wpin(i- 2,j)*Mask(16)+ &
                      &         wpin(i- 1,j)*Mask(17)+wpin(i   ,j)*Mask(18)+wpin(i+ 1,j)*Mask(19)+wpin(i+ 2,j)*Mask(20)+ &
                      &         wpin(i+ 3,j)*Mask(21)+wpin(i+ 4,j)*Mask(22)+wpin(i+ 5,j)*Mask(23)+wpin(i+ 6,j)*Mask(24)+ &
                      &         wpin(i+ 7,j)*Mask(25)+wpin(i+ 8,j)*Mask(26)+wpin(i+ 9,j)*Mask(27)+wpin(i+10,j)*Mask(28)+ &
                      &         wpin(i+11,j)*Mask(29)+wpin(i+12,j)*Mask(30)+wpin(i+13,j)*Mask(31)+wpin(i+14,j)*Mask(32)+ &
                      &         wpin(i+15,j)*Mask(33)+wpin(i+16,j)*Mask(34)+wpin(i+17,j)*Mask(35)
                   ENDDO !i=1,Nm(1)
                   FORALL (i=1:nnodes(1))
                      wpout(i,j)=tmp1_r(i)
                   END FORALL
                ENDDO !j=ldl(2),ldu(2)
             CASE (18)
                DO j=1,nnodes(2)
                   DO i=1,nnodes(1)
                      tmp1_r(i)=wpin(i-18,j)*Mask( 1)+wpin(i-17,j)*Mask( 2)+wpin(i-16,j)*Mask( 3)+wpin(i-15,j)*Mask( 4)+ &
                      &         wpin(i-14,j)*Mask( 5)+wpin(i-13,j)*Mask( 6)+wpin(i-12,j)*Mask( 7)+wpin(i-11,j)*Mask( 8)+ &
                      &         wpin(i-10,j)*Mask( 9)+wpin(i- 9,j)*Mask(10)+wpin(i- 8,j)*Mask(11)+wpin(i- 7,j)*Mask(12)+ &
                      &         wpin(i- 6,j)*Mask(13)+wpin(i- 5,j)*Mask(14)+wpin(i- 4,j)*Mask(15)+wpin(i- 3,j)*Mask(16)+ &
                      &         wpin(i- 2,j)*Mask(17)+wpin(i- 1,j)*Mask(18)+wpin(i   ,j)*Mask(19)+wpin(i+ 1,j)*Mask(20)+ &
                      &         wpin(i+ 2,j)*Mask(21)+wpin(i+ 3,j)*Mask(22)+wpin(i+ 4,j)*Mask(23)+wpin(i+ 5,j)*Mask(24)+ &
                      &         wpin(i+ 6,j)*Mask(25)+wpin(i+ 7,j)*Mask(26)+wpin(i+ 8,j)*Mask(27)+wpin(i+ 9,j)*Mask(28)+ &
                      &         wpin(i+10,j)*Mask(29)+wpin(i+11,j)*Mask(30)+wpin(i+12,j)*Mask(31)+wpin(i+13,j)*Mask(32)+ &
                      &         wpin(i+14,j)*Mask(33)+wpin(i+15,j)*Mask(34)+wpin(i+16,j)*Mask(35)+wpin(i+17,j)*Mask(36)+ &
                      &         wpin(i+18,j)*Mask(37)
                   ENDDO !i=1,Nm(1)
                   FORALL (i=1:nnodes(1))
                      wpout(i,j)=tmp1_r(i)
                   END FORALL
                ENDDO !j=ldl(2),ldu(2)
             CASE (19)
                DO j=1,nnodes(2)
                   DO i=1,nnodes(1)
                      tmp1_r(i)=wpin(i-19,j)*Mask( 1)+wpin(i-18,j)*Mask( 2)+wpin(i-17,j)*Mask( 3)+wpin(i-16,j)*Mask( 4)+ &
                      &         wpin(i-15,j)*Mask( 5)+wpin(i-14,j)*Mask( 6)+wpin(i-13,j)*Mask( 7)+wpin(i-12,j)*Mask( 8)+ &
                      &         wpin(i-11,j)*Mask( 9)+wpin(i-10,j)*Mask(10)+wpin(i- 9,j)*Mask(11)+wpin(i- 8,j)*Mask(12)+ &
                      &         wpin(i- 7,j)*Mask(13)+wpin(i- 6,j)*Mask(14)+wpin(i- 5,j)*Mask(15)+wpin(i- 4,j)*Mask(16)+ &
                      &         wpin(i- 3,j)*Mask(17)+wpin(i- 2,j)*Mask(18)+wpin(i- 1,j)*Mask(19)+wpin(i   ,j)*Mask(20)+ &
                      &         wpin(i+ 1,j)*Mask(21)+wpin(i+ 2,j)*Mask(22)+wpin(i+ 3,j)*Mask(23)+wpin(i+ 4,j)*Mask(24)+ &
                      &         wpin(i+ 5,j)*Mask(25)+wpin(i+ 6,j)*Mask(26)+wpin(i+ 7,j)*Mask(27)+wpin(i+ 8,j)*Mask(28)+ &
                      &         wpin(i+ 9,j)*Mask(29)+wpin(i+10,j)*Mask(30)+wpin(i+11,j)*Mask(31)+wpin(i+12,j)*Mask(32)+ &
                      &         wpin(i+13,j)*Mask(33)+wpin(i+14,j)*Mask(34)+wpin(i+15,j)*Mask(35)+wpin(i+16,j)*Mask(36)+ &
                      &         wpin(i+17,j)*Mask(37)+wpin(i+18,j)*Mask(38)+wpin(i+19,j)*Mask(39)
                   ENDDO !i=1,Nm(1)
                   FORALL (i=1:nnodes(1))
                      wpout(i,j)=tmp1_r(i)
                   END FORALL
                ENDDO !j=ldl(2),ldu(2)
             CASE (20)
                DO j=1,nnodes(2)
                   DO i=1,nnodes(1)
                      tmp1_r(i)=wpin(i-20,j)*Mask(1)+wpin(i-19,j)*Mask(2)+wpin(i-18,j)*Mask(3)+wpin(i-17,j)*Mask(4)+     &
                      &         wpin(i-16,j)*Mask(5)+wpin(i-15,j)*Mask(6)+wpin(i-14,j)*Mask(7)+wpin(i-13,j)*Mask(8)+     &
                      &         wpin(i-12,j)*Mask(9)+wpin(i-11,j)*Mask(10)+wpin(i-10,j)*Mask(11)+wpin(i-9,j)*Mask(12)+   &
                      &         wpin(i-8,j)*Mask(13)+wpin(i-7,j)*Mask(14)+wpin(i-6,j)*Mask(15)+wpin(i-5,j)*Mask(16)+     &
                      &         wpin(i-4,j)*Mask(17)+wpin(i-3,j)*Mask(18)+wpin(i-2,j)*Mask(19)+wpin(i-1,j)*Mask(20)+     &
                      &         wpin(i,j)*Mask(21)+wpin(i+1,j)*Mask(22)+wpin(i+2,j)*Mask(23)+wpin(i+3,j)*Mask(24)+       &
                      &         wpin(i+4,j)*Mask(25)+wpin(i+5,j)*Mask(26)+wpin(i+6,j)*Mask(27)+wpin(i+7,j)*Mask(28)+     &
                      &         wpin(i+8,j)*Mask(29)+wpin(i+9,j)*Mask(30)+wpin(i+10,j)*Mask(31)+wpin(i+11,j)*Mask(32)+   &
                      &         wpin(i+12,j)*Mask(33)+wpin(i+13,j)*Mask(34)+wpin(i+14,j)*Mask(35)+wpin(i+15,j)*Mask(36)+ &
                      &         wpin(i+16,j)*Mask(37)+wpin(i+17,j)*Mask(38)+wpin(i+18,j)*Mask(39)+wpin(i+19,j)*Mask(40)+ &
                      &         wpin(i+20,j)*Mask(41)
                   ENDDO !i=1,Nm(1)
                   FORALL (i=1:nnodes(1))
                      wpout(i,j)=tmp1_r(i)
                   END FORALL
                ENDDO !j=ldl(2),ldu(2)
             CASE (21)
                DO j=1,nnodes(2)
                   DO i=1,nnodes(1)
                      tmp1_r(i)=wpin(i-21,j)*Mask(1)+wpin(i-20,j)*Mask(2)+wpin(i-19,j)*Mask(3)+wpin(i-18,j)*Mask(4)+     &
                      &         wpin(i-17,j)*Mask(5)+wpin(i-16,j)*Mask(6)+wpin(i-15,j)*Mask(7)+wpin(i-14,j)*Mask(8)+     &
                      &         wpin(i-13,j)*Mask(9)+wpin(i-12,j)*Mask(10)+wpin(i-11,j)*Mask(11)+wpin(i-10,j)*Mask(12)+  &
                      &         wpin(i-9,j)*Mask(13)+wpin(i-8,j)*Mask(14)+wpin(i-7,j)*Mask(15)+wpin(i-6,j)*Mask(16)+     &
                      &         wpin(i-5,j)*Mask(17)+wpin(i-4,j)*Mask(18)+wpin(i-3,j)*Mask(19)+wpin(i-2,j)*Mask(20)+     &
                      &         wpin(i-1,j)*Mask(21)+wpin(i,j)*Mask(22)+wpin(i+1,j)*Mask(23)+wpin(i+2,j)*Mask(24)+       &
                      &         wpin(i+3,j)*Mask(25)+wpin(i+4,j)*Mask(26)+wpin(i+5,j)*Mask(27)+wpin(i+6,j)*Mask(28)+     &
                      &         wpin(i+7,j)*Mask(29)+wpin(i+8,j)*Mask(30)+wpin(i+9,j)*Mask(31)+wpin(i+10,j)*Mask(32)+    &
                      &         wpin(i+11,j)*Mask(33)+wpin(i+12,j)*Mask(34)+wpin(i+13,j)*Mask(35)+wpin(i+14,j)*Mask(36)+ &
                      &         wpin(i+15,j)*Mask(37)+wpin(i+16,j)*Mask(38)+wpin(i+17,j)*Mask(39)+wpin(i+18,j)*Mask(40)+ &
                      &         wpin(i+19,j)*Mask(41)+wpin(i+20,j)*Mask(42)+wpin(i+21,j)*Mask(43)
                   ENDDO !i=1,Nm(1)
                   FORALL (i=1:nnodes(1))
                      wpout(i,j)=tmp1_r(i)
                   END FORALL
                ENDDO !j=ldl(2),ldu(2)
             CASE (22)
                DO j=1,nnodes(2)
                   DO i=1,nnodes(1)
                      tmp1_r(i)=wpin(i-22,j)*Mask(1)+wpin(i-21,j)*Mask(2)+wpin(i-20,j)*Mask(3)+wpin(i-19,j)*Mask(4)+     &
                      &         wpin(i-18,j)*Mask(5)+wpin(i-17,j)*Mask(6)+wpin(i-16,j)*Mask(7)+wpin(i-15,j)*Mask(8)+     &
                      &         wpin(i-14,j)*Mask(9)+wpin(i-13,j)*Mask(10)+wpin(i-12,j)*Mask(11)+wpin(i-11,j)*Mask(12)+  &
                      &         wpin(i-10,j)*Mask(13)+wpin(i-9,j)*Mask(14)+wpin(i-8,j)*Mask(15)+wpin(i-7,j)*Mask(16)+    &
                      &         wpin(i-6,j)*Mask(17)+wpin(i-5,j)*Mask(18)+wpin(i-4,j)*Mask(19)+wpin(i-3,j)*Mask(20)+     &
                      &         wpin(i-2,j)*Mask(21)+wpin(i-1,j)*Mask(22)+wpin(i,j)*Mask(23)+wpin(i+1,j)*Mask(24)+       &
                      &         wpin(i+2,j)*Mask(25)+wpin(i+3,j)*Mask(26)+wpin(i+4,j)*Mask(27)+wpin(i+5,j)*Mask(28)+     &
                      &         wpin(i+6,j)*Mask(29)+wpin(i+7,j)*Mask(30)+wpin(i+8,j)*Mask(31)+wpin(i+9,j)*Mask(32)+     &
                      &         wpin(i+10,j)*Mask(33)+wpin(i+11,j)*Mask(34)+wpin(i+12,j)*Mask(35)+wpin(i+13,j)*Mask(36)+ &
                      &         wpin(i+14,j)*Mask(37)+wpin(i+15,j)*Mask(38)+wpin(i+16,j)*Mask(39)+wpin(i+17,j)*Mask(40)+ &
                      &         wpin(i+18,j)*Mask(41)+wpin(i+19,j)*Mask(42)+wpin(i+20,j)*Mask(43)+wpin(i+21,j)*Mask(44)+ &
                      &         wpin(i+22,j)*Mask(45)
                   ENDDO !i=1,Nm(1)
                   FORALL (i=1:nnodes(1))
                      wpout(i,j)=tmp1_r(i)
                   END FORALL
                ENDDO !j=ldl(2),ldu(2)
             CASE (23)
                DO j=1,nnodes(2)
                   DO i=1,nnodes(1)
                      tmp1_r(i)=wpin(i-23,j)*Mask(1)+wpin(i-22,j)*Mask(2)+wpin(i-21,j)*Mask(3)+wpin(i-20,j)*Mask(4)+     &
                      &         wpin(i-19,j)*Mask(5)+wpin(i-18,j)*Mask(6)+wpin(i-17,j)*Mask(7)+wpin(i-16,j)*Mask(8)+     &
                      &         wpin(i-15,j)*Mask(9)+wpin(i-14,j)*Mask(10)+wpin(i-13,j)*Mask(11)+wpin(i-12,j)*Mask(12)+  &
                      &         wpin(i-11,j)*Mask(13)+wpin(i-10,j)*Mask(14)+wpin(i-9,j)*Mask(15)+wpin(i-8,j)*Mask(16)+   &
                      &         wpin(i-7,j)*Mask(17)+wpin(i-6,j)*Mask(18)+wpin(i-5,j)*Mask(19)+wpin(i-4,j)*Mask(20)+     &
                      &         wpin(i-3,j)*Mask(21)+wpin(i-2,j)*Mask(22)+wpin(i-1,j)*Mask(23)+wpin(i,j)*Mask(24)+       &
                      &         wpin(i+1,j)*Mask(25)+wpin(i+2,j)*Mask(26)+wpin(i+3,j)*Mask(27)+wpin(i+4,j)*Mask(28)+     &
                      &         wpin(i+5,j)*Mask(29)+wpin(i+6,j)*Mask(30)+wpin(i+7,j)*Mask(31)+wpin(i+8,j)*Mask(32)+     &
                      &         wpin(i+9,j)*Mask(33)+wpin(i+10,j)*Mask(34)+wpin(i+11,j)*Mask(35)+wpin(i+12,j)*Mask(36)+  &
                      &         wpin(i+13,j)*Mask(37)+wpin(i+14,j)*Mask(38)+wpin(i+15,j)*Mask(39)+wpin(i+16,j)*Mask(40)+ &
                      &         wpin(i+17,j)*Mask(41)+wpin(i+18,j)*Mask(42)+wpin(i+19,j)*Mask(43)+wpin(i+20,j)*Mask(44)+ &
                      &         wpin(i+21,j)*Mask(45)+wpin(i+22,j)*Mask(46)+wpin(i+23,j)*Mask(47)
                   ENDDO !i=1,Nm(1)
                   FORALL (i=1:nnodes(1))
                      wpout(i,j)=tmp1_r(i)
                   END FORALL
                ENDDO !j=ldl(2),ldu(2)
             CASE (24)
                DO j=1,nnodes(2)
                   DO i=1,nnodes(1)
                      tmp1_r(i)=wpin(i-24,j)*Mask( 1)+wpin(i-23,j)*Mask( 2)+wpin(i-22,j)*Mask( 3)+wpin(i-21,j)*Mask( 4)+ &
                      &         wpin(i-20,j)*Mask( 5)+wpin(i-19,j)*Mask( 6)+wpin(i-18,j)*Mask( 7)+wpin(i-17,j)*Mask( 8)+ &
                      &         wpin(i-16,j)*Mask( 9)+wpin(i-15,j)*Mask(10)+wpin(i-14,j)*Mask(11)+wpin(i-13,j)*Mask(12)+ &
                      &         wpin(i-12,j)*Mask(13)+wpin(i-11,j)*Mask(14)+wpin(i-10,j)*Mask(15)+wpin(i- 9,j)*Mask(16)+ &
                      &         wpin(i- 8,j)*Mask(17)+wpin(i- 7,j)*Mask(18)+wpin(i- 6,j)*Mask(19)+wpin(i- 5,j)*Mask(20)+ &
                      &         wpin(i- 4,j)*Mask(21)+wpin(i- 3,j)*Mask(22)+wpin(i- 2,j)*Mask(23)+wpin(i- 1,j)*Mask(24)+ &
                      &         wpin(i   ,j)*Mask(25)+wpin(i+ 1,j)*Mask(26)+wpin(i+ 2,j)*Mask(27)+wpin(i+ 3,j)*Mask(28)+ &
                      &         wpin(i+ 4,j)*Mask(29)+wpin(i+ 5,j)*Mask(30)+wpin(i+ 6,j)*Mask(31)+wpin(i+ 7,j)*Mask(32)+ &
                      &         wpin(i+ 8,j)*Mask(33)+wpin(i+ 9,j)*Mask(34)+wpin(i+10,j)*Mask(35)+wpin(i+11,j)*Mask(36)+ &
                      &         wpin(i+12,j)*Mask(37)+wpin(i+13,j)*Mask(38)+wpin(i+14,j)*Mask(39)+wpin(i+15,j)*Mask(40)+ &
                      &         wpin(i+16,j)*Mask(41)+wpin(i+17,j)*Mask(42)+wpin(i+18,j)*Mask(43)+wpin(i+19,j)*Mask(44)+ &
                      &         wpin(i+20,j)*Mask(45)+wpin(i+21,j)*Mask(46)+wpin(i+22,j)*Mask(47)+wpin(i+23,j)*Mask(48)+ &
                      &         wpin(i+24,j)*Mask(49)
                   ENDDO !i=1,Nm(1)
                   FORALL (i=1:nnodes(1))
                      wpout(i,j)=tmp1_r(i)
                   END FORALL
                ENDDO !j=ldl(2),ldu(2)
             CASE DEFAULT
                m=MaskSize+1
                DO j=1,nnodes(2)
                   tmp1_r=__ZERO
                   DO i=1,nnodes(1)
                      DO l=-MaskSize,MaskSize
                         tmp1_r(i)=tmp1_r(i)+wpin(i+l,j)*Mask(l+m)
                      ENDDO
                   ENDDO !i=1,Nm(1)
                   FORALL (i=1:nnodes(1))
                      wpout(i,j)=tmp1_r(i)
                   END FORALL
                ENDDO !j=ldl(2),ldu(2)
             END SELECT

             DEALLOCATE(tmp1_r,STAT=info)
             or_fail_dealloc("Failed to deallocate tmp1_r!",exit_point=7777)
          CASE (2)
          !!! Y direction
             SELECT CASE (MaskSize)
             CASE (1)
                !At each pixel, I should keep three values of the first, second
                !and third row
                ALLOCATE(tmp2_r(nnodes(1),3),ll(3),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                ! This optimisation is done to avoid recomputing what has been computed once

                DO j=0,1
                   l=j+1
                   FORALL (i=1:nnodes(1))
                      tmp2_r(i,l)=wpin(i,j)
                   END FORALL
                ENDDO
                FORALL(i=1:3) ll(i)=i
                DO j=1,nnodes(2)
                   jp=j+1
                   l3=ll(3)
                   FORALL (i=1:nnodes(1))
                      tmp2_r(i,l3)=wpin(i,jp)
                   END FORALL
                   l1=ll(1)
                   l2=ll(2)
                   FORALL (i=1:nnodes(1))
                      wpout(i,j)=tmp2_r(i,l1)*Mask(1)+tmp2_r(i,l2)*Mask(2)+tmp2_r(i,l3)*Mask(3)
                   END FORALL
                   ll=CSHIFT(ll,SHIFT=1,DIM=1)
                ENDDO !j=1,nnodes(2)
             CASE (2)
                ALLOCATE(tmp2_r(nnodes(1),5),ll(5),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                DO j=-1,2
                   l=j+2
                   FORALL (i=1:nnodes(1))
                      tmp2_r(i,l)=wpin(i,j)
                   END FORALL
                ENDDO
                FORALL(i=1:5) ll(i)=i
                DO j=1,nnodes(2)
                   jp=j+2
                   l5=ll(5)
                   FORALL (i=1:nnodes(1))
                      tmp2_r(i,l5)=wpin(i,jp)
                   END FORALL
                   l1=ll(1)
                   l2=ll(2)
                   l3=ll(3)
                   l4=ll(4)
                   FORALL (i=1:nnodes(1))
                      wpout(i,j)=tmp2_r(i,l1)*Mask(1)+tmp2_r(i,l2)*Mask(2)+tmp2_r(i,l3)*Mask(3)+tmp2_r(i,l4)*Mask(4)+ &
                      &          tmp2_r(i,l5)*Mask(5)
                   END FORALL
                   ll=CSHIFT(ll,SHIFT=1,DIM=1)
                ENDDO !j=1,nnodes(2)
             CASE (3)
                ALLOCATE(tmp2_r(nnodes(1),7),ll(7),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                DO j=-2,3
                   l=j+3
                   FORALL (i=1:nnodes(1))
                      tmp2_r(i,l)=wpin(i,j)
                   END FORALL
                ENDDO
                FORALL(i=1:7) ll(i)=i
                DO j=1,nnodes(2)
                   jp=j+3
                   l7=ll(7)
                   FORALL (i=1:nnodes(1))
                      tmp2_r(i,l7)=wpin(i,jp)
                   END FORALL
                   l1=ll(1)
                   l2=ll(2)
                   l3=ll(3)
                   l4=ll(4)
                   l5=ll(5)
                   l6=ll(6)
                   FORALL (i=1:nnodes(1))
                      wpout(i,j)=tmp2_r(i,l1)*Mask(1)+tmp2_r(i,l2)*Mask(2)+tmp2_r(i,l3)*Mask(3)+tmp2_r(i,l4)*Mask(4)+ &
                      &          tmp2_r(i,l5)*Mask(5)+tmp2_r(i,l6)*Mask(6)+tmp2_r(i,l7)*Mask(7)
                   END FORALL
                   ll=CSHIFT(ll,SHIFT=1,DIM=1)
                ENDDO !j=1,nnodes(2)
             CASE (4)
                ALLOCATE(tmp2_r(nnodes(1),9),ll(9),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                DO j=-3,4
                   l=j+4
                   FORALL (i=1:nnodes(1))
                      tmp2_r(i,l)=wpin(i,j)
                   END FORALL
                ENDDO
                FORALL(i=1:9) ll(i)=i
                DO j=1,nnodes(2)
                   jp=j+4
                   l9=ll(9)
                   FORALL (i=1:nnodes(1))
                      tmp2_r(i,l9)=wpin(i,jp)
                   END FORALL
                   l1=ll(1)
                   l2=ll(2)
                   l3=ll(3)
                   l4=ll(4)
                   l5=ll(5)
                   l6=ll(6)
                   l7=ll(7)
                   l8=ll(8)
                   FORALL (i=1:nnodes(1))
                      wpout(i,j)=tmp2_r(i,l1)*Mask(1)+tmp2_r(i,l2)*Mask(2)+tmp2_r(i,l3)*Mask(3)+tmp2_r(i,l4)*Mask(4)+ &
                      &          tmp2_r(i,l5)*Mask(5)+tmp2_r(i,l6)*Mask(6)+tmp2_r(i,l7)*Mask(7)+tmp2_r(i,l8)*Mask(8)+ &
                      &          tmp2_r(i,l9)*Mask(9)
                   END FORALL
                   ll=CSHIFT(ll,SHIFT=1,DIM=1)
                ENDDO !j=1,nnodes(2)
             CASE (5)
                ALLOCATE(tmp2_r(nnodes(1),11),ll(11),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                DO j=-4,5
                   l=j+5
                   FORALL (i=1:nnodes(1))
                      tmp2_r(i,l)=wpin(i,j)
                   END FORALL
                ENDDO
                FORALL(i=1:11) ll(i)=i
                DO j=1,nnodes(2)
                   jp=j+5
                   l11=ll(11)
                   FORALL (i=1:nnodes(1))
                      tmp2_r(i,l11)=wpin(i,jp)
                   END FORALL
                   l1 =ll( 1)
                   l2 =ll( 2)
                   l3 =ll( 3)
                   l4 =ll( 4)
                   l5 =ll( 5)
                   l6 =ll( 6)
                   l7 =ll( 7)
                   l8 =ll( 8)
                   l9 =ll( 9)
                   l10=ll(10)
                   FORALL (i=1:nnodes(1))
                      wpout(i,j)=tmp2_r(i,l1)*Mask(1)+tmp2_r(i,l2 )*Mask( 2)+tmp2_r(i,l3 )*Mask( 3)+tmp2_r(i,l4)*Mask(4)+ &
                      &          tmp2_r(i,l5)*Mask(5)+tmp2_r(i,l6 )*Mask( 6)+tmp2_r(i,l7 )*Mask( 7)+tmp2_r(i,l8)*Mask(8)+ &
                      &          tmp2_r(i,l9)*Mask(9)+tmp2_r(i,l10)*Mask(10)+tmp2_r(i,l11)*Mask(11)
                   END FORALL
                   ll=CSHIFT(ll,SHIFT=1,DIM=1)
                ENDDO !j=1,nnodes(2)
             CASE (6)
                ALLOCATE(tmp2_r(nnodes(1),13),ll(13),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                DO j=-5,6
                   l=j+6
                   FORALL (i=1:nnodes(1))
                      tmp2_r(i,l)=wpin(i,j)
                   END FORALL
                ENDDO
                FORALL(i=1:13) ll(i)=i
                DO j=1,nnodes(2)
                   jp=j+6
                   l13=ll(13)
                   FORALL (i=1:nnodes(1))
                      tmp2_r(i,l13)=wpin(i,jp)
                   END FORALL
                   l1 =ll( 1)
                   l2 =ll( 2)
                   l3 =ll( 3)
                   l4 =ll( 4)
                   l5 =ll( 5)
                   l6 =ll( 6)
                   l7 =ll( 7)
                   l8 =ll( 8)
                   l9 =ll( 9)
                   l10=ll(10)
                   l11=ll(11)
                   l12=ll(12)
                   FORALL (i=1:nnodes(1))
                      wpout(i,j)=tmp2_r(i,l1 )*Mask( 1)+tmp2_r(i,l2 )*Mask( 2)+tmp2_r(i,l3 )*Mask( 3)+tmp2_r(i,l4 )*Mask( 4)+ &
                      &          tmp2_r(i,l5 )*Mask( 5)+tmp2_r(i,l6 )*Mask( 6)+tmp2_r(i,l7 )*Mask( 7)+tmp2_r(i,l8 )*Mask( 8)+ &
                      &          tmp2_r(i,l9 )*Mask( 9)+tmp2_r(i,l10)*Mask(10)+tmp2_r(i,l11)*Mask(11)+tmp2_r(i,l12)*Mask(12)+ &
                      &          tmp2_r(i,l13)*Mask(13)
                   END FORALL
                   ll=CSHIFT(ll,SHIFT=1,DIM=1)
                ENDDO !j=1,nnodes(2)
             CASE (7)
                ALLOCATE(tmp2_r(nnodes(1),15),ll(15),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                DO j=-6,7
                   l=j+7
                   FORALL (i=1:nnodes(1))
                      tmp2_r(i,l)=wpin(i,j)
                   END FORALL
                ENDDO
                FORALL(i=1:15) ll(i)=i
                DO j=1,nnodes(2)
                   jp=j+7
                   l15=ll(15)
                   FORALL (i=1:nnodes(1))
                      tmp2_r(i,l15)=wpin(i,jp)
                   END FORALL
                   l1 =ll( 1)
                   l2 =ll( 2)
                   l3 =ll( 3)
                   l4 =ll( 4)
                   l5 =ll( 5)
                   l6 =ll( 6)
                   l7 =ll( 7)
                   l8 =ll( 8)
                   l9 =ll( 9)
                   l10=ll(10)
                   l11=ll(11)
                   l12=ll(12)
                   l13=ll(13)
                   l14=ll(14)
                   FORALL (i=1:nnodes(1))
                      wpout(i,j)=tmp2_r(i,l1 )*Mask( 1)+tmp2_r(i,l2 )*Mask( 2)+tmp2_r(i,l3 )*Mask( 3)+tmp2_r(i,l4 )*Mask( 4)+ &
                      &          tmp2_r(i,l5 )*Mask( 5)+tmp2_r(i,l6 )*Mask( 6)+tmp2_r(i,l7 )*Mask( 7)+tmp2_r(i,l8 )*Mask( 8)+ &
                      &          tmp2_r(i,l9 )*Mask( 9)+tmp2_r(i,l10)*Mask(10)+tmp2_r(i,l11)*Mask(11)+tmp2_r(i,l12)*Mask(12)+ &
                      &          tmp2_r(i,l13)*Mask(13)+tmp2_r(i,l14)*Mask(14)+tmp2_r(i,l15)*Mask(15)
                   END FORALL
                   ll=CSHIFT(ll,SHIFT=1,DIM=1)
                ENDDO !j=1,nnodes(2)
             CASE (8)
                ALLOCATE(tmp2_r(nnodes(1),17),ll(17),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                DO j=-7,8
                   l=j+8
                   FORALL (i=1:nnodes(1))
                      tmp2_r(i,l)=wpin(i,j)
                   END FORALL
                ENDDO
                FORALL(i=1:17) ll(i)=i
                DO j=1,nnodes(2)
                   jp=j+8
                   l17=ll(17)
                   FORALL (i=1:nnodes(1))
                      tmp2_r(i,l17)=wpin(i,jp)
                   END FORALL
                   l1 =ll( 1)
                   l2 =ll( 2)
                   l3 =ll( 3)
                   l4 =ll( 4)
                   l5 =ll( 5)
                   l6 =ll( 6)
                   l7 =ll( 7)
                   l8 =ll( 8)
                   l9 =ll( 9)
                   l10=ll(10)
                   l11=ll(11)
                   l12=ll(12)
                   l13=ll(13)
                   l14=ll(14)
                   l15=ll(15)
                   l16=ll(16)
                   FORALL (i=1:nnodes(1))
                      wpout(i,j)=tmp2_r(i,l1 )*Mask( 1)+tmp2_r(i,l2 )*Mask( 2)+tmp2_r(i,l3 )*Mask( 3)+tmp2_r(i,l4 )*Mask( 4)+ &
                      &          tmp2_r(i,l5 )*Mask( 5)+tmp2_r(i,l6 )*Mask( 6)+tmp2_r(i,l7 )*Mask( 7)+tmp2_r(i,l8 )*Mask( 8)+ &
                      &          tmp2_r(i,l9 )*Mask( 9)+tmp2_r(i,l10)*Mask(10)+tmp2_r(i,l11)*Mask(11)+tmp2_r(i,l12)*Mask(12)+ &
                      &          tmp2_r(i,l13)*Mask(13)+tmp2_r(i,l14)*Mask(14)+tmp2_r(i,l15)*Mask(15)+tmp2_r(i,l16)*Mask(16)+ &
                      &          tmp2_r(i,l17)*Mask(17)
                   END FORALL
                   ll=CSHIFT(ll,SHIFT=1,DIM=1)
                ENDDO !j=1,nnodes(2)
             CASE (9)
                ALLOCATE(tmp2_r(nnodes(1),19),ll(19),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                DO j=-8,9
                   l=j+9
                   FORALL (i=1:nnodes(1))
                      tmp2_r(i,l)=wpin(i,j)
                   END FORALL
                ENDDO
                FORALL(i=1:19) ll(i)=i
                DO j=1,nnodes(2)
                   jp=j+9
                   l19=ll(19)
                   FORALL (i=1:nnodes(1))
                      tmp2_r(i,l19)=wpin(i,jp)
                   END FORALL
                   l1 =ll( 1)
                   l2 =ll( 2)
                   l3 =ll( 3)
                   l4 =ll( 4)
                   l5 =ll( 5)
                   l6 =ll( 6)
                   l7 =ll( 7)
                   l8 =ll( 8)
                   l9 =ll( 9)
                   l10=ll(10)
                   l11=ll(11)
                   l12=ll(12)
                   l13=ll(13)
                   l14=ll(14)
                   l15=ll(15)
                   l16=ll(16)
                   l17=ll(17)
                   l18=ll(18)
                   FORALL (i=1:nnodes(1))
                      wpout(i,j)=tmp2_r(i,l1)*Mask(1)+tmp2_r(i,l2)*Mask(2)+tmp2_r(i,l3)*Mask(3)+tmp2_r(i,l4)*Mask(4)+         &
                      &          tmp2_r(i,l5)*Mask(5)+tmp2_r(i,l6)*Mask(6)+tmp2_r(i,l7)*Mask(7)+tmp2_r(i,l8)*Mask(8)+         &
                      &          tmp2_r(i,l9)*Mask(9)+tmp2_r(i,l10)*Mask(10)+tmp2_r(i,l11)*Mask(11)+tmp2_r(i,l12)*Mask(12)+   &
                      &          tmp2_r(i,l13)*Mask(13)+tmp2_r(i,l14)*Mask(14)+tmp2_r(i,l15)*Mask(15)+tmp2_r(i,l16)*Mask(16)+ &
                      &          tmp2_r(i,l17)*Mask(17)+tmp2_r(i,l18)*Mask(18)+tmp2_r(i,l19)*Mask(19)
                   END FORALL
                   ll=CSHIFT(ll,SHIFT=1,DIM=1)
                ENDDO !j=1,nnodes(2)
             CASE (10)
                ALLOCATE(tmp2_r(nnodes(1),21),ll(21),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                DO j=-9,10
                   l=j+10
                   FORALL (i=1:nnodes(1))
                      tmp2_r(i,l)=wpin(i,j)
                   END FORALL
                ENDDO
                FORALL(i=1:21) ll(i)=i
                DO j=1,nnodes(2)
                   jp=j+10
                   l21=ll(21)
                   FORALL (i=1:nnodes(1))
                      tmp2_r(i,l21)=wpin(i,jp)
                   END FORALL
                   l1 =ll( 1)
                   l2 =ll( 2)
                   l3 =ll( 3)
                   l4 =ll( 4)
                   l5 =ll( 5)
                   l6 =ll( 6)
                   l7 =ll( 7)
                   l8 =ll( 8)
                   l9 =ll( 9)
                   l10=ll(10)
                   l11=ll(11)
                   l12=ll(12)
                   l13=ll(13)
                   l14=ll(14)
                   l15=ll(15)
                   l16=ll(16)
                   l17=ll(17)
                   l18=ll(18)
                   l19=ll(19)
                   l20=ll(20)
                   FORALL (i=1:nnodes(1))
                      wpout(i,j)=tmp2_r(i,l1 )*Mask( 1)+tmp2_r(i,l2 )*Mask( 2)+tmp2_r(i,l3 )*Mask( 3)+tmp2_r(i,l4 )*Mask( 4)+ &
                      &          tmp2_r(i,l5 )*Mask( 5)+tmp2_r(i,l6 )*Mask( 6)+tmp2_r(i,l7 )*Mask( 7)+tmp2_r(i,l8 )*Mask( 8)+ &
                      &          tmp2_r(i,l9 )*Mask( 9)+tmp2_r(i,l10)*Mask(10)+tmp2_r(i,l11)*Mask(11)+tmp2_r(i,l12)*Mask(12)+ &
                      &          tmp2_r(i,l13)*Mask(13)+tmp2_r(i,l14)*Mask(14)+tmp2_r(i,l15)*Mask(15)+tmp2_r(i,l16)*Mask(16)+ &
                      &          tmp2_r(i,l17)*Mask(17)+tmp2_r(i,l18)*Mask(18)+tmp2_r(i,l19)*Mask(19)+tmp2_r(i,l20)*Mask(20)+ &
                      &          tmp2_r(i,l21)*Mask(21)
                   END FORALL
                   ll=CSHIFT(ll,SHIFT=1,DIM=1)
                ENDDO !j=1,nnodes(2)

             CASE DEFAULT
                m=MaskSize*2+1

                ALLOCATE(tmp2_r(nnodes(1),m),ll(m),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                DO j=1,m-1
                   l=j-MaskSize
                   FORALL (i=1:nnodes(1))
                      tmp2_r(i,j)=wpin(i,l)
                   END FORALL
                ENDDO
                FORALL(i=1:m) ll(i)=i
                DO j=1,nnodes(2)
                   jp=j+MaskSize
                   l=ll(m)
                   FORALL (i=1:nnodes(1))
                      tmp2_r(i,l)=wpin(i,jp)
                   END FORALL
                   DO i=1,nnodes(1)
                      wpout(i,j)=tmp2_r(i,l)*Mask(m)
                      DO n=1,m-1
                         l1=ll(n)
                         wpout(i,j)=wpout(i,j)+tmp2_r(i,l1)*Mask(n)
                      ENDDO
                   ENDDO
                   ll=CSHIFT(ll,SHIFT=1,DIM=1)
                ENDDO !j=1,nnodes(2)
             END SELECT

             DEALLOCATE(tmp2_r,ll,STAT=info)
             or_fail_dealloc("tmp2",exit_point=7777)
          CASE DEFAULT
             fail("Wrong direction for 2D problem!!!",exit_point=7777)
          END SELECT
!           ELSE !IF (ASSOCIATED(wpout,wpin)) THEN
!           ENDIF !IF (ASSOCIATED(wpout,wpin)) THEN
        7777 CONTINUE
        END SUBROUTINE Update
#elif __DIME == __3D
        SUBROUTINE Update(Mask,MaskSize,wpin,wpout,nnodes,direction,info)
          IMPLICIT NONE
          REAL(MK),             DIMENSION(:),     INTENT(IN   ) :: Mask
          INTEGER,                                INTENT(IN   ) :: MaskSize
          REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: wpin
          REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: wpout
          INTEGER,              DIMENSION(:),     INTENT(IN   ) :: nnodes
          INTEGER,                                INTENT(IN   ) :: direction
          INTEGER,                                INTENT(  OUT) :: info


          INTEGER                            :: l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11
          INTEGER                            :: l12,l13,l14,l15,l16,l17,l18,l19,l20,l21
          INTEGER                            :: i,j,k,l,m,n
          INTEGER                            :: jp,kp
          INTEGER, DIMENSION(:), ALLOCATABLE :: ll

          SELECT CASE (direction)
          CASE (1)
          !!! X direction
             ALLOCATE(tmp1_r(nnodes(1)),STAT=info)
             or_fail_alloc("Failed to allocate tmp1_r!",exit_point=7777)

             SELECT CASE (MaskSize)
             CASE (1)
                DO k=1,nnodes(3)
                   DO j=1,nnodes(2)
                      DO i=1,nnodes(1)
                         tmp1_r(i)=wpin(i-1,j,k)*Mask(1)+wpin(i,j,k)*Mask(2)+wpin(i+1,j,k)*Mask(3)
                      ENDDO !i=1,nnodes(1)
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp1_r(i)
                      END FORALL
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (2)
                DO k=1,nnodes(3)
                   DO j=1,nnodes(2)
                      DO i=1,nnodes(1)
                         tmp1_r(i)=wpin(i-2,j,k)*Mask(1)+wpin(i-1,j,k)*Mask(2)+wpin(i,j,k)*Mask(3)+wpin(i+1,j,k)*Mask(4)+ &
                         &         wpin(i+2,j,k)*Mask(5)
                      ENDDO !i=1,nnodes(1)
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp1_r(i)
                      END FORALL
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (3)
                DO k=1,nnodes(3)
                   DO j=1,nnodes(2)
                      DO i=1,nnodes(1)
                         tmp1_r(i)=wpin(i-3,j,k)*Mask(1)+wpin(i-2,j,k)*Mask(2)+wpin(i-1,j,k)*Mask(3)+wpin(i,j,k)*Mask(4)+ &
                         &         wpin(i+1,j,k)*Mask(5)+wpin(i+2,j,k)*Mask(6)+wpin(i+3,j,k)*Mask(7)
                      ENDDO !i=1,nnodes(1)
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp1_r(i)
                      END FORALL
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (4)
                DO k=1,nnodes(3)
                   DO j=1,nnodes(2)
                      DO i=1,nnodes(1)
                         tmp1_r(i)=wpin(i-4,j,k)*Mask(1)+wpin(i-3,j,k)*Mask(2)+wpin(i-2,j,k)*Mask(3)+wpin(i-1,j,k)*Mask(4)+ &
                         &         wpin(i,j,k)*Mask(5)+wpin(i+1,j,k)*Mask(6)+wpin(i+2,j,k)*Mask(7)+wpin(i+3,j,k)*Mask(8)+   &
                         &         wpin(i+4,j,k)*Mask(9)
                      ENDDO !i=1,nnodes(1)
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp1_r(i)
                      END FORALL
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (5)
                DO k=1,nnodes(3)
                   DO j=1,nnodes(2)
                      DO i=1,nnodes(1)
                         tmp1_r(i)=wpin(i-5,j,k)*Mask(1)+wpin(i-4,j,k)*Mask(2)+wpin(i-3,j,k)*Mask(3)+wpin(i-2,j,k)*Mask(4)+ &
                         &         wpin(i-1,j,k)*Mask(5)+wpin(i,j,k)*Mask(6)+wpin(i+1,j,k)*Mask(7)+wpin(i+2,j,k)*Mask(8)+   &
                         &         wpin(i+3,j,k)*Mask(9)+wpin(i+4,j,k)*Mask(10)+wpin(i+5,j,k)*Mask(11)
                      ENDDO !i=1,nnodes(1)
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp1_r(i)
                      END FORALL
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (6)
                DO k=1,nnodes(3)
                   DO j=1,nnodes(2)
                      DO i=1,nnodes(1)
                         tmp1_r(i)=wpin(i-6,j,k)*Mask(1)+wpin(i-5,j,k)*Mask(2)+wpin(i-4,j,k)*Mask(3)+wpin(i-3,j,k)*Mask(4)+    &
                         &         wpin(i-2,j,k)*Mask(5)+wpin(i-1,j,k)*Mask(6)+wpin(i,j,k)*Mask(7)+wpin(i+1,j,k)*Mask(8)+      &
                         &         wpin(i+2,j,k)*Mask(9)+wpin(i+3,j,k)*Mask(10)+wpin(i+4,j,k)*Mask(11)+wpin(i+5,j,k)*Mask(12)+ &
                         &         wpin(i+6,j,k)*Mask(13)
                      ENDDO !i=1,nnodes(1)
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp1_r(i)
                      END FORALL
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (7)
                DO k=1,nnodes(3)
                   DO j=1,nnodes(2)
                      DO i=1,nnodes(1)
                         tmp1_r(i)=wpin(i-7,j,k)*Mask(1)+wpin(i-6,j,k)*Mask(2)+wpin(i-5,j,k)*Mask(3)+wpin(i-4,j,k)*Mask(4)+    &
                         &         wpin(i-3,j,k)*Mask(5)+wpin(i-2,j,k)*Mask(6)+wpin(i-1,j,k)*Mask(7)+wpin(i,j,k)*Mask(8)+      &
                         &         wpin(i+1,j,k)*Mask(9)+wpin(i+2,j,k)*Mask(10)+wpin(i+3,j,k)*Mask(11)+wpin(i+4,j,k)*Mask(12)+ &
                         &         wpin(i+5,j,k)*Mask(13)+wpin(i+6,j,k)*Mask(14)+wpin(i+7,j,k)*Mask(15)
                      ENDDO !i=1,nnodes(1)
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp1_r(i)
                      END FORALL
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (8)
                DO k=1,nnodes(3)
                   DO j=1,nnodes(2)
                      DO i=1,nnodes(1)
                         tmp1_r(i)=wpin(i-8,j,k)*Mask(1)+wpin(i-7,j,k)*Mask(2)+wpin(i-6,j,k)*Mask(3)+wpin(i-5,j,k)*Mask(4)+     &
                         &         wpin(i-4,j,k)*Mask(5)+wpin(i-3,j,k)*Mask(6)+wpin(i-2,j,k)*Mask(7)+wpin(i-1,j,k)*Mask(8)+     &
                         &         wpin(i,j,k)*Mask(9)+wpin(i+1,j,k)*Mask(10)+wpin(i+2,j,k)*Mask(11)+wpin(i+3,j,k)*Mask(12)+    &
                         &         wpin(i+4,j,k)*Mask(13)+wpin(i+5,j,k)*Mask(14)+wpin(i+6,j,k)*Mask(15)+wpin(i+7,j,k)*Mask(16)+ &
                         &         wpin(i+8,j,k)*Mask(17)
                      ENDDO !i=1,nnodes(1)
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp1_r(i)
                      END FORALL
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (9)
                DO k=1,nnodes(3)
                   DO j=1,nnodes(2)
                      DO i=1,nnodes(1)
                         tmp1_r(i)=wpin(i-9,j,k)*Mask(1)+wpin(i-8,j,k)*Mask(2)+wpin(i-7,j,k)*Mask(3)+wpin(i-6,j,k)*Mask(4)+     &
                         &         wpin(i-5,j,k)*Mask(5)+wpin(i-4,j,k)*Mask(6)+wpin(i-3,j,k)*Mask(7)+wpin(i-2,j,k)*Mask(8)+     &
                         &         wpin(i-1,j,k)*Mask(9)+wpin(i,j,k)*Mask(10)+wpin(i+1,j,k)*Mask(11)+wpin(i+2,j,k)*Mask(12)+    &
                         &         wpin(i+3,j,k)*Mask(13)+wpin(i+4,j,k)*Mask(14)+wpin(i+5,j,k)*Mask(15)+wpin(i+6,j,k)*Mask(16)+ &
                         &         wpin(i+7,j,k)*Mask(17)+wpin(i+8,j,k)*Mask(18)+wpin(i+9,j,k)*Mask(19)
                      ENDDO !i=1,nnodes(1)
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp1_r(i)
                      END FORALL
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (10)
                DO k=1,nnodes(3)
                   DO j=1,nnodes(2)
                      DO i=1,nnodes(1)
                         tmp1_r(i)=wpin(i-10,j,k)*Mask(1)+wpin(i-9,j,k)*Mask(2)+wpin(i-8,j,k)*Mask(3)+wpin(i-7,j,k)*Mask(4)+    &
                         &         wpin(i-6,j,k)*Mask(5)+wpin(i-5,j,k)*Mask(6)+wpin(i-4,j,k)*Mask(7)+wpin(i-3,j,k)*Mask(8)+     &
                         &         wpin(i-2,j,k)*Mask(9)+wpin(i-1,j,k)*Mask(10)+wpin(i,j,k)*Mask(11)+wpin(i+1,j,k)*Mask(12)+    &
                         &         wpin(i+2,j,k)*Mask(13)+wpin(i+3,j,k)*Mask(14)+wpin(i+4,j,k)*Mask(15)+wpin(i+5,j,k)*Mask(16)+ &
                         &         wpin(i+6,j,k)*Mask(17)+wpin(i+7,j,k)*Mask(18)+wpin(i+8,j,k)*Mask(19)+wpin(i+9,j,k)*Mask(20)+ &
                         &         wpin(i+10,j,k)*Mask(21)
                      ENDDO !i=1,nnodes(1)
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp1_r(i)
                      END FORALL
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (11)
                DO k=1,nnodes(3)
                   DO j=1,nnodes(2)
                      DO i=1,nnodes(1)
                         tmp1_r(i)=wpin(i-11,j,k)*Mask(1)+wpin(i-10,j,k)*Mask(2)+wpin(i-9,j,k)*Mask(3)+wpin(i-8,j,k)*Mask(4)+   &
                         &         wpin(i-7,j,k)*Mask(5)+wpin(i-6,j,k)*Mask(6)+wpin(i-5,j,k)*Mask(7)+wpin(i-4,j,k)*Mask(8)+     &
                         &         wpin(i-3,j,k)*Mask(9)+wpin(i-2,j,k)*Mask(10)+wpin(i-1,j,k)*Mask(11)+wpin(i,j,k)*Mask(12)+    &
                         &         wpin(i+1,j,k)*Mask(13)+wpin(i+2,j,k)*Mask(14)+wpin(i+3,j,k)*Mask(15)+wpin(i+4,j,k)*Mask(16)+ &
                         &         wpin(i+5,j,k)*Mask(17)+wpin(i+6,j,k)*Mask(18)+wpin(i+7,j,k)*Mask(19)+wpin(i+8,j,k)*Mask(20)+ &
                         &         wpin(i+9,j,k)*Mask(21)+wpin(i+10,j,k)*Mask(22)+wpin(i+11,j,k)*Mask(23)
                      ENDDO !i=1,nnodes(1)
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp1_r(i)
                      END FORALL
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (12)
                DO k=1,nnodes(3)
                   DO j=1,nnodes(2)
                      DO i=1,nnodes(1)
                         tmp1_r(i)=wpin(i-12,j,k)*Mask(1)+wpin(i-11,j,k)*Mask(2)+wpin(i-10,j,k)*Mask(3)+wpin(i-9,j,k)*Mask(4)+    &
                         &         wpin(i-8,j,k)*Mask(5)+wpin(i-7,j,k)*Mask(6)+wpin(i-6,j,k)*Mask(7)+wpin(i-5,j,k)*Mask(8)+       &
                         &         wpin(i-4,j,k)*Mask(9)+wpin(i-3,j,k)*Mask(10)+wpin(i-2,j,k)*Mask(11)+wpin(i-1,j,k)*Mask(12)+    &
                         &         wpin(i,j,k)*Mask(13)+wpin(i+1,j,k)*Mask(14)+wpin(i+2,j,k)*Mask(15)+wpin(i+3,j,k)*Mask(16)+     &
                         &         wpin(i+4,j,k)*Mask(17)+wpin(i+5,j,k)*Mask(18)+wpin(i+6,j,k)*Mask(19)+wpin(i+7,j,k)*Mask(20)+   &
                         &         wpin(i+8,j,k)*Mask(21)+wpin(i+9,j,k)*Mask(22)+wpin(i+10,j,k)*Mask(23)+wpin(i+11,j,k)*Mask(24)+ &
                         &         wpin(i+12,j,k)*Mask(25)
                      ENDDO !i=1,nnodes(1)
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp1_r(i)
                      END FORALL
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (13)
                DO k=1,nnodes(3)
                   DO j=1,nnodes(2)
                      DO i=1,nnodes(1)
                         tmp1_r(i)=wpin(i-13,j,k)*Mask(1)+wpin(i-12,j,k)*Mask(2)+wpin(i-11,j,k)*Mask(3)+wpin(i-10,j,k)*Mask(4)+   &
                         &         wpin(i-9,j,k)*Mask(5)+wpin(i-8,j,k)*Mask(6)+wpin(i-7,j,k)*Mask(7)+wpin(i-6,j,k)*Mask(8)+       &
                         &         wpin(i-5,j,k)*Mask(9)+wpin(i-4,j,k)*Mask(10)+wpin(i-3,j,k)*Mask(11)+wpin(i-2,j,k)*Mask(12)+    &
                         &         wpin(i-1,j,k)*Mask(13)+wpin(i,j,k)*Mask(14)+wpin(i+1,j,k)*Mask(15)+wpin(i+2,j,k)*Mask(16)+     &
                         &         wpin(i+3,j,k)*Mask(17)+wpin(i+4,j,k)*Mask(18)+wpin(i+5,j,k)*Mask(19)+wpin(i+6,j,k)*Mask(20)+   &
                         &         wpin(i+7,j,k)*Mask(21)+wpin(i+8,j,k)*Mask(22)+wpin(i+9,j,k)*Mask(23)+wpin(i+10,j,k)*Mask(24)+  &
                         &         wpin(i+11,j,k)*Mask(25)+wpin(i+12,j,k)*Mask(26)+wpin(i+13,j,k)*Mask(27)
                      ENDDO !i=1,nnodes(1)
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp1_r(i)
                      END FORALL
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (14)
                DO k=1,nnodes(3)
                   DO j=1,nnodes(2)
                      DO i=1,nnodes(1)
                         tmp1_r(i)=wpin(i-14,j,k)*Mask(1)+wpin(i-13,j,k)*Mask(2)+wpin(i-12,j,k)*Mask(3)+wpin(i-11,j,k)*Mask(4)+     &
                         &         wpin(i-10,j,k)*Mask(5)+wpin(i-9,j,k)*Mask(6)+wpin(i-8,j,k)*Mask(7)+wpin(i-7,j,k)*Mask(8)+        &
                         &         wpin(i-6,j,k)*Mask(9)+wpin(i-5,j,k)*Mask(10)+wpin(i-4,j,k)*Mask(11)+wpin(i-3,j,k)*Mask(12)+      &
                         &         wpin(i-2,j,k)*Mask(13)+wpin(i-1,j,k)*Mask(14)+wpin(i,j,k)*Mask(15)+wpin(i+1,j,k)*Mask(16)+       &
                         &         wpin(i+2,j,k)*Mask(17)+wpin(i+3,j,k)*Mask(18)+wpin(i+4,j,k)*Mask(19)+wpin(i+5,j,k)*Mask(20)+     &
                         &         wpin(i+6,j,k)*Mask(21)+wpin(i+7,j,k)*Mask(22)+wpin(i+8,j,k)*Mask(23)+wpin(i+9,j,k)*Mask(24)+     &
                         &         wpin(i+10,j,k)*Mask(25)+wpin(i+11,j,k)*Mask(26)+wpin(i+12,j,k)*Mask(27)+wpin(i+13,j,k)*Mask(28)+ &
                         &         wpin(i+14,j,k)*Mask(29)
                      ENDDO !i=1,nnodes(1)
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp1_r(i)
                      END FORALL
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (15)
                DO k=1,nnodes(3)
                   DO j=1,nnodes(2)
                      DO i=1,nnodes(1)
                         tmp1_r(i)=wpin(i-15,j,k)*Mask(1)+wpin(i-14,j,k)*Mask(2)+wpin(i-13,j,k)*Mask(3)+wpin(i-12,j,k)*Mask(4)+     &
                         &         wpin(i-11,j,k)*Mask(5)+wpin(i-10,j,k)*Mask(6)+wpin(i-9,j,k)*Mask(7)+wpin(i-8,j,k)*Mask(8)+       &
                         &         wpin(i-7,j,k)*Mask(9)+wpin(i-6,j,k)*Mask(10)+wpin(i-5,j,k)*Mask(11)+wpin(i-4,j,k)*Mask(12)+      &
                         &         wpin(i-3,j,k)*Mask(13)+wpin(i-2,j,k)*Mask(14)+wpin(i-1,j,k)*Mask(15)+wpin(i,j,k)*Mask(16)+       &
                         &         wpin(i+1,j,k)*Mask(17)+wpin(i+2,j,k)*Mask(18)+wpin(i+3,j,k)*Mask(19)+wpin(i+4,j,k)*Mask(20)+     &
                         &         wpin(i+5,j,k)*Mask(21)+wpin(i+6,j,k)*Mask(22)+wpin(i+7,j,k)*Mask(23)+wpin(i+8,j,k)*Mask(24)+     &
                         &         wpin(i+9,j,k)*Mask(25)+wpin(i+10,j,k)*Mask(26)+wpin(i+11,j,k)*Mask(27)+wpin(i+12,j,k)*Mask(28)+  &
                         &         wpin(i+13,j,k)*Mask(29)+wpin(i+14,j,k)*Mask(30)+wpin(i+15,j,k)*Mask(31)
                      ENDDO !i=1,nnodes(1)
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp1_r(i)
                      END FORALL
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (16)
                DO k=1,nnodes(3)
                   DO j=1,nnodes(2)
                      DO i=1,nnodes(1)
                         tmp1_r(i)=wpin(i-16,j,k)*Mask(1)+wpin(i-15,j,k)*Mask(2)+wpin(i-14,j,k)*Mask(3)+wpin(i-13,j,k)*Mask(4)+     &
                         &         wpin(i-12,j,k)*Mask(5)+wpin(i-11,j,k)*Mask(6)+wpin(i-10,j,k)*Mask(7)+wpin(i-9,j,k)*Mask(8)+      &
                         &         wpin(i-8,j,k)*Mask(9)+wpin(i-7,j,k)*Mask(10)+wpin(i-6,j,k)*Mask(11)+wpin(i-5,j,k)*Mask(12)+      &
                         &         wpin(i-4,j,k)*Mask(13)+wpin(i-3,j,k)*Mask(14)+wpin(i-2,j,k)*Mask(15)+wpin(i-1,j,k)*Mask(16)+     &
                         &         wpin(i,j,k)*Mask(17)+wpin(i+1,j,k)*Mask(18)+wpin(i+2,j,k)*Mask(19)+wpin(i+3,j,k)*Mask(20)+       &
                         &         wpin(i+4,j,k)*Mask(21)+wpin(i+5,j,k)*Mask(22)+wpin(i+6,j,k)*Mask(23)+wpin(i+7,j,k)*Mask(24)+     &
                         &         wpin(i+8,j,k)*Mask(25)+wpin(i+9,j,k)*Mask(26)+wpin(i+10,j,k)*Mask(27)+wpin(i+11,j,k)*Mask(28)+   &
                         &         wpin(i+12,j,k)*Mask(29)+wpin(i+13,j,k)*Mask(30)+wpin(i+14,j,k)*Mask(31)+wpin(i+15,j,k)*Mask(32)+ &
                         &         wpin(i+16,j,k)*Mask(33)
                      ENDDO !i=1,nnodes(1)
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp1_r(i)
                      END FORALL
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (17)
                DO k=1,nnodes(3)
                   DO j=1,nnodes(2)
                      DO i=1,nnodes(1)
                         tmp1_r(i)=wpin(i-17,j,k)*Mask(1)+wpin(i-16,j,k)*Mask(2)+wpin(i-15,j,k)*Mask(3)+wpin(i-14,j,k)*Mask(4)+     &
                         &         wpin(i-13,j,k)*Mask(5)+wpin(i-12,j,k)*Mask(6)+wpin(i-11,j,k)*Mask(7)+wpin(i-10,j,k)*Mask(8)+     &
                         &         wpin(i-9,j,k)*Mask(9)+wpin(i-8,j,k)*Mask(10)+wpin(i-7,j,k)*Mask(11)+wpin(i-6,j,k)*Mask(12)+      &
                         &         wpin(i-5,j,k)*Mask(13)+wpin(i-4,j,k)*Mask(14)+wpin(i-3,j,k)*Mask(15)+wpin(i-2,j,k)*Mask(16)+     &
                         &         wpin(i-1,j,k)*Mask(17)+wpin(i,j,k)*Mask(18)+wpin(i+1,j,k)*Mask(19)+wpin(i+2,j,k)*Mask(20)+       &
                         &         wpin(i+3,j,k)*Mask(21)+wpin(i+4,j,k)*Mask(22)+wpin(i+5,j,k)*Mask(23)+wpin(i+6,j,k)*Mask(24)+     &
                         &         wpin(i+7,j,k)*Mask(25)+wpin(i+8,j,k)*Mask(26)+wpin(i+9,j,k)*Mask(27)+wpin(i+10,j,k)*Mask(28)+    &
                         &         wpin(i+11,j,k)*Mask(29)+wpin(i+12,j,k)*Mask(30)+wpin(i+13,j,k)*Mask(31)+wpin(i+14,j,k)*Mask(32)+ &
                         &         wpin(i+15,j,k)*Mask(33)+wpin(i+16,j,k)*Mask(34)+wpin(i+17,j,k)*Mask(35)
                      ENDDO !i=1,nnodes(1)
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp1_r(i)
                      END FORALL
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (18)
                DO k=1,nnodes(3)
                   DO j=1,nnodes(2)
                      DO i=1,nnodes(1)
                         tmp1_r(i)=wpin(i-18,j,k)*Mask(1)+wpin(i-17,j,k)*Mask(2)+wpin(i-16,j,k)*Mask(3)+wpin(i-15,j,k)*Mask(4)+     &
                         &         wpin(i-14,j,k)*Mask(5)+wpin(i-13,j,k)*Mask(6)+wpin(i-12,j,k)*Mask(7)+wpin(i-11,j,k)*Mask(8)+     &
                         &         wpin(i-10,j,k)*Mask(9)+wpin(i-9,j,k)*Mask(10)+wpin(i-8,j,k)*Mask(11)+wpin(i-7,j,k)*Mask(12)+     &
                         &         wpin(i-6,j,k)*Mask(13)+wpin(i-5,j,k)*Mask(14)+wpin(i-4,j,k)*Mask(15)+wpin(i-3,j,k)*Mask(16)+     &
                         &         wpin(i-2,j,k)*Mask(17)+wpin(i-1,j,k)*Mask(18)+wpin(i,j,k)*Mask(19)+wpin(i+1,j,k)*Mask(20)+       &
                         &         wpin(i+2,j,k)*Mask(21)+wpin(i+3,j,k)*Mask(22)+wpin(i+4,j,k)*Mask(23)+wpin(i+5,j,k)*Mask(24)+     &
                         &         wpin(i+6,j,k)*Mask(25)+wpin(i+7,j,k)*Mask(26)+wpin(i+8,j,k)*Mask(27)+wpin(i+9,j,k)*Mask(28)+     &
                         &         wpin(i+10,j,k)*Mask(29)+wpin(i+11,j,k)*Mask(30)+wpin(i+12,j,k)*Mask(31)+wpin(i+13,j,k)*Mask(32)+ &
                         &         wpin(i+14,j,k)*Mask(33)+wpin(i+15,j,k)*Mask(34)+wpin(i+16,j,k)*Mask(35)+wpin(i+17,j,k)*Mask(36)+ &
                         &         wpin(i+18,j,k)*Mask(37)
                      ENDDO !i=1,nnodes(1)
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp1_r(i)
                      END FORALL
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (19)
                DO k=1,nnodes(3)
                   DO j=1,nnodes(2)
                      DO i=1,nnodes(1)
                         tmp1_r(i)=wpin(i-19,j,k)*Mask(1)+wpin(i-18,j,k)*Mask(2)+wpin(i-17,j,k)*Mask(3)+wpin(i-16,j,k)*Mask(4)+     &
                         &         wpin(i-15,j,k)*Mask(5)+wpin(i-14,j,k)*Mask(6)+wpin(i-13,j,k)*Mask(7)+wpin(i-12,j,k)*Mask(8)+     &
                         &         wpin(i-11,j,k)*Mask(9)+wpin(i-10,j,k)*Mask(10)+wpin(i-9,j,k)*Mask(11)+wpin(i-8,j,k)*Mask(12)+    &
                         &         wpin(i-7,j,k)*Mask(13)+wpin(i-6,j,k)*Mask(14)+wpin(i-5,j,k)*Mask(15)+wpin(i-4,j,k)*Mask(16)+     &
                         &         wpin(i-3,j,k)*Mask(17)+wpin(i-2,j,k)*Mask(18)+wpin(i-1,j,k)*Mask(19)+wpin(i,j,k)*Mask(20)+       &
                         &         wpin(i+1,j,k)*Mask(21)+wpin(i+2,j,k)*Mask(22)+wpin(i+3,j,k)*Mask(23)+wpin(i+4,j,k)*Mask(24)+     &
                         &         wpin(i+5,j,k)*Mask(25)+wpin(i+6,j,k)*Mask(26)+wpin(i+7,j,k)*Mask(27)+wpin(i+8,j,k)*Mask(28)+     &
                         &         wpin(i+9,j,k)*Mask(29)+wpin(i+10,j,k)*Mask(30)+wpin(i+11,j,k)*Mask(31)+wpin(i+12,j,k)*Mask(32)+  &
                         &         wpin(i+13,j,k)*Mask(33)+wpin(i+14,j,k)*Mask(34)+wpin(i+15,j,k)*Mask(35)+wpin(i+16,j,k)*Mask(36)+ &
                         &         wpin(i+17,j,k)*Mask(37)+wpin(i+18,j,k)*Mask(38)+wpin(i+19,j,k)*Mask(39)
                      ENDDO !i=1,nnodes(1)
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp1_r(i)
                      END FORALL
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (20)
                DO k=1,nnodes(3)
                   DO j=1,nnodes(2)
                      DO i=1,nnodes(1)
                         tmp1_r(i)=wpin(i-20,j,k)*Mask(1)+wpin(i-19,j,k)*Mask(2)+wpin(i-18,j,k)*Mask(3)+wpin(i-17,j,k)*Mask(4)+     &
                         &         wpin(i-16,j,k)*Mask(5)+wpin(i-15,j,k)*Mask(6)+wpin(i-14,j,k)*Mask(7)+wpin(i-13,j,k)*Mask(8)+     &
                         &         wpin(i-12,j,k)*Mask(9)+wpin(i-11,j,k)*Mask(10)+wpin(i-10,j,k)*Mask(11)+wpin(i-9,j,k)*Mask(12)+   &
                         &         wpin(i-8,j,k)*Mask(13)+wpin(i-7,j,k)*Mask(14)+wpin(i-6,j,k)*Mask(15)+wpin(i-5,j,k)*Mask(16)+     &
                         &         wpin(i-4,j,k)*Mask(17)+wpin(i-3,j,k)*Mask(18)+wpin(i-2,j,k)*Mask(19)+wpin(i-1,j,k)*Mask(20)+     &
                         &         wpin(i,j,k)*Mask(21)+wpin(i+1,j,k)*Mask(22)+wpin(i+2,j,k)*Mask(23)+wpin(i+3,j,k)*Mask(24)+       &
                         &         wpin(i+4,j,k)*Mask(25)+wpin(i+5,j,k)*Mask(26)+wpin(i+6,j,k)*Mask(27)+wpin(i+7,j,k)*Mask(28)+     &
                         &         wpin(i+8,j,k)*Mask(29)+wpin(i+9,j,k)*Mask(30)+wpin(i+10,j,k)*Mask(31)+wpin(i+11,j,k)*Mask(32)+   &
                         &         wpin(i+12,j,k)*Mask(33)+wpin(i+13,j,k)*Mask(34)+wpin(i+14,j,k)*Mask(35)+wpin(i+15,j,k)*Mask(36)+ &
                         &         wpin(i+16,j,k)*Mask(37)+wpin(i+17,j,k)*Mask(38)+wpin(i+18,j,k)*Mask(39)+wpin(i+19,j,k)*Mask(40)+ &
                         &         wpin(i+20,j,k)*Mask(41)
                      ENDDO !i=1,nnodes(1)
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp1_r(i)
                      END FORALL
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (21)
                DO k=1,nnodes(3)
                   DO j=1,nnodes(2)
                      DO i=1,nnodes(1)
                         tmp1_r(i)=wpin(i-21,j,k)*Mask(1)+wpin(i-20,j,k)*Mask(2)+wpin(i-19,j,k)*Mask(3)+wpin(i-18,j,k)*Mask(4)+     &
                         &         wpin(i-17,j,k)*Mask(5)+wpin(i-16,j,k)*Mask(6)+wpin(i-15,j,k)*Mask(7)+wpin(i-14,j,k)*Mask(8)+     &
                         &         wpin(i-13,j,k)*Mask(9)+wpin(i-12,j,k)*Mask(10)+wpin(i-11,j,k)*Mask(11)+wpin(i-10,j,k)*Mask(12)+  &
                         &         wpin(i-9,j,k)*Mask(13)+wpin(i-8,j,k)*Mask(14)+wpin(i-7,j,k)*Mask(15)+wpin(i-6,j,k)*Mask(16)+     &
                         &         wpin(i-5,j,k)*Mask(17)+wpin(i-4,j,k)*Mask(18)+wpin(i-3,j,k)*Mask(19)+wpin(i-2,j,k)*Mask(20)+     &
                         &         wpin(i-1,j,k)*Mask(21)+wpin(i,j,k)*Mask(22)+wpin(i+1,j,k)*Mask(23)+wpin(i+2,j,k)*Mask(24)+       &
                         &         wpin(i+3,j,k)*Mask(25)+wpin(i+4,j,k)*Mask(26)+wpin(i+5,j,k)*Mask(27)+wpin(i+6,j,k)*Mask(28)+     &
                         &         wpin(i+7,j,k)*Mask(29)+wpin(i+8,j,k)*Mask(30)+wpin(i+9,j,k)*Mask(31)+wpin(i+10,j,k)*Mask(32)+    &
                         &         wpin(i+11,j,k)*Mask(33)+wpin(i+12,j,k)*Mask(34)+wpin(i+13,j,k)*Mask(35)+wpin(i+14,j,k)*Mask(36)+ &
                         &         wpin(i+15,j,k)*Mask(37)+wpin(i+16,j,k)*Mask(38)+wpin(i+17,j,k)*Mask(39)+wpin(i+18,j,k)*Mask(40)+ &
                         &         wpin(i+19,j,k)*Mask(41)+wpin(i+20,j,k)*Mask(42)+wpin(i+21,j,k)*Mask(43)
                      ENDDO !i=1,nnodes(1)
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp1_r(i)
                      END FORALL
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (22)
                DO k=1,nnodes(3)
                   DO j=1,nnodes(2)
                      DO i=1,nnodes(1)
                         tmp1_r(i)=wpin(i-22,j,k)*Mask( 1)+wpin(i-21,j,k)*Mask( 2)+wpin(i-20,j,k)*Mask( 3)+wpin(i-19,j,k)*Mask( 4)+ &
                         &         wpin(i-18,j,k)*Mask( 5)+wpin(i-17,j,k)*Mask( 6)+wpin(i-16,j,k)*Mask( 7)+wpin(i-15,j,k)*Mask( 8)+ &
                         &         wpin(i-14,j,k)*Mask( 9)+wpin(i-13,j,k)*Mask(10)+wpin(i-12,j,k)*Mask(11)+wpin(i-11,j,k)*Mask(12)+ &
                         &         wpin(i-10,j,k)*Mask(13)+wpin(i- 9,j,k)*Mask(14)+wpin(i- 8,j,k)*Mask(15)+wpin(i- 7,j,k)*Mask(16)+ &
                         &         wpin(i- 6,j,k)*Mask(17)+wpin(i- 5,j,k)*Mask(18)+wpin(i- 4,j,k)*Mask(19)+wpin(i- 3,j,k)*Mask(20)+ &
                         &         wpin(i- 2,j,k)*Mask(21)+wpin(i- 1,j,k)*Mask(22)+wpin(i   ,j,k)*Mask(23)+wpin(i+ 1,j,k)*Mask(24)+ &
                         &         wpin(i+ 2,j,k)*Mask(25)+wpin(i+ 3,j,k)*Mask(26)+wpin(i+ 4,j,k)*Mask(27)+wpin(i+ 5,j,k)*Mask(28)+ &
                         &         wpin(i+ 6,j,k)*Mask(29)+wpin(i+ 7,j,k)*Mask(30)+wpin(i+ 8,j,k)*Mask(31)+wpin(i+ 9,j,k)*Mask(32)+ &
                         &         wpin(i+10,j,k)*Mask(33)+wpin(i+11,j,k)*Mask(34)+wpin(i+12,j,k)*Mask(35)+wpin(i+13,j,k)*Mask(36)+ &
                         &         wpin(i+14,j,k)*Mask(37)+wpin(i+15,j,k)*Mask(38)+wpin(i+16,j,k)*Mask(39)+wpin(i+17,j,k)*Mask(40)+ &
                         &         wpin(i+18,j,k)*Mask(41)+wpin(i+19,j,k)*Mask(42)+wpin(i+20,j,k)*Mask(43)+wpin(i+21,j,k)*Mask(44)+ &
                         &         wpin(i+22,j,k)*Mask(45)
                      ENDDO !i=1,nnodes(1)
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp1_r(i)
                      END FORALL
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (23)
                DO k=1,nnodes(3)
                   DO j=1,nnodes(2)
                      DO i=1,nnodes(1)
                         tmp1_r(i)=wpin(i-23,j,k)*Mask( 1)+wpin(i-22,j,k)*Mask( 2)+wpin(i-21,j,k)*Mask( 3)+wpin(i-20,j,k)*Mask( 4)+ &
                         &         wpin(i-19,j,k)*Mask( 5)+wpin(i-18,j,k)*Mask( 6)+wpin(i-17,j,k)*Mask( 7)+wpin(i-16,j,k)*Mask( 8)+ &
                         &         wpin(i-15,j,k)*Mask( 9)+wpin(i-14,j,k)*Mask(10)+wpin(i-13,j,k)*Mask(11)+wpin(i-12,j,k)*Mask(12)+ &
                         &         wpin(i-11,j,k)*Mask(13)+wpin(i-10,j,k)*Mask(14)+wpin(i- 9,j,k)*Mask(15)+wpin(i- 8,j,k)*Mask(16)+ &
                         &         wpin(i- 7,j,k)*Mask(17)+wpin(i- 6,j,k)*Mask(18)+wpin(i- 5,j,k)*Mask(19)+wpin(i- 4,j,k)*Mask(20)+ &
                         &         wpin(i- 3,j,k)*Mask(21)+wpin(i- 2,j,k)*Mask(22)+wpin(i- 1,j,k)*Mask(23)+wpin(i   ,j,k)*Mask(24)+ &
                         &         wpin(i+ 1,j,k)*Mask(25)+wpin(i+ 2,j,k)*Mask(26)+wpin(i+ 3,j,k)*Mask(27)+wpin(i+ 4,j,k)*Mask(28)+ &
                         &         wpin(i+ 5,j,k)*Mask(29)+wpin(i+ 6,j,k)*Mask(30)+wpin(i+ 7,j,k)*Mask(31)+wpin(i+ 8,j,k)*Mask(32)+ &
                         &         wpin(i+ 9,j,k)*Mask(33)+wpin(i+10,j,k)*Mask(34)+wpin(i+11,j,k)*Mask(35)+wpin(i+12,j,k)*Mask(36)+ &
                         &         wpin(i+13,j,k)*Mask(37)+wpin(i+14,j,k)*Mask(38)+wpin(i+15,j,k)*Mask(39)+wpin(i+16,j,k)*Mask(40)+ &
                         &         wpin(i+17,j,k)*Mask(41)+wpin(i+18,j,k)*Mask(42)+wpin(i+19,j,k)*Mask(43)+wpin(i+20,j,k)*Mask(44)+ &
                         &         wpin(i+21,j,k)*Mask(45)+wpin(i+22,j,k)*Mask(46)+wpin(i+23,j,k)*Mask(47)
                      ENDDO !i=1,nnodes(1)
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp1_r(i)
                      END FORALL
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (24)
                DO k=1,nnodes(3)
                   DO j=1,nnodes(2)
                      DO i=1,nnodes(1)
                         tmp1_r(i)=wpin(i-24,j,k)*Mask( 1)+wpin(i-23,j,k)*Mask( 2)+wpin(i-22,j,k)*Mask( 3)+wpin(i-21,j,k)*Mask( 4)+ &
                         &         wpin(i-20,j,k)*Mask( 5)+wpin(i-19,j,k)*Mask( 6)+wpin(i-18,j,k)*Mask( 7)+wpin(i-17,j,k)*Mask( 8)+ &
                         &         wpin(i-16,j,k)*Mask( 9)+wpin(i-15,j,k)*Mask(10)+wpin(i-14,j,k)*Mask(11)+wpin(i-13,j,k)*Mask(12)+ &
                         &         wpin(i-12,j,k)*Mask(13)+wpin(i-11,j,k)*Mask(14)+wpin(i-10,j,k)*Mask(15)+wpin(i- 9,j,k)*Mask(16)+ &
                         &         wpin(i- 8,j,k)*Mask(17)+wpin(i- 7,j,k)*Mask(18)+wpin(i- 6,j,k)*Mask(19)+wpin(i- 5,j,k)*Mask(20)+ &
                         &         wpin(i- 4,j,k)*Mask(21)+wpin(i- 3,j,k)*Mask(22)+wpin(i- 2,j,k)*Mask(23)+wpin(i- 1,j,k)*Mask(24)+ &
                         &         wpin(i   ,j,k)*Mask(25)+wpin(i+ 1,j,k)*Mask(26)+wpin(i+ 2,j,k)*Mask(27)+wpin(i+ 3,j,k)*Mask(28)+ &
                         &         wpin(i+ 4,j,k)*Mask(29)+wpin(i+ 5,j,k)*Mask(30)+wpin(i+ 6,j,k)*Mask(31)+wpin(i+ 7,j,k)*Mask(32)+ &
                         &         wpin(i+ 8,j,k)*Mask(33)+wpin(i+ 9,j,k)*Mask(34)+wpin(i+10,j,k)*Mask(35)+wpin(i+11,j,k)*Mask(36)+ &
                         &         wpin(i+12,j,k)*Mask(37)+wpin(i+13,j,k)*Mask(38)+wpin(i+14,j,k)*Mask(39)+wpin(i+15,j,k)*Mask(40)+ &
                         &         wpin(i+16,j,k)*Mask(41)+wpin(i+17,j,k)*Mask(42)+wpin(i+18,j,k)*Mask(43)+wpin(i+19,j,k)*Mask(44)+ &
                         &         wpin(i+20,j,k)*Mask(45)+wpin(i+21,j,k)*Mask(46)+wpin(i+22,j,k)*Mask(47)+wpin(i+23,j,k)*Mask(48)+ &
                         &         wpin(i+24,j,k)*Mask(49)
                      ENDDO !i=1,nnodes(1)
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp1_r(i)
                      END FORALL
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE DEFAULT
                m=MaskSize+1
                DO k=1,nnodes(3)
                   DO j=1,nnodes(2)
                      tmp1_r=__ZERO
                      DO i=1,nnodes(1)
                         DO l=-MaskSize,MaskSize
                            tmp1_r(i)=tmp1_r(i)+wpin(i+l,j,k)*Mask(l+m)
                         ENDDO
                      ENDDO !i=1,Nm(1)
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp1_r(i)
                      END FORALL
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             END SELECT

             DEALLOCATE(tmp1_r,STAT=info)
             or_fail_dealloc("Failed to deallocate tmp1_r!",exit_point=7777)
          CASE (2)
          !!! Y direction
             SELECT CASE (MaskSize)
             CASE (1)
                !At each pixel, I should keep three values of the first, second
                !and third row
                ALLOCATE(tmp2_r(nnodes(1),3),ll(3),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                ! This optimisation is done to avoid recomputing what has been computed once
                DO k=1,nnodes(3)
                   DO j=0,1
                      l=j+1
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,k)
                      END FORALL
                   ENDDO
                   FORALL(i=1:3) ll(i)=i
                   DO j=1,nnodes(2)
                      jp=j+1
                      l3=ll(3)
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l3)=wpin(i,jp,k)
                      END FORALL
                      l1=ll(1)
                      l2=ll(2)
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp2_r(i,l1)*Mask(1)+tmp2_r(i,l2)*Mask(2)+tmp2_r(i,l3)*Mask(3)
                      END FORALL
                      ll=CSHIFT(ll,SHIFT=1,DIM=1)
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (2)
                ALLOCATE(tmp2_r(nnodes(1),5),ll(5),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                DO k=1,nnodes(3)
                   DO j=-1,2
                      l=j+2
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,k)
                      END FORALL
                   ENDDO
                   FORALL(i=1:5) ll(i)=i
                   DO j=1,nnodes(2)
                      jp=j+2
                      l5=ll(5)
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l5)=wpin(i,jp,k)
                      END FORALL
                      l1=ll(1)
                      l2=ll(2)
                      l3=ll(3)
                      l4=ll(4)
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp2_r(i,l1)*Mask(1)+tmp2_r(i,l2)*Mask(2)+tmp2_r(i,l3)*Mask(3)+tmp2_r(i,l4)*Mask(4)+ &
                         &            tmp2_r(i,l5)*Mask(5)
                      END FORALL
                      ll=CSHIFT(ll,SHIFT=1,DIM=1)
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (3)
                ALLOCATE(tmp2_r(nnodes(1),7),ll(7),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                DO k=1,nnodes(3)
                   DO j=-2,3
                      l=j+3
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,k)
                      END FORALL
                   ENDDO
                   FORALL(i=1:7) ll(i)=i
                   DO j=1,nnodes(2)
                      jp=j+3
                      l7=ll(7)
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l7)=wpin(i,jp,k)
                      END FORALL
                      l1=ll(1)
                      l2=ll(2)
                      l3=ll(3)
                      l4=ll(4)
                      l5=ll(5)
                      l6=ll(6)
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp2_r(i,l1)*Mask(1)+tmp2_r(i,l2)*Mask(2)+tmp2_r(i,l3)*Mask(3)+tmp2_r(i,l4)*Mask(4)+ &
                         &            tmp2_r(i,l5)*Mask(5)+tmp2_r(i,l6)*Mask(6)+tmp2_r(i,l7)*Mask(7)
                      END FORALL
                      ll=CSHIFT(ll,SHIFT=1,DIM=1)
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (4)
                ALLOCATE(tmp2_r(nnodes(1),9),ll(9),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                DO k=1,nnodes(3)
                   DO j=-3,4
                      l=j+4
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,k)
                      END FORALL
                   ENDDO
                   FORALL(i=1:9) ll(i)=i
                   DO j=1,nnodes(2)
                      jp=j+4
                      l9=ll(9)
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l9)=wpin(i,jp,k)
                      END FORALL
                      l1=ll(1)
                      l2=ll(2)
                      l3=ll(3)
                      l4=ll(4)
                      l5=ll(5)
                      l6=ll(6)
                      l7=ll(7)
                      l8=ll(8)
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp2_r(i,l1)*Mask(1)+tmp2_r(i,l2)*Mask(2)+tmp2_r(i,l3)*Mask(3)+tmp2_r(i,l4)*Mask(4)+ &
                         &            tmp2_r(i,l5)*Mask(5)+tmp2_r(i,l6)*Mask(6)+tmp2_r(i,l7)*Mask(7)+tmp2_r(i,l8)*Mask(8)+ &
                         &            tmp2_r(i,l9)*Mask(9)
                      END FORALL
                      ll=CSHIFT(ll,SHIFT=1,DIM=1)
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (5)
                ALLOCATE(tmp2_r(nnodes(1),11),ll(11),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                DO k=1,nnodes(3)
                   DO j=-4,5
                      l=j+5
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,K)
                      END FORALL
                   ENDDO
                   FORALL(i=1:11) ll(i)=i
                   DO j=1,nnodes(2)
                      jp=j+5
                      l11=ll(11)
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l11)=wpin(i,jp,k)
                      END FORALL
                      l1 =ll( 1)
                      l2 =ll( 2)
                      l3 =ll( 3)
                      l4 =ll( 4)
                      l5 =ll( 5)
                      l6 =ll( 6)
                      l7 =ll( 7)
                      l8 =ll( 8)
                      l9 =ll( 9)
                      l10=ll(10)
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp2_r(i,l1 )*Mask( 1)+tmp2_r(i,l2 )*Mask( 2)+tmp2_r(i,l3 )*Mask( 3)+tmp2_r(i,l4 )*Mask( 4)+ &
                         &            tmp2_r(i,l5 )*Mask( 5)+tmp2_r(i,l6 )*Mask( 6)+tmp2_r(i,l7 )*Mask( 7)+tmp2_r(i,l8 )*Mask( 8)+ &
                         &            tmp2_r(i,l9 )*Mask( 9)+tmp2_r(i,l10)*Mask(10)+tmp2_r(i,l11)*Mask(11)
                      END FORALL
                      ll=CSHIFT(ll,SHIFT=1,DIM=1)
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (6)
                ALLOCATE(tmp2_r(nnodes(1),13),ll(13),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                DO k=1,nnodes(3)
                   DO j=-5,6
                      l=j+6
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,k)
                      END FORALL
                   ENDDO
                   FORALL(i=1:13) ll(i)=i
                   DO j=1,nnodes(2)
                      jp=j+6
                      l13=ll(13)
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l13)=wpin(i,jp,k)
                      END FORALL
                      l1 =ll( 1)
                      l2 =ll( 2)
                      l3 =ll( 3)
                      l4 =ll( 4)
                      l5 =ll( 5)
                      l6 =ll( 6)
                      l7 =ll( 7)
                      l8 =ll( 8)
                      l9 =ll( 9)
                      l10=ll(10)
                      l11=ll(11)
                      l12=ll(12)
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp2_r(i,l1 )*Mask( 1)+tmp2_r(i,l2 )*Mask( 2)+tmp2_r(i,l3 )*Mask( 3)+tmp2_r(i,l4 )*Mask( 4)+ &
                         &            tmp2_r(i,l5 )*Mask( 5)+tmp2_r(i,l6 )*Mask( 6)+tmp2_r(i,l7 )*Mask( 7)+tmp2_r(i,l8 )*Mask( 8)+ &
                         &            tmp2_r(i,l9 )*Mask( 9)+tmp2_r(i,l10)*Mask(10)+tmp2_r(i,l11)*Mask(11)+tmp2_r(i,l12)*Mask(12)+ &
                         &            tmp2_r(i,l13)*Mask(13)
                      END FORALL
                      ll=CSHIFT(ll,SHIFT=1,DIM=1)
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (7)
                ALLOCATE(tmp2_r(nnodes(1),15),ll(15),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                DO k=1,nnodes(3)
                   DO j=-6,7
                      l=j+7
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,k)
                      END FORALL
                   ENDDO
                   FORALL(i=1:15) ll(i)=i
                   DO j=1,nnodes(2)
                      jp=j+7
                      l15=ll(15)
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l15)=wpin(i,jp,k)
                      END FORALL
                      l1 =ll( 1)
                      l2 =ll( 2)
                      l3 =ll( 3)
                      l4 =ll( 4)
                      l5 =ll( 5)
                      l6 =ll( 6)
                      l7 =ll( 7)
                      l8 =ll( 8)
                      l9 =ll( 9)
                      l10=ll(10)
                      l11=ll(11)
                      l12=ll(12)
                      l13=ll(13)
                      l14=ll(14)
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp2_r(i,l1 )*Mask( 1)+tmp2_r(i,l2 )*Mask( 2)+tmp2_r(i,l3 )*Mask( 3)+tmp2_r(i,l4 )*Mask( 4)+ &
                         &            tmp2_r(i,l5 )*Mask( 5)+tmp2_r(i,l6 )*Mask( 6)+tmp2_r(i,l7 )*Mask( 7)+tmp2_r(i,l8 )*Mask( 8)+ &
                         &            tmp2_r(i,l9 )*Mask( 9)+tmp2_r(i,l10)*Mask(10)+tmp2_r(i,l11)*Mask(11)+tmp2_r(i,l12)*Mask(12)+ &
                         &            tmp2_r(i,l13)*Mask(13)+tmp2_r(i,l14)*Mask(14)+tmp2_r(i,l15)*Mask(15)
                      END FORALL
                      ll=CSHIFT(ll,SHIFT=1,DIM=1)
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (8)
                ALLOCATE(tmp2_r(nnodes(1),17),ll(17),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                DO k=1,nnodes(3)
                   DO j=-7,8
                      l=j+8
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,k)
                      END FORALL
                   ENDDO
                   FORALL(i=1:17) ll(i)=i
                   DO j=1,nnodes(2)
                      jp=j+8
                      l17=ll(17)
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l17)=wpin(i,jp,k)
                      END FORALL
                      l1 =ll( 1)
                      l2 =ll( 2)
                      l3 =ll( 3)
                      l4 =ll( 4)
                      l5 =ll( 5)
                      l6 =ll( 6)
                      l7 =ll( 7)
                      l8 =ll( 8)
                      l9 =ll( 9)
                      l10=ll(10)
                      l11=ll(11)
                      l12=ll(12)
                      l13=ll(13)
                      l14=ll(14)
                      l15=ll(15)
                      l16=ll(16)
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp2_r(i,l1 )*Mask( 1)+tmp2_r(i,l2 )*Mask( 2)+tmp2_r(i,l3 )*Mask( 3)+tmp2_r(i,l4 )*Mask( 4)+ &
                         &            tmp2_r(i,l5 )*Mask( 5)+tmp2_r(i,l6 )*Mask( 6)+tmp2_r(i,l7 )*Mask( 7)+tmp2_r(i,l8 )*Mask( 8)+ &
                         &            tmp2_r(i,l9 )*Mask( 9)+tmp2_r(i,l10)*Mask(10)+tmp2_r(i,l11)*Mask(11)+tmp2_r(i,l12)*Mask(12)+ &
                         &            tmp2_r(i,l13)*Mask(13)+tmp2_r(i,l14)*Mask(14)+tmp2_r(i,l15)*Mask(15)+tmp2_r(i,l16)*Mask(16)+ &
                         &            tmp2_r(i,l17)*Mask(17)
                      END FORALL
                      ll=CSHIFT(ll,SHIFT=1,DIM=1)
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (9)
                ALLOCATE(tmp2_r(nnodes(1),19),ll(19),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                DO k=1,nnodes(3)
                   DO j=-8,9
                      l=j+9
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,k)
                      END FORALL
                   ENDDO
                   FORALL(i=1:19) ll(i)=i
                   DO j=1,nnodes(2)
                      jp=j+9
                      l19=ll(19)
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l19)=wpin(i,jp,k)
                      END FORALL
                      l1 =ll( 1)
                      l2 =ll( 2)
                      l3 =ll( 3)
                      l4 =ll( 4)
                      l5 =ll( 5)
                      l6 =ll( 6)
                      l7 =ll( 7)
                      l8 =ll( 8)
                      l9 =ll( 9)
                      l10=ll(10)
                      l11=ll(11)
                      l12=ll(12)
                      l13=ll(13)
                      l14=ll(14)
                      l15=ll(15)
                      l16=ll(16)
                      l17=ll(17)
                      l18=ll(18)
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp2_r(i,l1 )*Mask( 1)+tmp2_r(i,l2 )*Mask( 2)+tmp2_r(i,l3 )*Mask( 3)+tmp2_r(i,l4 )*Mask( 4)+ &
                         &            tmp2_r(i,l5 )*Mask( 5)+tmp2_r(i,l6 )*Mask( 6)+tmp2_r(i,l7 )*Mask( 7)+tmp2_r(i,l8 )*Mask( 8)+ &
                         &            tmp2_r(i,l9 )*Mask( 9)+tmp2_r(i,l10)*Mask(10)+tmp2_r(i,l11)*Mask(11)+tmp2_r(i,l12)*Mask(12)+ &
                         &            tmp2_r(i,l13)*Mask(13)+tmp2_r(i,l14)*Mask(14)+tmp2_r(i,l15)*Mask(15)+tmp2_r(i,l16)*Mask(16)+ &
                         &            tmp2_r(i,l17)*Mask(17)+tmp2_r(i,l18)*Mask(18)+tmp2_r(i,l19)*Mask(19)
                      END FORALL
                      ll=CSHIFT(ll,SHIFT=1,DIM=1)
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (10)
                ALLOCATE(tmp2_r(nnodes(1),21),ll(21),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                DO k=1,nnodes(3)
                   DO j=-9,10
                      l=j+10
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,k)
                      END FORALL
                   ENDDO
                   FORALL(i=1:21) ll(i)=i
                   DO j=1,nnodes(2)
                      jp=j+10
                      l21=ll(21)
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l21)=wpin(i,jp,k)
                      END FORALL
                      l1 =ll( 1)
                      l2 =ll( 2)
                      l3 =ll( 3)
                      l4 =ll( 4)
                      l5 =ll( 5)
                      l6 =ll( 6)
                      l7 =ll( 7)
                      l8 =ll( 8)
                      l9 =ll( 9)
                      l10=ll(10)
                      l11=ll(11)
                      l12=ll(12)
                      l13=ll(13)
                      l14=ll(14)
                      l15=ll(15)
                      l16=ll(16)
                      l17=ll(17)
                      l18=ll(18)
                      l19=ll(19)
                      l20=ll(20)
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp2_r(i,l1 )*Mask( 1)+tmp2_r(i,l2 )*Mask( 2)+tmp2_r(i,l3 )*Mask( 3)+tmp2_r(i,l4 )*Mask( 4)+ &
                         &            tmp2_r(i,l5 )*Mask( 5)+tmp2_r(i,l6 )*Mask( 6)+tmp2_r(i,l7 )*Mask( 7)+tmp2_r(i,l8 )*Mask( 8)+ &
                         &            tmp2_r(i,l9 )*Mask( 9)+tmp2_r(i,l10)*Mask(10)+tmp2_r(i,l11)*Mask(11)+tmp2_r(i,l12)*Mask(12)+ &
                         &            tmp2_r(i,l13)*Mask(13)+tmp2_r(i,l14)*Mask(14)+tmp2_r(i,l15)*Mask(15)+tmp2_r(i,l16)*Mask(16)+ &
                         &            tmp2_r(i,l17)*Mask(17)+tmp2_r(i,l18)*Mask(18)+tmp2_r(i,l19)*Mask(19)+tmp2_r(i,l20)*Mask(20)+ &
                         &            tmp2_r(i,l21)*Mask(21)
                      END FORALL
                      ll=CSHIFT(ll,SHIFT=1,DIM=1)
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (11)
                ALLOCATE(tmp2_r(nnodes(1),23),ll(23),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                DO k=1,nnodes(3)
                   DO j=-10,11
                      l=j+11
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,k)
                      END FORALL
                   ENDDO
                   FORALL(i=1:23) ll(i)=i
                   DO j=1,nnodes(2)
                      jp=j+11
                      l=ll(23)
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,jp,k)
                      END FORALL
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp2_r(i,ll( 1))*Mask( 1)+tmp2_r(i,ll( 2))*Mask( 2)+tmp2_r(i,ll( 3))*Mask( 3)+tmp2_r(i,ll( 4))*Mask( 4)+ &
                         &            tmp2_r(i,ll( 5))*Mask( 5)+tmp2_r(i,ll( 6))*Mask( 6)+tmp2_r(i,ll( 7))*Mask( 7)+tmp2_r(i,ll( 8))*Mask( 8)+ &
                         &            tmp2_r(i,ll( 9))*Mask( 9)+tmp2_r(i,ll(10))*Mask(10)+tmp2_r(i,ll(11))*Mask(11)+tmp2_r(i,ll(12))*Mask(12)+ &
                         &            tmp2_r(i,ll(13))*Mask(13)+tmp2_r(i,ll(14))*Mask(14)+tmp2_r(i,ll(15))*Mask(15)+tmp2_r(i,ll(16))*Mask(16)+ &
                         &            tmp2_r(i,ll(17))*Mask(17)+tmp2_r(i,ll(18))*Mask(18)+tmp2_r(i,ll(19))*Mask(19)+tmp2_r(i,ll(20))*Mask(20)+ &
                         &            tmp2_r(i,ll(21))*Mask(21)+tmp2_r(i,ll(22))*Mask(22)+tmp2_r(i,ll(23))*Mask(23)
                      END FORALL
                      ll=CSHIFT(ll,SHIFT=1,DIM=1)
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (12)
                ALLOCATE(tmp2_r(nnodes(1),25),ll(25),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                DO k=1,nnodes(3)
                   DO j=-11,12
                      l=j+12
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,k)
                      END FORALL
                   ENDDO
                   FORALL(i=1:25) ll(i)=i
                   DO j=1,nnodes(2)
                      jp=j+12
                      l=ll(25)
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,jp,k)
                      END FORALL
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp2_r(i,ll( 1))*Mask( 1)+tmp2_r(i,ll( 2))*Mask( 2)+tmp2_r(i,ll( 3))*Mask( 3)+tmp2_r(i,ll( 4))*Mask( 4)+ &
                         &            tmp2_r(i,ll( 5))*Mask( 5)+tmp2_r(i,ll( 6))*Mask( 6)+tmp2_r(i,ll( 7))*Mask( 7)+tmp2_r(i,ll( 8))*Mask( 8)+ &
                         &            tmp2_r(i,ll( 9))*Mask( 9)+tmp2_r(i,ll(10))*Mask(10)+tmp2_r(i,ll(11))*Mask(11)+tmp2_r(i,ll(12))*Mask(12)+ &
                         &            tmp2_r(i,ll(13))*Mask(13)+tmp2_r(i,ll(14))*Mask(14)+tmp2_r(i,ll(15))*Mask(15)+tmp2_r(i,ll(16))*Mask(16)+ &
                         &            tmp2_r(i,ll(17))*Mask(17)+tmp2_r(i,ll(18))*Mask(18)+tmp2_r(i,ll(19))*Mask(19)+tmp2_r(i,ll(20))*Mask(20)+ &
                         &            tmp2_r(i,ll(21))*Mask(21)+tmp2_r(i,ll(22))*Mask(22)+tmp2_r(i,ll(23))*Mask(23)+tmp2_r(i,ll(24))*Mask(24)+ &
                         &            tmp2_r(i,ll(25))*Mask(25)
                      END FORALL
                      ll=CSHIFT(ll,SHIFT=1,DIM=1)
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (13)
                ALLOCATE(tmp2_r(nnodes(1),27),ll(27),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                DO k=1,nnodes(3)
                   DO j=-12,13
                      l=j+13
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,k)
                      END FORALL
                   ENDDO
                   FORALL(i=1:27) ll(i)=i
                   DO j=1,nnodes(2)
                      jp=j+13
                      l=ll(27)
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,jp,k)
                      END FORALL
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp2_r(i,ll( 1))*Mask( 1)+tmp2_r(i,ll( 2))*Mask( 2)+tmp2_r(i,ll( 3))*Mask( 3)+tmp2_r(i,ll( 4))*Mask( 4)+ &
                         &            tmp2_r(i,ll( 5))*Mask( 5)+tmp2_r(i,ll( 6))*Mask( 6)+tmp2_r(i,ll( 7))*Mask( 7)+tmp2_r(i,ll( 8))*Mask( 8)+ &
                         &            tmp2_r(i,ll( 9))*Mask( 9)+tmp2_r(i,ll(10))*Mask(10)+tmp2_r(i,ll(11))*Mask(11)+tmp2_r(i,ll(12))*Mask(12)+ &
                         &            tmp2_r(i,ll(13))*Mask(13)+tmp2_r(i,ll(14))*Mask(14)+tmp2_r(i,ll(15))*Mask(15)+tmp2_r(i,ll(16))*Mask(16)+ &
                         &            tmp2_r(i,ll(17))*Mask(17)+tmp2_r(i,ll(18))*Mask(18)+tmp2_r(i,ll(19))*Mask(19)+tmp2_r(i,ll(20))*Mask(20)+ &
                         &            tmp2_r(i,ll(21))*Mask(21)+tmp2_r(i,ll(22))*Mask(22)+tmp2_r(i,ll(23))*Mask(23)+tmp2_r(i,ll(24))*Mask(24)+ &
                         &            tmp2_r(i,ll(25))*Mask(25)+tmp2_r(i,ll(26))*Mask(26)+tmp2_r(i,ll(27))*Mask(27)
                      END FORALL
                      ll=CSHIFT(ll,SHIFT=1,DIM=1)
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (14)
                ALLOCATE(tmp2_r(nnodes(1),29),ll(29),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                DO k=1,nnodes(3)
                   DO j=-13,14
                      l=j+14
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,k)
                      END FORALL
                   ENDDO
                   FORALL(i=1:29) ll(i)=i
                   DO j=1,nnodes(2)
                      jp=j+14
                      l=ll(29)
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,jp,k)
                      END FORALL
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp2_r(i,ll( 1))*Mask( 1)+tmp2_r(i,ll( 2))*Mask( 2)+tmp2_r(i,ll( 3))*Mask( 3)+tmp2_r(i,ll( 4))*Mask( 4)+ &
                         &            tmp2_r(i,ll( 5))*Mask( 5)+tmp2_r(i,ll( 6))*Mask( 6)+tmp2_r(i,ll( 7))*Mask( 7)+tmp2_r(i,ll( 8))*Mask( 8)+ &
                         &            tmp2_r(i,ll( 9))*Mask( 9)+tmp2_r(i,ll(10))*Mask(10)+tmp2_r(i,ll(11))*Mask(11)+tmp2_r(i,ll(12))*Mask(12)+ &
                         &            tmp2_r(i,ll(13))*Mask(13)+tmp2_r(i,ll(14))*Mask(14)+tmp2_r(i,ll(15))*Mask(15)+tmp2_r(i,ll(16))*Mask(16)+ &
                         &            tmp2_r(i,ll(17))*Mask(17)+tmp2_r(i,ll(18))*Mask(18)+tmp2_r(i,ll(19))*Mask(19)+tmp2_r(i,ll(20))*Mask(20)+ &
                         &            tmp2_r(i,ll(21))*Mask(21)+tmp2_r(i,ll(22))*Mask(22)+tmp2_r(i,ll(23))*Mask(23)+tmp2_r(i,ll(24))*Mask(24)+ &
                         &            tmp2_r(i,ll(25))*Mask(25)+tmp2_r(i,ll(26))*Mask(26)+tmp2_r(i,ll(27))*Mask(27)+tmp2_r(i,ll(28))*Mask(28)+ &
                         &            tmp2_r(i,ll(29))*Mask(29)
                      END FORALL
                      ll=CSHIFT(ll,SHIFT=1,DIM=1)
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (15)
                ALLOCATE(tmp2_r(nnodes(1),31),ll(31),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                DO k=1,nnodes(3)
                   DO j=-14,15
                      l=j+15
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,k)
                      END FORALL
                   ENDDO
                   FORALL(i=1:31) ll(i)=i
                   DO j=1,nnodes(2)
                      jp=j+15
                      l=ll(31)
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,jp,k)
                      END FORALL
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp2_r(i,ll( 1))*Mask( 1)+tmp2_r(i,ll( 2))*Mask( 2)+tmp2_r(i,ll( 3))*Mask( 3)+tmp2_r(i,ll( 4))*Mask( 4)+ &
                         &            tmp2_r(i,ll( 5))*Mask( 5)+tmp2_r(i,ll( 6))*Mask( 6)+tmp2_r(i,ll( 7))*Mask( 7)+tmp2_r(i,ll( 8))*Mask( 8)+ &
                         &            tmp2_r(i,ll( 9))*Mask( 9)+tmp2_r(i,ll(10))*Mask(10)+tmp2_r(i,ll(11))*Mask(11)+tmp2_r(i,ll(12))*Mask(12)+ &
                         &            tmp2_r(i,ll(13))*Mask(13)+tmp2_r(i,ll(14))*Mask(14)+tmp2_r(i,ll(15))*Mask(15)+tmp2_r(i,ll(16))*Mask(16)+ &
                         &            tmp2_r(i,ll(17))*Mask(17)+tmp2_r(i,ll(18))*Mask(18)+tmp2_r(i,ll(19))*Mask(19)+tmp2_r(i,ll(20))*Mask(20)+ &
                         &            tmp2_r(i,ll(21))*Mask(21)+tmp2_r(i,ll(22))*Mask(22)+tmp2_r(i,ll(23))*Mask(23)+tmp2_r(i,ll(24))*Mask(24)+ &
                         &            tmp2_r(i,ll(25))*Mask(25)+tmp2_r(i,ll(26))*Mask(26)+tmp2_r(i,ll(27))*Mask(27)+tmp2_r(i,ll(28))*Mask(28)+ &
                         &            tmp2_r(i,ll(29))*Mask(29)+tmp2_r(i,ll(30))*Mask(30)+tmp2_r(i,ll(31))*Mask(31)
                      END FORALL
                      ll=CSHIFT(ll,SHIFT=1,DIM=1)
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (16)
                ALLOCATE(tmp2_r(nnodes(1),33),ll(33),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                DO k=1,nnodes(3)
                   DO j=-15,16
                      l=j+16
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,k)
                      END FORALL
                   ENDDO
                   FORALL(i=1:33) ll(i)=i
                   DO j=1,nnodes(2)
                      jp=j+16
                      l=ll(33)
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,jp,k)
                      END FORALL
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp2_r(i,ll( 1))*Mask( 1)+tmp2_r(i,ll( 2))*Mask( 2)+tmp2_r(i,ll( 3))*Mask( 3)+tmp2_r(i,ll( 4))*Mask( 4)+ &
                         &            tmp2_r(i,ll( 5))*Mask( 5)+tmp2_r(i,ll( 6))*Mask( 6)+tmp2_r(i,ll( 7))*Mask( 7)+tmp2_r(i,ll( 8))*Mask( 8)+ &
                         &            tmp2_r(i,ll( 9))*Mask( 9)+tmp2_r(i,ll(10))*Mask(10)+tmp2_r(i,ll(11))*Mask(11)+tmp2_r(i,ll(12))*Mask(12)+ &
                         &            tmp2_r(i,ll(13))*Mask(13)+tmp2_r(i,ll(14))*Mask(14)+tmp2_r(i,ll(15))*Mask(15)+tmp2_r(i,ll(16))*Mask(16)+ &
                         &            tmp2_r(i,ll(17))*Mask(17)+tmp2_r(i,ll(18))*Mask(18)+tmp2_r(i,ll(19))*Mask(19)+tmp2_r(i,ll(20))*Mask(20)+ &
                         &            tmp2_r(i,ll(21))*Mask(21)+tmp2_r(i,ll(22))*Mask(22)+tmp2_r(i,ll(23))*Mask(23)+tmp2_r(i,ll(24))*Mask(24)+ &
                         &            tmp2_r(i,ll(25))*Mask(25)+tmp2_r(i,ll(26))*Mask(26)+tmp2_r(i,ll(27))*Mask(27)+tmp2_r(i,ll(28))*Mask(28)+ &
                         &            tmp2_r(i,ll(29))*Mask(29)+tmp2_r(i,ll(30))*Mask(30)+tmp2_r(i,ll(31))*Mask(31)+tmp2_r(i,ll(32))*Mask(32)+ &
                         &            tmp2_r(i,ll(33))*Mask(33)
                      END FORALL
                      ll=CSHIFT(ll,SHIFT=1,DIM=1)
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (17)
                ALLOCATE(tmp2_r(nnodes(1),35),ll(35),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                DO k=1,nnodes(3)
                   DO j=-16,17
                      l=j+17
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,k)
                      END FORALL
                   ENDDO
                   FORALL(i=1:35) ll(i)=i
                   DO j=1,nnodes(2)
                      jp=j+17
                      l=ll(35)
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,jp,k)
                      END FORALL
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp2_r(i,ll( 1))*Mask( 1)+tmp2_r(i,ll( 2))*Mask( 2)+tmp2_r(i,ll( 3))*Mask( 3)+tmp2_r(i,ll( 4))*Mask( 4)+ &
                         &            tmp2_r(i,ll( 5))*Mask( 5)+tmp2_r(i,ll( 6))*Mask( 6)+tmp2_r(i,ll( 7))*Mask( 7)+tmp2_r(i,ll( 8))*Mask( 8)+ &
                         &            tmp2_r(i,ll( 9))*Mask( 9)+tmp2_r(i,ll(10))*Mask(10)+tmp2_r(i,ll(11))*Mask(11)+tmp2_r(i,ll(12))*Mask(12)+ &
                         &            tmp2_r(i,ll(13))*Mask(13)+tmp2_r(i,ll(14))*Mask(14)+tmp2_r(i,ll(15))*Mask(15)+tmp2_r(i,ll(16))*Mask(16)+ &
                         &            tmp2_r(i,ll(17))*Mask(17)+tmp2_r(i,ll(18))*Mask(18)+tmp2_r(i,ll(19))*Mask(19)+tmp2_r(i,ll(20))*Mask(20)+ &
                         &            tmp2_r(i,ll(21))*Mask(21)+tmp2_r(i,ll(22))*Mask(22)+tmp2_r(i,ll(23))*Mask(23)+tmp2_r(i,ll(24))*Mask(24)+ &
                         &            tmp2_r(i,ll(25))*Mask(25)+tmp2_r(i,ll(26))*Mask(26)+tmp2_r(i,ll(27))*Mask(27)+tmp2_r(i,ll(28))*Mask(28)+ &
                         &            tmp2_r(i,ll(29))*Mask(29)+tmp2_r(i,ll(30))*Mask(30)+tmp2_r(i,ll(31))*Mask(31)+tmp2_r(i,ll(32))*Mask(32)+ &
                         &            tmp2_r(i,ll(33))*Mask(33)+tmp2_r(i,ll(34))*Mask(34)+tmp2_r(i,ll(35))*Mask(35)
                      END FORALL
                      ll=CSHIFT(ll,SHIFT=1,DIM=1)
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (18)
                ALLOCATE(tmp2_r(nnodes(1),37),ll(37),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                DO k=1,nnodes(3)
                   DO j=-17,18
                      l=j+18
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,k)
                      END FORALL
                   ENDDO
                   FORALL(i=1:37) ll(i)=i
                   DO j=1,nnodes(2)
                      jp=j+18
                      l=ll(37)
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,jp,k)
                      END FORALL
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp2_r(i,ll( 1))*Mask( 1)+tmp2_r(i,ll( 2))*Mask( 2)+tmp2_r(i,ll( 3))*Mask( 3)+tmp2_r(i,ll( 4))*Mask( 4)+ &
                         &            tmp2_r(i,ll( 5))*Mask( 5)+tmp2_r(i,ll( 6))*Mask( 6)+tmp2_r(i,ll( 7))*Mask( 7)+tmp2_r(i,ll( 8))*Mask( 8)+ &
                         &            tmp2_r(i,ll( 9))*Mask( 9)+tmp2_r(i,ll(10))*Mask(10)+tmp2_r(i,ll(11))*Mask(11)+tmp2_r(i,ll(12))*Mask(12)+ &
                         &            tmp2_r(i,ll(13))*Mask(13)+tmp2_r(i,ll(14))*Mask(14)+tmp2_r(i,ll(15))*Mask(15)+tmp2_r(i,ll(16))*Mask(16)+ &
                         &            tmp2_r(i,ll(17))*Mask(17)+tmp2_r(i,ll(18))*Mask(18)+tmp2_r(i,ll(19))*Mask(19)+tmp2_r(i,ll(20))*Mask(20)+ &
                         &            tmp2_r(i,ll(21))*Mask(21)+tmp2_r(i,ll(22))*Mask(22)+tmp2_r(i,ll(23))*Mask(23)+tmp2_r(i,ll(24))*Mask(24)+ &
                         &            tmp2_r(i,ll(25))*Mask(25)+tmp2_r(i,ll(26))*Mask(26)+tmp2_r(i,ll(27))*Mask(27)+tmp2_r(i,ll(28))*Mask(28)+ &
                         &            tmp2_r(i,ll(29))*Mask(29)+tmp2_r(i,ll(30))*Mask(30)+tmp2_r(i,ll(31))*Mask(31)+tmp2_r(i,ll(32))*Mask(32)+ &
                         &            tmp2_r(i,ll(33))*Mask(33)+tmp2_r(i,ll(34))*Mask(34)+tmp2_r(i,ll(35))*Mask(35)+tmp2_r(i,ll(36))*Mask(36)+ &
                         &            tmp2_r(i,ll(37))*Mask(37)
                      END FORALL
                      ll=CSHIFT(ll,SHIFT=1,DIM=1)
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (19)
                ALLOCATE(tmp2_r(nnodes(1),39),ll(39),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                DO k=1,nnodes(3)
                   DO j=-18,19
                      l=j+19
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,k)
                      END FORALL
                   ENDDO
                   FORALL(i=1:39) ll(i)=i
                   DO j=1,nnodes(2)
                      jp=j+19
                      l=ll(39)
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,jp,k)
                      END FORALL
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp2_r(i,ll( 1))*Mask( 1)+tmp2_r(i,ll( 2))*Mask( 2)+tmp2_r(i,ll( 3))*Mask( 3)+tmp2_r(i,ll( 4))*Mask( 4)+ &
                         &            tmp2_r(i,ll( 5))*Mask( 5)+tmp2_r(i,ll( 6))*Mask( 6)+tmp2_r(i,ll( 7))*Mask( 7)+tmp2_r(i,ll( 8))*Mask( 8)+ &
                         &            tmp2_r(i,ll( 9))*Mask( 9)+tmp2_r(i,ll(10))*Mask(10)+tmp2_r(i,ll(11))*Mask(11)+tmp2_r(i,ll(12))*Mask(12)+ &
                         &            tmp2_r(i,ll(13))*Mask(13)+tmp2_r(i,ll(14))*Mask(14)+tmp2_r(i,ll(15))*Mask(15)+tmp2_r(i,ll(16))*Mask(16)+ &
                         &            tmp2_r(i,ll(17))*Mask(17)+tmp2_r(i,ll(18))*Mask(18)+tmp2_r(i,ll(19))*Mask(19)+tmp2_r(i,ll(20))*Mask(20)+ &
                         &            tmp2_r(i,ll(21))*Mask(21)+tmp2_r(i,ll(22))*Mask(22)+tmp2_r(i,ll(23))*Mask(23)+tmp2_r(i,ll(24))*Mask(24)+ &
                         &            tmp2_r(i,ll(25))*Mask(25)+tmp2_r(i,ll(26))*Mask(26)+tmp2_r(i,ll(27))*Mask(27)+tmp2_r(i,ll(28))*Mask(28)+ &
                         &            tmp2_r(i,ll(29))*Mask(29)+tmp2_r(i,ll(30))*Mask(30)+tmp2_r(i,ll(31))*Mask(31)+tmp2_r(i,ll(32))*Mask(32)+ &
                         &            tmp2_r(i,ll(33))*Mask(33)+tmp2_r(i,ll(34))*Mask(34)+tmp2_r(i,ll(35))*Mask(35)+tmp2_r(i,ll(36))*Mask(36)+ &
                         &            tmp2_r(i,ll(37))*Mask(37)+tmp2_r(i,ll(38))*Mask(38)+tmp2_r(i,ll(39))*Mask(39)
                      END FORALL
                      ll=CSHIFT(ll,SHIFT=1,DIM=1)
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (20)
                ALLOCATE(tmp2_r(nnodes(1),41),ll(41),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                DO k=1,nnodes(3)
                   DO j=-19,20
                      l=j+20
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,k)
                      END FORALL
                   ENDDO
                   FORALL(i=1:41) ll(i)=i
                   DO j=1,nnodes(2)
                      jp=j+20
                      l=ll(41)
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,jp,k)
                      END FORALL
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp2_r(i,ll( 1))*Mask( 1)+tmp2_r(i,ll( 2))*Mask( 2)+tmp2_r(i,ll( 3))*Mask( 3)+tmp2_r(i,ll( 4))*Mask( 4)+ &
                         &            tmp2_r(i,ll( 5))*Mask( 5)+tmp2_r(i,ll( 6))*Mask( 6)+tmp2_r(i,ll( 7))*Mask( 7)+tmp2_r(i,ll( 8))*Mask( 8)+ &
                         &            tmp2_r(i,ll( 9))*Mask( 9)+tmp2_r(i,ll(10))*Mask(10)+tmp2_r(i,ll(11))*Mask(11)+tmp2_r(i,ll(12))*Mask(12)+ &
                         &            tmp2_r(i,ll(13))*Mask(13)+tmp2_r(i,ll(14))*Mask(14)+tmp2_r(i,ll(15))*Mask(15)+tmp2_r(i,ll(16))*Mask(16)+ &
                         &            tmp2_r(i,ll(17))*Mask(17)+tmp2_r(i,ll(18))*Mask(18)+tmp2_r(i,ll(19))*Mask(19)+tmp2_r(i,ll(20))*Mask(20)+ &
                         &            tmp2_r(i,ll(21))*Mask(21)+tmp2_r(i,ll(22))*Mask(22)+tmp2_r(i,ll(23))*Mask(23)+tmp2_r(i,ll(24))*Mask(24)+ &
                         &            tmp2_r(i,ll(25))*Mask(25)+tmp2_r(i,ll(26))*Mask(26)+tmp2_r(i,ll(27))*Mask(27)+tmp2_r(i,ll(28))*Mask(28)+ &
                         &            tmp2_r(i,ll(29))*Mask(29)+tmp2_r(i,ll(30))*Mask(30)+tmp2_r(i,ll(31))*Mask(31)+tmp2_r(i,ll(32))*Mask(32)+ &
                         &            tmp2_r(i,ll(33))*Mask(33)+tmp2_r(i,ll(34))*Mask(34)+tmp2_r(i,ll(35))*Mask(35)+tmp2_r(i,ll(36))*Mask(36)+ &
                         &            tmp2_r(i,ll(37))*Mask(37)+tmp2_r(i,ll(38))*Mask(38)+tmp2_r(i,ll(39))*Mask(39)+tmp2_r(i,ll(40))*Mask(40)+ &
                         &            tmp2_r(i,ll(41))*Mask(41)
                      END FORALL
                      ll=CSHIFT(ll,SHIFT=1,DIM=1)
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (21)
                ALLOCATE(tmp2_r(nnodes(1),43),ll(43),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                DO k=1,nnodes(3)
                   DO j=-20,21
                      l=j+21
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,k)
                      END FORALL
                   ENDDO
                   FORALL(i=1:43) ll(i)=i
                   DO j=1,nnodes(2)
                      jp=j+21
                      l=ll(43)
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,jp,k)
                      END FORALL
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp2_r(i,ll( 1))*Mask( 1)+tmp2_r(i,ll( 2))*Mask( 2)+tmp2_r(i,ll( 3))*Mask( 3)+tmp2_r(i,ll( 4))*Mask( 4)+ &
                         &            tmp2_r(i,ll( 5))*Mask( 5)+tmp2_r(i,ll( 6))*Mask( 6)+tmp2_r(i,ll( 7))*Mask( 7)+tmp2_r(i,ll( 8))*Mask( 8)+ &
                         &            tmp2_r(i,ll( 9))*Mask( 9)+tmp2_r(i,ll(10))*Mask(10)+tmp2_r(i,ll(11))*Mask(11)+tmp2_r(i,ll(12))*Mask(12)+ &
                         &            tmp2_r(i,ll(13))*Mask(13)+tmp2_r(i,ll(14))*Mask(14)+tmp2_r(i,ll(15))*Mask(15)+tmp2_r(i,ll(16))*Mask(16)+ &
                         &            tmp2_r(i,ll(17))*Mask(17)+tmp2_r(i,ll(18))*Mask(18)+tmp2_r(i,ll(19))*Mask(19)+tmp2_r(i,ll(20))*Mask(20)+ &
                         &            tmp2_r(i,ll(21))*Mask(21)+tmp2_r(i,ll(22))*Mask(22)+tmp2_r(i,ll(23))*Mask(23)+tmp2_r(i,ll(24))*Mask(24)+ &
                         &            tmp2_r(i,ll(25))*Mask(25)+tmp2_r(i,ll(26))*Mask(26)+tmp2_r(i,ll(27))*Mask(27)+tmp2_r(i,ll(28))*Mask(28)+ &
                         &            tmp2_r(i,ll(29))*Mask(29)+tmp2_r(i,ll(30))*Mask(30)+tmp2_r(i,ll(31))*Mask(31)+tmp2_r(i,ll(32))*Mask(32)+ &
                         &            tmp2_r(i,ll(33))*Mask(33)+tmp2_r(i,ll(34))*Mask(34)+tmp2_r(i,ll(35))*Mask(35)+tmp2_r(i,ll(36))*Mask(36)+ &
                         &            tmp2_r(i,ll(37))*Mask(37)+tmp2_r(i,ll(38))*Mask(38)+tmp2_r(i,ll(39))*Mask(39)+tmp2_r(i,ll(40))*Mask(40)+ &
                         &            tmp2_r(i,ll(41))*Mask(41)+tmp2_r(i,ll(42))*Mask(42)+tmp2_r(i,ll(43))*Mask(43)
                      END FORALL
                      ll=CSHIFT(ll,SHIFT=1,DIM=1)
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (22)
                ALLOCATE(tmp2_r(nnodes(1),45),ll(45),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                DO k=1,nnodes(3)
                   DO j=-21,22
                      l=j+22
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,k)
                      END FORALL
                   ENDDO
                   FORALL(i=1:45) ll(i)=i
                   DO j=1,nnodes(2)
                      jp=j+22
                      l=ll(45)
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,jp,k)
                      END FORALL
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp2_r(i,ll( 1))*Mask( 1)+tmp2_r(i,ll( 2))*Mask( 2)+tmp2_r(i,ll( 3))*Mask( 3)+tmp2_r(i,ll( 4))*Mask( 4)+ &
                         &            tmp2_r(i,ll( 5))*Mask( 5)+tmp2_r(i,ll( 6))*Mask( 6)+tmp2_r(i,ll( 7))*Mask( 7)+tmp2_r(i,ll( 8))*Mask( 8)+ &
                         &            tmp2_r(i,ll( 9))*Mask( 9)+tmp2_r(i,ll(10))*Mask(10)+tmp2_r(i,ll(11))*Mask(11)+tmp2_r(i,ll(12))*Mask(12)+ &
                         &            tmp2_r(i,ll(13))*Mask(13)+tmp2_r(i,ll(14))*Mask(14)+tmp2_r(i,ll(15))*Mask(15)+tmp2_r(i,ll(16))*Mask(16)+ &
                         &            tmp2_r(i,ll(17))*Mask(17)+tmp2_r(i,ll(18))*Mask(18)+tmp2_r(i,ll(19))*Mask(19)+tmp2_r(i,ll(20))*Mask(20)+ &
                         &            tmp2_r(i,ll(21))*Mask(21)+tmp2_r(i,ll(22))*Mask(22)+tmp2_r(i,ll(23))*Mask(23)+tmp2_r(i,ll(24))*Mask(24)+ &
                         &            tmp2_r(i,ll(25))*Mask(25)+tmp2_r(i,ll(26))*Mask(26)+tmp2_r(i,ll(27))*Mask(27)+tmp2_r(i,ll(28))*Mask(28)+ &
                         &            tmp2_r(i,ll(29))*Mask(29)+tmp2_r(i,ll(30))*Mask(30)+tmp2_r(i,ll(31))*Mask(31)+tmp2_r(i,ll(32))*Mask(32)+ &
                         &            tmp2_r(i,ll(33))*Mask(33)+tmp2_r(i,ll(34))*Mask(34)+tmp2_r(i,ll(35))*Mask(35)+tmp2_r(i,ll(36))*Mask(36)+ &
                         &            tmp2_r(i,ll(37))*Mask(37)+tmp2_r(i,ll(38))*Mask(38)+tmp2_r(i,ll(39))*Mask(39)+tmp2_r(i,ll(40))*Mask(40)+ &
                         &            tmp2_r(i,ll(41))*Mask(41)+tmp2_r(i,ll(42))*Mask(42)+tmp2_r(i,ll(43))*Mask(43)+tmp2_r(i,ll(44))*Mask(44)+ &
                         &            tmp2_r(i,ll(45))*Mask(45)
                      END FORALL
                      ll=CSHIFT(ll,SHIFT=1,DIM=1)
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (23)
                ALLOCATE(tmp2_r(nnodes(1),47),ll(47),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                DO k=1,nnodes(3)
                   DO j=-22,23
                      l=j+23
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,k)
                      END FORALL
                   ENDDO
                   FORALL(i=1:47) ll(i)=i
                   DO j=1,nnodes(2)
                      jp=j+23
                      l=ll(47)
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,jp,k)
                      END FORALL
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp2_r(i,ll( 1))*Mask( 1)+tmp2_r(i,ll( 2))*Mask( 2)+tmp2_r(i,ll( 3))*Mask( 3)+tmp2_r(i,ll( 4))*Mask( 4)+ &
                         &            tmp2_r(i,ll( 5))*Mask( 5)+tmp2_r(i,ll( 6))*Mask( 6)+tmp2_r(i,ll( 7))*Mask( 7)+tmp2_r(i,ll( 8))*Mask( 8)+ &
                         &            tmp2_r(i,ll( 9))*Mask( 9)+tmp2_r(i,ll(10))*Mask(10)+tmp2_r(i,ll(11))*Mask(11)+tmp2_r(i,ll(12))*Mask(12)+ &
                         &            tmp2_r(i,ll(13))*Mask(13)+tmp2_r(i,ll(14))*Mask(14)+tmp2_r(i,ll(15))*Mask(15)+tmp2_r(i,ll(16))*Mask(16)+ &
                         &            tmp2_r(i,ll(17))*Mask(17)+tmp2_r(i,ll(18))*Mask(18)+tmp2_r(i,ll(19))*Mask(19)+tmp2_r(i,ll(20))*Mask(20)+ &
                         &            tmp2_r(i,ll(21))*Mask(21)+tmp2_r(i,ll(22))*Mask(22)+tmp2_r(i,ll(23))*Mask(23)+tmp2_r(i,ll(24))*Mask(24)+ &
                         &            tmp2_r(i,ll(25))*Mask(25)+tmp2_r(i,ll(26))*Mask(26)+tmp2_r(i,ll(27))*Mask(27)+tmp2_r(i,ll(28))*Mask(28)+ &
                         &            tmp2_r(i,ll(29))*Mask(29)+tmp2_r(i,ll(30))*Mask(30)+tmp2_r(i,ll(31))*Mask(31)+tmp2_r(i,ll(32))*Mask(32)+ &
                         &            tmp2_r(i,ll(33))*Mask(33)+tmp2_r(i,ll(34))*Mask(34)+tmp2_r(i,ll(35))*Mask(35)+tmp2_r(i,ll(36))*Mask(36)+ &
                         &            tmp2_r(i,ll(37))*Mask(37)+tmp2_r(i,ll(38))*Mask(38)+tmp2_r(i,ll(39))*Mask(39)+tmp2_r(i,ll(40))*Mask(40)+ &
                         &            tmp2_r(i,ll(41))*Mask(41)+tmp2_r(i,ll(42))*Mask(42)+tmp2_r(i,ll(43))*Mask(43)+tmp2_r(i,ll(44))*Mask(44)+ &
                         &            tmp2_r(i,ll(45))*Mask(45)+tmp2_r(i,ll(46))*Mask(46)+tmp2_r(i,ll(47))*Mask(47)
                      END FORALL
                      ll=CSHIFT(ll,SHIFT=1,DIM=1)
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (24)
                ALLOCATE(tmp2_r(nnodes(1),49),ll(49),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                DO k=1,nnodes(3)
                   DO j=-23,24
                      l=j+24
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,k)
                      END FORALL
                   ENDDO
                   FORALL(i=1:49) ll(i)=i
                   DO j=1,nnodes(2)
                      jp=j+24
                      l=ll(49)
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,jp,k)
                      END FORALL
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp2_r(i,ll( 1))*Mask( 1)+tmp2_r(i,ll( 2))*Mask( 2)+tmp2_r(i,ll( 3))*Mask( 3)+tmp2_r(i,ll( 4))*Mask( 4)+ &
                         &            tmp2_r(i,ll( 5))*Mask( 5)+tmp2_r(i,ll( 6))*Mask( 6)+tmp2_r(i,ll( 7))*Mask( 7)+tmp2_r(i,ll( 8))*Mask( 8)+ &
                         &            tmp2_r(i,ll( 9))*Mask( 9)+tmp2_r(i,ll(10))*Mask(10)+tmp2_r(i,ll(11))*Mask(11)+tmp2_r(i,ll(12))*Mask(12)+ &
                         &            tmp2_r(i,ll(13))*Mask(13)+tmp2_r(i,ll(14))*Mask(14)+tmp2_r(i,ll(15))*Mask(15)+tmp2_r(i,ll(16))*Mask(16)+ &
                         &            tmp2_r(i,ll(17))*Mask(17)+tmp2_r(i,ll(18))*Mask(18)+tmp2_r(i,ll(19))*Mask(19)+tmp2_r(i,ll(20))*Mask(20)+ &
                         &            tmp2_r(i,ll(21))*Mask(21)+tmp2_r(i,ll(22))*Mask(22)+tmp2_r(i,ll(23))*Mask(23)+tmp2_r(i,ll(24))*Mask(24)+ &
                         &            tmp2_r(i,ll(25))*Mask(25)+tmp2_r(i,ll(26))*Mask(26)+tmp2_r(i,ll(27))*Mask(27)+tmp2_r(i,ll(28))*Mask(28)+ &
                         &            tmp2_r(i,ll(29))*Mask(29)+tmp2_r(i,ll(30))*Mask(30)+tmp2_r(i,ll(31))*Mask(31)+tmp2_r(i,ll(32))*Mask(32)+ &
                         &            tmp2_r(i,ll(33))*Mask(33)+tmp2_r(i,ll(34))*Mask(34)+tmp2_r(i,ll(35))*Mask(35)+tmp2_r(i,ll(36))*Mask(36)+ &
                         &            tmp2_r(i,ll(37))*Mask(37)+tmp2_r(i,ll(38))*Mask(38)+tmp2_r(i,ll(39))*Mask(39)+tmp2_r(i,ll(40))*Mask(40)+ &
                         &            tmp2_r(i,ll(41))*Mask(41)+tmp2_r(i,ll(42))*Mask(42)+tmp2_r(i,ll(43))*Mask(43)+tmp2_r(i,ll(44))*Mask(44)+ &
                         &            tmp2_r(i,ll(45))*Mask(45)+tmp2_r(i,ll(46))*Mask(46)+tmp2_r(i,ll(47))*Mask(47)+tmp2_r(i,ll(48))*Mask(48)+ &
                         &            tmp2_r(i,ll(49))*Mask(49)
                      END FORALL
                      ll=CSHIFT(ll,SHIFT=1,DIM=1)
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE DEFAULT
                m=MaskSize*2+1

                ALLOCATE(tmp2_r(nnodes(1),m),ll(m),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                DO k=1,nnodes(3)
                   DO j=1-MaskSize,MaskSize
                      l=j+MaskSize
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,k)
                      END FORALL
                   ENDDO
                   FORALL(i=1:m) ll(i)=i
                   DO j=1,nnodes(2)
                      jp=j+MaskSize
                      l=ll(m)
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,jp,k)
                      END FORALL
                      DO i=1,nnodes(1)
                         wpout(i,j,k)=tmp2_r(i,l)*Mask(m)
                         DO n=1,m-1
                            l1=ll(n)
                            wpout(i,j,k)=wpout(i,j,k)+tmp2_r(i,l1)*Mask(n)
                         ENDDO
                      ENDDO
                      ll=CSHIFT(ll,SHIFT=1,DIM=1)
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             END SELECT

             DEALLOCATE(tmp2_r,ll,STAT=info)
             or_fail_dealloc("tmp2",exit_point=7777)
          CASE (3)
          !!! Z direction
             SELECT CASE (MaskSize)
             CASE (1)
                !At each pixel, I should keep three values of the first, second
                !and third row
                ALLOCATE(tmp2_r(nnodes(1),3),ll(3),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                ! This optimisation is done to avoid recomputing what has been computed once
                DO j=1,nnodes(2)
                   DO k=0,1
                      l=k+1
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,k)
                      END FORALL
                   ENDDO
                   FORALL(i=1:3) ll(i)=i
                   DO k=1,nnodes(3)
                      kp=k+1
                      l3=ll(3)
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l3)=wpin(i,j,kp)
                      END FORALL
                      l1=ll(1)
                      l2=ll(2)
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp2_r(i,l1)*Mask(1)+tmp2_r(i,l2 )*Mask(2)+tmp2_r(i,l3)*Mask(3)
                      END FORALL
                      ll=CSHIFT(ll,SHIFT=1,DIM=1)
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (2)
                ALLOCATE(tmp2_r(nnodes(1),5),ll(5),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                DO j=1,nnodes(2)
                   DO k=-1,2
                      l=k+2
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,k)
                      END FORALL
                   ENDDO
                   FORALL(i=1:5) ll(i)=i
                   DO k=1,nnodes(3)
                      kp=k+2
                      l5=ll(5)
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l5)=wpin(i,j,kp)
                      END FORALL
                      l1=ll(1)
                      l2=ll(2)
                      l3=ll(3)
                      l4=ll(4)
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp2_r(i,l1)*Mask(1)+tmp2_r(i,l2 )*Mask(2)+tmp2_r(i,l3)*Mask(3)+tmp2_r(i,l4)*Mask(4)+ &
                         &            tmp2_r(i,l5)*Mask(5)
                      END FORALL
                      ll=CSHIFT(ll,SHIFT=1,DIM=1)
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (3)
                ALLOCATE(tmp2_r(nnodes(1),7),ll(7),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                DO j=1,nnodes(2)
                   DO k=-2,3
                      l=k+3
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,k)
                      END FORALL
                   ENDDO
                   FORALL(i=1:7) ll(i)=i
                   DO k=1,nnodes(3)
                      kp=k+3
                      l7=ll(7)
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l7)=wpin(i,j,kp)
                      END FORALL
                      l1=ll(1)
                      l2=ll(2)
                      l3=ll(3)
                      l4=ll(4)
                      l5=ll(5)
                      l6=ll(6)
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp2_r(i,l1)*Mask(1)+tmp2_r(i,l2 )*Mask(2)+tmp2_r(i,l3)*Mask(3)+tmp2_r(i,l4)*Mask(4)+ &
                         &            tmp2_r(i,l5)*Mask(5)+tmp2_r(i,l6 )*Mask(6)+tmp2_r(i,l7)*Mask(7)
                      END FORALL
                      ll=CSHIFT(ll,SHIFT=1,DIM=1)
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (4)
                ALLOCATE(tmp2_r(nnodes(1),9),ll(9),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                DO j=1,nnodes(2)
                   DO k=-3,4
                      l=k+4
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,k)
                      END FORALL
                   ENDDO
                   FORALL(i=1:9) ll(i)=i
                   DO k=1,nnodes(3)
                      kp=k+4
                      l9=ll(9)
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l9)=wpin(i,j,kp)
                      END FORALL
                      l1=ll(1)
                      l2=ll(2)
                      l3=ll(3)
                      l4=ll(4)
                      l5=ll(5)
                      l6=ll(6)
                      l7=ll(7)
                      l8=ll(8)
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp2_r(i,l1)*Mask(1)+tmp2_r(i,l2 )*Mask(2)+tmp2_r(i,l3)*Mask(3)+tmp2_r(i,l4)*Mask(4)+ &
                         &            tmp2_r(i,l5)*Mask(5)+tmp2_r(i,l6 )*Mask(6)+tmp2_r(i,l7)*Mask(7)+tmp2_r(i,l8)*Mask(8)+ &
                         &            tmp2_r(i,l9)*Mask(9)
                      END FORALL
                      ll=CSHIFT(ll,SHIFT=1,DIM=1)
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (5)
                ALLOCATE(tmp2_r(nnodes(1),11),ll(11),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                DO j=1,nnodes(2)
                   DO k=-4,5
                      l=k+5
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,k)
                      END FORALL
                   ENDDO
                   FORALL(i=1:11) ll(i)=i
                   DO k=1,nnodes(3)
                      kp=k+5
                      l11=ll(11)
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l11)=wpin(i,j,kp)
                      END FORALL
                      l1 =ll( 1)
                      l2 =ll( 2)
                      l3 =ll( 3)
                      l4 =ll( 4)
                      l5 =ll( 5)
                      l6 =ll( 6)
                      l7 =ll( 7)
                      l8 =ll( 8)
                      l9 =ll( 9)
                      l10=ll(10)
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp2_r(i,l1 )*Mask( 1)+tmp2_r(i,l2 )*Mask( 2)+tmp2_r(i,l3 )*Mask( 3)+tmp2_r(i,l4 )*Mask( 4)+ &
                         &            tmp2_r(i,l5 )*Mask( 5)+tmp2_r(i,l6 )*Mask( 6)+tmp2_r(i,l7 )*Mask( 7)+tmp2_r(i,l8 )*Mask( 8)+ &
                         &            tmp2_r(i,l9 )*Mask( 9)+tmp2_r(i,l10)*Mask(10)+tmp2_r(i,l11)*Mask(11)
                      END FORALL
                      ll=CSHIFT(ll,SHIFT=1,DIM=1)
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (6)
                ALLOCATE(tmp2_r(nnodes(1),13),ll(13),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                DO j=1,nnodes(2)
                   DO k=-5,6
                      l=k+6
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,k)
                      END FORALL
                   ENDDO
                   FORALL(i=1:13) ll(i)=i
                   DO k=1,nnodes(3)
                      kp=k+6
                      l13=ll(13)
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l13)=wpin(i,j,kp)
                      END FORALL
                      l1 =ll( 1)
                      l2 =ll( 2)
                      l3 =ll( 3)
                      l4 =ll( 4)
                      l5 =ll( 5)
                      l6 =ll( 6)
                      l7 =ll( 7)
                      l8 =ll( 8)
                      l9 =ll( 9)
                      l10=ll(10)
                      l11=ll(11)
                      l12=ll(12)
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp2_r(i,l1 )*Mask( 1)+tmp2_r(i,l2 )*Mask( 2)+tmp2_r(i,l3 )*Mask( 3)+tmp2_r(i,l4 )*Mask( 4)+ &
                         &            tmp2_r(i,l5 )*Mask( 5)+tmp2_r(i,l6 )*Mask( 6)+tmp2_r(i,l7 )*Mask( 7)+tmp2_r(i,l8 )*Mask( 8)+ &
                         &            tmp2_r(i,l9 )*Mask( 9)+tmp2_r(i,l10)*Mask(10)+tmp2_r(i,l11)*Mask(11)+tmp2_r(i,l12)*Mask(12)+ &
                         &            tmp2_r(i,l13)*Mask(13)
                      END FORALL
                      ll=CSHIFT(ll,SHIFT=1,DIM=1)
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (7)
                ALLOCATE(tmp2_r(nnodes(1),15),ll(15),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                DO j=1,nnodes(2)
                   DO k=-6,7
                      l=k+7
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,k)
                      END FORALL
                   ENDDO
                   FORALL(i=1:15) ll(i)=i
                   DO k=1,nnodes(3)
                      kp=k+7
                      l15=ll(15)
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l15)=wpin(i,j,kp)
                      END FORALL
                      l1 =ll( 1)
                      l2 =ll( 2)
                      l3 =ll( 3)
                      l4 =ll( 4)
                      l5 =ll( 5)
                      l6 =ll( 6)
                      l7 =ll( 7)
                      l8 =ll( 8)
                      l9 =ll( 9)
                      l10=ll(10)
                      l11=ll(11)
                      l12=ll(12)
                      l13=ll(13)
                      l14=ll(14)
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp2_r(i,l1 )*Mask( 1)+tmp2_r(i,l2 )*Mask( 2)+tmp2_r(i,l3 )*Mask( 3)+tmp2_r(i,l4 )*Mask( 4)+ &
                         &            tmp2_r(i,l5 )*Mask( 5)+tmp2_r(i,l6 )*Mask( 6)+tmp2_r(i,l7 )*Mask( 7)+tmp2_r(i,l8 )*Mask( 8)+ &
                         &            tmp2_r(i,l9 )*Mask( 9)+tmp2_r(i,l10)*Mask(10)+tmp2_r(i,l11)*Mask(11)+tmp2_r(i,l12)*Mask(12)+ &
                         &            tmp2_r(i,l13)*Mask(13)+tmp2_r(i,l14)*Mask(14)+tmp2_r(i,l15)*Mask(15)
                      END FORALL
                      ll=CSHIFT(ll,SHIFT=1,DIM=1)
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (8)
                ALLOCATE(tmp2_r(nnodes(1),17),ll(17),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                DO j=1,nnodes(2)
                   DO k=-7,8
                      l=k+8
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,k)
                      END FORALL
                   ENDDO
                   FORALL(i=1:17) ll(i)=i
                   DO k=1,nnodes(3)
                      kp=k+8
                      l17=ll(17)
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l17)=wpin(i,j,kp)
                      END FORALL
                      l1 =ll( 1)
                      l2 =ll( 2)
                      l3 =ll( 3)
                      l4 =ll( 4)
                      l5 =ll( 5)
                      l6 =ll( 6)
                      l7 =ll( 7)
                      l8 =ll( 8)
                      l9 =ll( 9)
                      l10=ll(10)
                      l11=ll(11)
                      l12=ll(12)
                      l13=ll(13)
                      l14=ll(14)
                      l15=ll(15)
                      l16=ll(16)
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp2_r(i,l1 )*Mask( 1)+tmp2_r(i,l2 )*Mask( 2)+tmp2_r(i,l3 )*Mask( 3)+tmp2_r(i,l4 )*Mask( 4)+ &
                         &            tmp2_r(i,l5 )*Mask( 5)+tmp2_r(i,l6 )*Mask( 6)+tmp2_r(i,l7 )*Mask( 7)+tmp2_r(i,l8 )*Mask( 8)+ &
                         &            tmp2_r(i,l9 )*Mask( 9)+tmp2_r(i,l10)*Mask(10)+tmp2_r(i,l11)*Mask(11)+tmp2_r(i,l12)*Mask(12)+ &
                         &            tmp2_r(i,l13)*Mask(13)+tmp2_r(i,l14)*Mask(14)+tmp2_r(i,l15)*Mask(15)+tmp2_r(i,l16)*Mask(16)+ &
                         &            tmp2_r(i,l17)*Mask(17)
                      END FORALL
                      ll=CSHIFT(ll,SHIFT=1,DIM=1)
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (9)
                ALLOCATE(tmp2_r(nnodes(1),19),ll(19),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                DO j=1,nnodes(2)
                   DO k=-8,9
                      l=k+9
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,k)
                      END FORALL
                   ENDDO
                   FORALL(i=1:19) ll(i)=i
                   DO k=1,nnodes(3)
                      kp=k+9
                      l19=ll(19)
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l19)=wpin(i,j,kp)
                      END FORALL
                      l1 =ll( 1)
                      l2 =ll( 2)
                      l3 =ll( 3)
                      l4 =ll( 4)
                      l5 =ll( 5)
                      l6 =ll( 6)
                      l7 =ll( 7)
                      l8 =ll( 8)
                      l9 =ll( 9)
                      l10=ll(10)
                      l11=ll(11)
                      l12=ll(12)
                      l13=ll(13)
                      l14=ll(14)
                      l15=ll(15)
                      l16=ll(16)
                      l17=ll(17)
                      l18=ll(18)
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp2_r(i,l1 )*Mask( 1)+tmp2_r(i,l2 )*Mask( 2)+tmp2_r(i,l3 )*Mask( 3)+tmp2_r(i,l4 )*Mask( 4)+ &
                         &            tmp2_r(i,l5 )*Mask( 5)+tmp2_r(i,l6 )*Mask( 6)+tmp2_r(i,l7 )*Mask( 7)+tmp2_r(i,l8 )*Mask( 8)+ &
                         &            tmp2_r(i,l9 )*Mask( 9)+tmp2_r(i,l10)*Mask(10)+tmp2_r(i,l11)*Mask(11)+tmp2_r(i,l12)*Mask(12)+ &
                         &            tmp2_r(i,l13)*Mask(13)+tmp2_r(i,l14)*Mask(14)+tmp2_r(i,l15)*Mask(15)+tmp2_r(i,l16)*Mask(16)+ &
                         &            tmp2_r(i,l17)*Mask(17)+tmp2_r(i,l18)*Mask(18)+tmp2_r(i,l19)*Mask(19)
                      END FORALL
                      ll=CSHIFT(ll,SHIFT=1,DIM=1)
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (10)
                ALLOCATE(tmp2_r(nnodes(1),21),ll(21),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                DO j=1,nnodes(2)
                   DO k=-9,10
                      l=k+10
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,k)
                      END FORALL
                   ENDDO
                   FORALL(i=1:21) ll(i)=i
                   DO k=1,nnodes(3)
                      kp=k+10
                      l21=ll(21)
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l21)=wpin(i,j,kp)
                      END FORALL
                      l1 =ll( 1)
                      l2 =ll( 2)
                      l3 =ll( 3)
                      l4 =ll( 4)
                      l5 =ll( 5)
                      l6 =ll( 6)
                      l7 =ll( 7)
                      l8 =ll( 8)
                      l9 =ll( 9)
                      l10=ll(10)
                      l11=ll(11)
                      l12=ll(12)
                      l13=ll(13)
                      l14=ll(14)
                      l15=ll(15)
                      l16=ll(16)
                      l17=ll(17)
                      l18=ll(18)
                      l19=ll(19)
                      l20=ll(20)
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp2_r(i,l1 )*Mask( 1)+tmp2_r(i,l2 )*Mask( 2)+tmp2_r(i,l3 )*Mask( 3)+tmp2_r(i,l4 )*Mask( 4)+ &
                         &            tmp2_r(i,l5 )*Mask( 5)+tmp2_r(i,l6 )*Mask( 6)+tmp2_r(i,l7 )*Mask( 7)+tmp2_r(i,l8 )*Mask( 8)+ &
                         &            tmp2_r(i,l9 )*Mask( 9)+tmp2_r(i,l10)*Mask(10)+tmp2_r(i,l11)*Mask(11)+tmp2_r(i,l12)*Mask(12)+ &
                         &            tmp2_r(i,l13)*Mask(13)+tmp2_r(i,l14)*Mask(14)+tmp2_r(i,l15)*Mask(15)+tmp2_r(i,l16)*Mask(16)+ &
                         &            tmp2_r(i,l17)*Mask(17)+tmp2_r(i,l18)*Mask(18)+tmp2_r(i,l19)*Mask(19)+tmp2_r(i,l20)*Mask(20)+ &
                         &            tmp2_r(i,l21)*Mask(21)
                      END FORALL
                      ll=CSHIFT(ll,SHIFT=1,DIM=1)
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (11)
                ALLOCATE(tmp2_r(nnodes(1),23),ll(23),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                DO j=1,nnodes(2)
                   DO k=-10,11
                      l=k+11
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,k)
                      END FORALL
                   ENDDO
                   FORALL(i=1:23) ll(i)=i
                   DO k=1,nnodes(3)
                      kp=k+11
                      l=ll(23)
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,kp)
                      END FORALL
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp2_r(i,ll( 1))*Mask( 1)+tmp2_r(i,ll( 2))*Mask( 2)+tmp2_r(i,ll( 3))*Mask( 3)+tmp2_r(i,ll( 4))*Mask( 4)+ &
                         &            tmp2_r(i,ll( 5))*Mask( 5)+tmp2_r(i,ll( 6))*Mask( 6)+tmp2_r(i,ll( 7))*Mask( 7)+tmp2_r(i,ll( 8))*Mask( 8)+ &
                         &            tmp2_r(i,ll( 9))*Mask( 9)+tmp2_r(i,ll(10))*Mask(10)+tmp2_r(i,ll(11))*Mask(11)+tmp2_r(i,ll(12))*Mask(12)+ &
                         &            tmp2_r(i,ll(13))*Mask(13)+tmp2_r(i,ll(14))*Mask(14)+tmp2_r(i,ll(15))*Mask(15)+tmp2_r(i,ll(16))*Mask(16)+ &
                         &            tmp2_r(i,ll(17))*Mask(17)+tmp2_r(i,ll(18))*Mask(18)+tmp2_r(i,ll(19))*Mask(19)+tmp2_r(i,ll(20))*Mask(20)+ &
                         &            tmp2_r(i,ll(21))*Mask(21)+tmp2_r(i,ll(22))*Mask(22)+tmp2_r(i,ll(23))*Mask(23)
                      END FORALL
                      ll=CSHIFT(ll,SHIFT=1,DIM=1)
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (12)
                ALLOCATE(tmp2_r(nnodes(1),25),ll(25),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                DO j=1,nnodes(2)
                   DO k=-11,12
                      l=k+12
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,k)
                      END FORALL
                   ENDDO
                   FORALL(i=1:25) ll(i)=i
                   DO k=1,nnodes(3)
                      kp=k+12
                      l=ll(25)
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,kp)
                      END FORALL
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp2_r(i,ll( 1))*Mask( 1)+tmp2_r(i,ll( 2))*Mask( 2)+tmp2_r(i,ll( 3))*Mask( 3)+tmp2_r(i,ll( 4))*Mask( 4)+ &
                         &            tmp2_r(i,ll( 5))*Mask( 5)+tmp2_r(i,ll( 6))*Mask( 6)+tmp2_r(i,ll( 7))*Mask( 7)+tmp2_r(i,ll( 8))*Mask( 8)+ &
                         &            tmp2_r(i,ll( 9))*Mask( 9)+tmp2_r(i,ll(10))*Mask(10)+tmp2_r(i,ll(11))*Mask(11)+tmp2_r(i,ll(12))*Mask(12)+ &
                         &            tmp2_r(i,ll(13))*Mask(13)+tmp2_r(i,ll(14))*Mask(14)+tmp2_r(i,ll(15))*Mask(15)+tmp2_r(i,ll(16))*Mask(16)+ &
                         &            tmp2_r(i,ll(17))*Mask(17)+tmp2_r(i,ll(18))*Mask(18)+tmp2_r(i,ll(19))*Mask(19)+tmp2_r(i,ll(20))*Mask(20)+ &
                         &            tmp2_r(i,ll(21))*Mask(21)+tmp2_r(i,ll(22))*Mask(22)+tmp2_r(i,ll(23))*Mask(23)+tmp2_r(i,ll(24))*Mask(24)+ &
                         &            tmp2_r(i,ll(25))*Mask(25)
                      END FORALL
                      ll=CSHIFT(ll,SHIFT=1,DIM=1)
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (13)
                ALLOCATE(tmp2_r(nnodes(1),27),ll(27),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                DO j=1,nnodes(2)
                   DO k=-12,13
                      l=k+13
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,k)
                      END FORALL
                   ENDDO
                   FORALL(i=1:27) ll(i)=i
                   DO k=1,nnodes(3)
                      kp=k+13
                      l=ll(27)
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,kp)
                      END FORALL
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp2_r(i,ll( 1))*Mask( 1)+tmp2_r(i,ll( 2))*Mask( 2)+tmp2_r(i,ll( 3))*Mask( 3)+tmp2_r(i,ll( 4))*Mask( 4)+ &
                         &            tmp2_r(i,ll( 5))*Mask( 5)+tmp2_r(i,ll( 6))*Mask( 6)+tmp2_r(i,ll( 7))*Mask( 7)+tmp2_r(i,ll( 8))*Mask( 8)+ &
                         &            tmp2_r(i,ll( 9))*Mask( 9)+tmp2_r(i,ll(10))*Mask(10)+tmp2_r(i,ll(11))*Mask(11)+tmp2_r(i,ll(12))*Mask(12)+ &
                         &            tmp2_r(i,ll(13))*Mask(13)+tmp2_r(i,ll(14))*Mask(14)+tmp2_r(i,ll(15))*Mask(15)+tmp2_r(i,ll(16))*Mask(16)+ &
                         &            tmp2_r(i,ll(17))*Mask(17)+tmp2_r(i,ll(18))*Mask(18)+tmp2_r(i,ll(19))*Mask(19)+tmp2_r(i,ll(20))*Mask(20)+ &
                         &            tmp2_r(i,ll(21))*Mask(21)+tmp2_r(i,ll(22))*Mask(22)+tmp2_r(i,ll(23))*Mask(23)+tmp2_r(i,ll(24))*Mask(24)+ &
                         &            tmp2_r(i,ll(25))*Mask(25)+tmp2_r(i,ll(26))*Mask(26)+tmp2_r(i,ll(27))*Mask(27)
                      END FORALL
                      ll=CSHIFT(ll,SHIFT=1,DIM=1)
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (14)
                ALLOCATE(tmp2_r(nnodes(1),29),ll(29),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                DO j=1,nnodes(2)
                   DO k=-13,14
                      l=k+14
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,k)
                      END FORALL
                   ENDDO
                   FORALL(i=1:29) ll(i)=i
                   DO k=1,nnodes(3)
                      kp=k+14
                      l=ll(29)
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,kp)
                      END FORALL
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp2_r(i,ll( 1))*Mask( 1)+tmp2_r(i,ll( 2))*Mask( 2)+tmp2_r(i,ll( 3))*Mask( 3)+tmp2_r(i,ll( 4))*Mask( 4)+ &
                         &            tmp2_r(i,ll( 5))*Mask( 5)+tmp2_r(i,ll( 6))*Mask( 6)+tmp2_r(i,ll( 7))*Mask( 7)+tmp2_r(i,ll( 8))*Mask( 8)+ &
                         &            tmp2_r(i,ll( 9))*Mask( 9)+tmp2_r(i,ll(10))*Mask(10)+tmp2_r(i,ll(11))*Mask(11)+tmp2_r(i,ll(12))*Mask(12)+ &
                         &            tmp2_r(i,ll(13))*Mask(13)+tmp2_r(i,ll(14))*Mask(14)+tmp2_r(i,ll(15))*Mask(15)+tmp2_r(i,ll(16))*Mask(16)+ &
                         &            tmp2_r(i,ll(17))*Mask(17)+tmp2_r(i,ll(18))*Mask(18)+tmp2_r(i,ll(19))*Mask(19)+tmp2_r(i,ll(20))*Mask(20)+ &
                         &            tmp2_r(i,ll(21))*Mask(21)+tmp2_r(i,ll(22))*Mask(22)+tmp2_r(i,ll(23))*Mask(23)+tmp2_r(i,ll(24))*Mask(24)+ &
                         &            tmp2_r(i,ll(25))*Mask(25)+tmp2_r(i,ll(26))*Mask(26)+tmp2_r(i,ll(27))*Mask(27)+tmp2_r(i,ll(28))*Mask(28)+ &
                         &            tmp2_r(i,ll(29))*Mask(29)
                      END FORALL
                      ll=CSHIFT(ll,SHIFT=1,DIM=1)
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (15)
                ALLOCATE(tmp2_r(nnodes(1),31),ll(31),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                DO j=1,nnodes(2)
                   DO k=-14,15
                      l=k+15
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,k)
                      END FORALL
                   ENDDO
                   FORALL(i=1:31) ll(i)=i
                   DO k=1,nnodes(3)
                      kp=k+15
                      l=ll(31)
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,kp)
                      END FORALL
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp2_r(i,ll( 1))*Mask( 1)+tmp2_r(i,ll( 2))*Mask( 2)+tmp2_r(i,ll( 3))*Mask( 3)+tmp2_r(i,ll( 4))*Mask( 4)+ &
                         &            tmp2_r(i,ll( 5))*Mask( 5)+tmp2_r(i,ll( 6))*Mask( 6)+tmp2_r(i,ll( 7))*Mask( 7)+tmp2_r(i,ll( 8))*Mask( 8)+ &
                         &            tmp2_r(i,ll( 9))*Mask( 9)+tmp2_r(i,ll(10))*Mask(10)+tmp2_r(i,ll(11))*Mask(11)+tmp2_r(i,ll(12))*Mask(12)+ &
                         &            tmp2_r(i,ll(13))*Mask(13)+tmp2_r(i,ll(14))*Mask(14)+tmp2_r(i,ll(15))*Mask(15)+tmp2_r(i,ll(16))*Mask(16)+ &
                         &            tmp2_r(i,ll(17))*Mask(17)+tmp2_r(i,ll(18))*Mask(18)+tmp2_r(i,ll(19))*Mask(19)+tmp2_r(i,ll(20))*Mask(20)+ &
                         &            tmp2_r(i,ll(21))*Mask(21)+tmp2_r(i,ll(22))*Mask(22)+tmp2_r(i,ll(23))*Mask(23)+tmp2_r(i,ll(24))*Mask(24)+ &
                         &            tmp2_r(i,ll(25))*Mask(25)+tmp2_r(i,ll(26))*Mask(26)+tmp2_r(i,ll(27))*Mask(27)+tmp2_r(i,ll(28))*Mask(28)+ &
                         &            tmp2_r(i,ll(29))*Mask(29)+tmp2_r(i,ll(30))*Mask(30)+tmp2_r(i,ll(31))*Mask(31)
                      END FORALL
                      ll=CSHIFT(ll,SHIFT=1,DIM=1)
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (16)
                ALLOCATE(tmp2_r(nnodes(1),33),ll(33),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                DO j=1,nnodes(2)
                   DO k=-15,16
                      l=k+16
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,k)
                      END FORALL
                   ENDDO
                   FORALL(i=1:33) ll(i)=i
                   DO k=1,nnodes(3)
                      kp=k+16
                      l=ll(33)
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,kp)
                      END FORALL
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp2_r(i,ll( 1))*Mask( 1)+tmp2_r(i,ll( 2))*Mask( 2)+tmp2_r(i,ll( 3))*Mask( 3)+tmp2_r(i,ll( 4))*Mask( 4)+ &
                         &            tmp2_r(i,ll( 5))*Mask( 5)+tmp2_r(i,ll( 6))*Mask( 6)+tmp2_r(i,ll( 7))*Mask( 7)+tmp2_r(i,ll( 8))*Mask( 8)+ &
                         &            tmp2_r(i,ll( 9))*Mask( 9)+tmp2_r(i,ll(10))*Mask(10)+tmp2_r(i,ll(11))*Mask(11)+tmp2_r(i,ll(12))*Mask(12)+ &
                         &            tmp2_r(i,ll(13))*Mask(13)+tmp2_r(i,ll(14))*Mask(14)+tmp2_r(i,ll(15))*Mask(15)+tmp2_r(i,ll(16))*Mask(16)+ &
                         &            tmp2_r(i,ll(17))*Mask(17)+tmp2_r(i,ll(18))*Mask(18)+tmp2_r(i,ll(19))*Mask(19)+tmp2_r(i,ll(20))*Mask(20)+ &
                         &            tmp2_r(i,ll(21))*Mask(21)+tmp2_r(i,ll(22))*Mask(22)+tmp2_r(i,ll(23))*Mask(23)+tmp2_r(i,ll(24))*Mask(24)+ &
                         &            tmp2_r(i,ll(25))*Mask(25)+tmp2_r(i,ll(26))*Mask(26)+tmp2_r(i,ll(27))*Mask(27)+tmp2_r(i,ll(28))*Mask(28)+ &
                         &            tmp2_r(i,ll(29))*Mask(29)+tmp2_r(i,ll(30))*Mask(30)+tmp2_r(i,ll(31))*Mask(31)+tmp2_r(i,ll(32))*Mask(32)+ &
                         &            tmp2_r(i,ll(33))*Mask(33)
                      END FORALL
                      ll=CSHIFT(ll,SHIFT=1,DIM=1)
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (17)
                ALLOCATE(tmp2_r(nnodes(1),35),ll(35),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                DO j=1,nnodes(2)
                   DO k=-16,17
                      l=k+17
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,k)
                      END FORALL
                   ENDDO
                   FORALL(i=1:35) ll(i)=i
                   DO k=1,nnodes(3)
                      kp=k+17
                      l=ll(35)
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,kp)
                      END FORALL
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp2_r(i,ll( 1))*Mask( 1)+tmp2_r(i,ll( 2))*Mask( 2)+tmp2_r(i,ll( 3))*Mask( 3)+tmp2_r(i,ll( 4))*Mask( 4)+ &
                         &            tmp2_r(i,ll( 5))*Mask( 5)+tmp2_r(i,ll( 6))*Mask( 6)+tmp2_r(i,ll( 7))*Mask( 7)+tmp2_r(i,ll( 8))*Mask( 8)+ &
                         &            tmp2_r(i,ll( 9))*Mask( 9)+tmp2_r(i,ll(10))*Mask(10)+tmp2_r(i,ll(11))*Mask(11)+tmp2_r(i,ll(12))*Mask(12)+ &
                         &            tmp2_r(i,ll(13))*Mask(13)+tmp2_r(i,ll(14))*Mask(14)+tmp2_r(i,ll(15))*Mask(15)+tmp2_r(i,ll(16))*Mask(16)+ &
                         &            tmp2_r(i,ll(17))*Mask(17)+tmp2_r(i,ll(18))*Mask(18)+tmp2_r(i,ll(19))*Mask(19)+tmp2_r(i,ll(20))*Mask(20)+ &
                         &            tmp2_r(i,ll(21))*Mask(21)+tmp2_r(i,ll(22))*Mask(22)+tmp2_r(i,ll(23))*Mask(23)+tmp2_r(i,ll(24))*Mask(24)+ &
                         &            tmp2_r(i,ll(25))*Mask(25)+tmp2_r(i,ll(26))*Mask(26)+tmp2_r(i,ll(27))*Mask(27)+tmp2_r(i,ll(28))*Mask(28)+ &
                         &            tmp2_r(i,ll(29))*Mask(29)+tmp2_r(i,ll(30))*Mask(30)+tmp2_r(i,ll(31))*Mask(31)+tmp2_r(i,ll(32))*Mask(32)+ &
                         &            tmp2_r(i,ll(33))*Mask(33)+tmp2_r(i,ll(34))*Mask(34)+tmp2_r(i,ll(35))*Mask(35)
                      END FORALL
                      ll=CSHIFT(ll,SHIFT=1,DIM=1)
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (18)
                ALLOCATE(tmp2_r(nnodes(1),37),ll(37),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                DO j=1,nnodes(2)
                   DO k=-17,18
                      l=k+18
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,k)
                      END FORALL
                   ENDDO
                   FORALL(i=1:37) ll(i)=i
                   DO k=1,nnodes(3)
                      kp=k+18
                      l=ll(37)
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,kp)
                      END FORALL
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp2_r(i,ll( 1))*Mask( 1)+tmp2_r(i,ll( 2))*Mask( 2)+tmp2_r(i,ll( 3))*Mask( 3)+tmp2_r(i,ll( 4))*Mask( 4)+ &
                         &            tmp2_r(i,ll( 5))*Mask( 5)+tmp2_r(i,ll( 6))*Mask( 6)+tmp2_r(i,ll( 7))*Mask( 7)+tmp2_r(i,ll( 8))*Mask( 8)+ &
                         &            tmp2_r(i,ll( 9))*Mask( 9)+tmp2_r(i,ll(10))*Mask(10)+tmp2_r(i,ll(11))*Mask(11)+tmp2_r(i,ll(12))*Mask(12)+ &
                         &            tmp2_r(i,ll(13))*Mask(13)+tmp2_r(i,ll(14))*Mask(14)+tmp2_r(i,ll(15))*Mask(15)+tmp2_r(i,ll(16))*Mask(16)+ &
                         &            tmp2_r(i,ll(17))*Mask(17)+tmp2_r(i,ll(18))*Mask(18)+tmp2_r(i,ll(19))*Mask(19)+tmp2_r(i,ll(20))*Mask(20)+ &
                         &            tmp2_r(i,ll(21))*Mask(21)+tmp2_r(i,ll(22))*Mask(22)+tmp2_r(i,ll(23))*Mask(23)+tmp2_r(i,ll(24))*Mask(24)+ &
                         &            tmp2_r(i,ll(25))*Mask(25)+tmp2_r(i,ll(26))*Mask(26)+tmp2_r(i,ll(27))*Mask(27)+tmp2_r(i,ll(28))*Mask(28)+ &
                         &            tmp2_r(i,ll(29))*Mask(29)+tmp2_r(i,ll(30))*Mask(30)+tmp2_r(i,ll(31))*Mask(31)+tmp2_r(i,ll(32))*Mask(32)+ &
                         &            tmp2_r(i,ll(33))*Mask(33)+tmp2_r(i,ll(34))*Mask(34)+tmp2_r(i,ll(35))*Mask(35)+tmp2_r(i,ll(36))*Mask(36)+ &
                         &            tmp2_r(i,ll(37))*Mask(37)
                      END FORALL
                      ll=CSHIFT(ll,SHIFT=1,DIM=1)
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (19)
                ALLOCATE(tmp2_r(nnodes(1),39),ll(39),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                DO j=1,nnodes(2)
                   DO k=-18,19
                      l=k+19
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,k)
                      END FORALL
                   ENDDO
                   FORALL(i=1:39) ll(i)=i
                   DO k=1,nnodes(3)
                      kp=k+19
                      l=ll(39)
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,kp)
                      END FORALL
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp2_r(i,ll( 1))*Mask( 1)+tmp2_r(i,ll( 2))*Mask( 2)+tmp2_r(i,ll( 3))*Mask( 3)+tmp2_r(i,ll( 4))*Mask( 4)+ &
                         &            tmp2_r(i,ll( 5))*Mask( 5)+tmp2_r(i,ll( 6))*Mask( 6)+tmp2_r(i,ll( 7))*Mask( 7)+tmp2_r(i,ll( 8))*Mask( 8)+ &
                         &            tmp2_r(i,ll( 9))*Mask( 9)+tmp2_r(i,ll(10))*Mask(10)+tmp2_r(i,ll(11))*Mask(11)+tmp2_r(i,ll(12))*Mask(12)+ &
                         &            tmp2_r(i,ll(13))*Mask(13)+tmp2_r(i,ll(14))*Mask(14)+tmp2_r(i,ll(15))*Mask(15)+tmp2_r(i,ll(16))*Mask(16)+ &
                         &            tmp2_r(i,ll(17))*Mask(17)+tmp2_r(i,ll(18))*Mask(18)+tmp2_r(i,ll(19))*Mask(19)+tmp2_r(i,ll(20))*Mask(20)+ &
                         &            tmp2_r(i,ll(21))*Mask(21)+tmp2_r(i,ll(22))*Mask(22)+tmp2_r(i,ll(23))*Mask(23)+tmp2_r(i,ll(24))*Mask(24)+ &
                         &            tmp2_r(i,ll(25))*Mask(25)+tmp2_r(i,ll(26))*Mask(26)+tmp2_r(i,ll(27))*Mask(27)+tmp2_r(i,ll(28))*Mask(28)+ &
                         &            tmp2_r(i,ll(29))*Mask(29)+tmp2_r(i,ll(30))*Mask(30)+tmp2_r(i,ll(31))*Mask(31)+tmp2_r(i,ll(32))*Mask(32)+ &
                         &            tmp2_r(i,ll(33))*Mask(33)+tmp2_r(i,ll(34))*Mask(34)+tmp2_r(i,ll(35))*Mask(35)+tmp2_r(i,ll(36))*Mask(36)+ &
                         &            tmp2_r(i,ll(37))*Mask(37)+tmp2_r(i,ll(38))*Mask(38)+tmp2_r(i,ll(39))*Mask(39)
                      END FORALL
                      ll=CSHIFT(ll,SHIFT=1,DIM=1)
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (20)
                ALLOCATE(tmp2_r(nnodes(1),41),ll(41),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                DO j=1,nnodes(2)
                   DO k=-19,20
                      l=k+20
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,k)
                      END FORALL
                   ENDDO
                   FORALL(i=1:41) ll(i)=i
                   DO k=1,nnodes(3)
                      kp=k+20
                      l=ll(41)
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,kp)
                      END FORALL
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp2_r(i,ll( 1))*Mask( 1)+tmp2_r(i,ll( 2))*Mask( 2)+tmp2_r(i,ll( 3))*Mask( 3)+tmp2_r(i,ll( 4))*Mask( 4)+ &
                         &            tmp2_r(i,ll( 5))*Mask( 5)+tmp2_r(i,ll( 6))*Mask( 6)+tmp2_r(i,ll( 7))*Mask( 7)+tmp2_r(i,ll( 8))*Mask( 8)+ &
                         &            tmp2_r(i,ll( 9))*Mask( 9)+tmp2_r(i,ll(10))*Mask(10)+tmp2_r(i,ll(11))*Mask(11)+tmp2_r(i,ll(12))*Mask(12)+ &
                         &            tmp2_r(i,ll(13))*Mask(13)+tmp2_r(i,ll(14))*Mask(14)+tmp2_r(i,ll(15))*Mask(15)+tmp2_r(i,ll(16))*Mask(16)+ &
                         &            tmp2_r(i,ll(17))*Mask(17)+tmp2_r(i,ll(18))*Mask(18)+tmp2_r(i,ll(19))*Mask(19)+tmp2_r(i,ll(20))*Mask(20)+ &
                         &            tmp2_r(i,ll(21))*Mask(21)+tmp2_r(i,ll(22))*Mask(22)+tmp2_r(i,ll(23))*Mask(23)+tmp2_r(i,ll(24))*Mask(24)+ &
                         &            tmp2_r(i,ll(25))*Mask(25)+tmp2_r(i,ll(26))*Mask(26)+tmp2_r(i,ll(27))*Mask(27)+tmp2_r(i,ll(28))*Mask(28)+ &
                         &            tmp2_r(i,ll(29))*Mask(29)+tmp2_r(i,ll(30))*Mask(30)+tmp2_r(i,ll(31))*Mask(31)+tmp2_r(i,ll(32))*Mask(32)+ &
                         &            tmp2_r(i,ll(33))*Mask(33)+tmp2_r(i,ll(34))*Mask(34)+tmp2_r(i,ll(35))*Mask(35)+tmp2_r(i,ll(36))*Mask(36)+ &
                         &            tmp2_r(i,ll(37))*Mask(37)+tmp2_r(i,ll(38))*Mask(38)+tmp2_r(i,ll(39))*Mask(39)+tmp2_r(i,ll(40))*Mask(40)+ &
                         &            tmp2_r(i,ll(41))*Mask(41)
                      END FORALL
                      ll=CSHIFT(ll,SHIFT=1,DIM=1)
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (21)
                ALLOCATE(tmp2_r(nnodes(1),43),ll(43),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                DO j=1,nnodes(2)
                   DO k=-20,21
                      l=k+21
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,k)
                      END FORALL
                   ENDDO
                   FORALL(i=1:43) ll(i)=i
                   DO k=1,nnodes(3)
                      kp=k+21
                      l=ll(43)
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,kp)
                      END FORALL
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp2_r(i,ll( 1))*Mask( 1)+tmp2_r(i,ll( 2))*Mask( 2)+tmp2_r(i,ll( 3))*Mask( 3)+tmp2_r(i,ll( 4))*Mask( 4)+ &
                         &            tmp2_r(i,ll( 5))*Mask( 5)+tmp2_r(i,ll( 6))*Mask( 6)+tmp2_r(i,ll( 7))*Mask( 7)+tmp2_r(i,ll( 8))*Mask( 8)+ &
                         &            tmp2_r(i,ll( 9))*Mask( 9)+tmp2_r(i,ll(10))*Mask(10)+tmp2_r(i,ll(11))*Mask(11)+tmp2_r(i,ll(12))*Mask(12)+ &
                         &            tmp2_r(i,ll(13))*Mask(13)+tmp2_r(i,ll(14))*Mask(14)+tmp2_r(i,ll(15))*Mask(15)+tmp2_r(i,ll(16))*Mask(16)+ &
                         &            tmp2_r(i,ll(17))*Mask(17)+tmp2_r(i,ll(18))*Mask(18)+tmp2_r(i,ll(19))*Mask(19)+tmp2_r(i,ll(20))*Mask(20)+ &
                         &            tmp2_r(i,ll(21))*Mask(21)+tmp2_r(i,ll(22))*Mask(22)+tmp2_r(i,ll(23))*Mask(23)+tmp2_r(i,ll(24))*Mask(24)+ &
                         &            tmp2_r(i,ll(25))*Mask(25)+tmp2_r(i,ll(26))*Mask(26)+tmp2_r(i,ll(27))*Mask(27)+tmp2_r(i,ll(28))*Mask(28)+ &
                         &            tmp2_r(i,ll(29))*Mask(29)+tmp2_r(i,ll(30))*Mask(30)+tmp2_r(i,ll(31))*Mask(31)+tmp2_r(i,ll(32))*Mask(32)+ &
                         &            tmp2_r(i,ll(33))*Mask(33)+tmp2_r(i,ll(34))*Mask(34)+tmp2_r(i,ll(35))*Mask(35)+tmp2_r(i,ll(36))*Mask(36)+ &
                         &            tmp2_r(i,ll(37))*Mask(37)+tmp2_r(i,ll(38))*Mask(38)+tmp2_r(i,ll(39))*Mask(39)+tmp2_r(i,ll(40))*Mask(40)+ &
                         &            tmp2_r(i,ll(41))*Mask(41)+tmp2_r(i,ll(42))*Mask(42)+tmp2_r(i,ll(43))*Mask(43)
                      END FORALL
                      ll=CSHIFT(ll,SHIFT=1,DIM=1)
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (22)
                ALLOCATE(tmp2_r(nnodes(1),45),ll(45),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                DO j=1,nnodes(2)
                   DO k=-21,22
                      l=k+22
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,k)
                      END FORALL
                   ENDDO
                   FORALL(i=1:45) ll(i)=i
                   DO k=1,nnodes(3)
                      kp=k+22
                      l=ll(45)
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,kp)
                      END FORALL
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp2_r(i,ll( 1))*Mask( 1)+tmp2_r(i,ll( 2))*Mask( 2)+tmp2_r(i,ll( 3))*Mask( 3)+tmp2_r(i,ll( 4))*Mask( 4)+ &
                         &            tmp2_r(i,ll( 5))*Mask( 5)+tmp2_r(i,ll( 6))*Mask( 6)+tmp2_r(i,ll( 7))*Mask( 7)+tmp2_r(i,ll( 8))*Mask( 8)+ &
                         &            tmp2_r(i,ll( 9))*Mask( 9)+tmp2_r(i,ll(10))*Mask(10)+tmp2_r(i,ll(11))*Mask(11)+tmp2_r(i,ll(12))*Mask(12)+ &
                         &            tmp2_r(i,ll(13))*Mask(13)+tmp2_r(i,ll(14))*Mask(14)+tmp2_r(i,ll(15))*Mask(15)+tmp2_r(i,ll(16))*Mask(16)+ &
                         &            tmp2_r(i,ll(17))*Mask(17)+tmp2_r(i,ll(18))*Mask(18)+tmp2_r(i,ll(19))*Mask(19)+tmp2_r(i,ll(20))*Mask(20)+ &
                         &            tmp2_r(i,ll(21))*Mask(21)+tmp2_r(i,ll(22))*Mask(22)+tmp2_r(i,ll(23))*Mask(23)+tmp2_r(i,ll(24))*Mask(24)+ &
                         &            tmp2_r(i,ll(25))*Mask(25)+tmp2_r(i,ll(26))*Mask(26)+tmp2_r(i,ll(27))*Mask(27)+tmp2_r(i,ll(28))*Mask(28)+ &
                         &            tmp2_r(i,ll(29))*Mask(29)+tmp2_r(i,ll(30))*Mask(30)+tmp2_r(i,ll(31))*Mask(31)+tmp2_r(i,ll(32))*Mask(32)+ &
                         &            tmp2_r(i,ll(33))*Mask(33)+tmp2_r(i,ll(34))*Mask(34)+tmp2_r(i,ll(35))*Mask(35)+tmp2_r(i,ll(36))*Mask(36)+ &
                         &            tmp2_r(i,ll(37))*Mask(37)+tmp2_r(i,ll(38))*Mask(38)+tmp2_r(i,ll(39))*Mask(39)+tmp2_r(i,ll(40))*Mask(40)+ &
                         &            tmp2_r(i,ll(41))*Mask(41)+tmp2_r(i,ll(42))*Mask(42)+tmp2_r(i,ll(43))*Mask(43)+tmp2_r(i,ll(44))*Mask(44)+ &
                         &            tmp2_r(i,ll(45))*Mask(45)
                      END FORALL
                      ll=CSHIFT(ll,SHIFT=1,DIM=1)
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (23)
                ALLOCATE(tmp2_r(nnodes(1),47),ll(47),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                DO j=1,nnodes(2)
                   DO k=-22,23
                      l=k+23
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,k)
                      END FORALL
                   ENDDO
                   FORALL(i=1:47) ll(i)=i
                   DO k=1,nnodes(3)
                      kp=k+23
                      l=ll(47)
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,kp)
                      END FORALL
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp2_r(i,ll( 1))*Mask( 1)+tmp2_r(i,ll( 2))*Mask( 2)+tmp2_r(i,ll( 3))*Mask( 3)+tmp2_r(i,ll( 4))*Mask( 4)+ &
                         &            tmp2_r(i,ll( 5))*Mask( 5)+tmp2_r(i,ll( 6))*Mask( 6)+tmp2_r(i,ll( 7))*Mask( 7)+tmp2_r(i,ll( 8))*Mask( 8)+ &
                         &            tmp2_r(i,ll( 9))*Mask( 9)+tmp2_r(i,ll(10))*Mask(10)+tmp2_r(i,ll(11))*Mask(11)+tmp2_r(i,ll(12))*Mask(12)+ &
                         &            tmp2_r(i,ll(13))*Mask(13)+tmp2_r(i,ll(14))*Mask(14)+tmp2_r(i,ll(15))*Mask(15)+tmp2_r(i,ll(16))*Mask(16)+ &
                         &            tmp2_r(i,ll(17))*Mask(17)+tmp2_r(i,ll(18))*Mask(18)+tmp2_r(i,ll(19))*Mask(19)+tmp2_r(i,ll(20))*Mask(20)+ &
                         &            tmp2_r(i,ll(21))*Mask(21)+tmp2_r(i,ll(22))*Mask(22)+tmp2_r(i,ll(23))*Mask(23)+tmp2_r(i,ll(24))*Mask(24)+ &
                         &            tmp2_r(i,ll(25))*Mask(25)+tmp2_r(i,ll(26))*Mask(26)+tmp2_r(i,ll(27))*Mask(27)+tmp2_r(i,ll(28))*Mask(28)+ &
                         &            tmp2_r(i,ll(29))*Mask(29)+tmp2_r(i,ll(30))*Mask(30)+tmp2_r(i,ll(31))*Mask(31)+tmp2_r(i,ll(32))*Mask(32)+ &
                         &            tmp2_r(i,ll(33))*Mask(33)+tmp2_r(i,ll(34))*Mask(34)+tmp2_r(i,ll(35))*Mask(35)+tmp2_r(i,ll(36))*Mask(36)+ &
                         &            tmp2_r(i,ll(37))*Mask(37)+tmp2_r(i,ll(38))*Mask(38)+tmp2_r(i,ll(39))*Mask(39)+tmp2_r(i,ll(40))*Mask(40)+ &
                         &            tmp2_r(i,ll(41))*Mask(41)+tmp2_r(i,ll(42))*Mask(42)+tmp2_r(i,ll(43))*Mask(43)+tmp2_r(i,ll(44))*Mask(44)+ &
                         &            tmp2_r(i,ll(45))*Mask(45)+tmp2_r(i,ll(46))*Mask(46)+tmp2_r(i,ll(47))*Mask(47)
                      END FORALL
                      ll=CSHIFT(ll,SHIFT=1,DIM=1)
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE (24)
                ALLOCATE(tmp2_r(nnodes(1),49),ll(49),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                DO j=1,nnodes(2)
                   DO k=-23,24
                      l=k+24
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,k)
                      END FORALL
                   ENDDO
                   FORALL(i=1:49) ll(i)=i
                   DO k=1,nnodes(3)
                      kp=k+24
                      l=ll(49)
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,kp)
                      END FORALL
                      FORALL (i=1:nnodes(1))
                         wpout(i,j,k)=tmp2_r(i,ll( 1))*Mask( 1)+tmp2_r(i,ll( 2))*Mask( 2)+tmp2_r(i,ll( 3))*Mask( 3)+tmp2_r(i,ll( 4))*Mask( 4)+ &
                         &            tmp2_r(i,ll( 5))*Mask( 5)+tmp2_r(i,ll( 6))*Mask( 6)+tmp2_r(i,ll( 7))*Mask( 7)+tmp2_r(i,ll( 8))*Mask( 8)+ &
                         &            tmp2_r(i,ll( 9))*Mask( 9)+tmp2_r(i,ll(10))*Mask(10)+tmp2_r(i,ll(11))*Mask(11)+tmp2_r(i,ll(12))*Mask(12)+ &
                         &            tmp2_r(i,ll(13))*Mask(13)+tmp2_r(i,ll(14))*Mask(14)+tmp2_r(i,ll(15))*Mask(15)+tmp2_r(i,ll(16))*Mask(16)+ &
                         &            tmp2_r(i,ll(17))*Mask(17)+tmp2_r(i,ll(18))*Mask(18)+tmp2_r(i,ll(19))*Mask(19)+tmp2_r(i,ll(20))*Mask(20)+ &
                         &            tmp2_r(i,ll(21))*Mask(21)+tmp2_r(i,ll(22))*Mask(22)+tmp2_r(i,ll(23))*Mask(23)+tmp2_r(i,ll(24))*Mask(24)+ &
                         &            tmp2_r(i,ll(25))*Mask(25)+tmp2_r(i,ll(26))*Mask(26)+tmp2_r(i,ll(27))*Mask(27)+tmp2_r(i,ll(28))*Mask(28)+ &
                         &            tmp2_r(i,ll(29))*Mask(29)+tmp2_r(i,ll(30))*Mask(30)+tmp2_r(i,ll(31))*Mask(31)+tmp2_r(i,ll(32))*Mask(32)+ &
                         &            tmp2_r(i,ll(33))*Mask(33)+tmp2_r(i,ll(34))*Mask(34)+tmp2_r(i,ll(35))*Mask(35)+tmp2_r(i,ll(36))*Mask(36)+ &
                         &            tmp2_r(i,ll(37))*Mask(37)+tmp2_r(i,ll(38))*Mask(38)+tmp2_r(i,ll(39))*Mask(39)+tmp2_r(i,ll(40))*Mask(40)+ &
                         &            tmp2_r(i,ll(41))*Mask(41)+tmp2_r(i,ll(42))*Mask(42)+tmp2_r(i,ll(43))*Mask(43)+tmp2_r(i,ll(44))*Mask(44)+ &
                         &            tmp2_r(i,ll(45))*Mask(45)+tmp2_r(i,ll(46))*Mask(46)+tmp2_r(i,ll(47))*Mask(47)+tmp2_r(i,ll(48))*Mask(48)+ &
                         &            tmp2_r(i,ll(49))*Mask(49)
                      END FORALL
                      ll=CSHIFT(ll,SHIFT=1,DIM=1)
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             CASE DEFAULT
                m=MaskSize*2+1

                ALLOCATE(tmp2_r(nnodes(1),m),ll(m),STAT=info)
                or_fail_alloc("tmp2_r",exit_point=7777)

                DO j=1,nnodes(2)
                   DO k=1-MaskSize,MaskSize
                      l=k+MaskSize
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,k)
                      END FORALL
                   ENDDO
                   FORALL(i=1:m) ll(i)=i
                   DO k=1,nnodes(3)
                      kp=k+MaskSize
                      l=ll(m)
                      FORALL (i=1:nnodes(1))
                         tmp2_r(i,l)=wpin(i,j,kp)
                      END FORALL
                      DO i=1,nnodes(1)
                         wpout(i,j,k)=tmp2_r(i,l)*Mask(m)
                         DO n=1,m-1
                            l1=ll(n)
                            wpout(i,j,k)=wpout(i,j,k)+tmp2_r(i,l1)*Mask(n)
                         ENDDO
                      ENDDO
                      ll=CSHIFT(ll,SHIFT=1,DIM=1)
                   ENDDO !j=1,nnodes(2)
                ENDDO !k=1,nnodes(3)
             END SELECT

             DEALLOCATE(tmp2_r,ll,STAT=info)
             or_fail_dealloc("tmp2",exit_point=7777)

          CASE DEFAULT
             fail("Wrong direction for 3D problem!!!",exit_point=7777)
          END SELECT

        7777 CONTINUE
        END SUBROUTINE Update
#endif
        END SUBROUTINE DTYPE(ppm_rc_GaussianImageFilter)

