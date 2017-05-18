      !-------------------------------------------------------------------------
      !  Program      :                    ppm_rc_EdgeDetection
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
      !  Subroutine   :                    ppm_rc_EdgeDetection
      !-------------------------------------------------------------------------
      !
      !  Purpose      :  A 2D and 3D edge detection using the Sobel operator on a
      !                  blured image convolved with discrete gaussian kernels.
      !
      !                  We use the Sobel operator to calculate the image gradient and then
      !                  finds the magnitude of this gradient vector.
      !                  The Sobel gradient magnitude (square-root sum of squares) is
      !                  an indication of edge strength.
      !
      !  Input        :
      !
      !  Input/output :
      !
      !  Output       : info          (I) return status. 0 on success.
      !
      !  Routines     :
      !
      !  Remarks      : The original image should have at least 4 ghost layers
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      SUBROUTINE DTYPE(ppm_rc_EdgeDetection)(FieldIn,MeshIn,info,MCMCZe_)

        USE ppm_module_interfaces, ONLY : ppm_t_equi_mesh_,ppm_t_field_
        USE ppm_module_field_typedef, ONLY : ppm_t_field
        USE ppm_module_mpi

        USE ppm_rc_module_write, ONLY : DTYPE(ppm_rc_write_image)
        USE ppm_rc_module_rnd, ONLY : ppm_rc_GenerateImageDiscrDistr
        USE ppm_rc_module_util, ONLY : DTYPE(ppm_rc_CopyImageAndNormalize)
        IMPLICIT NONE

        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        CLASS(ppm_t_field_),             POINTER       :: FieldIn
        !!! the source filed

        CLASS(ppm_t_equi_mesh_),         POINTER       :: MeshIn
        !!! The FieldIn is discretized on this Mesh

        INTEGER,                         INTENT(  OUT) :: info
        !!! Status

        REAL(ppm_kind_double), OPTIONAL, INTENT(INOUT) :: MCMCZe_
        !!!
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
        CLASS(ppm_t_field_), POINTER :: Fieldtmp
        !!! the source filed

        CLASS(ppm_t_subpatch_), POINTER :: sbpitr

        REAL(ppm_kind_double)                               :: t0
#if   __DIME == __2D
        REAL(MK), CONTIGUOUS, DIMENSION(:,:),   POINTER     :: DTYPE(wpi)
#elif __DIME == __3D
        REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER     :: DTYPE(wpi)
#endif
        REAL(ppm_kind_double)                               :: MCMCZe
        REAL(MK)                                            :: FMinVal,FMaxVal
        REAL(MK)                                            :: Normalfac

#if   __DIME == __2D
        INTEGER, CONTIGUOUS, DIMENSION(:,:),   POINTER :: DTYPE(wpl)
#elif __DIME == __3D
        INTEGER, CONTIGUOUS, DIMENSION(:,:,:), POINTER :: DTYPE(wpl)
#endif
        INTEGER,             DIMENSION(:),     POINTER :: Nm
        INTEGER                                        :: i,j
#if   __DIME == __3D
        INTEGER                                        :: k
#endif
#ifdef __Linux
        INTEGER                                        :: memory
#endif

        LOGICAL :: tmpNormalize

        CHARACTER(LEN=ppm_char) :: caller = 'ppm_rc_EdgeDetection'
        !-------------------------------------------------------------------------
        !  Initialize
        !-------------------------------------------------------------------------
        CALL substart(caller,t0,info)

        CALL check

        MCMCZe=MERGE(MCMCZe_,zerod,PRESENT(MCMCZe_))

        ALLOCATE(ppm_t_field::Fieldtmp,STAT=info)
        or_fail_alloc('Failed to allocate Fieldtmp.')

        CALL Fieldtmp%create(1,info,dtype=FieldIn%data_type,name="edge_image")
        or_fail("Failed to create field!" )

        SELECT CASE (lNormalize)
        CASE (.TRUE.)
           Normalfac=ImageNormalfac

           ! First blur the image by a Gaussian kernel with sigma=SQRT(2) (!rull of thumb)
           ! then use the Sobel operator to calculate the image gradient
           CALL DTYPE(ppm_rc_GaussianImageFilter)((/SQRT(two)/),FieldIn,MeshIn,info,&
           &    Fieldtmp,UseImageSpacing=.FALSE.)
        CASE (.FALSE.)
           CALL Fieldtmp%discretize_on(MeshIn,info)
           or_fail("Failed to discretize Field on Mesh!")

           CALL DTYPE(ppm_rc_CopyImageAndNormalize)(FieldIn,Fieldtmp,MeshIn,info, &
           &    Normalize=.TRUE.,FieldMinVal=FMinVal,FieldMaxVal=FMaxVal)
           or_fail("Normalizing the FieldIn faield!")

           Normalfac=FMaxVal-FMinVal

           ! First blur the image by a Gaussian kernel with sigma=2 (!rull of thumb)
           ! then use the Sobel operator to calculate the image gradient
           CALL DTYPE(ppm_rc_GaussianImageFilter)((/SQRT(two)/),Fieldtmp,MeshIn,info,&
           &    UseImageSpacing=.FALSE.)
        END SELECT
        or_fail("ppm_rc_GaussianImageFilter!")

        IF (ppm_nproc.GT.1) THEN
           ghostsize=1
           !-------------------------------------------------------------------------
           !  Get field ghosts for image
           !-------------------------------------------------------------------------
           CALL MeshIn%map_ghost_get(info,ghostsize=ghostsize)
           or_fail("MeshIn%map_ghost_get")
           !-------------------------------------------------------------------------
           !  fill the buffer
           !-------------------------------------------------------------------------
           CALL Fieldtmp%map_ghost_push(MeshIn,info)
           or_fail("Fieldtmp%map_ghost_push")
           !-------------------------------------------------------------------------
           !  Non-blocking send & recieve
           !-------------------------------------------------------------------------
           CALL MeshIn%map_isend(info,sendrecv=.TRUE.)
           or_fail("MeshIn%map_isend")
        ENDIF

        !-------------------------------------------------------------------------
        !  Pads the input at the boarder of the global domain
        !  pad_size: 1
        !-------------------------------------------------------------------------
        sbpitr => MeshIn%subpatch%begin()
        DO WHILE (ASSOCIATED(sbpitr))
           Nm => sbpitr%nnodes

           CALL sbpitr%get_field(Fieldtmp,DTYPE(wpi),info)
           or_fail("Failed to get field wpi data.")

#if   __DIME == __2D
           !padding the left xz image face
           IF (sbpitr%istart(2).EQ.1) THEN
              j=0
              DO i=0,Nm(1)+1
                 DTYPE(wpi)(i,j)=DTYPE(wpi)(i,1)
              ENDDO !i=0,Nm(1)+1
           ENDIF !(sbpitr%istart(2).EQ.1)
           !padding the right xz image face
           IF (sbpitr%iend(2).EQ.Ngrid(2)) THEN
              j=Nm(2)+1
              DO i=0,Nm(1)+1
                 DTYPE(wpi)(i,j)= DTYPE(wpi)(i,Nm(2))
              ENDDO !i=0,Nm(1)+1
           ENDIF !(sbpitr%iend(2).EQ.Ngrid(2))
           !padding the left yz image face
           IF (sbpitr%istart(1).EQ.1) THEN
              i=0
              DO j=0,Nm(2)+1
                 DTYPE(wpi)(i,j)=DTYPE(wpi)(1,j)
              ENDDO !j=0,Nm(2)+1
           ENDIF !(sbpitr%istart(1).EQ.1)
           !padding the right yz image face
           IF (sbpitr%iend(1).EQ.Ngrid(1)) THEN
              i=Nm(1)+1
              DO j=0,Nm(2)+1
                 DTYPE(wpi)(i,j)=DTYPE(wpi)(Nm(1),j)
              ENDDO !j=0,Nm(2)+1
           ENDIF !(sbpitr%iend(1).EQ.Ngrid(1))
#elif __DIME == __3D
           !padding the bottom xy image face
           IF (sbpitr%istart(3).EQ.1) THEN
              k=0
              DO j=0,Nm(2)+1
                 DO i=0,Nm(1)+1
                    DTYPE(wpi)(i,j,k)=DTYPE(wpi)(i,j,1)
                 ENDDO !i=0,Nm(1)+1
              ENDDO !j=0,Nm(2)+1
           ENDIF !(sbpitr%istart(3).EQ.1)
           !padding the top xy image face
           IF (sbpitr%iend(3).EQ.Ngrid(3)) THEN
              k=Nm(3)+1
              DO j=0,Nm(2)+1
                 DO i=0,Nm(1)+1
                    DTYPE(wpi)(i,j,k)= DTYPE(wpi)(i,j,Nm(3))
                 ENDDO !i=0,Nm(1)+1
              ENDDO !j=0,Nm(2)+1
           ENDIF !(sbpitr%iend(3).EQ.Ngrid(3))

           !padding the left xz image face
           IF (sbpitr%istart(2).EQ.1) THEN
              j=0
              DO k=0,Nm(3)+1
                 DO i=0,Nm(1)+1
                    DTYPE(wpi)(i,j,k)=DTYPE(wpi)(i,1,k)
                 ENDDO !i=0,Nm(1)+1
              ENDDO !k=0,Nm(3)+1
           ENDIF !(sbpitr%istart(2).EQ.1)
           !padding the right xz image face
           IF (sbpitr%iend(2).EQ.Ngrid(2)) THEN
              j=Nm(2)+1
              DO k=0,Nm(3)+1
                 DO i=0,Nm(1)+1
                    DTYPE(wpi)(i,j,k)= DTYPE(wpi)(i,Nm(2),k)
                 ENDDO !i=0,Nm(1)+1
              ENDDO !k=0,Nm(3)+1
           ENDIF !(sbpitr%iend(2).EQ.Ngrid(2))

           !padding the left yz image face
           IF (sbpitr%istart(1).EQ.1) THEN
              i=0
              DO k=0,Nm(3)+1
                 DO j=0,Nm(2)+1
                    DTYPE(wpi)(i,j,k)= DTYPE(wpi)(1,j,k)
                 ENDDO !j=0,Nm(2)+1
              ENDDO !k=0,Nm(3)+1
           ENDIF !(sbpitr%istart(1).EQ.1)
           !padding the right yz image face
           IF (sbpitr%iend(1).EQ.Ngrid(1)) THEN
              i=Nm(1)+1
              DO k=0,Nm(3)+1
                 DO j=0,Nm(2)+1
                    DTYPE(wpi)(i,j,k)= DTYPE(wpi)(Nm(1),j,k)
                 ENDDO !j=0,Nm(2)+1
              ENDDO !k=0,Nm(3)+1
           ENDIF !(sbpitr%iend(1).EQ.Ngrid(1))
#endif
           sbpitr => MeshIn%subpatch%next()
        ENDDO
        NULLIFY(DTYPE(wpi))

        !!! allocate the helper field for computing gradient in different direction
        !!! Here I take advantage of non-blocking communication to allocate and
        !!! discretize the filed
        ALLOCATE(ppm_t_field::fld1,STAT=info)
        or_fail_alloc('Failed to allocate fld1.')

        CALL fld1%create(1,info,dtype=FieldIn%data_type,name="Ygradient")
        or_fail("Failed to create field!" )

        CALL fld1%discretize_on(MeshIn,info)
        or_fail("Failed to discretize Field on Mesh!")

#if   __DIME == __3D
        ALLOCATE(ppm_t_field::fld2,STAT=info)
        or_fail_alloc('Failed to allocate fld2.')

        CALL fld2%create(1,info,dtype=FieldIn%data_type,name="Zgradient")
        or_fail("Failed to create field!" )

        CALL fld2%discretize_on(MeshIn,info)
        or_fail("Failed to discretize Field on Mesh!")
#endif

        IF (ppm_nproc.GT.1) THEN
           !-------------------------------------------------------------------------
           !  wait to complete the recieve
           !-------------------------------------------------------------------------
           CALL MeshIn%map_isend(info,sendrecv=.FALSE.)
           or_fail("MeshIn%map_isend")
           !-------------------------------------------------------------------------
           !  empty the buffer
           !-------------------------------------------------------------------------
           CALL Fieldtmp%map_ghost_pop(MeshIn,info)
           or_fail("Fieldtmp%map_ghost_pop")
        ENDIF

#if   __DIME == __2D
        ! The Sobel operator in 2D
        ! Sx  ----  Gx=Sx*A
        ! Sy  ----  Gy=Sy*A
        ! the gradient magnitud G=sqrt(Gx^2+Gy^2)
#elif __DIME == __3D
        ! The Sobel operator in 3D
        ! Sx  ----  Gx=Sx*A
        ! Sy  ----  Gy=Sy*A
        ! Sz  ----  Gz=Sz*A
        ! the gradient magnitud G=sqrt(Gx^2+Gy^2+Gz^2)
#endif
        !-------------------------------------------------------------------------
        !  Use Sobel operator to calculate the image gradient and then finds the
        !  magnitude of this gradient vector. The Sobel gradient magnitude
        !  (square-root sum of squares) is an indication of edge strength.
        !-------------------------------------------------------------------------
        CALL DTYPE(ppm_rc_SobelImageFilter)(Fieldtmp,MeshIn,info)
        or_fail("ppm_rc_SobelImageFilter!")

        sbpitr => MeshIn%subpatch%begin()
        DO WHILE (ASSOCIATED(sbpitr))
           Nm => sbpitr%nnodes

           CALL sbpitr%get_field(Fieldtmp,DTYPE(wpi),info)
           or_fail("Failed to get field wpi data.")

#if   __DIME == __2D
           DO j=1,Nm(2)
              DO i=1,Nm(1)
                 MCMCZe=MCMCZe+REAL(DTYPE(wpi)(i,j),ppm_kind_double)
              ENDDO
           ENDDO
#elif __DIME == __3D
           DO k=1,Nm(3)
              DO j=1,Nm(2)
                 DO i=1,Nm(1)
                    MCMCZe=MCMCZe+REAL(DTYPE(wpi)(i,j,k),ppm_kind_double)
                 ENDDO
              ENDDO
           ENDDO
#endif
           sbpitr => MeshIn%subpatch%next()
        ENDDO
        NULLIFY(DTYPE(wpi))

#ifdef __MPI
        CALL MPI_Allreduce(MPI_IN_PLACE,MCMCZe,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,info)
        or_fail_MPI("MPI_Allreduce")
#endif

        sbpitr => MeshIn%subpatch%begin()
        DO WHILE (ASSOCIATED(sbpitr))
           Nm => sbpitr%nnodes

           CALL sbpitr%get_field(Fieldtmp,DTYPE(wpi),info)
           or_fail("Failed to get field wpi data.")

           IF (MCMCZe.LE.smalld*tend) THEN
#if   __DIME == __2D
              DO j=1,Nm(2)
                 DO i=1,Nm(1)
                    DTYPE(wpi)(i,j)=one
                 ENDDO
              ENDDO
              MCMCZe=REAL(Nm(1),ppm_kind_double)*REAL(Nm(2),ppm_kind_double)
#elif __DIME == __3D
              DO k=1,Nm(3)
                 DO j=1,Nm(2)
                    DO i=1,Nm(1)
                       DTYPE(wpi)(i,j,k)=one
                    ENDDO
                 ENDDO
              ENDDO
              MCMCZe=REAL(Nm(1),ppm_kind_double)*REAL(Nm(2),ppm_kind_double)*REAL(Nm(3),ppm_kind_double)
#endif
           ENDIF
           !-------------------------------------------------------------------------
           !  Constructs a discrete_distribution from an edge image.
           !  The values of the range represent weights for the possible values of the distribution.
           !-------------------------------------------------------------------------
#if   __DIME == __2D
           info=ppm_rc_GenerateImageDiscrDistr(DTYPE(wpi),Nm(1),Nm(2))
#elif __DIME == __3D
           info=ppm_rc_GenerateImageDiscrDistr(DTYPE(wpi),Nm(1),Nm(2),Nm(3))
#endif
           or_fail("ppm_rc_GenerateImageDiscrDistr")

           sbpitr => MeshIn%subpatch%next()
        ENDDO
        NULLIFY(DTYPE(wpi))

        IF (debug.GT.0) THEN
           ! We are shifting and scaling the Fieldtmp as it has negative value.
           CALL DTYPE(ppm_rc_CopyImageAndNormalize)(Fieldtmp,Fieldtmp,MeshIn,info, &
           &    Normalize=.TRUE.,FieldMinVal=FMinVal,FieldMaxVal=FMaxVal)
           or_fail("ppm_rc_CopyImageAndNormalize")


           IF (Normalfac.LE.255._MK) THEN
              Normalfac=255._MK
           ELSE
              Normalfac=65535._MK
           ENDIF

           tmpNormalize=lNormalize

           lNormalize=.TRUE.

           ! Fieldtmp is already shifted and normalized, so we Write a 8 BITS Edge image
           CALL DTYPE(ppm_rc_write_image)(Fieldtmp,MeshIn,"EdgeImage",info,Scalefac=Normalfac)
           or_fail("ppm_rc_write_image")

           lNormalize=tmpNormalize
        ENDIF

        !-------------------------------------------------------------------------
        !  Free memory
        !-------------------------------------------------------------------------
        CALL Fieldtmp%destroy(info)
        or_fail("Destroying field failed!")

        ! deallocate the pointer:
        dealloc_pointer("Fieldtmp")

        IF (PRESENT(MCMCZe_)) MCMCZe_=MCMCZe
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
      9999 CONTINUE
        CALL substop(caller,t0,info)
        RETURN
      CONTAINS
      SUBROUTINE check
      IMPLICIT NONE
        IF (.NOT.MeshIn%subpatch%nb.EQ.1) THEN
           fail("This algorithm is designed and optimized for one subdomain per node!!!", &
           & ppm_error=ppm_error_fatal,exit_point=8888)
        ENDIF
        IF (ANY(MeshIn%ghostsize.LT.4)) THEN
           fail("The current implementation is heavily depends on the FieldIn ghost size which should not be less than 4 pixels on each side!", &
           & ppm_error=ppm_error_fatal,exit_point=8888)
        ENDIF
      8888 CONTINUE
      END SUBROUTINE check
      END SUBROUTINE DTYPE(ppm_rc_EdgeDetection)

