        !-------------------------------------------------------------------------
        !  Subroutine   :                    ppm_rc_SobelImageFilter
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
        !  Author           - y.afshar           Dec       2015
        !-------------------------------------------------------------------------

        !-------------------------------------------------------------------------
        !  Subroutine   :                    ppm_rc_SobelImageFilter
        !-------------------------------------------------------------------------
        !
        !  Purpose      : Filter an image by a Sobel filter.
        !
        !
        !  Input        :
        !
        !  Input/output : FieldIn, MeshIn, FieldOut
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
        SUBROUTINE DTYPE(ppm_rc_SobelImageFilter)(FieldIn,MeshIn,info,FieldOut)

          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------

          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(ppm_t_field_),           POINTER       :: FieldIn
          !!! the source filed which is supposed to be

          CLASS(ppm_t_equi_mesh_),       POINTER       :: MeshIn
          !!! source mesh on which fields are descritized

          INTEGER,                       INTENT(  OUT) :: info

          CLASS(ppm_t_field_), OPTIONAL, POINTER       :: FieldOut
          !!! the output filed which is supposed to be

          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          CLASS(ppm_t_subpatch_), POINTER :: sbpitr

          REAL(ppm_kind_double)                               :: t0
#if   __DIME == __2D
          REAL(MK), CONTIGUOUS, DIMENSION(:,:),   POINTER     :: DTYPE(wpin)
          REAL(MK), CONTIGUOUS, DIMENSION(:,:),   POINTER     :: DTYPE(wpout)
          REAL(MK), CONTIGUOUS, DIMENSION(:,:),   POINTER     :: DTYPE(wpy)
#elif __DIME == __3D
          REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER     :: DTYPE(wpin)
          REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER     :: DTYPE(wpout)
          REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER     :: DTYPE(wpy)
          REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER     :: wpz
#endif

          REAL(MK),             DIMENSION(:,:),   ALLOCATABLE :: sob
          REAL(MK),             DIMENSION(:),     ALLOCATABLE :: MASK

          INTEGER, DIMENSION(:), POINTER :: Nm
          INTEGER                        :: i,j,jp,l,l1,l2,l3
#if   __DIME == __3D
          INTEGER                        :: k,m
#endif

          CHARACTER(LEN=ppm_char) :: caller = 'ppm_rc_SobelImageFilter'
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          CALL check()

#if   __DIME == __2D
          ! The Sobel operator in 2D
          ! Sx  ----  Gx=Sx*A
          ! -1  0  1
          ! -2  0  2
          ! -1  0  1
          !
          ! Sy  ----  Gy=Sy*A
          ! -1 -2 -1
          !  0  0  0
          !  1  2  1
          !
          ! the gradient magnitud G=sqrt(Gx^2+Gy^2)
#elif __DIME == __3D
          ! The Sobel operator in 3D
          ! Sx  ----  Gx=Sx*A
          ! -1  0  1    -3  0  3    -1  0  1
          ! -3  0  3    -6  0  6    -3  0  3
          ! -1  0  1    -3  0  3    -1  0  1
          !
          ! Sy  ----  Gy=Sy*A
          ! -1 -3 -1    -3 -6 -3    -1 -3 -1
          !  0  0  0     0  0  0     0  0  0
          !  1  3  1     3  6  3     1  3  1
          !
          ! Sz  ----  Gz=Sz*A
          ! -1 -3 -1     0  0  0     1  3  1
          ! -3 -6 -3     0  0  0     3  6  3
          ! -1 -3 -1     0  0  0     1  3  1
          !
          ! the gradient magnitud G=sqrt(Gx^2+Gy^2+Gz^2)
#endif

#if   __DIME == __2D
          NULLIFY(DTYPE(wpin),DTYPE(wpout),DTYPE(wpy))
#elif __DIME == __3D
          NULLIFY(DTYPE(wpin),DTYPE(wpout),DTYPE(wpy),wpz)
#endif
          sbpitr => MeshIn%subpatch%begin()
          DO WHILE (ASSOCIATED(sbpitr))
             Nm => sbpitr%nnodes

             CALL sbpitr%get_field(FieldIn,DTYPE(wpin),info)
             or_fail("Failed to get field data wpin.")

             IF (PRESENT(FieldOut)) THEN
                IF (ASSOCIATED(FieldIn,FieldOut)) THEN
                   DTYPE(wpout) => DTYPE(wpin)
                ELSE
                   CALL sbpitr%get_field(FieldOut,DTYPE(wpout),info)
                   or_fail("Failed to get field data wpout.")
                ENDIF
             ELSE
                DTYPE(wpout) => DTYPE(wpin)
             ENDIF

             CALL sbpitr%get_field(fld1,DTYPE(wpy),info)
             or_fail("Failed to get field data wpy.")

#if   __DIME == __3D
             CALL sbpitr%get_field(fld2,wpz,info)
             or_fail("Failed to get field data wpz.")
#endif

#if   __DIME == __2D
             !At each pixel, I should keep two values of the first
             !and third row mask multiply by the pixel row
             ALLOCATE(sob(Nm(1),6),MASK(6),STAT=info)
             or_fail_alloc("Failed to allocate sob & MASK!")
             ! The Sobel operator in 2D
             ! Sy  ----  Gy=Sy*A
             ! -1 -2 -1
             !  0  0  0
             !  1  2  1
             !  ---------------->
             !  1  2  1
             !  0  0  0
             ! -1 -2 -1
             MASK(1)= one
             MASK(2)= two
             MASK(3)= one
             MASK(4)=-one
             MASK(5)=-two
             MASK(6)=-one

             ! This optimisation is done to avoid recomputing what has been computed once
             j=0
             DO i=1,Nm(1)
!               sob(i,1)=DTYPE(wpin)(i-1,j)*MASK(1)+DTYPE(wpin)(i,j)*MASK(2)+DTYPE(wpin)(i+1,j)*MASK(3)
                sob(i,1)=DTYPE(wpin)(i-1,j)+DTYPE(wpin)(i,j)*MASK(2)+DTYPE(wpin)(i+1,j)
             ENDDO

             j=1
             DO i=1,Nm(1)
!               sob(i,3)=DTYPE(wpin)(i-1,j)*MASK(1)+DTYPE(wpin)(i,j)*MASK(2)+DTYPE(wpin)(i+1,j)*MASK(3)
                sob(i,3)=DTYPE(wpin)(i-1,j)+DTYPE(wpin)(i,j)*MASK(2)+DTYPE(wpin)(i+1,j)
             ENDDO

             DO j=1,Nm(2)
                jp=j+1
                l=MOD(jp*2,6)+1
                l2=l+1
                DO i=1,Nm(1)
!                  sob(i,l )=DTYPE(wpin)(i-1,jp)*MASK(1)+DTYPE(wpin)(i,jp)*MASK(2)+DTYPE(wpin)(i+1,jp)*MASK(3)
!                  sob(i,l2)=DTYPE(wpin)(i-1,jp)*MASK(4)+DTYPE(wpin)(i,jp)*MASK(5)+DTYPE(wpin)(i+1,jp)*MASK(6)
                   sob(i,l )= DTYPE(wpin)(i-1,jp)+DTYPE(wpin)(i,jp)*MASK(2)+DTYPE(wpin)(i+1,jp)
                   sob(i,l2)=-DTYPE(wpin)(i-1,jp)+DTYPE(wpin)(i,jp)*MASK(5)-DTYPE(wpin)(i+1,jp)
                ENDDO
                l1=MOD(l+2,6)

                FORALL (i=1:Nm(1))
                   DTYPE(wpy)(i,j)=sob(i,l1)
                END FORALL
                FORALL (i=1:Nm(1))
                   DTYPE(wpy)(i,j)=DTYPE(wpy)(i,j)+sob(i,l2)
                END FORALL
                !calculate the square Gy^2
                FORALL (i=1:Nm(1))
                   DTYPE(wpy)(i,j)=DTYPE(wpy)(i,j)*DTYPE(wpy)(i,j)
                END FORALL
             ENDDO !j=1,Nm(2)

             DEALLOCATE(sob,STAT=info)
             or_fail_dealloc("sob")

             ALLOCATE(sob(Nm(1),9),STAT=info)
             or_fail_alloc("sob & MASK")
             ! The Sobel operator in 2D
             ! Sx  ----  Gx=Sx*A
             ! -1  0  1
             ! -2  0  2
             ! -1  0  1
             !  ---------------->
             !  1  0 -1
             !  2  0 -2
             !  1  0 -1
             MASK(1)= one
             MASK(2)=-one
             MASK(3)= two
             MASK(4)=-two
             MASK(5)= one
             MASK(6)=-one
             ! This optimisation is done to avoid recomputing what has been computed once
             j=0
             DO i=1,Nm(1)
!               sob(i,1)=DTYPE(wpin)(i-1,j)*MASK(1)+DTYPE(wpin)(i+1,j)*MASK(2)
                sob(i,1)=DTYPE(wpin)(i-1,j)-DTYPE(wpin)(i+1,j)
             ENDDO

             j=1
             DO i=1,Nm(1)
!               sob(i,4)=DTYPE(wpin)(i-1,j)*MASK(1)+DTYPE(wpin)(i+1,j)*MASK(2)
!               sob(i,5)=DTYPE(wpin)(i-1,j)*MASK(3)+DTYPE(wpin)(i+1,j)*MASK(4)
                sob(i,4)=DTYPE(wpin)(i-1,j)        -DTYPE(wpin)(i+1,j)
                sob(i,5)=DTYPE(wpin)(i-1,j)*MASK(3)+DTYPE(wpin)(i+1,j)*MASK(4)
             ENDDO

             DO j=1,Nm(2)
                jp=j+1
                l=MOD(jp*3,9)+1
                l2=l+1
                l3=l+2
                DO i=1,Nm(1)
!                  sob(i,l )=DTYPE(wpin)(i-1,jp)*MASK(1)+DTYPE(wpin)(i+1,jp)*MASK(2)
!                  sob(i,l2)=DTYPE(wpin)(i-1,jp)*MASK(3)+DTYPE(wpin)(i+1,jp)*MASK(4)
!                  sob(i,l3)=DTYPE(wpin)(i-1,jp)*MASK(5)+DTYPE(wpin)(i+1,jp)*MASK(6)
                   sob(i,l )=DTYPE(wpin)(i-1,jp)        -DTYPE(wpin)(i+1,jp)
                   sob(i,l2)=DTYPE(wpin)(i-1,jp)*MASK(3)+DTYPE(wpin)(i+1,jp)*MASK(4)
                   sob(i,l3)=DTYPE(wpin)(i-1,jp)        -DTYPE(wpin)(i+1,jp)
                ENDDO
                l1=MOD(l+3,9)
                l2=MOD(l+7,9)
                FORALL (i=1:Nm(1))
                   DTYPE(wpout)(i,j)=sob(i,l1)
                END FORALL
                FORALL (i=1:Nm(1))
                   DTYPE(wpout)(i,j)=DTYPE(wpout)(i,j)+sob(i,l2)
                END FORALL
                FORALL (i=1:Nm(1))
                   DTYPE(wpout)(i,j)=DTYPE(wpout)(i,j)+sob(i,l3)
                END FORALL
                !calculate the square Gx^2
                FORALL (i=1:Nm(1))
                   DTYPE(wpout)(i,j)=DTYPE(wpout)(i,j)*DTYPE(wpout)(i,j)
                END FORALL
             ENDDO !j=1,Nm(2)

             DEALLOCATE(sob,MASK,STAT=info)
             or_fail_dealloc("sob & MASK")

             ! the gradient magnitud G=sqrt(Gx^2+Gy^2)
             FORALL (i=1:Nm(1),j=1:Nm(2))
                DTYPE(wpout)(i,j)=SQRT(DTYPE(wpout)(i,j)+DTYPE(wpy)(i,j))
             END FORALL
#elif __DIME == __3D
             !At each pixel, I should keep two values of the first
             !and third row mask multiply by the pixel row
             ALLOCATE(sob(Nm(1),9),MASK(18),STAT=info)
             or_fail_alloc("sob & MASK")
             ! The Sobel operator in 3D
             ! Sz  ----  Gz=Sz*A
             ! -1 -3 -1     0  0  0     1  3  1
             ! -3 -6 -3     0  0  0     3  6  3
             ! -1 -3 -1     0  0  0     1  3  1
             !  ---------------->
             !  1  3  1     0  0  0    -1 -3 -1
             !  3  6  3     0  0  0    -3 -6 -3
             !  1  3  1     0  0  0    -1 -3 -1
             MASK( 1)= one
             MASK( 2)= three
             MASK( 3)= one
             MASK( 4)=-one
             MASK( 5)=-three
             MASK( 6)=-one

             MASK( 7)= three
             MASK( 8)= 6.0_MK
             MASK( 9)= three
             MASK(10)=-three
             MASK(11)=-6.0_MK
             MASK(12)=-three

             MASK(13)= one
             MASK(14)= three
             MASK(15)= one
             MASK(16)=-one
             MASK(17)=-three
             MASK(18)=-one

             DO k=1,Nm(3)
                m=k-1
                j=0
                DO i=1,Nm(1)
!                  sob(i,1)=         DTYPE(wpin)(i-1,j,m)*MASK(1)+DTYPE(wpin)(i,j,m)*MASK(2)+DTYPE(wpin)(i+1,j,m)*MASK(3)
                   sob(i,1)=         DTYPE(wpin)(i-1,j,m)+DTYPE(wpin)(i,j,m)*MASK(2)+DTYPE(wpin)(i+1,j,m)
                ENDDO
                j=1
                DO i=1,Nm(1)
!                  sob(i,4)=         DTYPE(wpin)(i-1,j,m)*MASK(1)+DTYPE(wpin)(i,j,m)*MASK(2)+DTYPE(wpin)(i+1,j,m)*MASK(3)
!                  sob(i,5)=         DTYPE(wpin)(i-1,j,m)*MASK(7)+DTYPE(wpin)(i,j,m)*MASK(8)+DTYPE(wpin)(i+1,j,m)*MASK(9)
                   sob(i,4)=         DTYPE(wpin)(i-1,j,m)        +DTYPE(wpin)(i,j,m)*MASK(2)+DTYPE(wpin)(i+1,j,m)
                   sob(i,5)=         DTYPE(wpin)(i-1,j,m)*MASK(7)+DTYPE(wpin)(i,j,m)*MASK(8)+DTYPE(wpin)(i+1,j,m)*MASK(9)
                ENDDO

                m=k+1
                j=0
                DO i=1,Nm(1)
!                  sob(i,1)=sob(i,1)+DTYPE(wpin)(i-1,j,m)*MASK(4)+DTYPE(wpin)(i,j,m)*MASK(5)+DTYPE(wpin)(i+1,j,m)*MASK(6)
                   sob(i,1)=sob(i,1)-DTYPE(wpin)(i-1,j,m)+DTYPE(wpin)(i,j,m)*MASK(5)-DTYPE(wpin)(i+1,j,m)
                ENDDO
                j=1
                DO i=1,Nm(1)
!                  sob(i,4)=sob(i,4)+DTYPE(wpin)(i-1,j,m)*MASK( 4)+DTYPE(wpin)(i,j,m)*MASK( 5)+DTYPE(wpin)(i+1,j,m)*MASK( 6)
!                  sob(i,5)=sob(i,5)+DTYPE(wpin)(i-1,j,m)*MASK(10)+DTYPE(wpin)(i,j,m)*MASK(11)+DTYPE(wpin)(i+1,j,m)*MASK(12)
                   sob(i,4)=sob(i,4)-DTYPE(wpin)(i-1,j,m)         +DTYPE(wpin)(i,j,m)*MASK( 5)-DTYPE(wpin)(i+1,j,m)
                   sob(i,5)=sob(i,5)+DTYPE(wpin)(i-1,j,m)*MASK(10)+DTYPE(wpin)(i,j,m)*MASK(11)+DTYPE(wpin)(i+1,j,m)*MASK(12)
                ENDDO

                DO j=1,Nm(2)
                   jp=j+1
                   l=MOD(jp*3,9)+1

                   m=k-1
                   l2=l+1
                   l3=l+2
                   DO i=1,Nm(1)
!                     sob(i,l )=          DTYPE(wpin)(i-1,j,m)*MASK( 1)+DTYPE(wpin)(i,j,m)*MASK( 2)+DTYPE(wpin)(i+1,j,m)*MASK( 3)
!                     sob(i,l2)=          DTYPE(wpin)(i-1,j,m)*MASK( 7)+DTYPE(wpin)(i,j,m)*MASK( 8)+DTYPE(wpin)(i+1,j,m)*MASK( 9)
!                     sob(i,l3)=          DTYPE(wpin)(i-1,j,m)*MASK(13)+DTYPE(wpin)(i,j,m)*MASK(14)+DTYPE(wpin)(i+1,j,m)*MASK(15)
                      sob(i,l )=          DTYPE(wpin)(i-1,j,m)         +DTYPE(wpin)(i,j,m)*MASK( 2)+DTYPE(wpin)(i+1,j,m)
                      sob(i,l2)=          DTYPE(wpin)(i-1,j,m)*MASK( 7)+DTYPE(wpin)(i,j,m)*MASK( 8)+DTYPE(wpin)(i+1,j,m)*MASK( 9)
                      sob(i,l3)=          DTYPE(wpin)(i-1,j,m)         +DTYPE(wpin)(i,j,m)*MASK(14)+DTYPE(wpin)(i+1,j,m)
                   ENDDO
                   m=k+1
                   DO i=1,Nm(1)
!                     sob(i,l )=sob(i,l )+DTYPE(wpin)(i-1,j,m)*MASK( 4)+DTYPE(wpin)(i,j,m)*MASK( 5)+DTYPE(wpin)(i+1,j,m)*MASK( 6)
!                     sob(i,l2)=sob(i,l2)+DTYPE(wpin)(i-1,j,m)*MASK(10)+DTYPE(wpin)(i,j,m)*MASK(11)+DTYPE(wpin)(i+1,j,m)*MASK(12)
!                     sob(i,l3)=sob(i,l3)+DTYPE(wpin)(i-1,j,m)*MASK(16)+DTYPE(wpin)(i,j,m)*MASK(17)+DTYPE(wpin)(i+1,j,m)*MASK(18)
                      sob(i,l )=sob(i,l )-DTYPE(wpin)(i-1,j,m)         +DTYPE(wpin)(i,j,m)*MASK( 5)-DTYPE(wpin)(i+1,j,m)
                      sob(i,l2)=sob(i,l2)+DTYPE(wpin)(i-1,j,m)*MASK(10)+DTYPE(wpin)(i,j,m)*MASK(11)+DTYPE(wpin)(i+1,j,m)*MASK(12)
                      sob(i,l3)=sob(i,l3)-DTYPE(wpin)(i-1,j,m)         +DTYPE(wpin)(i,j,m)*MASK(17)-DTYPE(wpin)(i+1,j,m)
                   ENDDO

                   !index l-6
                   l1=MOD(l+3,9)
                   !index l-2
                   l2=MOD(l+7,9)

                   FORALL (i=1:Nm(1))
                      wpz(i,j,k)=sob(i,l1)
                   END FORALL
                   FORALL (i=1:Nm(1))
                      wpz(i,j,k)=wpz(i,j,k)+sob(i,l2)
                   END FORALL
                   FORALL (i=1:Nm(1))
                      wpz(i,j,k)=wpz(i,j,k)+sob(i,l3)
                   END FORALL

                   !calculate the square "Gz^2"
                   FORALL (i=1:Nm(1))
                      wpz(i,j,k)=wpz(i,j,k)*wpz(i,j,k)
                   END FORALL
                ENDDO !j=1,Nm(2)
             ENDDO !k=1,Nm(3)

             DEALLOCATE(sob,STAT=info)
             or_fail_dealloc("sob")

             !At each pixel, I should keep two values of the first
             !and third row mask multiply by the pixel row
             ALLOCATE(sob(Nm(1),6),STAT=info)
             or_fail_alloc("sob")
             ! The Sobel operator in 3D
             ! Sy  ----  Gy=Sy*A
             ! -1 -3 -1    -3 -6 -3    -1 -3 -1
             !  0  0  0     0  0  0     0  0  0
             !  1  3  1     3  6  3     1  3  1
             !  ---------------->
             !  1  3  1     3  6  3     1  3  1
             !  0  0  0     0  0  0     0  0  0
             ! -1 -3 -1    -3 -6 -3    -1 -3 -1

             MASK( 1)= one
             MASK( 2)= three
             MASK( 3)= one
             MASK( 4)= three
             MASK( 5)= 6.0_MK
             MASK( 6)= three
             MASK( 7)= one
             MASK( 8)= three
             MASK( 9)= one

             MASK(10)=-one
             MASK(11)=-three
             MASK(12)=-one
             MASK(13)=-three
             MASK(14)=-6.0_MK
             MASK(15)=-three
             MASK(16)=-one
             MASK(17)=-three
             MASK(18)=-one

             DO k=1,Nm(3)
                m=k-1
                j=0
                DO i=1,Nm(1)
!                  sob(i,1)=         DTYPE(wpin)(i-1,j,m)*MASK(1)+DTYPE(wpin)(i,j,m)*MASK(2)+DTYPE(wpin)(i+1,j,m)*MASK(3)
                   sob(i,1)=         DTYPE(wpin)(i-1,j,m)+DTYPE(wpin)(i,j,m)*MASK(2)+DTYPE(wpin)(i+1,j,m)
                ENDDO
                j=1
                DO i=1,Nm(1)
!                  sob(i,3)=         DTYPE(wpin)(i-1,j,m)*MASK( 1)+DTYPE(wpin)(i,j,m)*MASK( 2)+DTYPE(wpin)(i+1,j,m)*MASK( 3)
                   sob(i,3)=         DTYPE(wpin)(i-1,j,m)+DTYPE(wpin)(i,j,m)*MASK( 2)+DTYPE(wpin)(i+1,j,m)
                ENDDO
                !!!
                m=k
                j=0
                DO i=1,Nm(1)
                   sob(i,1)=sob(i,1)+DTYPE(wpin)(i-1,j,m)*MASK(4)+DTYPE(wpin)(i,j,m)*MASK(5)+DTYPE(wpin)(i+1,j,m)*MASK(6)
                ENDDO
                j=1
                DO i=1,Nm(1)
                   sob(i,3)=sob(i,3)+DTYPE(wpin)(i-1,j,m)*MASK( 4)+DTYPE(wpin)(i,j,m)*MASK( 5)+DTYPE(wpin)(i+1,j,m)*MASK( 6)
                ENDDO
                !!!
                m=k+1
                j=0
                DO i=1,Nm(1)
!                  sob(i,1)=sob(i,1)+DTYPE(wpin)(i-1,j,m)*MASK(7)+DTYPE(wpin)(i,j,m)*MASK(8)+DTYPE(wpin)(i+1,j,m)*MASK(9)
                   sob(i,1)=sob(i,1)+DTYPE(wpin)(i-1,j,m)+DTYPE(wpin)(i,j,m)*MASK(8)+DTYPE(wpin)(i+1,j,m)
                ENDDO
                j=1
                DO i=1,Nm(1)
!                  sob(i,3)=sob(i,3)+DTYPE(wpin)(i-1,j,m)*MASK( 7)+DTYPE(wpin)(i,j,m)*MASK( 8)+DTYPE(wpin)(i+1,j,m)*MASK( 9)
                   sob(i,3)=sob(i,3)+DTYPE(wpin)(i-1,j,m)+DTYPE(wpin)(i,j,m)*MASK( 8)+DTYPE(wpin)(i+1,j,m)
                ENDDO

                DO j=1,Nm(2)
                   jp=j+1
                   l=MOD(jp*2,6)+1

                   m=k-1
                   l2=l+1
                   DO i=1,Nm(1)
!                     sob(i,l )=          DTYPE(wpin)(i-1,j,m)*MASK( 1)+DTYPE(wpin)(i,j,m)*MASK( 2)+DTYPE(wpin)(i+1,j,m)*MASK( 3)
!                     sob(i,l2)=          DTYPE(wpin)(i-1,j,m)*MASK(10)+DTYPE(wpin)(i,j,m)*MASK(11)+DTYPE(wpin)(i+1,j,m)*MASK(12)
                      sob(i,l )=          DTYPE(wpin)(i-1,j,m)+DTYPE(wpin)(i,j,m)*MASK( 2)+DTYPE(wpin)(i+1,j,m)
                      sob(i,l2)=         -DTYPE(wpin)(i-1,j,m)+DTYPE(wpin)(i,j,m)*MASK(11)-DTYPE(wpin)(i+1,j,m)
                   ENDDO

                   DO i=1,Nm(1)
                      sob(i,l )=sob(i,l )+DTYPE(wpin)(i-1,j,k)*MASK( 4)+DTYPE(wpin)(i,j,k)*MASK( 5)+DTYPE(wpin)(i+1,j,k)*MASK( 6)
                      sob(i,l2)=sob(i,l2)+DTYPE(wpin)(i-1,j,k)*MASK(13)+DTYPE(wpin)(i,j,k)*MASK(14)+DTYPE(wpin)(i+1,j,k)*MASK(15)
                   ENDDO

                   m=k+1
                   DO i=1,Nm(1)
!                     sob(i,l )=sob(i,l )+DTYPE(wpin)(i-1,j,m)*MASK( 7)+DTYPE(wpin)(i,j,m)*MASK( 8)+DTYPE(wpin)(i+1,j,m)*MASK( 9)
!                     sob(i,l2)=sob(i,l2)+DTYPE(wpin)(i-1,j,m)*MASK(16)+DTYPE(wpin)(i,j,m)*MASK(17)+DTYPE(wpin)(i+1,j,m)*MASK(18)
                      sob(i,l )=sob(i,l )+DTYPE(wpin)(i-1,j,m)+DTYPE(wpin)(i,j,m)*MASK( 8)+DTYPE(wpin)(i+1,j,m)
                      sob(i,l2)=sob(i,l2)-DTYPE(wpin)(i-1,j,m)+DTYPE(wpin)(i,j,m)*MASK(17)-DTYPE(wpin)(i+1,j,m)
                   ENDDO

                   l1=MOD(l+2,6)

                   FORALL (i=1:Nm(1))
                      DTYPE(wpy)(i,j,k)=sob(i,l1)
                   END FORALL
                   FORALL (i=1:Nm(1))
                      DTYPE(wpy)(i,j,k)=DTYPE(wpy)(i,j,k)+sob(i,l2)
                   END FORALL

                   !calculate the square "Gy^2"
                   FORALL (i=1:Nm(1))
                      DTYPE(wpy)(i,j,k)=DTYPE(wpy)(i,j,k)*DTYPE(wpy)(i,j,k)
                   END FORALL
                ENDDO !j=1,Nm(2)
             ENDDO !k=1,Nm(3)

             DEALLOCATE(sob,STAT=info)
             or_fail_dealloc("sob")

             !calculate the sum "Gy^2+Gz^2" to use the memory
             FORALL (i=1:Nm(1),j=1:Nm(2),k=1:Nm(3))
                DTYPE(wpy)(i,j,k)=DTYPE(wpy)(i,j,k)+wpz(i,j,k)
             END FORALL

             !At each pixel, I should keep two values of the first
             !and third row mask multiply by the pixel row
             ALLOCATE(sob(Nm(1),6),STAT=info)
             or_fail_alloc("sob")
             ! The Sobel operator in 3D
             ! Sx  ----  Gx=Sx*A
             ! -1  0  1    -3  0  3    -1  0  1
             ! -3  0  3    -6  0  6    -3  0  3
             ! -1  0  1    -3  0  3    -1  0  1
             !  ---------------->
             !  1  0 -1     3  0 -3     1  0 -1
             !  3  0 -3     6  0 -6     3  0 -3
             !  1  0 -1     3  0 -3     1  0 -1

             MASK( 1)= one
             MASK( 2)=-one
             MASK( 3)= three
             MASK( 4)=-three
             MASK( 5)= one
             MASK( 6)=-one

             MASK( 7)= three
             MASK( 8)=-three
             MASK( 9)= 6.0_MK
             MASK(10)=-6.0_MK
             MASK(11)= three
             MASK(12)=-three

             MASK(13)= one
             MASK(14)=-one
             MASK(15)= three
             MASK(16)=-three
             MASK(17)= one
             MASK(18)=-one

             DO k=1,Nm(3)
                m=k-1
                j=0
                DO i=1,Nm(1)
!                  sob(i,1)=         DTYPE(wpin)(i-1,j,m)*MASK(1)+DTYPE(wpin)(i+1,j,m)*MASK(2)
                   sob(i,1)=         DTYPE(wpin)(i-1,j,m)-DTYPE(wpin)(i+1,j,m)
                ENDDO
                j=1
                DO i=1,Nm(1)
!                  sob(i,4)=         DTYPE(wpin)(i-1,j,m)*MASK( 1)+DTYPE(wpin)(i+1,j,m)*MASK( 2)
!                  sob(i,5)=         DTYPE(wpin)(i-1,j,m)*MASK( 7)+DTYPE(wpin)(i+1,j,m)*MASK( 8)
                   sob(i,4)=         DTYPE(wpin)(i-1,j,m)         -DTYPE(wpin)(i+1,j,m)
                   sob(i,5)=         DTYPE(wpin)(i-1,j,m)*MASK( 7)+DTYPE(wpin)(i+1,j,m)*MASK( 8)
                ENDDO
                !!!
                m=k
                j=0
                DO i=1,Nm(1)
                   sob(i,1)=sob(i,1)+DTYPE(wpin)(i-1,j,m)*MASK(3)+DTYPE(wpin)(i+1,j,m)*MASK(4)
                ENDDO
                j=1
                DO i=1,Nm(1)
                   sob(i,4)=sob(i,4)+DTYPE(wpin)(i-1,j,m)*MASK( 3)+DTYPE(wpin)(i+1,j,m)*MASK( 4)
                   sob(i,5)=sob(i,5)+DTYPE(wpin)(i-1,j,m)*MASK( 9)+DTYPE(wpin)(i+1,j,m)*MASK(10)
                ENDDO
                !!!
                m=k+1
                j=0
                DO i=1,Nm(1)
!                  sob(i,1)=sob(i,1)+DTYPE(wpin)(i-1,j,m)*MASK(5)+DTYPE(wpin)(i+1,j,m)*MASK(6)
                   sob(i,1)=sob(i,1)+DTYPE(wpin)(i-1,j,m)-DTYPE(wpin)(i+1,j,m)
                ENDDO
                j=1
                DO i=1,Nm(1)
!                  sob(i,4)=sob(i,4)+DTYPE(wpin)(i-1,j,m)*MASK( 5)+DTYPE(wpin)(i+1,j,m)*MASK( 6)
!                  sob(i,5)=sob(i,5)+DTYPE(wpin)(i-1,j,m)*MASK(11)+DTYPE(wpin)(i+1,j,m)*MASK(12)
                   sob(i,4)=sob(i,4)+DTYPE(wpin)(i-1,j,m)         -DTYPE(wpin)(i+1,j,m)
                   sob(i,5)=sob(i,5)+DTYPE(wpin)(i-1,j,m)*MASK(11)+DTYPE(wpin)(i+1,j,m)*MASK(12)
                ENDDO
                DO j=1,Nm(2)
                   jp=j+1
                   l=MOD(jp*3,9)+1

                   m=k-1
                   l2=l+1
                   l3=l+2
                   DO i=1,Nm(1)
!                     sob(i,l )=          DTYPE(wpin)(i-1,j,m)*MASK( 1)+DTYPE(wpin)(i+1,j,m)*MASK( 2)
!                     sob(i,l2)=          DTYPE(wpin)(i-1,j,m)*MASK( 7)+DTYPE(wpin)(i+1,j,m)*MASK( 8)
!                     sob(i,l3)=          DTYPE(wpin)(i-1,j,m)*MASK(13)+DTYPE(wpin)(i+1,j,m)*MASK(14)
                      sob(i,l )=          DTYPE(wpin)(i-1,j,m)         -DTYPE(wpin)(i+1,j,m)
                      sob(i,l2)=          DTYPE(wpin)(i-1,j,m)*MASK( 7)+DTYPE(wpin)(i+1,j,m)*MASK( 8)
                      sob(i,l3)=          DTYPE(wpin)(i-1,j,m)         -DTYPE(wpin)(i+1,j,m)
                   ENDDO
                   DO i=1,Nm(1)
                      sob(i,l )=sob(i,l )+DTYPE(wpin)(i-1,j,k)*MASK( 3)+DTYPE(wpin)(i+1,j,k)*MASK( 4)
                      sob(i,l2)=sob(i,l2)+DTYPE(wpin)(i-1,j,k)*MASK( 9)+DTYPE(wpin)(i+1,j,k)*MASK(10)
                      sob(i,l3)=sob(i,l3)+DTYPE(wpin)(i-1,j,k)*MASK(15)+DTYPE(wpin)(i+1,j,k)*MASK(16)
                   ENDDO
                   m=k+1
                   DO i=1,Nm(1)
!                     sob(i,l )=sob(i,l )+DTYPE(wpin)(i-1,j,m)*MASK( 5)+DTYPE(wpin)(i+1,j,m)*MASK( 6)
!                     sob(i,l2)=sob(i,l2)+DTYPE(wpin)(i-1,j,m)*MASK(11)+DTYPE(wpin)(i+1,j,m)*MASK(12)
!                     sob(i,l3)=sob(i,l3)+DTYPE(wpin)(i-1,j,m)*MASK(17)+DTYPE(wpin)(i+1,j,m)*MASK(18)
                      sob(i,l )=sob(i,l )+DTYPE(wpin)(i-1,j,m)         -DTYPE(wpin)(i+1,j,m)
                      sob(i,l2)=sob(i,l2)+DTYPE(wpin)(i-1,j,m)*MASK(11)+DTYPE(wpin)(i+1,j,m)*MASK(12)
                      sob(i,l3)=sob(i,l3)+DTYPE(wpin)(i-1,j,m)         -DTYPE(wpin)(i+1,j,m)
                   ENDDO

                   !index l-6
                   l1=MOD(l+3,9)
                   !index l-2
                   l2=MOD(l+7,9)

                   FORALL (i=1:Nm(1))
                      wpz(i,j,k)=sob(i,l1)
                   END FORALL
                   FORALL (i=1:Nm(1))
                      wpz(i,j,k)=wpz(i,j,k)+sob(i,l2)
                   END FORALL
                   FORALL (i=1:Nm(1))
                      wpz(i,j,k)=wpz(i,j,k)+sob(i,l3)
                   END FORALL

                   !calculate the square "Gx^2"
                   FORALL (i=1:Nm(1))
                      wpz(i,j,k)=wpz(i,j,k)*wpz(i,j,k)
                   END FORALL
                ENDDO !j=1,Nm(2)
             ENDDO !k=1,Nm(3)

             DEALLOCATE(sob,MASK,STAT=info)
             or_fail_dealloc("sob & MASK")

             ! the gradient magnitud G=sqrt(Gx^2+Gy^2+Gz^2)
             FORALL (i=1:Nm(1),j=1:Nm(2),k=1:Nm(3))
                DTYPE(wpout)(i,j,k)=SQRT(wpz(i,j,k)+DTYPE(wpy)(i,j,k))
             END FORALL
#endif

             sbpitr => MeshIn%subpatch%next()
          ENDDO
#if   __DIME == __2D
          NULLIFY(DTYPE(wpin),DTYPE(wpout),DTYPE(wpy))
#elif __DIME == __3D
          NULLIFY(DTYPE(wpin),DTYPE(wpout),DTYPE(wpy),wpz)
#endif

          CALL fld1%destroy(info)
          or_fail("Destroy fld1 field failed!")

          ! deallocate the pointer:
          dealloc_pointer("fld1")

#if   __DIME == __3D
          CALL fld2%destroy(info)
          or_fail("Destroy fld2 field failed!")

          ! deallocate the pointer:
          dealloc_pointer("fld2")
#endif

          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        CONTAINS
        SUBROUTINE check
        IMPLICIT NONE
          check_true(<#ALL(MeshIn%ghostsize.GE.1)#>, &
          & "There should be at least one layer of ghost at the borders!",exit_point=8888)

          IF (PRESENT(FieldOut)) THEN
             IF (.NOT.ASSOCIATED(FieldIn,FieldOut)) THEN
                check_true(<#FieldOut%data_type.EQ.FieldIn%data_type#>, &
                & "The Input and output fields have different data types!",exit_point=8888)
                IF (.NOT.FieldOut%is_discretized_on(MeshIn)) THEN
                   CALL FieldOut%discretize_on(MeshIn,info)
                   or_fail("Failed to discretize Field on Mesh!",exit_point=8888)
                ENDIF
             ENDIF
          ENDIF

          IF (.NOT.ASSOCIATED(fld1)) THEN
             ALLOCATE(ppm_t_field::fld1,STAT=info)
             or_fail_alloc('Failed to allocate fld1.',exit_point=8888)

             CALL fld1%create(1,info,dtype=FieldIn%data_type,name="Ygradient")
             or_fail("Failed to create field!",exit_point=8888)

             CALL fld1%discretize_on(MeshIn,info)
             or_fail("Failed to discretize Field on Mesh!",exit_point=8888)
          ENDIF
#if   __DIME == __3D
          IF (.NOT.ASSOCIATED(fld2)) THEN
             ALLOCATE(ppm_t_field::fld2,STAT=info)
             or_fail_alloc('Failed to allocate fld2.',exit_point=8888)

             CALL fld2%create(1,info,dtype=FieldIn%data_type,name="Zgradient")
             or_fail("Failed to create field!",exit_point=8888)

             CALL fld2%discretize_on(MeshIn,info)
             or_fail("Failed to discretize Field on Mesh!",exit_point=8888)
          ENDIF
#endif
          !-------------------------------------------------------------------------
          !  Return
          !-------------------------------------------------------------------------
        8888 CONTINUE
          RETURN
        END SUBROUTINE check
        END SUBROUTINE DTYPE(ppm_rc_SobelImageFilter)

