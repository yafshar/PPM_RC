        !-------------------------------------------------------------------------
        !  Subroutine   :                  ppm_rc_convolve
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
        !  Author           - y.afshar           Oct       2014
        !                                        Dec       2015
        !-------------------------------------------------------------------------

        !-------------------------------------------------------------------------
        !  Subroutine   :                    ppm_rc_convolve
        !-------------------------------------------------------------------------
        !
        !  Purpose      : Convolve image by a mask.
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
        !                inplace convolution in 3D does not work!!!!
        !
        !  References   :
        !
        !  Revisions    :
        !-------------------------------------------------------------------------
        SUBROUTINE DTYPE(CTYPE(ppm_rc_convolve))(MaskIn,FieldIn,MeshIn,info, &
        &          FieldOut,boundary_condition,is_normalized,Filter)

          IMPLICIT NONE
          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------

          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
#if   __DIME == __2D
          __TYPE, DIMENSION(:,:),        INTENT(IN   ) :: MaskIn
#elif __DIME == __3D
          __TYPE, DIMENSION(:,:,:),      INTENT(IN   ) :: MaskIn
#endif
          CLASS(ppm_t_field_),           POINTER       :: FieldIn
          !!! the source filed which is supposed to be

          CLASS(ppm_t_equi_mesh_),       POINTER       :: MeshIn
          !!! source mesh on which fields are descritized

          INTEGER,                       INTENT(  OUT) :: info

          CLASS(ppm_t_field_), OPTIONAL, POINTER       :: FieldOut
          !!! the output filed which is supposed to be

          INTEGER,             OPTIONAL, INTENT(IN   ) :: boundary_condition
          !!!not implemented

          LOGICAL,             OPTIONAL, INTENT(IN   ) :: is_normalized
          !!!not implemented
          LOGICAL,             OPTIONAL, INTENT(IN   ) :: Filter
          !!!Use the correlation or convolution

          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          CLASS(ppm_t_subpatch_), POINTER :: sbpitr

          REAL(ppm_kind_double) :: t0

#if   __DIME == __2D
          __TYPE, CONTIGUOUS, DIMENSION(:,:),   POINTER :: DTYPE(wpin)
          __TYPE, CONTIGUOUS, DIMENSION(:,:),   POINTER :: DTYPE(wpout)
#elif __DIME == __3D
          __TYPE, CONTIGUOUS, DIMENSION(:,:,:), POINTER :: DTYPE(wpin)
          __TYPE, CONTIGUOUS, DIMENSION(:,:,:), POINTER :: DTYPE(wpout)
#endif

#if   __DIME == __2D
          __TYPE, DIMENSION(SIZE(MaskIn,1),SIZE(MaskIn,2))                :: Mask_
#elif __DIME == __3D
          __TYPE, DIMENSION(SIZE(MaskIn,1),SIZE(MaskIn,2),SIZE(MaskIn,2)) :: Mask_
#endif

          INTEGER, DIMENSION(:), POINTER :: nnodes
          INTEGER, DIMENSION(__DIME)          :: msize
          INTEGER                        :: i,j,s1,s2
#if   __DIME == __3D
          INTEGER                        :: k,s3
#endif
          CHARACTER(LEN=ppm_char) :: caller = 'ppm_rc_convolve'

          LOGICAL :: Filter_
          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          Filter_=MERGE(Filter,.FALSE.,PRESENT(Filter))
          !Per default we are convolving the input with a mask

          s1=SIZE(MaskIn,1)
          s2=SIZE(MaskIn,2)
#if   __DIME == __3D
          s3=SIZE(MaskIn,3)
#endif

          SELECT CASE (Filter_)
          CASE (.FALSE.)
          !!! For convolution we need to reorder the input MASK array
#if   __DIME == __2D
             DO j=1,s2
                DO i=1,s1
                   Mask_(i,j)=MaskIn(s1-i+1,s2-j+1)
                ENDDO
             ENDDO
#elif __DIME == __3D
             Do k=1,s3
                DO j=1,s2
                   DO i=1,s1
                      Mask_(i,j,k)=MaskIn(s1-i+1,s2-j+1,s3-k+1)
                   ENDDO
                ENDDO
             ENDDO
#endif
         END SELECT

          msize(1)=(s1-1)/2
          msize(2)=(s2-1)/2
#if   __DIME == __3D
          msize(3)=(s3-1)/2
#endif

          CALL check()

          !currently we suppose that the field has been padded at the borders
          sbpitr => MeshIn%subpatch%begin()
          DO WHILE (ASSOCIATED(sbpitr))
             nnodes => sbpitr%nnodes

             NULLIFY(DTYPE(wpin))
             CALL sbpitr%get_field(FieldIn,DTYPE(wpin),info)
             or_fail("Failed to get field wpin data.")

             IF (PRESENT(FieldOut)) THEN
                NULLIFY(DTYPE(wpout))
                CALL sbpitr%get_field(FieldOut,DTYPE(wpout),info)
                or_fail("Failed to get field wpout data.")
             ELSE
                DTYPE(wpout) => DTYPE(wpin)
             ENDIF

#if   __DIME == __3D
             IF (ASSOCIATED(DTYPE(wpout),DTYPE(wpin))) THEN
                fail("For this version we can not do inplace convolution in 3D!",ppm_error=ppm_error_fatal)
             ENDIF
#endif

             SELECT CASE (Filter_)
             CASE (.FALSE.)
                CALL CTYPE(Update)(Mask_,DTYPE(wpin),DTYPE(wpout),nnodes,msize,info)
             CASE (.TRUE.)
                CALL CTYPE(Update)(MaskIn,DTYPE(wpin),DTYPE(wpout),nnodes,msize,info)
             END SELECT
             or_fail("Update")

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
          check_true(<#ALL(MeshIn%ghostsize(1:__DIME).GE.msize)#>, &
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
        SUBROUTINE CTYPE(Update)(Mask,wpin,wpout,nnodes_,msize_,info)
          IMPLICIT NONE
          __TYPE,             DIMENSION(:,:), INTENT(IN   ) :: Mask
          __TYPE, CONTIGUOUS, DIMENSION(:,:), POINTER       :: wpin
          __TYPE, CONTIGUOUS, DIMENSION(:,:), POINTER       :: wpout
          INTEGER,            DIMENSION(:),   INTENT(IN   ) :: nnodes_
          INTEGER,            DIMENSION(:),   INTENT(IN   ) :: msize_
          INTEGER,                            INTENT(  OUT) :: info

          INTEGER :: jp,l1,l2,l3,l4,l5,l6,l7
          INTEGER :: i,j,l,m,n

          IF (ALL(msize_.EQ.1)) THEN
             !At each pixel, I should keep three values of the first, second
             !and third row mask multiply by the pixel row
             ALLOCATE(CTYPE(tmp2)(nnodes_(1),9),STAT=info)
             or_fail_alloc("tmp2",exit_point=7777)

             ! This optimisation is done to avoid recomputing what has been computed once
             j=0
             DO i=1,nnodes_(1)
                CTYPE(tmp2)(i,1)=wpin(i-1,j)*Mask(1,1)+wpin(i,j)*Mask(2,1)+wpin(i+1,j)*Mask(3,1)
             ENDDO
             j=1
             DO i=1,nnodes_(1)
                CTYPE(tmp2)(i,4)=wpin(i-1,j)*Mask(1,1)+wpin(i,j)*Mask(2,1)+wpin(i+1,j)*Mask(3,1)
                CTYPE(tmp2)(i,5)=wpin(i-1,j)*Mask(1,2)+wpin(i,j)*Mask(2,2)+wpin(i+1,j)*Mask(3,2)
             ENDDO
             DO j=1,nnodes_(2)
                jp=j+1
                l=MOD(jp*3,9)+1
                l2=l+1
                l3=l+2
                DO i=1,nnodes_(1)
                   CTYPE(tmp2)(i,l )=wpin(i-1,jp)*Mask(1,1)+wpin(i,jp)*Mask(2,1)+wpin(i+1,jp)*Mask(3,1)
                   CTYPE(tmp2)(i,l2)=wpin(i-1,jp)*Mask(1,2)+wpin(i,jp)*Mask(2,2)+wpin(i+1,jp)*Mask(3,2)
                   CTYPE(tmp2)(i,l3)=wpin(i-1,jp)*Mask(1,3)+wpin(i,jp)*Mask(2,3)+wpin(i+1,jp)*Mask(3,3)
                ENDDO

                l1=MOD(l+3,9)
                l2=MOD(l+7,9)

                FORALL (i=1:nnodes_(1))
                   wpout(i,j)=CTYPE(tmp2)(i,l1)
                END FORALL
                FORALL (i=1:nnodes_(1))
                   wpout(i,j)=wpout(i,j)+CTYPE(tmp2)(i,l2)
                END FORALL
                FORALL (i=1:nnodes_(1))
                   wpout(i,j)=wpout(i,j)+CTYPE(tmp2)(i,l3)
                END FORALL
             ENDDO !j=1,nnodes_(2)

             DEALLOCATE(CTYPE(tmp2),STAT=info)
             or_fail_dealloc("tmp2",exit_point=7777)
          ELSE IF (ALL(msize_.EQ.2)) THEN
             ALLOCATE(CTYPE(tmp2)(nnodes_(1),25),STAT=info)
             or_fail_alloc("CTYPE(tmp2)",exit_point=7777)

             j=-1
             DO i=1,nnodes_(1)
                CTYPE(tmp2)(i,1)=wpin(i-2,j)*Mask(1,1)+wpin(i-1,j)*Mask(2,1)+wpin(i,j)*Mask(3,1)+wpin(i+1,j)*Mask(4,1)+wpin(i+2,j)*Mask(5,1)
             ENDDO
             j=0
             DO i=1,nnodes_(1)
                CTYPE(tmp2)(i,6)=wpin(i-2,j)*Mask(1,1)+wpin(i-1,j)*Mask(2,1)+wpin(i,j)*Mask(3,1)+wpin(i+1,j)*Mask(4,1)+wpin(i+2,j)*Mask(5,1)
                CTYPE(tmp2)(i,7)=wpin(i-2,j)*Mask(1,2)+wpin(i-1,j)*Mask(2,2)+wpin(i,j)*Mask(3,2)+wpin(i+1,j)*Mask(4,2)+wpin(i+2,j)*Mask(5,2)
             ENDDO
             j=1
             DO i=1,nnodes_(1)
                CTYPE(tmp2)(i,11)=wpin(i-2,j)*Mask(1,1)+wpin(i-1,j)*Mask(2,1)+wpin(i,j)*Mask(3,1)+wpin(i+1,j)*Mask(4,1)+wpin(i+2,j)*Mask(5,1)
                CTYPE(tmp2)(i,12)=wpin(i-2,j)*Mask(1,2)+wpin(i-1,j)*Mask(2,2)+wpin(i,j)*Mask(3,2)+wpin(i+1,j)*Mask(4,2)+wpin(i+2,j)*Mask(5,2)
                CTYPE(tmp2)(i,13)=wpin(i-2,j)*Mask(1,3)+wpin(i-1,j)*Mask(2,3)+wpin(i,j)*Mask(3,3)+wpin(i+1,j)*Mask(4,3)+wpin(i+2,j)*Mask(5,3)
             ENDDO
             j=2
             DO i=1,nnodes_(1)
                CTYPE(tmp2)(i,16)=wpin(i-2,j)*Mask(1,1)+wpin(i-1,j)*Mask(2,1)+wpin(i,j)*Mask(3,1)+wpin(i+1,j)*Mask(4,1)+wpin(i+2,j)*Mask(5,1)
                CTYPE(tmp2)(i,17)=wpin(i-2,j)*Mask(1,2)+wpin(i-1,j)*Mask(2,2)+wpin(i,j)*Mask(3,2)+wpin(i+1,j)*Mask(4,2)+wpin(i+2,j)*Mask(5,2)
                CTYPE(tmp2)(i,18)=wpin(i-2,j)*Mask(1,3)+wpin(i-1,j)*Mask(2,3)+wpin(i,j)*Mask(3,3)+wpin(i+1,j)*Mask(4,3)+wpin(i+2,j)*Mask(5,3)
                CTYPE(tmp2)(i,19)=wpin(i-2,j)*Mask(1,4)+wpin(i-1,j)*Mask(2,4)+wpin(i,j)*Mask(3,4)+wpin(i+1,j)*Mask(4,4)+wpin(i+2,j)*Mask(5,4)
             ENDDO

             DO j=1,nnodes_(2)
                jp=j+2
                l=MOD((jp+1)*5,25)+1
                l2=l+1
                l3=l+2
                l4=l+3
                l5=l+4
                DO i=1,nnodes_(1)
                   CTYPE(tmp2)(i,l )=wpin(i-2,jp)*Mask(1,1)+wpin(i-1,jp)*Mask(2,1)+wpin(i,jp)*Mask(3,1)+wpin(i+1,jp)*Mask(4,1)+wpin(i+2,jp)*Mask(5,1)
                   CTYPE(tmp2)(i,l2)=wpin(i-2,jp)*Mask(1,2)+wpin(i-1,jp)*Mask(2,2)+wpin(i,jp)*Mask(3,2)+wpin(i+1,jp)*Mask(4,2)+wpin(i+2,jp)*Mask(5,2)
                   CTYPE(tmp2)(i,l3)=wpin(i-2,jp)*Mask(1,3)+wpin(i-1,jp)*Mask(2,3)+wpin(i,jp)*Mask(3,3)+wpin(i+1,jp)*Mask(4,3)+wpin(i+2,jp)*Mask(5,3)
                   CTYPE(tmp2)(i,l4)=wpin(i-2,jp)*Mask(1,4)+wpin(i-1,jp)*Mask(2,4)+wpin(i,jp)*Mask(3,4)+wpin(i+1,jp)*Mask(4,4)+wpin(i+2,jp)*Mask(5,4)
                   CTYPE(tmp2)(i,l5)=wpin(i-2,jp)*Mask(1,5)+wpin(i-1,jp)*Mask(2,5)+wpin(i,jp)*Mask(3,5)+wpin(i+1,jp)*Mask(4,5)+wpin(i+2,jp)*Mask(5,5)
                ENDDO

                l1=MOD(l+ 5,25)
                l2=MOD(l+11,25)
                l3=MOD(l+17,25)
                l4=MOD(l+23,25)

                FORALL (i=1:nnodes_(1))
                   wpout(i,j)=CTYPE(tmp2)(i,l1)
                END FORALL
                FORALL (i=1:nnodes_(1))
                   wpout(i,j)=wpout(i,j)+CTYPE(tmp2)(i,l2)
                END FORALL
                FORALL (i=1:nnodes_(1))
                   wpout(i,j)=wpout(i,j)+CTYPE(tmp2)(i,l3)
                END FORALL
                FORALL (i=1:nnodes_(1))
                   wpout(i,j)=wpout(i,j)+CTYPE(tmp2)(i,l4)
                END FORALL
                FORALL (i=1:nnodes_(1))
                   wpout(i,j)=wpout(i,j)+CTYPE(tmp2)(i,l5)
                END FORALL
             ENDDO !j=1,nnodes_(2)

             DEALLOCATE(CTYPE(tmp2),STAT=info)
             or_fail_dealloc("tmp2",exit_point=7777)
          ELSE IF (ALL(msize_.EQ.3)) THEN
             ALLOCATE(CTYPE(tmp2)(nnodes_(1),49),STAT=info)
             or_fail_alloc("tmp2",exit_point=7777)

             j=-2
             DO i=1,nnodes_(1)
                CTYPE(tmp2)(i,1)=wpin(i-3,j)*Mask(1,1)+wpin(i-2,j)*Mask(2,1)+wpin(i-1,j)*Mask(3,1)+ &
                &                wpin(i,  j)*Mask(4,1)+wpin(i+1,j)*Mask(5,1)+wpin(i+2,j)*Mask(6,1)+ &
                &                wpin(i+3,j)*Mask(7,1)
             ENDDO
             j=-1
             DO i=1,nnodes_(1)
                CTYPE(tmp2)(i,8)=wpin(i-3,j)*Mask(1,1)+wpin(i-2,j)*Mask(2,1)+wpin(i-1,j)*Mask(3,1)+ &
                &                wpin(i,  j)*Mask(4,1)+wpin(i+1,j)*Mask(5,1)+wpin(i+2,j)*Mask(6,1)+ &
                &                wpin(i+3,j)*Mask(7,1)
                CTYPE(tmp2)(i,9)=wpin(i-3,j)*Mask(1,2)+wpin(i-2,j)*Mask(2,2)+wpin(i-1,j)*Mask(3,2)+ &
                &                wpin(i,  j)*Mask(4,2)+wpin(i+1,j)*Mask(5,2)+wpin(i+2,j)*Mask(6,2)+ &
                &                wpin(i+3,j)*Mask(7,2)
             ENDDO
             j=0
             DO i=1,nnodes_(1)
                CTYPE(tmp2)(i,15)=wpin(i-3,j)*Mask(1,1)+wpin(i-2,j)*Mask(2,1)+wpin(i-1,j)*Mask(3,1)+ &
                &                 wpin(i,  j)*Mask(4,1)+wpin(i+1,j)*Mask(5,1)+wpin(i+2,j)*Mask(6,1)+ &
                &                 wpin(i+3,j)*Mask(7,1)
                CTYPE(tmp2)(i,16)=wpin(i-3,j)*Mask(1,2)+wpin(i-2,j)*Mask(2,2)+wpin(i-1,j)*Mask(3,2)+ &
                &                 wpin(i,  j)*Mask(4,2)+wpin(i+1,j)*Mask(5,2)+wpin(i+2,j)*Mask(6,2)+ &
                &                 wpin(i+3,j)*Mask(7,2)
                CTYPE(tmp2)(i,17)=wpin(i-3,j)*Mask(1,3)+wpin(i-2,j)*Mask(2,3)+wpin(i-1,j)*Mask(3,3)+ &
                &                 wpin(i,  j)*Mask(4,3)+wpin(i+1,j)*Mask(5,3)+wpin(i+2,j)*Mask(6,3)+ &
                &                 wpin(i+3,j)*Mask(7,3)
             ENDDO
             j=1
             DO i=1,nnodes_(1)
                CTYPE(tmp2)(i,22)=wpin(i-3,j)*Mask(1,1)+wpin(i-2,j)*Mask(2,1)+wpin(i-1,j)*Mask(3,1)+ &
                &                 wpin(i,  j)*Mask(4,1)+wpin(i+1,j)*Mask(5,1)+wpin(i+2,j)*Mask(6,1)+ &
                &                 wpin(i+3,j)*Mask(7,1)
                CTYPE(tmp2)(i,23)=wpin(i-3,j)*Mask(1,2)+wpin(i-2,j)*Mask(2,2)+wpin(i-1,j)*Mask(3,2)+ &
                &                 wpin(i,  j)*Mask(4,2)+wpin(i+1,j)*Mask(5,2)+wpin(i+2,j)*Mask(6,2)+ &
                &                 wpin(i+3,j)*Mask(7,2)
                CTYPE(tmp2)(i,24)=wpin(i-3,j)*Mask(1,3)+wpin(i-2,j)*Mask(2,3)+wpin(i-1,j)*Mask(3,3)+ &
                &                 wpin(i,  j)*Mask(4,3)+wpin(i+1,j)*Mask(5,3)+wpin(i+2,j)*Mask(6,3)+ &
                &                 wpin(i+3,j)*Mask(7,3)
                CTYPE(tmp2)(i,25)=wpin(i-3,j)*Mask(1,4)+wpin(i-2,j)*Mask(2,4)+wpin(i-1,j)*Mask(3,4)+ &
                &                 wpin(i,  j)*Mask(4,4)+wpin(i+1,j)*Mask(5,4)+wpin(i+2,j)*Mask(6,4)+ &
                &                 wpin(i+3,j)*Mask(7,4)
             ENDDO
             j=2
             DO i=1,nnodes_(1)
                CTYPE(tmp2)(i,29)=wpin(i-3,j)*Mask(1,1)+wpin(i-2,j)*Mask(2,1)+wpin(i-1,j)*Mask(3,1)+ &
                &                 wpin(i,  j)*Mask(4,1)+wpin(i+1,j)*Mask(5,1)+wpin(i+2,j)*Mask(6,1)+ &
                &                 wpin(i+3,j)*Mask(7,1)
                CTYPE(tmp2)(i,30)=wpin(i-3,j)*Mask(1,2)+wpin(i-2,j)*Mask(2,2)+wpin(i-1,j)*Mask(3,2)+ &
                &                 wpin(i,  j)*Mask(4,2)+wpin(i+1,j)*Mask(5,2)+wpin(i+2,j)*Mask(6,2)+ &
                &                 wpin(i+3,j)*Mask(7,2)
                CTYPE(tmp2)(i,31)=wpin(i-3,j)*Mask(1,3)+wpin(i-2,j)*Mask(2,3)+wpin(i-1,j)*Mask(3,3)+ &
                &                 wpin(i,  j)*Mask(4,3)+wpin(i+1,j)*Mask(5,3)+wpin(i+2,j)*Mask(6,3)+ &
                &                 wpin(i+3,j)*Mask(7,3)
                CTYPE(tmp2)(i,32)=wpin(i-3,j)*Mask(1,4)+wpin(i-2,j)*Mask(2,4)+wpin(i-1,j)*Mask(3,4)+ &
                &                 wpin(i,  j)*Mask(4,4)+wpin(i+1,j)*Mask(5,4)+wpin(i+2,j)*Mask(6,4)+ &
                &                 wpin(i+3,j)*Mask(7,4)
                CTYPE(tmp2)(i,33)=wpin(i-3,j)*Mask(1,5)+wpin(i-2,j)*Mask(2,5)+wpin(i-1,j)*Mask(3,5)+ &
                &                 wpin(i,  j)*Mask(4,5)+wpin(i+1,j)*Mask(5,5)+wpin(i+2,j)*Mask(6,5)+ &
                &                 wpin(i+3,j)*Mask(7,5)
             ENDDO
             j=3
             DO i=1,nnodes_(1)
                CTYPE(tmp2)(i,36)=wpin(i-3,j)*Mask(1,1)+wpin(i-2,j)*Mask(2,1)+wpin(i-1,j)*Mask(3,1)+ &
                &                 wpin(i,  j)*Mask(4,1)+wpin(i+1,j)*Mask(5,1)+wpin(i+2,j)*Mask(6,1)+ &
                &                 wpin(i+3,j)*Mask(7,1)
                CTYPE(tmp2)(i,37)=wpin(i-3,j)*Mask(1,2)+wpin(i-2,j)*Mask(2,2)+wpin(i-1,j)*Mask(3,2)+ &
                &                 wpin(i,  j)*Mask(4,2)+wpin(i+1,j)*Mask(5,2)+wpin(i+2,j)*Mask(6,2)+ &
                &                 wpin(i+3,j)*Mask(7,2)
                CTYPE(tmp2)(i,38)=wpin(i-3,j)*Mask(1,3)+wpin(i-2,j)*Mask(2,3)+wpin(i-1,j)*Mask(3,3)+ &
                &                 wpin(i,  j)*Mask(4,3)+wpin(i+1,j)*Mask(5,3)+wpin(i+2,j)*Mask(6,3)+ &
                &                 wpin(i+3,j)*Mask(7,3)
                CTYPE(tmp2)(i,39)=wpin(i-3,j)*Mask(1,4)+wpin(i-2,j)*Mask(2,4)+wpin(i-1,j)*Mask(3,4)+ &
                &                 wpin(i,  j)*Mask(4,4)+wpin(i+1,j)*Mask(5,4)+wpin(i+2,j)*Mask(6,4)+ &
                &                 wpin(i+3,j)*Mask(7,4)
                CTYPE(tmp2)(i,40)=wpin(i-3,j)*Mask(1,5)+wpin(i-2,j)*Mask(2,5)+wpin(i-1,j)*Mask(3,5)+ &
                &                 wpin(i,  j)*Mask(4,5)+wpin(i+1,j)*Mask(5,5)+wpin(i+2,j)*Mask(6,5)+ &
                &                 wpin(i+3,j)*Mask(7,5)
                CTYPE(tmp2)(i,41)=wpin(i-3,j)*Mask(1,6)+wpin(i-2,j)*Mask(2,6)+wpin(i-1,j)*Mask(3,6)+ &
                &                 wpin(i,  j)*Mask(4,6)+wpin(i+1,j)*Mask(5,6)+wpin(i+2,j)*Mask(6,6)+ &
                &                 wpin(i+3,j)*Mask(7,6)
             ENDDO
             DO j=1,nnodes_(2)
                jp=j+3
                l=MOD((jp+2)*7,49)+1
                l2=l+1
                l3=l+2
                l4=l+3
                l5=l+4
                l6=l+5
                l7=l+6
                DO i=1,nnodes_(1)
                   CTYPE(tmp2)(i,l )=wpin(i-3,jp)*Mask(1,1)+wpin(i-2,jp)*Mask(2,1)+wpin(i-1,jp)*Mask(3,1)+ &
                   &                 wpin(i,  jp)*Mask(4,1)+wpin(i+1,jp)*Mask(5,1)+wpin(i+2,jp)*Mask(6,1)+ &
                   &                 wpin(i+3,jp)*Mask(7,1)
                   CTYPE(tmp2)(i,l2)=wpin(i-3,jp)*Mask(1,2)+wpin(i-2,jp)*Mask(2,2)+wpin(i-1,jp)*Mask(3,2)+ &
                   &                 wpin(i,  jp)*Mask(4,2)+wpin(i+1,jp)*Mask(5,2)+wpin(i+2,jp)*Mask(6,2)+ &
                   &                 wpin(i+3,jp)*Mask(7,2)
                   CTYPE(tmp2)(i,l3)=wpin(i-3,jp)*Mask(1,3)+wpin(i-2,jp)*Mask(2,3)+wpin(i-1,jp)*Mask(3,3)+ &
                   &                 wpin(i,  jp)*Mask(4,3)+wpin(i+1,jp)*Mask(5,3)+wpin(i+2,jp)*Mask(6,3)+ &
                   &                 wpin(i+3,jp)*Mask(7,3)
                   CTYPE(tmp2)(i,l4)=wpin(i-3,jp)*Mask(1,4)+wpin(i-2,jp)*Mask(2,4)+wpin(i-1,jp)*Mask(3,4)+ &
                   &                 wpin(i,  jp)*Mask(4,4)+wpin(i+1,jp)*Mask(5,4)+wpin(i+2,jp)*Mask(6,4)+ &
                   &                 wpin(i+3,jp)*Mask(7,4)
                   CTYPE(tmp2)(i,l5)=wpin(i-3,jp)*Mask(1,5)+wpin(i-2,jp)*Mask(2,5)+wpin(i-1,jp)*Mask(3,5)+ &
                   &                 wpin(i,  jp)*Mask(4,5)+wpin(i+1,jp)*Mask(5,5)+wpin(i+2,jp)*Mask(6,5)+ &
                   &                 wpin(i+3,jp)*Mask(7,5)
                   CTYPE(tmp2)(i,l6)=wpin(i-3,jp)*Mask(1,6)+wpin(i-2,jp)*Mask(2,6)+wpin(i-1,jp)*Mask(3,6)+ &
                   &                 wpin(i,  jp)*Mask(4,6)+wpin(i+1,jp)*Mask(5,6)+wpin(i+2,jp)*Mask(6,6)+ &
                   &                 wpin(i+3,jp)*Mask(7,6)
                   CTYPE(tmp2)(i,l7)=wpin(i-3,jp)*Mask(1,7)+wpin(i-2,jp)*Mask(2,7)+wpin(i-1,jp)*Mask(3,7)+ &
                   &                 wpin(i,  jp)*Mask(4,7)+wpin(i+1,jp)*Mask(5,7)+wpin(i+2,jp)*Mask(6,7)+ &
                   &                 wpin(i+3,jp)*Mask(7,7)
                ENDDO
                l1=MOD(l+ 7,49)
                l2=MOD(l+15,49)
                l3=MOD(l+23,49)
                l4=MOD(l+31,49)
                l5=MOD(l+39,49)
                l6=MOD(l+47,49)

                FORALL (i=1:nnodes_(1))
                   wpout(i,j)=CTYPE(tmp2)(i,l1)
                END FORALL
                FORALL (i=1:nnodes_(1))
                   wpout(i,j)=wpout(i,j)+CTYPE(tmp2)(i,l2)
                END FORALL
                FORALL (i=1:nnodes_(1))
                   wpout(i,j)=wpout(i,j)+CTYPE(tmp2)(i,l3)
                END FORALL
                FORALL (i=1:nnodes_(1))
                   wpout(i,j)=wpout(i,j)+CTYPE(tmp2)(i,l4)
                END FORALL
                FORALL (i=1:nnodes_(1))
                   wpout(i,j)=wpout(i,j)+CTYPE(tmp2)(i,l5)
                END FORALL
                FORALL (i=1:nnodes_(1))
                   wpout(i,j)=wpout(i,j)+CTYPE(tmp2)(i,l6)
                END FORALL
                FORALL (i=1:nnodes_(1))
                   wpout(i,j)=wpout(i,j)+CTYPE(tmp2)(i,l7)
                END FORALL
             ENDDO !j=1,nnodes_(2)

             DEALLOCATE(CTYPE(tmp2),STAT=info)
             or_fail_dealloc("tmp2",exit_point=7777)
          ELSE
             m=msize_(1)
             n=msize_(2)

             IF (n.EQ.0) THEN
                ALLOCATE(CTYPE(tmp1)(nnodes_(1)),STAT=info)
                or_fail_alloc("tmp1",exit_point=7777)

                DO j=1,nnodes_(2)
                   DO i=1,nnodes_(1)
                      CTYPE(tmp1)(i)=SUM(wpin(i-m:i+m,j)*Mask(:,1))
                   ENDDO

                   FORALL (i=1:nnodes_(1))
                      wpout(i,j)=CTYPE(tmp1)(i)
                   END FORALL
                ENDDO !j=1,nnodes_(2)

                DEALLOCATE(CTYPE(tmp1),STAT=info)
                or_fail_dealloc("tmp1",exit_point=7777)

                GOTO 7777
             ENDIF !n.EQ.0

             l1=2*m+1
             l2=2*n+1
             l3=l2*l2

             ALLOCATE(CTYPE(tmp2)(nnodes_(1),l3),STAT=info)
             or_fail_alloc("tmp2",exit_point=7777)

             l4=0
             DO j=1-n,n
                l5=l4*l2
                !counter
                l4=l4+1
                !counter
                DO l=1,l4
                   DO i=1,nnodes_(1)
                      CTYPE(tmp2)(i,l5+l)=SUM(wpin(i-m:i+m,j)*Mask(:,l))
                   ENDDO
                ENDDO
             ENDDO

             DO j=1,nnodes_(2)
                jp=j+n
                l=MOD((jp+n-1)*l2,l3)
                !Here l do not have +1 !
                DO l4=1,l1
                   l5=l+l4
                   DO i=1,nnodes_(1)
                      CTYPE(tmp2)(i,l5)=SUM(wpin(i-m:i+m,j)*Mask(:,l4))
                   ENDDO
                ENDDO

                l4=l2+1

                l5=1
                l6=MOD(l+l5*l4,l3)
                IF (l6.EQ.0) l6=l3
                FORALL (i=1:nnodes_(1))
                   wpout(i,j)=CTYPE(tmp2)(i,l6)
                END FORALL
                DO l5=2,l2
                   l6=MOD(l+l5*l4,l3)
                   IF (l6.EQ.0) l6=l3
                   FORALL (i=1:nnodes_(1))
                      wpout(i,j)=wpout(i,j)+CTYPE(tmp2)(i,l6)
                   END FORALL
                ENDDO
             ENDDO !j=1,nnodes_(2)

             DEALLOCATE(CTYPE(tmp2),STAT=info)
             or_fail_dealloc("tmp2",exit_point=7777)
          ENDIF
        7777 CONTINUE
        END SUBROUTINE CTYPE(Update)
#elif __DIME == __3D
        SUBROUTINE CTYPE(Update)(Mask,wpin,wpout,nnodes_,msize_,info)
          IMPLICIT NONE
          __TYPE,             DIMENSION(:,:,:), INTENT(IN   ) :: Mask
          __TYPE, CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: wpin
          __TYPE, CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: wpout
          INTEGER,            DIMENSION(:),     INTENT(IN   ) :: nnodes_
          INTEGER,            DIMENSION(:),     INTENT(IN   ) :: msize_
          INTEGER,                              INTENT(  OUT) :: info

          INTEGER :: i,j,k,l,m,n
          INTEGER :: jp,l1,l2,l3,l4,l5

          IF (ALL(msize_.EQ.1)) THEN
             !At each pixel, I should keep three values of the first, second
             !and third row mask multiply by the pixel row
             ALLOCATE(CTYPE(tmp2)(nnodes_(1),9),STAT=info)
             or_fail_alloc("tmp2",exit_point=7777)

             DO k=1,nnodes_(3)
                m=k-1
                j=0
                DO i=1,nnodes_(1)
                   CTYPE(tmp2)(i,1)=                 wpin(i-1,j,m)*Mask(1,1,1)+wpin(i,j,m)*Mask(2,1,1)+wpin(i+1,j,m)*Mask(3,1,1)
                ENDDO
                j=1
                DO i=1,nnodes_(1)
                   CTYPE(tmp2)(i,4)=                 wpin(i-1,j,m)*Mask(1,1,1)+wpin(i,j,m)*Mask(2,1,1)+wpin(i+1,j,m)*Mask(3,1,1)
                   CTYPE(tmp2)(i,5)=                 wpin(i-1,j,m)*Mask(1,2,1)+wpin(i,j,m)*Mask(2,2,1)+wpin(i+1,j,m)*Mask(3,2,1)
                ENDDO
                !!!
                j=0
                DO i=1,nnodes_(1)
                   CTYPE(tmp2)(i,1)=CTYPE(tmp2)(i,1)+wpin(i-1,j,k)*Mask(1,1,2)+wpin(i,j,k)*Mask(2,1,2)+wpin(i+1,j,k)*Mask(3,1,2)
                ENDDO
                j=1
                DO i=1,nnodes_(1)
                   CTYPE(tmp2)(i,4)=CTYPE(tmp2)(i,4)+wpin(i-1,j,k)*Mask(1,1,2)+wpin(i,j,k)*Mask(2,1,2)+wpin(i+1,j,k)*Mask(3,1,2)
                   CTYPE(tmp2)(i,5)=CTYPE(tmp2)(i,5)+wpin(i-1,j,k)*Mask(1,2,2)+wpin(i,j,k)*Mask(2,2,2)+wpin(i+1,j,k)*Mask(3,2,2)
                ENDDO
                !!!
                m=k+1
                j=0
                DO i=1,nnodes_(1)
                   CTYPE(tmp2)(i,1)=CTYPE(tmp2)(i,1)+wpin(i-1,j,m)*Mask(1,1,3)+wpin(i,j,m)*Mask(2,1,3)+wpin(i+1,j,m)*Mask(3,1,3)
                ENDDO
                j=1
                DO i=1,nnodes_(1)
                   CTYPE(tmp2)(i,4)=CTYPE(tmp2)(i,4)+wpin(i-1,j,m)*Mask(1,1,3)+wpin(i,j,m)*Mask(2,1,3)+wpin(i+1,j,m)*Mask(3,1,3)
                   CTYPE(tmp2)(i,5)=CTYPE(tmp2)(i,5)+wpin(i-1,j,m)*Mask(1,2,3)+wpin(i,j,m)*Mask(2,2,3)+wpin(i+1,j,m)*Mask(3,2,3)
                ENDDO

                DO j=1,nnodes_(2)
                   jp=j+1
                   l=MOD(jp*3,9)+1

                   m=k-1
                   l2=l+1
                   l3=l+2
                   DO i=1,nnodes_(1)
                      CTYPE(tmp2)(i,l )=                  wpin(i-1,j,m)*Mask(1,1,1)+wpin(i,j,m)*Mask(2,1,1)+wpin(i+1,j,m)*Mask(3,1,1)
                      CTYPE(tmp2)(i,l2)=                  wpin(i-1,j,m)*Mask(1,2,1)+wpin(i,j,m)*Mask(2,2,1)+wpin(i+1,j,m)*Mask(3,2,1)
                      CTYPE(tmp2)(i,l3)=                  wpin(i-1,j,m)*Mask(1,3,1)+wpin(i,j,m)*Mask(2,3,1)+wpin(i+1,j,m)*Mask(3,3,1)
                   ENDDO
                   DO i=1,nnodes_(1)
                      CTYPE(tmp2)(i,l )=CTYPE(tmp2)(i,l )+wpin(i-1,j,k)*Mask(1,1,2)+wpin(i,j,k)*Mask(2,1,2)+wpin(i+1,j,k)*Mask(3,1,2)
                      CTYPE(tmp2)(i,l2)=CTYPE(tmp2)(i,l2)+wpin(i-1,j,k)*Mask(1,2,2)+wpin(i,j,k)*Mask(2,2,2)+wpin(i+1,j,k)*Mask(3,2,2)
                      CTYPE(tmp2)(i,l3)=CTYPE(tmp2)(i,l3)+wpin(i-1,j,k)*Mask(1,3,2)+wpin(i,j,k)*Mask(2,3,2)+wpin(i+1,j,k)*Mask(3,3,2)
                   ENDDO
                   m=k+1
                   DO i=1,nnodes_(1)
                      CTYPE(tmp2)(i,l )=CTYPE(tmp2)(i,l )+wpin(i-1,j,m)*Mask(1,1,3)+wpin(i,j,m)*Mask(2,1,3)+wpin(i+1,j,m)*Mask(3,1,3)
                      CTYPE(tmp2)(i,l1)=CTYPE(tmp2)(i,l2)+wpin(i-1,j,m)*Mask(1,2,3)+wpin(i,j,m)*Mask(2,2,3)+wpin(i+1,j,m)*Mask(3,2,3)
                      CTYPE(tmp2)(i,l2)=CTYPE(tmp2)(i,l3)+wpin(i-1,j,m)*Mask(1,3,3)+wpin(i,j,m)*Mask(2,3,3)+wpin(i+1,j,m)*Mask(3,3,3)
                   ENDDO
                   !index l-6
                   l1=MOD(l+3,9)
                   !index l-2
                   l2=MOD(l+7,9)

                   FORALL (i=1:nnodes_(1))
                      wpout(i,j,k)=CTYPE(tmp2)(i,l1)
                   END FORALL
                   FORALL (i=1:nnodes_(1))
                      wpout(i,j,k)=wpout(i,j,k)+CTYPE(tmp2)(i,l2)
                   END FORALL
                   FORALL (i=1:nnodes_(1))
                      wpout(i,j,k)=wpout(i,j,k)+CTYPE(tmp2)(i,l3)
                   END FORALL
                ENDDO !j=1,nnodes_(2)
             ENDDO !k

             DEALLOCATE(CTYPE(tmp2),STAT=info)
             or_fail_dealloc("tmp2",exit_point=7777)
          ELSE IF (ALL(msize_.EQ.2)) THEN
             !At each pixel, I should keep three values of the first, second
             !and third row mask multiply by the pixel row
             ALLOCATE(CTYPE(tmp2)(nnodes_(1),25),STAT=info)
             or_fail_alloc("tmp2",exit_point=7777)
             !!!
             DO k=1,nnodes_(3)
                !!!
                m=k-2
                j=-1
                DO i=1,nnodes_(1)
                   CTYPE(tmp2)(i,1)=                 wpin(i-2,j,m)*Mask(1,1,1)+wpin(i-1,j,m)*Mask(2,1,1)+wpin(i,j,m)*Mask(3,1,1)+wpin(i+1,j,m)*Mask(4,1,1)+wpin(i+2,j,m)*Mask(5,1,1)
                ENDDO
                j=0
                DO i=1,nnodes_(1)
                   CTYPE(tmp2)(i,6)=                 wpin(i-2,j,m)*Mask(1,1,1)+wpin(i-1,j,m)*Mask(2,1,1)+wpin(i,j,m)*Mask(3,1,1)+wpin(i+1,j,m)*Mask(4,1,1)+wpin(i+2,j,m)*Mask(5,1,1)
                   CTYPE(tmp2)(i,7)=                 wpin(i-2,j,m)*Mask(1,2,1)+wpin(i-1,j,m)*Mask(2,2,1)+wpin(i,j,m)*Mask(3,2,1)+wpin(i+1,j,m)*Mask(4,2,1)+wpin(i+2,j,m)*Mask(5,2,1)
                ENDDO
                j=1
                DO i=1,nnodes_(1)
                   CTYPE(tmp2)(i,11)=                  wpin(i-2,j,m)*Mask(1,1,1)+wpin(i-1,j,m)*Mask(2,1,1)+wpin(i,j,m)*Mask(3,1,1)+wpin(i+1,j,m)*Mask(4,1,1)+wpin(i+2,j,m)*Mask(5,1,1)
                   CTYPE(tmp2)(i,12)=                  wpin(i-2,j,m)*Mask(1,2,1)+wpin(i-1,j,m)*Mask(2,2,1)+wpin(i,j,m)*Mask(3,2,1)+wpin(i+1,j,m)*Mask(4,2,1)+wpin(i+2,j,m)*Mask(5,2,1)
                   CTYPE(tmp2)(i,13)=                  wpin(i-2,j,m)*Mask(1,3,1)+wpin(i-1,j,m)*Mask(2,3,1)+wpin(i,j,m)*Mask(3,3,1)+wpin(i+1,j,m)*Mask(4,3,1)+wpin(i+2,j,m)*Mask(5,3,1)
                ENDDO
                j=2
                DO i=1,nnodes_(1)
                   CTYPE(tmp2)(i,16)=                  wpin(i-2,j,m)*Mask(1,1,1)+wpin(i-1,j,m)*Mask(2,1,1)+wpin(i,j,m)*Mask(3,1,1)+wpin(i+1,j,m)*Mask(4,1,1)+wpin(i+2,j,m)*Mask(5,1,1)
                   CTYPE(tmp2)(i,17)=                  wpin(i-2,j,m)*Mask(1,2,1)+wpin(i-1,j,m)*Mask(2,2,1)+wpin(i,j,m)*Mask(3,2,1)+wpin(i+1,j,m)*Mask(4,2,1)+wpin(i+2,j,m)*Mask(5,2,1)
                   CTYPE(tmp2)(i,18)=                  wpin(i-2,j,m)*Mask(1,3,1)+wpin(i-1,j,m)*Mask(2,3,1)+wpin(i,j,m)*Mask(3,3,1)+wpin(i+1,j,m)*Mask(4,3,1)+wpin(i+2,j,m)*Mask(5,3,1)
                   CTYPE(tmp2)(i,19)=                  wpin(i-2,j,m)*Mask(1,4,1)+wpin(i-1,j,m)*Mask(2,4,1)+wpin(i,j,m)*Mask(3,4,1)+wpin(i+1,j,m)*Mask(4,4,1)+wpin(i+2,j,m)*Mask(5,4,1)
                ENDDO
                !!!
                m=k-1
                j=-1
                DO i=1,nnodes_(1)
                   CTYPE(tmp2)(i,1)=CTYPE(tmp2)(i,1)+wpin(i-2,j,m)*Mask(1,1,2)+wpin(i-1,j,m)*Mask(2,1,2)+wpin(i,j,m)*Mask(3,1,2)+wpin(i+1,j,m)*Mask(4,1,2)+wpin(i+2,j,m)*Mask(5,1,2)
                ENDDO
                j=0
                DO i=1,nnodes_(1)
                   CTYPE(tmp2)(i,6)=CTYPE(tmp2)(i,6)+wpin(i-2,j,m)*Mask(1,1,2)+wpin(i-1,j,m)*Mask(2,1,2)+wpin(i,j,m)*Mask(3,1,2)+wpin(i+1,j,m)*Mask(4,1,2)+wpin(i+2,j,m)*Mask(5,1,2)
                   CTYPE(tmp2)(i,7)=CTYPE(tmp2)(i,7)+wpin(i-2,j,m)*Mask(1,2,2)+wpin(i-1,j,m)*Mask(2,2,2)+wpin(i,j,m)*Mask(3,2,2)+wpin(i+1,j,m)*Mask(4,2,2)+wpin(i+2,j,m)*Mask(5,2,2)
                ENDDO
                j=1
                DO i=1,nnodes_(1)
                   CTYPE(tmp2)(i,11)=CTYPE(tmp2)(i,11)+wpin(i-2,j,m)*Mask(1,1,2)+wpin(i-1,j,m)*Mask(2,1,2)+wpin(i,j,m)*Mask(3,1,2)+wpin(i+1,j,m)*Mask(4,1,2)+wpin(i+2,j,m)*Mask(5,1,2)
                   CTYPE(tmp2)(i,12)=CTYPE(tmp2)(i,12)+wpin(i-2,j,m)*Mask(1,2,2)+wpin(i-1,j,m)*Mask(2,2,2)+wpin(i,j,m)*Mask(3,2,2)+wpin(i+1,j,m)*Mask(4,2,2)+wpin(i+2,j,m)*Mask(5,2,2)
                   CTYPE(tmp2)(i,13)=CTYPE(tmp2)(i,13)+wpin(i-2,j,m)*Mask(1,3,2)+wpin(i-1,j,m)*Mask(2,3,2)+wpin(i,j,m)*Mask(3,3,2)+wpin(i+1,j,m)*Mask(4,3,2)+wpin(i+2,j,m)*Mask(5,3,2)
                ENDDO
                j=2
                DO i=1,nnodes_(1)
                   CTYPE(tmp2)(i,16)=CTYPE(tmp2)(i,16)+wpin(i-2,j,m)*Mask(1,1,2)+wpin(i-1,j,m)*Mask(2,1,2)+wpin(i,j,m)*Mask(3,1,2)+wpin(i+1,j,m)*Mask(4,1,2)+wpin(i+2,j,m)*Mask(5,1,2)
                   CTYPE(tmp2)(i,17)=CTYPE(tmp2)(i,17)+wpin(i-2,j,m)*Mask(1,2,2)+wpin(i-1,j,m)*Mask(2,2,2)+wpin(i,j,m)*Mask(3,2,2)+wpin(i+1,j,m)*Mask(4,2,2)+wpin(i+2,j,m)*Mask(5,2,2)
                   CTYPE(tmp2)(i,18)=CTYPE(tmp2)(i,18)+wpin(i-2,j,m)*Mask(1,3,2)+wpin(i-1,j,m)*Mask(2,3,2)+wpin(i,j,m)*Mask(3,3,2)+wpin(i+1,j,m)*Mask(4,3,2)+wpin(i+2,j,m)*Mask(5,3,2)
                   CTYPE(tmp2)(i,19)=CTYPE(tmp2)(i,19)+wpin(i-2,j,m)*Mask(1,4,2)+wpin(i-1,j,m)*Mask(2,4,2)+wpin(i,j,m)*Mask(3,4,2)+wpin(i+1,j,m)*Mask(4,4,2)+wpin(i+2,j,m)*Mask(5,4,2)
                ENDDO
                !!!
                j=-1
                DO i=1,nnodes_(1)
                   CTYPE(tmp2)(i,1)=CTYPE(tmp2)(i,1)+wpin(i-2,j,k)*Mask(1,1,3)+wpin(i-1,j,k)*Mask(2,1,3)+wpin(i,j,k)*Mask(3,1,3)+wpin(i+1,j,k)*Mask(4,1,3)+wpin(i+2,j,k)*Mask(5,1,3)
                ENDDO
                j=0
                DO i=1,nnodes_(1)
                   CTYPE(tmp2)(i,6)=CTYPE(tmp2)(i,6)+wpin(i-2,j,k)*Mask(1,1,3)+wpin(i-1,j,k)*Mask(2,1,3)+wpin(i,j,k)*Mask(3,1,3)+wpin(i+1,j,k)*Mask(4,1,3)+wpin(i+2,j,k)*Mask(5,1,3)
                   CTYPE(tmp2)(i,7)=CTYPE(tmp2)(i,7)+wpin(i-2,j,k)*Mask(1,2,3)+wpin(i-1,j,k)*Mask(2,2,3)+wpin(i,j,k)*Mask(3,2,3)+wpin(i+1,j,k)*Mask(4,2,3)+wpin(i+2,j,k)*Mask(5,2,3)
                ENDDO
                j=1
                DO i=1,nnodes_(1)
                   CTYPE(tmp2)(i,11)=CTYPE(tmp2)(i,11)+wpin(i-2,j,k)*Mask(1,1,3)+wpin(i-1,j,k)*Mask(2,1,3)+wpin(i,j,k)*Mask(3,1,3)+wpin(i+1,j,k)*Mask(4,1,3)+wpin(i+2,j,k)*Mask(5,1,3)
                   CTYPE(tmp2)(i,12)=CTYPE(tmp2)(i,12)+wpin(i-2,j,k)*Mask(1,2,3)+wpin(i-1,j,k)*Mask(2,2,3)+wpin(i,j,k)*Mask(3,2,3)+wpin(i+1,j,k)*Mask(4,2,3)+wpin(i+2,j,k)*Mask(5,2,3)
                   CTYPE(tmp2)(i,13)=CTYPE(tmp2)(i,13)+wpin(i-2,j,k)*Mask(1,3,3)+wpin(i-1,j,k)*Mask(2,3,3)+wpin(i,j,k)*Mask(3,3,3)+wpin(i+1,j,k)*Mask(4,3,3)+wpin(i+2,j,k)*Mask(5,3,3)
                ENDDO
                j=2
                DO i=1,nnodes_(1)
                   CTYPE(tmp2)(i,16)=CTYPE(tmp2)(i,16)+wpin(i-2,j,k)*Mask(1,1,3)+wpin(i-1,j,k)*Mask(2,1,3)+wpin(i,j,k)*Mask(3,1,3)+wpin(i+1,j,k)*Mask(4,1,3)+wpin(i+2,j,k)*Mask(5,1,3)
                   CTYPE(tmp2)(i,17)=CTYPE(tmp2)(i,17)+wpin(i-2,j,k)*Mask(1,2,3)+wpin(i-1,j,k)*Mask(2,2,3)+wpin(i,j,k)*Mask(3,2,3)+wpin(i+1,j,k)*Mask(4,2,3)+wpin(i+2,j,k)*Mask(5,2,3)
                   CTYPE(tmp2)(i,18)=CTYPE(tmp2)(i,18)+wpin(i-2,j,k)*Mask(1,3,3)+wpin(i-1,j,k)*Mask(2,3,3)+wpin(i,j,k)*Mask(3,3,3)+wpin(i+1,j,k)*Mask(4,3,3)+wpin(i+2,j,k)*Mask(5,3,3)
                   CTYPE(tmp2)(i,19)=CTYPE(tmp2)(i,19)+wpin(i-2,j,k)*Mask(1,4,3)+wpin(i-1,j,k)*Mask(2,4,3)+wpin(i,j,k)*Mask(3,4,3)+wpin(i+1,j,k)*Mask(4,4,3)+wpin(i+2,j,k)*Mask(5,4,3)
                ENDDO
                !!!
                m=k+1
                j=-1
                DO i=1,nnodes_(1)
                   CTYPE(tmp2)(i,1)=CTYPE(tmp2)(i,1)+wpin(i-2,j,m)*Mask(1,1,4)+wpin(i-1,j,m)*Mask(2,1,4)+wpin(i,j,m)*Mask(3,1,4)+wpin(i+1,j,m)*Mask(4,1,4)+wpin(i+2,j,m)*Mask(5,1,4)
                ENDDO
                j=0
                DO i=1,nnodes_(1)
                   CTYPE(tmp2)(i,6)=CTYPE(tmp2)(i,6)+wpin(i-2,j,m)*Mask(1,1,4)+wpin(i-1,j,m)*Mask(2,1,4)+wpin(i,j,m)*Mask(3,1,4)+wpin(i+1,j,m)*Mask(4,1,4)+wpin(i+2,j,m)*Mask(5,1,4)
                   CTYPE(tmp2)(i,7)=CTYPE(tmp2)(i,7)+wpin(i-2,j,m)*Mask(1,2,4)+wpin(i-1,j,m)*Mask(2,2,4)+wpin(i,j,m)*Mask(3,2,4)+wpin(i+1,j,m)*Mask(4,2,4)+wpin(i+2,j,m)*Mask(5,2,4)
                ENDDO
                j=1
                DO i=1,nnodes_(1)
                   CTYPE(tmp2)(i,11)=CTYPE(tmp2)(i,11)+wpin(i-2,j,m)*Mask(1,1,4)+wpin(i-1,j,m)*Mask(2,1,4)+wpin(i,j,m)*Mask(3,1,4)+wpin(i+1,j,m)*Mask(4,1,4)+wpin(i+2,j,m)*Mask(5,1,4)
                   CTYPE(tmp2)(i,12)=CTYPE(tmp2)(i,12)+wpin(i-2,j,m)*Mask(1,2,4)+wpin(i-1,j,m)*Mask(2,2,4)+wpin(i,j,m)*Mask(3,2,4)+wpin(i+1,j,m)*Mask(4,2,4)+wpin(i+2,j,m)*Mask(5,2,4)
                   CTYPE(tmp2)(i,13)=CTYPE(tmp2)(i,13)+wpin(i-2,j,m)*Mask(1,3,4)+wpin(i-1,j,m)*Mask(2,3,4)+wpin(i,j,m)*Mask(3,3,4)+wpin(i+1,j,m)*Mask(4,3,4)+wpin(i+2,j,m)*Mask(5,3,4)
                ENDDO
                j=2
                DO i=1,nnodes_(1)
                   CTYPE(tmp2)(i,16)=CTYPE(tmp2)(i,16)+wpin(i-2,j,m)*Mask(1,1,4)+wpin(i-1,j,m)*Mask(2,1,4)+wpin(i,j,m)*Mask(3,1,4)+wpin(i+1,j,m)*Mask(4,1,4)+wpin(i+2,j,m)*Mask(5,1,4)
                   CTYPE(tmp2)(i,17)=CTYPE(tmp2)(i,17)+wpin(i-2,j,m)*Mask(1,2,4)+wpin(i-1,j,m)*Mask(2,2,4)+wpin(i,j,m)*Mask(3,2,4)+wpin(i+1,j,m)*Mask(4,2,4)+wpin(i+2,j,m)*Mask(5,2,4)
                   CTYPE(tmp2)(i,18)=CTYPE(tmp2)(i,18)+wpin(i-2,j,m)*Mask(1,3,4)+wpin(i-1,j,m)*Mask(2,3,4)+wpin(i,j,m)*Mask(3,3,4)+wpin(i+1,j,m)*Mask(4,3,4)+wpin(i+2,j,m)*Mask(5,3,4)
                   CTYPE(tmp2)(i,19)=CTYPE(tmp2)(i,19)+wpin(i-2,j,m)*Mask(1,4,4)+wpin(i-1,j,m)*Mask(2,4,4)+wpin(i,j,m)*Mask(3,4,4)+wpin(i+1,j,m)*Mask(4,4,4)+wpin(i+2,j,m)*Mask(5,4,4)
                ENDDO
                !!!
                m=k+2
                j=-1
                DO i=1,nnodes_(1)
                   CTYPE(tmp2)(i,1)=CTYPE(tmp2)(i,1)+wpin(i-2,j,m)*Mask(1,1,5)+wpin(i-1,j,m)*Mask(2,1,5)+wpin(i,j,m)*Mask(3,1,5)+wpin(i+1,j,m)*Mask(4,1,5)+wpin(i+2,j,m)*Mask(5,1,5)
                ENDDO
                j=0
                DO i=1,nnodes_(1)
                   CTYPE(tmp2)(i,6)=CTYPE(tmp2)(i,6)+wpin(i-2,j,m)*Mask(1,1,5)+wpin(i-1,j,m)*Mask(2,1,5)+wpin(i,j,m)*Mask(3,1,5)+wpin(i+1,j,m)*Mask(4,1,5)+wpin(i+2,j,m)*Mask(5,1,5)
                   CTYPE(tmp2)(i,7)=CTYPE(tmp2)(i,7)+wpin(i-2,j,m)*Mask(1,2,5)+wpin(i-1,j,m)*Mask(2,2,5)+wpin(i,j,m)*Mask(3,2,5)+wpin(i+1,j,m)*Mask(4,2,5)+wpin(i+2,j,m)*Mask(5,2,5)
                ENDDO
                j=1
                DO i=1,nnodes_(1)
                   CTYPE(tmp2)(i,11)=CTYPE(tmp2)(i,11)+wpin(i-2,j,m)*Mask(1,1,5)+wpin(i-1,j,m)*Mask(2,1,5)+wpin(i,j,m)*Mask(3,1,5)+wpin(i+1,j,m)*Mask(4,1,5)+wpin(i+2,j,m)*Mask(5,1,5)
                   CTYPE(tmp2)(i,12)=CTYPE(tmp2)(i,12)+wpin(i-2,j,m)*Mask(1,2,5)+wpin(i-1,j,m)*Mask(2,2,5)+wpin(i,j,m)*Mask(3,2,5)+wpin(i+1,j,m)*Mask(4,2,5)+wpin(i+2,j,m)*Mask(5,2,5)
                   CTYPE(tmp2)(i,13)=CTYPE(tmp2)(i,13)+wpin(i-2,j,m)*Mask(1,3,5)+wpin(i-1,j,m)*Mask(2,3,5)+wpin(i,j,m)*Mask(3,3,5)+wpin(i+1,j,m)*Mask(4,3,5)+wpin(i+2,j,m)*Mask(5,3,5)
                ENDDO
                j=2
                DO i=1,nnodes_(1)
                   CTYPE(tmp2)(i,16)=CTYPE(tmp2)(i,16)+wpin(i-2,j,m)*Mask(1,1,5)+wpin(i-1,j,m)*Mask(2,1,5)+wpin(i,j,m)*Mask(3,1,5)+wpin(i+1,j,m)*Mask(4,1,5)+wpin(i+2,j,m)*Mask(5,1,5)
                   CTYPE(tmp2)(i,17)=CTYPE(tmp2)(i,17)+wpin(i-2,j,m)*Mask(1,2,5)+wpin(i-1,j,m)*Mask(2,2,5)+wpin(i,j,m)*Mask(3,2,5)+wpin(i+1,j,m)*Mask(4,2,5)+wpin(i+2,j,m)*Mask(5,2,5)
                   CTYPE(tmp2)(i,18)=CTYPE(tmp2)(i,18)+wpin(i-2,j,m)*Mask(1,3,5)+wpin(i-1,j,m)*Mask(2,3,5)+wpin(i,j,m)*Mask(3,3,5)+wpin(i+1,j,m)*Mask(4,3,5)+wpin(i+2,j,m)*Mask(5,3,5)
                   CTYPE(tmp2)(i,19)=CTYPE(tmp2)(i,19)+wpin(i-2,j,m)*Mask(1,4,5)+wpin(i-1,j,m)*Mask(2,4,5)+wpin(i,j,m)*Mask(3,4,5)+wpin(i+1,j,m)*Mask(4,4,5)+wpin(i+2,j,m)*Mask(5,4,5)
                ENDDO

                DO j=1,nnodes_(2)
                   jp=j+2
                   l=MOD((jp+1)*5,25)+1

                   m=k-2
                   l2=l+1
                   l3=l+2
                   l4=l+3
                   l5=l+4
                   DO i=1,nnodes_(1)
                      CTYPE(tmp2)(i,l )=                  wpin(i-2,j,m)*Mask(1,1,1)+wpin(i-1,j,m)*Mask(2,1,1)+wpin(i,j,m)*Mask(3,1,1)+wpin(i+1,j,m)*Mask(4,1,1)+wpin(i+2,j,m)*Mask(5,1,1)
                      CTYPE(tmp2)(i,l2)=                  wpin(i-2,j,m)*Mask(1,2,1)+wpin(i-1,j,m)*Mask(2,2,1)+wpin(i,j,m)*Mask(3,2,1)+wpin(i+1,j,m)*Mask(4,2,1)+wpin(i+2,j,m)*Mask(5,2,1)
                      CTYPE(tmp2)(i,l3)=                  wpin(i-2,j,m)*Mask(1,3,1)+wpin(i-1,j,m)*Mask(2,3,1)+wpin(i,j,m)*Mask(3,3,1)+wpin(i+1,j,m)*Mask(4,3,1)+wpin(i+2,j,m)*Mask(5,3,1)
                      CTYPE(tmp2)(i,l4)=                  wpin(i-2,j,m)*Mask(1,4,1)+wpin(i-1,j,m)*Mask(2,4,1)+wpin(i,j,m)*Mask(3,4,1)+wpin(i+1,j,m)*Mask(4,4,1)+wpin(i+2,j,m)*Mask(5,4,1)
                      CTYPE(tmp2)(i,l5)=                  wpin(i-2,j,m)*Mask(1,5,1)+wpin(i-1,j,m)*Mask(2,5,1)+wpin(i,j,m)*Mask(3,5,1)+wpin(i+1,j,m)*Mask(4,5,1)+wpin(i+2,j,m)*Mask(5,5,1)
                   ENDDO
                   m=k-1
                   DO i=1,nnodes_(1)
                      CTYPE(tmp2)(i,l )=CTYPE(tmp2)(i,l )+wpin(i-2,j,m)*Mask(1,1,2)+wpin(i-1,j,m)*Mask(2,1,2)+wpin(i,j,m)*Mask(3,1,2)+wpin(i+1,j,m)*Mask(4,1,2)+wpin(i+2,j,m)*Mask(5,1,2)
                      CTYPE(tmp2)(i,l2)=CTYPE(tmp2)(i,l2)+wpin(i-2,j,m)*Mask(1,2,2)+wpin(i-1,j,m)*Mask(2,2,2)+wpin(i,j,m)*Mask(3,2,2)+wpin(i+1,j,m)*Mask(4,2,2)+wpin(i+2,j,m)*Mask(5,2,2)
                      CTYPE(tmp2)(i,l3)=CTYPE(tmp2)(i,l3)+wpin(i-2,j,m)*Mask(1,3,2)+wpin(i-1,j,m)*Mask(2,3,2)+wpin(i,j,m)*Mask(3,3,2)+wpin(i+1,j,m)*Mask(4,3,2)+wpin(i+2,j,m)*Mask(5,3,2)
                      CTYPE(tmp2)(i,l4)=CTYPE(tmp2)(i,l4)+wpin(i-2,j,m)*Mask(1,4,2)+wpin(i-1,j,m)*Mask(2,4,2)+wpin(i,j,m)*Mask(3,4,2)+wpin(i+1,j,m)*Mask(4,4,2)+wpin(i+2,j,m)*Mask(5,4,2)
                      CTYPE(tmp2)(i,l5)=CTYPE(tmp2)(i,l5)+wpin(i-2,j,m)*Mask(1,5,2)+wpin(i-1,j,m)*Mask(2,5,2)+wpin(i,j,m)*Mask(3,5,2)+wpin(i+1,j,m)*Mask(4,5,2)+wpin(i+2,j,m)*Mask(5,5,2)
                   ENDDO
                   m=k
                   DO i=1,nnodes_(1)
                      CTYPE(tmp2)(i,l )=CTYPE(tmp2)(i,l )+wpin(i-2,j,m)*Mask(1,1,3)+wpin(i-1,j,m)*Mask(2,1,3)+wpin(i,j,m)*Mask(3,1,3)+wpin(i+1,j,m)*Mask(4,1,3)+wpin(i+2,j,m)*Mask(5,1,3)
                      CTYPE(tmp2)(i,l2)=CTYPE(tmp2)(i,l2)+wpin(i-2,j,m)*Mask(1,2,3)+wpin(i-1,j,m)*Mask(2,2,3)+wpin(i,j,m)*Mask(3,2,3)+wpin(i+1,j,m)*Mask(4,2,3)+wpin(i+2,j,m)*Mask(5,2,3)
                      CTYPE(tmp2)(i,l3)=CTYPE(tmp2)(i,l3)+wpin(i-2,j,m)*Mask(1,3,3)+wpin(i-1,j,m)*Mask(2,3,3)+wpin(i,j,m)*Mask(3,3,3)+wpin(i+1,j,m)*Mask(4,3,3)+wpin(i+2,j,m)*Mask(5,3,3)
                      CTYPE(tmp2)(i,l4)=CTYPE(tmp2)(i,l4)+wpin(i-2,j,m)*Mask(1,4,3)+wpin(i-1,j,m)*Mask(2,4,3)+wpin(i,j,m)*Mask(3,4,3)+wpin(i+1,j,m)*Mask(4,4,3)+wpin(i+2,j,m)*Mask(5,4,3)
                      CTYPE(tmp2)(i,l5)=CTYPE(tmp2)(i,l5)+wpin(i-2,j,m)*Mask(1,5,3)+wpin(i-1,j,m)*Mask(2,5,3)+wpin(i,j,m)*Mask(3,5,3)+wpin(i+1,j,m)*Mask(4,5,3)+wpin(i+2,j,m)*Mask(5,5,3)
                   ENDDO
                   m=k+1
                   DO i=1,nnodes_(1)
                      CTYPE(tmp2)(i,l )=CTYPE(tmp2)(i,l )+wpin(i-2,j,m)*Mask(1,1,4)+wpin(i-1,j,m)*Mask(2,1,4)+wpin(i,j,m)*Mask(3,1,4)+wpin(i+1,j,m)*Mask(4,1,4)+wpin(i+2,j,m)*Mask(5,1,4)
                      CTYPE(tmp2)(i,l2)=CTYPE(tmp2)(i,l2)+wpin(i-2,j,m)*Mask(1,2,4)+wpin(i-1,j,m)*Mask(2,2,4)+wpin(i,j,m)*Mask(3,2,4)+wpin(i+1,j,m)*Mask(4,2,4)+wpin(i+2,j,m)*Mask(5,2,4)
                      CTYPE(tmp2)(i,l3)=CTYPE(tmp2)(i,l3)+wpin(i-2,j,m)*Mask(1,3,4)+wpin(i-1,j,m)*Mask(2,3,4)+wpin(i,j,m)*Mask(3,3,4)+wpin(i+1,j,m)*Mask(4,3,4)+wpin(i+2,j,m)*Mask(5,3,4)
                      CTYPE(tmp2)(i,l4)=CTYPE(tmp2)(i,l4)+wpin(i-2,j,m)*Mask(1,4,4)+wpin(i-1,j,m)*Mask(2,4,4)+wpin(i,j,m)*Mask(3,4,4)+wpin(i+1,j,m)*Mask(4,4,4)+wpin(i+2,j,m)*Mask(5,4,4)
                      CTYPE(tmp2)(i,l5)=CTYPE(tmp2)(i,l5)+wpin(i-2,j,m)*Mask(1,5,4)+wpin(i-1,j,m)*Mask(2,5,4)+wpin(i,j,m)*Mask(3,5,4)+wpin(i+1,j,m)*Mask(4,5,4)+wpin(i+2,j,m)*Mask(5,5,4)
                   ENDDO
                   m=k+2
                   DO i=1,nnodes_(1)
                      CTYPE(tmp2)(i,l )=CTYPE(tmp2)(i,l )+wpin(i-2,j,m)*Mask(1,1,5)+wpin(i-1,j,m)*Mask(2,1,5)+wpin(i,j,m)*Mask(3,1,5)+wpin(i+1,j,m)*Mask(4,1,5)+wpin(i+2,j,m)*Mask(5,1,5)
                      CTYPE(tmp2)(i,l2)=CTYPE(tmp2)(i,l2)+wpin(i-2,j,m)*Mask(1,2,5)+wpin(i-1,j,m)*Mask(2,2,5)+wpin(i,j,m)*Mask(3,2,5)+wpin(i+1,j,m)*Mask(4,2,5)+wpin(i+2,j,m)*Mask(5,2,5)
                      CTYPE(tmp2)(i,l3)=CTYPE(tmp2)(i,l3)+wpin(i-2,j,m)*Mask(1,3,5)+wpin(i-1,j,m)*Mask(2,3,5)+wpin(i,j,m)*Mask(3,3,5)+wpin(i+1,j,m)*Mask(4,3,5)+wpin(i+2,j,m)*Mask(5,3,5)
                      CTYPE(tmp2)(i,l4)=CTYPE(tmp2)(i,l4)+wpin(i-2,j,m)*Mask(1,4,5)+wpin(i-1,j,m)*Mask(2,4,5)+wpin(i,j,m)*Mask(3,4,5)+wpin(i+1,j,m)*Mask(4,4,5)+wpin(i+2,j,m)*Mask(5,4,5)
                      CTYPE(tmp2)(i,l5)=CTYPE(tmp2)(i,l5)+wpin(i-2,j,m)*Mask(1,5,5)+wpin(i-1,j,m)*Mask(2,5,5)+wpin(i,j,m)*Mask(3,5,5)+wpin(i+1,j,m)*Mask(4,5,5)+wpin(i+2,j,m)*Mask(5,5,5)
                   ENDDO

                   l1=MOD(l+ 5,25)
                   l2=MOD(l+11,25)
                   l3=MOD(l+17,25)
                   l4=MOD(l+23,25)

                   FORALL (i=1:nnodes_(1))
                      wpout(i,j,k)=             CTYPE(tmp2)(i,l1)
                   END FORALL
                   FORALL (i=1:nnodes_(1))
                      wpout(i,j,k)=wpout(i,j,k)+CTYPE(tmp2)(i,l2)
                   END FORALL
                   FORALL (i=1:nnodes_(1))
                      wpout(i,j,k)=wpout(i,j,k)+CTYPE(tmp2)(i,l3)
                   END FORALL
                   FORALL (i=1:nnodes_(1))
                      wpout(i,j,k)=wpout(i,j,k)+CTYPE(tmp2)(i,l4)
                   END FORALL
                   FORALL (i=1:nnodes_(1))
                      wpout(i,j,k)=wpout(i,j,k)+CTYPE(tmp2)(i,l5)
                   END FORALL
                ENDDO !j=1,nnodes_(2)
             ENDDO !k

             DEALLOCATE(CTYPE(tmp2),STAT=info)
             or_fail_dealloc("tmp2",exit_point=7777)
          ELSE
             l=msize_(1)
             m=msize_(2)
             n=msize_(3)

             IF (m.EQ.0.AND.n.EQ.0) THEN
                ALLOCATE(CTYPE(tmp1)(1:nnodes_(1)),STAT=info)
                or_fail_alloc("tmp1",exit_point=7777)

                l=2*m+1

                DO k=1,nnodes_(3)
                   DO j=1,nnodes_(2)
                      DO i=1,nnodes_(1)
                         CTYPE(tmp1)(i)=SUM(wpin(i-m:i+m,j,k)*Mask(1:l,1,1))
                      ENDDO

                      DO i=1,nnodes_(1)
                         wpout(i,j,k)=CTYPE(tmp1)(i)
                      ENDDO
                   ENDDO !j=1,nnodes_(2)
                ENDDO !k=1,nnodes_(3)

                DEALLOCATE(CTYPE(tmp1),STAT=info)
                or_fail_dealloc("tmp1",exit_point=7777)

                GOTO 7777
             ENDIF !n.EQ.0

             l=msize_(1)
             m=msize_(2)
             n=msize_(3)

             IF (ASSOCIATED(wpin,wpout)) THEN
                ALLOCATE(CTYPE(tmp3)(1:nnodes_(1),1:nnodes_(2),1:nnodes_(3)),STAT=info)
                or_fail_alloc("tmp3",exit_point=7777)


                FORALL (i=1:nnodes_(1),j=1:nnodes_(2),k=1:nnodes_(3))
                   CTYPE(tmp3)(i,j,k)=SUM(wpin(i-l:i+l,j-m:j+m,k-n:k+n)*Mask(1:2*l+1,1:2*m+1,1:2*n+1))
                END FORALL

                FORALL (i=1:nnodes_(1),j=1:nnodes_(2),k=1:nnodes_(3)) wpout(i,j,k)=CTYPE(tmp3)(i,j,k)

                DEALLOCATE(CTYPE(tmp3),STAT=info)
                or_fail_dealloc("tmp3",exit_point=7777)
              ELSE
                FORALL (i=1:nnodes_(1),j=1:nnodes_(2),k=1:nnodes_(3))
                   wpout(i,j,k)=SUM(wpin(i-l:i+l,j-m:j+m,k-n:k+n)*Mask(1:2*l+1,1:2*m+1,1:2*n+1))
                END FORALL
              ENDIF

          ENDIF
        7777 CONTINUE
        END SUBROUTINE CTYPE(Update)
#endif
        END SUBROUTINE DTYPE(CTYPE(ppm_rc_convolve))

