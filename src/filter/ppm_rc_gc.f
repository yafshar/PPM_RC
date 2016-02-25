        !-------------------------------------------------------------------------
        !  Subroutine   :                    ppm_rc_gc
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
        !  Subroutine   :                    ppm_rc_gc
        !-------------------------------------------------------------------------
        !
        !  Purpose      : Gaussian curvature filter
        !
        !
        !  Input        : FieldIn, FieldOut, MeshIn
        !
        !  Input/output :
        !
        !  Output       : info  (I) return status. 0 on success.
        !
        !  Routines     :
        !
        !
        !  Remarks      : For this filter ghost size should be 2 layers, to
        !                 avoid any ghost update
        !
        !  References   :
        !
        !  Revisions    :
        !-------------------------------------------------------------------------
        SUBROUTINE DTYPE(ppm_rc_gc)(FieldIn,MeshIn,info)
          !-------------------------------------------------------------------------
          !  Modules
          !-------------------------------------------------------------------------
          IMPLICIT NONE

          !-------------------------------------------------------------------------
          !  Includes
          !-------------------------------------------------------------------------
          !-------------------------------------------------------------------------
          !  Arguments
          !-------------------------------------------------------------------------
          CLASS(ppm_t_field_),     POINTER       :: FieldIn
          !!! the source filed which is supposed to be

          CLASS(ppm_t_equi_mesh_), POINTER       :: MeshIn
          !!! source mesh on which fields are descritized

          INTEGER,                 INTENT(  OUT) :: info

          !-------------------------------------------------------------------------
          !  Local variables
          !-------------------------------------------------------------------------
          CLASS(ppm_t_subpatch_), POINTER :: sbpitr

#if   __DIME == __2D
          REAL(MK), CONTIGUOUS, DIMENSION(:,:),   POINTER :: DTYPE(wpin_r)
#elif __DIME == __3D
          REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER :: DTYPE(wpin_r)
#endif
          REAL(ppm_kind_double) :: t0

          INTEGER, DIMENSION(:), POINTER :: istart_
          INTEGER, DIMENSION(:), POINTER :: iend_
  !         INTEGER, DIMENSION(__DIME)          :: ghostsize__
          INTEGER, DIMENSION(__DIME)          :: istart
          INTEGER, DIMENSION(__DIME)          :: iend

          CHARACTER(LEN=ppm_char) :: caller = 'ppm_rc_gc'

          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          CALL check()

          sbpitr => MeshIn%subpatch%begin()
          DO WHILE (ASSOCIATED(sbpitr))
             istart_ => sbpitr%istart
             iend_   => sbpitr%iend

             NULLIFY(DTYPE(wpin_r))
             CALL sbpitr%get_field(FieldIn,DTYPE(wpin_r),info)
             or_fail("Failed to get field wpin_r data.")

             iend=MIN(iend_,MeshIn%Nm-1)-istart_+1
             iend(1:__DIME)=MERGE(iend(1:__DIME)+1,iend(1:__DIME),iend_(1:__DIME).NE.Ngrid(1:__DIME))

             istart=1+MOD(istart_,2)
             istart=MERGE(0,istart,istart_.NE.1.AND.istart.EQ.2)
             !start for the Black Triangle (BT)
             CALL Update(DTYPE(wpin_r),istart,iend)

             istart(1:2)=MERGE(3,2-MOD(istart_(1:2),2),istart_(1:2).EQ.1)
             istart(1:2)=MERGE(0,istart(1:2),istart(1:2).EQ.2)
             !start for the Black Circle (BC)
             CALL Update(DTYPE(wpin_r),istart,iend)

             istart(1)=1+MOD(istart_(1),2)
             istart(2)=MERGE(3,2-MOD(istart_(2),2),istart_(2).EQ.1)
#if   __DIME == __3D
             istart(3)=istart(3)+1
#endif
             istart=MERGE(0,istart,istart_.NE.1.AND.istart.EQ.2)
             !start for the White Triangle in 2D and Black Ellipse in 3D
             CALL Update(DTYPE(wpin_r),istart,iend)

             istart(1)=MERGE(3,2-MOD(istart_(1),2),istart_(1).EQ.1)
             istart(2)=1+MOD(istart_(2),2)
             istart=MERGE(0,istart,istart_.NE.1.AND.istart.EQ.2)
             !start for the White Circle in 2D and Black Rectangle in 3D
             CALL Update(DTYPE(wpin_r),istart,iend)

#if   __DIME == __3D
             istart(1)=1+MOD(istart_(1),2)
             istart(2)=MERGE(3,2-MOD(istart_(2),2),istart_(2).EQ.1)
             istart(3)=1+MOD(istart_(3),2)
             istart=MERGE(0,istart,istart_.NE.1.AND.istart.EQ.2)
             !start for the White Circle in 3D (WC)
             CALL Update(DTYPE(wpin_r),istart,iend)

             istart(1)=MERGE(3,2-MOD(istart_(1),2),istart_(1).EQ.1)
             istart(2)=1+MOD(istart_(2),2)
             istart=MERGE(0,istart,istart_.NE.1.AND.istart.EQ.2)
             !start for the White Triangle in 3D (WT)
             CALL Update(DTYPE(wpin_r),istart,iend)

             istart=1+MOD(istart_,2)
             istart(3)=istart(3)+1
             istart=MERGE(0,istart,istart_.NE.1.AND.istart.EQ.2)
             !start for the White Ellipse in 3D (WE)
             CALL Update(DTYPE(wpin_r),istart,iend)

             istart(1:2)=MERGE(3,2-MOD(istart_(1:2),2),istart_(1:2).EQ.1)
             istart(1:2)=MERGE(0,istart(1:2),istart(1:2).EQ.2)
             istart=MERGE(0,istart,istart_.NE.1.AND.istart.EQ.2)
             !start for the White Rectangle in 3D (WR)
             CALL Update(DTYPE(wpin_r),istart,iend)
#endif

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
        GOTO 8888
        8888 CONTINUE
        END SUBROUTINE check
#if   __DIME == __2D
        SUBROUTINE Update(wpin,istart__,iend__)

        IMPLICIT NONE

        REAL(MK), CONTIGUOUS, DIMENSION(:,:), POINTER       :: wpin

        INTEGER,              DIMENSION(:),   INTENT(IN   ) :: istart__
        INTEGER,              DIMENSION(:),   INTENT(IN   ) :: iend__

        REAL(MK), DIMENSION(8) :: d

        INTEGER :: i,j,n

        DO CONCURRENT (i=istart__(1):iend__(1):2,j=istart__(2):iend__(2):2)
           d(1)=(wpin(i-1,j  )+wpin(i+1,j  ))/two-wpin(i,j)
           d(2)=(wpin(i  ,j-1)+wpin(i  ,j+1))/two-wpin(i,j)
           d(3)=(wpin(i-1,j-1)+wpin(i+1,j+1))/two-wpin(i,j)
           d(4)=(wpin(i-1,j+1)+wpin(i+1,j-1))/two-wpin(i,j)

           d(5)=(wpin(i-1,j  )+wpin(i  ,j-1)+wpin(i-1,j-1))/three-wpin(i,j)
           d(6)=(wpin(i-1,j  )+wpin(i  ,j+1)+wpin(i-1,j+1))/three-wpin(i,j)
           d(7)=(wpin(i  ,j-1)+wpin(i+1,j  )+wpin(i+1,j-1))/three-wpin(i,j)
           d(8)=(wpin(i  ,j+1)+wpin(i+1,j  )+wpin(i+1,j+1))/three-wpin(i,j)

           n=MINLOC(ABS(d),DIM=1)

           wpin(i,j)=wpin(i,j)+d(n)
        ENDDO

        END SUBROUTINE Update
#elif __DIME == __3D
        SUBROUTINE Update(wpin,istart__,iend__)

        IMPLICIT NONE

        REAL(MK), CONTIGUOUS, DIMENSION(:,:,:), POINTER       :: wpin

        INTEGER,              DIMENSION(:),     INTENT(IN   ) :: istart__
        INTEGER,              DIMENSION(:),     INTENT(IN   ) :: iend__

        REAL(MK), DIMENSION(49) :: d

        INTEGER :: i,j,k,n

        DO CONCURRENT (i=istart__(1):iend__(1):2,j=istart__(2):iend__(2):2,k=istart__(3):iend__(3):2)
           d( 1)=(wpin(i  ,j  ,k-1)+wpin(i  ,j  ,k+1))/two-wpin(i,j,k)
           d( 2)=(wpin(i  ,j-1,k  )+wpin(i  ,j+1,k  ))/two-wpin(i,j,k)
           d( 3)=(wpin(i-1,j  ,k  )+wpin(i+1,j  ,k  ))/two-wpin(i,j,k)

           d( 4)=(wpin(i-1,j-1,k-1)+wpin(i+1,j+1,k+1))/two-wpin(i,j,k)
           d( 5)=(wpin(i  ,j-1,k-1)+wpin(i  ,j+1,k+1))/two-wpin(i,j,k)
           d( 6)=(wpin(i+1,j-1,k-1)+wpin(i-1,j+1,k+1))/two-wpin(i,j,k)
           d( 7)=(wpin(i-1,j  ,k-1)+wpin(i+1,j  ,k+1))/two-wpin(i,j,k)
           d( 8)=(wpin(i+1,j  ,k-1)+wpin(i-1,j  ,k+1))/two-wpin(i,j,k)
           d( 9)=(wpin(i-1,j+1,k-1)+wpin(i+1,j-1,k+1))/two-wpin(i,j,k)
           d(10)=(wpin(i  ,j+1,k-1)+wpin(i  ,j-1,k+1))/two-wpin(i,j,k)
           d(11)=(wpin(i+1,j+1,k-1)+wpin(i-1,j-1,k+1))/two-wpin(i,j,k)

           d(12)=(wpin(i-1,j-1,k  )+wpin(i+1,j+1,k  ))/two-wpin(i,j,k)
           d(13)=(wpin(i+1,j-1,k  )+wpin(i-1,j+1,k  ))/two-wpin(i,j,k)

           d(14)=(wpin(i  ,j  ,k+1)+wpin(i-1,j-1,k  )+wpin(i-1,j-1,k+1))/three-wpin(i,j,k)
           d(15)=(wpin(i  ,j  ,k+1)+wpin(i  ,j-1,k  )+wpin(i  ,j-1,k+1))/three-wpin(i,j,k)
           d(16)=(wpin(i  ,j  ,k+1)+wpin(i+1,j-1,k  )+wpin(i+1,j-1,k+1))/three-wpin(i,j,k)
           d(17)=(wpin(i  ,j  ,k+1)+wpin(i-1,j  ,k  )+wpin(i-1,j  ,k+1))/three-wpin(i,j,k)
           d(18)=(wpin(i  ,j  ,k+1)+wpin(i+1,j  ,k  )+wpin(i+1,j  ,k+1))/three-wpin(i,j,k)
           d(19)=(wpin(i  ,j  ,k+1)+wpin(i-1,j+1,k  )+wpin(i-1,j+1,k+1))/three-wpin(i,j,k)
           d(20)=(wpin(i  ,j  ,k+1)+wpin(i  ,j+1,k  )+wpin(i  ,j+1,k+1))/three-wpin(i,j,k)
           d(21)=(wpin(i  ,j  ,k+1)+wpin(i+1,j+1,k  )+wpin(i+1,j+1,k+1))/three-wpin(i,j,k)

           d(22)=(wpin(i  ,j  ,k-1)+wpin(i-1,j-1,k  )+wpin(i-1,j-1,k-1))/three-wpin(i,j,k)
           d(23)=(wpin(i  ,j  ,k-1)+wpin(i  ,j-1,k  )+wpin(i  ,j-1,k-1))/three-wpin(i,j,k)
           d(24)=(wpin(i  ,j  ,k-1)+wpin(i+1,j-1,k  )+wpin(i+1,j-1,k-1))/three-wpin(i,j,k)
           d(25)=(wpin(i  ,j  ,k-1)+wpin(i-1,j  ,k  )+wpin(i-1,j  ,k-1))/three-wpin(i,j,k)
           d(26)=(wpin(i  ,j  ,k-1)+wpin(i+1,j  ,k  )+wpin(i+1,j  ,k-1))/three-wpin(i,j,k)
           d(27)=(wpin(i  ,j  ,k-1)+wpin(i-1,j+1,k  )+wpin(i-1,j+1,k-1))/three-wpin(i,j,k)
           d(28)=(wpin(i  ,j  ,k-1)+wpin(i  ,j+1,k  )+wpin(i  ,j+1,k-1))/three-wpin(i,j,k)
           d(29)=(wpin(i  ,j  ,k-1)+wpin(i+1,j+1,k  )+wpin(i+1,j+1,k-1))/three-wpin(i,j,k)

           d(30)=(wpin(i-1,j  ,k+1)+wpin(i  ,j+1,k  )+wpin(i-1,j+1,k+1))/three-wpin(i,j,k)
           d(31)=(wpin(i-1,j  ,k+1)+wpin(i  ,j-1,k  )+wpin(i-1,j-1,k+1))/three-wpin(i,j,k)
           d(32)=(wpin(i+1,j  ,k+1)+wpin(i  ,j+1,k  )+wpin(i+1,j+1,k+1))/three-wpin(i,j,k)
           d(33)=(wpin(i+1,j  ,k+1)+wpin(i  ,j-1,k  )+wpin(i+1,j-1,k+1))/three-wpin(i,j,k)
           d(34)=(wpin(i-1,j  ,k-1)+wpin(i  ,j+1,k  )+wpin(i-1,j+1,k-1))/three-wpin(i,j,k)
           d(35)=(wpin(i-1,j  ,k-1)+wpin(i  ,j-1,k  )+wpin(i-1,j-1,k-1))/three-wpin(i,j,k)
           d(36)=(wpin(i+1,j  ,k-1)+wpin(i  ,j+1,k  )+wpin(i+1,j+1,k-1))/three-wpin(i,j,k)
           d(37)=(wpin(i+1,j  ,k-1)+wpin(i  ,j-1,k  )+wpin(i+1,j-1,k-1))/three-wpin(i,j,k)
           d(38)=(wpin(i-1,j  ,k  )+wpin(i  ,j-1,k  )+wpin(i-1,j-1,k  ))/three-wpin(i,j,k)
           d(39)=(wpin(i-1,j  ,k  )+wpin(i  ,j+1,k  )+wpin(i-1,j+1,k  ))/three-wpin(i,j,k)
           d(40)=(wpin(i+1,j  ,k  )+wpin(i  ,j-1,k  )+wpin(i+1,j-1,k  ))/three-wpin(i,j,k)
           d(41)=(wpin(i+1,j  ,k  )+wpin(i  ,j+1,k  )+wpin(i+1,j+1,k  ))/three-wpin(i,j,k)

           d(42)=(wpin(i  ,j-1,k+1)+wpin(i-1,j  ,k  )+wpin(i-1,j-1,k+1))/three-wpin(i,j,k)
           d(43)=(wpin(i  ,j-1,k+1)+wpin(i+1,j  ,k  )+wpin(i+1,j-1,k+1))/three-wpin(i,j,k)
           d(44)=(wpin(i  ,j+1,k+1)+wpin(i-1,j  ,k  )+wpin(i-1,j+1,k+1))/three-wpin(i,j,k)
           d(45)=(wpin(i  ,j+1,k+1)+wpin(i+1,j  ,k  )+wpin(i+1,j+1,k+1))/three-wpin(i,j,k)
           d(46)=(wpin(i  ,j-1,k-1)+wpin(i-1,j  ,k  )+wpin(i-1,j-1,k-1))/three-wpin(i,j,k)
           d(47)=(wpin(i  ,j-1,k-1)+wpin(i+1,j  ,k  )+wpin(i+1,j-1,k-1))/three-wpin(i,j,k)
           d(48)=(wpin(i  ,j+1,k-1)+wpin(i-1,j  ,k  )+wpin(i-1,j+1,k-1))/three-wpin(i,j,k)
           d(49)=(wpin(i  ,j+1,k-1)+wpin(i+1,j  ,k  )+wpin(i+1,j+1,k-1))/three-wpin(i,j,k)

           n=MINLOC(ABS(d),DIM=1)

           wpin(i,j,k)=wpin(i,j,k)+d(n)
        ENDDO
        END SUBROUTINE Update
#endif
        END SUBROUTINE DTYPE(ppm_rc_gc)

