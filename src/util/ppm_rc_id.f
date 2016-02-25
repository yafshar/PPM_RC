      !-------------------------------------------------------------------------
      !  Subroutine   :                    ppm_rc_id
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
      !  Subroutine   :                    ppm_rc_id
      !-------------------------------------------------------------------------
      !
      !  Purpose      : return the global id of each cell, giving cell index
      !               : giving cell number will return local id
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
      SUBROUTINE DTYPE(ppm_rc_id_LocalToGlobal_i)(istart,Nm,ld,ltg)
        !!! function for getting global id based on local mesh index

        IMPLICIT NONE
        !-------------------------------------------------------------------------
        !  Includes
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        INTEGER, DIMENSION(:), INTENT(IN   ) :: istart
        INTEGER, DIMENSION(:), INTENT(IN   ) :: Nm
        INTEGER, DIMENSION(:), INTENT(IN   ) :: ld
        INTEGER,               INTENT(  OUT) :: ltg
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
        INTEGER :: i,j
#if   __DIME == __3D
        INTEGER :: k
#endif
        i=ld(1)+istart(1)-1
        j=ld(2)+istart(2)-2
#if   __DIME == __2D
        ltg=i+Nm(1)*j
#elif __DIME == __3D
        k=ld(3)+istart(3)-2
        ltg=i+Nm(1)*(j+Nm(2)*k)
#endif
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        RETURN
      END SUBROUTINE DTYPE(ppm_rc_id_LocalToGlobal_i)

      SUBROUTINE DTYPE(ppm_rc_id_LocalToGlobal_li)(istart,Nm,ld,ltg)
        !!! function for getting global id based on local mesh index

        IMPLICIT NONE
        !-------------------------------------------------------------------------
        !  Includes
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        INTEGER, DIMENSION(:),   INTENT(IN   ) :: istart
        INTEGER, DIMENSION(:),   INTENT(IN   ) :: Nm
        INTEGER, DIMENSION(:),   INTENT(IN   ) :: ld
        INTEGER(ppm_kind_int64), INTENT(  OUT) :: ltg
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
        INTEGER(ppm_kind_int64) :: i,j
#if   __DIME == __3D
        INTEGER(ppm_kind_int64) :: k
#endif
        i=INT(ld(1)+istart(1)-1,ppm_kind_int64)
        j=INT(ld(2)+istart(2)-2,ppm_kind_int64)
#if   __DIME == __2D
        ltg=i+Nm(1)*j
#elif __DIME == __3D
        k=INT(ld(3)+istart(3)-2,ppm_kind_int64)
        ltg=i+Nm(1)*(j+Nm(2)*k)
#endif
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        RETURN
      END SUBROUTINE DTYPE(ppm_rc_id_LocalToGlobal_li)

      SUBROUTINE DTYPE(ppm_rc_id_GlobalToLocal_i)(istart,Nm,ld,gtl)
        !!! function for getting local mesh index based on global id

        IMPLICIT NONE
        !-------------------------------------------------------------------------
        !  Includes
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        INTEGER,   DIMENSION(:), INTENT(IN   ) :: istart
        INTEGER,   DIMENSION(:), INTENT(IN   ) :: Nm
        INTEGER,                 INTENT(IN   ) :: ld
        INTEGER,   DIMENSION(:), INTENT(  OUT) :: gtl
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
        INTEGER :: i,j

#if   __DIME == __3D
        INTEGER :: k
#endif
        ! because the number of cell has been increased by one
#if   __DIME == __2D
        j=(ld-1)/Nm(1)
        i=ld-Nm(1)*j

        gtl(1)=i+1-istart(1)
        gtl(2)=j+2-istart(2)
#elif __DIME == __3D
        k=(ld-1)/(Nm(1)*Nm(2))
        j=(ld-1-k*Nm(1)*Nm(2))/Nm(1)
        i=ld-Nm(1)*(j+Nm(2)*k)

        gtl(1)=i+1-istart(1)
        gtl(2)=j+2-istart(2)
        gtl(3)=k+2-istart(3)
#endif
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        RETURN
      END SUBROUTINE DTYPE(ppm_rc_id_GlobalToLocal_i)

      SUBROUTINE DTYPE(ppm_rc_id_GlobalToLocal_li)(istart,Nm,ld,gtl)
        !!! function for getting local mesh index based on global id

        IMPLICIT NONE
        !-------------------------------------------------------------------------
        !  Includes
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        INTEGER,   DIMENSION(:), INTENT(IN   ) :: istart
        INTEGER,   DIMENSION(:), INTENT(IN   ) :: Nm
        INTEGER(ppm_kind_int64), INTENT(IN   ) :: ld
        INTEGER,   DIMENSION(:), INTENT(  OUT) :: gtl
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
        INTEGER(ppm_kind_int64), DIMENSION(__DIME) :: Nm_
        INTEGER(ppm_kind_int64) :: i,j

#if   __DIME == __3D
        INTEGER(ppm_kind_int64) :: k
#endif
        Nm_=INT(Nm,ppm_kind_int64)
        ! because the number of cell has been increased by one
#if   __DIME == __2D
        j=(ld-1_ppm_kind_int64)/Nm_(1)
        i=ld-Nm_(1)*j

        gtl(1)=INT(i)+1-istart(1)
        gtl(2)=INT(j)+2-istart(2)
#elif __DIME == __3D
        k=(ld-1_ppm_kind_int64)/(Nm_(1)*Nm_(2))
        j=(ld-1_ppm_kind_int64-k*Nm_(1)*Nm_(2))/Nm_(1)
        i=ld-Nm_(1)*(j+Nm_(2)*k)

        gtl(1)=INT(i)+1-istart(1)
        gtl(2)=INT(j)+2-istart(2)
        gtl(3)=INT(k)+2-istart(3)
#endif
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        RETURN
      END SUBROUTINE DTYPE(ppm_rc_id_GlobalToLocal_li)

      SUBROUTINE DTYPE(ppm_rc_id_GlobalToLocal__i)(Nm,ld,gtl)
        !!! function for getting local mesh index based on global id

        IMPLICIT NONE
        !-------------------------------------------------------------------------
        !  Includes
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        INTEGER,   DIMENSION(:), INTENT(IN   ) :: Nm
        INTEGER,                 INTENT(IN   ) :: ld
        INTEGER,   DIMENSION(:), INTENT(  OUT) :: gtl
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
        INTEGER :: j
#if   __DIME == __3D
        INTEGER :: k
#endif

        ! because the number of cell has been increased by one
#if   __DIME == __2D
        j=(ld-1)/Nm(1)
        gtl(1)=ld-Nm(1)*j
        gtl(2)=j+1
#elif __DIME == __3D
        k=(ld-1)/(Nm(1)*Nm(2))
        j=(ld-1-k*Nm(1)*Nm(2))/Nm(1)
        gtl(1)=ld-Nm(1)*(j+Nm(2)*k)
        gtl(2)=j+1
        gtl(3)=k+1
#endif
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        RETURN
      END SUBROUTINE DTYPE(ppm_rc_id_GlobalToLocal__i)

      SUBROUTINE DTYPE(ppm_rc_id_GlobalToLocal__li)(Nm,ld,gtl)
        !!! function for getting local mesh index based on global id

        IMPLICIT NONE
        !-------------------------------------------------------------------------
        !  Includes
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        INTEGER,   DIMENSION(:), INTENT(IN   ) :: Nm
        INTEGER(ppm_kind_int64), INTENT(IN   ) :: ld
        INTEGER,   DIMENSION(:), INTENT(  OUT) :: gtl
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
        INTEGER(ppm_kind_int64), DIMENSION(__DIME) :: Nm_
        INTEGER(ppm_kind_int64) :: j

#if   __DIME == __3D
        INTEGER(ppm_kind_int64) :: i,k
#endif
        Nm_=INT(Nm,ppm_kind_int64)
        ! because the number of cell has been increased by one
#if   __DIME == __2D
        j=(ld-1_ppm_kind_int64)/Nm_(1)
        gtl(1)=INT(ld-Nm_(1)*j)
        gtl(2)=INT(j)+1
#elif __DIME == __3D
        k=(ld-1_ppm_kind_int64)/(Nm_(1)*Nm_(2))
        j=(ld-1_ppm_kind_int64-k*Nm_(1)*Nm_(2))/Nm_(1)
        i=ld-Nm_(1)*(j+Nm_(2)*k)
        gtl(1)=INT(i)
        gtl(2)=INT(j)+1
        gtl(3)=INT(k)+1
#endif
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        RETURN
      END SUBROUTINE DTYPE(ppm_rc_id_GlobalToLocal__li)

      INTEGER FUNCTION DTYPE(ppm_rc_ghost_id)(sbpitr,ld)
        !!! function for getting local ghost mesh index based on ghost id

        IMPLICIT NONE
        !-------------------------------------------------------------------------
        !  Includes
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        CLASS(ppm_t_subpatch_), POINTER :: sbpitr

        INTEGER, DIMENSION(:),  POINTER :: ld
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
        INTEGER, DIMENSION(:), POINTER :: Nm
        INTEGER                        :: nsize
        !-------------------------------------------------------------------------
        !  Initialize
        !-------------------------------------------------------------------------
        Nm => sbpitr%nnodes

#if   __DIME == __2D
        !-y
        IF (ld(2).EQ.0)       THEN
           DTYPE(ppm_rc_ghost_id)=ld(1)
        ELSE IF (ld(2).EQ.Nm(2)+1) THEN
           nsize=Nm(1)
           DTYPE(ppm_rc_ghost_id)=nsize+ld(1)
        ELSE IF (ld(1).EQ.0)       THEN
           nsize=2*Nm(1)
           DTYPE(ppm_rc_ghost_id)=nsize+ld(2)+1
        ELSE IF (ld(1).EQ.Nm(1)+1) THEN
           nsize=2*Nm(1)+(Nm(2)+2)
           DTYPE(ppm_rc_ghost_id)=nsize+ld(2)+1
        ENDIF
#elif __DIME == __3D
        !-z
        IF      (ld(3).EQ.0)       THEN
           DTYPE(ppm_rc_ghost_id)=ld(1)+(ld(2)-1)*Nm(1)
        ELSE IF (ld(3).EQ.Nm(3)+1) THEN
           nsize=Nm(1)*Nm(2)
           DTYPE(ppm_rc_ghost_id)=nsize+ld(1)+(ld(2)-1)*Nm(1)
        ELSE IF (ld(2).EQ.0)       THEN
           nsize=2*Nm(1)*Nm(2)
           DTYPE(ppm_rc_ghost_id)=nsize+ld(1)+ld(3)*Nm(1)
        ELSE IF (ld(2).EQ.Nm(2)+1) THEN
           nsize=2*Nm(1)*Nm(2)+Nm(1)*(Nm(3)+2)
           DTYPE(ppm_rc_ghost_id)=nsize+ld(1)+ld(3)*Nm(1)
        ELSE IF (ld(1).EQ.0)       THEN
           nsize=2*Nm(1)*(Nm(2)+Nm(3)+2)
           DTYPE(ppm_rc_ghost_id)=nsize+ld(2)+ld(3)*(Nm(2)+2)+1
        ELSE IF (ld(1).EQ.Nm(1)+1) THEN
           nsize=2*Nm(1)*Nm(2)+(2*Nm(1)+(Nm(2)+2))*(Nm(3)+2)
           DTYPE(ppm_rc_ghost_id)=nsize+ld(2)+ld(3)*(Nm(2)+2)+1
        ENDIF
#endif
      END FUNCTION DTYPE(ppm_rc_ghost_id)

      FUNCTION DTYPE(ppm_rc_ghost_lidx)(sbpitr,ld)
        !!! function for getting local ghost mesh index based on ghost id

        IMPLICIT NONE
        !-------------------------------------------------------------------------
        !  Includes
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        CLASS(ppm_t_subpatch_), POINTER :: sbpitr

        INTEGER,                POINTER :: ld
        INTEGER, DIMENSION(__DIME)           :: DTYPE(ppm_rc_ghost_lidx)
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
        INTEGER, DIMENSION(:), POINTER :: Nm
#if __DIME == __3D
        INTEGER                        :: z1,z2
#endif
        INTEGER                        :: y1,y2
        INTEGER                        :: x1,x2
        !-------------------------------------------------------------------------
        !  Initialize
        !-------------------------------------------------------------------------
        Nm => sbpitr%nnodes

#if __DIME == __2D
        !-y
        y1=Nm(1)
        !+y
        y2=2*y1
        IF      (ld.LE.y1) THEN
           DTYPE(ppm_rc_ghost_lidx)(2)=0
           DTYPE(ppm_rc_ghost_lidx)(1)=ld
        ELSE IF (ld.LE.y2) THEN
           y1=ld-y1
           DTYPE(ppm_rc_ghost_lidx)(2)=Nm(2)+1
           DTYPE(ppm_rc_ghost_lidx)(1)=y1
        ELSE
           y1=ld-y2
           !-x
           x1=Nm(2)+2
           !+x
           x2=2*x1
           IF      (y1.LE.x1) THEN
              DTYPE(ppm_rc_ghost_lidx)(2)=y1-1
              DTYPE(ppm_rc_ghost_lidx)(1)=0
           ELSE IF (y1.LE.x2) THEN
              x1=y1-x1
              DTYPE(ppm_rc_ghost_lidx)(2)=x1-1
              DTYPE(ppm_rc_ghost_lidx)(1)=Nm(1)+1
           ENDIF
        ENDIF
#elif __DIME == __3D
        !-z
        z1=Nm(1)*Nm(2)
        !+z
        z2=2*z1

        IF      (ld.LE.z1) THEN
           DTYPE(ppm_rc_ghost_lidx)(3)=0
           DTYPE(ppm_rc_ghost_lidx)(2)=(ld-1)/Nm(1)
           DTYPE(ppm_rc_ghost_lidx)(1)=ld-DTYPE(ppm_rc_ghost_lidx)(2)*Nm(1)
           DTYPE(ppm_rc_ghost_lidx)(2)=DTYPE(ppm_rc_ghost_lidx)(2)+1
        ELSE IF (ld.LE.z2) THEN
           z1=ld-z1
           DTYPE(ppm_rc_ghost_lidx)(3)=Nm(3)+1
           DTYPE(ppm_rc_ghost_lidx)(2)=(z1-1)/Nm(1)
           DTYPE(ppm_rc_ghost_lidx)(1)=z1-DTYPE(ppm_rc_ghost_lidx)(2)*Nm(1)
           DTYPE(ppm_rc_ghost_lidx)(2)=DTYPE(ppm_rc_ghost_lidx)(2)+1
        ELSE
           z1=ld-z2
           !-y
           y1=Nm(1)*(Nm(3)+2)
           !+y
           y2=2*y1
           IF      (z1.LE.y1) THEN
              DTYPE(ppm_rc_ghost_lidx)(3)=(z1-1)/Nm(1)
              DTYPE(ppm_rc_ghost_lidx)(2)=0
              DTYPE(ppm_rc_ghost_lidx)(1)=z1-DTYPE(ppm_rc_ghost_lidx)(3)*Nm(1)
           ELSE IF (z1.LE.y2) THEN
              y1=z1-y1
              DTYPE(ppm_rc_ghost_lidx)(3)=(y1-1)/Nm(1)
              DTYPE(ppm_rc_ghost_lidx)(2)=Nm(2)+1
              DTYPE(ppm_rc_ghost_lidx)(1)=y1-DTYPE(ppm_rc_ghost_lidx)(3)*Nm(1)
           ELSE
              y1=z1-y2
              !-x
              x1=(Nm(2)+2)*(Nm(3)+2)
              !+x
              x2=2*x1
              IF      (y1.LE.x1) THEN
                 DTYPE(ppm_rc_ghost_lidx)(3)=(y1-1)/(Nm(1)+2)
                 DTYPE(ppm_rc_ghost_lidx)(2)=(y1-1)-DTYPE(ppm_rc_ghost_lidx)(3)*(Nm(2)+2)
                 DTYPE(ppm_rc_ghost_lidx)(1)=0
              ELSE IF (y1.LE.x2) THEN
                 x1=y1-x1
                 DTYPE(ppm_rc_ghost_lidx)(3)=(x1-1)/(Nm(1)+2)
                 DTYPE(ppm_rc_ghost_lidx)(2)=(x1-1)-DTYPE(ppm_rc_ghost_lidx)(3)*(Nm(2)+2)
                 DTYPE(ppm_rc_ghost_lidx)(1)=Nm(1)+1
              ENDIF
           ENDIF
        ENDIF
#endif
      END FUNCTION DTYPE(ppm_rc_ghost_lidx)

      PURE FUNCTION DTYPE(IndexHashFunctor32)(index_) RESULT(key)
        !!!
        IMPLICIT NONE
        !-------------------------------------------------------------------------
        !  Includes
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        INTEGER, DIMENSION(:), INTENT(IN   ) :: index_
        INTEGER                              :: key
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
#if   __DIME == __2D
        key=IOR(SHIFTL(index_(2),16),index_(1))
#elif __DIME == __3D
        key=IOR(SHIFTL(IAND(index_(3),1023),22),IOR(SHIFTL(IAND(index_(2),2047),11),IAND(index_(1),2047)))
#endif
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        RETURN
      END FUNCTION DTYPE(IndexHashFunctor32)

      PURE FUNCTION DTYPE(IndexHashFunctor64)(index_) RESULT(key)
        !!!
        IMPLICIT NONE
        !-------------------------------------------------------------------------
        !  Includes
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        INTEGER, DIMENSION(:), INTENT(IN   ) :: index_
        INTEGER(ppm_kind_int64)              :: key
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
        INTEGER(ppm_kind_int64), DIMENSION(__DIME) :: index_64
        index_64=INT(index_,ppm_kind_int64)
#if   __DIME == __2D
        key=IOR(SHIFTL(index_64(2),16),index_64(1))
#elif __DIME == __3D
        key=IOR(SHIFTL(index_64(3),32),IOR(SHIFTL(index_64(2),16),index_64(1)))
#endif
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        RETURN
      END FUNCTION DTYPE(IndexHashFunctor64)

#if   __DIME == __2D
      ELEMENTAL FUNCTION DTYPE(IndexHashFunctor32_)(index_1,index_2) RESULT(key)
#elif __DIME == __3D
      ELEMENTAL FUNCTION DTYPE(IndexHashFunctor32_)(index_1,index_2,index_3) RESULT(key)
#endif
        !!!
        IMPLICIT NONE
        !-------------------------------------------------------------------------
        !  Includes
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        INTEGER, INTENT(IN   ) :: index_1
        INTEGER, INTENT(IN   ) :: index_2
#if   __DIME == __3D
        INTEGER, INTENT(IN   ) :: index_3
#endif
        INTEGER                :: key
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
#if   __DIME == __2D
        key=IOR(SHIFTL(index_2,16),IAND(index_1,65535))
#elif __DIME == __3D
        key=IOR(SHIFTL(IAND(index_3,1023),22),IOR(SHIFTL(IAND(index_2,2047),11),IAND(index_1,2047)))
#endif
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        RETURN
      END FUNCTION DTYPE(IndexHashFunctor32_)


#if   __DIME == __2D
      ELEMENTAL FUNCTION DTYPE(IndexHashFunctor64_)(index_1,index_2) RESULT(key)
#elif __DIME == __3D
      ELEMENTAL FUNCTION DTYPE(IndexHashFunctor64_)(index_1,index_2,index_3) RESULT(key)
#endif
        !!!
        IMPLICIT NONE
        !-------------------------------------------------------------------------
        !  Includes
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        INTEGER,  INTENT(IN   ) :: index_1
        INTEGER,  INTENT(IN   ) :: index_2
#if   __DIME == __3D
        INTEGER,  INTENT(IN   ) :: index_3
#endif
        INTEGER(ppm_kind_int64) :: key
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
        INTEGER(ppm_kind_int64) :: index1,index2
#if   __DIME == __2D
        index1=INT(index_1,ppm_kind_int64)
        index2=INT(index_2,ppm_kind_int64)
        key=IOR(SHIFTL(index2,16),index1)
#elif __DIME == __3D
        INTEGER(ppm_kind_int64) :: index3
        index1=INT(index_1,ppm_kind_int64)
        index2=INT(index_2,ppm_kind_int64)
        index3=INT(index_3,ppm_kind_int64)
        key=IOR(SHIFTL(index3,32),IOR(SHIFTL(index2,16),index1))
#endif
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        RETURN
      END FUNCTION DTYPE(IndexHashFunctor64_)

      SUBROUTINE DTYPE(HashIndexFunctor32)(key,index_)
        !!!
        IMPLICIT NONE
        !-------------------------------------------------------------------------
        !  Includes
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        INTEGER,               INTENT(IN   ) :: key
        INTEGER, DIMENSION(:), INTENT(  OUT) :: index_
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
#if   __DIME == __2D
        index_(1)=IBITS(key,0,16)
        index_(2)=IBITS(key,16,16)
#elif __DIME == __3D
        index_(1)=IBITS(key,0,11)
        index_(2)=IBITS(key,11,11)
        index_(3)=IBITS(key,22,10)
#endif
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        RETURN
      END SUBROUTINE DTYPE(HashIndexFunctor32)

#if   __DIME == __2D
      SUBROUTINE DTYPE(HashIndexFunctor32_)(key,index_1,index_2)
#elif __DIME == __3D
      SUBROUTINE DTYPE(HashIndexFunctor32_)(key,index_1,index_2,index_3)
#endif
        !!!
        IMPLICIT NONE
        !-------------------------------------------------------------------------
        !  Includes
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        INTEGER, INTENT(IN   ) :: key
        INTEGER, INTENT(  OUT) :: index_1
        INTEGER, INTENT(  OUT) :: index_2
#if   __DIME == __3D
        INTEGER, INTENT(  OUT) :: index_3
#endif
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
#if   __DIME == __2D
        index_1=IBITS(key,0,16)
        index_2=IBITS(key,16,16)
#elif __DIME == __3D
        index_1=IBITS(key,0,11)
        index_2=IBITS(key,11,11)
        index_3=IBITS(key,22,10)
#endif
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        RETURN
      END SUBROUTINE DTYPE(HashIndexFunctor32_)

      SUBROUTINE DTYPE(HashIndexFunctor64)(key,index_)
        !!!
        IMPLICIT NONE
        !-------------------------------------------------------------------------
        !  Includes
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        INTEGER(ppm_kind_int64), INTENT(IN   ) :: key
        INTEGER, DIMENSION(:),   INTENT(  OUT) :: index_
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
        index_(1)=IBITS(key,0,16)
        index_(2)=IBITS(key,16,16)
#if   __DIME == __3D
        index_(3)=IBITS(key,32,16)
#endif
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        RETURN
      END SUBROUTINE DTYPE(HashIndexFunctor64)

#if   __DIME == __2D
      SUBROUTINE DTYPE(HashIndexFunctor64_)(key,index_1,index_2)
#elif __DIME == __3D
      SUBROUTINE DTYPE(HashIndexFunctor64_)(key,index_1,index_2,index_3)
#endif
        !!!
        IMPLICIT NONE
        !-------------------------------------------------------------------------
        !  Includes
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        INTEGER(ppm_kind_int64), INTENT(IN   ) :: key
        INTEGER,                 INTENT(  OUT) :: index_1
        INTEGER,                 INTENT(  OUT) :: index_2
#if   __DIME == __3D
        INTEGER,                 INTENT(  OUT) :: index_3
#endif
        !-------------------------------------------------------------------------
        !  Local variables
        !-------------------------------------------------------------------------
        index_1=IBITS(key,0,16)
        index_2=IBITS(key,16,16)
#if   __DIME == __3D
        index_3=IBITS(key,32,16)
#endif
        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        RETURN
      END SUBROUTINE DTYPE(HashIndexFunctor64_)
