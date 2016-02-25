      !-------------------------------------------------------------------------
      !  Module       :                    ppm_rc_module_bgl
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
      !  Module       :                    ppm_rc_module_bgl
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Creating graph and subgraphs using BOOST library
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      MODULE ppm_rc_module_bgl
        !----------------------------------------------------------------------
        !  Modules
        !----------------------------------------------------------------------
        USE ppm_module_data, ONLY : ppm_kind_double

        USE ISO_C_BINDING
        IMPLICIT NONE

        PRIVATE
        !----------------------------------------------------------------------
        !  Some work memory on the heap
        !----------------------------------------------------------------------
        !----------------------------------------------------------------------
        !  Define module interfaces
        !----------------------------------------------------------------------
        INTERFACE
          FUNCTION creategraph2d(npart_,vertices_neigh,nvertices_neigh,childsize) BIND(C,NAME='creategraph2d')
            IMPORT         :: C_INT,C_PTR
            INTEGER(C_INT) :: npart_
            INTEGER(C_INT) :: vertices_neigh(4,npart_)
            INTEGER(C_INT) :: nvertices_neigh(npart_)
            INTEGER(C_INT) :: childsize
            TYPE(C_PTR)    :: creategraph2d
          END FUNCTION
        END INTERFACE

        INTERFACE
          FUNCTION creategraph3d(npart_,vertices_neigh,nvertices_neigh,childsize) BIND(C,NAME='creategraph3d')
            IMPORT         :: C_INT,C_PTR
            INTEGER(C_INT) :: npart_
            INTEGER(C_INT) :: vertices_neigh(6,npart_)
            INTEGER(C_INT) :: nvertices_neigh(npart_)
            INTEGER(C_INT) :: childsize
            TYPE(C_PTR)    :: creategraph3d
          END FUNCTION
        END INTERFACE

        INTERFACE
          FUNCTION destroygraph() BIND(C,NAME='destroygraph')
            IMPORT         :: C_INT
            INTEGER(C_INT) :: destroygraph
          END FUNCTION
        END INTERFACE

        !----------------------------------------------------------------------
        ! Public
        !----------------------------------------------------------------------
        PUBLIC :: create_graph
        PUBLIC :: destroy_graph

      CONTAINS

        SUBROUTINE create_graph(Npart_,neigh,nneigh,subgraph,info)
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          USE ppm_rc_module_global
          IMPLICIT NONE

          INTEGER,                 INTENT(IN   ) :: Npart_
          INTEGER, DIMENSION(:,:), INTENT(IN   ) :: neigh
          INTEGER, DIMENSION(:),   INTENT(IN   ) :: nneigh
          INTEGER, DIMENSION(:),   POINTER       :: subgraph
          INTEGER,                 INTENT(  OUT) :: info

          REAL(ppm_kind_double) :: t0

          INTEGER(C_INT) :: gsize

          CHARACTER(LEN=ppm_char) :: caller='create_graph'

          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          IF (Npart_.LE.0) GOTO 9999

          SELECT CASE (SIZE(neigh,DIM=1))
          CASE (6)
             CALL C_F_POINTER(creategraph3d(Npart_,neigh,nneigh,gsize),subgraph,(/gsize/))

          CASE (4)
             CALL C_F_POINTER(creategraph2d(Npart_,neigh,nneigh,gsize),subgraph,(/gsize/))

          CASE DEFAULT
             fail('Space dimension must be 2 or 3!',ppm_err_wrong_dim,ppm_error=ppm_error_fatal)

          END SELECT

        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        END SUBROUTINE create_graph

        SUBROUTINE destroy_graph(subgraph,info)
          !----------------------------------------------------------------------
          !  Modules
          !----------------------------------------------------------------------
          USE ppm_rc_module_global, ONLY : MK,ppm_char,substart,substop
          IMPLICIT NONE

          INTEGER, DIMENSION(:), POINTER       :: subgraph
          INTEGER,               INTENT(  OUT) :: info

          REAL(ppm_kind_double) :: t0

          CHARACTER(LEN=ppm_char) :: caller='destroy_graph'

          !-------------------------------------------------------------------------
          !  Initialize
          !-------------------------------------------------------------------------
          CALL substart(caller,t0,info)

          IF (ASSOCIATED(subgraph)) THEN
             info=destroygraph()
          ENDIF
          NULLIFY(subgraph)

        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
        9999 CONTINUE
          CALL substop(caller,t0,info)
          RETURN
        END SUBROUTINE destroy_graph

      END MODULE ppm_rc_module_bgl
