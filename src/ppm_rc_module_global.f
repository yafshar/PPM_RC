      !-------------------------------------------------------------------------
      !  Module       :                    ppm_rc_module_global
      !-------------------------------------------------------------------------
      !
      !  Purpose      :  Declare global types and variables for ppm_rc.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  MOSAIC Group
      !  Max Planck Institute of Molecular Cell Biology and Genetics
      !  Pfotenhauerstr. 108, 01307 Dresden, Germany
      !
      !  Adapted           - y.afshar           June   2014
      !-------------------------------------------------------------------------

      MODULE ppm_rc_module_global
        !----------------------------------------------------------------------
        !  Modules
        !----------------------------------------------------------------------
        USE ppm_module_data, ONLY : ppm_kind_double,ppm_char,&
        & ppm_kind_single,ppm_integer,ppm_myepss,            &
        & ppm_kind_int32,ppm_kind_int64,ppm_myepsd,          &
        & ppm_type_int,ppm_type_real,ppm_rank,ppm_nproc,     &
        & ppm_mpi_kind,ppm_type_longint,                     &
        & ppm_param_io_ascii,ppm_param_io_binary,            &
        & ppm_param_bcdef_freespace,ppm_error_error,         &
        & ppm_param_alloc_fit,ppm_param_alloc_grow,          &
        & ppm_param_alloc_grow_preserve,                     &
        & ppm_param_dealloc,ppm_param_alloc_fit_preserve,    &
        & ppm_param_pop_replace,ppm_error_fatal,             &
        & ppm_param_bcdef_periodic,ppm_type_real_single,     &
        & ppm_kind
        USE ppm_module_alloc, ONLY : ppm_alloc
        USE ppm_module_substart, ONLY : substart
        USE ppm_module_substop, ONLY : substop
        USE ppm_module_error, ONLY : ppm_error,ppm_err_sub_failed, &
        & ppm_err_alloc, ppm_err_mpi_fail, ppm_err_dealloc, &
        & ppm_err_argument,ppm_err_wrong_dim
        USE ppm_module_util_time, ONLY : ppm_util_time
        USE ppm_module_interfaces, ONLY : ppm_t_discr_kind, &
        & ppm_part_areinside,&
        & ppm_t_discr_info_,ppm_part_partial,ppm_part_ghosts, &
        & ppm_ppt_map_ghosts,ppm_ppt_map_parts,ppm_ppt_partial, &
        & ppm_t_particles_d_,ppm_t_part_prop_d_,                &
        & ppm_t_particles_s_,ppm_t_part_prop_s_
        USE ppm_module_field_typedef, ONLY : ppm_t_field_,ppm_t_field
        USE ppm_module_mesh_typedef, ONLY : ppm_t_equi_mesh_,ppm_t_equi_mesh,&
        & ppm_mesh,ppm_t_subpatch_,ppm_t_subpatch
        USE ppm_module_write, ONLY : ppm_write
        USE ppm_module_log, ONLY : ppm_log
        USE ppm_module_inl_hash, ONLY : ppm_htable,htable_null
        USE ppm_module_particles_typedef, ONLY : ppm_t_particles_d,ppm_t_particles_s
        USE ppm_module_ctrl
        USE ppm_module_io_vtk
        IMPLICIT NONE

        !----------------------------------------------------------------------
        !  Header file
        !----------------------------------------------------------------------

        !-------------------------------------------------------------------------
        !  Precision
        ! Note, To do single Precision segmentation, you need to change the
        ! Part datatype in all the routines, currently it is the ppm_t_particles_d_
        !-------------------------------------------------------------------------
        INTEGER,                  PARAMETER :: MK        = ppm_kind_single
        INTEGER,                  PARAMETER :: FORBIDDEN = HUGE(0)-1
        PUBLIC :: MK,FORBIDDEN

        !----------------------------------------------------------------------
        !  Define constant parameters
        !----------------------------------------------------------------------
        REAL(MK),                 PARAMETER :: pi = ACOS(-1.0_MK)
        REAL(ppm_kind_double),    PARAMETER :: pid= REAL(pi,ppm_kind_double)

        REAL(MK),                 PARAMETER :: small=EPSILON(1.0_MK)
        REAL(ppm_kind_double),    PARAMETER :: smalld=REAL(small,ppm_kind_double)
        REAL(MK),                 PARAMETER :: smallest=TINY(1.0_MK)
        REAL(ppm_kind_double),    PARAMETER :: smallestd=REAL(smallest,ppm_kind_double)
        ! This factor has been chosen based on tha maximum number of
        ! Particles on each processor which could be 2^32 times by 2^-32
        ! will create one move anything less than this is obsured

        REAL(MK),                 PARAMETER :: big=HUGE(1.0_MK)/REAL(HUGE(1),MK)
        REAL(MK),                 PARAMETER :: bigs=SQRT(REAL(HUGE(1.0)/REAL(HUGE(1)),MK))
        REAL(ppm_kind_double),    PARAMETER :: bigsd=REAL(bigs,ppm_kind_double)

        REAL(MK),                 PARAMETER :: zero=0.0_MK
        REAL(ppm_kind_double),    PARAMETER :: zerod=0.0_ppm_kind_double
        REAL(MK),                 PARAMETER :: one=1.0_MK
        REAL(ppm_kind_double),    PARAMETER :: oned=REAL(one,ppm_kind_double)
        REAL(MK),                 PARAMETER :: oneplus=NEAREST(1.0_MK,1.0_MK)
        REAL(ppm_kind_double),    PARAMETER :: oneplusd=REAL(oneplus,ppm_kind_double)
        REAL(MK),                 PARAMETER :: oneminus=NEAREST(1.0_MK,-1.0_MK)
        REAL(ppm_kind_double),    PARAMETER :: oneminusd=REAL(oneminus,ppm_kind_double)
        REAL(MK),                 PARAMETER :: two=2.0_MK
        REAL(ppm_kind_double),    PARAMETER :: twod=2.0_ppm_kind_double
        REAL(MK),                 PARAMETER :: three=3.0_MK
        REAL(ppm_kind_double),    PARAMETER :: threed=3.0_ppm_kind_double
        REAL(MK),                 PARAMETER :: half=0.5_MK
        REAL(ppm_kind_double),    PARAMETER :: halfd=0.5_ppm_kind_double
        REAL(MK),                 PARAMETER :: ten=10.0_MK
        REAL(ppm_kind_double),    PARAMETER :: tend=10.0_ppm_kind_double

        INTEGER,                  PARAMETER :: bigi=HUGE(0)-2

        PUBLIC :: pi,pid
        PUBLIC :: small,smalld,smallest,smallestd
        PUBLIC :: big,bigs,bigsd,bigi
        PUBLIC :: zero,zerod
        PUBLIC :: one,oned,oneplus,oneplusd,oneminus,oneminusd
        PUBLIC :: two,twod,three,threed
        PUBLIC :: half,halfd
        PUBLIC :: ten,tend

        !----------------------------------------------------------------------
        !  InitializationKindType
        !----------------------------------------------------------------------
        INTEGER,                  PARAMETER :: e_fromFile=1
        INTEGER,                  PARAMETER :: e_rect    =2
        INTEGER,                  PARAMETER :: e_sphere  =3
        INTEGER,                  PARAMETER :: e_otsu    =4
        INTEGER,                  PARAMETER :: e_blob_det=5
        INTEGER,                  PARAMETER :: e_localmax=6

        CHARACTER(LEN=ppm_char)             :: init_mode

        PUBLIC :: e_fromFile,e_rect,e_sphere,e_otsu,e_blob_det,e_localmax,init_mode

        !----------------------------------------------------------------------
        !  Global TYPEs
        !----------------------------------------------------------------------
        TYPE(ppm_htable),           POINTER :: htable => NULL()
        !!! hash table for the region labels

        PUBLIC :: htable
        PUBLIC :: htable_null
        !----------------------------------------------------------------------
        !  Topology
        !----------------------------------------------------------------------
        !Digital topology (control)
        !TODO
        LOGICAL                             :: AllowFusion
        !!! Allow fusion of regions?
        !!! If true, regions will compete and fuse.
        !!! If allow_fusion is false, fronts stop propagation if they meet.
        LOGICAL                             :: AllowFusionZ

        LOGICAL                             :: AllowFission
        !!! Allow fusion of regions?
        LOGICAL                             :: AllowHandles
        !!! Allow handels at regions

        PUBLIC :: AllowFusion,AllowFusionZ,AllowFission,AllowHandles

        REAL(MK)                            :: AcceptedPointsFactor
        !!! stores acceptance reduction factor
        REAL(MK),                 PARAMETER :: AcceptedPointsReductionFactor = half
        !!! how much to reduce the number of accepted points

        INTEGER,                  PARAMETER :: OscillationHistoryLength = 10

        REAL(ppm_kind_double)               :: OscillationThreshold

        INTEGER                             :: initopoid
        !!! topology id for initialization
        INTEGER                             :: initmeshid
        !!! mesh id for initialization

        PUBLIC :: AcceptedPointsFactor,AcceptedPointsReductionFactor
        PUBLIC :: OscillationHistoryLength,OscillationThreshold
        PUBLIC :: initopoid,initmeshid

        INTEGER, DIMENSION(:),   ALLOCATABLE :: inighostsize
        !!! size of the ghostlayers in pixels used during init process
        !!! This is used in create_particles and depends on the convolution
        !!! kernel
        CLASS(ppm_t_equi_mesh_), POINTER    :: imesh => NULL()
        !!! mesh structure that is used during init process

        INTEGER                             :: topoid
        !!! main topology id for segmentation
        CLASS(ppm_t_equi_mesh_), POINTER    :: mesh => NULL()
        !!! Main mesh structure during the segmentation

        INTEGER                             :: meshid
        !!! main Mesh id for segmentation

        PUBLIC :: inighostsize,imesh,topoid,mesh,meshid
        !----------------------------------------------------------------------
        !  MCMC
        !----------------------------------------------------------------------
        LOGICAL                             :: UseMCMC
        !!! Algorithm is intended to use in MCMC mode for segmentation?
        LOGICAL                             :: MCMCcontinue
        !!! Continue MCMC run (no relabeling)?
        LOGICAL                             :: MCMCuseBiasedProposal
        !!! Use biased proposal for MCMC sampling?
        LOGICAL                             :: MCMCusePairProposal
        !!! One move corresponds to a neighboring pair of particles
        INTEGER                             :: MCMCstepsize
        !!! The step size for Markov-chain Monte-Carlo method.
        !!! This determines the posterior pdf of a Gibbs-Boltzmann distribution.
        REAL(MK)                            :: MCMCtemperature
        REAL(MK)                            :: MCMCburnInFactor

        PUBLIC :: UseMCMC,MCMCcontinue,MCMCuseBiasedProposal,MCMCusePairProposal
        PUBLIC :: MCMCstepsize,MCMCtemperature,MCMCburnInFactor
        !----------------------------------------------------------------------
        !  Domain decomposition
        !----------------------------------------------------------------------
        REAL(MK), DIMENSION(:),  POINTER    :: proc_speed => NULL()
        REAL(MK), DIMENSION(:),  POINTER    :: cost => NULL()

        REAL(MK), DIMENSION(:), POINTER     :: min_phys
        REAL(MK), DIMENSION(:), POINTER     :: max_phys

        ! size of the ghost layers in pixels
        INTEGER,  DIMENSION(:), ALLOCATABLE :: ghostsize
        INTEGER,  DIMENSION(:), ALLOCATABLE :: bcdef

        INTEGER                             :: decomp
        INTEGER                             :: assig

        INTEGER, DIMENSION(:),  ALLOCATABLE :: ineighproc
        INTEGER                             :: nneighproc

        INTEGER, DIMENSION(:),  ALLOCATABLE :: procflag

        LOGICAL, DIMENSION(2)               :: ConvergenceMASK

        PUBLIC :: proc_speed,cost
        PUBLIC :: min_phys,max_phys
        PUBLIC :: ghostsize,bcdef,decomp,assig
        PUBLIC :: ineighproc,nneighproc
        PUBLIC :: procflag,ConvergenceMASK
        !----------------------------------------------------------------------
        !  Particle positions and properties
        !----------------------------------------------------------------------
        INTEGER                             :: ghost_on_fire

        CLASS(ppm_t_particles_s_), POINTER  :: Part => NULL()
        INTEGER                             :: Npart
        INTEGER                             :: Mpart
        INTEGER                             :: NpartNew
        INTEGER                             :: NpartNewtmp

        INTEGER, DIMENSION(:),  ALLOCATABLE :: list_del_parts
        INTEGER                             :: del_parts

        PUBLIC :: ghost_on_fire
        PUBLIC :: Part,Npart,Mpart,NpartNew,NpartNewtmp
        PUBLIC :: list_del_parts,del_parts
        !----------------------------------------------------------------------
        !  Fields
        !----------------------------------------------------------------------
        CLASS(ppm_t_field_),       POINTER  :: image => NULL()
        !!! pixel intensities of the image, this field
        !!! is defined with INTEGER datatype

        CLASS(ppm_t_field_),       POINTER  :: labels => NULL()
        !!! pixel labels of the segmentation

        CLASS(ppm_t_field_),       POINTER  :: Backuplabels => NULL()
        !!! Backup pixel labels of the segmentation MCMC
        !!! A copy of the label image for the reconstruction of the
        !!! final results (marginals). It is a backup of the inital state.

        CLASS(ppm_t_field_),       POINTER  :: plabels => NULL()
        !!! particle labels of the segmentation

        CLASS(ppm_t_field_),       POINTER  :: pind => NULL()
        !!! index of the particle living in the mesh cell

        PUBLIC :: image,labels,plabels,pind
        !----------------------------------------------------------------------
        !  Normalized Image intensity information
        !----------------------------------------------------------------------
!         REAL(ppm_kind_double)               :: subimageintensity
!         REAL(ppm_kind_double)               :: imageintensity

        !-----------------------------------------------------------------------
        !  Which energy term to use?
        !-----------------------------------------------------------------------

        !-----------------------------------------------------------------------
        !  Weights of the energy terms
        !-----------------------------------------------------------------------
        REAL(MK)                            :: energy_coeff_data
        REAL(MK)                            :: energy_coeff_length
        REAL(MK)                            :: energy_coeff_balloon
        REAL(MK)                            :: energy_coeff_outward_flow
        REAL(MK)                            :: energy_region_merge_ths
        REAL(MK)                            :: energy_local_window_radius
        REAL(MK)                            :: energy_curvature_mask_radius
        !!! hypersphere radius which is used in curvature regularization

        PUBLIC :: energy_coeff_data,energy_coeff_length,energy_coeff_balloon
        PUBLIC :: energy_coeff_outward_flow,energy_region_merge_ths
        PUBLIC :: energy_local_window_radius,energy_curvature_mask_radius
        !-----------------------------------------------------------------------
        !  External and Internal energy names
        !-----------------------------------------------------------------------
        CHARACTER(LEN=ppm_char)             :: energy_ext_name
        CHARACTER(LEN=ppm_char)             :: energy_int_name

        REAL(MK), DIMENSION(:), ALLOCATABLE :: energya
        !!! Delta-E of the particles

        PUBLIC :: energy_ext_name,energy_int_name,energya

        !-----------------------------------------------------------------------
        !
        !-----------------------------------------------------------------------
        INTEGER, DIMENSION(:),  ALLOCATABLE :: labela
        ! label for this particle (what it label it is)
        INTEGER, DIMENSION(:),  ALLOCATABLE :: candlabela
        ! candidate label (l') for this particle (what it could be changed to)
        INTEGER, DIMENSION(:),  ALLOCATABLE :: ccandlabela
        ! count of parents with labels of candidate label (l')
        INTEGER, DIMENSION(:),  ALLOCATABLE :: ndaughtersa
        INTEGER, DIMENSION(:),  ALLOCATABLE :: nmothersa
        INTEGER, DIMENSION(:,:),ALLOCATABLE :: mothersa
        ! index (local!) of the mother particles
        INTEGER, DIMENSION(:,:),ALLOCATABLE :: daughtersa
        ! index (local!) of the daughter particles
        INTEGER, DIMENSION(:),  ALLOCATABLE :: accepteda
        INTEGER, DIMENSION(:),  ALLOCATABLE :: processeda

        PUBLIC :: labela,candlabela,ccandlabela
        PUBLIC :: ndaughtersa,nmothersa,mothersa,daughtersa
        PUBLIC :: accepteda,processeda
        !----------------------------------------------------------------------
        !  Image preprocessing parameters
        !----------------------------------------------------------------------
        REAL(MK)                            :: intensity_ths

        REAL(MK)                            :: m_OutsideValue
        REAL(MK)                            :: m_lowerBound
        REAL(MK)                            :: m_upperBound

        INTEGER                             :: histSize

        PUBLIC :: intensity_ths,m_OutsideValue
        PUBLIC :: m_lowerBound,m_upperBound
        PUBLIC :: histSize
        !----------------------------------------------------------------------
        !  Image parameters
        !----------------------------------------------------------------------
        !----------------------------------------------------------------------
        !  Number of pixels in each direction
        !----------------------------------------------------------------------
        INTEGER, DIMENSION(3)               :: Ngrid
        PUBLIC :: Ngrid
        !-----------------------------------------------------------------------
        !  Label generation
        !-----------------------------------------------------------------------
        ! this counter is local and constantly decreases when generating a new
        ! label.
        INTEGER                             :: loc_label
        PUBLIC :: loc_label
        !-----------------------------------------------------------------------
        !  Local variables
        !-----------------------------------------------------------------------
        REAL(MK), DIMENSION(:), ALLOCATABLE :: radius
        REAL(MK), DIMENSION(:), ALLOCATABLE :: space
        REAL(MK), DIMENSION(:), ALLOCATABLE :: rect
        REAL(MK), DIMENSION(:), ALLOCATABLE :: pixel

        PUBLIC :: radius,space,rect
        PUBLIC :: pixel

        REAL(MK), DIMENSION(:), ALLOCATABLE :: init_rd
        !!! radius user estimaton for blob detection
        REAL(MK), DIMENSION(:), ALLOCATABLE :: init_sp
        !!! radius user estimaton for blob detection

        PUBLIC :: init_rd,init_sp


        LOGICAL                             :: lNormalize
        REAL(MK)                            :: Normalfac

        INTEGER                             :: Tolerance
        REAL(MK)                            :: Tolerance_
        !!!

        PUBLIC :: lNormalize,Normalfac
        PUBLIC :: Tolerance,Tolerance_

        REAL(MK)                            :: Sigma
        !!!Sigma of the Gaussian kernel. */

        PUBLIC :: Sigma

        REAL(ppm_kind_double)               :: RegionSizeThreshold

        PUBLIC :: RegionSizeThreshold
        !-----------------------------------------------------------------------
        !  Which frame to read from a 4D image?
        !-----------------------------------------------------------------------
        !TODO
        !TOCHECK
        !still, it is not implemented
        !INTEGER                             :: frame
        !-----------------------------------------------------------------------
        !  Number of input images which from one image?
        !-----------------------------------------------------------------------
        INTEGER                             :: ninputimage
        INTEGER                             :: ninitimage
        LOGICAL                             :: createoneimage

        PUBLIC :: ninputimage,ninitimage,createoneimage
        !----------------------------------------------------------------------
        !  Global iteration number
        !----------------------------------------------------------------------
        INTEGER                             :: istep
        INTEGER                             :: maxiter
        !!! maximum number of iterations to perform
        PUBLIC :: istep,maxiter
        !----------------------------------------------------------------------
        !  I/O and control
        !----------------------------------------------------------------------
        CLASS(ppm_t_equi_mesh_),    POINTER :: iomesh => NULL()
        !!! mesh for I/O
        INTEGER                             :: iotopoid
        !!! topology id for I/O
        INTEGER                             :: iomeshid
        !!! mesh id for I/O
        INTEGER                             :: freqoutput
        INTEGER                             :: freqdiag

        PUBLIC :: iomesh,iotopoid,iomeshid,freqoutput,freqdiag

        CHARACTER(LEN=ppm_char)             :: inputimage
        CHARACTER(LEN=ppm_char)             :: abortfile
        CHARACTER(LEN=ppm_char)             :: outputfile
        CHARACTER(LEN=ppm_char)             :: outputname
        CHARACTER(LEN=ppm_char)             :: diagfile
        CHARACTER(LEN=ppm_char)             :: diagfmt
        CHARACTER(LEN=ppm_char)             :: casename
        CHARACTER(LEN=ppm_char)             :: initimage

        PUBLIC :: inputimage,abortfile,outputfile,outputname
        PUBLIC :: diagfile,diagfmt,casename,initimage

        INTEGER                             :: debug
        INTEGER                             :: ppm_log_unit

        PUBLIC :: debug,ppm_log_unit

        !debug time
        REAL(ppm_kind_double)               :: tmove_Simple
        REAL(ppm_kind_double)               :: tmove_SimpleComm
        REAL(ppm_kind_double)               :: tmove_NotSimple
        REAL(ppm_kind_double)               :: tmove_Part
        REAL(ppm_kind_double)               :: tmove_ghostfire

        PUBLIC :: tmove_Simple,tmove_SimpleComm,tmove_NotSimple
        PUBLIC :: tmove_Part,tmove_ghostfire
        !----------------------------------------------------------------------
        !  Program control flags
        !----------------------------------------------------------------------
        LOGICAL                             :: probeproc
        PUBLIC :: probeproc
        !----------------------------------------------------------------------
        !  Float comparison tolerance
        !----------------------------------------------------------------------
        REAL(MK)                            :: lmyeps
        PUBLIC :: lmyeps
        !----------------------------------------------------------------------
        !  MPI stuff
        !----------------------------------------------------------------------
        INTEGER                             :: ppm_rc_dim
        INTEGER                             :: ghost_size
        INTEGER                             :: rank
        INTEGER                             :: comm
        INTEGER                             :: n_io_procs_read  ! Number of I/O processors

        PUBLIC :: ppm_rc_dim,ghost_size,rank,comm,n_io_procs_read
        !----------------------------------------------------------------------
        !  Load-balancing stuff
        !----------------------------------------------------------------------
!         REAL(MK)                           :: step_time
!         INTEGER                            :: rebal_nstep=0

        !----------------------------------------------------------------------
        ! Public
        !----------------------------------------------------------------------
        PUBLIC :: define_args

      CONTAINS

        SUBROUTINE define_args

           IMPLICIT NONE
           !!! arg(variable,name,flag,long_flag,ctrl_name, &
           !!!     default,min,max,vtype,default_func,validator,help)
           !!!
           !!!
           !!! variable     :: Global variable to bind to.
           !!! name         :: Name of the arg for use in the auto generated usage message/ctrl file.
           !!! flag         :: Single character flag (eg. _'-f'_).
           !!! long_flag    :: Long flag (eg. _'--flag'_). Has to start with _'--'_!
           !!! ctrl_name    :: Control file variable name.
           !!! default      :: Default value.
           !!! min          :: Minimum value of arg.
           !!! max          :: Maximum value of arg.
           !!! vtype        :: Type of flag. Logical flags require no value to be supplied.
           !!! default_func :: Default function.
           !!! validator    :: Validator function.
           !!! help         :: Help string for the auto generated usage message/ctrl file.

           !-------------------------------------------------------------------------
           !  Run name
           !-------------------------------------------------------------------------
           CALL arg(casename,'casename',default='ppm_rc',ctrl_name='casename',help="Run casename")
           !-------------------------------------------------------------------------
           !  System dimensions
           !-------------------------------------------------------------------------
           CALL arg(ppm_rc_dim,'casedim',default=3,ctrl_name='casedim',help="System dimensions")
           !-------------------------------------------------------------------------
           !  Program control flags
           !-------------------------------------------------------------------------

           !-------------------------------------------------------------------------
           !  Number of I/O processors
           ! The number of reading processors. The default number of processors to read the data is 1.
           !-------------------------------------------------------------------------
           CALL arg(n_io_procs_read,'n_io_procs_read',default=1, &
           &    min=1,ctrl_name='n_io_procs_read',help="The number of reading processors")

           !-------------------------------------------------------------------
           !  Probe processor speeds at runtime?
           !-------------------------------------------------------------------
           CALL arg(probeproc,'probeproc',default=.FALSE.,ctrl_name='probeproc', &
           &    help="Probe processor speeds at runtime")
           !-------------------------------------------------------------------
           !  Digital topology (control)
           !-------------------------------------------------------------------
           CALL arg(AllowFusion,'allow_fusion',default=.TRUE.,   &
           &    ctrl_name='allow_fusion',help="Fusion is allowed?")
           CALL arg(AllowFusionZ,'allow_fusionz',default=.FALSE.,   &
           &    ctrl_name='allow_fusionz',help="Fusion in Z is allowed?")
           CALL arg(AllowFission,'allow_fission',default=.TRUE., &
           &    ctrl_name='allow_fission',help="Fission is allowed?")
           CALL arg(AllowHandles,'allow_handels',default=.TRUE., &
           &    ctrl_name='allow_handels',help="Handel is allowed?")


           !-------------------------------------------------------------------
           !  MCMC (control)
           !-------------------------------------------------------------------
           CALL arg(UseMCMC,'usemcmc',default=.FALSE.,   &
           &    ctrl_name='usemcmc',help="Algorithm is intended to use in MCMC mode for segmentation.")
           CALL arg(MCMCcontinue,'mcmccontinue',default=.FALSE.,   &
           &    ctrl_name='mcmccontinue',help="Continue MCMC run (no relabeling).")
           CALL arg(MCMCuseBiasedProposal,'mcmcusebiasedproposal',default=.FALSE.,   &
           &    ctrl_name='mcmcusebiasedproposal',help="Use biased proposal for MCMC sampling.")
           CALL arg(MCMCusePairProposal,'mcmcusepairproposal',default=.FALSE.,   &
           &    ctrl_name='mcmcusepairproposal',help="One move corresponds to a neighboring pair of particles.")

           CALL arg(MCMCstepsize,'mcmcstepsize',default=1, &
           &    min=1,ctrl_name='mcmcstepsize',help="The step size for Markov-chain Monte-Carlo method")

           CALL arg(MCMCtemperature,'mcmctemperature',default=1.0_MK, &
           &    ctrl_name='mcmctemperature',help="The Boltzmann temperature.")

           CALL arg(MCMCburnInFactor,'mcmcburninfactor',default=0.5_MK, &
           &    min=0.0_MK,max=1.0_MK,ctrl_name='mcmcburninfactor', &
           &    help="MCMC burn in factor, Length (factor) of burn-in phase (in [0,1[)..")

!            CALL arg(MCMCsampleOffBoundaryPercentage,'mcmcsampleoffboundarypercentage',default=0.05_MK, &
!            &    min=0.0_MK,max=1.0_MK,ctrl_name='mcmcsampleoffboundarypercentage', &
!            &    help="MCMC off-boundary samples percentage.")

!            CALL arg(MCMCFloatingParticlesProposalNormalizer,'mcmcfloatingparticlesproposalnormalizer',default=0.0_MK, &
!            &    min=0.0_MK,max=1.0_MK,ctrl_name='mcmcfloatingparticlesproposalnormalizer', &
!            &    help="MCMC off-boundary samples percentage.")











           !-------------------------------------------------------------------------
           !  Input filename (default is empty, thus no default is offered)
           !-------------------------------------------------------------------------
           CALL arg(inputimage,'inputimage',default='', &
           &    ctrl_name='inputimage',help="Input filename")
           CALL arg(ninputimage,'ninputimage',default=1, &
           &    ctrl_name='ninputimage',help="Number of image files which create the 3D image")
           CALL arg(createoneimage,'createoneimage',default=.FALSE., &
           &    ctrl_name='createoneimage',help="Control FLAG to create one TIFF file from several inputs")
           CALL arg(ninitimage,'ninitimage',default=1, &
           &    ctrl_name='ninitimage',help="Number of init image files which is used for the 3D image initialization")
           !-------------------------------------------------------------------
           !  Frame number to read from 4D image
           !-------------------------------------------------------------------
!            CALL arg(frame,'frame',default=1,ctrl_name='frame',help="Frame number to read from 4D image")

           !-------------------------------------------------------------------------
           !  Output filename (default is empty, thus no default is offered)
           !-------------------------------------------------------------------------
           CALL arg(outputfile,'outputfile',default='',ctrl_name='outputfile',help="Output filename")
           !-------------------------------------------------------------------------
           !  Output frequency
           !-------------------------------------------------------------------------
           CALL arg(freqoutput,'freqoutput',default=1,ctrl_name='freqoutput',help="Output frequency")
           !-------------------------------------------------------------------------
           !  Abort file
           !-------------------------------------------------------------------------
           CALL arg(abortfile,'abortfile',default='ABORT',ctrl_name='abortfile',help="Abort filename")
           !-------------------------------------------------------------------------
           !  Do 300 iterations per default
           !-------------------------------------------------------------------------
           CALL arg(maxiter,'maxiter',default=300,ctrl_name='maxiter',help="max number of iterations")
           !-------------------------------------------------------------------------
           !  Debug value
           !-------------------------------------------------------------------------
           CALL arg(debug,'debug',flag='-d',default=0,ctrl_name='debug',help="Debug level value")
           !-------------------------------------------------------------------
           !  Get diagnostics file
           !-------------------------------------------------------------------
           CALL arg(diagfile,'diagnosticsfile',default='',ctrl_name='diagnosticsfile', &
           &    help="diagnostics file")
           !-------------------------------------------------------------------
           !  Diagnostics file format
           !-------------------------------------------------------------------
           CALL arg(diagfmt,'diagnosticsfmt',ctrl_name='diagnosticsfmt',help="Diagnostics file format")
           !-------------------------------------------------------------------------
           !  Frequency of diagnogstics output
           !-------------------------------------------------------------------------
           CALL arg(freqdiag,'freqdiag',default=1,ctrl_name='freqdiag', &
           &    help="Frequency of diagnogstics output")





           !-------------------------------------------------------------------
           !  Energy term to use
           !-------------------------------------------------------------------
           CALL arg(energy_ext_name,'energy_ext_name',default='', &
           &    ctrl_name='energy_ext_name',help="External Energy term to use")
           CALL arg(energy_int_name,'energy_int_name',default='', &
           &    ctrl_name='energy_int_name',help="Internal Energy term to use")
           !-------------------------------------------------------------------
           !  weight for the energy
           !-------------------------------------------------------------------
           CALL arg(energy_coeff_data,'energy_coeff_data',default=1.0_MK, &
           &    ctrl_name='energy_coeff_data',help="coefficient for the energy")
           !-------------------------------------------------------------------
           !  weight for the energy
           !-------------------------------------------------------------------
           CALL arg(energy_coeff_length,'energy_coeff_length',default=0.04_MK, &
           & ctrl_name='energy_coeff_length',help="coefficient for the energy")
           !-------------------------------------------------------------------
           !  weight for the energy
           !-------------------------------------------------------------------
           CALL arg(energy_coeff_balloon,'energy_coeff_balloon',default=0.0_MK, &
           & ctrl_name='energy_coeff_balloon',help="coefficient for the energy")
           !-------------------------------------------------------------------
           !  weight for the energy
           !-------------------------------------------------------------------
           CALL arg(energy_coeff_outward_flow,'energy_coeff_outward_flow',default=0.0_MK, &
           & ctrl_name='energy_coeff_outward_flow',help="coefficient for the energy")
           !-------------------------------------------------------------------
           !  energy merging threshold coefficient
           !-------------------------------------------------------------------
           CALL arg(energy_region_merge_ths,'energy_region_merge_ths',default=0.2_MK, &
           &    ctrl_name='energy_region_merge_ths',help="merging threshold coefficient")
           CALL arg(energy_local_window_radius,'energy_local_window_radius',default=3.0_MK, &
           &    ctrl_name='energy_local_window_radius',help="local window radius")
           CALL arg(energy_curvature_mask_radius,'energy_curvature_mask_radius',default=4.0_MK, &
           &    ctrl_name='energy_curvature_mask_radius',help="hypersphere radius")


           CALL arg(OscillationThreshold,'Oscillation_ths',default=0.02_ppm_kind_double, &
           &    ctrl_name='Oscillation_ths',help="Oscillation threshold in (Convergence)")



           !-------------------------------------------------------------------
           !  Object initialization
           !-------------------------------------------------------------------
           CALL arg(init_mode,'init_mode',default='e_sphere',ctrl_name='init_mode', &
           &    help="Init term to use")
           !-------------------------------------------------------------------
           !  Object radius estimation
           !-------------------------------------------------------------------
           ALLOCATE(radius(3))
           CALL arg(radius,'init_rd',default=(/-1.0_MK,-1.0_MK,-1.0_MK/), &
           &    ctrl_name='init_rd',help="Init radius estimation")
           !-------------------------------------------------------------------
           !  Object dimension estimation
           !-------------------------------------------------------------------
           ALLOCATE(rect(3))
           CALL arg(rect,'init_rect',default=(/-1.0_MK,-1.0_MK,-1.0_MK/), &
           &    ctrl_name='init_rect',help="Init size estimation")
           !-------------------------------------------------------------------
           !  Object radius estimation
           !-------------------------------------------------------------------
           ALLOCATE(space(3))
           CALL arg(space,'init_sp',default=(/-1.0_MK,-1.0_MK,-1.0_MK/), &
           &    ctrl_name='init_sp',help="Init spacial distance estimation between objects")


           !----------------------------------------------------------------------
           !  Image preprocessing parameters
           !----------------------------------------------------------------------
           !-------------------------------------------------------------------
           !  Intensity threshold
           !-------------------------------------------------------------------
           CALL arg(intensity_ths,'intensity_ths',default=0.0_MK, &
           &    ctrl_name='intensity_ths',help="Intensity threshold")

           CALL arg(m_lowerBound,'lowerBound',default=-1.0_MK, &
           &    ctrl_name='lowerBound',help="Lower Intensity threshold")

           CALL arg(m_upperBound,'upperBound',default=-1.0_MK, &
           &    ctrl_name='upperBound',help="Upper Intensity threshold")

           CALL arg(histSize,'histsize',default=0, &
           &    ctrl_name='histsize',help="Histogram size (Histsize)")

           !-------------------------------------------------------------------
           !  Normalize the input image to the [0,1]
           !-------------------------------------------------------------------
           CALL arg(lNormalize,'Normalize',default=.TRUE., &
           &    ctrl_name='Normalize',help="Normalize the input image to the [0,1]")

           !-------------------------------------------------------------------
           !  Sigma of the gaussian kernel
           !-------------------------------------------------------------------
           CALL arg(Sigma,'Sigma',default=0.0_MK,ctrl_name='Sigma', &
           &    help="Sigma (standard deviation of the Gaussian kernel)")
           CALL arg(Tolerance,'Tolerance',default=1,ctrl_name='Tolerance', &
           &    help="Intensity Tolerance of the image")

           !-------------------------------------------------------------------------
           !  Input filename (default is empty, thus no default is offered)
           !-------------------------------------------------------------------------
           CALL arg(initimage,'initimage',default='',ctrl_name='initimage',help="Init image filename.")

           !-------------------------------------------------------------------
           !  Pixelsize
           !-------------------------------------------------------------------
           ALLOCATE(pixel(3))
           CALL arg(pixel(1),'pixelsize_x',default=1.0_MK,ctrl_name='pixelsize_x',help="Pixelsize")
           CALL arg(pixel(2),'pixelsize_y',default=1.0_MK,ctrl_name='pixelsize_y',help="Pixelsize")
           CALL arg(pixel(3),'pixelsize_z',default=1.0_MK,ctrl_name='pixelsize_z',help="Pixelsize")

           !-------------------------------------------------------------------
           !  Region threshold for removing Not Significant Regions
           !-------------------------------------------------------------------
           CALL arg(RegionSizeThreshold,'RegionSizeThreshold',default=1.0_ppm_kind_double, &
           &    ctrl_name='RegionSizeThreshold',help="Region threshold for removing Not Significant Regions")

        END SUBROUTINE define_args

      END MODULE ppm_rc_module_global
