
      SUBROUTINE ppm_rc_init_arg(info)

        !-------------------------------------------------------------------------
        !  Modules
        !-------------------------------------------------------------------------
        USE ppm_module_data
        USE ppm_module_mpi
        USE ppm_module_init, ONLY : ppm_init
        USE ppm_module_loadbal, ONLY : ppm_estimate_procspeed, &
        &   ppm_set_proc_speed

#ifdef __Linux
        USE ppm_rc_module_util, ONLY : ppm_rc_uppercase,ppm_rc_mem_usage,valueRSS0
#else
        USE ppm_rc_module_util, ONLY : ppm_rc_uppercase
#endif
        USE ppm_rc_module_topologicalnumber, ONLY : FG_ConnectivityType, &
        &   BG_ConnectivityType,BackgroundConnectivity,TopologicalNumberFunction
        USE ppm_rc_module_energy, ONLY : E_Gamma,E_ContourLengthApprox,E_PC, &
        &   E_PCGaussian,E_PCPoisson,E_PS,E_PSGaussian,E_PSPoisson,e_data,e_length
        USE ppm_rc_module_read, ONLY : ppm_rc_read_image_info
        USE ppm_rc_module_rnd, ONLY : saru_random_init,mt_random_init
        USE ppm_rc_module_mcmc, ONLY : MCMCCheckParameters
        IMPLICIT NONE

        !-------------------------------------------------------------------------
        !  Includes
        !-------------------------------------------------------------------------
        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        INTEGER, INTENT(  OUT) :: info

        REAL(ppm_kind_double) :: t0

        CHARACTER(LEN=ppm_char) :: caller="ppm_rc_init_arg"

#ifdef __Linux
        INTEGER               :: memory
#endif
        INTEGER               :: tolexp
        INTEGER               :: iopt
        INTEGER, DIMENSION(1) :: ldl,ldu
        !-------------------------------------------------------------------------
        !  Initialize
        !-------------------------------------------------------------------------
#ifdef __MPI
        !-------------------------------------------------------------------------
        !  Initialize MPI
        !-------------------------------------------------------------------------
        CALL MPI_Init(info)
        or_fail_MPI("MPI_Init failed.")

        comm=MPI_COMM_WORLD

        !-------------------------------------------------------------------------
        !  Get rank of processor
        !-------------------------------------------------------------------------
        CALL MPI_Comm_rank(comm,rank,info)
        or_fail_MPI("MPI_Comm_Rank failed.")
#endif

#ifdef __Linux
        valueRSS0=0
        CALL ppm_rc_mem_usage(memory,info)
        IF (memory.GE.0) valueRSS0=memory
#endif

        CALL ppm_util_time(t0)

        CALL define_args
        CALL parse_args(info)
        SELECT CASE (info)
        CASE (0)
        CASE (exit_gracefully)
           GOTO 9999
        CASE DEFAULT
           or_fail("parse_args",ppm_error=ppm_error_fatal)
        END SELECT
        CALL reset()

        !-------------------------------------------------------------------------
        !  Define defaults
        !-------------------------------------------------------------------------
        CALL ppm_rc_defaults(info)
        or_fail('ppm_rc_defaults failed')

        !-------------------------------------------------------------------------
        !  Check parameter validity
        !-------------------------------------------------------------------------
        CALL ppm_rc_check_ctrl(info)
        or_fail('ppm_rc_check_ctrl failed')

        !-------------------------------------------------------------------------
        !  Initialize the ppm library
        !-------------------------------------------------------------------------
        !TODO
        !TOCHECK
        tolexp = INT(LOG10(EPSILON(lmyeps)))+3

        CALL ppm_init(ppm_rc_dim,MK,tolexp,comm,debug,info,ppm_log_unit)
        or_fail('Failed to initialize PPM library.')

        rank=ppm_rank
        comm=ppm_comm

        !-------------------------------------------------------------------------
        ! Set geometry tolerance
        !-------------------------------------------------------------------------
        lmyeps=MERGE(REAL(ppm_myepsd,MK),REAL(ppm_myepss,MK),MK.EQ.ppm_kind_double)

        !-------------------------------------------------------------------------
        !  Initialize local label counter
        !  This way, we can initialize local labels on every processor
        !  from a big number and will reduce them one by one.
        !-------------------------------------------------------------------------
        IF (ppm_nproc.LE.1024) THEN
           !Max number of regions in one processor is limited to 262144
           loc_label=ISHFT(rank+1,SHIFT=18)
        ELSE
           !Max number of regions in one processor is limited to 131072
           loc_label=ISHFT(rank+1,SHIFT=17)
        ENDIF

        !-------------------------------------------------------------------------
        !  Initialize topology related members
        !-------------------------------------------------------------------------
        ALLOCATE(FG_ConnectivityType,STAT=info)
        or_fail_alloc("FG_ConnectivityType")

        CALL FG_ConnectivityType%create(ppm_rc_dim,ppm_rc_dim*2,info)
        or_fail("FG_ConnectivityType%create")

        BG_ConnectivityType => BackgroundConnectivity(FG_ConnectivityType)

        CALL TopologicalNumberFunction%create(FG_ConnectivityType,info)
        or_fail("TopologicalNumberFunction%create")

        !TOCHECK
        !TODO FIXME
        TopologicalNumberFunction%ComputeBackgroundTN=.FALSE.

        !-------------------------------------------------------------------------
        !  Initialize acceptance reduction factor (for Oscillations)
        !  Set ConvergenceMASK to TRUE which is used for global synchronization
        !-------------------------------------------------------------------------
        AcceptedPointsFactor=one
        ConvergenceMASK=.TRUE.

        !-------------------------------------------------------------------------
        !  Initialize init kind
        !-------------------------------------------------------------------------
        CALL ppm_rc_uppercase(init_mode)

        SELECT CASE (TRIM(init_mode))
        CASE ("E_RECT","RECT")
           vInitKind=e_rect

        CASE ("E_SPHERE","SPHERE")
           vInitKind=e_sphere

        CASE ("E_OTSU","OTSU")
           vInitKind=e_otsu

        CASE ("E_LOCALMAX","LOCALMAX")
           vInitKind=e_localmax

        CASE DEFAULT
           vInitKind=e_localmax

        END SELECT

        !-------------------------------------------------------------------------
        !  Initialize energy related terms
        !-------------------------------------------------------------------------
        CALL ppm_rc_uppercase(energy_ext_name)
        CALL ppm_rc_uppercase(energy_int_name)
        !-------------------------------------------------------------------------
        !  Initialize external energy related terms
        !-------------------------------------------------------------------------
        SELECT CASE (TRIM(energy_ext_name))
        CASE ("PC")
           ALLOCATE(E_PC::e_data,STAT=info)
           or_fail_alloc("e_data")

           CALL e_data%create(1,energy_coeff_data,info,energy_region_merge_ths)
           or_fail("e_data%create")

        CASE ("PCGAUSSIAN")
           ALLOCATE(E_PCGaussian::e_data,STAT=info)
           or_fail_alloc("e_data")

           CALL e_data%create(2,energy_coeff_data,info,energy_region_merge_ths)
           or_fail("e_data%create")

        CASE ("PCPOISSON")
           ALLOCATE(E_PCPoisson::e_data,STAT=info)
           or_fail_alloc("e_data")

           CALL e_data%create(3,energy_coeff_data,info,energy_region_merge_ths)
           or_fail("e_data%create")

        CASE ("PS")
           ALLOCATE(E_PS::e_data,STAT=info)
           or_fail_alloc("e_data")

           CALL e_data%create(11,energy_coeff_data,info,energy_region_merge_ths)
           or_fail("e_data%create")

           SELECT TYPE (e_data)
           TYPE IS (E_PS)
              e_data%m_Radius=energy_local_window_radius

              SELECT CASE (ppm_rc_dim)
              CASE (2)
                 CALL e_data%PrepareEnergyCalculation_2d(info)
              CASE (3)
                 CALL e_data%PrepareEnergyCalculation_3d(info)
              END SELECT
              or_fail("e_data%PrepareEnergyCalculation")
           END SELECT

        CASE ("PSGAUSSIAN")
           ALLOCATE(E_PSGaussian::e_data,STAT=info)
           or_fail_alloc("e_data")

           CALL e_data%create(12,energy_coeff_data,info,energy_region_merge_ths)
           or_fail("e_data%create")

           SELECT TYPE (e_data)
           TYPE IS (E_PSGaussian)
              e_data%m_Radius=energy_local_window_radius

              SELECT CASE (ppm_rc_dim)
              CASE (2)
                 CALL e_data%PrepareEnergyCalculation_2d(info)
              CASE (3)
                 CALL e_data%PrepareEnergyCalculation_3d(info)
              END SELECT
              or_fail("e_data%PrepareEnergyCalculation")
           END SELECT

        CASE ("PSPOISSON")
           ALLOCATE(E_PSPoisson::e_data,STAT=info)
           or_fail_alloc("e_data")

           CALL e_data%create(13,energy_coeff_data,info,energy_region_merge_ths)
           or_fail("e_data%create")

           SELECT TYPE (e_data)
           TYPE IS (E_PSPoisson)
              e_data%m_Radius=energy_local_window_radius

              SELECT CASE (ppm_rc_dim)
              CASE (2)
                 CALL e_data%PrepareEnergyCalculation_2d(info)
              CASE (3)
                 CALL e_data%PrepareEnergyCalculation_3d(info)
              END SELECT
              or_fail("e_data%PrepareEnergyCalculation")
           END SELECT

        END SELECT

        IF (energy_coeff_balloon.GT.smallest.OR. &
        &   energy_coeff_balloon.LT.-smallest) THEN
           e_data%m_EnergyFunctional=e_data%m_EnergyFunctional+1000
        ENDIF

        !-------------------------------------------------------------------------
        !  Initialize internal energy related terms
        !-------------------------------------------------------------------------
        SELECT CASE (TRIM(energy_int_name))
        CASE ("GAMMA")
           ALLOCATE(E_Gamma::e_length,STAT=info)
           or_fail_alloc("e_length")

           CALL e_length%create(31,energy_coeff_length,info)
           or_fail("e_length%create")

        CASE ("CURV")
           ALLOCATE(E_ContourLengthApprox::e_length,STAT=info)
           or_fail_alloc("e_length")

           CALL e_length%create(32,energy_coeff_length,info)
           or_fail("e_length%create")

           SELECT TYPE (e_length)
           TYPE IS (E_ContourLengthApprox)
              e_length%m_Radius=energy_curvature_mask_radius

              SELECT CASE (ppm_rc_dim)
              CASE (2)
                 CALL e_length%PrepareEnergyCalculation_2d(info)
              CASE (3)
                 CALL e_length%PrepareEnergyCalculation_3d(info)
              END SELECT
              or_fail("e_length%PrepareEnergyCalculation")
           END SELECT

        END SELECT

        IF (energy_coeff_outward_flow.GT.smallest.OR. &
        &   energy_coeff_outward_flow.LT.-smallest) THEN
           e_length%m_EnergyFunctional=e_length%m_EnergyFunctional+1000
        ENDIF

        !-------------------------------------------------------------------------
        !  Get image information and store the pixel size
        !-------------------------------------------------------------------------
        CALL ppm_rc_read_image_info(inputimage,ninputimage,info)
        or_fail('Failed to read image information.')

        max_phys(1:ppm_rc_dim) = REAL(Ngrid(1:ppm_rc_dim)-1,MK)  !+smallest

        !-------------------------------------------------------------------------
        !  Write some of the parameters
        !-------------------------------------------------------------------------
        CALL ppm_rc_dump_parameters(info)
        or_fail("ppm_rc_dump_parameters")

        !-------------------------------------------------------------------------
        !  Probe the processor speeds if needed
        !-------------------------------------------------------------------------
        IF (probeproc) THEN
           iopt   = ppm_param_alloc_fit
           ldl(1) = 0
           ldu(1) = ppm_nproc-1
           CALL ppm_alloc(proc_speed,ldl,ldu,iopt,info)
           or_fail_alloc('Failed to allocate proc_speed.')

           CALL ppm_estimate_procspeed(proc_speed,info)
           or_fail("Failed in ppm_estimate_procspeed")

           CALL ppm_set_proc_speed(proc_speed,info)
           or_fail("ppm_set_proc_speed")
        ENDIF

        SELECT CASE (UseMCMC)
        CASE (.TRUE.)
           CALL MCMCCheckParameters(info)
           or_fail("MCMCCheckParameters")

           CALL saru_random_init(info)
           or_fail("saru_random_init")

           CALL mt_random_init(info)
           or_fail("mt_random_init")

           SELECT CASE (lNormalize)
           CASE (.TRUE.)
              IF (rank.EQ.0) THEN
                 stdout("Warning: MCMC on normalized image!")
              ENDIF
           END SELECT
        END SELECT

        !-------------------------------------------------------------------------
        !  Return
        !-------------------------------------------------------------------------
      9999 CONTINUE
        CALL substop(caller,t0,info)
        RETURN
      END SUBROUTINE ppm_rc_init_arg
