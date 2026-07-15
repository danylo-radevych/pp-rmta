  ! Copyright (C) 2024-2026 Danylo Radevych
  !                                                                            
  ! This file is distributed under the terms of the MIT Non-AI License. 
  ! See the file `LICENSE' in the root directory of the               
  ! present distribution, or 
  ! https://github.com/non-ai-licenses/non-ai-licenses/blob/main/NON-AI-MIT .
  !
  ! Please cite: DOI: https://doi.org/10.1038/s41524-026-02141-7
  !
  !
  !=============================================================================
  PROGRAM rmta_prog
  !=============================================================================
  !!
  !! Program for rigid muffin-tin approximation (RMTA) factors
  !!
  !! Namelist name is &rmta.
  !!
  !
  !  Danylo Radevych
  !  updated: 2026/06/04
  !  started: 2024/07/20
  !
    !
    USE command_line_options, ONLY : npool_
    USE io_global, ONLY: stdout, ionode, ionode_id
    USE io_files,         ONLY : prefix, tmp_dir
    USE constants,        ONLY : rytoev
    USE kinds,            ONLY : DP
    USE io_global,        ONLY : ionode, ionode_id
    USE environment,      ONLY : environment_start, environment_end
    USE mp_world,         ONLY : mpime
    USE mp,               ONLY : mp_size, mp_bcast
    USE mp_global,        ONLY : mp_startup
    USE mp_images,        ONLY : intra_image_comm
    ! USE paw_variables,    ONLY : okpaw
    ! USE noncollin_module, ONLY : noncolin, lforcet
    ! USE control_flags,    ONLY : gamma_only
    USE klist,            ONLY : degauss, ngauss
    ! following modules needed for generation of tetrahedra
    USE ktetra,          ONLY : tetra_init
    USE symm_base, ONLY: nsym, s, time_reversal, t_rev
    USE klist, ONLY: nks, xk
    USE start_k, ONLY: k1, k2, k3, nk1, nk2, nk3
    USE cell_base, ONLY: at, bg
    USE lsda_mod, ONLY: lsda
    USE uspp_param, ONLY: upf
    ! custom module
    USE muffin_tin, ONLY: rmta_init, rmta_compute, rmta_quit
    USE mt_var, ONLY: &
      lmpi_single_rank, &
      luse_ref_pot, luse_tot_dos, &
      rmta_code, rmta_routine, & ! variables
      atomic_type, lsemiloc, lhybrid, &
      lsemilocupf, lnonlocal, &
      lwrite_dat, mt_ngauss, mt_degauss, &
      irf_delta, lrmt, &
      ltetra, ldense_r_grid, &
      rmt, rmt_method
    USE mt_printing, ONLY: &
      print_welcome_message, rmta_write, &
      print_clocks, check_input
    USE compare, ONLY: set_ref_pot, delete_ref_pot
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=256), EXTERNAL :: trimcheck
    !! trims names
    !
    CHARACTER(LEN=256) :: outdir
    !! outdir of rmta computation
    CHARACTER(len = 200) :: program_name
    !! name of this program
    LOGICAL :: verbose = .TRUE.
    !! set verbosity
    LOGICAL :: needwf = .TRUE.
    !! if true, read info about the wfc file
    INTEGER :: ios = 0
    !! system integer iostat
    INTEGER :: nks_aux
    !! auxiliary number of k-points needed for tetra_init
    INTEGER :: dnr
    !! number of points below the MT-radius to print out
    ! REAL(DP) :: rmt(natmax)
    ! !! muffin-tin radius
    !
    EXTERNAL :: errore, input_from_file, read_file_new, stop_pp
    !! error messages
    !
    NAMELIST / rmta / outdir, prefix, &
      rmt, lwrite_dat, ngauss, degauss, &
      dnr, lrmt, rmt, ltetra, lhybrid, rmt_method
    !
    ! defaults
    !
    program_name = "rmta.x"
    lmpi_single_rank = .TRUE. ! error if multiple ranks
    !
    lwrite_dat = .FALSE.
    lrmt = .FALSE.
    rmt_method = "touching"
    rmt(:) = -1.0_dp
    ngauss = -99
    degauss = 0.001_dp
    dnr = 0
    !
    ! initialize parallelization levels
    !
#if defined(__MPI)
    ! CALL mp_startup(start_images = .TRUE.)
    ! CALL mp_startup()
    CALL mp_startup(start_images = .TRUE.)
    !
    IF (lmpi_single_rank) THEN
      IF (mp_size(intra_image_comm) > 1) THEN
        WRITE(stdout, '(/5x, "Error: full MPI support for ", A, &
          " is not implemented yet.")') &
          TRIM(program_name)
        WRITE(stdout, '(5x, "Use: mpirun -n 1 ", A)') &
          TRIM(program_name)
        CALL errore(program_name, 'Use: mpirun -n 1', 1)
      END IF
    END IF
#endif
    !
    ! TODO possible solution to run on a single core
    ! CALL mp_startup(my_world_comm = MPI_COMM_SELF)
    !
    CALL print_welcome_message()
    !
    ! start environment management
    CALL environment_start(rmta_code)
    !
    ios = 0
    !
    ! set defaults for some variables in namelist
    prefix = 'pwscf'
    CALL get_environment_variable('ESPRESSO_TMPDIR', outdir)
    IF (trim(outdir) == ' ') outdir = './'
    !
    !
    ! stable defaults
    !
    ltetra = .TRUE.
    lhybrid = .FALSE.
    lnonlocal = .FALSE. ! non-local form is not functional at this point
    lsemiloc = .TRUE.
    lsemilocupf = .FALSE. ! will reset to .TRUE. if pseudo has PP_SEMILOCAL
    luse_ref_pot = .FALSE. ! only for debug with FLAPW potential
    luse_tot_dos = .TRUE.
    atomic_type = "na" ! obsolete
    ldense_r_grid = .FALSE.
    !
    !
    ! work on ionode
    IF (ionode) THEN
      ! read input from file
      CALL input_from_file()
      ! read namelist from file
      READ (5, nml = rmta, iostat = ios)
      tmp_dir = trimcheck(outdir)
      !
      IF (verbose) THEN
        ! check input
        CALL check_input(outdir)
      ENDIF
      !
      ! variables from input
      IF (.NOT. lsemiloc) &
        lsemilocupf = .FALSE.
      irf_delta = dnr
      mt_ngauss = ngauss
      mt_degauss = degauss
      !
    ENDIF
    !
    !
    !
    CALL mp_bcast(ios, ionode_id, intra_image_comm)
    IF (ios /= 0) CALL errore(rmta_routine, 'reading', rmta_routine, &
      'namelist', ABS(ios))
    !
    ! broadcast variables
    CALL mp_bcast(tmp_dir, ionode_id, intra_image_comm)
    CALL mp_bcast(prefix, ionode_id, intra_image_comm)
    CALL mp_bcast(atomic_type, ionode_id, intra_image_comm)
    CALL mp_bcast(lsemiloc, ionode_id, intra_image_comm)
    CALL mp_bcast(lsemilocupf, ionode_id, intra_image_comm)
    CALL mp_bcast(lhybrid, ionode_id, intra_image_comm)
    CALL mp_bcast(lnonlocal, ionode_id, intra_image_comm)
    CALL mp_bcast(lrmt, ionode_id, intra_image_comm)
    CALL mp_bcast(rmt_method, ionode_id, intra_image_comm)
    CALL mp_bcast(rmt, ionode_id, intra_image_comm)
    CALL mp_bcast(lwrite_dat, ionode_id, intra_image_comm)
    CALL mp_bcast(mt_ngauss, ionode_id, intra_image_comm)
    CALL mp_bcast(mt_degauss, ionode_id, intra_image_comm)
    CALL mp_bcast(luse_ref_pot, ionode_id, intra_image_comm)
    CALL mp_bcast(luse_tot_dos, ionode_id, intra_image_comm)
    CALL mp_bcast(irf_delta, ionode_id, intra_image_comm)
    CALL mp_bcast(ltetra, ionode_id, intra_image_comm)
    CALL mp_bcast(ldense_r_grid, ionode_id, intra_image_comm)
    !
    ! read xml data file produced by pw.x or cp.x
    CALL read_file_new(needwf)
    !
    IF ( ALL(upf%typ(:) == "SL") ) THEN
      lsemilocupf = .TRUE.
      CALL mp_bcast(lsemilocupf, ionode_id, intra_image_comm)
    END IF
    !
    !
    ! info on tetrahedra is no longer saved to file and must be rebuilt
    !
    ! in the lsda case, only the first half of the k points
    ! are needed in the input of "tetrahedra"
    !
    !
    IF (ltetra) THEN
      !
      IF (nk1 * nk2 * nk3 == 0) THEN
        CALL errore(program_name, &
          'tetrahedra integration only with automatic ' // &
          'Monkhorst-Pack k_point meshes.', 1)
      END IF
      !
      IF (lsda) THEN
        nks_aux = nks / 2
      ELSE
        nks_aux = nks
      END IF
      !
      WRITE(stdout,'(/5x,"Tetrahedra used"/)')
      CALL tetra_init(nsym, s, time_reversal, t_rev, at, bg, nks, &
        k1, k2, k3, nk1, nk2, nk3, nks_aux, xk)
      !
    END IF ! ltetra
    !
    !
    IF (luse_ref_pot) &
      CALL set_ref_pot(rmt(1))
    !
    ! init RMTA
    CALL rmta_init()
    !
    ! compute RMTA
    CALL rmta_compute()
    !
    !
    IF (ionode .AND. lwrite_dat) THEN
      ! print RMTA quantities
      CALL rmta_write()
    END IF
    !
    ! deallocate rmta-specific arrays
    CALL rmta_quit()
    !
    IF (luse_ref_pot) &
      CALL delete_ref_pot()
    !
    ! print RMTA clocks
    CALL print_clocks()
    !
    ! end environment management
    CALL environment_end(rmta_code)
    !
    ! synchronize processes before stopping
    CALL stop_pp
    !
  !=============================================================================
  END PROGRAM rmta_prog
  !=============================================================================
