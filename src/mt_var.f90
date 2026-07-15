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
  MODULE mt_var
  !=============================================================================
  !!
  !! Module defining and (de)allocating RMTA variables
  !!
  !=============================================================================
  !
  !  Danylo Radevych
  !
    USE kinds, ONLY: DP
    USE radial_grids, ONLY: radial_grid_type
    USE constants, ONLY: bohr_radius_si
    !
    IMPLICIT NONE
    !
    PUBLIC :: &
      rmta_set_vars, rmta_delete_vars, rmt_default, &
      lmpi_single_rank, &
      mt_prec, natoms, nspins, mt_rmt, norbs, &
      rmta_lmax, orb_label, &
      tau_cart, n_chem_types, natoms_per_chem_type, &
      fermi_energy, &
      mt_nr, mt_nr_max, mt_r, chir, betar, vlocionr, vlocaer, &
      vlocionr3d, &
      vlocscr00rf, vlocscf00rf, &
      vlociong, vlocscrr3d, vlocscfr3d, &
      vlocscrg3d, mt_g, vlocscfg3d, &
      vsemilocr, vsemilocrf, urf, duderf, dudrrf, &
      d2udrderf, &
      mll1rf_label, mll1rf, &
      vfullrf, rvfullrf, &
      irf_min, irf_max, &
      dos_nlmrf, dos_nlrf, dos_nrf, dos_n, &
      etall1rf, &
      loglrf, dloglderf, &
      luse_ref_pot, luse_tot_dos, &
      rmta_routine, rmta_code, &
      atomic_type, &
      lsemiloc, lsemilocupf, &
      lhybrid, lnonlocal, &
      lwrite_dat, &
      mt_ngauss, mt_degauss, &
      mt_rab, &
      nchis, chi_label, &
      nbetas, beta_label, &
      mt_grid, &
      mt_nrf, mt_rf, mt_rfab, mt_dx, rmta_ng, lrmt, &
      ltetra, ldense_r_grid, rmt, rmt_method
    !
    !
    PRIVATE :: set_orbitals, set_tau_cart, set_rmt, set_chem_type, &
      set_fermi_energy, set_grids, set_upf_vars, set_scf_vars, set_empty_arrays
    !
    !
    INTEGER, PARAMETER :: natmax = 256
    !! max number of atoms
    REAL(DP), PARAMETER :: bohrtoang = bohr_radius_si * 1.0E10_DP
    !! Bohr radius in angstroms
    !
    CHARACTER(LEN=12), SAVE :: rmta_code = 'RMTA'
    !! Name of the code
    CHARACTER(LEN=12), SAVE :: rmta_routine = 'rmta'
    !! Name of the code (short)
    CHARACTER(len = 256) :: atomic_type
    !! TODO: obsolete
    CHARACTER(LEN=128) :: rmt_method
    !! Method for the selection of MT radii
    CHARACTER(LEN=1), ALLOCATABLE :: orb_label(:)
    !! orbital labels
    CHARACTER(LEN=5), ALLOCATABLE :: mll1rf_label(:, :, :)
    !! mll1rf_label(norbs, nspins, natoms)
    !! labels for mll1rf
    CHARACTER(LEN=2), ALLOCATABLE :: beta_label(:, :)
    !! beta_label(nbetas(n_chem_types), n_chem_types)
    !! label of beta(r) for each atomic wfc and type
    CHARACTER(LEN=2), ALLOCATABLE :: chi_label(:, :)
    !! chi_label(nchis(n_chem_types), n_chem_types)
    !! label of chi(r) for each atomic wfc and type
    LOGICAL :: ltetra
    !! if true, use tetrahedron method for integration
    LOGICAL :: ldense_r_grid
    !! if true, use high number of points on radial grid
    LOGICAL :: lwrite_dat
    !! if true, write data in *.dat files
    LOGICAL :: lnonlocal
    !! use nonlocal operator V_{NL}(r, r')
    LOGICAL :: luse_ref_pot
    !! use reference (AE) potential instead of PS
    INTEGER :: mt_ngauss
    !! type of the delta-function
    INTEGER :: norbs
    !! number of orbitals for RMTA purposes
    !! 4 <-> s, p, d, f
    LOGICAL :: lsemiloc
    !! Include V^l_SL(r) pieces
    LOGICAL :: lsemilocupf
    !! Read V^l_SL(r) pieces from UPF
    LOGICAL :: lhybrid
    !! use hybrid Pettifor formula
    ! REAL(DP), ALLOCATABLE :: rmta_pse(:)
    ! !! energy corresponding to V_{SL}
    LOGICAL :: luse_tot_dos
    !! if true, use total DOS per spin,
    !! not the one inside each atom
    LOGICAL :: lrmt
    !! set rmt from file
    LOGICAL :: lmpi_single_rank
    !! if true, give an error when multiple ranks are used
    INTEGER :: rmta_lmax
    !! maximum l for RMTA
    INTEGER :: natoms
    !! number of atoms
    INTEGER :: nspins
    !! number of spins
    INTEGER :: n_chem_types
    !! number of chemical types
    INTEGER :: mt_nr_max
    !! max size in the mt_nr array
    INTEGER, ALLOCATABLE :: natoms_per_chem_type(:)
    !! number of atoms per each chemical type
    INTEGER, ALLOCATABLE:: mt_nr(:)
    !! mt_nr(n_chem_types)
    !! number of points on radial mesh for each type
    INTEGER :: rmta_ng
    !! number of \bm{G} vectors
    INTEGER :: irf_min
    !! min irf
    INTEGER :: irf_max
    !! max irf
    INTEGER :: irf_delta
    !! number of rf points below rmt
    INTEGER :: mt_nrf
    !! number of points for each radial grid
    INTEGER, ALLOCATABLE :: nchis(:)
    !! nchis(n_chem_types)
    !! number of atomic orbitals chi(r) for each type
    INTEGER, ALLOCATABLE :: nbetas(:)
    !! nbetas(n_chem_types)
    !! number of of projectors beta(r) for each type
    REAL(DP) :: mt_degauss
    !! degauss value for the delta-function
    REAL(DP) :: mt_prec
    !! default rmta precision
    REAL(DP) :: rmt(natmax)
    !! muffin-tin radii for each atom
    REAL(DP), ALLOCATABLE :: mt_dx(:)
    !! dx parameter for RMTA grids
    REAL(DP), ALLOCATABLE :: mt_rmt(:)
    !! MT radii for each symmetry type
    REAL(DP), ALLOCATABLE :: vlocionr3d(:)
    !! vlocionr3d(dfftp%nnr)
    !! total local ion potential in 3D real space
    REAL(DP), ALLOCATABLE :: fermi_energy(:)
    !! Fermi energy for each spin
    REAL(DP), ALLOCATABLE :: dos_n(:)
    !! total densities n^i(E_F)
    !! dos_n(nspins), per spin
    REAL(DP), ALLOCATABLE :: tau_cart(:, :)
    !! atomic positions in cartesian coordinates, bohr
    REAL(DP), ALLOCATABLE :: mt_r(:, :)
    !! mt_r(mt_nr(n_chem_types), n_chem_types)
    !! radial grid for each type's pseudo
    REAL(DP), ALLOCATABLE :: mt_rab(:, :)
    !! mt_rab(mt_nr(n_chem_types), n_chem_types)
    !! radial dr / di for each type's pseudo
    REAL(DP), ALLOCATABLE :: mt_rf(:, :)
    !! RMTA log grids mt_rf(mt_nrf, nst)
    REAL(DP), ALLOCATABLE :: mt_rfab(:, :)
    !! RMTA log grids mt_rfab(mt_nrf, nst)
    REAL(DP), ALLOCATABLE :: vlociong(:, :)
    !! vlociong(ngl, ntyp)
    !! copy of local pseudopotential
    !! V_ps(G) = V_{loc}^{ion}(|G|) for all types of atoms
    REAL(DP), ALLOCATABLE :: vlocscrr3d(:, :)
    !! vlocscrr3d(dfftp%nnr, nspins)
    !! total local potential on real space (smooth grid)
    !! see PW/set_vrs.f90 for more details
    REAL(DP), ALLOCATABLE :: vlocscfr3d(:, :)
    !! vlocscfr3d(dfftp%nnr, nspins)
    !! total local SCF part of the potential [scf(H+xc)]
    !! on real space (smooth grid)
    !! see PW/set_vrs.f90 for more details
    REAL(DP), ALLOCATABLE :: mt_g(:, :)
    !! mt_g(4, ngm)
    !! coordinates (1 : 3) and magnitude (4) of G-vectors
    REAL(DP), ALLOCATABLE :: chir(:, :, :)
    !! chir(mt_nr, nchis, n_chem_types)
    !! chi(r) for each chemical type
    REAL(DP), ALLOCATABLE :: betar(:, :, :)
    !! betar(mt_nr, nbetas, n_chem_types)
    !! beta(r) for each chemical type
    REAL(DP), ALLOCATABLE :: vlocionr(:, :)
    !! vlocionr(mt_nr, n_chem_types)
    !! V_{loc}^{ion}(r) for each chemical type
    REAL(DP), ALLOCATABLE :: vlocaer(:, :)
    !! vlocaer(mt_nr, n_chem_types)
    !! V_{loc}^{AE}(r) for each type
    REAL(DP), ALLOCATABLE :: vlocscr00rf(:, :, :)
    !! vlocscr00rf(mt_nrf, nspins, natoms)
    !! V_{L=00}(|r|) on fine grid for each atom and spin
    REAL(DP), ALLOCATABLE :: vlocscf00rf(:, :, :)
    !! vlocscf00rf(mt_nrf, nspins, natoms)
    !! V_{loc}^{scf}_{L=00}(|r|) on fine grid for each spin
    REAL(DP), ALLOCATABLE :: vsemilocr(:, :, :)
    !! vsemilocr(mt_nr(n_chem_types), norbs, &
    !! n_chem_types)
    !! semilocal V_SL(r) = beta(r) / chi(r)
    REAL(DP), ALLOCATABLE :: vsemilocrf(:, :, :)
    !! vsemilocrf(mt_nrf, norbs, natoms)
    !! semilocal V_SL(rf) interpolated from V_SL(r)
    REAL(DP), ALLOCATABLE :: urf(:, :, :, :)
    !! urf(mt_nrf, norbs, nspins, natoms)
    !! computed radial functions u(r) on fine grid
    REAL(DP), ALLOCATABLE :: duderf(:, :, :, :)
    !! duderf(mt_nrf, norbs, nspins, natoms)
    !! computed energy derivatives of radial functions
    !! d u(r, e) / de on fine r grid
    REAL(DP), ALLOCATABLE :: dudrrf(:, :, :, :)
    !! dudrrf(mt_nrf, rnorbs, nspins, natoms)
    !! computed r derivatives of radial functions
    !! d u(r, e) / dr on fine r grid
    REAL(DP), ALLOCATABLE :: d2udrderf(:, :, :, :)
    !! d2udrderf(mt_nrf, norbs, nspins, natoms)
    !! computed second derivatives of radial functions
    !! d2 u(r, e) / (dr de) on fine r grid
    REAL(DP), ALLOCATABLE :: mll1rf(:, :, :, :)
    !! mll1rf(mt_nrf, norbs, nspins, natoms)
    !! computed matrix elements M_{l, l + 1},
    !! as defined by Pettifor,
    !! on fine r grid
    REAL(DP), ALLOCATABLE :: vfullrf(:, :, :, :)
    !! vfullrf(mt_nrf, norbs, nspins, natoms)
    !! total potential V(r) on fine grid
    REAL(DP), ALLOCATABLE :: rvfullrf(:, :, :, :)
    !! rvfullrf(mt_nrf, norbs, nspins, natoms)
    !! total potential r V(r) on fine grid
    REAL(DP), ALLOCATABLE :: dos_nlmrf(:, :, :, :)
    !! partial densities n^i_{lm}(r, E_F)
    !! dos_nlmrf(nr, (lmax + 1)**2, norbs, nspins, natoms)
    REAL(DP), ALLOCATABLE :: dos_nlrf(:, :, :, :)
    !! partial densities n^i_{l}(r, E_F)
    !! dos_nlrf(nr, norbs, nspins, natoms)
    REAL(DP), ALLOCATABLE :: dos_nrf(:, :, :)
    !! total partial densities n^i(r, E_F)
    !! dos_nrf(nr, nspins, natoms), per atom, per spin
    REAL(DP), ALLOCATABLE :: etall1rf(:, :, :, :)
    !! etall1rf(mt_nrf, norbs, nspins, natoms)
    !! McMillan-Hopfield \eta_l,
    !! on fine r grid
    !! note etall1rf(:, norbs, :, :) = rmta_etarf (total eta)
    REAL(DP), ALLOCATABLE :: loglrf(:, :, :, :)
    !! loglrf(mt_nrf, norbs, nspins, natoms)
    !! computed log derivatives of radial functions
    !! L_l(r) = r u_l'(r) / u_l(r) - 1 on fine r grid
    REAL(DP), ALLOCATABLE :: dloglderf(:, :, :, :)
    !! dloglderf(mt_nrf, norbs, nspins, natoms)
    !! computed energy derivatives of log derivatives of radial functions
    !! d L_l(r, e) / d e
    COMPLEX(DP), ALLOCATABLE :: vlocscrg3d(:, :)
    !! vlocscrg3d(dfftp%nng, nspins)
    !! total local potential on reciprocal space (smooth grid)
    !! \tilde{V}(\bm{G}) = \int d\bm{r} V(\bm{r}) e^{i \bm{G} \bm{r}}
    COMPLEX(DP), ALLOCATABLE :: vlocscfg3d(:, :)
    !! vlocscfg3d(dfftp%nng, nspins)
    !! total local SCF part of the potential [scf(H+xc)]
    !! on reciprocal space (smooth grid)
    !! \tilde{V}_{scf}(\bm{G}) = \int d\bm{r} V(\bm{r}) e^{i \bm{G} \bm{r}}
    TYPE(radial_grid_type) :: mt_grid
    !! RMTA log grid
    !
  !
  CONTAINS
  !
    !
    !---------------------------------------------------------------------------
    SUBROUTINE rmta_set_vars()
    !---------------------------------------------------------------------------
    !!
    !! Sets some RMTA variables.
    !!
    !---------------------------------------------------------------------------
      USE lsda_mod, ONLY: nspin
      USE ions_base, ONLY: nat
      USE gvect, ONLY: ngm
      !
      IMPLICIT NONE
      !
      CHARACTER(len=256) :: routine_name
      !! name of this subroutine
      !
      routine_name = "rmta_set_vars"
      rmta_routine = "rmta"
      !
      !
      natoms = nat
      nspins = nspin
      mt_prec = 1.0e-12_dp
      !
      IF (ldense_r_grid) THEN
        mt_nrf = 3001
      ELSE
        mt_nrf = 401
      END IF
      !
      rmta_ng = ngm
      !
      ! irf_delta = 0
      irf_max = mt_nrf
      irf_min = irf_max - irf_delta
      !
      ! MT radii for each symmetry type
      CALL set_rmt()
      !
      ! radial grid for each symmetry type
      CALL set_grids()
      !
      ! orbitals
      CALL set_orbitals()
      !
      ! atomic positions in bohrs
      CALL set_tau_cart()
      !
      ! set chemical type info
      CALL set_chem_type()
      !
      ! Fermi energy
      CALL set_fermi_energy()
      !
      ! set upf variables on radial grid from upf
      CALL set_upf_vars()
      !
      ! set (n)scf vars
      CALL set_scf_vars()
      !
      ! allocate other RMTA arrays
      CALL set_empty_arrays()
      !
      !
    !---------------------------------------------------------------------------
    END SUBROUTINE rmta_set_vars
    !---------------------------------------------------------------------------
    !
    !
    !---------------------------------------------------------------------------
    SUBROUTINE set_rmt()
    !---------------------------------------------------------------------------
    !!
    !! Sets MT radii for each symmetry type as a minimum of the 
    !! nearest-neighbor distance divided in ratio of internal MT radii 
    !! defaults, then ensuring that the spheres do touch...
    !! ... more in the code below
    !!
    !---------------------------------------------------------------------------
    !
    !  D. Radevych
    !
      !
      USE io_global, ONLY: stdout
      USE sym_type, ONLY: nst, ist_i, st_name, ist_ityp
      USE neighbor, ONLY: nneighbors, nn_dist, inn_i
      USE constants, ONLY: eps6
      USE uspp_param, ONLY: upf
      USE ions_base, ONLY: ityp
      !
      IMPLICIT NONE
      !
      CHARACTER(len=256) :: routine_name
      !! name of this subroutine
      LOGICAL :: ltouch
      !! if true, make additional pass to ensure touching spheres
      !! only if lrmt = .false.
      LOGICAL, ALLOCATABLE :: lrmt_fixed(:)
      !! if true, mt_rmt(ist) is already constrained 
      INTEGER :: ierr
      !! error code
      INTEGER :: ist, iat, inn !, jat
      !! iterators
      REAL(DP) :: rtmp, rtmp2
      !! real temporary vars
      REAL(DP) :: rmt_d_iat, rmt_d_iat_nn
      !! temprary default MT radii
      !
      EXTERNAL :: errore
      !
      routine_name = "set_rmt"
      !
      ltouch = .TRUE.
      !
      ALLOCATE(mt_rmt(nst), STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, 'Error allocating mt_rmt', 1)
      ALLOCATE(lrmt_fixed(nst), STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, 'Error allocating lrmt_fixed', 1)
      !
      mt_rmt(:) = -1.0_dp
      lrmt_fixed(:) = .FALSE.
      !
      IF (.NOT. lrmt) THEN
        !
        DO iat = 1, natoms
          !
          DO inn = 1, nneighbors(iat)
            !
            ! this symmetry type
            rmt_d_iat = rmt_default(st_name(ist_i(iat)))
            ! its nearest neighbor symmetry type
            rmt_d_iat_nn = &
              rmt_default(st_name(ist_i(inn_i(iat, inn))))
            !
            !
            rtmp = nn_dist(iat)  * & ! distance to the neighbor
              rmt_d_iat / & ! this atom
              (rmt_d_iat + rmt_d_iat_nn) ! this atom and its nearest neighbor
            !
            IF ((mt_rmt(ist_i(iat)) < 0.0_dp) .OR. &
              (mt_rmt(ist_i(iat)) > rtmp)) THEN
              !
              mt_rmt(ist_i(iat)) = rtmp
              !
            END IF
            !
          END DO ! inn
          !
          !
          !
          IF (mt_rmt(ist_i(iat)) < MAXVAL(upf(ityp(iat))%rcut(:))) THEN
            WRITE(stdout, '(6x, "symmetry type #", I4)') ist_i(iat)
            WRITE(stdout, '(6x, "MT radius: ", &
              F10.8, " bohr = ", F10.8, " A")') &
              mt_rmt(ist_i(iat)), mt_rmt(ist_i(iat)) * bohrtoang
            CALL errore(routine_name, &
              "First MT radius guess is too small.", 1)
          ELSE IF (mt_rmt(ist_i(iat)) > nn_dist(iat)) THEN
            WRITE(stdout, '(6x, "symmetry type #", I4)') ist_i(iat)
            WRITE(stdout, '(6x, "MT radius: ", &
              F10.8, " bohr = ", F10.8, " A")') &
              mt_rmt(ist_i(iat)), mt_rmt(ist_i(iat)) * bohrtoang
            CALL errore(routine_name, &
              "First MT radius guess is too high.", 1)
          ELSE IF (mt_rmt(ist_i(iat)) < 0.0_dp) THEN
            WRITE(stdout, '(6x, "symmetry type #", I4)') ist_i(iat)
            WRITE(stdout, '(6x, "MT radius: ", &
              F10.8, " bohr = ", F10.8, " A")') &
              mt_rmt(ist_i(iat)), mt_rmt(ist_i(iat)) * bohrtoang
            CALL errore(routine_name, &
              "First MT radius guess not assigned.", 1)
          END IF
          !
          !
        END DO ! iat
        !
        !
        IF (ltouch) THEN
          !
          DO iat = 1, natoms
            !
            DO inn = 1, nneighbors(iat)
              !
              IF (ABS(mt_rmt(ist_i(iat)) + mt_rmt(ist_i(inn_i(iat, inn))) - &
                nn_dist(iat)) < eps6) THEN
                lrmt_fixed(ist_i(iat)) = .TRUE.
                lrmt_fixed(ist_i(inn_i(iat, inn))) = .TRUE.
              END IF
              !
            END DO ! inn
            !
            IF (.NOT. lrmt_fixed(ist_i(iat))) THEN
              !
              rtmp = -1.0_dp
              !
              DO inn = 1, nneighbors(iat)
                !
                ! IF (lrmt_fixed(ist_i(inn_i(iat, inn)))) THEN
                  !
                  rtmp2 = nn_dist(iat) - mt_rmt(ist_i(inn_i(iat, inn)))
                  !
                  IF (((rtmp > 0.0_dp) .AND. (rtmp2 < rtmp)) .OR. &
                    (rtmp < 0.0_dp)) THEN
                    rtmp = rtmp2
                  END IF
                  !
                ! END IF
                !
              END DO ! inn
              !
              IF ((rtmp > 0.0_dp) .AND. (rtmp > mt_rmt(ist_i(iat)))) THEN
                mt_rmt(ist_i(iat)) = rtmp
              END IF
              !
              lrmt_fixed(ist_i(iat)) = .TRUE.
              !
              !
            END IF
            !
          END DO ! iat
          !
        END IF ! ltouch
        !
        !
        !
      ELSE
        !
        ! user-specified MT radii
        !
        DO iat = 1, natoms
          IF (mt_rmt(ist_i(iat)) < 0._dp .AND. rmt(iat) > 0._dp) THEN
            WRITE(stdout, '(6x, "MT-radius of atom ", I4, &
              " is used for the symmetry type ", I4)') iat, ist_i(iat)
            mt_rmt(ist_i(iat)) = rmt(iat)
          ELSE IF (ABS(mt_rmt(ist_i(iat)) - rmt(iat)) > eps6) THEN
            WRITE(stdout, '(6x, "Check MT-radius of atom ", I4, &
              " of the symmetry type ", I4)') iat, ist_i(iat)
            CALL errore(routine_name, &
              "MT-radii inconsistent for the symmetry type.", 1)
          END IF
        END DO
        !
      END IF ! lrmt
      !
      !
      ! check
      !
      DO ist = 1, nst
        !
        IF (mt_rmt(ist) < 0._dp) THEN
          !
          WRITE(stdout, '(6x, "MT radius for the type ", I0, &
            " is not defined")') ist
          CALL errore(routine_name, "Error in mt_rmt array", 1)
          !
        ELSE IF (mt_rmt(ist) < MAXVAL(upf(ist_ityp(ist))%rcut(:))) THEN
            !
            WRITE(stdout, '(6x, "symmetry type #", I4)') ist
            WRITE(stdout, '(6x, "MT radius: ", &
              F10.8, " bohr = ", F10.8, " A")') &
              mt_rmt(ist), mt_rmt(ist) * bohrtoang
            CALL errore(routine_name, "MT radius is too small.", 1)
            !
        END IF
        !
      END DO
      !
      !
      DO iat = 1, natoms
        !
        DO inn = 1, nneighbors(iat)
          !
          IF ((mt_rmt(ist_i(iat)) + mt_rmt(ist_i(inn_i(iat, inn))) - &
            nn_dist(iat)) > eps6) THEN
            !
            WRITE(stdout, '(/5x, "Spheres ", I0, " and ", I0, " overlap:")') &
              iat, inn_i(iat, inn)
            WRITE(stdout, '(/5x, "Check: ", F0.16, " + ", F0.16, " > " F0.16)') &
              mt_rmt(ist_i(iat)), mt_rmt(ist_i(inn_i(iat, inn))), nn_dist(iat)
            CALL errore(routine_name, "Error for overlapping spheres", 1)
            !
          END IF
          !
        END DO ! inn
        !
      END DO ! iat
      !
      !
      ! printing
      !
      WRITE(stdout, '(/5x, "MT radii for symmetry types")')
      DO ist = 1, nst
        WRITE(stdout, '(5x)')
        WRITE(stdout, '(6x, "symmetry type #", I4, "  ", A2)') ist, st_name(ist)
        WRITE(stdout, '(6x, "MT radius: ", F10.8, " bohr = ", F10.8, " A")') &
          mt_rmt(ist), mt_rmt(ist) * bohrtoang
      END DO ! ist
      WRITE(stdout, '(/5x, /5x)')
      !
      WRITE(stdout, '(/5x, "MT radii for atoms")')
      DO iat = 1, natoms
        WRITE(stdout, '(6x, "rmt(", I0, ") =  ", F0.16)') &
          iat, mt_rmt(ist_i(iat))
      END DO ! iat
      WRITE(stdout, '(/5x, /5x)')
      !
      DEALLOCATE(lrmt_fixed, STAT = ierr)
      IF (ierr /= 0) &
        CALL errore(routine_name, 'Error deallocating lrmt_fixed', 1)
      !
    !---------------------------------------------------------------------------
    END SUBROUTINE set_rmt
    !---------------------------------------------------------------------------
    !
    !
    !
    !---------------------------------------------------------------------------
    SUBROUTINE set_grids()
    !---------------------------------------------------------------------------
    !!
    !! Sets log grid for n MT radii
    !!
    !---------------------------------------------------------------------------
    !
    ! D. Radevych
    !
      USE kinds, ONLY: DP
      USE radial_grids, ONLY: radial_grid_type
      USE radial_grids, ONLY: do_mesh, check_mesh, deallocate_radial_grid
      USE io_global, ONLY: stdout
      USE sym_type, ONLY: nst
      USE uspp_param, ONLY: upf
      USE ions_base, ONLY: ityp
      !
      IMPLICIT NONE
      !
      EXTERNAL :: errore
      !
      INTEGER :: nr
      !! number of points on radial grid
      REAL(DP) :: rmin
      !! minimum radius of the mesh
      REAL(DP) :: rmax
      !! maximum radius of the mesh
      REAL(DP) :: zmesh
      !! the nuclear charge used for mesh
      REAL(DP) :: xmin
      !! minimum x of the linear mesh
      REAL(DP) :: dx
      !! the deltax of the linear mesh
      INTEGER :: ir, ist, itp
      !! iterators
      REAL(DP) :: tmp
      !! temporary variable
      LOGICAL :: lprint_grid = .FALSE.
      !! print grid (for testing)
      CHARACTER(len=256) :: routine_name
      !! name of this subroutine
      INTEGER :: ierr
      !! error code
      TYPE(radial_grid_type) :: tmp_grid
      !! temporary grid needed to bypass deallocation in do_mesh subroutine
      !
      !
      routine_name = "set_grids"
      !
      ! rmax = 100.0_dp ! bohr
      ! zmesh = 118.0_dp ! max Z in periodic table
      dx = 0.01_dp
      nr = mt_nrf
      zmesh = 1._dp
      rmin = 1.5094148E-05
      xmin = LOG(zmesh * rmin)
      !
      ALLOCATE(mt_rf(mt_nrf, nst), STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, 'Error allocating mt_rf', 1)
      !
      ALLOCATE(mt_rfab(mt_nrf, nst), STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, 'Error allocating mt_rfab', 1)
      !
      ALLOCATE(mt_dx(nst), STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, 'Error allocating mt_dx', 1)
      !
      DO ist = 1, nst
        !
        rmax = mt_rmt(ist)
        !
        dx = LOG(rmax / rmin) / (nr - 1)
        !
        IF (dx <= mt_prec) THEN
          WRITE(stdout, '(5x, "WARNING: Setting default dx.")')
          dx = 0.01_dp
        END IF
        !
        IF (rmax <= mt_prec) THEN
          WRITE(stdout, '(5x, "WARNING: Setting default rmax.")')
          rmax = 1E+12_dp
          DO itp = 1, SIZE(upf)
            tmp = MAXVAL(upf(itp)%r(:))
            IF (tmp < rmax) rmax = tmp
          END DO ! ityp
        END IF
        !
        ! log grid: note that tmp_grid is deallocated in do_mesh
        ! DO NOT deallocate tmp_grid in this module to avoid
        ! double deallocation!
        !
        CALL do_mesh(rmax, zmesh, xmin, dx, 0, tmp_grid)
        mt_rf(:, ist) = tmp_grid%r(:)
        mt_rfab(:, ist) = tmp_grid%rab(:)
        mt_dx(ist) = tmp_grid%dx
        !
        ! DO NOT uncomment this!
        ! CALL deallocate_radial_grid(tmp_grid)
        !
        !
        WRITE(stdout, '(/4x, "RMTA grid ", I0, " info")') ist
        !
        WRITE(stdout, '(6x, A, I3, A, I10)') &
          "grid(", ist, ")%mesh:", tmp_grid%mesh
        !
        IF (tmp_grid%mesh /= nr) THEN
          CALL errore(routine_name, "incorrect size of generated grid", 1)
        END IF
        !
        WRITE(stdout, '(6x, A, I3, A, F16.7)') &
          "grid(", ist, ")%rmax:", tmp_grid%rmax
        WRITE(stdout, '(6x, A, I3, A, F16.7)') &
          "grid(", ist, ")%zmesh:", tmp_grid%zmesh
        WRITE(stdout, '(6x, A, I3, A, F16.7)') &
          "grid(", ist, ")%xmin:", tmp_grid%xmin
        WRITE(stdout, '(6x, A, I3, A, F16.7)') &
          "grid(", ist, ")%dx:", tmp_grid%dx
        WRITE(stdout, '(6x, A, I3, A, F16.7)') &
          "grid(", ist, ")%r(1):", tmp_grid%r(1)
        WRITE(stdout, '(6x, A, I3, A, F16.7)') &
          "grid(", ist, ")%r(2):", tmp_grid%r(2)
        WRITE(stdout, '(6x, A, I3, A, F16.7)') &
          "grid(", ist, ")%r(3):", tmp_grid%r(3)
        WRITE(stdout, '(6x, A, I3, A, F16.7)') &
          "grid(", ist, ")%r(mesh):", tmp_grid%r(nr)
        !
        CALL check_mesh(tmp_grid)
        !
        CALL deallocate_radial_grid(tmp_grid)
        !
        IF (lprint_grid) THEN
          WRITE(stdout, '(/4x, "Generated grid ", I0, " :")')
          DO ir = 1, mt_nrf
            WRITE(stdout, '(6x, ES16.7)') mt_rf(ir, ist)
          END DO
          !
          ! CALL errore(routine_name, "TEST: grids printed", 1)
          !
        END IF
        !
      END DO ! ist
      !
      ! check
      !
      WRITE(stdout, '(/4x, "VERIFY GRID")')
      !
      DO ist = 1, nst
        !
        !
        WRITE(stdout, '(/4x, "RMTA grid ", I0, " info")') ist
        !
        WRITE(stdout, '(6x, A, I3, A, F16.7)') &
          "grid(", ist, ")%r(1):", mt_rf(1, ist)
        WRITE(stdout, '(6x, A, I3, A, F16.7)') &
          "grid(", ist, ")%r(2):", mt_rf(2, ist)
        WRITE(stdout, '(6x, A, I3, A, F16.7)') &
          "grid(", ist, ")%r(3):", mt_rf(3, ist)
        WRITE(stdout, '(6x, A, I3, A, F16.7)') &
          "grid(", ist, ")%r(mesh):", mt_rf(nr, ist)
        !
        !
      END DO ! ist
      !
      !
      ! CALL errore(routine_name, "Test DONE", 1)
      !
    !---------------------------------------------------------------------------
    END SUBROUTINE set_grids
    !---------------------------------------------------------------------------
    !
    !
    !---------------------------------------------------------------------------
    SUBROUTINE set_orbitals()
    !---------------------------------------------------------------------------
    !!
    !! Sets orbitals for RMTA purposes.
    !!
    !---------------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      CHARACTER(len=256) :: routine_name
      !! name of this subroutine
      INTEGER :: ierr
      !! error code
      INTEGER :: iorb
      !! iterators
      !
      EXTERNAL :: errore
      !
      routine_name = "set_orbitals"
      !
      ! RMTA orbitals
      ! s, p, d, f
      !
      norbs = 4
      rmta_lmax = norbs - 1
      !
      ALLOCATE(orb_label(norbs), STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, &
        'Error allocating orb_label', 1)
      !
      DO iorb = 1, norbs
        SELECT CASE (iorb)
          CASE (1)
            orb_label(iorb) = "s"
          CASE (2)
            orb_label(iorb) = "p"
          CASE (3)
            orb_label(iorb) = "d"
          CASE (4)
            orb_label(iorb) = "f"
          CASE DEFAULT
            orb_label(iorb) = "g"
        END SELECT
      END DO ! iorb
      !
      !
    !---------------------------------------------------------------------------
    END SUBROUTINE set_orbitals
    !---------------------------------------------------------------------------
    !
    !
    !---------------------------------------------------------------------------
    SUBROUTINE set_tau_cart()
    !---------------------------------------------------------------------------
    !!
    !! Sets tau in Cartisian coordinates (in bohr)
    !!
    !---------------------------------------------------------------------------
      USE io_global, ONLY : stdout
      USE ions_base, ONLY: tau
      USE cell_base, ONLY: alat
      !
      IMPLICIT NONE
      !
      CHARACTER(len=256) :: routine_name
      !! name of this subroutine
      INTEGER :: ierr
      !! error code
      INTEGER :: iat
      !! iterators
      !
      EXTERNAL :: errore
      !
      routine_name = "set_tau_cart"
      !
      ALLOCATE(tau_cart(3, natoms))
      IF (ierr /= 0) CALL errore(routine_name, &
        'Error allocating tau_cart', 1)
      !
      WRITE(stdout, '(/5x, A)') "tau (cartesian, bohr)"
      DO iat = 1, natoms
        tau_cart(:, iat) = tau(:, iat) * alat
        WRITE(stdout, '(6x, F10.4, "  ", F10.4, "  ", F10.4)') &
          tau_cart(:, iat)
      END DO
      WRITE(stdout, '(/5x)')
      !
    !---------------------------------------------------------------------------
    END SUBROUTINE set_tau_cart
    !---------------------------------------------------------------------------
    !
    !
    !---------------------------------------------------------------------------
    SUBROUTINE set_fermi_energy()
    !---------------------------------------------------------------------------
    !!
    !! Sets the Fermi energy
    !!
    !---------------------------------------------------------------------------
      ! USE io_global, ONLY : stdout, ionode
      USE lsda_mod, ONLY: nspin
      USE ener, ONLY: ef, ef_up, ef_dw
      USE compare, ONLY: ref_ef
      USE klist, ONLY: two_fermi_energies
      !
      IMPLICIT NONE
      !
      CHARACTER(len=256) :: routine_name
      !! name of this subroutine
      INTEGER :: ierr
      !! error code
      !
      EXTERNAL :: errore
      !
      routine_name = "set_fermi_energy"
      !
      ALLOCATE(fermi_energy(nspin), STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, &
        "Error allocating fermi_energy", 1)
      !
      ! Fermi energy
      fermi_energy(:) = ef
      IF (.NOT. luse_ref_pot) THEN
        !
        IF (two_fermi_energies) THEN
          fermi_energy(1) = ef_up
          fermi_energy(2) = ef_dw
        END IF
        !
      ELSE
        fermi_energy(:) = ref_ef
      END IF
      !
    !---------------------------------------------------------------------------
    END SUBROUTINE set_fermi_energy
    !---------------------------------------------------------------------------
    !
    !
    !---------------------------------------------------------------------------
    SUBROUTINE set_chem_type()
    !---------------------------------------------------------------------------
    !!
    !! Sets chemical types
    !!
    !---------------------------------------------------------------------------
      USE io_global, ONLY : stdout
      USE ions_base, ONLY: ityp
      USE uspp_param, ONLY: upf
      !
      IMPLICIT NONE
      !
      CHARACTER(len=256) :: routine_name
      !! name of this subroutine
      INTEGER :: ierr
      !! error code
      INTEGER :: iat, ict
      !! iterators
      !
      EXTERNAL :: errore
      !
      routine_name = "set_nat_per_chem_tp"
      !
      n_chem_types = SIZE(upf)
      WRITE(stdout, '(/6x, "n_chem_types: ", I0)') &
        n_chem_types
      !
      ! get number of atoms per type
      ALLOCATE(natoms_per_chem_type(n_chem_types), STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, &
        'Error allocating natoms_per_chem_type', 1)
      natoms_per_chem_type(:) = 0
      !
      DO ict = 1, n_chem_types
        DO iat = 1, natoms
          IF (ict == ityp(iat)) &
            natoms_per_chem_type(ict) = natoms_per_chem_type(ict) + 1
        END DO ! iat
      END DO ! ict
      !
    !---------------------------------------------------------------------------
    END SUBROUTINE set_chem_type
    !---------------------------------------------------------------------------
    !
    !
    !---------------------------------------------------------------------------
    SUBROUTINE set_upf_vars()
    !---------------------------------------------------------------------------
    !!
    !! Sets pseudo vars on corresponding radial grids from upf files
    !!
    !---------------------------------------------------------------------------
      USE uspp_param, ONLY: upf
      USE ions_base, ONLY: ityp
      !
      IMPLICIT NONE
      !
      CHARACTER(len=256) :: routine_name
      !! name of this subroutine
      INTEGER :: ierr
      !! error code
      INTEGER :: ict, ir, ichi, ibeta
      !! iterators
      !
      EXTERNAL :: errore
      !
      routine_name = "set_upf_vars"
      !
      ! pseudopotentials (upf)
      !
      ! number of points on r mesh from upf for each atom
      ALLOCATE(mt_nr(n_chem_types), STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, "Error allocating mt_nr", 1)
      DO ict = 1, n_chem_types
        mt_nr(ict) = upf(ict)%mesh
      END DO ! ict
      mt_nr_max = MAXVAL(mt_nr(:))
      !
      ! r values on r mesh for each type
      ALLOCATE(mt_r(mt_nr_max, n_chem_types), STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, "Error allocating mt_r", 1)
      mt_r(:, :) = 0.0_dp
      DO ict = 1, n_chem_types
        DO ir = 1, mt_nr(ict)
          mt_r(ir, ict) = upf(ict)%r(ir)
        END DO ! ir
      END DO ! ict
      !
      !
      ! dr / di values on r mesh for each type
      ALLOCATE(mt_rab(mt_nr_max, n_chem_types), STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, "Error allocating mt_rab", 1)
      mt_rab(:, :) = 0.0_dp
      DO ict = 1, n_chem_types
        DO ir = 1, mt_nr(ict)
          mt_rab(ir, ict) = upf(ict)%rab(ir)
        END DO ! ir
      END DO ! ict
      !
      ! radial functions chi(r) and beta-projectors beta(r)
      !
      ALLOCATE(nchis(n_chem_types), STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, "Error allocating nchis", 1)
      !
      ALLOCATE(nbetas(n_chem_types), STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, "Error allocating nbetas", 1)
      !
      DO ict = 1, n_chem_types
        nchis(ict) = upf(ict)%nwfc
        nbetas(ict) = upf(ict)%nbeta
      END DO ! ict
      !
      ALLOCATE(chi_label(MAXVAL(nchis), n_chem_types), STAT = ierr)
      IF (ierr /= 0) &
        CALL errore(routine_name, "Error allocating chi_label", 1)
      chi_label(:, :) = "??"
      !
      ALLOCATE(beta_label(MAXVAL(nbetas), n_chem_types), &
        STAT = ierr)
      IF (ierr /= 0) &
        CALL errore(routine_name, "Error allocating beta_label", 1)
      beta_label(:, :) = "??"
      DO ict = 1, n_chem_types
        !
        DO ichi = 1, nchis(ict)
          chi_label(ichi, ict) = upf(ict)%els(ichi)
        END DO ! ichi
        !
        DO ibeta = 1, nbetas(ict)
          beta_label(ibeta, ict) = upf(ict)%els_beta(ibeta)
        END DO ! ibeta
        !
      END DO ! ict
      !
      IF (lwrite_dat) THEN
        !
        ! chir: radial functions chi
        !
        ALLOCATE(chir(mt_nr_max, MAXVAL(nchis(:)), &
          n_chem_types), STAT = ierr)
        IF (ierr /= 0) &
          CALL errore(routine_name, "Error allocating chir", 1)
        chir(:, :, :) = 0.0_dp
        !
        DO ict = 1, n_chem_types
          DO ichi = 1, nchis(ict)
            DO ir = 1, mt_nr(ict)
                chir(ir, ichi, ict) = upf(ict)%chi(ir, ichi)
            END DO ! ir
          END DO ! ichi
        END DO ! ict
        !
        !
        ! betar: beta-projectors
        !
        ALLOCATE(betar(mt_nr_max, MAXVAL(nbetas(:)), &
          n_chem_types), STAT = ierr)
        IF (ierr /= 0) &
          CALL errore(routine_name, "Error allocating betar", 1)
        betar(:, :, :) = 0.0_dp
        !
        DO ict = 1, n_chem_types
          DO ibeta = 1, nbetas(ict)
            DO ir = 1, mt_nr(ict)
                betar(ir, ibeta, ict) = upf(ict)%beta(ir, ibeta)
            END DO ! ir
          END DO ! ibeta
        END DO ! ict
        !
      END IF
      !
      ! vlocionr: local ionic psuedopotential from upf
      ! vlocaer: local ionic all-electron potential from upf
      ! (non-zero only for PAW)
      !
      ALLOCATE(vlocionr(mt_nr_max, n_chem_types), STAT = ierr)
      IF (ierr /= 0) &
        CALL errore(routine_name, "Error allocating vlocionr", 1)
      vlocionr(:, :) = 0.0_dp
      !
      ALLOCATE(vlocaer(mt_nr_max, n_chem_types), STAT = ierr)
      IF (ierr /= 0) &
        CALL errore(routine_name, "Error allocating vlocaer", 1)
      vlocaer(:, :) = 0.0_dp
      !
      DO ict = 1, n_chem_types
        !
        DO ir = 1, mt_nr(ict)
          !
          ! copy pseudopotential = V_{loc}^{ion}(r)
          !
          vlocionr(ir, ict) = upf(ict)%vloc(ir)
          !
          !
        END DO ! ir
        !
        IF (upf(ict)%tpawp) THEN
          DO ir = 1, mt_nr(ict)
            !
            ! copy AE potential = V_{loc}^{AE}(r)
            !
            vlocaer(ir, ict) = upf(ict)%paw%ae_vloc(ir)
            !
            !
          END DO ! ir
        END IF ! tpawp
        !
      END DO ! ict
      !
      ! CALL errore(routine_name, "Test DONE", 1)
      !
    !---------------------------------------------------------------------------
    END SUBROUTINE set_upf_vars
    !---------------------------------------------------------------------------
    !
    !
    !---------------------------------------------------------------------------
    SUBROUTINE set_scf_vars()
    !---------------------------------------------------------------------------
    !!
    !! Sets vars from previous (N)SCF run
    !!
    !---------------------------------------------------------------------------
      USE io_global, ONLY: stdout
      USE vlocal, ONLY: vloc
      USE fft_base, ONLY: dfftp
      USE scf, ONLY: vltot, vrs, v, kedtau ! v_of_0
      USE gvecs, ONLY: doublegrid
      USE gvect, ONLY: ngl, gl, g, ngm ! ecutrho, mill
      USE fft_rho, ONLY: rho_r2g ! rho_g2r
      USE cell_base, ONLY: tpiba ! omega, at, bg, alat
      !
      IMPLICIT NONE
      !
      CHARACTER(len=256) :: routine_name
      !! name of this subroutine
      INTEGER :: ierr
      !! error code
      INTEGER :: ict, igm
      !! iterators
      REAl(DP) :: g3d(3)
      !! current g3d vector
      !
      EXTERNAL :: errore, set_vrs
      !
      routine_name = "set_scf_vars"
      !
      !
      ! G
      !
      ALLOCATE(mt_g(4, ngm), STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, "Error allocating mt_g", 1)
      !
      DO igm = 1, ngm
        !
        g3d(:) = g(:, igm)
        g3d(:) = tpiba * g3d(:)
        mt_g(1 : 3, igm) = g3d(:)
        mt_g(4, igm) = SQRT(SUM(g3d(:) * g3d(:)))
        !
      END DO ! igm
      !
      ! vlociong: V_ps(G) = V_{loc}^{ion}(|G|)
      !
      ALLOCATE(vlociong(ngl, n_chem_types), STAT = ierr)
      IF (ierr /= 0) &
        CALL errore(routine_name, "Error allocating vlociong", 1)
      DO ict = 1, n_chem_types
        vlociong(:, ict) = vloc(:, ict)
      END DO ! ict
      !
      !
      ! vlocionr3d: total local ion potential in 3D real space
      !
      ALLOCATE(vlocionr3d(dfftp%nnr), STAT = ierr)
      IF (ierr /= 0) &
        CALL errore(routine_name, "Error allocating vlocionr3d", 1)
      vlocionr3d(:) = vltot(:)
      ! vlocionr3d(:) = vrs(:, 1)
      !
      !
      ! vlocscrr3d: total local screened potential in 3D real space
      !
      CALL set_vrs(vrs, vltot, v%of_r, kedtau, v%kin_r, &
        dfftp%nnr, nspins, doublegrid)
      !
      ALLOCATE(vlocscrr3d(dfftp%nnr, nspins), STAT = ierr)
      IF (ierr /= 0) &
        CALL errore(routine_name, "Error allocating vlocscrr3d", 1)
      vlocscrr3d(:, :) = vrs(:, :)
      !
      ! vlocscfr3d: total local SCF part of the potential [scf(H+xc)]
      ALLOCATE(vlocscfr3d(dfftp%nnr, nspins), STAT = ierr)
      IF (ierr /= 0) &
        CALL errore(routine_name, "Error allocating vlocscfr3d", 1)
      vlocscfr3d(:, :) = v%of_r
      !
      !
      ! vlocscrg3d:
      ! \tilde{V}(\bm{G}) = \int d\bm{r} V(\bm{r}) e^{i \bm{G} \bm{r}}
      ALLOCATE(vlocscrg3d(ngm, nspins), STAT = ierr)
      IF (ierr /= 0) &
        CALL errore(routine_name, "Error allocating vlocscrg3d", 1)
      vlocscrg3d(:, :) = 0.0_dp
      !
      WRITE(stdout, '(/7x, A)') "FFT vrs to vgs"
      CALL rho_r2g(dfftp, v%of_r(:, :), vlocscrg3d(:, :), v=vltot)
      WRITE(stdout, '(7x, A)') "DONE FFT vrs to vgs"
      ! test scaling V(G = 0)
      vlocscrg3d(1, :) = vlocscrg3d(1, :) ! * scale_vg0
      !
      !
      ! vlocscfg3d:
      ! \tilde{V}_{scf}(\bm{G}) = \int d\bm{r} V(\bm{r}) e^{i \bm{G} \bm{r}}
      ALLOCATE(vlocscfg3d(ngm, nspins), STAT = ierr)
      IF (ierr /= 0) &
        CALL errore(routine_name, "Error allocating vlocscfg3d", 1)
      CALL rho_r2g(dfftp, vlocscfr3d(:, :), vlocscfg3d(:, :))
      !
    !---------------------------------------------------------------------------
    END SUBROUTINE set_scf_vars
    !---------------------------------------------------------------------------
    !
    !
    !---------------------------------------------------------------------------
    SUBROUTINE set_empty_arrays()
    !---------------------------------------------------------------------------
    !!
    !! Initialize empty arrays.
    !!
    !---------------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      CHARACTER(len=256) :: routine_name
      !! name of this subroutine
      INTEGER :: ierr
      !! error code
      !
      EXTERNAL :: errore
      !
      routine_name = "set_empty_arrays"
      !
      !
      ! allocate other arrays
      !
      ALLOCATE(loglrf(mt_nrf, norbs, nspins, natoms), &
        STAT = ierr)
      IF (ierr /= 0) &
         CALL errore(routine_name, "Error allocating loglrf", 1)
      loglrf(:, :, :, :) = 0._dp
      !
      ALLOCATE(dloglderf(mt_nrf, norbs, nspins, natoms), &
        STAT = ierr)
      IF (ierr /= 0) &
         CALL errore(routine_name, "Error allocating dloglderf", 1)
      dloglderf(:, :, :, :) = 0._dp
      !
      ALLOCATE(mll1rf(mt_nrf, norbs, &
        nspins, natoms), STAT = ierr)
      IF (ierr /= 0) &
        CALL errore(routine_name, "Error allocating mll1rf", 1)
      mll1rf(:, :, :, :) = 0.0_dp
      !
      ALLOCATE(mll1rf_label(norbs, &
        nspins, natoms), STAT = ierr)
      IF (ierr /= 0) &
        CALL errore(routine_name, "Error allocating mll1rf_label", 1)
      mll1rf_label(:, :, :) = "?????"
      !
      ALLOCATE(dos_nlmrf(1 : irf_delta + 1, &
        (rmta_lmax + 1) * (rmta_lmax + 1), nspins, natoms), STAT = ierr)
      IF (ierr /= 0) &
        CALL errore(routine_name, 'Error allocating dos_nlmrf', 1)
      dos_nlmrf(:, :, :, :) = 0._dp
      !
      ALLOCATE(dos_nlrf(1 : irf_delta + 1, &
        norbs, nspins, natoms), STAT = ierr)
      IF (ierr /= 0) &
        CALL errore(routine_name, 'Error allocating dos_nlrf', 1)
      dos_nlrf(:, :, :, :) = 0._dp
      !
      ALLOCATE(dos_nrf(1 : irf_delta + 1, &
        nspins, natoms), STAT = ierr)
      IF (ierr /= 0) &
        CALL errore(routine_name, 'Error allocating dos_nrf', 1)
      dos_nrf(:, :, :) = 0._dp
      !
      ALLOCATE(dos_n(nspins), STAT = ierr)
      IF (ierr /= 0) &
        CALL errore(routine_name, 'Error allocating dos_n', 1)
      dos_n(:) = 0._dp
      !
      ALLOCATE(etall1rf(1 : irf_delta + 1, norbs, &
        nspins, natoms), STAT = ierr)
      IF (ierr /= 0) &
        CALL errore(routine_name, "Error allocating etall1rf", 1)
      etall1rf(:, :, :, :) = 0.0_dp
      !
      ! V_SL
      !
      ! on upf grid
      !
      ALLOCATE(vsemilocr(mt_nr_max, norbs, n_chem_types), &
        STAT = ierr)
      IF (ierr /= 0) &
        CALL errore(routine_name, 'Error allocating vsemilocr', 1)
      vsemilocr(:, :, :) = 0.0_dp
      !
      ! on RMTA grid
      !
      ALLOCATE(vsemilocrf(mt_nrf, norbs, natoms), STAT = ierr)
      IF (ierr /= 0) &
        CALL errore(routine_name, 'Error allocating vsemilocrf', 1)
      vsemilocrf(:, :, :) = 0.0_dp
      !
      ! Vlocscr00
      !
      ALLOCATE(vlocscr00rf(mt_nrf, nspins, natoms))
      IF (ierr /= 0) &
        CALL errore(routine_name, 'Error allocating vlocscr00rf', 1)
      vlocscr00rf(:, :, :) = 0.0_dp
      !
      ! Vlocscf00
      ! TODO: not useful
      !
      ALLOCATE(vlocscf00rf(mt_nrf, nspins, natoms))
      IF (ierr /= 0) &
        CALL errore(routine_name, 'Error allocating vlocscf00rf', 1)
      vlocscf00rf(:, :, :) = 0.0_dp
      !
      ! CALL errore(routine_name, "Test DONE", 1)
      !
    !---------------------------------------------------------------------------
    END SUBROUTINE set_empty_arrays
    !---------------------------------------------------------------------------
    !
    !
    !---------------------------------------------------------------------------
    SUBROUTINE rmta_delete_vars()
    !---------------------------------------------------------------------------
    !!
    !! Deallocates RMTA variables.
    !!
    !---------------------------------------------------------------------------
      USE radial_grids, ONLY: deallocate_radial_grid
      !
      IMPLICIT NONE
      !
      CHARACTER(len=256) :: routine_name
      !! name of this subroutine
      INTEGER :: ierr
      !! error code
      !
      EXTERNAL :: errore
      !
      routine_name = "deallocate_rmta_vars"
      !
      DEALLOCATE(mt_rmt, STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, &
        'Error deallocating mt_rmt', 1)
      !
      DEALLOCATE(orb_label, STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, &
        'Error deallocating orb_label', 1)
      !
      DEALLOCATE(tau_cart, STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, &
        'Error deallocating tau_cart', 1)
      !
      DEALLOCATE(natoms_per_chem_type, STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, &
        'Error deallocating natoms_per_chem_type', 1)
      !
      DEALLOCATE(mt_nr, STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, &
        'Error deallocating mt_nr', 1)
      !
      DEALLOCATE(mt_r, STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, &
        'Error deallocating mt_r', 1)
      !
      DEALLOCATE(vlocionr, STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, &
        'Error deallocating vlocionr', 1)
      !
      DEALLOCATE(vlocaer, STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, &
        'Error deallocating vlocaer', 1)
      !
      DEALLOCATE(vlocionr3d, STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, &
        'Error deallocating vlocionr3d', 1)
      !
      DEALLOCATE(mt_rf, STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, &
        'Error deallocating mt_rf', 1)
      !
      !
      DEALLOCATE(vlocscr00rf, STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, &
        'Error deallocating vlocscr00rf', 1)
      !
      DEALLOCATE(vlocscf00rf, STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, &
        'Error deallocating vlocscf00rf', 1)
      !
      !
      DEALLOCATE(vlociong, STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, &
        'Error deallocating vlociong', 1)
      !
      DEALLOCATE(vlocscrr3d, STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, &
        'Error deallocating vlocscrr3d', 1)
      !
      DEALLOCATE(vlocscfr3d, STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, &
        'Error deallocating vlocscfr3d', 1)
      !
      DEALLOCATE(vlocscrg3d, STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, &
        'Error deallocating vlocscrg3d', 1)
      !
      DEALLOCATE(mt_g, STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, &
        'Error deallocating mt_g', 1)
      !
      DEALLOCATE(vlocscfg3d, STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, &
        'Error deallocating vlocscfg3d', 1)
      !
      IF (lwrite_dat) THEN
        !
        DEALLOCATE(chir, STAT = ierr)
        IF (ierr /= 0) CALL errore(routine_name, &
          'Error deallocating chir', 1)
        !
        DEALLOCATE(betar, STAT = ierr)
        IF (ierr /= 0) CALL errore(routine_name, &
          'Error deallocating betar', 1)
        !
      END IF
      !
      DEALLOCATE(vsemilocr, STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, &
        'Error deallocating vsemilocr', 1)
      !
      DEALLOCATE(vsemilocrf, STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, &
        'Error deallocating vsemilocrf', 1)
      !
      ! u(r) and derivatives
      DEALLOCATE(urf, STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, &
        'Error deallocating urf', 1)
      !
      DEALLOCATE(duderf, STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, &
        'Error deallocating duderf', 1)
      !
      DEALLOCATE(d2udrderf, STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, &
        'Error deallocating d2udrderf', 1)
      !
      ! computed RMTA quantities
      DEALLOCATE(mll1rf, STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, &
        'Error deallocating mll1rf', 1)
      !
      DEALLOCATE(mll1rf_label, STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, &
          'Error deallocating mll1rf_label', 1)
      !
      ! V(r)
      DEALLOCATE(vfullrf, STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, &
        'Error deallocating vfullrf', 1)
      !
      IF (lwrite_dat) THEN
        ! r V(r)
        DEALLOCATE(rvfullrf, STAT = ierr)
        IF (ierr /= 0) CALL errore(routine_name, &
          'Error deallocating rvfullrf', 1)
      END IF
      !
      !
      ! DOS
      !
      DEALLOCATE(dos_nlmrf, STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, &
        'Error deallocating dos_nlmrf', 1)
      !
      DEALLOCATE(dos_nlrf, STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, &
        'Error deallocating dos_nlrf', 1)
      !
      DEALLOCATE(dos_nrf, STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, &
        'Error deallocating dos_nrf', 1)
      !
      DEALLOCATE(dos_n, STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, &
        'Error deallocating dos_n', 1)
      !
      DEALLOCATE(etall1rf, STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, &
        'Error deallocating etall1rf', 1)
      !
      DEALLOCATE(loglrf, STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, &
        'Error deallocating loglrf', 1)
      !
      DEALLOCATE(dloglderf, STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, &
        'Error deallocating dloglderf', 1)
      !
      DEALLOCATE(mt_rab, STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, &
        'Error deallocating mt_rab', 1)
      !
      !
      ! chi(r), beta(r), and V_sl(r)
      DEALLOCATE(nchis, STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, &
        'Error deallocating nchis', 1)
      !
      DEALLOCATE(chi_label, STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, &
        'Error deallocating chi_label', 1)
      !
      DEALLOCATE(nbetas, STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, &
        'Error deallocating nbetas', 1)
      !
      DEALLOCATE(beta_label, STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, &
        'Error deallocating beta_label', 1)
      !
      DEALLOCATE(mt_rfab, STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, &
        'Error deallocating mt_rfab', 1)
      !
      DEALLOCATE(mt_dx, STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, &
        'Error deallocating mt_dx', 1)
      !
      DEALLOCATE(fermi_energy, STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, &
        'Error deallocating fermi_energy', 1)
      !
    !---------------------------------------------------------------------------
    END SUBROUTINE rmta_delete_vars
    !---------------------------------------------------------------------------
    !
    !
    !---------------------------------------------------------------------------
    REAL(DP) FUNCTION rmt_default(element_label)
    !---------------------------------------------------------------------------
    !!
    !! Get default MT-radius for the type
    !!
    !---------------------------------------------------------------------------
    !
    !  D. Radevych
    !
    !  Default values are
    !  courtesy of M. Weinert and flair: FLAPW code.
    !  https://sites.uwm.edu/weinert/flair/
    !
      !
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=2), INTENT(in) :: element_label
      REAL(DP) :: rmt_d
      !
      CHARACTER(len=256) :: routine_name
       !! name of this subroutine
      EXTERNAL :: errore
      !
      routine_name = "rmt_default"
      !
      SELECT CASE(element_label)
        !
        CASE(' ')  ! vacancy
          rmt_d = 2.50_dp
        !
        CASE('H')  ! H 1
          rmt_d = 0.60_dp
        !
        CASE('He')  ! He 2
          rmt_d = 0.70_dp
        !
        CASE('Li')  ! Li 3
          rmt_d = 1.70_dp
        !
        CASE('Be')  ! Be 4
          rmt_d = 1.70_dp
        !
        CASE('B')  ! B 5
          rmt_d = 1.55_dp
        !
        CASE('C')  ! C 6
          rmt_d = 1.30_dp
        !
        CASE('N')  ! N 7
          rmt_d = 1.40_dp
        !
        CASE('O')  ! O 8
          rmt_d = 1.50_dp
        !
        CASE('F')  ! F 9
          rmt_d = 2.20_dp
        !
        CASE('Ne') ! Ne 10
          rmt_d = 2.00_dp
        !
        CASE('Na') ! Na 11
          rmt_d = 2.20_dp
        !
        CASE('Mg') ! Mg 12
          rmt_d = 2.00_dp
        !
        CASE('Al') ! Al 13
          rmt_d = 2.00_dp
        !
        CASE('Si') ! Si 14
          rmt_d = 2.00_dp
        !
        CASE('P') ! P 15
          rmt_d = 2.00_dp
        !
        CASE('S') ! S 16
          rmt_d = 2.00_dp
        !
        CASE('Cl') ! Cl 17
          rmt_d = 2.80_dp
        !
        CASE('Ar') ! Ar 18
          rmt_d = 2.80_dp
        !
        CASE('K') ! K 19
          rmt_d = 2.80_dp
        !
        CASE('Ca') ! Ca 20
          rmt_d = 3.40_dp
        !
        CASE('Sc') ! Sc 21
          rmt_d = 2.50_dp
        !
        CASE('Ti') ! Ti 22
          rmt_d = 2.50_dp
        !
        CASE('V') ! V 23
          rmt_d = 2.50_dp
        !
        CASE('Cr') ! Cr 24
          rmt_d = 2.50_dp
        !
        CASE('Mn') ! Mn 25
          rmt_d = 2.50_dp
        !
        CASE('Fe') ! Fe 26
          rmt_d = 2.50_dp
        !
        CASE('Co') ! Co 27
          rmt_d = 2.50_dp
        !
        CASE('Ni') ! Ni 28
          rmt_d = 2.50_dp
        !
        CASE('Cu') ! Cu 29
          rmt_d = 2.50_dp
        !
        CASE('Zn') ! Zn 30
          rmt_d = 2.50_dp
        !
        CASE('Ga') ! Ga 31
          rmt_d = 2.50_dp
        !
        CASE('Ge') ! Ge 32
          rmt_d = 2.50_dp
        !
        CASE('As') ! As 33
          rmt_d = 2.50_dp
        !
        CASE('Se') ! Se 34
          rmt_d = 2.50_dp
        !
        CASE('Br') ! Br 35
          rmt_d = 2.50_dp
        !
        CASE('Kr') ! Kr 36
          rmt_d = 2.50_dp
        !
        CASE('Rb') ! Rb 37
          rmt_d = 2.50_dp
        !
        CASE('Sr') ! Sr 38
          rmt_d = 2.50_dp
        !
        CASE('Y') ! Y 39
          rmt_d = 2.50_dp
        !
        CASE('Zr') ! Zr 40
          rmt_d = 2.50_dp
        !
        CASE('Nb') ! Nb 41
          rmt_d = 2.50_dp
        !
        CASE('Mo') ! Mo 42
          rmt_d = 2.50_dp
        !
        CASE('Tc') ! Tc 43
          rmt_d = 2.50_dp
        !
        CASE('Ru') ! Ru 44
          rmt_d = 2.50_dp
        !
        CASE('Rh') ! Rh 45
          rmt_d = 2.50_dp
        !
        CASE('Pd') ! Pd 46
          rmt_d = 2.50_dp
        !
        CASE('Ag') ! Ag 47
          rmt_d = 2.50_dp
        !
        CASE('Cd') ! Cd 48
          rmt_d = 2.50_dp
        !
        CASE('In') ! In 49
          rmt_d = 2.50_dp
        !
        CASE('Sn') ! Sn 50
          rmt_d = 2.50_dp
        !
        CASE('Sb') ! Sb 51
          rmt_d = 2.50_dp
        !
        CASE('Te') ! Te 52
          rmt_d = 2.50_dp
        !
        CASE('I') ! I 53
          rmt_d = 2.50_dp
        !
        CASE('Xe') ! Xe 54
          rmt_d = 2.50_dp
        !
        CASE('Cs') ! Cs 55
          rmt_d = 2.80_dp
        !
        CASE('Ba') ! Ba 56
          rmt_d = 2.80_dp
        !
        CASE('La') ! La 57
          rmt_d = 2.80_dp
        !
        CASE('Ce') ! Ce 58
          rmt_d = 2.50_dp
        !
        CASE('Pr') ! Pr 59
          rmt_d = 2.50_dp
        !
        CASE('Nd') ! Nd 60
          rmt_d = 2.50_dp
        !
        CASE('Pm') ! Pm 61
          rmt_d = 2.50_dp
        !
        CASE('Sm') ! Sm 62
          rmt_d = 2.50_dp
        !
        CASE('Eu') ! Eu 63
          rmt_d = 2.50_dp
        !
        CASE('Gd') ! Gd 64
          rmt_d = 2.50_dp
        !
        CASE('Tb') ! Tb 65
          rmt_d = 2.50_dp
        !
        CASE('Dy') ! Dy 66
          rmt_d = 2.50_dp
        !
        CASE('Ho') ! Ho 67
          rmt_d = 2.50_dp
        !
        CASE('Er') ! Er 68
          rmt_d = 2.50_dp
        !
        CASE('Tm') ! Tm 69
          rmt_d = 2.50_dp
        !
        CASE('Yb') ! Yb 70
          rmt_d = 2.50_dp
        !
        CASE('Lu') ! Lu 71
          rmt_d = 2.50_dp
        !
        CASE('Hf') ! Hf 72
          rmt_d = 2.50_dp
        !
        CASE('Ta') ! Ta 73
          rmt_d = 2.50_dp
        !
        CASE('W') ! W 74
          rmt_d = 2.50_dp
        !
        CASE('Re') ! Re 75
          rmt_d = 2.50_dp
        !
        CASE('Os') ! Os 76
          rmt_d = 2.50_dp
        !
        CASE('Ir') ! Ir 77
          rmt_d = 2.50_dp
        !
        CASE('Pt') ! Pt 78
          rmt_d = 2.50_dp
        !
        CASE('Au') ! Au 79
          rmt_d = 2.50_dp
        !
        CASE('Hg') ! Hg 80
          rmt_d = 2.50_dp
        !
        CASE('Tl') ! Tl 81
          rmt_d = 2.50_dp
        !
        CASE('Pb') ! Pb 82
          rmt_d = 2.50_dp
        !
        CASE('Bi') ! Bi 83
          rmt_d = 2.50_dp
        !
        CASE('Po') ! Po 84
          rmt_d = 2.50_dp
        !
        CASE('At') ! At 85
          rmt_d = 2.50_dp
        !
        CASE('Rn') ! Rn 86
          rmt_d = 2.50_dp
        !
        CASE('Fr') ! Fr 87
          rmt_d = 2.50_dp
        !
        CASE('Ra') ! Ra 88
          rmt_d = 2.50_dp
        !
        CASE('Ac') ! Ac 89
          rmt_d = 2.50_dp
        !
        CASE('Th') ! Th 90
          rmt_d = 2.50_dp
        !
        CASE('Pa') ! Pa 91
          rmt_d = 2.50_dp
        !
        CASE('U') ! U 92
          rmt_d = 2.50_dp
        !
        CASE('Np') ! Np 93
          rmt_d = 2.50_dp
        !
        CASE('Pu') ! Pu 94
          rmt_d = 2.50_dp
        !
        CASE('Am') ! Am 95
          rmt_d = 2.50_dp
        !
        CASE('Cm') ! Cm 96
          rmt_d = 2.50_dp
        !
        CASE('Bk') ! Bk 97
          rmt_d = 2.50_dp
        !
        CASE('Cf') ! Cf 98
          rmt_d = 2.50_dp
        !
        CASE('Es') ! Es 99
          rmt_d = 2.50_dp
        !
        CASE('Fm') ! Fm 100
          rmt_d = 2.50_dp
        !
        CASE('Md') ! Md 101
          rmt_d = 2.50_dp
        !
        CASE('No') ! No 102
          rmt_d = 2.50_dp
        !
        CASE('Lr') ! Lr 103
          rmt_d = 2.50_dp
        !
        CASE('Rf') ! Rf 104
          rmt_d = 2.50_dp
        !
        CASE DEFAULT
          CALL errore(TRIM(routine_name), "element label is not defined", 1)
        !
      END SELECT
      !
      rmt_default = rmt_d
      !
    !---------------------------------------------------------------------------
    END FUNCTION rmt_default
    !---------------------------------------------------------------------------
    !
    !
  !=============================================================================
  END MODULE mt_var
  !=============================================================================
