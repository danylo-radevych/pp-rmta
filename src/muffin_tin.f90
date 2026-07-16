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
  MODULE muffin_tin
  !=============================================================================
  !!
  !! Main module for the RMTA program
  !!
  !=============================================================================
  !
  !  Danylo Radevych
  !  updated: 2026/06/06
  !  started: 2024/07/20
  !
    USE kinds, ONLY : DP
    USE io_global, ONLY : ionode, ionode_id, stdout
    !
    IMPLICIT NONE
    !
    ! shorthand name of the module: mt
    !
    PUBLIC :: & ! subroutines
      rmta_init, rmta_compute, vg3d_to_vloc00r, rmta_quit
    !
    ! private
    !
    !!
    !
    ! public
    !
    !
    INTERFACE spline_interpolation
      MODULE PROCEDURE spline_interpolation_rarr
      MODULE PROCEDURE spline_interpolation_rsca
      MODULE PROCEDURE spline_interpolation_rsca_range
    END INTERFACE
  !
  CONTAINS
  !
    !
    !
    !
    !---------------------------------------------------------------------------
    SUBROUTINE rmta_init()
    !---------------------------------------------------------------------------
    !!
    !! Init everything needed for RMTA run.
    !!
    !---------------------------------------------------------------------------
    !
    !  D. Radevych
    !
      USE kinds,             ONLY : DP
      USE ions_base,         ONLY : ityp ! tau
      ! USE scf,               ONLY : rho, v_of_0, vltot, vrs, v, kedtau
      USE gvect,             ONLY : ngl, gl, g ! mill, ecutrho
      ! USE cell_base,         ONLY : at, bg, alat, tpiba, omega
      ! USE vlocal,            ONLY : vloc
      ! USE gvecs,             ONLY : doublegrid
      ! USE fft_rho,           ONLY : rho_g2r, rho_r2g
      ! USE constants,         ONLY : fpi, pi, tpi
      USE compare, ONLY: ref_pot
      ! USE symm_base, ONLY: nosym
      USE sym_type, ONLY: allocate_st, ist_i
      USE neighbor, ONLY: allocate_nn
      USE mt_var, ONLY: &
        rmta_set_vars, &
        natoms, nspins, &
        tau_cart, &
        luse_ref_pot, &
        vlocscr00rf, vlocscf00rf, &
        vlocscrg3d, mt_g, vlocscfg3d, &
        nchis, nbetas, &
        rmta_ng, mt_nrf, mt_rf
      USE mt_printing, ONLY: check_upf
      USE constants, ONLY: eps12
      !
      IMPLICIT NONE
      !
      INTEGER :: iat, ispin
      !! iterators
      CHARACTER(len=256) :: routine_name
      !! name of this subroutine
      !
      EXTERNAL :: sph_bes
      EXTERNAL :: lschps
      EXTERNAL :: errore, start_clock, stop_clock
      !
      !
      routine_name = "rmta_init"
      CALL start_clock(routine_name)
      !
      !
      ! symmetry types
      ! IF (.NOT. nosym) &
      CALL allocate_st()
      !
      ! nearest neighbors
      CALL allocate_nn()
      !
      ! RMTA vars
      CALL rmta_set_vars()
      !
      !
      !
      !
      !
      !
      ! check upf data
      CALL check_upf()
      !
      !
      ! fine radial grid >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !
      !
      ! vlocscr00rf: local screened atomic potential,
      ! on fine radial grid
      ! DR: VERY IMPORTANT
      !
      IF (.NOT. luse_ref_pot) THEN
        !
        ! V(r) = \sum_{\bm{G}} e^{-i \bm{G} \cdot \bm{tau}}
        ! V(\bm{G}) sin(G r) / (G r)
        !
        CALL vg3d_to_vloc00r(natoms, nspins, mt_nrf, &
          ist_i, mt_rf, &
          rmta_ng, mt_g, tau_cart, vlocscrg3d, vlocscr00rf)
        !
      ELSE
        DO iat = 1, natoms
          DO ispin = 1, nspins
            vlocscr00rf(:, ispin, iat) = ref_pot(:, 2)
          END DO ! ispin
        END DO ! iat
      END IF
      !
      !
      !
      ! vlocscf00rf: scf[H+XC] part of the local screened atomic potential,
      ! on fine radial grid
      ! TODO: not useful
      !
      ! V(r) = \sum_{\bm{G}} e^{- i \bm{G} \cdot \bm{tau}}
      ! V(\bm{G}) sin(G r) / (G r)
      !
      CALL vg3d_to_vloc00r(natoms, nspins, mt_nrf, &
        ist_i, mt_rf, &
        rmta_ng, mt_g, tau_cart, vlocscfg3d, vlocscf00rf)
      !
      !
      !
      !
      !
      CALL stop_clock(routine_name)
      !
    !---------------------------------------------------------------------------
    END SUBROUTINE rmta_init
    !---------------------------------------------------------------------------
    !
    !
    !---------------------------------------------------------------------------
    SUBROUTINE vg3d_to_vloc00r(nat, nspin, nr, stp, r, &
      ng, g, tau, vg3d, vloc00r)
    !---------------------------------------------------------------------------
    !!
    !! Extracts local atomic potential V(r, spin, atom)
    !! from 3D periodic potential V(\bm{G}, spin, atom)
    !! by using properties of the Fourier transform and keeping
    !! only spherically-symmetric (L = 00) part.
    !!
    !---------------------------------------------------------------------------
    !
    !  Danylo Radevych
    !
      USE constants, ONLY: eps12
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(in) :: nat
      !! number of atoms
      INTEGER, INTENT(in) :: nspin
      !! number of spins
      INTEGER, INTENT(in) :: nr
      !! number of points on radial grids
      INTEGER, INTENT(in) :: stp(:)
      !! array converting atom indices into symmetry type indices, stp(atom)
      REAL(DP), INTENT(in) :: r(:, :)
      !! radial grid for each symmetry type, r(:, sym_type)
      INTEGER, INTENT(in) :: ng
      !! number of G-vectors
      REAL(DP), INTENT(in) :: g(:, :)
      !! array of G-vectors,
      !! where g(1:3, :) are G_i components and g(4, :) is |G|
      !! g(4, ng)
      REAL(DP), INTENT(in) :: tau(:, :)
      !! atomic coordinates (cartesian, bohr), tau(3, nat)
      COMPLEX(DP), INTENT(in) :: vg3d(:, :)
      !! total 3D potential in reciprocal space, V(\bm{G}, spin)
      REAL(DP), INTENT(inout) :: vloc00r(:, :, :)
      !! local atomic potential, V(r, spin, atom)
      !
      REAL(DP) :: gtau
      !! \bm{G} \cdot \bm{tau}
      REAL(DP) :: sinc
      !! value of the sinc function
      REAL(DP) :: exp_factor
      !! exponential factor
      INTEGER :: iat, ispin, ir, ig
      !! iterators
      CHARACTER(len=256) :: routine_name
      !! name of this subroutine
      !
      routine_name = "vg3d_to_vloc00r"
      !
      ! OLD
      ! V(r) = \sum_{\bm{G}} e^{- i \bm{G} \cdot \bm{tau}}
      ! V(\bm{G}) sin(G r) / (G r)
      !
      ! agreed with literature
      ! V(r) = \sum_{\bm{G}} e^{+ i \bm{G} \cdot \bm{tau}}
      ! V(\bm{G}) sin(G r) / (G r)
      !
      DO iat = 1, nat
        DO ispin = 1, nspin
          DO ir = 1, nr
            DO ig = 1, ng
              !
              ! \bm{G} \cdot \bm{tau}
              !
              gtau = DOT_PRODUCT(g(1 : 3, ig), tau(:, iat))
              !
              ! exponential factor
              !
              ! old
              ! exp_factor = CMPLX(COS(gtau), -SIN(gtau), KIND=DP)
              !
              ! agreed with literature
              exp_factor = CMPLX(COS(gtau), +SIN(gtau), KIND=DP)
              !
              ! value of the sinc function
              !
              sinc = g(4, ig) * r(ir, stp(iat))
              IF (ABS(sinc) > eps12) THEN
                ! definition of sinc
                sinc = SIN(sinc) / sinc
              ELSE
                ! zero-argument limit is 1
                sinc = 1.0_dp
              END IF
              !
              vloc00r(ir, ispin, iat) = &
                vloc00r(ir, ispin, iat) + &
                REAL(exp_factor * vg3d(ig, ispin)) * sinc
              !
            END DO ! ig
          END DO ! ir
        END DO ! ispin
      END DO ! iat
      !
      ! CALL errore(routine_name, "Test DONE", 1)
      !
    !---------------------------------------------------------------------------
    END SUBROUTINE vg3d_to_vloc00r
    !---------------------------------------------------------------------------
    !
    !
    !---------------------------------------------------------------------------
    SUBROUTINE rmta_compute()
    !---------------------------------------------------------------------------
    !!
    !! Executive part of RMTA.
    !!
    !---------------------------------------------------------------------------
    !
    !  D. Radevych
    !
      USE mt_var, ONLY: norbs, ltetra
      USE part_dos, ONLY: set_dos_n
      USE mh_eta, ONLY: set_eta
      USE mt_var, ONLY: &
        mt_nr, mt_r, vsemilocr, &
        vsemilocrf, &
        irf_min, irf_max, &
        urf, dudrrf, duderf, d2udrderf, &
        loglrf, dloglderf, &
        dos_nlmrf, dos_nlrf, dos_nrf, dos_n, &
        dos_nlmrf_nodloglde, dos_nlrf_nodloglde, dos_nrf_nodloglde, &
        luse_tot_dos, &
        mll1rf, mll1rf_nodloglde, &
        etall1rf, etall1rf_nodloglde, &
        natoms, tau_cart, n_chem_types, &
        nspins, fermi_energy, &
        lsemilocupf, mt_rab, &
        lhybrid, mt_ngauss, mt_degauss, fermi_energy, &
        mt_nrf, mt_rf, mt_dx
      USE mt_printing, ONLY: print_at_rmt
      USE sym_type, ONLY: ist_i
      !
      IMPLICIT NONE
      !
      ! semilocal potential V^l_{SL}(r)
      !
      CALL set_vsemiloc(lsemilocupf, &
        norbs, n_chem_types, &
        mt_nr, mt_r, mt_rab, vsemilocr, &
        mt_nrf, ist_i, mt_rf, vsemilocrf)
      !
      ! radial functions u_l(r, E_F)
      !
      CALL solve_radial()
      !
      !
      ! log derivatives L_l(r, E_F), d L_l / d e
      !
      CALL set_log_ders(irf_max, mt_dx, mt_rf, &
        natoms, ist_i, norbs, nspins, lhybrid, &
        urf, dudrrf, duderf, d2udrderf, &
        loglrf, dloglderf)
      !
      ! Pettifor's M_{l, l+1}
      !
      CALL set_pet_mll1()
      !
      ! partial DOS
      !
      CALL set_dos_n(ltetra, mt_nrf, irf_min, irf_max, &
        ist_i, &
        natoms, norbs, nspins, mt_ngauss, mt_rf, &
        tau_cart(1 : 3, 1 : natoms), &
        dloglderf(1 : irf_max, &
          1 : norbs, 1 : nspins, 1 : natoms), &
        mt_degauss, fermi_energy, &
        dos_nlmrf, dos_nlrf, dos_nrf, dos_n, &
        dos_nlmrf_nodloglde, dos_nlrf_nodloglde, dos_nrf_nodloglde)
      !
      ! McMillan-Hopfield \eta_l and \eta = \sum_l \eta_l
      !
      CALL set_eta(mt_nrf, irf_min, irf_max, &
        natoms, norbs, nspins, &
        dos_nlrf, dos_nrf, &
        dos_nlrf_nodloglde, dos_nrf_nodloglde, &
        dos_n, luse_tot_dos, &
        mll1rf, mll1rf_nodloglde, &
        etall1rf, etall1rf_nodloglde)
      !
      ! Interpolate and print quantities at specified MT-radius
      !
      CALL print_at_rmt()
      !
      !
    !---------------------------------------------------------------------------
    END SUBROUTINE rmta_compute
    !---------------------------------------------------------------------------
    !
    !
    !---------------------------------------------------------------------------
    SUBROUTINE set_vsemiloc(lread_semiloc_from_upf, &
      norb, n_chem_tp, &
      nr, r, rab, vsemilocr, &
      nrf, stp, rf, &
      vsemilocrf)
    !---------------------------------------------------------------------------
    !!
    !! Sets l-depndent semilocal part of the potential V^l_{SL}(r)
    !!
    !---------------------------------------------------------------------------
    !
    ! D. Radevych
    !
      !
      USE kinds, ONLY: DP
      USE uspp_param, ONLY: upf
      USE ions_base, ONLY: nat, ityp
      USE mt_var, ONLY: luse_ref_pot, mt_nr, mt_r, nchis, nbetas
      USE constants, ONLY: eps12, eps32
      !
      IMPLICIT NONE
      !
      LOGICAL, INTENT(in) :: lread_semiloc_from_upf
      !! if true, read V_SL(r) from upf file;
      !! if false, compute V_SL from V_NL(r)
      INTEGER, INTENT(in) :: norb
      !! number of orbitals
      INTEGER, INTENT(in) :: n_chem_tp
      !! number of types
      INTEGER, INTENT(in) :: nr(:)
      !! number of points on pseudo r grid for each type
      REAL(DP), INTENT(in) :: r(:, :)
      !! pseudo r grid for each type
      REAL(DP), INTENT(in) :: rab(:, :)
      !! pseudo rab for each type
      REAL(DP), INTENT(out) :: vsemilocr(:, :, :)
      !! V_SL(r) on pseudo r grid
      INTEGER, INTENT(in) :: nrf
      !! number of points on fine r grid
      INTEGER, INTENT(in) :: stp(:)
      !! array converting atom indices into symmetry type indices, stp(atom)
      REAL(DP), INTENT(in) :: rf(:, :)
      !! radial grid for each symmetry type, rf(:, sym_type)
      ! REAL(DP), INTENT(in) :: rf(:)
      ! !! fine r grid
      REAL(DP), INTENT(out) :: vsemilocrf(:, :, :)
      !! V_SL(r) on fine r grid
      !
      INTEGER :: ir, ichi, ibeta, jbeta, iorb, ict, iat
      !! iterators
      INTEGER :: l
      !! l number
      INTEGER :: lmax
      !! max l number
      REAL(DP) :: betachi
      !! local integral \int d r' beta(r') chi(r')
      REAL(DP) :: rmin = 0.01_dp ! bohr
      ! REAL(DP) :: rmin = 0.05_dp ! bohr
      !! values below rmin are "extrapolated"
      REAL(DP) :: v0
      !! set V_SL(r < rmin) value
      LOGICAL, ALLOCATABLE :: lvlfound(:)
      !! indicate if corresponding l is found
      INTEGER, ALLOCATABLE :: lvlchiindx(:)
      !! index of chi corresponding to l
      CHARACTER(len=256) :: routine_name
      !! name of this subroutine
      INTEGER :: ierr
      !! error code
      !
      EXTERNAL :: errore, start_clock, stop_clock
      !
      routine_name = "set_vsemiloc"
      CALL start_clock(routine_name)
      !
      !
      IF (lread_semiloc_from_upf) THEN
        !
        ! read semilocal from upf
        !
        DO ict = 1, n_chem_tp
          !
          lmax = 0
          DO ichi = 1, nchis(ict)
            l = rmta_get_lchi(ichi, ict)
            IF (l > lmax) lmax = l
          END DO
          ! WRITE(*, *) "lmax == ", lmax
          !
          DO iorb = 1, lmax + 1
            !
            vsemilocr(1 : mt_nr(ict), iorb, ict) = &
              upf(ict)%vnl(1 : mt_nr(ict), iorb - 1, 1)
            !
            DO ir = 1, mt_nr(ict)
              IF (ABS(vsemilocr(ir, iorb, ict)) < eps12) &
                vsemilocr(ir, iorb, ict) = 0.0_dp
            END DO ! ir
            !
          END DO ! iorb
        END DO ! ict
        !
      ELSE
        !
        ! recover semilocal from non-local
        !
        DO ict = 1, n_chem_tp
          !
          lmax = 0
          DO ichi = 1, nchis(ict)
            l = rmta_get_lchi(ichi, ict)
            IF (l > lmax) lmax = l
          END DO
          ! WRITE(*, *) "lmax == ", lmax
          !
          ! array with flags whether V_SL for particular l is already found
          ALLOCATE(lvlfound(lmax + 1), STAT = ierr)
          IF (ierr /= 0) &
            CALL errore(routine_name, 'Error allocating lvlfound', 1)
          lvlfound(:) = .FALSE.
          !
          ! array with indices of chi corresponding to each l
          ALLOCATE(lvlchiindx(lmax + 1), STAT = ierr)
          IF (ierr /= 0) &
            CALL errore(routine_name, 'Error allocating lvlchiindx', 1)
          lvlchiindx(:) = 0
          !
          !
          DO ichi = 1, nchis(ict)
            !
            l = rmta_get_lchi(ichi, ict)
            !
            IF (l <= lmax) THEN
              IF (.NOT. lvlfound(l + 1)) THEN
                DO jbeta = 1, nbetas(ict)
                  IF (l == upf(ict)%lll(jbeta)) THEN
                    DO ibeta = 1, nbetas(ict)
                      !
                      betachi = rmta_integrate(nr(ict), &
                        rab(:, ict), &
                        upf(ict)%beta(:, jbeta) * &
                        upf(ict)%chi(:, ichi))
                      !
                      DO ir = 1, mt_nr(ict)
                        IF (ABS(upf(ict)%chi(ir, ichi)) > eps12) THEN
                          ! IF (upf(ict)%tpawp) THEN
                          !   vsemilocr(ir, l + 1, ict) = &
                          !     vsemilocr(ir, l + 1, ict) + &
                          !     (upf(ict)%dion(ibeta, jbeta) - &
                          !     rmta_pse(jbeta) * &
                          !     upf(ict)%qqq(ibeta, jbeta)) * &
                          !     upf(ict)%beta(ir, ibeta) * betachi / &
                          !     upf(ict)%chi(ir, ichi)
                          ! ELSE
                            vsemilocr(ir, l + 1, ict) = &
                              vsemilocr(ir, l + 1, ict) + &
                              upf(ict)%dion(ibeta, jbeta) * &
                              upf(ict)%beta(ir, ibeta) * betachi / &
                              upf(ict)%chi(ir, ichi)
                          ! END IF
                        END IF
                      END DO ! ir
                      !
                    END DO ! ibeta
                  END IF ! l
                END DO ! jbeta
                !
                !
                ! clean-up of the region r -> 0
                !
                DO ir = mt_nr(ict), 1, - 1
                  IF (mt_r(ir, ict) >= rmin) THEN
                    v0 = vsemilocr(ir, l + 1, ict)
                  ELSE
                    vsemilocr(ir, l + 1, ict) = v0
                  END IF
                END DO ! ir
                !
                lvlfound(l + 1) = .TRUE.
                lvlchiindx(l + 1) = ichi
                !
              END IF ! not lvlfound
              !
            END IF ! l >= lmax
            !
          END DO ! ichi
          !
          DEALLOCATE(lvlfound, STAT = ierr)
          IF (ierr /= 0) &
            CALL errore(routine_name, 'Error deallocating lvlfound', 1)
          !
          DEALLOCATE(lvlchiindx, STAT = ierr)
          IF (ierr /= 0) &
            CALL errore(routine_name, 'Error deallocating lvlchiindx', 1)
          !
        END DO ! ict
        !
      END IF ! read from upf
      !
      ! on the fine rf grid >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !
      WRITE(stdout, '(5x, "nat ", I4)') nat
      ! WRITE(*, *) "n_chem_tp", n_chem_tp
      WRITE(stdout, '(5x, "norb ", I4)') norb
      !
      IF (.NOT. luse_ref_pot) THEN
        DO iat = 1, nat
          DO iorb = 1, norb
            ! WRITE(*, *) "mt_nr(ityp(iat))", nr(ityp(iat))
            CALL spline_interpolation(nr(ityp(iat)), &
              r(1 : nr(ityp(iat)), ityp(iat)), &
              vsemilocr(1 : nr(ityp(iat)), iorb, ityp(iat)), &
              nrf, &
              rf(1 : nrf, stp(iat)), vsemilocrf(1 : nrf, iorb, iat))
            !
            DO ir =  1, nrf
              IF (ABS(vsemilocrf(ir, iorb, iat)) < eps32) THEN
                vsemilocrf(ir, iorb, iat) = 0.0_dp
              END IF
            END DO
            !
          END DO ! iorb
        END DO ! ict
      END IF
      !
      CALL stop_clock(routine_name)
      !
    !---------------------------------------------------------------------------
    END SUBROUTINE set_vsemiloc
    !---------------------------------------------------------------------------
    !
    !
    !---------------------------------------------------------------------------
    REAL(DP) FUNCTION rmta_integrate(nr, rab, fr)
    !---------------------------------------------------------------------------
    !!
    !! Integrates f(r) on logarithmic mesh
    !!
    !---------------------------------------------------------------------------
    !
    ! D. Radevych
    !
      USE kinds, ONLY: DP
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(in) :: nr
      !! number of elements in arrays rab and f(r)
      REAL(DP), INTENT(in) :: rab(:)
      !! dr / di of the log mesh
      REAL(DP), INTENT(in) :: fr(:)
      !! array f(r), function of the mesh
      !
      REAL(DP) :: integral
      !! stores final result
      INTEGER :: ir
      !! iterator
      LOGICAL :: lsimpson
      !! if true, use Simpson's integration
      !
      EXTERNAL :: simpson
      !
      integral = 0.0_dp
      lsimpson = .TRUE.
      !
      IF (lsimpson) THEN
        !
        ! Simpson method from QE
        !
        CALL simpson(nr, fr, rab, integral)
      ELSE
        !
        ! DR: simple step method
        !
        DO ir = 1, nr
          integral = integral + rab(ir) * fr(ir)
        END DO ! ir
      END IF ! lsimpson
      !
      !
      rmta_integrate = integral
      !
    !---------------------------------------------------------------------------
    END FUNCTION rmta_integrate
    !---------------------------------------------------------------------------
    !
    !
    !---------------------------------------------------------------------------
    PURE REAL(DP) FUNCTION rmta_integrate_u2(nr, dx, r, ur, l)
    !---------------------------------------------------------------------------
    !!
    !! Integrates u^2(r) on logarithmic mesh.
    !!
    !! Inspired by Hamann's code.
    !!
    !---------------------------------------------------------------------------
    !
    ! D. Radevych
    !
      USE kinds, ONLY: DP
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(in) :: nr
      !! max index of integration
      REAL(DP), INTENT(in) :: dx
      !! dx parameter for the radial grid
      REAL(DP), INTENT(in) :: r(:)
      !! radial grid
      REAL(DP), INTENT(in) :: ur(:)
      !! array u(r), function on the log mesh
      INTEGER, INTENT(in) :: l
      !! angular momentum quantum number
      !
      INTEGER :: ir
      !! iterator
      REAL(DP) :: al
      !! dx
      REAL(DP) :: amesh
      !! e^dx
      REAL(DP) :: ro
      !! r1 / sqrt(e^dx)
      REAL(DP) :: sn
      !! result
      !
      al = dx
      amesh = EXP(al)
      ro = r(1) / SQRT(amesh)
      sn = ro**(2 * l + 3) / (2 * l + 3)
      !
      DO ir = 1, nr - 3
        sn = sn + al * r(ir) * ur(ir)**2
      END DO
      !
      sn = sn + al * (23.0_dp * r(nr - 2) * ur(nr - 2)**2 + &
        28.0_dp * r(nr - 1) * ur(nr - 1)**2 + &
        9.0_dp * r(nr) * ur(nr)**2) / 24.0_dp
      !
      rmta_integrate_u2 = sn
      !
    !---------------------------------------------------------------------------
    END FUNCTION rmta_integrate_u2
    !---------------------------------------------------------------------------
    !
    !
    !---------------------------------------------------------------------------
    PURE INTEGER FUNCTION rmta_get_lchi(ichi, ict)
    !---------------------------------------------------------------------------
    !!
    !! Returns l of the corresponding chi(r)
    !!
    !---------------------------------------------------------------------------
    !
    ! D. Radevych
    !
      USE uspp_param, ONLY: upf
      USE mt_var, ONLY: nchis
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(in) :: ichi
      !! index of chi
      INTEGER, INTENT(in) :: ict
      !! index of type
      !
      INTEGER :: l
      !! l number
      !
      IF (ichi <= nchis(ict)) THEN
        l = upf(ict)%lchi(ichi)
      ELSE
        l = 3 ! 4F
      END IF
      !
      rmta_get_lchi = l
      !
    !---------------------------------------------------------------------------
    END FUNCTION rmta_get_lchi
    !---------------------------------------------------------------------------
    !
    !
    !
    !---------------------------------------------------------------------------
    SUBROUTINE set_log_ders(nin, dx, rf, nat, stp, &
      norb, nspins, lhybrid, urf, dudrrf, duderf, d2udrderf, &
      logl, dloglde)
    !---------------------------------------------------------------------------
    !!
    !! Sets log derivatives of radial functions
    !!
    !---------------------------------------------------------------------------
    !
    ! D. Radevych
    !
       USE kinds, ONLY: DP
       USE constants, ONLY: eps12
       !
       IMPLICIT NONE
       !
       INTEGER, INTENT(in) :: nin
       !! max index used for integration
       REAL(DP), INTENT(in) :: dx(:)
       !! dx of each radial RMTA grid
       REAL(DP), INTENT(in) :: rf(:, :)
       !! radial fine grid
       INTEGER, INTENT(in) :: nat
       !! number of atoms
       INTEGER, INTENT(in) :: stp(:)
       !! array converting atom indices into symmetry type indices, stp(atom)
       INTEGER, INTENT(in) :: norb
       !! number of orbitals
       INTEGER, INTENT(in) :: nspins
       !! number of spins
       LOGICAL, INTENT(in) :: lhybrid
       !! if true, use hybrid formula
       REAL(DP), INTENT(in) :: urf(:, :, :, :)
       !! u(r) on fine radial grid
       REAL(DP), INTENT(in) :: dudrrf(:, :, :, :)
       !! d u(r) / d r on fine radial grid
       REAL(DP), INTENT(in) :: duderf(:, :, :, :)
       !! d u(r, e) / d e on fine radial grid
       REAL(DP), INTENT(in) :: d2udrderf(:, :, :, :)
       !! d2 u(r, e) / d r / d e on fine radial grid
       !
       REAL(DP), INTENT(inout) :: logl(:, :, :, :)
       !! L_l(r, E_F)
       REAL(DP), INTENT(inout) :: dloglde(:, :, :, :)
       !! d L_l(r, E_F) / d e
       !
       REAL(DP) :: rmtf
       !! current MT radius, on fine grid
       INTEGER :: ir, iat, ispin, iorb
       !! iterators
       CHARACTER(len=256) :: routine_name
       !! name of this subroutine
       !
       EXTERNAL :: start_clock, stop_clock
       !
       routine_name = "set_log_ders"
       !
       CALL start_clock(routine_name)
       !
       DO iat = 1, nat
         DO iorb = 1, norb
           DO ispin = 1, nspins
             DO ir = 1, nin
               !
               IF (ABS(urf(ir, iorb, ispin, iat)) > eps12) THEN
                 !
                 rmtf = rf(ir, stp(iat))
                 !
                 logl(ir, iorb, ispin, iat) = &
                   rmtf * dudrrf(ir, iorb, ispin, iat) / &
                   urf(ir, iorb, ispin, iat) - 1.0_dp
                 !
                 IF (.NOT. lhybrid) THEN
                   dloglde(ir, iorb, ispin, iat) = rmtf / &
                     (urf(ir, iorb, ispin, iat) * &
                     urf(ir, iorb, ispin, iat)) * &
                     (d2udrderf(ir, iorb, ispin, iat) * &
                     urf(ir, iorb, ispin, iat) - &
                     dudrrf(ir, iorb, ispin, iat) * &
                     duderf(ir, iorb, ispin, iat))
                 ELSE
                   !
                   dloglde(ir, iorb, ispin, iat) =  - rmtf / &
                     (urf(ir, iorb, ispin, iat) * &
                     urf(ir, iorb, ispin, iat)) * &
                     rmta_integrate_u2(ir, dx(stp(iat)), rf(:, stp(iat)), &
                     urf(:, iorb, ispin, iat), iorb - 1)
                   !
                 END IF ! lhybrid
                 !
               END IF ! urf > eps12
               !
             END DO ! ir
           END DO ! ispin
         END DO ! iorb
       END DO ! iat
       !
       CALL stop_clock(routine_name)
       !
    !---------------------------------------------------------------------------
    END SUBROUTINE set_log_ders
    !---------------------------------------------------------------------------
    !
    !
    !---------------------------------------------------------------------------
    SUBROUTINE solve_radial()
    !---------------------------------------------------------------------------
    !!
    !! Solve radial Schrodinger equation for all available l.
    !!
    !---------------------------------------------------------------------------
    !
    ! D. Radevych
    !
      USE kinds, ONLY: DP
      USE uspp_param, ONLY: upf ! see upflib/pseudo_types.f90
      USE io_global, ONLY: stdout
      USE hamann, ONLY: hamann_lschps, hamann_lschvkbs
      USE ions_base, ONLY: ityp
      USE mt_var, ONLY: lwrite_dat, &
        natoms, orb_label, norbs, nspins, fermi_energy, &
        urf, duderf, dudrrf, d2udrderf, &
        vfullrf, rvfullrf, vlocscr00rf, vsemilocrf, &
        irf_max, &
        lsemiloc, lnonlocal, &
        mt_nrf, mt_rf, mt_dx
      USE sym_type, ONLY: ist_i
      USE constants, ONLY: eps12
      !
      IMPLICIT NONE
      !
      EXTERNAL :: errore
      !
      CHARACTER(len=200) :: routine_name
      !! name of this subroutine
      INTEGER :: nvkb
      !! number of VKB projectors
      REAL(DP), ALLOCATABLE :: vkb(:, :)
      REAL(DP), ALLOCATABLE :: evkb(:)
      INTEGER :: nin
      !! integration up to r(nin) (mode=3,4,5)
      INTEGER :: ierr
      !! error code
      REAL(DP) :: z
      !! atomic number
      REAL(DP) :: eps
      !! convergence factor
      INTEGER :: n, l
      !! main and angular quantum numbers
      REAL(DP) :: e
      !! energy at which wavefunction is calculated
      INTEGER :: nstop
      !! output of lschps used for error decision
      REAL(DP), ALLOCATABLE :: v0(:)
      !! potential
      INTEGER :: iat, ispin, ir, iorb
      !! iterators
      LOGICAL :: lint_to_rmt = .TRUE.
      !! integrate to rcut only
      LOGICAL :: lnormalize = .FALSE.
      !! normalize energy derivatives at each r
      REAL(DP), ALLOCATABLE :: urf_at_de_up(:)
      !! u(r, e + de) for energy derivative
      REAL(DP), ALLOCATABLE :: urf_at_de_down(:)
      !! u(r, e - de) for energy derivative
      REAL(DP), ALLOCATABLE :: dudrrf_at_de_up(:)
      !! du(r, e + de) / dr for energy derivative
      REAL(DP), ALLOCATABLE :: dudrrf_at_de_down(:)
      !! du(r, e - de) / dr for energy derivative
      REAL(DP), ALLOCATABLE :: norm_urf_at_e(:)
      !! norm of u(r, e)
      REAL(DP), ALLOCATABLE :: norm_urf_at_de_up(:)
      !! norm of u(r, e + de) for energy derivative
      REAL(DP), ALLOCATABLE :: norm_urf_at_de_down(:)
      !! norm of u(r, e - de) for energy derivative
      REAL(DP) :: de = 1.0E-6_dp
      !! de step for energy derivative
      !
      ! pre-defined variables (QE lschps)
      !
      ! INTEGER :: mode = 1
      !! find energy and wavefunction of bound states,
      !! scalar-relativistic (all-electron)
      ! INTEGER :: mode = 2
      !! mode number for finding energy and wavefunction of bound state
      ! INTEGER :: mode = 3
      !! mode number for fixed-energy calculation
      ! INTEGER :: mode = 5
      !! for pseudopotential to produce wavefunction beyond
      !! radius used for pseudopotential construction
      !
      routine_name = "solve_radial"
      !
      !
      ! prepare input for lschps
      !
      !
      nin = mt_nrf
      eps = eps12
      !
      !
      ALLOCATE(v0(mt_nrf), STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, "Error allocating v0", 1)
      v0(:) = 0.0_dp
      !
      ALLOCATE(vkb(mt_nrf, 3), STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, "Error allocating vkb", 1)
      vkb(:, :) = 0.0_dp
      ALLOCATE(evkb(3), STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, "Error allocating evkb", 1)
      evkb(:) = 0.0_dp
      !
      !
      ! allocate u(r, e) and its derivatives
      !
      ALLOCATE(urf(mt_nrf, norbs, &
        nspins, natoms), STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, "Error allocating urf", 1)
      urf(:, :, :, :) = 0.0_dp
      !
      ALLOCATE(duderf(mt_nrf, norbs, &
        nspins, natoms), STAT = ierr)
      IF (ierr /= 0) &
        CALL errore(routine_name, "Error allocating duderf", 1)
      duderf(:, :, :, :) = 0.0_dp
      !
      ALLOCATE(dudrrf(mt_nrf, norbs, &
          nspins, natoms), STAT = ierr)
      IF (ierr /= 0) &
        CALL errore(routine_name, "Error allocating dudrrf", 1)
      dudrrf(:, :, :, :) = 0.0_dp
      !
      ALLOCATE(d2udrderf(mt_nrf, norbs, &
        nspins, natoms), STAT = ierr)
      IF (ierr /= 0) &
        CALL errore(routine_name, "Error allocating d2udrderf", 1)
      d2udrderf(:, :, :, :) = 0.0_dp
      !
      !
      !
      ALLOCATE(vfullrf(mt_nrf, norbs, &
        nspins, natoms), STAT = ierr)
      IF (ierr /= 0) &
        CALL errore(routine_name, "Error allocating vfullrf", 1)
      vfullrf(:, :, :, :) = 0.0_dp
      !
      IF (lwrite_dat) THEN
        ALLOCATE(rvfullrf(mt_nrf, norbs, &
          nspins, natoms), STAT = ierr)
        IF (ierr /= 0) &
          CALL errore(routine_name, "Error allocating rvfullrf", 1)
        rvfullrf(:, :, :, :) = 0.0_dp
      END IF
      !
      IF (lint_to_rmt) THEN
        nin = irf_max
      END IF
      ! WRITE(*, *) "mt_rf(ist_i(iat))%r(nin)", mt_rf(ist_i(iat))%r(nin)
      !
      ! local arrays for energy derivative computation
      !
      ALLOCATE(urf_at_de_up(mt_nrf), STAT = ierr)
      IF (ierr /= 0) &
        CALL errore(routine_name, "Error allocating urf_at_de_up", 1)
      !
      ALLOCATE(urf_at_de_down(mt_nrf), STAT = ierr)
      IF (ierr /= 0) &
        CALL errore(routine_name, "Error allocating urf_at_de_down", 1)
      !
      ALLOCATE(dudrrf_at_de_up(mt_nrf), STAT = ierr)
      IF (ierr /= 0) &
         CALL errore(routine_name, "Error allocating dudrrf_at_de_up", 1)
      !
      ALLOCATE(dudrrf_at_de_down(mt_nrf), STAT = ierr)
      IF (ierr /= 0) &
        CALL errore(routine_name, "Error allocating dudrrf_at_de_down", 1)
      !
      !
      DO iat = 1, natoms
        !
        !
        z = upf(ityp(iat))%zp
        !
        !
        !
        ! DO ichi = 1, nchis(iat)
        DO iorb = 1, norbs
          !
          !
          WRITE(stdout, '(5x, "rf(nin) ", F8.6)') mt_rf(nin, ist_i(iat))
          !
          ! n = upf(ityp(iat))%nchi(ichi)
          n = 0
          !
          ! char to int
          ! READ(chi_label(ichi, ityp(iat))(1:1),'(I)') n
          ! WRITE(*, *) "n", n
          ! l = rmta_get_lchi(ichi, ityp(iat))
          l = iorb - 1
          WRITE(stdout, '(5x, "l ", I4)') l
          !
          DO ispin = 1, nspins
            !
            ! zero arrays before updating them
            urf_at_de_up(:) = 0.0_dp
            urf_at_de_down(:) = 0.0_dp
            !
            dudrrf_at_de_up(:) = 0.0_dp
            dudrrf_at_de_down(:) = 0.0_dp
            !
            !
            ! TODO 2025/04/18: try to use vlocscrrf - DR: bad idea
            !
            ! vfullrf(:, iorb, ispin, iat) = vlocscrrf(:, ispin, iat)
            !
            vfullrf(:, iorb, ispin, iat) = vlocscr00rf(:, ispin, iat)
            !
            !
            IF (lsemiloc) THEN
              vfullrf(:, iorb, ispin, iat) = &
                vfullrf(:, iorb, ispin, iat) + &
                vsemilocrf(:, iorb, iat)
            END IF
            !
            IF (lwrite_dat) THEN
              rvfullrf(:, iorb, ispin, iat) = vfullrf(:, iorb, ispin, iat) * &
                mt_rf(:, ist_i(iat))
            END IF
            !
            ! solve Schrodinger equation
            WRITE(stdout, '(/5x, A)') ""
            WRITE(stdout, '(5x, A)') "Solving radial Schrodinger for"
            WRITE(stdout, '(5x, "n l spin", I3, I3, I3, A4)') n, l, ispin, &
                orb_label(iorb)
            !
            !
            ! at e
            !
            e = fermi_energy(ispin)
            !
            IF (lnonlocal) THEN
              !
              ! CALL rmta_set_vkb(l, iat, mt_nrf, nvkb, vkb, evkb)
              ! WRITE(stdout, '(/5x, A, 3es13.4)') "evkb(:) == ", evkb(:)
              !
              CALL hamann_lschvkbs(l, nvkb, e, mt_rf(:, ist_i(iat)), &
                vfullrf(:, iorb, ispin, iat), vkb, &
                evkb, &
                urf(:, iorb, ispin, iat), &
                dudrrf(:, iorb, ispin, iat), &
                mt_nrf, nin)
            ELSE
              CALL hamann_lschps(l, nstop, e, mt_rf(:, ist_i(iat)), &
                vfullrf(:, iorb, ispin, iat), &
                urf(:, iorb, ispin, iat), &
                dudrrf(:, iorb, ispin, iat), mt_nrf, nin)
            END IF
            !
            ! at e + de
            !
            e = fermi_energy(ispin) + de
            IF (lnonlocal) THEN
              CALL hamann_lschvkbs(l, nvkb, e, mt_rf(:, ist_i(iat)), &
                vfullrf(:, iorb, ispin, iat), vkb, &
                evkb, &
                urf_at_de_up(:), &
                dudrrf_at_de_up(:), &
                mt_nrf, nin)
            ELSE 
              CALL hamann_lschps(l, nstop, e, mt_rf(:, ist_i(iat)), &
                vfullrf(:, iorb, ispin, iat), &
                urf_at_de_up(:), &
                dudrrf_at_de_up(:), mt_nrf, nin)
            END IF
            !
            ! dudrrf_at_de_up(:) = dudrrf_at_de_up(:) / mt_rabf(:)
            !
            ! at e - de
            !
            e = fermi_energy(ispin) - de
            !
            IF (lnonlocal) THEN
              CALL hamann_lschvkbs(l, nvkb, e, mt_rf(:, ist_i(iat)), &
                vfullrf(:, iorb, ispin, iat), vkb, &
                evkb, &
                urf_at_de_down(:), &
                dudrrf_at_de_down(:), &
                mt_nrf, nin)
            ELSE
              CALL hamann_lschps(l, nstop, e, mt_rf(:, ist_i(iat)), &
                vfullrf(:, iorb, ispin, iat), &
                urf_at_de_down(:), &
                dudrrf_at_de_down(:), mt_nrf, nin)
            END IF
            !
            !
            IF (lnormalize) THEN
              ALLOCATE(norm_urf_at_e(mt_nrf), STAT = ierr)
              IF (ierr /= 0) &
                CALL errore(routine_name, &
                  "Error allocating norm_urf_at_e", 1)
              norm_urf_at_e(:) = 1.0_dp
              !
              ALLOCATE(norm_urf_at_de_up(mt_nrf), STAT = ierr)
              IF (ierr /= 0) &
                CALL errore(routine_name, &
                  "Error allocating norm_urf_at_de_up", 1)
              norm_urf_at_de_up(:) = 1.0_dp
              !
              ALLOCATE(norm_urf_at_de_down(mt_nrf), STAT = ierr)
              IF (ierr /= 0) &
                CALL errore(routine_name, &
                  "Error allocating norm_urf_at_de_down", 1)
              norm_urf_at_de_down(:) = 1.0_dp
              !
              DO ir = 1, nin
                !
                norm_urf_at_e(:) = 1.0_dp / &
                  SQRT(rmta_integrate_u2(ir, mt_dx(ist_i(iat)), &
                  mt_rf(:, ist_i(iat)), &
                  urf(:, iorb, ispin, iat), iorb - 1))
                !
                norm_urf_at_de_up(ir) = 1.0_dp / &
                  SQRT(rmta_integrate_u2(ir, mt_dx(ist_i(iat)), &
                  mt_rf(:, ist_i(iat)), &
                  urf_at_de_up(:), iorb - 1))
                !
                norm_urf_at_de_down(ir) = 1.0_dp / &
                  SQRT(rmta_integrate_u2(ir, mt_dx(ist_i(iat)), &
                  mt_rf(:, ist_i(iat)), &
                  urf_at_de_down(:), iorb - 1))
                !
              END DO ! ir
              !
              urf(1 : nin, iorb, ispin, iat) = &
                urf(1 : nin, iorb, ispin, iat) * &
                norm_urf_at_e(1 : nin)
              dudrrf(1 : nin, iorb, ispin, iat) = &
                dudrrf(1 : nin, iorb, ispin, iat) * &
                norm_urf_at_e(1 : nin)
              !
              urf_at_de_up(1 : nin) = &
                urf_at_de_up(1 : nin) * &
                norm_urf_at_de_up(1 : nin)
              dudrrf_at_de_up(1 : nin) = &
                dudrrf_at_de_up(1 : nin) * &
                norm_urf_at_de_up(1 : nin)
              !
              urf_at_de_down(1 : nin) = &
                urf_at_de_down(1 : nin) * &
                norm_urf_at_de_down(1 : nin)
              dudrrf_at_de_down(1 : nin) = &
                dudrrf_at_de_up(1 : nin) * &
                norm_urf_at_de_down(1 : nin)
              !
              DEALLOCATE(norm_urf_at_e, STAT = ierr)
              IF (ierr /= 0) &
                CALL errore(routine_name, &
                  "Error deallocating norm_urf_at_e", 1)
              DEALLOCATE(norm_urf_at_de_up, STAT = ierr)
              IF (ierr /= 0) &
                CALL errore(routine_name, &
                  "Error deallocating norm_urf_at_de_up", 1)
              DEALLOCATE(norm_urf_at_de_down, STAT = ierr)
              IF (ierr /= 0) &
                CALL errore(routine_name, &
                  "Error deallocating norm_urf_at_de_down", 1)
              !
            END IF
            !
            !
            !
            ! du(e) / de = [u(e + de) - u(e - de)] / [2 de]
            duderf(1 : nin, iorb, ispin, iat) = &
              (urf_at_de_up(1 : nin) - urf_at_de_down(1 : nin)) / &
              (2.0_dp * de)
            !
            ! du'(r, e) / de = [u'(r, e + de) - u'(r, e - de)] / [2 de]
            d2udrderf(1 : nin, iorb, ispin, iat) = &
              (dudrrf_at_de_up(1 : nin) - dudrrf_at_de_down(1 : nin)) / &
              (2.0_dp * de)
            !
            IF (nstop == 1) THEN
              WRITE(stdout, '("WARNING: irregular termination in lschps", &
                I3)') n + l
              ! CALL errore(rmta_routine, &
              !   "irregular termination in lschps", n + l)
            END IF
            WRITE(stdout, '(A)') ""
          END DO ! ispin
          !
        END DO ! iorb
        !
      END DO ! iat
      !
      ! clean-up
      !
      DEALLOCATE(dudrrf_at_de_up, STAT = ierr)
      IF (ierr /= 0) &
        CALL errore(routine_name, "Error deallocating dudrrf_at_de_up", 1)
      DEALLOCATE(dudrrf_at_de_down, STAT = ierr)
      IF (ierr /= 0) &
        CALL errore(routine_name, "Error deallocating dudrrf_at_de_down", 1)
      !
      DEALLOCATE(urf_at_de_up, STAT = ierr)
      IF (ierr /= 0) &
        CALL errore(routine_name, "Error deallocating urf_at_de_up", 1)
      DEALLOCATE(urf_at_de_down, STAT = ierr)
      IF (ierr /= 0) &
        CALL errore(routine_name, "Error deallocating urf_at_de_down", 1)
      !
      DEALLOCATE(v0, STAT = ierr)
      IF (ierr /= 0) &
        CALL errore(routine_name, "Error deallocating v0", 1)
      !
      DEALLOCATE(vkb, STAT = ierr)
      IF (ierr /= 0) &
        CALL errore(routine_name, "Error deallocating vkb", 1)
      DEALLOCATE(evkb, STAT = ierr)
      IF (ierr /= 0) &
        CALL errore(routine_name, "Error deallocating evkb", 1)
      !
    !---------------------------------------------------------------------------
    END SUBROUTINE solve_radial
    !---------------------------------------------------------------------------
    !
    !
    !---------------------------------------------------------------------------
    SUBROUTINE set_pet_mll1()
    !---------------------------------------------------------------------------
    !!
    !! Computes matrix elements M_{l, l + 1}, as defined by Pettifor
    !! in 10.1088/0305-4608/7/6/017,
    !! on fine r grid
    !!
    !---------------------------------------------------------------------------
    !
    !  D. Radevych
    !
      USE kinds, ONLY: DP
      USE uspp_param, ONLY: upf
      USE io_global, ONLY : stdout
      USE constants, ONLY : rytoev
      USE const, ONLY : bohrtoang
      USE ions_base, ONLY: ityp
      USE mt_var, ONLY: &
        natoms, norbs, orb_label, &
        nspins, fermi_energy, &
        mll1rf_label, mll1rf, mll1rf_nodloglde, &
        irf_max, &
        mt_nrf, mt_rf, mt_rmt, &
        vlocscr00rf, vsemilocrf, &
        urf, loglrf, dloglderf, &
        lhybrid
      USE sym_type, ONLY: ist_i
      USE constants, ONLY: eps6, eps12
      !
      IMPLICIT NONE
      !
      INTEGER :: nin
      !! max index used for integration of u(r)
      ! INTEGER :: imt
      ! !! index of fine-grid point closest to r_{MT}
      INTEGER :: iat, iorb, ispin, ir
      !! iterators
      ! INTEGER :: l1, l2
      ! !! l values
      ! INTEGER :: n1, n2
      ! !! n values
      REAL(DP) :: v0
      !! V_{L=00}(r)
      REAL(DP) :: veff
      !! V(r) - E_F
      REAL(DP) :: rmtf
      !! current r_{MT} radius, on fine grid
      REAL(DP) :: ul
      !! u_l(r)
      REAL(DP) :: ul1
      !! u_l+1(r)
      REAL(DP) :: rl
      !! R_l(r) = u_l(r) / r
      REAL(DP) :: rl1
      !! R_l+1(r) = u_l+1(r) / r
      REAL(DP) :: logl
      !! L_l(r) = r R' / R = r u' / u - 1
      REAL(DP) :: logl1
      !! L_{l + 1}(r)
      REAL(DP) :: dloglde
      !! d L_l(r) / d e = r / u_l^2 [(d2 u / de dr) u - (d u / d r) (d u / d e)]
      REAL(DP) :: dlogl1de
      !! d L_{l + 1}(r) / d e
      REAL(DP) :: mll1_at_rmt
      !! Value of M_{l, l+1} at r_mt
      REAL(DP) :: mll1_nodloglde_at_rmt
      !! Value of m_{l, l+1} at r_mt
      !
      CHARACTER(len=256) :: routine_name
      !! name of this subroutine
      !
      EXTERNAL :: start_clock, stop_clock
      !
      routine_name = "set_pet_mll1"
      CALL start_clock(routine_name)
      !
      mll1rf(:, :, :, :) = 0.0_dp
      mll1rf_nodloglde(:, :, :, :) = 0.0_dp
      mll1rf_label(:, :, :) = "?????"
      !
      nin = irf_max
      !
      IF (lhybrid) THEN
        WRITE(stdout, '(/5x, /5x, /5x, &
          ">>>>>>>>>>>>  HYBRID BEGIN  <<<<<<<<<<<<")')
      ELSE
        WRITE(stdout, '(/5x, /5x, /5x, &
          ">>>>>>>>>>>> PETTIFOR BEGIN <<<<<<<<<<<<")')
      END IF
      !
      !
      DO iat = 1, natoms
        !
        WRITE(stdout, '(/7x, "Atom # ", I0, A4, &
          " ==============================================================")') &
          iat, upf(ityp(iat))%psd
        rmtf = mt_rf(mt_nrf, ist_i(iat))
        WRITE(stdout, '(/7x, "Desired r_mt: ", F10.4)') mt_rmt(ist_i(iat))
        WRITE(stdout, '(/7x, "Actual r_mt: ", F10.4)') rmtf
        !
        DO ispin = 1, nspins
          !
          WRITE(stdout, '(/7x, A, I3, A, I3)') &
              "atom: ", iat, " spin: ", ispin
          !
          ! interpolate veff
          !
          v0 = vlocscr00rf(mt_nrf, ispin, iat)
          !
          ! TODO 2025/04/18: try to use rmta_vlocscrrf - DR: bad idea
          !
          veff = v0 - fermi_energy(ispin)
          !
          WRITE(stdout, '(/7x, A, F16.8, A)') &
            "E_F", fermi_energy(ispin), " (Ry)"
          WRITE(stdout, '(7x, A, F16.8, A)') &
            "V(r_mt)", v0, " (Ry)"
          WRITE(stdout, '(7x, A, F16.8, A)') &
           "V(r_mt) - E_F", veff, " (Ry)"
          !
          !
          DO iorb = 1, norbs - 1
            !
            mll1rf_label(iorb, ispin, iat) = &
              orb_label(iorb) // &
              orb_label(iorb + 1)
            !
            DO ir = 1, nin
              !
              v0 = vlocscr00rf(ir, ispin, iat)
              !
              ! TODO 2025/04/18: try to use rmta_vlocscrrf - DR: bad idea
              !
              ! v0 = rmta_vlocscrrf(ir, ispin, iat)
              !
              !
              ! need to add 1 / 2 [V_l + V_l+1],
              ! although it is zero at r_mt
              !
              v0 = v0 + &
                0.5_dp * (vsemilocrf(ir, iorb, iat) + &
                vsemilocrf(ir, iorb + 1, iat))
              !
              veff = v0 - fermi_energy(ispin)
              rmtf = mt_rf(ir, ist_i(iat))
              !
              logl = loglrf(ir, iorb, ispin, iat)
              dloglde = dloglderf(ir, iorb, ispin, iat)
              !
              !
              ! WRITE(*, *) "logl == ", logl
              ! WRITE(*, *) "dloglde == ", dloglde
              !
              logl1 = loglrf(ir, iorb + 1, ispin, iat)
              dlogl1de = dloglderf(ir, iorb + 1, ispin, iat)
              !
              !
              ! WRITE(*, *) "logl1 == ", logl1
              ! WRITE(*, *) "dlogl1de == ", dlogl1de
              !
              !
              IF (ABS(rmtf) > eps12) THEN
                !
                mll1rf_nodloglde(ir, iorb, ispin, iat) = &
                  veff * rmtf * rmtf - &
                  (logl - (iorb - 1)) * (logl1 + (iorb - 1) + 2)
                mll1rf_nodloglde(ir, iorb, ispin, iat) = &
                  mll1rf_nodloglde(ir, iorb, ispin, iat) / rmtf
                !
                IF (ABS(dloglde) > eps12 .AND. &
                  ABS(dlogl1de) > eps12 .AND. &
                  (dloglde * dlogl1de > 0.0_dp)) THEN
                  !
                  mll1rf(ir, iorb, ispin, iat) = &
                    mll1rf_nodloglde(ir, iorb, ispin, iat)
                  mll1rf(ir, iorb, ispin, iat) = &
                    mll1rf(ir, iorb, ispin, iat) / &
                    SQRT(ABS(dloglde * dlogl1de))
                  !
                ELSE
                  !
                  IF (ir == nin) THEN
                    IF (ABS(dloglde) <= eps12) THEN
                      CALL errore(routine_name, &
                        "dloglde is close to zero at rmt", 1)
                    ELSE IF (ABS(dlogl1de) <= eps12) THEN
                      CALL errore(routine_name, &
                        "dlogl1de is close to zero at rmt", 1)
                    ELSE IF (dloglde * dlogl1de <= 0.0_dp) THEN
                      WRITE(stdout, '(/5x, &
                        "WARNING: dloglde * dlogl1de <= 0 at rmt")')
                      WRITE(stdout, '(6x, "dloglde = ", F0.16)') dloglde
                      WRITE(stdout, '(6x, "dlogl1de = ", F0.16)') dlogl1de
                      CALL errore(routine_name, &
                        "dloglde * dlogl1de <= 0", 1)
                    ELSE
                      CALL errore(routine_name, &
                        "problem with dloglde or dlogl1de", 1)
                    END IF
                  END IF
                  !
                END IF
                !
              ELSE
                !
                IF (ir == nin) &
                  CALL errore(routine_name, "rmtf is close to zero", 1)
                !
              END IF
              !
              !
            END DO ! ir
            !
            !
            ! at rmtf
            !
            rmtf = mt_rf(mt_nrf, ist_i(iat))
            !
            ul = urf(mt_nrf, iorb, ispin, iat)
            ul1 = urf(mt_nrf, iorb + 1, ispin, iat)
            rl = ul / rmtf
            rl1 = ul1 / rmtf
            logl = loglrf(mt_nrf, iorb, ispin, iat)
            logl1 = loglrf(mt_nrf, iorb + 1, ispin, iat)
            dloglde = dloglderf(mt_nrf, iorb, ispin, iat)
            dlogl1de = dloglderf(mt_nrf, iorb + 1, ispin, iat)
            !
            WRITE(stdout, '(/7x, A, I1, A, F16.8)') &
              "u_", iorb - 1, "(r_mt) = ", ul
            WRITE(stdout, '(7x, A, I1, A, F16.8)') &
              "R_", iorb - 1, "(r_mt) = ", rl
            WRITE(stdout, '(7x, A, I1, A, F16.8)') &
              "L_", iorb - 1, "(r_mt) = ", logl
            WRITE(stdout, '(7x, A, I1, A, F16.8, A)') &
              "dL_", iorb - 1, "(r_mt) / de = ", dloglde, " (1 / Ry)"
            !
            IF (ABS(ul) < eps6) THEN
              WRITE(stdout, '(5x, "rmt is close to a node of u_", I0, "(r)")') & 
                iorb - 1
              CALL errore(routine_name, "Adjust rmt.")
            END IF
            !
            WRITE(stdout, '(/7x, A, I1, A, F16.8)') &
              "u_", iorb, "(r_mt) = ", ul1
            WRITE(stdout, '(7x, A, I1, A, F16.8)') &
              "R_", iorb, "(r_mt) = ", rl1
            WRITE(stdout, '(7x, A, I1, A, F16.8)') &
              "L_", iorb, "(r_mt) = ", logl1
            WRITE(stdout, '(7x, A, I1, A, F16.8, A)') &
              "dL_", iorb, "(r_mt) / de = ", dlogl1de, " (1 / Ry)"
            !
            IF (ABS(ul1) < eps6) THEN
              WRITE(stdout, '(5x, "rmt is close to a node of u_", I0, "(r)")') & 
                iorb
              CALL errore(routine_name, "Adjust rmt.")
            END IF
            !
            !
            ! interpolate M_{l, l+1} at r_mt (Ry / bohr)
            !
            ! CALL spline_interpolation(nin, mt_rf(:), &
            !   mll1rf(:, iorb, ispin, iat), &
            !   rmtf, mll1_at_rmt)
            !
            mll1_nodloglde_at_rmt = mll1rf_nodloglde(mt_nrf, iorb, ispin, iat)
            mll1_at_rmt = mll1rf(mt_nrf, iorb, ispin, iat)
            !
            WRITE(stdout, '(/7x, A10, A, A, F16.8, A, F16.8, A)') &
              " m_", mll1rf_label(iorb, ispin, iat), &
              "(",  rmtf, "):", &
              mll1_nodloglde_at_rmt, &
              " "
            !
            WRITE(stdout, '(7x, A10, A, A, F16.8, A, F16.8, A)') &
              " M_", mll1rf_label(iorb, ispin, iat), &
              "(",  rmtf, "):", &
              mll1_at_rmt, &
              " (Ry / bohr)"
            WRITE(stdout, '(7x, A10, A, A, F16.8, A, F16.8, A)') &
              " M^2_", mll1rf_label(iorb, ispin, iat), &
              "(",  rmtf, "):", &
              mll1_at_rmt * &
              mll1_at_rmt, &
              " (Ry / bohr)^2"
            !
            WRITE(stdout, '(7x, A10, A, A, F16.8, A, F16.8, A)') &
              " M_", mll1rf_label(iorb, ispin, iat), &
              "(",  rmtf, "):", &
              mll1_at_rmt * &
              (rytoev / bohrtoang), &
              " (eV / A)"
            WRITE(stdout, '(7x, A10, A, A, F16.8, A, F16.8, A)') &
              " M^2_", mll1rf_label(iorb, ispin, iat), &
              "(",  rmtf, "):", &
              mll1_at_rmt * &
              mll1_at_rmt * &
              (rytoev / bohrtoang)**2, &
              " (eV / A)^2"
            !
            ! WRITE(stdout, '(5x, "Scaled:")')
            ! WRITE(stdout, '(7x, A10, A, A, F16.8, A, F16.8, A)') &
            !   " M_", mll1rf_label(iorb, ispin, iat), &
            !   "(",  rmtf, "):", &
            !   mll1_at_rmt * &
            !   (rytoev / bohrtoang) / SQRT(rmta_wds4), &
            !   " (Wd^1/2 / S^2)"
            ! WRITE(stdout, '(7x, A10, A, A, F16.8, A, F16.8, A)') &
            !   " M^2_", mll1rf_label(iorb, ispin, iat), &
            !   "(",  rmtf, "):", &
            !   mll1_at_rmt * &
            !   mll1_at_rmt * &
            !   (rytoev / bohrtoang)**2 / rmta_wds4, &
            !   " (Wd / S^4)"
            !
            WRITE(stdout, '("")')
            !
          END DO ! iorb
        END DO ! ispin
      END DO ! iat
      !
      IF (lhybrid) THEN
        WRITE(stdout, &
        '(/5x, ">>>>>>>>>>>>   HYBRID END   <<<<<<<<<<<<", /5x, /5x, /5x)')
      ELSE
        WRITE(stdout, &
        '(/5x, ">>>>>>>>>>>>  PETTIFOR END  <<<<<<<<<<<<", /5x, /5x, /5x)')
      END IF
      !
      CALL stop_clock(routine_name)
      !
    !---------------------------------------------------------------------------
    END SUBROUTINE set_pet_mll1
    !---------------------------------------------------------------------------
    !
    !
    !---------------------------------------------------------------------------
    SUBROUTINE spline_interpolation_rarr(xdim, x, y, xfdim, xf, yf)
    !---------------------------------------------------------------------------
    !!
    !! Spline interpolation to real array.
    !!
    !---------------------------------------------------------------------------
      USE kinds,             ONLY: DP
      USE splinelib,         ONLY: spline, splint
      USE constants, ONLY: eps12
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(in) :: xdim
      !! coarse grid dimension
      REAL(DP), INTENT(in) :: x(:)
      !! coarse grid
      REAL(DP), INTENT(in) :: y(:)
      !! function values on coarse grid
      INTEGER, INTENT(in) :: xfdim
      !! fine grid dimension
      REAL(DP), INTENT(in) :: xf(:)
      !! fine grid
      !
      REAL(DP), INTENT(inout) :: yf(:)
      !! function values on fine grid
      !
      REAL(DP), ALLOCATABLE :: d2y(:)
      REAL(DP) :: startu, startd
      REAL(DP) :: delta
      INTEGER :: ix
      REAL(DP) :: xmin
      !! min value of x on coarse grid
      REAL(DP) :: xmax
      !! max value of x on coarse grid
      CHARACTER(len=256) :: routine_name
      !! name of this subroutine
      INTEGER :: ierr
      !! error code
      !
      EXTERNAL :: errore
      !
      routine_name = "spline_interpolation_rarr"
      !
      ! xmin
      xmin = MINVAL(x(:))
      ! xmax
      xmax = MAXVAL(x(:))
      !
      ! first derivative
      delta = x(2) - x(1)
      IF ((xdim >= 2) .AND. (ABS(delta) > eps12)) THEN
        startu = (y(2) - y(1)) / delta
      ELSE
        startu = 0.0_dp
      END IF
      !
      ! second derivative
      delta = x(3) - x(1)
      IF ((xdim >= 3) .AND. (ABS(delta) > eps12)) THEN
        startd = (y(3) - 2.0_dp * y(2) + y(1)) / &
          delta
      ELSE
        startd = 0.0_dp
      END IF
      !
      ALLOCATE(d2y(xdim), STAT = ierr)
      IF (ierr /= 0) &
        CALL errore(routine_name, "Error allocating d2y", 1)
      !
      ! prepare spline
      CALL spline(x, y, startu, startd, d2y)
      !
      DO ix = 1, xfdim
        !
        IF ((xf(ix) <= xmax) .AND. (xf(ix) >= xmin)) THEN
          yf(ix) = splint(x, y, d2y, xf(ix))
        ELSE IF (xf(ix) >= xmax) THEN
          yf(ix) = y(xdim)
        ELSE
          IF (yf(1) /= yf(1)) THEN
            CALL errore(routine_name, 'yf(1) value is NAN', 1)
          END IF
          !
          yf(ix) = y(1)
        END IF
        !
        IF (yf(ix) /= yf(ix)) THEN
          CALL errore(routine_name, 'yf value is NAN', 1)
        END IF
        !
        IF (ABS(yf(ix)) < eps12) THEN
          yf(ix) = 0.0_dp
        END IF
        !
      END DO ! ix
      !
      DEALLOCATE(d2y, STAT = ierr)
      IF (ierr /= 0) &
        CALL errore(routine_name, 'Error deallocating d2y', 1)
      !
    !---------------------------------------------------------------------------
    END SUBROUTINE spline_interpolation_rarr
    !---------------------------------------------------------------------------
    !
    !
    !---------------------------------------------------------------------------
    SUBROUTINE spline_interpolation_rsca(xdim, x, y, xf, yf)
    !---------------------------------------------------------------------------
    !!
    !! Spline interpolation to real scalar.
    !!
    !---------------------------------------------------------------------------
      USE kinds, ONLY: DP
      USE splinelib, ONLY: spline, splint
      USE constants, ONLY: eps12
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(in) :: xdim
      !! coarse grid dimension
      REAL(DP), INTENT(in) :: x(:)
      !! coarse grid
      REAL(DP), INTENT(in) :: y(:)
      !! function values on coarse grid
      REAL(DP), INTENT(in) :: xf
      !! x on fine grid
      !
      REAL(DP), INTENT(out) :: yf
      !! function value on fine grid
      !
      REAL(DP), ALLOCATABLE :: d2y(:)
      REAL(DP) :: startu, startd
      REAL(DP) :: delta
      REAL(DP) :: xmin
      !! min value of x on coarse grid
      REAL(DP) :: xmax
      !! max value of x on coarse grid
      CHARACTER(len=256) :: routine_name
      !! name of this subroutine
      INTEGER :: ierr
      !! error code
      !
      EXTERNAL :: errore
      !
      routine_name = "spline_interpolation_rsca"
      !
      ! xmin
      xmin = x(1)
      ! xmax
      xmax = x(xdim)
      !
      ! first derivative
      delta = x(2) - x(1)
      IF ((xdim >= 2) .AND. (ABS(delta) > eps12)) THEN
        startu = (y(2) - y(1)) / delta
      ELSE
        startu = 0.0_dp
      END IF
      !
      ! second derivative
      delta = x(3) - x(1)
      IF ((xdim >= 3) .AND. (ABS(delta) > eps12)) THEN
        startd = (y(3) - 2.0_dp * y(2) + y(1)) / &
          delta
      ELSE
        startd = 0.0_dp
      END IF
      !
      ALLOCATE(d2y(xdim), STAT = ierr)
      IF (ierr /= 0) &
        CALL errore(routine_name, "Error allocating d2y", 1)
      !
      ! prepare spline
      CALL spline(x(1 : xdim), y(1 : xdim), startu, startd, d2y)
      !
      ! IF ((xf <= xmax) .AND. (xf >= xmin)) THEN
      !   yf = splint(x(1 : xdim), y(1 : xdim), d2y, xf)
      ! ELSE IF (xf >= xmax) THEN
      !   yf = y(xdim)
      ! ELSE
      !   IF (yf /= yf) THEN
      !     CALL errore(routine_name, 'yf value is NAN', 1)
      !   END IF
      !   !
      !   yf = y(1)
      ! END IF
      !
      IF ((xf <= xmax) .AND. (xf >= xmin)) THEN
        yf = splint(x(1 : xdim), y(1 : xdim), d2y, xf)
      ELSE
        CALL errore(routine_name, "interpolation out of bounds", 1)
      END IF
      !
      IF (yf /= yf) THEN
        CALL errore(routine_name, 'yf value is NAN', 1)
      END IF
      !
      IF (ABS(yf) < eps12) THEN
        yf = 0.0_dp
      END IF
      !
      DEALLOCATE(d2y, STAT = ierr)
      IF (ierr /= 0) &
        CALL errore(routine_name, 'Error deallocating d2y', 1)
      !
    !---------------------------------------------------------------------------
    END SUBROUTINE spline_interpolation_rsca
    !---------------------------------------------------------------------------
    !
    !
    !---------------------------------------------------------------------------
    SUBROUTINE spline_interpolation_rsca_range(imin, imax, x, y, xf, yf)
    !---------------------------------------------------------------------------
    !!
    !! Spline interpolation to real scalar.
    !!
    !---------------------------------------------------------------------------
      USE kinds, ONLY: DP
      USE splinelib, ONLY: spline, splint
      USE constants, ONLY: eps12
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(in) :: imin
      !! min coarse grid index
      INTEGER, INTENT(in) :: imax
      !! max coarse grid index
      REAL(DP), INTENT(in) :: x(:)
      !! coarse grid
      REAL(DP), INTENT(in) :: y(:)
      !! function values on coarse grid
      REAL(DP), INTENT(in) :: xf
      !! x on fine grid
      !
      REAL(DP), INTENT(out) :: yf
      !! function value on fine grid
      !
      REAL(DP), ALLOCATABLE :: d2y(:)
      REAL(DP) :: startu, startd
      REAL(DP) :: delta
      REAL(DP) :: xmin
      !! min value of x on coarse grid
      REAL(DP) :: xmax
      !! max value of x on coarse grid
      CHARACTER(len=256) :: routine_name
      !! name of this subroutine
      INTEGER :: ierr
      !! error code
      !
      EXTERNAL :: errore
      !
      routine_name = "spline_interpolation_rsca_range"
      !
      ! xmin
      xmin = x(imin)
      ! xmax
      xmax = x(imax)
      !
      ! first derivative
      delta = x(imin + 1) - x(imin)
      IF ((imax >= imin + 1) .AND. (ABS(delta) > eps12)) THEN
        startu = (y(imin + 1) - y(imin)) / delta
      ELSE
        startu = 0.0_dp
      END IF
      !
      ! second derivative
      delta = x(imin + 2) - x(imin)
      IF ((imax >= imin + 2) .AND. (ABS(delta) > eps12)) THEN
        startd = (y(imin + 2) - 2.0_dp * y(imin + 1) + y(imin)) / &
          delta
      ELSE
        startd = 0.0_dp
      END IF
      !
      ALLOCATE(d2y(imax - imin + 1), STAT = ierr)
      IF (ierr /= 0) &
        CALL errore(routine_name, "Error allocating d2y", 1)
      !
      ! prepare spline
      CALL spline(x(imin : imax), y(imin : imax), startu, startd, d2y)
      !
      IF ((xf <= xmax) .AND. (xf >= xmin)) THEN
        yf = splint(x(imin : imax), y(imin : imax), d2y, xf)
      ELSE
        CALL errore(routine_name, "interpolation out of specified bounds")
      END IF
      !
      IF (yf /= yf) THEN
        CALL errore(routine_name, 'yf value is NAN', 1)
      END IF
      !
      IF (ABS(yf) < eps12) THEN
        yf = 0.0_dp
      END IF
      !
      DEALLOCATE(d2y, STAT = ierr)
      IF (ierr /= 0) &
        CALL errore(routine_name, 'Error deallocating d2y', 1)
      !
    !---------------------------------------------------------------------------
    END SUBROUTINE spline_interpolation_rsca_range
    !---------------------------------------------------------------------------
    !
    !
    !
    !---------------------------------------------------------------------------
    PURE INTEGER FUNCTION rmta_get_index(n, x, x0)
    !---------------------------------------------------------------------------
    !!
    !! Returns mesh index of array x corresponding to x0
    !!
    !---------------------------------------------------------------------------
    !
    !  D. Radevych
    !
      USE kinds, ONLY : DP
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(in) :: n
      !! length of array x
      REAL(DP), INTENT(in) :: x(:)
      !! array x
      REAL(DP), INTENT(in) :: x0
      !! point
      INTEGER :: indx
      !! resulting index
      INTEGER :: i
      !! iterator
      !
      indx = 0
      !
      IF (n > 1) THEN
        IF (x0 > (x(n - 1) + x(n)) / 2.0_dp ) THEN
          indx = n
        ELSE
          DO i = 2, n - 1
            IF (((x(i) + x(i - 1)) / 2.0_dp < x0) .AND. &
              ((x(i) + x(i + 1)) / 2.0_dp) >= x0) THEN
              indx = i
            END IF
          END DO ! i
        END IF
      END IF
      !
      rmta_get_index = indx
      !
    !---------------------------------------------------------------------------
    END FUNCTION rmta_get_index
    !---------------------------------------------------------------------------
    !
    !
    !---------------------------------------------------------------------------
    SUBROUTINE rmta_quit()
    !---------------------------------------------------------------------------
    !!
    !! Cleans allocated rmta-specific arrays.
    !!
    !---------------------------------------------------------------------------
    !
    !  D. Radevych
    !
      USE io_global, ONLY : stdout
      USE sym_type, ONLY: deallocate_st
      USE neighbor, ONLY: deallocate_nn
      USE mt_var, ONLY: &
        rmta_delete_vars, rmta_routine
      !
      IMPLICIT NONE
      !
      CHARACTER(len=256) :: routine_name
      !! name of this subroutine
      !
      EXTERNAL :: errore, start_clock, stop_clock
      !
      routine_name = "rmta_quit"
      rmta_routine = "rmta"
      !
      CALL start_clock(routine_name)
      !
      !
      CALL deallocate_st()
      !
      CALL deallocate_nn()
      !
      CALL rmta_delete_vars()
      !
      !
      WRITE(stdout, '(5x, A, ": module-specific arrays are deallocated.")') &
        rmta_routine
      !
      !
      CALL stop_clock(routine_name)
      !
    !---------------------------------------------------------------------------
    END SUBROUTINE rmta_quit
    !---------------------------------------------------------------------------
    !
    !
    !
  !=============================================================================
  END MODULE muffin_tin
  !=============================================================================
