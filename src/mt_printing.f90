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
  MODULE mt_printing
  !=============================================================================
  !!
  !! Module printing/writing RMTA output.
  !!
  !=============================================================================
  !
  !  Danylo Radevych
  !
    USE kinds, ONLY: DP
    !
    IMPLICIT NONE
    !
    PUBLIC :: print_welcome_message, rmta_write, print_clocks, &
      print_at_rmt, check_input, check_upf
    !
  !
  CONTAINS
  !
    !
    !
    !
    !---------------------------------------------------------------------------
    SUBROUTINE print_welcome_message()
    !---------------------------------------------------------------------------
    !!
    !! Print welcome message.
    !!
    !---------------------------------------------------------------------------
    !
    !  D. Radevych
    !
      USE io_global, ONLY : ionode, stdout
      !
      IMPLICIT NONE
      !
      IF (ionode) THEN
        WRITE(stdout,'(/5x, "==================================", &
          "========================")')
        WRITE(stdout,'(/5x, ">>>>>>>>>>>>>>>>>>>>>>>>>> RMTA ", &
          "<<<<<<<<<<<<<<<<<<<<<<<<<<")')
        WRITE(stdout,'(/5x, "==================================", &
          "========================")')
        WRITE(stdout,'()')
        WRITE(stdout,'(/5x, "Please cite:")')
        WRITE(stdout,'(5x, "D. Radevych, T. Shishidou, M. Weinert, ")')
        WRITE(stdout,'(5x, "E. R. Margine, A. N. Kolmogorov, I. I. Mazin, ")')
        WRITE(stdout,'(5x, "Rigid muffin-tin approximation in plane wave ")')
        WRITE(stdout,'(5x, "codes for fast modeling of phonon-mediated ")')
        WRITE(stdout,'(5x, "superconductors, npj Comput Mater (2026) ")')
        WRITE(stdout,'(5x, "DOI: https://doi.org/10.1038/s41524-026-02141-7")')
        WRITE(stdout,'()')
      ENDIF
      !
    !---------------------------------------------------------------------------
    END SUBROUTINE print_welcome_message
    !---------------------------------------------------------------------------
    !
    !
    !---------------------------------------------------------------------------
    SUBROUTINE rmta_write()
    !---------------------------------------------------------------------------
    !!
    !! Prints quantities from SCF data for further plotting.
    !!
    !
    !  D. Radevych
    !
    ! TODO move potential init to rmta_init
    !
      USE fft_base,          ONLY : dfftp
      USE fft_types,         ONLY : fft_index_to_3d
      USE fft_rho,           ONLY : rho_g2r, rho_r2g
      USE gvect,             ONLY : ngl, gl, g, ngm, mill
      USE uspp_param,        ONLY : upf
      USE cell_base,         ONLY : at, bg, alat, tpiba, omega
      USE io_global,         ONLY : stdout
      USE ions_base,         ONLY : ityp
      USE mt_var, ONLY: mt_prec, natoms, norbs, orb_label, &
        n_chem_types, &
        nchis, nbetas, nspins, &
        mt_nrf, mt_rf, &
        mt_nr_max, mt_r, chir, betar, vlocionr, vlocaer, &
        vlocionr3d, &
        vlocscr00rf, vlocscf00rf, &
        vlociong, vlocscrr3d, &
        vlocscfr3d, vlocscrg3d, mt_g, vlocscfg3d, &
        vsemilocr, vsemilocrf, urf, duderf, dudrrf, &
        d2udrderf, &
        mll1rf_label, mll1rf, &
        vfullrf, rvfullrf, &
        irf_min, &
        dos_nlmrf, dos_nlrf, dos_nrf, &
        etall1rf, loglrf, dloglderf
      USE sym_type, ONLY: ist_i
      !
      !
      IMPLICIT NONE
      !
      REAL(DP) :: r
      !!
      INTEGER :: i, igl, iat, ir, j, k, ispin, igm, ig, iorb, &
        ict, im, ichi, ibeta
      !! iterators
      INTEGER :: if_chir = 1299
      !! unit for chir
      INTEGER :: if_betar = 1300
      !! unit for betar
      INTEGER :: if_vlocionr = 1301
      !! unit for vlocionr
      INTEGER :: if_vlocionr3d = 1302
      !! unit for vloctot
      INTEGER :: if_vlociong = 1303
      !! unit for vlocg
      INTEGER :: if_vlocscrr3d = 1304
      !! unit for vrs
      INTEGER :: if_vlocscr00rf = 1308
      !! unit for vlocscr00rf
      INTEGER :: if_vlocscfr3d = 1311
      !! unit for vlocscfr3d
      INTEGER :: if_vlocscf00rf = 1312
      !! unit for vlocscf00rf
      INTEGER :: if_vlocscfg3d = 1316
      !! unit for vlocscfg3d
      INTEGER :: if_vsemilocr = 1318
      !! unit for vsemilocr
      INTEGER :: if_vsemilocrf = 1319
      !! unit for vsemilocrf
      INTEGER :: if_urf = 1320
      !! unit for urf
      INTEGER :: if_vfullrf = 1321
      !! unit for vfullrf
      INTEGER :: if_rvfullrf = 1322
      !! unit for rvfullrf
      INTEGER :: if_vlocaer = 1324
      !! unit for vlocaer
      INTEGER :: if_duderf = 1325
      !! unit for duderf
      INTEGER :: if_dudrrf = 1326
      !! unit for dudrrf
      INTEGER :: if_d2udrderf = 1327
      !! unit for d2udrderf
      INTEGER :: if_mll1rf = 1328
      !! unit for mll1rf
      INTEGER :: if_vlocscrg3d = 1329
      !! unit for vlocscrg3d
      INTEGER :: if_dos_nlmrf = 1330
      !! unit for dos_nlmrf
      INTEGER :: if_loglrf = 1331
      !! unit for loglrf
      INTEGER :: if_dloglderf = 1332
      !! unit for dloglderf
      INTEGER :: if_dos_nlrf = 1333
      !! unit for dos_nlrf
      INTEGER :: if_dos_nrf = 1334
      !! unit for dos_nrf
      INTEGER :: if_dos_nlnrf = 1335
      !! unit for dos_nlrf / dos_nrf
      INTEGER :: if_etall1rf = 1336
      !! unit for etall1rf
      REAl(DP) :: r3d(3)
      !! current r3d vector
      REAl(DP) :: g3d(3)
      !! current g3d vector
      CHARACTER(LEN=256) :: stra, strs, strc, strm
      !! temporary strings for conversion
      INTEGER :: ishift(3) = 0
      !! integer shift
      CHARACTER(len=256) :: routine_name
      !! name of this subroutine
      INTEGER :: lmax
      !! max value of angular momentum
      INTEGER :: l
      !! current value of l angular momentum number
      INTEGER :: l0
      !!
      INTEGER :: m
      !! current value of m quantum number
      LOGICAL :: lskip_nonessential = .TRUE.
      !! if true, do not write to some files
      REAL(DP), ALLOCATABLE :: rtmp
      !! real temporary variable
      !
      EXTERNAL :: sph_bes, errore, start_clock, stop_clock
      !
      routine_name = "rmta_write"
      CALL start_clock(routine_name)
      !
      lmax = norbs - 1
      !
      !
      ! chir
      !
      WRITE(stdout, '(5x, A)') "writing chir"
      OPEN(UNIT = if_chir, FILE = TRIM('chir.dat'), &
        FORM = 'formatted', STATUS = 'unknown')
      WRITE(if_chir, '(1x, "chi(r) for each chemical type ")')
      !
      DO ict = 1, n_chem_types
        DO ichi = 1, nchis(ict)
          WRITE(if_chir, '(1x, A13)', advance='no') "r"
          WRITE(strc, '(I0)') ichi
          WRITE(if_chir, '(A14)', advance='no') &
            TRIM(upf(ict)%psd) // TRIM(strc)
        END DO ! ichi
      END DO ! ict
      WRITE(if_chir, '()')
      !
      DO ir = 1, mt_nr_max
        DO ict = 1, n_chem_types
          r = mt_r(ir, ict)
          DO ichi = 1, nchis(ict)
            WRITE(if_chir, '(1x, es13.4)', advance='no') r
            WRITE(if_chir, '(1x, es13.4)', advance='no') chir(ir, ichi, ict)
          END DO ! ichi
        END DO ! ict
        WRITE(if_chir, '()')
      END DO ! ir
      CLOSE(if_chir)
      !
      !
      ! betar
      !
      WRITE(stdout, '(5x, A)') "writing betar"
      OPEN(UNIT = if_betar, FILE = TRIM('betar.dat'), &
        FORM = 'formatted', STATUS = 'unknown')
      WRITE(if_betar, '(1x, "beta(r) for each chemical type")')
      DO ict = 1, n_chem_types
        DO ibeta = 1, nbetas(ict)
          WRITE(if_betar, '(1x, A13)', advance='no') "r"
          WRITE(strc, '(I0)') ibeta
          WRITE(if_betar, '(A14)', advance='no') &
            TRIM(upf(ict)%psd) // TRIM(strc)
        END DO ! ibeta
      END DO ! ict
      WRITE(if_betar, '()')
      !
      DO ir = 1, mt_nr_max
        DO ict = 1, n_chem_types
          r = mt_r(ir, ict)
          DO ibeta = 1, nbetas(ict)
            WRITE(if_betar, '(1x, es13.4)', advance='no') r
            WRITE(if_betar, '(1x, es13.4)', advance='no') betar(ir, ibeta, ict)
          END DO ! ibeta
        END DO ! ict
        WRITE(if_betar, '()')
      END DO ! ir
      CLOSE(if_betar)
      !
      !
      !
      ! vlocionr
      !
      WRITE(stdout, '(5x, A)') "writing vlocionr"
      OPEN(UNIT = if_vlocionr, FILE = TRIM('vlocionr.dat'), &
        FORM = 'formatted', STATUS = 'unknown')
      WRITE(if_vlocionr, '(1x, "V_{loc}^{ion}(r) for each chemical type")')
      DO ict = 1, n_chem_types
        WRITE(if_vlocionr, '(1x, A13)', advance='no') "r"
        WRITE(if_vlocionr, '(A14)', advance='no') TRIM(upf(ict)%psd)
      END DO ! ict
      WRITE(if_vlocionr, '()')
      DO ir = 1, mt_nr_max
        DO ict = 1, n_chem_types
          r = mt_r(ir, ict)
          !
          WRITE(if_vlocionr, '(1x, es13.4)', advance='no') r
          WRITE(if_vlocionr, '(1x, es13.4)', advance='no') vlocionr(ir, ict)
          !
        END DO ! ict
        WRITE(if_vlocionr, '()')
      END DO ! ir
      CLOSE(if_vlocionr)
      !
      !
      IF (.NOT. lskip_nonessential) THEN
        !
        ! vlocaer
        !
        WRITE(stdout, '(5x, A)') "writing vlocaer"
        OPEN(UNIT = if_vlocaer, FILE = TRIM('vlocaer.dat'), &
          FORM = 'formatted', STATUS = 'unknown')
        WRITE(if_vlocaer, '(1x, "V_{loc}^{ion}(r) for each chemical type ")')
        DO ict = 1, n_chem_types
          WRITE(if_vlocaer, '(1x, A13)', advance='no') "r"
          WRITE(if_vlocaer, '(A14)', advance='no') TRIM(upf(ict)%psd)
        END DO ! ict
        WRITE(if_vlocaer, '()')
        DO ir = 1, mt_nr_max
          DO ict = 1, n_chem_types
            r = mt_r(ir, ict)
            !
            WRITE(if_vlocaer, '(1x, es13.4)', advance='no') r
            WRITE(if_vlocaer, '(1x, es13.4)', advance='no') vlocaer(ir, ict)
            !
          END DO ! ict
          WRITE(if_vlocaer, '()')
        END DO ! ir
        CLOSE(if_vlocaer)
      END IF ! lskip_nonessential
      !
      !
      !
      ! vsemilocr
      !
      WRITE(stdout, '(5x, A)') "writing vsemilocr"
      OPEN(UNIT = if_vsemilocr, FILE = TRIM('vsemilocr.dat'), &
        FORM = 'formatted', STATUS = 'unknown')
      WRITE(stdout, '(5x, A)') "writing vsemilocr"
      WRITE(if_vsemilocr, '(1x, "V_sl(r) each chemical type")')
      DO ict = 1, n_chem_types
        DO iorb = 1, norbs
          WRITE(if_vsemilocr, '(1x, A13)', advance='no') "r"
          WRITE(strc, '(I0)') iorb
          WRITE(if_vsemilocr, '(A14)', advance='no') &
            TRIM(upf(ict)%psd) // TRIM(strc)
        END DO ! iorb
      END DO ! ict
      WRITE(if_vsemilocr, '()')
      DO ir = 1, mt_nr_max
        DO ict = 1, n_chem_types
          r = mt_r(ir, ict)
          DO iorb = 1, norbs
            !
            WRITE(if_vsemilocr, '(1x, es13.4)', advance='no') r
            WRITE(if_vsemilocr, '(1x, es13.4)', advance='no') &
              vsemilocr(ir, iorb, ict)
            !
          END DO ! iorb
        END DO ! ict
        WRITE(if_vsemilocr, '()')
      END DO ! ir
      CLOSE(if_vsemilocr)
      !
      !
      !
      !
      ! vlocscr00rf
      !
      WRITE(stdout, '(5x, A)') "writing vlocscr00rf"
      OPEN(UNIT = if_vlocscr00rf, FILE = TRIM('vlocscr00rf.dat'), &
        FORM = 'formatted', STATUS = 'unknown')
      WRITE(if_vlocscr00rf, '(1x, "V_{loc}(r) for each atom and spin")')
      DO iat = 1, natoms
        DO ispin = 1, nspins
          WRITE(if_vlocscr00rf, '(1x, A13)', advance='no') "r"
          WRITE(stra, '(I0)') iat
          WRITE(strs, '(I0)') ispin
          WRITE(if_vlocscr00rf, '(A14)', advance='no') &
            TRIM(upf(ityp(iat))%psd) // TRIM(stra) // &
            "s" // TRIM(strs)
        END DO ! ispin
      END DO ! iat
      WRITE(if_vlocscr00rf, '()')
      DO ir = 1, mt_nrf
        DO iat = 1, natoms
          DO ispin = 1, nspins
            r = mt_rf(ir, ist_i(iat))
            WRITE(if_vlocscr00rf, '(1x, es13.4)', advance='no') r
            WRITE(if_vlocscr00rf, '(1x, es13.4)', advance='no') &
              vlocscr00rf(ir, ispin, iat)
          END DO ! ispin
        END DO ! iat
        WRITE(if_vlocscr00rf, '()')
      END DO ! ir
      CLOSE(if_vlocscr00rf)
      !
      !
      WRITE(stdout, '(5x, A)') "writing vlocscf00rf"
      OPEN(UNIT = if_vlocscf00rf, FILE = TRIM('vlocscf00rf.dat'), &
        FORM = 'formatted', STATUS = 'unknown')
      WRITE(if_vlocscf00rf, '(1x, "V_{loc}^{scf}(r) for each atom and spin")')
      WRITE(if_vlocscf00rf, '(1x, A13)', advance='no') "r"
      DO iat = 1, natoms
        DO ispin = 1, nspins
          write (stra, '(I0)') iat
          write (strs, '(I0)') ispin
          WRITE(if_vlocscf00rf, '(A14)', advance='no') &
            TRIM(upf(ityp(iat))%psd) // TRIM(stra) // &
            "s" // TRIM(strs)
        END DO ! ispin
      END DO ! iat
      WRITE(if_vlocscf00rf, '()')
      DO ir = 1, mt_nrf
        DO iat = 1, natoms
          DO ispin = 1, nspins
            r = mt_rf(ir, ist_i(iat))
            WRITE(if_vlocscf00rf, '(1x, es13.4)', advance='no') r
            WRITE(if_vlocscf00rf, '(1x, es13.4)', advance='no') &
              vlocscf00rf(ir, ispin, iat)
          END DO ! ispin
        END DO ! iat
        WRITE(if_vlocscf00rf, '()')
      END DO ! ir
      CLOSE(if_vlocscf00rf)
      !
      !
      !
      WRITE(stdout, '(5x, A)') "writing vsemilocrf"
      OPEN(UNIT = if_vsemilocrf, FILE = TRIM('vsemilocrf.dat'), &
        FORM = 'formatted', STATUS = 'unknown')
      WRITE(if_vsemilocrf, '(1x, "V_{SL}(r) for each atom and l")')
      DO iat = 1, natoms
        DO iorb = 1, norbs
          ! WRITE(strl, '(I0)') iorb - 1
          WRITE(if_vsemilocrf, '(1x, A13)', advance='no') "r"
          write (stra, '(I0)') iat
          WRITE(if_vsemilocrf, '(A14)', advance='no') &
            TRIM(upf(ityp(iat))%psd) // &
            TRIM(stra) // &
            TRIM(orb_label(iorb))
        END DO ! iorb
      END DO ! iat
      WRITE(if_vsemilocrf, '()')
      DO ir = 1, mt_nrf
        DO iat = 1, natoms
          DO iorb = 1, norbs
            r = mt_rf(ir, ist_i(iat))
            WRITE(if_vsemilocrf, '(1x, es13.4)', advance='no') r
            WRITE(if_vsemilocrf, '(1x, es13.4)', advance='no') &
              vsemilocrf(ir, iorb, iat)
          END DO ! iorb
        END DO ! iat
        WRITE(if_vsemilocrf, '()')
      END DO ! ir
      CLOSE(if_vsemilocrf)
      !
      !
      !
      ! urf
      !
      WRITE(stdout, '(5x, A)') "writing urf"
      OPEN(UNIT = if_urf, FILE = TRIM('urf.dat'), &
        FORM = 'formatted', STATUS = 'unknown')
      WRITE(if_urf, '(1x, "u_l(r, e_f) for each atom, spin, and l")')
      DO iat = 1, natoms
        DO ispin = 1, nspins
          DO iorb = 1, norbs
            WRITE(if_urf, '(1x, A13)', advance='no') "r"
            WRITE(stra, '(I0)') iat
            WRITE(strs, '(I0)') ispin
            WRITE(if_urf, '(A14)', advance='no') &
              TRIM(upf(ityp(iat))%psd) // TRIM(stra) // &
              "s" // TRIM(strs) // &
              TRIM(orb_label(iorb))
          END DO ! iorb
        END DO ! ispin
      END DO ! iat
      WRITE(if_urf, '()')
      DO ir = 1, mt_nrf
        DO iat = 1, natoms
          DO ispin = 1, nspins
            DO iorb = 1, norbs
              r = mt_rf(ir, ist_i(iat))
              WRITE(if_urf, '(1x, es13.4)', advance='no') r
              WRITE(if_urf, '(1x, es13.4)', advance='no') &
                urf(ir, iorb, ispin, iat)
            END DO ! iorb
          END DO ! ispin
        END DO ! iat
        WRITE(if_urf, '()')
      END DO ! ir
      CLOSE(if_urf)
      !
      !
      ! duderf
      !
      WRITE(stdout, '(5x, A)') "writing duderf"
      OPEN(UNIT = if_duderf, FILE = TRIM('duderf.dat'), &
        FORM = 'formatted', STATUS = 'unknown')
      WRITE(if_duderf, '(1x, "dude_l(r, e_f) for each atom, spin, and l")')
      DO iat = 1, natoms
        DO ispin = 1, nspins
          DO iorb = 1, norbs
            WRITE(if_duderf, '(1x, A13)', advance='no') "r"
            WRITE(stra, '(I0)') iat
            WRITE(strs, '(I0)') ispin
            WRITE(if_duderf, '(A14)', advance='no') &
              TRIM(upf(ityp(iat))%psd) // TRIM(stra) // &
              "s" // TRIM(strs) // &
              TRIM(orb_label(iorb))
          END DO ! iorb
        END DO ! ispin
      END DO ! iat
      WRITE(if_duderf, '()')
      DO ir = 1, mt_nrf
        DO iat = 1, natoms
          DO ispin = 1, nspins
            DO iorb = 1, norbs
              r = mt_rf(ir, ist_i(iat))
              WRITE(if_duderf, '(1x, es13.4)', advance='no') r
              WRITE(if_duderf, '(1x, es13.4)', advance='no') &
                duderf(ir, iorb, ispin, iat)
            END DO ! iorb
          END DO ! ispin
        END DO ! iat
        WRITE(if_duderf, '()')
      END DO ! ir
      CLOSE(if_duderf)
      !
      !
      ! dudrrf
      !
      WRITE(stdout, '(5x, A)') "writing dudrrf"
      OPEN(UNIT = if_dudrrf, FILE = TRIM('dudrrf.dat'), &
        FORM = 'formatted', STATUS = 'unknown')
      WRITE(if_dudrrf, '(1x, "dudr_l(r, e_f) for each atom, spin, and l")')
      DO iat = 1, natoms
        DO ispin = 1, nspins
          DO iorb = 1, norbs
            WRITE(if_dudrrf, '(1x, A13)', advance='no') "r"
            WRITE(stra, '(I0)') iat
            WRITE(strs, '(I0)') ispin
            WRITE(if_dudrrf, '(A14)', advance='no') &
              TRIM(upf(ityp(iat))%psd) // TRIM(stra) // &
              "s" // TRIM(strs) // &
              TRIM(orb_label(iorb))
          END DO ! iorb
        END DO ! ispin
      END DO ! iat
      WRITE(if_dudrrf, '()')
      DO ir = 1, mt_nrf
        DO iat = 1, natoms
          DO ispin = 1, nspins
            DO iorb = 1, norbs
              r = mt_rf(ir, ist_i(iat))
              WRITE(if_dudrrf, '(1x, es13.4)', advance='no') r
              WRITE(if_dudrrf, '(1x, es13.4)', advance='no') &
                dudrrf(ir, iorb, ispin, iat)
            END DO ! iorb
          END DO ! ispin
        END DO ! iat
        WRITE(if_dudrrf, '()')
      END DO ! ir
      CLOSE(if_dudrrf)
      !
      !
      ! d2udrderf
      !
      WRITE(stdout, '(5x, A)') "writing d2udrderf"
      OPEN(UNIT = if_d2udrderf, FILE = TRIM('d2udrderf.dat'), &
        FORM = 'formatted', STATUS = 'unknown')
      WRITE(if_d2udrderf, '(1x, &
        "d2udrde_l(r, e_f) for each atom, spin, and l")')
      DO iat = 1, natoms
        DO ispin = 1, nspins
          DO iorb = 1, norbs
            WRITE(if_d2udrderf, '(1x, A13)', advance='no') "r"
            WRITE(stra, '(I0)') iat
            WRITE(strs, '(I0)') ispin
            WRITE(if_d2udrderf, '(A14)', advance='no') &
              TRIM(upf(ityp(iat))%psd) // TRIM(stra) // &
              "s" // TRIM(strs) // &
              TRIM(orb_label(iorb))
          END DO ! iorb
        END DO ! ispin
      END DO ! iat
      WRITE(if_d2udrderf, '()')
      DO ir = 1, mt_nrf
        DO iat = 1, natoms
          DO ispin = 1, nspins
            DO iorb = 1, norbs
              r = mt_rf(ir, ist_i(iat))
              WRITE(if_d2udrderf, '(1x, es13.4)', advance='no') r
              WRITE(if_d2udrderf, '(1x, es13.4)', advance='no') &
                d2udrderf(ir, iorb, ispin, iat)
            END DO ! iorb
          END DO ! ispin
        END DO ! iat
        WRITE(if_d2udrderf, '()')
      END DO ! ir
      CLOSE(if_d2udrderf)
      !
      !
      ! IF (.NOT. lskip_nonessential) THEN
      !
      ! mll1rf
      !
      WRITE(stdout, '(5x, A)') "writing mll1rf"
      OPEN(UNIT = if_mll1rf, FILE = TRIM('mll1rf.dat'), &
        FORM = 'formatted', STATUS = 'unknown')
      WRITE(if_mll1rf, '(1x, "M_{l, l + 1}(r) for each atom, spin, and l")')
      !WRITE(if_mll1rf, '(1x)', advance='no')
      DO iat = 1, natoms
        DO ispin = 1, nspins
          DO iorb = 1, norbs - 1
            WRITE(if_mll1rf, '(1x, A13)', advance='no') "r"
            WRITE(stra, '(I0)') iat
            WRITE(strs, '(I0)') ispin
            WRITE(if_mll1rf, '(1x, A13)', advance='no') &
              TRIM(upf(ityp(iat))%psd) // TRIM(stra) // &
              "s" // TRIM(strs) // &
              TRIM(mll1rf_label(iorb, ispin, iat))
          END DO ! iorb
        END DO ! ispin
      END DO ! iat
      WRITE(if_mll1rf, '()')
      DO ir = irf_min, mt_nrf ! only M vs. r_MT, not r, makes sense!
        DO iat = 1, natoms
          DO ispin = 1, nspins
            DO iorb = 1, norbs - 1
              r = mt_rf(ir, ist_i(iat))
              WRITE(if_mll1rf, '(1x, es13.4)', advance='no') r
              WRITE(if_mll1rf, '(1x, es13.4)', advance='no') &
                mll1rf(ir, iorb, ispin, iat)
            END DO ! iorb
          END DO ! ispin
        END DO ! iat
        WRITE(if_mll1rf, '()')
      END DO ! ir
      CLOSE(if_mll1rf)
      ! END IF ! lskip_nonessential
      !
      !
      !
      !
      ! vfullrf
      !
      WRITE(stdout, '(5x, A)') "writing vfullrf"
      OPEN(UNIT = if_vfullrf, FILE = TRIM('vfullrf.dat'), &
        FORM = 'formatted', STATUS = 'unknown')
      WRITE(if_vfullrf, '(1x, "V_{full}(r) for each atom, spin, and l")')
      DO iat = 1, natoms
        DO ispin = 1, nspins
          DO iorb = 1, norbs
            WRITE(if_vfullrf, '(1x, A13)', advance='no') "r"
            WRITE(stra, '(I0)') iat
            WRITE(strs, '(I0)') ispin
            WRITE(if_vfullrf, '(A14)', advance='no') &
              TRIM(upf(ityp(iat))%psd) // TRIM(stra) // &
              "s" // TRIM(strs) // &
              orb_label(iorb)
          END DO ! iorb
        END DO ! ispin
      END DO ! iat
      WRITE(if_vfullrf, '()')
      DO ir = 1, mt_nrf
        DO iat = 1, natoms
          DO ispin = 1, nspins
            DO iorb = 1, norbs
              r = mt_rf(ir, ist_i(iat))
              WRITE(if_vfullrf, '(1x, es13.4)', advance='no') r
              WRITE(if_vfullrf, '(1x, es13.4)', advance='no') &
                vfullrf(ir, iorb, ispin, iat)
            END DO ! iorb
          END DO ! ispin
        END DO ! iat
        WRITE(if_vfullrf, '()')
      END DO ! ir
      CLOSE(if_vfullrf)
      !
      !
      ! rvfullrf
      !
      WRITE(stdout, '(5x, A)') "writing rvfullrf"
      OPEN(UNIT = if_rvfullrf, FILE = TRIM('rvfullrf.dat'), &
        FORM = 'formatted', STATUS = 'unknown')
      WRITE(stdout, '(5x, A)') "writing rvfullrf"
      WRITE(if_rvfullrf, '(1x, "rV_{full}(r) for each atom, spin, and l")')
      DO iat = 1, natoms
        DO ispin = 1, nspins
          DO iorb = 1, norbs
            WRITE(if_rvfullrf, '(1x, A13)', advance='no') "r"
            WRITE(stra, '(I0)') iat
            WRITE(strs, '(I0)') ispin
            WRITE(if_rvfullrf, '(A14)', advance='no') &
              TRIM(upf(ityp(iat))%psd) // TRIM(stra) // &
              "s" // TRIM(strs) // &
              orb_label(iorb)
          END DO ! iorb
        END DO ! ispin
      END DO ! iat
      WRITE(if_rvfullrf, '()')
      DO ir = 1, mt_nrf
        DO iat = 1, natoms
          DO ispin = 1, nspins
            DO iorb = 1, norbs
              r = mt_rf(ir, ist_i(iat))
              WRITE(if_rvfullrf, '(1x, es13.4)', advance='no') r
              WRITE(if_rvfullrf, '(1x, es13.4)', advance='no') &
                rvfullrf(ir, iorb, ispin, iat)
            END DO ! iorb
          END DO ! ispin
        END DO ! iat
        WRITE(if_rvfullrf, '()')
      END DO ! ir
      CLOSE(if_rvfullrf)
      !
      !
      !
      ! dos_nlmrf
      !
      WRITE(stdout, '(5x, A)') "writing dos_nlmrf"
      OPEN(UNIT = if_dos_nlmrf, FILE = TRIM('dos_nlmrf.dat'), &
        FORM = 'formatted', STATUS = 'unknown')
      WRITE(if_dos_nlmrf, '(1x, "n_islm for each atom, spin, l, and m")')
      DO iat = 1, natoms
        DO ispin = 1, nspins
          DO iorb = 1, norbs
            !
            l = iorb - 1
            !
            DO im = 1, 2 * l + 1
              !
              m = im - l - 1
              !
              WRITE(if_dos_nlmrf, '(1x, A13)', advance='no') "r"
              WRITE(stra, '(I0)') iat
              WRITE(strs, '(I0)') ispin
              WRITE(strm, '(I0)') m
              WRITE(if_dos_nlmrf, '(A14)', advance='no') &
                TRIM(upf(ityp(iat))%psd) // TRIM(stra) // "_" // &
                "spin" // TRIM(strs) // "_" // &
                TRIM(orb_label(iorb)) // &
                TRIM(strm)
              !
            END DO ! im
          END DO ! iorb
        END DO ! ispin
      END DO ! iat
      WRITE(if_dos_nlmrf, '()')
      DO ir = irf_min, mt_nrf
        DO iat = 1, natoms
          DO ispin = 1, nspins
            DO iorb = 1, norbs
              !
              l = iorb - 1
              l0 = l * (l + 1) + 1
              !
              DO im = 1, 2 * l + 1
                !
                m = im - l - 1
                !
                r = mt_rf(ir, ist_i(iat))
                WRITE(if_dos_nlmrf, '(1x, es13.4)', advance='no') r
                WRITE(if_dos_nlmrf, '(1x, es13.4)', advance='no') &
                  dos_nlmrf(ir - irf_min + 1, l0 + m, ispin, iat)
              END DO ! im
            END DO ! iorb
          END DO ! ispin
        END DO ! iat
        WRITE(if_dos_nlmrf, '()')
      END DO ! ir
      CLOSE(if_dos_nlmrf)
      !
      !
      ! dos_nlrf
      !
      WRITE(stdout, '(5x, A)') "writing dos_nlrf"
      OPEN(UNIT = if_dos_nlrf, FILE = TRIM('dos_nlrf.dat'), &
        FORM = 'formatted', STATUS = 'unknown')
      WRITE(if_dos_nlrf, '(1x, "n_isl for each atom, spin, and l")')
      DO iat = 1, natoms
        DO ispin = 1, nspins
          DO iorb = 1, norbs
            !
            WRITE(if_dos_nlrf, '(1x, A13)', advance='no') "r"
            WRITE(stra, '(I0)') iat
            WRITE(strs, '(I0)') ispin
            WRITE(if_dos_nlrf, '(A14)', advance='no') &
              TRIM(upf(ityp(iat))%psd) // TRIM(stra) // "_" // &
              "spin" // TRIM(strs) // "_" // &
              TRIM(orb_label(iorb))
            !
          END DO ! iorb
        END DO ! ispin
      END DO ! iat
      WRITE(if_dos_nlrf, '()')
      DO ir = irf_min, mt_nrf
        DO iat = 1, natoms
          DO ispin = 1, nspins
            DO iorb = 1, norbs
              r = mt_rf(ir, ist_i(iat))
              WRITE(if_dos_nlrf, '(1x, es13.4)', advance='no') r
              WRITE(if_dos_nlrf, '(1x, es13.4)', advance='no') &
                dos_nlrf(ir - irf_min + 1, iorb, ispin, iat)
            END DO ! iorb
          END DO ! ispin
        END DO ! iat
        WRITE(if_dos_nlrf, '()')
      END DO ! ir
      CLOSE(if_dos_nlrf)
      !
      !
      ! dos_nrf
      !
      WRITE(stdout, '(5x, A)') "writing dos_nrf"
      OPEN(UNIT = if_dos_nrf, FILE = TRIM('dos_nrf.dat'), &
        FORM = 'formatted', STATUS = 'unknown')
      WRITE(if_dos_nrf, '(1x, "n_is for each atom, spin")')
      DO iat = 1, natoms
        DO ispin = 1, nspins
          !
          WRITE(if_dos_nrf, '(1x, A13)', advance='no') "r"
          WRITE(stra, '(I0)') iat
          WRITE(strs, '(I0)') ispin
          WRITE(if_dos_nrf, '(A14)', advance='no') &
            TRIM(upf(ityp(iat))%psd) // TRIM(stra) // "_" // &
            "spin" // TRIM(strs)
          !
        END DO ! ispin
      END DO ! iat
      WRITE(if_dos_nrf, '()')
      DO ir = irf_min, mt_nrf
        DO iat = 1, natoms
          DO ispin = 1, nspins
            r = mt_rf(ir, ist_i(iat))
            WRITE(if_dos_nrf, '(1x, es13.4)', advance='no') r
            WRITE(if_dos_nrf, '(1x, es13.4)', advance='no') &
              dos_nrf(ir - irf_min + 1, ispin, iat)
          END DO ! ispin
        END DO ! iat
        WRITE(if_dos_nrf, '()')
      END DO ! ir
      CLOSE(if_dos_nrf)
      !
      !
      ! etall1rf
      !
      WRITE(stdout, '(5x, A)') "writing etall1rf"
      OPEN(UNIT = if_etall1rf, FILE = TRIM('etall1rf.dat'), &
        FORM = 'formatted', STATUS = 'unknown')
      WRITE(if_etall1rf, '(1x, "eta_l for each atom, spin, and l")')
      DO iat = 1, natoms
        DO ispin = 1, nspins
          DO iorb = 1, norbs
            WRITE(if_etall1rf, '(1x, A13)', advance='no') "r"
            !
            WRITE(stra, '(I0)') iat
            WRITE(strs, '(I0)') ispin
            !
            IF (iorb < norbs) THEN
              WRITE(if_etall1rf, '(1x, A13)', advance='no') &
                TRIM(upf(ityp(iat))%psd) // TRIM(stra) // "_" // &
                "spin" // TRIM(strs) // "_" // &
                TRIM(mll1rf_label(iorb, ispin, iat))
            ELSE
              WRITE(if_etall1rf, '(1x, A13)', advance='no') &
                TRIM(upf(ityp(iat))%psd) // TRIM(stra) // "_" // &
                "spin" // TRIM(strs) // "_tot"
            END IF
          END DO ! iorb
        END DO ! ispin
      END DO ! iat
      WRITE(if_etall1rf, '()')
      DO ir = irf_min, mt_nrf
        DO iat = 1, natoms
          DO ispin = 1, nspins
            DO iorb = 1, norbs
              r = mt_rf(ir, ist_i(iat))
              WRITE(if_etall1rf, '(1x, es13.4)', advance='no') r
              WRITE(if_etall1rf, '(1x, es13.4)', advance='no') &
                etall1rf(ir - irf_min + 1, iorb, ispin, iat)
            END DO ! iorb
          END DO ! ispin
        END DO ! iat
        WRITE(if_etall1rf, '()')
      END DO ! ir
      CLOSE(if_etall1rf)
      !
      !
      !
      IF (.NOT. lskip_nonessential) THEN
        !
        ! vlociong
        !
        WRITE(stdout, '(5x, A)') "writing vlociong"
        OPEN(UNIT = if_vlociong, FILE = TRIM('vlociong.dat'), &
          FORM = 'formatted', STATUS = 'unknown')
        DO ict = 1, n_chem_types
          WRITE(if_vlociong, '(1x, "VG for type ", A)') TRIM(upf(ict)%psd)
          DO igl = 1, ngl
            WRITE(if_vlociong, '(1x, F10.6, "  ", F10.6)') &
              SQRT(gl(igl)) * tpiba, vlociong(igl, ict)
          END DO ! igl
        END DO ! ict
        CLOSE(if_vlociong)
      END IF ! lskip_nonessential
      !
      !
      !
      IF (.NOT. lskip_nonessential) THEN
        ! vlocionr3d
        WRITE(stdout, '(5x, A)') "writing vlocionr3d"
        WRITE(stdout, '(A)') "at"
        WRITE(stdout, '(F10.4, "  ", F10.4, "  ", F10.4)') at(:, 1) * alat
        WRITE(stdout, '(F10.4, "  ", F10.4, "  ", F10.4)') at(:, 2) * alat
        WRITE(stdout, '(F10.4, "  ", F10.4, "  ", F10.4)') at(:, 3) * alat
        WRITE(stdout, '(A)') "omega"
        WRITE(stdout, '(F10.4)') omega
        OPEN(UNIT = if_vlocionr3d, FILE = TRIM('vlocionr3d.dat'), &
          FORM = 'formatted', STATUS = 'unknown')
        WRITE(if_vlocionr3d, '(1x, "Vloc_tot(r1, r2, r3)")')
        WRITE(if_vlocionr3d, '(1x, I0, "  ", I0, "  ", I0)') &
            dfftp%nr1, dfftp%nr2, dfftp%nr3
        DO i = 1, dfftp%nr1
          DO j = 1, dfftp%nr2
            DO k = 1, dfftp%nr3
              ir = i + (j - 1) * dfftp%nr1 + (k - 1) * dfftp%nr1 * dfftp%nr2
              r3d(:) = &
                (i - 1) * at(:, 1) / dfftp%nr1 + &
                (j - 1) * at(:, 2) / dfftp%nr2 + &
                (k - 1) * at(:, 3) / dfftp%nr3
              r3d(:) = alat * r3d(:)
              WRITE(if_vlocionr3d, &
                '(1x, I3, "  ", I3, "  ", I3, "  ", &
                F10.6, "  ", F10.6, "  ", F10.6, "  ", F10.6)') &
                i, j, k, &
                r3d(1), r3d(2), r3d(3), &
                vlocionr3d(ir)
            END DO
          END DO
        END DO
        CLOSE(if_vlocionr3d)
      END IF ! lskip_nonessential
      !
      !
      !
      !
      IF (.NOT. lskip_nonessential) THEN
        !
        ! vlocscrr3d
        !
        WRITE(stdout, '(5x, A)') "writing vlocscrr3d"
        OPEN(UNIT = if_vlocscrr3d, FILE = TRIM('vlocscrr3d.dat'), &
          FORM = 'formatted', STATUS = 'unknown')
        DO ispin = 1, nspins
          WRITE(if_vlocscrr3d, '(1x, "Vrs for spin ", I0)') ispin
          WRITE(if_vlocscrr3d, '(1x, I0, "  ", I0, "  ", I0)') &
            dfftp%nr1, dfftp%nr2, dfftp%nr3
          DO i = 1, dfftp%nr1
            DO j = 1, dfftp%nr2
              DO k = 1, dfftp%nr3
                ir = i + (j - 1) * dfftp%nr1x + (k - 1) * dfftp%nr1x * dfftp%nr2x
                r3d(:) = &
                  (i - 1) * at(:, 1) / dfftp%nr1 + &
                  (j - 1) * at(:, 2) / dfftp%nr2 + &
                  (k - 1) * at(:, 3) / dfftp%nr3
                r3d(:) = alat * r3d(:)
                WRITE(if_vlocscrr3d, &
                  '(1x, I3, "  ", I3, "  ", I3, "  ", &
                  F10.6, "  ", F10.6, "  ", F10.6, "  ", F10.6)') &
                  i, j, k, &
                  r3d(1), r3d(2), r3d(3), &
                  vlocscrr3d(ir, ispin)
              END DO ! k
            END DO ! j
          END DO ! i
        END DO ! ispin
        CLOSE(if_vlocscrr3d)
      END IF ! lskip_nonessential
      !
      !
      IF (.NOT. lskip_nonessential) THEN
        !
        ! vlocscfr3d
        !
        WRITE(stdout, '(5x, A)') "writing vlocscfr3d"
        OPEN(UNIT = if_vlocscfr3d, FILE = TRIM('vlocscfr3d.dat'), &
          FORM = 'formatted', STATUS = 'unknown')
        DO ispin = 1, nspins
          WRITE(if_vlocscfr3d, '(1x, "Vrs for spin ", I0)') ispin
          WRITE(if_vlocscfr3d, '(1x, I0, "  ", I0, "  ", I0)') &
            dfftp%nr1, dfftp%nr2, dfftp%nr3
          DO i = 1, dfftp%nr1
            DO j = 1, dfftp%nr2
              DO k = 1, dfftp%nr3
                ir = i + (j - 1) * dfftp%nr1x + (k - 1) * dfftp%nr1x * dfftp%nr2x
                r3d(:) = &
                  (i - 1) * at(:, 1) / dfftp%nr1 + &
                  (j - 1) * at(:, 2) / dfftp%nr2 + &
                  (k - 1) * at(:, 3) / dfftp%nr3
                r3d(:) = alat * r3d(:)
                WRITE(if_vlocscfr3d, &
                  '(1x, I3, "  ", I3, "  ", I3, "  ", &
                  F10.6, "  ", F10.6, "  ", F10.6, "  ", F10.6)') &
                  i, j, k, &
                  r3d(1), r3d(2), r3d(3), &
                  vlocscfr3d(ir, ispin)
              END DO ! k
            END DO ! j
          END DO ! i
        END DO ! ispin
        CLOSE(if_vlocscfr3d)
      END IF ! lskip_nonessential
      !
      !
      IF (.NOT. lskip_nonessential) THEN
        !
        ! vlocscrg3d
        !
        WRITE(stdout, '(5x, A)') "writing vlocscrg3d"
        OPEN(UNIT = if_vlocscrg3d, FILE = TRIM('vlocscrg3d.dat'), &
          FORM = 'formatted', STATUS = 'unknown')
        DO ispin = 1, nspins
          WRITE(if_vlocscrg3d, '(1x, "V_{loc-tot}^{scr}(G) for spin ", I0)') &
            ispin
          ishift(1) = MAXVAL(mill(1, :))
          ishift(2) = MAXVAL(mill(2, :))
          ishift(3) = MAXVAL(mill(3, :))
          WRITE(if_vlocscrg3d, '(1x, I0, "  ", I0, "  ", I0)') &
            2 * ishift(1) + 1, &
            2 * ishift(2) + 1, &
            2 * ishift(3) + 1
          ishift(:) = ishift(:) + 1
          DO ig = 1, ngm
            g3d(:) = g(:, ig)
            g3d(:) = tpiba * g3d(:)
            WRITE(if_vlocscrg3d, &
              '(1x, I3, "  ", I3, "  ", I3, "  ", &
              F10.6, "  ", F10.6, "  ", F10.6, "  ", F10.6)') &
              mill(1, ig) + ishift(1), &
              mill(2, ig) + ishift(2), &
              mill(3, ig) + ishift(3), &
              g3d(1), g3d(2), g3d(3), &
              REAL(vlocscrg3d(ig, ispin))
            ! WRITE(if_vlocscrg3d, &
            !   '(1x, I3, "  ", I3, "  ", I3, "  ", &
            !   F10.6, "  ", F10.6, "  ", F10.6, "  ", F10.6, "  ", F10.6)') &
            !   mill(1, ig), mill(2, ig), mill(3, ig), &
            !   g3d(1), g3d(2), g3d(3), &
            !   REAL(vlocscrg3d(ig, ispin)), &
            !   AIMAG(vlocscrg3d(ig, ispin))
          END DO ! ig
        END DO ! ispin
        CLOSE(if_vlocscrg3d)
      END IF ! lskip_nonessential
      !
      !
      !
      IF (.NOT. lskip_nonessential) THEN
        !
        ! vlocscfg3d
        !
        WRITE(stdout, '(5x, A)') "writing vlocscfg3d"
        OPEN(UNIT = if_vlocscfg3d, FILE = TRIM('vlocscfg3d.dat'), &
          FORM = 'formatted', STATUS = 'unknown')
        DO ispin = 1, nspins
          WRITE(if_vlocscfg3d, '(1x, "Vgs for spin ", I0)') ispin
          WRITE(if_vlocscfg3d, '(1x, I0)') ngm
          WRITE(if_vlocscfg3d, '("   h    k    l", &
            "          gx          gy          gz", &
            "           g", &
            "      Re(Vg)      Im(Vg)")')
          DO igm = 1, ngm
            g3d(:) = g(:, igm)
            g3d(:) = tpiba * g3d(:)
            WRITE(if_vlocscfg3d, &
              '(1x, I3, "  ", I3, "  ", I3, "  ", &
              F10.6, "  ", &
              F10.6, "  ", F10.6, "  ", F10.6, "  ", F10.6, "  ", F10.6)') &
              mill(1, igm), mill(2, igm), mill(3, igm), &
              g3d(1), g3d(2), g3d(3), &
              mt_g(4, igm), &
              REAL(vlocscfg3d(igm, ispin)), &
              AIMAG(vlocscfg3d(igm, ispin))
          END DO ! igm
        END DO ! ispin
        CLOSE(if_vlocscfg3d)
      END IF ! lskip_nonessential
      !
      !
      !
      !
      !
      IF (.NOT. lskip_nonessential) THEN
        !
        ! loglrf
        !
        OPEN(UNIT = if_loglrf, FILE = TRIM('loglrf.dat'), &
          FORM = 'formatted', STATUS = 'unknown')
        WRITE(if_loglrf, '(1x, "L_l(r, e_f) for each atom, spin, and l")')
        WRITE(if_loglrf, '(1x, A13)', advance='no') "r"
        DO iat = 1, natoms
          DO ispin = 1, nspins
            DO iorb = 1, norbs
              WRITE(stra, '(I0)') iat
              WRITE(strs, '(I0)') ispin
              WRITE(if_loglrf, '(A14)', advance='no') &
                TRIM(upf(ityp(iat))%psd) // TRIM(stra) // &
                "s" // TRIM(strs) // &
                TRIM(orb_label(iorb))
            END DO ! iorb
          END DO ! ispin
        END DO ! iat
        WRITE(if_loglrf, '()')
        DO ir = 1, mt_nrf
          DO iat = 1, natoms
            DO ispin = 1, nspins
              DO iorb = 1, norbs
                r = mt_rf(ir, ist_i(iat))
                WRITE(if_loglrf, '(1x, es13.4)', advance='no') r
                WRITE(if_loglrf, '(1x, es13.4)', advance='no') &
                  loglrf(ir, iorb, ispin, iat)
              END DO ! iorb
            END DO ! ispin
          END DO ! iat
          WRITE(if_loglrf, '()')
        END DO ! ir
        CLOSE(if_loglrf)
        !
        !
        ! dloglderf
        !
        OPEN(UNIT = if_dloglderf, FILE = TRIM('dloglderf.dat'), &
          FORM = 'formatted', STATUS = 'unknown')
        WRITE(if_dloglderf, &
          '(1x, "d L_l(r, e_f) / d e for each atom, spin, and l")')
        WRITE(if_dloglderf, '(1x, A13)', advance='no') "r"
        DO iat = 1, natoms
          DO ispin = 1, nspins
            DO iorb = 1, norbs
              WRITE(stra, '(I0)') iat
              WRITE(strs, '(I0)') ispin
              WRITE(if_dloglderf, '(A14)', advance='no') &
                TRIM(upf(ityp(iat))%psd) // TRIM(stra) // &
                "s" // TRIM(strs) // &
                TRIM(orb_label(iorb))
            END DO ! iorb
          END DO ! ispin
        END DO ! iat
        WRITE(if_dloglderf, '()')
        DO ir = 1, mt_nrf
          DO iat = 1, natoms
            DO ispin = 1, nspins
              DO iorb = 1, norbs
                r = mt_rf(ir, ist_i(iat))
                WRITE(if_dloglderf, '(1x, es13.4)', advance='no') r
                WRITE(if_dloglderf, '(1x, es13.4)', advance='no') &
                  dloglderf(ir, iorb, ispin, iat)
              END DO ! iorb
            END DO ! ispin
          END DO ! iat
          WRITE(if_dloglderf, '()')
        END DO ! ir
        CLOSE(if_dloglderf)
        !
        !
        ! dos_nlrf / dos_nrf
        !
        OPEN(UNIT = if_dos_nlnrf, FILE = TRIM('dos_nlnrf.dat'), &
          FORM = 'formatted', STATUS = 'unknown')
        WRITE(if_dos_nlnrf, '(1x, "n_isl / n_is for each atom, spin, and l")')
        WRITE(if_dos_nlnrf, '(1x, A13)', advance='no') "r"
        DO iat = 1, natoms
          DO ispin = 1, nspins
            DO iorb = 1, norbs
              !
              WRITE(stra, '(I0)') iat
              WRITE(strs, '(I0)') ispin
              WRITE(if_dos_nlnrf, '(A14)', advance='no') &
                TRIM(upf(ityp(iat))%psd) // TRIM(stra) // "_" // &
                "spin" // TRIM(strs) // "_" // &
                TRIM(orb_label(iorb))
              !
            END DO ! iorb
          END DO ! ispin
        END DO ! iat
        WRITE(if_dos_nlnrf, '()')
        DO ir = 1, mt_nrf
          DO iat = 1, natoms
            DO ispin = 1, nspins
              DO iorb = 1, norbs
                rtmp = dos_nlrf(ir - irf_min + 1, iorb, ispin, iat)
                IF (rtmp < mt_prec) THEN
                  WRITE(if_dos_nlnrf, '(1x, es13.4)', advance='no') &
                    0._dp
                ELSE
                  WRITE(if_dos_nlnrf, '(1x, es13.4)', advance='no') &
                    dos_nlrf(ir - irf_min + 1, iorb, ispin, iat) / &
                    dos_nrf(ir - irf_min + 1, ispin, iat)
                END IF
              END DO ! iorb
            END DO ! ispin
          END DO ! iat
          WRITE(if_dos_nlnrf, '()')
        END DO ! ir
        CLOSE(if_dos_nlnrf)
        !
      END IF ! lskip_nonessential
      !
      !
      CALL stop_clock(routine_name)
      !
    !---------------------------------------------------------------------------
    END SUBROUTINE rmta_write
    !---------------------------------------------------------------------------
    !
    !
    !---------------------------------------------------------------------------
    SUBROUTINE print_clocks()
    !---------------------------------------------------------------------------
    !!
    !! Prints RMTA clocks.
    !!
      USE io_global, ONLY: stdout
      !
      IMPLICIT NONE
      !
      EXTERNAL :: print_clock
      !
      WRITE(stdout, *)
      !
      CALL print_clock('rmta_init')
      CALL print_clock('set_vsemiloc')
      CALL print_clock('set_log_ders')
      CALL print_clock('set_pet_mll1')
      CALL print_clock('tetra_delta_weights')
      CALL print_clock('set_dos_n')
      CALL print_clock('set_eta')
      CALL print_clock('print_at_rmt')
      CALL print_clock('rmta_write')
      CALL print_clock('spline_interpolation')
      CALL print_clock('rmta_quit')
      !
    !---------------------------------------------------------------------------
    END SUBROUTINE print_clocks
    !---------------------------------------------------------------------------
    !
    !
    !
    !---------------------------------------------------------------------------
    SUBROUTINE print_at_rmt()
    !---------------------------------------------------------------------------
    !!
    !! Printes RMTA results interpolated exactly at MT-radii
    !!
      USE io_global, ONLY: stdout
      USE ions_base, ONLY: ityp
      USE uspp_param, ONLY: upf
      USE constants, ONLY: rytoev, eps12
      USE const, ONLY: bohrtoang
      USE sym_type, ONLY: nst, ist_nat, ist_i, st_name
      USE mt_var, ONLY: natoms, norbs, orb_label, &
        nspins, fermi_energy, &
        luse_tot_dos, &
        irf_min, irf_max, &
        dos_nrf, dos_n, dos_nlrf, &
        vlocscr00rf, &
        mll1rf_label, mll1rf, etall1rf, &
        mt_nrf, mt_rf
      !
      IMPLICIT NONE
      !
      INTEGER :: ierr
      !! error code
      CHARACTER(len=256) :: routine_name
      !! name of this subroutine
      REAL(DP) :: rmtf
      !! MT-radius on fine grid
      REAL(DP) :: mll1
      !! current value of M_{l, l+1}
      CHARACTER(len=256) :: m_label
      !! label of M_{l, l+1}
      REAL(DP) :: veff
      !! effective potential on MT-sphere
      REAL(DP) :: ntot
      !! total DOS inside MT-sphere
      REAL(DP) :: nl
      !! n_{l} DOS inside MT-sphere
      REAL(DP) :: nl1
      !! n_{l+1} DOS inside MT-sphere
      REAL(DP) :: etall1
      !! \eta_{l, l+1}
      REAL(DP) :: eta
      !! \eta = sum_l \eta_{l, l+1}
      INTEGER :: iat, ispin, iorb, ist
      !! iterators
      INTEGER :: imin
      !! minimum index of defined partial DOS
      INTEGER :: imax
      !! maximum index of defined partial DOS
      REAL(DP), ALLOCATABLE :: mll1_sym_tp(:, :, :)
      !! value of m_{l, l + 1} for each type
      REAL(DP), ALLOCATABLE :: etall1_sym_tp(:, :, :)
      !! value of eta_{l, l + 1} for each type,
      !! last element etall1_sym_tp(norbs, :) is eta_tot
      REAL(DP) :: sum_n_dos
      !! sum of all DOS in all atoms
      !
      EXTERNAL :: errore, start_clock, stop_clock
      !
      !
      !
      routine_name = "print_at_rmt"
      CALL start_clock(routine_name)
      !
      imax = irf_max
      imin = irf_min
      !
      ALLOCATE(mll1_sym_tp(norbs, nspins, nst), STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, &
        'Error allocating mll1_sym_tp', 1)
      mll1_sym_tp(:, :, :) = 0._dp
      !
      ALLOCATE(etall1_sym_tp(norbs, nspins, nst), STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, &
        'Error allocating etall1_sym_tp', 1)
      etall1_sym_tp(:, :, :) = 0._dp
      !
      !
      WRITE(stdout, '(/5x, /5x, /5x, &
        ">>>>>>>>>", &
        ">>>>>>>>>>>>", &
        ">>>>>>>>>>>> OUTPUT BEGIN <<<<<<<<<<<<", &
        "<<<<<<<<<<<<", &
        "<<<<<<<<<")')
      !
      ! per atom
      !
      sum_n_dos = 0._dp
      !
      WRITE(stdout, '(/6x, "################################### ATOMS ######", &
        "###############################")')
      !
      DO iat = 1, natoms
        !
        WRITE(stdout, '(/6x, "Atom # ", I4, A4, &
          " ==============================================================")') &
          iat, upf(ityp(iat))%psd
        !
        rmtf = mt_rf(mt_nrf, ist_i(iat))
        WRITE(stdout, '(/6x, "r_mt = ", F8.6, " (bohr)")') rmtf
        !
        !
        DO ispin = 1, nspins
          !
          WRITE(stdout, '(/7x, "Spin: ", I1, &
            " --------------------------------------------------------------", &
            "-------")') &
          ispin
          !
          IF (.NOT. luse_tot_dos) THEN
            ! CALL spline_interpolation(imin, imax, mt_rf(:), &
            !   dos_nrf(:, ispin, iat), &
            !   rmtf, ntot)
            ntot = dos_nrf(mt_nrf - irf_min + 1, ispin, iat)
          ELSE
            ntot = dos_n(ispin)
          END IF
          sum_n_dos = sum_n_dos + ntot
          !
          veff = vlocscr00rf(mt_nrf, ispin, iat)
          !
          WRITE(stdout, '(/8x, A16, F10.4, A8, F8.4, A6)') &
            "E_F = ", fermi_energy(ispin), " (Ry) = ", &
            fermi_energy(ispin) * rytoev, " (eV)"
          WRITE(stdout, '(8x, A5, F10.5, A4, F10.4, A8, F10.4, A8)') &
            "V(", rmtf, ") = ", veff, &
            " (Ry) = ", veff * rytoev, " (eV)"
          WRITE(stdout, '(8x, A5, F10.5, A10, F10.4, A8, F10.4, A8)') &
            "V(", rmtf, ") - E_F = ", veff - fermi_energy(ispin), &
            " (Ry) = ", (veff - fermi_energy(ispin)) * rytoev, " (eV)"
          !
          WRITE(stdout, '(/8x, A16, F10.4)') &
            "! E_F:", fermi_energy(ispin)
          !
          WRITE(stdout, '(8x, A7, F10.4)') &
            "! V:", veff
          !
          WRITE(stdout, '(8x, A12, es14.4)') &
            "! N(E_F):", dos_n(ispin) * natoms
          !
          !
          !
          DO iorb = 1, norbs
            !
            !
            IF (iorb < norbs) THEN
              !
              WRITE(stdout, '(/8x, "l: ", A2, &
                " ........................................................", &
                ".......", ".......")') &
                orb_label(iorb)
              !
              !   CALL spline_interpolation(mt_nrf, mt_rf(:), &
              !     vlocscr00rf(:, ispin, iat), &
              !     rmtf, veff)
              !
              !
              m_label = mll1rf_label(iorb, ispin, iat)
              !
              WRITE(stdout, '("")')
              !   !
              !   ! interpolate M_{l, l+1} at r_mt (Ry / bohr)
              !   !
              !   CALL spline_interpolation(imax, mt_rf(:), &
              !     mll1rf(:, iorb, ispin, iat), &
              !     rmtf, mll1)
              !
              mll1 = mll1rf(mt_nrf, iorb, ispin, iat)
              !
              mll1_sym_tp(iorb, ispin, ist_i(iat)) = &
                mll1_sym_tp(iorb, ispin, ist_i(iat)) + &
                mll1 / ist_nat(ist_i(iat))
              !
              WRITE(stdout, '(8x, A10, A, A, F10.5, A, F14.4, A16)') &
                " M_", TRIM(m_label), &
                "(",  rmtf, "):", &
                mll1, " (Ry / bohr)"
              WRITE(stdout, '(8x, A10, A, A, F10.5, A, F14.4, A16)') &
                " M^2_", TRIM(m_label), &
                "(",  rmtf, "):", &
                mll1 * mll1, " (Ry / bohr)^2"
              !
              WRITE(stdout, '(8x, A10, A, A, F10.5, A, F14.4, A16)') &
                " M_", TRIM(m_label), &
                "(",  rmtf, "):", &
                mll1 * (rytoev / bohrtoang), " (eV / A)"
              WRITE(stdout, '(8x, A10, A, A, F10.5, A, F14.4, A16)') &
                " M^2_", TRIM(m_label), &
                "(",  rmtf, "):", &
                mll1 * mll1 * (rytoev / bohrtoang)**2, " (eV / A)^2"
              !
              ! WRITE(stdout, '(8x, A10, A, A, F10.5, A, F14.4)') &
              !   "! M^2_", TRIM(m_label), &
              !   "(",  rmtf, ") (Ry / bohr)^2:", &
              !   mll1 * mll1
              
              ! WRITE(stdout, '(5x, "Scaled:")')
              ! WRITE(stdout, '(7x, A10, A, A, F10.4, A, F14.4, A)') &
              !   " M_", TRIM(m_label), &
              !   "(",  rmtf, "):", &
              !   mll1 * &
              !   (rytoev / bohrtoang) / SQRT(rmta_wds4), &
              !   " (Wd^1/2 / S^2)"
              ! WRITE(stdout, '(7x, A10, A, A, F10.4, A, F14.4, A)') &
              !   " M^2_", TRIM(m_label), &
              !   "(",  rmtf, "):", &
              !   mll1 * &
              !   mll1 * &
              !   (rytoev / bohrtoang)**2 / rmta_wds4, &
              !   " (Wd / S^4)"
              !
              WRITE(stdout, '("")')
              !
              ! CALL spline_interpolation(imin, imax, mt_rf(:), &
              !   dos_nlrf(:, iorb, ispin, iat), &
              !   rmtf, nl)
              ! CALL spline_interpolation(imin, imax, mt_rf(:), &
              !   dos_nlrf(:, iorb + 1, ispin, iat), &
              !   rmtf, nl1)
              !
              nl = dos_nlrf(mt_nrf - irf_min + 1, &
                iorb, ispin, iat)
              nl1 = dos_nlrf(mt_nrf - irf_min + 1, &
                iorb + 1, ispin, iat)
              !
              ! WRITE(stdout, '(8x, A21, I1, A2, es14.4, A12)') &
              !   "N(E_F, ", ispin, "):", &
              !   dos_n(ispin) * natoms, " (1 / Ry)"
              WRITE(stdout, '(8x, A11, A, F10.5, A, es14.4)') &
                "n", &
                "(",  rmtf, ") (1 / Ry):", &
                ntot
              WRITE(stdout, '(8x, A10, A, A, F10.5, A, es14.4)') &
                "n_", TRIM(orb_label(iorb)), &
                "(",  rmtf, ") (1 / Ry):", &
                nl
              WRITE(stdout, '(8x, A10, A, A, F10.5, A, es14.4)') &
                "n_", TRIM(orb_label(iorb + 1)), &
                "(",  rmtf, ") (1 / Ry):", &
                nl1
              WRITE(stdout, '(8x, A10, A, A, es14.4)') &
                "n_", TRIM(orb_label(iorb)), &
                " / n : ", &
                nl / ntot
              WRITE(stdout, '(8x, A10, A, A, es14.4)') &
                "n_", TRIM(orb_label(iorb + 1)), &
                " / n : ", &
                nl1 / ntot
              !
              WRITE(stdout, '(8x, A10, A, A, es14.4)') &
                "! n_", TRIM(orb_label(iorb)), &
                ":", nl
              !
              IF (iorb == norbs - 1) THEN
                WRITE(stdout, '(8x, A10, A, A, es14.4)') &
                  "! n_", TRIM(orb_label(iorb + 1)), &
                  ":", nl1
              END IF
              !
              WRITE(stdout, '("")')
              !
              ! CALL spline_interpolation(imin, imax, mt_rf(:), &
              !   etall1rf(:, iorb, ispin, iat), &
              !   rmtf, etall1)
              !
              etall1 = etall1rf(mt_nrf - irf_min + 1, &
               iorb, ispin, iat)
              !
              ! etall1_sym_tp(iorb, ispin, ityp(iat)) = &
              !   etall1_sym_tp(iorb, ispin, ityp(iat)) + &
              !   etall1 / nat_per_chem_tp(ityp(iat))
              !
              etall1_sym_tp(iorb, ispin, ist_i(iat)) = &
                etall1_sym_tp(iorb, ispin, ist_i(iat)) + &
                etall1 / natoms
              !
              WRITE(stdout, '(8x, A7, A, A, F10.5, A, es14.4, A16, &
                es14.4, A)') &
                "eta_", TRIM(m_label), &
                "(",  rmtf, "):", &
                etall1, " (Ry / bohr^2) = ", &
                etall1 * rytoev / bohrtoang**2, " (eV / A^2)"
              ! WRITE(stdout, '(8x, A7, A, A, F10.5, A, &
              !   es14.4)') &
              !   "! eta_", TRIM(m_label), &
              !   "(",  rmtf, ") (eV / A^2):", &
              !   etall1 * rytoev / bohr**2
              !
              WRITE(stdout, '("")')
              !
            ELSE
              WRITE(stdout, '(8x, "sum_l: ", A, &
                " ........................................................", &
                ".......")') &
                "total"
              !
              WRITE(stdout, '("")')
              !
              ! CALL spline_interpolation(imin, imax, mt_rf(:), &
              !   etall1rf(:, norbs, ispin, iat), &
              !   rmtf, eta)
              !
              eta = etall1rf(mt_nrf - irf_min + 1, &
                norbs, ispin, iat)
              !
              ! etall1_sym_tp(norbs, ispin, ityp(iat)) = &
              !   etall1_sym_tp(norbs, ispin, ityp(iat)) + &
              !   eta / nat_per_chem_tp(ityp(iat))
              !
              etall1_sym_tp(norbs, ispin, ist_i(iat)) = &
                etall1_sym_tp(norbs, ispin, ist_i(iat)) + &
                eta / natoms
              !
              WRITE(stdout, '(8x, A7, A, F10.5, A, es14.4, A16, &
                es14.4, A)') &
                "eta_tot", &
                "(",  rmtf, "):", &
                eta, " (Ry / bohr^2) = ", &
                eta * rytoev / bohrtoang**2, " (eV / A^2)"
              ! WRITE(stdout, '(8x, A9, A, F10.4, A, &
              !   es14.4)') &
              !   "! eta_tot", &
              !   "(",  rmtf, ") (eV / A^2):", &
              !   eta * rytoev / bohrtoang**2
              !
              WRITE(stdout, '("")')
              !
            END IF
            !
            WRITE(stdout, '("")')
            !
          END DO ! iorb
          !
        END DO ! ispin
      END DO ! iat
      !
      IF (nspins == 1) sum_n_dos = sum_n_dos * 2.
      !
      WRITE(stdout, '("")')
      WRITE(stdout, '(/6x, "sum_n_dos: ", ES16.8, " (1 / Ry)")') sum_n_dos
      WRITE(stdout, '("")')
      WRITE(stdout, '("")')
      WRITE(stdout, '("")')
      !
      ! type-specific
      !
      WRITE(stdout, '(/6x, "########################## SYMMETRY TYPES ######", &
        "##############################")')
      DO ist = 1, nst
        !
        WRITE(stdout, '(/6x, "Symmetry type # ", I4, A4, &
          " =====================================================")') &
          ist, st_name(ist)
        rmtf = mt_rf(mt_nrf, ist)
        !
        !
        DO ispin = 1, nspins
          !
          WRITE(stdout, '(/7x, "Spin: ", I1, &
            " --------------------------------------------------------------", &
            "-------")') &
          ispin
          !
          DO iorb = 1, norbs
            !
            !
            IF (iorb < norbs) THEN
              !
              WRITE(stdout, '(/8x, "l: ", A2, &
                " ........................................................", &
                ".......", ".......")') &
                orb_label(iorb)
              !
              m_label = orb_label(iorb) // orb_label(iorb + 1)
              !
              WRITE(stdout, '("")')
              !
              ! M_{l, l + 1}
              !
              mll1 = mll1_sym_tp(iorb, ispin, ist)
              !
              WRITE(stdout, '(8x, A10, A, A, F10.5, A, F14.4, A16)') &
                " M_", TRIM(m_label), &
                "(",  rmtf, "):", &
                mll1, " (Ry / bohr)"
              WRITE(stdout, '(8x, A10, A, A, F10.5, A, F14.4, A16)') &
                " M^2_", TRIM(m_label), &
                "(",  rmtf, "):", &
                mll1 * mll1, " (Ry / bohr)^2"
              !
              WRITE(stdout, '(8x, A10, A, A, F10.5, A, F14.4, A16)') &
                " M_", TRIM(m_label), &
                "(",  rmtf, "):", &
                mll1 * (rytoev / bohrtoang), " (eV / A)"
              WRITE(stdout, '(8x, A10, A, A, F10.5, A, F14.4, A16)') &
                " M^2_", TRIM(m_label), &
                "(",  rmtf, "):", &
                mll1 * mll1 * (rytoev / bohrtoang)**2, " (eV / A)^2"
              !
              WRITE(stdout, '(8x, A10, A, A, F14.4)') &
                "! M^2_", TRIM(m_label), &
                ":", mll1 * mll1
              !
              WRITE(stdout, '("")')
              !
              !
              etall1 = etall1_sym_tp(iorb, ispin, ist)
              !
              WRITE(stdout, '(8x, A7, A, A, F10.5, A, es14.4, A16, &
                es14.4, A)') &
                "eta_", TRIM(m_label), &
                "(",  rmtf, "):", &
                etall1, " (Ry / bohr^2) = ", &
                etall1 * rytoev / bohrtoang**2, " (eV / A^2)"
              !
              WRITE(stdout, '(8x, A7, A, A, es14.4)') &
                "! eta_", TRIM(m_label), &
                ":", etall1 * rytoev / bohrtoang**2
              !
              WRITE(stdout, '("")')
              !
            ELSE
              WRITE(stdout, '(8x, "sum_l: ", A, &
                " ........................................................", &
                ".......")') &
                "total"
              !
              !
              eta = etall1_sym_tp(norbs, ispin, ist)
              !
              WRITE(stdout, '(8x, A7, A, F10.5, A, es14.4, A16, &
                es14.4, A)') &
                "eta_tot", &
                "(",  rmtf, "):", &
                eta, " (Ry / bohr^2) = ", &
                eta * rytoev / bohrtoang**2, " (eV / A^2)"
              !
              WRITE(stdout, '(8x, A12, es14.4)') &
                "! eta_tot:", &
                eta * rytoev / bohrtoang**2
              !
              WRITE(stdout, '("")')
              !
            END IF ! iorb
          END DO ! iorb
        END DO ! ispin
        !
        WRITE(stdout, '("")')
        !
      END DO ! ict
      !
      WRITE(stdout, '(/5x, &
        ">>>>>>>>>", &
        ">>>>>>>>>>>>", &
        ">>>>>>>>>>>>  OUTPUT END  <<<<<<<<<<<<", &
        "<<<<<<<<<<<<", &
        "<<<<<<<<<", &
        /5x, /5x, /5x)')
      !
      DEALLOCATE(mll1_sym_tp, STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, &
        'Error deallocating mll1_sym_tp', 1)
      !
      DEALLOCATE(etall1_sym_tp, STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, &
        'Error deallocating etall1_sym_tp', 1)
      !
      CALL stop_clock(routine_name)
      !
    !---------------------------------------------------------------------------
    END SUBROUTINE print_at_rmt
    !---------------------------------------------------------------------------
    !
    !
    !
    !---------------------------------------------------------------------------
    SUBROUTINE check_input(outdir)
    !---------------------------------------------------------------------------
    !!
    !! Checks input from the namelist.
    !!
    !
    !  D. Radevych
    !
      USE io_global, ONLY: stdout, ionode
      USE io_files, ONLY: prefix, tmp_dir
      USE mt_var, ONLY: lnonlocal
      !
      IMPLICIT NONE
      !
      EXTERNAL :: errore
      !
      CHARACTER(len=200) :: routine_name
      !! name of this subroutine
      CHARACTER(LEN=256), INTENT(IN) :: outdir
      CHARACTER(LEN=256), EXTERNAL :: trimcheck
      !
      routine_name = "check_input"
      !
      IF (ionode) THEN
        !
        WRITE(stdout,'(/5x, "Check input:")')
        WRITE(stdout,'(7x, "prefix = ", a)') prefix
        WRITE(stdout,'(7x, "outdir = ", a)') outdir
        WRITE(stdout,'(7x, "tmp_dir = ", a)') tmp_dir
        !
        IF (lnonlocal) &
          CALL errore(routine_name, "lnonlocal is not supported", 1)
        !
      ENDIF
      !
    !---------------------------------------------------------------------------
    END SUBROUTINE check_input
    !---------------------------------------------------------------------------
    !
    !
    !---------------------------------------------------------------------------
    SUBROUTINE check_upf()
    !---------------------------------------------------------------------------
    !!
    !! Checks upf info.
    !!
    !
    !  D. Radevych
    !
      USE uspp_param, ONLY : upf ! see upflib/pseudo_types.f90
      USE io_global, ONLY : stdout
      USE basis, ONLY: natomwfc
!       USE projections, ONLY: nlmchi
      USE ions_base, ONLY : ityp, nat
      USE ener, ONLY : ef
      !
      IMPLICIT NONE
      !
      INTEGER :: ict, j
      !! iterators
      !
      WRITE(stdout, '(/5x, "PP info:")')
      WRITE(stdout, '(7x, "SIZE(upf): ", I0, /)') SIZE(upf)
      DO ict = 1, SIZE(upf)
        WRITE(stdout, '(7x, "element: ", A)') upf(ict)%psd
        WRITE(stdout, '(7x, "z_valence: ", F10.6)') upf(ict)%zp
        WRITE(stdout, '(7x, "type: ", A)') upf(ict)%typ
        WRITE(stdout, '(7x, "functional: ", A)') upf(ict)%dft
        WRITE(stdout, '(7x, "relativistic: ", A)') upf(ict)%rel
        WRITE(stdout, '(7x, "is_ultrasoft: ", L)') upf(ict)%tvanp
        WRITE(stdout, '(7x, "is_paw: ", L)') upf(ict)%tpawp
        WRITE(stdout, '(7x, "is_coulomb: ", L)') upf(ict)%tcoulombp
        WRITE(stdout, '(7x, "has_so: ", L)') upf(ict)%has_so
        WRITE(stdout, '(7x, "has_wfc: ", L)') upf(ict)%has_wfc
        WRITE(stdout, '(7x, "has_gipaw: ", L)') upf(ict)%has_gipaw
        WRITE(stdout, '(7x, "paw_as_gipaw: ", L)') upf(ict)%paw_as_gipaw
        WRITE(stdout, '(7x, "core_correction: ", L)') upf(ict)%nlcc
        WRITE(stdout, '(7x, "with_metagga_info: ", L)') &
          upf(ict)%with_metagga_info
        WRITE(stdout, '(7x, "total_psenergy: ", F0.20)') upf(ict)%etotps
        WRITE(stdout, '(7x, "wfc_cutoff: ", F10.6)') upf(ict)%ecutwfc
        WRITE(stdout, '(7x, "rho_cutoff: ", F10.6)') upf(ict)%ecutrho
        ! maximum l component in beta
        WRITE(stdout, '(7x, "l_max: ", I0)') upf(ict)%lmax
        ! 2 * lmax
        WRITE(stdout, '(7x, "l_max_rho: ", I0)') upf(ict)%lmax_rho
        ! L of channel used to generate local potential
        ! (if < 0 it was generated by smoothing AE potential)
        WRITE(stdout, '(7x, "l_local: ", I0)') upf(ict)%lloc
        ! the maximum radius of the mesh
        WRITE(stdout, '(7x, "r_max: ", F10.6)') upf(ict)%rmax
        WRITE(stdout, '(7x, "number_of_wfc: ", I0)') upf(ict)%nwfc
        DO j = 1, upf(ict)%nwfc
          WRITE(stdout, '(9x, "rcut(", I3, "): ", F6.3, &
            "  rcutus(", I3, "): ", F6.3, &
            "  n(", I3, "): ", I3, &
            "  l(", I3, "): ", I3)') &
            j, upf(ict)%rcut(j), &
            j, upf(ict)%rcutus(j), &
            j, upf(ict)%nchi(j), &
            j, upf(ict)%lchi(j)
        END DO ! j
        WRITE(stdout, '(7x)')
        WRITE(stdout, '(7x, "nat: ", I0)') nat
        WRITE(stdout, '(7x, "natomwfc: ", I0)') natomwfc
        ! DO j = 1, natomwfc
        !   WRITE(stdout, '(9x, &
        !     "  n(", I3, "): ", I3, &
        !     "  l(", I3, "): ", I3, &
        !     "  m(", I3, "): ", I3)') &
        !     j, nlmchi(j)%n, &
        !     j, nlmchi(j)%l, &
        !     j, nlmchi(j)%m
        ! END DO ! j
        WRITE(stdout, '(7x)')
        WRITE(stdout, '(7x, "number_of_proj: ", I0)') upf(ict)%nbeta
        ! number of points in the radial mesh
        WRITE(stdout, '(7x, "mesh_size: ", I0)') upf(ict)%mesh
        !
        WRITE(stdout, '(7x)')
      END DO ! ict
      !
      !
      WRITE(stdout, '(/5x, "ef:", F10.6)') ef
      !
    !---------------------------------------------------------------------------
    END SUBROUTINE check_upf
    !---------------------------------------------------------------------------
    !
    !
  !=============================================================================
  END MODULE mt_printing
  !=============================================================================
