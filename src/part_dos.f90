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
  MODULE part_dos
  !=============================================================================
  !!
  !! Module for calculation of partial DOS at Fermi level
  !! inside a muffin-tin sphere
  !!
  !=============================================================================
  !
  !  D. Radevych
  !
    USE kinds, ONLY : DP
    !
    !
    IMPLICIT NONE
    !
    PUBLIC :: gauss_points, ylm4, set_dos_n, get_sum_wk, get_sum_wdk, &
      tetra_delta_weights
    !
    PRIVATE :: grule
    !
  !
  CONTAINS
  !
    !
    !---------------------------------------------------------------------------
    SUBROUTINE gauss_points(vgauss, wt, lmax)
    !---------------------------------------------------------------------------
    !!
    !! Generates Gaussian points (unit vectors corrsponding to specific
    !! angles phi and theta) to integrate spherical harmonics
    !! up to lmax exactly (!).
    !! <lm|l'm'> for l, l' <= lmax, number of points = (2*lmax+1)*(lmax+1))
    !!
    !---------------------------------------------------------------------------
    !
    ! Courtesy of M. Weinert and flair: FLAPW code.
    ! https://sites.uwm.edu/weinert/flair/
    !
      !
      IMPLICIT NONE
      !
      REAL(DP), INTENT(out) :: vgauss(:, :)
      !! points
      REAL(DP), INTENT(out) :: wt(:)
      !! weights
      INTEGER, INTENT(in) :: lmax
      !! max angular momentum
      !
      REAL(DP) :: xx(lmax / 2 + 1), w(lmax / 2 + 1)
      REAL(DP) :: delphi, phi, rxy
      INTEGER :: ngpt, nphi, i, j, k
      !
      ! determine the number of points cos(theta)
      !
      ngpt = lmax + 1
      CALL grule(ngpt, xx, w) ! outputs (ngpt + 1) / 2 points
      !
      ! in phi, use nyquist frequency, i.e.,  2 * lmax + 1
      !
      nphi = 2 * lmax + 1
      delphi = 8._dp * ATAN(1._dp) / nphi
      !
      j = 0
      DO i = 1, ngpt / 2
        rxy = SQRT(1._dp - xx(i) * xx(i))
        DO k=1, nphi
          phi = k * delphi
          j = j + 1
          vgauss(1, j) = rxy * COS(phi)
          vgauss(2, j) = rxy * SIN(phi)
          vgauss(3, j) = xx(i)
          wt(j) = w(i) * delphi
          j = j + 1
          vgauss(1, j) = vgauss(1, j - 1)
          vgauss(2, j) = vgauss(2, j - 1)
          vgauss(3, j) = -xx(i)
          wt(j) = w(i) * delphi
        ENDDO
      ENDDO
      !
      IF (MODULO(ngpt, 2) == 1) THEN
        DO k = 1, nphi
          vgauss(1, j + k) = COS(k * delphi)
          vgauss(2, j + k) = SIN(k * delphi)
          vgauss(3, j + k) = 0._dp
          wt(j + k) = w((ngpt + 1) / 2) * delphi
        ENDDO
        j = j + nphi
      END IF
      !
    !---------------------------------------------------------------------------
    END SUBROUTINE gauss_points
    !---------------------------------------------------------------------------
    !
    !
    !---------------------------------------------------------------------------
    SUBROUTINE grule(n, x, w)
    !---------------------------------------------------------------------------
    !!
    !! Determines the (n + 1) / 2 nonnegative points x(i) and
    !! the corresponding weights w(i) of the n-point
    !! Gauss-Legendre integration rule, normalized to the
    !! interval [-1, 1]. The x(i) appear in descending order.
    !!
    !! This routine is from 'Methods of Numerical Integration',
    !! P.J. Davis and P. Rabinowitz, p. 369.
    !---------------------------------------------------------------------------
    !
    ! Courtesy of M. Weinert and flair: FLAPW code.
    ! https://sites.uwm.edu/weinert/flair/
    !
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(in) :: n
      REAL(DP), INTENT(out) :: x(*), w(*)
      !
      INTEGER :: i, it, k, m
      REAL(DP) :: e1, p, pi, pk, pkm1, pkp1, t, t1, u, v, x0
      REAL(DP) :: den, d1, d_p, dpn, d2pn, d3pn, d4pn, fx, h
      !
      !
      pi = 4._dp * ATAN(1._dp)
      m = (n + 1) / 2
      e1 = n * (n + 1)
      !
      DO i = 1, m
        t = (4*i - 1) * pi / (4 * n + 2)
        x0 = (1._dp - (1._dp - 1._dp / n) / (8._dp * n * n)) * COS(t)
        !
        ! iterate on the value (M. Weinert Jan. 1982)
        !
        DO it=1,2
          pkm1 = 1._dp
          pk = x0
          !
          DO k = 2, n
            t1 = x0 * pk
            pkp1 = t1 - pkm1 - (t1 - pkm1) / k + t1
            pkm1 = pk
            pk = pkp1
          ENDDO
          !
          den = 1._dp - x0 * x0
          d1 = n * (pkm1 - x0 * pk)
          dpn = d1 / den
          d2pn = (2._dp * x0 * dpn - e1 * pk) / den
          d3pn = (4._dp * x0 * d2pn + (2._dp - e1) * dpn) / den
          d4pn = (6._dp * x0 * d3pn + (6._dp - e1) * d2pn) / den
          u = pk / dpn
          v = d2pn / dpn
          h = -u * (1._dp + .5_dp * u * &
            (v + u * (v * v - u * d3pn /(3._dp * dpn))))
          p = pk + h * (dpn + .5_dp * h * &
            (d2pn + h / 3._dp * (d3pn + .25_dp * h * d4pn)))
          d_p = dpn + h * (d2pn + .5_dp * h * (d3pn + h * d4pn / 3._dp))
          h = h - p / d_p
          x0 = x0 + h
        ENDDO
        x(i) = x0
        fx = d1 - h * e1 * (pk + .5 * h * &
          (dpn + h / 3. * (d2pn + .25 * h * (d3pn + .2 * h * d4pn))))
        w(i) = 2. * (1. - x(i) * x(i)) / (fx * fx)
      ENDDO
      !
      IF (m + m .GT. n) x(m) = 0.
      !
    !---------------------------------------------------------------------------
    END SUBROUTINE grule
    !---------------------------------------------------------------------------
    !
    !
    !---------------------------------------------------------------------------
    SUBROUTINE ylm4(v, ylm, lmax)
    !---------------------------------------------------------------------------
    !!
    !! Generate the spherical harmonics for the vector v
    !! using a stable upward recursion in l. (see notes
    !! by M. Weinert.)
    !! M. Weinert   January 1982
    !! modified by R. Podloucky (added in ynorm); July 1989
    !! cleaned up    mw 1995
    !!
    !!  Modified to make use of f90 constructs. note that
    !!  the normalization is an internal subroutine and hence
    !!  can only be called from here. also, no need to dimension
    !!  arrays for ynorm, done dynamically.          mw 1999
    !!
    !---------------------------------------------------------------------------
    !
    ! Courtesy of M. Weinert and flair: FLAPW code.
    ! https://sites.uwm.edu/weinert/flair/
    !
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(in) :: lmax
      REAL(DP), INTENT(in) :: v(3)
      COMPLEX(DP), INTENT(out) :: ylm(*)
      !
      REAL(DP), PARAMETER :: small = 1.0e-12_dp
      !
      COMPLEX(DP) :: ylms
      REAL(DP) :: p(0:lmax,0:lmax), c(0:lmax), s(0:lmax)
      REAL(DP) :: fac, x, y, z, xy, r, rxy, cth, sth, cph, sph, cph2
      INTEGER :: l, lm0, m
      !
      ! check whether  or not normalizations are needed
      !
      REAL(DP), ALLOCATABLE, SAVE :: ynorm(:)
      INTEGER, SAVE :: lmaxd = -1  ! initial value
      !
      IF (lmax .gt. lmaxd) THEN
        ! write(6,*) ' calling ylmnorm, lmax=',lmax,lmaxd
        ! first deallocate the array if it exists
        IF (ALLOCATED(ynorm)) DEALLOCATE(ynorm)
        ALLOCATE(ynorm((lmax+1)**2)) ! allocate array
        lmaxd = lmax
        CALL ylmnorm
      END IF
      !
      ! calculate sin and cos of theta and phi
      !
      x = v(1)
      y = v(2)
      z = v(3)
      xy = x * x + y * y
      r = SQRT(xy + z * z)
      rxy = SQRT(xy)

      IF (r .gt. small) THEN
         cth = z / r
         sth = rxy / r
      ELSE
         sth = 0._dp
         cth = 1._dp
      END IF
      IF (rxy .gt. small) THEN
         cph = x / rxy
         sph = y / rxy
      ELSE
         cph = 1._dp
         sph = 0._dp
      END IF
      !
      ! generate associated legendre functions for m.ge.0
      !
      fac = 1._dp
      ! loop over m values
      DO m = 0, lmax - 1
         fac = -(m + m - 1) * fac
         p(m, m) = fac
         p(m + 1, m) = (m + m + 1) * cth * fac
         ! recurse upward in l
         DO l = m + 2, lmax
            p(l, m)=((l + l - 1)* cth * &
              p(l - 1, m) - (l + m - 1) *p(l - 2, m)) / (l - m)
         END DO
         fac = fac * sth
      END DO
      p(lmax, lmax) = -(lmax + lmax - 1) * fac
      !
      ! determine sin and cos of phi
      s(0) = 0._dp
      s(1) = sph
      c(0) = 1._dp
      c(1) = cph
      cph2 = cph + cph
      DO m = 2, lmax
         s(m) = cph2 * s(m - 1) - s(m - 2)
         c(m) = cph2 * c(m - 1) - c(m - 2)
      END DO
      !
      ! multiply in the normalization factors
      DO l = 0, lmax
         ylm(l * (l + 1) + 1) = ynorm(l * (l + 1) + 1) * &
           CMPLX(p(l, 0), 0._dp, KIND=DP)
      END DO
      DO m = 1, lmax
         DO l = m, lmax
            lm0 = l * (l + 1) + 1
            ylms = p(l, m) * CMPLX(c(m), s(m), KIND=DP)
            ylm(lm0 + m) = ynorm(lm0 + m) * ylms
            ylm(lm0 - m) = CONJG(ylms) * ynorm(lm0 - m)
         END DO
      END DO
      !
      RETURN
    !
    CONTAINS ! INTERNAL SUBROUTINE
    !
      !.........................................................................
      SUBROUTINE ylmnorm
      !.........................................................................
      !
      ! normalization constants for ylm (internal subroutine has access
      ! to lmax and ynorm from above)
      !
      !.........................................................................
        !
        USE constants, ONLY: fpi
        !
        IMPLICIT NONE
        !
        INTEGER :: l,lm0,m
        REAL(DP) :: a, cd
        !
        !
        DO l=0,lmax
          lm0 = l * (l + 1) + 1
          cd = 1._dp
          a = SQRT((2 * l + 1) / fpi)
          ynorm(lm0) = a
          DO m = 1, l
            cd = cd / ((l + 1 - m) * (l + m))
            ynorm(lm0 + m) = a * SQRT(cd)
            ynorm(lm0 - m) = ((-1._dp)**m) * ynorm(lm0 + m)
          END DO
        END DO
        !
      !.........................................................................
      END SUBROUTINE ylmnorm
      !.........................................................................
      !
    !---------------------------------------------------------------------------
    END SUBROUTINE ylm4
    !---------------------------------------------------------------------------
    !
    !
    !---------------------------------------------------------------------------
    SUBROUTINE set_dos_n(ltetra, nr, imin, imax, stp, nat, norb, &
      nspin, ngauss, r, &
      tau_cart, dloglde, degauss, efermi, &
      dos_nlmr, dos_nlr, dos_nr, dos_n, &
      dos_nlmr_nodloglde, dos_nlr_nodloglde)
    !---------------------------------------------------------------------------
    !!
    !! Computes partial densities n at Fermi level for each atom i,
    !! spin s, and angular L = l, m on a range of r-indices [nmin, nmax]
    !!
    !---------------------------------------------------------------------------
    !
    !  D. Radevych
    !
      USE io_global, ONLY: stdout
      USE constants, ONLY: tpi, eps32, eps6, eps4
      USE cell_base, ONLY: tpiba, omega
      USE parameters, ONLY: npk
      USE klist, ONLY: xk, nkstot, ngk, wk, igk_k, two_fermi_energies
      USE wvfct, ONLY: npwx, et, nbnd
      ! USE wavefunctions, ONLY: evc
      USE io_files, ONLY: restart_dir
      USE pw_restart_new, ONLY: read_collected_wfc
      USE gvect, ONLY: g ! mill, gl, ngl
      USE symm_base, ONLY: nrot, irt, nosym
      USE lsda_mod, ONLY: isk
      ! USE lsda_mod, ONLY: nspin, isk, current_spin ! TODO
      ! USE ener, ONLY : ef
      !
      IMPLICIT NONE
      !
      EXTERNAL :: errore, start_clock, stop_clock, cryst_to_cart
      REAL(DP), EXTERNAL :: w0gauss
      ! EXTERNAL :: ylmr2 ! no good: produces REAL spherical harmonics
      !
      LOGICAL, INTENT(in) :: ltetra
      !! if true, use tetrahedron method
      INTEGER, INTENT(in) :: nr
      !! number of points on radial grid
      INTEGER, INTENT(in) :: imin
      !! min index of the radial point for partial DOS
      !! partial DOS for ir < nmin will be left zero
      INTEGER, INTENT(in) :: imax
      !! max index of the radial point for partial DOS
      !! partial DOS for ir > nmax will be left zero
      INTEGER, INTENT(in) :: stp(:)
      !! array converting atom indices into symmetry type indices, stp(atom)
      INTEGER, INTENT(in) :: nat
      !! number of atoms
      INTEGER, INTENT(in) :: norb
      !! number of orbitals
      INTEGER, INTENT(in) :: nspin
      !! number of spins
      INTEGER, INTENT(in) :: ngauss
      !! type of delta-function
      REAL(DP), INTENT(in) :: r(:, :)
      !! radial grid for each symmetry type, r(sym_type)
      REAL(DP), INTENT(in) :: tau_cart(:, :)
      !! atomic Cartesian coordinates
      REAL(DP), INTENT(in) :: dloglde(:, :, :, :)
      !! d L_l(r, e) / d e
      REAL(DP), INTENT(in) :: degauss
      !! smearing value
      REAL(DP), INTENT(in) :: efermi(:)
      !! Fermi energy
      !
      REAL(DP), INTENT(inout) :: dos_nlmr(:, :, :, :)
      !! partial densities n**i_{lm}(r, E_F)
      !! dos_nlmr(nr, (lmax + 1)**2, nspin, nat)
      !! lmax = norb - 1
      REAL(DP), INTENT(inout) :: dos_nlr(:, :, :, :)
      !! partial densities n**i_{l}(r, E_F)
      !! dos_nlr(nr, norb, nspin, nat)
      !! lmax = norb - 1
      REAL(DP), INTENT(inout) :: dos_nr(:, :, :)
      !! partial densities n**i(r, E_F), per atom, per spin
      !! dos_nr(nr, nspin, nat)
      !! lmax = norb - 1
      REAL(DP), INTENT(inout) :: dos_n(:)
      !! total DOS at the Fermi level n(E_F), per atom per spin
      !! dos_n(nspin)
      REAL(DP), INTENT(inout) :: dos_nlmr_nodloglde(:, :, :, :)
      !! reduced by dloglde partial densities n**i_{lm}(r, E_F)
      !! dos_nlmr(nr, (lmax + 1)**2, nspin, nat)
      !! lmax = norb - 1
      REAL(DP), INTENT(inout) :: dos_nlr_nodloglde(:, :, :, :)
      !! reduced by dloglde partial densities n**i_{l}(r, E_F)
      !! dos_nlr(nr, norb, nspin, nat)
      !! lmax = norb - 1
      !
      ! local variables
      !
      CHARACTER(len=256) :: routine_name
      !! name of this subroutine
      LOGICAL :: lspeedup = .TRUE.
      !! if true, ignore bands far on the tails of the delta-function
      !! (should be done)
      LOGICAL :: lselect
      !! auxiliary flag to include specific bands in integration
      LOGICAL :: lfound
      !! flag, if found
      LOGICAL :: lsymmetrize = .TRUE.
      !! symmatrize partial DOS
      LOGICAL :: lcheck_sym = .TRUE.
      !! check symmetry relations
      INTEGER :: ierr
      !! error code
      INTEGER :: iat, jat, iorb, im, ir, ispin, igp, ik, ig, ibnd, irot
      !! iterators
      INTEGER :: lmax
      !! max angular momentum
      INTEGER :: gp_ntheta
      !! number of Gauss points for theta
      INTEGER :: gp_nphi
      !! number of Gauss points for phi
      INTEGER :: ngp
      !! number of all Gauss points
      INTEGER :: npw
      !! number of plane-waves
      INTEGER :: l
      !! current ang momentum number
      INTEGER :: l0
      !! l0
      INTEGER :: m
      !! current m quantum number
      INTEGER, ALLOCATABLE :: counters(:)
      !! counters on number of equivalent sites
      REAL(DP) :: deltaf
      !! current value of Dirac-delta function
      REAL(DP) :: prefactor_part_dos
      !! constant prefactor for partial DOS
      REAL(DP) :: prefactor
      !! common prefactor for given rmt, atom, and l
      REAL(DP) :: prefactor_tot_dos
      !! constant prefactor for total DOS
      REAL(DP) :: prefactor_nodloglde
      !! reduced by dloglde common prefactor for given rmt, atom, and l
      REAL(DP) :: rvec_cart(3)
      !! temporary MT-sphere vector for Gauss integration
      !! in Cartesian coordinates
      REAL(DP) :: kvec_cart(3)
      !! temporary k-point vector
      !! in Cartesian coordinates
      REAL(DP) :: gvec_cart(3)
      !! temporary G-vector
      !! in Cartesian coordinates
      REAL(DP) :: arg
      !! temporary argument of the exponent
      REAL(DP) :: sum_wk
      !! sum of all k-point weights
      REAL(DP) :: sum_wdk
      !! sum of all k-point and band weights (ltetra = .true.)
      REAL(DP) :: psi_kg_norm
      !! norm of psi(k + G) coefficients
      REAL(DP) :: avg_dos_per_spin
      !! average of the total DOS for both spins
      REAL(DP), ALLOCATABLE :: gp_vec(:, :)
      !! unit vectors for Gauss integration
      REAL(DP), ALLOCATABLE :: gp_wt(:)
      !! weights for Gauss integration
      REAL(DP), ALLOCATABLE :: wdk(:, :)
      !! tetrahedron weights for integration with the delta-function
      REAL(DP), ALLOCATABLE :: dos_nlmr_sym(:, :, :, :)
      !! symmetrized partial densities n**i_{lm}(r, E_F)
      !! for a specific atom
      !! dos_nlmr_sym(nr, 2 * lmax + 1, norb, nspin, nat)
      !! lmax = norb - 1
      REAL(DP), ALLOCATABLE :: dos_nlr_sym(:, :, :, :)
      !! symmetrized  partial densities n**i_{l}(r, E_F)
      !! for a specific atom
      !! dos_nlr_sym(nr, norb, nspin, nat)
      !! lmax = norb - 1
      REAL(DP), ALLOCATABLE :: dos_nlmr_nodloglde_sym(:, :, :, :)
      !! reduced by dloglde symmetrized partial densities n**i_{lm}(r, E_F)
      !! for a specific atom
      !! dos_nlmr_sym(nr, 2 * lmax + 1, norb, nspin, nat)
      !! lmax = norb - 1
      REAL(DP), ALLOCATABLE :: dos_nlr_nodloglde_sym(:, :, :, :)
      !! reduced by dloglde symmetrized  partial densities n**i_{l}(r, E_F)
      !! for a specific atom
      !! dos_nlr_sym(nr, norb, nspin, nat)
      !! lmax = norb - 1
      COMPLEX(DP) :: cnr_aux
      !! complex version of partial DOS on r grid
      COMPLEX(DP), ALLOCATABLE :: psi_kg(:, :, :)
      !! FT coefficients of SCF wavefunctions for all bands and kpoints
      !! psi_kg(npwx, nbnd, nks)
      COMPLEX(DP), ALLOCATABLE :: psi_krtau_aux(:, :, :)
      !! temporary value of
      !! \delta(\varepsilon_{\bm{k}, i} - E_F)
      !! \sum_{\bm{G}} ef**{i (\bm{k} + \bm{G}) \cdot (\bm{r} + \bm{tau})}
      !! psi_i(\bm{k} + \bm{G}) =
      !! \delta(\varepsilon_{\bm{k}, i} - E_F)
      !! psi_{\bm{k}, i}(\bm{r} + \bm{\tau})
      !! on the given MT sphere
      COMPLEX(DP), ALLOCATABLE :: ylm(:, :)
      !! temporary array containing spherical harmonics
      !! for specific unit vector and l
      !
      !
      !
      !
      routine_name = "set_dos_n"
      CALL start_clock(routine_name)
      !
      IF (ltetra) THEN
        !
        ALLOCATE(wdk(nbnd, nkstot), STAT = ierr)
        IF (ierr /= 0) CALL errore(routine_name, 'Error allocating wdk', 1)
        !
        DO ispin = 1, nspin
          !
          CALL tetra_delta_weights(nkstot, nspin, ispin, isk, nbnd, &
            efermi(ispin), et, wdk)
          !
        END DO ! ispin
        !
      END IF
      !
      IF (imax > nr) &
        CALL errore(routine_name, "imax > nr", 1)
      !
      lsymmetrize = (.NOT. nosym)
      lcheck_sym = (.NOT. nosym)
      !
      WRITE(stdout, '(/5x, ">>>>>>>>>>   PartDOS BEGIN   <<<<<<<<<<<")')
      !
      WRITE(stdout, '(/7x, "ltetra = ", L2)') ltetra
      WRITE(stdout, '(7x, "two_fermi_energies = ", L2)') two_fermi_energies
      WRITE(stdout, '(7x, "nspin = ", I0)') nspin
      IF (.NOT. ltetra) THEN
        WRITE(stdout, '(7x, "ngauss = ", I0)') ngauss
        WRITE(stdout, '(7x, "degauss = ", F10.4)') degauss
      END IF
      WRITE(stdout, '(7x, "nkstot = ", I0)') nkstot
      WRITE(stdout, '(7x, "nbnd = ", I0)') nbnd
      WRITE(stdout, '(7x, "omega = ", F10.4)') omega
      !
      DO ik = 1, nkstot, nkstot - 1
        WRITE(stdout, '(7x, "isk(", I6, ") = ", I1)') ik, isk(ik)
      END DO
      !
      sum_wk = get_sum_wk()
      WRITE(stdout, '(7x, "sum_wk = ", F10.4, /7x)') sum_wk
      !
      IF (ltetra) THEN
        sum_wdk = get_sum_wdk(wdk)
        WRITE(stdout, '(7x, "sum_wdk = ", F10.4, /7x)') sum_wdk
      END IF
      !
      lmax = norb - 1
      !
      !
      ! generate Gauss-integration points for integration with
      ! spherical harmonics
      !
      gp_ntheta = lmax + 1
      gp_nphi = 2 * lmax + 1
      ngp = gp_ntheta * gp_nphi
      !
      ALLOCATE(gp_vec(3, ngp), STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, 'Error allocating gp_vec', 1)
      !
      ALLOCATE(gp_wt(ngp), STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, 'Error allocating gp_wt', 1)
      !
      CALL gauss_points(gp_vec, gp_wt, lmax)
      !
      !
      ! prepare corresponding spherical harmonics for
      ! Gauss-point integration
      !
      ALLOCATE(ylm((lmax + 1) * (lmax + 1), ngp), STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, 'Error allocating ylm', 1)
      !
      ylm(:, :) = CMPLX(0._dp, 0._dp, KIND=DP)
      !
      DO igp = 1, ngp
        CALL ylm4(gp_vec(:, igp), ylm(:, igp), lmax)
      END DO
      !
      !
      !
      ! read and store all psi_kg (evc) coefficients of the wavefunctions
      !
      WRITE(stdout, &
        '(/6x, "Reading all stored psi_kg ", &
        & "coefficients of the wavefunctions...")')
      !
      ALLOCATE(psi_kg(npwx, nbnd, nkstot), STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, 'Error allocating psi_kg', 1)
      !
      !
      DO ik = 1, nkstot
        !
        ! WRITE(stdout, '(/5x, "k-point #", I)') ik
        !
        ! npw = ngk(ik)
        ! CALL read_collected_wfc(restart_dir(), ik, evc)
        CALL read_collected_wfc(restart_dir(), ik, psi_kg(:, :, ik))
      END DO ! ik
      !
      WRITE(stdout, '(6x, "Done reading all stored psi_kg ", &
        & "coefficients of the wavefunctions.", /6x, /6x)')
      !
      psi_kg_norm = 0._dp
      DO ik = 1, nkstot
        DO ibnd = 1, nbnd
          npw = ngk(ik)
          DO ig = 1, npw
            psi_kg_norm = psi_kg_norm + &
              REAL(psi_kg(ig, ibnd, ik) * CONJG(psi_kg(ig, ibnd, ik)))
          END DO ! ig
        END DO ! ibnd
      END DO ! ik
      psi_kg_norm = psi_kg_norm / nkstot / nbnd
      WRITE(stdout, '(7x, "psi_kg_norm = ", F10.4, /7x)') psi_kg_norm
      !
      !
      !
      ! partial DOS = f(r, E_F) integration
      !
      WRITE(stdout, &
        '(/6x, "Computing partial DOS per atom per spin = f(r, E_F)...")')
      !
      dos_nlmr(:, :, :, :) = 0._dp
      dos_nlr(:, :, :, :) = 0._dp
      dos_nr(:, :, :) = 0._dp
      !
      dos_nlmr_nodloglde(:, :, :, :) = 0._dp
      dos_nlr_nodloglde(:, :, :, :) = 0._dp
      !
      !
      ! spin!
      IF (ltetra) THEN
        prefactor_part_dos = 1._dp / omega
        prefactor_tot_dos = 1._dp / nat
        !
        If (nspin == 1) THEN
          prefactor_part_dos = prefactor_part_dos / 2._dp
          prefactor_tot_dos = prefactor_tot_dos / 2._dp
        END IF
        !
      ELSE
        ! times nspin to compensate sum_wk = 2
        prefactor_part_dos = 1._dp / omega / sum_wk * nspin
        prefactor_tot_dos = 1._dp / sum_wk / nat * nspin
      END IF
      !
      !
      ! prefactor0 = 1._dp
      ! prefactor0 = 1.0 / omega / (tpi * tpi * tpi)
      ! prefactor0 = prefactor0 * prefactor0 * prefactor0
      !
      !
      ALLOCATE(psi_krtau_aux(ngp, nbnd, nkstot), STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, &
        'Error allocating psi_krtau_aux', 1)
      !
      DO iat = 1, nat
        !
        WRITE(stdout, '(/6x, "atom #", I0, " out of ", I0)') iat, nat
        !
        DO ispin = 1, nspin
          !
          WRITE(stdout, '(/7x, "spin #", I0, " out of ", I0)') ispin, nspin
          !
          DO ir = imin, imax
            !
            WRITE(stdout, '(8x, "ir #", I0, " on [", I0, ", ", I0, "]")') &
             ir, imin, imax
            !
            ! precompute psi_krtau with delta function for all bands
            !
            psi_krtau_aux(:, :, :) = CMPLX(0._dp, 0._dp, KIND=DP)
            !
            WRITE(stdout, &
              '(/9x, "Precomputing ", &
              & "psi_krtau_aux...")')
            !
            DO ik = 1, nkstot
              !
              ! only k-points of the corresponding spin
              !
              IF (isk(ik) == ispin) THEN
                !
                npw = ngk(ik)
                !
                kvec_cart(:) = tpiba * xk(:, ik)
                ! WRITE(*, *) "xk(:, ik)", xk(:, ik)
                !
                ! loop over KS-states
                !
                DO ibnd = 1, nbnd
                  !
                  IF (ltetra) THEN
                    lselect = ABS(wdk(ibnd, ik)) > eps32
                    ! lselect = .true.
                  ELSE
                    deltaf = &
                      w0gauss((et(ibnd, ik) - efermi(ispin)) / &
                      degauss, ngauss) / degauss
                      ! w0gauss((et(ibnd, ik) - efermi) / degauss, ngauss)
                    lselect = ABS(deltaf) > eps32 .OR. (.NOT. lspeedup)
                  END IF
                  !
                  IF (lselect) THEN
                    !
                    DO igp = 1, ngp
                      !
                      rvec_cart(:) = r(ir, stp(iat)) * gp_vec(:, igp)
                      !
                      ! sum over G-vectors
                      !
                      DO ig = 1, npw
                        !
                        gvec_cart(:) = tpiba * g(:, igk_k(ig, ik))
                        !
                        ! WRITE(*, *) "g(:, igk_k(ig, ik))", g(:, igk_k(ig, ik))
                        !
                        arg = DOT_PRODUCT(kvec_cart(:) + gvec_cart(:), &
                          rvec_cart(:) + tau_cart(:, iat))
                        !
                        ! TEST (naive, against all equations)
                        ! arg = DOT_PRODUCT(kvec_cart(:) - gvec_cart(:), &
                        !   rvec_cart(:) + tau_cart(:, iat))
                        !
                        psi_krtau_aux(igp, ibnd, ik) = &
                          psi_krtau_aux(igp, ibnd, ik) + &
                          CMPLX(COS(arg), SIN(arg), KIND=DP) * &
                          psi_kg(ig, ibnd, ik)
                        !
                      END DO ! ig
                      !
                    END DO ! igp
                    !
                  END IF ! lselect
                  !
                END DO ! ibnd
                !
              END IF ! isk
              !
            END DO ! ik
            !
            !
            WRITE(stdout, '(9x, "Done precomputing ", &
              & "psi_krtau_aux.", /6x)')
            !
            DO iorb = 1, norb
              !
              l = iorb - 1
              l0 = l * (l + 1) + 1
              !
              ! WRITE(stdout, '(/8x, "orbit #", I0, " out of ", I0)') l, norb - 1
              !
              DO im = 1, 2 * l + 1
                !
                m = im - l - 1
                !
                ! WRITE(stdout, '(/9x, "m #", I)') m
                !
                !
                prefactor = prefactor_part_dos * r(ir, stp(iat)) * &
                  ABS(dloglde(ir, iorb, ispin, iat))
                prefactor_nodloglde = prefactor_part_dos * r(ir, stp(iat))
                !
                ! cnr_aux(:, :) = CMPLX(0._dp, 0._dp, KIND=DP)
                !
                ! integral over k-points
                !
                DO ik = 1, nkstot
                  !
                  IF (isk(ik) == ispin) THEN
                  !
                  npw = ngk(ik)
                  !
                  !
                  ! loop over KS-states
                  !
                  DO ibnd = 1, nbnd
                    !
                    deltaf = 0.0_dp
                    !
                    IF (ltetra) THEN
                      lselect = (ABS(wdk(ibnd, ik)) > eps32)
                      ! lselect = .true.
                    ELSE
                      deltaf = &
                        w0gauss((et(ibnd, ik) - efermi(ispin)) / &
                        degauss, ngauss) / degauss
                        ! w0gauss((et(ibnd, ik) - efermi) / degauss, ngauss)
                      lselect = (ABS(deltaf) > eps32) .OR. (.NOT. lspeedup)
                    END IF
                    !
                    IF (lselect) THEN
                      !
                      cnr_aux = CMPLX(0._dp, 0._dp, KIND=DP)
                      !
                      ! integral over r angle (Gauss points)
                      !
                      DO igp = 1, ngp
                        !
                        cnr_aux = cnr_aux + &
                          psi_krtau_aux(igp, ibnd, ik) * &
                          CONJG(ylm(l * (l + 1) + 1 + m, igp)) * &
                          gp_wt(igp)
                        !
                      END DO ! igp
                      !
                      IF (ltetra) THEN
                        dos_nlmr(ir - imin + 1, l0 + m, ispin, iat) = &
                          dos_nlmr(ir - imin + 1, l0 + m, ispin, iat) + &
                          prefactor * &
                          REAL(cnr_aux * CONJG(cnr_aux), KIND=DP) * &
                          wdk(ibnd, ik)
                        !
                        dos_nlmr_nodloglde(ir - imin + 1, l0 + m, ispin, iat) &
                          = &
                          dos_nlmr_nodloglde(ir - imin + 1, l0 + m, &
                            ispin, iat) + &
                          prefactor_nodloglde * &
                          REAL(cnr_aux * CONJG(cnr_aux), KIND=DP) * &
                          wdk(ibnd, ik)
                      ELSE
                        dos_nlmr(ir - imin + 1, l0 + m, ispin, iat) = &
                          dos_nlmr(ir - imin + 1, l0 + m, ispin, iat) + &
                          prefactor * &
                          REAL(cnr_aux * CONJG(cnr_aux), KIND=DP) * &
                          wk(ik) * deltaf
                        !
                        dos_nlmr_nodloglde(ir - imin + 1, l0 + m, ispin, iat) &
                          = &
                          dos_nlmr_nodloglde(ir - imin + 1, l0 + m, &
                            ispin, iat) + &
                          prefactor_nodloglde * &
                          REAL(cnr_aux * CONJG(cnr_aux), KIND=DP) * &
                          wk(ik) * deltaf
                      END IF
                      !
                    END IF ! lselect
                    !
                  END DO ! ibnd
                  !
                  END IF ! isk
                  !
                END DO ! ik
                !
                dos_nlr(ir - imin + 1, iorb, ispin, iat) = &
                  dos_nlr(ir - imin + 1, iorb, ispin, iat) + &
                  dos_nlmr(ir - imin + 1, l0 + m, ispin, iat)
                !
                dos_nlr_nodloglde(ir - imin + 1, iorb, ispin, iat) = &
                  dos_nlr_nodloglde(ir - imin + 1, iorb, ispin, iat) + &
                  dos_nlmr_nodloglde(ir - imin + 1, l0 + m, ispin, iat)
                !
                dos_nr(ir - imin + 1, ispin, iat) = &
                  dos_nr(ir - imin + 1, ispin, iat) + &
                  dos_nlmr(ir - imin + 1, l0 + m, ispin, iat)
                !
                !
              END DO ! im
              !
            END DO ! iorb
            !
          END DO ! ir
        END DO ! ispin
      END DO ! iat
      !
      !
      IF (lsymmetrize) THEN
        !
        WRITE(stdout, &
          '(7x, "Symmetrizing partial DOS = f(r, E_F)...")')
        !
        !
        ALLOCATE(dos_nlmr_sym(imax - imin + 1, (lmax + 1) * (lmax + 1), &
          nspin, nat), STAT = ierr)
        IF (ierr /= 0) &
          CALL errore(routine_name, 'Error allocating dos_nlmr_sym', 1)
        dos_nlmr_sym(:, :, :, :) = dos_nlmr(:, :, :, :)
        !
        ALLOCATE(dos_nlr_sym(imax - imin + 1, norb, nspin, nat), STAT = ierr)
        IF (ierr /= 0) &
          CALL errore(routine_name, 'Error allocating dos_nlr_sym', 1)
        dos_nlr_sym(:, :, :, :) = dos_nlr(:, :, :, :)
        !
        ALLOCATE(dos_nlmr_nodloglde_sym(imax - imin + 1, &
          (lmax + 1) * (lmax + 1), &
          nspin, nat), STAT = ierr)
        IF (ierr /= 0) &
          CALL errore(routine_name, &
            'Error allocating dos_nlmr_nodloglde_sym', 1)
        dos_nlmr_nodloglde_sym(:, :, :, :) = dos_nlmr_nodloglde(:, :, :, :)
        !
        ALLOCATE(dos_nlr_nodloglde_sym(imax - imin + 1, norb, nspin, nat), &
          STAT = ierr)
        IF (ierr /= 0) &
          CALL errore(routine_name, 'Error allocating dos_nlr_nodloglde_sym', 1)
        dos_nlr_nodloglde_sym(:, :, :, :) = dos_nlr_nodloglde(:, :, :, :)
        !
        ALLOCATE(counters(nat), STAT = ierr)
        IF (ierr /= 0) &
          CALL errore(routine_name, 'Error allocating counters', 1)
        counters(:) = 1
        !
        !
        DO iat = 1, nat
          !
          !
          DO jat = 1, nat
            !
            lfound = .FALSE.
            !
            irot_loop: DO irot = 1, nrot
              IF ((irt(irot, iat) == jat) .AND. (iat /= jat)) THEN
                lfound = .TRUE.
                counters(iat) = counters(iat) + 1
                EXIT irot_loop
              END IF
            END DO irot_loop ! irot
            !
            IF (lfound) THEN
              !
              DO ispin = 1, nspin
                DO ir = imin, imax
                  !
                  DO iorb = 1, norb
                    !
                    l = iorb - 1
                    l0 = l * (l + 1) + 1
                    !
                    DO im = 1, 2 * l + 1
                      !
                      m = im - l - 1
                      !
                      dos_nlmr_sym(ir - imin + 1, l0 + m, ispin, iat) = &
                        dos_nlmr_sym(ir - imin + 1, l0 + m, ispin, iat) + &
                        dos_nlmr(ir - imin + 1, l0 + m, ispin, jat)
                      !
                      dos_nlmr_nodloglde_sym(ir - imin + 1, l0 + m, &
                        ispin, iat) = &
                        dos_nlmr_nodloglde_sym(ir - imin + 1, l0 + m, &
                          ispin, iat) + &
                        dos_nlmr_nodloglde(ir - imin + 1, l0 + m, ispin, jat)
                      !
                    END DO ! im
                    !
                    dos_nlr_sym(ir - imin + 1, iorb, ispin, iat) = &
                      dos_nlr_sym(ir - imin + 1, iorb, ispin, iat) + &
                      dos_nlr(ir - imin + 1, iorb, ispin, jat)
                    !
                    dos_nlr_nodloglde_sym(ir - imin + 1, iorb, ispin, iat) = &
                      dos_nlr_nodloglde_sym(ir - imin + 1, iorb, ispin, iat) + &
                      dos_nlr_nodloglde(ir - imin + 1, iorb, ispin, jat)
                    !
                  END DO ! iorb
                  !
                END DO ! ir
              END DO ! ispin
              !
            END IF ! lfound
            !
          END DO ! jat
          !
        END DO ! iat
        !
        !
        DO iat = 1, nat
          !
          ! updating partial DOS, inside MT sphere(s)
          !
          dos_nlr(:, :, :, iat) = &
            dos_nlr_sym(:, :, :, iat) / counters(iat)
          dos_nlmr(:, :, :, iat) = &
            dos_nlmr_sym(:, :, :, iat) / counters(iat)
          !
          dos_nlr_nodloglde(:, :, :, iat) = &
            dos_nlr_nodloglde_sym(:, :, :, iat) / counters(iat)
          dos_nlmr_nodloglde(:, :, :, iat) = &
            dos_nlmr_nodloglde_sym(:, :, :, iat) / counters(iat)
          !
          ! updating total DOS, inside MT sphere(s)
          !
          dos_nr(:, :, iat) = 0._dp
          !
          DO ispin = 1, nspin
            DO ir = imin, imax
              DO iorb = 1, norb
                !
                dos_nr(ir - imin + 1, ispin, iat) = &
                  dos_nr(ir - imin + 1, ispin, iat) + &
                  dos_nlr(ir - imin + 1, iorb, ispin, iat)
                !
              END DO ! iorb
            END DO ! ir
          END DO ! ispin
          !
        END DO ! iat
        !
        !
        DEALLOCATE(counters, STAT = ierr)
        IF (ierr /= 0) &
          CALL errore(routine_name, 'Error deallocating counters', 1)
        !
        DEALLOCATE(dos_nlr_sym, STAT = ierr)
        IF (ierr /= 0) &
          CALL errore(routine_name, 'Error deallocating dos_nlr_sym', 1)
        !
        DEALLOCATE(dos_nlmr_sym, STAT = ierr)
        IF (ierr /= 0) &
          CALL errore(routine_name, 'Error deallocating dos_nlmr_sym', 1)
        !
        DEALLOCATE(dos_nlr_nodloglde_sym, STAT = ierr)
        IF (ierr /= 0) &
          CALL errore(routine_name, &
            'Error deallocating dos_nlr_nodloglde_sym', 1)
        !
        DEALLOCATE(dos_nlmr_nodloglde_sym, STAT = ierr)
        IF (ierr /= 0) &
          CALL errore(routine_name, &
            'Error deallocating dos_nlmr_nodloglde_sym', 1)
        !
        WRITE(stdout, &
          '(7x, "Done symmetrizing partial DOS = f(r, E_F).", /6x, /6x)')
        !
      END IF ! lsymmetrize
      !

      !
      !
      WRITE(stdout, '(6x, "Done computing partial DOS per atom per spin", &
        & " = f(r, E_F).", /6x, /6x)')
      !
      !
      ! clean-up
      !
      DEALLOCATE(psi_krtau_aux, STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, &
        'Error deallocating psi_krtau_aux', 1)
      !
      DEALLOCATE(ylm, STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, 'Error deallocating ylm', 1)
      !
      DEALLOCATE(psi_kg, STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, 'Error deallocating psi_kg', 1)
      !
      DEALLOCATE(gp_vec, STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, 'Error deallocating gp_vec', 1)
      !
      DEALLOCATE(gp_wt, STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, 'Error deallocating gp_wt', 1)
      !
      !
      WRITE(stdout, &
        '(/6x, "Computing total DOS per atom per spin = f(spin, E_F)...")')
      !
      dos_n(:) = 0._dp
      !
      DO ispin = 1, nspin
          !
          WRITE(stdout, '(/6x, "spin #", I0, " out of ", I0)') ispin, nspin
          !
          DO ik = 1, nkstot
            !
            IF (isk(ik) == ispin) THEN
              !
              DO ibnd = 1, nbnd
                !
                IF (ltetra) THEN
                  lselect = (ABS(wdk(ibnd, ik)) > eps32)
                  ! lselect = .true.
                ELSE
                  deltaf = &
                    w0gauss((et(ibnd, ik) - efermi(ispin)) / &
                    degauss, ngauss) / degauss
                    lselect = (ABS(deltaf) > eps32) .OR. (.NOT. lspeedup)
                END IF
                !
                IF (lselect) THEN
                  IF (ltetra) THEN
                    dos_n(ispin) = dos_n(ispin) + &
                      prefactor_tot_dos * wdk(ibnd, ik)
                  ELSE
                    dos_n(ispin) = dos_n(ispin) + &
                      prefactor_tot_dos * deltaf * wk(ik)
                  END IF
                END IF ! lselect
                !
              END DO ! ibnd
              !
            END IF ! isk
            !
          END DO ! ik
          !
          WRITE(stdout, &
            '(/7x, "n(E_F, ", I1, ") = ", ES16.8)') ispin, dos_n(ispin)
          !
      END DO ! ispin
      !
      ! averaging of total DOS per spin
      !
      IF (nspin == 2) THEN
        avg_dos_per_spin = 0.5_dp * (dos_n(1) + dos_n(2))
        dos_n(1) = avg_dos_per_spin
        dos_n(2) = avg_dos_per_spin
      END IF
      !
      WRITE(stdout, &
        '(6x, "Done computing total DOS", &
        & " per atom per spin = f(spin, E_F).", /6x, /6x)')
      !
      WRITE(stdout, '(/5x, ">>>>>>>>>>    PartDOS END    <<<<<<<<<<<", &
        & /5x, /5x, /5x)')
      !
      IF (ltetra) THEN
        !
        DEALLOCATE(wdk, STAT = ierr)
        IF (ierr /= 0) CALL errore(routine_name, 'Error deallocating wdk', 1)
        !
      END IF
      !
      CALL stop_clock(routine_name)
      !
      !
    !---------------------------------------------------------------------------
    END SUBROUTINE set_dos_n
    !---------------------------------------------------------------------------
    !
    !
    !---------------------------------------------------------------------------
    PURE REAL(DP) FUNCTION get_sum_wk()
    !---------------------------------------------------------------------------
    !!
    !! Checks if weights wk sum up to 1.
    !!
    !---------------------------------------------------------------------------
      USE klist, ONLY: nkstot, wk
      !
      IMPLICIT NONE
      !
      CHARACTER(len=256) :: routine_name
      !! name of this subroutine
      INTEGER :: ik
      !! iterators
      REAL(DP) :: sum_wk
      !! sum of all weights wk for all k-points
      !
      routine_name = "get_sum_wk"
      !
      sum_wk = 0._dp
      !
      DO ik = 1, nkstot
        sum_wk = sum_wk + wk(ik)
      END DO ! ik
      !
      get_sum_wk = sum_wk
      !
    !---------------------------------------------------------------------------
    END FUNCTION get_sum_wk
    !---------------------------------------------------------------------------
    !
    !
    !---------------------------------------------------------------------------
    PURE REAL(DP) FUNCTION get_sum_wdk(wdk)
    !---------------------------------------------------------------------------
    !!
    !! Checks if weights wdk sum up to 1
    !!
    !---------------------------------------------------------------------------
      USE klist, ONLY: nkstot
      USE wvfct, ONLY: nbnd
      !
      IMPLICIT NONE
      !
      REAL(DP), INTENT(in) :: wdk(:, :)
      !! weights wdk for all bands and k-points
      !
      ! local variables
      !
      CHARACTER(len=256) :: routine_name
      !! name of this subroutine
      INTEGER :: ik, ibnd
      !! iterators
      REAL(DP) :: sum_wdk
      !! sum of all weights wk for all k-points
      !
      routine_name = "get_sum_wdk"
      !
      sum_wdk = 0._dp
      !
      DO ik = 1, nkstot
        DO ibnd = 1, nbnd
        sum_wdk = sum_wdk + wdk(ibnd, ik)
        END DO
      END DO ! ik
      !
      get_sum_wdk = sum_wdk
      !
    !---------------------------------------------------------------------------
    END FUNCTION get_sum_wdk
    !---------------------------------------------------------------------------
    !
    !
    !---------------------------------------------------------------------------
    SUBROUTINE tetra_delta_weights(nks, nspin, is, isk, nbnd, ef, et, wdk)
    !---------------------------------------------------------------------------
    !!
    !! Tetrahedron weights for integration with the delta-function,
    !! as Fermi-energy derivatives of the weights published by
    !! Blochl et al. in 10.1103/PhysRevB.49.16223
    !!
    !! Not to be confused with regular tetrahedron weights for
    !! integration with the occupation function
    !!
    !! Adapted from tetra_weights_only of the PW package
    !!
    !---------------------------------------------------------------------------
    !
    ! D. Radevych
    !
      USE mp_images, ONLY: intra_image_comm
      USE mp, ONLY: mp_sum
      USE ktetra, ONLY: ntetra, tetra
      USE constants, ONLY: eps32
      !
      IMPLICIT NONE
      !
      EXTERNAL :: errore, divide, hpsort, start_clock, stop_clock
      !
      INTEGER, INTENT(IN) :: nks
      !! number of k in irreducible BZ
      INTEGER, INTENT(IN) :: nspin
      !! number of spins
      INTEGER, INTENT(IN) :: is
      !! spin label
      INTEGER, INTENT(IN) :: isk(nks)
      !! for each k-point: 1 = spin up, 2 = spin down
      INTEGER, INTENT(IN) :: nbnd
      !! number of bands
      REAL(DP), INTENT(IN) :: ef
      !! Fermi energy
      REAL(DP), INTENT(IN) :: et(nbnd, nks)
      !! eigenvalues of the hamiltonian
      REAL(DP), INTENT(INOUT) :: wdk(nbnd, nks)
      !! the weight of each k point and band
      ! wg must be (inout) and not (out) because if is /= 0 only terms for
      ! spin = is are initialized; the remaining terms should be kept, not lost
      !
      ! local variables
      !
      CHARACTER(len = 200) :: routine_name
      !! name of this subroutine
      INTEGER :: ik, ibnd, nt, nk, ns, i
      !! iterators
      INTEGER :: kp1, kp2, kp3, kp4
      !! irreductible k-points corresponding to the given tetrahedron
      INTEGER :: itetra(4)
      !! tetrahedron indices
      INTEGER :: nspin_lsda
      !! number of spinds
      INTEGER :: s_tetra, l_tetra
      !!
      REAL(DP) :: rfac
      !! real prefactor
      REAL(DP) :: e1, e2, e3, e4
      !! energies in tetrahedron's corners
      REAL(DP) :: etetra(4)
      !! tetrahedron energies
      REAL(DP) :: e21, e31, e41, e32, e42, e43
      !! energy differences
      REAL(DP) :: C, C1, C2, C3, DC, DC1, DC2, DC3, E
      !! tmp, aux
      !
      E = ef
      !
      routine_name = "tetra_delta_weights"
      CALL start_clock(routine_name)
      !
      rfac = 1._dp / ntetra
      !
      ! check if tetrahedrons were initialized
      IF (ntetra == 0) &
        CALL errore(routine_name, 'called before tetrahedron initialization', 1)
      !
      !
      ! set weights for spin "is" to zero
      DO ik = 1, nks
        !
        IF (is /= 0) THEN
          IF (isk(ik) /= is) CYCLE
        ENDIF
        !
        DO ibnd = 1, nbnd
          wdk(ibnd, ik) = 0._dp
        ENDDO ! ibnd
        !
      ENDDO ! ik
      !
      ! set local number of spins (DR: nspin could be used instead?)
      IF (nspin == 2) THEN
        nspin_lsda = 2
      ELSE
        nspin_lsda = 1
      ENDIF
      !
      CALL divide(intra_image_comm, ntetra, s_tetra, l_tetra)
      !
      !$OMP PARALLEL DEFAULT(NONE) &
      !$OMP & SHARED(nspin_lsda, is, nks, s_tetra, l_tetra, nbnd, et, tetra, &
      !$OMP & wdk, ef, ntetra) &
      !$OMP & PRIVATE(ns, nk, nt, ibnd, i, etetra, itetra, e1, e2, e3, e4, &
      !$OMP & kp1, kp2, kp3, kp4)
      !
      DO ns = 1, nspin_lsda
        !
        ! skip if not desired spin
        IF (is /= 0) THEN
          IF (ns /= is) CYCLE
        ENDIF
        !
        ! nk is used to select k-points with up (ns = 1) or down (ns = 2) spin
        !
        IF (ns == 1) THEN
          nk = 0
        ELSE
          nk = nks / 2
        ENDIF
        !
        DO nt = s_tetra, l_tetra
          !
          !$OMP DO
          !
          loop_ibnd: DO ibnd = 1, nbnd
            !
            ! etetra are the energies at the verteces of the nt-th tetrahedron
            !
            DO i = 1, 4
              etetra(i) = et(ibnd, tetra(i, nt) + nk)
            ENDDO
            itetra(1) = 0
            !
            ! sort in ascending order: e1 < e2 < e3 < e4
            !
            CALL hpsort(4, etetra, itetra)
            !
            e1 = etetra(1)
            e2 = etetra(2)
            e3 = etetra(3)
            e4 = etetra(4)
            !
            ! energy differences
            !
            e21 = e2 - e1
            e31 = e3 - e1
            e41 = e4 - e1
            e32 = e3 - e2
            e42 = e4 - e2
            e43 = e4 - e3
            !
            ! if one of the energy differences is zero, cycle loop_ibnd
            ! IF (ABS(e21) < eps32 .OR. ABS(e31) < eps32 .OR. &
            !   ABS(e41) < eps32 .OR. ABS(e32) < eps32 .OR. &
            !   ABS(e42) < eps32 .OR. ABS(e43) < eps32) THEN
            !   CYCLE loop_ibnd
            ! END IF
            !
            ! kp1 - kp4 are the irreducible k-points corresponding to e1 - e4
            !
            kp1 = tetra(itetra(1), nt) + nk
            kp2 = tetra(itetra(2), nt) + nk
            kp3 = tetra(itetra(3), nt) + nk
            kp4 = tetra(itetra(4), nt) + nk
            !
            ! calculate weights wg
            !
            ! IF (e1 > ef .OR. ef >= e4) THEN
            !   !
            !   ! redundancy
            !   !
            !   wdk(ibnd, kp1) = 0._dp
            !   wdk(ibnd, kp2) = 0._dp
            !   wdk(ibnd, kp3) = 0._dp
            !   wdk(ibnd, kp4) = 0._dp
            !   !
            ! ELSE
            IF (e1 <= ef .AND. ef < e2) THEN
              !
              ! wdk(ibnd, kp1) = wdk(ibnd, kp1) + rfac * (ef - e1)**2 / &
              !   (e21 * e31 * e41) * &
              !   (-(ef - e1) * (1._dp / e21 + 1._dp / e31 + 1._dp / e41) + &
              !   3._dp)
              ! !
              ! wdk(ibnd, kp2) = wdk(ibnd, kp2) + rfac * (ef - e1)**3 / &
              !   (e21 * e31 * e41) / &
              !   e21
              ! !
              ! wdk(ibnd, kp3) = wdk(ibnd, kp3) + rfac * (ef - e1)**3 / &
              !   (e21 * e31 * e41) / &
              !   e31
              ! !
              ! wdk(ibnd, kp4) = wdk(ibnd, kp4) + rfac * (ef - e1)**3 / &
              !   (e21 * e31 * e41) / &
              !   e41
              !
              C=rfac / 4.* (E-e1)**3/(e2-e1)/(e3-e1)/(e4-e1)
              DC=rfac / 4.* 3*(E-e1)**2/(e2-e1)/(e3-e1)/(e4-e1)
              wdk(ibnd, kp1) = wdk(ibnd, kp1) + &
                DC*(4-(E-e1)*(1/(e2-e1)+1/(e3-e1)+1/(e4-e1))) &
                -C*(1/(e2-e1)+1/(e3-e1)+1/(e4-e1))
              !
              wdk(ibnd, kp2) = wdk(ibnd, kp2) + &
                DC*(E-e1)/(e2-e1)+C/(e2-e1)
              !
              wdk(ibnd, kp3) = wdk(ibnd, kp3) + &
                DC*(E-e1)/(e3-e1)+C/(e3-e1)
              !
              wdk(ibnd, kp4) = wdk(ibnd, kp4) + &
                DC*(E-e1)/(e4-e1)+C/(e4-e1)
              !
            ELSEIF (e2 <= ef .AND. ef < e3) THEN
              !
              ! wdk(ibnd, kp1) = wdk(ibnd, kp1) + rfac * &
              !   (ef**3 * e1**2 - 2._dp * ef**3 * e1 * e3 - &
              !   2._dp * ef**3 * e1 * e4 + 2._dp * e2 * ef**3 * e1 + &
              !   ef**3 * e3**2 + ef**3 * e3 * e4 - e2 * ef**3 * e3 + &
              !   ef**3 * e4**2 - e2 * ef**3 * e4 - 3._dp * e2 * ef**2 * e1**2 + &
              !   6._dp * ef**2 * e1 * e3 * e4 - 3._dp * ef**2 * e3**2 * e4 - &
              !   3._dp * ef**2 * e3 * e4**2 + 3._dp * e2 * ef**2 * e3 * e4 - &
              !   3._dp * ef * e1**2 * e3 * e4 + 3._dp * e2 * ef * e1**2 * e3 + &
              !   3._dp * e2 * ef * e1**2 * e4 - &
              !   6._dp * e2 * ef * e1 * e3 * e4 + &
              !   3._dp * ef* e3**2 * e4**2 + e1**2 * e3**2 * e4 - &
              !   e2 * e1**2 * e3**2 + e1**2 * e3 * e4**2 - &
              !   e2 * e1**2 * e3 * e4 - e2 * e1**2 * e4**2 - &
              !   2._dp * e1 * e3**2 * e4**2 + 2._dp * e2 * e1 * e3**2 * e4 + &
              !   2._dp * e2 * e1 * e3 * e4**2 - e2 * e3**2 * e4**2) / &
              !   (e31 * e32 * e41 * e42) / (e31 * e41)
              ! !
              ! wdk(ibnd, kp2) = wdk(ibnd, kp2) + rfac * &
              !   (ef**3 * e2**2 - 2._dp * ef**3 * e2 * e3 - &
              !   2._dp * ef**3 * e2 * e4 + 2._dp * e1 * ef**3 * e2 + &
              !   ef**3 * e3**2 + ef**3 * e3 * e4 - e1 * ef**3 * e3 + &
              !   ef**3 * e4**2 - e1 * ef**3 * e4 - 3._dp * e1 * ef**2 * e2**2 + &
              !   6._dp * ef**2 * e2 * e3 * e4 - 3._dp * ef**2 * e3**2 * e4 - &
              !   3._dp * ef**2 * e3 * e4**2 + 3._dp * e1 * ef**2 * e3 * e4 - &
              !   3._dp * ef * e2**2 * e3 * e4 + 3._dp * e1 * ef * e2**2 * e3 + &
              !   3._dp * e1 * ef * e2**2 * e4 - &
              !   6._dp * e1 * ef * e2 * e3 * e4 + &
              !   3._dp * ef * e3**2 * e4**2 + e2**2 * e3**2 * e4 - &
              !   e1 * e2**2 * e3**2 + e2**2 * e3 * e4**2 - &
              !   e1 * e2**2 * e3 * e4 - e1 * e2**2 * e4**2 - &
              !   2._dp * e2 * e3**2 * e4**2 + 2._dp * e1 * e2 * e3**2 * e4 + &
              !   2._dp * e1 * e2 * e3 * e4**2 - e1 * e3**2 * e4**2) / &
              !   (e31 * e32 * e41 * e42) / (e32 * e42)
              ! !
              ! wdk(ibnd, kp3) = wdk(ibnd, kp3) - rfac * &
              !   (ef**3 * e1**2 + ef**3 * e1 * e2 - 2._dp * ef**3 * e1 * e3 - &
              !   e4 * ef**3 * e1 + ef**3 * e2**2 - 2._dp * ef**3 * e2 * e3 - &
              !   e4 * ef**3 * e2 + ef**3 * e3**2 + 2._dp * e4 * ef**3 * e3 - &
              !   3._dp * ef**2 * e1**2 * e2 - 3._dp * ef**2 * e1 * e2**2 + &
              !   6._dp * ef**2 * e1 * e2 * e3 + 3._dp * e4 * ef**2 * e1 * e2 - &
              !   3._dp * e4 * ef**2 * e3**2 + 3._dp * ef * e1**2 * e2**2 - &
              !   3._dp * ef * e1 * e2 * e3**2 - &
              !   6._dp * e4 * ef * e1 * e2 * e3 + &
              !   3._dp * e4 * ef * e1 * e3**2 + 3._dp * e4 * ef * e2 * e3**2 - &
              !   2._dp * e1**2 * e2**2 * e3 - e4 * e1**2 * e2**2 + &
              !   e1**2 * e2 * e3**2 + 2._dp * e4 * e1**2 * e2 * e3 - &
              !   e4 * e1**2 * e3**2 + e1 * e2**2 * e3**2 + &
              !   2._dp * e4 * e1 * e2**2 * e3 - e4 * e1 * e2 * e3**2 - &
              !   e4 * e2**2 * e3**2) / &
              !   (e31 * e32 * e41 * e42) / (e31 * e32)
              ! !
              ! wdk(ibnd, kp4) = wdk(ibnd, kp4) - rfac * &
              !   (ef**3 * e1**2 + ef**3 * e1 * e2 - 2._dp * ef**3 * e1 * e4 - &
              !   e3 * ef**3 * e1 + ef**3 * e2**2 - 2._dp * ef**3 * e2 * e4 - &
              !   e3 * ef**3 * e2 + ef**3 * e4**2 + 2._dp * e3 * ef**3 * e4 - &
              !   3._dp * ef**2 * e1**2 * e2 - 3._dp * ef**2 * e1 * e2**2 + &
              !   6._dp * ef**2 * e1 * e2 * e4 + 3._dp * e3 * ef**2 * e1 * e2 - &
              !   3._dp * e3 * ef**2 * e4**2 + 3._dp * ef * e1**2 * e2**2 - &
              !   3._dp * ef * e1 * e2 * e4**2 - &
              !   6._dp * e3 * ef * e1 * e2 * e4 + &
              !   3._dp * e3 * ef * e1 * e4**2 + 3._dp * e3 * ef * e2 * e4**2 - &
              !   2._dp * e1**2 * e2**2 * e4 - e3 * e1**2 * e2**2 + &
              !   e1**2 * e2 * e4**2 + 2._dp * e3 * e1**2 * e2 * e4 - &
              !   e3 * e1**2 * e4**2 + e1 * e2**2 * e4**2 + &
              !   2._dp * e3 * e1 * e2**2 * e4 - e3 * e1 * e2 * e4**2 - &
              !   e3 * e2**2 * e4**2) / &
              !   (e31 * e32 * e41 * e42) / (e41 * e42)
              !
              C1=rfac / 4.* (E-e1)**2/(e4-e1)/(e3-e1)
              DC1=rfac / 4. * 2*(E-e1)/(e4-e1)/(e3-e1)
              C2=rfac / 4. * (E-e1)*(E-e2)*(e3-E)/(e4-e1)/(e3-e2)/(e3-e1)
              DC2=rfac / 4. * ((E-e2)*(e3-E)+(E-e1)*(e3-E)-(E-e1)*(E-e2)) &
                 /(e4-e1)/(e3-e2)/(e3-e1)
              C3=rfac / 4. * (E-e2)**2*(e4-E)/(e4-e2)/(e3-e2)/(e4-e1)
              DC3=rfac / 4. * (2*(E-e2)*(e4-E)-(E-e2)**2) &
              /(e4-e2)/(e3-e2)/(e4-e1)
              !
              wdk(ibnd, kp1) = wdk(ibnd, kp1) + &
                DC1+(DC1+DC2)*(e3-E)/(e3-e1)+(DC1+DC2+DC3)*(e4-E)/(e4-e1) &
                -(C1+C2)/(e3-e1)-(C1+C2+C3)/(e4-e1)
              !
              wdk(ibnd, kp2) = wdk(ibnd, kp2) + &
                DC1+DC2+DC3+(DC2+DC3)*(e3-E)/(e3-e2)+DC3*(e4-E)/(e4-e2) &
                       -(C2+C3)/(e3-e2)-C3/(e4-e2)
              !
              wdk(ibnd, kp3) = wdk(ibnd, kp3) + &
                (DC1+DC2)*(E-e1)/(e3-e1)+(DC2+DC3)*(E-e2)/(e3-e2) &
                +(C1+C2)/(e3-e1)+(C2+C3)/(e3-e2)
              !
              wdk(ibnd, kp4) = wdk(ibnd, kp4) + &
                (DC1+DC2+DC3)*(E-e1)/(e4-e1)+DC3*(E-e2)/(e4-e2) &
                +(C1+C2+C3)/(e4-e1)+C3/(e4-e2)
              !
            ELSEIF (e3 <= ef .AND. ef < e4) THEN
              !
              ! wdk(ibnd, kp1) = wdk(ibnd, kp1) - rfac * (ef - e4)**3 / &
              !   (e41 * e42 * e43) / e41
              ! !
              ! wdk(ibnd, kp2) = wdk(ibnd, kp2) - rfac * (ef - e4)**3 / &
              !   (e41 * e42 * e43) / e42
              ! !
              ! wdk(ibnd, kp3) = wdk(ibnd, kp3) - rfac * (ef - e4)**3 / &
              !   (e41 * e42 * e43) / e43
              ! !
              ! wdk(ibnd, kp4) = wdk(ibnd, kp4) + rfac * (ef - e4)**2 / &
              !   (e41 * e42 * e43) * &
              !   ((ef - e4) * (1._dp / e41 + 1._dp / e42 + 1._dp / e43) + 3._dp)
              !
              C = rfac / 4. * (e4-E)**3 / (e4-e1) / (e4-e2) / (e4-e3)
              DC = rfac / 4. * -3*(e4-E)**2 / (e4-e1) / (e4-e2) / (e4-e3)
              !
              wdk(ibnd, kp1) = wdk(ibnd, kp1) &
                -DC*(e4-E)/(e4-e1)+C/(e4-e1)
              !
              wdk(ibnd, kp2) = wdk(ibnd, kp2) &
                -DC*(e4-E)/(e4-e2)+C/(e4-e2)
              !
              wdk(ibnd, kp3) = wdk(ibnd, kp3) &
                -DC*(e4-E)/(e4-e3)+C/(e4-e3)
              !
              wdk(ibnd, kp4) = wdk(ibnd, kp4) &
                -DC*(4-(e4-E)*(1/(e4-e1)+1/(e4-e2)+1/(e4-e3))) &
                -C*(1/(e4-e1)+1/(e4-e2)+1/(e4-e3))
              !
            ENDIF
            !
          END DO loop_ibnd ! ibnd
          !
          !$OMP END DO NOWAIT
          !
        END DO ! nt
        !
      END DO ! ns
      !
      !$OMP END PARALLEL
      !
      ! Each process in a "image" communicator has contributions from
      ! only a part of tetrahedra (s_tetra - l_tetra).
      ! Contribution from full BZ is computed by combining them.
      !
      CALL mp_sum(wdk, intra_image_comm)
      !
      ! add correct spin normalization (2 for LDA, 1 for all other cases)
      IF (nspin == 1) THEN
        wdk(1 : nbnd, 1 : nks) = wdk(1 : nbnd, 1 : nks) * 2._dp
      END IF
      !
      CALL stop_clock(routine_name)
      !
    !---------------------------------------------------------------------------
    END SUBROUTINE tetra_delta_weights
    !---------------------------------------------------------------------------
    !
    !
  !=============================================================================
  END MODULE part_dos
  !=============================================================================
