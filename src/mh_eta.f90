  ! Copyright (C) 2024-2026 Danylo Radevych
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  ! Please cite: DOI: https://doi.org/10.1038/s41524-026-02141-7
  !
  !
  !=============================================================================
  MODULE mh_eta
  !=============================================================================
  !!
  !! Module for computation of l-resolved McMillan-Hopfield
  !! paremeters eta_l.
  !!
  !
  !  D. Radevych
  !
    USE kinds, ONLY : DP
    !
    !
    IMPLICIT NONE
    !
    PUBLIC :: set_eta
    !
  !
  CONTAINS
  !
    !
    !---------------------------------------------------------------------------
    SUBROUTINE set_eta(nr, imin, imax, nat, norb, nspin, &
      dos_nlr, dos_nr, dos_n, luse_tot_dos, mll1, etall1)
    !---------------------------------------------------------------------------
    !!
    !! Computes McMillan-Hopfield eta_l
    !!
      USE constants, ONLY: eps12
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(in) :: nr
      !! number of points on radial grid
      INTEGER, INTENT(in) :: imin
      !! min index of the radial point for eta
      !! eta for ir < nmin will be left zero
      INTEGER, INTENT(in) :: imax
      !! max index of the radial point for eta
      !! eta for ir > nmax will be left zero
      INTEGER, INTENT(in) :: nat
      !! number of atoms
      INTEGER, INTENT(in) :: norb
      !! number of orbitals
      INTEGER, INTENT(in) :: nspin
      !! number of spins
      REAL(DP), INTENT(in) :: dos_nlr(:, :, :, :)
      !! partial densities n^i_{l}(r, E_F)
      !! dos_nlr(nr, norb, nspin, nat)
      !! lmax = norb - 1
      REAL(DP), INTENT(in) :: dos_nr(:, :, :)
      !! partial densities n^i(r, E_F), per atom, per spin
      !! dos_nr(nr, nspin, nat)
      !! lmax = norb - 1
      REAL(DP), INTENT(in) :: dos_n(:)
      !! total densities n(E_F), per spin
      !! dos_nr(nspin)
      LOGICAL, INTENT(in) :: luse_tot_dos
      !! if true, use total DOS per spin,
      !! not the one inside each atom
      REAL(DP), INTENT(in) :: mll1(:, :, :, :)
      !! M_{l, l+1}(r, E_F)
      !! mll1(nr, norb, nspin, nat)
      REAL(DP), INTENT(inout) :: etall1(:, :, :, :)
      !! \eta_{l}(r, E_F)
      !! etall1(nr, norb, nspin, nat)
      !! note etall1(:, norb, :, :) = etar (total eta)
      !
      CHARACTER(len=256) :: routine_name
      !! name of this subroutine
      INTEGER :: iat, ispin, iorb, ir
      !! iterators
      INTEGER :: l
      !! angular momentum l
      REAL(DP) :: rtmp_mll1, rtmp_n, rtmp_nl, rtmp_nl1
      !! temporary reals
      !
      EXTERNAL :: errore, start_clock, stop_clock
      !
      routine_name = "set_eta"
      !
      CALL start_clock(routine_name)
      !
      IF (imax > nr) &
        CALL errore(routine_name, "imax > nr", 1)
      !
      etall1(:, :, :, :) = 0._dp
      !
      DO iat = 1, nat
        DO ispin = 1, nspin
          DO iorb = 1, norb - 1
            !
            l = iorb - 1
            !
            DO ir = imin, imax
              !
              IF (.NOT. luse_tot_dos) THEN
                rtmp_n = dos_nr(ir - imin + 1, ispin, iat)
              ELSE
                rtmp_n = dos_n(ispin)
              END IF
              !
              IF (rtmp_n < 0._dp) &
                CALL errore(routine_name, "n < 0", 1)
              !
              rtmp_nl = dos_nlr(ir - imin + 1, iorb, ispin, iat)
              IF (rtmp_nl < 0._dp) &
                CALL errore(routine_name, "n_l < 0", 1)
              !
              rtmp_nl1 = dos_nlr(ir - imin + 1, iorb + 1, ispin, iat)
              IF (rtmp_nl1 < 0._dp) &
                CALL errore(routine_name, "n_{l+1} < 0", 1)
              !
              rtmp_mll1 = mll1(ir, iorb, ispin, iat)
              !
              IF (ABS(rtmp_n) > eps12) THEN
                !
                ! \eta_{l, l+1}
                !
                etall1(ir - imin + 1, iorb, ispin, iat) = &
                  2._dp * (l + 1._dp) / &
                  (2._dp * l + 1._dp) / (2._dp * l + 3._dp) * &
                  rtmp_mll1 * rtmp_mll1 * &
                  (rtmp_nl / SQRT(rtmp_n)) * &
                  (rtmp_nl1 / SQRT(rtmp_n))
                !
                ! \eta = \sum_l \eta_{l, l+1}
                !
                etall1(ir - imin + 1, norb, ispin, iat) = &
                  etall1(ir - imin + 1, norb, ispin, iat) + &
                  etall1(ir - imin + 1, iorb, ispin, iat)
                !
              END IF
              !
            END DO ! ir
            !
          END DO ! iorb
        END DO ! ispin
      END DO ! iat
      !
      CALL stop_clock(routine_name)
      !
    !---------------------------------------------------------------------------
    END SUBROUTINE set_eta
    !---------------------------------------------------------------------------
    !
    !
  !=============================================================================
  END MODULE mh_eta
  !=============================================================================
