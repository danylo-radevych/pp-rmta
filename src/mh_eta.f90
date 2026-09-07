  !---------------------------------------------------------------------------!
  ! Copyright (C) 2024-2026 Danylo Radevych                                   !
  !                                                                           !
  ! PP-RMTA is free software: you can redistribute it and/or modify           !
  ! it under the terms of the GNU General Public License as published by the  !
  ! Free Software Foundation, either version 3 of the License, or             !
  ! (at your option) any later version.                                       !
  !                                                                           !
  ! PP-RMTA is distributed in the hope that it will be useful,                !
  ! but WITHOUT ANY WARRANTY; without even the implied warranty of            !
  ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                      !
  ! See the GNU General Public License for more details.                      !
  !                                                                           !
  ! You should have received a copy of the GNU General Public License along   !
  ! with PP-RMTA. If not, see <https://www.gnu.org/licenses/>.                !
  !                                                                           !
  ! Please cite: DOI: https://doi.org/10.1038/s41524-026-02141-7              !
  !---------------------------------------------------------------------------!
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
      dos_nlr, dos_nr, &
      dos_n, luse_tot_dos, &
      mll1, etall1)
    !---------------------------------------------------------------------------
    !!
    !! Computes McMillan-Hopfield eta_l
    !!
      USE constants, ONLY: eps12
      USE const, ONLY: zero, one, two, three
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
      etall1(:, :, :, :) = zero
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
              rtmp_n = rtmp_n * nat ! turn into total DOS, per whole cell
              !
              IF (rtmp_n < zero) &
                CALL errore(routine_name, "n < 0", 1)
              !
              rtmp_nl = dos_nlr(ir - imin + 1, iorb, ispin, iat)
              IF (rtmp_nl < zero) &
                CALL errore(routine_name, "n_l < 0", 1)
              !
              rtmp_nl1 = dos_nlr(ir - imin + 1, iorb + 1, ispin, iat)
              IF (rtmp_nl1 < zero) &
                CALL errore(routine_name, "n_{l+1} < 0", 1)
              !
              rtmp_mll1 = mll1(ir, iorb, ispin, iat)
              !
              IF (ABS(rtmp_n) > eps12) THEN
                !
                ! \eta_{l, l+1}
                !
                etall1(ir - imin + 1, iorb, ispin, iat) = &
                  two * (l + one) / &
                  (two * l + one) / (two * l + three) * &
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
              ELSE
                !
                CALL errore(routine_name, "problem with total DOS n", 1)
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
