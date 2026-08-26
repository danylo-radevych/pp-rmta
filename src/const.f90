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
  !=============================================================================
  MODULE const
  !=============================================================================
  !!
  !! RMTA constants
  !!
  !
    USE kinds, ONLY: DP
    USE constants, ONLY: bohr_radius_si
    !
    IMPLICIT NONE
    !
    REAL(DP), PARAMETER :: bohrtoang = bohr_radius_si * 1.0E+10_DP
    !! Bohr radius in angstroms
    REAL(DP), PARAMETER :: zero = 0.0_DP
    !! real zero
    REAL(DP), PARAMETER :: one = 1.0_DP
    !! real one
    REAL(DP), PARAMETER :: two = 2.0_DP
    !! real two
    REAL(DP), PARAMETER :: three = 3.0_DP
    !! real three
    REAL(DP), PARAMETER :: four = 4.0_DP
    !! real four
    REAL(DP), PARAMETER :: five = 5.0_DP
    !! real five
    COMPLEX(DP), PARAMETER :: czero = CMPLX(0.0_DP, 0.0_DP, KIND = DP)
    !! complex zero
    !
  !=============================================================================
  END MODULE const
  !=============================================================================
