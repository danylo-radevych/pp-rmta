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
    REAL(DP), PARAMETER :: bohrtoang = bohr_radius_si * 1.0E10_DP
    !! Bohr radius in angstroms
    !
  !=============================================================================
  END MODULE const
  !=============================================================================
