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
  MODULE neighbor
  !=============================================================================
  !!
  !! Module for finding nearest neighbors and compute r_mt radii based on found
  !! distances
  !!
  !=============================================================================
  !
  !  Danylo Radevych
  !
    USE kinds, ONLY: DP
    !
    IMPLICIT NONE
    !
    ! shorthand name of the module: nn
    ! integer nn: inn_
    !
    PUBLIC :: &
      allocate_nn, deallocate_nn, &
      nneighbors, inn_i, nn_dist
    !
    INTEGER, ALLOCATABLE :: nneighbors(:)
    !! number of the nearest neighbors of each atom
    INTEGER, ALLOCATABLE :: inn_i(:, :)
    !! indeces of the nearest neighbors of each atom
    REAL(DP), ALLOCATABLE :: nn_dist(:)
    !! distances to the nearest neighbors of each atom
  !
  CONTAINS
  !
    !
    !---------------------------------------------------------------------------
    SUBROUTINE allocate_nn
    !---------------------------------------------------------------------------
    !!
    !! Allocate and set nearest-neighbor variables
    !!
    !---------------------------------------------------------------------------
      USE io_global, ONLY: stdout
      USE ions_base, ONLY: nat, tau, ityp
      USE uspp_param, ONLY: upf
      USE cell_base, ONLY: at, bg, alat, omega
      USE const, ONLY: bohrtoang
      USE constants, ONLY: eps6, eps12
      !
      IMPLICIT NONE
      !
      EXTERNAL :: errore
      !
      CHARACTER(len=256) :: routine_name
      !! name of this subroutine
      INTEGER :: ierr
      !! error code
      INTEGER :: iat, jat, nc1, nc2, nc3, inn
      !! iterators
      INTEGER :: nc_max
      !! maximum index of the supercell
      REAL(DP) :: jat_dist
      !! current distance to atom jat
      REAL(DP) :: tmp_dist
      !! temporary distance
      REAL(DP) :: T(3)
      !! translation vector
      REAL(DP) :: dist_vec(3)
      !! distance vector
      !
      routine_name = "allocate_nn"
      !
      ! checking only translational vectors going to neighboring cells
      !
      nc_max = 1
      !
      ALLOCATE(nneighbors(nat), STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, 'Error allocating nneighbors', 1)
      nneighbors(:) = 0
      !
      ALLOCATE(inn_i(nat, nat), STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, 'Error allocating inn_i', 1)
      inn_i(:, :) = -1
      !
      ALLOCATE(nn_dist(nat), STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, 'Error allocating nn_dist', 1)
      nn_dist(:) = -1.0_dp
      !
      !
      ! minimize distance
      !
      DO iat = 1, nat
        DO jat = 1, nat
          !
          !
          jat_dist = -1.0_dp
          !
          DO nc3 = - nc_max, nc_max
            DO nc2 = - nc_max, nc_max
              DO nc1 = - nc_max, nc_max
                !
                T(:) = (nc1 * at(:, 1) + &
                  nc2 * at(:, 2) + &
                  nc3 * at(:, 3))
                !
                dist_vec(:) = ((tau(:, jat) + T(:)) - tau(:, iat)) * &
                  alat
                tmp_dist = dist_vec(1) * dist_vec(1)
                tmp_dist = tmp_dist + &
                  dist_vec(2) * dist_vec(2)
                tmp_dist = tmp_dist + &
                  dist_vec(3) * dist_vec(3)
                !
                tmp_dist = SQRT(tmp_dist)
                !
                IF (((jat_dist < 0._dp) .OR. &
                  (tmp_dist < jat_dist)) .AND. (tmp_dist > eps12)) THEN
                  jat_dist = tmp_dist
                END IF
                !
              END DO ! nc1
            END DO ! nc2
          END DO ! nc3
          !
          IF ((nn_dist(iat) < 0.0_dp) .OR. &
            (jat_dist < nn_dist(iat))) THEN
            nn_dist(iat) = jat_dist
            ! inn_i(iat) = jat
          END IF
          !
          !
        END DO ! jat
        !
        IF (nn_dist(iat) < 0.0_dp) THEN
          WRITE(stdout, '("Distance to nearest-neighbor of atom ", I0, &
            " not found")')
          CALL errore(routine_name, "Error in nn_dist", 1)
        END IF
        !
      END DO ! iat
      !
      !
      ! record all atoms with the minimal distance
      !
      DO iat = 1, nat
        DO jat = 1, nat
          !
          !
          jat_dist = -1.0_dp
          !
          DO nc3 = - nc_max, nc_max
            DO nc2 = - nc_max, nc_max
              DO nc1 = - nc_max, nc_max
                !
                T(:) = (nc1 * at(:, 1) + &
                  nc2 * at(:, 2) + &
                  nc3 * at(:, 3))
                !
                dist_vec(:) = ((tau(:, jat) + T(:)) - tau(:, iat)) * &
                  alat
                tmp_dist = dist_vec(1) * dist_vec(1)
                tmp_dist = tmp_dist + &
                  dist_vec(2) * dist_vec(2)
                tmp_dist = tmp_dist + &
                  dist_vec(3) * dist_vec(3)
                !
                tmp_dist = SQRT(tmp_dist)
                !
                IF (((jat_dist < 0._dp) .OR. &
                  (tmp_dist < jat_dist)) .AND. (tmp_dist > eps12)) THEN
                  jat_dist = tmp_dist
                END IF
                !
              END DO ! nc1
            END DO ! nc2
          END DO ! nc3
          !
          IF ((nn_dist(iat) >= 0.0_dp) .AND. &
            (ABS(jat_dist - nn_dist(iat)) < eps6)) THEN
            nneighbors(iat) = nneighbors(iat) + 1
            inn_i(iat, nneighbors(iat)) = jat
          END IF
          !
        END DO ! jat
        !
        IF (ALL(inn_i(iat, :) < 0.0_dp)) THEN
          WRITE(stdout, '("Nearest-neighbor indices of atom ", I0, &
            " not found")')
          CALL errore(routine_name, "Error in inn_i", 1)
        END IF
        !
      END DO ! iat
      !
      ! printing
      !
      WRITE(stdout, '(/5x, "Nearest neighbors")')
      DO iat = 1, nat
        WRITE(stdout, '(5x)')
        WRITE(stdout, '(6x, "atom #", I4, " out of ", I4, ": ", A4)') &
          iat, nat, upf(ityp(iat))%psd
        WRITE(stdout, '(7x, "neighbor distance: ", F10.8, " bohr = ", &
          F10.8, " A")') nn_dist(iat), nn_dist(iat) * bohrtoang
        !
        WRITE(stdout, '(7x, "neighbor indices:", /8x)', advance='no')
        DO inn = 1, nneighbors(iat)
          WRITE(stdout, '(I3, " ")', advance='no') &
            inn_i(iat, inn)
        END DO
        WRITE(stdout, '()')
        !
        WRITE(stdout, '(7x, "neighbor types:", /8x)', advance='no')
        DO inn = 1, nneighbors(iat)
          WRITE(stdout, '(A3, " ")', advance='no') &
            upf(ityp(inn_i(iat, inn)))%psd
        END DO
        WRITE(stdout, '()')
        !
      END DO ! iat
      WRITE(stdout, '(/5x, /5x)')
      !
      !
      ! CALL errore(routine_name, "Test DONE", 1)
      !
      !
    !---------------------------------------------------------------------------
    END SUBROUTINE allocate_nn
    !---------------------------------------------------------------------------
    !
    !
    !---------------------------------------------------------------------------
    SUBROUTINE deallocate_nn
    !---------------------------------------------------------------------------
    !!
    !! Deallocate nearset-neighbor variables
    !!
    !---------------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      EXTERNAL :: errore
      !
      CHARACTER(len=256) :: routine_name
      !! name of this subroutine
      INTEGER :: ierr
      !! error code
      !
      routine_name = "deallocate_nn"
      !
      DEALLOCATE(nneighbors, STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, &
        'Error deallocating nneighbors', 1)
      !
      DEALLOCATE(inn_i, STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, 'Error deallocating inn_i', 1)
      !
      DEALLOCATE(nn_dist, STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, 'Error deallocating nn_dist', 1)
      !
    !---------------------------------------------------------------------------
    END SUBROUTINE deallocate_nn
    !---------------------------------------------------------------------------
    !
  !=============================================================================
  END MODULE neighbor
  !=============================================================================
