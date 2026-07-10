  ! Copyright (C) 2024-2026 Danylo Radevych
  !                                                                            
  ! This file is distributed under the terms of the MIT Non-AI License. 
  ! See the file `LICENSE' in the root directory of the               
  ! present distribution, or 
  ! https://github.com/non-ai-licenses/non-ai-licenses/blob/main/NON-AI-MIT .
  !
  ! Please cite: DOI: https://doi.org/10.1038/s41524-026-02141-7
  !
  !=============================================================================
  MODULE compare
  !=============================================================================
  !!
  !! Module to compare RMTA results with FLAPW.
  !!
  !
  !  D. Radevych
  !
    USE kinds, ONLY : DP
    !
    !
    IMPLICIT NONE
    !
    PRIVATE :: compare_spline_interpolation_rsca
    !
    PUBLIC :: set_ref_pot, &
      ref_nr, ref_pot, ref_ef
    !
    ! REAL(DP), ALLOCATABLE :: ref_r(:)
    REAL(DP), ALLOCATABLE :: ref_pot(:, :)
    INTEGER :: ref_nr
    REAL(DP) :: ref_ef
    !
  !
  CONTAINS
  !
    !---------------------------------------------------------------------------
    SUBROUTINE set_ref_pot(rmt)
    !---------------------------------------------------------------------------
    !!
    !! Read reference FLAPW potential from file.
    !!
      USE io_global, ONLY : stdout
      USE constants, ONLY : eps12
      !
      IMPLICIT NONE
      !
      EXTERNAL :: errore
      !
      REAL(DP), INTENT(in) :: rmt
      !! selected muffin-tin radius
      !
      CHARACTER(len=256) :: routine_name
      !! name of this subroutine
      CHARACTER(len=256) :: tmp_char
      !! temporary character
      INTEGER :: if_data = 9000
      !! unit for data file
      INTEGER :: i, irow, l, ir
      !! iterators
      INTEGER :: ierr
      !! error code
      INTEGER :: nrow, ncol
      !! number of rows and columns
      INTEGER :: lmax = 4
      !! max l for the data from FLAPW
      !
      REAL(DP), ALLOCATABLE :: data(:, :)
      ! REAL(DP), ALLOCATABLE :: ref_r(:)
      ! REAL(DP), ALLOCATABLE :: ref_pot(:)
      REAL(DP) :: rtmp
      !! real temporary variable
      REAL(DP) :: e_fermi
      !! Fermi energy
      REAL(DP) :: vmt
      !! value of the potential at rmt
      !
      routine_name = "set_ref_pot"
      !
      WRITE(stdout, '(/1x, A, &
        ":: START --------------------------------------------------------")') &
        TRIM(routine_name)
      !
      OPEN(UNIT = if_data, FILE = 'pot.txt', STATUS = 'old', ACTION = 'read')
      DO i = 1, 3
        READ(if_data, *), tmp_char
!         WRITE(stdout, '(1x, A256)') tmp_char
      END DO
      !
      READ(if_data, *), e_fermi
      ref_ef = e_fermi
      READ(if_data, *), nrow
      READ(if_data, *), ncol
      WRITE(stdout, '(1x, A, ":: nrow = ", I3, ", ncol = ", I3)') &
        TRIM(routine_name), nrow, ncol
      !
      IF (ALLOCATED(data)) THEN
        DEALLOCATE(data, STAT = ierr)
        IF (ierr /= 0) &
          CALL errore(routine_name, "Error deallocating data", 1)
      END IF
      IF (.NOT. ALLOCATED(data)) THEN
        ALLOCATE(data(nrow, ncol), STAT = ierr)
        IF (ierr /= 0) &
          CALL errore(routine_name, "Error allocating data", 1)
      END IF
      !
      DO irow = 1, nrow
        READ(if_data, *) data(irow, :)
!         DO icol = 1, ncol
!           WRITE(stdout, '(2x, "data(", I3, ", ", I3, ") == ", ES16.8)') &
!             irow, icol, data(irow, icol)
!         END DO ! icol
      END DO ! irow
      !
      CLOSE(if_data)
      !
      ! print reference quantities at rmt
      !
      ! v
      !
      WRITE(stdout, '(2x, A10, ES16.8)') "rmt = ", rmt
      WRITE(stdout, '(2x, A10, ES16.8)') "rmin = ", MINVAL(data(:, 1))
      WRITE(stdout, '(2x, A10, ES16.8)') "rmax = ", MAXVAL(data(:, 1))
      CALL compare_spline_interpolation_rsca(nrow, &
        data(:, 1), data(:, 2), rmt, vmt)
      WRITE(stdout, '(2x, A10, ES16.8)') "r V = ", vmt
      WRITE(stdout, '(2x, A10, ES16.8)') "V = ", vmt / rmt
      WRITE(stdout, '(2x, A10, ES16.8)') "E_F = ", e_fermi
      WRITE(stdout, '(2x, A10, ES16.8)') "V - E_F = ", (vmt / rmt - e_fermi)
      !
      DO l = 1, lmax
        CALL compare_spline_interpolation_rsca(nrow, &
          data(:, 1), data(:, 3 + 2 * (l - 1)), rmt, rtmp)
        WRITE(stdout, '(3x, A2, I1, A, ES16.8)') "L_", l - 1, " = ", rtmp
        CALL compare_spline_interpolation_rsca(nrow, &
          data(:, 1), data(:, 4 + 2 * (l - 1)), rmt, rtmp)
        WRITE(stdout, '(3x, A3, I1, A, ES16.8)') "dL_", l - 1, " / de = ", rtmp
      END DO ! l
      !
      ! copy to ref_pot
      !
      ref_nr = nrow
      !
      IF (ALLOCATED(ref_pot)) THEN
        DEALLOCATE(ref_pot, STAT = ierr)
        IF (ierr /= 0) &
          CALL errore(routine_name, "Error deallocating ref_pot", 1)
      END IF
      IF (.NOT. ALLOCATED(ref_pot)) THEN
        ALLOCATE(ref_pot(ref_nr, 2), STAT = ierr)
        IF (ierr /= 0) &
          CALL errore(routine_name, "Error allocating ref_pot", 1)
      END IF
      !
      ref_pot(:, 1 : 2) = data(:, 1 : 2)
      !
      ! conver r V to V
      !
      DO ir = 1, ref_nr
        IF (ABS(ref_pot(ir, 1)) > eps12) &
          ref_pot(ir, 2) = ref_pot(ir, 2) / ref_pot(ir, 1)
      END DO
      !
      IF (ALLOCATED(data)) THEN
        DEALLOCATE(data, STAT = ierr)
        IF (ierr /= 0) &
          CALL errore(routine_name, "Error deallocating data", 1)
      END IF
      !
      WRITE(stdout, '(1x, A, &
        ":: END ----------------------------------------------------------")') &
        TRIM(routine_name)
      WRITE(stdout, '(1x, "\n")')
      !
      !
      ! CALL errore(routine_name, "TEST DONE", 1)
      !
    !---------------------------------------------------------------------------
    END SUBROUTINE set_ref_pot
    !---------------------------------------------------------------------------
    !
    !
    !---------------------------------------------------------------------------
    SUBROUTINE delete_ref_pot()
    !---------------------------------------------------------------------------
    !!
    !! Destructor
    !!
      !
      IMPLICIT NONE
      !
      CHARACTER(len=256) :: routine_name
      !! name of this subroutine
      INTEGER :: ierr
      !! error code
      EXTERNAL :: errore
      !
      routine_name = "delete_ref_pot"
      !
      DEALLOCATE(ref_pot, STAT = ierr)
      IF (ierr /= 0) &
        CALL errore(routine_name, "Error deallocating ref_pot", 1)
      !
    !---------------------------------------------------------------------------
    END SUBROUTINE delete_ref_pot
    !---------------------------------------------------------------------------
    !
    !
    !---------------------------------------------------------------------------
    SUBROUTINE compare_spline_interpolation_rsca(xdim, x, y, xf, yf)
    !---------------------------------------------------------------------------
    !!
    !! Spline interpolation to real scalar.
    !!
      USE kinds, ONLY: DP
      USE splinelib, ONLY: spline, splint
      USE constants, ONLY: eps12
      !
      IMPLICIT NONE
      !
      EXTERNAL :: errore
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
      routine_name = "compare_spline_interpolation_rsca"
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
      IF (ALLOCATED(d2y)) DEALLOCATE(d2y)
      IF (.NOT. ALLOCATED(d2y)) THEN
        ALLOCATE(d2y(xdim), STAT = ierr)
        IF (ierr /= 0) &
          CALL errore(routine_name, "Error allocating d2y", 1)
      END IF
      !
      ! prepare spline
      CALL spline(x(1 : xdim), y(1 : xdim), startu, startd, d2y)
      !
      !
      IF ((xf <= xmax) .AND. (xf >= xmin)) THEN
        yf = splint(x(1 : xdim), y(1 : xdim), d2y, xf)
      ELSE
        CALL errore(routine_name, "interpolation out of bounds", 1)
      END IF
      !
      !
      IF (yf /= yf) THEN
        CALL errore(routine_name, 'yf value is NAN', 1)
      END IF
      !
      IF (ABS(yf) < eps12) THEN
        yf = 0.0_dp
      END IF
      !
      IF (ALLOCATED(d2y)) THEN
        DEALLOCATE(d2y, STAT = ierr)
        IF (ierr /= 0) &
          CALL errore(routine_name, 'Error deallocating d2y', 1)
      END IF
      !
    !---------------------------------------------------------------------------
    END SUBROUTINE compare_spline_interpolation_rsca
    !---------------------------------------------------------------------------
    !
    !
  !=============================================================================
  END MODULE compare
  !=============================================================================
