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
  MODULE sym_type
  !=============================================================================
  !!
  !! Module to group atoms of the same chemical element into
  !! symmetry types related by the site-symmetry operations.
  !!
  !=============================================================================
  !
  !  Danylo Radevych
  !
    !
    USE kinds, ONLY : DP
    !
    IMPLICIT NONE
    !
    ! shorthand name of the module: st
    ! integer st: ist_
    !
    PUBLIC :: allocate_st, deallocate_st, &
              nst, st_wt, ist_i, ist_nat, st_name, ist_ityp
    !
    CHARACTER(len=2), ALLOCATABLE :: st_name(:)
    !! name of each symmetry type
    INTEGER :: nst
    !! number of symmetry types
    INTEGER, ALLOCATABLE :: ist_nat(:)
    !! number of atoms in each symmetry type
    INTEGER, ALLOCATABLE :: ist_i(:)
    !! returns symmetry type index of each atom
    INTEGER, ALLOCATABLE :: ist_ityp(:)
    !! corresponding chemical type of each symmetry type
    REAL(DP), ALLOCATABLE :: st_wt(:)
    !! weight of the atom in its symmetry type
    ! TODO weights are wrong right now
    !
  !
  CONTAINS
  !
    !
    !
    !---------------------------------------------------------------------------
    SUBROUTINE allocate_st
    !---------------------------------------------------------------------------
    !!
    !! Allocate and set symmetry type variables
    !!
    !---------------------------------------------------------------------------
      USE io_global, ONLY: stdout
      USE ions_base, ONLY: nat, ityp
      USE symm_base, ONLY: nrot, irt
      USE uspp_param, ONLY: upf
      !
      IMPLICIT NONE
      !
      EXTERNAL :: errore
      !
      CHARACTER(len=256) :: routine_name
      !! name of this subroutine
      INTEGER :: irot
      !! symmatry index
      INTEGER :: iat
      !! atomic index
      INTEGER :: counter
      !! counts number of symmetry partners
      INTEGER :: ierr
      !! error code
      INTEGER :: ist
      !! index of the symmetry type
      REAL(DP) :: sum_st_wt
      !! sum of the symmetry weights
      !
      routine_name = "allocate_st"
      ! CALL start_clock(routine_name)
      !
      ! allocate symmetry weights for each atom and set them to 1
      !
      IF (ALLOCATED(st_wt)) &
        DEALLOCATE(st_wt)
      IF (.NOT. ALLOCATED(st_wt)) THEN
        ALLOCATE(st_wt(nat), STAT = ierr)
        IF (ierr /= 0) CALL errore(routine_name, &
          'Error allocating st_wt', 1)
      END IF
      st_wt(:) = 1.0_dp
      nst = 0
      sum_st_wt = 0.0_dp
      !
      ALLOCATE(ist_i(nat), STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, 'Error allocating ist_i', 1)
      ist_i(:) = 0
      !
      WRITE(stdout, '(/5x, "Checking irt, whether symmetry operations ", &
        & "relate different atoms:")')
      DO iat = 1, nat
        WRITE(stdout, '(5x)')
        WRITE(stdout, '(6x, "atom #", I4, " out of ", I4, ": ", A4)') &
          iat, nat, upf(ityp(iat))%psd
        counter = 1
        !
        IF (ist_i(iat) == 0) THEN
          IF (iat == 1) THEN
            ist_i(iat) = 1
          ELSE
            ist_i(iat) = MAXVAL(ist_i(:)) + 1
          END IF
        END IF
        !
        DO irot = 1, nrot
          IF ((irt(irot, iat) > 0) .AND. (irt(irot, iat) /= iat) .AND. &
            (irt(irot, iat) <= nat)) THEN
            WRITE(stdout, '(7x, "irt(", I2, ", ", I2, ") = ", I0, ": ", A)') &
              irot, iat, irt(irot, iat), upf(ityp(irt(irot, iat)))%psd
            counter = counter + 1
            !
            IF (ist_i(irt(irot, iat)) == 0) THEN
              ist_i(irt(irot, iat)) = ist_i(iat)
            END IF
            !
          END IF
        END DO ! irot
        WRITE(stdout, '(6x, "symmetry type: ", I4)') ist_i(iat)
        WRITE(stdout, '(6x, "# of partners: ", I4)') counter
        ! st_wt(iat) = st_wt(iat) / counter
        ! sum_st_wt = sum_st_wt + st_wt(iat)
        ! WRITE(stdout, '(6x, "symmetry weight: ", F10.8)') st_wt(iat)
      END DO ! iat
      WRITE(stdout, '(/5x)')
      !
      ! WRITE(stdout, '(5x, "sum of symmetry weights: ", F10.8)') sum_st_wt
      !
      ! ! TODO
      ! nst = NINT(sum_st_wt)
      nst = MAXVAL(ist_i(:))
      !
      IF (nst /= MAXVAL(ist_i(:))) THEN
        CALL errore(routine_name, &
          "Number of the symmetry types is not equal to max symmetry type index", 1)
      END IF
      !
      WRITE(stdout, '(5x, "# of symmetry types: ", I4)') nst
      WRITE(stdout, '(/5x)')
      !
      ! number of atoms in each symmetry type
      !
      ALLOCATE(ist_nat(nst), STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, 'Error allocating ist_nat', 1)
      ist_nat(:) = 0
      !
      ALLOCATE(st_name(nst), STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, 'Error allocating st_name', 1)
      st_name(:) = ""
      !
      ALLOCATE(ist_ityp(nst), STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, 'Error allocating ist_ityp', 1)
      ist_ityp(:) = 0
      !
      DO iat = 1, nat
        ist_nat(ist_i(iat)) = ist_nat(ist_i(iat)) + 1
        !
        IF (st_name(ist_i(iat)) == "") THEN
          st_name(ist_i(iat)) = upf(ityp(iat))%psd
        END IF
        !
      END DO ! iat
      !
      DO ist = 1, nst
        WRITE(stdout, '(5x, "type #: ", I4, ": ", A, " with ", I4, &
          & " members:")') &
          ist, st_name(ist), ist_nat(ist)
        !
        WRITE(stdout, '(7x)', ADVANCE = "no")
        DO iat = 1, nat
          IF (ist_i(iat) == ist) THEN
            !
            WRITE(stdout, '(I4, " ")', ADVANCE = "no") iat
            !
            st_wt(iat) = 1.0_dp / ist_nat(ist)
            ist_ityp(ist) = ityp(iat)
            !
          END IF
        END DO ! iat
        WRITE(stdout, '(7x)')
        !
      END DO ! ist
      WRITE(stdout, '(/7x, /7x)')
      !
      ! CALL errore(routine_name, "Test DONE", 1)
      !
      ! CALL stop_clock(routine_name)
      !
      !
    !---------------------------------------------------------------------------
    END SUBROUTINE allocate_st
    !---------------------------------------------------------------------------
    !
    !
    !---------------------------------------------------------------------------
    SUBROUTINE deallocate_st
    !---------------------------------------------------------------------------
    !!
    !! Deallocate symmetry type variables
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
      routine_name = "deallocate_st"
      !
      DEALLOCATE(st_wt, STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, 'Error deallocating st_wt', 1)
      !
      DEALLOCATE(ist_nat, STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, 'Error deallocating ist_nat', 1)
      !
      DEALLOCATE(ist_i, STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, 'Error deallocating ist_i', 1)
      !
      DEALLOCATE(st_name, STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, 'Error deallocating st_name', 1)
      !
      DEALLOCATE(ist_ityp, STAT = ierr)
      IF (ierr /= 0) CALL errore(routine_name, 'Error deallocating ist_ityp', 1)
      !
    !---------------------------------------------------------------------------
    END SUBROUTINE deallocate_st
    !---------------------------------------------------------------------------
    !
    !
  !=============================================================================
  END MODULE sym_type
  !=============================================================================
