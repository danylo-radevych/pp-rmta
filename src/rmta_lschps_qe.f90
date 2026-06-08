    !
    !---------------------------------------------------------------------------
    SUBROUTINE rmta_lschps_qe(mode, z, eps, grid, nin, n, l, e, v, u, up, nstop)
    !---------------------------------------------------------------------------
    !!
    !! Modified radial Schrodinger equation solver from atomic/src/lschps.f90
    !! tailored for RMTA needs.
    !!
    !
    ! D. Radevych
    !
      !
      ! integrates radial pauli-type scalar-relativistic equation
      ! on a logarithmic grid
      ! modified routine to be used in finding norm-conserving
      ! pseudopotential
      !
      ! on input:
      !   mode = 1 find energy and wavefunction of bound states,
      !            scalar-relativistic (all-electron)
      !   mode = 2 find energy and wavefunction of bound state,
      !            nonrelativistic (pseudopotentials)
      !   mode = 3 fixed-energy calculation, for logarithmic derivatives
      !   mode = 4 find energy which produces a specified logarithmic
      !            derivative (nonrelativistic, pseudopotentials)
      !   mode = 5 is for pseudopotential to produce wavefunction beyond
      !            radius used for pseudopotential construction
      !   z    = atomic number
      !   eps  = convergence factor: eiganvalue is considered converged if
      !          the correction to eigenvalue is smaller in magnitude than
      !          eps times the magnitude of the current guess
      !   grid = structure containing radial grid information
      !   l, n = main and angular quantum numbers
      !   e    = starting estimate of the energy (mode=1,2)
      !          fixed energy at which the wavefctn is calculated (mode=3,4)
      !   v(i) = self-consistent potential
      !   nin  = integration up to r(nin) (mode=3,4,5)
      !
      ! on output:
      !   e    = final energy (mode=1,2)
      !   u(i) = radial wavefunction (defined as the radial part of the wavefct
      !          multiplied by r)
      !   nstop= 0 if regular termination, 1 otherwise
      !   nin  = last grid point for which the wavefct is calculated (mode=1,2)
      !
      !
      !
      ! In the future, to use pure lschps from atomic/src/lschps.f90:
      ! - set cau_fact from ld1inc to c_au from constants
      !
      !
      !
      USE kinds, ONLY : DP
      USE radial_grids, ONLY: radial_grid_type
      ! USE ld1inc, ONLY : cau_fact
      USE constants, ONLY: c_au
      IMPLICIT NONE
      !
      ! I/O variables
      !
      INTEGER, INTENT (in) :: mode, n, l
      REAL(DP), INTENT(in) :: z, eps
      TYPE (radial_grid_type), INTENT(in) :: grid
      REAL(DP), INTENT(in) :: v(grid%mesh)
      INTEGER, INTENT(inout) :: nin
      REAL(DP), INTENT(inout) :: e
      INTEGER, INTENT(out) :: nstop
      REAL(DP), INTENT(inout) :: u(grid%mesh)
      REAL(DP), INTENT(inout) :: up(grid%mesh)
      !
      ! local variables
      !
      INTEGER, PARAMETER :: maxter=1
      REAL(DP), EXTERNAL:: aei, aeo, aii, aio
      ! arrays  used as work space
      REAL(DP),ALLOCATABLE :: upp(:),cf(:),dv(:),fr(:),frp(:)
      REAL(DP):: al, als, cn
      REAL(DP):: de, emax, emin
      REAL(DP):: fss, gamma, ro, sc
      REAL(DP):: sls, sn, uld, uout,  upin, upout
      REAL(DP):: xkap
      INTEGER:: i, it, mmax, n_it, node, mch, ierr
      !
      LOGICAL :: lnormalize = .TRUE.
      !
      ! cau_fact = c_au
      u(:) = 0.0_dp
      !
      nstop=0
      al = grid%dx
      mmax = grid%mesh
      !
!       ALLOCATE(up(mmax), stat=ierr)
      ALLOCATE(upp(mmax), stat=ierr)
      ALLOCATE(cf(mmax), stat=ierr)
      ALLOCATE(dv(mmax), stat=ierr)
      ALLOCATE(fr(mmax), stat=ierr)
      ALLOCATE(frp(mmax), stat=ierr)
      !
      uld = 0.0_dp
      !
      !
      IF (mode == 1 .OR. mode == 3) THEN
         !
         ! relativistic calculation
         !
         ! fss=(1.0_dp/137.036_dp)**2
         fss = (1.0_dp / c_au)**2
         IF (l == 0) THEN
            gamma = SQRT(1.0_dp - fss * z**2)
         ELSE
            gamma = (l * SQRT(l**2 - fss * z**2) + &
             (l+1) * SQRT((l + 1)**2 - fss * z**2)) / (2 * l + 1)
         ENDIF
      ELSE
         !
         ! non-relativistic calculation
         !
         fss = 1.0e-20_dp
         gamma = l + 1
      ENDIF
      WRITE(*, *) fss
      !
      sls = l * (l+1)
      !
      ! emin, emax = estimated bounds for e
      !
      IF (mode == 1 .or. mode == 2) THEN
         emax=v(mmax)+sls/grid%r(mmax)**2
         emin=0.0_dp
         DO i=1,mmax
            emin=min(emin,v(i)+sls/grid%r(i)**2)
         ENDDO
         IF(e > emax) e=1.25_dp*emax
         IF(e < emin) e=0.75_dp*emin
         IF(e > emax) e=0.5_dp*(emax+emin)
      ELSE IF(mode == 4) THEN
         emax=e + 10.0_dp
         emin=e - 10.0_dp
      ELSE IF (mode == 3) THEN
         emax = e
         emin = e
      ENDIF
      !
      DO i = 1, 4
         u(i) = 0.0_dp
         up(i) = 0.0_dp
         upp(i) = 0.0_dp
      ENDDO
      als = al**2
      !
      ! calculate dv/dr for darwin correction
      !
      CALL derv(mmax, al, grid%r, v, dv)
      !
      ! starting of loop on energy for bound state
      !
      DO n_it = 1, maxter
         !
         ! coefficient array for u in differential eq.
         DO i=1,mmax
            cf(i) = als * (sls + (v(i) - e) * grid%r(i)**2)
         ENDDO
         !
         ! find classical turning point for matching
         !
         IF (mode == 1 .OR. mode == 2) THEN
            DO i = mmax, 2, -1
               IF (cf(i-1) <= 0.0_dp .AND. cf(i) > 0.0_dp) THEN
                  mch = i
                  GOTO 40
               ENDIF
            ENDDO
            PRINT '('' warning: wfc '',2i2,'' no turning point'')', n, l
            e = 0.0_dp
            DO i = 1, mmax
               u(i) = 0.0_dp
            ENDDO
            nstop=1
            GOTO 999
         ELSE
            mch = nin
         ENDIF
    40 CONTINUE
         !
         !  relativistic coefficient arrays for u(fr) and up(frp).
         DO i = 1, mmax
            fr(i) = als * (grid%r(i)**2) * 0.25_dp * (-fss * (v(i) - e)**2 + &
                 fss * dv(i) / &
                 (grid%r(i) * (1.0_dp + 0.25_dp * fss * (e - v(i)))))
            frp(i) = - al * grid%r(i) * 0.25_dp * fss * dv(i) / &
              (1.0_dp + 0.25_dp * fss * (e - v(i)))
         ENDDO
         !
         ! start wavefunction with series
         !
         DO i = 1,4
            u(i) = grid%r(i)**gamma
!             WRITE(*, *) gamma
            up(i) = al * gamma * grid%r(i)**gamma
            upp(i) = (al + frp(i)) * up(i) + (cf(i) + fr(i)) * u(i)
         ENDDO
         !
         ! outward integration using predictor once, corrector
         ! twice
         node = 0
         !
         DO i = 4, mch - 1
            u(i + 1) = u(i) + aeo(up, i)
            up(i + 1) = up(i) + aeo(upp, i)
            DO it = 1, 2
               upp(i + 1) = (al + frp(i + 1)) * up(i + 1) + &
                 (cf(i + 1) + fr(i + 1)) * u(i + 1)
               up(i + 1) = up(i) + aio(upp, i)
               u(i + 1) = u(i) + aio(up, i)
            ENDDO
            IF(u(i+1) * u(i) <= 0.0_dp) node = node + 1
         ENDDO
         !
         uout = u(mch)
         upout = up(mch)
         !
         IF(node - n + l + 1 == 0 .OR. mode == 3 .OR. mode == 5) THEN
            !
            IF(mode == 1 .or. mode == 2) THEN
               !
               ! start inward integration at 10*classical turning
               ! point with simple exponential
               nin=mch+2.3_dp/al
               IF(nin+4 > mmax) nin=mmax-4
               xkap=sqrt(sls/grid%r(nin)**2 + 2.0_dp*(v(nin)-e))
               !
               DO i=nin,nin+4
                  u(i)=exp(-xkap*(grid%r(i)-grid%r(nin)))
                  up(i)=-grid%r(i)*al*xkap*u(i)
                  upp(i)=(al+frp(i))*up(i)+(cf(i)+fr(i))*u(i)
               ENDDO
               !
               ! integrate inward
               !
               DO i=nin,mch+1,-1
                  u(i-1)=u(i)+aei(up,i)
                  up(i-1)=up(i)+aei(upp,i)
                  DO it=1,2
                     upp(i-1)=(al+frp(i-1))*up(i-1)+(cf(i-1)+fr(i-1))*u(i-1)
                     up(i-1)=up(i)+aii(upp,i)
                     u(i-1)=u(i)+aii(up,i)
                  ENDDO
               ENDDO
               !
               ! scale outside wf for continuity
               sc=uout/u(mch)
               !
               DO i=mch,nin
                  up(i)=sc*up(i)
                  u (i)=sc*u (i)
               ENDDO
               !
               upin=up(mch)
               !
            ELSE
               !
               upin = uld * uout
               !
            ENDIF
            !
            ! perform normalization sum
            !
            ro = grid%r(1) * exp(-0.5_dp * grid%dx)
            sn = ro**(2.0_dp * gamma + 1.0_dp) / (2.0_dp * gamma + 1.0_dp)
            !
            DO i = 1, nin - 3
               sn = sn + al * grid%r(i) * u(i)**2
            ENDDO
            !
            sn=sn + al*(23.0_dp*grid%r(nin-2)*u(nin-2)**2 &
                 + 28.0_dp*grid%r(nin-1)*u(nin-1)**2 &
                 +  9.0_dp*grid%r(nin  )*u(nin  )**2)/24.0_dp
            !
            ! normalize u
            IF (lnormalize) THEN
              cn = 1.0_dp / SQRT(sn)
              uout =cn * uout
              upout = cn * upout
              upin = cn * upin
              !
              DO i = 1, nin
                 up(i) = cn * up(i)
                 u(i) = cn * u(i)
              ENDDO
              !
              ! DR: make up true du / dr
              DO i = 1, nin
                up(i) = up(i) / grid%r(i) / al
              END DO
              !
            END IF
            !
            DO i = nin + 1, mmax
               u(i) = 0.0_dp
            ENDDO
            !
            ! exit for fixed-energy calculation
            !
            IF(mode == 3 .or. mode == 5) GOTO 999
            !
            ! perturbation theory for energy shift
            de = uout * (upout - upin) / (al * grid%r(mch))
            !
            ! convergence test and possible exit
            !
            IF ( abs(de) < max(abs(e), 0.2_dp)*eps) GOTO 999
            !
            IF(de > 0.0_dp) THEN
               emin=e
            ELSE
               emax=e
            ENDIF
            e=e+de
            IF(e > emax .or. e < emin) e=0.5_dp*(emax+emin)
            !
         ELSEIF(node-n+l+1 < 0) THEN
            ! too few nodes
            emin=e
            e=0.5_dp*(emin+emax)
            !
         ELSE
            ! too many nodes
            emax=e
            e=0.5_dp*(emin+emax)
         ENDIF
      ENDDO ! n_it

      PRINT '('' warning: wfc '',2i2,'' not converged'')', n, l
      ! u=0.0_dp
      nstop=1
      !
      ! deallocate arrays and exit
      !
    999 CONTINUE
      WRITE(*, *) "emin", emin
      WRITE(*, *) "emax", emax
      WRITE(*, *) "e", e
      DEALLOCATE(frp)
      DEALLOCATE(fr)
      DEALLOCATE(dv)
      DEALLOCATE(cf)
      DEALLOCATE(upp)
!       DEALLOCATE(up)
      RETURN
      !
    !---------------------------------------------------------------------------
    END SUBROUTINE rmta_lschps_qe
    !---------------------------------------------------------------------------
    !
