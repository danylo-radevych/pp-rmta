  ! Copyright (C) 2024-2026 Danylo Radevych
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  ! Please cite: DOI: https://doi.org/10.1038/s41524-026-02141-7
  !
  !=============================================================================
  MODULE hamann
  !=============================================================================
  !!
  !! Module with subroutines for radial Schrodinger equation from
  !! Hamann's ONCVPSP code.
  !!
  !=============================================================================
  !
  !  Modified by D. Radevych
  !
    !
    USE kinds, ONLY : DP
    !
    IMPLICIT NONE
    !
    !
    PUBLIC :: hamann_lschps, hamann_lschvkbs
    !
    PRIVATE :: hmn_vpinteg, hmn_lschkb, hmn_lschps, hmn_vkboutwf
    !
  !
  CONTAINS
  !
    !---------------------------------------------------------------------------
    SUBROUTINE hamann_lschps(ll, ierr, ee, rr, vv, uu, up, mmax, mch)
    !---------------------------------------------------------------------------
    !!
    !! outward integration of Srcroedinger equation for
    !! semi-local pseudopotential
    !! on a logarithmic mesh
    !!
    !
    !nn  principal quantum number (not used actually)
    !ll  angular-momentum quantum number
    !ierr  non-zero return if error
    !ee  bound-state energy, input guess and output calculated value
    !rr  log radial mesh
    !vv  semi-local atomic pseudopotential
    !uu  output radial wave function (*rr)
    !up  d(uu)/dr
    !mmax  size of log grid
    !mch matching mesh point for inward-outward integrations
    !
      !
      IMPLICIT NONE
      !
      !Input variables
      INTEGER, INTENT(in) :: mmax,mch
      REAL(DP), INTENT(in) :: rr(mmax),vv(mmax)
      INTEGER, INTENT(in) :: ll
      !
      !Output variables
      REAL(DP), INTENT(inout) :: uu(mmax), up(mmax)
      REAL(DP), INTENT(inout) :: ee
      INTEGER, INTENT(inout) :: ierr
      !
      !Local variables
      REAL(DP) :: amesh,al
      REAL(DP) :: aeo, aio, als, cn
      REAL(DP) :: ro
      REAL(DP) :: sls, sn, uout, upout
      INTEGER :: ii, it
      !
      REAL(DP), allocatable :: upp(:),cf(:)
      !
      ALLOCATE(upp(mmax),cf(mmax))
      !
      al = 0.01d0 * dlog(rr(101) / rr(1))
      amesh = dexp(al)
      !
      ierr = 0
      !
      sls=ll*(ll+1)
      !
      ! null arrays to remove leftover garbage
      !
      uu(:)=0.0d0
      up(:)=0.0d0
      upp(:)=0.0d0
      !
      als=al**2
      !
      ! coefficient array for u in differential eq.
      DO ii=1,mmax
        ! cf(ii)=als*sls + 2.0d0*als*(vv(ii)-ee)*rr(ii)**2
        cf(ii)=als*sls + als*(vv(ii)-ee)*rr(ii)**2 ! DR, vv and ee are in Ry already
      END DO
      !
      ! start wavefunction with series
      !
      DO ii=1,4
        uu(ii)=rr(ii)**(ll+1)
        up(ii)=al*(ll+1)*rr(ii)**(ll+1)
        upp(ii)=al*up(ii)+cf(ii)*uu(ii)
      END DO
      !
      ! outward integration using predictor once, corrector
      ! twice
      !
      DO ii=4,mch-1
        uu(ii+1)=uu(ii)+aeo(up,ii)
        up(ii+1)=up(ii)+aeo(upp,ii)
        DO it=1,2
          upp(ii+1)=al*up(ii+1)+cf(ii+1)*uu(ii+1)
          up(ii+1)=up(ii)+aio(upp,ii)
          uu(ii+1)=uu(ii)+aio(up,ii)
        END DO
      END DO
      !
      uout=uu(mch)
      upout=up(mch)
      !
      !perform normalization sum
      !
      ro=rr(1)/dsqrt(amesh)
      sn=ro**(2*ll+3)/(2*ll+3)
      !
      DO ii=1,mch-3
        sn=sn+al*rr(ii)*uu(ii)**2
      END DO
      !
      sn =sn + al*(23.0d0*rr(mch-2)*uu(mch-2)**2 + &
        28.0d0*rr(mch-1)*uu(mch-1)**2 + &
        9.0d0*rr(mch  )*uu(mch  )**2)/24.0d0
      !
      !normalize u
      !
      cn=1.0d0/dsqrt(sn)
      uout=cn*uout
      upout=cn*upout
      !
      DO ii=1,mch
        up(ii)=cn*up(ii)
        uu(ii)=cn*uu(ii)
      END DO
      DO ii=mch+1,mmax
        uu(ii)=0.0d0
      END DO
      !
      ! DR: turn up = du into true du / dr
      !
      DO ii=1, mch
        up(ii) = up(ii) / rr(ii) / al
      END DO
      !
      DEALLOCATE(upp,cf)
      !
      RETURN
      !
    !---------------------------------------------------------------------------
    END subroutine hamann_lschps
    !---------------------------------------------------------------------------
    !
    !
    !---------------------------------------------------------------------------
    SUBROUTINE hamann_lschvkbs(ll, ivkb, ee, rr, vloc, vkb, &
      evkb, uu, up, mmax, mch)
    !---------------------------------------------------------------------------
      ! Normalized scattering state for pseudopotential, fully non-local unless
      !
      ! ivkb=0, which should only happen if ll=lloc
      !ll  angular-momentum quantum number
      !ivkb  = 0, 1 or 2 VKB proectors to be used
      !ee  scattering state energy
      !rr  log radial mesh
      !vloc  local pseudopotential
      !vkb  VKB projectors
      !evkb  coefficients of VKB projectors
      !uu  output radial wave function (*rr)
      !up  d(uu)/dr
      !mmax  size of log grid
      !mch  index of radius to which uu is computed

    IMPLICIT NONE
    INTEGER, PARAMETER :: dp=kind(1.0d0)

    !Input variables
    INTEGER :: mmax,mch
    REAL(DP) :: rr(mmax),vloc(mmax),vkb(mmax,2),evkb(2)
    REAL(DP) :: ee
    INTEGER :: ivkb,ll

    !Output variables
    REAL(DP) :: uu(mmax),up(mmax)


    !Local variables
    REAL(DP) :: amesh,al
    REAL(DP) :: cn,ro,sn
    INTEGER :: ii,node


    al = 0.01d0 * dlog(rr(101) / rr(1))
    amesh = dexp(al)

    ! null arrays to remove leftover garbage

    uu(:)=0.0d0
    up(:)=0.0d0

    CALL hmn_vkboutwf(ll, ivkb, ee, vkb, evkb, rr, vloc, uu, up, node, mmax, mch)

    !perform normalization sum

    ro=rr(1)/dsqrt(amesh)
    sn=ro**(2*ll+3)/(2*ll+3)

    DO ii=1,mch-3
    sn=sn+al*rr(ii)*uu(ii)**2
    END DO

    sn =sn + al*(23.0d0*rr(mch-2)*uu(mch-2)**2 &
    &       + 28.0d0*rr(mch-1)*uu(mch-1)**2 &
    &      +  9.0d0*rr(mch  )*uu(mch  )**2)/24.0d0

    !normalize u

    cn=1.0d0/dsqrt(sn)

    DO ii=1,mch
    up(ii)=cn*up(ii)
    uu(ii)=cn*uu(ii)
    END DO
    DO ii=mch+1,mmax
    uu(ii)=0.0d0
    END DO
    !
    ! DR: make up true du / dr
    DO ii=1, mch
    up(ii) = up(ii) / rr(ii) / al
    END DO
    !

    RETURN
  !-----------------------------------------------------------------------------
    END SUBROUTINE hamann_lschvkbs
    !---------------------------------------------------------------------------
    !
    !
    !---------------------------------------------------------------------------
    SUBROUTINE hmn_vkboutwf(ll,nvkb,ep,vkb,evkb,rr,vloc,uu,up,node,mmax,mch)
    !---------------------------------------------------------------------------
    !!
    !! computes Vanderbilt / Kleinman-Bylander outward-integrated wave functions
    !!
    !ll  angular momentum
    !nvkb  switch for 1 or 2 projedtors
    !ep  energy at which wave function is to be calculated
    !vkb  Vanderbilt-Kleinman-Bylander projectors for this l
    !evkb  projector coefficients
    !rr  log radial mesh
    !vloc  local pseudopotential
    !uu  wave function
    !up  1st derivative of uu
    !node  count of number of nodes from 0 to rr(mch)
    !mmax  dimension of log mesh
    !mch  index of radius to which wave function is to be integrated

    IMPLICIT NONE
    INTEGER, PARAMETER :: dp=kind(1.0d0)

    !Input variables
    REAL(DP) :: rr(mmax),vloc(mmax),vkb(mmax,nvkb),evkb(nvkb)
    REAL(DP) :: ep
    INTEGER nvkb,ll,mmax,mch

    !Output variables
    REAL(DP) :: uu(mmax),up(mmax)
    INTEGER node

    !Local variables
    REAL(DP), allocatable ::  phi(:,:), phip(:,:)
    REAL(DP), allocatable ::  gg0(:), gg(:,:)
    INTEGER, allocatable :: ipiv(:)

    REAL(DP) :: rcut
    INTEGER ii,jj,krc,ierr,info

    uu(:)=0.0d0
    up(:)=0.0d0

    ! homogeneous solution
    CALL hmn_lschps(ll,ierr,ep,rr,vloc,uu,up,mmax,mch)

    ! default lower bound for node counting when nvkb==0
    rcut=0.1d0

    IF(nvkb/=0) THEN

    ALLOCATE(phi(mmax,nvkb),phip(mmax,nvkb))
    ALLOCATE(gg(nvkb,nvkb),gg0(nvkb))
    ALLOCATE(ipiv(nvkb))

    phi(:,:)=0.0d0
    phip(:,:)=0.0d0
    gg(:,:)=0.0d0
    gg0(:)=0.0d0

    ! inhomogeneous solutions
    DO ii=1,nvkb
    CALL hmn_lschkb(ll,ierr,ep,vkb(1,ii),rr,vloc,phi(1,ii),phip(1,ii),mmax,mch)
    END DO


    ! projector matrix elements and coefficient matrix
    DO jj=1,nvkb
    CALL hmn_vpinteg(uu,vkb(1,jj),mch,2*ll+2,gg0(jj),rr)
    gg0(jj)=evkb(jj)*gg0(jj)
    !    gg0(jj)=evkb(jj)*gg0(jj) / 2.0_dp ! DR
    DO ii=1,nvkb
    CALL hmn_vpinteg(phi(1,ii),vkb(1,jj),mch,2*ll+2,gg(jj,ii),rr)
    gg(jj,ii)=-evkb(jj)*gg(jj,ii)
    !     gg(jj,ii)=-evkb(jj)*gg(jj,ii) / 2.0_dp ! DR
    END DO
    gg(jj,jj)=1.0d0+gg(jj,jj)
    END DO

    ! solve linear equations for coefficients

    !    SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )

    CALL dgesv(nvkb, 1, gg, nvkb, ipiv, gg0, nvkb, info)
    IF(info/=0) THEN
    write(6,'(/a,i4)') 'vkboutwf: dgesv ERROR, stopping info =',info
    STOP
    END IF

    ! output wave functions
    DO jj=1,nvkb
    uu(:)=uu(:)+gg0(jj)*phi(:,jj)
    up(:)=up(:)+gg0(jj)*phip(:,jj)
    END DO

    DEALLOCATE(phi,phip)
    DEALLOCATE(gg,gg0)
    DEALLOCATE(ipiv)


    ! rcut is lower cutoff for node counting to avoid small-r noise
    ! this method of "finding" rc is cumbersome but it avoids a lot of re-coding
    ! to simply pass rc or irc along.
    DO ii=mch,1,-1
    IF(dabs(vkb(ii,1))>0.0d0) THEN
    krc=ii+1
    exit
    END IF
    END DO
    ! the constants below might need future adjustment
    rcut=dmin1(0.5d0, 0.25d0*rr(krc))

    END IF !nvkb>0

    node=0
    DO ii=6,mch
    IF(rr(ii)>rcut .and. uu(ii-1)*uu(ii)<0.0d0) THEN
    node=node+1
    END IF
    END DO

    RETURN
    !---------------------------------------------------------------------------
    END SUBROUTINE hmn_vkboutwf
    !---------------------------------------------------------------------------
    !
    !
    !---------------------------------------------------------------------------
    SUBROUTINE hmn_lschps(ll,ierr,ee,rr,vv,uu,up,mmax,mch)
    !---------------------------------------------------------------------------
    !!
    !! outward integration of Srcroedinger equation for semi-local pseudopotential
    !! on a logarithmic mesh
    !!
    !
    !nn  principal quantum number (not used actually)
    !ll  angular-momentum quantum number
    !ierr  non-zero return if error
    !ee  bound-state energy, input guess and output calculated value
    !rr  log radial mesh
    !vv  semi-local atomic pseudopotential
    !uu  output radial wave function (*rr)
    !up  d(uu)/dr
    !mmax  size of log grid
    !mch matching mesh point for inward-outward integrations

    IMPLICIT NONE
    INTEGER, PARAMETER :: dp=kind(1.0d0)

    !Input variables
    INTEGER :: mmax,mch
    REAL(DP) :: rr(mmax),vv(mmax)
    INTEGER :: ll

    !Output variables
    REAL(DP) :: uu(mmax),up(mmax)
    REAL(DP) :: ee
    INTEGER :: ierr


    !Local variables
    REAL(DP) :: amesh,al
    REAL(DP) :: aeo, aio, als, cn
    REAL(DP) :: ro
    REAL(DP) :: sls, sn, uout, upout
    INTEGER :: ii, it

    REAL(DP), allocatable :: upp(:),cf(:)

    ALLOCATE(upp(mmax),cf(mmax))


    al = 0.01d0 * dlog(rr(101) / rr(1))
    amesh = dexp(al)

    ierr = 0

    sls=ll*(ll+1)

    ! null arrays to remove leftover garbage

    uu(:)=0.0d0
    up(:)=0.0d0
    upp(:)=0.0d0

    als=al**2

    ! coefficient array for u in differential eq.
    DO ii=1,mmax
    ! cf(ii)=als*sls + 2.0d0*als*(vv(ii)-ee)*rr(ii)**2
    cf(ii)=als*sls + als*(vv(ii)-ee)*rr(ii)**2 ! DR
    END DO

    ! start wavefunction with series

    DO ii=1,4
    uu(ii)=rr(ii)**(ll+1)
    up(ii)=al*(ll+1)*rr(ii)**(ll+1)
    upp(ii)=al*up(ii)+cf(ii)*uu(ii)
    END DO

    ! outward integration using predictor once, corrector
    ! twice

    DO ii=4,mch-1
    uu(ii+1)=uu(ii)+aeo(up,ii)
    up(ii+1)=up(ii)+aeo(upp,ii)
    DO it=1,2
    upp(ii+1)=al*up(ii+1)+cf(ii+1)*uu(ii+1)
    up(ii+1)=up(ii)+aio(upp,ii)
    uu(ii+1)=uu(ii)+aio(up,ii)
    END DO
    END DO

    uout=uu(mch)
    upout=up(mch)

    !perform normalization sum

    ro=rr(1)/dsqrt(amesh)
    sn=ro**(2*ll+3)/(2*ll+3)

    DO ii=1,mch-3
    sn=sn+al*rr(ii)*uu(ii)**2
    END DO

    sn =sn + al*(23.0d0*rr(mch-2)*uu(mch-2)**2 &
    &       + 28.0d0*rr(mch-1)*uu(mch-1)**2 &
    &      +  9.0d0*rr(mch  )*uu(mch  )**2)/24.0d0

    !normalize u

    cn=1.0d0/dsqrt(sn)
    uout=cn*uout
    upout=cn*upout

    DO ii=1,mch
    up(ii)=cn*up(ii)
    uu(ii)=cn*uu(ii)
    END DO
    DO ii=mch+1,mmax
    uu(ii)=0.0d0
    END DO

    DEALLOCATE(upp,cf)

    RETURN
    !---------------------------------------------------------------------------
    END SUBROUTINE hmn_lschps
    !---------------------------------------------------------------------------
    !
    !
    !---------------------------------------------------------------------------
    SUBROUTINE hmn_lschkb(ll,ierr,ee,vkb,rr,vv,uu,up,mmax,mch)
    !---------------------------------------------------------------------------
    !!
    !! outward integration of the inhomogeneous radial Schroedinger equation
    !! on a logarithmic mesh with local potential and one proector term
    !!
    !nn  principal quantum number (not used)
    !ll  angular-momentum quantum number
    !ierr  non-zero return if error
    !ee  bound-state energy, input guess and output calculated value
    !vkb Vanderbilt-Kleinman-bylander projector
    !rr  log radial mesh
    !vv  local pseudopotential
    !uu  output radial wave function (*rr)
    !up  d(uu)/dr
    !zz  atomic number
    !mmax  size of log grid
    !mch matching mesh point for inward-outward integrations

    IMPLICIT NONE
    INTEGER, PARAMETER :: dp=kind(1.0d0)

    !Input variables
    INTEGER :: mmax,mch
    REAL(DP) :: rr(mmax),vv(mmax),vkb(mmax)
    INTEGER :: ll

    !Output variables
    REAL(DP) :: uu(mmax),up(mmax)
    REAL(DP) :: ee
    INTEGER :: ierr


    !Local variables
    REAL(DP) :: amesh,al
    REAL(DP) :: aeo, aio, als
    REAL(DP) :: sls, uout, upout
    REAL(DP) :: akb,ckb
    INTEGER :: ii, it

    REAL(DP), allocatable :: upp(:),cf(:)

    ALLOCATE(upp(mmax),cf(mmax))


    al = 0.01d0 * dlog(rr(101) / rr(1))
    amesh = dexp(al)

    ierr = 0

    sls=ll*(ll+1)

    ! null arrays to remove leftover garbage

    uu(:)=0.0d0
    up(:)=0.0d0
    upp(:)=0.0d0

    als=al**2

    ! coefficient array for uu in differential eq.
    DO ii=1,mmax
    ! cf(ii)=als*sls + 2.0d0*als*(vv(ii)-ee)*rr(ii)**2
    cf(ii)=als*sls + als*(vv(ii)-ee)*rr(ii)**2 ! DR
    END DO

    ! start wavefunction with series based on projector

    !  ckb = 2.0d0*vkb(1)/rr(1)**(ll+1)
    ckb = vkb(1)/rr(1)**(ll+1) ! DR
    akb = ckb/(6.0d0+4.0d0*ll)
    DO ii=1,4
    uu(ii)=akb*rr(ii)**(ll+3)
    up(ii)= al*(ll+3)*uu(ii)
    upp(ii)= als*(ll+3)**2*uu(ii)
    END DO

    ! outward integration using predictor once, corrector
    ! twice

    DO ii=4,mch-1
    uu(ii+1)=uu(ii)+aeo(up,ii)
    up(ii+1)=up(ii)+aeo(upp,ii)
    DO it=1,2
    upp(ii+1)=al*up(ii+1)+cf(ii+1)*uu(ii+1) &
    ! &        + 2.0d0*als*vkb(ii+1)*rr(ii+1)**2
    &        + als*vkb(ii+1)*rr(ii+1)**2 ! DR
    up(ii+1)=up(ii)+aio(upp,ii)
    uu(ii+1)=uu(ii)+aio(up,ii)
    END DO
    END DO

    uout=uu(mch)
    upout=up(mch)

    DEALLOCATE(upp,cf)

      RETURN
      !
    !---------------------------------------------------------------------------
    END SUBROUTINE hmn_lschkb
    !---------------------------------------------------------------------------
    !
    !
    !---------------------------------------------------------------------------
    SUBROUTINE hmn_vpinteg(gg, hh, nn, mm, ss, rr)
    !---------------------------------------------------------------------------
    !!
    !! integrals that go into construction of Vanderbilt separable pseudopotential
    !!
      !
      IMPLICIT NONE
      !
      ! product of functions gg*hh goes like rr**mm at rr -> 0
      ! integral on usual log mesh from 1 to nn
      !
      !Input variables
      REAL(DP) :: gg(nn), hh(nn), rr(nn)
      INTEGER :: nn, mm
      !
      !Output variable
      REAL(DP) :: ss
      !
      !Local variables
      REAL(DP) :: r0, amesh, al
      INTEGER :: ii
      !
      al = 0.01d0 * LOG(rr(101) / rr(1))
      amesh = EXP(al)
      !
      r0 = rr(1) / SQRT(amesh)
      ss = r0**(mm + 1) * (gg(1) * hh(1) / rr(1)**mm) / FLOAT(mm + 1)
      !
      DO ii = 4, nn - 3
      ss =  ss + al * gg(ii) * hh(ii) * rr(ii)
      END DO
      !
      ss=ss + al * (23.d0 * rr(nn - 2) * gg(nn - 2) * hh(nn - 2) + &
        + 28.d0 * rr(nn-1) * gg(nn - 1) * hh(nn - 1) + &
        +  9.d0 * rr(nn) * gg(nn) * hh(nn)) / 24.d0
      !
      RETURN
      !
    !---------------------------------------------------------------------------
    END SUBROUTINE hmn_vpinteg
    !---------------------------------------------------------------------------
    !
    !
  !=============================================================================
  END MODULE hamann
  !=============================================================================
