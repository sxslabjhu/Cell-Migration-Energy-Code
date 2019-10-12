!!-----------------------------------------------------------------
!!
!! 2D code for two phase model for cell migration, by Lingxing Yao
!! 
!! For simulations prepared for article
!! "On the Energy Efficiency of Cell Migration in Diverse Physical Environments"
!! by Yizeng Li, Lingxing Yao, Yoichiro Mori and Sean X. Sun
!!
!!-----------------------------------------------------------------
!!
!! Copyright (c) 2018-2019, Lingxing Yao
!! All rights reserved.
!!
!! Redistribution and use of the code, with or without modifications, are
!! permitted provided that the following conditions are met:
!!
!!    * Redistributions of source code must retain the above copyright notice,
!!      this list of conditions and the following disclaimer.
!!
!!    * Redistributions in binary form must reproduce the above copyright
!!      notice, this list of conditions and the following disclaimer in the
!!      documentation and/or other materials provided with the distribution.
!!
!! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
!! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
!! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
!! ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
!! LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!! CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!! SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!! CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!! POSSIBILITY OF SUCH DAMAGE.

module energy
  use, intrinsic :: iso_c_binding
  use parameters
  use geometry
  use chemical_mod, only : cc_pump
  use network_mod, only : nt_theta, Jactin
  use linsys, only : setbc
  use ibforce, only : calcplateletforce

  implicit none

private

  integer, parameter :: nxb = nx
  integer, parameter :: nyb = ny
!
  integer, dimension(1:nx*ny) :: idx
  double precision, dimension(-1:nx+2,-1:ny+2) :: omega, muc

public :: getGs, getIflow, getEn, getInet, getImem, getJmem

contains
!========================================================================
!
!========================================================================
subroutine getIflow(Iflow, nuf, difC, u, v, cc, iaary, llen)
! in entire domain, use v_c @MAC grid and c @cell centers, simplified complete
  integer, intent(in) :: llen(nent-2)
  type(iapt), dimension(:), intent(in) :: iaary(nent-2)
  double precision, intent(in out):: Iflow
  double precision, intent(in) :: nuf, difC
  double precision, intent(in) :: u(-1:nx+1,-1:ny+1)   ! x velo of solvent, MAC grid
  double precision, intent(in) :: v(-1:nx+1,-1:ny+1)   ! y velo of solvent, MAC grid
  double precision, intent(in) :: cc(-1:nx+2,-1:ny+2)  ! chemical @ cell center
!
  double precision :: tmp, tpc, tpv, ux, uy, vx, vy, gx, gy, tp1, tp2

!!  double precision :: uc(-1:nx+1,-1:ny+1)   ! copy of u
!!  double precision :: vc(-1:nx+1,-1:ny+1)   ! copy of v
  double precision, dimension(-1:nx+1,4) :: uvbc !ubt, utp, vbt, vtp
  integer :: i, j
  type(cc_augvar), pointer :: cp_lp
!
!!  uc = u; vc = v
!!  uvbc = zero;
!!  call setbc(uc,uvbc,0)
!!  call setbc(vc,uvbc,1)
!
  tmp = zero
  do i = 2, nx-2
  do j = 2, ny-2 ! interior points only
    if (id(i,j) .eq. 0) then
      tp1 = (cc(i+1,j)-cc(i-1,j)); tp2 = (cc(i,j+1)-cc(i,j-1));
      tmp = tmp + half*half*(tp1*tp1+tp2*tp2)/hg/hg
    else if (id(i,j) > 0) then
      tmp = tmp
    else if (id(i,j) < 0) then
      tmp = tmp
    endif
  enddo
  enddo
  tpc = tmp*hg*hg*difC
!
  do i = 2, nx-2
  do j = 2, ny-2! interior points only
    ux = (u(i,j)-u(i-1,j)); vy = (v(i,j)-v(i,j-1))
    uy = half*((u(i,j+1)+u(i-1,j+1)) - (u(i,j-1)+u(i-1,j-1)))
    vx = half*((v(i+1,j)+v(i+1,j-1)) - (v(i-1,j)+v(i-1,j-1)))
    tmp= tmp + ux*ux + vy*vy + (uy+vx)*(uy+vx)*half
  enddo
  enddo
  tpv = two*nuf*tmp
!
  Iflow = tpc + tpv
!
  return
end subroutine getIflow
!!
!========================================================================
!
!========================================================================
!!subroutine getInet(Inet, etafl, etast, sig, u,v, iaary, llen)
subroutine getInet(Inet, sig, u,v, iaary, llen)
!I_n evaluated inside membrane only
!
  double precision, intent(out):: Inet
  integer, intent(in) :: llen(nent-2)
  type(iapt), dimension(:), intent(in) :: iaary(nent-2)
  !!double precision, intent(in) :: etafl, etast
  double precision :: etafl, etast
  double precision, intent(in) :: u(-1:nx+1,-1:ny+1)   ! x velocity for solvent
  double precision, intent(in) :: v(-1:nx+1,-1:ny+1)   ! y velocity for solvent
  double precision, intent(in) :: sig(-1:nx+2,-1:ny+2) ! network stress, extended on input
  !double precision, intent(in) :: theta(-1:nx+2,-1:ny+2) !
  !double precision :: theta(-1:nx+2,-1:ny+2) ! actin volume, extended
!
  integer :: i, j
  double precision :: tpx, tpy, theta
  double precision, dimension(-1:nx+1,-1:ny+1) :: uc, vc
  double precision :: tmp, tp1, tp2, coe1, coe2, coe3, tuh
!
  etast = etas; etafl = eta
  tuh = two*hg
  uc(1:nx,1:ny) = half*(u(0:nx-1,1:ny) + u(1:nx,1:ny))
  vc(1:nx,1:ny) = half*(v(1:nx,0:ny-1) + v(1:nx,1:ny))
!
  tmp = zero; tp1 = zero; tpx = zero; tpy = zero
!
  coe1 = one/(etafl+etast); coe2 = etafl*coe1; coe3 = etast*coe1

!!  do i = 2, nx-2
!!  do j = 2, ny-2
!!    if (idf(i,j)
!!    theta(i,j) = nt_theta(sig(i,j)) ! compute theta from sigma
!!  enddo
!!  enddo
!
  do i = 2, nx-2
  do j = 2, ny-2
    !!if (idf(i,j) .eq. 1 .or. id(i,j) .eq. -1) then
    if ( idf(i,j) .eq. 1 ) then
      theta = nt_theta(sig(i,j))
      ! u-component
      tp1 = coe3*theta*uc(i,j)+coe1*(sig(i+1,j)-sig(i-1,j))/tuh
      tp1 = tp1*tp1*etafl
      tp2 = coe2*theta*uc(i,j)-coe1*(sig(i+1,j)-sig(i-1,j))/tuh
      tp2 = tp2*tp2*etast
      tmp = tmp+ (tp1 + tp2)/theta
      ! v-component
      tp1 = coe3*theta*vc(i,j)+coe1*(sig(i,j+1)-sig(i,j-1))/tuh
      tp1 = tp1*tp1*etafl
      tp2 = coe2*theta*vc(i,j)-coe1*(sig(i,j+1)-sig(i,j-1))/tuh
      tp2 = tp2*tp2*etast
      tmp = tmp+ (tp1 + tp2)/theta
    endif
  enddo
  enddo
  Inet = tmp*hg*hg
!
  return
end subroutine getInet
!!
!========================================================================

!========================================================================
subroutine getJnet
! J_n = 0 analytically as it evaluate \grad\cdot v_n

  return
end subroutine getJnet
!========================================================================

!========================================================================
subroutine getEsol

  return
end subroutine getEsol
!========================================================================

!========================================================================
subroutine getEn(En,sig,iflg)
! only in membrane, completed without counting partial computational cells
  integer, intent(in) :: iflg
  double precision :: En, sig(-1:nx+1,-1:ny+1)
!!
  integer :: i, j
  double precision :: tmp, tp1
!
  tmp = zero
  do i = 1, nx
  do j = 1, ny
    if (idf(i,j) > 0) then
      tp1 = nt_theta(sig(i,j))
      tmp = tmp +fen(tp1)
    endif
  enddo
  enddo
  En = tmp*hg*hg
!
  return
end subroutine getEn

!========================================================================

!========================================================================
subroutine getEmem

  return
end subroutine getEmem

!========================================================================
subroutine getGs(gs, cc,iflg)
! complete without counting on the partial computational cells near membrane
  integer :: iflg
  double precision :: gs, cc(-1:nx+1,-1:nx+1)
!
  integer :: i, j
  double precision :: tmp, tp1
!
  tmp = zero
  do i = 2, nx-2
  do j = 2, ny-2
    tp1 = cc(i,j)
    if (tp1 < zero) then
      print *,'concentration is less than 0!'
      stop
    endif
    tmp = tmp +tp1*(dlog(tp1)-one)
    !!tmp = tmp +tp1*((tp1)-one)
    !tmp = tmp +log(tp1)
  enddo
  enddo
  gs = tmp*hg*hg
!
  return
end subroutine getGs

!========================================================================
subroutine getImem(Imem, ckc, fkw, mkc, mkm, iaary, llen,iflg)
!
  integer :: iflg
  double precision, intent(in out):: Imem
  double precision, intent(in) :: ckc, fkw
  double precision, intent(in) :: mkc(7,nring,cent), mkm(7,nring,cent)! for chemical and stress
  type(iapt), dimension(:), intent(in) :: iaary(nent-2)
  integer, intent(in) :: llen(nent-2)
!
  integer :: i, il, info
  double precision :: cj, muj, tmp, tp1, tp2, s0, ds, psi, sig, frc
  double precision, dimension(2*nring) :: f0, cf0
  double precision :: gp(mpts)
!
  
  ds = two*cpi/dble(nring)

  call calcplateletforce(xpt,ypt,f0,cf0,gp,1)
  il = 1
  tmp = zero
  do i = 1, nring
    !!s0 = ds*dble(i-1)
    s0 = mkc(1,i,il)
    tp1= paraval(s0, nring, mkc, 1, 0, 0, info) ! on - side old time chemical
    tp2= paraval(s0, nring, mkc, 1, 1, 0, info) ! on + side old time chemical
    cj = tp1-tp2
    muj= dlog(tp1)-dlog(tp2)
    sig= paraval(s0, nring, mkm, 1, 1, 0, info) ! on - side old time chemical
    frc= f0(2*i-1) *ndv(nring-i+1,1,il) + f0(2*i)*ndv(nring-i+1,2,il)
    psi=(cj+(sig+frc)/rtc)
!
    tmp = tmp+ ckc*cj*muj+fkw*psi*psi
  enddo

  Imem = tmp*ds
!
  return
end subroutine getImem
!========================================================================

!========================================================================
subroutine getJmem(Jmem,mkc,mkm,iaary,llen,iflg,time)
! along membrane
  integer, intent(in) :: iflg
  double precision, intent(in) :: time
  double precision, intent(in out) :: Jmem
  double precision, intent(in) ::  mkc(7,nring,cent) ! hold chemical 
  double precision, intent(in) ::  mkm(7,nring,cent) ! hold stress
  type(iapt), dimension(:), intent(in) :: iaary(nent-2)
  integer, intent(in) :: llen(nent-2)
!
  integer :: i, il, info
  double precision :: tmp, muj, det, jpmp, jact, s0, ds, tp1, tp2

  ds = two*cpi/dble(nring)

  tmp = zero
  il = 1
  do i = 1, nring
    s0 = mkc(1,i,il)
    tp1= paraval(s0, nring, mkc, 1, 0, 0, info) ! on - side chemical
    tp2= paraval(s0, nring, mkc, 1, 1, 0, info) ! on + side chemical
    muj = dlog(tp1)-dlog(tp2)
    jpmp = cc_pump(s0,time)
    jact = Jactin(s0,time)
    tp1= paraval(s0, nring, mkm, 1, 1, 0, info) ! on + side stress
    tp2= nt_theta(tp1)
    det= dfen(tp2)
    tmp = tmp + muj*jpmp+jact*det
  enddo

  Jmem = tmp*ds

!
  return
end subroutine getJmem
!========================================================================

double precision function fen(theta)
! this is e_n(theta)
  double precision, intent(in) :: theta

  fen = half*(theta0/theta+ theta*dlog(theta)/theta0)

end function fen

double precision function dfen(theta)
! this is de_n(theta)/d\theta
  double precision, intent(in) :: theta

  dfen = half*(-theta0/theta/theta+ (dlog(theta)+one)/theta0)

end function dfen


end module energy
