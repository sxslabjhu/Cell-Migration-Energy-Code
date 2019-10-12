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

module network_mod
  use, intrinsic :: iso_c_binding
  use parameters
  use geometry
  use solver_mod
  use chemical_mod, only : vc_cfs

  implicit none

  private

  integer, parameter :: nxb = nx
  integer, parameter :: nyb = ny
!
  double precision, dimension(-1:nx+2,-1:ny+2) :: xavgx, pavg
  double precision :: pnu(1:nx*ny)
  double precision :: nsigma
  integer :: igd(1:nx*ny)
!
  
  public :: nt_init, nt_IBbdy
  public :: nt_sigma, nt_dsigma, nt_theta, Jactin, nt_sga
contains

subroutine nt_init(iaary, llen) ! chemical concentration initialization
  implicit none
  type(iapt), dimension(:) :: iaary(nent-2)
  integer :: llen(nent-2)
!
  double precision :: nmvec(1:nring), npvec(nring)
  type(cc_augvar), pointer :: cp_lp
  double precision :: time = 0., x, y, dx, dy, sigma
  integer :: isel, i, j, info, il
!
  dx = hg; dy = hg

  sigma = nt_sigma(theta0) ! get sigma_n
  print *, sigma, npsi
  !!sigma = half
  nsigma = sigma; 
  nt_pc = sigma
  nt_pn = sigma
!
  do il = 1, nib
    cp_lp => iaary(il)%p
    do i = 1, llen(il)
      cp_lp%nm = sigma !!set initial sigma_n
      cp_lp%np = sigma !!set initial sigma_n
      cp_lp => cp_lp%next
    enddo
    nmvec = sigma
    npvec = sigma
    mkq(1,:,il) = mkn(1,:,il)
    call smoothcrossing(qlen, nib,il,mkq,iaary,llen, nmvec,npvec, info)
  enddo
!!  do il = 1, nib
!!    !!print *, 'ck here', nl
!!  enddo
!!prt  do i = 1, nring
!!prt    dx = paraval(mkq(1,i,nib), nring, mkq, 1, 0, 0, info) ! on - side old time chemical
!!prt    dy = paraval(mkq(1,i,nib), nring, mkq, 1, 1, 0, info) ! on + side old time chemical
!!prt!!    print '(1x, " why ", i5, 20(e14.6,1x))', i, mkq(1,i,nib), dx, dy, npvec(i), nmvec(i), mkn(1,i,nib), npsi, one
!!prt  enddo
!
  return
end subroutine nt_init
!=======================================================================
subroutine f2m(x, M, nl, ux, uy, iaary, llen, vnt, info, dt, mydif)
! f2m subroutine will project value at grid crossing to corresponding
! im or ip index. Note variable x should have exact order as in iaary
! since we will need delta info. The order if "-" first then "+" side
! pointer will continue on the "next" direction
!
  double precision :: dt, mydif
  double precision :: M(1:nx*ny), vnt(1:nx*ny)
  double precision :: ux(1:nx*ny), uy(1:nx*(ny+1))
  type(iapt), dimension(:) :: iaary(nent-2)
  integer :: nl, info, llen(nent-2)
  double precision, dimension(nl) :: x
!
  integer :: i, k, il, ic, im(2), ip(2), gm(2)
  double precision :: tmp, coef, delta, utx, vtx, fm(2), uca, vca, theta, dthe, s0, ets
  type(cc_augvar), pointer :: cp_lp
!
  !!coef =-dt*dif(1)/hg/hg  ! Euler
  ets = one/(eta+etas)
  coef =-ets*dt/hg/hg  ! for network volume
  utx = 0.5d0*dt/hg; vtx = utx
!
  M = 0.
  do il = 1, nib
    cp_lp => iaary(il)%p
    ic = 0
    do i = 1, llen(il)
      delta = (cp_lp%delm); im = cp_lp%im; ip = cp_lp%ip; gm = ip-im; fm = dble(gm)
      !!tmp = 6./(1.+delta)/(2.+delta)*coef
!twoside      ic = ic + 1
!twoside      call sub2ind(nxb,im(1),im(2),k)
!twoside      !!M(im(1),im(2)) = tmp*x(ic)
!twoside      !!M(im(1),im(2)) = M(im(1),im(2))+ coef*cp_lp%st_m(1)*x(ic)
!twoside      !!old M(k) = M(k)+ coef*cp_lp%st_m(1)*x(ic)
!twoside      !!M(k) = M(k)+ (delta*(utx*ux(k)*fm(1)+vtx*uy(k)*fm(2))+coef)*cp_lp%st_m(1)*x(ic)
!twoside      !!M(k) = M(k)+ coef*cp_lp%st_m(1)*cmsi !
!twoside      !!M(k) = M(k)+ (delta*(utx*ux(k)*fm(1)+vtx*uy(k)*fm(2))+coef)*cp_lp%st_m(1)*x(ic) ! here ux is @ cell center
!twoside!!deb      uca = half*(ux(k)+ux(k+1)); vca = half*(uy(k)+uy(k+nx))
!twoside!!deb      utx = dlt/hg/(delta+half);  vtx = dlt/hg/(delta+half)
!twoside!!deb
!twoside!!deb      M(k) = M(k)+ ((utx*uca*fm(1)+vtx*vca*fm(2))+coef*pnu(k)*cp_lp%st_m(1))*x(ic) ! here ux is @ cell center
!twoside      M(k) = zero ! debugging mode

      !!delta = (cp_lp%delp); im = cp_lp%im; ip = cp_lp%ip; gm = im-ip; fm = dble(gm)
      delta = (cp_lp%delp);                               gm = im-ip; fm = dble(gm)
      !!tmp = 6./(1.+delta)/(2.+delta)*coef
      ic = ic + 1
      call sub2ind(nxb,ip(1),ip(2),k)
      !!M(ip(1),ip(2)) = tmp*x(ic)
      !!M(ip(1),ip(2)) = M(ip(1),ip(2))+ coef*cp_lp%st_p(1)*x(ic)
      !!old M(k) = M(k)+ coef*cp_lp%st_p(1)*x(ic)
      !!M(k) = M(k)+ (delta*(utx*ux(k)*fm(1)+vtx*uy(k)*fm(2))+coef)*cp_lp%st_p(1)*x(ic)
      !!M(k) = M(k)+ coef*cp_lp%st_p(1)*cpsi ! 
      !!uca = half*(ux(k)+ux(k+1)); vca = half*(uy(k)+uy(k+nx))
      uca = half*(ux(k)+ux(k+gm(1))); vca = half*(uy(k)+uy(k+gm(2)*nx))
      utx = dt/hg/(delta+half);  vtx = dt/hg/(delta+half)
      M(k) = M(k)+ ((utx*uca*fm(1)+vtx*vca*fm(2))+coef*pnu(k)*cp_lp%st_p(1))*x(ic) ! here ux is @ cell center
      !!M(k) = M(k)+ (delta*(utx*uca*fm(1)+vtx*vca*fm(2))+coef*pnu(k)*cp_lp%st_p(1))*x(ic) ! here ux is @ cell center
      !!M(k) = M(k)+ (                           coef*pnu(k)*cp_lp%st_p(1))*x(ic) ! here ux is @ cell center
!
      cp_lp => cp_lp%next
    enddo
  enddo
  if (ic .ne. nl) then
    print *, 'error in data processing between linked list and array'
    stop
  endif
!R  print *, 'IC= ', ic
!
  return
end subroutine f2m

!=======================================================================
subroutine EM(Y, x, nl, iaary, llen, info, mydif)
! Y is the input vector defined for chemicals @ cell centers
!
  !!double precision :: Y(-1:nxb+2,-1:nyb+2)  ! output vector
  double precision, intent(in out) :: Y(1:nx*ny)  ! output vector
  type(iapt), dimension(:) :: iaary(nent-2)
  integer :: nl, info, llen(nent-2)
  double precision, dimension(nl) :: x
  double precision, intent(in) :: mydif
!
  integer :: il, i, j, k, ic, im(2), ip(2), g_ali(4,2), g_off(4,2)
  integer :: cpl, gpl, opl
  double precision :: tmp, coef, delta, c_ali(4), c_off(4), s0, theta, dthe, ets
  type(cc_augvar), pointer :: cp_lp
!
  ets = one/(eta+etas)
  do il = 1, nib
    cp_lp => iaary(il)%p
    ic = 0
    do i = 1, llen(il)
      !!g_ali = cp_lp%i_ali_grd_m; c_ali = cp_lp%c_ali_grd_m; 
      !!g_off = cp_lp%i_off_grd_m; c_off = cp_lp%c_off_grd_m; 
      !!cpl   = cp_lp%iml; gpl = g_ali(4,1); opl = g_off(4,1)
!
      !!tmp = 0.
      !!!!do j = 1, cpl
      !!do j = 1, gpl ! note c_ali(j+1) = 0 if g_ali(j,:) is not a stencil pt
      !!  call sub2ind(nxb,g_ali(j,1),g_ali(j,2),k)
      !!  tmp = tmp + Y(k)*c_ali(j+1)
      !!  !!tmp = tmp + Y(g_ali(j,1),g_ali(j,2))*c_ali(j+1)
      !!enddo
      !!!!do j = 1, g_off(4,2)
      !!do j = 1, opl
      !!  call sub2ind(nxb,g_off(j,1),g_off(j,2),k)
      !!  tmp = tmp + Y(k)*c_off(j+1)
      !!  !!tmp = tmp + Y(g_off(j,1),g_off(j,2))*c_off(j+1)
      !!enddo
!!twoside      ic = ic + 1
!!twoside      !!x(ic) = tmp
!!twoside      x(ic) = 0.d0 ! debugging mode
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
      g_ali = 0; g_off = 0; c_ali = 0.; c_off = 0.;
      g_ali = cp_lp%i_ali_grd_p; c_ali = cp_lp%c_ali_grd_p; 
      g_off = cp_lp%i_off_grd_p; c_off = cp_lp%c_off_grd_p; 
      !!cpl   = cp_lp%ipl;
      cpl   = cp_lp%ipl; gpl = g_ali(4,1); opl = g_off(4,1)

      tmp = 0.
      !!do j = 1, cpl
      do j = 1, gpl
        call sub2ind(nxb,g_ali(j,1),g_ali(j,2),k)
        tmp = tmp + Y(k)*c_ali(j+1)
        !!tmp = tmp + Y(g_ali(j,1),g_ali(j,2))*c_ali(j+1)
      enddo
      !!do j = 1, g_off(4,2)
      do j = 1, opl
        call sub2ind(nxb,g_off(j,1),g_off(j,2),k)
        tmp = tmp + Y(k)*c_off(j+1)
        !!tmp = tmp + Y(g_off(j,1),g_off(j,2))*c_off(j+1)
      enddo
      ic = ic + 1
      x(ic) = tmp*rho0*ets
!
      cp_lp => cp_lp%next
    enddo
  enddo
!    
  if (ic .ne. nl) then
    print *, 'error in data processing between linked list and array'
    stop
  endif
!R  print *, 'IC= ', ic

!
  return
end subroutine EM
!
!=======================================================================
subroutine Qmat(x, vmx, nl, iaary, llen, info, time, mydif)
! return I*x on output and stowd in vmx
  type(iapt), dimension(:) :: iaary(nent-2)
  integer :: nl, info, llen(nent-2)
  double precision, dimension(nl) :: vmx, x
  double precision, intent(in) :: time, mydif
!
  integer :: il, i, j, ic, im(2), ip(2), g_ali(4,2), g_off(4,2), cpl
  double precision :: tmp, tp1, coef, delta, c_ali(4), c_off(4), c_alim(4), c_offm(4)
  double precision :: xt, yt, tp0, theta, dthe, s0, ets
  type(cc_augvar), pointer :: cp_lp
!
  vmx = zero
  ets = one/(eta+etas)
  do il = 1, nib
    cp_lp => iaary(il)%p
    ic = 0
    do i = 1, llen(il)
      !!tp0 = cc_gate(cp_lp%s,time)/cp_lp%dxs ! pump is concentration dependent
      !!tp0 = 0. ! constpmp
      !!s0=cp_lp%s
      !!theta = paraval(s0, nring, mkq, il, 1, 0, info) ! on + side old network volume
      !!dthe = nt_dsigma(theta)

!!twoside      ic = ic + 1
!!twoside      !!debc_ali = cp_lp%c_ali_grd_m; c_off = cp_lp%c_off_grd_m
!!twoside      !!deb!!prt print '(1x, "- side ", 1024(e16.8,1x))', c_ali, c_off, tmp, x(ic), vmx(ic)
!!twoside      !!deb!oldtmp =-kc(il) + (c_off(1) + c_ali(1)) + cp_lp%ncm !should be same in solver for prem 
!!twoside      !!deb!oldvmx(ic) = tmp*x(ic) + (kc(il)-cp_lp%ncm)*x(ic+1) +max(tp0,0.)*x(ic+1) + min(tp0,0.)*x(ic)
!!twoside      !!debtp1 =  (c_off(1) + c_ali(1)) ! be same in solver for prem 
!!twoside      !!debvmx(ic) = tp1*x(ic) 
!!twoside      vmx(ic) = x(ic) !debugging mode

      ic = ic + 1
      c_ali = cp_lp%c_ali_grd_p; c_off = cp_lp%c_off_grd_p
      !!prt print '(1x, "+ side ", 1024(e16.8,1x))', c_ali, c_off, tmp, x(ic), vmx(ic)
      !!print '(1x, " coefficient ", 1024(e16.8,1x))', tmp, tp1
      !oldtp1 = kc(il) + (c_off(1) + c_ali(1)) - cp_lp%ncp !+tp0!should be same in solver for prem 
      !oldvmx(ic) = tp1*x(ic) -(kc(il)-cp_lp%ncp)*x(ic-1)+min(tp0,0.)*x(ic-1) + max(tp0,0.)*x(ic)
      tp1 = (c_off(1)+c_ali(1))*rho0*ets + rtr*cp_lp%udn !+tp0!should be same in solver for prem 
      vmx(ic) = tp1*x(ic)
!
      cp_lp => cp_lp%next
    enddo
  enddo
  !!print *, 'finish one'
!
  return
end subroutine Qmat
!
!=======================================================================
!=======================================================================
!=======================================================================
subroutine nmatv(xx, Ax, ux, uy, iaary, llen, time, info, dt, mydif)
! here the vector is reshaped to 1d, in the order of storing each 
! row, which is
! {1(1:...:nxb) ... j(1:.i.:nxb) ... nyb(1:...:nxb)}
  integer :: llen(nent-2), info, nl
  double precision :: time, dt, mydif
  double precision, dimension(1:nx*ny) :: x, xx, Ax
  double precision :: ux(1:nx*ny), uy(1:nx*(ny+1))
  type(iapt), dimension(:) :: iaary(nent-2)
!
  integer :: il, i, j, k, l, ic, g_ali(4,2), g_off(4,2), cpl
  integer :: ik, jk, jm(2), jp(2), is(2), jc(2), iloc(4,2)
  integer :: imv(-1:nx+2,-1:ny+2)
  double precision :: tmp, coef, con, cox, coy, delta, c_ali(4), c_off(4), st(4)
  type(cc_augvar), pointer :: cp_lp
  logical :: ichk
!
  double precision :: x1, y1, utx, vtx, uca, vca, fm(2), tp1, tp2, ets
  integer :: gm(2), km(2)
! 
  ets = one/(eta+etas)
  con = ets*dt/hg/hg; cox = con; coy = con;
  utx = 0.5d0*dt/hg; vtx = utx
  Ax = zero
  x = xx
!
!========================================================================
  !Ax(cip)=(1.+4.*con+utx*(ux(cip+1)-ux(cip))+vtx*(uy(cip+nx)-uy(cip)))*x(cip)&
  !    + ( utx*ux(cip+1 )-cox)*x(cip+1 ) + (-utx*ux(cip)-cox)*x(cip-1 )  &
  !    + ( vtx*uy(cip+nx)-coy)*x(cip+nx) + (-vtx*uy(cip)-coy)*x(cip-nx)
  Ax(cip)=(1.d0+4.d0*con*pnu(cip)+utx*(ux(cip+1)-ux(cip))+vtx*(uy(cip+nx)-uy(cip)))*x(cip)&
    + ( utx*ux(cip+1 )-cox*pnu(cip))*x(cip+1 ) + (-utx*ux(cip)-cox*pnu(cip))*x(cip-1 )&
    + ( vtx*uy(cip+nx)-coy*pnu(cip))*x(cip+nx) + (-vtx*uy(cip)-coy*pnu(cip))*x(cip-nx)

!========================================================================
! bottom: Dirichlet BC (need change corresponding RHS part)
  !Ax(cbt)=(1.+4.*con+utx*(ux(cbt+1)-ux(cbt))+vtx*(uy(cbt+nx)-uy(cbt)))*x(cbt)&
  !  + ( utx*ux(cbt+1 )-cox)*x(cbt+1 )  + (-utx*ux(cbt)-cox)*x(cbt-1)      &
  !  + ( vtx*uy(cbt+nx)-coy)*x(cbt+nx)  - (-vtx*uy(cbt)-coy)*x(cbt  )
  Ax(cbt)=(1.+4.*con*pnu(cbt)+utx*(ux(cbt+1)-ux(cbt))+vtx*(uy(cbt+nx)-uy(cbt)))*x(cbt)&
    + ( utx*ux(cbt+1 )-cox*pnu(cbt))*x(cbt+1 )  + (-utx*ux(cbt)-cox*pnu(cbt))*x(cbt-1)&
    + ( vtx*uy(cbt+nx)-coy*pnu(cbt))*x(cbt+nx)  - (-vtx*uy(cbt)-coy*pnu(cbt))*x(cbt  )
!========================================================================
! top Dirichlet BC (need change corresponding RHS part) 
  !Ax(ctp)=(1.+4.*con+utx*(ux(ctp+1)-ux(ctp))+vtx*(uy(ctp+nx)-uy(ctp)))*x(ctp)&
  !    + ( utx*ux(ctp+1 )-cox)*x(ctp+1)  + (-utx*ux(ctp)-cox)*x(ctp-1)   &
  !    - ( vtx*uy(ctp+nx)-coy)*x(ctp  )  + (-vtx*uy(ctp)-coy)*x(ctp-nx)
  Ax(ctp)=(1.+4.*con*pnu(ctp)+utx*(ux(ctp+1)-ux(ctp))+vtx*(uy(ctp+nx)-uy(ctp)))*x(ctp)&
    + ( utx*ux(ctp+1 )-cox*pnu(ctp))*x(ctp+1)  + (-utx*ux(ctp)-cox*pnu(ctp))*x(ctp- 1)&
    - ( vtx*uy(ctp+nx)-coy*pnu(ctp))*x(ctp  )  + (-vtx*uy(ctp)-coy*pnu(ctp))*x(ctp-nx)
!========================================================================
! left periodic
  !Ax(clf)=(1.+4.*con+utx*(ux(clf+1)-ux(clf))+vtx*(uy(clf+nx)-uy(clf)))*x(clf)&
  !    + ( utx*ux(clf+1 )-cox)*x(clf+1 ) + (-utx*ux(clf)-cox)*x(crt   )&
  !    + ( vtx*uy(clf+nx)-coy)*x(clf+nx) + (-vtx*uy(clf)-coy)*x(clf-nx)
  Ax(clf)=(1.+4.*con*pnu(clf)+utx*(ux(clf+1)-ux(clf))+vtx*(uy(clf+nx)-uy(clf)))*x(clf)&
    + ( utx*ux(clf+1 )-cox*pnu(clf))*x(clf+1 ) + (-utx*ux(clf)-cox*pnu(clf))*x(clf   )&
    + ( vtx*uy(clf+nx)-coy*pnu(clf))*x(clf+nx) + (-vtx*uy(clf)-coy*pnu(clf))*x(clf-nx)
!========================================================================
! right periodic
  !Ax(crt)=(1.+4.*con+utx*(ux(clf  )-ux(crt))+vtx*(uy(crt+nx)-uy(crt)))*x(crt)&
  !    + ( utx*ux(clf   )-cox)*x(clf   ) + (-utx*ux(crt)-cox)*x(crt-1 )  &
  !    + ( vtx*uy(crt+nx)-coy)*x(crt+nx) + (-vtx*uy(crt)-coy)*x(crt-nx)
  Ax(crt)=(1.+4.*con*pnu(crt)+utx*(ux(clf  )-ux(crt))+vtx*(uy(crt+nx)-uy(crt)))*x(crt)&
    + ( utx*ux(clf   )-cox*pnu(crt))*x(crt   ) + (-utx*ux(crt)-cox*pnu(crt))*x(crt-1 )&
    + ( vtx*uy(crt+nx)-coy*pnu(crt))*x(crt+nx) + (-vtx*uy(crt)-coy*pnu(crt))*x(crt-nx)
                                
!========================================================================
! bottom left corner: Dirichlet BC (need change corresponding RHS part)
  !Ax(cbl)=(1.+4.*con+utx*(ux(cbl+1)-ux(cbl))+vtx*(uy(cbl+nx)-uy(cbl)))*x(cbl)&
  !  + ( utx*ux(cbl+1 )-cox)*x(cbl+1 ) + (-utx*ux(cbl)-cox)*x(cbr   )&
  !  + ( vtx*uy(cbl+nx)-coy)*x(cbl+nx) - (-vtx*uy(cbl)-coy)*x(cbl   )
  Ax(cbl)=(1.+4.*con*pnu(cbl)+utx*(ux(cbl+1)-ux(cbl))+vtx*(uy(cbl+nx)-uy(cbl)))*x(cbl)&
    + ( utx*ux(cbl+1 )-cox*pnu(cbl))*x(cbl+1 ) + (-utx*ux(cbl)-cox*pnu(cbl))*x(cbr   )&
    + ( vtx*uy(cbl+nx)-coy*pnu(cbl))*x(cbl+nx) - (-vtx*uy(cbl)-coy*pnu(cbl))*x(cbl   )
!========================================================================
!!!right, bot: Dirichlet
  !Ax(cbr)=(1.+4.*con+utx*(ux(cbl  )-ux(cbr))+vtx*(uy(cbr+nx)-uy(cbr)))*x(cbr)&
  !    + ( utx*ux(cbl   )-cox)*x(cbl   ) + (-utx*ux(cbr)-cox)*x(cbr-1 )  &
  !    + ( vtx*uy(cbr+nx)-coy)*x(cbr+nx) - (-vtx*uy(cbr)-coy)*x(cbr   )
  Ax(cbr)=(1.+4.*con*pnu(cbr)+utx*(ux(cbl  )-ux(cbr))+vtx*(uy(cbr+nx)-uy(cbr)))*x(cbr)&
    + ( utx*ux(cbl   )-cox*pnu(cbr))*x(cbl   ) + (-utx*ux(cbr)-cox*pnu(cbr))*x(cbr-1 )&
    + ( vtx*uy(cbr+nx)-coy*pnu(cbr))*x(cbr+nx) - (-vtx*uy(cbr)-coy*pnu(cbr))*x(cbr   )
!========================================================================
!! left, top: Dirichlet
  !Ax(ctl)=(1.+4.*con+utx*(ux(ctl+1)-ux(ctl))+vtx*(uy(ctl+nx)-uy(ctl)))*x(ctl)&
  !    + ( utx*ux(ctl+1 )-cox)*x(ctl+1 ) + (-utx*ux(ctl)-cox)*x(ctr   )  &
  !    - ( vtx*uy(ctl+nx)-coy)*x(ctl   ) + (-vtx*uy(ctl)-coy)*x(ctl-nx)
  Ax(ctl)=(1.+4.*con*pnu(ctl)+utx*(ux(ctl+1)-ux(ctl))+vtx*(uy(ctl+nx)-uy(ctl)))*x(ctl)&
    + ( utx*ux(ctl+1 )-cox*pnu(ctl))*x(ctl+1 ) + (-utx*ux(ctl)-cox*pnu(ctl))*x(ctr   )&
    - ( vtx*uy(ctl+nx)-coy*pnu(ctl))*x(ctl   ) + (-vtx*uy(ctl)-coy*pnu(ctl))*x(ctl-nx)
!========================================================================
!!!right, top: Dirichlet
  !Ax(ctr)=(1.+4.*con+utx*(ux(ctl  )-ux(ctr))+vtx*(uy(ctr+nx)-uy(ctr)))*x(ctr)&
  !    + ( utx*ux(ctl   )-cox)*x(ctl   ) + (-utx*ux(ctr)-cox)*x(ctr-1 )  &
  !    - ( vtx*uy(ctr+nx)-coy)*x(ctr   ) + (-vtx*uy(ctr)-coy)*x(ctr-nx)
  Ax(ctr)=(1.+4.*con*pnu(ctr)+utx*(ux(ctl  )-ux(ctr))+vtx*(uy(ctr+nx)-uy(ctr)))*x(ctr)&
    + ( utx*ux(ctl   )-cox*pnu(ctr))*x(ctl   ) + (-utx*ux(ctr)-cox*pnu(ctr))*x(ctr-1 )&
    - ( vtx*uy(ctr+nx)-coy*pnu(ctr))*x(ctr   ) + (-vtx*uy(ctr)-coy*pnu(ctr))*x(ctr-nx)
!=======================================================================
  imv = 0
  ichk = .false.
  do il = 1, nib
    cp_lp => iaary(il)%p
    do i = 1, llen(il) ! assume irregular cells are 1 grid away from ctb/cbt, clf/crt
    !!do i = 1, 2
!!twoside      jm = 0; jp = 0; iloc = 0; cpl = 0; st = 0.; fm = 0.
!!twoside      jm = cp_lp%im; jp = cp_lp%ip; cpl = cp_lp%iml; gm = jp-jm; fm = dble(gm)
!!twoside      delta = cp_lp%delm; iloc = cp_lp%ilocm; st = cp_lp%st_m
!!twoside      !!prt print '(1x,1024(e16.8,1x))', cp_lp%st_m
!!twoside      call sub2ind(nx, jm(1), jm(2), k) 
!!twoside      ic = k
!!twoside      call sub2ind(nx, jp(1), jp(2), j) 
!!twoside      !!prt print *, 'Ax(ic) ', Ax(ic), ' x(ic) ', x(ic), ic, con
!!twoside      !!diffusion only Ax(k) = Ax(k) + con*x(j) ! add contribution from ghost point back
!!twoside!========================================================================
!!twoside!!ref  Ax(cip)=(1.+4.*con+utx*(ux(cip+1)-ux(cip))+vtx*(uy(cip+nx)-uy(cip)))*x(cip)&
!!twoside!!ref      + ( utx*ux(cip+1 )-cox)*x(cip+1 ) + (-utx*ux(cip)-cox)*x(cip-1 )  &
!!twoside!!ref      + ( vtx*uy(cip+nx)-coy)*x(cip+nx) + (-vtx*uy(cip)-coy)*x(cip-nx)
!!twoside!========================================================================
!!twoside      !!Ax(k) = Ax(k) - (utx*ux(k+ik)*fm(1)+vtx*uy(k+jk)*fm(2)-con)*x(j) ! add values from ghost point back
!!twoside      !!tp1 = utx*(ux(ic+1 )*(x(ic+1 )+x(ic))-ux(ic)*(x(ic)+x(ic-1 )))
!!twoside      !!tp2 = vtx*(uy(ic+nx)*(x(ic+nx)+x(ic))-uy(ic)*(x(ic)+x(ic-nx)))
!!twoside      !!Ax(k) = Ax(k) - ((tp1*abs(fm(1))+tp2*abs(fm(2)))-con*pnu(k)*x(j)) ! add values from ghost point back
!!twoside      !!Ax(k) = Ax(k) - (utx*ux(k+ik)*fm(1)+vtx*uy(k+jk)*fm(2)-con)*x(j) !new, add values from ghost point back
!!twoside      !!Ax(k) = Ax(k) - (utx*ux(k)*fm(1)+vtx*uy(k)*fm(2)-con)*x(j) ! add values from ghost point back
!!twoside      !!prt print *, 'Ax(ic) ', Ax(ic), ' diff', con*x(j), x(j), j
!!twoside
!!twoside      tmp = zero
!!twoside      Ax(ic) = Ax(ic) + tmp
      !!prt print *, 'Ax(ic) ', Ax(ic), ' x(ic) ', x(ic), ic
!
      jm = 0; jp = 0; iloc = 0; cpl = 0; st = 0.; fm = 0.
      jm = cp_lp%im; jp = cp_lp%ip; cpl = cp_lp%ipl; gm = jm-jp; fm = dble(gm)
      delta = cp_lp%delp; iloc = cp_lp%ilocp; st = cp_lp%st_p
      call sub2ind(nxb, jp(1), jp(2), k) 
      ic = k
      call sub2ind(nxb, jm(1), jm(2), j) ! on the "-" side !!
      !!diffusion only Ax(k) = Ax(k) + con*x(j) ! add contribution from ghost point back
      tp1 = utx*(ux(ic+1 )*(x(ic+1 )+x(ic))-ux(ic)*(x(ic)+x(ic-1 )))
      tp2 = vtx*(uy(ic+nx)*(x(ic+nx)+x(ic))-uy(ic)*(x(ic)+x(ic-nx)))
      Ax(k) = Ax(k) - ((tp1*abs(fm(1))+tp2*abs(fm(2)))-con*pnu(k)*x(j)) ! add values from ghost point back
      !!Ax(k) = 
      !!oldAx(k) = Ax(k) - (utx*ux(k)*fm(1)+vtx*uy(k)*fm(2)-con)*x(j) ! add values from ghost point back
      tmp= 0.
      do j = 1, cpl
        call sub2ind(nxb, iloc(j,1), iloc(j,2), k) 
        tmp = tmp + x(k)*st(j+1)
      enddo
      !!diffusion only tmp =-tmp*con
      !!tmp =tmp*(utx*ux(ic+ik)*fm(1)+vtx*uy(ic+jk)*fm(2)-con) ! conservative form
      tmp =tmp*(-con*pnu(ic))
      !!uca = half*(ux(ic)+ux(ic+1)); vca = half*(uy(ic)+uy(ic+nx))
      uca = half*(ux(ic)+ux(ic+gm(1))); vca = half*(uy(ic)+uy(ic+gm(2)*nx))
      tp1 = -fm(1)*uca*half*(x(ic-gm(1))+x(ic))
      tp2 = -fm(2)*vca*half*(x(ic-gm(2)*nx)+x(ic))
      tp1 = tp1*dt/hg/(delta+half)
      tp2 = tp2*dt/hg/(delta+half)
      tmp = tmp + tp1 + tp2
      !!if (abs(idn(jp(1),jp(2))).eq. 3) then !attention here
      if ((ivf(jp(1),jp(2))).eq. 1.and. imv(jp(1),jp(2)).eq. 0) then !attention here
!!deb        print '("pnu(ic)", e12.5,1x,100(i4,1x))', nt_pc(jp(1),jp(2)), jp, jm, idn(jp(1),jp(2)), &
!!deb           idf(jp(1),jp(2)), imv(jp(1),jp(2))
!   print *, maxval(pnu), minval(pnu)
        if (cp_lp%pi_st .eq. 1) then
          call sub2ind(nxb, cp_lp%p_is(1), cp_lp%p_is(2), k)
          tmp = tmp + x(k)*cp_lp%stp_is(2)
          tmp = tmp + x(ic)*cp_lp%stp_is(1)
        endif
        if (cp_lp%pj_st .eq. 1) then
          call sub2ind(nxb, cp_lp%p_js(1), cp_lp%p_js(2), k)
          tmp = tmp + x(k)*cp_lp%stp_js(2)
          tmp = tmp + x(ic)*cp_lp%stp_js(1)
        endif
!!	ichk = .true.
        imv(jp(1),jp(2)) = 1
!!prt        print '("after  ", e12.5,1x,100(i4,1x))', pnu(ic), jp, jm, idn(jp(1),jp(2)), &
!!prt           idf(jp(1),jp(2)), imv(jp(1),jp(2))
      endif
     

      Ax(ic) = Ax(ic) + tmp
!
      cp_lp => cp_lp%next
    enddo
  enddo
!
  Ax = Ax*dble(igd) + (one-dble(igd))*xx ! debugging mode
!!  Ax = x
!
!========================================================================
!
  return
end subroutine nmatv
!
!=======================================================================
subroutine amatv(x, Ax, nl, ux, uy, iaary, llen, rhsv, time, info, dt, mydif)
! the whole vector has length nxb*nyb+nl, will be output
! the vector x and Ax have length nxb*nyb+nl
!
  integer :: llen(nent-2), info, nl
  double precision, dimension(1:nx*ny) :: rhsv
  double precision :: ux(1:nx*ny), uy(1:nx*(ny+1))
  double precision :: time, dt, mydif
  double precision, dimension(1:nx*ny+nl) :: x, Ax
  double precision, dimension(-1:nxb+2,-1:nyb+2,2) :: vmat
  type(iapt), dimension(:) :: iaary(nent-2)
!
  integer :: i, j, il, ierr, lenj, lenw, ic
  double precision, dimension(1:nxb*ny) :: cg, cgv, pca
  double precision, dimension(1:nl) :: ca, ecg, qca
!
  lenj = nx*ny; lenw = lenj + nl
  cg = x(1:lenj); ca = x(lenj+1:lenw)
!!  print *, 'input ', maxval(cg), maxval(ca)
  call nmatv(cg, cgv, ux, uy, iaary, llen, time, info, dt, mydif)
!!  print *, 'one substep ', maxval(cg), maxval(ca)
!!!=======================================================================
!!  pca = 0.
  call f2m(ca, pca, nl, ux, uy, iaary, llen, rhsv, info, dt, mydif)! pca(1:nxb*nyb): output; ca(nl): input
  Ax(1:lenj) = cgv+pca
!!  print *, 'how large', maxval(pca), maxloc(pca)
!!  Ax(1:lenj) =pca
!!!=======================================================================
  call EM(cg, ecg, nl, iaary, llen, info, mydif)! cg(1:nxb*nyb); ecg(nl): output
!!  print *, 'after EM', maxval(ecg), maxval(cg)
  call Qmat(ca,qca, nl, iaary, llen, info, time,mydif)! qca(nl): output; ca(nl): input
!!  print *, 'after Qmat', maxval(qca), maxval(ca)
! use the 4 above to get
  Ax(lenj+1:lenw) =ecg+qca
!!!=======================================================================
!=======================================================================
!
  return
end subroutine amatv

!=======================================================================
!
!!oldsubroutine nt_solver(rg,ra, cg, rb, nl, vmat, iaary, llen,info, isel,time,dt,mydif)
subroutine nt_solver(rg,ra, cg, rb, nl, iaary, llen,info, isel,time,dt,mydif)
  implicit none
! solve the linear system for A*xc = rs, where xc is chemical variable only
! vmat is for velocity vector (it should be umat & vmat)
  !double precision, dimension(nxb*nyb) :: xc, rs, vmat
  double precision :: time, dt, mydif
!!  double precision :: vmat(-1:nxb+2,-1:nyb+2,2)
  double precision, dimension(-1:nx+2,-1:ny+2) :: rg, cg
  double precision, dimension(1:nl) :: ra, rb
  integer :: nl, info, llen(nent-2), isel
  type(iapt), dimension(:) :: iaary(nent-2), itary(nent-2)
!
  integer, parameter :: mm=50
  integer :: ierr, ipl, ip(2), im(2), iloc(4,2), ict
  integer :: i, j, k, l, il, iter, it, mit, infon, mi, j1, lenj, lenw
  double precision :: s0, theta, dthe, tp1, tp2, delt, st_p(4)
  double precision :: con, cox, coy, tmp, normb, dr, omega, beta, eps, resid
  double precision :: ser0, ser1, serr, dr1, hit, st, error, tol, dx, dy
  double precision, dimension(nx*ny+nl) :: s, vk, x11, bf, r, yy, w, xx
!!  double precision, dimension(nx*ny+nl) :: ph, vv, rt, sv, pv, tv
!
  double precision :: vv(nxb*ny+nl,mm),  hj(nxb*ny+nl,2), hq(nxb*ny+nl,mm)
!========================================================================
  double precision :: rgg(nx*ny), rff(nx*ny)!!, idz(nxb*nyb+nl)
  double precision :: utx, vtx, ctop, cbot, etap, ets
  double precision :: ux(1:nx*ny), uy(1:nx*(ny+1))
!!  integer, dimension(nx-2) :: cbt, ctp
!========================================================================
! 
  !!double precision, dimension(1:nxb*nyb+nl) :: rhs, xc
  double precision, dimension(:), allocatable :: rhsv, xcv, prem
!
  type(cc_augvar), pointer :: cp_lp
  !!character(40) ibfile  , ffile
!
  !!print *, 'In Solver....'
  !!nl = 0
  !!do il = 1, nib
  !!  nl = nl + 2*llen(il) 
  !!enddo
  !!stop
!
  !!prt print *, 'Are we here? even'
 !!dt ; ! now dt is an input
  dx = hg; dy = hg
  mit = nx*ny
  l   = mit; lenj = nx*ny; lenw = lenj + nl
  etap= eta/(eta+etas)


  allocate(rhsv(lenw), xcv(lenw), prem(lenw), stat=ierr)
  if (ierr.ne. 0) then
    print *, 'error in allocating variable in solver! stop'
    stop
  endif

  !!deb print *, 'Are we here?'

!========================================================================
! set RHS for Dirichlet BC 
!========================================================================
  !!cox = dif(1)*dlt/dx/dx; coy = dif(1)*dlt/dy/dy ! backward Euler
  con = dt/dx/dy; cox = con; coy = con 
  utx = 0.5d0*dt/hg; vtx = utx
!
  ux(1:nx*ny) = reshape(unc(0:nx-1,  1:ny  ), (/ nx*ny /)) 
  uy(1:nx*(ny+1)) = reshape(vnc(1:nx,  0:ny  ), (/ nx*(ny+1) /))
  ux = ux*etap
  uy = uy*etap
!========================================================================
!  setting all - side pts (including auxiliary variables) to 0 in index idz
!========================================================================
!========================================================================
! modify top/bottom & 4 corners for Dirichlet BCs
!========================================================================
  rgg = zero
  !!ctop = npsi; cbot = npsi;
!!!!  ctop = zero; cbot = zero;
!!!!  rgg(cbt)=-2.*cbot*(-vtx*uy(cbt)-coy)
!!!!  rgg(cbl)=-2.*cbot*(-vtx*uy(cbl)-coy)
!!!!  rgg(cbr)=-2.*cbot*(-vtx*uy(cbr)-coy)
!!!!  rgg(ctp)=-2.*ctop*( vtx*uy(ctp+nx)-coy) ! add nx to reach the top for
!!!!  rgg(ctl)=-2.*ctop*( vtx*uy(ctl+nx)-coy)
!!!!  rgg(ctr)=-2.*ctop*( vtx*uy(ctr+nx)-coy)
!========================================================================
  rhsv(1:lenj) = reshape(rg(1:nx,1:ny), (/ nx*ny /))
  rhsv(1:lenj) = rhsv(1:lenj) + rgg(1:lenj)
  rhsv(lenj+1:lenw) = ra(1:nl)
  rgg = rhsv(1:lenj)
!
!!prt  print '("max rhs ", 10(e14.6,1x))', maxval(abs(rhsv)), maxval(abs(nt_pc)), &
!!prt    maxval(abs(ra))
  normb=dsqrt(dot_product(rhsv,rhsv))
!!prt  print *, 'RHS norm is ', normb, lenj
!!  normb=sqrt(innerless(lenw,rhsv,rhsv,idz))
!!  print *, 'RHS norm is ', normb
!
  !!eps = 1.d-11*normb/2.**(2.*dble(l2ny-5)); resid=eps; tol = eps;
  eps = 1.d-10*normb/dble(nring/80); resid=eps; tol = eps;
  if (abs(normb)< eps) then
    xx=0.d0
    iter = 1
    resid = 0.
    cg = rg
    rb = ra
    print *, 'NO GMRES, iter', iter, normb
    return
  endif
!========================================================================
!========================================================================
!=======================================================================
! construct Jacobian predconditioner for GMRS
!=======================================================================
  tp2 = 1.d0/(1.d0+2.d0*(cox+coy))
!
!  prem = 1.
  !!prem(1:lenj) = tp2
  prem(1:lenj) = one/(one+two*(pnu*cox+pnu*coy))

  prem(1:lenj) = prem(1:lenj)*dble(igd)+(one-dble(igd))
!!! 
  it = 0
  do il = 1, nib
    cp_lp => iaary(il)%p
    do i = 1, llen(il)
      s0 = cp_lp%s
!!      theta = paraval(s0, nring, mkq,il,1,0,info);
!!twoside      it = it + 1
!!twoside!!      call sub2ind(nxb, cp_lp%im(1), cp_lp%im(2), k)
!!twoside!!      tp2 =(cp_lp%c_off_grd_m(1)+cp_lp%c_ali_grd_m(1))
!!twoside
!!twoside      !prem(it+lenj) = 1./tp2
!!twoside      prem(it+lenj) = one
!
      it = it + 1
      call sub2ind(nxb, cp_lp%ip(1), cp_lp%ip(2), k)
      tp2 = (cp_lp%c_off_grd_p(1)+cp_lp%c_ali_grd_p(1))*rho0/(eta+etas) +rtr*cp_lp%udn !
      prem(it+lenj) = 1.d0/tp2
!
      cp_lp => cp_lp%next
    enddo
  enddo
!
!=======================================================================
! GMRES
!=======================================================================
  !xc = 1.+0.001*tp2
  xcv = .0
!!  xcv(1:lenj) = reshape(cg(1:nxb,1:nyb), (/ nxb*nyb /))
!!  xcv(lenj+1:lenw) = rb
  iter =0
  vv = 0.d0; hq = 0.d0; hj = 0.d0; serr = 0.d0
  do j = 1, mit
    ser0 = 0.
    !!call matv(xcv, yy, iblst, fqlst,ilq, fplst, ilp, cox, coy, isel)
    !!deb xcv(:) = 1.
    !!call amatv(xcv, yy, nl, ux, uy, iaary, llen, rhsv(1:lenj), time, info, dt, mydif)
    call amatv(xcv, yy, nl, ux, uy, iaary, llen, rgg, time, info, dt, mydif)
!!deb    print *, 'check 1',maxval(yy(1:lenj)), maxval(yy(lenj+1:lenw)), lenw, lenj
    !!deb Y(1:nxb,1:nyb)=(reshape(yy, (/ nxb, nyb /)))
    !!deb print *, 'debugging'

    r =  rhsv - yy
    r = r*prem !PRECONDITIONER
    dr1 = dsqrt(dot_product(r,r))
    !!dr1 = sqrt(innerless(lenw,r,r,idz))
    if(dr1 .le. 1.d-14) then
      print *, 'guess solved EQN already'
      return
    endif
    
    vv(:,1) = r/dr1
    s = 0.
    s(1) = dr1
      do i = 1, mm-1
        iter = iter + 1
        vk = vv(:,i)

        !!call amatv(vk,w, nl, ux, uy, iaary, llen, rhsv(1:lenj), time, info, dt, mydif)
        call amatv(vk,w, nl, ux, uy, iaary, llen, rgg, time, info, dt, mydif)
!!deb        print *, 'check 2',maxval(w(1:lenj)), maxval(w(lenj+1:lenw)), lenw, lenj
        w = w*prem !PRECONDITIONER
        do k = 1, i
          vk = vv(:,k)
          hq(k,i) = dot_product(w,vk)
          !!hq(k,i) = innerless(lenw,w,vk,idz)
          w = w - hq(k,i)*vk
        enddo

        hq(i+1,i) = sqrt(dot_product(w,w))
        !!hq(i+1,i) = sqrt(innerless(lenw,w,w,idz))
        if (abs(hq(i+1,i)) < 1d-14) then
          infon = 1
        else
          infon = 0
          vv(:, i+1) = w/hq(i+1,i)
        endif

        do k = 1, i-1
          hit = hj(k,1)*hq(k,i) + hj(k,2)*hq(k+1,i)
          hq(k+1,i) = -hj(k,2)*hq(k,i) + hj(k,1)*hq(k+1,i)
          hq(k,i) = hit
        enddo
        if(infon .eq.0) then
          dthe = sqrt(hq(i,i)*hq(i,i)+hq(i+1,i)*hq(i+1,i))
          hj(i,1) = hq(i,i)/dthe
          hj(i,2) = hq(i+1,i)/dthe
          st = hj(i,1)*s(i) + hj(i,2)*s(i+1)
          s(i+1) = -hj(i,2)*s(i) + hj(i,1)*s(i+1)
          s(i) = st

          hit = hj(i,1)* hq(i,i) + hj(i,2)*hq(i+1,i)
          hq(i+1,i) = -hj(k,2)*hq(i,i) + hj(i,1)*hq(i+1,i)
          hq(i,i) = hit
        endif

        ser1 = abs(s(i+1))
        if (ser1 .gt. 1.0d-15) then
           serr = abs((ser1-ser0)/ser1)
        else
           serr = 0.0
        endif
        if(ser1 .le. tol .or. serr .le. tol) then
          serr = tol - 1.0d-15
          mi = i
          goto 100
        endif

        ser0 = ser1
      enddo
      mi = mm - 1
 100  yy(mi) = s(mi)/(hq(mi,mi)+1.0d-14)
      do k=mi-1,1,-1
         yy(k) = s(k)
         do j1 = k+1,mi
           yy(k) = yy(k) - hq(k,j1)*yy(j1)
         enddo
         yy(k) = yy(k)/hq(k,k)
      enddo        
      x11 = xcv
      do k=1,mi
         x11 = x11 + yy(k)*vv(:,k)
      enddo
      !!call amatv(x11,w, nl, ux, uy, iaary, llen, rhsv(1:lenj), time, info, dt, mydif)
      call amatv(x11,w, nl, ux, uy, iaary, llen, rgg, time, info, dt, mydif)
      r = rhsv - w
      r = r*prem !PRECONDITIONER
      s(mi+1) = sqrt(dot_product(r,r))
      !!s(mi+1) = sqrt(innerless(lenw,r,r,idz))
      xcv = x11
      if( abs(s(mi+1)) .lt. tol .or. serr .lt. tol ) then
        ! this should be the return point 
        print '(1x,"IN   net GMRES, iter", i5,1x, e14.6,5(i5,1x))', iter, serr, nl, llen
!!	print *, mi, s(mi), s(mi+1), maxval(ra), maxval(rg), maxval(prem), s(1)
        cg(1:nx,1:ny) = (reshape(xcv(1:nx*ny), (/ nx, ny /)))
        rb = xcv(lenj+1:lenw)
        deallocate(rhsv,xcv,prem)
        return
      endif
      error = s(mi+1)
  enddo
  print *, 'ENDGMRES w/ MAX iter', iter, serr
!=======================================================================
! end of GMRES
!=======================================================================
!
  deallocate(rhsv,xcv,prem)
!
  return
end subroutine nt_solver
!
! Advance adv-diff equation, which is solved by a fractional step
! method. The advection step assumes the variable p(si) is defined ever-
! ywhere (even inside the moving cell)
!=======================================================================
subroutine nt_advDiff(ntpn, ntpc, idn, isel, time, mkk, iaary, llen, info,dt,mydif)
  implicit none
  type(iapt), dimension(:) :: iaary(nent-2)
  integer :: idn(-1:nxb+2,-1:nyb+2)  ! tagger matrix for each grid point
  double precision, dimension(:,:,:) :: mkk(7,nring,nent-2)
!!  type(ibpt), dimension(:) :: iblst(nent-2)
  integer :: isel, llen(nent-2), info
  type(cc_bp), pointer :: curr
  double precision :: time, dt, mydif
  character(40) :: efile, ffile
  integer :: nt, strlen
  double precision :: ntpn(-1:nxb+2,-1:nyb+2)  ! chemical concentration new
  double precision :: ntpc(-1:nxb+2,-1:nyb+2)  ! chemical concentration old(current)
  double precision :: ntpp(-1:nxb+2,-1:nyb+2)  ! temparary
!
  integer :: i, j, k, il, ierr, iter, it, jt, ik, jk, im(2), ip(2)
  type(cc_augvar), pointer :: cp_lp
  double precision, dimension(-1:nxb+2,-1:nyb+2) :: rg, xx
  double precision, dimension(-1:nxb+2,-1:nyb+2,2) :: vmat
  double precision :: dx, dy, eps, resid, normb, tp1, tp2, tp3
  integer :: nl
  double precision :: nmvec(nring), npvec(nring)
!
  double precision :: ra(nxb*nyb/2), rb(nxb*nyb/2) !may have problem in definition
!
!!!  if ( isel .eq. -1) then
!!!    goto 9999
!!!  endif

  nl = 0
  do il = 1, nib ! # of aug variables @ interface
!!twosides    nl = nl + 2*llen(il)
    nl = nl + llen(il)
  enddo

  ierr = 0
  rb = 0.; ra = 0.
 
!  print '(1x, "TIME IS: ", e13.6)', time
  eps = 1.d-12; resid = eps
  dx = hg; dy = hg

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! start the advancing 
!========================================================================
  do i = 1, nxb! set ghost cell @top/bottom w/ no flux, not used in linear solve
!    x = xmin+(dble(i)-xshift)*dx; 
!    y = ymin bottom, no flux BC
!!rmtemp    ccpo(i,0)  = ccpo(i,1)!! + h*cc_btt(i,1)
!!    ntpc(i,0)  = ntpc(i,1)!! + h*cc_bto(i,1)
!!    ntpn(i,0)  = ntpn(i,1)!! + h*cc_btn(i,1)
!    y = ymax top, no flux BC
!!rmtemp    ccpo(i,nyb+1)= ccpo(i,nyb)!! + h*cc_btt(i,2)
!!    ntpc(i,nyb+1)= ntpc(i,nyb)!! + h*cc_bto(i,2)
!!    ntpn(i,nyb+1)= ntpn(i,nyb)!! + h*cc_btn(i,2)
  enddo

  do j = 1, nyb! set periodic @left/right ghost, not used in linear solve
!    x = xmin left, periodic BC
    i = 0
!!rmtemp    ccpo(i,j) = ccpo(nxb,j)
    ntpc(i,j) = ntpc(nxb,j)
    ntpn(i,j) = ntpn(nxb,j)
!    x = xmax right, periodic BC
    i = nxb+1
!!rmtemp    ccpo(i,j) = ccpo(1,j)
    ntpc(i,j) = ntpc(1,j)
    ntpn(i,j) = ntpn(1,j)
  enddo
!========================================================================
   
!!  vmat(-1:nx+1,-1:ny+1,1) = eta/(eta+etas)*unc(-1:nx+1,-1:ny+1)
!!  vmat(-1:nx+1,-1:ny+1,2) = eta/(eta+etas)*vnc(-1:nx+1,-1:ny+1)
!!  !!rg(1:nxb, 1:nyb) = ntpc(1:nxb, 1:nyb) + dlt*ccps(1:nxb, 1:nyb) 
!========================================================================
!
  ivf = 0
  do il = 1, nib
    cp_lp =>iaary(il)%p
    do i = 1, llen(il)
!!      call nt_setRHS(ntpc, cp_lp, llen, il, ivf, idn, mkk, -1, iaary, info, dt)
!
      call nt_setRHS(ntpc, cp_lp, llen, il, ivf, idn, mkk,  1, iaary, info, dt)
!
      cp_lp => cp_lp%next
    enddo
  enddo
!========================================================================
! set the extra-cellular part to always constant
  igd = reshape(iid(1:nx,1:ny), (/ nx*ny /))
  !!ntpp(1:nx,1:ny) = ntpc(1:nx,1:ny)*dble(iid(1:nx,1:ny)) + (one-dble(iid(1:nx,1:ny)))*nsigma
  ntpp(1:nx,1:ny) = ntpc(1:nx,1:ny)*dble(iid(1:nx,1:ny)) + (one-dble(iid(1:nx,1:ny)))
  ntpc(1:nx,1:ny) = ntpp(1:nx,1:ny) ! debugging mode
!========================================================================
  do j = 1, ny
  do i = 1, nx
    pavg(i,j) = nt_rDthetaDsigma(ntpc(i,j))
  enddo
  enddo
  pnu = reshape(pavg(1:nx,1:ny), (/ nx*ny /))
  rg(1:nxb, 1:nyb) = ntpc(1:nxb, 1:nyb)!! + dlt*ccps(1:nxb, 1:nyb) 
!!deb  print *, 'max', maxval(pnu), minval(pnu)
!=======================================================================
  
!!  nt = nint(time/dlt)
!!  strlen = len_trim(runname)
  normb = sqrt(mydot(rg,rg))
  if (abs(normb) < eps) then
    xx = 0.d0
    iter = 1
    resid = 0.d0
    ntpn = xx
    print *, 'GET ZERO here!!'
    return
  else
    j = 0
    do il = 1, nib
      ! set RHS for aux variabls defined @ grid crossing
      cp_lp => iaary(il)%p
      do i = 1, llen(il)
        !!tp1 = cp_lp%s
        !!tp2 = paraval(tp1, nring, mkq, il, 1, 0, info) ! on - side old time chemical
        !!tp3 = nt_dthsig(tp2)
        !!tmp = cc_gate(cp_lp%s,time)/cp_lp%dxs! constpmp
        !!tmp = 0.d0 ! pump is concentration dependent

!!twoside        ! on "-" side
!!twoside        j = j + 1
!!twoside        !!tmp = 0.d0 ! no flux interface cond!
!!twoside        ra(j) = cp_lp%ncm
!!        cp_lp%rs_m = zero
!
        ! on "+" side
        j = j + 1
        !!tmp = 0.
!        ra(j) = rho0/(eta+etas)*(cp_lp%nx)
        ra(j) = cp_lp%ncp
!!        print '(1x, " test data ", 10(e14.6,1x))', tp1, tp2, tp3, cp_lp%udn, ra(j)
        !!ra(j) = kp(il)
        !!cp_lp%rs_p = kp(il)
!
        cp_lp => cp_lp%next
      enddo
    enddo
!=======================================================================
    !!xx = ntpc
    !!rb = npsi
    rb = one
    xx = rg
!    xx = zero
!!deb    cp_lp => iaary(1)%p
!!deb    j = 0
!!deb    do i = 1, nl/2
!!deb      j = j+1
!!deb      print '("rhs ",i5, 1024(e14.6,1x))', i, ra(2*i-1), ra(2*i), ndv(i,1,1), ndv(i,2,1), &
!!deb        cp_lp%ncp, cp_lp%ncm
!!deb    enddo
!!deb    print *, 'Are we here? before solver'
    call nt_solver(rg,ra(1:nl), xx, rb(1:nl), nl, iaary, llen,info, isel,time,dt,mydif)
    ntpn = xx
!!    print '(1024(e14.6,1x))', rb(1:nl)
    do il = 1, nib
      !save aux chemicals @ grid crossing
      j = 0
      cp_lp => iaary(il)%p
      do i = 1, llen(il)
!twoside        j = j+1
!twoside!twoside        cp_lp%nm = rb(j) ! for getting mkq
!twoside        nmvec(i) = rb(j)
!
        j = j+1
        cp_lp%np = rb(j) ! for getting mkq
        mkq(1,i,il) = cp_lp%s
        npvec(i) = rb(j)
        nmvec(i) = rb(j)
!
        cp_lp => cp_lp%next
      enddo
    enddo
!========================================================================
  endif
!
!!!9999 continue
!=======================================================================
  do il = 1, nib
    cp_lp => iaary(il)%p
    !nmvec & npvec are set abovek
    mkq(1,:,il) = mkn(1,:,il)
    !!print *, 'ck here', nl
    call smoothcrossing(olen, nib,il,mkq,iaary,llen, nmvec,npvec, info)
  enddo
!=======================================================================
! on output, mkq store interpolation of network volume @ + & - side of grid crossing
!=======================================================================
!
  return
end subroutine nt_advDiff
!=======================================================================
!=======================================================================
! We will use this subroutine to initiate all the computations. So this
! one will call nt_advDiff to advance one time step in chemical
!=======================================================================
subroutine nt_IBbdy(ntpn, ntpc, ibary, iaary, llen, mk, isel, time, dt, mydif)
  implicit none
  type(ibpt), dimension(:), intent(in) :: ibary(nent-2)
  type(iapt), dimension(:), intent(in) :: iaary(nent-2)
  double precision, dimension(:,:,:) :: mk(7,nring,nent-2), nb(nring,2, nent-2)
  integer :: llen(nent-2)
  integer :: isel, im(2), iloc(4,2), ipl, ict
  double precision :: time, dt, mydif, st_p(4)
  double precision :: ntpn(-1:nxb+2,-1:nyb+2)  ! chemical concentration function
  double precision :: ntpc(-1:nxb+2,-1:nyb+2)  ! chemical concentration function
!!  double precision :: ccer(-1:nxb+2,-1:nyb+2)  ! chemical concentration function
!!  double precision :: ccrs(-1:nxb+2,-1:nyb+2)  ! chemical concentration function
!
  integer :: i, j, k, ik, jk, il, nt, info, strlen
  type(cc_bp), pointer :: curr
  type(cc_augvar), pointer :: cp_lp, cp_lq
  double precision ::  dx, dy, eps, rh, tmp, tp1, tp2, tp3, tp4, tpp, s0, &
     nx, ny, vx, vy
!!  integer :: ip, jp, im, jm, is, js, nw, ns, it, jt, i1, i2
  double precision :: cbq(nring)
  character(40) :: efile, ffile
!!  logical :: fl(5), as, alll, allr
!!  logical :: al(nent-2)
!!  character :: equed
  double precision :: cmp(nring,6,nent-2)
!
  eps = 1.d-8
  dx = hg; dy = hg
  tupi = 2.d0*cpi
!
!=======================================================================
! Read positions of IB points from (xpt,ypt), which will define the cell
! region and the boundary is interpolated from those pts. Then the index
! (i,j) is calculated according to the offset to (xmin,ymin) and h, and
! the index is always the left bottom corner of the square formed by the
! cell centers (computational grids starting from (1,1) to (m,n))
!=======================================================================
  xamin = xmin; yamin = ymin
  rh= one/hg
!================

!///////////////////////////////////////////////////////////////////////

  if (isel .eq. -1) then
    goto 5555
  endif
!
  call vc_cfs(xpt,ypt, vc_fe) ! this can't be used in two phase model
  do il = 1, nib ! 
    ! stores the speed of -(dx/dt-eta/(eta+etas)u)
!!    tmp = eta/(eta+etas)
!!exp     curr => ibary(il)%p
!!exp     do i = 1, nring
!!exp       cbq(i)=  (curr%vx)*curr%nx+(curr%vy)*curr%ny
!!exp !!      print '(1x,"u-dxdt", 6(e14.8,1x))', curr%vx, curr%vy, curr%ux, curr%uy, curr%nx, curr%ny
!!exp       curr => curr%prev
!!exp     enddo
    cp_lp => iaary(il)%p
    do i = 1, llen(il)
      s0 = cp_lp%s; 
      !!exp call brs(nring,nib,il,mk,cbq  ,s0,tmp, cp_lp%k,time,0)
      !!exp nx = cp_lp%nx; ny = cp_lp%ny;
      !!exp call cc_VelatIB(cp_lp%xb, cp_lp%yb, vx, vy, unc, vnc, info)
      !!exp tmp = tmp-eta/(etas+eta)*(nx*vx+ny*vy)
      !!exp tp1 = paraval(s0, nring, mkq, il, 1, 0, info) ! on + side old time chemical
      !!exp tp2 = nt_rDthetaDsigma(tp1)
      !!exp tp3 = nt_theta(tp1)
      !!exp cp_lp%vdn = tp2 ! we change the %vdn value, which is used in chemical_mod
      !!exp cp_lp%udn = tp3*kw(il) ! (theta(\sigma_n))
      !!exp tp2 = paraval(s0, nring, mkr, il, 0, 0, info) ! on + side sigma_a
      !!exp cp_lp%ncp =-tp3*tmp-tp2*nx/(eta+etas)       +Jactin(s0,time) !
      !!exp cp_lp%ncm = one
      call cc_VelatIB(cp_lp%xb, cp_lp%yb, vx, vy, unc, vnc, info)
      nx = cp_lp%nx; ny = cp_lp%ny;
      tmp =etas/(etas+eta)*(nx*vx+ny*vy)
      tp1 = paraval(s0, nring, mkq, il, 1, 0, info) ! on + side old sigma_n
      tp2 = nt_theta(tp1) ! theta(sigma_n)
      cp_lp%udn = tp2*kw(il) ! (theta(\sigma_n))

      tp3 = paraval(s0, nring, mkr, il, 0, 0, info) ! on + side gradient sigma_a
      tp4 = paraval(s0, nring, mkr, il, 1, 0, info) ! on + side sigma_a
      call brs(nring,nib,il,mkn,vc_fe,s0,tpp, cp_lp%k,time,0)
      !!cp_lp%ncp =-tp3*nx/(eta+etas)-tp2*tmp+Jactin(s0,time)/cp_lp%dxs -tp2*kw(il)* ( &
      !!oldcp_lp%ncp =-rho0*tp3*nx/(eta+etas)-tp2*tmp+Jactin(s0,time) -tp2*kw(il)* ( &
      !!old  tp4*rtr +tpp/rtc) ! missing jump in chemicals
      tp1 = paraval(s0, nring, mkp, il, 1, 0, info) ! on + side old chemical
      tp1 = tp1 - paraval(s0, nring, mkp, il, 0, 0, info) ! on + side old chemical
      cp_lp%ncp =-rho0*tp3*nx/(eta+etas)-tp2*(tmp-tp1*kw(il))+Jactin(s0,time) -tp2*kw(il)* ( &
        tp4*rtr +tpp/rtc) ! 
      !!cp_lp%ncp =-rho0*tp3*nx/(eta+etas)-tp2*tmp+Jactin(cp_lp%sl,time)/cp_lp%dxs -tp2*kw(il)* ( &
      !!  tp4*rtr +tpp/rtc) ! missing jump in chemicals
      cp_lp%ncm = one*1d-8
      
!!deb      print '("mark ", (i5,1x), 1024(e14.6,1x))', i,s0,tmp,cp_lp%xb,cp_lp%yb, &
!!deb         vx,vy,nx,ny, tp3, tpp, tp1, tp2, cp_lp%ncp, cp_lp%udn, cp_lp%vdn, cp_lp%dxs
      !
      cp_lp => cp_lp%next
    enddo
  enddo

!///////////////////////////////////////////////////////////////////////
!=======================================================================
! Advance advection diffusion system one step
  5555 continue
  info = 0
  call nt_advDiff(ntpn,ntpc, idn, isel, time, mko, iaary,llen, info,dt,mydif)

!!  tp1 = mydot(ntpn,ntpn); tp2 = mydot(ntpc,ntpc)
!!  print *, 'norm ', tp1, tp2, maxval(abs(ntpn))
!=======================================================================
  cmp = 0
    do il = 1, nib ! save lagrangian chemicals on the new time level
      curr => ibary(il)%p
      !!cp_lp=>iaary(il)%p !save @ grid crossing
      !!do i = 1, llen(il) ! save grid crossings
      do i = 1, nring ! save @ IB points
        !!s0 = cp_lp%s; ! @ grid crossing
        s0 = curr%s; ! @IB points
        tp1 = paraval(s0, nring, mkq,il,1,0,info);
        tp2 = nt_theta(tp1)
        cmp(i,1,il) = tp1
        cmp(i,2,il) = tp2

        !!call cc_VelatIB(cp_lp%xb, cp_lp%yb, vx, vy, unc, vnc, info)
        !!nx = cp_lp%nx; ny = cp_lp%ny;
        !!tp1 = (vx*nx+vy*ny)/(eta+etas)
        !!cmp(i,1,il) = tp1
        !!cmp(i,2,il) = cp_lp%udn/kw(il) 
        cmp(i,3,il) = curr%x
        cmp(i,4,il) = curr%y
        !!cmp(i,4,il) = cp_lp%s
        !!cmp(i,5,il) = curr%uy
        tp1 = paraval(s0, nring, mkq,il,1,0,info);
        call nrmgrad(tp2, cp_lp, nt_pn, tp1)
        cmp(i,5,il) = s0
        cmp(i,6,il) = tp2
!!        if (idf(im(1),im(2)) >0) then
!!          print *, 'wrong pt stop'
!!          stop
!!        endif
!   
        curr => curr%prev
        !!cp_lp => cp_lp%next !save @ grid crossing
      enddo
    enddo
!!
!!if (isel .ne. -1) then
    nt = nint(time/dt)
    strlen = len_trim(runname)
  if (nt-nfreq*(nt/nfreq) .eq. 0) then
    write(efile,'(2a,i4.4)') runname(1:strlen),'.nib.',outcount+1
    open(77,file=efile,access='stream',action='write')
    write(77)cmp
    close(77)
    write(ffile,'(2a,i4.4)') runname(1:strlen),'.ns.',outcount+1
    open(67,file=ffile,access='stream',action='write')
    write(67)ntpn
    close(67)
    write(efile,'(2a,i4.4)') runname(1:strlen),'.chk.',outcount+1
    open(88,file=efile,form='formatted',action='write')
    do i = -1, nxb+2
      write(88,'(1x,1024(i5,1x))')(idf(i,j),j=-1,nyb+2)
    enddo
    close(88)
    write(efile,'(2a,i4.4)') runname(1:strlen),'.cib.',outcount+1
    open(77,file=efile,access='stream',action='write')
    write(77)cmp
    close(77)
  endif
!=======================================================================
!!deb  tp1 = sqrt(mydot(ntpc,ntpc))*hg; tp2 = sqrt(mydot(ntpn,ntpn))*hg
!!  tp1 = sqrt(mydot(ntpc,ntpc))*hg; tp2 = sqrt(mydot(ntpn,ntpn))*hg
!!  tp3 = sum(xpt)/dble(nring); tp4 = maxval(xpt); tpp = minval(xpt)
!!  eps = sum(ypt)/dble(nring)
  
!!  print '(1x,"In network: ", 2(i5,1x), 10(e16.8,1x))', im, idf(im(1),im(2)), &
!!  print *, im
!!deb  im = maxloc(ntpn)
!!  print *, im
  !!print '(1x,"In network: ", 2(i5,1x), 10(e16.8,1x))', im, &
!!deb  print '(1x,"In network: ", 10(e16.8,1x))',  &
!!deb     tp1, tp2, maxval(ntpn), minval(ntpn), maxval(ntpn), minval(ntpn)!, &
!        ntpn(im(1), im(2)), ntpc(im(1),im(2))
!=======================================================================

!=======================================================================
! deallocate memorys here
 4001 format(1x,1024(e18.6,1x))
!TT
  return
end subroutine nt_IBbdy
!
!
subroutine nt_setRHS(rg, cp_lp, llen, il, ivf, idn, mkk, isel, iaary, info, dt)
  ! for freshly cleared pt
  implicit none
  double precision :: dt
  type(iapt), dimension(:) :: iaary(nent-2)
  integer :: ivf(-1:nxb+2,-1:nyb+2)  ! 
  integer :: idn(-1:nxb+2,-1:nyb+2)  ! 
  double precision :: rg(-1:nxb+2,-1:nyb+2)  ! chemical concentration function
  type(cc_augvar), pointer :: cp_lp
  double precision, dimension(:,:,:) :: mkk(7,nring,nent-2)
  integer :: isel, info, llen(nent-2), il
!
  integer :: i1(2), k, it, jt, i_is(2), i_js(2)
  double precision :: tp(2), x, y, xm, ym, s0, s1, vx, vy, dx, dy
  double precision :: st_is(2), st_js(2), i_st, j_st
  double precision :: tp1, tp2, tp3, tp4, tmp
!
  dx = hg; dy =hg;

  if (isel .eq. -1) then
    i1 = cp_lp%im; !cm = cp_lp%cm
  else if (isel .eq. 1) then
    i1 = cp_lp%ip; !cm = cp_lp%cp, intracellular space
  endif
  if ((idn(i1(1),i1(2))) .eq. 2.and. ivf(i1(1),i1(2)).eq.0) then ! freshly cleared pt
  !!if ((idn(i1(1),i1(2))) .eq. 2*isel.and. ivf(i1(1),i1(2)).eq.0) then
    !!rg(i1(1),i1(2)) = cm
    x = xmin+(dble(i1(1))-xshift)*dx; y = ymin+(dble(i1(2))-yshift)*dy; 
   print '(1x,"freshly cleared(net): side", 4(i3,1x),1024(e16.8,1x))', isel, idn(i1(1),i1(2)), &
      idf(i1(1),i1(2)), ivf(i1(1),i1(2)),x, y, nt_pc(i1(1),i1(2)), rg(i1(1),i1(2))
    call FindBdyPtForFP(s0,s1,il,x,y,xm,ym,mkk,k,1) ! (xm,ym) along old boundary
    if (isel .eq. 1) then
!!deb      print '("New value ", 2(e12.4,1x),1024(i4,1x))', rg(i1(1),i1(2)), nt_pc(i1(1),i1(2)), &
!!deb      i1, idf(i1(1),i1(2)), idn(i1(1),i1(2)), ivf(i1(1),i1(2))
      tmp = paraval(s1, nring, mkq, il, 1, 0, info)
      rg(i1(1),i1(2)) = tmp
!!deb      print '(" New value ", 1024(e12.5,1x))',x,y,xm,ym
      print '(" New value ", 1024(e12.5,1x))', tmp, dist(x,y,xm,ym)
!!    else
!!      tmp = paraval(s1, nring, mkq, il, 1, 0, info)
    endif
!
    !!vx = (xm-x)/dt; vy = (ym-y)/dt
    vx = (xm-x); vy = (ym-y)
    it = dsign(1.d0,dble(vx)); jt = dsign(1.d0,dble(vy))
    i_is(1) = i1(1) + it; i_is(2) = i1(2)
    i_js(1) = i1(1)     ; i_js(2) = i1(2) + jt
    !if (idf(cp_lp%m_is(1),cp_lp%m_is(2)).ne. idf(i1(1),i1(2)) .or. idf(cp_lp%m_js(1),cp_lp%m_js(2)).ne. idf(i1(1),i1(2))) then
    if (idf(i_is(1),i_is(2)).eq. idf(i1(1),i1(2))) then
      i_st = 1.d0
    else
      i_st = 0.d0
    endif
    if (idf(i_js(1),i_js(2)).eq. idf(i1(1),i1(2))) then
      j_st = 1.d0
    else
      j_st = 0.d0
    endif
    !if (idf(i_is(1),i_is(2)).ne. idf(i1(1),i1(2)) .or. idf(i_js(1),i_js(2)).ne. idf(i1(1),i1(2))) then
    if (i_st+j_st <0.5d0 ) then ! no enough resolution
      print *,  i_st, j_st, it, jt
      print *, 'there has to at least one point!'
      print '(1x,1024(i3,1x))', i_is, idf(i_is(1),i_is(2)), i_js, idf(i_js(1),i_js(2)),&
        i1,idf(i1(1),i1(2)),  isel, idn(i1(1),i1(2))
      print '(1x,1024(i3,1x))', idf(i_is(1),i_is(2)),idf(i1(1),i1(2)), idf(i_js(1),i_js(2))
      print '(1x,1024(i3,1x))', chker0(i_is(1),i_is(2)),chker0(i1(1),i1(2)), chker0(i_js(1),i_js(2))
      print '(1x,1024(i3,1x))', chker1(i_is(1),i_is(2)),chker1(i1(1),i1(2)), chker1(i_js(1),i_js(2))
      print '(1x,1024(i3,1x))', k, it, jt, i1, chker0(i1(1),i1(2)), chker1(i1(1),i1(2))!!, chker2(i1(1),i1(2))
      print '(1x,1024(e16.8,1x))', s0, s1, dist(x,y,xm,ym)/hg, dist(x,y,xm,ym)
!
      print *, ''
      print '(1x,1024(e16.8,1x))', x, y
      print '(1x,1024(e16.8,1x))', xm, ym

      print *, ''
      do it = 1, llen(1)
        print '(1x,1024(e16.8,1x))', cp_lp%xb, cp_lp%yb
        cp_lp => cp_lp%next
      enddo
      print *, ' finish intercept'
      !!cm = cpi*2./dble(nring)
      do it = 1, nring
        !!tp(1) =  paraval(cm*dble(it), nring, mkk, 1, 0, 0, info)
        !!tp(2) =  paraval(cm*dble(it), nring, mkk, 1, 1, 0, info)
        tp(1) =  paraval(mko(1,it,1), nring, mko, 1, 0, 0, info)
        tp(2) =  paraval(mko(1,it,1), nring, mko, 1, 1, 0, info)
        print '(1x,1024(e16.8,1x))', tp
      enddo
!!      print *, 'older'
!!      do it = 1, nring
!!        !!tp(1) =  paraval(cm*dble(it), mkk, 1, 0, 0, info)
!!        !!tp(2) =  paraval(cm*dble(it), mkk, 1, 1, 0, info)
!!        print '(1x,1024(e16.8,1x))', tp
!!      enddo

      stop
    endif
    !!tp(1) = vx*dt/hg; tp(2) =-vx*dt/hg !should there be 0.5 factor in front?
    tp(1) = vx/hg; tp(2) =-vx/hg !should there be 0.5 factor in front?
    if (it .eq. 1) then ! m_is(1) > i1(1)
      st_is(1) = tp(1) ! coef @ i1
      st_is(2) = tp(2) ! coef @ neighbor pt
    else if (it .eq. -1) then ! i1(1) > m_is(1)
      st_is(1) = tp(2) ! coef @ i1
      st_is(2) = tp(1) ! coef @ ineighbor pt
    else
      print *, 'no such case, we need to stop to check!'
      print *, 'it ', it, vx, vy, i1
      stop
    endif
    st_is = st_is*i_st
    !!tp(1) = vy*dt/hg; tp(2) =-vy*dt/hg
    tp(1) = vy/hg; tp(2) =-vy/hg
    if (jt .eq. 1) then
      st_js(1) = tp(1) ! coef @ i1
      st_js(2) = tp(2) ! coef @ neighbor pt
    else if (jt .eq. -1) then
      st_js(1) = tp(2) ! coef @ i1
      st_js(2) = tp(1) ! coef @ neighbor pt
    else
      print *, 'no such case, we need to stop to check 2!'
      stop
    endif
    st_js = st_js*j_st
    if (isel .eq. -1) then
      cp_lp%m_is(:) = i_is(:);
      cp_lp%m_js(:) = i_js(:)
      cp_lp%stm_is  = st_is
      cp_lp%stm_js  = st_js
      cp_lp%sm = s1; cp_lp%xm = xm;; cp_lp%ym = ym
      cp_lp%mi_st = i_st
      cp_lp%mj_st = j_st
    else if (isel .eq. 1) then
      cp_lp%p_is(:) = i_is(:);
      cp_lp%p_js(:) = i_js(:)
      cp_lp%stp_is  = st_is
      cp_lp%stp_js  = st_js
      cp_lp%sp = s1; cp_lp%xp = xm;; cp_lp%yp = ym
      cp_lp%pi_st = i_st
      cp_lp%pj_st = j_st
    else
      print *, 'can not be here!'
      stop
    endif
    ivf(i1(1),i1(2)) = isel
  endif
!
  return
end subroutine nt_setRHS
!!
!!contains
double precision function nt_sga(theta)
  double precision theta
!
  nt_sga = -0.00
!
  return
end function nt_sga
!
!
double precision function nt_sigma(theta)
  !return \sigma_n(\theta)
  double precision theta
!
  !!!!nt_sigma = -rho0*(                           -theta*half/theta0)! free
  !!nt_sigma = -rho0*((theta0/theta)**(one/three)-theta*half/theta0)! free
  !!nt_sigma = -rho0*((theta0/theta)-theta*half/theta0)! uniaxial
  nt_sigma = -((theta0/theta)-theta*half/theta0)! uniaxial
!
  return
end function nt_sigma
!
double precision function nt_dsigma(theta)
  !return d\sigma_n/d\theta
  double precision theta
!
  !!!!nt_dsigma = rho0*(                                        half/theta0)! free
  !!nt_dsigma = rho0*((theta0/theta)**(one/three)/three/theta+half/theta0)! free
  nt_dsigma = rho0*((theta0/(theta*theta))+half/theta0) ! uniaxial
!
  return
end function nt_dsigma
!
double precision function nt_theta(sig)
  double precision :: sig
! The diffusion coefficient for sigma_n equation
  !!nt_theta = theta0*(sig+sqrt(two*rho0*rho0+sig*sig))/rho0
  nt_theta = theta0*(sig+sqrt(two+sig*sig)) ! scaled
!
  return
end function nt_theta
!
double precision function nt_rDthetaDsigma(sig)
  double precision :: sig
! The diffusion coefficient for sigma_n equation
  !!nt_rDthetaDsigma = rho0/theta0/(one+sig/sqrt(two*rho0*rho0+sig*sig))
  nt_rDthetaDsigma = rho0/theta0/(one+sig/sqrt(two+sig*sig)) ! scaled
!
  return
end function nt_rDthetaDsigma

double precision function Jactin(s, time)
implicit none
  double precision :: s, time
  integer :: iside, isel
  ! we have to make sure s \in (0,2pi)
  double precision :: shl, shw, shs ! location, width, & strength leading gaussian
  double precision :: stl, stw, sts ! location, width, & strength back gaussian
  double precision :: sul, suw, sus ! location, width, & strength back gaussian
  double precision :: sbl, sbw, sbs ! location, width, & strength back gaussian

  ! sth*exp(-(s-s?l)^2/s?w^2)
  double precision :: tmp, su, sl, s0, tupi
  double precision :: pah, paw, pwt
!
  tupi = two*cpi
  tmp = zero
!
  s0 = mod(mod(s, tupi)+tupi, tupi) ! archlength we will used

! if both shs and sts are -, we have shrinking
!  shs=-.10;      sts= 1.0*shs
! if both shs and sts are +, we have swelling
!  shs= .10;      sts= 1.0*shs
! if  shs and sts have opposite -/+, we move the cell
  shs= .005d0;      sts=-1.0d0*shs
  sus=-shs*0.0d0;  sbs= 1.0d0*sus
  shl=cpi*.0d0 ;   stl=cpi*1.0d0
  sul=cpi*.80d0;   sbl=cpi*1.20d0
  shw=cpi*.210d0;  stw=shw*1.0d0;
  suw=shw*.5d0 ;   sbw=suw*1.0d0

!!o  tmp=shs*exp(-.5*(s-shl)**2./shw**2.)+0. *exp(-.5*(s-shl-tupi)**2./shw**2.)+&
!!o      sts*exp(-.5*(s-stl)**2./stw**2.)
!working  tmp=shs*exp(-.5*(s-shl)**2./shw**2.)+shs*exp(-.5*(s-shl-tupi)**2./shw**2.)+&
!working      sts*exp(-.5*(s-stl)**2./stw**2.)
  tmp=shs*exp(-.5d0*(s-shl)*(s-shl)/(shw*shw))+ &
      shs*exp(-.5d0*(s-shl-tupi)*(s-shl-tupi)/(shw*shw))+&
      sts*exp(-.5d0*(s-stl)*(s-stl)/(stw*stw))+ &
      sus*exp(-.5d0*(s-sul)*(s-sul)/(suw*suw))+ &
      sbs*exp(-.5d0*(s-sbl)*(s-sbl)/(sbw*sbw))
!!2  tmp=shs*exp(-.5*(s-shl)**2./shw**2.)+sts*exp(-.5*(s-stl-tupi)**2./shw**2.)
!!  tmp=shs*exp(-(s-shl)**2./shw**2.)/sqrt(2.*cpi)/shw+sts*exp(-(s-stl)**2./stw**2.)/sqrt(2.*cpi)/stw
!!  tmp=shs*exp(-(s-shl)**2./shw**2.)/sqrt(2.*cpi)/shw+sts*exp(-(s-stl)**2./stw**2.)/sqrt(2.*cpi)/stw+shs*exp(-(s-shl-tupi)**2./shw**2.)/sqrt(2.*cpi)/shw
!!!
!!steps  tmp= (pah-cosh(sin(4.*(s-cpi)**4.)))*myhv(s-paw)*myhv(tupi-paw-s)*sts+&
!!steps       (pah-cosh(sin(4.*(s-tupi)**4.)))*myhv(s-cpi-paw)*shs         +   &
!!steps       (pah-cosh(sin(4.*(s)**4.)))*myhv(cpi-paw-s)*shs
!!  if (time > 12.and. time < 30 .or. time > 45) then
!!    tmp = -tmp
!!  endif
  !!Jactin = dble(isel)*tmp! don't use this anymore
  !!tmp = zero
!!  if ( zero <= s .and. s <= shw .or. tupi - shw <= s .and. s <= tupi) then
!!    tmp = shs
!!  endif
!!  if ( shw < s .and. s < cpi - stw .or. cpi+ stw < s .and. s < tupi - shw) then
!!    tmp = zero
!!  endif
!!  if ( cpi - stw <= s .and. s <= cpi + stw ) then
!!    tmp = sts
!!  endif

  Jactin =tmp
  !!cc_gate = 0.1
!
  return
end function Jactin

subroutine nrmgrad(nrm, cp_lp, cc, cbval)
  double precision :: nrm, cbval
  double precision :: cc(-1:nx+1,-1:ny+1)
  type(cc_augvar), pointer :: cp_lp
  !
  integer :: j, gpl, opl, im(2), ip(2), g_ali(4,2), g_off(4,2)
  double precision :: c_ali(4), c_off(4), tmp
  

      g_ali = 0; g_off = 0; c_ali = 0.; c_off = 0.;
      g_ali = cp_lp%i_ali_grd_p; c_ali = cp_lp%c_ali_grd_p; 
      g_off = cp_lp%i_off_grd_p; c_off = cp_lp%c_off_grd_p; 
      !!cpl   = cp_lp%ipl;
      gpl = g_ali(4,1); opl = g_off(4,1)

      tmp = 0.d0
      !!do j = 1, cpl
      do j = 1, gpl
        !!call sub2ind(nxb,g_ali(j,1),g_ali(j,2),k)
        !!tmp = tmp + Y(k)*c_ali(j+1)
        tmp = tmp + cc(g_ali(j,1),g_ali(j,2))*c_ali(j+1)
      enddo
      !!do j = 1, g_off(4,2)
      do j = 1, opl
        !!call sub2ind(nxb,g_off(j,1),g_off(j,2),k)
        !!tmp = tmp + Y(k)*c_off(j+1)
        tmp = tmp + cc(g_off(j,1),g_off(j,2))*c_off(j+1)
      enddo
      tmp = tmp + (c_ali(1)+c_off(1))*cbval

  nrm = tmp
  !
  return
end subroutine nrmgrad
!
end module network_mod
