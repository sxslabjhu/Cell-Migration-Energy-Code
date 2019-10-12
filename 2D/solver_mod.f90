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

module solver_mod
! solving Laplace equation with Dirichlet BC over the cell membrane
  use, intrinsic :: iso_c_binding
  use parameters
  use geometry
!
  implicit none
!
  private
  integer, parameter :: nxb = nx, nyb = ny
  integer :: igd(1:nx*ny)
  !!double precision, parameter :: sik =-2.d-4
  double precision, parameter :: sik = 0.d0

  public :: LapSolver, sg_IBbdy, siga, sik

!!  abstract interface 
!!    subroutine matv(x, Ax, nl, ux, uy, iaary, llen, dt, mydif, info)
!!      integer :: llen(1), info, nl
!!      double precision, dimension(:) :: ux, uy, rhsv
!!      double precision :: time, dt, mydif
!!      double precision, dimension(:) :: x, Ax
!!      type(iapt), dimension(:) :: iaary(1)
!!    end subroutine matv
!!  end interface


contains

!========================================================================
subroutine sg_IBbdy(sga, iaary,llen,time,dt,info)
  type(iapt), dimension(:) :: iaary(nent-2)
  integer :: nl, info, llen(nent-2), isel
  double precision, dimension(-1:nx+2,-1:ny+2) :: sga
!!  double precision, dimension(:,:,:) :: mk(7,nring,nent-2)
  double precision, intent(in) :: time, dt
!
  integer :: i, j, il, im(2), ip(2), imv(-1:nx+2,-1:ny+2), gpl, opl, nt, strlen
  character*40 :: efile
  type(cc_augvar), pointer :: cp_lp
  double precision ::  dx, dy, tmp, tp1, tp2, s0, xbc, st(4)
  integer :: g_ali(4,2), g_off(4,2)
  double precision :: c_ali(4), c_off(4)
  double precision, dimension(-1:nxb+2,-1:nyb+2) :: rhs, xx
  double precision :: smvec(1:nring), spvec(nring)
!
  tmp = one/hg/hg
!!  tmp = one
  rhs(1:nx,1:ny) = one-dble(iid(1:nx,1:ny))

!!deb  smvec = zero; spvec = zero
  do il = 1, nib
    mkr(1,:,il) = mkn(1,:,il)
    cp_lp => iaary(il)%p
    do i = 1, llen(il)
      s0 = cp_lp%s; ip = cp_lp%ip; st = cp_lp%st_p
      xbc = siga(s0,cp_lp%xb,cp_lp%yb,1)
      cp_lp%sgm = xbc
      spvec(i) = xbc
!!deb      spvec(i) = xbc
      rhs(ip(1),ip(2)) = rhs(ip(1),ip(2)) + st(1)*xbc*tmp
!
      cp_lp => cp_lp%next
    enddo
    smvec = spvec
!!deb    call smoothcrossing(olen, nib,il,mkr,iaary,llen, smvec,spvec, info)
!!deb    tp1 = paraval(half, nring, mkr, il, 0, 0, info) ! on + side old network volume
!!deb    tp2 = paraval(half, nring, mkr, il, 1, 0, info) ! on + side old network volume
!!deb    print '("really", 1024(e14.6,1x))', tp1, tp2, s0, siga(s0,half,half,0), mkr(:,1,il)
!!    print *, spvec(1:llen(il))
  enddo

  igd = reshape(iid(1:nx,1:ny), (/ nx*ny /))
!
!!simplified  call LapSolver(rhs,xx, nl, iaary, llen, info)
!
!========================================================================
! debug
!!  tp1 = siga(0.5)
!!  xx = tp1
!========================================================================
  sga = xx
!
  do il = 1, nib
    cp_lp => iaary(il)%p
    do i = 1, llen(il)
      g_ali = 0; g_off = 0; c_ali = 0.; c_off = 0.;
      g_ali = cp_lp%i_ali_grd_p; c_ali = cp_lp%c_ali_grd_p; 
      g_off = cp_lp%i_off_grd_p; c_off = cp_lp%c_off_grd_p; 
      gpl = g_ali(4,1); opl = g_off(4,1)
!
      tp1 = zero
      do j = 1, gpl
        tp1 = tp1 + sga(g_ali(j,1),g_ali(j,2))*c_ali(j+1)
      enddo
      tp1 = tp1 + cp_lp%sgm*c_ali(1)
      tp2 = zero
      do j = 1, opl
        tp2 = tp2 + sga(g_off(j,1),g_off(j,2))*c_off(j+1)
      enddo
      tp2 = tp2 + cp_lp%sgm*c_off(1)
      cp_lp%gsgm = tp1+tp2
!
!!simplified      smvec(i) = tp1+tp2
      smvec(i) = sik
!
      cp_lp => cp_lp%next
    enddo
    call smoothcrossing(olen, nib,il,mkr,iaary,llen, smvec,spvec, info)
!!prt    tp1 = paraval(half, nring, mkr, il, 0, 0, info) ! on + side old network volume
!!prt    tp2 = paraval(half, nring, mkr, il, 1, 0, info) ! on + side old network volume
!!prt    print '("REAL  ", 1024(e14.6,1x))', tp1, tp2, s0, siga(s0,half,half,0), mkr(:,1,il)
  enddo
!========================================================================
!!  strlen = len_trim(runname)
!!  write(ibfile,'(2a,i4.4)') runname(1:strlen),'.ib.',outcount
!!  open(67,file=ibfile,form='formatted',action='write')
!!  do i = 1, nring
!!    write(67,'(2(e22.14,1x))')xpt(i),ypt(i)
!!  enddo
!!  write(67,'(2(e22.14,1x))')xpt(1),ypt(1)
!!  close(67)
    nt = nint(time/dt)
    strlen = len_trim(runname)
  if (nt-nfreq*(nt/nfreq) .eq. 0) then
!!notused    write(efile,'(2a,i4.4)') runname(1:strlen),'.sa.',outcount+1
!!notused    open(77,file=efile,access='stream',action='write')
!!notused    write(77)sga
!!notused    close(77)
!!    write(ffile,'(2a,i4.4)') runname(1:strlen),'.ns.',outcount+1
!!    open(67,file=ffile,access='stream',action='write')
!!    write(67)ntpn
!!    close(67)
  endif
  
!========================================================================
!
  return
end subroutine sg_IBbdy
!p
!========================================================================

subroutine LapSolver(rg,cg, nl, iaary, llen, info)
!!  double precision, dimension(-1:nx+2,-1:ny+2) :: rg
  implicit none
! solve the linear system for A*xc = rs, where xc is chemical variable only
  double precision :: time, dt, mydif
  double precision, dimension(-1:nxb+2,-1:nyb+2) :: rg, cg
  integer :: nl, info, llen(nent-2), isel
  type(iapt), dimension(:) :: iaary(nent-2)
!
  integer, parameter :: mm=100
  integer :: ierr, ipl, ip(2), im(2), iloc(4,2), ict
  integer :: i, j, k, l, il, iter, it, mit, infon, mi, j1, lenj, lenw
  double precision :: s0, theta, dthe, tp1, tp2, delt, st_p(4)
  double precision :: con, cox, coy, tmp, normb, dr, omega, beta, eps, resid
  double precision :: ser0, ser1, serr, dr1, hit, st, error, tol, dx, dy
  double precision, dimension(nxb*nyb) :: s, vk, x11, bf, r, yy, w, xx
  double precision, dimension(nxb*nyb) :: ph, vv, rt, sv, pv, tv
!
  double precision :: v(nxb*nyb,mm),  hj(nxb*nyb,2), hq(nxb*nyb,mm)
!========================================================================
  double precision :: utx, vtx
!!  integer :: ctl, ctr, cbl, cbr, idx(1:nxb*nyb)
!!  integer, dimension(nx-2) :: cbt, ctp
!========================================================================
! 
  double precision, dimension(nx*ny) :: rhsv, xcv, prem
!
  type(cc_augvar), pointer :: cp_lp
  !!character(40) ibfile  , ffile
!
  dx = hg; dy = hg
  mit = nxb*nyb
  l   = mit; lenj = nxb*nyb; lenw = lenj
!
  v = 0.; hj = 0.; hq = 0.; serr = 0.
!========================================================================
! set RHS for Dirichlet BC 
!========================================================================
  !!cox = dif(1)*dlt/dx/dx; coy = dif(1)*dlt/dy/dy ! backward Euler
!!  con = one/dx/dy; cox = con; coy = con 
!
!========================================================================
!  setting all - side pts (including auxiliary variables) to 0 in index idz
!========================================================================
!!  jgd = reshape(iid(1:nx,1:ny), (/ nx*ny /))
!========================================================================
! modify top/bottom & 4 corners for Dirichlet BCs
!========================================================================
!========================================================================
  rhsv(1:lenj) = reshape(rg(1:nxb,1:nyb), (/ nxb*nyb /))
!!
  normb=sqrt(dot_product(rhsv,rhsv))
!!  print *, 'RHS norm is ', normb
!
  !!eps = 1.d-11*normb/2.**(2.*dble(l2ny-5)); resid=eps; tol = eps;
  eps = 1.d-9*normb; resid=eps; tol = eps;
  if (abs(normb)< eps) then
    xx=0.d0
    iter = 1
    resid = 0.
!!    cg = rg
    print *, 'NO GMRES, iter', iter, normb
    return
  endif
!========================================================================
!!  rg(1:nx,1:ny) = rg(1:nx,1:ny)*dble(iid(1:nx,1:ny))
!========================================================================
!=======================================================================
! construct Jacobian predconditioner for GMRS
!=======================================================================
  !!tp2 = 1./(2.*(cox+coy))
  !!prem(1:lenj) = tp2
!
  prem(1:lenj) = 0.25d0
!=======================================================================
! GMRES
!=======================================================================
  !xc = 1.+0.001*tp2
  xcv = .0
  xcv(1:lenj) = reshape(rg(1:nxb,1:nyb), (/ nxb*nyb /))
  iter =0
  v = 0.; hq = 0.; hj = 0.; serr = 0.
  do j = 1, mit
    ser0 = 0.
!!    call amatv(xcv, yy, nl, iaary, llen, info)
    call nmatv(xcv, yy, iaary, llen, info)
!!        print *, 'check 1',maxval(yy(1:lenj)), minval(yy), maxval(xcv), minval(xcv)
!
    r =  rhsv - yy
    r = r*prem !PRECONDITIONER
    dr1 = sqrt(dot_product(r,r))
    if(dr1 .le. 1d-14) then
      print *, 'guess solved EQN already'
      return
    endif
    
    v(:,1) = r/dr1
    s = 0.
    s(1) = dr1
      do i = 1, mm-1
        iter = iter + 1
        vk = v(:,i)

!!        call amatv(vk,w, nl, iaary, llen, info)
        call nmatv(vk, w, iaary, llen, info)
!!        print *, 'check 2',maxval(w(1:lenj)), minval(w), maxval(vk), minval(vk)
        w = w*prem !PRECONDITIONER
        do k = 1, i
          vk = v(:,k)
          hq(k,i) = dot_product(w,vk)
          w = w - hq(k,i)*vk
        enddo

        hq(i+1,i) = sqrt(dot_product(w,w))
        if (abs(hq(i+1,i)) < 1d-14) then
          infon = 1
        else
          infon = 0
          v(:, i+1) = w/hq(i+1,i)
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
         x11 = x11 + yy(k)*v(:,k)
      enddo
!!      call amatv(x11,w, nl, iaary, llen, info)
      call nmatv(x11,w, iaary, llen, info)
      r = rhsv - w
      r = r*prem !PRECONDITIONER
      s(mi+1) = sqrt(dot_product(r,r))
      xcv = x11
      if( abs(s(mi+1)) .lt. tol .or. serr .lt. tol ) then
        ! this should be the return point 
        print '(1x,"IN   poisson  , iter", i5,1x, 3(e14.6,1x),3(i5,1x))', iter, &    
            serr,normb, s(mi+1),qlen, mi
!!	print *, mi, s(mi), s(mi+1), maxval(ra), maxval(rg), maxval(prem), s(1)
        cg(1:nxb,1:nyb) = (reshape(xcv(1:nxb*nyb), (/ nxb, nyb /)))
!!	print *, 'result ', maxval(cg), maxloc(cg), ijd(45,22), iid(45,22), id(45,22), idf(45,22)
!!        deallocate(rhsv,xcv,prem)
        return
      endif
      error = s(mi+1)
  enddo
  print *, 'ENDGMRES w/ MAX iter', iter, serr
!=======================================================================
! end of GMRES
!=======================================================================
!
!!  deallocate(rhsv,xcv,prem)
!
  return
end subroutine LapSolver
!
!=======================================================================
!=======================================================================
subroutine nmatv(xx, Ax, iaary, llen, info)
! here the vector is reshaped to 1d, in the order of storing each 
! row, which is
! {1(1:...:nxb) ... j(1:.i.:nxb) ... nyb(1:...:nxb)}
  integer :: llen(nent-2), info, nl
  double precision :: time, dt, mydif
  double precision, dimension(1:nxb*nyb) :: x, xx, Ax
  type(iapt), dimension(:) :: iaary(nent-2)
!
  integer :: il, i, j, k, l, ic, g_ali(4,2), g_off(4,2), cpl
  integer :: ik, jk, jm(2), jp(2), is(2), jc(2), iloc(4,2), gm(2)
  integer :: imv(-1:nxb+2,-1:nyb+2)
  double precision :: tmp, coef, con, cox, coy, delta, c_ali(4), c_off(4), st(4)
  type(cc_augvar), pointer :: cp_lp
!
  double precision :: x1, y1, fm(2)
  double precision :: ptop(nxb)
! 
  con = one/hg/hg; cox = con; coy = con;
!!  con = one; cox = con; coy = con;
  Ax = 0.
!!  x = xx*jgd
  x = xx
!
!!  x = x*jgd
!========================================================================
  Ax(cip) = 4.*con*x(cip)-cox*x(cip+1)-cox*x(cip-1)-coy*x(cip+nxb)-coy*x(cip-nxb)
!========================================================================
!!! bottom: no flux
  Ax(cbt) = (4.*con)*x(cbt)-cox*x(cbt+1)-cox*x(cbt-1)-coy*x(cbt+nxb)-coy*x(cbt)
!!!!! top no flux
  Ax(ctp) = (4.*con)*x(ctp)-cox*x(ctp+1)-cox*x(ctp-1)-coy*x(ctp)-coy*x(ctp-nxb)
!!! left periodic
!!  Ax(clf) = (4.*con)*x(clf)-cox*x(clf+1)-cox*x(crt  )-coy*x(clf+nxb)-coy*x(clf-nxb)
! left noflux
  Ax(clf) = (4.*con)*x(clf)-cox*x(clf+1)-cox*x(clf  )-coy*x(clf+nxb)-coy*x(clf-nxb)
!!! right periodic
!!  Ax(crt) = (4.*con)*x(crt)-cox*x(clf  )-cox*x(crt-1)-coy*x(crt+nxb)-coy*x(crt-nxb)
! right noflux
  Ax(crt) = (4.*con)*x(crt)-cox*x(crt  )-cox*x(crt-1)-coy*x(crt+nxb)-coy*x(crt-nxb)
!========================================================================
! bottom left corner: noflux/periodic BC
  Ax(cbl) = (4.*con)*x(cbl)-cox*x(cbl+1)-cox*x(cbr  )-coy*x(cbl+nxb)-coy*x(cbl)
!!!!!right, bot: Dirichlet
  Ax(cbr) = (4.*con)*x(cbr)-cox*x(cbl  )-cox*x(cbr-1)-coy*x(cbr+nxb)-coy*x(cbr)
!!!!top left conner: noflux/periodic
  Ax(ctl) = (4.*con)*x(ctl)-cox*x(ctl+1)-cox*x(ctr  )-coy*x(ctl)-coy*x(ctl-nxb)
!!!!!right, top: no flux
  Ax(ctr) = (4.*con)*x(ctr)-cox*x(ctl  )-cox*x(ctr-1)-coy*x(ctr)-coy*x(ctr-nxb)
!=======================================================================
!
  imv = 0
  do il = 1, nib
    cp_lp => iaary(il)%p
    do i = 1, llen(il)
!========================================================================
!========================================================================
      jm = 0; jp = 0; iloc = 0; cpl = 0; st = 0.; fm = 0.
      jm = cp_lp%im; jp = cp_lp%ip; cpl = cp_lp%ipl; gm = jm-jp; fm = dble(gm)
      iloc = cp_lp%ilocp; st = cp_lp%st_p
      call sub2ind(nxb, jp(1), jp(2), k) 
      ic = k
      call sub2ind(nxb, jm(1), jm(2), j) ! on the "-" side !!
      Ax(k) = Ax(k) + con*x(j) ! add contribution from ghost point back
      tmp= zero
      do j = 1, cpl
        call sub2ind(nxb, iloc(j,1), iloc(j,2), k) 
        tmp = tmp + x(k)*st(j+1)
      enddo
      tmp =-tmp*con

      Ax(ic) = Ax(ic) + tmp
!========================================================================
!
      cp_lp => cp_lp%next
    enddo
  enddo

  Ax = Ax*dble(igd) + (one-dble(igd))*xx ! debugging mode
!!  Ax = Ax*jgd
!
!=======================================================================
!
  return
end subroutine nmatv
!
!!!=======================================================================
!!subroutine f2m(x, M, nl, iaary, llen, info)
!! here this should go to the RHS of the Laplace Eqn solver
!!!
!!  double precision :: dt, mydif
!!  double precision :: M(1:nxb*nyb), vnt(1:nx*ny)
!!  type(iapt), dimension(:) :: iaary(nent-2)
!!  integer :: nl, info, llen(nent-2)
!!  double precision, dimension(nl) :: x
!!!
!!  integer :: i, j, k, il, ic, im(2), ip(2), cpl, iloc(4,2), gm(2)
!!  double precision :: tmp, tp1, coef, delta, theta, dthe, s0, st(4), fm(2)
!!  type(cc_augvar), pointer :: cp_lp
!!!
!!  coef =-one/hg/hg  ! for network volume
!!!
!!  M = 0.
!!  do il = 1, nib
!!    cp_lp => iaary(il)%p
!!    ic = 0
!!    do i = 1, llen(il)
!!      ic = ic + 1
!!      ip = cp_lp%ip
!!!
!!      call sub2ind(nxb,ip(1),ip(2),k) ! find the irregular pt for the auxiliary variable
!!      
!!      M(k) = M(k)+ coef*cp_lp%st_p(1)*x(ic)
!!!!      print '(1x, "M(-)= ", e16.8, "k: ", 2(i5,1x), 10(e16.8,1x))', M(k), k, ic, s0, theta, tp1, vnt(k), tmp, dthe, x(ic), tp1
!!      cp_lp => cp_lp%next
!!    enddo
!!  enddo
!!  if (ic .ne. nl) then
!!    print *, 'error in data processing between linked list and array'
!!    stop
!!  endif
!!!R  print *, 'IC= ', ic
!!!
!!  return
!!end subroutine f2m
!!
!!!=======================================================================

double precision function siga(s,x,y,isel)
! assign sigma_a along boundary
  double precision :: s,x,y
  integer :: isel
!
  !!siga = -rho0*theta0*0.01
  !!siga = zero
!!implicit none
!!  double precision :: s, time
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

  s0 = mod(mod(s, tupi)+tupi, tupi) ! archlength we will used

! if both shs and sts are -, we have shrinking
!  shs=-.10;      sts= 1.0*shs
! if both shs and sts are +, we have swelling
!  shs= .10;      sts= 1.0*shs
! if  shs and sts have opposite -/+, we move the cell
  shs=-.00001;      sts=-2.0*shs
  sus=-shs*0.0;  sbs= 1.0*sus
  shl=cpi*.0 ;   stl=cpi*1.0
  sul=cpi*.80;   sbl=cpi*1.20
  shw=cpi*.125;  stw=shw*0.5;
  suw=shw*.5 ;   sbw=shw*.5

!!o  tmp=shs*exp(-.5*(s-shl)**2./shw**2.)+0. *exp(-.5*(s-shl-tupi)**2./shw**2.)+&
!!o      sts*exp(-.5*(s-stl)**2./stw**2.)
!working  tmp=shs*exp(-.5*(s-shl)**2./shw**2.)+shs*exp(-.5*(s-shl-tupi)**2./shw**2.)+&
!working      sts*exp(-.5*(s-stl)**2./stw**2.)
  tmp=shs*exp(-.5*(s-shl)**2./shw**2.)+shs*exp(-.5*(s-shl-tupi)**2./shw**2.)+&
      sts*exp(-.5*(s-stl)**2./stw**2.)+ &
      sus*exp(-.5*(s-sul)**2./suw**2.)+sbs*exp(-.5*(s-sbl)**2./sbw**2.)
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
  !!tmp = 0.001
!!  siga =tmp
  !!cc_gate = 0.1
!
!
  if (isel .eq. 1) then
    !!siga = 0.001d0*x
    siga = sik*x
  else if (isel .eq. 0) then
    !!siga = -rho0*theta0*0.1d0
    siga = -rho0*theta0*0.1d0
  !siga = -rho0*0.1d0*(one+dcos(s)*dcos(s))
  else 
    print *, 'need choose 0 or 1 for sigma_a along cell membrane!'
    stop
  endif
!
  return
end function siga

!!deb!========================================================================
!!deb! gmres solver with mat-vec product done with procedure pointer(suppoted
!!deb! in Fortran 2003 (gfortran 4.5 and newer)
!!deb!========================================================================
!!debsubroutine gmres(x, Ax, nl, ux, uy, iaary, llen, dt, mydif, info, mvp)
!!deb  integer :: llen(nent-2), info, nl
!!deb  double precision, dimension(1:nx*ny) :: ux, uy, rhsv
!!deb  double precision :: time, dt, mydif
!!deb  double precision, dimension(1:nx*ny+nl) :: x, Ax
!!deb  type(iapt), dimension(:) :: iaary(nent-2)
!!deb
!!deb  procedure(matv), pointer :: mvp
!!deb
!!deb
!!deb
!!deb  return
!!debend subroutine gmres
!========================================================================

end module solver_mod
