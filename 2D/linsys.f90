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

!include subroutines for solving Stokes equation, reading parameters and 
! save (u,v,p) data
module linsys
  use, intrinsic :: iso_c_binding
  use parameters
!!  use IBmod
!!  use IBforce
  use myfft

  implicit none
  public :: getRHS, linmatk, sollinsys, rsoltouvp, asol, inituvp, cmpuvp, &
       readpar, outuvpx, setbc
!
contains
!--------------------------------------------------------------         
!
  subroutine readpar(CFL,nu,bk,vco,ntmax, nfreq, dlt)
  integer :: ntmax, nfreq, isel
  double precision :: CFL,nu,bk,vco,dlt
  namelist/PARAM/CFL,nu,bk,vco,ntmax,nfreq,dlt
!
  open(unit=99, file = 'input.par',status='old',form='formatted')
  read(99,nml=PARAM)
  close(unit=99)
!
  return
  end subroutine readpar
!
!========================================================================
!========================================================================
subroutine setbc(rho,uvbc,isel)
! periodic on left/right, Dirichlet on top/bottom
!
  integer :: isel
  double precision, dimension(-1:nx+1,4) :: uvbc !ubt, utp, vbt, vtp
  double precision, dimension(-1:nx+1,-1:ny+1), intent(in out) :: rho
!
  integer :: i, j, k
  double precision :: tmp
  double precision, dimension(-1:nx+1) :: top, bot
!
  select case (isel)
  case (1,-1) !u
    if (isel .eq. 1) then
      bot(0:nx-1) = uvbc(0:nx-1,1); top(0:nx-1) = uvbc(0:nx-1,2)
    else ! homogeneous Dirichlet
      bot(0:nx-1) = 0.d0     ; top(0:nx-1) = 0.d0     
    endif
    j = 0
    rho(0:nx-1,j) = two*bot(0:nx-1) - rho(0:nx-1,j+1)
    j = -1
    rho(0:nx-1,j) = two*bot(0:nx-1) - rho(0:nx-1,j+2)
    j = ny+1
    rho(0:nx-1,j) = two*top(0:nx-1) - rho(0:nx-1,j-1)
    !
    i = -1
    rho(i,0:ny+1) = rho(nx+i,0:ny+1)
    i = nx
    rho(i,0:ny+1) = rho(i-nx,0:ny+1)
!!    i = nx+1
!!    rho(i,0:ny+1) = rho(i-nx,0:ny+1)
  case (2,-2) !v
    if (isel .eq. 2) then
      bot(1:nx) = uvbc(1:nx,3); top(1:nx) = uvbc(1:nx,4)
    else ! homogeneous Dirichlet
      bot(1:nx) = 0.d0     ; top(1:nx) = 0.d0     
    endif
    j = 0 ! bottom Bdry
    rho(1:nx,j) = bot(1:nx);
    j = -1
    rho(1:nx,j) = two*bot(1:nx) - rho(1:nx,j+2)
    j = ny ! top Bdry
    rho(1:nx,j) = top(1:nx)
    j = ny+1
    rho(1:nx,j) = two*top(1:nx) - rho(1:nx,j-2)
    !
    i = 0
    rho(i,-1:ny+1) = rho(nx+i,-1:ny+1)
!!    i = -1
!!    rho(i,-1:ny+1) = rho(nx+i,-1:ny+1)
    i = nx+1
    rho(i,-1:ny+1) = rho(i-nx,-1:ny+1)
  case (3) !p
    i = 0
    rho(i, 1:ny) = rho(nx+i, 1:ny)
    i = nx+1
    rho(i, 1:ny) = rho(i-nx, 1:ny)
  case default
    print *,'Can only choose (-)1 (u), or (-)2(v)! stop'
    stop
  end select 
!
  return
end subroutine setbc
!========================================================================
subroutine lapop(lap, rho, isel)
! periodic BC on left/right, Homogeneous Dirichlet on top/bottom 
! need advu(0:nx-1,1:ny), advv(1:nx,1:ny-1)
!
  integer :: isel
  double precision, dimension(-1:nx+1,-1:ny+1) :: lap, rho
  double precision, dimension(-1:nx+1,4) :: uvbc !ubt, utp, vbt, vtp
!
  integer :: i, j, istr, iend, jstr, jend
  double precision :: tp1, tp2, tmp, h, dx, dy, rhh
!
  rhh = one/(hg*hg) 
  lap = 0.d0
!
  select case (isel) 
   ! u component                                            
  case (0) 
    istr = 0 
    iend = nx-1 
    jstr = 1 
    jend = ny
   ! v component                                            
  case (1) 
    istr = 1 
    iend = nx
    jstr = 1 
    jend = ny-1 
   ! rho @ cell centers                                     
  case (2) 
    istr = 1 
    iend = nx
    jstr = 1 
    jend = ny
  case default 
    print *, isel 
    print *, 'no such option in lap, stop' 
    stop 
  end select 
  do j = jstr, jend 
  do i = istr, iend 
    lap(i,j) = (rho(i+1,j)+rho(i-1,j)-two*rho(i,j))          &
     &       + (rho(i,j+1)+rho(i,j-1)-two*rho(i,j))          
  enddo 
  enddo 
  lap = rhh*lap
!
  return
end subroutine lapop
!========================================================================
subroutine getAdv(advu, advv, u, v, uvbc)
! periodic BC on left/right, Homogeneous Dirichlet on top/bottom 
! need advu(0:nx-1,1:ny), advv(1:nx,1:ny-1)
!
  double precision, dimension(-1:nx+1,-1:ny+1) :: u, v
  double precision, dimension(-1:nx+1,-1:ny+1), intent(out) :: advu, advv
  double precision, dimension(-1:nx+1,4) :: uvbc !ubt, utp, vbt, vtp
!
  integer :: i, j, k
  double precision :: tp1, tp2, tmp, h, dx, dy
  double precision, dimension(-1:nx+1,-1:ny+1) :: ua, va
!
  h = hg; dx = h; dy = h
  advu = 0.d0; advv = 0.d0;
! 
!!  call setbc(u,uvbc,1) ! ghost of u CHANGED here
!!  call setbc(v,uvbc,2) ! ghost of v CHANGED here
!
  do i = 0,nx-1 ! average of v @ u locations
    do j = 1, ny
      va(i,j) = 0.25d0*(v(i,j-1) + v(i+1,j-1) + v(i,j) + v(i+1,j))
    enddo
  enddo
  call setbc(va,uvbc,1) ! need to be modified to account for the shifted BC loc
  do i = 1, nx ! average of u @ v locations
    do j = 1, ny-1
      ua(i,j) = 0.25d0*(u(i-1,j) + u(i,j) + u(i-1,j+1) + u(i,j+1))
    enddo
  enddo
  call setbc(ua,uvbc,2) ! 
!
  do i = 0, nx-1
    do j = 1, ny
      advu(i,j) = u(i,j)*(u(i+1,j)-u(i-1,j))+va(i,j)*(u(i,j+1)-u(i,j-1)) + & 
                  (u(i+1,j)*u(i+1,j)  - u(i-1,j)*u(i-1,j) ) +  &
                  (u(i,j+1)*va(i,j+1) - u(i,j-1)*va(i,j-1))
    enddo
  enddo
  advu = advu*half*half/h
  do i = 1, nx
    do j = 1, ny-1
      advv(i,j) = ua(i,j)*(v(i+1,j)-v(i-1,j))+v(i,j)*(v(i,j+1)-v(i,j-1)) + &
                  (ua(i+1,j)*v(i+1,j) - ua(i-1,j)*v(i-1,j)) +  &
                  (v(i,j+1)*v(i,j+1)  - v(i,j-1)*v(i,j-1))
    enddo
  enddo
  advv = advv*half*half/h
!
  return
end subroutine getAdv
!========================================================================
subroutine getRHS(rsu, rsv, u, v, p, f, g, uvbc, uold, vold, ins, iadv, dt)
! get RHS for NS solver
  implicit none
!!  include '/Users/yaol/local/fftw/include/fftw3.f'
! do FFT on u & v to get RHS.
!
  complex*16, dimension(-1:nx+1,-1:ny+1), intent(in out) :: rsu, rsv
  double precision, dimension(-1:nx+1,-1:ny+1) :: u, v, p, f, g, uold, vold
  double precision, dimension(-1:nx+1,4) :: uvbc !ubt, utp, vbt, vtp
  integer :: ins, iadv
  double precision :: dt
!!deb  double precision, dimension(-1:nx+1,-1:ny+1) :: p0, ur, vr, fr, gr
!
  double precision :: dx, dy, x, y, eu, mu, gam, beta
  double precision, dimension(-1:nx+1,-1:ny+1) :: advu, advv, lapu, lapv
  complex*16, dimension(1:nx) :: din, dout
  integer :: i, j, nl, iflg
  integer*8 :: plan
!!  double precision :: sol
!!  external sol
! 
!!  h = hg; 
  dx = hg; dy = hg; nl = nx; 
!
  advu = 0.d0;  advv = 0.d0
  if ( iadv .eq. 1) then
    call getAdv(advu, advv, u, v, uvbc) 
  endif
  select case (ins)
  case (0)
    advu = dt*f;
    advv = dt*g;
  case (1)! Implicit Euler
    advu = -dt*advu+ u + dt*f;
    advv = -dt*advv+ v + dt*g;
    beta = dt*nu/hg/hg
  case (2)! CN
    call lapop(lapu,u,0)
    call lapop(lapv,v,1)
    !
    beta = half*dt*nu/hg/hg
    gam=half*dt/hg
    !
    advu = -dt*advu+ u + dt*f +nu*half*dt*lapu;
    advu(0:nx-1,1:ny) = advu(0:nx-1,1:ny)-gam*(p(1:nx,1:ny)-p(0:nx-1,1:ny))
    advv = -dt*advv+ v + dt*g +nu*half*dt*lapv;
    advv(1:nx,1:ny-1) = advv(1:nx,1:ny-1)-gam*(p(1:nx,2:ny)-p(1:nx,1:ny-1))
  case default
    print *, 'No such option for RHS of linear system. Stop!'
    stop
  end select
!
  do j = 1, ny
    din = advu(0:nx-1,j)
    call dfftw_plan_dft_1d(plan, nl, din,dout,FFTW_FORWARD,FFTW_ESTIMATE)
    call dfftw_execute_dft(plan,din,dout)
    rsu(0:nx-1,j) = dout(1:nx)
    call dfftw_destroy_plan(plan)
  enddo
!========================================================================
! modify the RHS due to Dirichlet BC from u-component @ bottom
!!  din = uvbc(0:nx-1,1)
!!  call dfftw_plan_dft_1d(plan, nl, din,dout,FFTW_FORWARD,FFTW_ESTIMATE)
!!  call dfftw_execute_dft(plan,din,dout)
!!  rsu(0:nx-1,1) = rsu(0:nx-1,1)+two*beta*dout(1:nx)
!!  call dfftw_destroy_plan(plan)
!!! modify the RHS due to Dirichlet BC from u-component @ top
!!  din = uvbc(0:nx-1,2)
!!  call dfftw_plan_dft_1d(plan, nl, din,dout,FFTW_FORWARD,FFTW_ESTIMATE)
!!  call dfftw_execute_dft(plan,din,dout)
!!  rsu(0:nx-1,ny) = rsu(0:nx-1,ny)+two*beta*dout(1:nx)
!!  call dfftw_destroy_plan(plan)
!========================================================================
  ! FFT along each row (j=1:ny-1) for v-component
  do j = 1, ny-1
    !!din = 0.d0; 
    din = advv(1:nx,j)
    call dfftw_plan_dft_1d(plan, nl, din,dout,FFTW_FORWARD,FFTW_ESTIMATE)
    call dfftw_execute_dft(plan,din,dout)
    rsv(1:nx,j) = dout(1:nx)
    call dfftw_destroy_plan(plan)
  enddo

!!  do i = 1, ny
!!    print'(1x,i3, 1024("(",e11.3,1x,e11.3,")"))', i, rsu(0:nx-1,i)
!!  enddo
!!  do i = 1, ny-1
!!    print'(1x,i3, 1024("(",e11.3,1x,e11.3,")"))', i, rsv(1:nx,i)
!!  enddo
!!  print *, 'FFT of data'
!
  return
end subroutine getRHS

!--------------------------------------------------------------         
!--------------------------------------------------------------         
subroutine linmatk(AB, k, alpha, beta, gam)
  implicit none
  !!include '/Users/yaol/local/fftw/include/fftw3.f'
  integer :: isel, k ! k is the wave number
  double precision :: alpha, beta, gam
  ! LDAB = 2*3+3+1 = 10
  complex*16, dimension(10,3*ny-1) :: AB
!
  integer :: i, j, jj, ik, jk, kl, ku, nl
  double precision :: tpikh, pikh, cstpikh, sitpikh, cspikh, sipikh
  double precision :: h, alpbet, alpbet1
  complex*16 :: Epikh, Etpikh, Etpikh1, Etpikhp
  complex*16, dimension(3) :: l1, l2, l3, u1, u2, u3
  complex*16, dimension(7) :: c1, c2, c3
!
  nl = 3*ny-1;  h = hg
  tpikh = 2.d0*cpi*dble(k)*hg;
  cstpikh = dcos(tpikh); sitpikh = dsin(tpikh)
!
  alpbet = alpha-2.d0*beta*cstpikh
  alpbet1= alpha-2.d0*beta*cstpikh+beta
!
  Etpikh1= (cmplx(cstpikh-1.d0,sitpikh))*gam
  Etpikhp= (cmplx(1.d0-cstpikh,sitpikh))*gam
!
  AB = cmplx(zero,zero)
!
  u1(1) = 0.d0;     u1(2) = 0.d0;  u1(3) =-gam
  l1(1) = Etpikhp;  l1(2) = gam;   l1(3) = 0.d0
  u2(1) =-beta;     u2(2) = 0.d0;  u2(3) = Etpikh1
  l2(1) = 0.d0;     l2(2) = 0.d0;  l2(3) =-beta
  u3(1) =-beta;     u3(2) =-gam  ; u3(3) = 0.d0
  l3(1) = gam  ;    l3(2) = 0.d0;  l3(3) =-beta
  c1(1:3) = u1; c1(4) = cmplx(0.d0,0.d0);   c1(5:7) = l1
  c2(1:3) = u2; c2(4) = cmplx(alpbet,0.d0); c2(5:7) = l2
  c3(1:3) = u3; c3(4) = cmplx(alpbet,0.d0); c3(5:7) = l3
!
  kl = 3; ku = 3
!
  jj = 1;
  j  = 3*(jj-1)+1; i = 1; ik = kl + ku + 1 + i -j
  AB(ik:ik+3,j) = c1(4:7)
  if (k .eq. 0) then
    AB(ik,j) = cmplx(1.0d0,0.d0)
  endif
  j  = 3*(jj-1)+2; i = 1; ik = kl + ku + 1 + i -j
  AB(ik:ik+4,j) = c2(3:7)
  AB(ik+1,j) = cmplx(alpbet1,0.d0)
  j  = 3*(jj-1)+3; i = 1; ik = kl + ku + 1 + i -j
  AB(ik:ik+5,j)= c3(2:7)
  do jj = 2, ny-2
    j = 3*(jj-1)+1; i = j-3; ik = kl + ku + 1 + i -j
    AB(ik:ik+6,j) = c1
    j = 3*(jj-1)+2; i = j-3; ik = kl + ku + 1 + i -j
    AB(ik:ik+6,j) = c2
    j = 3*(jj-1)+3; i = j-3; ik = kl + ku + 1 + i -j
    AB(ik:ik+6,j) = c3
  enddo
  jj= ny-1
  j = 3*(jj-1)+1; i = j-3; ik = kl + ku + 1 + i -j
  AB(ik:ik+6,j) = c1
  j = 3*(jj-1)+2; i = j-3; ik = kl + ku + 1 + i -j
  AB(ik:ik+6,j) = c2
  j = 3*(jj-1)+3; i = j-3; ik = kl + ku + 1 + i -j
  AB(ik:ik+5,j) = c3(1:6)
!
  jj= ny
  j = 3*(jj-1)+1; i = j-3; ik = kl + ku + 1 + i -j 
  AB(ik:ik+4,j) = c1(1:5)
  j = 3*(jj-1)+2; i = j-3; ik = kl + ku + 1 + i -j 
  AB(ik:ik+2,j) = c2(1:3)
  AB(ik+3,j) = cmplx(alpbet1,0.d0)
!
  !!prt if (k .eq. 0) then
  !!prt print *, 'Coefficient Matrix start'
  !!prt do i = 1, nl
  !!prt   !!print '(i3, 1024("(",e8.1,1x,e8.1,")"))', i, AB(4:10,i)
  !!prt   print '(i3, 1024(e9.2,1x))', i, real(AB(4:10,i))
  !!prt enddo
  !!prt print *, 'Coefficient Matrix end'
  !!prt endif
!
  return
end subroutine linmatk
!--------------------------------------------------------------         
subroutine rsoltouvp(rsol, u, v, p)
  double precision, dimension(3*ny-1,nx) :: rsol
  double precision, dimension(-1:nx+1,-1:ny+1) :: u,v,p
  integer ::j, jj, k
!
!!  u = 0.; v= 0.; p = 0.
!
  !!do k = 1, nx
  !!  do jj = 1, ny-1
  !!    j = 3*(jj-1)+1
  !!    p(k,  ny-jj+1) = rsol(j,  k)
  !!    u(k-1,ny-jj+1) = rsol(j+1,k)
  !!    v(k,  ny-jj) = rsol(j+2,k)
  !!  enddo
  !!  jj = ny
  !!  j = 3*(jj-1)+1
  !!  p(k,  ny-jj+1) = rsol(j,  k)
  !!  u(k-1,ny-jj+1) = rsol(j+1,k)
  !!enddo
  do k = 1, nx
    do jj = 1, ny-1
      j = 3*(jj-1)+1
      p(k,  ny-jj+1) = rsol(j,  k)
      u(k-1,ny-jj+1) = rsol(j+1,k)
      v(k,  ny-jj) = rsol(j+2,k)
    enddo
    jj = ny
    j = 3*(jj-1)+1
    p(k,  ny-jj+1) = rsol(j,  k)
    u(k-1,ny-jj+1) = rsol(j+1,k)
  enddo
!
  return
end subroutine rsoltouvp


subroutine sollinsys(rsu,rsv,resul,rsol, alpha, beta, gam, dt)
  implicit none
  complex*16, dimension(-1:nx+1,-1:ny+1), intent(in) :: rsu, rsv
  double precision, dimension(3*ny-1,nx), intent(out) :: rsol
  double precision :: dt
!
  integer :: i, j, k, iflg, info, nl
  double precision :: dx, tupih, pih, h
  double precision :: alpha, beta, gam
  complex*16, dimension(10,3*ny-1) :: AB, AFB
  complex*16, dimension(3*ny-1) :: rhs, solv, solw
  complex*16, dimension(3*ny-1,nx) :: resul, csol!!, rhsm
  integer :: ipiv(3*ny-1), nlen
  character :: equed
  double precision :: mvR(3*ny-1), mvC(3*ny-1), ferr(1), berr(1), rwork(6*ny-2)
  double precision :: rcond
  complex*16 :: work(6*ny-2)
  integer*8 :: plan
  complex*16, dimension(nx) :: din, dout

!
!!  h = hg
!!prt print *, '====='
!!prt   print*,alpha, beta, gam
!!prt print *, '====='
!
  nl = 3*ny-1
!
  !!do k = nx, 1, -1
  do k = 1, nx
    rhs = cmplx(zero,zero)
    do i = 1, ny-1
      j = 3*(i-1) + 1
      rhs(j  ) = cmplx(0.d0,0.d0)
      rhs(j+1) = rsu(k-1,ny-i+1)
      rhs(j+2) = rsv(k,  ny-i)
    enddo
    i = ny
    j = 3*(i-1)+1
    rhs(j) = cmplx(0.d0,0.d0)
    rhs(j+1) = rsu(k-1,ny-i+1)
!!      print '(1x,"about??? ", i5, 1024(e14.6,1x))', k, maxval(real(rhs)), &
!!      minval(real(rhs))
!
!!    rhsm(:,k) = rhs
!
    AB = zero
!
    !!call linmatk(AB, k-nx/2+1, alpha, beta, gam, 0)
    call linmatk(AB, k-1, alpha, beta, gam)
!!    print '(1x,"howabout AB ", 1024(e14.6,1x))', maxval(real(AB)), &
!!    minval(real(AB))
    info = 0
    ipiv = 0

    solw = rhs
    call zgbsv(nl, 3, 3, 1, AB, 10, ipiv, rhs, nl, info)
!!deb    AFB = AB
!!deb    equed='N'
!!deb    call ZGBTRF(nl, nl, 3, 3, AFB, 10, ipiv, info )
!!deb!!DEBUG
!!deb    IF( INFO.NE.0 ) THEN
!!deb      print '(1x,1024(e14.6,1x))', AFB
!!deb      stop
!!deb    ENDIF
!!deb!!DEBUG
!!deb    mvR = zero; mvC = zero; rhs = solw
!!deb    call zgbsvx('F','N',nl,3,3,1,AB,10,AFB,nl,ipiv,&
!!deb       equed,mvR,mvC,rhs,nl,solv,nl,rcond,ferr,berr,work,rwork,info)
!!deb    rhs = solv
!!    if (k .eq. 2) then
!!      do i = 1, 3*ny-1
!!      print *, rhs(i)
!!      enddo
!!    endif
    if (info .eq. 0) then
      resul(1:nl,k) = rhs
!!      print '(1x,"get rhs? ", i5, 1024(e14.6,1x))', k, maxval(real(rhs)), &
!!      minval(real(rhs)), maxval(dimag(rhs)), maxval(dimag(rhs))
    else
      print *, nl
      print '("Error in wave stop ", 2(i5,1x) 10(e22.15,1x))', k, info, rcond
!!      print'(1x, 1024("(",e11.3,1x,e11.3,")"))', sum(AB(4:10,:),1)
      print'(1x, 1024(e11.3,1x))', (sum(AB(4:10,:),1))
      stop
    endif
  enddo
  rhs = zero
!

  rsol = 0.d0; csol = 0.d0
  nlen = nx
  do i = 1, 3*ny-1
    dout = resul(i,1:nx); 
!!    print '("dout ", i5,1024(e14.6,1x))',i, maxval(real(dout)), minval(real(dout))
    call dfftw_plan_dft_1d(plan, nlen, dout,din,FFTW_BACKWARD,FFTW_ESTIMATE)
    call dfftw_execute_dft(plan,dout,din)
    csol(i,:) = din/dble(nlen)
    rsol(i,:) = real(din/dble(nlen))
!!    print '("din  ", i5,1024(e14.6,1x))',i, maxval(real(din)), minval(real(din)),&
!!       maxval(dimag(csol)), minval(dimag(csol))
    call dfftw_destroy_plan(plan)
  enddo
!!prt  print *, norm2(real(csol)), norm2(imag(csol))
!!prt  do i = 1, 3*ny-1
!!prt  !!do i = 1, 2
!!prt    !!print'(1x, i3, 1024("(",e11.3,1x,e11.3,")"))', i, resul(i,:)
!!prt    !!print'(1x, i3, 1024("(",e11.3,1x,e11.3,")"))', i, csol(i,:)
!!prt    print'(1x, 1024(e11.3,1x))', rsol(i,:)
!!prt  enddo
!!  print *, 'rhs'
!!  do i = 1, 3*ny-1
!!  !!do i = 1, 2
!!    !!print'(1x, i3, 1024("(",e11.3,1x,e11.3,")"))', i, resul(i,:)
!!    print'(1x, i3, 1024(e11.3,1x))', i, real(resul(i,:))
!!    !!print'(1x, i3, 1024("(",e11.3,1x,e11.3,")"))', i, rhsm(i,:)
!!    !!print'(1x, i3, 1024(e11.3,1x))', i, real(rhsm(i,:))
!!  enddo
!
  return
end subroutine sollinsys
!--------------------------------------------------------------         
!!double precision function asol(x,y,t,isel)
subroutine asol(mysol, x,y,t,isel)
!!  use parameters
  implicit none
!
  double precision :: x, y, t, mysol
  integer :: isel
!
  tupi = two*cpi; pi = cpi
  select case (isel) ! first two rows are for steady Stokes
  case (0) !u
    mysol =-sin(tupi*y)*sin(pi*x)*sin(pi*x)*exp(-t)
    !!case4 mysol = cos(tupi*x)*(three*y*y-two*y)*exp(-t)
    !!case2 mysol =-sin(tupi*y)*sin(pi*x)*sin(pi*x)*exp(-t)
    !!case1 mysol = cos(tupi*x)*(three*y*y-two*y)*exp(-t)
    !!case3 mysol = 0.d0
  case (1) !v
    mysol = sin(tupi*x)*sin(pi*y)*sin(pi*y)*exp(-t)
    !!case4 mysol = tupi*sin(tupi*x)*y*y*(y-one)*exp(-t)
    !!case2 mysol = sin(tupi*x)*sin(pi*y)*sin(pi*y)*exp(-t)
    !!case1 mysol = tupi*sin(tupi*x)*y*y*(y-one)*exp(-t)
    !!case3 mysol = 0.d0
!!  case (-2) !px
!!    !!mysol =cos(tupi*x)*cos(tupi*y)
!!    mysol = tupi*pi*cos(tupi*x)*sin(tupi*y)
!!    !!mysol = 0.d0
!!  case (-3) !py
!!    !!mysol =cos(tupi*x)*cos(tupi*y)
!!    mysol = tupi*pi*sin(tupi*x)*cos(tupi*y)
!!    !!mysol = 0.d0
  case (2) !p
    mysol = pi*sin(tupi*x)*sin(tupi*y)*exp(-t)
    !!case4 mysol = nu*sin(tupi*x)*sin(tupi*y)*exp(-t)
    !!case2 mysol = pi*sin(tupi*x)*sin(tupi*y)*exp(-t)
    !!case1 mysol = nu*sin(tupi*x)*sin(tupi*y)*exp(-t)
    !!case3 mysol = 0.d0
  case (3) ! forcing on x
    !!mysol = two*pi*pi*(-nu+( one +two*nu)*cos(tupi*x))*sin(tupi*y)*exp(-t)
    mysol = half*(one-four*pi*pi*nu+(-one +pi*pi*four*( one+two*nu))*cos(tupi*x))*sin(tupi*y)*exp(-t)
    !!case4 mysol = two*nu*cos(tupi*x)*(-three+two*pi*pi*y*(-two +three*y)+pi*sin(tupi*y))*exp(-t)
    !!case2 mysol = half*(one-four*pi*pi*nu+(-one +pi*pi*four*( one+two*nu))*cos(tupi*x))*sin(tupi*y)*exp(-t)
    !!case1 mysol = exp(-t)*cos(tupi*x)*(-6.d0*nu+y*two*(one-four*pi*pi*nu)+three*y*y*(-one+tupi*tupi*nu)+tupi*nu*sin(tupi*y))
    !!case3 mysol = 0.d0
  case (4) ! forcing on y
    !!mysol =-two*pi*pi*(-nu+(-one +two*nu)*cos(tupi*y))*sin(tupi*x)*exp(-t)
    mysol =-half*(one-four*pi*pi*nu+(-one +pi*pi*four*(-one+two*nu))*cos(tupi*y))*sin(tupi*x)*exp(-t)
    !!case4 mysol = two*pi*nu*(two-6.d0*y-four*pi*pi*y*y+four*pi*pi*y*y*y+cos(tupi*y))*sin(tupi*x)*exp(-t)
    !!case2 mysol =-half*(one-four*pi*pi*nu+(-one +pi*pi*four*(-one+two*nu))*cos(tupi*y))*sin(tupi*x)*exp(-t)
    !!case1 mysol = exp(-t)*tupi*sin(tupi*x)*(two*nu-6.d0*y*nu+y*y*(one-four*pi*pi*nu)+y*y*y*(-one+tupi*tupi*nu)+nu*cos(tupi*y))
    !!case3 mysol = 0.d0
  case default
    print *, 'no this thing'
    stop
  end select
!
  return
END subroutine asol
!
subroutine inituvp(u0, v0, p0, f, g, uvbc, isel, dt, time)
  use parameters
!!  use linsys
  !use linsys
  implicit none
!
  double precision, dimension(-1:nx+1,-1:ny+1), intent(in out) :: u0, v0, p0, f, g
  double precision, dimension(-1:nx+1,4), intent(in out) :: uvbc !ubt, utp, vbt, vtp
  integer :: isel
  double precision :: time, dt
!
  integer :: i, j
  double precision :: x, y, dx, dy, realt
!!  double precision asol
!!  external asol
!
  dx = hg; dy = hg
  tupi = two*cpi; pi = cpi
!
!!  if (isel .eq. 0) then
  do j = 1, ny
    y = (dble(j)-0.5d0)*dy
    do i = 0, nx-1
      x = (dble(i))*dx
      !!u0(i,j) = asol(x,y,time,0)
      call asol(u0(i,j),x,y,time,0)
    enddo
  enddo
  do i = 0,nx-1
    x = (dble(i))*dx
    y = ymin
    !!uvbc(i,1) = asol(x,y,time,0)
    call asol(uvbc(i,1),x,y,time,0)
    y = ymax
    !!uvbc(i,2) = asol(x,y,time,0)
    call asol(uvbc(i,2),x,y,time,0)
  enddo
  do j = 1, ny-1
    y = (dble(j))*dy
    do i = 1, nx
      x = (dble(i)-0.5)*dx
      !!v0(i,j) = asol(x,y,time,1)
      call asol(v0(i,j),x,y,time,1)
    enddo
  enddo
  do i = 1, nx
    x = (dble(i)-0.5)*dx
    y = ymin
    call asol(uvbc(i,3),x,y,time,1)
    y = ymax
    call asol(uvbc(i,4),x,y,time,1)
  enddo
  do j = 1, ny
    y = (dble(j)-.5)*dy
    do i = 1, nx
      x = (dble(i)-0.5)*dx
      !!p0(i,j) = asol(x,y,time,2)
      call asol(p0(i,j),x,y,time,2)
    enddo
  enddo
!!  endif
  select case (isel)
  case (0)
    realt = time
  case (1)
    realt = time - 0.5*dt
!!    realt = time
  end select

    do j = 1, ny
      y = (dble(j)-0.5)*dy
      do i = 0, nx-1
        x = (dble(i)) *dx
        !!f(i,j) = asol(x,y,realt,3)
        call asol(f(i,j),x,y,realt,3)
      enddo
    enddo
    do j = 1, ny-1
      y = (dble(j))*dy
      do i = 1, nx
        x = (dble(i)-0.5)*dx
        !!g(i,j) = asol(x,y,realt,4)
        call asol(g(i,j),x,y,realt,4)
      enddo
    enddo
!!  endif
!!  call setbc(v0, p0, p0, time, isel) 
!
  return
end subroutine inituvp
!
subroutine cmpuvp(u,v,p,time)
  implicit none
!
  double precision, dimension(-1:nx+1,-1:ny+1) :: u,v,p,px, py
  double precision, dimension(-1:nx+1,-1:ny+1) :: ru,rv,rp, rpx, rpy, rqx, rqy
  double precision :: time
  !!double precision sol
  !!external sol
!
  double precision :: eu, ev, ep, dx, dy, x, y, mu, mv, mp, mpx, mpy, epx, epy
  integer :: i, j, k, isel
!
  dx = hg; dy = hg
!
  do j = 1, ny
    y = (dble(j)-0.5)*dy
    do i = 0, nx-1
      x = (dble(i)) *dx
      !!ru(i,j) = sol(x,y,time,0)
      call asol(ru(i,j),x,y,time,0)
    enddo
  enddo
  do j = 1, ny
    y = (dble(j))*dy
    do i = 1, nx
      x = (dble(i)-0.5)*dx
      !!rv(i,j) = sol(x,y,time,1)
      call asol(rv(i,j),x,y,time,1)
    enddo
  enddo
  do j = 1, ny
    y = dble(j-0.5)*dy
    do i = 1, nx
      x = dble(i-0.5)*dx
      !!rp(i,j) = sol(x,y,time,2)
      call asol(rp(i,j),x,y,time,2)
    enddo
  enddo
  eu = norm2(     ru(0:nx-1,1:ny)-u(0:nx-1,1:ny))
  mu = maxval(abs(ru(0:nx-1,1:ny)-u(0:nx-1,1:ny)))
  ev = norm2(     rv(1:nx,1:ny-1)-v(1:nx,1:ny-1))
  mv = maxval(abs(rv(1:nx,1:ny-1)-v(1:nx,1:ny-1)))
  ep = norm2(rp(1:nx,1:ny)-p(1:nx,1:ny))
  mp = maxval(abs(rp(1:nx,1:ny)-p(1:nx,1:ny)))

!!deb  rqx = 0.
!!deb  do j = 1, ny
!!deb    do i = 1, nx-1
!!deb      rqx(i,j) = p(i+1,j)-p(i,j)
!!deb    enddo
!!deb    i = 0
!!deb    rqx(i,j) = p(i+1,j)-p(nx,j)
!!deb  enddo
!!deb  rqx = rqx/dx
!!deb  do j = 1, ny
!!deb    y = (dble(j)-.5)*dy
!!deb    do i = 0, nx-1
!!deb      x = dble(i)*dx
!!deb      rpx(i,j) = sol(x,y,time,-2)
!!deb    enddo
!!deb  enddo
!!deb  rqy = 0.
!!deb  do j = 1, ny-1
!!deb    do i = 1, nx
!!deb      rqy(i,j) = p(i,j+1) - p(i,j)
!!deb    enddo
!!deb  enddo
!!deb  rqy = rqy/dy
!!deb  do j = 1, ny
!!deb    y = (dble(j))*dy
!!deb    do i = 1, nx
!!deb      x = (dble(i)-0.5)*dx
!!deb      rpy(i,j) = sol(x,y,time,-3)
!!deb    enddo
!!deb  enddo
!!deb  epx= norm2(rpx(0:nx-1,1:ny)-rqx(0:nx-1,1:ny))
!!deb  mpx= maxval(abs(rpx(0:nx-1,1:ny)-rqx(0:nx-1,1:ny)))
!!deb  epy= norm2(rpy(1:nx,1:ny-1)-rqy(1:nx,1:ny-1))
!!deb  mpy= maxval(abs(rpy(1:nx,1:ny-1)-rqy(1:nx,1:ny-1)))

  print '("Time@", e11.3, " e(L_2)u,v,p ",3(e14.6,1x), " e(L_inf) ", 3(e14.6,1x))', &
    & time, eu*hg, ev*hg, ep*hg, mu, mv, mp
  !!print '(1x,"Time @", e11.3, " error",1024(e14.6,1x))', time, eu*hg, ev*hg, ep*hg, &
  !!  & epx*hg, epx*hg, mu, mv, mp, mpx, mpy
  !!prt print *, '===================='
!!  print *, ''
!!!!  print *, 'Uc=['
!!  do i = 1,ny
!!    print '(1x,1024(e14.6,1x))', u(0:nx-1,i)
!!  enddo
!!!!  print *, ']'
!!!!  print *, 'Vc=['
!!  do i = 1,ny-1
!!    print '(1x,1024(e14.6,1x))', v(1:nx,i)
!!  enddo
!!  do i = 1,ny
!!    print '(1x,1024(e14.6,1x))', p(1:nx,i)
!!  enddo
!!!!  print *, ']'
!!!!  print *, 'U=['
!!  do i = 1,ny
!!    print '(1x,1024(e14.6,1x))', ru(0:nx-1,i)
!!  enddo
!!!!  print *, ']'
!!!!  print *, 'V=['
!!  do i = 1,ny-1
!!    print '(1x,1024(e14.6,1x))', rv(1:nx,i)
!!  enddo
!!  do i = 1,ny
!!    print '(1x,1024(e14.6,1x))', rp(1:nx,i)
!!  enddo
!!  print *, ']'

!
  return
end subroutine cmpuvp
!========================================================================
subroutine outuvpx(n,u,v,p,xpt,ypt,nct)
  integer :: n, nct
  double precision, dimension(-1:nx+1,-1:ny+1), intent(in) :: u,v,p
  double precision, dimension(n), intent(in) :: xpt, ypt
!
  integer :: i, strlen
  character(40) efile
!
  strlen = len_trim(runname)
  write(efile,'(2a,i4.4)') runname(1:strlen),'.ib.',nct
  open(67,file=efile,form='formatted',action='write')
  do i = 1, n
    write(67,'(2(e22.14,1x))')xpt(i),ypt(i)
  enddo
  write(67,'(2(e22.14,1x))')xpt(1),ypt(1)
  close(67)
  write(efile,'(2a,i4.4)') runname(1:strlen),'.u.',nct
  open(67,file=efile,access='stream',action='write')
  write(67)u
  close(67)
  write(efile,'(2a,i4.4)') runname(1:strlen),'.v.',nct
  open(67,file=efile,access='stream',action='write')
  write(67)v
  close(67)
  write(efile,'(2a,i4.4)') runname(1:strlen),'.p.',nct
  open(67,file=efile,access='stream',action='write')
  write(67)p
  close(67)
! 
  return
end subroutine outuvpx
!
!========================================================================

end module linsys
