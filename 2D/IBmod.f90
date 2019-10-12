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

module IBmod
  use, intrinsic :: iso_c_binding
  use parameters
  use IBforce
  use geometry, only :paraval
  use solver_mod, only : sik

  implicit none
!
private
  double precision, parameter :: h = hg
  double precision :: vc_nn((1)*nring*2)
  double precision :: vc_sm((1)*nring*2)
  double precision :: vc_cu((1)*nring*2)
  double precision :: vc_nt((1)*nring*2)
  double precision :: vc_us((1)*nring*2)
  double precision :: vc_xo((1)*nring*2)
  double precision :: vc_xs((1)*nring*2)
  double precision :: vc_xt((1)*nring*2)

  integer, parameter :: lensav=4*nring+15
  integer, parameter :: lenwrk=2*nring
  double precision :: vc_sv(lensav), fftwk(lenwrk)
  double precision :: vc_ft((1)*nring)
  integer :: ier, inc
!
!!public :: InitLocCell, xmovemac, spreadmac, IBspreadS, IBinterpS
!!public :: InitLocCell, xmovemac, IBspreadS, IBinterpS, fullspread
public :: InitLocCell, xmovemac, IBspreadS, IBinterpS

contains
!=======================================================================
subroutine InitLocCell(xpt, ypt, fs, cfs, gp)
  double precision, dimension(mpts) :: xpt, ypt
  double precision :: fs(mcoor), cfs(mcoor), gp(mpts)
!
  double precision :: r, rb, ra, ang, ds, xcmp, ycmp, xcma, ycma
  integer :: i, j, ient, jent
  double precision :: pa, pb, pc, theta
!
  xcma=0.35d0
  !!xcma=half
  ycma=half
  !!ycma=0.25d0
  !!ycma=0.125d0
!
  pi = cpi
  ang= two*pi/dble(nring)
  ds = dl/sin(half*ang)
  r  = half*ds
  rb = r*0.8
  ra = r/0.8
  !ra = r
  !rb = r
  pa = 6.0d0
  pb = 100.d0
  pc = 3.0d0
  pa = 4.0d0
  pb = 64.d0
  pc = 3.0d0

  ctild=-dl*dl*sin(ang)
  !!ctild=0. ! no preferential curvature
  cs(1) = ctild
!
  xcmp = xcma; ycmp = ycma
  do i = 1, nring
    theta = ang*dble(i-1) 
    !!xpt(i) = xcmp + ra*cos(dble(i-1)*ang)
    !!ypt(i) = ycmp + rb*sin(dble(i-1)*ang)
    xpt(i) = xcmp + cpi*(pa+cos(pc*theta))*cos(theta)/pb
    ypt(i) = ycmp + cpi*(pa+cos(pc*theta))*sin(theta)/pb
  enddo
  oxpt = xpt; oypt = ypt;
  xbp = xpt; ybp = ypt;
  fs = 0.d0; cfs = 0.d0
  call calcplateletforce(xpt,ypt,fs,cfs,gp,-1)
  !!fs = zero
  !!cfs = zero
!!prt  print *, 'the force?', maxval(fs), maxval(cfs), maxval(gp)
!!  call spreadmac(xpt,ypt,fs,ff)
!!  fs = cfs
!
  return
end subroutine InitLocCell
!========================================================================
double precision function ddel(r,isel)
  double precision :: r
  integer :: isel
!
!
  select case (isel)
  case (0) ! 0<r<1
    ddel = 0.125D0*(3.d0-2.d0*dabs(r)+dsqrt(1.d0+4.d0*dabs(r)-4.d0*r*r))
  case (1) ! 1<r<2
    ddel = 0.125D0*(5.d0-2.d0*dabs(r)-dsqrt(-7.d0+12.d0*dabs(r)-4.d0*r*r))
  case default
    print *, 'No such choice! stop'
    stop
  end select
!
!
  return
end function ddel
!========================================================================


!
!=======================================================================
!=======================================================================
!!!subroutine fullspread(xpt,ypt,mkp,mkq,mks,gp,fs,ff)
!!!  double precision :: gp(mpts)
!!!  double precision :: fs(mcoor), ff(-1:nx+1,-1:ny+1,2)
!!!  double precision :: xpt(mpts), ypt(mpts)
!!!  double precision :: mkp(7,nring,cent), mkq(7,nring,cent), mks(7,nring,cent)
!!!!
!!!!
!!!  return
!!!end subroutine
!=======================================================================
!========================================================================
subroutine IBspreadS(xpt,ypt,gp,fs,ff)
  implicit none
  double precision :: fs(mcoor)
  double precision :: xpt(mpts), ypt(mpts)
  double precision :: f(-1:nxp2,-1:nyp2,2), ff(-1:nx+1,-1:ny+1,2)
  double precision :: g(-1:nxp2,-1:nyp2,2)
  double precision :: gp(mpts)
!
! for mac grid 
!
! The spreading is done to an extended domain with 2 extra cells above
! and below the tube walls.  This force will be mapped (with reflection)
! to a smaller force array before being passed to fluidstep.
!
! spreads to 4x4 portion of grid directly
!
! lower left of the spreading grid  is (xamin,yamin)
! upper right of the spreading grid is (xamax,yamax)
!
! ib point should be located between
!    xamin + 2h and xamax-2h
!    yamin + 2h and yamax-2h
! 
  integer :: i,j,k, n, il, nlp, nil, njl, info
  integer :: jmarrLR,kmarrLR
  integer :: jmarrTB,kmarrTB
  integer :: iflg
  integer :: jmp,kmp
  integer :: ix, jx
!
  double precision :: xxsft, yysft, rh, ds, f1, f2, xper100, yper100
  double precision :: x, y, x0, y0, x1, y1, x2, y2
  double precision :: dxl, dyl, rxl, ryl, dxt, dyt, rxt, ryt
  double precision :: uwtLR(4), vwtLR(4), uwtTB(4), vwtTB(4)
  double precision :: s0, tmp, tp1, tp2
!
  xamin = xmin
  yamin = ymin
  xper100 = (xmax-xmin)*100.0d0
  yper100 = (ymax-ymin)*100.0d0
!
  rh   = 1.0d0/hg
!========================================================================
  f(0:nx,0:ny,1:2)=zero
  g(0:nx,0:ny,1:2)=zero
!
  ds=cpi*two/dble(nring)
  do n=1,mpts
        x=xpt(n)+xper100
        y=ypt(n)+yper100
        x0=xpt(n); y0=ypt(n); 
    il  = (n-2*nfil-1)/nring + 1
    nlp = (n-2*nfil) - (n-2*nfil-1)/nring*nring
    nil = nlp + (il-1)*2*nring
    njl = nil + nring
    s0 = mkq(1,nlp,il)
    tmp= paraval(s0,nring,mkq,il,1,0,info)
    !!tp1=(tmp+sik)*ndv(nring-nlp+1,1,il)
    !!tp2=(tmp+sik)*ndv(nring-nlp+1,2,il)
    tp1=(tmp+sik)*mtd(nring-nlp+1,1,il)
    tp2=(tmp+sik)*mtd(nring-nlp+1,2,il)

!        y=ypt(n)

        !!oldForce f1= dl*fs(2*n-1)
        !!oldForce f2= dl*fs(2*n  )
!	tp1 = zero !!DEBUG
!	tp2 = zero !!DEBUG
        f1= ds*(fs(2*n-1)+tp1)
        f2= ds*(fs(2*n  )+tp2)
        !!f1= ds*(fs(2*n-1)+tp1/gp(n))
        !!f2= ds*(fs(2*n  )+tp2/gp(n))

! LR edges
        xxsft = 0.0d0
        yysft = 0.5d0
        jmarrLR=int((x-xamin)*rh+xxsft)
        kmarrLR=int((y-yamin)*rh+yysft)
!
! compute weights
!
        ix = int((x0-xamin)*rh+xxsft); jx = int((y0-yamin)*rh+yysft)
        x1 = xamin+h*(dble(ix)-xxsft); y1 = yamin+h*(dble(jx)-yysft);
        dxl= x1-x0; dyl = y1-y0; rxl= dxl*rh; ryl = dyl*rh
!
        uwtLR(1)=ddel(rxl-1.d0,1); vwtLR(1)=ddel(ryl-1.d0,1);
        uwtLR(2)=ddel(rxl     ,0); vwtLR(2)=ddel(ryl     ,0);
        uwtLR(3)=ddel(rxl+1.d0,0); vwtLR(3)=ddel(ryl+1.d0,0);
        uwtLR(4)=ddel(rxl+2.d0,1); vwtLR(4)=ddel(ryl+2.d0,1);
        uwtLR = uwtLR*rh;          vwtLR = vwtLR*rh;

! for TB edges
        xxsft = 0.5
        yysft = 0.0
        jmarrTB=int((x-xamin)*rh+xxsft)
        kmarrTB=int((y-yamin)*rh+yysft)
!
        ix = int((x0-xamin)*rh+xxsft); jx = int((y0-yamin)*rh+yysft)
        x2 = xamin+h*(dble(ix)-xxsft); y2 = yamin+h*(dble(jx)-yysft);
        dxt= x2-x0; dyt = y2-y0; rxt= dxt*rh; ryt = dyt*rh
!
        uwtTB(1)=ddel(rxt-1.d0,1); vwtTB(1)=ddel(ryt-1.d0,1);
        uwtTB(2)=ddel(rxt     ,0); vwtTB(2)=ddel(ryt     ,0);
        uwtTB(3)=ddel(rxt+1.d0,0); vwtTB(3)=ddel(ryt+1.d0,0);
        uwtTB(4)=ddel(rxt+2.d0,1); vwtTB(4)=ddel(ryt+2.d0,1);
        uwtTB = uwtTB*rh;          vwtTB = vwtTB*rh;
!
!!deb      print '(1x,"which part ", 2(i5,1x),1024(e14.6,1x))', jmarrTB, kmarrTB, &
!!deb        x,y,xpt(n),ypt(n), xper100, yper100, xxsft,yysft, uwtLR(2), vwtLR(2), f1, f2
!  spread these forces
!
!  forces spread to grid points with indices jmarr-1,..,jmarr+2 and kmarr-1,..,kmarr+2.
!
    do k=1,4
    do j=1,4
      jmp=mod((jmarrLR-2+j)+100*nx,nx)
      kmp=mod((kmarrLR-2+k)+100*ny,ny)
      g(jmp,kmp,1) = g(jmp,kmp,1) + uwtLR(j)*vwtLR(k)*f1
      jmp=mod((jmarrTB-2+j)+100*nx,nx)
      kmp=mod((kmarrTB-2+k)+100*ny,ny)
      g(jmp,kmp,2) = g(jmp,kmp,2) + uwtTB(j)*vwtTB(k)*f2
    enddo
    enddo
  enddo        ! end do n=1,npts
!
!!  print '(1x,"Difference in force: ", 6(e16.8,1x))', &
!!   & norm2(f(0:nx,0:ny,1)-g(0:nx,0:ny,1)), norm2(f(0:nx,0:ny,2)-g(0:nx,0:ny,2)), &
!!   & norm2(f(0:nx,0:ny,1)), norm2(g(0:nx,0:ny,1)),  &
!!   & norm2(f(0:nx,0:ny,2)), norm2(g(0:nx,0:ny,2)) 
! extend periodically
!
  do j=1,ny
    g(nx  ,j,1)=g(0,j,1)
    g(nx+1,j,1)=g(1,j,1)
    g(nx  ,j,2)=g(0,j,2)
    g(nx+1,j,2)=g(1,j,2)
  enddo
  g(nx  ,ny+1,2)=g(0,ny+1,2)
  g(nx+1,ny+1,2)=g(1,ny+1,2)
!
  ff(-1:nx+1,-1:ny+1,1:2) = g(-1:nx+1,-1:ny+1,1:2)
!========================================================================
! Delta function is new
!========================================================================
!
  return
end subroutine IBspreadS
!========================================================================
!========================================================================
subroutine IBinterpS(xpt,ypt,u0,v0,us)
!  interpolate Eulerian variable (u0,v0) to the variable us, which is
!  defined @ IB locations (xpt,ypt). 
  implicit none
  double precision :: xpt(mpts), ypt(mpts)
  double precision, dimension(-1:nx+1,-1:ny+1) :: u0, v0
  double precision :: us(mpts,2) !<-this is the output
!
  integer :: jmarrLR,kmarrLR
  integer :: jmarrTB,kmarrTB
  integer :: n,ngrid, nlp, ncp, nil, njl, il, info, iflag
  integer :: j,k, iflg
  integer :: jmp,kmp
  double precision :: xper100,yper100
  double precision :: uwtLR(4),vwtLR(4)
  double precision :: uwtTB(4),vwtTB(4)
  double precision :: x,y, x0,y0, x2,y2, x1, y1, xp, yp, xxsft, yysft
  double precision :: h, rh, dxl, dyl, rxl, ryl, dxt, dyt, rxt, ryt
  integer :: ix, jx
!
  xamin = xmin; yamin = ymin
  xper100 = (xmax-xmin)*100.0d0
  yper100 = (ymax-ymin)*100.0d0
  rh   = one/h
!
  do n=2*nfil+1,npts
    x0 = xpt(n); y0 = ypt(n)
    x  = x0 + xper100; y = y0 + yper100
    !
! for LR edges
    xxsft = 0.0d0
    yysft = 0.5d0
    jmarrLR=int((x-xamin)*rh+xxsft)
    kmarrLR=int((y-yamin)*rh+yysft)
!
    ix = int((x0-xamin)*rh+xxsft); jx = int((y0-yamin)*rh+yysft)
    x1 = xamin+h*(dble(ix)-xxsft); y1 = yamin+h*(dble(jx)-yysft);
    dxl= x1-x0; dyl = y1-y0; rxl= dxl*rh; ryl = dyl*rh
!
    uwtLR(1)=ddel(rxl-1.d0,1); vwtLR(1)=ddel(ryl-1.d0,1);
    uwtLR(2)=ddel(rxl     ,0); vwtLR(2)=ddel(ryl     ,0);
    uwtLR(3)=ddel(rxl+1.d0,0); vwtLR(3)=ddel(ryl+1.d0,0);
    uwtLR(4)=ddel(rxl+2.d0,1); vwtLR(4)=ddel(ryl+2.d0,1);
    uwtLR = uwtLR*rh;          vwtLR = vwtLR*rh;

! for TB edges
    xxsft = 0.5d0
    yysft = 0.0d0
    jmarrTB=int((x-xamin)*rh+xxsft)
    kmarrTB=int((y-yamin)*rh+yysft)
!
    ix = int((x0-xamin)*rh+xxsft); jx = int((y0-yamin)*rh+yysft)
    x2 = xamin+h*(dble(ix)-xxsft); y2 = yamin+h*(dble(jx)-yysft);
    dxt= x2-x0; dyt = y2-y0; rxt= dxt*rh; ryt = dyt*rh
!
    uwtTB(1)=ddel(rxt-1.d0,1); vwtTB(1)=ddel(ryt-1.d0,1);
    uwtTB(2)=ddel(rxt     ,0); vwtTB(2)=ddel(ryt     ,0);
    uwtTB(3)=ddel(rxt+1.d0,0); vwtTB(3)=ddel(ryt+1.d0,0);
    uwtTB(4)=ddel(rxt+2.d0,1); vwtTB(4)=ddel(ryt+2.d0,1);
    uwtTB = uwtTB*rh;          vwtTB = vwtTB*rh;
!
!compute weights
!
!  interpolate these Lagrangian quantity by computing
!  the weighted sum that gives the values at (xpt,ypt).
!
    us(n,1)=0.d0
    us(n,2)=0.d0
!
    do k=1,4
    do j=1,4
      jmp=mod((jmarrLR-2+j)+100*nx,nx)
      kmp=mod((kmarrLR-2+k)+100*ny,ny)
      !!CosDelt us(n,1)=us(n,1)+xwtLR(j)*ywtLR(k)*u0(jmp,kmp) 
      us(n,1)=us(n,1)+uwtLR(j)*vwtLR(k)*u0(jmp,kmp) 
      jmp=mod((jmarrTB-2+j)+100*nx,nx)
      kmp=mod((kmarrTB-2+k)+100*ny,ny)
      !!CosDelt us(n,2)=us(n,2)+xwtTB(j)*ywtTB(k)*v0(jmp,kmp) 
      us(n,2)=us(n,2)+uwtTB(j)*vwtTB(k)*v0(jmp,kmp) 
    enddo
    enddo
    us(n,:)=us(n,:)*h*h
    
  enddo

!
  return
end subroutine IBinterpS

!========================================================================
!========================================================================
subroutine Jacob(n, xp0, xpw, f0, fdir,iflag)
  !! here f0 = F(xp0) is input
  implicit none
  integer :: n, iflag
  double precision :: xp0(n), xpw(n), fdir(n), f0(n)
  double precision :: xt(mpts), yt(mpts), xw(mpts), yw(mpts), gp(mpts)
  double precision :: vc_cf(n), vc_f(n), del(n)
!
  double precision :: tmp, epsnew, normw, xs
!

  fdir = zero
  epsnew = 1.d-8
  normw = sqrt(dot_product(xpw,xpw))

  if (normw < 1.d-14) then
    fdir = zero
    print '(1x,"direction should not be 0!!", e13.6)',  normw
    goto 123
  endif
  epsnew = epsnew/normw;

  xs = sqrt(dot_product(xp0,xp0))
  if (abs(xs) > 1.d-12) then
    epsnew=epsnew*xs
  endif
  del = xp0 + epsnew*xpw
  
  call vc_xr(del,xt,yt)
!
  call calcplateletforce(xt,yt,vc_f,vc_cf,gp,1)
  !!call vc_fn(n, del, vc_cf,iflag)
  fdir = (vc_cf - f0)/epsnew
! 
  123 continue
!
  return
end subroutine Jacob
!========================================================================
!

!!old subroutine xmovemac(xpt,ypt, us,u,fs,cfs)
subroutine xmovemac(xpt,ypt,us,u0,v0,fs,cfs,dt)
  implicit none
  double precision :: dt
  double precision :: xpt(mpts), ypt(mpts), gp(mpts)
  double precision, dimension(-1:nx+1,-1:ny+1) :: u0, v0
!
! for mac grid
!
! The interpolation is done from an extended domain with 2 extra cells above
! and below the tube walls.  This velocity was mapped (with reflection)
! from a smaller velocity array that was passed from fluidstep.
!
!   interpolates from 4x4 portion of grid directly
!
! lower left  of the spreading grid is (xamin,yamin)
! upper right of the spreading grid is (xamax,yamax)
!
! ib point should be located between
!    xamin + 2h and xamax-2h
!    yamin + 2h and yamax-2h
!
  integer :: jmarrLR,kmarrLR
  integer :: jmarrTB,kmarrTB
  integer :: n,ngrid, nlp, ncp, nil, njl, il, info, iflag
  integer :: j,k, iflg
  integer :: jmp,kmp
!
  double precision :: tol, res
  double precision :: fs(mcoor)
  double precision :: cfs(mcoor),vc_f(mcoor), vc_cf(mcoor)
  double precision :: vc_xpt(mpts), vc_ypt(mpts)
  double precision :: ccxt, ccyt, ccx1, ccy1, ccx2, ccy2, tp1, tp2, tmp
  double precision :: vc_xx((1)*nring*2)
  double precision :: vc_x0((1)*nring*2)
  double precision ::  vc_b((1)*nring*2)
  double precision :: vc_bh((1)*nring*2)


  double precision :: h2,x,y
  double precision :: wt1x,wt1y,wt2x,wt2y
  double precision :: rh
  double precision :: argxLR,argyLR
  double precision :: argxTB,argyTB
  double precision :: xxsft,yysft
  double precision :: cy1,cy2,sy1,sy2
  double precision :: cx1,cx2,sx1,sx2
  double precision :: cx,cy,sx,sy,xper100,yper100
  double precision :: us(mpts,2)
  double precision :: xwtLR(4),ywtLR(4)
  double precision :: xwtTB(4),ywtTB(4)
  double precision :: u(-1:nxp2,-1:nyp2,2)
!
!========================================================================
  double precision :: uwtLR(4),vwtLR(4)
  double precision :: uwtTB(4),vwtTB(4)
  double precision :: x0,y0, x2,y2, x1, y1, xp, yp
  double precision :: dxl, dyl, rxl, ryl, dxt, dyt, rxt, ryt
  integer :: ix, jx
!========================================================================
!
!!prt    print '("in sub ", 1024(e14.6,1x))', maxval(u0), minval(u0),&
!!prt    maxval(v0), minval(v0), u0(88,63), v0(88,63)
! once only ccccccccccccccccccccccccccccccccccccccccc
!
  xamin = xmin; yamin = ymin
  xper100 = (xmax-xmin)*100.0d0
  yper100 = (ymax-ymin)*100.0d0
  rh   = one/hg
  wt1x = two*atan(1.0d0)*rh
  wt1y = two*atan(1.0d0)*rh
  wt2x = .25d0*rh
  wt2y = .25d0*rh
  h2   = h*h
  cy1  = cos(wt1y*h)
  sy1  = sin(wt1y*h)
  cy2  = cos(2.0d0*wt1y*h)
  sy2  = sin(2.0d0*wt1y*h)
  cx1  = cos(wt1x*h)
  sx1  = sin(wt1x*h)
  cx2  = cos(2.0d0*wt1x*h)
  sx2  = sin(2.0d0*wt1x*h)
!
!!!!!      do n=1,npts
!      write(6,*)' xmovemac'
!      write(6,*)' 2*nfil+1 = ',2*nfil+1
!      write(6,*)' npts = ',npts
      vc_xpt = xpt; vc_ypt = ypt
      call calcplateletforce(vc_xpt,vc_ypt,vc_f,vc_cf,gp,1)
  do n=2*nfil+1,npts
    il  = (n-2*nfil-1)/nring + 1
    nlp = (n-2*nfil) - (n-2*nfil-1)/nring*nring
    ncp = n - 2*nfil
    nil = nlp + (il-1)*2*nring
    njl = nil + nring
    !!print '(1x,1024(i4,1x))', n, il, nlp, ncp, nil, njl
!
! add large multiple of period to x,y so x>0 and y>0
! to avoid testing
!
    x0 = xpt(n); y0 = ypt(n)
    x  = x0 + xper100; y = y0 + yper100
!!    x=xpt(n)+xper100
!!!!!            y=ypt(n)
!!    y=ypt(n)+yper100
!
! for LR edges
    xxsft = 0.0d0
    yysft = 0.5d0
    jmarrLR=int((x-xamin)*rh+xxsft)
    kmarrLR=int((y-yamin)*rh+yysft)
    argxLR=wt1x*((x-xamin)-dble(jmarrLR)*h+xxsft*h)
    argyLR=wt1y*((y-yamin)-dble(kmarrLR)*h+yysft*h)
!
!compute weights
!
    !!cy = cos(argyLR)
    !!sy = sin(argyLR)
    !!ywtLR(1)=wt2y*(1.0 + cy*cy1 - sy*sy1)
    !!ywtLR(2)=wt2y*(1.0 + cy)
    !!ywtLR(3)=wt2y*(1.0 + cy*cy1 + sy*sy1)
    !!ywtLR(4)=wt2y*(1.0 + cy*cy2 + sy*sy2)
    !!cx = cos(argxLR)
    !!sx = sin(argxLR)
    !!xwtLR(1)=wt2x*(1.0 + cx*cx1 - sx*sx1)
    !!xwtLR(2)=wt2x*(1.0 + cx)
    !!xwtLR(3)=wt2x*(1.0 + cx*cx1 + sx*sx1)
    !!xwtLR(4)=wt2x*(1.0 + cx*cx2 + sx*sx2)
!
    ix = int((x0-xamin)*rh+xxsft); jx = int((y0-yamin)*rh+yysft)
    x1 = xamin+h*(dble(ix)-xxsft); y1 = yamin+h*(dble(jx)-yysft);
    dxl= x1-x0; dyl = y1-y0; rxl= dxl*rh; ryl = dyl*rh
!
    uwtLR(1)=ddel(rxl-1.d0,1); vwtLR(1)=ddel(ryl-1.d0,1);
    uwtLR(2)=ddel(rxl     ,0); vwtLR(2)=ddel(ryl     ,0);
    uwtLR(3)=ddel(rxl+1.d0,0); vwtLR(3)=ddel(ryl+1.d0,0);
    uwtLR(4)=ddel(rxl+2.d0,1); vwtLR(4)=ddel(ryl+2.d0,1);
    uwtLR = uwtLR*rh;          vwtLR = vwtLR*rh;

! for TB edges
    xxsft = 0.5d0
    yysft = 0.0d0
    jmarrTB=int((x-xamin)*rh+xxsft)
    kmarrTB=int((y-yamin)*rh+yysft)
    argxTB=wt1x*((x-xamin)-dble(jmarrTB)*h+xxsft*h)
    argyTB=wt1y*((y-yamin)-dble(kmarrTB)*h+yysft*h)
!
!compute weights
!
    !!cy = cos(argyTB)
    !!sy = sin(argyTB)
    !!ywtTB(1)=wt2y*(1.0 + cy*cy1 - sy*sy1)
    !!ywtTB(2)=wt2y*(1.0 + cy)
    !!ywtTB(3)=wt2y*(1.0 + cy*cy1 + sy*sy1)
    !!ywtTB(4)=wt2y*(1.0 + cy*cy2 + sy*sy2)
    !!cx = cos(argxTB)
    !!sx = sin(argxTB)
    !!xwtTB(1)=wt2x*(1.0 + cx*cx1 - sx*sx1)
    !!xwtTB(2)=wt2x*(1.0 + cx)
    !!xwtTB(3)=wt2x*(1.0 + cx*cx1 + sx*sx1)
    !!xwtTB(4)=wt2x*(1.0 + cx*cx2 + sx*sx2)
!
    ix = int((x0-xamin)*rh+xxsft); jx = int((y0-yamin)*rh+yysft)
    x2 = xamin+h*(dble(ix)-xxsft); y2 = yamin+h*(dble(jx)-yysft);
    dxt= x2-x0; dyt = y2-y0; rxt= dxt*rh; ryt = dyt*rh
!
    uwtTB(1)=ddel(rxt-1.d0,1); vwtTB(1)=ddel(ryt-1.d0,1);
    uwtTB(2)=ddel(rxt     ,0); vwtTB(2)=ddel(ryt     ,0);
    uwtTB(3)=ddel(rxt+1.d0,0); vwtTB(3)=ddel(ryt+1.d0,0);
    uwtTB(4)=ddel(rxt+2.d0,1); vwtTB(4)=ddel(ryt+2.d0,1);
    uwtTB = uwtTB*rh;          vwtTB = vwtTB*rh;
!
!  interpolate these velocities by computing
!  the weighted sum that gives the velocity at (xpt,ypt).
!
    us(n,1)=0.d0
    us(n,2)=0.d0
!
    do k=1,4
    do j=1,4
      jmp=mod((jmarrLR-2+j)+100*nx,nx)
      kmp=mod((kmarrLR-2+k)+100*ny,ny)
      !!old us(n,1)=us(n,1)+xwtLR(j)*ywtLR(k)*u(jmp,kmp,1) 
      !!CosDelt us(n,1)=us(n,1)+xwtLR(j)*ywtLR(k)*u0(jmp,kmp) 
      us(n,1)=us(n,1)+uwtLR(j)*vwtLR(k)*u0(jmp,kmp) 
      jmp=mod((jmarrTB-2+j)+100*nx,nx)
      kmp=mod((kmarrTB-2+k)+100*ny,ny)
      !!old us(n,2)=us(n,2)+xwtTB(j)*ywtTB(k)*u(jmp,kmp,2) 
      !!CosDelt us(n,2)=us(n,2)+xwtTB(j)*ywtTB(k)*v0(jmp,kmp) 
      us(n,2)=us(n,2)+uwtTB(j)*vwtTB(k)*v0(jmp,kmp) 
!!deb      print '("useinfo ", 7(i3,1x), 1024(e14.6,1x))', n, jmp, kmp, &
!!deb       jmp,kmp, j,k, uwtLR(j),vwtLR(k),u0(jmp,kmp),uwtTB(j),vwtTB(k),v0(jmp,kmp)
    enddo
    enddo
!=======================================================================
! set geometric data, mkn(1:nring) is for xpt(2*nfil+1:2*nfil+nring):backwards
!=======================================================================
!!debif (n .eq. 2*nfil+1) then
!!deb  !!il = 1
!!deb  ccxt = paraval(mkn(1,1,il), nring, mkn, il, 0, 0, info)
!!deb  ccyt = paraval(mkn(1,1,il), nring, mkn, il, 1, 0, info)
!!deb  ccx1 = paraval(mkn(1,2,il), nring, mkn, il, 0, 0, info)
!!deb  ccy1 = paraval(mkn(1,2,il), nring, mkn, il, 1, 0, info)
!!deb  ccx2 = paraval(mkn(1,nring,il), nring, mkn, il, 0, 0, info)
!!deb  ccy2 = paraval(mkn(1,nring,il), nring, mkn, il, 1, 0, info)
!!deb  print '(1x,"P:", 1024(e15.8,1x))', ccxt, ccyt, &
!!deb    xpt(npts), ypt(npts), ccx1, ccy1, xpt(npts-1), ypt(npts-1), ccx2, ccy2, xpt(n), ypt(n)
!!debendif
    oxpt(n)=xpt(n)
    oypt(n)=ypt(n)

    vc_nn(nil)=ndv(nring-nlp+1,1,il)
    vc_nn(njl)=ndv(nring-nlp+1,2,il)

    vc_cu(nil) = paraval(mkn(1,nring-nlp+1,il), nring, mkp, il, 0, 0, info)
    vc_cu(njl) = paraval(mkn(1,nring-nlp+1,il), nring, mkp, il, 1, 0, info)
    !!vc_cu(nil) = 0.d0
    !!vc_cu(njl) = 0.d0
    tmp = paraval(mkn(1,nring-nlp+1,il), nring, mkq, il, 1, 0, info)
    tp1 = paraval(mkn(1,nring-nlp+1,il), nring, mkr, il, 1, 0, info)
!!    print *, 'mk', tmp, tp1
    tmp = tp1 + tmp
    !!vc_nt(nil) = nt_sigma(tmp) ! tmp is theta
    !!vc_nt(njl) = vc_nt(nil)
    vc_nt(nil) = tmp ! tmp is sigma_n
    vc_nt(njl) = vc_nt(nil)
    !!if (nil .eq. 1) print *, vc_cu(nil), vc_cu(njl)
    !!prt print *, vc_cu(nil), vc_cu(njl)
    vc_us(nil) = us(n,1)*h2
    vc_us(njl) = us(n,2)*h2
    vc_xo(nil)=xpt(n);     
    vc_xo(njl)=ypt(n)
    vc_xx(nil)=xpt(n);     
    vc_xx(njl)=ypt(n)
!
    !!tmp = vc_cu(njl)-vc_cu(nil) + (vc_cf(2*n-1)*vc_nn(nil)+vc_cf(2*n)*vc_nn(njl))
    !!!!tmp = vc_cu(njl)-vc_cu(nil) + (cfs(2*n-1)*vc_nn(nil)+cfs(2*n)*vc_nn(njl))
    tp1 = vc_cu(njl)-vc_cu(nil)
    vc_cu(nil) = (tp1+tmp*rtr)*vc_nn(nil)*kw(il) 
    vc_cu(njl) = (tp1+tmp*rtr)*vc_nn(njl)*kw(il) 
!!old    vc_cu(nil) = (tp1+tmp/gp(n))*vc_nn(nil)*kw(il) 
!!old    vc_cu(njl) = (tp1+tmp/gp(n))*vc_nn(njl)*kw(il) 
    !!debvc_cu(nil) = (tp1          )*vc_nn(nil)*kw(il) 
    !!debvc_cu(njl) = (tp1          )*vc_nn(njl)*kw(il) 
    !!vc_cu(nil) = tp1*vc_nn(nil)*kw(il) 
    !!vc_cu(njl) = tp1*vc_nn(njl)*kw(il)
    !!vc_cu(nil) = tp1*vc_nn(nil)*kw(il) + (vc_nt(nil)+rho0/(eta+etas)*vc_nn(nil))*vc_nn(nil)*gp(n)*kw(il)
    !!vc_cu(njl) = tp1*vc_nn(njl)*kw(il) + vc_nt(njl)*vc_nn(njl)*gp(n)*kw(il)
    !!vc_cu(nil) = tp1*vc_nn(nil)*kw(il)*rtc
    !!vc_cu(njl) = tp1*vc_nn(njl)*kw(il)*rtc
    !! print *, 'kw=', il, nil, njl, kw(il)

    !!explicitupdate xpt(n)=xpt(n)+dt*(us(n,1)*h2+kw(il)*tmp*ndv(nring-nlp+1,1,il) )
    !!explicitupdate ypt(n)=ypt(n)+dt*(us(n,2)*h2+kw(il)*tmp*ndv(nring-nlp+1,2,il) )

    !!original_explicit xpt(n)=xpt(n)+dt*us(n,1)*h2 ! original IB loc updates:  x
    !!original_explicit ypt(n)=ypt(n)+dt*us(n,2)*h2 ! original IB loc updates:  y
!=======================================================================
! new method
!=======================================================================
!=======================================================================
  enddo    ! end do n=2*nfil+1,npts
!!prt      write(*,*)' max us = ',maxval(abs(us)*h2)
  !!call calcplateletforce(vc_xpt,vc_ypt,vc_f,vc_cf)
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
! solving x^{n+1} with Newton's method
   !!vc_x0 = 0.
  if (dot_product(vc_xo,vc_xo) < 1d-14) then
    print *, 'zero vector!'
    stop
  endif
  vc_xs = vc_xo
!
  !!vc_sm = 0.
  vc_sm = vc_cu+vc_us
        
  call vc_rh(nring*2,vc_xo,vc_bh,iflag,dt)
  !!call vc_rb(nring*2,vc_bh,iflag,dt)
  !!vc_x0 = 0.
  vc_x0 = vc_xo
  res=1d-7
!!deb  do n = 2*nfil+1, npts
!!deb    il  = (n-2*nfil-1)/nring + 1
!!deb    nlp = (n-2*nfil) - (n-2*nfil-1)/nring*nring 
!!deb    nil = nlp + (il-1)*2*nring
!!deb    njl = nil + nring
!!deb    print '(1x,"checking ", 1024(e14.6,1x))', vc_xo(nil), vc_xo(njl), &
!!deb     vc_bh(nil), vc_bh(njl), vc_cu(nil), vc_cu(njl), vc_us(nil), &
!!deb     vc_us(njl), us(n,1), us(n,2)
!!deb  enddo
  
  !!call vc_gmrs(nring*2, 10, vc_x0, vc_b, vc_xo, res)
  call vc_gmrs(nring*2, 10, vc_x0, vc_bh, vc_xo, res, dt)
!!Newton !!  print *, 'partition 2'
  vc_xt = vc_x0
  do n = 2*nfil+1, npts
    il  = (n-2*nfil-1)/nring + 1
    nlp = (n-2*nfil) - (n-2*nfil-1)/nring*nring 
    nil = nlp + (il-1)*2*nring
    njl = nil + nring
!!    print '(1x,10(e12.6,1x))', vc_xx(nil), vc_xx(njl), xpt(n), ypt(n),  &
!!              vc_b(nil), vc_b(njl)
!!    xpt(n) = vc_xx(nil)
!!    ypt(n) = vc_xx(njl)
    xpt(n) = vc_xt(nil)
    ypt(n) = vc_xt(njl)
  enddo
  !!call vc_rx(xpt,ypt,vc_xx)
  !!call vc_cfs(vc_xx, vc_fe)

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!========================================================================
!
      return
end subroutine xmovemac
!
!=======================================================================
!

!!rmdup double precision function paraval(s, n, mk, il, idir, icase, info)
!!rmdup   integer :: lda, n
!!rmdup   double precision :: s, mk(7,n, 1)
!!rmdup   integer :: il, idir, icase, info
!!rmdup !
!!rmdup   double precision :: ck, eps
!!rmdup   integer :: i, k, inc, ijump
!!rmdup ! 
!!rmdup   !tupi = 2.*3.14159265358979323846d0
!!rmdup   tupi = two*cpi
!!rmdup   select case (idir)
!!rmdup   case (0)
!!rmdup     inc=2
!!rmdup   case (1)
!!rmdup     inc=5
!!rmdup   case default
!!rmdup     print *, 'No such option, check carefully'
!!rmdup     stop
!!rmdup   end select
!!rmdup !
!!rmdup   eps = 1e-10
!!rmdup   info = -1
!!rmdup   ! locate k
!!rmdup   ck = mod(mod(s, tupi)+tupi, tupi) ! archlength we will used
!!rmdup !  ck = s0
!!rmdup !PRT  print *, 's0', s0, 'ck', ck
!!rmdup   k = 1
!!rmdup   ijump = 0
!!rmdup   if ( ck > eps ) then
!!rmdup     if (ck < mk(1,k,il)) then
!!rmdup       ck = ck + tupi
!!rmdup       ijump = 1
!!rmdup       goto 654
!!rmdup     endif
!!rmdup     do while (k < n .and. ijump .eq. 0)
!!rmdup       !!if ( ck > mk(1,k,il) ) then 
!!rmdup       !!if ( ck > mk(1,k,il) ) then 
!!rmdup       !!  k = k + 1
!!rmdup       !!else  
!!rmdup       !!  ijump = 1
!!rmdup       !!endif  
!!rmdup       if ( ck .ge. mk(1,k,il) .and. ck < mk(1,k+1,il)) then
!!rmdup         ijump = 1
!!rmdup       else
!!rmdup         k = k + 1
!!rmdup       endif
!!rmdup     enddo
!!rmdup   else
!!rmdup     !!ck = ck + tupi
!!rmdup     k = 1
!!rmdup     ijump = 1
!!rmdup   endif
!!rmdup   if (k .eq. n) then 
!!rmdup     if (ijump .eq. 0) then
!!rmdup       if ( ck .ge. mk(1,k,il) .and. ck < tupi) then
!!rmdup         ijump = 1
!!rmdup       else
!!rmdup         print *, "really?", ck, s, mk(1,1,il)
!!rmdup         stop
!!rmdup       endif
!!rmdup     else
!!rmdup       print *, 'this could not happen!'
!!rmdup       stop
!!rmdup     endif
!!rmdup   endif
!!rmdup   if (k .eq. 1 .and. ijump .eq. 1) then
!!rmdup     !!print *, 'small ck', ck
!!rmdup     ck = tupi + ck
!!rmdup     !!ck = mod(ck, tupi)+tupi ! archlength we will used
!!rmdup   endif
!!rmdup   654 continue
!!rmdup   paraval = mk(inc,k,il) + mk(inc+1,k,il)*ck + mk(inc+2,k,il)*ck*ck
!!rmdup   if (icase .eq. 1) then
!!rmdup     !!paraval = mk(inc,k,il) + mk(inc+1,k,il)*ck + mk(inc+2,k,il)*ck*ck
!!rmdup     paraval =  mk(inc+1,k,il) + 2.*mk(inc+2,k,il)*ck
!!rmdup   else if (icase .eq. -1) then ! linear interpolation
!!rmdup   endif
!!rmdup   !!prt print *, 'ck = ', ck, ' k = ', k, 'inc = ', inc, 'idir = ', idir, 's= ', s
!!rmdup   if (ijump .eq. 1) then
!!rmdup     info = 0
!!rmdup   else
!!rmdup     print '(1x,"why here?", 4(i3,1x), 1024(e16.8,1x))',k, info, n, ijump, s, ck, paraval
!!rmdup     stop
!!rmdup   endif
!!rmdup !
!!rmdup   return
!!rmdup end function paraval
!
!
subroutine vc_matv(n,vc_xx,vc_xn,vc_rs,iflag,dt)
implicit none
! vc_xn gives the old time value x^n
  integer :: n, iflag
  double precision :: dt
  double precision :: vc_xx(n), vc_xn(n), vc_rs(n)
  double precision :: vc_f0(n), vc_r(n)
!
  integer :: ncp, nn, nlp, il, nil, njl
  double precision :: tmp
!
  call vc_fn(n,vc_xn,vc_f0,iflag)

  call vc_dirder(n,vc_xn,vc_xx,vc_r,vc_f0,iflag)
!!pr  print '(1x, "mtv dir", 10(e22.14,1x))', sqrt(dot_product(vc_r,vc_r))

  vc_rs = vc_xx - dt*vc_r
!
  return
end subroutine vc_matv
!
!
subroutine vc_rh(n,vc_xn,vc_rs,iflag,dt)
! vc_xn gives the old time value x^n
  integer :: n, iflag
  double precision :: dt
  double precision :: vc_xn(n), vc_rs(n)
  double precision :: vc_f0(n), vc_r(n)

  integer :: il, nn, nlp, nil, njl
!
  call vc_fn(n,vc_xn,vc_f0,iflag)
!
  call vc_dirder(n,vc_xn,vc_xn,vc_r,vc_f0,iflag)
!!pr  print '(1x, "RHS dir", 10(e22.14,1x))', sqrt(dot_product(vc_r,vc_r))
!
  vc_rs = vc_xn + dt*(vc_f0-vc_r)

!!deb  print *, 'what is dt', dt
!!deb  do nn = 2*nfil+1, npts
!!deb    il  = (nn-2*nfil-1)/nring + 1
!!deb    nlp = (nn-2*nfil) - (nn-2*nfil-1)/nring*nring 
!!deb    nil = nlp + (il-1)*2*nring
!!deb    njl = nil + nring
!!deb    print '(1x,"wrong parts ", 1024(e14.6,1x))', vc_xn(nil), vc_xn(njl), &
!!deb     vc_r(nil), vc_r(njl), vc_f0(nil), vc_f0(njl)
!!deb  enddo
  
!
  return
end subroutine vc_rh
!
!
subroutine vc_rb(n,vc_rs,iflag, dt)
! vc_xn gives the old time value x^n
  integer :: n, iflag
  double precision :: dt
  double precision :: vc_rs(n)
!
  vc_rs = vc_xs + dt*vc_sm
!
  return
end subroutine vc_rb
!
!
subroutine vc_fn(n,vc_xx,vc_fs,iflag)
implicit none
  integer :: n, iflag
  double precision :: vc_xx(n), vc_fs(n)
  double precision :: vc_r(n)
  double precision :: vc_f(mcoor)
  double precision :: vc_cf(mcoor)
  double precision :: vc_xpt(mpts), vc_ypt(mpts), gp(mpts)
  double precision :: vc_ff((1)*nring*2)
!
  integer :: ncp, nn, nlp, il, nil, njl
  double precision :: tmp
!
  vc_xpt = 0.; vc_ypt = 0.
  call vc_xr(vc_xx,vc_xpt,vc_ypt)
  vc_f = 0.; vc_cf = 0.
  call calcplateletforce(vc_xpt,vc_ypt,vc_f,vc_cf,gp,1)
  vc_ff = 0.
  do nn = 2*nfil+1, npts
    il  = (nn-2*nfil-1)/nring + 1
    nlp = (nn-2*nfil) - (nn-2*nfil-1)/nring*nring 
    ncp = nn - 2*nfil
    nil = nlp + (il-1)*2*nring
    njl = nil + nring
    tmp = zero
    tmp =  tmp + vc_cf(2*nn-1)*vc_nn(nil)*kw(il)/rtc
    tmp =  tmp + vc_cf(2*nn  )*vc_nn(njl)*kw(il)/rtc
    !!tmp =  tmp + vc_cf(2*nn-1)*vc_nn(nil)*kw(il)
    !!tmp =  tmp + vc_cf(2*nn  )*vc_nn(njl)*kw(il)
    vc_ff(nil)= tmp*vc_nn(nil)
    vc_ff(njl)= tmp*vc_nn(njl)
  enddo
!!  vc_fs= (vc_us+vc_sm+vc_ff)
  vc_fs= (vc_sm+vc_ff)
!!  vc_fs= (vc_ff)
!!deb  do nn = 2*nfil+1, npts
!!deb    il  = (nn-2*nfil-1)/nring + 1
!!deb    nlp = (nn-2*nfil) - (nn-2*nfil-1)/nring*nring 
!!deb    ncp = nn - 2*nfil
!!deb    nil = nlp + (il-1)*2*nring
!!deb    njl = nil + nring
!!deb    print '("needprt ", 1024(e14.6,1x))', vc_ff(nil), vc_ff(njl), vc_sm(nil), vc_sm(njl), vc_cf(2*nn-1), vc_sm(2*nn)
!!deb  enddo
!
  return
end subroutine vc_fn
!
!
subroutine vc_xr(vc_xx,vc_xpt,vc_ypt)
  implicit none
  double precision :: vc_xx(nring*2)
  double precision :: vc_xpt(mpts), vc_ypt(mpts)
!
  integer :: ncp, n, nlp, il, nil, njl
 
!
  vc_xpt = xpt; vc_ypt = ypt
  do n=2*nfil+1,npts
    il  = (n-2*nfil-1)/nring + 1
    nlp = (n-2*nfil) - (n-2*nfil-1)/nring*nring 
    ncp = n - 2*nfil
    nil = nlp + (il-1)*2*nring
    njl = nil + nring
    vc_xpt(n) = vc_xx(nil)
    vc_ypt(n) = vc_xx(njl)
  enddo
!
  return
end subroutine vc_xr
!
!
subroutine vc_dirder(n,vc_xx,vc_ww,vc_f,vc_f0,iflag)
  implicit none
  integer :: n, iflag
  double precision :: vc_xx(n), vc_ww(n), vc_f0(n), vc_f(n)
  double precision :: vc_cf(n), vc_df(n), vc_tf(n), del(n)
!
  double precision :: tmp, epsnew, normw, xs
  integer :: i, j, ncp, nn, nlp, il, nil, njl
!
  epsnew = 1.d-8
  normw = sqrt(mydot(n,vc_ww,vc_ww))
  if (normw < 1.d-14) then
    vc_f = 0.
    print '(1x,"direction should not be 0!!", e13.6)',  normw
    goto 123
  endif
  xs = mydot(n,vc_xx,vc_ww)/normw
  if (abs(xs) > 1.d-10) then
    epsnew=epsnew*max(abs(xs),1.d0)*sign(1.d0,xs)
  endif
  epsnew=epsnew/normw
!
  del = vc_xx + epsnew*vc_ww
  call vc_fn(n, del, vc_cf,iflag)
!!!  call vc_fn(n, vc_xx, vc_df,iflag)
  vc_f = (vc_cf - vc_f0)/epsnew
!!  print '(1x, "dirder ",10(e22.14,1x))', (sum(abs(vc_cf-vc_df))),   &
!!   sqrt(sum(abs(vc_cf-vc_f0))), epsnew, normw, xs, (sum(abs(vc_xx-del)))/dble(nring)
!!
!!!
!!  do nn = 2*nfil+1, npts
!!    il  = (n-2*nfil-1)/nring + 1
!!    nlp = (n-2*nfil) - (n-2*nfil-1)/nring*nring 
!!    ncp = n - 2*nfil
!!    nil = nlp + (il-1)*2*nring
!!    njl = nil + nring
!!  enddo
!!
! 
  123 continue
  return
end subroutine vc_dirder
!
double precision function mydot(n,x,y)
  integer :: n
  double precision, dimension(n) :: x, y
!
  integer :: i
  double precision :: tsum 
  tsum = 0.d0

  do i = 1, n
    tsum = tsum + x(i)*y(i)
  enddo
  mydot = tsum
! 
  return
end function mydot
!
!=======================================================================
!
!
subroutine vc_gmrs(n, mm, xc, rs, xo, error, dt)
  implicit none
  integer :: n, mm
  double precision :: xc(n), rs(n), xo(n), error, dt
!
  integer :: infon, iter, imax, iflag
  double precision :: ser0, ser1, serr, tol, dr1, hit, dht, st
  double precision :: ss(n), v(n,mm), vk(n),x11(n), hj(n,2), hq(n,mm), bf(n), r(n),  &
          yy(n), w(n), xx(n)
!
  integer :: i, j, k, j1, j2, mi
  imax = n
  !!tol = 1e-9*norm2(rs)/dble(n)
  tol = 1d-9*norm2(rs)/dble(n)
  !!tol = 1e-6
  bf = rs
  
  iter = 0
! 
  hj = 0.; v = 0.; hq = 0.
  do 10 j = 1, imax
    ser0 = zero
    call vc_matv(n,xc,xo,yy,iflag,dt)
    !!call vc_mv(n,xc,xo,yy,iflag)
    r =  bf - yy
    dr1 = sqrt(mydot(n,r,r))
    if(dr1 .le. 1d-13) then
      print *, '1 iter in IB loc GMRS'
      return
    endif
    
    v(:,1) = r/dr1
    ss = 0.
    ss(1) = dr1
      do i = 1, mm-1
        iter = iter + 1
        vk = v(:,i)

        call vc_matv(n,vk,xo,w,iflag,dt)
        !!call vc_mv(n,vk,xo,w,iflag)
        do k = 1, i
          vk = v(:,k)
          hq(k,i) = mydot(n,w,vk)
          w = w - hq(k,i)*vk
        enddo

        hq(i+1,i) = sqrt(mydot(n,w,w))
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
          dht = sqrt(hq(i,i)*hq(i,i)+hq(i+1,i)*hq(i+1,i))
          hj(i,1) = hq(i,i)/dht
          hj(i,2) = hq(i+1,i)/dht
          st = hj(i,1)*ss(i) + hj(i,2)*ss(i+1)
          ss(i+1) = -hj(i,2)*ss(i) + hj(i,1)*ss(i+1)
          ss(i) = st

          hit = hj(i,1)* hq(i,i) + hj(i,2)*hq(i+1,i)
          hq(i+1,i) = -hj(k,2)*hq(i,i) + hj(i,1)*hq(i+1,i)
          hq(i,i) = hit
        endif

        ser1 = abs(ss(i+1))
        if (ser1 .gt. 1.0d-14) then
           serr = abs((ser1-ser0)/ser1)
        else
           serr = zero
        endif
        if(ser1 .le. tol .or. serr .le. tol) then
          serr = tol - 1.0d-14
          mi = i
          goto 100
        endif

        ser0 = ser1
      enddo
      mi = mm - 1
 100  yy(mi) = ss(mi)/(hq(mi,mi)+1.0d-14)
      do k=mi-1,1,-1
         yy(k) = ss(k)
         do j1 = k+1,mi
           yy(k) = yy(k) - hq(k,j1)*yy(j1)
         enddo
         yy(k) = yy(k)/hq(k,k)
      enddo        
      x11 = xc
      do k=1,mi
         x11 = x11 + yy(k)*v(:,k)
      enddo
      call vc_matv(n,x11,xo,w,iflag,dt)
      !!call vc_mv(n,x11,xo,w,iflag)
      r = bf - w
      ss(mi+1) = sqrt(mydot(n,r,r))
      xc = x11
      if( abs(ss(mi+1)) .lt. tol .or. serr .lt. tol ) then
        print '(1x,"IN force gmres, iter", i5,1x,4(e14.6,1x))', iter, serr, ss(mi+1)
        return
      endif
      error = ss(mi+1)
  10 continue
!  xc = xx

  print *, 'IBmod IN FORCE GMRES, iter', iter, serr

  return
end subroutine vc_gmrs
!
!
subroutine vc_fx(n,vc_xx,vc_fs,iflag)
implicit none
  integer :: n, iflag
  double precision :: vc_xx(n), vc_fs(n)
  double precision :: vc_r(n)
  double precision :: vc_f(mcoor)
  double precision :: vc_cf(mcoor)
  double precision :: vc_xpt(mpts), vc_ypt(mpts), gp(mpts)
  double precision :: vc_ff((1)*nring*2)
!
  integer :: ncp, nn, nlp, il, nil, njl
  double precision :: tmp
!
  vc_xpt = 0.; vc_ypt = 0.
  call vc_xr(vc_xx,vc_xpt,vc_ypt)
  vc_f = 0.; vc_cf = 0.
  call calcplateletforce(vc_xpt,vc_ypt,vc_f,vc_cf,gp,1)
  vc_ff = 0.
  do nn = 2*nfil+1, npts
  !!do nn = 2*nfil+1, 2*nfil+2*nring
    il  = (nn-2*nfil-1)/nring + 1
    !!il  = 1
    nlp = (nn-2*nfil) - (nn-2*nfil-1)/nring*nring 
    nil = nlp + (il-1)*2*nring
    njl = nil + nring
    !!if (il .ne. 1) then 
    !!if (nlp .ne. 0) then 
    !!  print '(1x,"this is il",1024(i5,1x))', il, nn, nn-2*nfil-1, nlp, nil, njl
    !!  print '(1x,1024(e16.8,1x))', vc_cf(nil), vc_cf(njl), vc_cf(2*nn-1), vc_cf(2*nn)
    !!endif
    
    tmp = 0.
    !!oldMaybeRong tmp =  tmp + vc_cf(nil)*vc_nn(nil)*kw(il)
    !!oldMaybeRong tmp =  tmp + vc_cf(njl)*vc_nn(njl)*kw(il)
    tmp =  tmp + vc_cf(2*nn-1)*vc_nn(nil)*kw(il)
    tmp =  tmp + vc_cf(2*nn  )*vc_nn(njl)*kw(il)
    vc_ff(nil)= tmp*vc_nn(nil)
    vc_ff(njl)= tmp*vc_nn(njl)
  enddo
!!  vc_fs= (vc_us+vc_sm+vc_ff)
!!  vc_fs= (vc_sm+vc_ff)
  vc_fs= (vc_ff)
!
  return
end subroutine vc_fx
!
!
subroutine vc_mv(n,vc_xx,vc_xn,vc_rs,iflag,dt)
implicit none
! vc_xn gives the old time value x^n
  integer :: n, iflag
  double precision :: dt
  double precision :: vc_xx(n), vc_xn(n), vc_rs(n)
  double precision :: vc_r(n)
!
!
  call vc_fx(n,vc_xx,vc_r,iflag)

  vc_rs = vc_xx - dt*vc_r
!
  return
end subroutine vc_mv
!
!=======================================================================
end module IBmod
