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

module IBforce
  use, intrinsic :: iso_c_binding
  use parameters


implicit none

!!integer, parameter :: mpts=ment*nring
!!integer, parameter :: npts=mpts
double precision, parameter :: h=hg

private

public :: calcplateletforce

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calcplateletforce(xpt,ypt,fs,cfs,gp,isel)

integer, intent(in) :: isel
integer      :: j
double precision  :: xpt(mpts)
double precision  :: ypt(mpts)
double precision  :: fs(mcoor)
double precision  :: cfs(mcoor)
double precision  :: g(mpts,2)
double precision  :: gp(mpts) !jacobian
!
! now we evaluate force density at these locations:
!
!      write(*,*)' calcplateletforce'
      !!oldForce call grad(xpt,ypt,g,gp)
      
      call ngrad(xpt,ypt,g,gp)
      !!  
!      write(*,*)' after grad max g=',maxval(abs(g))
!
!  contributions to gradient from periodicity conditions (closed IB list)
!
!
! calculate entity force density (force/length)
!
      cfs = 0.
      do j=1,npts
         fs(2*j-1)=-g(j,1)
         fs(2*j)  =-g(j,2)
      enddo
      !!oldIdx do j=2*nfil+1,npts
      do j=1,npts ! compute f*|dx/ds|
         cfs(2*j-1)=-g(j,1)*gp(j)
         cfs(2*j)  =-g(j,2)*gp(j)
      enddo
!      write(*,*)' plateletforce max fs=',maxval(abs(fs))

!!prt      print *, 'force max ', maxval(gp), isel

end subroutine calcplateletforce
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!123456789212345678931234567894123456789512345678961234567897123456789
subroutine ngrad(x,y,g,gp)
!
  integer  :: isel
  integer :: i,im,ip,iend,ione,jent
  double precision  :: ecv,g1,g2,gp1,gp2
  double precision  :: se,sm,xd1,xd2,xm1,xp,xp1, sp
  double precision  :: yd1,yd2,ym1,yp,yp1
  double precision  :: xs1,xq,xa1, ys1,yq,ya1
  double precision  :: g(mpts,2)
  double precision, dimension(mpts) :: gp, x,y
!
  double precision  :: tmp, tp1, tp2, tpx, tpy, ds,rds,rds2, ds2
  double precision, dimension(mpts) :: xx, yy, cx, cy
!
  g = zero
  gp= zero
!
!=======================================================================
!!  find difference of current IB pts and preferred shape
!!  do i = 1, nring
!!    print '(1x,10(e22.16,1x))', x(i), y(i), xbp(i), ybp(i), x(i)-xbp(i), y(i)-ybp(i)
!!  enddo
  
!!prt  print '(1x,10(e22.16,1x))', x(1), y(1), xbp(1), ybp(1), x(1)-xbp(1), y(1)-ybp(1)
!!prt  print *, ''
!!  if (abs(sum(x-xbp))+abs(sum(y-ybp))>1d-10 ) then
    xx= x-xbp; yy= y-ybp
!!  else
!!    xx= x; yy= y
!!  endif
!=======================================================================
!
  ds=cpi*two/dble(nring); rds=one/ds; rds2=rds*rds
!=======================================================================
  do i = 1, nring 
    if (i>1 .and. i< nring) then
      im = i-1; ip = i+1
    elseif (i.eq.1) then
      im = nring; ip = i+1
    elseif (i.eq.nring)then
      im = i-1; ip = 1
    endif
    cx(i)=xx(ip)-two*xx(i)+xx(im)
    cy(i)=yy(ip)-two*yy(i)+yy(im)
  enddo
  cx=cx*rds2; cy=cy*rds2

  tmp = rsl*ds ! rest lenth = rsl
  do i = 1, nring 
    if (i>1 .and. i< nring) then
      im = i-1; ip = i+1
    elseif (i.eq.1) then
      im = nring; ip = i+1
    elseif (i.eq.nring)then
      im = i-1; ip = 1
    endif
    xp=xx(i); xp1=xx(ip); xm1=xx(im)
    yp=yy(i); yp1=yy(ip); ym1=yy(im)
    xq=x(i); xa1=x(ip); xs1=x(im)
    yq=y(i); ya1=y(ip); ys1=y(im)
!
    sp=sqrt((xa1-xq)*(xa1-xq) + (ya1-yq)*(ya1-yq))
    sm=sqrt((xs1-xq)*(xs1-xq) + (ys1-yq)*(ys1-yq))
!
!=======================================================================
! with rest length
    g(i,1)=-rds*((one       )*rds*(xp1-xp)-(one       )*rds*(xp-xm1))
    g(i,2)=-rds*((one       )*rds*(yp1-yp)-(one       )*rds*(yp-ym1))
    g(i,:)=sw(1)*g(i,:)
!=======================================================================
! with bending force
    g(i,1)=g(i,1)+sb(1)*rds2*(cx(im)-two*cx(i)+cx(ip))
    g(i,2)=g(i,2)+sb(1)*rds2*(cy(im)-two*cy(i)+cy(ip))
!=======================================================================
    gp(i)=ds*two/(sqrt(sm)+sqrt(sp))
  enddo
!=======================================================================

!
  return
end subroutine ngrad
!=======================================================================
!123456789212345678931234567894123456789512345678961234567897123456789
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
end module IBforce
