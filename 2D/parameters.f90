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

! all parameters set here except: kc, and active pump for chemical and J_actin
!
module parameters
  implicit none
  double precision, parameter :: cpi=3.141592653589793238D0
  double precision, parameter :: zero=0.d0
  double precision, parameter :: half=0.5d0
  double precision, parameter :: one=1.d0
  double precision, parameter :: two=2.d0
  double precision, parameter :: three=3.d0
  double precision, parameter :: four=4.d0
  double precision, parameter :: eight=8.d0
  double precision :: pi, tupi
!
  integer :: outcount
  integer,parameter:: ment=1                    !max number of IB entities
  integer,parameter:: nfil=0                    !number of points/wall
  integer,parameter:: cent=1
!  integer,parameter:: nring=200                    !number of points/platelet
!  integer,parameter:: nring=80                    !number of points/platelet
!  integer,parameter:: nring=160                    !number of points/platelet
  integer,parameter:: nring=320                    !number of points/platelet
!  integer,parameter:: nring=640                    !number of points/platelet
  integer,parameter:: mpts=ment*nring             !number of IB points
  integer,parameter:: npts=ment*nring             !number of IB points
  integer,parameter:: mcoor=2*mpts                !number of IB coordinates
!=======================================================================
! for Cartesian grid!
  integer,parameter:: l2ny=7                               !log of fluid grid in y-direction
  integer,parameter:: maspect=1                            !aspect ratio
  integer,parameter:: ny=2**l2ny                           !fluid grid sizes
  integer,parameter:: nx=maspect*ny                        !fluid grid sizes
  !!integer,parameter:: ny=32
  !!integer,parameter:: nx=128
  integer,parameter:: nxp2=nx+2                        !fluid grid sizes
  integer,parameter:: nyp2=ny+2                        !fluid grid sizes
!
!=======================================================================
  double precision :: CFL
  double precision :: nu
  double precision :: bk
  double precision :: vco
  !!double precision, parameter :: rtc = 0.36 ! for c(Na+)=145mol/m^3 & 300K
  !!double precision, parameter :: rtc = 1.d0
  double precision, parameter :: rtc = 8.4d3
  double precision, parameter :: rtr = .005d0, rho0=rtc*rtr ! R*T*rho0,scaled with RTC
  ! here rtr = {R*T*rho0}/{R*T*c(Na+)}
!=======================================================================
  double precision,parameter :: xmin = 0.0    ! domain extent
  double precision,parameter :: xmax = 1.0    !,ymin,ymax,xlength  ! 
  double precision,parameter :: ymin = 0.0
  double precision,parameter :: ymax = 1.0
  !!double precision,parameter :: ymax = 0.25
  double precision,parameter :: xlength = xmax-xmin
  double precision,parameter :: ylength = ymax-ymin
  double precision,parameter :: ytop = ymax
  double precision,parameter :: ybot = ymin
  double precision,parameter :: ymid = 0.5*(ybot+ytop)
  double precision,parameter :: hg=xlength/dble(nx)        ! spacestep
  !!double precision :: dlt = CFL*hg
  double precision :: dlt
  !double precision,parameter :: dlt = 0.02
  integer :: ntmax
  integer :: nfreq ! print frequency
!
!=======================================================================
! IB links parameters
  double precision,parameter :: cbstiff=.000000 ! bending stiffness
  double precision,parameter :: sb(1)=cbstiff   ! bending stiffness
  double precision,parameter :: clstiff=1.0     ! elastic modulus
  double precision,parameter :: sw(1)=clstiff   ! elastic modulus
  double precision,parameter :: dl =0.5d0*hg                 ! IB point separation
  !!double precision,parameter :: rsl=4.0*dble(nring/80)*dl
  double precision,parameter :: rsl=zero
  !!double precision,parameter :: rsl=.16
  double precision :: ctild                  !
  double precision :: cs(1)   ! bending angle?
  double precision :: xpt(mpts),ypt(mpts)  !IB coordinates (current coor)
  double precision :: oxpt(mpts),oypt(mpts)!IB coordinates ("old" coor)
  double precision :: xbp(mpts),ybp(mpts)! preferred shape of cell
!=======================================================================
  character(40)   :: runname  = './Data/frun'
!=======================================================================
  double precision, parameter :: xshift=0.5, yshift=0.5 ! cell center for chemicals
  double precision, parameter :: cmsi = 1.0, cpsi = 1.0 ! initial chemicals on two sides
  double precision, parameter :: nmsi = 1.0, npsi = 1.0 ! initial network volumes
  integer, parameter :: nent=3
  integer, parameter :: nib=1
  double precision :: kw(4), kc(4), dif(4)
  integer, parameter :: mglevel=l2ny/2
  integer, dimension(-1:nx+2,-1:ny+2) :: chker0, chker1, kdf, oid, ijd
  integer :: jgd(1:nx*ny)
  integer, dimension(-1:nx+2,-1:ny+2) :: idu,idv, jdu,jdv, idf,id, iid, idn, isf, ivf
  Data kc/1.00, 1.0, 0.0, 0.0/
  Data kw/.100, 0.0, 0.0, 0.0/ !here kw is scaled by RTC to become dimensionless
  Data dif/.01, 1.0, 0.0, 0.0/
  double precision, parameter :: eta=1.d0 , etas=5.d2 , theta0=.02d0
  double precision :: mkp(7,nring,cent), mkq(7,nring,cent), mkr(7,nring,cent)
  double precision :: mko(7,nring,cent), mkn(7,nring,cent)
!!  double precision :: cbo(nring,cent), cbn(nring,cent)
  double precision :: ndv(nring,2,cent) ! normal direction at each IB point
  double precision :: tdv(nring,2,cent) ! tangent direction at each IB point
  double precision :: mtd(nring,2,cent) ! normal direction at each IB point
  double precision :: vc_fe(nring)
  double precision :: cc_pn(-1:nx+2,-1:ny+2)  ! chemical concentration function
  double precision :: cc_pc(-1:nx+2,-1:ny+2)  ! chemical concentration function
  !!double precision :: cc_po(-1:nx+2,-1:ny+2)  ! chemical concentration function
  double precision :: nt_pn(-1:nx+2,-1:ny+2)  ! network volume fraction
  double precision :: nt_pc(-1:nx+2,-1:ny+2)  ! network volume fraction
  logical :: iskip = .false.
  double precision ::  sgma(-1:nx+2,-1:ny+2)  ! network volume fraction
!=======================================================================
  integer :: olen(cent), qlen(cent)
  double precision :: xamin,xamax,yamin,yamax      ! domaina extent
  double precision :: unc(-1:nx+1,-1:ny+1)   ! x velocity for fluid & chemical
  double precision :: vnc(-1:nx+1,-1:ny+1)   ! x velocity for fluid & chemical
!!  double precision :: umc(-1:nx+1,-1:ny+1)   ! x velocity for network
!!  double precision :: vmc(-1:nx+1,-1:ny+1)   ! x velocity for network
!========================================================================

end module parameters
