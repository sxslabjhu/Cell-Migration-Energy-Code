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

MODULE  chemical_mod
  use, intrinsic :: iso_c_binding
! GMRS + Backward Euler + cell center + No flux top/bot & Periodic left/right
  use parameters
  use geometry
  use IBforce
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!========================================================================
  implicit none
!
  PRIVATE 
!
  integer, parameter :: nxb = nx
  integer, parameter :: nyb = ny
!!  integer, parameter :: nst = 8
!!  double precision :: conserv=1.
!!  double precision, parameter :: cmsi = 1.0, cpsi = 1.0 ! initial chemicals on two sides
  double precision, parameter :: fco=half
!!  double precision :: pco = 1.d3, toc, vs= 0.d0, ht=1.00d20

!
  double precision :: vc_xpt(mpts), vc_ypt(mpts), gp(mpts)
  double precision :: vc_cf(mcoor), vc_f(mcoor)
!!      double precision :: vc_fe(nring)
!##
  public :: cc_init, cc_IBbdy, cc_advDiff, vc_cfs, cc_pump

!DDinterface
!DDend interface
!-----B--1----+----2----+----3----+----4----+----5----+----6----+----7-E--+----|
!=======================================================================&
CONTAINS

!##
subroutine cc_init(iaary, llen) ! chemical concentration initialization
  type(iapt), dimension(:) :: iaary(nent-2)
  integer :: llen(nent-2)
  integer :: isel, i, j, jp, i1, i2, k, ic, info, il
  integer     :: strlen,scount, nt
  double precision :: time = 0., tpx, tpy, x, y, dx, dy
  character(40) ibfile  , ffile
!
  double precision :: cmvec(nring), cpvec(nring)
  type(cc_augvar), pointer :: cp_lp
!
  dx = hg; dy = hg
!=======================================================================
!=======================================================================
! Initialize chemical field
  time = 0.d0; 
  do i = 1, nxb ! include corner
  do j = 1, nyb
    x = xmin+(dble(i)-xshift)*dx; y = ymin+(dble(j)-yshift)*dy
!    cc_pc(i,j) = cc_p(x, y, time)
    cc_pc(i,j) = cpsi
  enddo
  enddo
!
!!  cc_po = cc_pc
  cc_pn = cc_pc
!
!!  isel = -1
!!  xpt0 = xpt; ypt0 = ypt;
  oid = 0
  oxpt = xpt; oypt = ypt
!!  call cc_IBbdy(cc_pn, cc_pc, cc_po, isel, time, dt, mydif)
!
  do i = 1, nxb
  do j = 1, nyb
    if (idf(i,j) .eq. -1) then
      cc_pc(i,j) = cmsi
    endif
  enddo
  enddo
!!  cc_po = cc_pc
  cc_pn = cc_pc
!
!!!!  if (isel .eq. -1) then
!!    j = 0
!!    do il = 1, nib
!!      cp_lp => iaary(il)%p
!!      do i = 1, llen(il)
!!        j = j+1
!!        cp_lp%cm = cmsi
!!        !!cp_lp%cm = tmp
!!        j = j+1
!!        cp_lp%cp = cpsi
!!        !!cp_lp%cp = tmp
!!!
!!        cp_lp => cp_lp%next
!!      enddo
!!    enddo
!!!!  endif
!=======================================================================
  do il = 1, nib
    cp_lp => iaary(il)%p
    do i = 1, llen(il)
      cp_lp%cm=cmsi
      cp_lp%cp=cpsi
      cmvec(i) = cmsi
      cpvec(i) = cpsi
      cp_lp => cp_lp%next
    enddo
!!    cmvec = cmsi !noneed, only cp_lp%cm/cp values are used
!!    cpvec = cpsi
    mkp(1,:,il) = mkn(1,:,il)
    !!print *, 'ck here', nring
    !!call smoothcrossing(llen(il),olen(il), nib,il,mkp,iaary,llen,iblst, info)
    call smoothcrossing(olen, nib,il,mkp,iaary,llen, cmvec,cpvec, info)
  enddo
!=======================================================================
! on output, mkp store interpolation of + & - side of grid crossing
  nt = 0
  strlen = len_trim(runname)
  write(ffile,'(2a,i4.4)') runname(1:strlen),'.ps.',nt
  print *, ffile
  open(67,file=ffile,access='stream',action='write')
  write(67)cc_pc
  close(67)
  write(ffile,'(2a,i4.4)') runname(1:strlen),'.chk.',nt
  open(88,file=ffile,form='formatted',action='write')
  do i = -1, nxb+2
    write(88,'(1x,1024(i5,1x))')(idf(i,j),j=-1,nyb+2)
  enddo
  close(88)
  print '(1x,"Done initializing chemical module")'
!=======================================================================
! after calling cc_init, mkp and olen should be set for initial conditions
! olen is not used anymore, so we only need mkp
!=======================================================================
  return
end subroutine cc_init
!
!=======================================================================
! We will use this subroutine to initiate all the computations. So this
! one will call cc_advDiff to advance one time step in chemical
!=======================================================================
subroutine cc_IBbdy(ccpn, ccpc, ibary, iaary, llen, uin,isel, time, dt, mydif)
  implicit none
  type(ibpt), dimension(:), intent(in) :: ibary(nent-2)
  type(iapt), dimension(:), intent(in) :: iaary(nent-2)
  integer :: llen(nent-2)
  double precision, dimension(:,:,:) :: uin(2,nring,nent-2)
  integer :: isel
  double precision :: time, dt, mydif
  double precision :: ccpn(-1:nxb+2,-1:nyb+2)  ! chemical concentration function
  double precision :: ccpc(-1:nxb+2,-1:nyb+2)  ! chemical concentration function
!!  double precision :: ccpo(-1:nxb+2,-1:nyb+2)  ! chemical concentration function
!!  double precision :: ccer(-1:nxb+2,-1:nyb+2)  ! chemical concentration function
!!  double precision :: ccrs(-1:nxb+2,-1:nyb+2)  ! chemical concentration function
!
  integer :: i, j, k, ik, jk, il, nt, info, strlen
  type(cc_bp), pointer :: curr
  type(cc_augvar), pointer :: cp_lp, cp_lq
  double precision ::  dx, dy, dt2, eps, rh, tmp, tp1, tp2, tp3, tp4, tpp, s0, convexity
!!  integer :: ip, jp, im, jm, is, js, nw, ns, it, jt, i1, i2
  double precision :: cbq(nring)
  character(40) :: efile, ffile
  double precision, dimension(:,:,:) :: mk(7,nring,nent-2), nb(nring,2, nent-2)
!!  logical :: fl(5), as, alll, allr
!!  logical :: al(nent-2)
!!  character :: equed
  double precision :: cmp(nring,5,nent-2)
!
  eps = 1.d-8
  dt2 = 0.5*dt; dx = hg; dy = hg
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
  mk = mkn
  vc_xpt = xpt; vc_ypt = ypt
  call vc_cfs(vc_xpt,vc_ypt, vc_fe) ! this can't be used in two phase model

! uin stores dx/dt-u at each IB point
!///////////////////////////////////////////////////////////////////////

  call curvatureInt(mkn(:,:,1:nent-2),convexity)

  if (isel .eq. -1) then
    goto 5555
  endif
!
  do il = 1, nib ! 
    cbq = uin(1,:,il); ! uin stores the speed of -(dx/dt-u)
    cp_lp => iaary(il)%p
    do i = 1, llen(il)
      s0 = cp_lp%s; k = cp_lp%k
      call brs(nring,nib,il,mk,cbq  ,s0,tmp, cp_lp%k,time,0)
      cp_lp%vdn = tmp ! interpolate to grid crossing
      tp1 = paraval(s0 , nring, mkq, il, 1, 0, info) 
      tp2 = paraval(s0 , nring, mkr, il, 1, 0, info) 
!!      cp_lp%udn = (tp1+tp2)*kw(il)/cp_lp%dxs
      call brs(nring,nib,il,mk,vc_fe,s0,tp1, cp_lp%k,time,0)
      cp_lp%fwF = tp1/rtc ! interpolate F to grid crossing
      !!cp_lp%fwF = 0. ! interpolate F to grid crossing, assuming no elastic force
      !!print '(1x, 4(e13.6,1x), 2(i3,1x))', tmp, cp_lp%s, s0, cbq(k), cp_lp%k, i
      tmp = paraval(s0, nring, mkp, il, 0, 0, info) ! on - side old time chemical
      cp_lp%ocm = tmp*kw(il) ! store c^-*kw(t_n)
      !!cp_lp%ocm = rtc*tmp*kw(il) ! store c^-*kw(t_n)
      tmp = paraval(s0, nring, mkp, il, 1, 0, info) ! on + side old time chemical
      cp_lp%ocp = tmp*kw(il) ! store c^+*kw(t_n)
      !!cp_lp%ocp = rtc*tmp*kw(il) ! store c^+*kw(t_n)
      !
      cp_lp => cp_lp%next
    enddo
  enddo

!///////////////////////////////////////////////////////////////////////
!=======================================================================
! Advance advection diffusion system one step
  5555 continue
  info = 0
  call cc_advDiff(ccpn,ccpc, idn, isel, time, mko, iaary,llen, info,dt,mydif)

!=======================================================================
!--+-B--1----+----2----+----3----+----4----+----5----+----6----+----7-E

    do il = 1, nib ! save lagrangian chemicals on the new time level
      curr => ibary(il)%p
      cp_lp=>iaary(il)%p
      do i = 1, nring
        tmp = curr%s
        tp1 = paraval(tmp, nring, mkp,il,0,0,info);
        tp2 = paraval(tmp, nring, mkp,il,1,0,info);
        cmp(i,1,il) = tp1
        cmp(i,2,il) = tp2
        !cmp(i,1,il) = cc_cm(i,il)
        !cmp(i,2,il) = cc_cp(i,il)
        cmp(i,3,il) = curr%x
        cmp(i,4,il) = curr%y
        !cmp(i,3,il) = mkp(2,i,il)
        !cmp(i,4,il) = mkp(5,i,il)
        !!cmp(i,1,il) = curr%vx+curr%ux
        !!cmp(i,2,il) = uin(1,i,il)
        !!cmp(i,3,il) = tmp
        !!cmp(i,4,il) = cp_lp%s
        !!cmp(i,5,il) = curr%uy
        cmp(i,5,il) = curr%s
        !!if (i > llen(il)) cmp(i,5,il) =0.
        !!cp_lp =>cp_lp%next
!   
        curr => curr%prev
      enddo
    enddo
!!
if (isel .ne. -1) then
    nt = nint(time/dt)
    strlen = len_trim(runname)
  if (nt-nfreq*(nt/nfreq) .eq. 0) then
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
    write(ffile,'(2a,i4.4)') runname(1:strlen),'.ps.',outcount+1
    open(67,file=ffile,access='stream',action='write')
    write(67)ccpc
    close(67)
  endif
else
    strlen = len_trim(runname)
    write(efile,'(2a,i4.4)') runname(1:strlen),'.chk.',outcount+1
    open(88,file=efile,form='formatted',action='write')
    do i = -1, nxb+2
      write(88,'(1x,1024(i5,1x))')(idf(i,j),j=-1,nyb+2)
    enddo
    close(88)
    write(efile,'(2a,i4.4)') runname(1:strlen),'.cib.', 0
    open(77,file=efile,access='stream',action='write')
    write(77)cmp
    close(77)
endif
!=======================================================================
  tp1 = sqrt(mydot(ccpc,ccpc))*hg; tp2 = sqrt(mydot(ccpn,ccpn))*hg
  tp3 = sum(xpt)/dble(nring); !tp4 = maxval(xpt); tpp = minval(xpt)
  tp4 = sum(ypt)/dble(nring)
  print '(1x,"Time: ", 10(e16.8,1x))', time, tp2, tp3, tp4, convexity
!=======================================================================

!=======================================================================
! deallocate memorys here
 4001 format(1x,1024(e18.6,1x))
!TT
  return
end subroutine cc_IBbdy
!
!=======================================================================
! Advance adv-diff equation, which is solved by a fractional step
! method. The advection step assumes the variable p(si) is defined ever-
! ywhere (even inside the moving cell)
!=======================================================================
subroutine cc_advDiff(ccpn, ccpc, idn, isel, time, mkk, iaary, llen, info,dt,mydif)
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
  double precision :: ccpn(-1:nxb+2,-1:nyb+2)  ! chemical concentration new
  double precision :: ccpc(-1:nxb+2,-1:nyb+2)  ! chemical concentration old(current)
  double precision :: ccps(-1:nxb+2,-1:nyb+2)  ! source term
!
  integer :: i, j, k, il, ierr, iter, it, jt, ik, jk, im(2), ip(2)
  type(cc_augvar), pointer :: cp_lp
  double precision, dimension(-1:nxb+2,-1:nyb+2) :: rg, xx
  double precision :: dx, dy, eps, resid, normb, tmp
  integer :: nl
  double precision :: cmvec(nring), cpvec(nring)
!
  double precision :: ra(nxb*nyb/2), rb(nxb*nyb/2) !may have problem in definition
!
  if ( isel .eq. -1) then
    goto 9999
  endif

  nl = 0
  do il = 1, nib ! # of aug variables @ interface
    nl = nl + 2*llen(il)
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
    ccpc(i,0)  = ccpc(i,1)!! + h*cc_bto(i,1)
    ccpn(i,0)  = ccpn(i,1)!! + h*cc_btn(i,1)
!    y = ymax top, no flux BC
!!rmtemp    ccpo(i,nyb+1)= ccpo(i,nyb)!! + h*cc_btt(i,2)
    ccpc(i,nyb+1)= ccpc(i,nyb)!! + h*cc_bto(i,2)
    ccpn(i,nyb+1)= ccpn(i,nyb)!! + h*cc_btn(i,2)
  enddo

  do j = 1, nyb! set periodic @left/right ghost, not used in linear solve
!    x = xmin left, periodic BC
    i = 0
!!rmtemp    ccpo(i,j) = ccpo(nxb,j)
    ccpc(i,j) = ccpc(nxb,j)
    ccpn(i,j) = ccpn(nxb,j)
!    x = xmax right, periodic BC
    i = nxb+1
!!rmtemp    ccpo(i,j) = ccpo(1,j)
    ccpc(i,j) = ccpc(1,j)
    ccpn(i,j) = ccpn(1,j)
  enddo
!========================================================================
   
  !!rg(1:nxb, 1:nyb) = ccpc(1:nxb, 1:nyb) + dlt*ccps(1:nxb, 1:nyb) 
  rg(1:nxb, 1:nyb) = ccpc(1:nxb, 1:nyb)!! + dlt*ccps(1:nxb, 1:nyb) 
!========================================================================
!
  isf = 0
  do il = 1, nib
    cp_lp =>iaary(il)%p
    do i = 1, llen(il)
      call cc_setRHS(rg, cp_lp, llen, il, isf, idn, mkk, -1, iaary, info, dt)
!
      call cc_setRHS(rg, cp_lp, llen, il, isf, idn, mkk,  1, iaary, info, dt)
!
      cp_lp => cp_lp%next
    enddo
  enddo
!=======================================================================
  
!!  nt = nint(time/dlt)
!!  strlen = len_trim(runname)
  normb = sqrt(mydot(rg,rg))
  if (abs(normb) < eps) then
    xx = 0.d0
    iter = 1
    resid = 0.
    ccpn = xx
    print *, 'GET ZERO here!!'
    return
  else
    j = 0
    do il = 1, nib
      ! set RHS for aux variabls defined @ grid crossing
      cp_lp => iaary(il)%p
      do i = 1, llen(il)
        tmp = cc_pump(cp_lp%sl,time)/cp_lp%dxs! constpmp
        !!tmp = 0. ! pump is concentration dependent

        ! on "-" side
        j = j + 1
        !!tmp = 0. ! no flux interface cond!
        ra(j) = tmp
        !!velo ra(j) = tmp + cp_lp%fwF*cp_lp%ocm ! explicit treatment of force in fw
        !!cp_lp%rs_m = tmp + cp_lp%fwF*cp_lp%ocm ! explicit treatment of force in fw
        !!ra(j) =-kp(il) ! pump is constant
        !!cp_lp%rs_m =-kp(il) ! pump is constant
!
        ! on "+" side
        j = j + 1
        !!tmp = 0.
        !!ra(j) = tmp
        ra(j) = tmp ! constpmp
        !!velo ra(j) = tmp + cp_lp%fwF*cp_lp%ocp ! explicit treatment of force in fw
        !!cp_lp%rs_p = tmp + cp_lp%fwF*cp_lp%ocp ! explicit treatment of force in fw
        !!ra(j) = kp(il)
        !!cp_lp%rs_p = kp(il)
!
        cp_lp => cp_lp%next
      enddo
    enddo
!=======================================================================
    xx = ccpc
    !!xx = 0.
    !!prt print *, 'Are we here? before solver'
    call cc_solver(rg,ra(1:nl), xx, rb(1:nl), nl, iaary, llen,info, isel,time,dt,mydif)
    ccpn = xx
    j = 0
    do il = 1, nib
      !save aux chemicals @ grid crossing
      cp_lp => iaary(il)%p
      do i = 1, llen(il)
        j = j+1
        cp_lp%cm = rb(j)
        j = j+1
        cp_lp%cp = rb(j)
!
        cp_lp => cp_lp%next
      enddo
    enddo
!========================================================================
  endif
!
9999 continue
!!noneed  if (isel .eq. -1) then
!!noneed    j = 0
!!noneed    do il = 1, nib
!!noneed      cp_lp => iaary(il)%p
!!noneed      do i = 1, llen(il)
!!noneed!!        tmp = cc_p(cp_lp%xb,cp_lp%yb, time)
!!noneed        j = j+1
!!noneed        cp_lp%cm = cmsi
!!noneed        !!cp_lp%cm = tmp
!!noneed        j = j+1
!!noneed        cp_lp%cp = cpsi
!!noneed        !!cp_lp%cp = tmp
!!noneed!
!!noneed        cp_lp => cp_lp%next
!!noneed      enddo
!!noneed    enddo
!!noneed  endif
!=======================================================================
  do il = 1, nib
    cp_lp => iaary(il)%p
    do i = 1, llen(il)
      !!tmp = cc_p(cp_lp%xb,cp_lp%yb, time)
      !!cmvec(i) = tmp
      !!cpvec(i) = tmp
      cmvec(i) = cp_lp%cm
      cpvec(i) = cp_lp%cp
      cp_lp => cp_lp%next
    enddo
    mkp(1,:,il) = mkn(1,:,il)
    !!print *, 'ck here', nl
    !!call smoothcrossing(llen(il),olen(il), nib,il,mkp,iaary,llen,iblst, info)
    call smoothcrossing(olen, nib,il,mkp,iaary,llen, cmvec,cpvec, info)
  enddo
!=======================================================================
! on output, mkp store interpolation of + & - side of grid crossing
!=======================================================================
!
  return
end subroutine cc_advDiff
!
!=======================================================================
subroutine f2m(x, M, nl, ux, uy, iaary, llen, info, dt, mydif)
! f2m subroutine will project value at grid crossing to corresponding
! im or ip index. Note variable x should have exact order as in iaary
! since we will need delta info. The order if "-" first then "+" side
! pointer will continue on the "next" direction
!
  !!double precision :: M(-1:nxb+2,-1:nyb+2)  ! output vector
  double precision :: dt, mydif
  double precision :: M(1:nx*ny)
  double precision :: ux(1:nx*ny), uy(1:nx*(ny+1))
  type(iapt), dimension(:) :: iaary(nent-2)
  integer :: nl, info, llen(nent-2)
  double precision, dimension(nl) :: x
!
  integer :: i, k, il, ic, im(2), ip(2), gm(2)
  double precision :: tmp, coef, delta, utx, vtx, fm(2), uca, vca
  type(cc_augvar), pointer :: cp_lp
!
  !!coef =-dt*dif(1)/hg/hg  ! Euler
  coef =-dt*mydif/hg/hg  ! Euler
  utx = 0.5d0*dt/hg; vtx = utx
!
  M = 0.
  do il = 1, nib
    cp_lp => iaary(il)%p
    ic = 0
    do i = 1, llen(il)
      delta = (cp_lp%delm); im = cp_lp%im; ip = cp_lp%ip; gm = ip-im; fm = dble(gm)
      !!tmp = 6./(1.+delta)/(2.+delta)*coef
      ic = ic + 1
      call sub2ind(nxb,im(1),im(2),k)
      !!M(im(1),im(2)) = tmp*x(ic)
      !!M(im(1),im(2)) = M(im(1),im(2))+ coef*cp_lp%st_m(1)*x(ic)
      !!old M(k) = M(k)+ coef*cp_lp%st_m(1)*x(ic)
      !!M(k) = M(k)+ (delta*(utx*ux(k)*fm(1)+vtx*uy(k)*fm(2))+coef)*cp_lp%st_m(1)*x(ic)
      !!M(k) = M(k)+ coef*cp_lp%st_m(1)*cmsi ! debugging mode
      !!M(k) = M(k)+ (delta*(utx*ux(k)*fm(1)+vtx*uy(k)*fm(2))+coef)*cp_lp%st_m(1)*x(ic) ! here ux is @ cell center
      uca = half*(ux(k)+ux(k+1)); vca = half*(uy(k)+uy(k+nx))
      utx = dt/hg/(delta+half);  vtx = dt/hg/(delta+half)

      M(k) = M(k)+ ((utx*uca*fm(1)+vtx*vca*fm(2))+coef*cp_lp%st_m(1))*x(ic) ! here ux is @ cell center

      !!delta = (cp_lp%delp); im = cp_lp%im; ip = cp_lp%ip; gm = im-ip; fm = dble(gm)
      delta = (cp_lp%delp);                               gm = im-ip; fm = dble(gm)
      !!tmp = 6./(1.+delta)/(2.+delta)*coef
      ic = ic + 1
      call sub2ind(nxb,ip(1),ip(2),k)
      !!M(ip(1),ip(2)) = tmp*x(ic)
      !!M(ip(1),ip(2)) = M(ip(1),ip(2))+ coef*cp_lp%st_p(1)*x(ic)
      !!old M(k) = M(k)+ coef*cp_lp%st_p(1)*x(ic)
      !!M(k) = M(k)+ (delta*(utx*ux(k)*fm(1)+vtx*uy(k)*fm(2))+coef)*cp_lp%st_p(1)*x(ic)
      !!M(k) = M(k)+ coef*cp_lp%st_p(1)*cpsi ! debugging mode
      uca = half*(ux(k)+ux(k+1)); vca = half*(uy(k)+uy(k+nx))
      utx = dt/hg/(delta+half);  vtx = dt/hg/(delta+half)
      M(k) = M(k)+ ((utx*uca*fm(1)+vtx*vca*fm(2))+coef*cp_lp%st_p(1))*x(ic) ! here ux is @ cell center
      
      !!print '(1x, "M(-)= ", e16.8, " M(+) ", 10(e16.8,1x))', M(im(1),im(2)), M(ip(1),ip(2)), tmp, x(ic)
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
  double precision, intent(in out) :: Y(1:nxb*nyb)  ! output vector
  type(iapt), dimension(:) :: iaary(nent-2)
  integer :: nl, info, llen(nent-2)
  double precision, dimension(nl) :: x
  double precision :: mydif
!
  integer :: il, i, j, k, ic, im(2), ip(2), g_ali(4,2), g_off(4,2)
  integer :: cpl, gpl, opl
  double precision :: tmp, coef, delta, c_ali(4), c_off(4)
  type(cc_augvar), pointer :: cp_lp
!
  do il = 1, nib
    cp_lp => iaary(il)%p
    ic = 0
    do i = 1, llen(il)
      g_ali = cp_lp%i_ali_grd_m; c_ali = mydif*cp_lp%c_ali_grd_m; 
      g_off = cp_lp%i_off_grd_m; c_off = mydif*cp_lp%c_off_grd_m; 
      cpl   = cp_lp%iml; gpl = g_ali(4,1); opl = g_off(4,1)
!
      tmp = 0.
      !!do j = 1, cpl
      do j = 1, gpl ! note c_ali(j+1) = 0 if g_ali(j,:) is not a stencil pt
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
      x(ic) = tmp
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
      g_ali = 0; g_off = 0; c_ali = 0.; c_off = 0.;
      g_ali = cp_lp%i_ali_grd_p; c_ali = mydif*cp_lp%c_ali_grd_p; 
      g_off = cp_lp%i_off_grd_p; c_off = mydif*cp_lp%c_off_grd_p; 
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
      x(ic) = tmp
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
  double precision::time, mydif
!
  integer :: il, i, j, ic, im(2), ip(2), g_ali(4,2), g_off(4,2), cpl
  double precision :: tmp, tp1, coef, delta, c_ali(4), c_off(4), c_alim(4), c_offm(4)
  double precision :: xt, yt, tp0, s0, sn, sa
  type(cc_augvar), pointer :: cp_lp
!
  vmx = 0.
  do il = 1, nib
    cp_lp => iaary(il)%p
    ic = 0
    do i = 1, llen(il)
      s0 = cp_lp%s; 
      sn = paraval(s0, nring, mkq, il, 1, 0, info)
      sa = paraval(s0, nring, mkr, il, 1, 0, info)
      tp0 = cc_pump(cp_lp%sl,time)/cp_lp%dxs ! pump is concentration dependent
      tp0 = 0. ! constpmp

      ic = ic + 1
      c_ali = mydif*cp_lp%c_ali_grd_m; c_off = mydif*cp_lp%c_off_grd_m
      !!constpmp tmp = kc(il) + c_off(1) + c_ali(1) + cp_lp%vdn !should be same in solver for prem 
      !!constpmp vmx(ic) = tmp*x(ic) - kc(il)*x(ic+1) 
      tmp =-kc(il) + (c_off(1) + c_ali(1)) - cp_lp%vdn !should be same in solver for prem 
      vmx(ic) = tmp*x(ic) + kc(il)*x(ic+1) +max(tp0,0.)*x(ic+1) + min(tp0,0.)*x(ic)
      !!prt print '(1x, "- side ", 1024(e16.8,1x))', c_ali, c_off, tmp, x(ic), vmx(ic)
      !!velo tmp =-kc(il) + (c_off(1) + c_ali(1))+cp_lp%udn+ cp_lp%ocm !should be same in solver for prem 
      !!velo vmx(ic) = tmp*x(ic) + (kc(il)-cp_lp%ocm)*x(ic+1) +max(tp0,0.d0)*x(ic+1) + min(tp0,0.d0)*x(ic)

      ic = ic + 1
      c_ali = mydif*cp_lp%c_ali_grd_p; c_off = mydif*cp_lp%c_off_grd_p
      !!constpmp tp1 = kc(il) + c_off(1) + c_ali(1) - cp_lp%vdn !should be same in solver for prem 
      !!constpmp vmx(ic) = tp1*x(ic) - kc(il)*x(ic-1)
      tp1 = kc(il) + (c_off(1) + c_ali(1)) - cp_lp%vdn !should be same in solver for prem 
      vmx(ic) = tp1*x(ic) - kc(il)*x(ic-1)+min(tp0,0.)*x(ic-1) + max(tp0,0.)*x(ic)
      !!prt print '(1x, "+ side ", 1024(e16.8,1x))', c_ali, c_off, tmp, x(ic), vmx(ic)
      !!print '(1x, " coefficient ", 1024(e16.8,1x))', tmp, tp1
      !!velo tp1 = kc(il) + (c_off(1) + c_ali(1))+cp_lp%udn - cp_lp%ocp !+tp0!should be same in solver for prem 
      !!velo vmx(ic) = tp1*x(ic) -(kc(il)-cp_lp%ocp)*x(ic-1)+min(tp0,0.)*x(ic-1) + max(tp0,0.)*x(ic)
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
subroutine nmatv(x, Ax, ux, uy, iaary, llen, time, info, dt, mydif)
! here the vector is reshaped to 1d, in the order of storing each 
! row, which is
! {1(1:...:nxb) ... j(1:.i.:nxb) ... nyb(1:...:nxb)}
  integer :: llen(nent-2), info, nl
  double precision :: time, dt, mydif
  double precision, dimension(1:nxb*nyb) :: x, Ax
  double precision :: ux(1:nx*ny), uy(1:nx*(ny+1))
  double precision :: uu((nxb+1)*nyb)
  double precision :: vv((nyb+1)*nxb)
  type(iapt), dimension(:) :: iaary(nent-2)
!
  integer :: il, i, j, k, l, ic, g_ali(4,2), g_off(4,2), cpl
  integer :: ik, jk, jm(2), jp(2), is(2), jc(2), iloc(4,2)
  integer :: imv(-1:nxb+2,-1:nyb+2)
  double precision :: tmp, coef, con, cox, coy, delta, c_ali(4), c_off(4), st(4)
  type(cc_augvar), pointer :: cp_lp
  double precision :: x1, y1, utx, vtx, uca, vca, fm(2), tp1, tp2
  integer :: gm(2), km(2)
! 
  con = mydif*dt/hg/hg; cox = con; coy = con;
  utx = 0.5d0*dt/hg; vtx = utx
  Ax = 0.
!
  l = nxb*nyb
!
!
!========================================================================
  Ax(cip)=(1.+4.*con+utx*(ux(cip+1)-ux(cip))+vtx*(uy(cip+nx)-uy(cip)))*x(cip)&
      + ( utx*ux(cip+1 )-cox)*x(cip+1 ) + (-utx*ux(cip)-cox)*x(cip-1 )  &
      + ( vtx*uy(cip+nx)-coy)*x(cip+nx) + (-vtx*uy(cip)-coy)*x(cip-nx)

!========================================================================
! bottom: Dirichlet BC (need change corresponding RHS part)
  Ax(cbt)=(1.+4.*con+utx*(ux(cbt+1)-ux(cbt))+vtx*(uy(cbt+nx)-uy(cbt)))*x(cbt)&
    + ( utx*ux(cbt+1 )-cox)*x(cbt+1 )  + (-utx*ux(cbt)-cox)*x(cbt-1)      &
    + ( vtx*uy(cbt+nx)-coy)*x(cbt+nx)  - (-vtx*uy(cbt)-coy)*x(cbt  )
!========================================================================
! top Dirichlet BC (need change corresponding RHS part) 
  Ax(ctp)=(1.+4.*con+utx*(ux(ctp+1)-ux(ctp))+vtx*(uy(ctp+nx)-uy(ctp)))*x(ctp)&
      + ( utx*ux(ctp+1 )-cox)*x(ctp+1)  + (-utx*ux(ctp)-cox)*x(ctp-1)   &
      - ( vtx*uy(ctp+nx)-coy)*x(ctp  )  + (-vtx*uy(ctp)-coy)*x(ctp-nx)
!========================================================================
! left periodic
  Ax(clf)=(1.+4.*con+utx*(ux(clf+1)-ux(clf))+vtx*(uy(clf+nx)-uy(clf)))*x(clf)&
      + ( utx*ux(clf+1 )-cox)*x(clf+1 ) + (-utx*ux(clf)-cox)*x(crt   )  &
      + ( vtx*uy(clf+nx)-coy)*x(clf+nx) + (-vtx*uy(clf)-coy)*x(clf-nx)
!========================================================================
! right periodic
  Ax(crt)=(1.+4.*con+utx*(ux(clf  )-ux(crt))+vtx*(uy(crt+nx)-uy(crt)))*x(crt)&
      + ( utx*ux(clf   )-cox)*x(clf   ) + (-utx*ux(crt)-cox)*x(crt-1 )  &
      + ( vtx*uy(crt+nx)-coy)*x(crt+nx) + (-vtx*uy(crt)-coy)*x(crt-nx)
                                
!========================================================================
! bottom left corner: Dirichlet BC (need change corresponding RHS part)
  Ax(cbl)=(1.+4.*con+utx*(ux(cbl+1)-ux(cbl))+vtx*(uy(cbl+nx)-uy(cbl)))*x(cbl)&
      + ( utx*ux(cbl+1 )-cox)*x(cbl+1 ) + (-utx*ux(cbl)-cox)*x(cbr   )  &
      + ( vtx*uy(cbl+nx)-coy)*x(cbl+nx) - (-vtx*uy(cbl)-coy)*x(cbl   )
!========================================================================
!!!right, bot: Dirichlet
  Ax(cbr)=(1.+4.*con+utx*(ux(cbl  )-ux(cbr))+vtx*(uy(cbr+nx)-uy(cbr)))*x(cbr)&
      + ( utx*ux(cbl   )-cox)*x(cbl   ) + (-utx*ux(cbr)-cox)*x(cbr-1 )  &
      + ( vtx*uy(cbr+nx)-coy)*x(cbr+nx) - (-vtx*uy(cbr)-coy)*x(cbr   )
!========================================================================
!! left, top: Dirichlet
  Ax(ctl)=(1.+4.*con+utx*(ux(ctl+1)-ux(ctl))+vtx*(uy(ctl+nx)-uy(ctl)))*x(ctl)&
      + ( utx*ux(ctl+1 )-cox)*x(ctl+1 ) + (-utx*ux(ctl)-cox)*x(ctr   )  &
      - ( vtx*uy(ctl+nx)-coy)*x(ctl   ) + (-vtx*uy(ctl)-coy)*x(ctl-nx)
!========================================================================
!!!right, top: Dirichlet
  Ax(ctr)=(1.+4.*con+utx*(ux(ctl  )-ux(ctr))+vtx*(uy(ctr+nx)-uy(ctr)))*x(ctr)&
      + ( utx*ux(ctl   )-cox)*x(ctl   ) + (-utx*ux(ctr)-cox)*x(ctr-1 )  &
      - ( vtx*uy(ctr+nx)-coy)*x(ctr   ) + (-vtx*uy(ctr)-coy)*x(ctr-nx)
!=======================================================================
  imv = 0
  do il = 1, nib
    cp_lp => iaary(il)%p
    do i = 1, llen(il) ! assume irregular cells are 1 grid away from ctb/cbt, clf/crt
    !!do i = 1, 2
      jm = 0; jp = 0; iloc = 0; cpl = 0; st = 0.; fm = 0.
      jm = cp_lp%im; jp = cp_lp%ip; cpl = cp_lp%iml; gm = jp-jm; fm = dble(gm)
      delta = cp_lp%delm; iloc = cp_lp%ilocm; st = cp_lp%st_m
      !!prt print '(1x,1024(e16.8,1x))', cp_lp%st_m
      call sub2ind(nx, jm(1), jm(2), k) 
      ic = k
      call sub2ind(nx, jp(1), jp(2), j) 
      !!prt print *, 'Ax(ic) ', Ax(ic), ' x(ic) ', x(ic), ic, con
      !!diffusion only Ax(k) = Ax(k) + con*x(j) ! add contribution from ghost point back
!========================================================================
!!ref  Ax(cip)=(1.+4.*con+utx*(ux(cip+1)-ux(cip))+vtx*(uy(cip+nx)-uy(cip)))*x(cip)&
!!ref      + ( utx*ux(cip+1 )-cox)*x(cip+1 ) + (-utx*ux(cip)-cox)*x(cip-1 )  &
!!ref      + ( vtx*uy(cip+nx)-coy)*x(cip+nx) + (-vtx*uy(cip)-coy)*x(cip-nx)
!========================================================================
      !!Ax(k) = Ax(k) - (utx*ux(k+ik)*fm(1)+vtx*uy(k+jk)*fm(2)-con)*x(j) ! add values from ghost point back
      tp1 = utx*(ux(ic+1 )*(x(ic+1 )+x(ic))-ux(ic)*(x(ic)+x(ic-1 )))
      tp2 = vtx*(uy(ic+nx)*(x(ic+nx)+x(ic))-uy(ic)*(x(ic)+x(ic-nx)))
      Ax(k) = Ax(k) - ((tp1*abs(fm(1))+tp2*abs(fm(2)))-con*x(j)) ! add values from ghost point back
      !!Ax(k) = Ax(k) - (utx*ux(k+ik)*fm(1)+vtx*uy(k+jk)*fm(2)-con)*x(j) !new, add values from ghost point back
      !!Ax(k) = Ax(k) - (utx*ux(k)*fm(1)+vtx*uy(k)*fm(2)-con)*x(j) ! add values from ghost point back
      !!prt print *, 'Ax(ic) ', Ax(ic), ' diff', con*x(j), x(j), j
      tmp= 0.
      do j = 1, cpl
        call sub2ind(nxb, iloc(j,1), iloc(j,2), k) 
        tmp = tmp + x(k)*st(j+1)
      enddo
      !!diffusion only tmp =-tmp*con
      !!tmp =tmp*(utx*ux(ic)*fm(1)+vtx*uy(ic)*fm(2)-con)
      tmp =tmp*(-con)
      uca = half*(ux(ic)+ux(ic+1)); vca = half*(uy(ic)+uy(ic+nx))
      tp1 = -fm(1)*uca*half*(x(ic-gm(1))+x(ic))
      tp2 = -fm(2)*vca*half*(x(ic-gm(2)*nx)+x(ic))
      tp1 = tp1*dt/hg/(delta+half)
      tp2 = tp2*dt/hg/(delta+half)
      tmp = tmp +  tp1 + tp2 
      !!prt print *, 'tmp ', tmp, con
      !!if (abs(idn(jm(1),jm(2))) .eq. 3) then ! attention here
      if ((isf(jm(1),jm(2))) .eq. -1.and. imv(jm(1),jm(2)).eq. 0) then ! attention here
        if (cp_lp%mi_st .eq. 1) then
          call sub2ind(nxb, cp_lp%m_is(1), cp_lp%m_is(2), k)
          tmp = tmp + x(k)*cp_lp%stm_is(2)
          tmp = tmp + x(ic)*cp_lp%stm_is(1)
        endif
        if (cp_lp%mj_st .eq. 1) then
          call sub2ind(nxb, cp_lp%m_js(1), cp_lp%m_js(2), k)
          tmp = tmp + x(k)*cp_lp%stm_js(2)
          tmp = tmp + x(ic)*cp_lp%stm_js(1)
        endif
        imv(jm(1),jm(2)) = -1
      endif

      Ax(ic) = Ax(ic) + tmp
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
      Ax(k) = Ax(k) - ((tp1*abs(fm(1))+tp2*abs(fm(2)))-con*x(j)) ! add values from ghost point back
      !!Ax(k) = 
      !!oldAx(k) = Ax(k) - (utx*ux(k)*fm(1)+vtx*uy(k)*fm(2)-con)*x(j) ! add values from ghost point back
      tmp= 0.
      do j = 1, cpl
        call sub2ind(nxb, iloc(j,1), iloc(j,2), k) 
        tmp = tmp + x(k)*st(j+1)
      enddo
      !!diffusion only tmp =-tmp*con
      !!tmp =tmp*(utx*ux(ic+ik)*fm(1)+vtx*uy(ic+jk)*fm(2)-con) ! conservative form
      tmp =tmp*(-con)
      uca = half*(ux(ic)+ux(ic+1)); vca = half*(uy(ic)+uy(ic+nx))
      tp1 = -fm(1)*uca*half*(x(ic-gm(1))+x(ic))
      tp2 = -fm(2)*vca*half*(x(ic-gm(2)*nx)+x(ic))
      tp1 = tp1*dt/hg/(delta+half)
      tp2 = tp2*dt/hg/(delta+half)
      tmp = tmp + tp1 + tp2
      !!if (abs(idn(jp(1),jp(2))).eq. 3) then !attention here
      if ((isf(jp(1),jp(2))).eq. 1.and. imv(jp(1),jp(2)).eq. 0) then !attention here
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
        imv(jp(1),jp(2)) = 1
      endif

      Ax(ic) = Ax(ic) + tmp
!
      cp_lp => cp_lp%next
    enddo
  enddo
!
!=======================================================================
!
  return
end subroutine nmatv
!
!=======================================================================
subroutine amatv(x, Ax, nl, ux, uy, iaary, llen, time, info, dt, mydif)
! the whole vector has length nxb*nyb+nl, will be output
! the vector x and Ax have length nxb*nyb+nl
!
  integer :: llen(nent-2), info, nl
  double precision :: time, dt, mydif
  double precision, dimension(1:nxb*nyb+nl) :: x, Ax
  double precision :: ux(1:nx*ny), uy(1:nx*(ny+1))
  type(iapt), dimension(:) :: iaary(nent-2)
!
  integer :: i, j, il, ierr, lenj, lenw, ic
  double precision, dimension(1:nxb*nyb) :: cg, cgv, pca
  double precision, dimension(1:nl) :: ca, ecg, qca
!
  lenj = nxb*nyb; lenw = lenj + nl
  cg = x(1:lenj); ca = x(lenj+1:lenw)
  !!call nmatv(cg, cgv, ux, uy, iaary, llen, time, info, dt, mydif)
  call nmatv(cg, cgv, ux, uy, iaary, llen, time, info, dt, mydif)
!!!=======================================================================
  call f2m(ca, pca, nl, ux, uy, iaary, llen, info, dt, mydif)! pca(1:nxb*nyb): output; ca(nl): input
  Ax(1:lenj) = cgv+pca
!!!=======================================================================
  call EM(cg, ecg, nl, iaary, llen, info, mydif)! cg(1:nxb*nyb); ecg(nl): output
  call Qmat(ca,qca, nl, iaary, llen, info, time,mydif)! qca(nl): output; ca(nl): input
! use the 4 above to get
  Ax(lenj+1:lenw) =ecg+qca
!!!=======================================================================
!!! start debugging mode
!!!=======================================================================
!!  !do i = 1, nl/2
!!  ic = 0
!!  do i = 1, llen(1)
!!    ic = ic + 1
!!    qca(ic) = ca(ic)*cmsi
!!    ic = ic + 1
!!    qca(ic) = ca(ic)*cpsi
!!  enddo
!!  ax(lenj+1:lenw) = qca
!!! finish debugging
!=======================================================================
!
  return
end subroutine amatv

!=======================================================================
!
subroutine cc_solver(rg,ra, cg, rb, nl, iaary, llen,info, isel,time,dt,mydif)
  implicit none
! solve the linear system for A*xc = rs, where xc is chemical variable only
! vmat is for velocity vector (it should be umat & vmat)
  !double precision, dimension(nxb*nyb) :: xc, rs, vmat
  double precision :: time, dt, mydif
  double precision, dimension(-1:nxb+2,-1:nyb+2) :: rg, cg
  double precision, dimension(1:nl) :: ra, rb
  integer :: nl, info, llen(nent-2), isel
  type(iapt), dimension(:) :: iaary(nent-2), itary(nent-2)
!
  integer, parameter :: mm=100
  integer :: ierr
  integer :: i, j, k, l, il, iter, it, mit, infon, mi, j1, lenj, lenw
  double precision :: con, cox, coy, tmp, normb, dr, omega, beta, tp2, eps, resid
  double precision :: ser0, ser1, serr, dr1, hit, st, error, tol, dx, dy, dthe
  double precision, dimension(nxb*nyb+nl) :: s, vk, x11, bf, r, yy, w, xx
  double precision, dimension(nxb*nyb+nl) :: ph, vv, rt, sv, pv, tv

  double precision :: v(nxb*nyb+nl,mm),  hj(nxb*nyb+nl,2), hq(nxb*nyb+nl,mm)
!========================================================================
  double precision :: rgg(nxb*nyb)
  double precision :: utx, vtx, ctop, cbot
  double precision :: ux(1:nx*ny), uy(1:nx*(ny+1))
!========================================================================
  !!double precision :: prem(nxb*nyb+nl)
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
 !!dt = dlt; ! now dt is an input
  dx = hg; dy = hg
  mit = nxb*nyb
  l   = mit; lenj = nxb*nyb; lenw = lenj + nl

  !!allocate(rhsv(nxb*nyb+nl), xc(nxb*nyb+nl), stat=ierr)
  allocate(rhsv(lenw), xcv(lenw), prem(lenw), stat=ierr)
  if (ierr.ne. 0) then
    print *, 'error in allocating variable in solver! stop'
    stop
  endif

  !!deb print *, 'Are we here?'

  v = 0.; hj = 0.; hq = 0.; serr = 0.
!========================================================================
! set RHS for Dirichlet BC 
!========================================================================
  !!cox = dif(1)*dlt/dx/dx; coy = dif(1)*dlt/dy/dy ! backward Euler
  con = mydif*dt/dx/dy; cox = con; coy = con 
  utx = 0.5d0*dt/hg; vtx = utx
!
  ux(1:nx*ny) = reshape(unc(0:nx-1,  1:ny  ), (/ nx*ny /)) 
  uy(1:nx*(ny+1)) = reshape(vnc(1:nx,  0:ny  ), (/ nx*(ny+1) /))
  rgg = 0.
!========================================================================
! modify top/bottom & 4 corners for Dirichlet BCs
!========================================================================
  ctop = cmsi; cbot = cmsi;
  rgg(cbt)=-2.*cbot*(-vtx*uy(cbt)-coy)
  rgg(cbl)=-2.*cbot*(-vtx*uy(cbl)-coy)
  rgg(cbr)=-2.*cbot*(-vtx*uy(cbr)-coy)
  rgg(ctp)=-2.*ctop*( vtx*uy(ctp+nx)-coy) ! add nx to reach the top for 
  rgg(ctl)=-2.*ctop*( vtx*uy(ctl+nx)-coy)
  rgg(ctr)=-2.*ctop*( vtx*uy(ctr+nx)-coy)
!========================================================================
  rhsv(1:lenj) = reshape(rg(1:nxb,1:nyb), (/ nxb*nyb /))
  rhsv(1:lenj) = rhsv(1:lenj) + rgg(1:lenj)
  rhsv(lenj+1:lenw) = ra(1:nl)
!
  normb=sqrt(dot_product(rhsv,rhsv))
!
  !!eps = 1.d-11*normb/2.**(2.*dble(l2ny-5)); resid=eps; tol = eps;
  eps = 1.d-9*normb; resid=eps; tol = eps;
  if (abs(normb)< eps) then
    xx=0.d0
    iter = 1
    resid = 0.
    cg = rg
    rb = ra
    print *, 'NO GMRES, iter', iter, normb
    return
  endif
!=======================================================================
! construct Jacobian predconditioner for GMRS
!=======================================================================
  tp2 = 1./(1.+2.*(cox+coy))
!
  prem = 1.
  prem(1:lenj) = tp2
  it = 0
  do il = 1, nib
    cp_lp => iaary(il)%p
    do i = 1, llen(il)
      it = it + 1
      tmp= cc_pump(cp_lp%sl,time)/cp_lp%dxs !pump is concentration dependent
      tmp= 0.; ! constant pump 
      call sub2ind(nxb, cp_lp%im(1), cp_lp%im(2), k)
      tp2 =-kc(il) + mydif*(cp_lp%c_off_grd_m(1) + cp_lp%c_ali_grd_m(1)) - cp_lp%vdn + min(tmp,0.)
      !!velo tp2 =-kc(il)+mydif*(cp_lp%c_off_grd_m(1)+cp_lp%c_ali_grd_m(1))+cp_lp%udn+ cp_lp%ocm + min(tmp,0.)
      ! last term above should have same sign as in Qmat
      prem(it+lenj) = 1./tp2

      it = it + 1
      call sub2ind(nxb, cp_lp%ip(1), cp_lp%ip(2), k)
      !!constpmp tp2 = kc(il) + cp_lp%c_off_grd_p(1) + cp_lp%c_ali_grd_p(1) - cp_lp%vdn 
      tp2 = kc(il) + mydif*(cp_lp%c_off_grd_p(1) + cp_lp%c_ali_grd_p(1)) - cp_lp%vdn + max(tmp,0.)
      !!velo tp2 = kc(il)+mydif*(cp_lp%c_off_grd_p(1)+cp_lp%c_ali_grd_p(1))+cp_lp%udn- cp_lp%ocp + max(tmp,0.)
      ! last term above should have same sign as in Qmat
      prem(it+lenj) = 1./tp2

      cp_lp => cp_lp%next
    enddo
  enddo
!
!=======================================================================
! GMRES
!=======================================================================
  !xc = 1.+0.001*tp2
  xcv = .0
  xcv(1:lenj) = reshape(cg(1:nxb,1:nyb), (/ nxb*nyb /))
  xcv(lenj+1:lenw) = rb



  iter =0
  v = 0.; hq = 0.; hj = 0.; serr = 0.
  do    j = 1, mit
    ser0 = 0.
    !!call matv(xcv, yy, iblst, fqlst,ilq, fplst, ilp, cox, coy, isel)
    !!deb xcv(:) = 1.
    call amatv(xcv, yy, nl, ux, uy, iaary, llen, time, info, dt, mydif)
    !!deb Y(1:nxb,1:nyb)=(reshape(yy, (/ nxb, nyb /)))
    !!deb print *, 'debugging'

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

        call amatv(vk,w, nl, ux, uy, iaary, llen, time, info, dt, mydif)
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
      call amatv(x11,w, nl, ux, uy, iaary, llen, time, info, dt, mydif)
      r = rhsv - w
      r = r*prem !PRECONDITIONER
      s(mi+1) = sqrt(dot_product(r,r))
      xcv = x11
      if( abs(s(mi+1)) .lt. tol .or. serr .lt. tol ) then
        ! this should be the return point 
        print '(1x,"IN  chem GMRES, iter", i5,1x, e14.6,5(i5,1x))', iter, serr, nl, llen
        cg(1:nxb,1:nyb) = (reshape(xcv(1:nxb*nyb), (/ nxb, nyb /)))
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
end subroutine cc_solver
!
!
!!!rmrdtdouble precision function mydot(x,y)
!!!rmrdt  double precision, dimension(-1:nxb+2,-1:nyb+2) :: x, y
!!!rmrdt!
!!!rmrdt  integer :: i, j
!!!rmrdt  double precision :: tsum 
!!!rmrdt  tsum = 0.d0
!!!rmrdt
!!!rmrdt  do i = 1, nxb
!!!rmrdt  do j = 1, nyb
!!!rmrdt    tsum = tsum + x(i,j)*y(i,j)
!!!rmrdt  enddo
!!!rmrdt  enddo
!!!rmrdt  mydot = tsum
!!!rmrdt! 
!!!rmrdt  return
!!!rmrdtend function mydot
! 
double precision function f(ax,bx,cx,ay,by,cy,xi,yj,s)
  implicit none
  double precision :: ax,bx,cx,ay,by,cy,xi,yj,s
  
!!old  f = 2.d0*(ax*ax+ay*ay)*s*s*s + 3.d0*(ax*bx+ay*by)*s*s + (2.d0*ax*cx+  &
!!old      2.d0*ay*cy+bx*bx+by*by-2.d0*(ax*xi+ay*yj))*s + (bx*cx+by*cy-bx*xi &
!!old      -by*yj)
  f = ((2.d0*(ax*ax+ay*ay)*s + 3.d0*(ax*bx+ay*by))*s + &
    (2.d0*ax*cx+2.d0*ay*cy+bx*bx+by*by-2.d0*(ax*xi+ay*yj)))*s &
    + (bx*cx+by*cy-bx*xi -by*yj)
  return 
end function f

double precision function fd(ax,bx,cx,ay,by,cy,xi,yj,s)
  implicit none
  double precision :: ax,bx,cx,ay,by,cy,xi,yj,s
  
!!old  fd = 6.d0*(ax*ax+ay*ay)*s*s + 6.d0*(ax*bx+ay*by)*s + (2.d0*ax*cx+     &
!!old      2.d0*ay*cy+bx*bx+by*by-2.d0*(ax*xi+ay*yj))
  fd = (6.d0*(ax*ax+ay*ay)*s + 6.d0*(ax*bx+ay*by))*s + (2.d0*ax*cx+     &
      2.d0*ay*cy+bx*bx+by*by-2.d0*(ax*xi+ay*yj))
  return 
end function fd

double precision function prdin(x1,y1,x2,y2)
  double precision :: x1,y1,x2,y2

  prdin = x1*x2+y1*y2

  return
end function prdin

!---+-B--1----+----2----+----3----+----4----+----5----+----6----+----7-E
!========================================================================
subroutine cc_vel(x, y, vx, vy, t)
  double precision :: x, y, vx, vy, t, r, rt, cs, si
  double precision :: x0, y0, xn, yn, tht


!!  x0 = .5; y0 = 0.5;
!!  xn = x - x0; yn = y - y0
!!  r = sqrt(xn*xn + yn*yn)
!!  rt= r*exp(-cpi*t)
!!  tht = atan2(y-y0,x-x0);
!!
!!  vx = rt*cos(tht)
!!  vy = rt*sin(tht)
  
  vx = vco*(.25d0-(y-half)*(y-half))
!  vx = 1.25
  vy = zero
!
  return 
end subroutine cc_vel

!!double precision function cc_pump(s,iside, isel)
double precision function cc_pump(s,time)
implicit none
  double precision :: s, time
  integer :: iside, isel
  ! we have to make sure s \in (0,2pi)
  double precision :: shl, shw, shs ! location, width, & strength leading gaussian
  double precision :: stl, stw, sts ! location, width, & strength back gaussian
  double precision :: sul, suw, sus ! location, width, & strength top/side gaussian
  double precision :: sbl, sbw, sbs ! location, width, & strength bot/side gaussian
!
  double precision :: tmp, su, sl, s0, tupi
  double precision :: pah, paw, pwt
!
  tupi = two*cpi
  tmp = 0.

  pah = 1.54308
  paw = 2.349975910046
  pwt = 1.655

  s0 = mod(mod(s, tupi)+tupi, tupi) ! archlength we will used

! if both shs and sts are -, we have shrinking
!  shs=-.10;      sts= 1.0*shs
! if both shs and sts are +, we have swelling
!  shs= .10;      sts= 1.0*shs
! if  shs and sts have opposite -/+, we move the cell
  shs= .01000;   sts=-1.0*shs
  sus=-shs*0.0;  sbs= 1.0*sus
  shl=cpi*.0 ;   stl=cpi*1.0
  sul=cpi*.80;   sbl=cpi*1.20
  shw=cpi*.21;   stw=shw*1.0;
  suw=shw*.5 ;   sbw=suw*1.

!!o  tmp=shs*exp(-.5*(s-shl)**2./shw**2.)+0. *exp(-.5*(s-shl-tupi)**2./shw**2.)+&
!!o      sts*exp(-.5*(s-stl)**2./stw**2.)
!working  tmp=shs*exp(-.5*(s-shl)**2./shw**2.)+shs*exp(-.5*(s-shl-tupi)**2./shw**2.)+&
!working      sts*exp(-.5*(s-stl)**2./stw**2.)
  tmp=shs*exp(-.5*(s-shl)**2./shw**2.)+shs*exp(-.5*(s-shl-tupi)**2./shw**2.)+&
      sts*exp(-.5*(s-stl)**2./stw**2.)+ &
      sus*exp(-.5*(s-sul)**2./suw**2.)+sbs*exp(-.5*(s-sbl)**2./sbw**2.)
!!  tmp=shs*(1.+cos(-(s-shl)))+sts*(1.+cos(-(s-stl)));
!!
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
  !!cc_pump = dble(isel)*tmp! don't use this anymore
  !!tmp = zero
  cc_pump =tmp
  !!cc_pump = 0.1
!
  return
end function cc_pump
!
double precision function cc_p(x, y, t)
  double precision :: x, y, t


!t  cc_p = sin(cpi*x)*sin(cpi*(y))*exp(-cpi*t)
!*  cc_p = (1.+sin(cpi*(x-0.5)))*(1.+sin(cpi*(y)))*exp(-cpi*t)
  cc_p = (1.+sin(2.*cpi*(x-0.25)))*(1.+sin(2.*cpi*(y-.25)))*exp(-cpi*t)
!o  cc_p = (1.+sin(2.*cpi*(x-0.25)))*(1.+sin(2.*cpi*(y-0.25)))*exp(-2.*cpi*t)
!!  cc_p = .8*exp(-(x*x+y*y)*.2/(t+1.))/cpi/(t+1.)
!
!=  cc_p = .5*pa*pa/cpi*exp(-(pa*(x-1.))**2.)*exp(-(pa*2.*(y-.5))**2.)*exp(-cpi*t)
  
  !!cc_p = 2.0
  cc_p = cc_p*fco
!
  return
end function cc_p

!

subroutine cc_cg(x, y, cx, cy, t)
  double precision :: x, y, t, cx, cy, c

  c = cc_p(x, y, t)

!t  cx = cpi*cos(cpi*x)*(sin(cpi*(y)))*exp(-cpi*t)
!t  cy = (sin(cpi*x))*cpi*cos(cpi*(y))*exp(-cpi*t)
!*  cx = cpi*(cos(cpi*(x-0.5)))*(1.+sin(cpi*(y)))*exp(-cpi*t)
!*  cy = (1.+sin(cpi*(x-0.5)))*cpi*(cos(cpi*(y)))*exp(-cpi*t)
  cx = cpi*2.*(cos(2.*cpi*(x-0.5)))*(1.+sin(2.*cpi*(y-.25)))*exp(-cpi*t)
  cy = (1.25+sin(2.*cpi*(x-0.5)))*2.*cpi*(cos(2.*cpi*(y-.25)))*exp(-cpi*t)
!=  cx = -2.*pa*pa*(x-1.0)*c
!=  cy = -4.*pa*pa*(2.*y-1.)*c
!o  cx = 2.*cpi*(cos(2.*cpi*(x-0.25)))*(1.+sin(2.*cpi*(y-0.25)))*exp(-2.*cpi*t)
!o  cy = (1.+sin(2.*cpi*(x-0.25)))*2.*cpi*(cos(2.*cpi*(y-0.25)))*exp(-2.*cpi*t)
!!  cx =-c*x*.4/(t+1.)
!!  cy =-c*y*.4/(t+1.)

  cx = cx*fco;
  cy = cy*fco;
!
  return
end subroutine cc_cg

double precision function cc_s(x, y, t)
  double precision :: x, y, c, t, vx, vy, ct, cx, cy, lapc

  call cc_vel(x, y, vx, vy, t)
!
  c  = cc_p(x, y, t);
  call cc_cg(x, y, cx, cy, t)
!*  ct = -cpi*c
!*  lapc = -cpi*cpi*(sin(cpi*(x-0.5))*(1.+sin(cpi*(y))) &
!*     &  +(1.+sin(cpi*(x-0.5)))*sin(cpi*(y)))*exp(-cpi*t)
  ct = -cpi*c
  lapc = -cpi*cpi*(4.*sin(2.*cpi*(x-0.5))*(1.+sin(2.*cpi*(y-.25))) &
     &  +(1.25+sin(2.*cpi*(x-0.5)))*4.*sin(2.*cpi*(y-.25)))*exp(-cpi*t)
  lapc = lapc*fco
!?  ct = -cpi*c
!?  lapc = -cpi*cpi*(sin(cpi*x)*(sin(cpi*(y))) &
!?     &  +(sin(cpi*x))*sin(cpi*(y)))*exp(-cpi*t)
!t  lapc = -cpi*cpi*(4.*sin(2.*cpi*(x-0.25))*(1.+sin(2.*cpi*(y-0.25))) &
!t     &  +4.*(1.+sin(2.*cpi*(x-0.25)))*sin(2.*cpi*(y-0.25)))*exp(-2.*cpi*t)
!!  ct = c*(x*x+y*y-5.*(t+1.))*.2/(t+1.)/(t+1.)
!!  lapc = 0.16*c/(t+1.)/(t+1)*(x*x+y*y-5.*(t+1.))
!!
!=  ct = -cpi*c
!=  lapc = pa*pa*((-2.+4.*pa*pa*(x-1.)**2.) + 8.*(2.*pa*pa*(2.*y-1.)**2.-1.))*c
!=
!
  vx = 0.; vy = 0.
  !!ct = 0.; 
  !!lapc = 0.;
  cc_s = ct + vx*cx+vy*cy - dif(1)*lapc
!
  return
end function cc_s

subroutine BdyCubic(iblst,mk,ccpc,nb,uin,isel)
  implicit none
  type(ibpt), dimension(:) :: iblst(nent-2)
  type(cc_bp), pointer :: curr, cb_bc
  integer :: isel, iflg
  double precision :: ccpc(-1:nxb+2,-1:nyb+2)  ! chemical concentration function
  double precision, dimension(:,:,:) :: mk(7,nring,nent-2), nb(nring,2, nent-2)
  double precision, dimension(:,:,:) :: uin(2,nring,nent-2)
!
  integer :: ima, jma, is, js, im, jm, cor(4,2), icor(4), ic(2)
  character :: equed
!
  double precision, dimension(nring+1) :: hs1, hs2, cox1, coy1, xp, yp
  double precision, dimension(nring  ) :: cv, dv, ev
  double precision, dimension(4*nring) :: work1, work2, work3
!
  integer i, j, k, il, info
  double precision :: tpx, tpy, tp1, tp2, tp3, tp4, sm, tmp
! 
  hs1 = 0.; sm = 0.
  do il = 1, nib
    curr => iblst(il)%p
    do i = 1, nring
      xp(i) = curr%x; yp(i) = curr%y
      tp1 = curr%prev%x; tp2 = curr%prev%y
      tp3 = curr%x;      tp4 = curr%y
      tmp = sqrt((tp1-tp3)*(tp1-tp3)+(tp2-tp4)*(tp2-tp4))
      hs1(i+1) = hs1(i) + tmp
      hs2(i+1) = hs1(i+1)-hs1(i)
      sm  = sm + tmp
      curr=> curr%prev
    enddo
  enddo
  hs1=hs1*two*cpi/sm
  hs2=hs2*two*cpi/sm
  xp(nring+1) = xp(1)
  yp(nring+1) = yp(1)
!  
  call splint(nring,xp,hs2,cox1)
  call splint(nring,yp,hs2,coy1)
!
!!deb  print *, 'test spline'
!!deb  
!!deb  tmp = cpi/dble(nring)
!!deb  do i = 1,2*nring
!!deb    sm = dble(i-0.5)*tmp
!!deb    call splval(nring,sm,xp,cox1,hs1,hs2,tpx,info)
!!deb    call splval(nring,sm,yp,coy1,hs1,hs2,tpy,info)
!!deb    print '(1x,1024(e16.7,1x))', tpx, tpy
!!deb  enddo
!!deb  print *, 'show original'
!!deb  do i = 1,nring+1
!!deb    print '(1x,1024(e16.7,1x))', xp(i),yp(i)
!!deb  enddo
!!deb  stop
!
  return
end subroutine BdyCubic
!
subroutine splint(n, x, hs, cx)
  integer :: n
  double precision, dimension(n+1) :: x, hs, cx
  double precision, dimension(1:n) :: c, d, e
  double precision, dimension(4*nring) :: work1, work2, work3
!
  integer :: i, j
  double precision :: tmp, tp1, tp2
!
  do i = 2, n-1
    tmp  = one/(hs(i)+hs(i+1))
    cx(i)= 6.0*tmp*((x(i+1)-x(i))/hs(i+1)-(x(i)-x(i-1))/hs(i))
    d(i) = two
    c(i) = hs(i)*tmp
    e(i) = hs(i+1)*tmp
  enddo
  d(n) = two
  tmp  = one/(hs(n)+hs(n+1))
  cx(n)= 6.0*tmp*((x(2)-x(1))/hs(2)-(x(1)-x(n))/hs(n+1))
  e(n) = hs(n+1)*tmp
  c(n) = hs(n)*tmp
  d(1) = two
  tmp  = one/(hs(2)+hs(n+1))
  cx(1)= 6.0*tmp*((x(2)-x(1))/hs(2)-(x(1)-x(n))/hs(n+1))
  e(1) = hs(2)*tmp
  c(1) = hs(n+1)*tmp
!
  call dgltsl(n,c,d,e,cx,work1,work2,work3)
!
  cx(n+1) = cx(1)
!
  return
end subroutine splint
!
!!subroutine splval(n,s,x,y,cox1,coy1,hs1,hs2,x1,y1,info)
subroutine splval(n,s,x,cx,hs1,hs2,x1,isel,info)
  integer :: n, info, isel
  double precision :: s, x1
  double precision, dimension(n+1) :: x, cx, hs1, hs2
!
  integer :: i
  double precision :: xt, a1, a2, a3, t
!
  do i = 1, n
    if (s .ge. hs1(i) .and. s .le. hs1(i+1)) then
      xt = x(i+1)
      a1 = (xt-x(i))/hs2(i+1)-(two*cx(i)+cx(i+1))*hs2(i+1)/6.
      a2 = cx(i)/two
      a3 = (cx(i+1)-cx(i))/(hs2(i+1)*6.)
      t  = s - hs1(i)
      x1 = x(i) + t*(a1 + t*(a2 + a3*t))
    endif
  enddo
!origSlow  do i = 1, n
!origSlow    if (s .ge. hs1(i) .and. s .le. hs1(i+1)) then
!origSlow      xt = x(i+1)
!origSlow      a1 = (xt-x(i))/hs2(i+1)-(two*cx(i)+cx(i+1))*hs2(i+1)/6.
!origSlow      a2 = cx(i)/two
!origSlow      a3 = (cx(i+1)-cx(i))/(hs2(i+1)*6.)
!origSlow      t  = s - hs1(i)
!origSlow      x1 = x(i) + t*(a1 + t*(a2 + a3*t))
!origSlow    endif
!origSlow  enddo
!
  return
end subroutine splval
!
subroutine dgltsl(n,c,d,e,b,work1,work2,work3)
!*******************************************************************************
!
!  dgltsl DECOMPOSE and solve the almost tridiagonal linear system with
!  the extra elements on (1,n) and (n,1). This routine calls cmlib dgtsl.
!
!   n   on entry, the dimension of the system.
!   d   double d(n) the diagonal  elements. On output d is destroyed.
!   c   double c(2),...,c(n) are the subdiagonal elements and c(1) is
!       the elements at (1,n). On output c is destroyed.
!   e   double e(1),e(2),..., e(n-1) are the superdiagonal elements and
!       e(n) is the elements at (n,1). On output e is destroyed.
!   b   double the right hand vector. On return it is the solution.
!   work1, work2, work3    double working spaces.
!
!*******************************************************************************

!!        implicit double precision (a-h,o-z)
  integer :: n
  double precision, dimension(n) :: c, d, e, b, work1, work2, work3
  double precision :: alf
  integer :: info, i

  work1(1:n) = c(1:n)
  work2(1:n) = d(1:n)
  work3(1:n) = e(1:n)

!  First we solve the n-1 system:
!  c(i)*y(i-1) + d(i)*y(i) + e(i)*y(i+1) = b(i), i = 1,2, ..., n-1
!  with y(0), y(n) = 0.

  !call dgtsl(n-1,c,d,e,b,info)
  info = 0
  call DGTSV(n-1, 1, c, D, e, B, n-1, INFO )
  if (info .ne. 0) then
    print *, 'Error in finding cubic spline, part0', info
    stop
  endif

!  Next we solve the n-1 system:
!  c(i)*z(i-1) + d(i)*z(i) + e(i)*z(i+1) = 0, i = 1,2, ..., n-1
!  with z(0)=y(n) = 1.

  e(1:n-1) = 0.d0
  e(1)   = -work1(1)
  e(n-1) = -work3(n-1)

  !call dgtsl(n-1,work1,work2,work3,e,info)
  info = 0
  call DGTSV(n-1, 1, work1, work2, work3, e, n-1, INFO )
  if (info .ne. 0) then
    print *, 'Error in finding cubic spline, part1', info
    stop
  endif

!  Then the solution x can be expressed
!  x = y + alf z.  alf is determined by the last equation.

  alf = b(n) - e(n)*b(1) - c(n)*b(n-1)
  b(n) = alf/(work2(n)+c(n)*e(n-1)+e(n)*e(1) )

  !!do i=1,n-1
  !!  b(i) = b(i) + b(n)*e(i)
  !!enddo
  b(1:n-1) = b(1:n-1) + b(n)*e(1:n-1)
!
  return
end subroutine dgltsl

!

subroutine ind2sub(m, i, j, l) 
implicit none 
integer :: m, n, i, j, l 
! 
 
i=l-((l-1)/m)*m 
j=(l-1)/m+1 
 
return 
end subroutine ind2sub 
!
subroutine cc_setRHS(rg, cp_lp, llen, il, isf, idn, mkk, isel, iaary, info, dt)
  implicit none
  double precision, intent(in) :: dt
  type(iapt), dimension(:) :: iaary(nent-2)
  integer :: isf(-1:nxb+2,-1:nyb+2)  ! chemical concentration function
  integer :: idn(-1:nxb+2,-1:nyb+2)  ! chemical concentration function
  double precision :: rg(-1:nxb+2,-1:nyb+2)  ! chemical concentration function
  type(cc_augvar), pointer :: cp_lp
  double precision, dimension(:,:,:) :: mkk(7,nring,nent-2)
  integer :: isel, info, llen(nent-2), il
!
  integer :: i1(2), k, it, jt, i_is(2), i_js(2)
  double precision :: tp(2), x, y, xm, ym, s0, s1, vx, vy, dx, dy, cm
  double precision :: st_is(2), st_js(2), i_st, j_st
  double precision :: tp1, tp2, tp3, tp4, tmp
!
  dx = hg; dy =hg;

  if (isel .eq. -1) then
    i1 = cp_lp%im; cm = cp_lp%cm
  else if (isel .eq. 1) then
    i1 = cp_lp%ip; cm = cp_lp%cp
  endif
  if (abs(idn(i1(1),i1(2))) .eq. 2.and. isf(i1(1),i1(2)).eq.0) then ! freshly cleared pt
  !!if ((idn(i1(1),i1(2))) .eq. 2*isel.and. isf(i1(1),i1(2)).eq.0) then
    !!rg(i1(1),i1(2)) = cm
    x = xmin+(dble(i1(1))-xshift)*dx; y = ymin+(dble(i1(2))-yshift)*dy; 
    !!print '(1x,"freshly cleared: side", 4(i3,1x),1024(e16.8,1x))', isel, i1, idf(i1(1),i1(2)), x, y, cm
    print '(1x,"freshly cleared(chem): side", 5(i3,1x),1024(e16.8,1x))', isel, i1,&
      idf(i1(1),i1(2)), idn(i1(1),i1(2))
    call FindBdyPtForFP(s0,s1,il,x,y,xm,ym,mkk,k,1) ! (xm,ym) along old boundary
    if (isel .eq. -1) then
      tmp = paraval(s1, nring, mkp, il, 0, 0, info)
    else
      tmp = paraval(s1, nring, mkp, il, 1, 0, info)
    endif
    rg(i1(1),i1(2)) = tmp
!
    vx = (xm-x)/dt; vy = (ym-y)/dt
    it = dsign(1.d0,dble(vx)); jt = dsign(1.d0,dble(vy))
    i_is(1) = i1(1) + it; i_is(2) = i1(2)
    i_js(1) = i1(1)     ; i_js(2) = i1(2) + jt
    !if (idf(cp_lp%m_is(1),cp_lp%m_is(2)).ne. idf(i1(1),i1(2)) .or. idf(cp_lp%m_js(1),cp_lp%m_js(2)).ne. idf(i1(1),i1(2))) then
    if (idf(i_is(1),i_is(2)).eq. idf(i1(1),i1(2))) then
      i_st = 1.
    else
      i_st = 0.
    endif
    if (idf(i_js(1),i_js(2)).eq. idf(i1(1),i1(2))) then
      j_st = 1.
    else
      j_st = 0.
    endif
    !if (idf(i_is(1),i_is(2)).ne. idf(i1(1),i1(2)) .or. idf(i_js(1),i_js(2)).ne. idf(i1(1),i1(2))) then
    if (i_st+j_st <0.5 ) then
      print *,  i_st, j_st, it, jt
      print *, 'there has to at least one point!'
      print '(1x,1024(i3,1x))', i_is, idf(i_is(1),i_is(2)), i_js, idf(i_js(1),i_js(2)), i1,idf(i1(1),i1(2)),  isel, idn(i1(1),i1(2))
      print '(1x,1024(i3,1x))', idf(i_is(1),i_is(2)),idf(i1(1),i1(2)), idf(i_js(1),i_js(2))
      print '(1x,1024(i3,1x))', chker0(i_is(1),i_is(2)),chker0(i1(1),i1(2)), chker0(i_js(1),i_js(2))
      print '(1x,1024(i3,1x))', chker1(i_is(1),i_is(2)),chker1(i1(1),i1(2)), chker1(i_js(1),i_js(2))
      print '(1x,1024(i3,1x))', k, it, jt, i1, chker0(i1(1),i1(2)), chker1(i1(1),i1(2))!!, chker2(i1(1),i1(2))
      print '(1x,1024(e16.8,1x))', s0, s1, cm, dist(x,y,xm,ym)/hg, dist(x,y,xm,ym)
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
      cm = cpi*2./dble(nring)
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
!!        tp(1) =  paraval(mkt(1,it,1), nring, mkt, 1, 0, 0, info)
!!        tp(2) =  paraval(mkt(1,it,1), nring, mkt, 1, 1, 0, info)
!!        print '(1x,1024(e16.8,1x))', tp
!!      enddo

      stop
    endif
    tp(1) = vx*dt/hg; tp(2) =-vx*dt/hg
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
    tp(1) = vy*dt/hg; tp(2) =-vy*dt/hg
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
    isf(i1(1),i1(2)) = isel
  endif
!
  return
end subroutine cc_setRHS
!
!
double precision function myhv(x)
  double precision :: x

  if (x >0) then
    myhv = 1.
  else
    myhv = 0.
  endif
!
  return
end function myhv
!
subroutine vc_cfs(vc_xpt,vc_ypt, fe)
implicit none
  double precision :: fe(nring)
  double precision :: vc_xpt(mpts), vc_ypt(mpts)
  double precision :: vc_cf(mcoor), vc_f(mcoor)
!
  integer :: ncp, nn, n, nlp, il, nil, njl
  double precision :: vc_nm((1)*nring*2), vc_mk(7,nring), tmp
!
  !!vc_nm = vc_nn
!!  do nn = 2*nfil+1, npts
!!    il  = (n-2*nfil-1)/nring + 1
!!    nlp = (n-2*nfil) - (n-2*nfil-1)/nring*nring
!!    ncp = n - 2*nfil
!!    nil = nlp + (il-1)*2*nring
!!    njl = nil + nring
!!    vc_nm(nil) = ndv(nring-nlp+1,1,il)
!!    vc_nm(njl) = ndv(nring-nlp+1,2,il)
!!  enddo
  !!vc_xpt = 0.; vc_ypt = 0.
  !!call vc_xr(vc_xx,vc_xpt,vc_ypt)
  vc_f = 0.; vc_cf = 0.
  call calcplateletforce(vc_xpt,vc_ypt,vc_f,vc_cf,gp,1)
  do nn = 2*nfil+1, npts
    il  = (nn-2*nfil-1)/nring + 1
    nlp = (nn-2*nfil) - (nn-2*nfil-1)/nring*nring 
    ncp = nn - 2*nfil
    nil = nlp + (il-1)*2*nring
    njl = nil + nring
    tmp = 0.
!!    tmp =  tmp + vc_cf(nil)*vc_nm(nil)
!!    tmp =  tmp + vc_cf(njl)*vc_nm(njl)
    tmp =  tmp + vc_cf(2*nn-1)*ndv(nring-nlp+1,1,il)
    tmp =  tmp + vc_cf(2*nn  )*ndv(nring-nlp+1,2,il)
    fe(nil)= tmp
  enddo
!
  return
end subroutine vc_cfs
!
subroutine cc_dxdt(x, y, vx, vy, t)
  double precision :: x, y, vx, vy, t, r, rt, cs, si
  double precision :: x0, y0, xn, yn, tht


!!  x0 = .5; y0 = 0.5;
!!  xn = x - x0; yn = y - y0
!!  r = sqrt(xn*xn + yn*yn)
!!  rt= r*exp(-cpi*t)
!!  tht = atan2(y-y0,x-x0);
!!
!!  vx = rt*cos(tht)
!!  vy = rt*sin(tht)
  
!!  vx =-vco*0.05
!!  vx =-0.05* vco*(.25-(y-0.5)**2.)
!  vx = 1.25
!!  vx = 0.
  call cc_vel(x,y,vx,vy,t)
  !!vx = vco*(.25-(y-0.5)**2.)
  !!vy = 0.d0

  !!vx = vx-0.5*vco*cos(t)
  !vx = 0.d0
  !!vy = 0.d0
!
  return 
end subroutine cc_dxdt
!---+-B--1----+----2----+----3----+----4----+----5----+----6----+----7-E
END MODULE chemical_mod
