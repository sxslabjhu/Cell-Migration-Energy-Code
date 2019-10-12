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

!Main program of calling initialization, and 1) FSI 2) Osmosis 3) Actin network
!
program TwoPhase2D
  use, intrinsic :: iso_c_binding
  use parameters
  use geometry
  use linsys
  use myfft
  use IBmod
  use IBforce
  use chemical_mod
  use network_mod
  use solver_mod
  use energy

  implicit none

!
!========================================================================
  double precision :: t0, tf, dt, dx, dy, h
  double precision :: Esol, Iflow, En, Inet, Imem, Jmem
!
  integer :: it, Nsteps
!
  double precision, dimension(-1:nx+1,-1:ny+1) :: mb1, ma1, &
    u1, u0, v1, v0, uold, vold, ua, va, pa, p0, p1
  double precision, dimension(-1:nx+2,-1:ny+2) :: thetu,thetv, sigx,sigy,&
          thetx,thety, grsgx,grsgy,grx, gry, sga
  double precision, dimension(-1:nx+1,-1:ny+1,2) :: ibfrc
  double precision, dimension(mcoor) :: fs, cfs
  double precision :: time, time0, time1, time2, time3, us(mpts,2), gp(mpts)
  integer :: nit, jt, ik, jk, i, j, k, isel, ins, iadv, iter, N, M, numvol
  double precision :: x, y, tp1, tp2, tp3, tp4, tmp, totVol
!========================================================================
  complex*16, dimension(-1:nx+1,-1:ny+1) :: rsu, rsv
  double precision :: alphab, alpha, beta, gam, totalvol, meanvol
  double precision, dimension(-1:nx+1,4) :: uvbc0, uvbc1, uvbca !ubt, utp, vbt, vtp
!========================================================================
  integer :: iflg, strlen
  character*40 :: ibfile
  complex*16, dimension(3*ny-1,nx) :: resul
  double precision, dimension(3*ny-1,nx) :: rsol
  double precision, dimension(2*nring) :: rhs, fold, f0, cf0, fdir, vxi, vxo
  double precision, dimension(nring) :: xn,yn, xc,yc
  double precision :: smvec(1:nring), spvec(nring)
!=======================================================================
!  new variables
!========================================================================
!!  integer, parameter :: nent = 3
  type(ibpt), dimension(:) :: ibary(nent-2)
  type(iapt), dimension(:) :: iaary(nent-2)
  type(cc_augvar), pointer :: cp_lp
!!  type(cc_bp), pointer :: curr
  integer :: info, llen(nent-2), im(2), ict, iloc(4,2), ipl, il
  double precision :: st_p(4), theta, s0, etap
  double precision :: uin(2,nring,nent-2)
!=======================================================================
!
  print '("Before reading", 2(i7,1x), 10(e14.6,1x))', ntmax, nfreq, dlt, CFL, nu, bk, vco
  call readpar(CFL,nu,bk,vco,ntmax,nfreq,dlt)
  print '("After reading ", 2(i7,1x), 10(e14.6,1x))', ntmax, nfreq, dlt, CFL, nu, bk, vco

  !!fname='real1.dat'
  N = nx; M = ny;
  dt = dlt;
  dx = hg; dy = hg;  h = hg
  etap=eta/(eta+etas)
!
  print '(1x,"Spacial-x increment: ", 1(e14.6,1x), " Grid nx: ", i5)', dx, nx
  print '(1x,"Spacial-y increment: ", 1(e14.6,1x), " Grid ny: ", i5)', dy, ny
  print '(1x,"viscosity: ", 1(e14.6,1x))', nu
!
  t0 = zero
  tf = .1d0
  !Nsteps = nint(tf/dt) ! # of time steps determined by dt and tfinal (tf)
  Nsteps = ntmax! # of time steps determined in input file
!
  print '(1x,"Temporal increment: ", 1(e14.6,1x))', dt
  print '(1x,"Total # of time steps: ", 1(i8,1x))', Nsteps
  print *,'============================================================'
  print *,'Starting....'
  print *,'============================================================'
!
!=======================================================================
  call cpu_time(time0) ! initialize time stamp
  time3 = time0
!=======================================================================
! initialize IB locations, and calculate elastic force
  call InitLocCell(xpt, ypt, fs, cfs, gp)
!
!!!========================================================================
!!! initialize geometry and set initial chemical and actin network 
!!!========================================================================
!!prt  print *, 'step -1'
  call getGeometry(ibary,iaary,llen,-1,uin, dt)
!!  do il = 1, nib
!!    cp_lp => iaary(il)%p
!!    do i = 1, llen(il)
!!      print *, cp_lp%s, cp_lp%im, 
!!      cp_lp => cp_lp%next
!!    enddo
!!  enddo
  call cc_init(iaary,llen)
  call nt_init(iaary,llen)

! spread the force to ibfrc
! note fs is force @ IB pts, and ibfrc is over Eulerian grids
  call IBspreadS(xpt,ypt,gp,fs,ibfrc) !<-F(X0)
!!!========================================================================
!
  strlen = len_trim(runname)
  write(ibfile,'(2a,i4.4)') runname(1:strlen),'.ib.',outcount
  open(67,file=ibfile,form='formatted',action='write')
  do i = 1, nring
    write(67,'(2(e22.14,1x))')xpt(i),ypt(i)
  enddo
  write(67,'(2(e22.14,1x))')xpt(1),ypt(1)
  close(67)
!
  ins = 0
  select case (ins)
  case (0)
    alpha=4.d0*dt*nu/h/h + bk*dt
    beta=dt*nu/h/h
    gam=dt/h
  case (1) !Euler
    alpha=1.d0+4.d0*dt*nu/h/h ! 
    beta=dt*nu/h/h
    gam=dt/h
    isel = 0 ! =1 CN ; =0 Euler
  case (2) !CN
    alpha=1.d0+2.d0*dt*nu/h/h ! 
    beta=0.5d0*dt*nu/h/h
    gam=0.5d0*dt/h
    !!gam=dt/h
    isel = 1 ! =1 CN ; =0 Euler
  case default
    print *, 'wrong isel for linear solve. stop'
    stop
  end select
!
! initialize velocity, pressure and source terms
!
  time = 0.
  iadv = 0 ! =1 w/ advection; =0 w/o advection
  call inituvp(u0, v0, p0, ma1, mb1, uvbc0, 0, dt, time)  ! initialize velocity, uvbc0 is assigned here
  !!call setbc(u0,uvbc0,1)
  !!call setbc(v0,uvbc0,2)
  !!call setbc(p0,uvbc0,3)
  u0 = 0.; v0 = 0.; p0 = 0.; ! set all to initial 0
  uold = u0; vold = v0;
  uvbc0 = zero
!=======================================================================
! save initial IB locations and u,v, p
  outcount = 0
  call outuvpx(nring, u0,v0,p0,xpt,ypt,outcount)
!=======================================================================
! Advance in time
!=======================================================================
  do it = 1, Nsteps
    time = dble(it)*dt
!!deb    call inituvp(ua, va, pa, ma1, mb1, uvbc1, isel, dt, time)  ! get RHS source term
!!deb    ua=0; va=0; pa = 0;
    oxpt = xpt; oypt = ypt ! save old IB locations
    fold = fs; !save old force
!!
    ma1 = ibfrc(:,:,1)
    mb1 = ibfrc(:,:,2)
    uvbc0 = 0.d0;
    uvbc1 = 0.d0;
!
!!deb    smvec=zero
!!deb    spvec=zero
!!deb    il = 1
!!deb    mkr(1,:,il) = mkn(1,:,il)
!!deb    cp_lp => iaary(il)%p
!!deb    do i = 1, llen(il)
!!deb      s0 = cp_lp%s
!!deb      tmp= siga(s0,cp_lp%xb,cp_lp%yb,1)
!!deb      spvec(i) = tmp
!!deb !      mkr(1,i,il) = s0
!!deb      cp_lp => cp_lp%next
!!deb    enddo
!!deb    smvec = spvec
!!deb    call smoothcrossing(olen, nib,il,mkr,iaary,llen, smvec,spvec, info)
!!deb    tp1 = paraval(half , nring, mkr, il, 0, 0, info) ! on + side old network volume
!!deb    tp2 = paraval(half , nring, mkr, il, 1, 0, info) ! on + side old network volume
!!deb    print '("TEST h", 1024(e14.6,1x))', tp1, tp2, s0, siga(s0,half,half,0), mkr(:,1,il)

    call sg_IBbdy(sgma, iaary, llen, time, dt, info)
!!prt    il = 1
!!prt    tp1 = paraval(half , nring, mkr, il, 0, 0, info) ! on + side old network volume
!!prt    tp2 = paraval(half , nring, mkr, il, 1, 0, info) ! on + side old network volume
!!prt    print '("WHY   ", 1024(e14.6,1x))', tp1, tp2, s0, siga(s0,half,half,0), mkr(:,1,il)
!!simplifed    call extend(iaary,llen,mkr,sgma)
!!!!deb!!deb    print *,'do we get it?'

!!deb    sgma = zero
    sga = (nt_pc)*dble(iid)
    totalvol = zero
    numvol = 0
    do j = 2, ny-1
    do i = 2, nx-1
      if (idf(i,j) > 0) then
      tmp = nt_theta(sga(i,j))
      totalvol = totalvol + tmp
      numvol = numvol + 1
      endif
    enddo
    enddo
    nit = sum(iid(1:nx,1:ny))
    meanvol= totalvol/dble(nit)
!!    meanvol= totalvol/dble(numvol)
    totalvol=totalvol*hg*hg
    call extend(iaary,llen,mkq,sga)
!!    sga = (sga)*iid
    !!if (i> 1) call extend(iaary,llen,mkq,sga)
!!    call extend(iaary,llen,mkq,sga)
    sigx = zero; thetx = zero; thetu = zero; grx = zero; grsgx = zero; 
    do i = 1, nx-1
      sigx(i,1:ny) = half*( sga(i,1:ny)+sga(i+1,1:ny) )
    enddo
    do j = 2, ny-2
    do i = 2, nx-2
      if ( jdu(i,j) > 0 ) then
        thetx(i,j) = nt_theta(sigx(i,j))
        grx(i,j) = sga(i+1,j) - sga(i,j)
      endif
    enddo
    enddo
    continue
    grsgx = eta*(grx/dx+sik)/(eta+etas)
    !!grsgx = grsgx/dx
    grx(0:nx,0:ny) = eta*(etas/(eta+etas))*thetx(0:nx,0:ny)*u0(0:nx,0:ny) +grsgx(0:nx,0:ny)
    !!thetu(0:nx,0:ny) = grsgx(0:nx,0:ny)
    thetu = grx*dble(jdu)
!
    sigy = zero; thety = zero; thetv = zero; gry = zero; grsgy = zero; 
    do j = 2, ny-1
      sigy(1:nx,j) = half*( sga(1:nx,j)+sga(1:nx,j+1) )
    enddo
    do j = 2, ny-2
    do i = 2, nx-2
      if ( jdv(i,j) > 0 ) then
        thety(i,j) = nt_theta(sigy(i,j))
        gry(i,j) = sga(i,j+1) - sga(i,j)
      endif
    enddo
    enddo
    continue
    grsgy = eta*gry/dy/(eta+etas)
    !!grsgy = grsgy/dy
    gry(0:nx,0:ny) = eta*(etas/(eta+etas))*thety(0:nx,0:ny)*v0(0:nx,0:ny) +grsgy(0:nx,0:ny)
    !!thetv(0:nx,0:ny) = grsgy(0:nx,0:ny)
    thetv = gry*dble(jdv)
!
    grx(-1:nx+1,-1:ny+1) = ma1(-1:nx+1,-1:ny+1);
    gry(-1:nx+1,-1:ny+1) = mb1(-1:nx+1,-1:ny+1);
!========================================================================
!========================================================================
!!needchange    alpha=4.d0*dt*nu/h/h + bk*dt
!!needchange    beta=dt*nu/h/h
!!needchange    gam=dt/h
!========================================================================
!========================================================================
    ma1(-1:nx+1,-1:ny+1) = grx(-1:nx+1,-1:ny+1) + thetu(-1:nx+1,-1:ny+1)
    mb1(-1:nx+1,-1:ny+1) = gry(-1:nx+1,-1:ny+1) + thetv(-1:nx+1,-1:ny+1)
!========================================================================
!!prt    print '(" rhs 4 uv ", 1024(e14.6,1x))', maxval(ma1), maxval(mb1), &
!!prt      maxval(thetu), maxval(thetv)
    call releasemem(ibary,iaary,llen)
!
    uold = u0; vold = v0;

!!deb    print '("before ", 1024(e14.6,1x))', maxval(u0), minval(u0),&
!!deb    maxval(v0), minval(v0), u0(88,63), v0(88,63)
!
    call getRHS(rsu, rsv, u0, v0, p0, ma1, mb1, uvbc0, uold, vold, ins, iadv, dt) 
!!deb    print '("results ", 1024(e14.6,1x))',maxval(real(rsu)),rsol(1,1),&
!!deb   rsol(100,20),maxval(real(rsv)), maxval(dimag(rsv)), minval(dimag(rsv))
    call sollinsys(rsu,rsv,resul,rsol, alpha, beta, gam, dt)
!!deb    print '("middle ", 1024(e14.6,1x))', maxval(rsol), minval(rsol), rsol(1,1), &
!!deb    rsol(100,nx/2), maxval(real(rsu)), maxval(real(rsv))
    call rsoltouvp(rsol,ua,va,pa)
!!deb    print '("after ", 1024(e14.6,1x))', maxval(ua), minval(ua),&
!!deb    maxval(va), minval(va), ua(88,63), va(88,63)
    !4) prepare for new step
    !!call calcplateletforce(xpt,ypt,f0,cf0)
    !!call IBspreadS(xpt,ypt,f0,ibfrc) !
    !
    u1 = ua; v1 = va; p1 = pa ! 
    fs = f0
    u0 = ua; v0 = va; p0 = pa
    call setbc(u0,uvbc1,-1)
    call setbc(v0,uvbc1,-2)
    call setbc(p0,uvbc1,3)
!
!!prt    print *, '(x,y)', xpt(1),ypt(1), maxval(u1),maxval(v1)
    call xmovemac(xpt,ypt,us,u1,v1,fs,cfs,dt)


!!prt    print '("indeed ", 1024(e14.6,1x))', maxval(u1), minval(u1),&
!!prt    maxval(v1), minval(v1), u1(88,63), v1(88,63)
    fs = 0.; cfs = 0.
!!prt    print *, 'done (x,y)', xpt(1),ypt(1), it
    call calcplateletforce(xpt,ypt,fs,cfs,gp,1)

!========================================================================
! need to update geometry here 
!========================================================================
!
!!prt    print *, 'step 0', it
    unc = u0; vnc = v0; 
    call getGeometry(ibary,iaary,llen,1,uin, dt)
    
!!prt    print *, 'step 1', it
!
    call cc_IBbdy(cc_pn, cc_pc, ibary, iaary, llen, uin, 1, time, dt, dif(1))
    cc_pc = cc_pn
    call getIflow(Iflow, nu, dif(1), unc, vnc, cc_pc, iaary, llen)
    call getGs(Esol, cc_pc(-1:nx+1,-1:ny+1),iflg)
!!    print *, "max val", maxval(log(cc_pc(1:nx,1:ny)))
!
!!prt    print *, 'step 2', it
    call nt_IBbdy(nt_pn, nt_pc, ibary, iaary, llen, mkn, 1, time, dt, dif(2))
    nt_pc = nt_pn
    call getInet(Inet, nt_pc, unc,vnc, iaary, llen)
    call getEn(En,nt_pc,iflg)
    call getImem(Imem, kc(1), kw(1), mkp, mkq, iaary, llen,iflg)
    call getJmem(Jmem,mkp,mkq,iaary,llen,iflg,time)
!
    ibfrc = 0.
    call IBspreadS(xpt,ypt,gp,fs,ibfrc)

!!prtdeb    print *, 'Lag variable '
!!prtdeb    do il = 1, nib
!!prtdeb      cp_lp => iaary(il)%p
!!prtdeb      do i = 1, llen(il)
!!prtdeb        tmp = cp_lp%s
!!prtdeb        tp1 = paraval(tmp, nring, mkq, il, 1, 0, info) ! on + side old network volume
!!prtdeb 
!!prtdeb        print '(1x,1024(e14.8,1x))', cp_lp%xb,  cp_lp%yb, tmp, tp1, cp_lp%cp, cp_lp%cm
!!prtdeb
!!prtdeb        cp_lp => cp_lp%next
!!prtdeb      enddo
!!prtdeb    enddo
!!prtdeb    print *, 'Ending variable '
!!
!!Did it before, do not releasing here   call releasemem(ibary,iaary,llen)
!
    if (it -nfreq*(it/nfreq) .eq. 0) then ! save data
      !!prt print '(1x, "Step ", i5, " Saving Data ....")', it
      outcount = outcount + 1
      call outuvpx(nring, u0,v0,p0,xpt,ypt,outcount)
!!notused      write(ibfile,'(2a,i4.4)') runname(1:strlen),'.sg.',outcount
!!notused      open(67,file=ibfile,access='stream',action='write')
!!notused      write(67)sga
!!notused      close(67)
      !
!!       do il = 1, nib
!!         cp_lp => iaary(il)%p
!!         do i = 1, llen(il)
!!           print '(1x," check data ", i5, 10(e14.8,1x))', i, cp_lp%sgm, cp_lp%gsgm, cp_lp%s,&
!!             cp_lp%nx, cp_lp%ny, cp_lp%xb, cp_lp%yb
!!           cp_lp => cp_lp%next
!!         enddo
!!       enddo
    endif
!
    do j = 1, ny
    do i = 1, nx
      if (idf(i,j) > 0) then
        grx(i,j) = nt_theta(nt_pn(i,j))
      endif
    enddo
    enddo
    tp1=sum(abs(grx(1:nx,1:ny))*dble(iid(1:nx,1:ny)))*hg
    tp2=sum(nt_pn(1:nx,1:ny)*dble(iid(1:nx,1:ny)))
!!    tp1=sum(grx*dble(iid))
!!    tp2=sum(nt_pn*dble(iid))
    !!tp3=maxval((u0(1:nx-1,1:ny-1)))
    !!tp4=maxval((v0(1:nx-1,1:ny-1)))
    tp3 = sum(xpt)/dble(nring); !!tp4 = maxval(xpt);
    !!tmp=maxval(abs(vc_fe))
    tmp=sum(ypt)/dble(nring)
    tp2=dot_product((xpt-minval(xpt))-(xbp-minval(xbp)),(xpt-minval(xpt))-(xbp-minval(xbp)))
    tp2=sqrt(tp2)
    tp4=dot_product((ypt-minval(ypt))-(ybp-minval(ybp)),(ypt-minval(ypt))-(ybp-minval(ybp)))
    tp4=sqrt(tp4)
    print '(1x,"Iter: ", 2(i10,1x), 20(e16.9,1x))', it, nit,time, tp1, tp3, &
            tmp, tp2+tp4,totalvol, Jmem-(Iflow+Inet+Imem), &
            Esol, En, Iflow, Inet, Imem, Jmem
!
  enddo ! <-main time evolution loop, it=+1

!=======================================================================
!
!=======================================================================
! check time spent
!
  call cpu_time(time2)
  !!print *, 'Time:', time
  write(*,200)it-1, time2-time3
!
!=======================================================================
200 format(1x,'<<-->>Cputime in ', i7,' step(s):'f16.5,'s')
!=======================================================================
!
end program TwoPhase2D
!
!=======================================================================
