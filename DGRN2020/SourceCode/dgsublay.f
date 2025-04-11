	subroutine dgsublay(ierr)
      use dgalloc
	implicit none
c
	integer*4 ierr
c
c	work space
c
	integer*4 i,l
	real*8 dh,dla,dmu,drho,z,dz
c
      integer*4, allocatable:: i0(:)
c
      allocate(i0(l0),stat=ierr)
      if(ierr.ne.0)stop 'Error in dgsublay: i0 not allocated!'
c
	n0=1
      do l=1,l0-1
	  dz=z2(l)-z1(l)
	  dla=2.d0*dabs(la2(l)-la1(l))/(la2(l)+la1(l))
	  dmu=2.d0*dabs(mu2(l)-mu1(l))/(mu2(l)+mu1(l))
        if(rho2(l)+rho1(l).gt.0.d0)then
	  drho=2.d0*dabs(rho2(l)-rho1(l))/(rho2(l)+rho1(l))
        else
          drho=0.d0
        endif
        i0(l)=idnint(dmax1(1.d0,dla/reslm,dmu/reslm,drho/resld))
        n0=n0+i0(l)
      enddo
c
      allocate(h(n0),stat=ierr)
      if(ierr.ne.0)stop 'Error in dgsublay: h not allocated!'
      allocate(la(n0),stat=ierr)
      if(ierr.ne.0)stop 'Error in dgsublay: la not allocated!'
      allocate(mu(n0),stat=ierr)
      if(ierr.ne.0)stop 'Error in dgsublay: mu not allocated!'
      allocate(rho(n0),stat=ierr)
      if(ierr.ne.0)stop 'Error in dgsublay: rho not allocated!'
c
      lp=n0+2
c
      allocate(xh(lp),stat=ierr)
      if(ierr.ne.0)stop 'Error in dgsublay: xh not allocated!'
      allocate(xh0(lp),stat=ierr)
      if(ierr.ne.0)stop 'Error in dgsublay: xh0 not allocated!'
c
      allocate(hp(lp),stat=ierr)
      if(ierr.ne.0)stop 'Error in dgsublay: hp not allocated!'
      allocate(nno(lp),stat=ierr)
      if(ierr.ne.0)stop 'Error in dgsublay: nno not allocated!'
c
      allocate(maup(6,6,lp),stat=ierr)
      if(ierr.ne.0)stop 'Error in dgsublay: maup not allocated!'
      allocate(maiup(6,6,lp),stat=ierr)
      if(ierr.ne.0)stop 'Error in dgsublay: maiup not allocated!'
      allocate(malw(6,6,lp),stat=ierr)
      if(ierr.ne.0)stop 'Error in dgsublay: malw not allocated!'
      allocate(mailw(6,6,lp),stat=ierr)
      if(ierr.ne.0)stop 'Error in dgsublay: mailw not allocated!'
      allocate(hkup(2,2,lp),stat=ierr)
      if(ierr.ne.0)stop 'Error in dgsublay: hkup not allocated!'
      allocate(hklw(2,2,lp),stat=ierr)
      if(ierr.ne.0)stop 'Error in dgsublay: hklw not allocated!'
c
      n0=0
c
	do l=1,l0-1
	  dz=z2(l)-z1(l)
	  dla=(la2(l)-la1(l))/dz
	  dmu=(mu2(l)-mu1(l))/dz
	  drho=(rho2(l)-rho1(l))/dz
	  dh=dz/dble(i0(l))
	  do i=1,i0(l)
	    n0=n0+1
	    h(n0)=dh
	    z=(dble(i)-0.5d0)*dh
	    la(n0)=la1(l)+dla*z
	    mu(n0)=mu1(l)+dmu*z
	    rho(n0)=rho1(l)+drho*z
	  enddo
	enddo
c
c	last layer is half-space
c
	n0=n0+1
	h(n0)=0.d0
	la(n0)=la1(l0)
	mu(n0)=mu1(l0)
	rho(n0)=rho1(l0)
c
	write(*,'(7a)')'  no',' thick(m)    ','  la(Pa)    ',
     &    '  mu(Pa)    ','rho(kg/m^3) '
	do i=1,n0
	  write(*,1001)i,h(i),la(i),mu(i),rho(i)
	enddo
1001	format(i4,f11.4,5E12.4)
	ierr=0
	return
	end
