	subroutine dgsublay(ierr)
	implicit none
c
	integer ierr
c
	include 'dgglob.h'
c
c	work space
c
	integer i,i0,l
	double precision dh,dla,dmu,drho,z,dz
c
	n0=0
c
	do l=1,l0-1
	  dz=z2(l)-z1(l)
	  dla=2.d0*dabs(la2(l)-la1(l))/(la2(l)+la1(l))
	  dmu=2.d0*dabs(mu2(l)-mu1(l))/(mu2(l)+mu1(l))
          if(rho2(l)+rho1(l).gt.0.d0)then
	    drho=2.d0*dabs(rho2(l)-rho1(l))/(rho2(l)+rho1(l))
          else
            drho=0.d0
          endif
	  i0=idnint(dmax1(1.d0,dla/reslm,dmu/reslm,drho/resld))
	  dla=(la2(l)-la1(l))/dz
	  dmu=(mu2(l)-mu1(l))/dz
	  drho=(rho2(l)-rho1(l))/dz
	  dh=dz/dble(i0)
	  do i=1,i0
	    n0=n0+1
	    if(n0.ge.lmax)then
	      ierr=1
	      return
	    endif
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
