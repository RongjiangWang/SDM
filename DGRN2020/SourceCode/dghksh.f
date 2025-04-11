      subroutine dghksh(hk,m,k,z,n)
      use dgalloc
      implicit none
c
      integer*4 m,n
      real*8 k,z
      real*8 hk(m,m)
c
	real*8 c2x,cem,cch,csh
c
	c2x=2.d0*k*z
c
	if(m.eq.2)then
c
c	  haskell propagator matrix for SH waves
c
	if(z.gt.0.d0)then
	  cem=dexp(-c2x)
	  cch=0.5d0*(1.d0+cem)
	  csh=0.5d0*(1.d0-cem)
	else
	  cem=dexp(c2x)
	  cch=0.5d0*(1.d0+cem)
	  csh=-0.5d0*(1.d0-cem)
	endif
c
c	propagator matrix for SH waves
c
	hk(1,1)=cch
	hk(1,2)=csh/(mu(n)*k)
	hk(2,1)=csh*mu(n)*k
	hk(2,2)=cch
	else
	  print *,'error in peghask: m schould be 2!'
	  return
	endif
c
	return
	end
