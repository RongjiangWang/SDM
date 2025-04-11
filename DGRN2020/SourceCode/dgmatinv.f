      subroutine dgmatinv(a,k,z,n)
      use dgalloc
      implicit none
c
      integer*4 n
      real*8 k,z
      real*8 a(6,6)
c
      real*8 k2,cx,cp,cm,cxi,cet,cgr
c
      k2=k*k
      cgr=gamma*rho(n)
c
      cx=k*z
      cxi=la(n)+2.d0*mu(n)
      cet=la(n)+mu(n)
      cp=1.d0+cx
      cm=1.d0-cx
c
      a(1,1)=cet*cx/(2.d0*cxi)
      a(1,2)=(mu(n)+cet*cx)/(4.d0*k*mu(n)*cxi)
      a(1,3)=cet*cm/(2.d0*cxi)
      a(1,4)=(cxi-cet*cx)/(4.d0*k*mu(n)*cxi)
      a(1,5)=0.d0
      a(1,6)=0.d0
c
      a(2,1)=-cet*cx/(2.d0*cxi)
      a(2,2)=(-mu(n)+cet*cx)/(4.d0*k*mu(n)*cxi)
      a(2,3)=-cet*cp/(2.d0*cxi)
      a(2,4)=(cxi+cet*cx)/(4.d0*k*mu(n)*cxi)
      a(2,5)=0.d0
      a(2,6)=0.d0
c
      a(3,1)=mu(n)/(2.d0*cxi)
      a(3,2)=1.d0/(4.d0*k*cxi)
      a(3,3)=-mu(n)/(2.d0*cxi)
      a(3,4)=-1.d0/(4.d0*k*cxi)
      a(3,5)=0.d0
      a(3,6)=0.d0
c
      a(4,1)=mu(n)/(2.d0*cxi)
      a(4,2)=-1.d0/(4.d0*k*cxi)
      a(4,3)=mu(n)/(2.d0*cxi)
      a(4,4)=-1.d0/(4.d0*k*cxi)
      a(4,5)=0.d0
      a(4,6)=0.d0
c
      a(5,1)=cgr*(cet-mu(n)*cx)/(2.d0*k*cxi)
      a(5,2)=-cgr*cx/(4.d0*k2*cxi)
      a(5,3)=cgr*mu(n)*cx/(2.d0*k*cxi)
      a(5,4)=cgr*cp/(4.d0*k2*cxi)
      a(5,5)=1.d0/2.d0
      a(5,6)=1.d0/(2.d0*k)
c
      a(6,1)=-cgr*(cet+mu(n)*cx)/(2.d0*k*cxi)
      a(6,2)=cgr*cx/(4.d0*k2*cxi)
      a(6,3)=-cgr*mu(n)*cx/(2.d0*k*cxi)
      a(6,4)=-cgr*cm/(4.d0*k2*cxi)
      a(6,5)=1.d0/2.d0
      a(6,6)=-1.d0/(2.d0*k)
c
      return
      end
