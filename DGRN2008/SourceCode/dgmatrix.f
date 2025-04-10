      subroutine dgmatrix(a,k,z,n)
      implicit none
c
      include 'dgglob.h'
c
      integer n
      double precision k,z
      double precision a(6,6)
c
      double precision cx,cp,cm
      double precision cet,cgar
c
      cx=k*z
      cgar=gamma*rho(n)
      cet=la(n)+mu(n)
      cp=1.d0+cx
      cm=1.d0-cx
c
      a(1,1)=1.d0
      a(2,1)=2.d0*mu(n)*k
      a(3,1)=1.d0
      a(4,1)=2.d0*mu(n)*k
      a(5,1)=0.d0
      a(6,1)=-cgar
c
      a(1,2)=1.d0
      a(2,2)=-2.d0*mu(n)*k
      a(3,2)=-1.d0
      a(4,2)=2.d0*mu(n)*k
      a(5,2)=0.d0
      a(6,2)=-cgar
c
      a(1,3)=1.d0+cet*cm/mu(n)
      a(2,3)=2.d0*cet*cm*k
      a(3,3)=-1.d0-cet*cx/mu(n)
      a(4,3)=-2.d0*cet*cx*k
      a(5,3)=cgar*z
      a(6,3)=cgar*(cx-cet*cm/mu(n))
c
      a(1,4)=1.d0+cet*cp/mu(n)
      a(2,4)=-2.d0*cet*cp*k
      a(3,4)=1.d0-cet*cx/mu(n)
      a(4,4)=2.d0*cet*cx*k
      a(5,4)=cgar*z
      a(6,4)=cgar*(-cx-cet*cp/mu(n))
c
      a(1,5)=0.d0
      a(2,5)=0.d0
      a(3,5)=0.d0
      a(4,5)=0.d0
      a(5,5)=1.d0
      a(6,5)=k
c
      a(1,6)=0.d0
      a(2,6)=0.d0
      a(3,6)=0.d0
      a(4,6)=0.d0
      a(5,6)=1.d0
      a(6,6)=-k
c
      return
      end
