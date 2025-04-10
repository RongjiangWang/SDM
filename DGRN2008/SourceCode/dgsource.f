      subroutine dgsource(dislocation)
      implicit none
c
      double precision dislocation
c
      include 'dgglob.h'
c
      integer i,n,istp
      double precision dam
c
      double precision pi
      data pi/3.14159265358979d0/
c
      do istp=1,4
        do i=1,8
          sfct0(i,istp)=0.d0
          sfct1(i,istp)=0.d0
        enddo
      enddo
c
      dam=dislocation/(2.d0*pi)
      n=nno(ls)
c
c     explosion (m11=m22=m33=1*kappa)
c
      ms(1)=0
      ics(1)=1
      sfct0(1,1)=-dam*(la(n)+2.d0*mu(n)/3.d0)/(la(n)+2.d0*mu(n))
      sfct1(4,1)=2.d0*mu(n)*sfct0(1,1)
c
c     strike-slip (m12=m21=1*mue)
c
      ms(2)=2
      ics(2)=-1
      sfct1(4,2)=dam*mu(n)
      sfct1(6,2)=-sfct1(4,2)
c
c     dip-slip (m13=m31=1*mue)
c
      ms(3)=1
      ics(3)=1
      sfct0(3,3)=-dam
      sfct0(5,3)=sfct0(3,3)
c
c     compensated linear vector dipole (CLVD) (m11=m22=-1*mue/2, m33=1*mue)
c
      ms(4)=0
      ics(4)=1
      sfct0(1,4)=-dam*mu(n)/(la(n)+2.d0*mu(n))
      sfct1(4,4)=dam*mu(n)*(3.d0-4.d0*mu(n)/(la(n)+2.d0*mu(n)))/2.d0
c
      return
      end
