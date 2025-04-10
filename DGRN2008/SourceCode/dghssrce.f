      subroutine dghssrce(dislocation,lahs,muhs)
      implicit none
c
      double precision dislocation
      double precision lahs,muhs
c
      include 'dgglob.h'
c
      integer i,istp
      double precision dam
c
      double precision pi
      double precision d2,d3,d4
c
      data pi/3.14159265358979d0/
      data d2,d3,d4/2.d0,3.d0,4.d0/
c
      do istp=1,4
        do i=1,8
          sfcths0(i,istp)=0.d0
          sfcths1(i,istp)=0.d0
        enddo
      enddo
c
      dam=dislocation/(2.d0*pi)
c
c     explosion (m11=m22=m33=M0)
c
      sfcths0(1,1)=-dam*(lahs+d2*muhs/d3)/(lahs+d2*muhs)
      sfcths1(4,1)=d2*muhs*sfcths0(1,1)
c
c     strike-slip (m12=m21=M0)
c
      sfcths1(4,2)=dam*muhs
      sfcths1(6,2)=-sfcths1(4,2)
c
c     dip-slip (m13=m31=M0)
c
      sfcths0(3,3)=-dam
      sfcths0(5,3)=sfcths0(3,3)
c
c     compensated linear vector dipole (CLVD) (m11=m22=-M0/2, M33=M0)
c
      sfcths0(1,4)=-dam*muhs/(lahs+d2*muhs)
      sfcths1(4,4)=dam*muhs*(d3-d4*muhs/(lahs+d2*muhs))/d2
c
      return
      end
