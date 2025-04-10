      subroutine hsmatinv(a,k,z,lahs,muhs)
      implicit none
c
      double precision k,z
      double precision lahs,muhs,a(6,6)
c
      double precision cx,cp,cm,xi,et
c
      cx=k*z
      xi=lahs+2.d0*muhs
      et=lahs+muhs
      cp=1.d0+cx
      cm=1.d0-cx
c
      a(1,1)=0.5d0*et*cx/xi
      a(1,2)=(muhs+et*cx)/(4.d0*k*muhs*xi)
      a(1,3)=0.5d0*et*cm/xi
      a(1,4)=0.25d0*(xi-et*cx)/(k*muhs*xi)
      a(1,5)=0.d0
      a(1,6)=0.d0
c
      a(2,1)=-0.5d0*et*cx/xi
      a(2,2)=0.25d0*(-muhs+et*cx)/(k*muhs*xi)
      a(2,3)=-et*cp/(2.d0*xi)
      a(2,4)=0.25d0*(xi+et*cx)/(k*muhs*xi)
      a(2,5)=0.d0
      a(2,6)=0.d0
c
      a(3,1)=0.5d0*muhs/xi
      a(3,2)=0.25d0/(k*xi)
      a(3,3)=-muhs/(2.d0*xi)
      a(3,4)=-0.25d0/(k*xi)
      a(3,5)=0.d0
      a(3,6)=0.d0
c
      a(4,1)=muhs/(2.d0*xi)
      a(4,2)=-0.25d0/(k*xi)
      a(4,3)=muhs/(2.d0*xi)
      a(4,4)=-0.25d0/(k*xi)
      a(4,5)=0.d0
      a(4,6)=0.d0
c
      a(5,1)=0.d0
      a(5,2)=0.d0
      a(5,3)=0.d0
      a(5,4)=0.d0
      a(5,5)=0.5d0
      a(5,6)=0.5d0/(k*muhs)
c
      a(6,1)=0.d0
      a(6,2)=0.d0
      a(6,3)=0.d0
      a(6,4)=0.d0
      a(6,5)=0.5d0
      a(6,6)=-0.5d0/(k*muhs)
c
      return
      end
