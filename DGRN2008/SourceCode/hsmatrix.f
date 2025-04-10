      subroutine hsmatrix(a,k,z,lahs,muhs)
      implicit none
c
      double precision k,z
      double precision lahs,muhs,a(6,6)
c
      double precision cx,cp,cm,et
c
      cx=k*z
      et=lahs+muhs
      cp=1.d0+cx
      cm=1.d0-cx
c
      a(1,1)=1.d0
      a(2,1)=2.d0*muhs*k
      a(3,1)=1.d0
      a(4,1)=2.d0*muhs*k
      a(5,1)=0.d0
      a(6,1)=0.d0
c
      a(1,2)=1.d0
      a(2,2)=-2.d0*muhs*k
      a(3,2)=-1.d0
      a(4,2)=2.d0*muhs*k
      a(5,2)=0.d0
      a(6,2)=0.d0
c
      a(1,3)=1.d0+et*cm/muhs
      a(2,3)=2.d0*et*cm*k
      a(3,3)=-1.d0-et*cx/muhs
      a(4,3)=-2.d0*et*cx*k
      a(5,3)=0.d0
      a(6,3)=0.d0
c
      a(1,4)=1.d0+et*cp/muhs
      a(2,4)=-2.d0*et*cp*k
      a(3,4)=1.d0-et*cx/muhs
      a(4,4)=2.d0*et*cx*k
      a(5,4)=0.d0
      a(6,4)=0.d0
c
      a(1,5)=0.d0
      a(2,5)=0.d0
      a(3,5)=0.d0
      a(4,5)=0.d0
      a(5,5)=1.d0
      a(6,5)=k*muhs
c
      a(1,6)=0.d0
      a(2,6)=0.d0
      a(3,6)=0.d0
      a(4,6)=0.d0
      a(5,6)=1.d0
      a(6,6)=-k*muhs
c
      return
      end
