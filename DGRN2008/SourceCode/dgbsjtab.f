      subroutine dgbsjtab(ierr)
      implicit none
      integer ierr
c
      include 'dgglob.h'
c
      integer i,j
      double precision x,xsqrt,pi2
      double precision bessj0,bessj1,bessj
c
      pi2=8.d0*datan(1.d0)
      do j=-1,3
        bsjfct(0,j)=0.d0
      enddo
      dxbsj=pi2/dble(ndbsj)
      do i=1,nnbsj1
        x=dxbsj*dble(i)
        xsqrt=dsqrt(x)
        bsjfct(i,0)=xsqrt*bessj0(x)
        bsjfct(i,1)=xsqrt*bessj1(x)
        bsjfct(i,2)=xsqrt*bessj(2,x)
        bsjfct(i,3)=xsqrt*bessj(3,x)
        bsjfct(i,-1)=-bsjfct(i,1)
      enddo
c
      ierr=0
      return
      end
