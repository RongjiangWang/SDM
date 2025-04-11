      subroutine dgkern(y,k,lahs,muhs)
      use dgalloc
      implicit none
c
c     calculation of response function in Laplace domain
c     y(8,4): solution vector
c     k0: wave number (input)
c     cs: precision Laplace variable
c
      real*8 k
      real*8 lahs,muhs,y(8,4)
c
      integer*4 i,istp
      real*8 yhs(8,4)
c
      do istp=1,4
        do i=1,8
          y(i,istp)=0.d0
          yhs(i,istp)=0.d0
        enddo
      enddo
c
      call dgpsv(y,k)
      call dgsh(y,k)
c
c     subtract the halfspace solution
c
      call dghskern(yhs,k,lahs,muhs)
      do istp=1,4
        do i=1,8
          y(i,istp)=y(i,istp)-yhs(i,istp)
        enddo
      enddo
c
      return
      end	  
