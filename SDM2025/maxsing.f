      real*8 function maxsing(mat,vec,swp,n,eps)
      implicit none
      integer*4 n
      real*8 eps
      real*8 mat(n,n),vec(n),swp(n)
c
c     dominant eigenvalue of a positive defined symmetric matrix
c     by the power method
c
      integer*4 i,j,iter
      real*8 a,b
c
      integer*4 niter
      data niter/10000/
c
      do i=1,n
        swp(i)=1.d0
      enddo
c
      b=0.d0
      do iter=1,niter
        a=0.d0
        do i=1,n
          vec(i)=0.d0
          do j=1,n
            vec(i)=vec(i)+mat(i,j)*swp(j) 
          enddo
          a=a+vec(i)**2
        enddo
        a=dsqrt(a/dble(n))
        do i=1,n
          swp(i)=vec(i)/a
        enddo
        if(dabs(a-b).le.eps*dabs(a))then
          goto 100
        else
          b=a
        endif
      enddo
      print *,' Warning in maxsing: convergence not achieved!'
100   continue
c
c     the dominant eigenvalue sigma2 of sysmat found
c
      maxsing=a
c
      return
      end