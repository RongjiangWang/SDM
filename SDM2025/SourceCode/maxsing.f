      real*8 function maxsing(mat,n,eps,vecini,ierr)
      implicit none
      integer*4 n,ierr
      real*8 eps
      real*8 mat(n,n),vecini(n)
c
c     dominant eigenvalue of a positive defined symmetric matrix
c     by the power method
c
      integer*4 i,j,itr
      real*8 a,b
      real*8, allocatable:: vec(:),swp(:)
c
      integer*4 nitr
      data nitr/10000/
c
      allocate(vec(n),stat=ierr)
      if(ierr.ne.0)stop ' Error in maxsing: vec not allocated!'
      allocate(swp(n),stat=ierr)
      if(ierr.ne.0)stop ' Error in maxsing: vec not allocated!'
c
      do i=1,n
        swp(i)=vecini(i)
      enddo
c
      b=0.d0
      do itr=1,nitr
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
      deallocate(vec,swp)
c
      return
      end