      subroutine svd(a,b,n,m,eps,key)
      implicit none
c-------------------------------------------------------------------------------
c     Solve linear equation system by single-value decomposiztion              i
c     method (modified from ludcmp and lubksb in the <Numerical Recipies>      i
c     a: coefficient matrix(n,n);                                              i
c     b: right-hand matrix(n,m) by input,                                      i
c         solution matrix(n,m) by return;                                      i
c     eps: control constant;                                                   i
c     key: if the main term of a column is                                     i
c          smaller than eps, key=0: anormal return,                            i
c          else key=1: normal return.                                          i
c                                                                              i
c     Note: n <= 500 will be NOT CHECKED!                                      i
c-------------------------------------------------------------------------------
      integer*4 n,m,key
      real*8 eps
      real*8 a(n,n),b(n,m)
c
      integer*4 i,ii,imax,j,k,ll
      real*8 aamax,dum,sum
c
      integer*4, allocatable:: indx(:)
      real*8, allocatable:: vv(:)
c
      key=0
      allocate(indx(n),stat=key)
      if(key.ne.0)stop ' Error in svd: indx not allocated!'
      allocate(vv(n),stat=key)
      if(key.ne.0)stop ' Error in svd: vv not allocated!'
c
      do i=1,n
        aamax=0.d0
        do j=1,n
          aamax=dmax1(aamax,dabs(a(i,j)))
        enddo
        if(aamax.le.eps)then
          key=0
          goto 100
        endif
        vv(i)=1.d0/aamax
      enddo
      do j=1,n
        do i=1,j-1
          sum=a(i,j)
          do k=1,i-1
            sum=sum-a(i,k)*a(k,j)
          enddo
          a(i,j)=sum
        enddo
        aamax=0.d0
        do i=j,n
          sum=a(i,j)
          do k=1,j-1
            sum=sum-a(i,k)*a(k,j)
          enddo
          a(i,j)=sum
          dum=vv(i)*dabs(sum)
          if(dum.ge.aamax)then
            imax=i
            aamax=dum
          endif
        enddo
        if(j.ne.imax)then
          do k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
          enddo
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(dabs(a(j,j)).le.eps)then
          key=0
          goto 100
        endif
        if(j.ne.n)then
          dum=1.d0/a(j,j)
          do i=j+1,n
            a(i,j)=a(i,j)*dum
          enddo
        endif
      enddo
c
      do k=1,m
        ii=0
        do i=1,n
          ll=indx(i)
          sum=b(ll,k)
          b(ll,k)=b(i,k)
          if(ii.ne.0)then
            do j=ii,i-1
              sum=sum-a(i,j)*b(j,k)
            enddo
          else if(dabs(sum).ne.0.d0)then
            ii=i
          endif
          b(i,k)=sum
        enddo
        do i=n,1,-1
          sum=b(i,k)
          do j=i+1,n
            sum=sum-a(i,j)*b(j,k)
          enddo
          b(i,k)=sum/a(i,i)
        enddo
      enddo
      key=1
100   deallocate(indx,vv)
      return
      end