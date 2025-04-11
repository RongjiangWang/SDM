      subroutine dgpsv(y,k)
      use dgalloc
      implicit none
c
c     calculation of response to p-sv source
c     y(8,4): solution vector
c     k: wave number
c
      real*8 k
      real*8 y(8,4)
c
c     work space
c
      integer*4 i,istp,j,l,n,lup,llw,key
      real*8 b(6,4),cy(6,4)
      real*8 yup(6,3),ylw(6,3)
      real*8 ma0(6,6),c0(6,3),coef(6,6)
      real*8 yup0(6,3),ylw0(6,3)
c
c===============================================================================
c
      lup=1
      llw=lp
c
      do l=lup,ls-1
        call dgmatrix(maup(1,1,l),k,hp(l),nno(l))
        call dgmatinv(maiup(1,1,l),k,0.d0,nno(l))
      enddo
      do l=ls,llw-1
        call dgmatrix(malw(1,1,l),k,-hp(l),nno(l))
        call dgmatinv(mailw(1,1,l),k,0.d0,nno(l))
      enddo
c
c     matrix propagation from surface to source
c
      do j=1,3
        do i=1,6
          c0(i,j)=0.d0
          yup(i,j)=0.d0
        enddo
      enddo
c
      yup(1,1)=1.d0
      yup(3,2)=1.d0
      yup(5,3)=1.d0
      yup(6,3)=k*yup(5,3)
      if(lup.eq.lzrec)call memcpy(yup,yup0,18)
c
      call dgproppsv(lup,ls,yup,yup0,k)
c
c===============================================================================
c
c     matrix propagation from half-space to source
c
      do j=1,3
        do i=1,6
          c0(i,j)=0.d0
          ylw(i,j)=0.d0
        enddo
      enddo
c
c     coefficient vectors in the half-space
c
      c0(2,1)=1.d0
      c0(4,2)=1.d0
      c0(6,3)=1.d0
      n=nno(llw)
      call dgmatrix(ma0,k,0.d0,n)
      call axb(ma0,c0,6,6,3,ylw)
c
      if(llw.eq.lzrec)call memcpy(ylw,ylw0,18)
c
      call dgproppsv(llw,ls,ylw,ylw0,k)
c
c===============================================================================
c
c     conditions on the source surface
c
c
c     source function
c
      do istp=1,4
        do i=1,4
          b(i,istp)=sfct0(i,istp)+sfct1(i,istp)*k
        enddo
        do i=5,6
          b(i,istp)=sfct0(i+2,istp)+sfct1(i+2,istp)*k
        enddo
      enddo
      do i=1,6
        do j=1,3
          coef(i,j)=yup(i,j)
          coef(i,j+3)=-ylw(i,j)
        enddo
      enddo
      key=0
      call svd(coef,b,6,4,0.d0,key)
      if(key.eq.0)then
        print *,'warning in dgpsv: anormal exit from svd!'
        return
      endif
      if(lzrec.lt.ls)then
        do istp=1,4
          do i=1,6
            cy(i,istp)=0.d0
            do j=1,3
              cy(i,istp)=cy(i,istp)+b(j,istp)*yup0(i,j)
            enddo
          enddo
        enddo
      else if(lzrec.gt.ls)then
        do istp=1,4
          do i=1,6
            cy(i,istp)=0.d0
            do j=1,3
              cy(i,istp)=cy(i,istp)+b(j+3,istp)*ylw0(i,j)
            enddo
          enddo
        enddo
      else
        do istp=1,4
          do i=1,6
            cy(i,istp)=0.d0
            do j=1,3
              cy(i,istp)=cy(i,istp)
     &             +0.5d0*(b(j,istp)*yup0(i,j)+b(j+3,istp)*ylw0(i,j))
            enddo
          enddo
        enddo
      endif
c
c     transform y6 -> gravity
c               y2(Eulerian) -> y2(Lagrangian)
c
      if(zrec.gt.0.d0)then
        n=nno(lzrec)
        do istp=1,4
          cy(6,istp)=cy(6,istp)+gamma*rho(n)*cy(1,istp)
          cy(2,istp)=cy(2,istp)-rho(n)*g0*cy(1,istp)
        enddo
      endif
      do istp=1,4
        do i=1,4
          y(i,istp)=cy(i,istp)
        enddo
        y(7,istp)=cy(5,istp)
        y(8,istp)=cy(6,istp)
      enddo
c
      return
      end
