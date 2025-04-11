      subroutine dghskern(yhs,k,lahs,muhs)
      use dgalloc
      implicit none
c
c     calculation of response of the half-space
c     yhs(8,4): solution vector
c     k: wave number
c
      real*8 k
      real*8 lahs,muhs,yhs(8,4)
c
c     work space
c
      integer*4 i,istp,j,l,lup,llw,key
      real*8 hl
      real*8 wave,wave2,cdet
      real*8 b(6,4),orth(3,3)
      real*8 yup(6,3),ylw(6,3),y1(6,3)
      real*8 c0(6,3),c1(6,3),yup0(6,3),ylw0(6,3)
      real*8 mat(6,6),mai(6,6),coef(6,6)
c
      call hsmatinv(mai,k,0.d0,lahs,muhs)
c===============================================================================
      lup=1
      llw=max0(ls,lzrec)
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
      if(lup.eq.lzrec)call memcpy(yup,yup0,18)
c
      do l=lup+1,ls
        hl=hp(l-1)
        wave=dexp(-k*hl)
        wave2=wave*wave
c
        call axb(mai,yup,6,6,3,c0)
c
c       orthogonalization
c
        cdet=c0(1,1)*c0(3,2)*c0(5,3)
     &      +c0(3,1)*c0(5,2)*c0(1,3)
     &      +c0(5,1)*c0(1,2)*c0(3,3)
     &      -c0(5,1)*c0(3,2)*c0(1,3)
     &      -c0(3,1)*c0(1,2)*c0(5,3)
     &      -c0(1,1)*c0(5,2)*c0(3,3)
        orth(1,1)=(c0(3,2)*c0(5,3)-c0(3,3)*c0(5,2))/cdet
        orth(2,1)=(c0(3,3)*c0(5,1)-c0(3,1)*c0(5,3))/cdet
        orth(3,1)=(c0(3,1)*c0(5,2)-c0(3,2)*c0(5,1))/cdet
        orth(1,2)=(c0(1,3)*c0(5,2)-c0(1,2)*c0(5,3))/cdet
        orth(2,2)=(c0(1,1)*c0(5,3)-c0(1,3)*c0(5,1))/cdet
        orth(3,2)=(c0(1,2)*c0(5,1)-c0(1,1)*c0(5,2))/cdet
        orth(1,3)=(c0(1,2)*c0(3,3)-c0(1,3)*c0(3,2))/cdet
        orth(2,3)=(c0(1,3)*c0(3,1)-c0(1,1)*c0(3,3))/cdet
        orth(3,3)=(c0(1,1)*c0(3,2)-c0(1,2)*c0(3,1))/cdet
c
        call axb(c0,orth,6,3,3,c1)
        if(l.gt.lzrec)then
c
c	  orthonormalization of the receiver vectors
c
          call axb(yup0,orth,6,3,3,y1)
          call memcpy(y1,yup0,18)
          do j=1,3
            do i=1,6
              yup0(i,j)=yup0(i,j)*wave
            enddo
          enddo
        endif
        c1(1,1)=1.d0
        c1(2,1)=c1(2,1)*wave2
        c1(3,1)=0.d0
        c1(4,1)=c1(4,1)*wave2
        c1(5,1)=0.d0
        c1(6,1)=c1(6,1)*wave2
c
        c1(1,2)=0.d0
        c1(2,2)=c1(2,2)*wave2
        c1(3,2)=1.d0
        c1(4,2)=c1(4,2)*wave2
        c1(5,2)=0.d0
        c1(6,2)=c1(6,2)*wave2
c
        c1(1,3)=0.d0
        c1(2,3)=c1(2,3)*wave2
        c1(3,3)=0.d0
        c1(4,3)=c1(4,3)*wave2
        c1(5,3)=1.d0
        c1(6,3)=c1(6,3)*wave2
c
        call hsmatrix(mat,k,hl,lahs,muhs)
        call axb(mat,c1,6,6,3,yup)
        if(l.eq.lzrec)call memcpy(yup,yup0,18)
      enddo
c
c===============================================================================
c
c     matrix propagation from halfspace to source
c
      do i=1,6
        do j=1,3
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
      call hsmatrix(mat,k,0.d0,lahs,muhs)
      call axb(mat,c0,6,6,3,ylw)
c
      if(llw.eq.lzrec)call memcpy(ylw,ylw0,18)
c
      do l=llw-1,ls,-1
        hl=hp(l)
        wave=dexp(-k*hl)
        wave2=wave*wave
c
        call axb(mai,ylw,6,6,3,c0)
c
c       orthogonalization
c
        cdet=c0(2,1)*c0(4,2)*c0(6,3)
     &      +c0(4,1)*c0(6,2)*c0(2,3)
     &      +c0(6,1)*c0(2,2)*c0(4,3)
     &      -c0(6,1)*c0(4,2)*c0(2,3)
     &      -c0(4,1)*c0(2,2)*c0(6,3)
     &      -c0(2,1)*c0(6,2)*c0(4,3)
        orth(1,1)=(c0(4,2)*c0(6,3)-c0(4,3)*c0(6,2))/cdet
        orth(2,1)=(c0(4,3)*c0(6,1)-c0(4,1)*c0(6,3))/cdet
        orth(3,1)=(c0(4,1)*c0(6,2)-c0(4,2)*c0(6,1))/cdet
        orth(1,2)=(c0(2,3)*c0(6,2)-c0(2,2)*c0(6,3))/cdet
        orth(2,2)=(c0(2,1)*c0(6,3)-c0(2,3)*c0(6,1))/cdet
        orth(3,2)=(c0(2,2)*c0(6,1)-c0(2,1)*c0(6,2))/cdet
        orth(1,3)=(c0(2,2)*c0(4,3)-c0(2,3)*c0(4,2))/cdet
        orth(2,3)=(c0(2,3)*c0(4,1)-c0(2,1)*c0(4,3))/cdet
        orth(3,3)=(c0(2,1)*c0(4,2)-c0(2,2)*c0(4,1))/cdet
        call axb(c0,orth,6,3,3,c1)
        if(l.lt.lzrec)then
c
c	  orthonormalization of the receiver vectors
c
          call axb(ylw0,orth,6,3,3,y1)
          call memcpy(y1,ylw0,18)
          do j=1,3
            do i=1,6
              ylw0(i,j)=ylw0(i,j)*wave
            enddo
          enddo
        endif
c
        c1(1,1)=c1(1,1)*wave2
        c1(2,1)=1.d0
        c1(3,1)=c1(3,1)*wave2
        c1(4,1)=0.d0
        c1(5,1)=c1(5,1)*wave2
        c1(6,1)=0.d0
c
        c1(1,2)=c1(1,2)*wave2
        c1(2,2)=0.d0
        c1(3,2)=c1(3,2)*wave2
        c1(4,2)=1.d0
        c1(5,2)=c1(5,2)*wave2
        c1(6,2)=0.d0
c
        c1(1,3)=c1(1,3)*wave2
        c1(2,3)=0.d0
        c1(3,3)=c1(3,3)*wave2
        c1(4,3)=0.d0
        c1(5,3)=c1(5,3)*wave2
        c1(6,3)=1.d0
c
        call hsmatrix(mat,k,-hl,lahs,muhs)
        call axb(mat,c1,6,6,3,ylw)
        if(l.eq.lzrec)call memcpy(ylw,ylw0,18)
      enddo
c
c===============================================================================
c
c     conditions on the source surface
c
c
c     source function
c
      do istp=1,4
        do i=1,6
          b(i,istp)=sfcths0(i,istp)+sfcths1(i,istp)*k
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
        print *,'warning in pegpsv: anormal exit from svd!'
        return
      endif
      if(lzrec.lt.ls)then
        do istp=1,4
          do i=1,6
            yhs(i,istp)=0.d0
            do j=1,3
              yhs(i,istp)=yhs(i,istp)+b(j,istp)*yup0(i,j)
            enddo
          enddo
        enddo
      else if(lzrec.gt.ls)then
        do istp=1,4
          do i=1,6
            yhs(i,istp)=0.d0
            do j=1,3
              yhs(i,istp)=yhs(i,istp)+b(j+3,istp)*ylw0(i,j)
            enddo
          enddo
        enddo
      else
        do istp=1,4
          do i=1,6
            yhs(i,istp)=0.d0
            do j=1,3
              yhs(i,istp)=yhs(i,istp)+0.5d0*(b(j,istp)*yup0(i,j)
     &                  +b(j+3,istp)*ylw0(i,j))
            enddo
          enddo
        enddo
      endif
c
      return
      end
