      subroutine dgsh(y,k)
      use dgalloc
      implicit none
c
c     calculation of response to sh source
c     y(8,4): solution vector
c     k: wave number
c
      real*8 k
      real*8 y(8,4)
c
c     work space
c
      integer*4 i,istp,l,n,lup,llw,key
      real*8 b(2,4)
      real*8 yup(2),ylw(2),yup0(2),ylw0(2),coef(2,2)
c
c===============================================================================
c
c     matrix propagation from surface to source
c
      do i=1,2
        yup(i)=0.d0
      enddo
c
c     yup: the starting solution vector
c
      lup=1
c
c     determination of starting sublayer for half-space
c
      llw=lp
c
      do l=lup,ls-1
        call dghksh(hkup(1,1,l),2,k,hp(l),nno(l))
      enddo
      do l=ls,llw-1
        call dghksh(hklw(1,1,l),2,k,-hp(l),nno(l))
      enddo
c
      yup(1)=1.d0
      if(lup.gt.1)then
        n=nno(lup-1)
        yup(2)=mu(n)*k
      endif
      if(lup.eq.lzrec)call memcpy(yup,yup0,2)
c
      call dgpropsh(lup,ls,k,yup,yup0)
c
c===============================================================================
c
c     matrix propagation from half-space to source
c
      do i=1,2
        ylw(i)=0.d0
      enddo
c
c     ylw: the starting solution vector
c
      n=nno(llw)
      ylw(1)=1.d0
      ylw(2)=-mu(n)*k
      if(llw.eq.lzrec)call memcpy(ylw,ylw0,2)
c
      call dgpropsh(llw,ls,k,ylw,ylw0)
c
c===============================================================================
c
c     conditions on the source surface
c
c
c     point source function
c
      do istp=1,4
        do i=1,2
          b(i,istp)=sfct0(i+4,istp)+sfct1(i+4,istp)*k
        enddo
      enddo
      do i=1,2
        coef(i,1)=yup(i)
        coef(i,2)=-ylw(i)
      enddo
      key=0
      call svd(coef,b,2,4,0.d0,key)
      if(key.eq.0)then
        print *,'warning in dgsh: anormal exit from svd!'
        return
      endif
      if(lzrec.lt.ls)then
        do istp=1,4
          do i=1,2
            y(i+4,istp)=b(1,istp)*yup0(i)
          enddo
        enddo
      else if(lzrec.gt.ls)then
        do istp=1,4
          do i=1,2
            y(i+4,istp)=b(2,istp)*ylw0(i)
          enddo
        enddo
      else
        do istp=1,4
          do i=1,2
            y(i+4,istp)=0.5d0*(b(1,istp)*yup0(i)+b(2,istp)*ylw0(i))
          enddo
        enddo
      endif
      return
      end
