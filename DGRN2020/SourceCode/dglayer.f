      subroutine dglayer(ierr)
      use dgalloc
      implicit none
      integer*4 ierr
c
      integer*4 l,n,li,lp0
      real*8 x1
c
      lp0=1
      xh0(lp0)=0.d0
      do n=1,n0-1
        lp0=lp0+1
        xh0(lp0)=xh0(lp0-1)+h(n)
      enddo
      lp0=lp0+1
      xh0(lp0)=zrec
      lp0=lp0+1
      xh0(lp0)=zs
c
c     sort the z0-profile
c
      do l=1,lp0-1
        do li=l+1,lp0
          if(xh0(li).lt.xh0(l))then
            x1=xh0(l)
            xh0(l)=xh0(li)
            xh0(li)=x1
          endif
        enddo
      enddo
c
c     delete duplicates
c
      lp=1
      xh(lp)=0.d0
      do l=2,lp0
        if(xh0(l).gt.xh(lp))then
          hp(lp)=xh0(l)-xh(lp)
          lp=lp+1
          xh(lp)=xh0(l)
        endif
      enddo
      hp(lp)=0.d0
c
c     determine ls,lzrec
c
      do l=1,lp
        if(xh(l).eq.zs)ls=l
        if(xh(l).eq.zrec)lzrec=l
      enddo
c
c     determine layer no of each depth
c
      li=1
      x1=h(1)
      nno(1)=1
      do l=2,lp
        if(xh(l).ge.x1.and.li.lt.n0)then
          li=li+1
          x1=x1+h(li)
        endif
        nno(l)=li
      enddo
c
      ierr=0
c
      return
      end
