      subroutine reweft(nmax,n,nd,lat1,lon1,lat2,lon2,dip,step,z)
      implicit none
c
      integer*4 nmax,n,nd
      real*8 step,z
      real*8 lat1(nmax),lon1(nmax)
      real*8 lat2(nmax),lon2(nmax)
      real*8 dip(nmax)
c
c     Last modified: Potsdam, March, 2013, by R. Wang
c
      integer*4 nn
      parameter(nn=10000)
c
      integer*4 i,j,k,l,nl
      real*8 x,y,length,delta,width,widthmax,diptop,dipbtm
      real*8 x1(nn),y1(nn),x2(nn),y2(nn),d(nn)
      real*8 lenk(nn)
c
      real*8 fdip
c
      real*8 DEG2RAD
      data DEG2RAD/1.745329252d-02/
c
      delta=0.1d0*step
c
      k=1
      x1(1)=lat1(1)
      y1(1)=lon1(1)
      x2(1)=lat2(1)
      y2(1)=lon2(1)
      d(1)=dip(1)
      lenk(1)=0.d0
      do i=2,n
        call disazi(6.371d+06,lat1(i-1),lon1(i-1),lat1(i),lon1(i),x,y)
        length=dsqrt(x**2+y**2)
        nl=1+idint(length/delta)
        do l=1,nl
          k=k+1
          if(k.gt.nn)then
            stop 'Error in reweft: nn too small!'
          endif
          x1(k)=lat1(i-1)+(lat1(i)-lat1(i-1))*dble(l)/dble(nl)
          y1(k)=lon1(i-1)+(lon1(i)-lon1(i-1))*dble(l)/dble(nl)
          x2(k)=lat2(i-1)+(lat2(i)-lat2(i-1))*dble(l)/dble(nl)
          y2(k)=lon2(i-1)+(lon2(i)-lon2(i-1))*dble(l)/dble(nl)
          d(k)=dip(i-1)+(dip(i)-dip(i-1))*dble(l)/dble(nl)
          lenk(k)=lenk(k-1)+length/dble(nl)
        enddo
      enddo
      n=2+idint(lenk(k)/step)
      delta=lenk(k)/dble(n-1)
      if(n.gt.nmax)then
        stop 'Error in reweft: nmax too small!'
      endif
c
      do i=2,n-1
        do j=2,k
          if(lenk(j).ge.dble(i-1)*delta)then
            lat1(i)=x1(j)
            lon1(i)=y1(j)
            lat2(i)=x2(j)
            lon2(i)=y2(j)
            dip(i)=d(j)
            goto 10
          endif
        enddo
10      continue
      enddo
      lat1(n)=x1(k)
      lon1(n)=y1(k)
      lat2(n)=x2(k)
      lon2(n)=y2(k)
      dip(n)=d(k)
      if(nd.eq.1)then
        widthmax=0.d0
        do i=1,n
          call disazi(6.371d+06,lat1(i),lon1(i),lat2(i),lon2(i),x,y)
          x=dsqrt(x*x+y*y)
          diptop=dip(i)
          dipbtm=fdip(dip(i),x,z)
          width=(z/((dipbtm-diptop)*DEG2RAD))
     &       *dlog(dtan(0.5d0*dipbtm*DEG2RAD)
     &            /dtan(0.5d0*diptop*DEG2RAD))
          widthmax=dmax1(widthmax,width)
        enddo
        nd=1+idint(widthmax/step)
      else
        nd=0
      endif
      return
      end
