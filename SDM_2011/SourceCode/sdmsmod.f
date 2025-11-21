      real*8 function sdmsmod(ns,nps,ismooth)
      implicit none
c
c     calculate the slip gradient part of the cost function
c
c     Last modified: Potsdam, Oct, 2008, by R. Wang
c
      include 'sdmglob.h'
      integer*4 ns,nps,ismooth
c
      integer*4 i,ips,jps
      real*8 area,sd(3)
c
      if(ismooth.eq.1)then
        sdmsmod=0.d0
        area=0.d0
        do ips=1,nps
          area=area+parea(ips)
          sdmsmod=sdmsmod+parea(ips)*(slpmdl(ips,1)**2+slpmdl(ips,2)**2)
        enddo
        sdmsmod=sdmsmod/area
      else
        sdmsmod=0.d0
        area=0.d0
        do ips=1,nps
          do i=1,3
            sd(i)=0.d0
            do jps=1,nps
              sd(i)=sd(i)
     &          +slpmdl(jps,1)*strgrn(jps,1,ips,i)
     &          +slpmdl(jps,2)*strgrn(jps,2,ips,i)
            enddo
          enddo
          area=area+parea(ips)
          sdmsmod=sdmsmod+parea(ips)*(sd(1)**2+sd(2)**2+sd(3)**2)
        enddo
        sdmsmod=sdmsmod/area
      endif
c
      return
      end
