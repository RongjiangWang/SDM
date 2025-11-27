      real*8 function sdmscost(nps,slip)
      implicit none
c
c     calculate the slip gradient part of the cost function
c
c     Last modified: Potsdam, Oct, 2008, by R. Wang
c
      include 'sdmglob.h'
      integer*4 nps
      real*8 slip(NPSMAX,2)
c
      integer*4 i,ips,jps
      real*8 scost,sd,sdsum
c
      scost=0.d0
      do jps=1,nps
        sdsum=0.d0
        do i=1,6
          sd=0.d0
          do ips=1,nps
            sd=sd+slip(ips,1)*dcgrn(ips,1,jps,i)
     &           +slip(ips,2)*dcgrn(ips,2,jps,i)
          enddo
          strdc(jps,i)=sd
          sdsum=sdsum+zhy(ips)*sd*sd
        enddo
        scost=scost+parea(jps)*sdsum
      enddo
c
      sdmscost=0.5d0*scost
      return
      end
