      double precision function sdmscost(nps,slip)
      implicit none
c
c     calculate the slip gradient part of the cost function
c
c     Last modified: Potsdam, Oct, 2008, by R. Wang
c
      include 'sdmglob.h'
      integer nps
      double precision slip(NPSMAX,2)
c
      integer i,ips,jps
      double precision scost,sd,sdsum
c
      scost=0.d0
      do ips=1,nps
        sdsum=0.d0
        do i=1,6
          sd=0.d0
          do jps=1,nps
            sd=sd+slip(jps,1)*dcgrn(jps,1,ips,i)
     &           +slip(jps,2)*dcgrn(jps,2,ips,i)
          enddo
          strdc(ips,i)=sd
          sdsum=sdsum+sd*sd
        enddo
        scost=scost+parea(ips)*sdsum
      enddo
c
      sdmscost=0.5d0*scost
      return
      end
