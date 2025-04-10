      subroutine sdmddsp(nobs,nps,slip)
      implicit none
c
c     calculate the slip gradient part of the cost function
c
c     Last modified: Potsdam, Oct, 2008, by R. Wang
c
      include 'sdmglob.h'
      integer nobs,nps
      double precision slip(NPSMAX,2)
c
      integer ips,iobs
c
      do iobs=1,nobs
        ddsp(iobs)=0.d0
        do ips=1,nps
          ddsp(iobs)=ddsp(iobs)+slip(ips,1)*dspmdl(ips,iobs,1)
     &                         +slip(ips,2)*dspmdl(ips,iobs,2)
        enddo
      enddo
      return
      end
