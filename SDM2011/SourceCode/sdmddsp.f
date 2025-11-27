      subroutine sdmddsp(nobs,nps,slip)
      implicit none
c
c     Last modified: Potsdam, Oct, 2008, by R. Wang
c
      include 'sdmglob.h'
      integer*4 nobs,nps
      real*8 slip(NPSMAX,2)
c
      integer*4 ips,iobs
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
