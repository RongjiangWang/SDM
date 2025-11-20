      subroutine sdmzhy(nps,nobs)
      implicit none
c
c     calculate location-dependent weights of roughness smoothing suggested by Zhang (2025)
c
c     Last modified: Potsdam, Oct, 2025, by R. Wang
c
      include 'sdmglob.h'
      integer*4 nps,nobs
c
      integer*4 ips,iobs
      real*8 zhysum
c
      if(izhy.eq.1)then
        zhysum=0.d0
        do ips=1,nps
          zhy(ips)=0.d0
          do iobs=1,nobs
            zhy(ips)=zhy(ips)
     &              +dspmdl(ips,iobs,1)**2+dspmdl(ips,iobs,2)**2
          enddo
          zhysum=zhysum+zhy(ips)
        enddo
        zhysum=zhysum/dble(nps)
        do ips=1,nps
          zhy(ips)=zhy(ips)/zhysum
        enddo
      else
        do ips=1,nps
          zhy(ips)=1.d0
        enddo
      endif
      return
      end
