      subroutine sdmzhy(nps,nobs,ismooth)
      implicit none
c
c     calculate location-dependent weights of roughness smoothing suggested by Zhang (2025)
c
c     Last modified: Potsdam, Oct, 2025, by R. Wang
c
      include 'sdmglob.h'
      integer*4 nps,nobs,ismooth
c
      integer*4 i,ira,ips,jps,iobs
      real*8 zhysum
c
      if(ismooth.eq.3)then
        zhysum=0.d0
        do ips=1,nps
          zhy(ips)=0.d0
          do iobs=1,nobs
            zhy(ips)=zhy(ips)
     &              +dspmdl(ips,iobs,1)**2+dspmdl(ips,iobs,2)**2
          enddo
          zhy(ips)=dsqrt(zhy(ips))
          zhysum=zhysum+zhy(ips)**2
        enddo
        zhysum=dsqrt(zhysum/dble(nps))
        do ips=1,nps
          zhy(ips)=zhy(ips)/zhysum
        enddo
        do ips=1,nps
          do ira=1,2
            do jps=1,nps
              do i=1,6
                dcgrn(ips,ira,jps,i)=dcgrn(ips,ira,jps,i)*zhy(jps)
              enddo
            enddo
          enddo
        enddo
      else
        do ips=1,nps
          zhy(ips)=1.d0
        enddo
      endif
      return
      end
