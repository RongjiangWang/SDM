      subroutine sdmzhy(ierr)
      use sdmalloc
      implicit none
      integer*4 ierr
c
c     calculate location-dependent weights of roughness smoothing suggested by Zhang (2025)
c
c     Last modified: Potsdam, Oct, 2025, by R. Wang
c
      integer*4 i,ira,ips,jps,iobs
      real*8 zhysum
c
      allocate(zhy(nps),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmzhy: zhy not allocated!'
c
      if(izhy.eq.1)then
        zhysum=0.d0
        do ips=1,nps
          zhy(ips)=0.d0
          do iobs=1,nobs
            zhy(ips)=zhy(ips)
     &              +datgrn(1,ips,iobs)**2+datgrn(2,ips,iobs)**2
          enddo
          zhy(ips)=dsqrt(zhy(ips)/dble(nobs))
          zhysum=zhysum+zhy(ips)**2
        enddo
        zhysum=dsqrt(zhysum/dble(nps))
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
