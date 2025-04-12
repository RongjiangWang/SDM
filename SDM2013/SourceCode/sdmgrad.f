      subroutine sdmgrad(nps,nobs,sdgrad,wgrad)
      implicit none
c
c     calculate slip gradient at the ips-th patch - derivative
c     of the gradient part of the cost function
c
c     Last modified: Potsdam, Oct, 2008, by R. Wang
c
      include 'sdmglob.h'
c
      integer*4 nobs,nps
      real*8 wgrad
      real*8 sdgrad(NPSMAX,2)
c
      integer*4 i,iobs,ips,jps
      real*8 og,sg
c
      do ips=1,nps
        do i=1,2
c
c         derivative of misfit part of cost function
c
          og=0.d0
          do iobs=1,nobs
            og=og+wf(iobs)*dspmdl(ips,iobs,i)*dspres(iobs)
          enddo
c
c         derivative of the weighting factor for stress drop
c
          sg=0.d0
          do jps=1,nps
            sg=sg+parea(jps)
     &        *(dcgrn(ips,i,jps,1)*strdc(jps,1)
     &         +dcgrn(ips,i,jps,2)*strdc(jps,2)
     &         +dcgrn(ips,i,jps,3)*strdc(jps,3)
     &         +dcgrn(ips,i,jps,4)*strdc(jps,4)
     &         +dcgrn(ips,i,jps,5)*strdc(jps,5)
     &         +dcgrn(ips,i,jps,6)*strdc(jps,6))
          enddo
          sdgrad(ips,i)=og-wgrad*sg
        enddo
      enddo
c
      return
      end