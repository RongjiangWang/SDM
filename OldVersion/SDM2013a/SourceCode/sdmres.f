      subroutine sdmres(ngd,nps,doffset,rmsres,maxres,minres,meanres)
      implicit none
c
c     calculate rms residual
c
c     Last modified: Potsdam, Oct, 2008, by R. Wang
c
      include 'sdmglob.h'
c
      integer*4 ngd,nps
      real*8 rmsres
      real*8 doffset(NGDMAX)
      real*8 maxres(NGDMAX),minres(NGDMAX),meanres(NGDMAX)
c
      integer*4 igd,ips,iobs,iobs1,iobs2
c
      rmsres=0.d0
      iobs1=0
      do igd=1,ngd
        maxres(igd)=0.d0
        minres(igd)=0.d0
        meanres(igd)=0.d0
        iobs2=iobs1+nobsj(igd)
        do iobs=iobs1+1,iobs2
          dspres(iobs)=dspres(iobs)-ddsp(iobs)-doffset(igd)
          rmsres=rmsres+wf(iobs)*dspres(iobs)**2
          maxres(igd)=dmax1(maxres(igd),dspres(iobs))
          minres(igd)=dmin1(minres(igd),dspres(iobs))
          meanres(igd)=meanres(igd)+wf(iobs)*dspres(iobs)**2
        enddo
        meanres(igd)=dsqrt(meanres(igd)/wfmsum(igd))
        iobs1=iobs2
      enddo
      rmsres=dsqrt(rmsres/wfsum)
c
      return
      end
