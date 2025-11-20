      real*8 function sdmcorr(ngd,nobs,nps)
      implicit none
c
c     Last modified: 2025 in Zhuhai by R. Wang
c
      include 'sdmglob.h'
      integer*4 ngd,nobs,nps
c
      integer*4 igd,iobs,iobs1,iobs2,ips
      real*8 b2,m2,bm,dm
c
      iobs1=0
      do igd=1,ngd
        iobs2=iobs1+nobsj(igd)
        do0(igd)=0.d0
        dm0(igd)=0.d0
        do iobs=iobs1+1,iobs2
          do0(igd)=do0(igd)+wf(iobs)*(dspobs(iobs)-offset(igd))
          dm=0.d0
          do ips=1,nps
            dm=dm+slpmdl(ips,1)*dspmdl(ips,iobs,1)
     &           +slpmdl(ips,2)*dspmdl(ips,iobs,2)
          enddo
          dm0(igd)=dm0(igd)+wf(iobs)*dm
        enddo
        iobs1=iobs2
      enddo
c
      b2=0.d0
      m2=0.d0
      bm=0.d0
      iobs1=0
      do igd=1,ngd
        iobs2=iobs1+nobsj(igd)
        do iobs=iobs1+1,iobs2
          b2=b2+wf(iobs)*(dspobs(iobs)-offset(igd)-do0(igd))**2
          dm=0.d0
          do ips=1,nps
            dm=dm+slpmdl(ips,1)*dspmdl(ips,iobs,1)
     &           +slpmdl(ips,2)*dspmdl(ips,iobs,2)
          enddo
          m2=m2+wf(iobs)*(dm-dm0(igd))**2
          bm=bm
     &      +wf(iobs)*(dspobs(iobs)-offset(igd)-do0(igd))*(dm-dm0(igd))
        enddo
        iobs1=iobs2
      enddo
      if(m2.gt.0.d0)then
        sdmcorr=bm/dsqrt(b2*m2)
      else
        sdmcorr=0.d0
      endif
      return
      end
