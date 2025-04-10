      double precision function sdmcorr(ngd,nobs,nps)
      implicit none
c
c     calculate the slip gradient part of the cost function
c
c     Last modified: Potsdam, Oct, 2008, by R. Wang
c
      include 'sdmglob.h'
      integer ngd,nobs,nps
c
      integer igd,iobs,iobs1,iobs2,ips
      double precision wei,b0,m0,b2,m2,bm,dm
c
      b2=0.d0
      m2=0.d0
      bm=0.d0
      iobs1=0
      do igd=1,ngd
        iobs2=iobs1+nobsj(igd)
        b0=0.d0
        m0=0.d0
        wei=0.d0
        do iobs=iobs1+1,iobs2
          b0=wf(iobs)*(dspobs(iobs)-offset(igd))
          dm=0.d0
          do ips=1,nps
            dm=dm+slpmdl(ips,1)*dspmdl(ips,iobs,1)
     &           +slpmdl(ips,2)*dspmdl(ips,iobs,2)
          enddo
          m0=m0+wf(iobs)*dm
          wei=wei+wf(iobs)
        enddo
        m0=m0/wei
        b0=b0/wei
        do iobs=iobs1+1,iobs2
          b2=b2+wf(iobs)*(dspobs(iobs)-offset(igd)-b0)**2
          dm=0.d0
          do ips=1,nps
            dm=dm+slpmdl(ips,1)*dspmdl(ips,iobs,1)
     &           +slpmdl(ips,2)*dspmdl(ips,iobs,2)
          enddo
          m2=m2+wf(iobs)*(dm-m0)**2
          bm=bm+wf(iobs)*(dspobs(iobs)-offset(igd)-b0)*(dm-m0)
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
