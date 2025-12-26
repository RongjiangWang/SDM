      real*8 function sdmcorl(ierr)
      use sdmalloc
      implicit none
c
c     Last modified: 2025 in Zhuhai by R. Wang
c
      integer*4 ierr
c
      integer*4 igd,iobs,iobs1,iobs2,ips
      real*8 b2,m2,bm,dm
c
      allocate(do0(ngd),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmcorl: do0 not allocated!'
      allocate(dm0(ngd),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmcorl: dm0 not allocated!'
c
      iobs1=0
      do igd=1,ngd
        iobs2=iobs1+nobsj(igd)
        do0(igd)=0.d0
        dm0(igd)=0.d0
        do iobs=iobs1+1,iobs2
          do0(igd)=do0(igd)+wf(iobs)*(datobs(iobs)-corrmdl(iobs))
          dm=0.d0
          do ips=1,nps
            dm=dm+slpmdl(1,ips)*datgrn(1,ips,iobs)
     &           +slpmdl(2,ips)*datgrn(2,ips,iobs)
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
          b2=b2+wf(iobs)*(datobs(iobs)-corrmdl(iobs)-do0(igd))**2
          dm=0.d0
          do ips=1,nps
            dm=dm+slpmdl(1,ips)*datgrn(1,ips,iobs)
     &           +slpmdl(2,ips)*datgrn(2,ips,iobs)
          enddo
          m2=m2+wf(iobs)*(dm-dm0(igd))**2
          bm=bm
     &      +wf(iobs)*(datobs(iobs)-corrmdl(iobs)-do0(igd))
     &               *(dm-dm0(igd))
        enddo
        iobs1=iobs2
      enddo
      if(m2.gt.0.d0)then
        sdmcorl=bm/dsqrt(b2*m2)
      else
        sdmcorl=0.d0
      endif
c
      deallocate(do0,dm0)
c
      return
      end
