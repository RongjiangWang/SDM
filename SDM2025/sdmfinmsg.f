      subroutine sdmfinmsg(ierr)
      use sdmalloc
      implicit none
      integer*4 ierr
c
c     for final information messages
c
c     Last modified: Zhuhai, Nov. 2025, by R. Wang
c
      integer*4 is,ips
      real*8 slp
c
      allocate(slpp(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmfinmsg: slpp not allocated!'
      allocate(ipsp(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmfinmsg: ipsp not allocated!'
      allocate(slpm(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmfinmsg: slpm not allocated!'
      allocate(ssm(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmfinmsg: ssm not allocated!'
      allocate(dsm(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmfinmsg: dsm not allocated!'
      allocate(ram(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmfinmsg: ram not allocated!'
      allocate(rap(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmfinmsg: rap not allocated!'
c
      do is=1,ns
        slpp(is)=0.d0
        ipsp(is)=0
        slpm(is)=0.d0
        ssm(is)=0.d0
        dsm(is)=0.d0
        do ips=nps1(is),nps2(is)
          slp=dsqrt(slpmdl(1,ips)**2+slpmdl(2,ips)**2)
          if(slpp(is).lt.slp)then
            slpp(is)=slp
            ipsp(is)=ips
          endif
          slpm(is)=slpm(is)+slp
          ssm(is)=ssm(is)+slpmdl(1,ips)
          dsm(is)=dsm(is)+slpmdl(2,ips)
        enddo
        slpm(is)=slpm(is)/dble(1+nps2(is)-nps1(is))
        ssm(is)=ssm(is)/dble(1+nps2(is)-nps1(is))
        dsm(is)=dsm(is)/dble(1+nps2(is)-nps1(is))
        ram(is)=dmod(datan2(dsm(is),ssm(is))/DEG2RAD+rake360(is),360.d0)
        ips=ipsp(is)
        rap(is)=dmod(datan2(slpmdl(2,ips),slpmdl(1,ips))/DEG2RAD
     &         +rake360(is),360.d0)
      enddo
c
      return
      end
