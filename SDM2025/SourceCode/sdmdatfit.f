      subroutine sdmdatfit(ierr)
      use sdmalloc
      implicit none
      integer*4 ierr
c
c     first step to prepare SDM iteration
c     Last modified: Zhuhai, Nov. 2025, by R. Wang
c
      integer*4 i,j,m,n,ira,ips,jps,igd,iobs,iusrp
      real*8 dswp,res,resvar,wfsum,wfgd,sdmod
      real*8 sd(6)
c
      do iobs=1,nobs
        datmdl(iobs)=0.d0
        do ips=1,nps
          do ira=1,2
            datmdl(iobs)=datmdl(iobs)
     &                  +slpmdl(ira,ips)*datgrn(ira,ips,iobs)
          enddo
        enddo
      enddo
      do iobs=1,nobs
        corrmdl(iobs)=0.d0
        do iusrp=1,nusrp
          corrmdl(iobs)=corrmdl(iobs)
     &                 +corrusrp(iusrp)*corrgrn(iusrp,iobs)
        enddo
      enddo
c
      rmsresall=0.d0
      wfsum=0.d0
      do igd=1,ngd
        rmsres(igd)=0.d0
        resmin(igd)=0.d0
        resmax(igd)=0.d0
        wfgd=0.d0
        do iobs=nobs1(igd),nobs2(igd)
          res=datobs(iobs)-datmdl(iobs)-corrmdl(iobs)
c
          resmin(igd)=dmax1(resmin(igd),-res)
          resmax(igd)=dmax1(resmax(igd),+res)
c
          resvar=(wf(iobs)*res)**2
          wfgd=wfgd+wf(iobs)**2
c
          rmsresall=rmsresall+resvar
          rmsres(igd)=rmsres(igd)+resvar
        enddo
        resmin(igd)=-resmin(igd)
        wfsum=wfsum+wfgd
        rmsres(igd)=dsqrt(rmsres(igd)/wfgd)
      enddo
      rmsresall=dsqrt(rmsresall/wfsum)/rmsdatall
c
      roughness=0.d0
      do jps=1,nps
        sdmod=0.d0
        do i=1,nsmocmp
          sd(i)=0.d0
          do ips=1,nps
            do ira=1,2
              m=(jps-1)*nsmocmp+i
              n=(ips-1)*2+ira
              sd(i)=sd(i)
     &             +slpmdl(ira,ips)*smogrnmat(m,n)*zhy(ips)/zhy(jps)
            enddo
          enddo
          sdmod=sdmod+sd(i)*sd(i)
        enddo
        roughness=roughness+sdmod
      enddo
c
      if(ismooth.eq.2)then
        sdmod=0.d0
        do jps=1,nps
          do i=1,3
            sd(i)=0.d0
            do ips=1,nps
              sd(i)=sd(i)+slpmdl(1,ips)*strgrn(1,ips,i,jps)
     &                   +slpmdl(2,ips)*strgrn(2,ips,i,jps)
            enddo
          enddo
          sdmod=sdmod+parea(jps)*(sd(1)**2+sd(2)**2+sd(3)**2)
        enddo
      else
        sdmod=0.d0
        do ips=1,nps
          sdmod=sdmod+parea(ips)*(slpmdl(1,ips)**2+slpmdl(2,ips)**2)
        enddo
      endif
c
      roughness=roughness/sdmod
c
      return
      end