      subroutine sdmscale(ngd,nobs,ns,nps,wgrad,
     &                    dslpmdl,doffset,converge)
      implicit none
c
c     scale the downhill step by least-squares fitting to the data
c
c     Last modified: Potsdam, Oct, 2008, by R. Wang
c
      include 'sdmglob.h'
c
      integer ngd,nobs,ns,nps
      double precision wgrad
      double precision dslpmdl(NPSMAX,2),doffset(NGDMAX)
      logical converge
c
      integer i,imin,is,ips,igd,iobs,iobs1,iobs2
      double precision a,b,slp0,dslp,slpabs,dsp,fac
      double precision x,x0,dx,y,y0,dy,ra,cs,ss,cg,cg1,cg2
      double precision a1,a1s,a2,a2s,b1,b2,ab,db1,db2,dab,det
      double precision alfa,beta,scostp,scostm,scost0
      double precision xp(0:5),yp(0:5)
      double precision swapslp(NPSMAX,2),swapoffs(NGDMAX)
c
      double precision sdmocost,sdmscost
c
      double precision eps,alfamin
      data eps,alfamin/1.0d-02,1.0d-06/
c
      converge=.false.
      slp0=0.d0
      dslp=0.d0
      do ips=1,nps
        slp0=slp0+slpmdl(ips,1)**2+slpmdl(ips,2)**2
        dslp=dslp+dslpmdl(ips,1)**2+dslpmdl(ips,2)**2
      enddo
      slp0=dsqrt(slp0/dble(nps))
      dslp=dsqrt(dslp/dble(nps))
c
      if(dslp.le.0.d0)then
        do igd=1,ngd
          doffset(igd)=0.d0
        enddo
        converge=.true.
        return
      else if(slp0.gt.0.d0)then
c
c       normalize the gradient of the cost function before the scaling
c
        fac=slp0/dslp
        do ips=1,nps
          dslpmdl(ips,1)=dslpmdl(ips,1)*fac
          dslpmdl(ips,2)=dslpmdl(ips,2)*fac
        enddo
      endif
c
c     cost = f0 - a1*alfa + 0.5*a2*alfa**2 - b1*beta + 0.5*b2*beta**2 + ab*alfa*beta
c               - a1s*alfa + 0.5*a2s*alfa**2
c
      call sdmddsp(nobs,nps,dslpmdl)
      do ips=1,nps
        swapslp(ips,1)=slpmdl(ips,1)
        swapslp(ips,2)=slpmdl(ips,2)
      enddo
      scost0=sdmscost(nps,swapslp)
      do ips=1,nps
        swapslp(ips,1)=slpmdl(ips,1)+dslpmdl(ips,1)
        swapslp(ips,2)=slpmdl(ips,2)+dslpmdl(ips,2)
      enddo
      scostp=sdmscost(nps,swapslp)
      do ips=1,nps
        swapslp(ips,1)=slpmdl(ips,1)-dslpmdl(ips,1)
        swapslp(ips,2)=slpmdl(ips,2)-dslpmdl(ips,2)
      enddo
      scostm=sdmscost(nps,swapslp)
      a2s=scostp+scostm-2.d0*scost0
      a1s=0.5d0*(scostm-scostp)
c
      if(caloffset.and.slp0.gt.0.d0)then
        a1=0.d0
        a2=0.d0
        b1=0.d0
        b2=0.d0
        ab=0.d0
c
        iobs1=0
        do igd=1,ngd
          iobs2=iobs1+nobsj(igd)
          db1=0.d0
          dab=0.d0
          do iobs=iobs1+1,iobs2
            a1=a1+wf(iobs)*dspres(iobs)*ddsp(iobs)
            a2=a2+wf(iobs)*ddsp(iobs)**2
            db1=db1+wf(iobs)*dspres(iobs)
            dab=dab+wf(iobs)*ddsp(iobs)
          enddo
          b1=b1+db1*doffset(igd)
          b2=b2+wfmsum(igd)*doffset(igd)**2
          ab=ab+dab*doffset(igd)
          iobs1=iobs2
        enddo
c
        a1=a1+wgrad*a1s
        a2=a2+wgrad*a2s
c
        det=a2*b2-ab**2
        if(dabs(det).gt.0.d0)then
          alfa=(a1*b2-b1*ab)/det
          beta=(b1*a2-a1*ab)/det
        else
          alfa=0.d0
          beta=0.d0
        endif
      else
        a1=0.d0
        a2=0.d0
c
        iobs1=0
        do igd=1,ngd
          iobs2=iobs1+nobsj(igd)
          do iobs=iobs1+1,iobs2
            a1=a1+wf(iobs)*dspres(iobs)*ddsp(iobs)
            a2=a2+wf(iobs)*ddsp(iobs)**2
          enddo
          iobs1=iobs2
        enddo
c
        a1=a1+wgrad*a1s
        a2=a2+wgrad*a2s
c
        if(dabs(a2).gt.0.d0)then
          alfa=a1/a2
        else
          beta=0.d0
        endif
        beta=0.d0
      endif
c
      if(slp0.gt.0.d0.and.alfa.le.alfamin)then
        do ips=1,nps
          dslpmdl(ips,1)=0.d0
          dslpmdl(ips,2)=0.d0
        enddo
        do igd=1,ngd
          doffset(igd)=0.d0
        enddo
        converge=.true.
        return
      endif
c
      dslp=0.d0
      do ips=1,nps
        dslpmdl(ips,1)=dslpmdl(ips,1)*alfa
        dslpmdl(ips,2)=dslpmdl(ips,2)*alfa
        dslp=dslp+dslpmdl(ips,1)**2+dslpmdl(ips,2)**2
      enddo
      dslp=eps*dsqrt(dslp/dble(nps))
c
      if(slp0.eq.0.d0)scostref=sdmscost(nps,dslpmdl)
c
      do igd=1,ngd
        doffset(igd)=doffset(igd)*beta
      enddo
c
c     use a-priori conditions for slip amplitude and rake
c
      do is=1,ns
        do ips=nps1(is),nps2(is)
          dx=dslpmdl(ips,1)
          dy=dslpmdl(ips,2)
          if(dx**2+dy**2.gt.0.d0)then
            x0=slpmdl(ips,1)
            y0=slpmdl(ips,2)
            x=x0+dx
            y=y0+dy
            slpabs=dsqrt(x**2+y**2)
            ra=dmod(datan2(y,x)/DEG2RAD+rake360(is),360.d0)
            cs=dcos(ra*DEG2RAD)
            ss=dsin(ra*DEG2RAD)
            if(rake1(is).eq.rake2(is))then
              if(cs*cs1(is)+ss*ss1(is).le.0.d0)then
                x=0.d0
                y=0.d0
              else if(slpabs.ge.maxslip(is))then
                x=maxslip(is)*cs
                y=maxslip(is)*ss
              endif
c
c             classify slip verctors:
c               0 = zero;
c               1 = max. amplitude;
c               2 = within amplitude range.
c
              slpabs=dsqrt(x**2+y**2)
              if(slpabs.le.dslp)then
                slppos(ips)=0
              else if(slpabs.ge.maxslip(is)-dslp)then
                slppos(ips)=1
              else
                slppos(ips)=2
              endif
            else
              if(rake2(is).ge.rake1(is)+359.999999d0.or.
     &           ra.gt.rake1(is).and.ra.lt.rake2(is))then
                slpabs=dmin1(slpabs,maxslip(is))
                x=slpabs*cs
                y=slpabs*ss
              else
                cg1=cs*cs1(is)+ss*ss1(is)
                cg2=cs*cs2(is)+ss*ss2(is)
                cg=dmax1(cg1,cg2)
                if(cg.le.0.d0)then
                  x=0.d0
                  y=0.d0
                else
                  slpabs=dmin1(slpabs*cg,maxslip(is))
                  if(cg1.ge.cg2)then
                    x=slpabs*cs1(is)
                    y=slpabs*ss1(is)
                  else
                    x=slpabs*cs2(is)
                    y=slpabs*ss2(is)
                  endif
                endif
              endif
c
c             classify slip verctors:
c               0 = zero;
c               1 = max. amplitude and min. rake;
c               2 = max. amplitude and max. rake;
c               3 = within amplitude range and min. rake;
c               4 = within amplitude range and max. rake;
c               5 = max. amplitude and within rake range;
c               6 = within amplitude and rake range.
c
              slppos(ips)=6
              slpabs=dsqrt(x**2+y**2)
              xp(0)=0.d0
              yp(0)=0.d0
              xp(1)=maxslip(is)*cs1(is)
              yp(1)=maxslip(is)*ss1(is)
              xp(2)=maxslip(is)*cs2(is)
              yp(2)=maxslip(is)*ss2(is)
              xp(3)=slpabs*cs1(is)
              yp(3)=slpabs*ss1(is)
              xp(4)=slpabs*cs2(is)
              yp(4)=slpabs*ss2(is)
              if(slpabs.gt.0.d0)then
                xp(5)=maxslip(is)*x/slpabs
                yp(5)=maxslip(is)*y/slpabs
              endif
              do i=0,5
                if(dsqrt((x-xp(i))**2+(y-yp(i))**2).le.dslp)then
                  slppos(ips)=i
                  goto 400
                endif
              enddo
400           continue
              if(rake2(is).ge.rake1(is)+359.999999d0.and.
     &           slppos(ips).ne.5)slppos(ips)=6
            endif
            dslpmdl(ips,1)=x-x0
            dslpmdl(ips,2)=y-y0
          endif
        enddo
      enddo
      return
      end
