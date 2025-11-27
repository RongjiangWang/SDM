      subroutine sdmdslip(ngd,nobs,ns,nps,wgrad,dslpmdl,doffset)
      implicit none
c
c     calculate maximum gradient of cost function within the variable range
c
c     Last modified: Potsdam, Oct, 2008, by R. Wang
c
      include 'sdmglob.h'
c
      integer*4 ngd,nobs,ns,nps
      real*8 wgrad
      real*8 dslpmdl(NPSMAX,2),doffset(NGDMAX)
c
      integer*4 is,ips,igd,iobs,iobs1,iobs2
      real*8 slpabs,fac,ra,rag
      real*8 cs,ss,csin,ssin,cslw,sslw,csup,ssup
      real*8 g1,g2,gc,gclw,gcup,gcin
      real*8 sdgrad(NPSMAX,2)
c
      call sdmgrad(nps,nobs,sdgrad,wgrad)
      do is=1,ns
        do ips=nps1(is),nps2(is)
          g1=sdgrad(ips,1)
          g2=sdgrad(ips,2)
          slpabs=dsqrt(slpmdl(ips,1)**2+slpmdl(ips,2)**2)
          ra=dmod(datan2(slpmdl(ips,2),slpmdl(ips,1))/DEG2RAD
     &               +rake360(is),360.d0)
c
          if(rake2(is).eq.rake1(is))then
c
c           uniform rake angle
c
            gc=g1*cs1(is)+g2*ss2(is)
            if(slppos(ips).eq.0)then
              gc=dmax1(gc,0.d0)
            else if(slppos(ips).eq.1)then
              gc=dmin1(gc,0.d0)
            endif
            dslpmdl(ips,1)=gc*cs1(is)
            dslpmdl(ips,2)=gc*ss1(is)
          else if(slppos(ips).eq.0)then
c
c           start slip = 0
c
            rag=dmod(datan2(g2,g1)/DEG2RAD+rake360(is),360.d0)
            if(rag.ge.rake1(is).and.rag.le.rake2(is))then
              dslpmdl(ips,1)=g1
              dslpmdl(ips,2)=g2
            else
              gclw=g1*cs1(is)+g2*ss1(is)
              gcup=g1*cs2(is)+g2*ss2(is)
              gc=dmax1(gclw,gcup)
              if(gc.le.0.d0)then
                dslpmdl(ips,1)=0.d0
                dslpmdl(ips,2)=0.d0
			else if(gclw.ge.gcup)then
                dslpmdl(ips,1)=gc*cs1(is)
                dslpmdl(ips,2)=gc*ss1(is)
              else
                dslpmdl(ips,1)=gc*cs2(is)
                dslpmdl(ips,2)=gc*ss2(is)
              endif
            endif
         else if(slppos(ips).eq.1)then
c
c           start slip with the maximum amplitude and with the minimum rake
c
            csin=dcos((ra+90.d0)*DEG2RAD)
            ssin=dsin((ra+90.d0)*DEG2RAD)
            gcin=g1*csin+g2*ssin
            gclw=g1*cs1(is)+g2*ss1(is)
            if(gcin.ge.0.d0.and.gclw.le.0.d0)then
              dslpmdl(ips,1)=g1
              dslpmdl(ips,2)=g2
            else
              gc=dmax1(gcin,-gclw)
              if(gc.le.0.d0)then
                dslpmdl(ips,1)=0.d0
                dslpmdl(ips,2)=0.d0
              else if(gcin.ge.-gclw)then
                dslpmdl(ips,1)=gcin*csin
                dslpmdl(ips,2)=gcin*ssin
              else
                dslpmdl(ips,1)=gclw*cs1(is)
                dslpmdl(ips,2)=gclw*ss1(is)
              endif
            endif
         else if(slppos(ips).eq.2)then
c
c           start slip with the maximum amplitude and with the maximum rake
c
            csin=dcos((ra-90.d0)*DEG2RAD)
            ssin=dsin((ra-90.d0)*DEG2RAD)
            gcin=g1*csin+g2*ssin
            gclw=g1*cs2(is)+g2*ss2(is)
            if(gcin.ge.0.d0.and.gclw.le.0.d0)then
              dslpmdl(ips,1)=g1
              dslpmdl(ips,2)=g2
            else
              gc=dmax1(gcin,-gclw)
              if(gc.le.0.d0)then
                dslpmdl(ips,1)=0.d0
                dslpmdl(ips,2)=0.d0
              else if(gcin.ge.-gclw)then
                dslpmdl(ips,1)=gcin*csin
                dslpmdl(ips,2)=gcin*ssin
              else
                dslpmdl(ips,1)=gclw*cs2(is)
                dslpmdl(ips,2)=gclw*ss2(is)
              endif
            endif
          else if(slppos(ips).eq.3)then
c
c           start slip within the amplitude range and with the minimum rake
c
            csin=dcos((ra+90.d0)*DEG2RAD)
            ssin=dsin((ra+90.d0)*DEG2RAD)
            gcin=g1*csin+g2*ssin
            gclw=g1*cs1(is)+g2*ss1(is)
            if(gcin.ge.0.d0)then
              dslpmdl(ips,1)=g1
              dslpmdl(ips,2)=g2
            else
              dslpmdl(ips,1)=gclw*cs1(is)
              dslpmdl(ips,2)=gclw*ss1(is)
            endif
          else if(slppos(ips).eq.4)then
c
c           start slip within the amplitude range and with the maximum rake
c
            csin=dcos((ra-90.d0)*DEG2RAD)
            ssin=dsin((ra-90.d0)*DEG2RAD)
            gcin=g1*csin+g2*ssin
            gclw=g1*cs2(is)+g2*ss2(is)
            if(gcin.ge.0.d0)then
              dslpmdl(ips,1)=g1
              dslpmdl(ips,2)=g2
            else
              dslpmdl(ips,1)=gclw*cs2(is)
              dslpmdl(ips,2)=gclw*ss2(is)
            endif
          else if(slppos(ips).eq.5)then
c
c           start slip with the maximum amplitude and within the rake range
c
            cs=dcos(ra*DEG2RAD)
            ss=dsin(ra*DEG2RAD)
            if(g1*cs+g2*ss.le.0.d0)then
              dslpmdl(ips,1)=g1
              dslpmdl(ips,2)=g2
            else
              cslw=dcos((ra-90.d0)*DEG2RAD)
              sslw=dsin((ra-90.d0)*DEG2RAD)
              gclw=g1*cslw+g2*sslw
              csup=dcos((ra+90.d0)*DEG2RAD)
              ssup=dsin((ra+90.d0)*DEG2RAD)
              gcup=g1*csup+g2*ssup
              gc=dmax1(gclw,gcup)
              if(gc.le.0.d0)then
                dslpmdl(ips,1)=0.d0
                dslpmdl(ips,2)=0.d0
              else if(gclw.ge.gcup)then
                dslpmdl(ips,1)=gc*cslw
                dslpmdl(ips,2)=gc*sslw
              else
                dslpmdl(ips,1)=gc*csup
                dslpmdl(ips,2)=gc*ssup
              endif
            endif
          else
c
c           start slip within the amplitude and rake range
c
            dslpmdl(ips,1)=g1
            dslpmdl(ips,2)=g2
          endif
        enddo
      enddo
c
      iobs1=0
      do igd=1,ngd
        doffset(igd)=0.d0
        iobs2=iobs1+nobsj(igd)
        if(seloffset(igd).eq.1)then
          do iobs=iobs1+1,iobs2
            doffset(igd)=doffset(igd)+wf(iobs)*dspres(iobs)
          enddo
        endif
        iobs1=iobs2
      enddo
      return
      end
