      subroutine sdminv(ngd,ns,nps,nobs,niter,wg0,ismooth)
      implicit none
c
c     Last modified: Potsdam, Oct, 2008, by R. Wang
c
      integer*4 ngd,ns,nps,nobs,niter,ismooth
      real*8 wg0
c
      include 'sdmglob.h'
c
      integer*4 i,j,is,igd,iobs,iobs1,iobs2,ips,ira,iter
      integer*4 ipsmax(NSMAX)
      real*8 ra,rac,st,di,sdam,slp,rake,farea
      real*8 wgrad,misfit,roughness
      real*8 meanslip,slpabs,ms,mw,me,mep,fac,corr
      real*8 rmsres,cost,ocost,scost
      real*8 sm(3,3),smp(3,3),swap(100)
      real*8 ssm(NSMAX),dsm(NSMAX),slpm(NSMAX),slpmax(NSMAX)
      real*8 maxres(NGDMAX),minres(NGDMAX),meanres(NGDMAX)
      real*8 swapslp(NPSMAX,2),swapoffs(NGDMAX)
      character logfile*80,text*80
      real*8 sdmocost,sdmscost,sdmcorr,sdmsmod
      logical*2 converge
c
      integer*4 icf,ncf
      parameter(ncf=10)
      real*8 cf1,cf2,cf(2*ncf)
c
      real*8 eps
      data eps/1.0d-06/
c
      wgrad=0.d0
      if(niter.gt.0)then
        do ips=1,nps
          slpmdl(ips,1)=0.d0
          slpmdl(ips,2)=0.d0
          swapslp(ips,1)=0.d0
          swapslp(ips,2)=0.d0
          slppos(ips)=0
        enddo
c
        do icf=1,2*ncf
          cf(icf)=0.d0
        enddo
      else
        open(20,file=slipout,status='old')
        read(20,'(a)')text
        do ips=1,nps
          read(20,*)(swap(i),i=1,9)
          slpmdl(ips,1)=swap(8)
          slpmdl(ips,2)=-swap(9)
          swapslp(ips,1)=slpmdl(ips,1)
          swapslp(ips,2)=slpmdl(ips,2)
        enddo
        close(20)
      endif
c
      do igd=1,ngd
        offset(igd)=0.d0
        swapoffs(igd)=0.d0
      enddo
c
      do iobs=1,nobs
        dspres(iobs)=dspobs(iobs)
      enddo
c
      call sdmddsp(nobs,nps,swapslp)
      call sdmres(ngd,nps,swapoffs,rmsres,maxres,minres,meanres)
      costref=sdmocost(ngd,nobs,swapoffs)
c
      mw=0.d0
c
      iter=0
      open(30,file='converge.dat',status='unknown')
      write(30,'(a)')' iter.       Mw      misfit   roughness'
      logfile='log_'//slipout
      open(32,file=logfile,status='unknown')
      write(*,'(a)') ' iter.       Mw      misfit   roughness'
      write(32,'(a)')' iter.       Mw      misfit   roughness'
      misfit=1.d0
      roughness=0.d0
      write(30,1000)iter,mw,misfit,roughness
      write(*, 1000)iter,mw,misfit,roughness
      write(32,1000)iter,mw,misfit,roughness
      do igd=1,ngd
        write(*,'(i6,a,4f9.4)')igd,'. rmsres/max/min/offset: ',
     &        meanres(igd),maxres(igd),minres(igd),offset(igd)
        write(32,'(i6,a,4f9.4)')igd,'. rmsres/max/min/offset: ',
     &        meanres(igd),maxres(igd),minres(igd),offset(igd)
      enddo
c
      do iter=1,niter
        call sdmdslip(ngd,nobs,ns,nps,wgrad,swapslp,swapoffs)
        call sdmscale(ngd,nobs,ns,nps,wgrad,swapslp,swapoffs,converge)
c
c       evaluate the new slip model (max. min. and rms residuals)
c
        call sdmddsp(nobs,nps,swapslp)
        call sdmres(ngd,nps,swapoffs,rmsres,maxres,minres,meanres)
        do ips=1,nps
          slpmdl(ips,1)=slpmdl(ips,1)+swapslp(ips,1)
          slpmdl(ips,2)=slpmdl(ips,2)+swapslp(ips,2)
          swapslp(ips,1)=slpmdl(ips,1)
          swapslp(ips,2)=slpmdl(ips,2)
        enddo
c
        do igd=1,ngd
          offset(igd)=offset(igd)+swapoffs(igd)
          swapoffs(igd)=0.d0
        enddo
c
        if(converge)then
          write(*,'(a)')' Convergence arrived!'
          goto 200
        endif
c
c       statistics
c
        do i=1,3
          do j=1,3
            sm(i,j)=0.d0
          enddo
        enddo
c
        ms=0.d0
        me=0.d0
        do is=1,ns
          do ips=nps1(is),nps2(is)
            st=strike(ips)*DEG2RAD
            di=dip(ips)*DEG2RAD
            ra=datan2(slpmdl(ips,2),slpmdl(ips,1))
            fac=parea(ips)*dsqrt(slpmdl(ips,2)**2+slpmdl(ips,1)**2)
            if(hsmodel)then
              fac=fac*3.d+10
            else
              fac=fac*muz(iz(ips))
            endif
            ms=ms+fac
c
            smp(1,1)=-dsin(di)*dcos(ra)*dsin(2.d0*st)
     &              -dsin(2.d0*di)*dsin(ra)*(dsin(st))**2
            smp(1,2)=dsin(di)*dcos(ra)*dcos(2.d0*st)
     &              +0.5d0*dsin(2.d0*di)*dsin(ra)*dsin(2.d0*st)
            smp(1,3)=-dcos(di)*dcos(ra)*dcos(st)
     &             -dcos(2.d0*di)*dsin(ra)*dsin(st)
            smp(2,2)=dsin(di)*dcos(ra)*sin(2.d0*st)
     &              -dsin(2.d0*di)*dsin(ra)*(dcos(st))**2
            smp(2,3)=-dcos(di)*dcos(ra)*dsin(st)
     &              +dcos(2.d0*di)*dsin(ra)*dcos(st)
            smp(3,3)=dsin(2.d0*di)*dsin(ra)
c
            do i=1,3
              do j=i,3
                sm(i,j)=sm(i,j)+smp(i,j)*fac
              enddo
            enddo
c
            mep=0.d0
            do i=1,3
              mep=mep+0.5d0*smp(i,i)**2
              do j=i+1,3
                mep=mep+smp(i,j)**2
              enddo
            enddo
            me=me+dsqrt(mep)*fac
          enddo
        enddo
        if(me.gt.0.d0)then
          me=(dlog10(me)-9.1d0)/1.5d0
        else
          me=0.d0
        endif
        if(ms.gt.0.d0)then
          ms=(dlog10(ms)-9.1d0)/1.5d0
        else
          ms=0.d0
        endif
c
        mw=0.d0
        do i=1,3
          mw=mw+0.5d0*sm(i,i)**2
          do j=i+1,3
            mw=mw+sm(i,j)**2
          enddo
        enddo
        if(mw.gt.0.d0)then
          mw=(dlog10(dsqrt(mw))-9.1d0)/1.5d0
        else
          mw=0.d0
        endif
c
        ocost=0.5d0*wfsum*rmsres**2
        misfit=dsqrt(ocost/costref)
        scost=sdmscost(nps,swapslp)
        if(scostref.gt.0.d0.and.wgrad.le.0.d0)then
c
c         fixing normalized smoothing factor
c         scostref = roughness after the first iteration
c         costref = data variance
c
          wgrad=wg0*costref/scostref
        endif
        roughness=scost/sdmsmod(ns,nps,ismooth)
        write(30,1000)iter,mw,misfit,roughness
        write(*, 1000)iter,mw,misfit,roughness
        write(32,1000)iter,mw,misfit,roughness
        do igd=1,ngd
          write(*,'(i6,a,4f9.4)')igd,'. rmsres/max/min/offset: ',
     &          meanres(igd),maxres(igd),minres(igd),offset(igd)
          write(32,'(i6,a,4f9.4)')igd,'. rmsres/max/min/offset: ',
     &          meanres(igd),maxres(igd),minres(igd),offset(igd)
        enddo
c
        cost=ocost+wgrad*scost
c
        cf1=cf(1)
        do icf=1,ncf
          cf(icf)=cf(icf+1)
          cf1=dmin1(cf1,cf(icf))
        enddo
        cf2=cost
        do icf=ncf+1,2*ncf-1
          cf(icf)=cf(icf+1)
          cf2=dmin1(cf2,cf(icf))
        enddo
        cf(2*ncf)=cost
c
        if(cost.eq.cf2.and.dabs(cf2-cf1).le.eps*cf2)then
          write(*,'(a)')' Convergence arrived!'
          goto 200
        endif
      enddo
200   continue
c
      do is=1,ns
        ssm(is)=0.d0
        dsm(is)=0.d0
        slpm(is)=0.d0
        slpmax(is)=0.d0
        ipsmax(is)=nps1(is)
        farea=0.d0
        do ips=nps1(is),nps2(is)
          ssm(is)=ssm(is)+slpmdl(ips,1)*parea(ips)
          dsm(is)=dsm(is)+slpmdl(ips,2)*parea(ips)
          slpabs=dsqrt(slpmdl(ips,1)**2+slpmdl(ips,2)**2)
          slpm(is)=slpm(is)+slpabs*parea(ips)
          if(slpmax(is).lt.slpabs)then
            slpmax(is)=slpabs
            ipsmax(is)=ips
          endif
          farea=farea+parea(ips)
        enddo
        slpm(is)=slpm(is)/farea
      enddo
c
      do igd=1,ngd
        if(seloffset(igd).eq.1)then
          write(*,'(a,i2,a,f12.6)')'  offset found for ',igd,
     &                           '. data set: ',offset(igd)
          write(32,'(a,i2,a,f12.6)')'  offset found for ',igd,
     &                           '. data set: ',offset(igd)
        endif
      enddo
      write(*,'(a)')'  seg  mean_slp mean_rake   max_slp'
     &      //'      rake   pos_lat   pos_lon     pos_z'
      write(32,'(a)')'  seg  mean_slp mean_rake   max_slp'
     &      //'      rake   pos_lat   pos_lon     pos_z'
      do is=1,ns
        ramean(is)=dmod(datan2(dsm(is),ssm(is))/DEG2RAD
     &            +rake360(is),360.d0)
        rac=dmod(datan2(slpmdl(ipsmax(is),2),
     &             slpmdl(ipsmax(is),1))/DEG2RAD+rake360(is),360.d0)
        write(*,'(i5,7f10.2)')is,slpm(is),ramean(is),slpmax(is),rac,
     &      plat(ipsmax(is)),plon(ipsmax(is)),pz(ipsmax(is))/KM2M
        write(32,'(i5,7f10.2)')is,slpm(is),ramean(is),slpmax(is),rac,
     &      plat(ipsmax(is)),plon(ipsmax(is)),pz(ipsmax(is))/KM2M
      enddo
      write(*,'(a)')' =================================='
      write(32,'(a)')' =================================='
      write(*,'(3(a,f8.4),a)')' Derived moment magnitude Mw = ',
     &                          mw,' - ', me,' (',ms,')'
      write(32,'(3(a,f8.4),a)')' Derived moment magnitude Mw = ',
     &                           mw,' - ',me,' (',ms,')'
c
      call sdmsdrop(ns,nps)
c
      write(*,'(a)')' Ave/Std/Max stress drop [MPa]: '
      write(32,'(a)')' Ave/Std/Max stress drop [MPa]: '
      do is=1,ns
        write(*,'(i4,3f12.4)')is,measdrop(is)*1.0d-06,
     &      stdsdrop(is)*1.0d-06,maxsdrop(is)*1.0d-06
        write(32,'(i4,3f12.4)')is,measdrop(is)*1.0d-06,
     &      stdsdrop(is)*1.0d-06,maxsdrop(is)*1.0d-06
      enddo
      corr=sdmcorr(ngd,nobs,nps)
      write(*,'(a,f8.4)')' Data-model correlation: ',corr
      write(32,'(a,f8.4)')' Data-model correlation: ',corr
      close(30)
      close(32)
1000  format(i5,f10.4,f12.6,E14.6)
c
      return
      end
