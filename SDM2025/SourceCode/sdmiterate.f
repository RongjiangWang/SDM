      subroutine sdmiterate(ierr)
      use sdmalloc
      implicit none
      integer*4 ierr
c
c     SDM iteration
c     Last modified: Beijing, July 2026, by R. Wang
c
      integer*4 i,j,ira,ips,igd,is,ipar
      integer*4 irelax,jter,nrelax
      real*8 misfit,corl,vecvar,dvcvar,fmodi
      character*1 text
c
      real*8 sdmcorl
      real*8 relax(nrelaxmax)
      logical*2 convergence
c
      real*8 eps
      data eps/1.0d-12/
c
      sysmis=1.d0
      fmodi=1.d0
      logfile='log_'//slipout
      open(32,file=logfile,status='unknown')
      if(niter.gt.0)then
        write(*,'(a)') '   iter.      cost_function'
     &               //'     iteration_step'
        write(32,'(a)')'   iter.      cost_function'
        write(*, '(i8,2f19.15)')iter,sysmis
        write(32,'(i8,2f19.15)')iter,sysmis
        sysmis0=1.d0
c
c       Modified Landweber iteration
c
        open(20,file=zhyrelax,status='old')
        nrelax=0
        do irelax=1,nrelaxmax
          read(20,*,end=10)relax(irelax)
          nrelax=nrelax+1
        enddo
10      close(20)
c
        irelax=0
        mweq=0.d0
        mwssum=0.d0
        mwpsum=0.d0
c
        jter=0
        do i=1,nsys
          vecswp(i)=sysvec(i)
        enddo
        call sdmproj(ierr)
        call sdmresbat(ierr)
        step=1.d0/sig2max
      endif
c
      do iter=1,niter
c
20      continue
c
        do i=1,nsys
          sysvec(i)=vecswp(i)-step*resbat(i)
        enddo
c
        call sdmproj(ierr)
        call sdmresbat(ierr)
c
        vecvar=0.d0
        dvcvar=0.d0
        do i=1,nsys
          vecvar=vecvar+sysvec(i)**2
          dvcvar=dvcvar+(sysvec(i)-vecswp(i))**2
        enddo
c
c       ckeck convergence
c
        convergence=dabs(sysmis-sysmis0).le.eps*sysmis.and.
     &              dvcvar.le.eps*vecvar
c
        if(sysmis.lt.sysmis0.or.iter.le.1)then
          if(sysmis.le.1.d0)then
            write( *,'(i8,f19.15,f19.6)')iter,sysmis,step*sig2max
            write(32,'(i8,f19.15)')iter,sysmis
          else
            write( *,'(i8,E19.11,f19.6)')iter,sysmis,step*sig2max
            write(32,'(i8,E19.11)')iter,sysmis
          endif
          if(step*sig2max.gt.1000.d0)then
            sig2max=0.5d0*sig2max
            fmodi=fmodi*2.d0
            write(*,'(a,E19.11)')' Landweber step doubled to ',fmodi
          endif
        else
          jter=jter+1
          do i=1,nsys
            sysvec(i)=vecswp(i)
          enddo
          call sdmproj(ierr)
          call sdmresbat(ierr)
          if(step*sig2max.le.10.d0)then
            sig2max=2.d0*sig2max
            fmodi=fmodi*0.5d0
            write(*,'(a,E19.11)')' Landweber step halved to ',fmodi
            if(fmodi.le.eps)then
              convergence=.true.
              goto 50
            endif
          endif
          step=1.d0/sig2max
          goto 20
        endif
c
        irelax=irelax+1
        if(irelax.gt.nrelax)irelax=1
        step=relax(irelax)/sig2max
c
        sysmis0=sysmis
        do i=1,nsys
          vecswp(i)=sysvec(i)
        enddo
c
50      continue
c
        if(convergence)then
          write(*,'(a,i6,a)')' Convergence achieved by ',
     &                  iter,' successful iterations!'
          goto 100
        endif
      enddo
      if(niter.gt.0)then
        write(*,'(a,i6,a)')' Convergence not achieved by ',
     &                  niter,' iterations!'
      endif
100   continue
      if(niter.gt.0)then
        write(*,'(a,i6)')' --- failed interations: ',jter
      endif
c
c     end of Landwerber iteration
c
      if(niter.eq.0)then
c
c       forward model
c
        open(20,file=slipout,status='old')
        read(20,'(a)')text
        do ips=1,nps
          read(20,*)plat(ips),plon(ips),pz(ips),pl(ips),pw(ips),
     &              dlen(ips),dwid(ips),slpmdl(1,ips),slpmdl(2,ips)
          pz(ips)=pz(ips)*KM2M
          pl(ips)=pl(ips)*KM2M
          pw(ips)=pw(ips)*KM2M
          dlen(ips)=dlen(ips)*KM2M
          dwid(ips)=dwid(ips)*KM2M
          slpmdl(2,ips)=-slpmdl(2,ips)
        enddo
        close(20)
        do ipar=1,ndpar
          corrpar(ipar)=0.d0
        enddo
        i=0
        do ips=1,nps
          do ira=1,2
            i=i+1
            sysvec(i)=slpmdl(ira,ips)*zhy(ips)
          enddo
        enddo
c
        do i=1,nsys
          vecswp(i)=sysvec(i)
        enddo
c
        call sdmresbat(ierr)
c
        sysmis=0.d0
        do i=1,nsys
          sysmis=sysmis+sysvec(i)*(resbat(i)-sysbat(i))
        enddo
        sysmis=1+sysmis/datnrm
      endif
      write(*,'(a)') '      Mw data_misfit cost_function mdl_roughness'
      write(32,'(a)')'      Mw data_misfit cost_function mdl_roughness'
      call sdmcalmw(ierr)
      call sdmdatfit(ierr)
      misfit=rmsresall
      write(*, 1000)mwssum,misfit,sysmis,roughness
      write(32,1000)mwssum,misfit,sysmis,roughness
c
      do igd=1,ngd
        write(*,'(i6,a,4f9.4)')igd,'. rmsres/max/min: ',
     &      rmsres(igd),resmax(igd),resmin(igd)
        write(32,'(i6,a,4f9.4)')igd,'. rmsres/max/min: ',
     &      rmsres(igd),resmax(igd),resmin(igd)
      enddo
c
200   continue
c
      corl=sdmcorl(ierr)
      call sdmsdrop(ierr)
      call sdmfinmsg(ierr)
c
c     final message
c
      do ipar=1,ndpar
        write(*,'(a,i2,a,E14.6)')'  data-correction parameter ',ipar,
     &                         ' found: ',corrpar(ipar)
        write(32,'(a,i2,a,E14.6)')'  data-correction parameter ',ipar,
     &                         ' found: ',corrpar(ipar)
      enddo
      write(*,'(a)')'  seg  mean_slp mean_rake   max_slp'
     &      //'      rake   pos_lat   pos_lon     pos_z'
      write(32,'(a)')'  seg  mean_slp mean_rake   max_slp'
     &      //'      rake   pos_lat   pos_lon     pos_z'
      do is=1,ns
        write(*,'(i5,7f10.2)')is,slpm(is),ram(is),slpp(is),rap(is),
     &      plat(ipsp(is)),plon(ipsp(is)),pz(ipsp(is))/KM2M
        write(32,'(i5,7f10.2)')is,slpm(is),ram(is),slpp(is),rap(is),
     &      plat(ipsp(is)),plon(ipsp(is)),pz(ipsp(is))/KM2M
      enddo
      write(*,'(a)')' =================================='
      write(32,'(a)')' =================================='
      write(*,'(3(a,f8.4),a)')' Derived moment magnitude Mw = ',
     &                          mweq,' - ', mwpsum,' (',mwssum,')'
      write(32,'(3(a,f8.4),a)')' Derived moment magnitude Mw = ',
     &                          mweq,' - ',mwpsum,' (',mwssum,')'
c
      write(*,'(a)')' Ave/Std/Max stress change [MPa]: '
      write(32,'(a)')' Ave/Std/Max stress change [MPa]: '
      do is=1,ns
        write(*,'(i4,3f12.4)')is,measdrop(is)/MEGA,
     &      stdsdrop(is)/MEGA,maxsdrop(is)/MEGA
        write(32,'(i4,3f12.4)')is,measdrop(is)/MEGA,
     &      stdsdrop(is)/MEGA,maxsdrop(is)/MEGA
      enddo
c
      write(*,'(a,f8.4)')' Data-model correlation: ',corl
      write(32,'(a,f8.4)')' Data-model correlation: ',corl
      close(32)
1000  format(f8.4,f12.6,2E14.6)
      return
      end