      subroutine sdmiterate(ierr)
      use sdmalloc
      implicit none
      integer*4 ierr
c
c     SDM iteration
c     Last modified: Zhuhai, Nov. 2025, by R. Wang
c
      integer*4 i,j,ira,ips,igd,is,ipar
      integer*4 irelax,jter,nrelax
      real*8 misfit,corl,a,b,c
      character*1 text
c
      real*8 sdmcorl
      real*8 relax(nrelaxmax)
      logical*2 convergence,landweber
c
      real*8 eps
      data eps/1.0d-12/
c
      iter=0
      sysmis=1.d0
      open(30,file='converge.dat',status='unknown')
      logfile='log_'//slipout
      open(32,file=logfile,status='unknown')
      write(*,'(a)') '   iter.      cost_function'
      write(30,'(a)')'   iter.      cost_function'
      write(32,'(a)')'   iter.      cost_function'
      write(*, '(i8,f19.15)')iter,sysmis
      write(30,'(i8,f19.15)')iter,sysmis
      write(32,'(i8,f19.15)')iter,sysmis
      if(niter.gt.0)then
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
        step=step0
c
        landweber=.true.
        jter=0
        do i=1,nsys
          vecswp(i)=sysvec(i)
        enddo
        do i=1,nsys
          resbat(i)=-sysbat(i)
          do j=1,nsys
            resbat(i)=resbat(i)+sysmat(i,j)*vecswp(j)
          enddo
        enddo
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
c
        sysmis=0.d0
        do i=1,nsys
          resbat(i)=-sysbat(i)
          do j=1,nsys
            resbat(i)=resbat(i)+sysmat(i,j)*sysvec(j)
          enddo
          sysmis=sysmis+sysvec(i)*(resbat(i)-sysbat(i))
        enddo
        sysmis=1+sysmis/datnrm
c
c       ckeck convergence
c
        convergence=dabs(sysmis-sysmis0).le.eps*sysmis
c
        if(sysmis.le.sysmis0)then
          write(*, '(i8,f19.15)')iter,sysmis
          write(30,'(i8,f19.15)')iter,sysmis
          write(32,'(i8,f19.15)')iter,sysmis
        endif
c
        if(sysmis.gt.sysmis0.and..not.landweber)then
          step=1.d0
          landweber=.true.
          jter=jter+1
          do i=1,nsys
            resbat(i)=-sysbat(i)
            do j=1,nsys
              resbat(i)=resbat(i)+sysmat(i,j)*vecswp(j)
            enddo
          enddo
          goto 20
        else
          irelax=irelax+1
          if(irelax.gt.nrelax)irelax=1
          step=relax(irelax)
          landweber=.false.
        endif
        sysmis0=sysmis
        do i=1,nsys
          vecswp(i)=sysvec(i)
        enddo
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
          slpmdl(2,ips)=-slpmdl(2,ips)
        enddo
        close(20)
        call sdmcalmw(ierr)
        call sdmdatfit(ierr)
c
        do ipar=1,npar
          corrpar(ipar)=0.d0
        enddo
      endif
      write(*,'(a)') '      Mw data_misfit cost_function mdl_roughness'
      write(30,'(a)')'      Mw data_misfit cost_function mdl_roughness'
      write(32,'(a)')'      Mw data_misfit cost_function mdl_roughness'
      call sdmcalmw(ierr)
      call sdmdatfit(ierr)
      misfit=rmsresall
      write(*, 1000)mwssum,misfit,sysmis,roughness
      write(30,1000)mwssum,misfit,sysmis,roughness
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
      do ipar=1,npar
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
      write(*,'(a)')' Ave/Std/Max stress drop [MPa]: '
      write(32,'(a)')' Ave/Std/Max stress drop [MPa]: '
      do is=1,ns
        write(*,'(i4,3f12.4)')is,measdrop(is)/MEGA,
     &      stdsdrop(is)*1.0d-06,maxsdrop(is)/MEGA
        write(32,'(i4,3f12.4)')is,measdrop(is)/MEGA,
     &      stdsdrop(is)*1.0d-06,maxsdrop(is)/MEGA
      enddo
c
      write(*,'(a,f8.4)')' Data-model correlation: ',corl
      write(32,'(a,f8.4)')' Data-model correlation: ',corl
      close(30)
      close(32)
1000  format(f8.4,f12.6,2E14.6)
      return
      end