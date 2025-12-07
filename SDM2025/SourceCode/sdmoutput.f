      subroutine sdmoutput(ierr)
      use sdmalloc
      implicit none
      integer*4 ierr
c
c     Last modified: 7 Dec. 2025 by R. Wang
c
      integer*4 i,j,is,ira,ips,jps,igd,iobs,ipar
      real*8 slp,rake,sdam
      character*1 text,txtis*2,slipfile*80
c
c     calculate elastic stress changes
c
      do jps=1,nps
        do i=1,3
          strdrop(i,jps)=0.d0
          do ips=1,nps
            do ira=1,2
              strdrop(i,jps)=strdrop(i,jps)
     &                      +slpmdl(ira,ips)*strgrn(ira,ips,i,jps)
            enddo
          enddo
        enddo
      enddo
c
      open(31,file=slipout,status='unknown')
      write(31,'(a)')'     lat_deg     lon_deg    depth_km'
     &             //'  x_local_km  y_local_km'
     &             //'   length_km    width_km'
     &             //'  slp_strk_m  slp_ddip_m    slp_am_m'
     &             //'  strike_deg     dip_deg    rake_deg'
     &             //' sig_stk_MPa sig_ddi_MPa sig_nrm_MPa'
      do is=1,ns
        i=is/10
        txtis(1:1)=char(ichar('0')+i)
        i=mod(is,10)
        txtis(2:2)=char(ichar('0')+i)
c
        slipfile='s'//txtis//'_'//slipout
        open(32,file=slipfile,status='unknown')
        write(32,'(a)')'     lat_deg     lon_deg    depth_km'
     &             //'  x_local_km  y_local_km'
     &             //'   length_km    width_km'
     &             //'  slp_strk_m  slp_ddip_m    slp_am_m'
     &             //'  strike_deg     dip_deg    rake_deg'
     &             //' sig_stk_MPa sig_ddi_MPa sig_nrm_MPa'
        do ips=nps1(is),nps2(is)
          slp=dsqrt(slpmdl(1,ips)**2+slpmdl(2,ips)**2)
          rake=dmod(datan2(slpmdl(2,ips),slpmdl(1,ips))/DEG2RAD
     &              +rake360(is),360.d0)
          sdam=dsqrt(strdrop(1,ips)**2+strdrop(2,ips)**2)
          write(31,'(16f12.4)')plat(ips),plon(ips),
     &       pz(ips)/KM2M,pl(ips)/KM2M,pw(ips)/KM2M,
     &       dlen(ips)/KM2M,dwid(ips)/KM2M,
     &       slpmdl(1,ips),-slpmdl(2,ips),
     &       slp,strike(ips),dip(ips),rake,
     &       strdrop(1,ips)/MEGA,-strdrop(2,ips)/MEGA,
     &       strdrop(3,ips)/MEGA
          write(32,'(16f12.4)')plat(ips),plon(ips),
     &       pz(ips)/KM2M,pl(ips)/KM2M,pw(ips)/KM2M,
     &       dlen(ips)/KM2M,dwid(ips)/KM2M,
     &       slpmdl(1,ips),-slpmdl(2,ips),
     &       slp,strike(ips),dip(ips),rake,
     &       strdrop(1,ips)/MEGA,-strdrop(2,ips)/MEGA,
     &       strdrop(3,ips)/MEGA
        enddo
        close(32)
        if(is.lt.ns)write(31,'(a)')'          '
      enddo
      close(31)
c
      do igd=1,ngd
        open(33,file=gdout(igd),status='unknown')
        write(33,'(a)')'       latitude      longitude'
     &               //' observation     resiual'
     &               //'  prediction  correction'
        do iobs=nobs1(igd),nobs2(igd)
          write(33,'(2f15.5,4f12.5)')latobs(iobs),lonobs(iobs),
     &        datobs(iobs)/datunit,(datobs(iobs)-datmdl(iobs))/datunit,
     &        datmdl(iobs)/datunit,corrmdl(iobs)/datunit
        enddo
        close(33)
      enddo
c
      if(npar.gt.0)then
        open(34,file=parout,status='unknown')
        write(34,'(a)')'  no data_correction_parameter'
        do ipar=1,npar
          write(34,'(i4,E16.6)')ipar,corrpar(ipar)
        enddo
        close(34)
      endif
      return
      end