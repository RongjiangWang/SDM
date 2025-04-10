      subroutine sdmout(ns,nps,ngd)
      implicit none
c
c     Last modified: Potsdam, Oct, 2008, by R. Wang
c
      integer ns,nps,ngd
c
      include 'sdmglob.h'
c
      integer i,j,is,ira,ips,igd,iobs,nobs0
      double precision latobs0,lonobs0,dspobs0,slp,sd,rake,st,di
      double precision dinc0,dazi0,xcs0,ycs0,zcs0,cs,ss,sdam,delta0
      double precision premean,obsmean,dspres0,dspres1,dspres2,dspres3
      double precision dsp0(4),ra(2),sm(3,3)
      character*1 text,txtis*2,slipfile*80
c
      integer nobs0max
      data nobs0max/10000000/
c
      write(*,'(a)')' ...output slip model ...'
      open(31,file=slipout,status='unknown')
      write(31,'(a)')'     lat_deg     lon_deg    depth_km'
     &             //'   length_km    width_km'
     &             //'  slp_strk_m  slp_ddip_m    slp_am_m'
     &             //'  strike_deg     dip_deg    rake_deg'
     &             //' sig_stk_MPa sig_ddi_MPa sig_nrm_MPa'
      do is=1,ns
        do ips=nps1(is),nps2(is)
          slp=dsqrt(slpmdl(ips,1)**2+slpmdl(ips,2)**2)
          rake=dmod(datan2(slpmdl(ips,2),slpmdl(ips,1))/DEG2RAD
     &              +rake360(is),360.d0)
          sdam=dsqrt(strdrop(ips,1)**2+strdrop(ips,2)**2)
          write(31,'(14f12.4)')plat(ips),plon(ips),pz(ips)/KM2M,
     &       dlen(ips)/KM2M,dwid(ips)/KM2M,
     &       slpmdl(ips,1),-slpmdl(ips,2),
     &       slp,strike(ips),dip(ips),rake,
     &       strdrop(ips,1)*1.d-6,-strdrop(ips,2)*1.d-6,
     &       strdrop(ips,3)*1.d-6
        enddo
        if(is.lt.ns)write(31,'(a)')'          '
      enddo
      close(31)
      do is=1,ns
        i=is/10
        txtis(1:1)=char(ichar('0')+i)
        i=mod(is,10)
        txtis(2:2)=char(ichar('0')+i)
c
        slipfile='s'//txtis//'_'//slipout
        open(31,file=slipfile,status='unknown')
        write(31,'(a)')'     lat_deg     lon_deg    depth_km'
     &             //'   length_km    width_km'
     &             //'  slp_strk_m  slp_ddip_m    slp_am_m'
     &             //'  strike_deg     dip_deg    rake_deg'
     &             //' sig_stk_MPa sig_ddi_MPa sig_nrm_MPa'
        do ips=nps1(is),nps2(is)
          slp=dsqrt(slpmdl(ips,1)**2+slpmdl(ips,2)**2)
          rake=dmod(datan2(slpmdl(ips,2),slpmdl(ips,1))/DEG2RAD
     &              +rake360(is),360.d0)
          sdam=dsqrt(strdrop(ips,1)**2+strdrop(ips,2)**2)
          write(31,'(14f12.4)')plat(ips),plon(ips),pz(ips)/KM2M,
     &       dlen(ips)/KM2M,dwid(ips)/KM2M,
     &       slpmdl(ips,1),-slpmdl(ips,2),
     &       slp,strike(ips),dip(ips),rake,
     &       strdrop(ips,1)*1.d-6,-strdrop(ips,2)*1.d-6,
     &       strdrop(ips,3)*1.d-6
        enddo
        close(31)
      enddo
c
      do igd=1,ngd
        premean=0.d0
        obsmean=0.d0
        dspres1=0.d0
        dspres2=0.d0
        dspres3=0.d0
        open(30,file=gddata0(igd),status='old')
        do i=1,nheader
          read(30,'(a)')text
        enddo
        open(31,file=gdout0(igd),status='unknown')
        write(31,'(a)')'       latitude      longitude'
     &               //' observation     resiual'
     &               //'  prediction     u_north'
     &               //'      u_east        u_up'
        if(csconst(igd))then
          xcs0=dsin(dinc_const(igd)*DEG2RAD)
     &        *dcos(dazi_const(igd)*DEG2RAD)
          ycs0=dsin(dinc_const(igd)*DEG2RAD)
     &        *dsin(dazi_const(igd)*DEG2RAD)
          zcs0=-dcos(dinc_const(igd)*DEG2RAD)
        endif
        nobs0=0
        do iobs=1,nobs0max
          if(csconst(igd))then
            read(30,*,end=200)latobs0,lonobs0,dspobs0
          else
            read(30,*,end=200)latobs0,lonobs0,dspobs0,delta0,dinc0,dazi0
            xcs0=dsin(dinc0*DEG2RAD)*dcos(dazi0*DEG2RAD)
            ycs0=dsin(dinc0*DEG2RAD)*dsin(dazi0*DEG2RAD)
            zcs0=-dcos(dinc0*DEG2RAD)
          endif
          dspobs0=dspobs0*dspunit
          call sdmfit(ns,xcs0,ycs0,zcs0,latobs0,lonobs0,dsp0)
          dspres0=dspobs0-offset(igd)-dsp0(4)
          write(31,'(2f15.5,6f12.5)')latobs0,lonobs0,
     &        dspobs0/dspunit,dspres0/dspunit,
     &        dsp0(4)/dspunit,dsp0(1)/dspunit,
     &        dsp0(2)/dspunit,-dsp0(3)/dspunit
          premean=premean+dsp0(4)
          obsmean=obsmean+dspobs0
          dspres1=dspres1+dspres0**2
          dspres2=dmax1(dspres2,dspres0)
          dspres3=dmax1(dspres3,-dspres0)
          nobs0=nobs0+1
        enddo
        print *,' Warning in eqsout: nobs0max defined too small!'
200     close(30)
        close(31)
        premean=premean/dble(nobs0)
        obsmean=obsmean/dble(nobs0)
        dspres1=dsqrt(dspres1/dble(nobs0))
        write(*,'(a,i3)')' Data set: ',igd
        write(*,'(a,2f10.3)') ' ...predicted/observed mean values[m]: ',
     &                        premean,obsmean
        write(*,'(a,3f10.3)')' ...msr/max/min[m] residuals: ',
     &                            dspres1,dspres2,-dspres3
      enddo
c
1001  format(i3,a,f8.2,a,f8.2,a,f8.2,2f12.2)
      return
      end