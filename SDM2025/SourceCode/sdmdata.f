      subroutine sdmdata(ierr)
      use sdmalloc
      implicit none
      integer*4 ierr
c
c     read in geodetic data
c
      integer*4 i,iobs,iobs0,igd,unit
      real*8 lat,lon,obs,d4,d5,d6,obsmn,wfsum
      character*1 text
      character*180 dataline
c
      unit=21
c
c     test data format
c
      do igd=1,ngd
        open(unit,file=gddata(igd),status='old')
        do i=1,nheader
          read(unit,'(a)')text
        enddo
        if(csconst(igd))then
          read(unit,'(a)')dataline
          read(dataline,*,end=101)lat,lon,obs,d4
          goto 102
        else
          read(unit,'(a)')dataline
          read(dataline,*,end=101)lat,lon,obs,d4,d5,d6
          goto 102
        endif
101     close(unit)
        write(*,'(a,i2)')' Wrong format of dataset ',igd
        stop
102     close(unit)
      enddo
c
c     count the number of data
c
      nobs=0
      do igd=1,ngd
        open(unit,file=gddata(igd),status='old')
        nobsj(igd)=0
        do i=1,nheader
          read(unit,'(a)')text
        enddo
        if(csconst(igd))then
103       continue
          read(unit,*,end=105)lat,lon,obs,d4
          nobsj(igd)=nobsj(igd)+1
          goto 103
        else
104       continue
          read(unit,'(a)')dataline
          read(unit,*,end=105)lat,lon,obs,d4,d5,d6
          nobsj(igd)=nobsj(igd)+1
          goto 104
        endif
105     close(unit)
        nobs=nobs+nobsj(igd)
      enddo
c
      allocate(latobs(nobs),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmdata: latobs not allocated!'
      allocate(lonobs(nobs),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmdata: lonobs not allocated!'
      allocate(datobs(nobs),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmdata: datobs not allocated!'
      allocate(wf(nobs),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmdata: wf not allocated!'
      allocate(xcs(nobs),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmdata: xcs not allocated!'
      allocate(ycs(nobs),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmdata: ycs not allocated!'
      allocate(zcs(nobs),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmdata: zcs not allocated!'
c
c     read data
c
      iobs0=0
      do igd=1,ngd
        open(unit,file=gddata(igd),status='old')
        do i=1,nheader
          read(21,'(a)')text
        enddo
        obsmn=0.d0
        wfsum=0.d0
        do iobs=iobs0+1,iobs0+nobsj(igd)
          if(csconst(igd))then
            read(unit,*)lat,lon,obs,d4
            d5=dinc_const(igd)*DEG2RAD
            d6=dazi_const(igd)*DEG2RAD
          else
            read(unit,*)lat,lon,obs,d4,d5,d6
            d5=d5*DEG2RAD
            d6=d6*DEG2RAD
          endif
          xcs(iobs)=dsin(d5)*dcos(d6)
          ycs(iobs)=dsin(d5)*dsin(d6)
          zcs(iobs)=-dcos(d5)
          latobs(iobs)=lat
          lonobs(iobs)=lon
          if(d4.le.0.d0)then
            stop 'Error in sdmdata: measurement error <= 0!'
          endif
          wf(iobs)=1.d0/d4**2
          wfsum=wfsum+wf(iobs)
          datobs(iobs)=obs*datunit
          obsmn=obsmn+wf(iobs)*datobs(iobs)
        enddo
        close(unit)
        obsmn=obsmn/wfsum
        do iobs=iobs0+1,iobs0+nobsj(igd)
          wf(iobs)=dsqrt(wfm(igd)*wf(iobs)/wfsum)
        enddo
        iobs0=iobs0+nobsj(igd)
c
        write(*,'(i7,a,i2,a,f10.4)')nobsj(igd),
     +       ' data read from data set: ',
     +       igd,', mean value: ',obsmn
c
        close(unit)
      enddo
c
      nobs1(1)=1
      nobs2(1)=nobsj(1)
      do igd=2,ngd
        nobs1(igd)=nobs2(igd-1)+1
        nobs2(igd)=nobs1(igd)+nobsj(igd)-1
      enddo
c
      wfsum=0.d0
      do iobs=1,nobs
        wfsum=wfsum+wf(iobs)**2
      enddo
c
      do iobs=1,nobs
        wf(iobs)=wf(iobs)/dsqrt(wfsum)
      enddo
c
      wfsum=0.d0
      rmsdatall=0.d0
      do iobs=1,nobs
        rmsdatall=rmsdatall+(wf(iobs)*datobs(iobs))**2
        wfsum=wfsum+wf(iobs)**2
      enddo
      rmsdatall=dsqrt(rmsdatall/wfsum)
c
      return
      end
