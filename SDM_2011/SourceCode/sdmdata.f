      subroutine sdmdata(ngd,nobs)
      implicit none
c
c     read in geodetic data
c
      integer ngd,nobs
c
      include 'sdmglob.h'
c
      integer i,igd
      double precision lat,lon,obs,d4,d5,d6,obsmn
      character*1 text
      character*180 dataline
c
      nobs=0
      wfsum=0.d0
      do igd=1,ngd
        open(21,file=gddata(igd),status='old')
        do i=1,nheader
          read(21,'(a)')text
        enddo
c
c       test data format
c
        if(csconst(igd))then
          call getdata(21,dataline)
          read(dataline,*,end=100)lat,lon,obs,d4
          goto 101
        else
          call getdata(21,dataline)
          read(dataline,*,end=100)lat,lon,obs,d4,d5,d6
          goto 101
        endif
100     close(21)
        write(*,'(a,i2)')' Wrong format of dataset ',igd
        stop
101     close(21)
        open(21,file=gddata(igd),status='old')
        do i=1,nheader
          read(21,'(a)')text
        enddo
        nobsj(igd)=0
        obsmn=0.d0
        wfmsum(igd)=0.d0
        do i=1,nobsmax
          if(csconst(igd))then
            read(21,*,end=200)lat,lon,obs,d4
            if(d4.le.0.d0)then
              stop 'Error in sdmdata: measurement error <= 0!'
            endif
            d5=dinc_const(igd)*DEG2RAD
            d6=dazi_const(igd)*DEG2RAD
          else
            read(21,*,end=200)lat,lon,obs,d4,d5,d6
            if(d4.le.0.d0)then
              stop 'Error in sdmdata: measurement error <= 0!'
            endif
            d5=d5*DEG2RAD
            d6=d6*DEG2RAD
          endif
          nobs=nobs+1
          nobsj(igd)=nobsj(igd)+1
          if(nobs.gt.nobsmax)then
            stop ' Error in sdmdata: nobsmax defined too small!'
          endif
          xcs(nobs)=dsin(d5)*dcos(d6)
          ycs(nobs)=dsin(d5)*dsin(d6)
          zcs(nobs)=-dcos(d5)
          latobs(nobs)=lat
          lonobs(nobs)=lon
          wf(nobs)=wfm(igd)/d4**2
          wfmsum(igd)=wfmsum(igd)+wf(nobs)
          dspobs(nobs)=obs*dspunit
          obsmn=obsmn+dspobs(nobs)
        enddo
        read(21,*,end=200)lat,lon,obs
        stop ' Error in sdmdata: nobsmax too small!'
200     close(21)
        wfsum=wfsum+wfmsum(igd)
        obsmn=obsmn/dble(nobsj(igd))
        write(*,'(i7,a,i2,a,f10.4)')nobsj(igd),
     +       ' data read from (filtered) data set: ',
     +       igd,', mean value: ',obsmn
      enddo
      do i=1,nobs
        wf(i)=wf(i)/wfsum
      enddo
      do igd=1,ngd
        wfmsum(igd)=wfmsum(igd)/wfsum
      enddo
      wfsum=1.d0
      return
      end
