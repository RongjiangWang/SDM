      subroutine sdmgetinp(ierr)
      use sdmalloc
      implicit none
c
c     Last modified: 2025 in Zhuhai by R. Wang
c
      integer*4 ierr
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     LOCAL WORK SPACES
c     =================
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer*4 i,j,is,igd,unit
      real*8 wfsum
      real*8 dswap(10)
      logical*2 testread
      character*180 dataline
c
      unit=10
c
      testread=.true.
      nftmax=0
10    open(unit,file=infile,status='old')
      call skipdoc(unit)
      read(unit,*)iearth
      if(iearth.eq.0)then
        hsmodel=.true.
        call skipdoc(unit)
        read(unit,*)poisson
        if(poisson.lt.0.d0.or.poisson.ge.0.5d0)then
          stop ' Error in sdmgetinp: Wrong Poisson ratio!'
        endif
        zobs=0.d0
        muehs=MUEREF
        lamhs=MUEREF*2.d0*poisson/(1.d0-2.d0*poisson)
      else if(iearth.eq.1)then
        hsmodel=.false.
        muehs=MUEREF
        lamhs=MUEREF
        call skipdoc(unit)
        read(unit,*)grndir,(green(i),i=1,3)
      else if(iearth.eq.2)then
        hsmodel=.true.
        call skipdoc(unit)
        read(unit,*)grndir,usr3dgrn(1),usr3dgrn(2)
        poisson=0.25d0
        muehs=MUEREF
        lamhs=MUEREF*2.d0*poisson/(1.d0-2.d0*poisson)
      else
        stop ' Error in sdmgetinp: Wrong selection of earth model!'
      endif
c00000000000000000000000000000000000000000000000000000000000000000000000
c     READ IN PARAMETERS FOR RECTANGULAR SOURCES
c     ==========================================
c00000000000000000000000000000000000000000000000000000000000000000000000
      call skipdoc(unit)
      read(unit,*)idisc
      if(idisc.lt.0.or.idisc.gt.1)then
        print *,'Error in sdminp: bad switch for subfaults input form'
        stop
      endif
      call skipdoc(unit)
      read(unit,*)ns
      if(testread)then
        allocate(rake1(ns),stat=ierr)
        if(ierr.ne.0)stop ' Error in sdmgetinp: rake1 not allocated!'
        allocate(rake2(ns),stat=ierr)
        if(ierr.ne.0)stop ' Error in sdmgetinp: rake2 not allocated!'
        allocate(maxslip(ns),stat=ierr)
        if(ierr.ne.0)stop ' Error in sdmgetinp: maxslip not allocated!'
        allocate(inpatches(ns),stat=ierr)
        if(ierr.ne.0)stop ' Error in sdmgetinp: inpatches'
     &                  //'not allocated!'
        allocate(isref(ns),stat=ierr)
        if(ierr.ne.0)stop ' Error in sdmgetinp: isref not allocated!'
        allocate(topdep(ns),stat=ierr)
        if(ierr.ne.0)stop ' Error in sdmgetinp: topdep not allocated!'
        allocate(width(ns),stat=ierr)
        if(ierr.ne.0)stop ' Error in sdmgetinp: width not allocated!'
        allocate(mstrike(ns),stat=ierr)
        if(ierr.ne.0)stop ' Error in sdmgetinp: mstrike not allocated!'
        allocate(patchsize(ns),stat=ierr)
        if(ierr.ne.0)then
          stop ' Error in sdmgetinp: patchsize not allocated!'
        endif
        allocate(nft(ns),stat=ierr)
        if(ierr.ne.0)stop ' Error in sdmgetinp: nft not allocated!'
        allocate(cs1(ns),stat=ierr)
        if(ierr.ne.0)stop ' Error in sdmgetinp: cs1 not allocated!'
        allocate(ss1(ns),stat=ierr)
        if(ierr.ne.0)stop ' Error in sdmgetinp: ss1 not allocated!'
        allocate(cs2(ns),stat=ierr)
        if(ierr.ne.0)stop ' Error in sdmgetinp: cs2 not allocated!'
        allocate(ss2(ns),stat=ierr)
        if(ierr.ne.0)stop ' Error in sdmgetinp: ss2 not allocated!'
        allocate(rake360(ns),stat=ierr)
        if(ierr.ne.0)stop ' Error in sdmgetinp: rake360 not allocated!'
      endif
      if(idisc.eq.1)goto 20
      do is=1,ns
        call skipdoc(unit)
        read(unit,*)i,rake1(is),rake2(is),
     &              maxslip(is),inpatches(is),isref(is)
      enddo
      if(idisc.eq.0)goto 30
20    continue
      do is=1,ns
        call skipdoc(unit)
        read(unit,*)i,topdep(is),width(is),
     &                  mstrike(is),patchsize(is)
        call skipdoc(unit)
        read(unit,*)rake1(is),rake2(is),maxslip(is)
        if(rake1(is).gt.rake2(is).or.
     &     rake1(is).lt.-180.d0.or.rake2(is).gt.360.d0)then
          print *, 'Error in sdmgetinp:'
     &           //'bad rake range for fault segment ',is
          stop
        endif
        if(maxslip(is).le.0.d0)then
          print *, 'Error in sdmgetinp:'
     &           //'bad maximum slip for fault segment ',is
          stop
        endif
        call skipdoc(unit)
        read(unit,*)nft(is)
c
        if(testread)nftmax=max0(nftmax,nft(is))
c
        if(.not.testread)then
          allocate(latft(nftmax,ns),stat=ierr)
          if(ierr.ne.0)stop ' Error in sdmgetinp: latft not allocated!'
          allocate(lonft(nftmax,ns),stat=ierr)
          if(ierr.ne.0)stop ' Error in sdmgetinp: lonft not allocated!'
          allocate(topdip(nftmax,ns),stat=ierr)
          if(ierr.ne.0)stop ' Error in sdmgetinp: topdip not allocated!'
          allocate(botdip(nftmax,ns),stat=ierr)
          if(ierr.ne.0)stop ' Error in sdmgetinp: botdip not allocated!'
        endif
c
        do i=1,nft(is)
          call skipdoc(unit)
          read(unit,*)(dswap(j),j=1,4)
          if(dswap(3).le.0.d0.or.dswap(3).ge.180.d0.or.
     &       dswap(4).le.0.d0.or.dswap(4).ge.180.d0)then
            print *,is,'. fault segment: bad dip angle!'
            stop
          endif
          if(.not.testread)then
            latft(i,is)=dswap(1)
            lonft(i,is)=dswap(2)
            topdip(i,is)=dswap(3)
            botdip(i,is)=dswap(4)
          endif
        enddo
c
        if(testread)then
          testread=.false.
          close(unit)
          goto 10
        endif
c
        topdep(is)=KM2M*topdep(is)
        width(is)=KM2M*width(is)
        patchsize(is)=KM2M*patchsize(is)
      enddo
30    continue
      do is=1,ns
        if(maxslip(is).le.0.d0)then
          stop ' Bad parameter max_slip!'
        endif
c
        if(rake1(is).lt.-180.d0.or.rake1(is).gt.rake2(is))then
          stop ' Wrong rake range!'
        endif
c
        if(rake1(is).ge.0.d0.and.rake2(is).ge.0.d0)then
          rake360(is)=360.d0
        else
          rake360(is)=0.d0
        endif
        cs1(is)=dcos(rake1(is)*DEG2RAD)
        ss1(is)=dsin(rake1(is)*DEG2RAD)
        cs2(is)=dcos(rake2(is)*DEG2RAD)
        ss2(is)=dsin(rake2(is)*DEG2RAD)
      enddo
c00000000000000000000000000000000000000000000000000000000000000000000000
c     READ IN PARAMETERS FOR OBSERVATION DATA
c     =======================================
c00000000000000000000000000000000000000000000000000000000000000000000000
      call skipdoc(unit)
      read(unit,*)ngd,datunit,nheader
c
      allocate(gddata(ngd),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmgetinp: gddata not allocated!'
      allocate(wfm(ngd),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmgetinp: wfm not allocated!'
      allocate(dinc_const(ngd),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmgetinp: dinc_const not allocated!'
      allocate(dazi_const(ngd),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmgetinp: dazi_const not allocated!'
      allocate(csconst(ngd),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmgetinp: csconst not allocated!'
      allocate(gdout(ngd),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmgetinp: gdout not allocated!'
      allocate(nobsj(ngd),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmgetinp: nobsj not allocated!'
      allocate(nobs1(ngd),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmgetinp: nobs1 not allocated!'
      allocate(nobs2(ngd),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmgetinp: nobs2 not allocated!'
c
      wfsum=0.d0
      do igd=1,ngd
        call skipdoc(unit)
        read(unit,*)gddata(igd)
        call skipdoc(unit)
        read(unit,'(a)')dataline
        read(dataline,*)wfm(igd),j
        if(wfm(igd).le.0.d0)then
          stop ' Error in sdmgetinp: wrong weighting factor!'
        endif
        wfsum=wfsum+wfm(igd)
        if(j.eq.1)then
          csconst(igd)=.true.
          read(dataline,*)wfm(igd),j,
     &                dinc_const(igd),dazi_const(igd)
        else
          csconst(igd)=.false.
        endif
      enddo
      if(wfsum.gt.0.d0)then
        do igd=1,ngd
          wfm(igd)=wfm(igd)/wfsum
        enddo
      endif
c00000000000000000000000000000000000000000000000000000000000000000000000
c     READ IN DATA-CORRECTION PARAMETERS
c     ==================================
c00000000000000000000000000000000000000000000000000000000000000000000000
      call skipdoc(unit)
      read(unit,'(a)')dataline
      read(dataline,*)npar
      if(npar.gt.0)then
        read(dataline,*)npar,corrgrnfile
        allocate(corrpar(npar),stat=ierr)
        if(ierr.ne.0)stop ' Error in sdmgetinp: corrpar not allocated!'
        allocate(parmin(npar),stat=ierr)
        if(ierr.ne.0)stop ' Error in sdmgetinp: parmin not allocated!'
        allocate(parmax(npar),stat=ierr)
        if(ierr.ne.0)stop ' Error in sdmgetinp: parmax not allocated!'
        do i=1,npar
          call skipdoc(unit)
          read(unit,*)parmin(i),parmax(i)
        enddo
      endif
c00000000000000000000000000000000000000000000000000000000000000000000000
c     READ IN PARAMETERS FOR INVERSION REQUIREMENTS
c     =============================================
c00000000000000000000000000000000000000000000000000000000000000000000000
      call skipdoc(unit)
      read(unit,*)niter,zhyrelax
      call skipdoc(unit)
      read(unit,*)ismooth,wei2smo0,izhy
c
c     ismooth = 1: smoothing slip
c               2: smoothing stress drop
c     izhy = 1: using Zhang's approach to increase deep slip resolution
c
      if(ismooth.eq.1)then
        nsmocmp=4
      else if(ismooth.eq.2)then
        nsmocmp=6
      else
        stop ' Error in sdmmain: wrong ´smoothing selection!'
      endif
          
      if(wei2smo0.lt.0.d0)then
        wei2smo0=0.d0
      endif
c00000000000000000000000000000000000000000000000000000000000000000000000
c     READ IN PARAMETERS FOR OUTPUTS
c     ==============================
c00000000000000000000000000000000000000000000000000000000000000000000000
      call skipdoc(unit)
      read(unit,*)slipout
      call skipdoc(unit)
      if(npar.gt.0)then
        read(unit,*)(gdout(i),i=1,ngd),parout
      else
        read(unit,*)(gdout(i),i=1,ngd)
      endif
c
      close(unit)
c
      return
      end
