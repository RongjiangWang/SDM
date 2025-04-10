      program sdmmain
      implicit none
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     this program synthesizes seismograms due to a number of          c
c     rectanglar rupture planes using the Green's function approach.   c
c     The input data will be read from an input file                   c
c                                                                      c
c     Last modified: Potsdam, Oct, 2008, by R. Wang                    c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     BEGIN DECLARATIONS
c     ==================
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     GLOBAL CONSTANTS
c     ================
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      include 'sdmglob.h'
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     LOCAL WORK SPACES
c     =================
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer i,j,is,iobs,nobs,ns,ngd,nps,ierr
      integer niter,ismooth
      double precision wgrad,swap,tuser
      character*180 dataline
      integer time
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     END DECLARATIONS
c     ================
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c00000000000000000000000000000000000000000000000000000000000000000000000
c     BEGIN READ IN INPUT PARAMETERS
c     ==============================
c00000000000000000000000000000000000000000000000000000000000000000000000
c
      nwarn=0
c
      print *,' #######################################################'
      print *,' #                                                     #'
      print *,' #                   Welcome to                        #'
      print *,' #                                                     #'
      print *,' #                                                     #'
      print *,' #              SSSS     DDDD      M   M               #'
      print *,' #             S         D   D     M M M               #'
      print *,' #              SSS      D   D     M M M               #'
      print *,' #                 S     D   D     M   M               #'
      print *,' #             SSSS      DDDD      M   M               #'
      print *,' #                                                     #'
      print *,' #                   * * * * * * *                     #'
      print *,' # --------------------------------------------------- #'
      print *,' #                 for inverting                       #'
      print *,' #        slip distribution on fault plane             #'
      print *,' #         from surface deformation data               #'
      print *,' #     by the Steepst Descent (or Gradient) Method     #'
      print *,' # --------------------------------------------------- #'
      print *,' #                   * * * * * * *                     #'
      print *,' #                                                     #'
      print *,' #                   Version 2011                      #'
      print *,' #                                                     #'
      print *,' #                       by                            #'
      print *,' #                                                     #'
      print *,' #                  Rongjiang Wang                     #'
      print *,' #               (wang@gfz-potsdam.de)                 #'
      print *,' #                                                     #'
      print *,' #            GeoForschungsZentrum Potsdam             #'
      print *,' #                     June 2011                       #'
      print *,' #######################################################'
      print *,'                                                       '
c      tuser=dble(time()/24/3600)-38.d0*365.25d0
c      if(tuser.gt.180.d0)then
c        pause ' Time for demo applications is over.'
c        stop
c      endif
      write(*,'(a,$)')' Please type the file name of input data: '
      read(*,'(a)')infile
      open(10,file=infile,status='old')
c00000000000000000000000000000000000000000000000000000000000000000000000
c     READ IN PARAMETERS FOR EARTH MODEL CHOICE
c     =========================================
c00000000000000000000000000000000000000000000000000000000000000000000000
      call getdata(10,dataline)
      read(dataline,*)i
      if(i.eq.1)then
        hsmodel=.true.
        call getdata(10,dataline)
        read(dataline,*)poisson
        if(poisson.lt.0.d0.or.poisson.ge.0.5d0)then
          stop ' Wrong Poisson ratio!'
        endif
        zobs=0.d0
      else
        hsmodel=.false.
        call getdata(10,dataline)
        read(dataline,*)grndir,(green(i),i=1,3)
      endif
c00000000000000000000000000000000000000000000000000000000000000000000000
c     READ IN PARAMETERS FOR RECTANGULAR SOURCES
c     ==========================================
c00000000000000000000000000000000000000000000000000000000000000000000000
      call getdata(10,dataline)
      read(dataline,*)ns
      if(ns.gt.NSMAX)then
        print *,'Error in input: too max. number of segments exceeded'
        stop
      endif
      do is=1,ns
        call getdata(10,dataline)
        read(dataline,*)i,topdep(is),width(is),
     &                  mstrike(is),patchsize(is)
        call getdata(10,dataline)
        read(dataline,*)rake1(is),rake2(is),maxslip(is)
        call getdata(10,dataline)
        read(dataline,*)nft(is)
        if(nft(is).lt.2.or.nft(is).gt.NFTMAX)then
          print *,is,'. fault segment:"
     &             //" bad number of fault trace locations!'
          stop
        endif
        do i=1,nft(is)
          call getdata(10,dataline)
          read(dataline,*)latft(i,is),lonft(i,is),
     &                    topdip(i,is),botdip(i,is)
          if(topdip(i,is).le.0.d0.or.topdip(i,is).ge.180.d0.or.
     &       botdip(i,is).le.0.d0.or.botdip(i,is).ge.180.d0)then
            print *,is,'. fault segment: bad dip angle!'
            stop
          endif
        enddo
c
        topdep(is)=KM2M*topdep(is)
        width(is)=KM2M*width(is)
        patchsize(is)=KM2M*patchsize(is)
c
        if(maxslip(is).le.0.d0)then
          stop ' Bad parameter max_slip!'
        endif
c
        if(rake1(is).gt.rake2(is))then
          stop ' Wrong rake range!'
        endif
c
        if(rake1(is).lt.-180.d0)then
          stop ' Wrong rake difinition used!'
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
      call getdata(10,dataline)
      read(dataline,*)ngd,dspunit,nheader
      if(ngd.gt.NGDMAX)then
        stop ' Error: max. no of data sets too small defined!'
      endif
      caloffset=.false.
      do i=1,ngd
        call getdata(10,dataline)
        read(dataline,*)gddata0(i),gddata(i)
        call getdata(10,dataline)
        read(dataline,*)wfm(i),seloffset(i),j
        if(wfm(i).le.0.d0)then
          stop ' Error in sdmmain: wrong weighting factor!'
        endif
        if(j.eq.1)then
          csconst(i)=.true.
          read(dataline,*)wfm(i),seloffset(i),j,
     &                    dinc_const(i),dazi_const(i)
        else
          csconst(i)=.false.
        endif
        if(seloffset(i).eq.1)caloffset=.true.
      enddo
c00000000000000000000000000000000000000000000000000000000000000000000000
c     READ IN PARAMETERS FOR INVERSION REQUIREMENTS
c     =============================================
c00000000000000000000000000000000000000000000000000000000000000000000000
      call getdata(10,dataline)
      read(dataline,*)niter
      call getdata(10,dataline)
      read(dataline,*)ismooth,wgrad
c
c     ismooth = 1: smoothing slip
c               2: smoothing stress drop
c
      if(ismooth.lt.1.or.ismooth.gt.2)then
        stop ' Error in sdmmain: wrong ´smoothing selection!'
      endif
      if(wgrad.lt.0.d0)wgrad=0.d0
      wgrad=wgrad**2
c00000000000000000000000000000000000000000000000000000000000000000000000
c     READ IN PARAMETERS FOR OUTPUTS
c     ==============================
c00000000000000000000000000000000000000000000000000000000000000000000000
      call getdata(10,dataline)
      read(dataline,*)slipout
      read(10,*)(gdout0(i),i=1,ngd)
c
      close(10)
c00000000000000000000000000000000000000000000000000000000000000000000000
c      END READ IN INPUT PARAMETERS
c      ============================
c00000000000000000000000000000000000000000000000000000000000000000000000
      write(*,'(a)')' ... the input file has been read successfully...'
c00000000000000000000000000000000000000000000000000000000000000000000000
c      BEGIN PROCESSING
c      ================
c00000000000000000000000000000000000000000000000000000000000000000000000
      nwarn=0
      write(*,'(a)')' ... read the observation data sets ...'
      call sdmdata(ngd,nobs)
c
      write(*,'(a)')' ... discretise the rupture area ...'
      call sdmdisc(ns,nps)
c
      if(.not.hsmodel)then
        write(*,'(a)')' ... read the differential Green functions ...'
        call sdmdgrn(ierr)
      endif
c
      write(*,'(a)')' ... calculate the complete Green functions ...'
      call sdmgrn(ns,nps,nobs,ismooth)
c
      write(*,'(a)')' ... derive the slip distribution ...'
      call sdminv(ngd,ns,nps,nobs,niter,wgrad)
c
      call sdmout(ns,nps,ngd)
c
c00000000000000000000000000000000000000000000000000000000000000000000000
c      END OF STANDARD PROCESSING
c      ==========================
c00000000000000000000000000000000000000000000000000000000000000000000000
      if(nwarn.eq.0)then
        print *,' ####################################################'
        print *,' #                                                  #'
        print *,' #        End of computations with SDM2011          #'
        print *,' #                                                  #'
        print *,' ####################################################'
      else
        print *,' ####################################################'
        print *,'      Sorry, there have been',nwarn,' warnings.      '
        print *,'              Results may be inaccurate!             '
        print *,' ####################################################'
      endif
c
      stop
      end
