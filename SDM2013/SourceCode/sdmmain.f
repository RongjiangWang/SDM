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
      integer*4 i,j,is,iobs,nobs,ns,ngd,nps,ierr
      integer*4 niter,ismooth,idip
      real*8 wgrad,swap,tuser,x,y,z,btmdip
      character*180 dataline
      integer*4 time
c
      real*8 fdip
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
      print *,' #                   Version 2013a                     #'
      print *,' #                                                     #'
      print *,' #                       by                            #'
      print *,' #                                                     #'
      print *,' #                  Rongjiang Wang                     #'
      print *,' #               (wang@gfz-potsdam.de)                 #'
      print *,' #                                                     #'
      print *,' #            GeoForschungsZentrum Potsdam             #'
      print *,' #             Last update Nov 20, 2025                #'
      print *,' #######################################################'
      print *,'                                                       '
c
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
      read(dataline,*)idisc
      if(idisc.lt.0.or.idisc.gt.1)then
        print *,'Error in input: bad switch for subfaults input form'
        stop
      endif
      if(idisc.eq.0)then
        call getdata(10,dataline)
        read(dataline,*)ns
        do is=1,ns
          call getdata(10,dataline)
          read(dataline,*)i,rake1(is),rake2(is),
     &                    maxslip(is),inpatches(is),iref(is)
        enddo
      else
        call getdata(10,dataline)
        read(dataline,*)ns
        do is=1,ns
          call getdata(10,dataline)
          read(dataline,*)i,topdep(is),btmdep(is),
     &         patchsize(is)
          if(topdep(is).ge.btmdep(is))then
            print *,is,'. subfault: bad top/bottom depth!'
            stop
          endif
          topdep(is)=KM2M*topdep(is)
          btmdep(is)=KM2M*btmdep(is)
          patchsize(is)=KM2M*patchsize(is)
c
          call getdata(10,dataline)
          read(dataline,*)rake1(is),rake2(is),maxslip(is)
          call getdata(10,dataline)
          read(dataline,*)nwft(is),idip,ndgrid(is)
          if(nwft(is).lt.2.or.nwft(is).gt.NWFTMAX)then
            print *,is,'. subfault: bad number of fault wefts!'
            stop
          endif
          do i=1,nwft(is)
            call getdata(10,dataline)
            if(idip.eq.1)then
              read(dataline,*)toplat(i,is),toplon(i,is),
     &                        btmlat(i,is),btmlon(i,is),topdip(i,is)
            else if(idip.eq.2)then
              read(dataline,*)toplat(i,is),toplon(i,is),
     &                        btmlat(i,is),btmlon(i,is),btmdip
              if(toplat(i,is).eq.btmlat(i,is).and.
     &           toplon(i,is).eq.btmlon(i,is))then
                topdip(i,is)=90.d0
              else if(btmdip.le.0.d0.or.btmdip.ge.180.d0.or.
     &                btmdip.eq.90.d0)then
                print *,is,'. subfault: bad dip angle!'
                stop
              else
                call disazi(REARTH,toplat(i,is),toplon(i,is),
     &                             btmlat(i,is),btmlon(i,is),x,y)
                x=dsqrt(x**2+y**2)
                z=btmdep(is)-topdep(is)
                if(btmdip.le.90.d0)then
                  topdip(i,is)=fdip(btmdip,x,z)
                else
                  topdip(i,is)=180.d0-fdip(180.d0-btmdip,x,z)
                endif
              endif
            else
              read(dataline,*)toplat(i,is),toplon(i,is),
     &                        btmlat(i,is),btmlon(i,is)
              if(toplat(i,is).eq.btmlat(i,is).and.
     &           toplon(i,is).eq.btmlon(i,is))then
                topdip(i,is)=90.d0
              else
                call disazi(REARTH,toplat(i,is),toplon(i,is),
     &                             btmlat(i,is),btmlon(i,is),x,y)
                x=dsqrt(x**2+y**2)
                z=btmdep(is)-topdep(is)
                topdip(i,is)=datan(z/x)/DEG2RAD
              endif
            endif
            if(i.gt.1)then
              if(topdip(i,is).gt.90.d0.and.topdip(i-1,is).lt.90.d0.or.
     &           topdip(i,is).lt.90.d0.and.topdip(i-1,is).gt.90.d0)then
                print *, 'Error: dip changes from < 90° to > 90°'
     &                 //'within one of weft intervals.'
                stop
              endif
            endif
          enddo
          call reweft(NWFTMAX,nwft(is),ndgrid(is),
     &                toplat(1,is),toplon(1,is),
     &                btmlat(1,is),btmlon(1,is),
     &                topdip(1,is),patchsize(is),btmdep(is)-topdep(is))
        enddo
      endif
      do is=1,ns
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
c               3: smoothing slip by Zhangs approach
c
      if(ismooth.lt.1.or.ismooth.gt.3)then
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
      write(*,'(a)')' ... calculate smoothing weights...'
      call sdmzhy(nps,nobs,ismooth)
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
        print *,' #        End of computations with SDM2013          #'
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
