      program sdmmain
      use sdmalloc
      implicit none
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     this program synthesizes seismograms due to a number of          c
c     rectanglar rupture planes using the Green's function approach.   c
c     The input data will be read from an input file                   c
c                                                                      c
c     Last modified: Zhuhai, Nov, 2025, by R. Wang                    c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     BEGIN DECLARATIONS
c     ==================
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     GLOBAL CONSTANTS
c     ================
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     LOCAL WORK SPACES
c     =================
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer*4 ierr
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
      print *,' #                  to invert                          #'
      print *,' #             surface deformation data                #'
      print *,' #         for slip distribution on fault plane        #'
      print *,' #     by the Steepst Descent (or Gradient) Method     #'
      print *,' # --------------------------------------------------- #'
      print *,' #                   * * * * * * *                     #'
      print *,' #                                                     #'
      print *,' #                   Version 2025                      #'
      print *,' #                                                     #'
      print *,' #                       by                            #'
      print *,' #                                                     #'
      print *,' #                  Rongjiang Wang                     #'
      print *,' #               (wang@gfz-potsdam.de)                 #'
      print *,' #                                                     #'
      print *,' #        GFZ Helmholtz Centre for Geosciences         #'
      print *,' #             Last modified Nov 28, 2025              #'
      print *,' #######################################################'
      print *,'                                                       '
c
      write(*,'(a,$)')' Please type the file name of input data: '
      read(*,'(a)')infile
      call sdmgetinp(ierr)
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
      call sdmdata(ierr)
c
      write(*,'(a)')' ... discretise the rupture area ...'
      call sdmdisc(ierr)
c
      if(.not.hsmodel)then
        write(*,'(a)')' ... read the differential Green functions ...'
        call sdmdgrn(ierr)
      endif
c
      write(*,'(a)')' ... calculate the complete Green functions ...'
      call sdmgrn(ierr)
c
      if(niter.gt.0)then
        write(*,'(a)')' ... calculate weight of roughness term ...'
        call sdmzhy(ierr)
      else
        write(*,'(a)')' ... forward modeling ...'
      endif
c
      write(*,'(a)')' ... prepare SDM iteration ...'
      call sdmpreinv(ierr)
      write(*,'(a)')' ... SDM iteration ...'
      call sdmiterate(ierr)
c
      write(*,'(a)')' ... output ...'
      call sdmoutput(ierr)
c
c00000000000000000000000000000000000000000000000000000000000000000000000
c      END OF STANDARD PROCESSING
c      ==========================
c00000000000000000000000000000000000000000000000000000000000000000000000
      if(nwarn.eq.0)then
        print *,' ####################################################'
        print *,' #                                                  #'
        print *,' #        End of computations with SDM2025          #'
        print *,' #                                                  #'
        print *,' ####################################################'
      else
        print *,' ####################################################'
        print *,'      There have been',nwarn,' warnings.      '
        print *,'              Results may be inaccurate!             '
        print *,' ####################################################'
      endif
c
      stop
      end
