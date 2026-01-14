      program sdmmain
      use sdmalloc
      implicit none
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     this program synthesizes seismograms due to a number of          c
c     rectanglar rupture planes using the Green's function approach.   c
c     The input data will be read from an input file                   c
c                                                                      c
c     Last modified: Berlin, Jan, 2026, by R. Wang                     c
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
      integer*4 ierr,time,it1,it2
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
      it1=time()
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
      write(*,'(a)')' ... calculate weight of roughness term ...'
      call sdmzhy(ierr)
c
      write(*,'(a)')' ... prepare SDM iteration ...'
      call sdmpreinv(ierr)
      write(*,'(a)')' ... SDM iteration ...'
      call sdmiterate(ierr)
c
      write(*,'(a)')' ... output ...'
      call sdmoutput(ierr)
      it2=time()
c
c00000000000000000000000000000000000000000000000000000000000000000000000
c      END OF STANDARD PROCESSING
c      ==========================
c00000000000000000000000000000000000000000000000000000000000000000000000
      print *,' ####################################################'
      print *,' #                                                  #'
      print *,' #        End of computations with SDM2025          #'
      print *,' #                                                  #'
      print *,' ####################################################'
      print *,' run time: ',it2-it1,'s.'
c
      stop
      end
