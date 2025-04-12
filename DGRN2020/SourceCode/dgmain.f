      program dgmain
      use dgalloc
      implicit none
c
c     work space
c
      integer*4 i,j,l,ir,nr,izs,ierr
      integer*4 istp,isp,nr1,nr2,nzs,nls,nlr
      integer*4 lend,lenf,leninp,iunit
      integer*4 unit(3,4)
      real*8 r1,r2,dr,rratio,dract
      real*8 pi,accuracy,zdis
      real*8 zs1,zs2,dzs,zratio,dzact,vp,vs
      real*8 swap
      logical homog,select(3,4)
      character*35 stype(4)
      character*35 comptxt(3)
      character*80 inputfile,fname(3),outdir
      character*163 green(3,4)
c
      pi=4.d0*datan(1.d0)
c
c     read input file file
c
      print *,' ######################################################'
      print *,' #                                                    #'
      print *,' #                  Welcome to                        #'
      print *,' #                                                    #'
      print *,' #                                                    #'
      print *,' #        DDDD       GGGG     RRRR      N   N         #'
      print *,' #        D   D     G         R   R     NN  N         #'
      print *,' #        D   D     G GGG     RRRR      N N N         #'
      print *,' #        D   D     G   G     R R       N  NN         #'
      print *,' #        DDDD       GGGG     R  R      N   N         #'
      print *,' #                                                    #'
      print *,' #                                                    #'
      print *,' #                   * * * * * * *                    #'
      print *,' # -------------------------------------------------- #'
      print *,' #                 for calculatng                     #'
      print *,' #        the differential Green functions            #'
      print *,' #      of a multi-layered elastic halfspace          #'
      print *,' # (in reference to a homogeneous elastic halfspace)  #'
      print *,' # -------------------------------------------------- #'
      print *,' #                   * * * * * * *                    #'
      print *,' #                                                    #'
      print *,' #                  Version 2020                      #'
      print *,' #                                                    #'
      print *,' #                      by                            #'
      print *,' #                 Rongjiang Wang                     #'
      print *,' #              (wang@gfz-potsdam.de)                 #'
      print *,' #                                                    #'
      print *,' #     GFZ German Research Centre for Geosciences     #'
      print *,' #                   April 2020                       #'
      print *,' ######################################################'
      print *,'                                                      '
c
      nwarn=0
      write(*,'(a,$)')' Please type the file name of input data: '
      read(*,'(a)')inputfile
      open(10,file=inputfile,status='old')
c
c     parameters for source-observation array
c     =======================================
c
      call skipdoc(10)
      read(10,*)zrec
      zrec=zrec*km2m
      call skipdoc(10)
      read(10,*)nr,r1,r2,rratio
      r1=r1*km2m
      r2=r2*km2m
      if(r1.lt.0.d0.or.r2.lt.r1.or.nr.le.1)then
        stop 'Error: wrong no of distance samples!'
      else if(nr.eq.1.or.r1.eq.r2)then
        stop 'Error: bad parameters for distance samples!'
      endif
c
      allocate(r(nr),stat=ierr)
      if(ierr.ne.0)stop 'Error in dgmain: r not allocated!'
      allocate(u(nr,3,4),stat=ierr)
      if(ierr.ne.0)stop 'Error in dgmain: u not allocated!'
      allocate(r0(nr),stat=ierr)
      if(ierr.ne.0)stop 'Error in dgmain: r0 not allocated!'
      allocate(u0(nr,3,4),stat=ierr)
      if(ierr.ne.0)stop 'Error in dgmain: u0 not allocated!'
c
      if(nr.eq.2)then
        dr=r2-r1
        r(1)=r1
        r(2)=r2
      else
        dr=2.d0*(r2-r1)/dble(nr-1)/(1.d0+rratio)
        r(1)=r1
        do ir=2,nr
          dract=dr*(1.d0+(rratio-1.d0)*dble(ir-2)/dble(nr-2))
          r(ir)=r(ir-1)+dract
        enddo
      endif
c
      call skipdoc(10)
      read(10,*)nzs,zs1,zs2,zratio
      zs1=zs1*km2m
      zs2=zs2*km2m
      if(zs1.le.0.d0.or.zs2.lt.zs1.or.nzs.le.1)then
        stop 'Error: wrong no of depth samples!'
      else if(nzs.eq.1.or.zs1.eq.zs2)then
        stop 'Error: bad parameters for depth samples!'
      endif
c
      allocate(zs0(nzs),stat=ierr)
      if(ierr.ne.0)stop 'Error in dgmain: zs0 not allocated!'
c
      if(nzs.eq.2)then
        dzs=zs2-zs1
        zs0(1)=zs1
        zs0(2)=zs2
      else
        dzs=2.d0*(zs2-zs1)/dble(nzs-1)/(1.d0+zratio)
        zs0(1)=zs1
        do izs=2,nzs
          dzact=dzs*(1.d0+(zratio-1.d0)*dble(izs-2)/dble(nzs-2))
          zs0(izs)=zs0(izs-1)+dzact
        enddo
      endif
c
c     wavenumber integration parameters
c     =================================
c
      call skipdoc(10)
      read(10,*)accuracy
      if(accuracy.le.0.d0.or.accuracy.ge.1.d0)accuracy=0.1d0
c
c     parameters for output files
c     ===========================
c
      call skipdoc(10)
      read(10,*)outdir
c
      do lend=80,1,-1
        if(outdir(lend:lend).ne.' ')goto 100
      enddo
100   continue
c
      if(lend.lt.1)then
        stop 'Error: wrong format for output directory!'
      endif
c
      call skipdoc(10)
      read(10,*)(fname(i),i=1,3)
      do i=1,3
        do lenf=80,1,-1
          if(fname(i)(lenf:lenf).ne.' ')goto 110
        enddo
110     continue
        green(i,1)=outdir(1:lend)//fname(i)(1:lenf)//'.ep'
        green(i,2)=outdir(1:lend)//fname(i)(1:lenf)//'.ss'
        green(i,3)=outdir(1:lend)//fname(i)(1:lenf)//'.ds'
        green(i,4)=outdir(1:lend)//fname(i)(1:lenf)//'.cl'
        do istp=1,4
          select(i,istp)=.true.
        enddo
      enddo
c
c     no tangential components for clvd sources
c
      select(3,1)=.false.
      select(3,4)=.false.
c
c     global model parameters
c     =======================
c
      call skipdoc(10)
      read(10,*)l
c
      allocate(h0(l),stat=ierr)
      if(ierr.ne.0)stop 'Error in dgmain: h0 not allocated!'
      allocate(la0(l),stat=ierr)
      if(ierr.ne.0)stop 'Error in dgmain: la0 not allocated!'
      allocate(mu0(l),stat=ierr)
      if(ierr.ne.0)stop 'Error in dgmain: mu0 not allocated!'
      allocate(rho0(l),stat=ierr)
      if(ierr.ne.0)stop 'Error in dgmain: rho0 not allocated!'
c
      allocate(z1(l),stat=ierr)
      if(ierr.ne.0)stop 'Error in dgmain: z1 not allocated!'
      allocate(la1(l),stat=ierr)
      if(ierr.ne.0)stop 'Error in dgmain: la1 not allocated!'
      allocate(mu1(l),stat=ierr)
      if(ierr.ne.0)stop 'Error in dgmain: mu1 not allocated!'
      allocate(rho1(l),stat=ierr)
      if(ierr.ne.0)stop 'Error in dgmain: rho1 not allocated!'
c
      allocate(z2(l),stat=ierr)
      if(ierr.ne.0)stop 'Error in dgmain: z2 not allocated!'
      allocate(la2(l),stat=ierr)
      if(ierr.ne.0)stop 'Error in dgmain: la2 not allocated!'
      allocate(mu2(l),stat=ierr)
      if(ierr.ne.0)stop 'Error in dgmain: mu2 not allocated!'
      allocate(rho2(l),stat=ierr)
      if(ierr.ne.0)stop 'Error in dgmain: rho2 not allocated!'
c
c      multilayered model parameters
c      =============================
c
      do i=1,l
        call skipdoc(10)
        read(10,*)j,h0(i),vp,vs,rho0(i)
c
c       convert to SI units
c
        h0(i)=h0(i)*km2m
        vp=vp*km2m
        vs=vs*km2m
        rho0(i)=rho0(i)*km2m
        mu0(i)=rho0(i)*vs*vs
        la0(i)=rho0(i)*vp*vp-2.d0*mu0(i)
        if(la0(i).le.0.d0)then
          stop 'inconsistent Vp/Vs ratio!'
        endif
      enddo
      if(l.eq.1)h0(l)=0.d0
      close(10)
c
c     end of inputs
c     =============
c
      homog=.true.
      do i=2,l
        if(la0(i).ne.la0(i-1).or.
     &     mu0(i).ne.mu0(i-1).or.
     &     rho0(i).ne.rho0(i-1))homog=.false.
      enddo
      if(homog)then
        stop ' Error: no layered earth model is given!'
      endif
c
      comptxt(1)='Uz (vertical displacement)'
      comptxt(2)='Ur (radial displacement)'
      comptxt(3)='Ut (tangential displacement)'
c
c     determine upper und lower parameter values of each layer
c
      l0=1
      z1(l0)=0.d0
      do i=2,l
        if(h0(i).gt.h0(i-1))then
          z1(l0)=h0(i-1)
          la1(l0)=la0(i-1)
          mu1(l0)=mu0(i-1)
          rho1(l0)=rho0(i-1)
c
          z2(l0)=h0(i)
          la2(l0)=la0(i)
          mu2(l0)=mu0(i)
          rho2(l0)=rho0(i)
          l0=l0+1
        else
          z1(l0)=h0(i)
          la1(l0)=la0(i)
          mu1(l0)=mu0(i)
          rho1(l0)=rho0(i)
        endif
      enddo
      z1(l0)=h0(l)
      la1(l0)=la0(l)
      mu1(l0)=mu0(l)
      rho1(l0)=rho0(l)
c
c     construction of sublayers
c
      write(*,*)'the multi-layered poroelastic model:'
c
      call dgsublay(ierr)
      if(ierr.eq.1)then
        stop 'the max. no of layers (lmax) too small defined!'
      endif
c
      zs=0.d0
      call dglayer(ierr)
      nlr=nno(lzrec)
c
      leninp=index(inputfile,' ')-1
c
      stype(1)='explosion (M11=M22=M33=1*kappa)'
      stype(2)='strike-slip (M12=M21=1*mue)'
      stype(3)='dip-slip (M13=M31=1*mue)'
      stype(4)='clvd (M33=1*mue, M11=M22=-M33/2)'
c
      call dgbsjtab(ierr)
c
      iunit=10
      do istp=1,4
        do i=1,3
          if(select(i,istp))then
            iunit=iunit+1
            unit(i,istp)=iunit
            open(unit(i,istp),file=green(i,istp),status='unknown')
	    write(unit(i,istp),'(a)')'################################'
	    write(unit(i,istp),'(a)')'# The input file used: '
     &                        //inputfile(1:leninp)
	    write(unit(i,istp),'(a)')'################################3'
	    write(unit(i,istp),'(a)')'# Greens function component: '
     &                        //comptxt(i)
	    write(unit(i,istp),'(a)')'# Source type: '//stype(istp)
	    write(unit(i,istp),'(a)')'# Observation distance sampling:'
	    write(unit(i,istp),'(a)')'#    nr        r1[m]        r2[m]'
     &                             //'       rratio'
	    write(unit(i,istp),'(i7,3E14.6)')nr,r1,r2,rratio
	    write(unit(i,istp),'(a)')'# Uniform obs. site parameters:'
	    write(unit(i,istp),'(a)')'#    depth[m]       la[Pa]       '
     &                             //'mu[Pa]  rho[kg/m^3]'
	    write(unit(i,istp),'(6d13.6)')zrec,la(nlr),mu(nlr),rho(nlr)
	    write(unit(i,istp),'(a)')'# Source depth sampling:'
	    write(unit(i,istp),'(a)')'#   nzs       zs1[m]       zs2[m]'
     &                             //'       zratio'
	    write(unit(i,istp),'(i7,3E14.6)')nzs,zs1,zs2,zratio
	    write(unit(i,istp),'(a)')'# Data in each source depth block'
	    write(unit(i,istp),'(a)')'# ==============================='
	    write(unit(i,istp),'(a)')'# 1. line: source layer parameters'
	    write(unit(i,istp),'(a)')'#  s_depth, la, mu, rho'
	    write(unit(i,istp),'(a)')'# 2. line: (f(ir),ir=1,nr)'
	  endif
        enddo
      enddo
      do izs=1,nzs
        zs=zs0(izs)
        write(*,'(a,E13.4,a)')' Processing for the '
     &                 //'source at depth:',zs,' m.'
c
        call dglayer(ierr)
        nls=nno(ls)
c
c       parameters for hankel integrations
c
        zdis=dabs(zrec-zs)
        if(zdis.eq.0.d0.and.r1.le.0.d0)then
          stop ' Error in dgmain: observation at source point!'
        endif
c
        isp=0
        nr1=1
        nr2=0
        do ir=nr1,nr
          if(r(ir).le.2.5d0*dsqrt(zdis**2+r(1)**2))nr2=ir
        enddo
        if(nr2.ge.nr1)then
          isp=1
          write(*,'(i3,a,F12.4,a,f12.4,a)')isp,'. sub-profile: ',
     &             r(nr1),' -> ',r(nr2),' m.'
c
          call dgwvint(nr1,nr2,accuracy)
        else
          isp=0
          nr2=0
        endif
200     isp=isp+1
        nr1=nr2+1
        nr2=nr1
        do ir=nr1+1,nr
          if(r(ir).le.2.5d0*r(nr1))nr2=ir
        enddo
        if(nr2.le.nr)then
          write(*,'(i3,a,F12.4,a,f12.4,a)')isp,'. sub-profile: ',
     &           r(nr1),' -> ',r(nr2),' m'
          call dgwvint(nr1,nr2,accuracy)
        endif
        if(nr2.lt.nr)goto 200
c
        do istp=1,4
          do i=1,3
            if(.not.select(i,istp))goto 400
	    write(unit(i,istp),'(a)')'##################################'
	    write(unit(i,istp),'(a,i2,a,E16.6,a)')
     &        '# the ',izs,'. source depth: ',zs,' m'
	    write(unit(i,istp),'(a)')'##################################'
	    write(unit(i,istp),'(6E13.6)')zs,la(nls),mu(nls),
     &                             rho(nls)
	    do ir=1,nr-1
	      write(unit(i,istp),'(E13.5,$)')u(ir,i,istp)
            enddo
	      write(unit(i,istp),'(E13.5)')u(nr,i,istp)
400         continue
          enddo
        enddo
      enddo
c
c     end of izs loop
c
800   continue
      if(nwarn.eq.0)then
        print *,'####################################################'
        print *,'#                                                  #'
        print *,'#          End of computations with DGRN2020       #'
        print *,'#                                                  #'
        print *,'####################################################'
      else
        print *,'####################################################'
        print *,'     Sorry, there have been',nwarn,' warnings.      '
        print *,'             Results may be inaccurate!             '
        print *,'####################################################'
      endif
      stop
      end
