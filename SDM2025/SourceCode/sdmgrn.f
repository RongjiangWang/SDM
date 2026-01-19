      subroutine sdmgrn(ierr)
      use sdmalloc
      implicit none
      integer*4 ierr
c
c     Last modified: Jan 2026 in Berlin by R. Wang
c
      integer*4 is,ips,jps,iobs,nlarge,ira,il,nl,lend,lenf
      integer*4 i,j,k,l,ir,ir1,ir2,izs,izs1,izs2,idiv
      integer*4 np10,npercent
      real*8 dr,ddl,st,di,dux,duy,duz,dur,dut
      real*8 xobs,yobs,xobs0,yobs0,xps,yps,dobs,dal,daw,wfsum
      real*8 eii,exx,eyy,ezz,exy,eyz,ezx
      real*8 si,si0,bsi,co,co0,bco,si2,co2,dis,dis0,azi,azi0,bazi
      real*8 w1,w2,psss,psds,pscl,shss,shds,strst,strdi,strnn
      real*8 dz1,dz2,dpz,dwei,pz1,pz2,xp,yp,px,py
      real*8 sig(3,3),dsp(3),stress(3),dux0(2),duy0(2),duz0(2)
      character*10 header
c
c     for calling Okada's subroutine DC3D0
c
      integer*4 IRET
      REAL*4 ALPHA,X,Y,Z,DEPTH,DIPS,DX,DY,DZ,
     &       UX,UY,UZ,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ
c
c     more details see Okada's subroutine DC3D
c
      REAL*4 AL1,AL2,AW1,AW2,DISL1,DISL2,DISL3
      REAL*4 SS(2),DS(2)
c
      real*8, allocatable:: csst(:),ssst(:)
      real*8, allocatable:: cs2st(:),ss2st(:)
      real*8, allocatable:: csdi(:),ssdi(:)
      real*8, allocatable:: rnn(:,:),rst(:,:),rdi(:,:)
c
      nlarge=0
c
      allocate(csst(nps),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmgrn: csst not allocated!'
      allocate(ssst(nps),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmgrn: ssst not allocated!'
      allocate(cs2st(nps),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmgrn: cs2st not allocated!'
      allocate(ss2st(nps),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmgrn: ss2st not allocated!'
      allocate(csdi(nps),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmgrn: csdi not allocated!'
      allocate(ssdi(nps),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmgrn: ssdi not allocated!'
      allocate(rnn(3,nps),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmgrn: rnn not allocated!'
      allocate(rst(3,nps),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmgrn: rst not allocated!'
      allocate(rdi(3,nps),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmgrn: rdi not allocated!'
c
      allocate(datgrn(2,nps,nobs),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmgrn: datgrn not allocated!'
      allocate(strgrn(2,nps,3,nps),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmgrn: strgrn not allocated!'
c
      allocate(corrmdl(nobs),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmgrn: corrmdl not allocated!'
      if(nusrp.gt.0)then
        allocate(corrgrn(nusrp,nobs),stat=ierr)
        if(ierr.ne.0)stop ' Error in sdmgrn: corrgrn not allocated!'
        allocate(usrpunit(nusrp),stat=ierr)
        if(ierr.ne.0)stop ' Error in sdmgrn: usrpunit not allocated!'
        open(20,file=corrgrnfile,status='old')
        read(20,'(a)')header
        do j=1,nusrp
          usrpunit(j)=0.d0
        enddo
        wfsum=0.d0
        do iobs=1,nobs
          read(20,*)(corrgrn(j,iobs),j=1,nusrp)
          do j=1,nusrp
            usrpunit(j)=usrpunit(j)+(wf(iobs)*corrgrn(j,iobs))**2
          enddo
          wfsum=wfsum+wf(iobs)**2
        enddo
        do j=1,nusrp
          usrpunit(j)=dsqrt(usrpunit(j)/wfsum)
        enddo
        close(20)
      endif
c
c     SUPERPOSITION OF ALL DISCRETE POINT SOURCES
c     ===========================================
c
      do iobs=1,nobs
        do ips=1,nps
          do ira=1,2
            datgrn(ira,ips,iobs)=0.d0
          enddo
        enddo
      enddo
      if(iearth.eq.2)then
        do lend=len(grndir),1,-1
          if(grndir(lend:lend).ne.' ')goto 100
        enddo
100     continue
        if(grndir(lend:lend).ne.'\'.and.grndir(lend:lend).ne.'/')then
          grndir=grndir(1:lend)//'/'
          lend=lend+1
        endif
c
c       read user's Green's functions
c
        do ira=1,2
          do lenf=len(usr3dgrn(ira)),1,-1
            if(usr3dgrn(ira)(lenf:lenf).ne.' ')goto 110
          enddo
110       continue
          usr3dgrn(ira)=grndir(1:lend)//usr3dgrn(ira)(1:lenf)
          open(20,file=usr3dgrn(ira),status='old')
          read(20,'(a)')header
          do iobs=1,nobs
            read(20,*,end=120)(datgrn(ira,ips,iobs),ips=1,nps)
          enddo
          iobs=iobs-1
120       close(20)
          if(iobs.lt.nobs)then
            print *,' Error in sdmgrn: a problem arises when reading '
     &            //'data from '//usr3dgrn(ira)
            stop
          endif
        enddo
      endif
c
      do ips=1,nps
        st=strike(ips)*DEG2RAD
	  di=dip(ips)*DEG2RAD
c
        csst(ips)=dcos(st)
        ssst(ips)=dsin(st)
        cs2st(ips)=dcos(2.d0*st)
        ss2st(ips)=dsin(2.d0*st)
c
        csdi(ips)=dcos(di)
        ssdi(ips)=dsin(di)
c
c       cosine vectors of normal of hanging plate
c
c        rnn(1,ips)= ssdi(ips)*ssst(ips)
c        rnn(2,ips)=-ssdi(ips)*csst(ips)
c        rnn(3,ips)= csdi(ips)
c
c       modified at Jan 19, 2026 by R. Wang:
c
c       cosine vectors of normal of foot plate
c
        rnn(1,ips)=-ssdi(ips)*ssst(ips)
        rnn(2,ips)= ssdi(ips)*csst(ips)
        rnn(3,ips)=-csdi(ips)
c
c       cosine vectors of strike direction
c
        rst(1,ips)=csst(ips)
        rst(2,ips)=ssst(ips)
        rst(3,ips)=0.d0
c
c       cosine vectors of up-dip (thrust) direction
c
        rdi(1,ips)= csdi(ips)*ssst(ips)
        rdi(2,ips)=-csdi(ips)*csst(ips)
        rdi(3,ips)=-ssdi(ips)
      enddo
c
      if(hsmodel)then
        ALPHA=sngl(0.5d0/(1.d0-poisson))
      else
        ALPHA=sngl((laobs+muobs)/(laobs+2.d0*muobs))
      endif
      SS(1)=1.0
      DS(1)=0.0
      SS(2)=0.0
      DS(2)=1.0
      DISL3=0.0
c
      np10=max0(1,nps/10)
c
      write(*,'(a)')' analytical Okada solutions - please wait: '
      do ips=1,nps
c
        DEPTH=sngl(pz(ips))
        DIPS=sngl(dip(ips))
c
        AL1=sngl(-0.5d0*dlen(ips))
        AL2=sngl(+0.5d0*dlen(ips))
        AW1=sngl(-0.4999999d0*dwid(ips))
        AW2=sngl(+0.4999999d0*dwid(ips))
c
c       displacement Green functions
c
        if(iearth.ne.2)then
          do iobs=1,nobs
c
c           transform from Aki's to Okada's system
c
            call disazi(REARTH,plat(ips),plon(ips),
     &                  latobs(iobs),lonobs(iobs),xobs,yobs)
c
            dis=dsqrt(xobs**2+yobs**2)
            if(dis.gt.dsqrt(dlen(ips)*dwid(ips)))then
              azi=datan2(yobs,xobs)
              call disazi(REARTH,latobs(iobs),lonobs(iobs),
     &                    plat(ips),plon(ips),xps,yps)
              bazi=datan2(yps,xps)-PI
            else
              azi=0.d0
              bazi=0.d0
            endif
            co=dcos(azi)
            si=dsin(azi)
            bco=dcos(bazi)
            bsi=dsin(bazi)
c
            X=sngl(xobs*csst(ips)+yobs*ssst(ips))
            Y=sngl(xobs*ssst(ips)-yobs*csst(ips))
            Z=-sngl(zobs)
c
            do ira=1,2
              DISL1=SS(ira)
              DISL2=DS(ira)
              IRET=1
              call DC3D(ALPHA,X,Y,Z,DEPTH,DIPS,AL1,AL2,AW1,AW2,
     &              DISL1,DISL2,DISL3,UX,UY,UZ,
     &              UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ,IRET)
c
c             transform from Okada's to Aki's system
c
              dux=dble(UX)*csst(ips)+dble(UY)*ssst(ips)
              duy=dble(UX)*ssst(ips)-dble(UY)*csst(ips)
              dur=dux*co+duy*si
              dut=duy*co-dux*si
              dux=dur*bco-dut*bsi
              duy=dur*bsi+dut*bco
              duz=-dble(UZ)
              datgrn(ira,ips,iobs)=datgrn(ira,ips,iobs)
     &              +dux*xcs(iobs)+duy*ycs(iobs)+duz*zcs(iobs)
            enddo
          enddo
        endif
c
c       stress drop Green's functions
c       strgrn(i,ips,j,jps) = contribution of i-th slip component
c       at patch ips to j-th stress-drop component at patch jps
c
        do jps=1,nps
c
c         transform from Aki's to Okada's system
c
          call disazi(REARTH,plat(ips),plon(ips),
     &                       plat(jps),plon(jps),xobs,yobs)
          X=sngl(xobs*csst(ips)+yobs*ssst(ips))
          Y=sngl(xobs*ssst(ips)-yobs*csst(ips))
          Z=-sngl(pz(jps))
c
          do ira=1,2
            DISL1=SS(ira)
            DISL2=DS(ira)
            IRET=1
            call DC3D(ALPHA,X,Y,Z,DEPTH,DIPS,AL1,AL2,AW1,AW2,
     &            DISL1,DISL2,DISL3,UX,UY,UZ,
     &            UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ,IRET)
c
c           transform from Okada's to Aki's system
c
            exx=dble(UXX)*csst(ips)*csst(ips)
     &         +dble(UYY)*ssst(ips)*ssst(ips)
     &         +0.5d0*dble(UXY+UYX)*ss2st(ips)
            eyy=dble(UXX)*ssst(ips)*ssst(ips)
     &         +dble(UYY)*csst(ips)*csst(ips)
     &         -0.5d0*dble(UXY+UYX)*ss2st(ips)
            ezz=dble(UZZ)
            exy=0.5d0*dble(UXX-UYY)*ss2st(ips)
     &         -0.5d0*dble(UXY+UYX)*cs2st(ips)
            eyz=-0.5d0*dble(UZX+UXZ)*ssst(ips)
     &          +0.5d0*dble(UYZ+UZY)*csst(ips)
            ezx=-0.5d0*dble(UZX+UXZ)*csst(ips)
     &          -0.5d0*dble(UYZ+UZY)*ssst(ips)
            eii=exx+eyy+ezz
c
            sig(1,1)=lamhs*eii+2.d0*muehs*exx
            sig(2,2)=lamhs*eii+2.d0*muehs*eyy
            sig(3,3)=lamhs*eii+2.d0*muehs*ezz
            sig(1,2)=2.d0*muehs*exy
            sig(2,3)=2.d0*muehs*eyz
            sig(3,1)=2.d0*muehs*ezx
            sig(2,1)=sig(1,2)
            sig(3,2)=sig(2,3)
            sig(1,3)=sig(3,1)
c
c           stress drop at jps-th patch
c
            do i=1,3
              stress(i)=0.d0
              do j=1,3
                stress(i)=stress(i)+sig(i,j)*rnn(j,jps)
              enddo
            enddo
c
            strst=0.d0
            strdi=0.d0
            strnn=0.d0
            do i=1,3
              strst=strst+stress(i)*rst(i,jps)
              strdi=strdi+stress(i)*rdi(i,jps)
              strnn=strnn+stress(i)*rnn(i,jps)
            enddo
            strgrn(ira,ips,1,jps)=strst
            strgrn(ira,ips,2,jps)=strdi
            strgrn(ira,ips,3,jps)=strnn
          enddo
        enddo
c
        if(mod(ips,np10).eq.0)then
          npercent=10*(ips/np10)
          write(*,'(a,i3,a)')'    - - - ',npercent,
     &                       '% calculated ...'
        endif
      enddo
c
      if(hsmodel)goto 500
c
      dr=r(2)-r(1)
      write(*,'(a)')' numerical layering effect - please wait: '
c
      allocate(iz(nps),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmgrn: iz not allocated!'
c
      np10=max0(1,nobs/10)
c
      do iobs=1,nobs
        do ips=1,nps
          do ira=1,2
            dux0(ira)=0.d0
            duy0(ira)=0.d0
            duz0(ira)=0.d0
          enddo
          call disazi(REARTH,plat(ips),plon(ips),
     &                latobs(iobs),lonobs(iobs),xobs0,yobs0)
          dis0=dsqrt(xobs0**2+yobs0**2)
          if(dis0.gt.dsqrt(dlen(ips)*dwid(ips)))then
            azi0=datan2(yobs0,xobs0)
            call disazi(REARTH,latobs(iobs),lonobs(iobs),
     &                  plat(ips),plon(ips),xps,yps)
            bazi=datan2(yps,xps)-PI
          else
            azi0=0.d0
            bazi=0.d0
          endif
          co0=dcos(azi0)
          si0=dsin(azi0)
          bco=dcos(bazi)
          bsi=dsin(bazi)
c
          pz1=pz(ips)-0.5d0*dwid(ips)*ssdi(ips)
          pz2=pz(ips)+0.5d0*dwid(ips)*ssdi(ips)
c
          izs1=1
          do izs=1,nzs
            if(zs(izs).le.pz1)then
              izs1=izs
            endif
          enddo
c
          izs2=nzs
          do izs=nzs,1,-1
            if(zs(izs).ge.pz2)then
              izs2=izs
            endif
          enddo
c
          iz(ips)=(izs1+izs2)/2
c
          nl=1+idint(dlen(ips)/dr)
          ddl=dlen(ips)/dble(nl)
c
          do izs=izs1,izs2
            if(izs.eq.izs1)then
              if(izs.lt.nzs)then
                dpz=dmax1(0.d0,dmin1(pz2,0.5d0*(zs(izs)+zs(izs+1)))-pz1)
              else
                dpz=pz2-pz1
              endif
            else if(izs.eq.izs2)then
              dpz=dmax1(0.d0,pz2-dmax1(pz1,0.5d0*(zs(izs)+zs(izs-1))))
            else
              dz1=zs(izs)-dmax1(pz1,0.5d0*(zs(izs)+zs(izs-1)))
              dz2=dmin1(pz2,0.5d0*(zs(izs)+zs(izs+1)))-zs(izs)
              dpz=dmax1(0.d0,dz1)+dmax1(0.d0,dz2)
            endif
            yp=(zs(izs)-pz(ips))/ssdi(ips)
            dwei=(ddl*dpz/ssdi(ips))/parea(ips)
            do il=1,nl
              xp=-0.5d0*dwid(ips)+(dble(il)-0.5d0)*ddl
              px=xp*csst(ips)-ssst(ips)*yp*csdi(ips)
              py=xp*ssst(ips)+csst(ips)*yp*csdi(ips)
c
              xobs=xobs0-px
              yobs=yobs0-py
c
              dis=dsqrt(xobs**2+yobs**2)
              if(dis.gt.0.d0)then
                azi=datan2(yobs,xobs)
              else
                azi=0.d0
              endif
              ir1=0
              do ir=1,nr
                if(dis.ge.r(ir))then
                  ir1=ir
                endif
              enddo
              ir2=ir1+1
              if(ir1.lt.1.or.ir2.gt.nr)then
                stop ' sdmgrn: distance range exceeded!'
              endif
c
              w1=dwei*(r(ir2)-dis)/(r(ir2)-r(ir1))
              w2=dwei*(dis-r(ir1))/(r(ir2)-r(ir1))
c
              co=dcos(azi)
              si=dsin(azi)
              co2=dcos(2.d0*azi)
              si2=dsin(2.d0*azi)
c
              do ira=1,2
c
c               add differential Green functions
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c               pmwei(1-5):
c               1 = weight for strike-slip: m12=m21=1;
c                   poloidal *sin(2*theta), toroidal *cos(2*theta)
c
c               2 = weight for dip-slip: m13=m31=1
c                   poloidal *cos(theta), toroidal *sin(theta)
c
c               3 = weight for clvd: m33=-m11=-m22=1
c                   axisymmetric
c
c               4 = weight for 45 deg strike-slip: m11=-m22=1
c                   greenfct4(theta) = green1(theta + 45 deg)
c
c               5 = weight for 45 deg dip-slip: m23=m32=1
c                   greenfct5(theta) = green2(theta - 90 deg)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                psss=pmwei(1,ips,ira)*si2+pmwei(4,ips,ira)*co2
                shss=pmwei(1,ips,ira)*co2-pmwei(4,ips,ira)*si2
                psds=pmwei(2,ips,ira)*co+pmwei(5,ips,ira)*si
                shds=pmwei(2,ips,ira)*si-pmwei(5,ips,ira)*co
                pscl=pmwei(3,ips,ira)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c               contributions from the strike-slip components
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                duz=psss
     &             *(w1*dgrns(ir1,izs,1,1)+w2*dgrns(ir2,izs,1,1))
                dur=psss
     &             *(w1*dgrns(ir1,izs,2,1)+w2*dgrns(ir2,izs,2,1))
                dut=shss
     &             *(w1*dgrns(ir1,izs,3,1)+w2*dgrns(ir2,izs,3,1))
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c               contributions from the dip-slip components
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                duz=duz+psds
     &             *(w1*dgrns(ir1,izs,1,2)+w2*dgrns(ir2,izs,1,2))
                dur=dur+psds
     &             *(w1*dgrns(ir1,izs,2,2)+w2*dgrns(ir2,izs,2,2))
                dut=dut+shds
     &             *(w1*dgrns(ir1,izs,3,2)+w2*dgrns(ir2,izs,3,2))
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c               contributions from the clvd components
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                duz=duz+pscl
     &             *(w1*dgrns(ir1,izs,1,3)+w2*dgrns(ir2,izs,1,3))
                dur=dur+pscl
     &             *(w1*dgrns(ir1,izs,2,3)+w2*dgrns(ir2,izs,2,3))
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                dux0(ira)=dux0(ira)+dur*co-dut*si
                duy0(ira)=duy0(ira)+dur*si+dut*co
                duz0(ira)=duz0(ira)+duz
              enddo
            enddo
          enddo
          do ira=1,2
            dur=dux0(ira)*co0+duy0(ira)*si0
            dut=duy0(ira)*co0-dux0(ira)*si0
            duz=duz0(ira)
c
            dux=dur*bco-dut*bsi
            duy=dur*bsi+dut*bco
c
            datgrn(ira,ips,iobs)=datgrn(ira,ips,iobs)
     &                          +dux*xcs(iobs)+duy*ycs(iobs)
     &                          +duz*zcs(iobs)
          enddo
        enddo
        if(mod(iobs,np10).eq.0)then
          npercent=10*(iobs/np10)
          write(*,'(a,i3,a)')'    - - - ',npercent,
     &                   '% of the diff. Green functions calculated ...'
        endif
      enddo
c
500   continue
c
      if(nlarge.gt.0)then
        nwarn=nwarn+nlarge
        write(*,'(a,i5,a)')' Warning: ',nlarge,' too large distances'
     &                   //' exceed the Green-function coverage!'
      endif
c
      deallocate(csst,ssst,cs2st,ss2st,csdi,ssdi,rnn,rst,rdi)
c
      return
      end
