      subroutine sdmgrn(ns,nps,nobs,ismooth)
      implicit none
c
c     Last modified: Potsdam, Oct, 2008, by R. Wang
c
      integer*4 ns,nps,nobs,ismooth
c
      include 'sdmglob.h'
c
      integer*4 is,ips,jps,iobs,nlarge,ira,il,nl
      integer*4 i,j,ir,ir1,ir2,izs,izs1,izs2,idiv
      integer*4 np5,ipsum,npercent
      real*8 dr,ddl,st,di,dux,duy,duz,dur,dut
      real*8 xobs,yobs,dobs,dal,daw
      real*8 eii,exx,eyy,ezz,exy,eyz,ezx
      real*8 si,co,si2,co2,dis,azi
      real*8 w1,w2,psss,psds,pscl,shss,shds
      real*8 dz1,dz2,dpz,dwei,pz1,pz2,xp,yp,px,py
      real*8 strst(0:4),strdi(0:4),strnn(0:4),dl(4),dw(4)
      real*8 csst(NPSMAX),ssst(NPSMAX)
      real*8 cs2st(NPSMAX),ss2st(NPSMAX)
      real*8 csdi(NPSMAX),ssdi(NPSMAX)
      real*8 rnn(3,NPSMAX),rst(3,NPSMAX),rdi(3,NPSMAX)
      real*8 sig(3,3),stress(3)
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
      nlarge=0
c
c     SUPERPOSITION OF ALL DISCRETE POINT SOURCES
c     ===========================================
c
      do ira=1,2
        do iobs=1,nobs
          do ips=1,nps
            dspmdl(ips,iobs,ira)=0.d0
          enddo
        enddo
      enddo
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
c       cosine vectors of normal of hanging wall
c
        rnn(1,ips)= ssdi(ips)*ssst(ips)
        rnn(2,ips)=-ssdi(ips)*csst(ips)
        rnn(3,ips)= csdi(ips)
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
      np5=nps*5
c
      ipsum=0
      write(*,'(a)')' analytical Okada solutions - please wait: '
      do ips=1,nps
c
        DEPTH=sngl(pz(ips))
        DIPS=sngl(dip(ips))
c
        AL1=sngl(-0.5d0*dlen(ips))
        AL2=sngl(+0.5d0*dlen(ips))
        AW1=sngl(-0.4999d0*dwid(ips))
        AW2=sngl(+0.4999d0*dwid(ips))
c
c       displacement Green functions
c
        do iobs=1,nobs
c
c         transform from Aki's to Okada's system
c
          call disazi(REARTH,plat(ips),plon(ips),
     &                latobs(iobs),lonobs(iobs),xobs,yobs)
          X=sngl(xobs*csst(ips)+yobs*ssst(ips))
          Y=sngl(xobs*ssst(ips)-yobs*csst(ips))
          Z=-sngl(zobs)
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
            dux=dble(UX)*csst(ips)+dble(UY)*ssst(ips)
            duy=dble(UX)*ssst(ips)-dble(UY)*csst(ips)
            duz=-dble(UZ)
c
            dspmdl(ips,iobs,ira)=dspmdl(ips,iobs,ira)
     &            +dux*xcs(iobs)+duy*ycs(iobs)+duz*zcs(iobs)
          enddo
        enddo
c
c       stress drop Green's functions
c       strgrn(ips,i,jps,j) = contribution of i-th slip component
c       at patch ips to j-th stress-drop component at patch jps
c
        do jps=1,nps
c
c         transform from Aki's to Okada's system
c
          call disazi(REARTH,plat(ips),plon(ips),
     &                  plat(jps),plon(jps),xobs,yobs)
          X=sngl(xobs*csst(ips)+yobs*ssst(ips))
          Y=sngl(xobs*ssst(ips)-yobs*csst(ips))
          Z=-sngl(pz(jps))
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
            sig(1,1)=LAMREF*eii+2.d0*MUEREF*exx
            sig(2,2)=LAMREF*eii+2.d0*MUEREF*eyy
            sig(3,3)=LAMREF*eii+2.d0*MUEREF*ezz
            sig(1,2)=2.d0*MUEREF*exy
            sig(2,3)=2.d0*MUEREF*eyz
            sig(3,1)=2.d0*MUEREF*ezx
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
            strst(1)=0.d0
            strdi(1)=0.d0
            strnn(1)=0.d0
            do i=1,3
              strst(1)=strst(1)+stress(i)*rst(i,jps)
              strdi(1)=strdi(1)+stress(i)*rdi(i,jps)
              strnn(1)=strnn(1)+stress(i)*rnn(i,jps)
            enddo
            strgrn(ips,ira,jps,1)=strst(1)
            strgrn(ips,ira,jps,2)=strdi(1)
            strgrn(ips,ira,jps,3)=strnn(1)
          enddo
        enddo
c
        if(ismooth.ne.2)then
          do jps=1,nps
            dal=1.d0/dlen(jps)**2
            daw=1.d0/dwid(jps)**2
            do i=0,4
              strst(i)=0.d0
              strdi(i)=0.d0
            enddo
            if(ips.eq.jps.and.
     &         ipsu(jps).gt.0.and.ipsd(jps).gt.0.and.
     &         ipsl(jps).gt.0.and.ipsr(jps).gt.0)then
              strst(0)=1.d0
              strdi(0)=1.d0
            endif
            if(ips.eq.ipsl(jps).and.ipsr(jps).gt.0)then
              strst(1)=1.d0
              strdi(1)=1.d0
            endif
            if(ips.eq.ipsr(jps).and.ipsl(jps).gt.0)then
              strst(3)=1.d0
              strdi(3)=1.d0
            endif
            if(ips.eq.ipsu(jps).and.ipsd(jps).gt.0)then
              strst(2)=1.d0
              strdi(2)=1.d0
            endif
            if(ips.eq.ipsd(jps).and.ipsu(jps).gt.0)then
              strst(4)=1.d0
              strdi(4)=1.d0
            endif
            dcgrn(ips,1,jps,1)=(strst(1)-2.d0*strst(0)+strst(3))*dal
            dcgrn(ips,1,jps,2)=(strst(2)-2.d0*strst(0)+strst(4))*daw
            dcgrn(ips,1,jps,3)=0.d0
            dcgrn(ips,1,jps,4)=0.d0
            dcgrn(ips,1,jps,5)=0.d0
            dcgrn(ips,1,jps,6)=0.d0
c
            dcgrn(ips,2,jps,1)=0.d0
            dcgrn(ips,2,jps,2)=0.d0
            dcgrn(ips,2,jps,3)=(strdi(1)-2.d0*strdi(0)+strdi(3))*dal
            dcgrn(ips,2,jps,4)=(strdi(2)-2.d0*strdi(0)+strdi(4))*daw
            dcgrn(ips,2,jps,5)=0.d0
            dcgrn(ips,2,jps,6)=0.d0
          enddo
        else
          do jps=1,nps
            dal=16.d0/dlen(jps)**2
            daw=16.d0/dwid(jps)**2
c
            dl(1)=0.25d0*dlen(jps)
            dw(1)=0.d0
c
            dl(2)=0.d0
            dw(2)=0.25d0*dwid(jps)
c
            dl(3)=-0.25d0*dlen(jps)
            dw(3)=0.d0
c
            dl(4)=0.d0
            dw(4)=-0.25d0*dwid(jps)
            do ira=1,2
              DISL1=SS(ira)
              DISL2=DS(ira)
              strst(0)=strgrn(ips,ira,jps,1)
              strdi(0)=strgrn(ips,ira,jps,2)
              strnn(0)=strgrn(ips,ira,jps,3)
              do idiv=1,4
                dobs=pz(jps)+dw(idiv)*ssdi(jps)
                if(dobs.gt.0.d0)then
                  call disazi(REARTH,plat(ips),plon(ips),
     &                        plat(jps),plon(jps),xobs,yobs)
                  xobs=xobs
     &                +dl(idiv)*csst(jps)-dw(idiv)*csdi(jps)*ssst(jps)
                  yobs=yobs
     &                +dl(idiv)*ssst(jps)+dw(idiv)*csdi(jps)*csst(jps)
                  X=sngl(xobs*csst(ips)+yobs*ssst(ips))
                  Y=sngl(xobs*ssst(ips)-yobs*csst(ips))
                  Z=-sngl(dobs)
                  IRET=1
                  call DC3D(ALPHA,X,Y,Z,DEPTH,DIPS,AL1,AL2,AW1,AW2,
     &                      DISL1,DISL2,DISL3,UX,UY,UZ,
     &                      UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ,IRET)
c
c                 transform from Okada's to Aki's system
c
                  exx=dble(UXX)*csst(ips)*csst(ips)
     &               +dble(UYY)*ssst(ips)*ssst(ips)
     &               +0.5d0*dble(UXY+UYX)*ss2st(ips)
                  eyy=dble(UXX)*ssst(ips)*ssst(ips)
     &               +dble(UYY)*csst(ips)*csst(ips)
     &               -0.5d0*dble(UXY+UYX)*ss2st(ips)
                  ezz=dble(UZZ)
                  exy=0.5d0*dble(UXX-UYY)*ss2st(ips)
     &               -0.5d0*dble(UXY+UYX)*cs2st(ips)
                  eyz=-0.5d0*dble(UZX+UXZ)*ssst(ips)
     &                +0.5d0*dble(UYZ+UZY)*csst(ips)
                  ezx=-0.5d0*dble(UZX+UXZ)*csst(ips)
     &                  -0.5d0*dble(UYZ+UZY)*ssst(ips)
                  eii=exx+eyy+ezz
c
                  sig(1,1)=LAMREF*eii+2.d0*MUEREF*exx
                  sig(2,2)=LAMREF*eii+2.d0*MUEREF*eyy
                  sig(3,3)=LAMREF*eii+2.d0*MUEREF*ezz
                  sig(1,2)=2.d0*MUEREF*exy
                  sig(2,3)=2.d0*MUEREF*eyz
                  sig(3,1)=2.d0*MUEREF*ezx
                  sig(2,1)=sig(1,2)
                  sig(3,2)=sig(2,3)
                  sig(1,3)=sig(3,1)
c
c                 stress drop at jps-th patch
c
                  do i=1,3
                    stress(i)=0.d0
                    do j=1,3
                      stress(i)=stress(i)+sig(i,j)*rnn(j,jps)
                    enddo
                  enddo
c
                  strst(idiv)=0.d0
                  strdi(idiv)=0.d0
                  strnn(idiv)=0.d0
                  do i=1,3
                    strst(idiv)=strst(idiv)+stress(i)*rst(i,jps)
                    strdi(idiv)=strdi(idiv)+stress(i)*rdi(i,jps)
                    strnn(idiv)=strnn(idiv)+stress(i)*rnn(i,jps)
                  enddo
                else
                  strst(idiv)=0.d0
                  strdi(idiv)=0.d0
                endif
              enddo
              dcgrn(ips,ira,jps,1)=(strst(1)-2.d0*strst(0)+strst(3))*dal
              dcgrn(ips,ira,jps,2)=(strst(2)-2.d0*strst(0)+strst(4))*daw
              dcgrn(ips,ira,jps,3)=(strdi(1)-2.d0*strdi(0)+strdi(3))*dal
              dcgrn(ips,ira,jps,4)=(strdi(2)-2.d0*strdi(0)+strdi(4))*daw
              dcgrn(ips,ira,jps,5)=(strnn(1)-2.d0*strnn(0)+strnn(3))*dal
              dcgrn(ips,ira,jps,6)=(strnn(2)-2.d0*strnn(0)+strnn(4))*daw
            enddo
          enddo
        endif
c
        if(ips.gt.0.and.nps.ge.20.and.mod(100*ips,np5).eq.0)then
          npercent=500*ips/np5
          write(*,'(a,i3,a)')'    - - - ',npercent,
     &                       '% calculated ...'
        endif
      enddo
c
      if(hsmodel)goto 500
c
      dr=r(2)-r(1)
      ipsum=0
      write(*,'(a)')' numerical layering effect - please wait: '
      do ips=1,nps
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
            do iobs=1,nobs
c
c             transform from Aki's to Okada's system
c
              call disazi(REARTH,plat(ips),plon(ips),
     &                latobs(iobs),lonobs(iobs),xobs,yobs)
              xobs=xobs-px
              yobs=yobs-py
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
     &             *(w1*grns(izs,ir1,1,1)+w2*grns(izs,ir2,1,1))
                dur=psss
     &             *(w1*grns(izs,ir1,2,1)+w2*grns(izs,ir2,2,1))
                dut=shss
     &             *(w1*grns(izs,ir1,3,1)+w2*grns(izs,ir2,3,1))
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c               contributions from the dip-slip components
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                duz=duz+psds
     &             *(w1*grns(izs,ir1,1,2)+w2*grns(izs,ir2,1,2))
                dur=dur+psds
     &             *(w1*grns(izs,ir1,2,2)+w2*grns(izs,ir2,2,2))
                dut=dut+shds
     &             *(w1*grns(izs,ir1,3,2)+w2*grns(izs,ir2,3,2))
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c               contributions from the clvd components
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                duz=duz+pscl
     &             *(w1*grns(izs,ir1,1,3)+w2*grns(izs,ir2,1,3))
                dur=dur+pscl
     &             *(w1*grns(izs,ir1,2,3)+w2*grns(izs,ir2,2,3))
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                dspmdl(ips,iobs,ira)=dspmdl(ips,iobs,ira)
     &                          +(dur*co-dut*si)*xcs(iobs)
     &                          +(dur*si+dut*co)*ycs(iobs)
     &                          +duz*zcs(iobs)
              enddo
            enddo
          enddo
        enddo
        if(ips.gt.0.and.nps.ge.20.and.mod(100*ips,np5).eq.0)then
          npercent=500*ips/np5
          write(*,'(a,i3,a)')'    - - - ',npercent,
     &                       '% of the Green functions calculated ...'
        endif
      enddo
c
c     for using Zhang's weighted smoothing of roughness
c
500   continue
c
      if(nlarge.gt.0)then
        nwarn=nwarn+nlarge
        write(*,'(a,i5,a)')' Warning: ',nlarge,' too large distances'
     &                   //' exceed the Green-function coverage!'
      endif
      return
      end
