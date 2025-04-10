      subroutine sdmdisc(ns,nps)
      implicit none
c
c     Last modified: Potsdam, Oct, 2008, by R. Wang
c
      integer ns,nps
c
      include 'sdmglob.h'
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     LOCAL WORK SPACES
c     =================
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer i,is,iw,il,ips,jps,ira,iwft,nlen,nwid
      double precision lat0,lon0,st,st0,di,dl,dw,pn,pe,wid
      double precision d0,dp0,dp,xp,yp,zp,dx,dy,lenpp
      double precision dx1,dx2,dy1,dy2,dz1,dz2,dnx,dny,dnz,hdis,vdis
      double precision x0,y0,bga,bgc,sma,smb,smc,alf,beta,d1,d2,dd
      double precision ra(2),sm(3,3)
      double precision xptop(2),yptop(2),xpbtm(2),ypbtm(2)
      double precision lenp(2),widp(2),depp(2),ditop(2),dibtm(2)
      double precision x(2,NPSMAX+1),y(2,NPSMAX+1),z(2,NPSMAX+1)
      character*180 header
      logical readok,bigdip
c
      double precision fdip
c
      double precision PI,DIPEPS
      data PI,DIPEPS/3.14159265358979d0,1.0d0/
c
      ra(1)=0.d0
      ra(2)=90.d0*DEG2RAD
c
      if(idisc.eq.0)then
        nps1(1)=1
        nps2(1)=0
        nps=0
        do is=1,ns
          if(is.gt.1)then
            nps1(is)=nps2(is-1)+1
            nps2(is)=nps1(is)-1
          endif
          open(20,file=inpatches(is),status='old')
          read(20,'(a)')header
          do i=1,NPSMAX
            nps2(is)=nps2(is)+1
            nps=nps+1
            if(nps.gt.NPSMAX)then
              stop 'NPSMAX defined too small'
            endif
            readok=.false.
            read(20,*,end=10)plat(nps),plon(nps),pz(nps),
     &                       dlen(nps),dwid(nps),strike(nps),dip(nps)
            pz(nps)=pz(nps)*KM2M
            dlen(nps)=dlen(nps)*KM2M
            dwid(nps)=dwid(nps)*KM2M
            parea(nps)=dlen(nps)*dwid(nps)
            if(iref(is).ge.1.and.iref(is).le.4)then
              st=strike(nps)*DEG2RAD
              di=dip(nps)*DEG2RAD
              if(iref(is).eq.1)then
                pz(nps)=pz(nps)+0.5d0*dwid(nps)*dsin(di)
                pn=0.5d0*dlen(nps)*dcos(st)
     &            -0.5d0*dwid(nps)*dcos(di)*dsin(st)
                pe=0.5d0*dlen(nps)*dsin(st)
     &            +0.5d0*dwid(nps)*dcos(di)*dcos(st)
              else if(iref(is).eq.2)then
                pz(nps)=pz(nps)+0.5d0*dwid(nps)*dsin(di)
                pn=-0.5d0*dlen(nps)*dcos(st)
     &             -0.5d0*dwid(nps)*dcos(di)*dsin(st)
                pe=-0.5d0*dlen(nps)*dsin(st)
     &             +0.5d0*dwid(nps)*dcos(di)*dcos(st)
              else if(iref(is).eq.3)then
                pz(nps)=pz(nps)-0.5d0*dwid(nps)*dsin(di)
                pn=0.5d0*dlen(nps)*dcos(st)
     &            +0.5d0*dwid(nps)*dcos(di)*dsin(st)
                pe=0.5d0*dlen(nps)*dsin(st)
     &            -0.5d0*dwid(nps)*dcos(di)*dcos(st)
              else
                pz(nps)=pz(nps)-0.5d0*dwid(nps)*dsin(di)
                pn=-0.5d0*dlen(nps)*dcos(st)
     &             +0.5d0*dwid(nps)*dcos(di)*dsin(st)
                pe=-0.5d0*dlen(nps)*dsin(st)
     &             -0.5d0*dwid(nps)*dcos(di)*dcos(st)
              endif
c
c             determine central point of the patch
c
c             spherical triangle:
c             A = pole, B = source position, C = reference position
c
              sma=dsqrt(pn**2+pe**2)/REARTH
              smb=0.5d0*PI-plat(nps)*DEG2RAD
              bgc=datan2(pe,pn)
              smc=dacos(dcos(sma)*dcos(smb)
     &           +dsin(sma)*dsin(smb)*dcos(bgc))
              bga=dasin(dsin(sma)*dsin(bgc)/dsin(smc))
c
c             geographic coordinate of the equivalent point source
c
              plat(nps)=90.d0-smc/DEG2RAD
              plon(nps)=dmod(plon(nps)+bga/DEG2RAD,360.d0)
c
            endif
c
            do ira=1,2
              sm(1,1)=-dsin(di)*dcos(ra(ira))*dsin(2.d0*st)
     &               -dsin(2.d0*di)*dsin(ra(ira))*(dsin(st))**2
              sm(2,2)= dsin(di)*dcos(ra(ira))*sin(2.d0*st)
     &               -dsin(2.d0*di)*dsin(ra(ira))*(dcos(st))**2
              sm(3,3)=-(sm(1,1)+sm(2,2))
              sm(1,2)= dsin(di)*dcos(ra(ira))*dcos(2.d0*st)
     &               +0.5d0*dsin(2.d0*di)*dsin(ra(ira))*dsin(2.d0*st)
              sm(2,1)=sm(1,2)
              sm(1,3)=-dcos(di)*dcos(ra(ira))*dcos(st)
     &              -dcos(2.d0*di)*dsin(ra(ira))*dsin(st)
              sm(3,1)=sm(1,3)
              sm(2,3)=-dcos(di)*dcos(ra(ira))*dsin(st)
     &               +dcos(2.d0*di)*dsin(ra(ira))*dcos(st)
              sm(3,2)=sm(2,3)
c
c             1 = weight for strike-slip: m12=m21=1;
c             2 = weight for dip-slip: m13=m31=1
c             3 = weight for clvd: m33=-m11=-m22=1
c             4 = weight for 45 deg strike-slip: m11=-m22=1
c             5 = weight for 45 deg dip-slip: m23=m32=1
c
              pmwei(1,nps,ira)=sm(1,2)*parea(nps)
              pmwei(2,nps,ira)=sm(1,3)*parea(nps)
              pmwei(3,nps,ira)=sm(3,3)*parea(nps)
              pmwei(4,nps,ira)=0.5d0*(sm(1,1)-sm(2,2))*parea(nps)
              pmwei(5,nps,ira)=sm(2,3)*parea(nps)
            enddo
            readok=.true.
          enddo
10        continue
          if(.not.readok)then
            nps2(is)=nps2(is)-1
            nps=nps-1
          endif
          write(*,'(a,i2,a,i6,a)')' the ',is,'. fault segment => ',
     &                         1+nps2(is)-nps1(is),' slip patches.'
          close(20)
        enddo
        goto 500
      endif
c
      nps1(1)=1
      nps2(1)=0
      nps=0
      do is=1,ns
        if(is.gt.1)then
          nps1(is)=nps2(is-1)+1
          nps2(is)=nps1(is)-1
        endif
        do iwft=1,nwft(is)-1
          lat0=toplat(iwft,is)
          lon0=toplon(iwft,is)
          if(topdip(iwft,is).gt.90.d0.or.
     &       topdip(iwft+1,is).gt.90.d0)then
            ditop(1)=(180.d0-topdip(iwft,is))*DEG2RAD
            ditop(2)=(180.d0-topdip(iwft+1,is))*DEG2RAD
            bigdip=.true.
          else
            ditop(1)=topdip(iwft,is)*DEG2RAD
            ditop(2)=topdip(iwft+1,is)*DEG2RAD
            bigdip=.false.
          endif
c
          xptop(1)=0.d0
          yptop(1)=0.d0
c
          call disazi(REARTH,lat0,lon0,
     &         toplat(iwft+1,is),toplon(iwft+1,is),xptop(2),yptop(2))
c
          call disazi(REARTH,lat0,lon0,
     &         btmlat(iwft,is),btmlon(iwft,is),xpbtm(1),ypbtm(1))
c
          call disazi(REARTH,lat0,lon0,
     &         btmlat(iwft+1,is),btmlon(iwft+1,is),xpbtm(2),ypbtm(2))
          do i=1,2
            lenp(i)=dsqrt((xpbtm(i)-xptop(i))**2
     &                   +(ypbtm(i)-yptop(i))**2)
            depp(i)=btmdep(is)-topdep(is)
            dibtm(i)=fdip(ditop(i)/DEG2RAD,lenp(i),depp(i))*DEG2RAD
c
            if(dabs(dibtm(i)-ditop(i)).le.DIPEPS)then
              widp(i)=depp(i)/dsin(0.5d0*(dibtm(i)+ditop(i)))
              ditop(i)=0.5d0*(dibtm(i)+ditop(i))
              dibtm(i)=ditop(i)
            else
              widp(i)=(depp(i)/(dibtm(i)-ditop(i)))
     &             *dlog(dtan(0.5d0*dibtm(i))/dtan(0.5d0*ditop(i)))
            endif
          enddo
c
          if(ndgrid(is).le.0)then
            nwid=2+idint(dmax1(widp(1),widp(2))/patchsize(is))
          else
            nwid=1+ndgrid(is)
          endif
c
          do i=1,2
            dw=widp(i)/dble(nwid-1)
            if(dabs(dibtm(i)-ditop(i)).le.DIPEPS)then
              dx=(xpbtm(i)-xptop(i))/dble(nwid-1)
              dy=(ypbtm(i)-yptop(i))/dble(nwid-1)
              do iw=1,nwid
                x(i,iw)=xptop(i)+dble(iw-1)*dx
                y(i,iw)=yptop(i)+dble(iw-1)*dy
                z(i,iw)=topdep(is)
     &                 +dble(iw-1)*dw*dsin(0.5d0*(dibtm(i)+ditop(i)))
              enddo
            else
              do iw=1,nwid
                wid=dble(iw-1)*dw
                zp=(2.d0*datan(dexp(wid*(dibtm(i)-ditop(i))/depp(i))
     &                      *dtan(0.5d0*ditop(i)))-ditop(i))
     &            *depp(i)/(dibtm(i)-ditop(i))
                lenpp=dlog(dsin(ditop(i)+zp*(dibtm(i)-ditop(i))/depp(i))
     &                /dsin(ditop(i)))*depp(i)/(dibtm(i)-ditop(i))
                x(i,iw)=xptop(i)+(xpbtm(i)-xptop(i))*lenpp/lenp(i)
                y(i,iw)=yptop(i)+(ypbtm(i)-yptop(i))*lenpp/lenp(i)
                z(i,iw)=topdep(is)+zp
              enddo
            endif
          enddo
c
c         determine patch parameters
c
          do iw=2,nwid
            nps=nps+1
            nps2(is)=nps2(is)+1
            if(nps.gt.NPSMAX)then
              stop ' Max. no. of point sources exceeded!'
            endif
            pn=0.25d0*(x(2,iw)+x(1,iw)+x(2,iw-1)+x(1,iw-1))
            pe=0.25d0*(y(2,iw)+y(1,iw)+y(2,iw-1)+y(1,iw-1))
            pz(nps)=0.25d0*(z(2,iw)+z(1,iw)+z(2,iw-1)+z(1,iw-1))
c
            strike(nps)=0.5d0*(datan2(y(2,iw-1)-y(1,iw-1),
     &                                x(2,iw-1)-x(1,iw-1))
     &                        +datan2(y(2,iw)-y(1,iw),
     &                                x(2,iw)-x(1,iw)))/DEG2RAD
            if(strike(nps).lt.0.d0)strike(nps)=strike(nps)+360.d0
c
c           determine two diagonal vectors
c
            dx1=x(2,iw-1)-x(1,iw)
            dy1=y(2,iw-1)-y(1,iw)
            dz1=z(2,iw-1)-z(1,iw)
            dx2=x(2,iw)-x(1,iw-1)
            dy2=y(2,iw)-y(1,iw-1)
            dz2=z(2,iw)-z(1,iw-1)
c
c           calculate cross product of the two diagonal vectors
c           amplitude = 2 times the area
c
            dnx=dy1*dz2-dy2*dz1
            dny=dx2*dz1-dx1*dz2
            dnz=dx1*dy2-dx2*dy1
c
            parea(nps)=0.5d0*dsqrt(dnx**2+dny**2+dnz**2)
c
            dip(nps)=dacos(0.5d0*dnz/parea(nps))/DEG2RAD
            if(bigdip)then
              dip(nps)=180.d0-dip(nps)
            endif
c
            dlen(nps)=0.5d0*(dsqrt((x(2,iw)-x(1,iw))**2
     &                            +(y(2,iw)-y(1,iw))**2)
     &                      +dsqrt((x(2,iw-1)-x(1,iw-1))**2
     &                            +(y(2,iw-1)-y(1,iw-1))**2))
            dwid(nps)=parea(nps)/dlen(nps)
c
c           convert to geographic coordinates
c
c           spherical triangle:
c             A = pole, B = source position, C = reference position
c
            sma=dsqrt(pn**2+pe**2)/REARTH
            smb=0.5d0*PI-lat0*DEG2RAD
            bgc=datan2(pe,pn)
            smc=dacos(dcos(sma)*dcos(smb)
     &         +dsin(sma)*dsin(smb)*dcos(bgc))
            bga=dasin(dsin(sma)*dsin(bgc)/dsin(smc))
c
c           geographic coordinate of the equivalent point source
c
            plat(nps)=90.d0-smc/DEG2RAD
            plon(nps)=dmod(lon0+bga/DEG2RAD,360.d0)
c
            st=strike(nps)*DEG2RAD
            di=dip(nps)*DEG2RAD
            if(pz(nps)-0.5d0*dwid(nps)*dsin(di).le.0.d0)then
c
c             upper patch edge exceeds the surface
c
              pz(nps)=0.5001d0*dwid(nps)*dsin(di)
            endif
            do ira=1,2
              sm(1,1)=-dsin(di)*dcos(ra(ira))*dsin(2.d0*st)
     &               -dsin(2.d0*di)*dsin(ra(ira))*(dsin(st))**2
              sm(2,2)= dsin(di)*dcos(ra(ira))*sin(2.d0*st)
     &               -dsin(2.d0*di)*dsin(ra(ira))*(dcos(st))**2
              sm(3,3)=-(sm(1,1)+sm(2,2))
              sm(1,2)= dsin(di)*dcos(ra(ira))*dcos(2.d0*st)
     &               +0.5d0*dsin(2.d0*di)*dsin(ra(ira))*dsin(2.d0*st)
              sm(2,1)=sm(1,2)
              sm(1,3)=-dcos(di)*dcos(ra(ira))*dcos(st)
     &              -dcos(2.d0*di)*dsin(ra(ira))*dsin(st)
              sm(3,1)=sm(1,3)
              sm(2,3)=-dcos(di)*dcos(ra(ira))*dsin(st)
     &               +dcos(2.d0*di)*dsin(ra(ira))*dcos(st)
              sm(3,2)=sm(2,3)
c
c             1 = weight for strike-slip: m12=m21=1;
c             2 = weight for dip-slip: m13=m31=1
c             3 = weight for clvd: m33=-m11=-m22=1
c             4 = weight for 45 deg strike-slip: m11=-m22=1
c             5 = weight for 45 deg dip-slip: m23=m32=1
c
              pmwei(1,nps,ira)=sm(1,2)*parea(nps)
              pmwei(2,nps,ira)=sm(1,3)*parea(nps)
              pmwei(3,nps,ira)=sm(3,3)*parea(nps)
              pmwei(4,nps,ira)=0.5d0*(sm(1,1)-sm(2,2))*parea(nps)
              pmwei(5,nps,ira)=sm(2,3)*parea(nps)
            enddo
          enddo
        enddo
        write(*,'(a,i2,a,i6,a)')' the ',is,'. fault segment => ',
     &                         1+nps2(is)-nps1(is),' slip patches.'
      enddo
500   continue
      write(*,*)'------------------------------------------------'
      write(*,'(a,i7)')' total number of point sources: ',nps
c
      sarea=0.d0
      do ips=1,nps
        sarea=sarea+parea(ips)
      enddo
c
      do ips=1,nps
        st=strike(ips)*DEG2RAD
        di=dip(ips)*DEG2RAD
        d0=dsqrt(dlen(ips)**2+dwid(ips)**2)
c
c       search left neighboring patch
c
        dp=d0
        ipsl(ips)=0
        xp=-dlen(ips)*dcos(st)
        yp=-dlen(ips)*dsin(st)
        zp=pz(ips)
c
        do jps=1,nps
          if(jps.ne.ips.and.dabs(pz(jps)-pz(ips)).le.
     &       dwid(ips)*dsin(di)+dwid(jps)*dsin(dip(jps)*DEG2RAD))then
            call disazi(rearth,plat(ips),plon(ips),
     &                         plat(jps),plon(jps),dx,dy)
            dp0=dsqrt((dx-xp)**2+(dy-yp)**2+(pz(jps)-zp)**2)
            if(dp.ge.dp0)then
              dp=dp0
              ipsl(ips)=jps
            endif
          endif
        enddo
c
c       search right neighboring patch
c
        dp=d0
        ipsr(ips)=0
        xp=dlen(ips)*dcos(st)
        yp=dlen(ips)*dsin(st)
        zp=pz(ips)
c
        do jps=1,nps
          if(jps.ne.ips.and.dabs(pz(jps)-pz(ips)).le.
     &       dwid(ips)*dsin(di)+dwid(jps)*dsin(dip(jps)*DEG2RAD))then
            call disazi(rearth,plat(ips),plon(ips),
     &                         plat(jps),plon(jps),dx,dy)
            dp0=dsqrt((dx-xp)**2+(dy-yp)**2+(pz(jps)-zp)**2)
            if(dp.ge.dp0)then
              dp=dp0
              ipsr(ips)=jps
            endif
          endif
        enddo
c
c       search upper neighboring patch
c
        if(pz(ips).le.dwid(ips)*dsin(di))then
          ipsu(ips)=-1
        else
          dp=d0
          ipsu(ips)=0
          xp= dwid(ips)*dcos(di)*dsin(st)
          yp=-dwid(ips)*dcos(di)*dcos(st)
          zp=pz(ips)-dwid(ips)*dsin(di)
c
          do jps=1,nps
            if(jps.ne.ips.and.pz(jps).lt.pz(ips))then
              call disazi(rearth,plat(ips),plon(ips),
     &                           plat(jps),plon(jps),dx,dy)
              dp0=dsqrt((dx-xp)**2+(dy-yp)**2+(pz(jps)-zp)**2)
              if(dp.ge.dp0)then
                dp=dp0
                ipsu(ips)=jps
              endif
            endif
          enddo
        endif
c
c       search lower neighboring patch
c
        dp=d0
        ipsd(ips)=0
        xp=-dwid(ips)*dcos(di)*dsin(st)
        yp= dwid(ips)*dcos(di)*dcos(st)
        zp=pz(ips)+dwid(ips)*dsin(di)
c
        do jps=1,nps
          if(jps.ne.ips.and.pz(jps).gt.pz(ips))then
            call disazi(rearth,plat(ips),plon(ips),
     &                         plat(jps),plon(jps),dx,dy)
            dp0=dsqrt((dx-xp)**2+(dy-yp)**2+(pz(jps)-zp)**2)
            if(dp.ge.dp0)then
              dp=dp0
              ipsd(ips)=jps
            endif
          endif
        enddo
      enddo
      return
      end
