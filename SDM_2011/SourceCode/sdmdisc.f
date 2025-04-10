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
      integer i,is,iw,il,ips,jps,ira
      double precision lat0,lon0,st,st0,di,dl,dw,pn,pe
      double precision dx,dy,xp,yp,zp,d0,dp0,dp
      double precision dx1,dx2,dy1,dy2,dz1,dz2,dnx,dny,dnz,hdis,vdis
      double precision x0,y0,bga,bgc,sma,smb,smc,alf,beta,d1,d2,dd
      double precision ra(2),sm(3,3)
      double precision x(NPSMAX+1,NPSMAX+1),y(NPSMAX+1,NPSMAX+1)
      double precision z(NPSMAX+1,NPSMAX+1)
      double precision dip1(NPSMAX+1),dip2(NPSMAX+1)
      double precision lft(NFTMAX),xft(NFTMAX),yft(NFTMAX)
c
      double precision PI,EPSDIP
      data PI,EPSDIP/3.14159265358979d0,0.1d0/
c
      ra(1)=0.d0
      ra(2)=90.d0*DEG2RAD
c
      nps=0
      do is=1,ns
        lat0=0.5d0*(latft(1,is)+latft(nft(is),is))
        lon0=0.5d0*(lonft(1,is)+lonft(nft(is),is))
c
c       local cartesian coordinates of surface ruptures
c
        do i=1,nft(is)
          call disazi(REARTH,lat0,lon0,latft(i,is),lonft(i,is),
     &                xft(i),yft(i))
        enddo
c
c       accumulated length of top fault trace
c
        lft(1)=0.d0
        do i=2,nft(is)
          lft(i)=lft(i-1)+dsqrt((xft(i)-xft(i-1))**2
     &                         +(yft(i)-yft(i-1))**2)
        enddo
        write(*,'(a,i4,a,f6.2,a)')' Length along the curved strike ',
     &       is,'. fault segment: ',
     &       lft(nft(is))/KM2M,' km'
c
c       length sampling parameters
c
        nlength(is)=max0(2,1+idnint(lft(nft(is))/patchsize(is)))
        if(nlength(is).ge.NPSMAX)then
          stop 'Error in sdmdisc: NPSMAX too small!'
        endif
        dl=lft(nft(is))/dble(nlength(is)-1)
c
c       discretise top fault trace
c       x(il,iw), y(il,iw) = cartesian coordinates
c       il = index increasing with length (strike)
c       iw = index increasing with width (down dip)
c
        x(1,1)=xft(1)
        y(1,1)=yft(1)
        z(1,1)=topdep(is)
        dip1(1)=topdip(1,is)
        dip2(1)=botdip(1,is)
        x(nlength(is),1)=xft(nft(is))
        y(nlength(is),1)=yft(nft(is))
        z(nlength(is),1)=topdep(is)
        dip1(nlength(is))=topdip(nft(is),is)
        dip2(nlength(is))=botdip(nft(is),is)
c
        do il=2,nlength(is)-1
          z(il,1)=topdep(is)
          do i=2,nft(is)
            if(lft(i).ge.dble(il-1)*dl)then
              beta=(dble(il-1)*dl-lft(i-1))/(lft(i)-lft(i-1))
              x(il,1)=xft(i-1)+(xft(i)-xft(i-1))*beta
              y(il,1)=yft(i-1)+(yft(i)-yft(i-1))*beta
              dip1(il)=topdip(i-1,is)+(topdip(i,is)-topdip(i-1,is))*beta
              dip2(il)=botdip(i-1,is)+(botdip(i,is)-botdip(i-1,is))*beta
              goto 100
            endif
          enddo
100       continue
        enddo
c
c       discretise whole fault plane
c
        nwidth(is)=max0(2,1+idnint(width(is)/patchsize(is)))
        if(nwidth(is).ge.NPSMAX)then
          stop 'Error in sdmdisc: NPSMAX too small!'
        endif
        dw=width(is)/dble(nwidth(is)-1)
c
c       --- dipping axis
c
        st0=mstrike(is)*DEG2RAD
        if(mstrike(is).lt.0.d0.or.mstrike(is).gt.360.d0)then
          x0=0.d0
          y0=0.d0
          do il=1,nlength(is)
            x0=x0+x(il,1)
            y0=y0+y(il,1)
          enddo
          x0=x0/dble(nlength(is))
          y0=y0/dble(nlength(is))
          d1=0.d0
          d2=0.d0
          dd=0.d0
          do il=1,nlength(is)
            d1=d1+(x(il,1)-x0)**2
            d2=d2+(y(il,1)-y0)**2
            dd=dd+(x(il,1)-x0)*(y(il,1)-y0)
          enddo
          if(d1.ge.d2)then
            alf=dd/d1
            st0=datan2(alf*(x(nlength(is),1)-x(1,1)),
     &                 x(nlength(is),1)-x(1,1))
          else
            alf=dd/d2
            st0=datan2(y(nlength(is),1)-y(1,1),
     &                 alf*(y(nlength(is),1)-y(1,1)))
          endif
          write(*,'(a,i4,a,f6.2,a)')' Average strike of ',
     &         is,'. fault segment: ',
     &         dmod(st0+2.d0*PI,2.d0*PI)/DEG2RAD,' deg'
        endif
c
c       --- cartesian coordinates of nodes at curved fault plane
c
        do il=1,nlength(is)
          d1=dip1(il)*DEG2RAD
          d2=dip2(il)*DEG2RAD
          if(dabs(d2-d1).le.EPSDIP)then
            do iw=2,nwidth(is)
              hdis=dble(iw-1)*dw*dcos(0.5d0*(d1+d2))
              vdis=dble(iw-1)*dw*dsin(0.5d0*(d1+d2))
              x(il,iw)=x(il,1)-hdis*dsin(st0)
              y(il,iw)=y(il,1)+hdis*dcos(st0)
              z(il,iw)=z(il,1)+vdis
            enddo
          else
            dd=(dip2(il)-dip1(il))*DEG2RAD/width(is)
            do iw=2,nwidth(is)
              hdis=(dsin(d1+dd*dble(iw-1)*dw)-dsin(d1))/dd
              vdis=(dcos(d1)-dcos(d1+dd*dble(iw-1)*dw))/dd
              x(il,iw)=x(il,1)-hdis*dsin(st0)
              y(il,iw)=y(il,1)+hdis*dcos(st0)
              z(il,iw)=z(il,1)+vdis
            enddo
          endif
        enddo
c
c       determine patch parameters
c
        ips=0
        do il=2,nlength(is)
          do iw=2,nwidth(is)
            ips=ips+1
            nps=nps+1
            if(nps.gt.NPSMAX)then
              stop ' Max. no. of point sources exceeded!'
            endif
            pl(nps)=(dble(il)-1.5d0)*dl
            pw(nps)=(dble(iw)-1.5d0)*dw
            pn=0.25d0*(x(il,iw)+x(il-1,iw)+x(il,iw-1)+x(il-1,iw-1))
            pe=0.25d0*(y(il,iw)+y(il-1,iw)+y(il,iw-1)+y(il-1,iw-1))
            pz(nps)=0.25d0*(z(il,iw)+z(il-1,iw)+z(il,iw-1)+z(il-1,iw-1))
c
            strike(nps)=0.5d0*(datan2(y(il,iw-1)-y(il-1,iw-1),
     &                                x(il,iw-1)-x(il-1,iw-1))
     &                        +datan2(y(il,iw)-y(il-1,iw),
     &                                x(il,iw)-x(il-1,iw)))/DEG2RAD
            if(strike(nps).lt.0.d0)strike(nps)=strike(nps)+360.d0
c
c           determine two diagonal vectors
c
            dx1=x(il,iw-1)-x(il-1,iw)
            dy1=y(il,iw-1)-y(il-1,iw)
            dz1=z(il,iw-1)-z(il-1,iw)
            dx2=x(il,iw)-x(il-1,iw-1)
            dy2=y(il,iw)-y(il-1,iw-1)
            dz2=z(il,iw)-z(il-1,iw-1)
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
c
            dlen(nps)=0.5d0*(dsqrt((x(il,iw)-x(il-1,iw))**2
     &                            +(y(il,iw)-y(il-1,iw))**2)
     &                      +dsqrt((x(il,iw-1)-x(il-1,iw-1))**2
     &                            +(y(il,iw-1)-y(il-1,iw-1))**2))
            dwid(nps)=parea(nps)/dlen(nps)
c
c           convert to geographic coordinates
c
c           spherical triangle:
c           A = pole, B = source position, C = reference position
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
c
        write(*,'(a,i2,a,i6,a)')' the ',is,'. fault segment => ',
     &                         ips,' slip patches.'
        nlength(is)=nlength(is)-1
        nwidth(is)=nwidth(is)-1
      enddo
      write(*,*)'------------------------------------------------'
      write(*,'(a,i7)')' total number of point sources: ',nps
c
      nps1(1)=1
      nps2(1)=nlength(1)*nwidth(1)
      do is=2,ns
        nps1(is)=nps2(is-1)+1
        nps2(is)=nps1(is)+nlength(is)*nwidth(is)-1
      enddo
c
      sarea=0.d0
      do ips=1,nps
        sarea=sarea+parea(ips)
      enddo
c
      do ips=1,nps
        st=strike(ips)*DEG2RAD
        di=dip(ips)*DEG2RAD
        d0=dmax1(d0,dsqrt(dlen(ips)**2+dwid(ips)**2
     &                  +(dwid(ips)*dsin(di))**2))
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
c
      enddo
      return
      end
