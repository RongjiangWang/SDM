      subroutine sdmdisc(ierr)
      use sdmalloc
      implicit none
      integer*4 ierr
c
c     Last modified: 2025 in Zhuhai by R. Wang
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     LOCAL WORK SPACES
c     =================
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer*4 i,is,iw,il,ips,jps,ira,nlmax,nwmax
      real*8 lat0,lon0,st,st0,di,dl,dw,pn,pe
      real*8 dx,dy,xp,yp,zp,d0,dp0,dp
      real*8 dx1,dx2,dy1,dy2,dz1,dz2,dnx,dny,dnz,hdis,vdis
      real*8 x0,y0,bga,bgc,sma,smb,smc,alf,beta,d1,d2,dd
      real*8 ra(2),sm(3,3)
      character*10 header
c
      real*8, allocatable:: lft(:,:),xft(:,:),yft(:,:)
      real*8, allocatable:: x(:,:),y(:,:),z(:,:),dip1(:),dip2(:)
c
      real*8 EPSDIP
      data EPSDIP/0.1d0/
c
      allocate(nps1(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmdisc: nps1 not allocated!'
      allocate(nps2(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmdisc: nps2 not allocated!'
c
      ra(1)=0.d0
      ra(2)=90.d0*DEG2RAD
c
      if(idisc.eq.0)then
c
c       test read
c
        nps=0
        do is=1,ns
          nps1(is)=0
          nps2(is)=0
          open(20,file=inpatches(is),status='old')
          read(20,'(a)')header
10        continue
          read(20,*,end=20)alf
          nps=nps+1
          if(nps.gt.0)goto 10
20        close(20)
          if(is.eq.1)then
            nps1(is)=1
            nps2(is)=nps
          else
            nps1(is)=nps2(is-1)+1
            nps2(is)=nps
          endif
        enddo
c
c       final read
c
        call sdmpatches(ierr)
c
        do is=1,ns
          open(20,file=inpatches(is),status='old')
          read(20,'(a)')header
          do ips=nps1(is),nps2(is)
          read(20,*)plat(ips),plon(ips),pz(ips),
     &              dlen(ips),dwid(ips),strike(ips),dip(ips)
          pz(ips)=pz(ips)*KM2M
          dlen(ips)=dlen(ips)*KM2M
          dwid(ips)=dwid(ips)*KM2M
          parea(ips)=dlen(ips)*dwid(ips)
c
c         added on May 18, 2024
c
          st=strike(ips)*DEG2RAD
          di=dip(ips)*DEG2RAD
c
          if(isref(is).ge.1.and.isref(is).le.4)then
c
c             deleted on May 18, 2024:
c
c             st=strike(ips)*DEG2RAD
c             di=dip(ips)*DEG2RAD
c
              if(isref(is).eq.1)then
                pz(ips)=pz(ips)+0.5d0*dwid(ips)*dsin(di)
                pn=0.5d0*dlen(ips)*dcos(st)
     &            -0.5d0*dwid(ips)*dcos(di)*dsin(st)
                pe=0.5d0*dlen(ips)*dsin(st)
     &            +0.5d0*dwid(ips)*dcos(di)*dcos(st)
              else if(isref(is).eq.2)then
                pz(ips)=pz(ips)+0.5d0*dwid(ips)*dsin(di)
                pn=-0.5d0*dlen(ips)*dcos(st)
     &             -0.5d0*dwid(ips)*dcos(di)*dsin(st)
                pe=-0.5d0*dlen(ips)*dsin(st)
     &             +0.5d0*dwid(ips)*dcos(di)*dcos(st)
              else if(isref(is).eq.3)then
                pz(ips)=pz(ips)-0.5d0*dwid(ips)*dsin(di)
                pn=0.5d0*dlen(ips)*dcos(st)
     &            +0.5d0*dwid(ips)*dcos(di)*dsin(st)
                pe=0.5d0*dlen(ips)*dsin(st)
     &            -0.5d0*dwid(ips)*dcos(di)*dcos(st)
              else
                pz(ips)=pz(ips)-0.5d0*dwid(ips)*dsin(di)
                pn=-0.5d0*dlen(ips)*dcos(st)
     &             +0.5d0*dwid(ips)*dcos(di)*dsin(st)
                pe=-0.5d0*dlen(ips)*dsin(st)
     &             -0.5d0*dwid(ips)*dcos(di)*dcos(st)
              endif
c
c             determine central point of the patch
c
c             spherical triangle:
c             A = pole, B = source position, C = reference position
c
              sma=dsqrt(pn**2+pe**2)/REARTH
              smb=0.5d0*PI-plat(ips)*DEG2RAD
              bgc=datan2(pe,pn)
              smc=dacos(dcos(sma)*dcos(smb)
     &           +dsin(sma)*dsin(smb)*dcos(bgc))
              bga=dasin(dsin(sma)*dsin(bgc)/dsin(smc))
c
c             geographic coordinate of the equivalent point source
c
              plat(ips)=90.d0-smc/DEG2RAD
              plon(ips)=dmod(plon(ips)+bga/DEG2RAD,360.d0)
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
              pmwei(1,ips,ira)=sm(1,2)*parea(ips)
              pmwei(2,ips,ira)=sm(1,3)*parea(ips)
              pmwei(3,ips,ira)=sm(3,3)*parea(ips)
              pmwei(4,ips,ira)=0.5d0*(sm(1,1)-sm(2,2))*parea(ips)
              pmwei(5,ips,ira)=sm(2,3)*parea(ips)
            enddo
          enddo
          write(*,'(a,i2,a,i6,a)')' the ',is,'. fault segment => ',
     &                         1+nps2(is)-nps1(is),' slip patches.'
          close(20)
        enddo
        goto 500
      endif
c
      allocate(nlength(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmdisc: nlength not allocated!'
      allocate(nwidth(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmdisc: nwidth not allocated!'
      allocate(lft(nftmax,ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmdisc: lft not allocated!'
      allocate(xft(nftmax,ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmdisc: xft not allocated!'
      allocate(yft(nftmax,ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmdisc: yft not allocated!'
c
      nps=0
      nlmax=0
      nwmax=0
      do is=1,ns
        lat0=0.5d0*(latft(1,is)+latft(nft(is),is))
        lon0=0.5d0*(lonft(1,is)+lonft(nft(is),is))
c
c       local cartesian coordinates of surface ruptures
c
        do i=1,nft(is)
          call disazi(REARTH,lat0,lon0,latft(i,is),lonft(i,is),
     &                xft(i,is),yft(i,is))
        enddo
c
c       accumulated length of top fault trace
c
        lft(1,is)=0.d0
        do i=2,nft(is)
          lft(i,is)=lft(i-1,is)+dsqrt((xft(i,is)-xft(i-1,is))**2
     &                         +(yft(i,is)-yft(i-1,is))**2)
        enddo
        write(*,'(a,i4,a,f6.2,a)')' Length along the curved strike ',
     &       is,'. fault segment: ',
     &       lft(nft(is),is)/KM2M,' km'
c
c       length sampling parameters
c
        nlength(is)=max0(2,1+idnint(lft(nft(is),is)/patchsize(is)))
        nwidth(is)=max0(2,1+idnint(width(is)/patchsize(is)))
        nps=nps+(nlength(is)-1)*(nwidth(is)-1)
        nlmax=max0(nlmax,nlength(is))
        nwmax=max0(nwmax,nwidth(is))
      enddo
c
      allocate(dip1(nlmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmdisc: dip1 not allocated!'
      allocate(dip2(nlmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmdisc: dip2 not allocated!'
      allocate(x(nlmax,nwmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmdisc: x not allocated!'
      allocate(y(nlmax,nwmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmdisc: y not allocated!'
      allocate(z(nlmax,nwmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmdisc: z not allocated!'
c
      call sdmpatches(ierr)
c
c     discretise whole fault plane
c
      ips=0
      do is=1,ns
        dl=lft(nft(is),is)/dble(nlength(is)-1)
c
c       discretise top fault trace
c       x(il,iw), y(il,iw) = cartesian coordinates
c       il = index increasing with length (strike)
c       iw = index increasing with width (down dip)
c
        x(1,1)=xft(1,is)
        y(1,1)=yft(1,is)
        z(1,1)=topdep(is)
        dip1(1)=topdip(1,is)
        dip2(1)=botdip(1,is)
        x(nlength(is),1)=xft(nft(is),is)
        y(nlength(is),1)=yft(nft(is),is)
        z(nlength(is),1)=topdep(is)
        dip1(nlength(is))=topdip(nft(is),is)
        dip2(nlength(is))=botdip(nft(is),is)
c
        do il=2,nlength(is)-1
          z(il,1)=topdep(is)
          do i=2,nft(is)
            if(lft(i,is).ge.dble(il-1)*dl)then
              beta=(dble(il-1)*dl-lft(i-1,is))/(lft(i,is)-lft(i-1,is))
              x(il,1)=xft(i-1,is)+(xft(i,is)-xft(i-1,is))*beta
              y(il,1)=yft(i-1,is)+(yft(i,is)-yft(i-1,is))*beta
              dip1(il)=topdip(i-1,is)+(topdip(i,is)-topdip(i-1,is))*beta
              dip2(il)=botdip(i-1,is)+(botdip(i,is)-botdip(i-1,is))*beta
              goto 400
            endif
          enddo
400       continue
        enddo
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
        do il=2,nlength(is)
          do iw=2,nwidth(is)
            ips=ips+1
            pl(ips)=(dble(il)-1.5d0)*dl
            pw(ips)=(dble(iw)-1.5d0)*dw
            pn=0.25d0*(x(il,iw)+x(il-1,iw)+x(il,iw-1)+x(il-1,iw-1))
            pe=0.25d0*(y(il,iw)+y(il-1,iw)+y(il,iw-1)+y(il-1,iw-1))
            pz(ips)=0.25d0*(z(il,iw)+z(il-1,iw)+z(il,iw-1)+z(il-1,iw-1))
c
            strike(ips)=0.5d0*(datan2(y(il,iw-1)-y(il-1,iw-1),
     &                                x(il,iw-1)-x(il-1,iw-1))
     &                        +datan2(y(il,iw)-y(il-1,iw),
     &                                x(il,iw)-x(il-1,iw)))/DEG2RAD
            if(strike(ips).lt.0.d0)strike(ips)=strike(ips)+360.d0
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
            parea(ips)=0.5d0*dsqrt(dnx**2+dny**2+dnz**2)
c
            dip(ips)=dacos(0.5d0*dnz/parea(ips))/DEG2RAD
c
            dlen(ips)=0.5d0*(dsqrt((x(il,iw)-x(il-1,iw))**2
     &                            +(y(il,iw)-y(il-1,iw))**2)
     &                      +dsqrt((x(il,iw-1)-x(il-1,iw-1))**2
     &                            +(y(il,iw-1)-y(il-1,iw-1))**2))
            dwid(ips)=parea(ips)/dlen(ips)
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
            plat(ips)=90.d0-smc/DEG2RAD
            plon(ips)=dmod(lon0+bga/DEG2RAD,360.d0)
c
            st=strike(ips)*DEG2RAD
            di=dip(ips)*DEG2RAD
            if(pz(ips)-0.5d0*dwid(ips)*dsin(di).le.0.d0)then
c
c             upper patch edge exceeds the surface
c
              pz(ips)=0.5001d0*dwid(ips)*dsin(di)
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
              pmwei(1,ips,ira)=sm(1,2)*parea(ips)
              pmwei(2,ips,ira)=sm(1,3)*parea(ips)
              pmwei(3,ips,ira)=sm(3,3)*parea(ips)
              pmwei(4,ips,ira)=0.5d0*(sm(1,1)-sm(2,2))*parea(ips)
              pmwei(5,ips,ira)=sm(2,3)*parea(ips)
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
500   continue
c
      do is=1,ns
        do ips=nps1(is),nps2(is)
          st=strike(ips)*DEG2RAD
          di=dip(ips)*DEG2RAD
          d0=dmax1(d0,dsqrt(dlen(ips)**2+dwid(ips)**2
     &                    +(dwid(ips)*dsin(di))**2))
c
c         search left neighboring patch
c
          dp=d0
          ipsl(ips)=0
          xp=-dlen(ips)*dcos(st)
          yp=-dlen(ips)*dsin(st)
          zp=pz(ips)
c
          do jps=nps1(is),nps2(is)
            if(jps.ne.ips.and.dabs(pz(jps)-pz(ips)).le.
     &         dwid(ips)*dsin(di)+dwid(jps)*dsin(dip(jps)*DEG2RAD))then
              call disazi(rearth,plat(ips),plon(ips),
     &                           plat(jps),plon(jps),dx,dy)
              dp0=dsqrt((dx-xp)**2+(dy-yp)**2+(pz(jps)-zp)**2)
              if(dp.ge.dp0)then
                dp=dp0
                ipsl(ips)=jps
              endif
            endif
          enddo
c
c         search right neighboring patch
c
          dp=d0
          ipsr(ips)=0
          xp=dlen(ips)*dcos(st)
          yp=dlen(ips)*dsin(st)
          zp=pz(ips)
c
          do jps=nps1(is),nps2(is)
            if(jps.ne.ips.and.dabs(pz(jps)-pz(ips)).le.
     &         dwid(ips)*dsin(di)+dwid(jps)*dsin(dip(jps)*DEG2RAD))then
              call disazi(rearth,plat(ips),plon(ips),
     &                           plat(jps),plon(jps),dx,dy)
              dp0=dsqrt((dx-xp)**2+(dy-yp)**2+(pz(jps)-zp)**2)
              if(dp.ge.dp0)then
                dp=dp0
                ipsr(ips)=jps
              endif
            endif
          enddo
c
c         search upper neighboring patch
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
            do jps=nps1(is),nps2(is)
              if(jps.ne.ips.and.pz(jps).lt.pz(ips))then
                call disazi(rearth,plat(ips),plon(ips),
     &                             plat(jps),plon(jps),dx,dy)
                dp0=dsqrt((dx-xp)**2+(dy-yp)**2+(pz(jps)-zp)**2)
                if(dp.ge.dp0)then
                  dp=dp0
                  ipsu(ips)=jps
                endif
              endif
            enddo
          endif
c
c         search lower neighboring patch
c
          dp=d0
          ipsd(ips)=0
          xp=-dwid(ips)*dcos(di)*dsin(st)
          yp= dwid(ips)*dcos(di)*dcos(st)
          zp=pz(ips)+dwid(ips)*dsin(di)
c
          do jps=nps1(is),nps2(is)
            if(jps.ne.ips.and.pz(jps).gt.pz(ips))then
              call disazi(rearth,plat(ips),plon(ips),
     &                           plat(jps),plon(jps),dx,dy)
              dp0=dsqrt((dx-xp)**2+(dy-yp)**2+(pz(jps)-zp)**2)
              if(dp.ge.dp0)then
                dp=dp0
                ipsd(ips)=jps
              endif
            endif
          enddo
        enddo
      enddo
c
c     rupture area
c
      slparea=0.d0
      do ips=1,nps
        slparea=slparea+parea(ips)
      enddo
c
      if(idisc.eq.1)then
        deallocate(lft,xft,yft,x,y,z,dip1,dip2)
      endif
c
      return
      end
