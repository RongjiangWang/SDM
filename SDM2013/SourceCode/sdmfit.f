      subroutine sdmfit(ns,xcs0,ycs0,zcs0,latobs0,lonobs0,disp0)
      implicit none
c
c     Last modified: Potsdam, Oct, 2008, by R. Wang
c
      integer ns
      double precision xcs0,ycs0,zcs0,latobs0,lonobs0
      double precision disp0(4)
c
      include 'sdmglob.h'
c
      integer i,is,ips,ira,il,nl
      integer ir,ir1,ir2,izs,izs1,izs2
      double precision dr,dl,st,di,dux,duy,duz,dur,dut
      double precision dz1,dz2,pz1,pz2,dpz,px,py,xp,yp
      double precision xobs0,yobs0
      double precision csst,ssst,csra,ssra,csdi,ssdi
      double precision si,co,si2,co2,dis,azi,area,slp2
      double precision w1,w2,dwei,psss,psds,pscl,shss,shds
c
c     from Okada's subroutine DC3D0:
c
      INTEGER IRET
      REAL*4 ALPHA,X,Y,Z,DEPTH,DIPS,DX,DY,DZ,
     &       UX,UY,UZ,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ
c
c     more from Okada's subroutine DC3D:
c
      REAL*4 AL1,AL2,AW1,AW2,DISL1,DISL2,DISL3
c
c     SUPERPOSITION OF ALL DISCRETE POINT SOURCES
c     ===========================================
      do i=1,4
        disp0(i)=0.d0
      enddo
      if(hsmodel)then
        ALPHA=sngl(0.5d0/(1.d0-poisson))
      else
        ALPHA=sngl((laobs+muobs)/(laobs+2.d0*muobs))
      endif
      Z=-sngl(zobs)
      DISL3=0.0
c
      do is=1,ns
        do ips=nps1(is),nps2(is)
          slp2=slpmdl(ips,1)**2+slpmdl(ips,2)**2
          if(slp2.gt.0.d0)then
            st=strike(ips)*DEG2RAD
            csst=dcos(st)
            ssst=dsin(st)
c
	      di=dip(ips)*DEG2RAD
            csdi=dcos(di)
            ssdi=dsin(di)
c
            DEPTH=sngl(pz(ips))
            DIPS=sngl(dip(ips))
c
            DISL1=slpmdl(ips,1)
            DISL2=slpmdl(ips,2)
c
            AL1=sngl(-0.5d0*dlen(ips))
            AL2=sngl(+0.5d0*dlen(ips))
            AW1=sngl(-0.5d0*dwid(ips))
            AW2=sngl(+0.5d0*dwid(ips))
c
c           transform from Aki's to Okada's system
c
            call disazi(REARTH,plat(ips),plon(ips),latobs0,lonobs0,
     &                  xobs0,yobs0)
            X=sngl(xobs0*csst+yobs0*ssst)
            Y=sngl(xobs0*ssst-yobs0*csst)
            IRET=1
            call DC3D(ALPHA,X,Y,Z,DEPTH,DIPS,AL1,AL2,AW1,AW2,
     &                DISL1,DISL2,DISL3,UX,UY,UZ,
     &                UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ,IRET)
c
c           transform from Okada's to Aki's system
c
            dux=dble(UX)*csst+dble(UY)*ssst
            duy=dble(UX)*ssst-dble(UY)*csst
            duz=-dble(UZ)
            disp0(1)=disp0(1)+dux
            disp0(2)=disp0(2)+duy
            disp0(3)=disp0(3)+duz
            disp0(4)=disp0(4)+dux*xcs0+duy*ycs0+duz*zcs0
          endif
        enddo
      enddo
c
      if(hsmodel)return
      dr=r(2)-r(1)
      do is=1,ns
        do ips=nps1(is),nps2(is)
          slp2=slpmdl(ips,1)**2+slpmdl(ips,2)**2
          if(slp2.gt.0.d0)then
            st=strike(ips)*DEG2RAD
            csst=dcos(st)
            ssst=dsin(st)
c
	      di=dip(ips)*DEG2RAD
            csdi=dcos(di)
            ssdi=dsin(di)
c
            area=dlen(ips)*dwid(ips)
            pz1=pz(ips)-0.5d0*dwid(ips)*ssdi
            pz2=pz(ips)+0.5d0*dwid(ips)*ssdi
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
            nl=1+idint(dlen(ips)/dr)
            dl=dlen(ips)/dble(nl)
c
            do izs=izs1,izs2
              if(izs.eq.izs1)then
                if(izs.lt.nzs)then
                  dpz=dmax1(0.d0,dmin1(pz2,
     &                      0.5d0*(zs(izs)+zs(izs+1)))-pz1)
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
              yp=(zs(izs)-pz(ips))/ssdi
              dwei=(dl*dpz/ssdi)/area
              do il=1,nl
                xp=-0.5d0*dwid(ips)+(dble(il)-0.5d0)*dl
                px=xp*csst-ssst*yp*csdi
                py=xp*ssst+csst*yp*csdi
c
c               transform from Aki's to Okada's system
c
                call disazi(REARTH,plat(ips),plon(ips),latobs0,lonobs0,
     &                  xobs0,yobs0)
                xobs0=xobs0-px
                yobs0=yobs0-py
c
                dis=dsqrt(xobs0**2+yobs0**2)
                if(dis.gt.0.d0)then
                  azi=datan2(yobs0,xobs0)
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
c                 add differential Green functions
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                 pmwei(1-5):
c                 1 = weight for strike-slip: m12=m21=1;
c                     poloidal *sin(2*theta), toroidal *cos(2*theta)
c
c                 2 = weight for dip-slip: m13=m31=1
c                     poloidal *cos(theta), toroidal *sin(theta)
c
c                 3 = weight for clvd: m33=-m11=-m22=1
c                     axisymmetric
c
c                 4 = weight for 45 deg strike-slip: m11=-m22=1
c                     greenfct4(theta) = green1(theta + 45 deg)
c
c                 5 = weight for 45 deg dip-slip: m23=m32=1
c                     greenfct5(theta) = green2(theta - 90 deg)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                  psss=pmwei(1,ips,ira)*si2+pmwei(4,ips,ira)*co2
                  shss=pmwei(1,ips,ira)*co2-pmwei(4,ips,ira)*si2
                  psds=pmwei(2,ips,ira)*co+pmwei(5,ips,ira)*si
                  shds=pmwei(2,ips,ira)*si-pmwei(5,ips,ira)*co
                  pscl=pmwei(3,ips,ira)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                 contributions from the strike-slip components
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                  duz=psss
     &               *(w1*grns(izs,ir1,1,1)+w2*grns(izs,ir2,1,1))
                  dur=psss
     &               *(w1*grns(izs,ir1,2,1)+w2*grns(izs,ir2,2,1))
                  dut=shss
     &               *(w1*grns(izs,ir1,3,1)+w2*grns(izs,ir2,3,1))
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                 contributions from the dip-slip components
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                  duz=duz+psds
     &               *(w1*grns(izs,ir1,1,2)+w2*grns(izs,ir2,1,2))
                  dur=dur+psds
     &               *(w1*grns(izs,ir1,2,2)+w2*grns(izs,ir2,2,2))
                  dut=dut+shds
     &               *(w1*grns(izs,ir1,3,2)+w2*grns(izs,ir2,3,2))
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                 contributions from the clvd components
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                  duz=duz+pscl
     &               *(w1*grns(izs,ir1,1,3)+w2*grns(izs,ir2,1,3))
                  dur=dur+pscl
     &               *(w1*grns(izs,ir1,2,3)+w2*grns(izs,ir2,2,3))
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                  dux=(dur*co-dut*si)*slpmdl(ips,ira)
                  duy=(dur*si+dut*co)*slpmdl(ips,ira)
                  duz=duz*slpmdl(ips,ira)
                  disp0(1)=disp0(1)+dux
                  disp0(2)=disp0(2)+duy
                  disp0(3)=disp0(3)+duz
                  disp0(4)=disp0(4)+dux*xcs0+duy*ycs0+duz*zcs0
                enddo
              enddo
            enddo
          endif
        enddo
      enddo
c
      return
      end
