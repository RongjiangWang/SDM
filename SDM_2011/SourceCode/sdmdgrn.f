      subroutine sdmdgrn(ierr)
      implicit none
c
c     Last modified: Potsdam, Oct, 2008, by R. Wang
c
      integer ierr
c
      include 'sdmglob.h'
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     GREEN'S FUNNCTIONN PARAMETERS
c     =============================
c
c     Green's function source types:
c       1 = strike-slip (m12=m21=1)
c       2 = dip-slip (m13=m31=1)
c       3 = compensated linear vector dipole (CLVD)
c           (m11=m22=-1/2, m33=1) (no tangential component)
c     Green's function coordinate system:
c       (z,r,t) = cylindrical with z being downward(!)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer unit(3,3)
      character*163 greens(3,3)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     LOCAL WORK SPACES
c     =================
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer i,n,ir,izs
      integer lend,lenf,istp
      double precision r1,r2,dr,dract,rratio
      double precision zs0,zs1,zs2,dzs,dzact,zratio
      double precision las,mus,rhos,rhoobs
      double precision d1,d2,d3,d4
      character*180 dataline
      logical clsh
c
c     OPEN GREEN'S FUNCTIONS FILES
c     ============================
c
      do lend=len(grndir),1,-1
        if(grndir(lend:lend).ne.' ')goto 100
      enddo
100   continue
      if(grndir(lend:lend).ne.'\\'.or.grndir(lend:lend).ne.'/')then
        grndir=grndir(1:lend)//'/'
        lend=lend+1
      endif
      do i=1,3
        do lenf=len(green(i)),1,-1
          if(green(i)(lenf:lenf).ne.' ')goto 110
        enddo
110     continue
        greens(i,1)=grndir(1:lend)//green(i)(1:lenf)//'.ss'
        greens(i,2)=grndir(1:lend)//green(i)(1:lenf)//'.ds'
        greens(i,3)=grndir(1:lend)//green(i)(1:lenf)//'.cl'
      enddo
c
      do istp=1,3
        do i=1,3
          clsh=istp.eq.3.and.i.eq.3
          if(clsh)goto 300
          unit(i,istp)=10+3*(istp-1)+i
          open(unit(i,istp),file=greens(i,istp),status='old')
          if(i*istp.eq.1)then
            call getdata(unit(i,istp),dataline)
            read(dataline,*)nr,r1,r2,rratio
            if(r1.lt.0.d0.or.r2.lt.r1.or.nr.le.1)then
              stop 'Error: wrong no of distance samples!'
            else if(nr.gt.NRMAX)then
              stop 'Error: max. no of distance samples exceeded!'
            else if(nr.eq.2)then
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
            call getdata(unit(i,istp),dataline)
            read(dataline,*)zobs,laobs,muobs,rhoobs
            call getdata(unit(i,istp),dataline)
            read(dataline,*)nzs,zs1,zs2,zratio
            if(zs1.le.0.d0.or.zs2.lt.zs1.or.nzs.le.1)then
              stop 'Error: wrong no of depth samples!'
            else if(nzs.gt.NZSMAX)then
              stop 'Error: max. no of depth samples exceeded!'
            else if(nzs.eq.2)then
              dzs=zs2-zs1
              zs(1)=zs1
              zs(2)=zs2
            else
              dzs=2.d0*(zs2-zs1)/dble(nzs-1)/(1.d0+zratio)
              zs(1)=zs1
              do izs=2,nzs
                dzact=dzs*(1.d0+(zratio-1.d0)*dble(izs-2)/dble(nzs-2))
                zs(izs)=zs(izs-1)+dzact
              enddo
            endif
          else
            call getdata(unit(i,istp),dataline)
            read(dataline,*)n,d1,d2,d3
            if(n.ne.nr.or.d1.ne.r1.or.d2.ne.r2.or.d3.ne.rratio)then
              stop 'Error: different observation sampling in Greens!'
            endif
            call getdata(unit(i,istp),dataline)
            read(dataline,*)d1,d2,d3,d4
            if(d1.ne.zobs.or.d2.ne.laobs.or.
     &         d3.ne.muobs.or.d4.ne.rhoobs)then
              stop 'Error: diff. observation site parameters in Greens!'
            endif
            call getdata(unit(i,istp),dataline)
            read(dataline,*)n,d1,d2,d3
            if(n.ne.nzs.or.d1.ne.zs1.or.d2.ne.zs2.or.d3.ne.zratio)then
              stop 'Error: different source sampling in Greens!'
            endif
          endif
300       continue
        enddo
      enddo
c
c     all Green's function files have been opened
c     ===========================================
c
      do izs=1,nzs
c
c       read in Green's functions
c
        do istp=1,3
          do i=1,3
            clsh=istp.eq.3.and.i.eq.3
            if(clsh)goto 400
            if(i.eq.1)then
              call getdata(unit(i,istp),dataline)
              read(dataline,*)zs0,las,mus,rhos
              muz(izs)=mus
            else
              call getdata(unit(i,istp),dataline)
              read(dataline,*)d1,d2,d3,d4
              if(d1.ne.zs0.or.d2.ne.las.or.
     &           d3.ne.mus.or.d4.ne.rhos)then
                stop 'Error: different s. layer parameters in greens!'
              endif
            endif
            read(unit(i,istp),*)(grns(izs,ir,i,istp),ir=1,nr)
400         continue
          enddo
        enddo
      enddo
c
      do istp=1,3
        do i=1,3
          clsh=istp.eq.3.and.i.eq.3
          if(.not.clsh)close(unit(i,istp))
        enddo
      enddo
c
      ierr=0
c
      return
      end
