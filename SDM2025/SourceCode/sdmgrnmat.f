      subroutine sdmgrnmat(ierr)
      use sdmalloc
      implicit none
      integer*4 ierr
c
c     Green's functions matrix element Mnm
c
c     Last modified: Zhuhai, Dec 7, 2025, by R. Wang
c
      integer*4 i,j,m,n,igd,ira,ips,iobs,jps,ipar
      real*8 sd,paru,wfsum,dal,daw
      real*8 dcg(6),flr(-1:1),fud(-1:1)
c
c     observation matrix
c     m(1:nobs): raw, n(1:nps*2+npar): column
c
      paru=0.d0
      wfsum=0.d0
      do iobs=1,nobs
        m=iobs
        n=0
        do ips=1,nps
          do ira=1,2
            n=n+1
            obsgrnmat(m,n)=wf(iobs)*datgrn(ira,ips,iobs)/zhy(ips)
            paru=paru+obsgrnmat(m,n)**2
            wfsum=wfsum+wf(iobs)**2
          enddo
        enddo
      enddo
c
c     offunit: a very important parameter to avoid loss-of-precision problem!
c
      paru=dsqrt(paru/wfsum)
      do ipar=1,npar
        parunit(ipar)=paru/parunit(ipar)
      enddo
c
      do iobs=1,nobs
        m=iobs
        do ipar=1,npar
          n=nps*2+ipar
          obsgrnmat(m,n)=wf(iobs)*corrgrn(ipar,iobs)*parunit(ipar)
        enddo
      enddo
c
c     smoothing matrix
c     m(1:nps*nsmocmp): raw, n(1:nps*2+npar): column
c
      do m=1,nps*nsmocmp
        do n=1,nsys
          smogrnmat(m,n)=0.d0
        enddo
      enddo
      if(ismooth.eq.1)then
        do ips=1,nps
          do jps=1,nps
            dal=parea(jps)/dlen(jps)**2
            daw=parea(jps)/dwid(jps)**2
            do ira=1,2
              do i=1,4
                dcg(i)=0.d0
              enddo
              do i=-1,1
                flr(i)=0.d0
                fud(i)=0.d0
              enddo
              if(ips.eq.jps)then
                if(ipsl(jps).le.0.or.ipsr(jps).le.0)then
                  flr(0)=1.d0
                else
                  flr(0)=2.d0
                endif
                if(ipsu(jps).le.0.or.ipsd(jps).le.0)then
                  fud(0)=1.d0
                else
                  fud(0)=2.d0
                endif
              else if(ips.eq.ipsl(jps))then
                flr(-1)=1.d0
              else if(ips.eq.ipsr(jps))then
                flr(1)=1.d0
              else if(ips.eq.ipsu(jps))then
                fud(-1)=1.d0
              else if(ips.eq.ipsd(jps))then
                fud(1)=1.d0
              endif
              dcg(2*ira-1)=(flr(-1)-flr(0)+flr(1))*dal
              dcg(2*ira  )=(fud(-1)-fud(0)+fud(1))*daw
c
              do i=1,nsmocmp
                m=(jps-1)*nsmocmp+i
                n=(ips-1)*2+ira
                smogrnmat(m,n)=dsqrt(parea(jps))
     &                        *dcg(i)*zhy(jps)/zhy(ips)
              enddo
            enddo
          enddo
        enddo
      else
        do ips=1,nps
          do jps=1,nps
            dal=parea(jps)/dlen(jps)**2
            daw=parea(jps)/dwid(jps)**2
            do ira=1,2
              do i=1,6
                dcg(i)=0.d0
              enddo
              if(ipsl(jps).gt.0)then
                do i=1,3
                  dcg(i)=dcg(i)
     &              +(strgrn(ira,ips,i,jps)-strgrn(ira,ips,i,ipsl(jps)))
     &              *dal
                enddo
              endif
              if(ipsr(jps).gt.0)then
                do i=1,3
                  dcg(i)=dcg(i)
     &              +(strgrn(ira,ips,i,jps)-strgrn(ira,ips,i,ipsr(jps)))
     &              *dal
                enddo
              endif
              if(ipsu(jps).gt.0)then
                do i=1,3
                  dcg(i+3)=dcg(i+3)
     &              +(strgrn(ira,ips,i,jps)-strgrn(ira,ips,i,ipsu(jps)))
     &              *daw
                enddo
              endif
              if(ipsd(jps).gt.0)then
                do i=1,3
                  dcg(i+3)=dcg(i+3)
     &              +(strgrn(ira,ips,i,jps)-strgrn(ira,ips,i,ipsd(jps)))
     &              *daw
                enddo
              endif
c
              do i=1,nsmocmp
                m=(jps-1)*nsmocmp+i
                n=(ips-1)*2+ira
                smogrnmat(m,n)=dsqrt(parea(jps))
     &                        *dcg(i)*zhy(jps)/zhy(ips)
              enddo
            enddo
          enddo
        enddo
      endif
c
      if(izhy.eq.0)return

      do n=1,nps*2
        sd=0.d0
        do m=1,nps*nsmocmp
          sd=sd+smogrnmat(m,n)**2
        enddo
        sd=dsqrt(sd)
        do m=1,nps*nsmocmp
          smogrnmat(m,n)=smogrnmat(m,n)/sd
        enddo
      enddo
c
      return
      end