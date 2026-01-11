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
      real*8 sd,paru,wfsum
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
            obsmat(m,n)=wf(iobs)*datgrn(ira,ips,iobs)/zhy(ips)
            paru=paru+obsmat(m,n)**2
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
          obsmat(m,n)=wf(iobs)*corrgrn(ipar,iobs)*parunit(ipar)
        enddo
      enddo
c
c     smoothing matrix
c     m(1:nps*nsmocmp): raw, n(1:nps*2+npar): column
c
      m=0
      do jps=1,nps
        do i=1,nsmocmp
          m=m+1
          n=0
          do ips=1,nps
            do ira=1,2
              n=n+1
              smomat(m,n)=dsqrt(parea(jps))*dcgrn(ira,ips,i,jps)
     &                   *zhy(jps)/zhy(ips)
            enddo
          enddo
          do ipar=1,npar
            n=n+1
            smomat(m,n)=0.d0
          enddo
        enddo
      enddo
c
      if(izhy.eq.0)return

      do n=1,nps*2
        sd=0.d0
        do m=1,nps*nsmocmp
          sd=sd+smomat(m,n)**2
        enddo
        sd=dsqrt(sd)
        do m=1,nps*nsmocmp
          smomat(m,n)=smomat(m,n)/sd
        enddo
      enddo
c
      return
      end