      subroutine sdmpreinv(ierr)
      use sdmalloc
      implicit none
      integer*4 ierr
c
c     first step to prepare SDM iteration
c     Last modified: Zhuhai, Nov. 2025, by R. Wang
c
      integer*4 i,j,k,m,n,ira,ips,jps,igd,iobs,ipar
      real*8 a,b
      real*8 datvar,smovar
      real*8 maxsing
c
      real*8 eps
      data eps/1.0d-08/
c
      nsys=nps*2+npar
c
      allocate(slpmdl(2,nps),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmpreinv: slpmdl not allocated!'
c
      allocate(datmdl(nobs),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmpreinv: datmdl not allocated!'
      allocate(rmsres(ngd),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmpreinv: rmsres not allocated!'
      allocate(resmin(ngd),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmpreinv: resmin not allocated!'
      allocate(resmax(ngd),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmpreinv: resmax not allocated!'
c
c      if(niter.le.0)return
c
      allocate(sysmat(nsys,nsys),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmpreinv: sysmat not allocated!'
      allocate(matswp(nsys,nsys),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmpreinv: matswp not allocated!'
      allocate(sysbat(nsys),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmpreinv: sysbat not allocated!'
      allocate(batswp(nsys),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmpreinv: batswp not allocated!'
      allocate(sysvec(nsys),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmpreinv: sysvec not allocated!'
      allocate(vecswp(nsys),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmpreinv: vecswp not allocated!'
c
      allocate(obsmat(nobs,nsys),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmpreinv: obsmat not allocated!'
      allocate(smomat(nps*nsmocmp,nsys),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmpreinv: smomat not allocated!'
c
      allocate(obsswp(max0(nps*nsmocmp,nobs)),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmpreinv: obsswp not allocated!'
c
      allocate(resbat(nsys),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmpreinv: resbat not allocated!'
c
      call sdmgrnmat(ierr)
c
      wei2smo=0.d0
c
      do m=1,nsys
        do n=1,m
          sysmat(m,n)=0.d0
          do k=1,nobs
            sysmat(m,n)=sysmat(m,n)+obsmat(k,m)*obsmat(k,n)
          enddo
        enddo
      enddo
c
      do m=1,nsys
        do n=m+1,nsys
          sysmat(m,n)=sysmat(n,m)
        enddo
      enddo
c
      do m=1,nsys
        sysbat(m)=0.d0
        do k=1,nobs
          sysbat(m)=sysbat(m)+obsmat(k,m)*wf(k)*datobs(k)
        enddo
      enddo
c
      do m=1,nsys
        do n=1,m
          matswp(m,n)=0.d0
          do k=1,nps*nsmocmp
            matswp(m,n)=matswp(m,n)+smomat(k,m)*smomat(k,n)
          enddo
        enddo
      enddo
c
      do m=1,nsys
        do n=m+1,nsys
          matswp(m,n)=matswp(n,m)
        enddo
      enddo
c
      a=0.d0
      b=0.d0
      do i=1,nsys
        sysvec(i)=0.d0
        vecswp(i)=0.d0
        do j=1,nsys
          sysvec(i)=sysvec(i)+sysmat(i,j)*sysbat(j)
          vecswp(i)=vecswp(i)+matswp(i,j)*sysbat(j)
        enddo
        a=a+sysbat(i)*sysvec(i)
        b=b+sysbat(i)*vecswp(i)
      enddo
c
      if(b.le.0.d0)then
        stop ' Error in sdmpreinv: no slip model can be found!'
      endif
c
      wei2smo=wei2smo0*a/b
c
      do m=1,nsys
        do n=1,nsys
          sysmat(m,n)=sysmat(m,n)+wei2smo*matswp(m,n)
        enddo
      enddo
c
      sig2max=maxsing(sysmat,sysvec,vecswp,nsys,eps)
c
      do m=1,nsys
        do n=1,nsys
          sysmat(m,n)=sysmat(m,n)/sig2max
        enddo
        sysbat(m)=sysbat(m)/sig2max
      enddo
c
      do i=1,nsys
        vecswp(i)=0.d0
        do j=1,nsys
          vecswp(i)=vecswp(i)+sysmat(i,j)*sysbat(j)
        enddo
      enddo
      a=0.d0
      b=0.d0
      do i=1,nsys
        a=a+vecswp(i)*vecswp(i)
        b=b+vecswp(i)*sysbat(i)
      enddo
c
      step0=b/a
c
c     initialization
c
      do i=1,nsys
        sysvec(i)=0.d0
        vecswp(i)=0.d0
      enddo
c
      do ips=1,nps
        do ira=1,2
          slpmdl(ira,ips)=0.d0
        enddo
      enddo
      do ipar=1,npar
        corrpar(ipar)=0.d0
      enddo
c
      datnrm=0.d0
      do iobs=1,nobs
        datnrm=datnrm+(wf(iobs)*datobs(iobs))**2
      enddo
      datnrm=datnrm/sig2max
c
      deallocate(matswp)
c
      return
      end