      subroutine sdmpreinv(ierr)
      use sdmalloc
      implicit none
      integer*4 ierr
c
c     first step to prepare SDM iteration
c     Last modified: Zhuhai, Nov. 2025, by R. Wang
c
      integer*4 i,j,k,m,n,ira,ips,jps,igd,iobs,ipar
      real*8 a,b,sum
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
      allocate(smosysmat(nsys,nsys),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmpreinv: smosysmat not allocated!'
      allocate(obssysmat(nsys,nsys),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmpreinv: obssysmat not allocated!'
      allocate(sysbat(nsys),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmpreinv: sysbat not allocated!'
      allocate(batswp(nsys),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmpreinv: batswp not allocated!'
      allocate(sysvec(nsys),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmpreinv: sysvec not allocated!'
      allocate(vecswp(nsys),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmpreinv: vecswp not allocated!'
      allocate(vecini(nsys),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmpreinv: vecini not allocated!'
c
      allocate(obsgrnmat(nobs,nsys),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmpreinv: obsgrnmat not allocated!'
      allocate(smogrnmat(nps*nsmocmp,nsys),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmpreinv: smogrnmat not allocated!'
c
      allocate(resbat(nsys),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmpreinv: resbat not allocated!'
c
      call sdmgrnmat(ierr)
      do i=1,nsys
        vecini(i)=1.d0
      enddo
c
      wei2smo=0.d0
c
      do m=1,nsys
        do n=1,m
          obssysmat(m,n)=0.d0
          do k=1,nobs
            obssysmat(m,n)=obssysmat(m,n)+obsgrnmat(k,m)*obsgrnmat(k,n)
          enddo
        enddo
      enddo
c
      do m=1,nsys
        do n=m+1,nsys
          obssysmat(m,n)=obssysmat(n,m)
        enddo
      enddo
c
      do m=1,nsys
        sysbat(m)=0.d0
        do k=1,nobs
          sysbat(m)=sysbat(m)+obsgrnmat(k,m)*wf(k)*datobs(k)
        enddo
      enddo
c
      do m=1,nsys
        do n=1,m
          smosysmat(m,n)=0.d0
          do k=1,nps*nsmocmp
            smosysmat(m,n)=smosysmat(m,n)
     &                    +smogrnmat(k,m)*smogrnmat(k,n)
          enddo
        enddo
      enddo
c
      do m=1,nsys
        do n=m+1,nsys
          smosysmat(m,n)=smosysmat(n,m)
        enddo
      enddo
c
      a=0.d0
      b=0.d0
      do i=1,nsys
        sysvec(i)=0.d0
        vecswp(i)=0.d0
        do j=1,nsys
          sysvec(i)=sysvec(i)+obssysmat(i,j)*sysbat(j)
          vecswp(i)=vecswp(i)+smosysmat(i,j)*sysbat(j)
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
          sysmat(m,n)=obssysmat(m,n)+wei2smo*smosysmat(m,n)
        enddo
      enddo
c
      sig2max=maxsing(sysmat,nsys,eps,vecini,ierr)
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
      step=b/a
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
      return
      end