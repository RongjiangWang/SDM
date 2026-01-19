      subroutine sdmsdrop(ierr)
      use sdmalloc
      implicit none
      integer*4 ierr
c
c     calculate stress chnage on fault plane
c
c     Last modified: Zhuhai, Nov. 2025, by R. Wang
c
      integer*4 i,is,ips,jps
      real*8 sd,sdss,sdds,sdnn,weisum,smsum,smmax,sm90
c
      real*8, allocatable:: wei(:),sm(:)
c
      allocate(wei(nps),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmsdrop: wei not allocated!'
      allocate(sm(nps),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmsdrop: sd not allocated!'
c
      allocate(costress(3,nps),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmsdrop: costress not allocated!'
      allocate(measdrop(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmsdrop: measdrop not allocated!'
      allocate(stdsdrop(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmsdrop: stdsdrop not allocated!'
      allocate(maxsdrop(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmsdrop: maxsdrop not allocated!'
c
      do ips=1,nps
        do i=1,3
          sd=0.d0
          do jps=1,nps
            sd=sd
     &        +slpmdl(1,jps)*strgrn(1,jps,i,ips)
     &        +slpmdl(2,jps)*strgrn(2,jps,i,ips)
          enddo
          costress(i,ips)=sd
        enddo
      enddo
c
      do is=1,ns
        smsum=0.d0
        do ips=nps1(is),nps2(is)
          wei(ips)=0.d0
          sm(ips)=dsqrt(slpmdl(1,ips)**2+slpmdl(2,ips)**2)
          smsum=smsum+parea(ips)*sm(ips)
        enddo
c
        sm90=0.d0
        weisum=0.d0
10      smmax=0.d0
        jps=0
        do ips=nps1(is),nps2(is)
          if(wei(ips).eq.0.d0.and.sm(ips).gt.smmax)then
            smmax=sm(ips)
            jps=ips
          endif
        enddo
        if(jps.gt.0)then
          wei(jps)=parea(jps)
          weisum=weisum+wei(jps)
          sm90=sm90+parea(jps)*smmax
          if(sm90.lt.0.9d0*smsum)goto 10
        endif
c
        if(weisum.gt.0.d0)then
          do ips=nps1(is),nps2(is)
            wei(ips)=wei(ips)/weisum
          enddo
        endif
c
        sdss=0.d0
        sdds=0.d0
        sdnn=0.d0
        maxsdrop(is)=0.d0
        do ips=nps1(is),nps2(is)
          sdss=sdss+wei(ips)*costress(1,ips)
          sdds=sdds+wei(ips)*costress(2,ips)
          sdnn=sdnn+wei(ips)*costress(3,ips)
          maxsdrop(is)=dmax1(maxsdrop(is),
     &       costress(1,ips)**2+costress(2,ips)**2+costress(3,ips)**2)
        enddo
        maxsdrop(is)=dsqrt(maxsdrop(is))
        measdrop(is)=dsqrt(sdss**2+sdds**2)
c
        stdsdrop(is)=0.d0
        do ips=nps1(is),nps2(is)
          stdsdrop(is)=stdsdrop(is)+wei(ips)
     &       *((costress(1,ips)-sdss)**2
     &        +(costress(2,ips)-sdds)**2
     &        +(costress(3,ips)-sdnn)**2)
        enddo
        stdsdrop(is)=dsqrt(stdsdrop(is))
      enddo
c
      return
      end