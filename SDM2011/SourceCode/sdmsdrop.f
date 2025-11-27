      subroutine sdmsdrop(ns,nps)
      implicit none
c
c     calculate the slip gradient part of the cost function
c
c     Last modified: Potsdam, Oct, 2008, by R. Wang
c
      include 'sdmglob.h'
      integer*4 ns,nps
c
      integer*4 i,is,ips,jps
      real*8 sd,sdss,sdds,sdnn,weisum,smsum,smmax,sm90
      real*8 wei(NPSMAX),sm(NPSMAX)
c
      do ips=1,nps
        do i=1,3
          sd=0.d0
          do jps=1,nps
            sd=sd
     &        +slpmdl(jps,1)*strgrn(jps,1,ips,i)
     &        +slpmdl(jps,2)*strgrn(jps,2,ips,i)
          enddo
          strdrop(ips,i)=sd
        enddo
      enddo
c
      do is=1,ns
        smsum=0.d0
        do ips=nps1(is),nps2(is)
          wei(ips)=0.d0
          sm(ips)=parea(ips)*dsqrt(slpmdl(ips,1)**2+slpmdl(ips,2)**2)
          smsum=smsum+sm(ips)
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
          wei(jps)=1.d0
          weisum=weisum+1.d0
          sm90=sm90+smmax
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
          sdss=sdss+wei(ips)*strdrop(ips,1)
          sdds=sdds+wei(ips)*strdrop(ips,2)
          sdnn=sdnn+wei(ips)*strdrop(ips,3)
          maxsdrop(is)=dmax1(maxsdrop(is),
     &       strdrop(ips,1)**2+strdrop(ips,2)**2+strdrop(ips,3)**2)
        enddo
        maxsdrop(is)=dsqrt(maxsdrop(is))
        measdrop(is)=dsqrt(sdss**2+sdds**2)
c
        stdsdrop(is)=0.d0
        do ips=nps1(is),nps2(is)
          stdsdrop(is)=stdsdrop(is)+wei(ips)
     &       *((strdrop(ips,1)-sdss)**2
     &        +(strdrop(ips,2)-sdds)**2
     &        +(strdrop(ips,3)-sdnn)**2)
        enddo
        stdsdrop(is)=dsqrt(stdsdrop(is))
      enddo
c
      return
      end
