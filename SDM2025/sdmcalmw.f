      subroutine sdmcalmw(ierr)
      use sdmalloc
      implicit none
      integer*4 ierr
c
c     calculate moment magnitude
c
c     Last modified: 29 Nov 2025 by R. Wang
c
      integer*4 i,j,is,ira,ips
      real*8 st,di,ra,pm,sm
      real*8 smp(3,3),sms(3,3),sme(3,3)
c
      do i=1,3
        do j=1,3
          sme(i,j)=0.d0
        enddo
      enddo
c
c     Mwpsum = sum of moment release by each fault patch
c     Mwssum = sum of moment release by each fault segment
c     Mweq = moment of whole earthquake
c
      mwpsum=0.d0
      mwssum=0.d0
      mweq=0.d0
      do is=1,ns
        do i=1,3
          do j=1,3
            sms(i,j)=0.d0
          enddo
        enddo
        do ips=nps1(is),nps2(is)
          st=strike(ips)*DEG2RAD
          di=dip(ips)*DEG2RAD
          ra=datan2(slpmdl(2,ips),slpmdl(1,ips))
          pm=parea(ips)*dsqrt(slpmdl(1,ips)**2+slpmdl(2,ips)**2)*muehs
          mwpsum=mwpsum+pm
c
          smp(1,1)=-dsin(di)*dcos(ra)*dsin(2.d0*st)
     &            -dsin(2.d0*di)*dsin(ra)*(dsin(st))**2
          smp(1,2)=dsin(di)*dcos(ra)*dcos(2.d0*st)
     &            +0.5d0*dsin(2.d0*di)*dsin(ra)*dsin(2.d0*st)
          smp(1,3)=-dcos(di)*dcos(ra)*dcos(st)
     &           -dcos(2.d0*di)*dsin(ra)*dsin(st)
          smp(2,2)=dsin(di)*dcos(ra)*sin(2.d0*st)
     &            -dsin(2.d0*di)*dsin(ra)*(dcos(st))**2
          smp(2,3)=-dcos(di)*dcos(ra)*dsin(st)
     &            +dcos(2.d0*di)*dsin(ra)*dcos(st)
          smp(3,3)=dsin(2.d0*di)*dsin(ra)
          do i=1,3
            do j=i,3
              sme(i,j)=sme(i,j)+smp(i,j)*pm
              sms(i,j)=sms(i,j)+smp(i,j)*pm
            enddo
          enddo
        enddo
        sm=0.d0
        do i=1,3
          sm=sm+0.5d0*sms(i,i)**2
          do j=i+1,3
            sm=sm+sms(i,j)**2
          enddo
        enddo
        mwssum=mwssum+dsqrt(sm)
      enddo
      sm=0.d0
      do i=1,3
        sm=sm+0.5d0*sme(i,i)**2
        do j=i+1,3
          sm=sm+sme(i,j)**2
        enddo
      enddo
      mweq=dsqrt(sm)
c
      if(mwpsum.gt.0.d0)then
        mwpsum=(dlog10(mwpsum)-9.1d0)/1.5d0
      else
        mwpsum=0.d0
      endif
      if(mwssum.gt.0.d0)then
        mwssum=(dlog10(mwssum)-9.1d0)/1.5d0
      else
        mwssum=0.d0
      endif
      if(mweq.gt.0.d0)then
        mweq=(dlog10(mweq)-9.1d0)/1.5d0
      else
        mweq=0.d0
      endif
c
      return
      end