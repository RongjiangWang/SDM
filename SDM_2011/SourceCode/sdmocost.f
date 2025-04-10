      double precision function sdmocost(ngd,nobs,doffset)
      implicit none
c
c     calculate the slip gradient part of the cost function
c
c     Last modified: Potsdam, Oct, 2008, by R. Wang
c
      include 'sdmglob.h'
      integer ngd,nobs
      double precision doffset(NGDMAX)
c
      integer igd,iobs,iobs1,iobs2
      double precision ocost
c
      ocost=0.d0
      iobs1=0
      do igd=1,ngd
        iobs2=iobs1+nobsj(igd)
        do iobs=iobs1+1,iobs2
          ocost=ocost+wf(iobs)*(dspres(iobs)-doffset(igd)-ddsp(iobs))**2
        enddo
        iobs1=iobs2
      enddo
      sdmocost=0.5d0*ocost
      return
      end
