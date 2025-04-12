      real*8 function fdip(dipup,x,z)
      implicit none
c
      real*8 dipup,x,z
c
c     calculate dip angle [deg] at the bottom edge of a
c     curved fault by assuming that the dip varies linearly
c     with depth
c     input parameters:
c     dipup = dip at upper edge [deg]
c     x = horizontal length of the fault
c     z = vertical length of the fault
c
c     Last modified: Potsdam, March, 2013, by R. Wang
c
      integer*4 i
      real*8 deg2rad,di0,di1,dil,dir,dim,fl,fr,fm
c
      integer*4 imax
      real*8 eps
      data imax,eps/10000,1.0d-06/
c
      if(x.le.0.d0.and.z.le.0.d0)then
        stop 'Error in fdip: bad parameters x and z!'
      endif
c
      deg2rad=datan(1.d0)/45.d0
c
      if(x.le.0.d0)then
        di0=90.d0*deg2rad
      else
        di0=datan(z/x)
      endif
      di1=dipup*deg2rad
      if(dabs(di1-di0).le.eps)then
        fdip=di1/deg2rad
        return
      endif
c
      if(di1.ge.di0)then
c
c       listric fault
c
        dil=0.d0
        dir=di0
        fl=1.d+06
        fr=dir-di1-(z/x)*dlog(dsin(dir)/dsin(di1))
      else
c
c       opposite listric fault
c
        dil=di0
        dir=2.d0*datan(1.d0)
        fl=dil-di1-(z/x)*dlog(dsin(dil)/dsin(di1))
        fr=dir-di1-(z/x)*dlog(dsin(dir)/dsin(di1))
      endif
      if(fl*fr.gt.0.d0)then
        stop 'Error in fdip: a problem appears.'
      endif
      dim=0.5d0*(dil+dir)
      fm=dim-di1-(z/x)*dlog(dsin(dim)/dsin(di1))
c
      i=0
c
10    i=i+1
      if(fm*fl.le.0.d0)then
        dir=dim
        fr=fm
      else
        dil=dim
        fl=fm
      endif
      dim=0.5d0*(dil+dir)   !(dil*fr-dir*fl)/(fr-fl)
      fm=dim-di1-(z/x)*dlog(dsin(dim)/dsin(di1))
      if(dabs(fm).gt.eps.and.i.lt.imax)goto 10
c
      fdip=dim/deg2rad
      return
      end
