      subroutine dgwvint(nr1,nr2,accuracy)
      use dgalloc
      implicit none
c
      integer*4 nr1,nr2
      real*8 accuracy
c
c     u: 1=uz, 2=ur, 3=ut,
c     NOTE: uz, ur
c           have the same azimuth-factor as the poloidal mode (p-sv);
c           ut
c           have the same azimuth-factor as the
c           toroidal mode (sh);
c
      integer*4 i,istp,ir
      integer*4 ik,nk,nx,nr0
      real*8 x,k,k0,dk,dk0,dkmax
      real*8 wl,wr,zdis,udif,uabs,fac
      real*8 umax(4),usum
      real*8 lahs,muhs
      real*8 cics(4),cms(4),cbs(3),uk(3)
      real*8 y(8,4),cy0(3,4)
      logical*2 again
c
      integer*4 iret
      real*8 xokada,yokada,dip
      real*8 alpha,pot1,pot2,pot3,pot4,cs45,ss45
      real*8 ux,uy,uz,uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz
      real*8 uza,ura,uta
c
      integer*4 nkmin,nkmax
      real*8 eps,pi2
      data nkmin,nkmax/1024,65536/
      data eps,pi2/1.0d-03,6.28318530717959d0/
c
      cs45=1.d0/dsqrt(2.d0)
      ss45=1.d0/dsqrt(2.d0)
c
      zdis=dabs(zs-zrec)
      if(r(nr1).eq.0.d0)then
        nr0=nr1+1
      else
        nr0=nr1
      endif
      do ir=nr1,nr2
        r0(ir)=0.01d0*dsqrt(r(ir)*r(ir)+zdis*zdis)
      enddo
c
c     source functions
c
      call dgsource(1.d0)
c
c     ics = 1  when the azmuth-factor is cos(ms*theta) for poloidal mode
c             (psv) and sin(ms*theta) for the toroidal mode (sh);
c     ics = -1 otherwise.
c
      do istp=1,4
        cics(istp)=dble(ics(istp))
        cms(istp)=dble(ms(istp))
      enddo
c
      dkmax=pi2/dsqrt(r(nr2)**2+zdis**2)
c
      do istp=1,4
        do i=1,3
          do ir=nr1,nr2
            u(ir,i,istp)=0.d0
          enddo
        enddo
      enddo
c
      lahs=la(nno(lzrec))
      muhs=mu(nno(lzrec))
      call dgsource(1.d0)
      call dghssrce(1.d0,lahs,muhs)
c
      dk0=eps*dkmax
      k0=dk0
      do istp=1,4
        umax(istp)=0.d0
      enddo
50    again=.false.
      call dgkern(y,k0,lahs,muhs)
      fac=k0*dsqrt(k0)*dexp(-0.5d0*(k0*r0(nr1))**2)
      do istp=1,4
        uabs=0.d0
        do i=1,5,2
          uabs=uabs+dabs(y(i,istp))
        enddo
        uabs=uabs*fac
        umax(istp)=dmax1(umax(istp),uabs)
        if(uabs.gt.eps*umax(istp))then
          again=.true.
        endif
      enddo
      if(again)then
        dk0=1.1d0*dk0
        k0=k0+dk0
        goto 50
      endif
      usum=0.d0
      do istp=1,4
        usum=usum+umax(istp)
      enddo
      if(usum.le.0.d0)then
        nk=0
        k0=0.d0
        goto 400
      endif
c
      nk=1
      dk=k0
100   nk=2*nk
      dk=0.5d0*dk
c
      do istp=1,4
        do i=1,3
          do ir=nr1,nr2
            u0(ir,i,istp)=u(ir,i,istp)
          enddo
        enddo
      enddo
c
      do ik=1,nk-1,2
        k=dble(ik)*dk
        call dgkern(y,k,lahs,muhs)
        do istp=1,4
c
c         for displacement components
c
          cy0(1,istp)=y(1,istp)
          cy0(2,istp)=0.5d0*(y(3,istp)+cics(istp)*y(5,istp))
          cy0(3,istp)=0.5d0*(y(3,istp)-cics(istp)*y(5,istp))
        enddo
c
c       u1-3 are displacement components:
c
        if(nr0.eq.nr1+1)then
c
c         for distance r = 0
c
          ir=nr1
          fac=k*dexp(-0.25d0*(k*r0(ir))**2)
          do istp=1,4
            do i=1,3
              uk(i)=cy0(i,istp)*fac
            enddo
            if(ms(istp).eq.0)then
              u(ir,1,istp)=u(ir,1,istp)+uk(1)
            else if(ms(istp).eq.1)then
              u(ir,2,istp)=u(ir,2,istp)+uk(2)
              u(ir,3,istp)=u(ir,3,istp)-cics(istp)*uk(2)
            endif
          enddo
        endif
        do ir=nr0,nr2
	    x=k*r(ir)
          fac=dsqrt(k)*dexp(-0.5d0*(k*r0(ir))**2)
c
c	    bessels functions from pre-calculated tables
c
	    nx=idint(x/dxbsj)
	    wr=dmod(x/dxbsj,1.d0)
	    wl=1.d0-wr
	    if(nx.gt.nnbsj)nx=nnbsj+mod(nx-nnbsj,ndbsj)
          do istp=1,4
            cbs(1)=(wl*bsjfct(nx,ms(istp)-1)
     &             +wr*bsjfct(nx+1,ms(istp)-1))*fac
            cbs(2)=(wl*bsjfct(nx,ms(istp))
     &             +wr*bsjfct(nx+1,ms(istp)))*fac
            cbs(3)=(wl*bsjfct(nx,ms(istp)+1)
     &             +wr*bsjfct(nx+1,ms(istp)+1))*fac
c
            u(ir,1,istp)=u(ir,1,istp)+cy0(1,istp)*cbs(2)
            u(ir,2,istp)=u(ir,2,istp)
     &        +cy0(2,istp)*cbs(1)-cy0(3,istp)*cbs(3)
            u(ir,3,istp)=u(ir,3,istp)
     &        -cics(istp)*(cy0(2,istp)*cbs(1)+cy0(3,istp)*cbs(3))
          enddo
        enddo
      enddo
c
      if(nk.lt.nkmin)then
        goto 100
      else
        again=.false.
        do i=1,3
          udif=0.d0
          uabs=0.d0
          do istp=1,4
            do ir=nr1,nr2
              udif=udif+(zdis**2+r(ir)**2)
     &             *dabs(u(ir,i,istp)-2.d0*u0(ir,i,istp))
              uabs=uabs+(zdis**2+r(ir)**2)
     &             *dabs(u(ir,i,istp))
            enddo
          enddo
          again=again.or.udif.gt.accuracy*uabs
        enddo
        if(again.and.nk.lt.nkmax)then
          goto 100
        else if(again.and.nk.ge.nkmax)then
          nwarn=nwarn+1
        endif
      endif
c
      do ir=nr1,nr2
        if(r(ir).eq.0.d0)then
          fac=dk
        else
          fac=dk/dsqrt(r(ir))
        endif
        do istp=1,4
          do i=1,3
            u(ir,i,istp)=u(ir,i,istp)*fac
          enddo
        enddo
      enddo
400   continue
c=============================================================================
c     end of wavenumber integration
c=============================================================================
      write(*,'(a,i8,a,2f8.1)')' Wavenumber samples: ',nk,
     &                         ' x1, x2: ',k0*r(nr1),k0*r(nr2)
c
      return
      end
