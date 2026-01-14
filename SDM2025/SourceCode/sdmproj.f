      subroutine sdmproj(ierr)
      use sdmalloc
      implicit none
      integer*4 ierr
c
c     Project on to subspace with positivity constraint
c
c     Last modified: Berlin, Jan. 2026, by R. Wang
c
      integer*4 i,j,k,is,ira,ips,iusrp
      real*8 a,b,slp,racs1,racs2,rass1,rass2,ra,racs,rass
c
      i=0
      do ips=1,nps
        do ira=1,2
          i=i+1
          slpmdl(ira,ips)=sysvec(i)/zhy(ips)
        enddo
      enddo
c
      do iusrp=1,nusrp
        corrusrp(iusrp)=sysvec(2*nps+iusrp)*usrpunit(iusrp)
      enddo
c
      do is=1,ns
        racs1=dcos(rake1(is)*DEG2RAD)
        rass1=dsin(rake1(is)*DEG2RAD)
        racs2=dcos(rake2(is)*DEG2RAD)
        rass2=dsin(rake2(is)*DEG2RAD)
        do ips=nps1(is),nps2(is)
          slp=dsqrt(slpmdl(1,ips)**2+slpmdl(2,ips)**2)
          if(rake1(is).eq.rake2(is))then
            slp=dmax1(0.d0,slpmdl(1,ips)*racs1+slpmdl(2,ips)*rass1)
            slpmdl(1,ips)=slp*racs1
            slpmdl(2,ips)=slp*rass1
          else if(slp.gt.0.d0.and.rake2(is)-rake1(is).lt.360.d0)then
            ra=datan2(slpmdl(2,ips),slpmdl(1,ips))
            if(rake360(is).gt.0.d0)then
              ra=dmod(ra+2.d0*PI,2.d0*PI)
            else
              ra=dmod(ra-2.d0*PI,2.d0*PI)
            endif
            if(ra/DEG2RAD.lt.rake1(is).or.ra/DEG2RAD.gt.rake2(is))then
              racs=dcos(ra)
              rass=dsin(ra)
              a=racs*racs1+rass*rass1
              b=racs*racs2+rass*rass2
              if(a.gt.0.d0.and.a.gt.b)then
                slpmdl(1,ips)=slp*a*racs1
                slpmdl(2,ips)=slp*a*rass1
              else if(b.gt.0.d0.and.b.gt.a)then
                slpmdl(1,ips)=slp*b*racs2
                slpmdl(2,ips)=slp*b*rass2
              else
                slpmdl(1,ips)=0.d0
                slpmdl(2,ips)=0.d0
              endif
            endif
          endif
          if(slp.gt.maxslip(is))then
            slpmdl(1,ips)=slpmdl(1,ips)*maxslip(is)/slp
            slpmdl(2,ips)=slpmdl(2,ips)*maxslip(is)/slp
          endif
        enddo
      enddo
      do iusrp=1,nusrp
        corrusrp(iusrp)=dmin1(usrpmax(iusrp),
     &                dmax1(usrpmin(iusrp),corrusrp(iusrp)))
      enddo
c
      i=0
      do ips=1,nps
        do ira=1,2
          i=i+1
          sysvec(i)=slpmdl(ira,ips)*zhy(ips)
        enddo
      enddo
      do iusrp=1,nusrp
        i=i+1
        sysvec(i)=corrusrp(iusrp)/usrpunit(iusrp)
      enddo
c
      sysmis=0.d0
      do i=1,nsys
        resbat(i)=-sysbat(i)
        do j=1,nsys
          resbat(i)=resbat(i)+sysmat(i,j)*sysvec(j)
        enddo
        sysmis=sysmis+sysvec(i)*(resbat(i)-sysbat(i))
      enddo
      sysmis=1+sysmis/datnrm
c
      return
      end