      subroutine sdmresbat(ierr)
      use sdmalloc
      implicit none
      integer*4 ierr
c
c     calculate weighted smoothing matrix
c     Last modified: Berlin, Jan. 2026, by R. Wang
c
      integer*4 i,j,k,ips,ira,iobs
      real*8 slpabs,sig2,sigsum
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