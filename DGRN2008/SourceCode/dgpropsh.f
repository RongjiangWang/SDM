      subroutine dgpropsh(l1,l2,k,ysh,ysh0)
      implicit none
c
c	propagation of sh vectors
c
      integer l1,l2
      double precision k
      double precision ysh(2),ysh0(2)
c
c
	include 'dgglob.h'
c
c     work space
c
      integer l
      double precision cnorm,yswab(2)
c
	if(l1.eq.l2)then
	  return
	else if(l1.lt.l2)then
        do l=l1+1,l2
c
c         determination of propagation matrix
c
	    call axb(hkup(1,1,l-1),ysh,2,2,1,yswab)
	    call memcpy(yswab,ysh,2)
          if(l.gt.lzrec)then
            cnorm=dexp(-k*hp(l-1))
            ysh0(1)=ysh0(1)*cnorm
            ysh0(2)=ysh0(2)*cnorm
	    else if(l.eq.lzrec)then
            call memcpy(ysh,ysh0,2)
          endif
        enddo
	else
        do l=l1-1,l2,-1
c
c         determination of propagation matrix
c
	    call axb(hklw(1,1,l),ysh,2,2,1,yswab)
	    call memcpy(yswab,ysh,2)
          if(l.lt.lzrec)then
            cnorm=dexp(-k*hp(l))
            ysh0(1)=ysh0(1)*cnorm
            ysh0(2)=ysh0(2)*cnorm
	    else if(l.eq.lzrec)then
            call memcpy(ysh,ysh0,2)
          endif
        enddo
	endif
	return
	end
