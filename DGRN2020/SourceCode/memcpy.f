c*******************************************************************************
c*******************************************************************************
        subroutine memcpy(a,b,n)
        implicit none
c
c       copy real array a to b
c
        integer*4 n
        real*8 a(n),b(n)
c
        integer*4 i
c
        do i=1,n
          b(i)=a(i)
        enddo
        return
        end
