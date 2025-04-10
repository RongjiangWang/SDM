c*******************************************************************************
c*******************************************************************************
        subroutine gemp(a,b,n,n1,eps,key)
        implicit none
c-------------------------------------------------------------------------------
c       subroutine gemp to solve linear equation system                        i
c       a: coefficient matrix(n,n);                                            i
c       b: right-hand matrix(n,n1) by input,                                   i
c          solution matrix(n,n1) by return;                                    i
c       eps: control constant;                                                 i
c       key: the main term of a column                                         i
c            smaller than eps, key=0: anormal return,                          i
c            else key=1: normal return.                                        i
c-------------------------------------------------------------------------------
        integer n,n1,key
        double precision eps
        double precision a(n,n),b(n,n1)
c
        integer i,j,k,l,m
        double precision p
        double precision q
c
        do m=1,n
          p=0.d0
          do i=m,n
            do j=1,n
              if(dabs(a(i,j)).le.p)goto 10
              p=dabs(a(i,j))
              k=i
              l=j
10          continue
            enddo
          enddo
          key=0
          if(p.le.eps)return
          do j=1,n
            q=a(k,j)
            a(k,j)=a(m,j)
            a(m,j)=q
          enddo
          do j=1,n1
            q=b(k,j)
            b(k,j)=b(m,j)
            b(m,j)=q
          enddo
          do i=1,n
            q=-a(i,l)/a(m,l)
            if(i.eq.m)goto 20
            do j=1,n
              a(i,j)=a(i,j)+a(m,j)*q
            enddo
            do j=1,n1
              b(i,j)=b(i,j)+b(m,j)*q
            enddo
20          continue
          enddo
          q=a(m,l)
          do j=1,n
            a(m,j)=a(m,j)/q
          enddo
          do j=1,n1
            b(m,j)=b(m,j)/q
          enddo
        enddo
        do j=1,n
          do i=1,n
            if(a(i,j).lt.0.5d0)goto 30
            do k=1,n
              q=a(i,k)
              a(i,k)=a(j,k)
              a(j,k)=q
            enddo
            do k=1,n1
              q=b(i,k)
              b(i,k)=b(j,k)
              b(j,k)=q
            enddo
30          continue
          enddo
        enddo
        key=1
        return
        end
