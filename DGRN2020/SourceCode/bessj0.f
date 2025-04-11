      real*8 function bessj0(x)
      implicit none
c
c     J_0(x)
c
      real*8 x,ax,x1,x2,pi,theta,fct
      real*8 a0,a1,a2,a3,a4,a5,a6
      real*8 b0,b1,b2,b3,b4,b5,b6
      real*8 c0,c1,c2,c3,c4,c5,c6
c
      data a0,a1,a2,a3,a4,a5,a6/1.00000000d0,
     + -2.24999970d0, 1.26562080d0,-0.31638660d0,
     +  0.04444790d0,-0.00394440d0, 0.00021000d0/
      data b0,b1,b2,b3,b4,b5,b6/0.79788456d0,
     + -0.00000077d0,-0.00552740d0,-0.00009512d0,
     +  0.00137237d0,-0.00072805d0, 0.00014476d0/
      data c0,c1,c2,c3,c4,c5,c6/-0.78539816d0,
     + -0.04166397d0,-0.00003954d0, 0.00262573d0,
     + -0.00054125d0,-0.00029333d0, 0.00013558d0/
c
      ax=dabs(x)
      if(ax.le.3.d0)then
        x2=(ax/3.d0)**2
        bessj0=a0+x2*(a1+x2*(a2+x2*(a3+x2*(a4+x2*(a5+x2*a6)))))
      else
        pi=4.d0*datan(1.d0)
        x1=3.d0/ax
        fct=b0+x1*(b1+x1*(b2+x1*(b3+x1*(b4+x1*(b5+x1*b6)))))
        theta=ax+c0+x1*(c1+x1*(c2+x1*(c3+x1*(c4+x1*(c5+x1*c6)))))
        theta=dmod(theta,2.d0*pi)
        bessj0=fct*dcos(theta)/dsqrt(ax)
      endif
c
      return
      end
