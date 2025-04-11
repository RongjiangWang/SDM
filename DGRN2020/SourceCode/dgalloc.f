      module dgalloc
c
c     INDEX PARAMETERS FOR BESSEL FUNCTION TABLES
c     ===========================================
c
      integer*4 nbsj,ndbsj
      parameter(nbsj=512,ndbsj=512)
      integer*4 nnbsj,nnbsj1
      parameter(nnbsj=nbsj*ndbsj,nnbsj1=nnbsj+ndbsj)
c
c     table of J_n(x), n = -1, 0, 1, 2, 3
c     all multiplied by sqrt(x)
c
      real*8 bsjfct(0:nnbsj1,-1:3)
c
c     GLOBAL CONSTANTS
c     ================
c
      real*8 km2m,day2sec,relaxmin
      parameter(km2m=1.0d+03,day2sec=8.64d+04,relaxmin=1.0d-06)
c
c     GRAVITY, GRAVITATIONAL CONSTANT AND EARTH RADIUS
c     ================================================
c     gamma = 4*pi*G
c
      real*8 g0,gamma,rearth
      parameter(g0=9.82d+00,gamma=8.38579d-10,rearth=6.371d+06)
c
c     DISCRETISATION ACCURACY FOR LAYERS WITH CONSTANT GRADIENT
c     =========================================================
c     reslm: for moduli
c     resld: for density
c
      real*8 reslm,resld
      parameter(reslm=0.05d0,resld=0.05d0)
c
      integer*4 nwarn,lzrec,ls,l0,lp,n0
      real*8 zrec,zs,dxbsj
c
      integer*4 ms(4),ics(4)
      real*8 sfct0(8,4),sfct1(8,4),sfcths0(8,4),sfcths1(8,4)
c
c     COMMON BLOCKS
c     =============
      integer*4, allocatable:: nno(:)
      real*8, allocatable:: h0(:),la0(:),mu0(:),rho0(:)
      real*8, allocatable:: z1(:),la1(:),mu1(:),rho1(:)
      real*8, allocatable:: z2(:),la2(:),mu2(:),rho2(:)
      real*8, allocatable:: h(:),la(:),mu(:),rho(:)
      real*8, allocatable:: xh(:),xh0(:),hp(:)
      real*8, allocatable:: maup(:,:,:),maiup(:,:,:)
      real*8, allocatable:: malw(:,:,:),mailw(:,:,:)
      real*8, allocatable:: hkup(:,:,:),hklw(:,:,:)
      real*8, allocatable:: r(:),r0(:),u(:,:,:),u0(:,:,:)
      real*8, allocatable:: zs0(:)
c
      end module
