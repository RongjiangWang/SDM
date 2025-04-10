c     GLOBAL INDEX PARAMETERS FOR DEFINING ARRAYS
c     ===========================================
c      nzmax: max. interface index;
c      lmax: max. no of total homogeneous layers (lmax <= nzmax-2);
c      nrmax: max. no of traces;
c      nzsmax: max. number of source depths
c
      integer lmax,nzmax,nrmax,nzsmax
      parameter(lmax=101)
      parameter(nzmax=lmax+3)
      parameter(nrmax=501)
      parameter(nzsmax=100)
c
c     INDEX PARAMETERS FOR BESSEL FUNCTION TABLES
c     ===========================================
c
      integer nbsj,ndbsj
      parameter(nbsj=512,ndbsj=512)
      integer nnbsj,nnbsj1
      parameter(nnbsj=nbsj*ndbsj,nnbsj1=nnbsj+ndbsj)
c
c     GLOBAL CONSTANTS
c     ================
c
      double precision km2m,day2sec,relaxmin
      parameter(km2m=1.0d+03,day2sec=8.64d+04,relaxmin=1.0d-06)
c
c     GRAVITY, GRAVITATIONAL CONSTANT AND EARTH RADIUS
c     ================================================
c     gamma = 4*pi*G
c
      double precision g0,gamma,rearth
      parameter(g0=9.82d+00,gamma=8.38579d-10,rearth=6.371d+06)
c
      integer nwarn
      common /warnings/nwarn
c
c     DISCRETISATION ACCURACY FOR LAYERS WITH CONSTANT GRADIENT
c     =========================================================
c     reslm: for moduli
c     resld: for density
c
      double precision reslm,resld
      parameter(reslm=0.05d0,resld=0.05d0)
c
c     COMMON BLOCKS
c     =============
      integer lp,nno(nzmax)
      double precision hp(nzmax)
      common /sublayer/ hp,lp,nno
c
c     zrec: receiver depth
c     lzrec: sublayer no of receiver
c
      integer lzrec
      double precision zrec
      common /receiver/ zrec,lzrec
c
c     original model parameters
c
      integer l0
      double precision z1(lmax),z2(lmax)
      double precision la1(lmax),la2(lmax),mu1(lmax),mu2(lmax)
      double precision rho1(lmax),rho2(lmax)
      common /model0/z1,z2,la1,la2,mu1,mu2,rho1,rho2,l0
c       
c     model parameter:
c     n0: number of homogeneous layers
c
      integer n0
      double precision h(lmax),la(lmax),mu(lmax),rho(lmax)
      common /model/ h,la,mu,rho,n0
c
c     source parameters
c
      integer ls,ms(4),ics(4)
      double precision zs
      double precision sfct0(8,4),sfct1(8,4)
      common /source/ zs,sfct0,sfct1,ls,ms,ics
c
c     half-space source parameters
c
      double precision sfcths0(8,4),sfcths1(8,4)
      common /sourcehs/ sfcths0,sfcths1
c
c     table of J_n(x), n = -1, 0, 1, 2, 3
c     all multiplied by sqrt(x)
c
      double precision dxbsj,bsjfct(0:nnbsj1,-1:3)
      common /bessels/ dxbsj,bsjfct
c
c     psv layer matrics
c
      double precision maup(6,6,nzmax),maiup(6,6,nzmax)
      double precision malw(6,6,nzmax),mailw(6,6,nzmax)
      common /psvlayma/ maup,maiup,malw,mailw
c
c     sh layer hask matrices
c
      double precision hkup(2,2,nzmax),hklw(2,2,nzmax)
      common /shhask/ hkup,hklw
c
c     output data
c
      double precision r(nrmax)
      double precision u(nrmax,3,4)
      common /outdata/ r,u

