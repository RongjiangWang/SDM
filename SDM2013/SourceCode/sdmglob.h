c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     GLOBAL CONSTANTS
c     ================
c
c     NOBSMAX = max. number of (filtered) data points
c     NSMAX = max. number of discontinuous fault segments
c     NZSMAX = max. number of the discrete source depths
c     NRMAX = max. number of the discrete radial diatances
c     NPSMAX = max. number of discrete point sources per source depth
c     NGDMAX = max. number of geodetic data sets
c     NWFTMAX = max. number of reference wefts of a curved fault segment
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer*4 NSMAX,NWFTMAX,NPSMAX,NXYMAX,NZSMAX,NRMAX
      integer*4 NOBSMAX,NGDMAX
      parameter(NSMAX=10,NWFTMAX=100,NPSMAX=2000)
      parameter(NZSMAX=101,NRMAX=251)
      parameter(NOBSMAX=5000)
      parameter(NGDMAX=20)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     CURVED FAULT SURFACES
c     =====================
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer*4 idisc
	  integer*4 nwft(NSMAX),ndgrid(NSMAX)
	  integer*4 nps1(NSMAX),nps2(NSMAX),iref(NSMAX)
	  real*8 rake1(NSMAX),rake2(NSMAX),maxslip(NSMAX)
      real*8 topdep(NSMAX),btmdep(NSMAX),patchsize(NSMAX)
      real*8 toplat(NWFTMAX,NSMAX),btmlat(NWFTMAX,NSMAX)
	  real*8 toplon(NWFTMAX,NSMAX),btmlon(NWFTMAX,NSMAX)
      real*8 topdip(NWFTMAX,NSMAX)
c
      real*8 ramean(NSMAX),rake360(NSMAX)
      real*8 cs1(NSMAX),ss1(NSMAX),cs2(NSMAX),ss2(NSMAX)
      real*8 measdrop(NSMAX),stdsdrop(NSMAX),maxsdrop(NSMAX)
c
      common/irects/idisc,nwft,ndgrid,nps1,nps2,iref
c
      common/drects/patchsize,topdep,btmdep,rake1,rake2,maxslip,
     &              toplat,btmlat,toplon,btmlon,topdip,
     &              ramean,rake360,
     &              cs1,ss1,cs2,ss2,measdrop,stdsdrop,maxsdrop
c
      real*8 poisson,dspunit,sarea,costref,scostref
c
      common/dglob/poisson,dspunit,sarea,costref,scostref
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     DISTRETE POINT SOURCES
c     ======================
c
c     (xs,ys,zs) = coordinates of the discrete point sources
c     with x = north, y = east, z = downward
c     angles in degree.
c     strgrn(ips,i,jps,j): j-th component of stress drop at patch jps
c                      induced by i-th component of slip at patch ips
c     divgrn(ips,i,jps): stress drop divergence at patch jps
c                      induced by i-th component of slip at patch ips
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer*4 iz(NPSMAX)
      integer*4 ipsl(NPSMAX),ipsr(NPSMAX),ipsu(NPSMAX),ipsd(NPSMAX)
      real*8 plat(NPSMAX),plon(NPSMAX),pz(NPSMAX)
      real*8 d2l(NPSMAX),d2r(NPSMAX)
      real*8 d2t(NPSMAX),d2b(NPSMAX)
      real*8 dlen(NPSMAX),dwid(NPSMAX),parea(NPSMAX)
      real*8 strike(NPSMAX),dip(NPSMAX)
      real*8 pmwei(5,NPSMAX,2)
      real*8 strgrn(NPSMAX,2,NPSMAX,3)
      real*8 dcgrn(NPSMAX,2,NPSMAX,6)
      real*8 strdrop(NPSMAX,3),strdc(NPSMAX,6)
c
      common/ipoints/iz,ipsl,ipsr,ipsu,ipsd
      common/dpoints/plat,plon,pz,d2l,d2r,d2t,d2b,
     &               dlen,dwid,parea,strike,dip,pmwei,
     &               strgrn,dcgrn,strdrop,strdc
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     GREEN'S FUNCTION INFO
c     =====================
c
c     zs = depth samples
c     r = distance samples
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      logical*2 hsmodel
      integer*4 nr,nzs
      real*8 zobs,laobs,muobs
      real*8 r(NRMAX),zs(NZSMAX),muz(NZSMAX)
c
      common /lgreeninfo/hsmodel
      common /greeninfo/zobs,laobs,muobs,r,zs,muz
      common /igreeninfo/nr,nzs
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     GREEN'S FUNNCTIONN PARAMETERS
c     =============================
c
c     Green's function source types:
c       1 = strike-slip (m12=m21=1)
c       2 = dip-slip (m13=m31=1)
c       3 = compensated linear vector dipole (CLVD)
c           (m11=m22=-1/2, m33=1) (no tangential component)
c     Green's function coordinate system:
c       (z,r,t) = cylindrical with z being downward(!)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      real*8 grns(NZSMAX,NRMAX,3,3)
c
      common/greenfcts/grns
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     OBSERVATION POSITIONS AND GEODETIC DATA
c     =======================================
c
c     (xobs(i),yobs(i))=coordinates of the observation positions
c     the 3 displcement components: ux,uy,uz
c     NOBSMAX = the max. number of (filtered) data points
c     NOBS0MAX = the max. number of (unfiltered) data points
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      logical*2 caloffset,csconst(NGDMAX)
      integer*4 nheader
      integer*4 nobsj(NGDMAX),seloffset(NGDMAX)
      real*8 dinc_const(NGDMAX),dazi_const(NGDMAX)
      real*8 xcs(NOBSMAX),ycs(NOBSMAX),zcs(NOBSMAX)
      real*8 offset(NGDMAX),do0(NGDMAX),dm0(NGDMAX)
      real*8 latobs(NOBSMAX),lonobs(NOBSMAX)
      real*8 dspobs(NOBSMAX),dspres(NOBSMAX)
      real*8 wf(NOBSMAX)
      real*8 wfsum
      real*8 wfm(NGDMAX),wfmsum(NGDMAX)
c
      common/incidence/caloffset,csconst
      common/gdarno/nheader,nobsj,seloffset
      common/projection/dinc_const,dazi_const,xcs,ycs,zcs
      common/positions/offset,latobs,lonobs
      common/datasets/dspobs,dspres,do0,dm0
      common/weighting/wfm,wfmsum,wf,wfsum
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     SLIP MODEL
c     ==========
c
c     slpmdl(2,..) = slip components in the strike and dip directions
c     zhang(ips) = weighted smoothing of roughness (Yong Zhang, 2025)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer*4 slppos(NPSMAX)
      real*8 slpmdl(NPSMAX,2)
      real*8 dspmdl(NPSMAX,NOBSMAX,2),ddsp(NOBSMAX)
      real*8 zhy(NPSMAX)
c
      common/islipm/slppos
      common/dslipm/slpmdl,dspmdl,ddsp,zhy
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     WARNING STATISTICS
c     ==================
c     nwarn = total number of warnings
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer*4 nwarn
c
      common/warnings/nwarn
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     FILE NAMES
c     ==========
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      character*80 grndir,infile,inpatches(NSMAX)
      character*80 green(3),gddata0(NGDMAX)
      character*80 gddata(NGDMAX),gdout0(NGDMAX),slipout
c
      common/memories/grndir,infile,inpatches,green,
     +                gddata0,gddata,gdout0,slipout
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     LOCAL CONSTANTS
c     ===============
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer*4 izhy
	  parameter(izhy=1)
      real*8 DEG2RAD,KM2M,DAY2SEC,REARTH,G0
      parameter(DEG2RAD=1.745329252d-02,KM2M=1.0d+03)
      parameter(DAY2SEC=8.64d+04,REARTH=6.371d+06,G0=9.82d+00)
      real*8 LAMREF,MUEREF
      parameter(LAMREF=3.0d+10,MUEREF=3.0d+10)
