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
      integer NSMAX,NWFTMAX,NPSMAX,NXYMAX,NZSMAX,NRMAX
      integer NOBSMAX,NGDMAX
      parameter(NSMAX=10,NWFTMAX=100,NPSMAX=2000)
      parameter(NZSMAX=101,NRMAX=251)
      parameter(NOBSMAX=5000)
      parameter(NGDMAX=20)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     CURVED FAULT SURFACES
c     =====================
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer idisc
	  integer nwft(NSMAX),ndgrid(NSMAX)
	  integer nps1(NSMAX),nps2(NSMAX),iref(NSMAX)
	  double precision rake1(NSMAX),rake2(NSMAX),maxslip(NSMAX)
      double precision topdep(NSMAX),btmdep(NSMAX),patchsize(NSMAX)
      double precision toplat(NWFTMAX,NSMAX),btmlat(NWFTMAX,NSMAX)
	  double precision toplon(NWFTMAX,NSMAX),btmlon(NWFTMAX,NSMAX)
      double precision topdip(NWFTMAX,NSMAX)
c
      double precision ramean(NSMAX),rake360(NSMAX)
      double precision cs1(NSMAX),ss1(NSMAX),cs2(NSMAX),ss2(NSMAX)
      double precision measdrop(NSMAX),stdsdrop(NSMAX),maxsdrop(NSMAX)
c
      common/irects/idisc,nwft,ndgrid,nps1,nps2,iref
c
      common/drects/patchsize,topdep,btmdep,rake1,rake2,maxslip,
     &              toplat,btmlat,toplon,btmlon,topdip,
     &              ramean,rake360,
     &              cs1,ss1,cs2,ss2,measdrop,stdsdrop,maxsdrop
c
      double precision poisson,dspunit,sarea,costref,scostref
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
      integer iz(NPSMAX)
      integer ipsl(NPSMAX),ipsr(NPSMAX),ipsu(NPSMAX),ipsd(NPSMAX)
      double precision plat(NPSMAX),plon(NPSMAX),pz(NPSMAX)
      double precision d2l(NPSMAX),d2r(NPSMAX)
      double precision d2t(NPSMAX),d2b(NPSMAX)
      double precision dlen(NPSMAX),dwid(NPSMAX),parea(NPSMAX)
      double precision strike(NPSMAX),dip(NPSMAX)
      double precision pmwei(5,NPSMAX,2)
      double precision strgrn(NPSMAX,2,NPSMAX,3)
      double precision dcgrn(NPSMAX,2,NPSMAX,6)
      double precision strdrop(NPSMAX,3),strdc(NPSMAX,6)
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
      logical hsmodel
      integer nr,nzs
      double precision zobs,laobs,muobs
      double precision r(NRMAX),zs(NZSMAX),muz(NZSMAX)
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
      double precision grns(NZSMAX,NRMAX,3,3)
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
      logical caloffset,csconst(NGDMAX)
      integer nheader
      integer nobsj(NGDMAX),seloffset(NGDMAX)
      double precision dinc_const(NGDMAX),dazi_const(NGDMAX)
      double precision xcs(NOBSMAX),ycs(NOBSMAX),zcs(NOBSMAX)
      double precision offset(NGDMAX)
      double precision latobs(NOBSMAX),lonobs(NOBSMAX)
      double precision dspobs(NOBSMAX),dspres(NOBSMAX)
      double precision wf(NOBSMAX)
      double precision wfsum
      double precision wfm(NGDMAX),wfmsum(NGDMAX)
c
      common/incidence/caloffset,csconst
      common/gdarno/nheader,nobsj,seloffset
      common/projection/dinc_const,dazi_const,xcs,ycs,zcs
      common/positions/offset,latobs,lonobs
      common/datasets/dspobs,dspres
      common/weighting/wfm,wfmsum,wf,wfsum
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     SLIP MODEL
c     ==========
c
c     slpmdl(2,..) = slip components in the strike and dip directions
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer slppos(NPSMAX)
      double precision slpmdl(NPSMAX,2)
      double precision dspmdl(NPSMAX,NOBSMAX,2),ddsp(NOBSMAX)
c
      common/islipm/slppos
      common/dslipm/slpmdl,dspmdl,ddsp
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     WARNING STATISTICS
c     ==================
c     nwarn = total number of warnings
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer nwarn
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
      double precision DEG2RAD,KM2M,DAY2SEC,REARTH,G0
      parameter(DEG2RAD=1.745329252d-02,KM2M=1.0d+03)
      parameter(DAY2SEC=8.64d+04,REARTH=6.371d+06,G0=9.82d+00)
      double precision LAMREF,MUEREF
      parameter(LAMREF=3.0d+10,MUEREF=3.0d+10)
