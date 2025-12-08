      module sdmalloc
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     GLOBAL CONSTANTS
c     ================
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      real*8 DEG2RAD,PI,KM2M,DAY2SEC,REARTH,G0,MUEREF,MEGA
      parameter(DEG2RAD=1.745329252d-02,PI=3.14159265358979d0)
      parameter(KM2M=1.0d+03,DAY2SEC=8.64d+04,REARTH=6.371d+06)
      parameter(G0=9.82d+00,MUEREF=3.0d+10,MEGA=1.0d+06)
      integer*4 nrelaxmax
      parameter(nrelaxmax=10000)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     OBSTANGULAR SOURCE PLANES
c     =========================
c     with x = north, y = east, z = downward.
c     all angles in degree.
c     NS = number of source rectangles
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer*4 ns,nftsum
      integer*4, allocatable:: nft(:),nlength(:),nwidth(:)
      integer*4, allocatable:: nps1(:),nps2(:),ipsp(:)
      real*8, allocatable:: topdep(:),width(:)
      real*8, allocatable:: mstrike(:),patchsize(:)
      real*8, allocatable:: latft(:,:),lonft(:,:)
      real*8, allocatable:: topdip(:,:),botdip(:,:)
      real*8, allocatable:: rake1(:),rake2(:),maxslip(:)
      real*8, allocatable:: ram(:),rake360(:),ssm(:),dsm(:),slpm(:)
      real*8, allocatable:: slpp(:),rap(:),cs1(:),ss1(:),cs2(:),ss2(:)
      real*8, allocatable:: measdrop(:),stdsdrop(:),maxsdrop(:)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     DISTRETE POINT SOURCES
c     ======================
c     (xs,ys,zs) = coordinates of the discrete point sources
c     with x = north, y = east, z = downward
c     angles in degree.
c     strgrn(ips,i,jps,j): j-th component of stress drop at patch jps
c                      induced by i-th component of slip at patch ips
c     divgrn(ips,i,jps): stress drop divergence at patch jps
c                      induced by i-th component of slip at patch ips
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer*4 nps,nftmax,idisc
      real*8 poisson,lamhs,muehs,slparea,roughness
      integer*4, allocatable:: iz(:),isref(:)
      integer*4, allocatable:: ipsl(:),ipsr(:),ipsu(:),ipsd(:)
      real*8, allocatable:: plat(:),plon(:)
      real*8, allocatable:: pl(:),pw(:),pz(:)
      real*8, allocatable:: d2l(:),d2r(:)
      real*8, allocatable:: d2t(:),d2b(:)
      real*8, allocatable:: dlen(:),dwid(:),parea(:)
      real*8, allocatable:: strike(:),dip(:)
      real*8, allocatable:: pmwei(:,:,:)
      real*8, allocatable:: strgrn(:,:,:,:)
      real*8, allocatable:: dcgrn(:,:,:,:)
      real*8, allocatable:: strdrop(:,:),strdc(:,:)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     GREEN'S FUNCTION INFO
c     =====================
c     zs = depth samples
c     r = distance samples
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      logical*2 hsmodel
      integer*4 nr,nzs
      real*8 zobs,laobs,muobs
      real*8, allocatable:: r(:),zs(:),muz(:)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     GREEN'S FUNNCTIONN PARAMETERS
c     =============================
c     Green's function source types:
c       1 = strike-slip (m12=m21=1)
c       2 = dip-slip (m13=m31=1)
c       3 = compensated linear vector dipole (CLVD)
c           (m11=m22=-1/2, m33=1) (no tangential component)
c     Green's function coordinate system:
c       (z,r,t) = cylindrical with z being downward(!)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      real*8, allocatable:: dgrns(:,:,:,:)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     OBSERVATION POSITIONS AND GEODETIC DATA
c     =======================================
c     (xobs(i),yobs(i))=coordinates of the observation positions
c     the 3 displcement components: ux,uy,uz
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer*4 nheader,ngd,nobs,npar
      real*8 datunit,offunit
      logical*2, allocatable:: csconst(:)
      integer*4, allocatable:: nobs1(:),nobs2(:),nobsj(:)
      real*8, allocatable:: dinc_const(:),dazi_const(:)
      real*8, allocatable:: xcs(:),ycs(:),zcs(:)
      real*8, allocatable:: do0(:),dm0(:),corrgrn(:,:)
      real*8, allocatable:: latobs(:),lonobs(:),parunit(:)
      real*8, allocatable:: datobs(:),datres(:),corrpar(:)
      real*8, allocatable:: wf(:),wfm(:),parmin(:),parmax(:)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     SLIP MODEL
c     ==========
c     slpmdl(2,..) = slip components in the strike and dip directions
c     zhang(ips) = weighted smoothing of roughness (Yong Zhang, 2025)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer*4 iearth,iter,niter,nsys,ismooth,izhy,nsmocmp
      real*8 wei2smo0,wei2smo,sig2sys,step,step0,mwpsum,mwssum,mweq
      real*8 rmsresall,rmsdatall,sysmis,sysmis0
      real*8, allocatable:: slpmdl(:,:),vecswp(:),matswp(:,:)
      real*8, allocatable:: datgrn(:,:,:),datmdl(:),batswp(:)
      real*8, allocatable:: sysmat(:,:),sysbat(:),sysvec(:),zhy(:)
      real*8, allocatable:: obsmat(:,:),smomat(:,:),obsswp(:)
      real*8, allocatable:: rmsres(:),resmin(:),resmax(:),corrmdl(:)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     WARNING STATISTICS
c     ==================
c     nwarn = total number of warnings
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer*4 nwarn
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     FILE NAMES
c     ==========
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      character*80 grndir,infile,slipout,logfile,zhyrelax
      character*80 usr3dgrn(2),green(3),corrgrnfile,parout
      character*80, allocatable:: gddata(:),gdout(:)
      character*80, allocatable:: inpatches(:)
c
      end module
