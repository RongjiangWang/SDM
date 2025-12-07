      subroutine sdmpatches(ierr)
      use sdmalloc
      implicit none
      integer*4 ierr
c
c     allocate arrays for slip patches
c     Last modified: Zhuhai, Nov. 2025, by R. Wang
c      
      allocate(plat(nps),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmdisc: plat not allocated!'
      allocate(plon(nps),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmdisc: plon not allocated!'
      allocate(pz(nps),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmdisc: pz not allocated!'
      allocate(pl(nps),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmdisc: pl not allocated!'
      allocate(pw(nps),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmdisc: dwid not allocated!'
      allocate(dlen(nps),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmdisc: pw not allocated!'
      allocate(dwid(nps),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmdisc: dwid not allocated!'
      allocate(parea(nps),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmdisc: parea not allocated!'
      allocate(strike(nps),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmdisc: strike not allocated!'
      allocate(dip(nps),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmdisc: dip not allocated!'
      allocate(pmwei(5,nps,2),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmdisc: pmwei not allocated!'
      allocate(ipsl(nps),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmdisc: ipsl not allocated!'
      allocate(ipsr(nps),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmdisc: ipsr not allocated!'
      allocate(ipsu(nps),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmdisc: ipsu not allocated!'
      allocate(ipsd(nps),stat=ierr)
      if(ierr.ne.0)stop ' Error in sdmdisc: ipsd not allocated!'
c
      return
      end