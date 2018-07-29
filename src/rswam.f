!-------------------------------------------------------------------------------
! RHYDRA: Version of hydra hydrological model to be called from R
!
! This uses the simpler SWAM description for monthly averages and an
! equilibrium run
!
! Ver. 0.1 Simple template to read in a DEM
! Ver. 0.2 Transferred code from hydra.f excluding the HDF/Netcdf parts
! Ver. 0.3 
!-------------------------------------------------------------------------------

      subroutine rswam( nyrs, startyear, converg, laket, spin,
     *                     normal, leap, irrig,
     *                     outnewi, outnewj, basin, 
     *                     dem, outdir, sillh,
     *                     prcpi, evapi, runin, drainin,
     *                     outelv, lakem, 
     *                     lakevolm, lakevola )
              
      ! Input variables
      integer sampi,sampj
      integer nc,nr,nrf,ncf ! Input grid dimensions coarse/fine
      integer istart,iend,jstart,jend ! Grid limits (coarse)
      integer i2start,i2end,j2start,j2end ! Fine limits
      integer i3start,i3end,j3start,j3end
      integer inc,inr,incf,inrf ! Number of rows/cols (region)
      integer nmons,spin ! Number of months/spinup

      double precision n,mean,sum
      double precision dgper ! Degrees per coarse grid
      double precision scale2 ! Number of fine cells per coarse
      double precision totarea ! Total area of study region
      double precision evapm ! Used to calculate basin evaporation
      double precision evapt ! Used to calculate total evaporation
      double precision lareat ! Used to calculate total lake area
 
      double precision circ,dy,dx,pi,rad,phi
      double precision delt,res,ic,io,ioo ! time step, velocity variable
      double precision grideps,dveps,gridif

      ! Parameters
      parameter (grideps = 1.e-10,dveps = 1.e-10)
c     parameter (grideps = 0.01,dveps = 1.e-06)
      parameter (circ=4.0024E+7, pi = 3.1415926536)
      parameter (rad = pi/180., dy = circ*((5./60.)/360.))
c
c Set up the full grid array for the 5' resolution (ncf,nrf) and
c the 1/2 degree input resolution (nc,nr)
c i2start, j2start etc. are the array bounds on 1/2 resolution
c i3start j3end etc. convert the 1/2 to 5' resolution
c
c The first two parameters here will need to be changed 
      parameter (dgper = 0.5)  !number of degrees/input grid cell (1/2x1/2 here)
      parameter (scale2 = 12*dgper) !number of 5' cells per grid cell
      parameter (i2start = (0+1),i2end = 720)    ! Globe
      parameter (j2start = (0+1),j2end = 360) ! for 0.5 deg
      parameter (i3start = (scale2*(i2start-1)+1),
     *           i3end = (scale2*i2end))                
      parameter (j3start = (scale2*(j2start-1)+1),
     *           j3end = (scale2*j2end))              
      parameter (nr = (j2end+1)-j2start, nc = (i2end+1)-i2start)  ! set array
      parameter (nrf = (j3end+1)-j3start, ncf = (i3end+1)-i3start)  ! set array

c -------------------------------------------------------------------------
c Sub region of main DEM
c For example, to run over the Lake Chad basin
c (about 2.5 million  km2) in north Africa:
      parameter (istart = 2245,iend = 2460) 
      parameter (jstart = 799,jend = 1014)
c Global    
c      parameter (istart = 1,iend = 4320) 
c      parameter (jstart = 1,jend = 2160)

c number of months to run the model (nmons) and the year the
c simulation will start (inity). These are entirely dependent
c on the input data and wishes of the modeler. The example below
c are for the Lake Chad data.

      !parameter (nmons = 708,inity = 1936) !CRU05 data
      parameter (nmons = 12,inity = 1) !Basic SWAM setup

c Define array size of 0.5 deg input data. This must match the sub-region 
c scale defined above.

      parameter (inc = int((iend-(istart-1))/scale2), 
     *           inr = int((jend-(jstart-1))/scale2)) 

c Calculate the fine grid regional array size.

      parameter (incf = iend-(istart-1), inrf = jend-(jstart-1))
c
c  area  = surface area of each grid cell in m**2
c  sfluxin   = flux into each cell at that timestep
c  fluxout   = flux out of a cell at that timestep
c  sfluxout  = mean monthly flux out of a cell (last year of run)
c  voll      = volume of water in river flow reservoir
c  dem       = land elevation
c  larea     = prescribed surface water area (unitless, 1/0)
c  prcpl     = precip at 5' resolution m**3/s
c  evapl     = evap at 5' resolution m**3/s
c  volb      = volume of water in baseflow reservoir
c  rin       = runoff into runoff reservoir m**3/s
c  rout      = flow out of runoff reservoir m**3/s
c  bin       = drainage into baseflow reservoir m**3/s
c  bout      = drainage out of baseflow reservoir m**3/s
c  irrout    = irrigation within the basin m**3/s
c  timed     = residence time of baseflow reservoir
c  timer     = residence time of runoff reservoir
c  effref    = reference effective velocity
c
      double precision tmph,effvel,dist,nspday,prcpl,evapl,
     *     timed,timer,timeg,effref,
     *     outelvp,volchk
      double precision bin,bout,rin,rout,irrout
      double precision drainin(inc,inr,nmons),runin(inc,inr,nmons),
     *     prcpi(inc,inr,nmons),evapi(inc,inr,nmons) 
      double precision area(nrf)
      double precision sfluxin(incf,inrf),
     *     fluxout(incf,inrf),laream(incf,inrf),
     *     volb(incf,inrf),volr(incf,inrf),areat(incf,inrf),
     *     tempdl(incf,inrf),tempdr(incf,inrf),volt(incf,inrf),
     *     basin2(incf,inrf),temp(incf,inrf)
      double precision dem(ncf,nrf),larea(ncf,nrf),basin(ncf,nrf),
     *     outdir(ncf,nrf),sillh(ncf,nrf),
     *     outnewi(ncf,nrf),outnewj(ncf,nrf)
      double precision dvoll(incf,inrf),outelv(incf,inrf),
     *     voll(incf,inrf),fin(incf,inrf),fout(incf,inrf)
c 
      double precision sflux(incf,inrf),elevm(incf,inrf),   !vollm(incf,inrf),
     *     deptm(incf,inrf),  !dvollm(incf,inrf),
     *     lakem(incf,inrf),aream(incf,inrf)
c
      double precision sfluxout(incf,inrf,nmons)   !,laream(ncf,nrf,12)
c
      integer i,j,k,k2,km,kt,i1,j1,i2,j2,ktstep,x,nyrs,
     *        tmpdir,tmpdir2,i3,j3,ix,jx,startyear,ii,jj,iii,jjj
      integer ioff(8),joff(8) 
      integer ndaypm(nmons)

c
c--------------------------------------------------------------
c set values for searching all gridcells around a given gridcell
c 8 possible directions
c
      data ioff /0,1,1,1,0,-1,-1,-1/ ! traditional way
      data joff /-1,-1,0,1,1,1,0,-1/
c--------------------------------------------------------------
c calculate the area of each gridcell phi starts at latitude 
c corresponding to area chosen. Therefore, 90-number of gridcells
c one is starting from top of array. dgper is the number of degrees
c per unit of coarse resolution array.
c
      phi = ((90.-((j2start-1)*dgper))-(2.5/60.))*rad
c
      write(*,*) ncf,nrf
      write(*,*) inc,incf,inr,inrf
      do j = 1,nrf
       dx = cos(phi-((j-1.)*(5./60.)*rad))*dy
       area(j) = dx*dy
      enddo
c
      totarea = 0.
      !chadarea = 0.
      do j = 1,nrf
       do i = 1,ncf
        !write(*,*) dem(i,j)
        if(dem(i,j) .ge. 1.)totarea = totarea + area(j)
        !if(basin(i,j) .eq. 40225.)chadarea = chadarea + area(j)
       enddo
      enddo
      write(*,*)'total area of Grid (km2) = ',totarea/1.e+06,phi 
c      write(*,*)'test cell dem = ',dem(108+2245-1,1+799-1)
c      write(*,*)'test cell ppt = ',prcpi(1,36,12)
c      write(*,*)'total area of Basin (km2) = ',chadarea/1.e+06 
c
c--------------------------------------------------------------
c initialize variables and constants
c nyr is the number of years to run the model, delt is the timestep.
c
      ioo  = 2.5e-03 
c     delt = 1800. !shorter timestep for poleward locations
c      delt = 3600. !good for tropics and mid-latitudes
      delt = 86400. ! One day integration for monthly run
      timed = 13.0e+05 ! Residence time for base flow
      timer = 7200.   ! Residence time for surface 
c     timed = timer     !for IBIS runs, which already calculates residence time
      effref = 0.8      !average of channel and floodplain flow
      ndaypm = (/31,28,31,30,31,30,31,31,30,31,30,31/)
      iday = 0
      imon = 1
      iyear = 1
      itime = 0
      ktstep = 0
      nspday = 24./(delt/3600.) ! Should be 1
      write(*,*) "Steps per day", nspday
      volchk = 0.

c--------------------------------------------------------------
c ***** DEM ADJUSTMENTS GO HERE *****
c THESE ARE CORRECTIONS TO THE DEM SPECIFIC TO THE LAKE CHAD BASIN
c The DEM is in error in many locations therefore, it is often
c necessary to consult accurate charts and correct the data in the
c DEM manually. The corrections below can be a guide.
c
c Correct the Chad basin at Grande Barierre. Currently the north 
c portion of the mid-level lake is too deep. There are elevations 
c which are from 220-260 m when they probably should be around 279 m.
c The south portion of the basin is too shallow. Elevations are
c about 281m and should be perhaps 275.
c
      do 910 j = 291,304
       do 911 i = 454,464
        if(dem(i,j) .le. 282.) then
         dem(i,j) = dem(i,j)-3.
        endif
        if(dem(i,j) .eq. 283.) then
         dem(i,j) = dem(i,j)-2.
        endif
        if(dem(i,j) .gt. 283.) then 
         dem(i,j) = dem(i,j)-1.
        endif
911    continue
910   continue
 
      do 912 j = 278,309
       do 913 i = 430,453
        dem(i,j) = max(dem(i,j), 280.)
913    continue
912   continue
c
c--------------------------------------------------------------
c initialize variables
      sampi=90
      sampj=135
c
      do 205 j = 1,inrf
       do 206 i = 1,incf
        ii = i + (istart-1) ! ii is the global scale index
        jj = j + (jstart-1) ! jj is the global scale index
        outelv(i,j) = dem(ii,jj) ! Local DEM
        elevm(i,j) = 0. ! Water surface elevation
        sflux(i,j) = 0. ! Net water flux / cell
        lakem(i,j) = 0. ! is 'm' the outlet?
        deptm(i,j) = 0.
c       vollm(i,j) = 0.
c       dvollm(i,j) = 0.
        fin(i,j) = 0.
        fout(i,j) = 0.
        voll(i,j) = 0.
        volb(i,j) = 0.
        volr(i,j) = 0.
        volt(i,j) = 0. ! Total (potential) volume for each PWA
        dvoll(i,j) = 0.
        tempdl(i,j) = 0.
        tempdr(i,j) = 0.
        temp(i,j)  = 0.
        fluxout(i,j) = 0.
        sfluxin(i,j) = 0.
        areat(i,j) = 0.
        basin2(i,j) = 0.
        do 207 k = 1,12
         sfluxout(i,j,k) = 0.
207     continue
206    continue
205   continue
c
      do 210 j = jstart,jend
       do 220 i = istart,iend
        if(laket .eq. 0)then
         larea(i,j) = 0.
        else
         larea(i,j) = min(larea(i,j),1.) ! Prescribed lakes
        endif
        if(sillh(i,j) .eq. 0.)then
         outnewi(i,j) = 0.
         outnewj(i,j) = 0.
c        else
c         write(*,*) "sillh",i,j,sillh(i,j)
        endif
        if(outnewi(i,j) .eq. 0.)sillh(i,j) = 0.
 220   continue
 210  continue
c

c--------------------------------------------------------------
c convert input from mm/day to m/s (could be pre-processed?)
c
      do 334 k = 1,nmons
       do 335 j = 1,inr
        do 336 i = 1,inc

         if(runin(i,j,k) .ge. 1.e+20) then
          runin(i,j,k) = 0.
         endif
         if(drainin(i,j,k) .ge. 1.e+20) then 
          drainin(i,j,k) = 0.
         endif
         if(prcpi(i,j,k) .ge. 1.e+20) then 
          prcpi(i,j,k) = 0.
         endif
         if(evapi(i,j,k) .ge. 1.e+20) then 
          evapi(i,j,k) = 0.
         endif
c
         runin(i,j,k)   = max((runin(i,j,k))/0.864e+8,0.)     !IBIS data
         drainin(i,j,k) = max((drainin(i,j,k))/0.864e+8,0.)   !IBIS data
         prcpi(i,j,k)   = max((prcpi(i,j,k))/0.864e+8,0.)     !IBIS data
         evapi(i,j,k)   = max((evapi(i,j,k))/0.864e+8,0.)     !IBIS data
c
c Test to see response of open lake
c
c       evapi(i,j,k)   = 0.
c       runin(i,j,k)   = runin(i,j,k)*100.
c
 336    continue
 335   continue
 334  continue
c

c
c--------------------------------------------------------------
c Calculate the total volume within each potential lake area (volt).
c First set the volume of the outlet location to initial value.
c volt is used to determine when the lake is full and outflow
c will occur to the river downstream of the sill.
c
      do 914 j = jstart,jend
       do 915 i = istart,iend
c
        if(sillh(i,j) .gt. 0.)then
         ! ii and jj give the outlet coordinates
         ii = min(max(outnewi(i,j)-(istart-1),1.),REAL(incf)) 
         jj = min(max(outnewj(i,j)-(jstart-1),1.),REAL(inrf))
         if ((ii .gt. 0).and. (ii .le. iend))then
          if((jj .gt. 0).and. (jj .le. jend))then
           volt(ii,jj) = volt(ii,jj) + 
     *        max(sillh(i,j)-dem(i,j),0.1)*area(j)
           !write(*,*) i,ii,j,jj,volt2(i,j),dem(i,j)
          endif
         endif
        endif
c
915    continue
914   continue
c
      write(*,*)'through lake volume calc'
c--------------------------------------------------------------
c This is used if you would like to specify the location that
c a lake will start filling from. This is useful for closed 
c basin lakes that have multiple sub-basins. If nothing is
c is changed here than the lake will begin to fill at the low
c point nearest to the sill location.
c
c Find the lake kernal location, the location that a lake will
c start to fill from. For now the kernal will be the first pit
c encountered. Later it can be modified to provide a more 
c physically realistic kernal location.
c
c set basin2 as a mask for lake kernal location. Currently am
c setting as the outlet location. Corrections to that will
c be made outside this loop
c
      do j = jstart,jend
       do i = istart,iend
        ii = i-(istart-1)
        jj = j-(jstart-1)
        if(outnewi(i,j) .gt. 0.)then
         if((i .eq. outnewi(i,j)) .and. (j .eq. outnewj(i,j)))then
          basin2(ii,jj) = 1.
         endif
        endif
       enddo
      enddo
c

c--------------------------------------------------------------
c CONVERGENCE (LAKE FILLING) GOES HERE

c--------------------------------------------------------------
c Main loop
      write(*,*)'to main loop'

c loop 130 is the total number of years that the model will be run.
c It is from the startyear to the to the number of years (nyrs) plus
c the number of spin-up years (spin). The spin-up years are there
c to get the model in equilibrium before printing out results.
c
      do 130 iyear = 1,spin+1 ! This can be replaced by convergence
      if (iyear.le.spin) then
       write(*,*)"spin up yr",iyear
      else 
       write(*,*)"simulation",iyear
      endif 
c
c monthly loop
c
      do 131 imon = 1,12
c
c write a few diagnostics
c
      write(*,*)'begin mon ',imon
c
c re-initialize some variables each month.
c
       do j = 1,inrf
        do i = 1,incf
         sfluxout(i,j,imon) = 0.
         laream(i,j) = 0.
        enddo
       enddo
c
c start the daily loop
c
      do 132 iday = 1,ndaypm(imon)
c
c start first spatial loop: linear reservoir model
c model is:
c dWV/dt = rin*(1-WA) + (P-E)*WA + sum(Fin) - Fout
c
c Fout = max(((WEl - WEd)*A),0) * u/D
c u = 0.003 m s-1; D is distance between centers
c
      do 120 j = jstart+1,jend-1
      do 110 i = istart+1,iend-1

      if(basin(i,j) .gt. 0.)then ! are we in a basin?
c
c----------------------------------------------------
c Initialize hourly climate input variables.
c
       bin     = 0. !sub-surface runoff
       rin     = 0. !surface runoff
       bout    = 0. !sub-surface flow to river
       rout    = 0. !surface flow to river
       irrout  = 0. !irrigation rate (prescribed below)
       prcpl   = 0. !precipitation rate
       evapl   = 0. !evaporation rate
c 
c----------------------------------------------------
c Calculate the daily runoff, precip, and evap input for a given cell
c from the monthly means.
c Derive 1/2 degree cell which corresponds to current 5' location
c Could be replaced with daily function from splash
c
       i2 = int((i-1)/scale2)-int((((istart-1)/scale2)-1))!longitude at 1/2 deg
       j2 = int((j-1)/scale2)-int((((jstart-1)/scale2)-1))!latitude at 1/2 deg
c
       if(iday .le. 15) then ! First part of month
c
        km = imon-1 ! Interpolate from preceding month
c
        if(imon .eq. 1) km = 12 ! Unless in January
c
        rin   = max((runin(i2,j2,km)*area(j)) + ((iday+15)/30.)*
     *     ((runin(i2,j2,imon)*area(j)) -
     *     (runin(i2,j2,km)*area(j))),0.)
        bin   = max((drainin(i2,j2,km)*area(j)) + ((iday+15)/30.)*
     *     ((drainin(i2,j2,imon)*area(j)) -
     *     (drainin(i2,j2,km)*area(j))),0.)
        prcpl =  (prcpi(i2,j2,km)*area(j)) + ((iday+15)/30.)*
     *     ((prcpi(i2,j2,imon)*area(j)) -
     *     (prcpi(i2,j2,km)*area(j)))
        evapl = (evapi(i2,j2,km)*area(j)) + ((iday+15)/30.)*
     *     ((evapi(i2,j2,imon)*area(j)) -
     *     (evapi(i2,j2,km)*area(j)))
c
       else ! Second part of month
c
        km = imon+1
c
        if(imon .eq. nmons) km = 697
c
        rin   = max((runin(i2,j2,imon)*area(j)) + ((iday-15)/30.)*
     *     ((runin(i2,j2,km)*area(j)) -
     *     (runin(i2,j2,imon)*area(j))),0.)
        bin   = max((drainin(i2,j2,imon)*area(j)) + ((iday-15)/30.)*
     *     ((drainin(i2,j2,km)*area(j)) -
     *     (drainin(i2,j2,imon)*area(j))),0.)
        prcpl   = (prcpi(i2,j2,imon)*area(j)) + ((iday-15)/30.)*
     *     ((prcpi(i2,j2,km)*area(j)) -
     *     (prcpi(i2,j2,imon)*area(j)))
        evapl   = (evapi(i2,j2,imon)*area(j)) + ((iday-15)/30.)*
     *     ((evapi(i2,j2,km)*area(j)) -
     *     (evapi(i2,j2,imon)*area(j)))
c
       endif ! End daily climate loop
c       if (runin(i2,j2,km).gt.0.) then
c        !write(*,*) rin,prcpl,evapl
c        write(*,*) runin(i2,j2,km),rin
c       endif
       if ((ii.eq.sampi).and.(jj.eq.sampj)) then
               write(*,*) icmon,iday
               write(*,*) prcpi(i2,j2,icmon),prcpi(i2,j2,km),prcpl
               write(*,*) evapi(i2,j2,icmon),evapi(i2,j2,km),evapl
       endif
 
c
c --------------------------------------------------------------------
c Irrigation adjustments go here
c

c --------------------------------------------------------------------
c set i and j locations at smaller grid over which model is being run (ii,jj)
c
       ii = i - (istart-1)
       jj = j - (jstart-1)
c
c calculate volume in runoff reservoir for land or lake
c
       rout = volr(ii,jj)/timer
       volr(ii,jj) = max(volr(ii,jj) + (rin-rout)*delt,0.)
       !rout = rin * delt
       !volr(ii,jj) = rout
       
c       if (rin.gt.0.) then
c        write(*,*) rin,rout,volr(ii,jj)
c       endif
c
c calculate volume in baseflow reservoir, land only
c
c      bout = 0.75*volb(ii,jj)/timed + 0.25*volb(ii,jj)/timeg
       bout = volb(ii,jj)/timed
       volb(ii,jj) = max(volb(ii,jj) + (bin-bout)*delt,0.)
       !bout = bin * delt
       !volb(ii,jj) = bout
c
c----------------------------------------------------
c calculate volume in transport reservoir
c The idea is to calculate the volume
c of each non-grid lake cell as usual. For the lakes
c cells add any water volume to the volume of the
c sill grid cell not to the local cell and don't 
c include the (fluxin-fluxout). It will all be rectified
c in the next loop. It is done this way to make sure that the
c lake calculation conserves mass and is stable.
c
c calculate the location of the outlet for this lake (i2,j2).
c
       i2 = min(max(int(outnewi(i,j)-(istart-1)),0),incf)
       j2 = min(max(int(outnewj(i,j)-(jstart-1)),0),inrf)
c
c lake water balance for this timestep
c
       temp(ii,jj) = (rout+bout)*(1.-larea(i,j))
c
c land water balance for this timestep.
c
       tempdr(ii,jj) = ((rout+
     *            bout)*(1.-larea(i,j))-irrout
     *           + (sfluxin(ii,jj) - fluxout(ii,jj)))*delt
c
c subtract any evaporation from the lake from the outlet
c location. The outlet is the accountant for the entire lake.
c
      if(larea(i,j) .gt. 0.)then
       !write(*,*) larea(i,j)
       tempdl(i2,j2) = tempdl(i2,j2) 
     *          + ((prcpl-evapl)*larea(i,j))*delt
      endif
c
      else                 !basin .ne. requested number
       voll(ii,jj) = 0.
      endif ! End if basin check
c
 110    continue
 120   continue ! End of linear resevoir model
c
c----------------------------------------------------------------
c Calculate the change in volume (dvoll) as the sum of the P-E (tempdl)
c and the river flow (tempdr). Set the minimum value of dvoll.
c Add dvoll to existing reservoir volume (voll). Set variables
c to 0. for next timestep.
c
       do 140 j = jstart+1,jend-1
        do 141 i = istart+1,iend-1
         ii = i - (istart-1)
         jj = j - (jstart-1)
         dvoll(ii,jj)  = tempdr(ii,jj) + tempdl(ii,jj)
         !if ((ii.eq.sampi).and.(jj.eq.sampj)) then
         !       write(*,*) iday,i,j,dvoll(ii,jj)
         !endif
         if(abs(dvoll(ii,jj))/delt/area(j) .lt. dveps) then
          dvoll(ii,jj) = 0.
         endif
         voll(ii,jj)  = max(voll(ii,jj) + dvoll(ii,jj),0.)
         tempdl(ii,jj) = 0.
         tempdr(ii,jj) = 0.
         !if ((ii.eq.sampi).and.(jj.eq.sampj)) then
         !       write(*,*) dvoll(ii,jj),voll(ii,jj)
         !endif
141     continue
140    continue
       
c
c----------------------------------------------------------------
c Distribute the dvoll of this timestep roughly into the existing lake area
c It will be evened out in loop 121
c
       do 122 j = jstart+1,jend-1
        do 112 i = istart+1,iend-1
        if(basin(i,j) .gt. 0.)then
c
         ii = i - (istart-1)
         jj = j - (jstart-1)
         i2 = min(max(outnewi(i,j)-(istart-1),0.),REAL(incf))
         j2 = min(max(outnewj(i,j)-(jstart-1),0.),REAL(inrf))
c
         if(((i2 .gt. 0).and. (j2 .gt. 0)) 
     *        .and. (laket .eq. 0))then
c
c          if ((ii.eq.sampi).and.(jj.eq.sampj).and.iday.eq.15) then
                !write(*,*) i,j,ii,jj,i2,j2,laket
                !write(*,*) outnewi(i,j),outnewj(i,j),basin(i,j)
                !write(*,*) voll(i2,j2),"gt",volt(i2,j2)
c                if (voll(i2,j2) .ge. volt(i2,j2))then
c                        write(*,*) "YES"
c                else
c                        write(*,*) "NO"
c                endif
c          endif
c if volume in basin > lake volume larea = 1. everywhere
c outelv = sill height
c
          if(voll(i2,j2) .ge. volt(i2,j2))then
           outelv(ii,jj) = max(outelv(ii,jj) + larea(i,j)*
     *            dvoll(i2,j2)/areat(i2,j2),dem(i,j))  !0.0 if larea = 0.
c          outelv(ii,jj) = sillh(ii,jj)  !simpler way of handling it
           larea(i,j) = max(min(outelv(ii,jj)-dem(i,j),1.),0.)
c
c if there is no volume then larea = 0.
c
          elseif(voll(i2,j2) .eq. 0.)then 
           larea(i,j)  = 0.
           outelv(ii,jj) = dem(i,j)
           areat(i2,j2) = 0.   !probably not necessary
c
          elseif((voll(i2,j2).gt.0.).and.
     *           (voll(i2,j2).lt.volt(i2,j2)))then
c
c if some lake already exists in closed basin distribute dvoll
c evenly to those existing cells 
c
           if(areat(i2,j2) .gt. 0.)then
c
c set outelv if larea > 0. Add depth of water if positive
c or negative. 
c
            outelv(ii,jj) = max(outelv(ii,jj) + larea(i,j)*
     *            dvoll(i2,j2)/areat(i2,j2),dem(i,j))  !0.0 if larea = 0.
c
c If no existing lake; set larea of outlet location = 1.
c Ideally would like to choose a kernal location in a 
c realistic location within a lake. This would be stored
c in array basin2.
c
           else   !if(areat(i2,j2) .eq. 0.))then 
            if(basin2(ii,jj) .eq. 1.)then
             outelv(ii,jj) = max(outelv(ii,jj) +
     *           dvoll(ii,jj)/area(j2+(jstart-1)),dem(i,j))
             larea(i,j) = max(min(outelv(ii,jj)-dem(i,j),1.),0.)
             areat(i2,j2) = max(area(j)*larea(i,j),0.)
c
            endif !basin2 if loop
c
           endif ! areat if loop
c
          endif !voll if loop
         endif ! i2/j2/laket loop
c
        else                 !basin .ne. requested number
         voll(ii,jj) = 0.
        endif
c
        fluxout(ii,jj) = 0.
        sfluxin(ii,jj) = 0.
c
 112    continue
 122   continue
c
c--------------------------------------------------------------
c In this loop calculate the flux out of the cell, to which direction,
c and the sum of the fluxes into cells. Also calculate the
c water depth and lake area (larea) for each cell. 
c
      do 121 j = jstart+1,jend-1
       do 111 i = istart+1,iend-1
c
       if(basin(i,j) .gt. 0.)then ! Check if we are in a basin
c
       ii = i - (istart-1) ! Fine grid index
       jj = j - (jstart-1)
       iii = min(max(outnewi(i,j)-(istart-1),1.),REAL(incf)) ! Outlet
       jjj = min(max(outnewj(i,j)-(jstart-1),1.),REAL(inrf))
c
c Find flow direction 
c need to replace this with GRASS derived directions
c
       tmpdir = outdir(i,j) 
c
       j2 = j + joff(tmpdir)
       i2 = i + ioff(tmpdir)
       dx = (area(j)/dy+area(j2)/dy)/2.
       dist = sqrt((dx*dx*abs(i2-i)*abs(i2-i))
     *        + (dy*dy*abs(j2-j)*abs(j2-j)))
       dist = max(dx,dist)
c
c set effective velocity of each cell. It is dependent on the
c lake volume of the cell. For large lakes the value is quite
c slow, for small lakes it approaches the reference value,
c effref. For non-lake points it is a function of the
c gradient as in Miller et al. 1994. This is fairly rough
c and could be improved with considerations of sinuosity or
c stream order for example.
c
c Invoke the top half of this if statement if you want the volume
c of the lake to impact the river discharge velocity. It is unique
c to each lake and should be studied before use
c
c      if((larea(i,j) .gt. 0.) .and.
c    *    (voll(ii,jj) .gt. 0.))then
c       vollocal = 1.*area(j)
c       volref = sqrt(vollocal/volt(iii,jjj))
c       effvel = min(effref,0.1*effref*volref)
c       effvel = min(effref,0.08*effref*volref)
c       effvel = max(effvel,1.0e-02)
c      else                         !non-lake
        ic = max(dem(i,j)-dem(i2,j2),1.)/dist
        io = ioo/(dx/dy)   !scale reference gradient to latitude
        effvel = effref*sqrt(ic/io)
        effvel = min(effvel,3.0)
c      endif
c
c calculate fluxout of each cell and send it as fluxin to
c either the cell downstream if it is not a lake downstream
c or to the outlet location. The fluxout of the cell which
c corresponds to the sill is calculated for only that water 
c volume in excess of the volume required to fill the lake.
c
c       effvel = 0.5 !alternatively could set velocity to a  constant
c
        fluxout(ii,jj) = max((voll(ii,jj)-volt(ii,jj))*
     *                  (effvel/dist),0.)
        fluxout(ii,jj) = max(min(fluxout(ii,jj),
     *           sfluxin(ii,jj) + temp(ii,jj) +
     *           ((voll(ii,jj)-volt(ii,jj))/(delt*2.))),0.)
c
c Truncate fluxout if too small for computation. 
c
        if(fluxout(ii,jj)/area(j) .lt. dveps) then 
         fluxout(ii,jj) = 0. 
        endif
        i3 = i2-(istart-1)
        j3 = j2-(jstart-1)
        if((i3.gt.0).and.(i3.le.incf).and.
     *      (j3.gt.0).and.(j3.le.inrf)) then
         sfluxin(i3,j3) =
     *    sfluxin(i3,j3) + fluxout(ii,jj)
        endif
c
c--------------------------------------------------------------
c Adjust height of water column in each cell using the cellular
c automata. Distribute water height within a lake basin only
c This flattens the lake surface so that there are no hills or
c valleys due to differences in the local water budget.
c
       if((laket .eq. 0) .and. ((iii .ne. 0).and. (jjj.ne.0)))then
c
        if((outelv(ii,jj) .gt. dem(i,j)) .and.
     *     (voll(iii,jjj) .lt.
     *      volt(iii,jjj)))then
c
         if((voll(iii,jjj) .ge. 0.).and.
     *      (sillh(i,j) .gt. 0.))then
c
          tmpdir2 = 0
          tmph = outelv(ii,jj)
          do 501 k2 = 1,8
           j3 = (j + joff(k2))-(jstart-1)
           i3 = (i + ioff(k2))-(istart-1)
           if((outelv(i3,j3) .lt. outelv(ii,jj)) .and.
     *        (sillh(i3+(istart-1),j3+(jstart-1)) .eq. sillh(i,j)))then
            if(outelv(i3,j3) .lt. tmph)then
             tmpdir2 = k2
             tmph = outelv(i3,j3)
            endif
           endif
501       continue
c
          if(tmpdir2 .ne. 0.)then
c
           j3 = (j + joff(tmpdir2))-(jstart-1)
           i3 = (i + ioff(tmpdir2))-(istart-1)
c
           gridif = min(max((outelv(ii,jj) -
     *               outelv(i3,j3)),0.),
     *               max((outelv(ii,jj) - dem(i,j)),0.))
           if(abs(gridif) .lt. grideps) then
            gridif = 0.
           endif
           gridif = gridif*area(j)*0.5  !move only 0.5 of difference
c
           outelv(i3,j3) = outelv(i3,j3) +
     *           gridif/area(j3+(jstart-1))
           outelv(ii,jj) = outelv(ii,jj) -
     *           gridif/area(j)
c
          endif
c
         else
          larea(i,j) = 0.
          outelv(ii,jj) = dem(i,j)
          areat(iii,jjj) = 0.
c 
         endif                !voll(i2,j2) .ge. 0.
c
        endif                !outelv .gt. dem
c
c Set lake area = either; 1 if depth is greater than 1 meter, to
c a % of the lake cell equal to the depth of water, or to 0. if
c depth = 0.  depth = outelv-grid
c Adjust the total lake area (areat) to represent current larea.
c Do this by first subtracting current larea from total then
c adding new larea to total. If larea = 0. area added = 0.
c
        if((outnewi(i,j) .gt. 0.) .and.
     *     (voll(iii,jjj) .lt. volt(iii,jjj)))then
c
         if((voll(iii,jjj) .gt. 0.) .and.
     *       (sillh(i,j) .gt. 0.))then
c
          ix = min(max(int(outnewi(i,j)-(istart-1)),1),incf)
          jx = min(max(int(outnewj(i,j)-(jstart-1)),1),inrf)
          areat(ix,jx) = max(areat(ix,jx) - area(j)*
     *         larea(i,j),0.)
          larea(i,j) = max(min(outelv(ii,jj)-
     *         dem(i,j),1.),0.)
          areat(ix,jx) = max(areat(ix,jx) + area(j)*
     *         larea(i,j),0.) !previous
         else
          larea(i,j) = 0.
          outelv(ii,jj) = dem(i,j)
          areat(iii,jjj) = 0.
         endif
        endif
c
        laream(ii,jj) = (laream(ii,jj)+((outelv(ii,jj)-dem(i,j))
     *                        /(ndaypm(imon)*nspday)))
        if(outelv(ii,jj)-dem(i,j) .lt. 0.1)then
         laream(ii,jj) = 0.
        endif
c
       endif        !laket = 0
c
       endif              !endif for limiting calculation to a basin
c
111    continue
121   continue          !end fluxout loop
c

132   continue
c
c This is a check to make sure that the volume of water stored in
c the lakes (volchk) is equal to the amount of water that is in
c voll. This tells me mass conservation is taking place. It is
c possible that it may not be conserved if I make a mistake in
c the way water is spread acros the surface.
c
      do 511 j = jstart+1,jend-1
       do 521 i = istart+1,iend-1
        ii = i - (istart-1)
        jj = j - (jstart-1)
        volchk = volchk + ((outelv(ii,jj)-dem(i,j))
     *           *area(j))*larea(i,j)
c        if (iyear.eq.spin+1) then
c                write(*,*) outelv(ii,jj),dem(i,j)
c        endif
521    continue
511   continue
c
c write thee out if you are interested in checking the mass conservation
c
c      write(308,*)'month = ',imon
c      write(308,*)'voll at step ',ktstep,'  = ',
c    *              voll(2343-(istart-1),969-(jstart-1))
c      write(308,*)'volchk at step ',ktstep,'  = ',volchk
c      write(*,*)'volchk = ',volchk
c      write(308,*)
       volchk = 0.
c
c Write monthly flux to output file. This can be modified to write out other
c variables if one wishes.
c
c      write(35,*)iyear,iwmon,' ','Chad areat = ',
c    *            areat(2343-(istart-1),969-(jstart-1))/1.e+6
c
       evapm = 0.    !variable to sum evap over entire lake basin
       evapt = 0.    !variable to sum evap over entire lake basin
       lareat = 0.    !variable to sum evap over entire lake basin
 
        do j = jstart,jend
         do i = istart,iend
         !if (sillh(i,j).eq.305.) then
         !        write(*,*) sillh(i,j),basin(i,j)
         !endif
          ii = i - (istart-1)
          jj = j - (jstart-1)
          elevm(ii,jj) = 0.
          lakem(ii,jj) = 0.
          deptm(ii,jj) = 0.
          sflux(ii,jj) = 0.
c         vollm(ii,jj) = 0.
          aream(ii,jj) = 0.
         if(laream(ii,jj) .gt. 0.)then
          elevm(ii,jj) = laream(ii,jj)+dem(i,j)
          lakem(ii,jj) = larea(i,j)
          deptm(ii,jj) = laream(ii,jj)
          evapm = evapm + evapl*larea(i,j)
          lareat = lareat + larea(i,j)
         endif
         evapt = evapt + evapl*area(j)
         sflux(ii,jj) = sfluxout(ii,jj,imon)
c        vollm(ii,jj) = vollout(ii,jj,imon)/1.e+09
c        dvollm(ii,jj) = dvoll(ii,jj)
         aream(ii,jj) = max(areat(ii,jj)/1.e+06,0.)
c        aream(ii,jj) = volt(ii,jj)
c        deptm(ii,jj) = basin2(ii,jj)
         laream(ii,jj) = 0.
        enddo
       enddo
c
       write(*,*)'evapm  = ',evapm     !sum of basin lake evaporation
       write(*,*)'evapt  = ',evapt     !sum of basin lake evaporation
       write(*,*)'lareat  = ',lareat    !total lake area
c
131   continue
130   continue
      write(*,*)'************'
      write(*,*)'Run Complete'
      write(*,*)'************'
c
      end
