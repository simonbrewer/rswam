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

      subroutine rswam( dem, basin, sillh, outnewi, outnewj,
     >                  prcpi,
     >                  laket, spin,
     >                  volt2)
              
      ! Input variables
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
      double precision volt2(incf,inrf)
      double precision dem(ncf,nrf),larea(ncf,nrf),basin(ncf,nrf),
     *     outdir(ncf,nrf),sillh(ncf,nrf),
     *     outnewi(ncf,nrf),outnewj(ncf,nrf)
      double precision dvoll(incf,inrf),outelv(incf,inrf),
     *     voll(incf,inrf)
c 
      double precision sflux(incf,inrf),elevm(incf,inrf),   !vollm(incf,inrf),
     *     deptm(incf,inrf),  !dvollm(incf,inrf),
     *     lakem(incf,inrf),aream(incf,inrf)
c
      double precision sfluxout(incf,inrf,nmons)   !,laream(ncf,nrf,12)
c
      integer i,j,k,k2,km,kt,i1,j1,i2,j2,ktstep,
     *        x,nyrs,tmpdir,i3,j3,ix,jx,startyear,ii,jj,iii,jjj
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
      write(*,*)'total area of Grid (km2) = ',totarea/1.e+06 
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
c
c--------------------------------------------------------------
c initialize variables
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
        voll(i,j) = 0.
        volb(i,j) = 0.
        volr(i,j) = 0.
        volt(i,j) = 0. ! Total (potential) volume for each PWA
        volt2(i,j) = 0. ! Total (potential) volume for each cell
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
      write(*,*) "laket",laket
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
c      do 334 k = 1,nmons
c       do 335 j = 1,inr
c        do 336 i = 1,inc
c
c         if(runin(i,j,k) .ge. 1.e+20) then
c          runin(i,j,k) = 0.
c         endif
c         if(drainin(i,j,k) .ge. 1.e+20) then 
c          drainin(i,j,k) = 0.
c         endif
c         if(prcpi(i,j,k) .ge. 1.e+20) then 
c          prcpi(i,j,k) = 0.
c         endif
c         if(evapi(i,j,k) .ge. 1.e+20) then 
c          evapi(i,j,k) = 0.
c         endif
c
c         runin(i,j,k)   = max((runin(i,j,k))/0.864e+8,0.)     !IBIS data
c         drainin(i,j,k) = max((drainin(i,j,k))/0.864e+8,0.)   !IBIS data
c         prcpi(i,j,k)   = max((prcpi(i,j,k))/0.864e+8,0.)     !IBIS data
c         evapi(i,j,k)   = max((evapi(i,j,k))/0.864e+8,0.)     !IBIS data
c
c Test to see response of open lake
c
c       evapi(i,j,k)   = 0.
c       runin(i,j,k)   = runin(i,j,k)*100.
c
c 336    continue
c 335   continue
c 334  continue
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
           volt2((i-istart+1),(j-jstart+1))=volt(ii,jj)
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
c re-initialize some variables each month.
c
      do j = 1,inrf
       do i = 1,incf
        sfluxout(i,j,imon) = 0.
        laream(i,j) = 0.
       enddo
      enddo
c
c write a few diagnostics
c
      write(*,*)'begin mon ',imon
c
c start the daily loop
c
      do 132 iday = 1,ndaypm(imon)
c
c start first spatial loop.
c note that the subdaily loop has been removed
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
        if(icmon .eq. 1) km = 12 ! Unless in January
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
        if(icmon .eq. nmons) km = 697
c
        rin   = max((runin(i2,j2,imon)*area(j)) + ((iday-15)/30.)*
     *     ((runin(i2,j2,km)*area(j)) -
     *     (runin(i2,j2,icmon)*area(j))),0.)
        bin   = max((drainin(i2,j2,imon)*area(j)) + ((iday-15)/30.)*
     *     ((drainin(i2,j2,km)*area(j)) -
     *     (drainin(i2,j2,icmon)*area(j))),0.)
        prcpl   = (prcpi(i2,j2,imon)*area(j)) + ((iday-15)/30.)*
     *     ((prcpi(i2,j2,km)*area(j)) -
     *     (prcpi(i2,j2,icmon)*area(j)))
        evapl   = (evapi(i2,j2,imon)*area(j)) + ((iday-15)/30.)*
     *     ((evapi(i2,j2,km)*area(j)) -
     *     (evapi(i2,j2,icmon)*area(j)))
c
       endif ! End daily climate loop
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
c
c calculate volume in baseflow reservoir, land only
c
c      bout = 0.75*volb(ii,jj)/timed + 0.25*volb(ii,jj)/timeg
       bout = volb(ii,jj)/timed
       volb(ii,jj) = max(volb(ii,jj) + (bin-bout)*delt,0.)
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

132   continue
131   continue
130   continue
      write(*,*)'************'
      write(*,*)'Run Complete'
      write(*,*)'************'
c
      end
