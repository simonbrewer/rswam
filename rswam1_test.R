## RSWAM1: uses original fortran code
## Test case (runs Chad example)

require(raster)
dyn.load("./src/rswam1.so")
rswam1 <- function(nyrs, startyear, converg = 1, laket = 0, spin = 1,
                  normal = 1, leap = 1, irrig = 1,
                  outnewi, outnewj, basin, dem, rivdir, mflac,
                  prcpi, evapi, runin, drainin) {
  
  simcf = .Fortran("rswam1",
                   nyrs = as.integer(nyrs),
                   startyear = as.integer(startyear), 
                   converg = as.integer(converg), 
                   laket = as.integer(laket), 
                   spin = as.integer(spin),
                   normal = as.integer(normal), 
                   leap = as.integer(leap), 
                   irrig = as.integer(irrig),
                   outnewi = as.double(outnewi), outnewj = as.double(outnewj),
                   basin = as.double(basin), dem = as.double(dem),
                   rivdir = as.double(rivdir), mflac = as.double(mflac),
                   prcpi = as.double(prec), evapi = as.double(evap),
                   runin = as.double(runoff), drainin = as.double(drain))
  return(simcf)
  
}

## Parameters
nyrs = 1 # of years to run the model after January 1937
startyear = 1 # year to start run (from January 1937) 
converg = 1 # 0 if set convergence helper, 1 if not
laket = 0 # 0 if predict lakes, 1 if parameterized with Cogley obs, 2 if with previously simulated
spin = 10  # number of spin-up years to run before reading transient data
normal = 1 # 0 if use normalization, 1 if not
leap = 2 # of years to the first leap year
irrig = 1 # 0 if use irrigation, function 1 if not

load("./data/chadExample.RData")

# Geomorphology data
outnewi = t(as.matrix(outnewi.r))
outnewj = t(as.matrix(outnewj.r))
basin = t(as.matrix(basin.r))
dem = t(as.matrix(dem.r))
rivdir = t(as.matrix(rivdir.r))
mflac = t(as.matrix(mflac.r))

## Commented just to speed things up
# Climate data
prec = aperm(as.array(prec.r.m), c(2,1,3))
evap = aperm(as.array(evap.r.m), c(2,1,3))
drain = aperm(as.array(drain.r.m), c(2,1,3))
runoff = aperm(as.array(runoff.r.m), c(2,1,3))

demf.r = crop(dem.r, extent(prec.r.y))

## Grid dimensions
# gridx = dim(dem)[1]
# gridy = dim(dem)[2]

hydro.out = rswam1(nyrs, startyear, converg, laket, spin,
                  normal, leap, irrig,
                  outnewi, outnewj, basin, 
                  dem, rivdir, mflac,
                  prcpi=prec, evapi=evap, 
                  runin=runoff, drainin=drain)
