## RSWAM1: uses original fortran code
## Test case (runs Chad example)

require(raster)
dyn.load("./src/rswam1_testing.so")
rswam1 <- function(nyrs, startyear, converg = 0, laket = 0, spin = 1,
                  normal = 1, leap = 1, irrig = 1,
                  outnewi, outnewj, basin, dem, rivdir, mflac,
                  prcpi, evapi, runin, drainin,
                  gridxf, gridyf) {
  
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
                   rivdir = as.double(rivdir), sillh = as.double(mflac),
                   prcpi = as.double(prec), evapi = as.double(evap),
                   runin = as.double(runoff), drainin = as.double(drain),
                   outelv = double(gridxf*gridyf), lakem = double(gridxf*gridyf), 
                   lakevolm = double(12*(spin+nyrs)), lakevola = double(spin+nyrs))
  return(simcf)
  
}

## Parameters
nyrs = 1 # of years to run the model after January 1937
startyear = 1 # year to start run (from January 1937) 
converg = 1 # 0 if set convergence helper, 1 if not
laket = 0 # 0 if predict lakes, 1 if parameterized with Cogley obs, 2 if with previously simulated
spin = 0  # number of spin-up years to run before reading transient data
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
basinf.r = crop(basin.r, extent(prec.r.y))
mflacf.r = crop(mflac.r, extent(prec.r.y))

## Grid dimensions
gridxf = 216
gridyf = 216

hydro.out = rswam1(nyrs, startyear, converg, laket, spin,
                  normal, leap, irrig,
                  outnewi, outnewj, basin, 
                  dem, rivdir, mflac,
                  prcpi=prec, evapi=evap, 
                  runin=runoff, drainin=drain,
                  gridxf, gridyf)

stop()
save(hydro.out, file=paste0("rswam1_out_",spin,".RData"))
outelv.r = setValues(demf.r, hydro.out$outelv)
plot(outelv.r)
plot(outelv.r-demf.r)

lakem.r = setValues(demf.r, hydro.out$lakem)
plot(lakem.r)
