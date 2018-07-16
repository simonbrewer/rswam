## RSWAM1: uses original fortran code
## Test case (runs Chad example)

require(raster)
# dyn.load("./src/rswam.so")
# rswam <- function(nyrs, startyear, converg = 1, laket = 0, spin = 1,
#                   normal = 1, leap = 1, irrig = 1,
#                   outnewi, outnewj, basin, dem, rivdir, mflac,
#                   prcpi, evapi, runin, drainin) {
#   
#   simcf = .Fortran("rswam",
#                    nyrs = as.integer(nyrs),
#                    startyear = as.integer(startyear), 
#                    converg = as.integer(converg), 
#                    laket = as.integer(laket), 
#                    spin = as.integer(spin),
#                    normal = as.integer(normal), 
#                    leap = as.integer(leap), 
#                    irrig = as.integer(irrig),
#                    outnewi = as.double(outnewi), outnewj = as.double(outnewj),
#                    basin = as.double(basin), dem = as.double(dem),
#                    rivdir = as.double(rivdir), mflac = as.double(mflac),
#                    prcpi = as.double(prec), evapi = as.double(evap),
#                    runin = as.double(runoff), drainin = as.double(drain))
#   return(simcf)
#   
# }

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

## Note that all matrices need to be transposed
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
gridxc = dim(dem.r)[2]
gridyc = dim(dem.r)[1]
gridxf = dim(demf.r)[2]
gridyf = dim(demf.r)[1]

dyn.load("./src/rswam.so")
simcf = .Fortran("rswam", 
                 dem = as.double(dem),
                 basin = as.double(basin),
                 sillh = as.double(mflac),
                 outnewi = as.double(outnewi),
                 outnewj = as.double(outnewj),
                 prcpi = as.double(prec), evapi = as.double(evap),
                 runin = as.double(runoff), drainin = as.double(drain),
                 laket = as.integer(laket),
                 spin = as.integer(spin),
                 volt2 = double(gridxf*gridyf),
                 voll = double(gridxf*gridyf))

  stop()

## Volume check
volt.r = setValues(demf.r, simcf$volt2)
plot(crop(volt.r, extent(mflac.r, 799, 1014, 2245, 2460)))
plot(crop(volt.r, extent(c(-180,0, 85,90))))
## Coordinate checks
dem.r[1,1]
dem.r[799,2245]
dem.r[1014,2460]
x.r = crop(dem.r, extent(mflac.r, 799, 1014, 2245, 2460))
x = t(as.matrix(x.r))
x[1,1]
x[216,216]

##
x.s = stack(x.r, x.r*-1)
x.a = aperm(as.array(x.s), c(2,1,3))


p.r = raster(prec.r.m,12)
plot(p.r)
prec.r.m[36,1]
p = aperm(as.array(prec.r.m), c(2,1,3))
p[1,36,12]

## Sill checks
x.r = crop(mflac.r, extent(mflac.r, 799, 1014, 2245, 2460))
plot(x.r==305)
image(t(as.matrix(x.r))==305)
