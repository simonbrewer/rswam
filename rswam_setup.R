## RSWAM: Uses rewritten code
## Test case (runs Chad example)

require(raster)
# Geomorphology data
outnewi.r = raster("../geomorphNetCDFin/outnewi.nc")
outnewj.r = raster("../geomorphNetCDFin/outnewj.nc")
basin.r = raster("../geomorphNetCDFin/basin.nc")
dem.r = raster("../geomorphNetCDFin/tbase.nc")
rivdir.r = raster("../geomorphNetCDFin/rivdir.nc")
mflac.r = raster("../geomorphNetCDFin/mflac.nc")

outnewi = as.matrix(outnewi.r)
outnewj = as.matrix(outnewj.r)
basin = as.matrix(basin.r)
dem = as.matrix(dem.r)
rivdir = as.matrix(rivdir.r)
mflac = as.matrix(mflac.r)

## Commented just to speed things up
## Load climate data
nyrs = 59

## Precip
prec.r = stack("../geomorphNetCDFin/rain.nc")
dem2.r = crop(dem.r, extent(prec.r))
prec.r.m = stackApply(prec.r, rep(seq(1,12), nyrs), fun=mean)
#prec.r.m = resample(prec.r.m, dem2.r, method="ngb")
prec.r.y = stackApply(prec.r.m, rep(1, 12), fun=mean)

## Evap
evap.r = stack("../geomorphNetCDFin/evap.chadf.nc")
evap.r.m = stackApply(evap.r, rep(seq(1,12), nyrs), fun=mean)
#evap.r.m = resample(evap.r.m, dem2.r, method="ngb")
evap.r.y = stackApply(evap.r.m, rep(1, 12), fun=mean)

## Runoff
runoff.r = stack("../geomorphNetCDFin/srunoff.ex4.nc")
runoff.r.m = stackApply(runoff.r, rep(seq(1,12), nyrs), fun=mean)
#runoff.r.m = resample(runoff.r.m, dem2.r, method="ngb")
runoff.r.y = stackApply(runoff.r.m, rep(1, 12), fun=mean)

## Runoff
drain.r = stack("../geomorphNetCDFin/drainage.ex4.nc")
drain.r.m = stackApply(drain.r, rep(seq(1,12), nyrs), fun=mean)
#drain.r.m = resample(drain.r.m, dem2.r, method="ngb")
drain.r.y = stackApply(drain.r.m, rep(1, 12), fun=mean)

save(outnewi.r, outnewj.r, 
     basin.r, dem.r, rivdir.r, mflac.r,
     prec.r.m, prec.r.y, 
     evap.r.m, evap.r.y, 
     runoff.r.m, runoff.r.y, 
     drain.r.m, drain.r.y, 
     file="chadExample.RData")
