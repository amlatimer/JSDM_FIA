######## Make layers of projected climate change for Sierras
# Using CMIP5 downscaled by worldclim

# Method: Calculate long-term annual means for 1970-1999 and for 30-yr period centered at 2050, subtract historical normals to get change. Add change to the PRISM long-term normals to get future projection. 
# Use 2.5 minute (~4.5 km) data which is about the resolution of PRISM

# Start with single model x projection: GFDL RCP 4.5 

setwd("/users/latimer/googledriveucd/samsi/fia")

library(raster)

d_site = read.csv("./Derek/FIA_non-WB_predictors.csv")
plotlocs=SpatialPoints(d_site[,c("LON", "LAT")])

### Precip
# future
r = raster("/Volumes/LaCie/Worldclim/gf45bi50/gf45bi5012.tif", pattern="grd", full.names=TRUE) 
#e = extent(c(-122, -117.5, 35, 41)) # sierras
e = extent(c(-125, -116, 32, 42)) # california
ppt.future = crop(r, e)
plot(ppt.future)

# historical
r = raster("./climatedata/bio_2-5m_bil/bio12.bil", pattern="grd", full.names=TRUE) 
#e = extent(c(-122, -117.5, 35, 41))
e = extent(c(-125, -116, 32, 42))
ppt.hist = crop(r, e)
par(mar=c(6,6,1,1))
plot(ppt.hist, xlab="LON", ylab="LAT", cex.axis=1.2, cex.lab=1.2)

# change
plot(ppt.future - ppt.hist)
points(plotlocs)

### Temperature
# future
r = raster("/Volumes/LaCie/Worldclim/gf45bi50/gf45bi501.tif", pattern="grd", full.names=TRUE) 
#e = extent(c(-122, -117.5, 35, 41))
e = extent(c(-125, -116, 32, 42))
tmp.future = crop(r, e)/10
plot(tmp.future)

# historical
r = raster("./climatedata/bio_2-5m_bil/bio1.bil", pattern="grd", full.names=TRUE) 
#e = extent(c(-122, -117.5, 35, 41))
e = extent(c(-125, -116, 32, 42))
tmp.hist = crop(r, e)/10
plot(tmp.hist)

# change
plot(tmp.future - tmp.hist)
#points(plotlocs)


########################################
# Next do ensemble average (all models for RCP45)
modelnames = c("ac", "bc", "cc", "ce", "cn", "gd", "gf", "gs", "hd", "hg", "he", "in", "ip", "mi", "mr", "mc", "mp", "mg", "no")

wc_ensemble <- function(modelnames, path, variable, extent) { # gets bioclim variable listed as an integer in "variables" from CMIP5 model outputs listed in "models" that is in folders at "path". Crops them to "extent" and returns them as a raster brick. 
  r.list = list()
  for (i in 1:length(modelnames)){
    f = paste(path, "/", modelnames[i], "45bi50/", modelnames[i], "45bi50", variable, ".tif", sep="")
    r = raster(f, pattern="grd", full.names=TRUE)
    r.list[i] = crop(r, e)
  }
  return(brick(r.list))
}

#e = extent(c(-122, -117.5, 35, 41))
e = extent(c(-125, -116, 32, 42)) # I'm cropping this to the Sierra region, but you could use a different extent

# Precipitation -- mean, change and variance 
ppt.ensemble = wc_ensemble(modelnames, path="/Volumes/LaCie/Worldclim", variable=12, extent=e)

ppt.2050.mean = calc(ppt.ensemble, mean)

# historical
r = raster("./climatedata/bio_2-5m_bil/bio12.bil", pattern="grd", full.names=TRUE) 
ppt.hist = crop(r, e)

# Change calculated for each model 
ppt.ensemble.change = ppt.ensemble - ppt.hist
ppt.change = calc(ppt.ensemble.change, mean)
ppt.change.sd = calc(ppt.ensemble.change, sd)
plot(ppt.change)

plot(ppt.change.sd) # variance among models in precip change


# display
par(mfrow=c(1,2))
plot(ppt.change)
map("county", "CA", col="black", add=T)
points(plotlocs, pch=".")
plot(ppt.change.sd)
map("county", "CA", col="black", add=T)
points(plotlocs, pch=".")
# Ensemble mean says precip will increase slightly in lowlands and decrease slightly at higher elevations. But the uncertainty is huge, and overwhelms the mean prediction. So, even for mid-century, model variance is the dominant source of uncertainty for precipitation. 


# Seems like very high variance in change. Check individual models
par(mfrow=c(3, 4))
plot(ppt.ensemble.change) # wildly different predictions
# not clear that the mean really means much, given the uncertainty


# Temperature
tmp.ensemble = wc_ensemble(modelnames, path="/Volumes/LaCie/Worldclim", variable=1, extent=e)

tmp.2050.mean = mean(tmp.ensemble)/10
plot(tmp.2050.mean)

tmp.2050.var = calc(tmp.ensemble/10, var)
plot(tmp.2050.var)

tmp.ensemble.change = (tmp.ensemble/10-tmp.hist)
tmp.change = calc(tmp.ensemble.change, mean)
tmp.change.sd = calc(tmp.ensemble.change, sd)
par(mfrow=c(1,2))
plot(tmp.change, col=heat.colors(128))
map("county", "CA", col="black", add=T)
points(plotlocs, pch=".")
plot(tmp.change.sd, col=heat.colors(128))

map("county", "CA", col="black", add=T)
points(plotlocs, pch=".")
# Interesting -- sites further inland expected to have greater temperature change and higher uncertainty. Change gets pretty high in upper elevations on west side as well as east side. But the uncertainty doesn't get very high until the east side and on into the great basin.




###### NEXT: Merge predicted change with the PRISM data in the FIA data files. 

load("sierra_prism_rasters.Rdata")

tmp.change.resamp = crop(tmp.change,extent(pptnorm.sierra))
tmp.change.resamp = resample(tmp.change.resamp,pptnorm.sierra)
newvals = getValues(tmp.change.resamp)
newvals[is.na(getValues(tmeannorm.sierra))] = NA
tmp.change.sierra = setValues(tmp.change.resamp, newvals)
tmp.future.sierra = tmeannorm.sierra + tmp.change.sierra
tmp.sd.resamp = crop(tmp.change.sd, extent(pptnorm.sierra))
tmp.sd.resamp = resample(tmp.sd.resamp,pptnorm.sierra)
newvals = getValues(tmp.sd.resamp)
newvals[is.na(getValues(tmeannorm.sierra))] = NA
tmp.sd.sierra = setValues(tmp.sd.resamp, newvals)

plot(tmp.change.sierra)
plot(tmp.sd.sierra)
plot(tmp.future.sierra)


ppt.change.resamp = crop(ppt.change,extent(pptnorm.sierra))
ppt.change.resamp = resample(ppt.change,pptnorm.sierra)
newvals = getValues(ppt.change.resamp)
newvals[is.na(getValues(pptnorm.sierra))] = NA
ppt.change.sierra = setValues(ppt.change.resamp, newvals)
ppt.future.sierra = pptnorm.sierra + ppt.change.sierra
ppt.sd.resamp = crop(ppt.change.sd, extent(pptnorm.sierra))
ppt.sd.resamp = resample(ppt.sd.resamp,pptnorm.sierra)
newvals = getValues(ppt.sd.resamp)
newvals[is.na(getValues(tmeannorm.sierra))] = NA
ppt.sd.sierra = setValues(ppt.sd.resamp, newvals)


plot(ppt.change.sierra)
plot(ppt.sd.sierra)

plot(pptnorm.sierra)
plot(ppt.future.sierra) # changes are really small


save(tmp.change.sierra, tmp.sd.sierra, tmp.future.sierra, ppt.change.sierra, ppt.sd.sierra, ppt.future.sierra, file="sierra_change_rasters.Rdata")



# Make picture of ensemble projections for presentation
par(mfrow=c(1,2), mar=rep(4, 4))
plot(ppt.future)
plot(tmp.future, yaxt="n")
