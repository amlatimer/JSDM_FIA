# Look at and organize Sierra forest FIA data 


plotpath = "../plots/" # place to put plots for easy access from Powerpoint

library(sp)
library(maps)
library(fields)
library(vegan)
library(gjam)
library(raster)
library(FNN)
library(mvtnorm)
library(tmvtnorm)

source("latimer_gjam_functions.R")

#############################
# Read in the data
#############################

#d = read.csv("FIA_WB_model_full_data_forAML.csv")
d_fireinfo = read.csv("../data/FinalSierraCWB_FIA_withCN_RadHours_fuzzedCoords.csv")
d_site = read.csv("../data/FIA_non-WB_predictors.csv")
d_ba = read.csv("../data/FIA_SPxBA.csv")

# compare coverage of the data sets
plot(d_fireinfo[,c("long", "lat")])
nrow(d_fireinfo)
points(d_fireinfo[d_fireinfo$FIA_CN2 %in% d_ba$CN,c("long", "lat")], col=2)
nrow(d_ba)

plot(d_site[,c("LON", "LAT")], col="blue")
points(d_site[d_site$CN %in% d_ba$CN,c("LON", "LAT")], col=2)


# It looks like the best data set for modeling sierra forest in general (not including fire info other than what's in FIA) is the combination of FIA_nonWB_predictors with the ba data set: 1139 plots. 

# Merge the info together
d = merge(d_site, d_ba, by="CN") # combine basal area and site info
plot(d[,c("LON", "LAT")])
d_fire = merge(d, d_fireinfo[,c("FIA_CN2", "YrLastFire")], by.x="CN", by.y="FIA_CN2", all.x=TRUE, all.y=FALSE)

names(d_fire)
plotlocs= SpatialPoints(d_fire[,c("LON", "LAT")])

# Load prism raster layers
load("../working-files/sierra_prism_rasters.Rdata")
load("../working-files/sierra_change_rasters.Rdata")

# Alternatively, load prism historical normals as raster first and then crop
#pptnorm = raster("./climatedata/PRISM_ppt_30yr_normal_4kmM2_annual_asc/PRISM_ppt_30yr_normal_4kmM2_annual_asc.asc")
#tmeannorm = raster("./ClimateData/PRISM_tmean_30yr_normal_4kmM2_annual_asc/PRISM_tmean_30yr_normal_4kmM2_annual_asc.asc")
# Select part of raster that matches area covered by Sierra FIA plots
#bb = extent(c(-122, -117, 34.5, 41))
#pptnorm = crop(pptnorm, bb)
#tmeannorm = crop(tmeannorm, bb)
#plot(tmeannorm)
#plot(plotlocs, pch=16, add=T, cex=0.5)
#snoutline = drawPoly()
#save(snoutline, file= "sierras.outline.Rdata")
#load("sierras.outline.Rdata")
#plot(snoutline,add=T)
# next get cell numbers for all cells in pptnorm inside the polygon, then set all other cell values (seValues) to NA
# And use the index of cell numbers to assign values and display model outputs. 
# can also then convert the raster to polygons using rasterToPolygons -- leaving behind the NA values. 

# Make precip raster layer
#cells.inside = extract(pptnorm, snoutline, cellnumbers=TRUE)[[1]]
#head(cells.inside)
#pptnorm.sierra = pptnorm
#cellvals = rep(NA, length(pptnorm.sierra))
#cellvals[cells.inside[,1]] = cells.inside[,2]
#pptnorm.sierra = setValues(pptnorm.sierra, cellvals)
#plot(pptnorm.sierra, col=heat.colors(100))

# make temperature raster layer
#cells.inside = extract(tmeannorm, snoutline, cellnumbers=TRUE)[[1]]
#cell.index = cells.inside[,1]
#tmeannorm.sierra = tmeannorm
#cellvals = rep(NA, length(tmeannorm.sierra))
#cellvals[cells.inside[,1]] = cells.inside[,2]
#tmeannorm.sierra = setValues(tmeannorm.sierra, cellvals)
#plot(tmeannorm.sierra, col=heat.colors(100))

# Write out the rasters and index of cells inside Sierra domain
#save(tmeannorm.sierra, pptnorm.sierra, cell.index,  file="sierra_prism_rasters.Rdata")

# Check the resulting rasters
plot(pptnorm.sierra, xlim=c(-122, -117), ylim=c(35, 41), xlab="LON", ylab="LAT", cex.axis=1.2, cex.lab=1.2)
plot(plotlocs, add=T)
plotmap(plotlocs, d_fire$ppt.tot, rain.colors)
d_fire$pptnew = extract(pptnorm.sierra, plotlocs)

plot(pptnew~ppt.tot, d_fire)

plot(tmeannorm.sierra, xlim=c(-122, -117), ylim=c(35, 41), xlab="LON", ylab="LAT", cex.axis=1.2, cex.lab=1.2)
#plot(plotlocs, add=T)
plotmap(plotlocs, d_fire$tmean.mean, rain.colors)
d_fire$tmeannew = extract(tmeannorm.sierra, plotlocs)
plot(tmeannew~tmean.mean, d_fire)
abline(c(0,1))
# Values for plots are close to, but not exactly same as in Derek's data set, because of fuzzing. So: use Derek's data for fitting (extracted from the actual locations), but map projections to the PRISM grid. 

# Make pictures of the climate layers (historical normals)
#png(filename=paste(plotpath, "pptnorm_PRISM.png", sep=""))
#plot(pptnorm, xlim=c(-122, -117), ylim=c(35, 41), xlab="LON", ylab="LAT", cex.axis=1.2, cex.lab=1.2)
#dev.off()
#png(filename=paste(plotpath, "tmeannorm_PRISM.png", sep=""))
#plot(tmeannorm, xlim=c(-122, -117), ylim=c(35, 41), xlab="LON", ylab="LAT", cex.axis=1.2, cex.lab=1.2)
#dev.off()

# Create data set for projecting model results to sierra raster
# Includes precip and temp
x.rast = data.frame(ppt.tot=pptnorm.sierra[cell.index], tmean.mean = tmeannorm.sierra[cell.index])
head(x.rast); dim(x.rast)

# Select a subset of common species to model
names(d_fire) 
spdata = d[,239:268]
prev = apply(spdata, 2, f<-function(x) {return(sum(x>0))})
rev(sort(prev)) 
# Select the most abundant species
spp = names(prev)[prev>=50] # take species that occur in at least 50 plots

# If it's still in the data set, remove mountain mahogany (doesn't converge, and it's not really a tree)
#spp = spp[spp != "sp.CELE3"]

# Set up the species basal area as the response matrix
y = spdata[,spp]
nsp = dim(y)[2]

# look at empty rows in y
ba_obs = apply(y, 1, sum)
sum(ba_obs==0)
ynonzero = which(ba_obs>0)
hist(ba_obs)
plotmap(plotlocs, ba_obs, heat.colors)
points(plotlocs[ba_obs==0,], pch=16, col="black")
# these are mostly at the edge of forested areas

# Remove empty rows in y and from d
y = y[ynonzero,]
d = d[ynonzero,]
d_fire = d_fire[ynonzero,]
plotlocs = plotlocs[ynonzero]

# From here on, let's call the data frame d_fire (including time of last fire) "d" for simplicity. 

d <- d_fire

# Add more variables to the analysis data frame 
d$sin.slope = d$SlopeFIA/180*pi
d$sin.slope.sin.aspect = sin(d$SlopeFIA/180*pi) * sin(d$AspectFIA/180*pi)
d$sin.slope.cos.aspect = sin(d$SlopeFIA/180*pi) * cos(d$AspectFIA/180*pi)

# check these correlations
cor(d[,c("ppt.tot", "tmean.mean", "sin.slope.sin.aspect", "sin.slope.cos.aspect")])
plot(sin.slope.sin.aspect~rad.07, data=d)
plot(sin.slope.cos.aspect~rad.07, data=d) # weirdly, not much correlation here



###################################
### READY TO RUN GJAM MODEL.  

# Data type for tree basal area is "CA" continuous abundance. 
# Per manual the default allows for a point mass at 0 (i.e. deals with zero inflation)
#   by applying a Tobit model.

# Effort for california FIA plots is all the same by design. 

# Model specification
ml = list(ng=2000, burnin=1000, typeNames="CA")#, reductList = list(r=4, N=8))
# note including dimension reduction may be necessary to avoid overfitting; see plots of inverse-predicted temp vs observed temp below

# Test run on all data 
out <- gjam(~ pptnew + tmeannew + I(pptnew^2) + I(tmeannew^2) + pptnew:tmeannew + rad.tot + sin.slope.cos.aspect+ sin.slope.sin.aspect, xdata=d, ydata=y, modelList=ml)

#~ ppt.tot + I(ppt.tot^2) + tmean.mean + I(tmean.mean^2) + ppt.tot:tmean.mean 

# One issue is whether to include topographic information. The potential problem is for projecting forward onto a coarse raster, there is no meaningful slope or aspect at that scale.

# But we could alternatively project onto a fine raster, or onto the FIA point locations. 

pl <- list(GRIDPLOTS = T, SMALLPLOTS =F, CLUSTERPLOTS=T)
gjamPlot(out, plotPars=pl)


###### Mapping imputed values of temperature and precipitation, and comparing those to observed. Are there any systematic over- or underestimates of temperature (or precip)?

plotlocs <- SpatialPoints(coord=d[,c("Lon", "Lat")])

plot(d$tmeannew, out$prediction$xpredMu[,3], pch=16, cex=0.5, col="darkgray")
abline(0,1)
points(d$tmeannew[d$YrLastFire>2005], out$prediction$xpredMu[d$YrLastFire>2005,3], pch=16, cex=0.5, col="red")
# no obvious association with fire. Maybe a hint that recently burned areas are "warmer" in the mid-elevation band.
tmean_errors <- out$prediction$xpredMu[,3] - d$tmeannew
tmean_errors[is.na(tmean_errors)] <- 0 # set the couple of NA's to 0 for plotting purposes

hist(tmean_errors)
# for plotting purposes, squeeze the negative outliers back to normal range
tmean_errors[tmean_errors < (-max(tmean_errors))] <- -max(tmean_errors)
plotmap(plotlocs, tmean_errors, rain.colors)
recent_fire <- which(d$YrLastFire>1996)
points(plotlocs[recent_fire], pch=1 )
# arbitrarily, color the map by whether there's a positive or negative error >0.2
tmean_errors_cat <- tmean_errors
tmean_errors_cat[abs(tmean_errors)<0.2] <- 0
tmean_errors_cat[tmean_errors > 0.2] <- 1
tmean_errors_cat[tmean_errors < -0.2] <- -1
plotmap(plotlocs, tmean_errors_cat, rain.colors)

plot(tmean_errors, d$tmeannew)
plot(tmean_errors, d$pptnew)
plot(tmean_errors, d$rad.tot)
plot(tmean_errors, d$cloudiness.07)
plot(tmean_errors, d$Deficit.BCM)
plot(tmean_errors~ d$elev)
abline(0,0, lwd=2, col="gold")

plot(d$pptnew, out$prediction$xpredMu[,2], pch=16, cex=0.5, col="darkgray")
pptmean_errors <- out$prediction$xpredMu[,2] - d$pptnew
hist(pptmean_errors)
pptmean_errors[is.na(pptmean_errors)] <- 0 # set the couple of NA's to 0 for plotting purposes
plotmap(plotlocs, pptmean_errors, rain.colors)
# interestingly, the plots that appear "wetter" than normal climate, i.e. have a wetter mix of trees in them, have a distinct spatial pattern and in the south might be associated with higher mortality in 2015-16?
# This might be because there are more than expected ABCO there? 

# compare errors in precip vs temperature inverse prediction
plot(tmean_errors, pptmean_errors) # moderate positive association
cor(tmean_errors, pptmean_errors)

# particular species associated with under-prediction? 
plot(pptmean_errors~y$sp.ABCO) # low ABCO cover associated with underprediction of precipitation
plot(pptmean_errors~y$sp.ABMA) 
plot(pptmean_errors~y$sp.JUOC) 
plot(pptmean_errors~y$sp.CADE27)
plot(pptmean_errors~y$sp.PICO) # high PICO cover associated with underpredictio of precip
plot(pptmean_errors~y$sp.PIJE)
plot(pptmean_errors~y$sp.PILA)
plot(pptmean_errors~y$sp.PIPO)
plot(pptmean_errors~y$sp.PISA2)
plot(pptmean_errors~y$sp.PIMO)
plot(pptmean_errors~y$sp.PSME)
plot(pptmean_errors~y$sp.QUCH2)
plot(pptmean_errors~y$sp.QUKE)
plot(pptmean_errors~y$sp.QUDO)
plot(pptmean_errors~y$sp.QUWI2)


### NEXT SHOULD MAP UNCERTAINTY IN INVERSE PREDICTION (xpredSD)
### AND COMPARE ERRORS TO MORTALITY IN DROUGHT








###### TESTING OUT OF SAMPLE IMPUTATION OF Y
# Simply deleting a row of the y data
y.test = y
holdouts <- seq(10,1100, by=10)
y.test[holdouts,] = NA
ml = list(ng=1000, burnin=500, typeNames="CA")
out.test <- gjam(~ ppt.tot + I(ppt.tot^2) + tmean.mean + I(tmean.mean^2) + ppt.tot:tmean.mean , xdata=d, ydata=y.test, modelList=ml)

# look at results
par(mfrow=c(4,4), mar=rep(2,4))
for (i in 1:16) plot(y[holdouts,i], out.test$prediction$ypredMu[holdouts,i], main=names(y)[i])
# of course, way fewer zeros in predictions than in data, but some species not predicted badly. mixed bag. 







#######################
# Explore model summary

# Model DIC and Gneiting-Raftery scores
out.test$modelSummary$DIC
out$modelSummary$yscore

# Fitted values of total basal area
dev.off(); plot.new()
ba.fit = apply(out$prediction$ypredMu, 1, sum)
ba.obs = apply(y,1,sum)
plot(ba.fit~ba.obs, xlim=c(0, 130), ylim=c(0,130))
abline(0,1)
# underspecified model -- variation in prediction much less than in observations, but not obviously over or under predicting.  

# Obs vs fit for species
plot.new()
par(mfrow=c(3,3))
for (i in 1:9) plot(y[,i], out$prediction$ypredMu[,i], ylab="BA model fit", xlab = "BA observed")
for (i in 10:18) plot(y[,i], out$prediction$ypredMu[,i], ylab="BA model fit", xlab = "BA observed")

# holdouts (if any data were held out in model fitting)
par(mfrow=c(1,2))
xMu  <- out$prediction$xpredMu
xSd  <- out$prediction$xpredSd
yMu  <- out$prediction$ypredMu
hold <- out$holdoutIndex

plot(out$x[hold,-1],xMu[hold,-1], cex=.2)
title('holdouts in x'); abline(0,1)
plot(out$y[hold,], yMu[hold,], cex=.2)
title('holdouts in y'); abline(0,1)
plot(apply(out$y[hold,],1, sum), apply(yMu[hold,], 1, sum), cex=.2)
title('holdouts in y, total by plot'); abline(0,1)inflation in the plots, which is largely because of their small size (?)

# KEY Q: How much should we worry about this? Can we still take model projections about basal area seriously? I would guess that for community change it should still be OK, but not sure about basal area. 


##### Display model coefficients and chains

# Ideas: 
# 1) Using community to predict environmental data. Hold out and predict environmental variables from the community data. (possibly: include this in interpolation? I.e. use species composition to improve climate surface creation from weather stations? probabl too big a task)

# Hold out and predict community values at sites to check accuracy. 

# 2) Use inverse prediction to predict vulnerability to drought. Identify areas on the map where sensitivity to CWD is high -- then see if those are the same areas where aerial mortality survey identifies high mortality. 

# 2) Use community response across space to predict community response to drought (roughly, info contained in E matrix).  Us fitted model to characterize how community abundances change when CWD increases. Does this qualitatively match the set of species that are dying in 2015 in the aerial mortality survey? 


##########################################
## Work on Idea 1: Using JSDM to predict environmental variables
##  
##  First to assess how well this works.  (q do with other data sets -- which? protea atlas for fynbos weather stations?)
## Second to assess sensitivity

### Map residuals of predictions of variables

v = d$ppt.tot - out$prediction$xpredMu[,2]
plot(v~d$ppt.tot)
pcols = (v*10-min(v*10))+1
rain.colors = colorRampPalette(c("red", "gray", "blue"))
palette(rain.colors(max(pcols)))
hist(v)
plot(plotlocs, pch=16, cex=0.7,col=pcols)
# no real clear pattern of under/over pred. 

# temperature
v = d$tmean.mean - out$prediction$xpredMu[,3]
plot(v~d$tmean.mean)
pcols = (v*10-min(v*10))+1
palette(rain.colors(max(pcols)))
plot(plotlocs, pch=16, cex=0.7, col=pcols)
# no real clear pattern of under/over pred. 


# Using community information to predict environment

# 1) First, try all community data across whole region 

# Leave out 10% of environmental data and predict
holdout = sample(1:nrow(d[ynonzero,]), 113, replace=FALSE)

# Map to show where the holdouts are
plot(plotlocs, col="gray")
points(plotlocs[holdout], col="red", pch=16)
map("county", "CA", add=T, fill=FALSE, col="black")

# fit model and predict held-out values
ml = list(ng=2000, burnin=1000, typeNames="CA", holdoutIndex=holdout)
outpred <- gjam(~ppt.tot + I(ppt.tot^2) + tmean.mean + I(tmean.mean^2) + ppt.tot:tmean.mean, d, y, modelList=ml)

# scatterplots for precip
#png(filename=paste(plotpath, "pred_x_fullmodel_ppt.tot.png", sep=""))
plot(outpred$inputs$x[outpred$modelList$holdoutIndex,"ppt.tot"],outpred$prediction$xpredMu[outpred$modelList$holdoutIndex,"ppt.tot"])
abline(0,1)
#dev.off()

# scatterplot for temp 
#png(filename=paste(plotpath, "pred_x_fullmodel_tmean.png", sep=""))
plot(outpred$inputs$x[outpred$modelList$holdoutIndex,"tmean.mean"],outpred$prediction$xpredMu[outpred$modelList$holdoutIndex,"tmean.mean"], xlab="", ylab=""); abline(0,1)
#dev.off()

# The holdout predictions are quite good for temperature, worse but still correlated for precipitation.


# 2) To gauge how important the joint information is, compare how well X can be predicted when fitting and predicting with only one common species (or a couple). I.e. a more standard-issue distrib model. 

y_sub = y[,c("sp.ABMA", "sp.PIPO")] # Use only one warm-temp pine and one cool-temp fir
ml_sub = list(ng=2000, burnin=1000, typeNames="CA", holdoutIndex=holdout)
outpred_sub <- gjam(~ppt.tot + I(ppt.tot^2) + tmean.mean + I(tmean.mean^2) + ppt.tot:tmean.mean, d, y_sub, modelList=ml)

plot(outpred_sub$inputs$x[outpred_sub$modelList$holdoutIndex,"tmean.mean"], outpred_sub$prediction$xpredMu[outpred_sub$modelList$holdoutIndex,"tmean.mean"])
abline(0,1)
# Not terrible, but much more scatter, and less discrimination at higher temperatures. 

plot(outpred_sub$inputs$x[outpred_sub$modelList$holdoutIndex,"ppt.tot"],outpred_sub$prediction$xpredMu[outpred_sub$modelList$holdoutIndex,"ppt.tot"])
abline(0,1)
# Precipitation is a mess, however.  Shows that having more species helps, as expected. 


# 3)  Use only pres-abs info for the whole community
y_presabs = apply(y, 2, f<-function(x){return(as.integer(x>0))})
ml_presabs = list(ng=2000, burnin=1000, typeNames="PA", holdoutIndex=holdout)
outpred_presabs <- gjam(~ppt.tot + I(ppt.tot^2) + tmean.mean + I(tmean.mean^2) + ppt.tot:tmean.mean, d, y_presabs, modelList=ml)

plot(outpred_presabs$inputs$x[outpred_presabs$modelList$holdoutIndex,"tmean.mean"],outpred_presabs$prediction$xpredMu[outpred_presabs$modelList$holdoutIndex,"tmean.mean"])
abline(0,1)
plot(outpred_presabs$inputs$x[outpred_presabs$modelList$holdoutIndex,"ppt.tot"],outpred_presabs$prediction$xpredMu[outpred_presabs$modelList$holdoutIndex,"ppt.tot"])
abline(0,1)
# Presabs data somewhat worse for inverse prediction of climate, but it's far better to have presabs data for most of the community than abundance data for only a small part of it


# IDEA: Predict relative suitablility of site (or predicted abundance) for trees Derek has sampled. Does the predicted abundance match growth rates or environmental sensitivity of those species in those grid cells (use his plot-level averages?)

# 5) Success across the whole Sierra mountains, including huge temperature and precip contrasts, is not so surprising. How well can this do with prediction within a certain forest type or elevation band? Q How to do this?? Restrict to only sites with mixed conifer species in them? 
mixcon = c("sp.PILA", "sp.CADE27", "sp.PILA", "sp.PIPO", "sp.PSME", "sp.QUKE")#, "sp.PIJE")
mixcon.ind = apply(y, 1, f<-function(x, spnames) {return(sum(spnames[x>0] %in% mixcon))}, spnames=colnames(y))
plot(plotlocs, col="gray", axes=T)
points(plotlocs[mixcon.ind>0], col="darkgreen")
map("county", "CA", col="black", add=T)
hist(mixcon.ind)
sum(mixcon.ind>0)

d_mixcon = d[mixcon.ind>0,]
y_mixcon = y[mixcon.ind>0,]
prev_mixcon = apply(y_mixcon, 2, f<-function(x){return(sum(x>0))})
y_mixcon = y_mixcon[,prev_mixcon>=30] # remove species with fewer than 30 occurrences in this narrower data set 
colnames(y_mixcon); nrow(y_mixcon)
# This leaves 9 species in 489 plots

ml_mixcon = list(ng=2000, burnin=1000, typeNames="CA", holdout=50)
outpred_mixcon <- gjam(~ppt.tot + I(ppt.tot^2) + tmean.mean + I(tmean.mean^2) + ppt.tot:tmean.mean, d_mixcon, y_mixcon, modelList=ml_mixcon)

plot(outpred_mixcon$x[outpred_mixcon$holdoutIndex,"tmean.mean"], outpred_mixcon$modelSummary$xpredMu[outpred_mixcon$holdoutIndex,"tmean.mean"], xlab="observed", ylab="predicted", cex.lab=1.5, cex.axis=1.2) 

MSPE = mean((outpred_mixcon$x[outpred_mixcon$holdoutIndex,"tmean.mean"]- outpred_mixcon$modelSummary$xpredMu[outpred_mixcon$holdoutIndex,"tmean.mean"])^2)
title(main="Mean annual temperature -- Mixed conifer forest", sub=paste("MSPE =", round(MSPE, 3)))
abline(0,1)

# Much worse prediction, temperature not no longer super strong predictor as it is across the whole elevation profile. 

summary(d$tmean.mean)
summary(d$tmean.mean[mixcon.ind>0]) # this is about 2/3 the range in the total data set. 
summary(d$ppt.tot)
summary(d$ppt.tot[mixcon.ind>0]) # almost all the range in precip, minus the driest areas




#########################
# Issue of fire disturbance -- obviously not all plots at "equilibrium" -- can we detect differences between recently burned plots (burned since 1980, i.e. within about 30 years of the most recent FIA survey data for most plots) and not burned in more than 60 years (not since 1950)? 

y_fire = d_fire[,239:268] 
y_fire = y_fire[,spp] 
prev_fire = apply(y_fire, 2, f<-function(x){return(sum(x>0))})
y_fire = y_fire[,prev_fire>=100] # remove species with fewer than 100 occurrences 
colnames(y_fire); nrow(y_fire)
# this leaves 10 species in 716 plots

# Add binary variable for "recent fire" 
d_fire$recentfire = 0
d_fire$recentfire[d_fire$YrLastFire>0] = 1

par(mfrow=c(1,1))
plot(SpatialPoints(d[,c("LON", 'LAT')]), col="gray", pch=1, cex.axis=1.2, cex.lab=1.2, axes=T, ylab="LAT", xlab="LON")
points(d_fire[,c("LON", 'LAT')])
points(d_fire[d_fire$YrLastFire>1980,c("LON", 'LAT')], col="red", pch=16)
map("county", "CA", col="black", add=T)
legend(-119, 40, c("Fire data", "Recent fire"), col=c("black", "red"), pch=c(1, 16),box.col= "white", cex=1.4)

ml_fire = list(ng=2000, burnin=1000, typeNames="CA",holdoutN=50)
out_fire <- gjam(~ppt.tot + I(ppt.tot^2) + tmean.mean + I(tmean.mean^2) + ppt.tot:tmean.mean + recentfire, d_fire, y_fire, modelList=ml_fire)

pl <- list(GRIDPLOTS = T, SMALLPLOTS =TRUE,CLUSTERPLOTS=T)
gjamPlot(out_fire, plotPars=pl)

plot(out_fire$inputs$x[out_fire$modelList$holdoutIndex,"tmean.mean"],out_fire$prediction$xpredMu[out_fire$modelList$holdoutIndex,"tmean.mean"]); abline(0,1)
plot(out_fire$inputs$x[out_fire$modelList$holdoutIndex,"ppt.tot"],out_fire$prediction$xpredMu[out_fire$modelList$holdoutIndex,"ppt.tot"]);  abline(0,1)

# Overall poorer fit because smaller data set and covers subset of environmental range. 
# inverse pred of temp is still good, inverse pred of fire is pretty good
# several species "significantly" affected by fire, red fire with most positive relationship, oddly 
# But fire in general is a weak predictor. Probably because lumping all fire severities together?


#################################################################




#################################################################
# Work on Idea 2 -- environmental sensitivity and projected change


# For these questions, the main task is to project the community composition forward in time using the climate model outputs 

# First, try projecting the community composition out to the full raster area using the 30-year historical normal data

# Model specification

# Q do we need to re-run the model using unstandardized values, or will the model recognize that new data are not yet standardized?

ml = list(ng=2000, burnin=1000, typeNames="CA")
out.current = gjam(~ ppt.tot + I(ppt.tot^2) + tmean.mean + I(tmean.mean^2) + ppt.tot:tmean.mean, xdata=d, ydata=y, modelList=ml)

gjamPredict(out.current, y2plot = "sp.ABCO")


# Try fitting model with missing y values for raster cells
index.sierra = which(!is.na(getValues(pptnorm.sierra)))
x.norm = data.frame(pptnew = getValues(pptnorm.sierra)[index.sierra], tmeannew = getValues(tmeannorm.sierra)[index.sierra])

x.norm$ppt.tot = x.norm$pptnew # (x.norm$pptnew-mean(d$pptnew, na.rm=T))/sd(d$pptnew, na.rm=T)
x.norm$tmean.mean = x.norm$tmeannew # (x.norm$tmeannew-mean(d$tmeannew, na.rm=T))/sd(d$tmeannew, na.rm=T)

pred.current = gjamPredict(out.current, newdata = list(xdata=x.norm))

# visualize result for one species
rastnew = pptnorm.sierra
rastnew[index.sierra] =  pred.current$sdList$yMu[,7]
plot(rastnew) 

hist(d$sp.ABCO)
hist(predtest$sdList$yMu[,1])

# Looks good. Now predict forward to future climate and look at change 


x.future = data.frame(ppt.tot = getValues(ppt.future.sierra)[index.sierra], tmean.mean = getValues(tmp.future.sierra)[index.sierra])

pred.future = gjamPredict(out.current, newdata = list(xdata=x.future))

# visualize result for one species
rastnew = pptnorm.sierra
par(mfrow=c(1, 2))
rastnew[index.sierra] =  pred.current$sdList$yMu[,1]
plot(rastnew) 
rastnew[index.sierra] =  pred.future$sdList$yMu[,1]
plot(rastnew) 


# NEXT: Make projections for future conditions for full community and calculate change metrics . . . . community composition, BA, etc. 

ba_present <- apply(pred.current$sdList$yMu, 1, sum)
ba_future <- apply(pred.future$sdList$yMu, 1, sum)
ba_change <- ba_future - ba_present
par(mfrow=c(1, 3))
rastnew[index.sierra] =  ba_present; plot(rastnew) 
rastnew[index.sierra] =  ba_future; plot(rastnew) 
rastnew[index.sierra] =  ba_change; plot(rastnew, col=rain.colors(64)) 

sum(ba_change) / sum(ba_present) # 3% decline in ba, but 
hist(ba_change/ba_present) # much stronger localized effects ranging from -50% to +100%

# Calculate jaccard dissimilarity between present and future communities
ncells <- length(index.sierra)
comm_jaccard <- rep(NA, ncells)
for (i in 1:ncells) {
  comm_jaccard[i] <- vegdist(rbind(pred.current$sdList$yMu[i,], pred.future$sdList$yMu[i,]), method = "jaccard")
}
rastnew[index.sierra] =  comm_jaccard; plot(rastnew)
# most change at lower elevation edge and in (apparently) subalpine areas  
# little change in some mid-elevation areas that will remain MEC, and some alpine areas that will remain alpine. 













##########################################################
#### Build a prediction function that will simulate values of y from fitted model for given values of x

# First, fit the model if not done already
#ml = list(ng=2000, burnin=1000, typeNames="CA")
#out2 <- gjamGibbs(~ ppt + ppt2 + tmean + tmean2 +pptxtemp, x, ytest, modelList=ml)
#pl <- list(width = 3,height = 2, GRIDPLOTS = T, SMALLPLOTS =T  ,CLUSTERPLOTS=T)
#gjamPlot(out2, plotPars=pl)
out2=out_full

# NOTE PIMO and PISA2 aren't converging well -- maybe have to remove (only present on fringes anyway)

# Predict all values of y (species abundances) for new values of x, drawing from posterior distribution of the model variables for each simulation. Assumes that the data are continuous abundance data. 

# First select a standard set of posterior samples to use (for consistency across data sets)
samples = sample(1:20000, size=100, replace=FALSE)


# Predict to Sierra raster cells for number of posterior simulations, which later can be used for inference about community composition, richness, etc. 
ypred.sierra = gjamPredictY.region(x.rast, out_full, samples=samples)

ypred.future = gjamPredictY.region(newx.future, out_full, samples=samples)

ypred.change = ypred.future - ypred.sierra





# Calculate mean and sd for species abundance -- historical 
ypred.sierra.mean = matrix(NA, nrow=nsp, ncol=nrow(x.rast))
for (i in 1:nsp) for (j in 1:nrow(x.rast)) ypred.sierra.mean[i,j] = mean(ypred.sierra[i,j,])
ypred.sierra.sd = matrix(NA, nrow=nsp, ncol=nrow(x.rast))
for (i in 1:nsp) for (j in 1:nrow(x.rast)) ypred.sierra.sd[i,j] = sd(ypred.sierra[i,j,])
ypred.sierra.cv = ypred.sierra.sd/ypred.sierra.mean

# Calculate mean and sd for species abundance -- future
ypred.future.mean = matrix(NA, nrow=nsp, ncol=nrow(x.rast))
for (i in 1:nsp) for (j in 1:nrow(x.rast)) ypred.future.mean[i,j] = mean(ypred.future[i,j,])
ypred.future.sd = matrix(NA, nrow=nsp, ncol=nrow(x.rast))
for (i in 1:nsp) for (j in 1:nrow(x.rast)) ypred.future.sd[i,j] = sd(ypred.future[i,j,])
ypred.future.cv = ypred.future.sd/ypred.future.mean



# Calculate mean and sd species richness
n.samples = dim(ypred.sierra)[3]
richness.sierra = matrix(NA, nrow=nrow(x.rast), ncol=n.samples)
for (j in 1:nrow(x.rast)) for (k in 1:n.samples) richness.sierra[j,k] = sum(ypred.sierra[,j,k]>0)
richness.sierra.mean = apply(richness.sierra, 1, mean)
richness.sierra.sd = apply(richness.sierra, 1, sd)
plot.layer.raster(pptnorm.sierra, p=richness.sierra.mean)
plot.layer.raster(pptnorm.sierra, p=richness.sierra.sd)

# Calculate mean and sd basal area
ba.sierra = matrix(NA, nrow=nrow(x.rast), ncol=n.samples)
for (j in 1:nrow(x.rast)) for (k in 1:n.samples) ba.sierra[j,k] = sum(ypred.sierra[,j,k])
ba.sierra.mean = apply(ba.sierra, 1, mean)
ba.sierra.sd = apply(ba.sierra, 1, sd)
plot.layer.raster(pptnorm.sierra, p=ba.sierra.mean)
plot.layer.raster(pptnorm.sierra, p=ba.sierra.sd)

# Calculate change in composition -- Jaccard dissimilarity
jaccard.sierra.future = gjamCalcChange.jaccard(ypred.sierra, ypred.future)
plot.layer.raster(pptnorm.sierra, p=jaccard.sierra.future[[1]])
plot.layer.raster(pptnorm.sierra, p=jaccard.sierra.future[[2]])

# Calculate change in composition -- Euclidean distance
euclid.sierra.future = gjamCalcChange.vegdist(ypred.sierra, ypred.future, "euclidean")
par(mfrow=c(1,2));plot.layer.raster(pptnorm.sierra, p=euclid.sierra.future[[1]])
plot.layer.raster(pptnorm.sierra, p=euclid.sierra.future[[2]])

# Calculate change in composition -- summed absolute ba change
abschange.sierra.future = gjamCalcChange.abschange(ypred.sierra, ypred.future)
plot.layer.raster(pptnorm.sierra, p=abschange.sierra.future[[1]])
plot.layer.raster(pptnorm.sierra, p=abschange.sierra.future[[2]])

# Calculate change in basal area
ba.sierra.change = gjamCalcChange.basalarea(ypred.sierra, ypred.future)
plot.layer.raster(pptnorm.sierra, p=ba.sierra.change[[1]])
plot.layer.raster(pptnorm.sierra, p=ba.sierra.change[[2]])
propchange = ba.sierra.change[[1]]/ba.sierra.mean
propchange[propchange>2] = 2
plot.layer.raster(pptnorm.sierra, p=propchange)


# Calculate change in hardwood basal area 
hw.sierra.change = gjamCalcChange.hardwood(ypred.sierra, ypred.future, spp)
plot.layer.raster(pptnorm.sierra, p=ba.sierra.change[[1]])
plot.layer.raster(pptnorm.sierra, p=ba.sierra.change[[2]])

# Calculate change in species richness
richness.sierra.change = gjamCalcChange.richness(ypred.sierra, ypred.future)
plot.layer.raster(pptnorm.sierra, p=richness.sierra.change[[1]])
plot.layer.raster(pptnorm.sierra, p=richness.sierra.change[[2]])

# WHich species are overall increasing/declining?
plot(apply(ypred.sierra.mean, 1, mean)~apply(ypred.future.mean, 1, mean), col="darkgreen", cex=0.8)

ychange.species = data.frame(species = spp, change=apply(ypred.future.mean, 1, mean)-apply(ypred.sierra.mean, 1, mean))
ychange.species = ychange.species[order(ychange.species$change),]
barplot(ychange.species$change, names.arg=ychange.species$species, horiz=T, las=1, col="darkgreen", cex.names=0.8)







# map one species at a time and compare visually to distribution data 

cell.index = !is.na(getValues(pptnorm.sierra))

plot.spdist.raster(pptnorm.sierra, ypred.sierra.mean[grep("PIPO", spp),], cell.index, "sp.PIPO")
plot.spdist.raster(pptnorm.sierra, ypred.sierra.sd[grep("PIPO", spp),], cell.index, "sp.PIPO")
plot.spdist.raster(pptnorm.sierra, ypred.sierra.cv[grep("PIPO", spp),], cell.index, "sp.PIPO")

plot.spdist.raster(pptnorm.sierra, ypred.sierra.mean[grep("ABMA", spp),], cell.index, "sp.ABMA")
plot.spdist.raster(pptnorm.sierra, ypred.sierra.sd[grep("ABMA", spp),], cell.index, "sp.ABMA")
plot.spdist.raster(pptnorm.sierra, ypred.sierra.cv[grep("ABMA", spp),], cell.index, "sp.ABMA")

# Comparing present to future
# Declining Species
par(mfrow=c(1, 2))
plot.spdist.raster(pptnorm.sierra, ypred.sierra.mean[grep("ABMA", spp),], cell.index, "sp.ABMA")
plot.spdist.raster(pptnorm.sierra, ypred.future.mean[grep("ABMA", spp),], cell.index, "sp.ABMA")
# Declining Species
par(mfrow=c(1, 2))
plot.spdist.raster(pptnorm.sierra, ypred.sierra.mean[grep("PICO", spp),], cell.index, "sp.PICO")
plot.spdist.raster(pptnorm.sierra, ypred.future.mean[grep("PICO", spp),], cell.index, "sp.PICO")


# increasing Species
par(mfrow=c(1, 2))
plot.spdist.raster(pptnorm.sierra, ypred.sierra.mean[grep("PSME", spp),], cell.index, "sp.PSME")
plot.spdist.raster(pptnorm.sierra, ypred.future.mean[grep("PSME", spp),], cell.index, "sp.PSME")

# increasing Species
par(mfrow=c(1, 2))
plot.spdist.raster(pptnorm.sierra, ypred.sierra.mean[grep("QUDO", spp),], cell.index, "sp.QUDO")
plot.spdist.raster(pptnorm.sierra, ypred.future.mean[grep("QUDO", spp),], cell.index, "sp.QUDO")


plot.spdist.raster(pptnorm.sierra, ypred.sierra.mean[grep("QUKE", spp),], cell.index, "sp.QUKE")
plot.spdist.raster(pptnorm.sierra, ypred.sierra.sd[grep("QUKE", spp),], cell.index, "sp.QUKE")
plot.spdist.raster(pptnorm.sierra, ypred.sierra.cv[grep("QUKE", spp),], cell.index, "sp.QUKE")

plot.spdist.raster(pptnorm.sierra, ypred.sierra.mean[grep("PILA", spp),], cell.index, "sp.PILA")

plot.spdist.raster(pptnorm.sierra, ypred.sierra.mean[grep("PIJE", spp),], cell.index, "sp.PIJE")

plot.spdist.raster(pptnorm.sierra, ypred.sierra.mean[grep("JUOC", spp),], cell.index, "sp.JUOC")
plot.spdist.raster(pptnorm.sierra, ypred.sierra.mean[spp=="sp.PIMO3",], cell.index, "sp.PIMO3")


plot.spdist.raster(pptnorm.sierra, ypred.sierra.mean[grep("CADE", spp),], cell.index, "sp.CADE27")
plot.spdist.raster(pptnorm.sierra, ypred.sierra.mean[grep("ABCO", spp),], cell.index, "sp.ABCO")
plot.spdist.raster(pptnorm.sierra, ypred.sierra.mean[grep("PISA", spp),], cell.index, "sp.PISA2")

#### NEED function to calculate skill of species distribution prediction -- since continuous abundance lots of options but try Root mean squared error (RMSE)
# see http://www.sciencedirect.com/science/article/pii/S030438000600247X 

spdist.RMSE <- function(out, sp, spp) {
  return(sqrt(sum((out$modelSummary$yMu[,grep(sp, spp)] - out$y[,grep(sp, spp)])^2)/nrow(out$y)))
}

# use gjam with only 2 species
ml = list(ng=20000, burnin=10000, typeNames="CA")
out_null <- gjamGibbs(~ ppt + ppt2 + tmean + tmean2 +pptxtemp, x, y[,1:2], modelList=ml)

testdata = cbind(y, x)
sp = "ABCO"

# Common species 
spdist.RMSE(out_full, "ABCO", spp) # 11.78
spdist.RMSE(out_null, "ABCO", c("ABCO", "ABMA")) # 11.79
# single-species linear model
m0 = lm(sp.ABCO~ppt+ppt2+tmean+tmean2+pptxtemp, data=testdata)
sqrt(sum((m0$fitted.values-y[,1])^2)/nrow(y)) # 12.56

# Common species 
spdist.RMSE(out_full, "PIPO", spp) # 4.87
out_null <- gjamGibbs(~ ppt + ppt2 + tmean + tmean2 +pptxtemp, x, y[,c(9,2)], modelList=ml)
spdist.RMSE(out_null, "PIPO", c("PIPO", "ABMA")) # 4.84
# single-species linear model
m0 = lm(sp.PIPO~ppt+ppt2+tmean+tmean2+pptxtemp, data=testdata)
sqrt(sum((m0$fitted.values-y[,1])^2)/nrow(y)) # 14.07


# Rarer species
spdist.RMSE(out_full, "JUOC", spp) # 2.39
out_null <- gjamGibbs(~ ppt + ppt2 + tmean + tmean2 +pptxtemp, x, y[,c(3,7)], modelList=ml)
spdist.RMSE(out_null, "JUOC", c("JUOC", "PILA")) #
# single-species linear model
m0 = lm(sp.JUOC~ppt+ppt2+tmean+tmean2+pptxtemp, data=testdata)
sqrt(sum((m0$fitted.values-y[,1])^2)/nrow(y)) # 




###########################
# 3 things left to do: 
# 1) compare multispecies vs single/two species distribution models
# 2) calculate sensitivity of communities to small constant delta temp and ppt (and associated uncertainty)
# 3) calculate projected change in community composition, functional composition (hw), productivity (ba) and diversity (species richness) under future climate (and attribute uncertainty in at least some of these to climate uncertainty and species distribution uncertainty). 




#### Now generate predictions with a perturbation to temperature 
# increase temperature by 1C everywhere
deltatemp = 0.5
deltatemp.scaled =  deltatemp/sd(d$tmean.mean)
xnew.rast = x.rast
xnew.rast$tmean = x.rast$tmean + deltatemp.scaled
xnew.rast$tmean2 = xnew.rast$tmean^2
xnew.rast$pptxtemp = xnew.rast$tmean*xnew.rast$ppt

ypred.deltatemp = gjamPredictY.region(xnew.rast, out=out_full, samples=samples)

n.cells = nrow(x.rast)
ypred.deltatemp.mean = matrix(NA, nsp, n.cells)
for (i in 1:nsp) for (j in 1:n.cells) ypred.deltatemp.mean[i,j] = mean(ypred.deltatemp[i,j,])

# NEXT: CHECK OUT THE SENSITIVITY -- MAPS, SPECIES EFFECTS

# Calculate change in composition -- Jacc
ard dissimilarity
jaccard.sierra.deltatemp = gjamCalcChange.jaccard(ypred.sierra, ypred.deltatemp)
plot.layer.raster(pptnorm.sierra, p=jaccard.sierra.deltatemp[[1]])
plot.layer.raster(pptnorm.sierra, p=jaccard.sierra.deltatemp[[2]])

# Calculate change in composition -- Euclidean distance
euclid.sierra.deltatemp = gjamCalcChange.vegdist(ypred.sierra, ypred.deltatemp, "euclidean")
plot.layer.raster(pptnorm.sierra, p=euclid.sierra.deltatemp[[1]])
plot.layer.raster(pptnorm.sierra, p=euclid.sierra.deltatemp[[2]])

# Calculate change in composition -- summed absolute ba change
abschange.sierra.deltatemp = gjamCalcChange.abschange(ypred.sierra, ypred.deltatemp)
plot.layer.raster(pptnorm.sierra, p=abschange.sierra.deltatemp[[1]])
plot.layer.raster(pptnorm.sierra, p=abschange.sierra.deltatemp[[2]])

# Calculate change in basal area
ba.sierra.deltatemp = gjamCalcChange.basalarea(ypred.sierra, ypred.deltatemp)
plot.layer.raster(pptnorm.sierra, p=ba.sierra.deltatemp[[1]])
plot.layer.raster(pptnorm.sierra, p=ba.sierra.deltatemp[[2]])

# Calculate change in hardwood basal area 
hw.sierra.deltatemp = gjamCalcChange.hardwood(ypred.sierra, ypred.deltatemp, spp)
plot.layer.raster(pptnorm.sierra, p=hw.sierra.deltatemp[[1]])
plot.layer.raster(pptnorm.sierra, p=hw.sierra.deltatemp[[2]])

# Calculate change in species richness
richness.sierra.deltatemp = gjamCalcChange.richness(ypred.sierra, ypred.deltatemp)
plot.layer.raster(pptnorm.sierra, p=richness.sierra.deltatemp[[1]])
plot.layer.raster(pptnorm.sierra, p=richness.sierra.deltatemp[[2]])

# WHich species are overall increasing/declining?

ychange.species.deltatemp = data.frame(species = spp, change=apply(ypred.deltatemp.mean, 1, mean)-apply(ypred.sierra.mean, 1, mean))
ychange.species.deltatemp = ychange.species.deltatemp[order(ychange.species.deltatemp$change),]
barplot(ychange.species.deltatemp$change, names.arg=ychange.species.deltatemp$species, horiz=T, las=1, col="darkgreen", cex.names=0.8)


###### Now evaluate sensitivity to precip
hist(d$ppt.tot)
deltappt = -100 
deltappt.scaled =  deltappt/sd(d$ppt.tot)
xnew.rast = x.rast
xnew.rast$ppt = x.rast$ppt + deltappt.scaled
xnew.rast$ppt2 = xnew.rast$ppt^2
xnew.rast$pptxtemp = xnew.rast$tmean*xnew.rast$ppt

ypred.deltappt = gjamPredictY.region(xnew.rast, out=out_full, samples=samples)

n.cells = nrow(x.rast)
ypred.deltappt.mean = matrix(NA, nsp, n.cells)
for (i in 1:nsp) for (j in 1:n.cells) ypred.deltappt.mean[i,j] = mean(ypred.deltappt[i,j,])

# NEXT: CHECK OUT THE SENSITIVITY -- MAPS, SPECIES EFFECTS

# Calculate change in composition -- Jaccard dissimilarity
jaccard.sierra.deltappt = gjamCalcChange.jaccard(ypred.sierra, ypred.deltappt)
plot.layer.raster(pptnorm.sierra, p=jaccard.sierra.deltappt[[1]])
plot.layer.raster(pptnorm.sierra, p=jaccard.sierra.deltappt[[2]])

# Calculate change in composition -- Euclidean distance
euclid.sierra.deltappt = gjamCalcChange.vegdist(ypred.sierra, ypred.deltappt, "euclidean")
plot.layer.raster(pptnorm.sierra, p=euclid.sierra.deltappt[[1]])
plot.layer.raster(pptnorm.sierra, p=euclid.sierra.deltappt[[2]])

# Calculate change in composition -- summed absolute ba change
abschange.sierra.deltappt = gjamCalcChange.abschange(ypred.sierra, ypred.deltappt)
plot.layer.raster(pptnorm.sierra, p=abschange.sierra.deltappt[[1]])
plot.layer.raster(pptnorm.sierra, p=abschange.sierra.deltappt[[2]])

# compare change in response to temp to change in response to precip
plot(abschange.sierra.deltappt[[1]]~abschange.sierra.deltatemp[[1]])

# Calculate change in basal area
ba.sierra.deltappt = gjamCalcChange.basalarea(ypred.sierra, ypred.deltappt)
plot.layer.raster(pptnorm.sierra, p=ba.sierra.deltappt[[1]])
plot.layer.raster(pptnorm.sierra, p=ba.sierra.deltappt[[2]])

# Calculate change in hardwood basal area 
hw.sierra.deltappt = gjamCalcChange.hardwood(ypred.sierra, ypred.deltappt, spp)
plot.layer.raster(pptnorm.sierra, p=hw.sierra.deltappt[[1]])
plot.layer.raster(pptnorm.sierra, p=hw.sierra.deltappt[[2]])

# Calculate change in species richness
richness.sierra.deltappt = gjamCalcChange.richness(ypred.sierra, ypred.deltappt)
plot.layer.raster(pptnorm.sierra, p=richness.sierra.deltappt[[1]])
plot.layer.raster(pptnorm.sierra, p=richness.sierra.deltappt[[2]])

# WHich species are overall increasing/declining?

ychange.species.deltappt = data.frame(species = spp, change=apply(ypred.deltappt.mean, 1, mean)-apply(ypred.sierra.mean, 1, mean))
ychange.species.deltappt = ychange.species.deltappt[order(ychange.species.deltappt$change),]
barplot(ychange.species.deltappt$change, names.arg=ychange.species.deltappt$species, horiz=T, las=1, col="darkgreen", cex.names=0.8)


##############################################################
# Combining variation in models AND variation in projected climate
#   to quantify and attribute uncertainty in future projections 

# 1) Run future projection for 100 posterior samples x 19 climate models
#     Q do for both temp and precip, or for temperature only also? 

# 2) Calculate community change for each combination in each grid cell
#     To keep it simple, calculate a) Eucl distance, b) basal area change, c) functional change (hardwood basal area) and d) species richness. 

# 3) Display standard deviation of these metrics in each cell. 

# 4) Do an anova or similar to partition variance in each grid cell -- display a map showing proportion variance attributed to model uncert vs climate uncert. 


# outline of steps: 
# get anomalies of T and P for each climate model
# standardize and add to x.rast climate variables data frame to get 19 future scenarios
# loop through the 19 scenarios
# loop through the 100 posterior samples
# save results as 4-D array [species, sites, samples, scenarios]
# For each metric of change: 
#   contrast the 1900 realizations of future community with historical fitted mean ypred.sierra
#   Save the change results as a 3-D array [sites, samples, scenarios]
# For each metric of change: 
#   calculate a total variance for each site across samples and scenarios
#   calculate model uncertainty by estimating variance across scenarios
#   this could be done by loop over sites: creating data frame with sample and scenario factor variables, and doing an ANOVA on it 
#   Map proportion variance associated with model and with climate (temp and precip) and climate (temp only)

# Get anomalies of T and P for each climate model
# historical temp
r = raster("./climatedata/bio_2-5m_bil/bio1.bil", pattern="grd", full.names=TRUE) 
e = extent(c(-122, -117.5, 35, 41))
tmp.hist = crop(r, e)/10
# future temp
path="/users/latimer/googledriveucd/samsi/fia/climatedata"
modelnames = c("ac", "bc", "cc", "ce", "cn", "gd", "gf", "gs", "hd", "hg", "he", "in", "ip", "mi", "mr", "mc", "mp", "mg", "no")
n.models = length(modelnames)
variable = 1
r.list = list()
for (i in 1:length(modelnames)){
    f = paste(path, "/", modelnames[i], "45bi50/", modelnames[i], "45bi50", variable, ".tif", sep="")
    r = raster(f, pattern="grd", full.names=TRUE)
    r.list[i] = crop(r, e)
}
temp.allmodels = brick(r.list)/10
# convert to anomalies
tmp.anom = temp.allmodels - tmp.hist
plot(tmp.anom)
# add to historical data in x.rast from PRISM
tmp.change.resamp = crop(tmp.anom,extent(pptnorm.sierra))
tmp.change.resamp = resample(tmp.change.resamp,pptnorm.sierra)
plot(tmp.change.resamp)
newvals = getValues(tmp.change.resamp)
newvals[is.na(getValues(tmeannorm.sierra))] = NA
tmp.change.sierra = setValues(tmp.change.resamp, newvals)

###  historical Precipitation
r = raster("./climatedata/bio_2-5m_bil/bio12.bil", pattern="grd", full.names=TRUE) 
e = extent(c(-122, -117.5, 35, 41))
ppt.hist = crop(r, e)
# future temp
path="/users/latimer/googledriveucd/samsi/fia/climatedata"
modelnames = c("ac", "bc", "cc", "ce", "cn", "gd", "gf", "gs", "hd", "hg", "he", "in", "ip", "mi", "mr", "mc", "mp", "mg", "no")
n.models = length(modelnames)
variable = 12
r.list = list()
for (i in 1:length(modelnames)){
    f = paste(path, "/", modelnames[i], "45bi50/", modelnames[i], "45bi50", variable, ".tif", sep="")
    r = raster(f, pattern="grd", full.names=TRUE)
    r.list[i] = crop(r, e)
}
ppt.allmodels = brick(r.list)
# convert to anomalies
ppt.anom = ppt.allmodels - ppt.hist
plot(ppt.anom)
# add to historical data in x.rast from PRISM
ppt.change.resamp = crop(ppt.anom,extent(pptnorm.sierra))
ppt.change.resamp = resample(ppt.change.resamp,pptnorm.sierra)
plot(ppt.change.resamp)
newvals = getValues(ppt.change.resamp)
newvals[is.na(getValues(pptnorm.sierra))] = NA
ppt.change.sierra = setValues(ppt.change.resamp, newvals)
plot(ppt.change.sierra)


# create list of future climate scenarios (changing both temp and precip)
x.future.allmodels = list()
for (i in 1:n.models) {
  newvals.tmp = getValues(tmp.change.sierra[[i]])
  newvals.tmp = newvals.tmp[!is.na(newvals.tmp)]
  newvals.tmp.std = newvals.tmp/sd(d$tmean.mean) # standardize the anomalies to same scale as x.rast
  newvals.ppt = getValues(ppt.change.sierra[[i]])
  newvals.ppt = newvals.ppt[!is.na(newvals.ppt)]
  newvals.ppt.std = newvals.ppt/sd(d$ppt.tot) # standardize the anomalies to same scale as x.rast
  xtemp = x.rast
  xtemp$tmean = x.rast$tmean + newvals.tmp.std
  xtemp$tmean2 = (x.rast$tmean + newvals.tmp.std)^2
  xtemp$pptxtemp = (x.rast$tmean + newvals.tmp.std)*(x.rast$ppt + newvals.ppt.std)
  xtemp$ppt = x.rast$ppt + newvals.ppt.std
  xtemp$ppt2 = (x.rast$ppt + newvals.ppt.std)^2
  x.future.allmodels[[i]] = xtemp
}

# create list of future climate scenarios (changing  temp ONLY)
x.future.tmp.allmodels = list()
for (i in 1:n.models) {
  newvals = getValues(tmp.change.sierra[[i]])
  newvals = newvals[!is.na(newvals)]
  newvals.std = newvals/sd(d$tmean.mean) # standardize the anomalies to same scale as x.rast
  xtemp = x.rast
  xtemp$tmean = x.rast$tmean + newvals.std
  xtemp$tmean2 = (x.rast$tmean + newvals.std)^2
  xtemp$pptxtemp = (x.rast$tmean + newvals.std)*x.rast$ppt
  x.future.tmp.allmodels[[i]] = xtemp
}


##### Project future species abundances for all scenarios for all posterior samples

ypred.allmodels = list()
for (i in 1:n.models) {
  print(i)
  ypred.allmodels[[i]] = gjamPredictY.region(x.future.allmodels[[i]], out_full, samples)
}

save(ypred.allmodels, file="ypred.allmodels.Rdata")

# Re-project historical data to make sure using the same posterior samples for direct comparison
ypred.sierra = gjamPredictY.region(x.rast, out_full, samples=samples)


###### Calculate change metrics for all these scenarios x posterior samples

ypred.euclid.allmodels = gjamCalcChange.vegdist.allmodels(ypred.allmodels, ypred.sierra, distmethod="euclidean")

save(ypred.euclid.allmodels, file="ypred.euclid.allmodels.Rdata")

###### Calculate total variance for each grid cell across posterior samples and scenarios
n.cells = nrow(x.rast)
n.samples = length(samples)
n.models = length(modelnames)
totalvar_euclid = rep(NA, n.cells)
for (i in 1:n.cells) totalvar_euclid[i] = var(as.vector(ypred.euclid.allmodels[i,,]), na.rm=T)
plot.layer.raster(pptnorm.sierra, totalvar_euclid)

###### Calculate variance attributable to climate
# Evaluate predictions for all climate models at the posterior mean of the fitted JSDM

postmeanpreds = list()
for (i in 1:19) postmeanpreds[[i]] = gjamPredictY.postmean(x.future.allmodels[[i]], out_full)

n.models = length(postmeanpreds)
ncells = ncol(ypred.sierra.mean)
distmethod = "euclidean"
vegdist.mat= array(NA, list(ncells, n.models))
  for (i in 1:n.models) {
    ypred.future = postmeanpreds[[i]]
    for (j in 1:ncells) {
      vegmat = rbind(ypred.future[,j], ypred.hist[j])
      #if (sum(vegmat[1,])==0 | sum(vegmat[2,])==0) vegdist.mat[j,i] = NA else 
      vegdist.mat[j,i] = vegdist(vegmat, method=distmethod, na.rm=T)
      }
    }

euclid_at_postmean_all_models = vegdist.mat
euclid_at_postmean_all_models.var = apply(vegdist.mat, 1, var, na.rm=T)

plot.layer.raster(pptnorm.sierra, euclid_at_postmean_all_models.var)


# Evaluate predictions for 100 posterior samples at the ensemble mean of climate models
ypred.ensemblemean = gjamPredictY.region(newx.future, out_full, samples)

euclid_at_ensemblemean_across_samples = gjamCalcChange.vegdist(ypred.sierra, ypred.ensemblemean, distmethod = "euclidean")

plot.layer.raster(pptnorm.sierra, euclid_at_ensemblemean_across_samples[[1]])

plot.layer.raster(pptnorm.sierra, euclid_at_ensemblemean_across_samples[[2]]^2)

# Relative importance of climate variance
plot.layer.raster(pptnorm.sierra, euclid_at_postmean_all_models.var / (euclid_at_ense
mblemean_across_samples[[2]]^2 + euclid_at_postmean_all_models.var))


plot.layer.raster.colorfun(pptnorm.sierra, euclid_at_ensemblemean_across_samples[[2]]^2 / (euclid_at_ensemblemean_across_samples[[2]]^2 + euclid_at_postmean_all_models.var), var.colors)

#var.colors = colorRampPalette(c("blue4","green2", "gold", "yellow"))
var.colors = colorRampPalette(c("blue", "lightgray", "red"))

plot.layer.raster.colorfun(pptnorm.sierra, euclid_at_postmean_all_models.var, var.colors)


plot.layer.raster.colorfun(pptnorm.sierra, euclid_at_ensemblemean_across_samples[[2]]^2, var.colors)

min(euclid_at_ensemblemean_across_samples[[2]]^2, na.rm=T)

euclid_at_postmean_all_models.sd = apply(euclid_at_postmean_all_models, 1, sd)

x.rast.tmean.unstd = (x.rast$tmean*sd(d$tmean.mean)+mean(d$tmean.mean))

plot(euclid_at_postmean_all_models.sd~x.rast.tmean.unstd)



#############################################################
### OLD STUFF BELOW HERE

######### Check for relationships between sensitivity and mortality
plot(log(d$mort.tph.nn+0.001)~d$tmean.mean)
plot(log(d$mort.tph.nn+0.001)~d$Defnorm)

plot(log(d$mort.tph.nn+0.001)~sens.vpd)
plot(log(d$mort.tph.nn[d$mort.bin==1]+0.001)~sens.vpd[d$mort.bin==1])

summary(lm(log(d$mort.tph.nn+0.001)~poly(sens.vpd, 2)))
summary(lm(log(d$mort.tph.nn[d$mort.bin==1]+0.001)~poly(sens.vpd[d$mort.bin==1], 2)))
# Not much of anything there. . . 







##############################################
# Using covariance information to improve species abundance predictions using partial data 
# If you know one (or a few species) how well does this improve predictions of the others? 
# How the hell to do this -- essentially imputation of missing values in the multivariate response. Google is no help. 

