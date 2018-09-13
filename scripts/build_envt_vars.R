# Build and check environmental variables for FIA plots
# Andrew Latimer
# April 24, 2015


plots = read.csv("/users/latimer/google drive/samsi/fia/workingdata/plotdata_western_states_PRISM.csv")



# A few checks on the plot-level data to see factors affecting them

# Planting -- overwhelmingly a PNW phenomenon 
table(plots$STDORGCD) # >1200 with clear evid of artificial regen
palette("default")
plot(LAT~LON, plots, col=plots$STDORGCD+1, pch=16,cex=0.4)
plot(LAT~LON, plots[plots$STDORGCD>0,], pch=16,cex=0.4)

# Cutting -- everywhere, though seems more evident in MT, AZ
table(plots$TRTCD1) # ~1000 with evidence of cutting
plot(LAT~LON, plots, col=as.integer(plots$TRTCD1==10)+1, pch=16, cex=0.4)
plot(LAT~LON, plots[plots$TRTCD1==10,], pch=16, cex=0.4)

# Natural regen "treatment" -- almost never noted -- not sure what that means, apparently nothing. 
plot(LAT~LON, plots, col=as.integer(plots$TRTCD1==40)+1, pch=16, cex=0.4)
plot(LAT~LON, plots[plots$TRTCD1==40,], pch=16, cex=0.4)

# Where are plots not according to standard design?
plot(LAT~LON, plots, col=as.integer(plots$DESIGNCD>1)+1, pch=16, cex=0.4)
plot(LAT~LON, plots[plots$DESIGNCD>1,], pch=16, cex=0.4)

# Q SHOULD WE TOSS THESE? THEY"RE All IN THE 400's (RMRS) and 500's (PNW)
# ??????
# what years? 
table(plots$DESIGNCD)
table(plots$INVYR.x[plots$DESIGNCD==410])
table(plots$INVYR.x[plots$DESIGNCD==552])
sum(plots$INVYR.x>=2000)
sum(plots$INVYR.x<2000)

table(plotdata.all$OWNGRPCD[plotdata.all$TRTCD1==10])
table(plotdata.all$TRTCD1[plotdata.all$OWNGRPCD==40]) # about 40% of private plots show evid of cutting
table(plotdata.all$TRTCD1[plotdata.all$OWNGRPCD==10]) # about 15% of USFS plots "" ""
table(plotdata.all$DSTRBCD1) # maybe 10% with fire evid
table(plots$CHAINING_CD) # >20 plots with evid of chaining, most are NA so I assume no chaining
table(plots$DESIGNCD) 
table(plotdata.all$DESIGNCD[plotdata.all$MANUAL>=2])

sum(is.na(plots$LON))
plots = plots[!is.na(plots$LON),]

plots.map = SpatialPointsDataFrame(plots[,c("LON", "LAT")], data=plots)

# more specifically check if survey type is the current, annual method
palette("default")
plot(plots.map, col=(as.integer(plots.map$MANUAL>=2)+1), pch=15, cex=0.2) # WY is hugely out of date, to a lesser extent all the Rockies states
# and further, regardless of when sampled, whether it uses the current national plot design
points(plots.map[plots.map$MANUAL<2 & plots.map$DESIGNCD==1,], col="blue",pch=15, cex=0.2) # This fortunately rescues a lot of WY, AZ, NM, ID and CO plots 


#########################################################################
# TRIM DATA

#### To start with, use fairly stringent inclusion criteria for analysis
#### Screen out plots that have evidence of artificial planting 
#### Remove plots that are not done according to the standard protocol (DESIGNCD == 1)

# first remove all the  plots that are not using current plot design. 
# this makes things easier but ultimately we may need to include some of these to get decent coverage in the Rockies
plots = plots[which(plots$INVYR.x >= 2000),]
plots = plots[which(plots$DESIGNCD==1),]
# Gives us ~29,000 plots to work with

# Remove plots with artificial regen
plots = plots[which(plots$STDORGCD == 0),]

# Remove plots with evidence of cutting or other forestry disturbance (codes 10, 20, 30, 50)
cutplots = which(plots$TRTCD1 %in% c(10,20,30,50) | plots$TRTCD2 %in% c(10,20,30,50) | plots$TRTCD3 %in% c(10,20,30,50))
plots = plots[-cutplots,]

# Remove plots on private land if not already removed
pvtplots = which(plots$OWNGRPCD == 40)
if (length(pvtplots)>0) plots = plots[-pvtplots,]

# Optionally, remove plots with permanent water on them 
# NOTE this could be revisited later or turned into a covariate; not doing that now because the number of plots involved is relatively small. 
# One alternative is to retain code 2 (small permanent water body) as a "riparian" dummy variable or something like that. 
#table(plots$PHYSCLCD[plots$WATERCD==2])
#waterplots = which(plots$WATERCD %in% c(1, 2, 3))
#plots = plots[-waterplots,]
#dim(plots)



######################################
#### ASSEMBLE COVARIATES AT PLOT LEVEL

# Extract the plot-level covariates other than physiographic class codes and disturbance codes since those have to be condensed separately. 
vars_to_include = c("CN.x", "PLT_CN", "plot_id.x", "STATECD.x", "ECOSUBCD", "MACRO_BREAKPOINT_DIA", "INVYR.x", "MEASYEAR", "MEASMON", "MEASDAY", "LON", "LAT", "ELEV", "SLOPE", "ASPECT",  "ppt", "tmax", "tmin", "vpdmax")

envtdata = plots[,vars_to_include]

# Create disturbance dummy variables
fire = which(plots$DSTRBCD1 %in% c(30, 31, 32) | plots$DSTRBCD2 %in% c(30, 31, 32) | plots$DSTRBCD3 %in% c(30, 31, 32))
envtdata$fire = rep(0, nrow(envtdata))
envtdata$fire[fire] = 1

crownfire = which(plots$DSTRBCD1 == 32 | plots$DSTRBCD2 == 32 | plots$DSTRBCD3 == 32) # only 310 crown fire plots !!!!
envtdata$crownfire = rep(0, nrow(envtdata))
envtdata$crownfire[crownfire] = 1

insect = which(plots$DSTRBCD1 %in% c(10, 12) | plots$DSTRBCD2 %in% c(10, 12) | plots$DSTRBCD3 %in% c(10, 12))
envtdata$insect = rep(0, nrow(envtdata))
envtdata$insect[insect] = 1

disease = which(plots$DSTRBCD1 %in% c(20, 22) | plots$DSTRBCD2 %in% c(20, 22) | plots$DSTRBCD3 %in% c(20, 22))
envtdata$disease = rep(0, nrow(envtdata))
envtdata$disease[disease] = 1

# Create physiographic class variable that simplifies into fewer groups
# 
envtdata$xeric = envtdata$mesic = envtdata$hydric = envtdata$dryslope = envtdata$otherdry = envtdata$moistslope = envtdata$othermesic = rep(0, nrow(envtdata))
envtdata$xeric[plots$PHYSCLCD %in% c(11, 12, 13, 19)] = 1
envtdata$mesic[plots$PHYSCLCD %in% c(21, 22, 23, 24, 25, 29)] = 1
envtdata$hydric[plots$PHYSCLCD %in% c(31, 32, 33, 34, 35, 39)] = 1
envtdata$dryslope[plots$PHYSCLCD == 12] = 1
envtdata$otherdry[plots$PHYSCLCD %in% c(11, 13, 19)] = 1
envtdata$othermesic[plots$PHYSCLCD %in% c(21, 22, 24, 25, 29)] = 1
envtdata$moistslope[plots$PHYSCLCD == 23] = 1


# Check for state discontinuities in physiographic classification
plot(LON~LAT, envtdata, pch=16, cex=0.4)
points(LON~LAT, envtdata[envtdata$xeric==1,], pch=16, cex=0.4, col="red")
# Bizarrely, the whole of eastern WA and OR are "xeric"" ?? 
points(LON~LAT, envtdata[envtdata$mesic==1,], pch=16, cex=0.4, col="green")
# Plus, WY and MT are in some cases unclassified (NA values)
points(LON~LAT, envtdata[envtdata$hydric==1,], pch=16, cex=0.4, col="blue")

# Put these on Google Earth to check
# Clearly many have been "fuzzed" over small ridgelines so they're classified as xeric but on DEM face north. 
# Cursory look at fire plots in CA and WA suggests they are almost all within real fires.
library(maptools)
xericpoints = SpatialPointsDataFrame(envtdata[envtdata$xeric==1,c("LON", "LAT")], data=envtdata[envtdata$xeric==1,c("CN.x", "LON", "LAT")])
kmlPoints(xericpoints, kmlfile="xericpoints.kml", kmlname="xeric.fia.plots", kmldescription="Locations of FIA plots classified as xeric in PHYSCLCD field",
          icon="http://google.com/mapfiles/kml/paddle/red-diamond.png")
mesicpoints = SpatialPointsDataFrame(envtdata[envtdata$mesic==1,c("LON", "LAT")], data=envtdata[envtdata$mesic==1,c("CN.x", "LON", "LAT")])
kmlPoints(mesicpoints, kmlfile="mesicpoints.kml", kmlname="xeric.fia.plots", kmldescription="Locations of FIA plots classified as mesic in PHYSCLCD field",
          icon="http://google.com/mapfiles/kml/paddle/blu-diamond.png")
hydricpoints = SpatialPointsDataFrame(envtdata[envtdata$hydric==1,c("LON", "LAT")], data=envtdata[envtdata$hydric==1,c("CN.x", "LON", "LAT")])
kmlPoints(hydricpoints, kmlfile="hydricpoints.kml", kmlname="hydric.fia.plots", kmldescription="Locations of FIA plots classified as hydric in PHYSCLCD field",
          icon="http://google.com/mapfiles/kml/paddle/wht-diamond.png")

firepoints = SpatialPointsDataFrame(envtdata[envtdata$fire==1,c("LON", "LAT")], data=envtdata[envtdata$fire==1,c("CN.x", "LON", "LAT")])
kmlPoints(firepoints, kmlfile="firepoints.kml", kmlname="fire.fia.plots", kmldescription="Locations of FIA plots classified as disturbed by fire in DSTRBCD field",
          icon="http://google.com/mapfiles/kml/paddle/wht-diamond.png")




head(envtdata)
dim(envtdata)

write.csv(envtdata, "./workingdata/envtdata_western_states.csv")




###################################

# FIRE DISTURBANCE LOCATIONS
plot(plots.map, col=plotdata.all$DSTRBCD1[z]==30, pch=15, cex=0.2)
plot(plots.map, col=plotdata.all$DSTRBCD2[z]==30, pch=15, cex=0.2, add=FALSE) # Weirdly fire seems to appear as DSTRBCD2 only in the 

fireplots = plots.map$DSTRBCD1 %in% c(30, 31, 32) | plots.map$DSTRBCD2 %in% c(30, 31, 32) | plots.map$DSTRBCD3 %in% c(30, 31, 32)
genericfireplots = which(plots.map$DSTRBCD1 == 30 | plots.map$DSTRBCD2 == 30 | plots.map$DSTRBCD3 == 30)
groundfireplots = which(plots.map$DSTRBCD1 == 31 | plots.map$DSTRBCD2 == 31 | plots.map$DSTRBCD3 == 31)
crownfireplots = which(plots.map$DSTRBCD1 ==32 | plots.map$DSTRBCD2 ==32 | plots.map$DSTRBCD3 ==32)

palette(terrain.colors(120))
plot(plots.map, col=plots.map$ELEV/100, pch=15, cex=0.4)
#points(plots.map[fireplots,], col="red", pch=17, cex=0.2)
points(plots.map[groundfireplots,], col="blue", pch=17, cex=0.3)
points(plots.map[crownfireplots,], col="red", pch=17, cex=0.3)
points(plots.map[genericfireplots,], col="black", pch=17, cex=0.3)

# Timing of fire
# Note not exact because sometimes fire is disturbance code 2 or 3
hist(plots.map$DSTRBYR1[genericfireplots])
hist(plots.map$DSTRBYR1[groundfireplots])
hist(plots.map$DSTRBYR1[crownfireplots])

# Get year of fire disturbance


# Time between fire and sampling
hist(plots.map$INVYR.x[genericfireplots] - plots.map$DSTRBYR1[genericfireplots])
hist(plots.map$INVYR.x[groundfireplots] - plots.map$DSTRBYR1[groundfireplots])
hist(plots.map$INVYR.x[crownfireplots] - plots.map$DSTRBYR1[crownfireplots])


# Does fire affect presence of seedlings?
plotstemdata = merge(plots, stemcounts.all, by.x="CN.x", by.y="PLT_CN")
head(plotstemdata)


plot(log(seedlingcount+1)~DSTRBYR1, plotstemdata[crownfireplots,]

boxplot(log(seedlingcount+1)~DSTRBCD1, plotstemdata[plotstemdata$SPECIES_SYMBOL=="ABCO",])

boxplot(log(treecount_1_5+1)~DSTRBCD1, plotstemdata[plotstemdata$SPECIES_SYMBOL=="ABCO",])
boxplot(log(treecount_5_MPBPD+1)~DSTRBCD1, plotstemdata[plotstemdata$SPECIES_SYMBOL=="ABCO",], notch=T)
boxplot(log(treecount_5_MPBPD+1)~DSTRBCD1, plotstemdata[plotstemdata$SPECIES_SYMBOL=="PIJE",], notch=T)
boxplot(log(treecount_5_MPBPD+1)~DSTRBCD1, plotstemdata[plotstemdata$SPECIES_SYMBOL=="PIPO",], notch=T)


summary(glm(seedlingcount~factor(DSTRBCD1), plotstemdata[plotstemdata$SPECIES_SYMBOL=="ABCO",], family="poisson"))
summary(glm(seedlingcount~factor(DSTRBCD1), plotstemdata[plotstemdata$SPECIES_SYMBOL=="PIPO",], family="poisson"))
summary(glm(seedlingcount~factor(DSTRBCD1), plotstemdata[plotstemdata$SPECIES_SYMBOL=="PICO",], family="poisson"))

summary(glm(treecount_1_5~factor(DSTRBCD1), plotstemdata[plotstemdata$SPECIES_SYMBOL=="ABCO",], family="poisson"))
summary(glm(treecount_1_5~factor(DSTRBCD1), plotstemdata[plotstemdata$SPECIES_SYMBOL=="PIPO",], family="poisson"))
summary(glm(treecount_1_5~factor(DSTRBCD1), plotstemdata[plotstemdata$SPECIES_SYMBOL=="PICO",], family="poisson"))

summary(glm(treecount_5_MPBPD~factor(DSTRBCD1), plotstemdata[plotstemdata$SPECIES_SYMBOL=="ABCO",], family="poisson"))
summary(glm(treecount_5_MPBPD~factor(DSTRBCD1), plotstemdata[plotstemdata$SPECIES_SYMBOL=="PIPO",], family="poisson"))
summary(glm(treecount_5_MPBPD~factor(DSTRBCD1), plotstemdata[plotstemdata$SPECIES_SYMBOL=="PICO",], family="poisson"))
# Species differences clearly seen in adult trees. Q whether this just indicates environmental association with fire. 

# Add climate covariates
summary(glm(treecount_5_MPBPD~ppt+tmax+vpdmax +factor(DSTRBCD1), plotstemdata[plotstemdata$SPECIES_SYMBOL=="ABCO",], family="poisson"))
summary(glm(treecount_5_MPBPD~ppt+tmax+vpdmax +factor(DSTRBCD1), plotstemdata[plotstemdata$SPECIES_SYMBOL=="PIPO",], family="poisson"))
summary(glm(treecount_5_MPBPD~ppt+tmax+vpdmax +factor(DSTRBCD1), plotstemdata[plotstemdata$SPECIES_SYMBOL=="PICO",], family="poisson"))
summary(glm(treecount_5_MPBPD~ppt+tmax+vpdmax +factor(DSTRBCD1), plotstemdata[plotstemdata$SPECIES_SYMBOL=="PSME",], family="poisson"))

plot(glm(treecount_5_MPBPD~ppt+tmax+vpdmax +factor(DSTRBCD1), plotstemdata[plotstemdata$SPECIES_SYMBOL=="PSME",], family="poisson"))


#points(plots.map[which(plots.map$INVYR.x==2013),], col="blue", pch=17, cex=0.3)



par(mfrow=c(5,5), mar=rep(0.2, 4))
for (i in sort(unique(plots.map$INVYR.x))) {
  plot(plots.map[which(plots.map$INVYR.x==i),], col="darkgreen", pch=15, cex=0.3, main=i, cex.main=0.8,
       xlim=range(plots.map$LON), ylim=range(plots.map$LAT))
  map(database="usa",xlim=range(plots.map$LON), ylim=range(plots.map$LAT), col="black", add=T)
}


#plot(plots.map[plots.map$INVYR.x == 1999,], col=plots.map$DESIGNCD[plots.map$INVYR.x == 1999])


# OTHER DISTURBANCE LOCATIONS
fireplots = plots.map$DSTRBCD1 %in% c(30, 31, 32) | plots.map$DSTRBCD2 %in% c(30, 31, 32) | plots.map$DSTRBCD3 %in% c(30, 31, 32)
logplots = plots.map$TRTCD1 %in% c(10) | plots.map$TRTCD2 %in% c(10) | plots.map$TRTCD3 %in% c(10)
insectplots = plots.map$DSTRBCD1 %in% c(10, 11, 12) | plots.map$DSTRBCD2 %in% c(10, 11, 12) | plots.map$DSTRBCD3 %in% c(10, 11, 12)
diseaseplots = plots.map$DSTRBCD1 %in% c(20, 21, 22) | plots.map$DSTRBCD2 %in% c(20, 21, 22) | plots.map$DSTRBCD3 %in% c(20, 21, 22)
grayscale <- colorRampPalette(c("lightgray", "darkgray", "black"))
palette(grayscale(120))
plot(plots.map, col=plots.map$ELEV/100, pch=15, cex=0.4, main="Disturbance records in the last 10 years of FIA data")
#points(plots.map[diseaseplots,], col="violet", pch=17, cex=0.4)
points(plots.map[logplots,], col="blue", pch=17, cex=0.4)
points(plots.map[insectplots,], col="yellow", pch=17, cex=0.4)
points(plots.map[fireplots,], col="red", pch=17, cex=0.4)
map(database="usa",xlim=range(plots.map$LON), ylim=range(plots.map$LAT), col="darkgray",lwd=2, add=T)
legend("bottomleft", c("Logging", "Insects", "Fire"), pch=rep(17, 3), col=c("blue", "yellow","red"))
