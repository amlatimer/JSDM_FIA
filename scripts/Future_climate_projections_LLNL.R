######## Make layers of projected climate change for Sierras
# Using CMIP5 downscaled by LLNL

# Method: Calculate long-term annual means for 1970-1999 and for 2040-2069, subtract to get change. Add change to the PRISM long-term normals to get future projection. 

# Start with single model x projection: GFDL RCP 4.5 

setwd("/users/latimer/googledriveucd/samsi/fia")

# create a raster to hold the period mean values. 


colvals = scan("./Climatedata/bcsd5/COLS_PeriodStat.txt")
rowvals = scan("./Climatedata/bcsd5/ROWS_PeriodStat.txt")


#r = raster(nrows=length(rowvals), ncols=length(colvals), ymn=min(rowvals), ymx=max(rowvals), xmn=min(colvals), xmx=max(colvals), crs='+proj=longlat +datum=WGS84')
#r = raster(xmn=237.875, xmx=243.375, ymn=34.5, ymx=41, res=1/8, crs='+proj=longlat +datum=WGS84')
r = raster(xmn=237.75, xmx=243.5, ymn=34.375, ymx=41.125, res=1/8, crs='+proj=longlat +datum=WGS84')

rr <- rotate(r)

#ppt.future = read.csv("./Climatedata/bcsd5_2040-69_gfdl_rcp45/pr_PeriodStat_mean.gfdl-cm3.1.rcp45.csv", header=FALSE)
ppt.future = read.csv("./Climatedata/bcsd5_2040-2069_gfdl_rcp60/pr_PeriodStat_mean.gfdl-cm3.1.rcp60.csv", header=FALSE)
ppt.future[ppt.future=="--"] = NA
ppt.future.mat = data.matrix(ppt.future)
image(ppt.future.mat)
ppt.future.raster = flip(raster(ppt.future.mat, template=rr), direction="y")
plot(ppt.future.raster)

tmp.future = read.csv("./Climatedata/bcsd5_2040-69_gfdl_rcp45/tas_PeriodStat_mean.gfdl-cm3.1.rcp45.csv", header=FALSE)
tmp.future.mat = data.matrix(tmp.future)
image(tmp.future.mat)
tmp.future.raster = flip(raster(tmp.future.mat, template=rr), direction="y")
plot(tmp.future.raster)
hist(tmp.future.raster)



ppt.hist = read.csv("./Climatedata/bcsd5_1970-1999/pr_PeriodStat_mean.gfdl-cm3.1.rcp45.csv", header=FALSE)
ppt.hist[ppt.hist=="--"] = NA
#ppt.hist = read.csv("./Climatedata/bcsd5_1970-1999/pr_PeriodStat_mean.fgoals-g2.1.rcp45.csv", header=FALSE)

ppt.hist.mat = data.matrix(ppt.hist)
image(ppt.hist.mat)
ppt.hist.raster = flip(raster(ppt.hist.mat, template=rr), direction="y")
plot(ppt.hist.raster)

ppt.change = raster(r)
ppt.change = setValues(ppt.change, ppt.future.mat-ppt.hist.mat)
plot(ppt.change)
 #Somethign corrupted here. 

image(ppt.future.mat-ppt.hist.mat) # Here too -- so it's not the raster that's the issue
image(ppt.hist.mat)
head(ppt.future)
dim(ppt.future)

