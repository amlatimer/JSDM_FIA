# Extract climate data from PRISM for the western FIA forest plots. 
# Uses plot locations obtained using "get_FIA_plot_data.R"
#   and saved to the file "
# Andrew Latimer
# April 24, 2015

library(raster)

# names including any subdirectory of files containing environmental data layers
layernames = c("PRISM_ppt_30yr_normal_800mM2_annual_bil/PRISM_ppt_30yr_normal_800mM2_annual_bil.bil", 
              "PRISM_tmax_30yr_normal_800mM2_07_bil/PRISM_tmax_30yr_normal_800mM2_07_bil.bil",
              "PRISM_tmin_30yr_normal_800mM2_01_bil/PRISM_tmin_30yr_normal_800mM2_01_bil.bil",
              "PRISM_vpdmax_30yr_normal_800mM2_annual_bil/PRISM_vpdmax_30yr_normal_800mM2_annual_bil.bil")
# names of those layers
envtvars = c("ppt", "tmax", "tmin", "vpdmax")

# Paths to environmental and plots data
PRISM_path = "~/documents/CAclimate/PRISM/"
plots_path = "~/google drive/SAMSI/FIA/workingdata/"
plots_filename = "plotdata_western_states.csv"

# Read in and clean the plots data
plots = read.csv(paste(plots_path, plots_filename, sep=""))
# remove plots that are on private land because they are more "fuzzed" and some are swapped. 
plots <- plots[plots$OWNGRPCD != 40,]
# remove plots with no lat/lon coordinates
plots = plots[which(!is.na(plots$LON) & !is.na(plots$LAT)),]

# Convert location info to spatial points object
plotlocs = SpatialPoints(plots[,c("LON", "LAT")])

# Optionally look at the data
#plot(r)
#points(plotlocs, pch=17, cex=0.5)

# Loop through the data layers and extract the data to the data frame envtdata
envtdata = data.frame(plot_id=plots$plot_id, PLT_CN=plots$CN.x)

for (i in 1:length(layernames)) {
  PRISM_filename = "PRISM_ppt_30yr_normal_800mM2_annual_bil/PRISM_ppt_30yr_normal_800mM2_annual_bil.bil"
  r = raster(paste(PRISM_path, layernames[i], sep=""))
  # Extract environmental variable
  envtdata[,2+i] = extract(r, plotlocs)
  names(envtdata)[2+i] = envtvars[i]
}

plots = merge(plots, envtdata, by.x="CN.x", by.y="PLT_CN")
head(plots)

# Write out the product
write.csv(plots, "/users/latimer/google drive/samsi/fia/workingdata/plotdata_western_states_PRISM.csv")


# quick look
palette(heat.colors(45))
plot(LAT~LON, plots, col=plots$vpd, pch=16, cex=0.5)


