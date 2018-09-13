
# Say where the files are and which state to summarize

data_path <- "/users/latimer/google drive/samsi/fia/rawdata/"


states_to_get = c("CA", "NV", "ID", "MT", "OR", "WA", "WY", "CO", "NM", "AZ", "UT")

# loop through and summarize basal area for trees >5 inches in all states separately
for (j in 1:length(states_to_get)) ba.summary(data_path, state=states_to_get[j], public.only=FALSE, annual.only=FALSE)

##############################################################
# Function that summarized BA by species by plot for one state
ba.summary <- function(data_path, state, public.only, annual.only) {
	# Note macroplot breakpoint diameter is different in different states and plots! Can be 21, 24 or 30 inches
	# and is found in the PLOT table. Often there is no macroplot breakpoint which means no macroplot was sampled. 
	
	plots <- read.csv(paste(data_path, state, "/", state, "_PLOT.CSV", sep=""),header=TRUE)
	trees <- read.csv(paste(data_path, state, "/", state, "_TREE.CSV", sep=""),header=TRUE)
	conditions <- read.csv(paste(data_path, state, "/", state, "_COND.CSV", sep=""),header=TRUE)
	ref.species <- read.csv(paste(data_path, "FIADB_REFERENCE/REF_SPECIES.CSV", sep=""),header=TRUE)

	# remove non-forested plots
  # Q should we remove all plots with any nonforest subplot?
  #   If not, this requires accounting for the area of the remaining subplots. 
	plot.cond <- aggregate(COND_STATUS_CD ~ PLT_CN,data=conditions,function(x) {all(x==1)})
	plots <- merge(plots,plot.cond,by.x="CN",by.y="PLT_CN")
	plots <- plots[plots$COND_STATUS_CD == TRUE,]

	# optionally, remove private land plots
	if (public.only) {
		plot.owner <- aggregate(OWNGRPCD ~ PLT_CN,data=conditions,function(x) {all(x != 40)})
		plots <- merge(plots,plot.owner,by.x="CN",by.y="PLT_CN")
		plots <- plots[plots$OWNGRPCD != 40,]
		trees <- trees[trees$plot_id %in% plots$plot_id,]
	}

	# optionally, eliminate all the surveys from before the annual surveys started
	if (annual.only) { 
		annual_survey_start = switch(state, AZ=2001, CA=2001, CO=2002, ID=2004, OR=2001,  MT=2003, NV=2004, NM=2005, UT=2000, WY=2011,  WA=2002)
		plots = plots[which(plots$INVYR>=annual_survey_start),]
		trees <- trees[trees$plot_id %in% plots$plot_id,]
	}
	
	
	# create unique code for each plot in both plots and trees tables
	plots$plot_id = apply(plots[,c("STATECD", "COUNTYCD", "PLOT")], 1, paste, collapse=".")
	trees$plot_id = apply(trees[,c("STATECD", "COUNTYCD", "PLOT")], 1, paste, collapse=".")
	
	# keep only most recent survey measurements
	plotids = unique(plots$plot_id)
	CN_recent = rep(NA, length(plotids))
	for (i in 1:length(plotids)) {
		plotrows = which(plots$plot_id==plotids[i])
		mostrecent = which.max(plots$MEASYEAR[plotrows])
		CN_recent[i] = plots$CN[plotrows[mostrecent]]
	}
	plots = plots[plots$CN %in% CN_recent,]
	trees = trees[trees$PLT_CN %in% CN_recent,]

	# thin to only the focal columns
	focal.cols <- c("CN","PLT_CN","SUBP","TREE","CONDID","STATUSCD","SPCD","SPGRPCD","DIA","HT","ACTUALHT","plot_id")
	trees <- trees[,focal.cols]
	# thin to only focal plots
	trees <- trees[trees$plot_id %in% plots$plot_id,]
	# thin to only live trees
	trees <- trees[trees$STATUSCD == 1,]
	# optionally, thin to only trees >= 5 inches DBH
	trees <- trees[trees$DIA >= 5,]
	#### Calc BA per Ha for each tree ####
	#MACRO_BREAKPOINT_DIA = switch(state, AZ=21, CA=24, CO=21, ID=24, OR=24,  MT=24, NV=21, NM=21, UT=21, WY=21,  WA=24)
  trees$macro_breakpoint = plots$MACRO_BREAKPOINT_DIA[match(trees$plot_id, plots$plot_id)] # assign macroplot diam to each tree
  trees$dia.cm <- trees$DIA * 2.54    # convert inches to cm
	trees$ba.ha <- rep(0,length(trees$dia.cm))

  
  # Calculate BA for the plots that don't have any macroplot
  has_macroplot = !is.na(trees$macro_breakpoint)
	trees$ba.ha[!has_macroplot & trees$dia.cm < 12.7] <- ((trees$dia.cm[!has_macroplot & trees$dia.cm < 12.7]/2)^2) * pi / 0.005398316 # plot expansion factor for smaller trees
	trees$ba.ha[!has_macroplot & trees$dia.cm >= 12.7] <- ((trees$dia.cm[!has_macroplot & trees$dia.cm >= 12.7]/2)^2) * pi / 0.06724546 # plot expansion factor for larger trees
	
  # Same for plots that DO have a macroplot  
  trees$ba.ha[has_macroplot & trees$dia.cm < 12.7] <- ((trees$dia.cm[has_macroplot & trees$dia.cm < 12.7]/2)^2) * pi / 0.005398316 # plot expansion factor
	trees$ba.ha[has_macroplot & trees$dia.cm >= 12.7 & trees$dia.cm <= trees$macro_breakpoint] <- ((trees$dia.cm[has_macroplot & trees$dia.cm >= 12.7 & trees$dia.cm <= trees$macro_breakpoint]/2)^2) * pi / 0.06724546 # plot expansion factor
	trees$ba.ha[has_macroplot & trees$dia.cm > trees$macro_breakpoint]  <- ((trees$dia.cm[has_macroplot & trees$dia.cm > trees$macro_breakpoint]/2)^2) * pi / 0.4050149   # plot expansion factor

  
  #### Sum the BA per Ha by plot and species ####
	trees$ba.ha <- trees$ba.ha/100^2 # convert to m2/ha	
  plot.ba.spcs.all <- aggregate(ba.ha ~ SPCD * plot_id, data=trees, sum)
	plot.ba.spcs.all <- merge(plot.ba.spcs.all,ref.species[,c("SPCD","SPECIES_SYMBOL")],by="SPCD",all.x=TRUE)
	plot.ba.spcs.all$SPECIES_SYMBOL <- droplevels(plot.ba.spcs.all$SPECIES_SYMBOL)
  plot.ba.spcs.all$PLT_CN = plots$CN[match(plot.ba.spcs.all$plot_id, plots$plot_id)] # append the FIA plot code number

	write.csv(plot.ba.spcs.all, paste("/users/latimer/google drive/samsi/fia/workingdata/", state, "_ba_summary_long.csv", sep=""))
}
#########################################



###########################
# From here on, this is to check and explore the data collected above. 
# Then combine into a large basal area matrix for all states together. 

plot.ba.spcs.all_AZ = read.csv("./workingdata/AZ_ba_summary_long.csv")
plot.ba.spcs.all_CA = read.csv("./workingdata/CA_ba_summary_long.csv")
plot.ba.spcs.all_CO = read.csv("./workingdata/CO_ba_summary_long.csv")
plot.ba.spcs.all_ID = read.csv("./workingdata/ID_ba_summary_long.csv")
plot.ba.spcs.all_MT = read.csv("./workingdata/MT_ba_summary_long.csv")
plot.ba.spcs.all_NM = read.csv("./workingdata/NM_ba_summary_long.csv")
plot.ba.spcs.all_NV = read.csv("./workingdata/NV_ba_summary_long.csv")
plot.ba.spcs.all_OR = read.csv("./workingdata/OR_ba_summary_long.csv")
plot.ba.spcs.all_UT = read.csv("./workingdata/UT_ba_summary_long.csv")
plot.ba.spcs.all_WA = read.csv("./workingdata/WA_ba_summary_long.csv")
plot.ba.spcs.all_WY = read.csv("./workingdata/WY_ba_summary_long.csv")

# check if overlap in PLT_CN numbers for different states -- apparently they're unique
z = unlist(sapply(plot.ba.spcs.all_OR$PLT_CN, FUN=f<-function(x) {return(x %in% plot.ba.spcs.all_WA$PLT_CN)}))
sum(z)

ba.all = rbind(plot.ba.spcs.all_AZ, plot.ba.spcs.all_CA, plot.ba.spcs.all_CO , plot.ba.spcs.all_ID, plot.ba.spcs.all_MT, plot.ba.spcs.all_NM, plot.ba.spcs.all_NV, plot.ba.spcs.all_OR, plot.ba.spcs.all_UT, plot.ba.spcs.all_WA, plot.ba.spcs.all_WY)


write.csv(ba.all, "/users/latimer/google drive/samsi/fia/workingdata/basal_area_western_states_long.csv")


#### Convert to wide format and save ####
ba.all.wide <- xtabs(ba.ha ~ PLT_CN + SPECIES_SYMBOL, data=ba.all.thinned)
write.csv(ba.all.wide, "/users/latimer/google drive/samsi/fia/workingdata/basal_area_western_states_wide.csv")






#########################
# Explore the data 


# decide which species to keep then thin out the rest
#### Get mean BA per Ha by plot across all species ####
ba.total <- aggregate(ba.ha ~ SPECIES_SYMBOL, data=ba.all, mean)
plot(rev(sort(ba.total$ba.ha)), col=1)
ba.total$SPECIES_SYMBOL[order(ba.total$ba.ha, decreasing=T)]

#### Get number of plots occupied across all species ####
prev.total <- aggregate(ba.ha ~ SPECIES_SYMBOL, data=ba.all, f<-function(x){return(sum(x>0))})
plot(rev(sort(prev.total$ba.ha)), col=1)
prev.total$SPECIES_SYMBOL[order(prev.total$ba.ha, decreasing=T)]
prev.total$ba.ha[order(prev.total$ba.ha, decreasing=T)]

# for now, remove all species that occur in fewer than 100 plots
species.to.keep = prev.total$SPECIES_SYMBOL[prev.total$ba.ha>= 100]

ba.all.thinned = ba.all[ba.all$SPECIES_SYMBOL %in% species.to.keep,]
ba.all.thinned$SPECIES_SYMBOL = as.character(ba.all.thinned$SPECIES_SYMBOL)

#### Convert to wide format ####
ba.all.wide <- xtabs(ba.ha ~ PLT_CN + SPECIES_SYMBOL, data=ba.all.thinned)

head(ba.all.wide)
dim(ba.all.wide)

# quick look at species abundances
abund=rev(sort(apply(ba.all.wide, 2, sum)))
plot(abund, col=0)
text(abund, names(abund), cex=0.75)
prev=rev(sort(apply(ba.all.wide, 2, f<-function(x){return(sum(x>0))})))
plot(prev, col=0)
text(prev, names(prev), cex=0.75)
plot(log(prev)~log(abund))

