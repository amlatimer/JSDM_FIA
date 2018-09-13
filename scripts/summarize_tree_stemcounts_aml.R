setwd("~/google drive/samsi/fia")

#subsec = read.dbf("./S_USA.EcomapSubsections/S_USA.EcomapSubsections.dbf")

# Say where the files are and which state to summarize

data_path <- "/users/latimer/google drive/samsi/fia/rawdata/"

states_to_get = c("CA", "NV", "ID", "MT", "OR", "WA", "WY", "CO", "NM", "AZ", "UT")

# loop through and summarize basal area for trees >5 inches in all states separately
for (j in 1:length(states_to_get)) stemcount.summary(data_path, state=states_to_get[j])


stemcount.summary <- function(data_path, state) {

	plots <- read.csv(paste(data_path, state, "/", state, "_PLOT.CSV", sep=""),header=TRUE)
	trees <- read.csv(paste(data_path, state, "/", state, "_TREE.CSV", sep=""),header=TRUE)
	conditions <- read.csv(paste(data_path, state, "/", state, "_COND.CSV", sep=""),header=TRUE)
	ref.species <- read.csv(paste(data_path, "FIADB_REFERENCE/REF_SPECIES.CSV", sep=""),header=TRUE)

	# keep only most recent survey measurements
	plotids = unique(plots$PLOT)
	CN_recent = rep(NA, length(plotids))
	for (i in 1:length(plotids)) {
		plotrows = which(plots$PLOT==plotids[i])
		mostrecent = which.max(plots$MEASYEAR[plotrows])
		CN_recent[i] = plots$CN[plotrows[mostrecent]]
	}
	plots = plots[plots$CN %in% CN_recent,]
	trees = trees[trees$PLT_CN %in% CN_recent,]
	
	# keep only "forested" plots
	plot.cond <- aggregate(COND_STATUS_CD ~ PLT_CN,data=conditions,function(x) {all(x==1)})
	plots <- merge(plots,plot.cond,by.x="CN",by.y="PLT_CN")
	plots <- plots[plots$COND_STATUS_CD == TRUE,]

	# thin to only the focal columns
	focal.cols <- c("CN","PLT_CN","SUBP","TREE","CONDID","STATUSCD","SPCD","SPGRPCD","DIA","HT","ACTUALHT")
	trees <- trees[,focal.cols]
	# thin to only focal plots
	trees <- trees[trees$PLT_CN %in% plots$CN,]
	# thin to only live trees
	trees <- trees[trees$STATUSCD == 1,]
	
	#### NOTE FROME HERE ON NEEDS TO BE UPDATED TO GET STEM INFO 
	# Here I believe we just need the counts because the area adjustment should happen in the statistical model. 
		
	# thin to only trees >= 5 inches DBH
	trees <- trees[trees$DIA >= 5,]
	#### Count stems per Ha for each plots and species ####
	plot.stems.spcs.all <- aggregate(DIA ~ SPCD * PLT_CN, data=trees, length)
	plot.stems.spcs.all <- merge(plot.stems.spcs.all,ref.species[,c("SPCD","SPECIES_SYMBOL")],by="SPCD",all.x=TRUE)
	plot.stems.spcs.all$SPECIES_SYMBOL <- droplevels(plot.stems.spcs.all$SPECIES_SYMBOL)

	write.csv(plot.ba.spcs.all, paste("/users/latimer/google drive/samsi/fia/workingdata/", state, "_stemcount_summary_long.csv", sep=""))
}





###########################
### Read in and combine all the data sets, filter them, then convert to wide format. 

plot.stemcount.spcs.all_AZ = read.csv("./workingdata/AZ_stemcount_summary_long.csv")
plot.stemcount.spcs.all_CA = read.csv("./workingdata/CA_stemcount_summary_long.csv")
plot.stemcount.spcs.all_CO = read.csv("./workingdata/CO_stemcount_summary_long.csv")
plot.stemcount.spcs.all_ID = read.csv("./workingdata/ID_stemcount_summary_long.csv")
plot.stemcount.spcs.all_MT = read.csv("./workingdata/MT_stemcount_summary_long.csv")
plot.stemcount.spcs.all_NM = read.csv("./workingdata/NM_stemcount_summary_long.csv")
plot.stemcount.spcs.all_NV = read.csv("./workingdata/NV_stemcount_summary_long.csv")
plot.stemcount.spcs.all_OR = read.csv("./workingdata/OR_stemcount_summary_long.csv")
plot.stemcount.spcs.all_UT = read.csv("./workingdata/UT_stemcount_summary_long.csv")
plot.stemcount.spcs.all_WA = read.csv("./workingdata/WA_stemcount_summary_long.csv")
plot.stemcount.spcs.all_WY = read.csv("./workingdata/WY_stemcount_summary_long.csv")

# check if overlap in PLT_CN numbers for different states -- apparently they're unique
z = unlist(sapply(plot.stemcount.spcs.all_OR$PLT_CN, FUN=f<-function(x) {return(x %in% plot.stemcount.spcs.all_WA$PLT_CN)}))
sum(z)

#### WEIRD??? HERE THERE SEEM TO BE OVERLAPPING PLOT IDS -- MAYBE A MISTAKE IN THE SCRIPT THAT MADE THESE FILES -- CHECK!!

stemcount.all = rbind(plot.stemcount.spcs.all_AZ, plot.stemcount.spcs.all_CA, plot.stemcount.spcs.all_CO , plot.stemcount.spcs.all_ID, plot.stemcount.spcs.all_MT, plot.stemcount.spcs.all_NM, plot.stemcount.spcs.all_NV, plot.stemcount.spcs.all_OR, plot.stemcount.spcs.all_UT, plot.stemcount.spcs.all_WA, plot.stemcount.spcs.all_WY)

# decide which species to keep then thin out the rest
#### Get mean stemcount per Ha by plot across all species ####
stemcount.total <- aggregate(stemcount.ha ~ SPECIES_SYMBOL, data=stemcount.all, mean)
plot(rev(sort(stemcount.total$stemcount.ha)), col=1)
stemcount.total$SPECIES_SYMBOL[order(stemcount.total$stemcount.ha)]

#### Get number of plots occupied across all species ####
prev.total <- aggregate(stemcount.ha ~ SPECIES_SYMBOL, data=stemcount.all, f<-function(x){return(sum(x>0))})
plot(rev(sort(prev.total$stemcount.ha)), col=1)
prev.total$SPECIES_SYMBOL[order(prev.total$stemcount.ha, decreasing=T)]
prev.total$stemcount.ha[order(prev.total$stemcount.ha, decreasing=T)]

# for now, remove all species that occur in fewer than 100 plots
species.to.keep = prev.total$SPECIES_SYMBOL[prev.total$stemcount.ha>= 100]

stemcount.all.thinned = stemcount.all[stemcount.all$SPECIES_SYMBOL %in% species.to.keep,]
stemcount.all.thinned$SPECIES_SYMBOL = as.character(stemcount.all.thinned$SPECIES_SYMBOL)

#### Convert to wide format ####
stemcount.all.wide <- xtabs(stemcount.ha ~ PLT_CN + SPECIES_SYMBOL, data=stemcount.all.thinned)

head(stemcount.all.wide)
dim(stemcount.all.wide)

# quick look at species abundances
abund=rev(sort(apply(stemcount.all.wide, 2, sum)))
plot(abund, col=0)
text(abund, names(abund), cex=0.75)
prev=rev(sort(apply(stemcount.all.wide, 2, f<-function(x){return(sum(x>0))})))
plot(prev, col=0)
text(prev, names(prev), cex=0.75)
plot(log(prev)~log(abund))

#### write to file ####
#write.csv(plot.stemcount, "./workingdata/FIA_SPxstemcount_CA.csv", row.names=FALSE)



