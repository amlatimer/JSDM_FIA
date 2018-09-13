plots <- read.csv("CA_PLOT.CSV",header=TRUE)
trees <- read.csv("CA_TREE.CSV",header=TRUE)
conditions <- read.csv("CA_COND.CSV",header=TRUE)
ref.species <- read.csv("REF_SPECIES.CSV",header=TRUE)

#### Thin plot data ####
# thin plots based on condition codes
# there are many options, but this one thins only to plots for which all conditions are 1 (forest)
plot.cond <- aggregate(COND_STATUS_CD ~ PLT_CN,data=conditions,function(x) {all(x==1)})
plots <- merge(plots,plot.cond,by.x="CN",by.y="PLT_CN")
plots <- plots[plots$COND_STATUS_CD == TRUE,]


#### Thin tree data ####
# thin to only the focal columns
focal.cols <- c("CN","PLT_CN","SUBP","TREE","CONDID","STATUSCD","SPCD","SPGRPCD","DIA","HT","ACTUALHT")
trees <- trees[,focal.cols]

# thin to only focal plots
trees <- trees[trees$PLT_CN %in% plots$CN,]

# thin to only live trees
trees <- trees[trees$STATUSCD == 1,]

# thin to only trees >= 5 inches DBH
trees <- trees[trees$DIA >= 5,]

#### Calc BA per Ha for each tree ####
trees$dia.cm <- trees$DIA * 2.54    # convert inches to cm
trees$ba.ha <- rep(0,length(trees$dia.cm))
trees$ba.ha[trees$dia.cm <= 60.96] <- ((trees$dia.cm[trees$dia.cm <= 60.96]/2)^2) * pi / 0.0672 # plot expansion factor
trees$ba.ha[trees$dia.cm > 60.96]  <- ((trees$dia.cm[trees$dia.cm > 60.96]/2)^2) * pi / 0.405   # plot expansion factor
trees$ba.ha <- trees$ba.ha/100^2 # convert to m2/ha

#### Sum the BA per Ha by plot and species ####
plot.ba.spcs.all <- aggregate(ba.ha ~ SPCD * PLT_CN, data=trees, sum)

#### Get total BA per Ha by plot across all species ####
plot.ba.total <- aggregate(ba.ha ~ PLT_CN, data=trees, sum)
colnames(plot.ba.total) <- c("PLT_CN","TOT_ALL")

#### Look up species abbreviations ####
plot.ba.spcs.all <- merge(plot.ba.spcs.all,ref.species[,c("SPCD","SPECIES_SYMBOL")],by="SPCD",all.x=TRUE)
plot.ba.spcs.all$SPECIES_SYMBOL <- droplevels(plot.ba.spcs.all$SPECIES_SYMBOL)

#### Convert to wide format ####
plot.ba.spcs.all.wide <- xtabs(ba.ha ~ PLT_CN + SPECIES_SYMBOL, data=plot.ba.spcs.all)


#### Narrow down species by removing species with little total BA across all plots ####
# I made cuttoff at 100 total BA/ha
sps.ba <- sort(round(apply(plot.ba.spcs.all.wide, 2, sum), 1), decreasing=T)
sps.ba <- sps.ba[sps.ba > 100]

plot.ba.spcs.wide <- plot.ba.spcs.all.wide[, which(colnames(plot.ba.spcs.all.wide) %in% names(sps.ba))]  
plot.ba.spcs.wide <- as.data.frame(plot.ba.spcs.wide)

#### Add CN column for plot ID, OT column for other species, total plot BA ####
plot.ba.spcs.wide$TOT_FOCAL <- apply(plot.ba.spcs.wide, 1, sum)
plot.ba.spcs.wide$PLT_CN <- rownames(plot.ba.spcs.wide)
plot.ba <- merge(plot.ba.spcs.wide,plot.ba.total,by="PLT_CN",all.x=TRUE)
plot.ba$OTHER <- plot.ba$TOT_ALL - plot.ba$TOT_FOCAL
plot.ba[,-1] <- round(plot.ba[,-1], 3)

#### write to file ####
write.csv(plot.ba, "FIA_SPxBA.csv", row.names=FALSE)
