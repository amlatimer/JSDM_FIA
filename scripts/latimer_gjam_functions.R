# Functions for using gjam library to analyze species distributions of trees in Sierra FIA plots. Most are more general and can operate on any gjam outputs. 

rain.colors = colorRampPalette(c("red", "lightgray", "blue"))

plotBetaChains <- function(out, variable, thinby) {
  betas = out$chains$bgibbs
  betas = betas[(out$burnin+1):nrow(betas),]
  z = sapply(colnames(betas), strsplit, split="_")
  z2= list()
  for (i in 1:length(z)) z2[[i]] = (z[[i]][2])
  columns = which(z2==variable)
  thinx = as.integer((1:(nrow(betas)/thinby))*thinby)
  plotcols=rep(c("brown","tan", "gold", "black", "darkgray", "lightgray"), 4)
  plot(betas[thinx,columns[1]], type="l", col=plotcols[1], ylim=c(min(betas[,columns]), max(betas[,columns])), xlab="Sample", ylab="Parameter value", cex.lab=1.5)
  for (i in 2:length(columns)) lines(betas[thinx,columns[i]], type="l", col=plotcols[i])
}

plotBetaBoxplot <- function(out, variable) {
  betas = out$chains$bgibbs
  betas = betas[(out$burnin+1):nrow(betas),]
  z = sapply(colnames(betas), strsplit, split="_")
  z2= list()
  labs = vector("character", length= ncol(betas))
  for (i in 1:length(z)) {
    z2[[i]] = (z[[i]][2])
    labs[i] = substr(z[[i]][1], start=4, stop=7)
  }
  columns = which(z2==variable)
  colnames(betas) = labs
  betasub = betas[,columns]
  parmean = apply(betasub, 2, mean)
  betasub = betasub[,order(parmean, decreasing=TRUE)]
  boxplot(betasub, main=paste(variable), col="orange3", xlab="", ylab="parameter value", cex.lab=1.2, pars=list(las=2))
  abline(0,0,col="darkgray", lty=2)
}

DIPscore <- function(x, xhat, predsigma) {
	scores = -((x-xhat) / predsigma)^2 - log(predsigma^2)
	return(scores)
}

scale.colors <- function(x, n.colors) { # take a vector of real numbers and convert to colors on scale 1 to 100
  if (is.null(n.colors)) colors=128
  temp1 = (x-min(x))
  temp2 = temp1 * (n.colors-1)/max(temp1)
  return(round(temp2) + 1)
}

sens.colors <- colorRampPalette(c("black", "black", "darkgray", "gray", "gray", "lightgray", "cyan", "yellow"))

plotmap <- function(plotlocs, x, colorfun) {
	x.col = scale.colors(x, 128)
  palette(colorfun(max(x.col)))
	#map("county", "CA", col="brown")
	#points(plotlocs, pch=16, cex=0.5, col="lightgray")
	plot(plotlocs, pch=16, cex=0.8, col=x.col)
}

# Function to plot results of out-of-sample prediction
pred_vs_obs_scatter <- function(x, xpred) {
  nplots = ncol(x)
  par(mfrow=c(sqrt(nplots), sqrt(nplots)+1))
  for (i in 2:nplots) { 
    plot(x[,i], xpred[,i], xlab="observed", ylab="predicted", cex.lab=1.5, cex.axis=1.2) 
    MSPE = mean((xpred[,i] - x[,i])^2)
    title(main=(colnames(x)[i]), sub=paste("MSPE =", round(MSPE, 3)))
    abline(0,1)
  }
}

sensitivityPlot <- function(model, varname, colorfun, showquantiles) { # model is an output from gjamGibbs(), variable name is the explanatory variable name (in the x matrix of the model). Colorfun is a palette function that is used to plot the variable, showquantiles is T/F saying whether colors should assigned by quantile instead
  xcolumn = match(varname, colnames(model$x))
  z = DIPscore(x=model$x[,xcolumn], xhat=model$modelSummary$xpredMu[,xcolumn], predsigma = model$modelSummary$xpredSd[,xcolumn])
  #hist(z); mean(z)
  if (showquantiles) { 
    x.col = rep("", length(z))
    z_quant = quantile(z, probs=c(0.5, 0.75, 0.9, 0.95))
    x.col[z<z_quant[1]] = "lightgray"
    x.col[z>=z_quant[1] & z<z_quant[2]] = "darkgray"
    x.col[z>=z_quant[2] & z<z_quant[3]] = "blue3"
    x.col[z>=z_quant[3] & z<z_quant[4]] = "cyan3"
    x.col[z>=z_quant[4]] = "gold"
    plot(plotlocs, pch=16, cex=0.8, col=x.col)
  } else plotmap(plotlocs, x=z, colorfun)
}

gjamPredictY <- function(newx, out, n.samples) {
  nsp = ncol(out$modelSummary$betaMu)
  n = nrow(out$chains$rgibbs)
  sampleindex = sample(1:n, n.samples, replace=TRUE)
  wpred = matrix(NA, n.samples, nsp)
  for (i in 1:n.samples) {
    betasim = matrix(out$chains$bgibbs[sampleindex[i],], nrow=nsp, byrow=TRUE)
    wpredmu = betasim %*% newx
    ssim = matrix(out$chains$sgibbs[sampleindex[i],], nrow=nsp, byrow=TRUE)
    wpred[i,] = rmvnorm(1, mean=wpredmu, sigma=ssim) # simulate from multivariate normal 
  }
  wpred[wpred<0] = 0 # truncate mvn draws at 0 (applies for continuous abundance only)
  ypred = apply(wpred, 2, mean) # take average of truncated normal draws as predicted abundance 
  return(ypred)
}

gjamPredictY.region <- function(newx, out, samples) { # Difference from gjamPredictY() is that this takes newx for a whole region, so newx is a matrix, where each row contains X values for one location/grid cell. Output is a list in which each element is one posterior simulation of Y for the whole region (an nsp x ncells matrix).
  # function to sample y for an individual cell
  sample.cell <- function(x, out, sampleindex) {
    nsp = ncol(out$modelSummary$betaMu)
    n = nrow(out$chains$rgibbs)
    betasim = matrix(out$chains$bgibbs[sampleindex,], nrow=nsp, byrow=TRUE)
    wpredmu = betasim %*% x
    ssim = matrix(out$chains$sgibbs[sampleindex,], nrow=nsp, byrow=TRUE)
    wpred = rmvnorm(1, mean=wpredmu, sigma=ssim) # simulate from multivariate normal
    wpred[wpred<0] = 0 # truncate at 0 (for continuous data only)
    return(wpred)
  }
  
  n.samples = length(samples)
  nsp = ncol(out$modelSummary$betaMu)
  ncells = nrow(newx)
  n = nrow(out$chains$rgibbs)
  ypred = array(NA, dim=c(nsp, ncells, n.samples))
    
  for (i in 1:n.samples) {
    #sampleindex = sample(1:n, 1, replace=TRUE)
    sampleindex = samples[i]
    ypred[,,i] = apply(newx, 1, sample.cell, out=out, sampleindex=sampleindex)
  }

  return(ypred)
}

gjamPredictY.postmean <- function(newx, out) { # Difference from gjamPredictY.region() is that this predicts for the posterior mean values of the model out only
  # function to sample y for an individual cell
  sample.cell <- function(x, out) {
    nsp = ncol(out$modelSummary$betaMu)
    beta_hat = matrix(out$modelSummary$betaMu, nrow=nsp, byrow=TRUE)
    wpredmu = beta_hat %*% x
    s_hat = matrix(out$modelSummary$sigMu, nrow=nsp, byrow=TRUE)
    wpred = rmvnorm(1, mean=wpredmu, sigma=s_hat) # simulate from multivariate normal
    wpred[wpred<0] = 0 # truncate at 0 (for continuous data only)
    return(wpred)
  }
  
  ypred = apply(newx, 1, sample.cell, out=out)
  return(ypred)
}

gjamCalcChange.jaccard <- function(ypred.hist, ypred.future) { # Takes current and future projections of species abundance (ypred.hist for historical and ypred.future for a future scenario) and calculates mean and standard deviation of change for each grid cell between the two environments
  # Calculate change in composition -- Jaccard dissimilarity
  ncells = ncol(ypred.hist)
  n.samples = dim(ypred.hist)[3]
  jaccard.mat = matrix(NA, nrow=ncells, ncol=n.samples)
  for (j in 1:ncells) {
    for (k in 1:n.samples) {
      vegmat = rbind(ypred.future[,j,k], ypred.hist[,j,k])
      if (sum(vegmat[1,])==0 | sum(vegmat[2,])==0) jaccard.mat[j,k] = NA else jaccard.mat[j,k] = vegdist(vegmat, method="jaccard", binary=FALSE, na.rm=T)
    }
  }
  jaccard.mean = apply(jaccard.mat, 1, mean, na.rm=T)
  jaccard.sd = apply(jaccard.mat, 1, sd, na.rm=T)
  return(list(jaccard.mean, jaccard.sd))
}

gjamCalcChange.vegdist <- function(ypred.hist, ypred.future, distmethod) { # Takes current and future projections of species abundance (ypred.hist for historical and ypred.future for a future scenario) and calculates mean and standard deviation of change for each grid cell between the two environments
  # Calculate change in composition -- apply vegdist using method=distmethod
  ncells = ncol(ypred.hist)
  n.samples = dim(ypred.hist)[3]
  vegdist.mat = matrix(NA, nrow=ncells, ncol=n.samples)
  for (j in 1:ncells) {
    for (k in 1:n.samples) {
      vegmat = rbind(ypred.future[,j,k], ypred.hist[,j,k])
      if (sum(vegmat[1,])==0 | sum(vegmat[2,])==0) vegdist.mat[j,k] = NA else vegdist.mat[j,k] = vegdist(vegmat, method=distmethod, na.rm=T)
    }
  }
  vegdist.mean = apply(vegdist.mat, 1, mean, na.rm=T)
  vegdist.sd = apply(vegdist.mat, 1, sd, na.rm=T)
  return(list(vegdist.mean, vegdist.sd))
}

gjamCalcChange.abschange <- function(ypred.hist, ypred.future) { # Takes current and future projections of species abundance (ypred.hist for historical and ypred.future for a future scenario) and calculates mean and standard deviation of summed absolute change for each grid cell across all species
  ncells = ncol(ypred.hist)
  n.samples = dim(ypred.hist)[3]
  change.mat = matrix(NA, nrow=ncells, ncol=n.samples)
  for (j in 1:ncells) {
    for (k in 1:n.samples) {
       change.mat[j,k] = sum(abs(ypred.future[,j,k]-ypred.hist[,j,k]))
    }
  }
  change.mean = apply(change.mat, 1, mean, na.rm=T)
  change.sd = apply(change.mat, 1, sd, na.rm=T)
  return(list(change.mean, change.sd))
}

gjamCalcChange.basalarea <- function(ypred.hist, ypred.future) { # Takes current and future projections of species abundance (ypred.hist for historical and ypred.future for a future scenario) and calculates mean and standard deviation of basal area change for each grid cell between the two environments
  ncells = ncol(ypred.hist)
  n.samples = dim(ypred.hist)[3]
  ba.mat = matrix(NA, nrow=ncells, ncol=n.samples)
  for (j in 1:ncells) for (k in 1:n.samples) ba.mat[j,k] = sum(ypred.future[,j,k])-sum( ypred.hist[,j,k])
  ba.mean = apply(ba.mat, 1, mean)
  ba.sd = apply(ba.mat, 1, sd)
  return(list(ba.mean, ba.sd))
}
  
gjamCalcChange.hardwood <- function(ypred.hist, ypred.future, spp) { # Takes current and future projections of species abundance (ypred.hist for historical and ypred.future for a future scenario) and calculates mean and standard deviation of basal area change for each grid cell between the two environments
  ncells = ncol(ypred.hist)
  n.samples = dim(ypred.hist)[3]
  hw.mat = matrix(NA, nrow=ncells, ncol=n.samples)
  hw.index = grep("QU", spp)
  for (j in 1:ncells) for (k in 1:n.samples) hw.mat[j,k] = sum(ypred.future[hw.index,j,k]) -sum(ypred.hist[hw.index,j,k])
  hw.mean = apply(hw.mat, 1, mean)
  hw.sd = apply(hw.mat, 1, sd)
  return(list(hw.mean, hw.sd))
}

gjamCalcChange.richness <- function(ypred.hist, ypred.future) { # Takes current and future projections of species abundance (ypred.hist for historical and ypred.future for a future scenario) and calculates mean and standard deviation of basal area change for each grid cell between the two environments
  ncells = ncol(ypred.hist)
  n.samples = dim(ypred.hist)[3]
  richness.mat = matrix(NA, nrow=ncells, ncol=n.samples)
  for (j in 1:ncells) for (k in 1:n.samples) richness.mat[j,k] = sum(ypred.future[,j,k]>0) - sum(ypred.hist[,j,k]>0)
  richness.mean = apply(richness.mat, 1, mean)
  richness.sd = apply(richness.mat, 1, sd)
  return(list(richness.mean, richness.sd))
}


plot.spdist.raster <- function(r, p, cell.index, sp){ # plot preducted species distributions and compare to observed data. r is the raster, cell.index tells which cells of the raster are data (not NAs), and p is the value to plot, and sp is the name of the species. 
  vals = rep(NA, ncell(r)); vals[cell.index]=p
  r = setValues(r, vals)
  plot(r)
  plot(plotlocs[y[,spp==sp]==0], pch=16, cex=0.5, col="darkgray", add=T)
  plot(plotlocs[y[,spp==sp]>0], pch=16, cex=0.5, col="black", add=T)
}

plot.layer.raster <- function(r, p){ # plot a set of model predictions to raster values, where r is the raster and p is the set of new values for the cells of the raster that are not NA. 
  cell.index = which(!is.na(getValues(r)))
  vals = rep(NA, ncell(r)); vals[cell.index]=p
  r = setValues(r, vals)
  plot(r)
}

plot.layer.raster.colorfun <- function(r, p, colorfun){ # plot a set of model predictions to raster values, where r is the raster and p is the set of new values for the cells of the raster that are not NA. 
  cell.index = which(!is.na(getValues(r)))
  vals = rep(NA, ncell(r)); vals[cell.index]=p
  r = setValues(r, vals)
  plot(r, col=colorfun(128))
}

make.sierra.raster <- function(r, p) {
  #return a raster layer that has NA's in the grid cells that are NA in the sierra raster layer, and the new cell values p in the other cells. 
  cell.index = which(!is.na(getValues(r)))
  vals = rep(NA, ncell(r)); vals[cell.index]=p
  rnew = setValues(r, vals)
  return(rnew)
}


hardwood_change <- function(ynew, yold, spp) {
  # yold is predicted or fitted abundances under current x
  # ynew is projected abundances under perturbed x (environment)
  # spp is the vector of species names in y
  #spp = colnames(yold)
  hw.index = grep("QU", spp) # hardwoods. Everythign else is conifer
  hw_old = apply(yold[hw.index,], 2, sum) 
  hw_new = apply(ynew[hw.index,], 2, sum) 
  return(hw_new-hw_old)
}


hardwood_conifer_change <- function(ynew, yold, spp) {
  # yold is predicted or fitted abundances under current x
  # ynew is projected abundances under perturbed x (environment)
  # spp is the vector of species names in y
  spp = colnames(yold)
  hw.index = grep("QU", spp) # hardwoods. Everythign else is conifer
  hw_old = apply(yold[hw.index,], 2, sum) / apply(yold, 2, sum)
  hw_new = apply(ynew[hw.index,], 2, sum) / apply(ynew, 2, sum)
  return(hw_new-hw_old)
}


### Functions to calculate change for every scenario x posterior sample combination

gjamCalcChange.vegdist.allmodels <- function(ypred.allmodels, ypred.hist, distmethod) { # Takes current and future projections of species abundance (ypred.hist for historical and ypred.future for a future scenario) and calculates mean and standard deviation of change for each grid cell between the two environments
  # Calculate change in composition -- apply vegdist using method=distmethod
  n.models = length(ypred.allmodels)
  ncells = ncol(ypred.hist)
  n.samples = dim(ypred.hist)[3]
  vegdist.array = array(NA, list(ncells, n.samples, n.models))
  for (i in 1:n.models) {
    ypred.future = ypred.allmodels[[i]]
    for (j in 1:ncells) {
      for (k in 1:n.samples) {
        vegmat = rbind(ypred.future[,j,k], ypred.hist[,j,k])
        if (sum(vegmat[1,])==0 | sum(vegmat[2,])==0) vegdist.array[j,k,i] = NA else vegdist.array[j,k,i] = vegdist(vegmat, method=distmethod, na.rm=T)
      }
    }
  }
  return(vegdist.array)
}