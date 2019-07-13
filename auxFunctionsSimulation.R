### Some functions to use later on:
#  Function to begin a new cohort: delete the last element and add a new rec
toNewYear = function(vec, firstVal = 5){
  
  newvec = vec[-length(vec)]
  newvec = c(firstVal, newvec)
  return(newvec)
  
}

toLenAgeMatSampled = function(mat, vec){

	tmpmat = matrix(NA, ncol = ncol(mat), nrow = nrow(mat))
	for(lk in 1:nrow(mat)){
	
		if(vec[lk] > 0){
			selcol = sample(x = 1:ncol(mat), size = vec[lk], replace = TRUE, prob = mat[lk,])
			tmp2 = table(selcol)
			tmpmat[lk ,as.numeric(names(tmp2))] = as.vector(tmp2)
		}
	
	}
	tmpmat[is.na(tmpmat)] = 0
	
	return(tmpmat)

}

sampleVec = function(x, ...) x[sample(length(x), ...)]


ageSample = function(mat, vec){

	newlist = apply(X = mat, MARGIN = 1, FUN = rep, x = allAges)
	newvec2 = mapply(sampleVec, newlist, size = vec, replace = FALSE) # sample of ages should be with a probability of age abundance?
	newvec3 = unlist(newvec2)
	return(newvec3)

}

NindAgeSample = function(vec, maxSam = 50){
	
	nsam = sum(vec)
	nrandom = nsam - maxSam
	while(nrandom > 0){
		vec[which.max(vec)] = vec[which.max(vec)] - 1
		nsam = sum(vec)
		nrandom = nsam - maxSam
	}

	return(vec)

}


  agesCRL = function(ages, freq, thr = 1){
	
	firstZero = which(freq <= thr)[1]
	cutAge2 = ages[firstZero - 1]
	
	return(cutAge2)
	
  }




map.heatmap <- function (lat, lon, data, 
                         color_low="white",color_high="darkred", 
						 color_na=gray(0.9),zeroiswhite=FALSE,
                         xlim=NULL, ylim=NULL, zlim=NULL,
                         mainTitle="", legendTitle="") {

  # Store the base data of the underlying map
  baseData <- map_data("state")

  # Combine the data into a dataframe
  dfMap           <- as.data.frame(cbind(lon, lat, data))
  colnames(dfMap) <- c("lon", "lat", "Value")

  # Set limits for x, y, z if not specified as parameters
  if (is.null(xlim)) { xlim <- range( lon,na.rm=TRUE) }
  if (is.null(ylim)) { ylim <- range( lat,na.rm=TRUE) }
  if (is.null(zlim)) { zlim <- range(data,na.rm=TRUE) }

  # Create the plot
  p <- ggplot(dfMap, aes(x=lon, y=lat, fill=Value)) + theme_bw()
  p <- p + geom_tile()
  p <- p + geom_polygon(data=baseData, aes(x=long, y=lat, group=group), 
                        colour="black", fill=8) 
  p <- p + labs(title=paste(mainTitle,"\n",sep=""), x="", y="")
  #p <- p + theme(plot.title = element_text(size = rel(1.5))) 
  p <- p + coord_fixed(ratio=1, xlim=xlim, ylim=ylim)

  if(zeroiswhite){
    p <- p + scale_fill_gradient2(low=color_low,
                                  high=color_high,
                                  na.value=color_na,
                                  limits=zlim,
                                  name=legendTitle) 
  }
  if(!zeroiswhite){
    p <- p + scale_fill_gradient(low=color_low, 
                                 high=color_high,
                                 na.value=color_na,
                                 limits=zlim,
                                 name=legendTitle) 
  }

return(p)}


# addalpha()
addalpha <- function(colors, alpha=1.0) {
  r <- col2rgb(colors, alpha=T)
  # Apply alpha
  r[4,] <- alpha*255
  r <- r/255.0
  return(rgb(r[1,], r[2,], r[3,], r[4,]))
}

# colorRampPaletteAlpha()
colorRampPaletteAlpha <- function(colors, n=32, interpolate='linear') {
  # Create the color ramp normally
  cr <- colorRampPalette(colors, interpolate=interpolate)(n)
  # Find the alpha channel
  a <- col2rgb(colors, alpha=T)[4,]
  # Interpolate
	if (interpolate=='linear') {
		l <- approx(a, n=n)
	} else {
		l <- spline(a, n=n)
	}
	l$y[l$y > 255] <- 255 # Clamp if spline is > 255
	cr <- addalpha(cr, l$y/255.0)
	return(cr)
}


normDF = function(dat, fac1 = 'AGE', fac2 = 'YEAR', normcol = 'LENGTH', maxAge = 8){

	outdat = NULL
	dat = dat[dat[,fac1] <= maxAge, ]
	dat[, fac1] = as.character(dat[, fac1])
	fctrs2 = sort(unique(dat[, fac2]))
	for(k in seq_along(fctrs2)){
	tmp2 = dat[dat[, fac2] == fctrs2[k], ]
	fctrs1 = sort(unique(tmp2[, fac1]))
		for(l in seq_along(fctrs1)){
			tmp = tmp2[tmp2[, fac1] == fctrs1[l], ]
			tmp[, normcol] = normalize(x = tmp[, normcol], method = 'range', range = c(-1,1))
			outdat = rbind(outdat, tmp)
		}
	}
	return(outdat)

}



