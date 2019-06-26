#  ------------------------------------------------------------------------
# Read again data:

# data2 = read.csv("simData/paccod_catch_Sim.csv")
# data3 = read.csv("simData/paccod_len_Sim.csv")
# data4 = read.csv("simData/paccod_age_Sim.csv")
# outdat2 = read.csv('simData/abunlen_stratum2Sim.csv')
# NAgeYearMatrix = read.csv('simData/NAgeYearMat.csv')

#  ------------------------------------------------------------------------
# Create matrix to save all prop per year by method
allMethodsPropYear = NULL

#  ------------------------------------------------------------------------
# CALCULATE REAL AGE PROPORTIONS:

# HERE I SHOULD MULTIPLY THIS MATRIX BY SELECTIVITY AT AGE OF THE SURVEY:
NAgeYearMatrix2 = sweep(NAgeYearMatrix, MARGIN=2, SelecSurv, `*`)
PropAgeYearMatrix = sweep(NAgeYearMatrix2, MARGIN=1, rowSums(NAgeYearMatrix2), `/`) # get prop matrix per year

# convert Matrix to DF:
MatNames = list(YEAR = allYears, AGE = allAges)
dimnames(PropAgeYearMatrix) = MatNames

# apply age plus group:
PropAgeYearMatrix[ ,which(as.numeric(colnames(PropAgeYearMatrix)) == agePlus)] = rowSums(PropAgeYearMatrix[ ,which(as.numeric(colnames(PropAgeYearMatrix)) == agePlus):ncol(PropAgeYearMatrix)])
PropAgeYearMatrix2 = PropAgeYearMatrix[,1:which(as.numeric(colnames(PropAgeYearMatrix)) == agePlus)]

PropAgeYearDF = as.data.frame(as.table(as.matrix(PropAgeYearMatrix2)))
names(PropAgeYearDF) = c('YEAR', 'AGE', 'FREQUENCY') # FREQUENCY is p^_a
PropAgeYearDF$AGE = as.numeric(as.character(PropAgeYearDF$AGE))
PropAgeYearDF$METHOD = 'True'

# Save data frame:
allMethodsPropYear = rbind(allMethodsPropYear, PropAgeYearDF)

#  ------------------------------------------------------------------------
# METHOD 1: UNIQUE ALK FOR ALL YEARS. (TO DO: select 0.5stations to simulate a 'data poor' case)

nages = length(allAges) # fix number
nlen = (maxLen - minLen) + 1 # use data3 because length data has more len bins
fakeLen = seq(minLen, maxLen, by = 1)

tmp = data4
tmpMat = matrix(NA, ncol = nages, nrow = nlen)
  
agefac = sort(unique(tmp$AGE)) # THE FIRST AGE MUST BE 1 !!!!!(OK FOR NOW)
for(j in seq_along(agefac)){
    
  tmp2 = tmp[tmp$AGE == agefac[j], ]
  tmp3 = table(tmp2$LENGTH)
  agelen = as.numeric(names(tmp3))
  tmpMat[match(agelen, fakeLen), j] = as.numeric(tmp3)
    
}

colnames(tmpMat) = allAges 
rownames(tmpMat) = fakeLen
tmpMat[which(is.na(tmpMat))] = 0

# Apply age plus group:
tmpMat[ ,which(as.numeric(colnames(tmpMat)) == agePlus)] = rowSums(tmpMat[ ,which(as.numeric(colnames(tmpMat)) == agePlus):ncol(tmpMat)])
tmpMat2 = tmpMat[,1:which(as.numeric(colnames(tmpMat)) == agePlus)]
# Calculating ALK:
outAL = calc_ALK(x = tmpMat2)

# calculate p^_a: (no consider stratum here)
tempo = data3$C_L_I * outAL[match(data3$LENGTH, rownames(outAL)), ]
tempo2 = rowsum(x = tempo, group = data3$ID_HAUL)
yearsfactor = data2$YEAR[match(rownames(tempo2), data2$ID_HAUL)]
tempo3 = rowsum(x = tempo2, group = yearsfactor)
catchtot = aggregate(data2$NUMBER_FISH, list(YEAR = data2$YEAR), sum)
finalmat = sweep(tempo3, MARGIN=1, catchtot$x, `/`)
finalmat = as.data.frame(finalmat)
finalmat$YEAR = rownames(finalmat)
met1df = melt(finalmat, id = c('YEAR'))
names(met1df) = c('YEAR', 'AGE', 'FREQUENCY')

met1df$METHOD = 'Method1'

# Save data frame:
allMethodsPropYear = rbind(allMethodsPropYear, met1df)

#  ------------------------------------------------------------------------
# METHOD 2: UNIQUE ALK FOR EACH YEAR.

nages = length(allAges) # fix number
nlen = (maxLen - minLen) + 1 # use data3 because length data has more len bins
fakeLen = seq(minLen, maxLen, by = 1)

met2df = NULL
for(k in seq_along(allYears)){

	tmp = data4[data4$YEAR == allYears[k], ]
	tmpMat = matrix(NA, ncol = nages, nrow = nlen)
	  
	agefac = sort(unique(tmp$AGE)) # THE FIRST AGE MUST BE 1 !!!!!(OK FOR NOW)
	for(j in seq_along(agefac)){
		
	  tmp2 = tmp[tmp$AGE == agefac[j], ]
	  tmp3 = table(tmp2$LENGTH)
	  agelen = as.numeric(names(tmp3))
	  tmpMat[match(agelen, fakeLen), j] = as.numeric(tmp3)
		
	}

	colnames(tmpMat) = allAges 
	rownames(tmpMat) = fakeLen
	tmpMat[which(is.na(tmpMat))] = 0

	# Apply age plus group:
	tmpMat[ ,which(as.numeric(colnames(tmpMat)) == agePlus)] = rowSums(tmpMat[ ,which(as.numeric(colnames(tmpMat)) == agePlus):ncol(tmpMat)])
	tmpMat2 = tmpMat[,1:which(as.numeric(colnames(tmpMat)) == agePlus)]
	# Calculating ALK:
	outAL = calc_ALK(x = tmpMat2)

	# calculate p^_a: (no consider stratum here)
	data3tmp = data3[data3$YEAR == allYears[k], ]
	data2tmp = data2[data2$YEAR == allYears[k], ]
	
	tempo = data3tmp$C_L_I * outAL[match(data3tmp$LENGTH, rownames(outAL)), ]
	tempo2 = rowsum(x = tempo, group = data3tmp$ID_HAUL)
	tempo3 = colSums(tempo2)
	catchtot = sum(data2tmp$NUMBER_FISH)
	finalmat = tempo3/catchtot
	finalmat = data.frame(YEAR = allYears[k], AGE = names(finalmat), FREQUENCY = as.vector(finalmat))
	met2df = rbind(met2df, finalmat)
	
}

met2df$METHOD = 'Method2'

# Save data frame:
allMethodsPropYear = rbind(allMethodsPropYear, met2df)


#  ------------------------------------------------------------------------
# 3) GAM (Lorenzo's approach) per year: AGE ~ s(LENGTH) + s(LON, LAT)

# continue script:
mydatagam = data4

# start loop over years
yearsfac = sort(unique(mydatagam$YEAR))
#saveModInd = NULL
dwriteAll = NULL
for(j in seq_along(yearsfac)){
  
  subdata = mydatagam[mydatagam$YEAR == yearsfac[j], ]
  data3tmp = data3[data3$YEAR == yearsfac[j], ]
  
  # run the model GAM:
  age_gam = gam(AGE~s(LENGTH)+s(LON,LAT,k=10), data=subdata, family = tw, 
			method = 'REML')


  # predict data
  data3tmp$AGE = as.vector(predict(age_gam,newdata=data3tmp,type='response'))

  #Round ages (makes sense?): 
  data3tmp$AGEROUND = round(data3tmp$AGE,0)
  data3tmp$AGEROUND = ifelse(test = data3tmp$AGEROUND > agePlus, agePlus, data3tmp$AGEROUND)
	
  dwriteAll = rbind(dwriteAll, data3tmp)
  
}

dwriteAll2 = aggregate(dwriteAll$FREQUENCY, list(YEAR = dwriteAll$YEAR, AGE = dwriteAll$AGEROUND), sum)
catchtot = aggregate(data2$NUMBER_FISH, list(YEAR = data2$YEAR), sum)
catchinorder = catchtot$x[match(dwriteAll2$YEAR, catchtot$YEAR)]
dwriteAll2$x = dwriteAll2$x/catchinorder
met3df = dwriteAll2 

names(met3df) = c('YEAR', 'AGE', 'FREQUENCY')

met3df$METHOD = 'Method3'

# Save data frame:
allMethodsPropYear = rbind(allMethodsPropYear, met3df)


#  ------------------------------------------------------------------------
# 4) GAM (CRL approach) per year: PROPAGE ~ LENGTH + s(LON, LAT)

mydatagam = data4
mydatagam = mydatagam[mydatagam$AGE >= firstAgeCRL, ] # CUT IN firstAgeCRL ALL DATA!
# Its better if I set the initial age to model and estimate the final age for each year.
# Its important to discuss this in the manuscript as a drawback of this method.

# define age plus and len vector:
lenvec = allLens

AgeYearTab = table(mydatagam$AGE, mydatagam$YEAR)

# start loop over years
yearsfac = sort(unique(mydatagam$YEAR))
fac = yearsfac

# Estimate age prop per year
# create age prop matrix:
met4df = NULL
#colnames(PropAgeMat) = allages

# Add age-abundance matrix:
#AbunAgeMat = matrix(NA, ncol = length(allages), nrow = length(yearsfac))
#colnames(AbunAgeMat) = sort(allages) 
for(j in seq_along(yearsfac)){

  tmptab = AgeYearTab[ ,colnames(AgeYearTab) == yearsfac[j]]
  tmpages = as.numeric(names(tmptab))
  tmpfreq = as.vector(tmptab)

  plusAgeCRL = agesCRL(ages = tmpages, freq = tmpfreq, thr = 5) # 5 ind as a thr?
	#print(plusAgeCRL)
	
  subdata = mydatagam[mydatagam$YEAR == yearsfac[j], ]
  data3sub = data3[data3$YEAR == yearsfac[j], ]
  
  subdata$AGE = ifelse(test = subdata$AGE > plusAgeCRL, plusAgeCRL, subdata$AGE) # plus group for that specific year

  # all ages for that year
  allages = sort(unique(subdata$AGE))

  # for each year
	PropAgeMat = matrix(NA, ncol = length(allages), nrow = 1) 
	colnames(PropAgeMat) = allages

  # run the model GAM:
  matPreds = matrix(NA, ncol = length(allages), nrow = nrow(data3sub))
  for(ii in seq_along(allages)){
    
    if(ii == length(allages)){
      predvals = rep(1, times = nrow(data3sub))
    } else {
      
      subdata$AGEFAC = ifelse(test = subdata$AGE > allages[ii], 0, 1)
      modtmp = gam(AGEFAC ~ LENGTH + s(LON,LAT,k=10),family = binomial,
                   data = subdata,method = "REML")
      predtmp = predict(modtmp, newdata = data3sub, type = 'response')
      predvals = as.vector(predtmp)
	  elimina = which(subdata$AGEFAC == 1)
	  if(length(elimina) > 0) {
		subdata = subdata[-which(subdata$AGEFAC == 1), ]
	  } else {
		subdata = subdata
	  }
      
    }
    
    matPreds[,ii] = predvals
    
  }
  
  matPreds2 = matrix(NA, ncol = length(allages), nrow = nrow(data3sub))
  for(kk in seq_along(allages)){
    
    if(kk == 1){
      matPreds2[,kk] = matPreds[,kk]
    } else {
      mattmp = 1 - as.matrix(matPreds[,(kk-1):1])
      matPreds2[,kk] =  matPreds[,kk]*apply(X = mattmp, MARGIN = 1, FUN = prod) 
    }
    
  }
    
  
    abunage = t(matPreds2) %*% as.matrix(data3sub$C_L_I)
    
    propage = abunage/sum(abunage)
    propage2 = as.data.frame(propage)
	names(propage2) = 'FREQUENCY'	
	propage2$YEAR = yearsfac[j]
	propage2$AGE = allages
	propage3 = propage2[ , c('YEAR', 'AGE', 'FREQUENCY')]
	
	met4df = rbind(met4df, propage3) 
    
}

met4df$METHOD = 'Method4'

# Save data frame:
allMethodsPropYear = rbind(allMethodsPropYear, met4df)

#  ------------------------------------------------------------------------
# SAVE FILE:

allMethodsPropYear$replicate = ix # add replicate column
if(!simulation){
	write.csv(allMethodsPropYear, paste0('simData/AllPropData', scenarioName, '.csv'), row.names = FALSE) # It is better to calculate RMSEage and RMSEyear and RMSEtot in Excel to avoid confusion.
}

#  ------------------------------------------------------------------------
# Plot all methods: Figure 1
if(!simulation){

	png(paste0('compareProportions_', scenarioName, '.png'), height = 700, width = 900, units = 'px', res = 130)
	print(ggplot(allMethodsPropYear, aes(AGE, FREQUENCY)) +
			geom_line(aes(color = factor(METHOD))) +
			facet_wrap( ~ factor(YEAR), ncol = 5) +
			theme_bw() +
			scale_x_discrete(limits = 0:agePlus))
	dev.off()

	bitmap(paste0('compareProportions_', scenarioName, '.tiff'), height = 190, width = 190, units = 'mm', res = 900)
	print(ggplot(allMethodsPropYear, aes(AGE, FREQUENCY)) +
			geom_line(aes(color = factor(METHOD))) +
			xlab('Age') +
			ylab('Proportion') +
			facet_wrap( ~ factor(YEAR), ncol = 5) +
			theme_bw() +
			theme(legend.position = c(0.75, 0.05), legend.title = element_blank()) +
			scale_x_discrete(limits = 0:agePlus)) 
	dev.off() 

}
#  ------------------------------------------------------------------------
# Plot all methods: Figure 2

nameallmethods = unique(allMethodsPropYear$METHOD)

allMethods = allMethodsPropYear[-which(allMethodsPropYear$METHOD == 'True'), ] # True just to be sure
TrueMethod = allMethodsPropYear[which(allMethodsPropYear$METHOD == 'True'), ]

allMethodsPropYear2 = NULL # THIS DF HAS (Pest - Ptrue)^2
FacYear = unique(TrueMethod$YEAR)
for(k in seq_along(FacYear)){
	
	tmp1 = TrueMethod[TrueMethod$YEAR == FacYear[k], ]
	tmp2 = allMethods[allMethods$YEAR == FacYear[k], ]
	
	FacMethods = unique(tmp2$METHOD)
	for(j in seq_along(FacMethods)){
	
		tmp3 = tmp2[tmp2$METHOD == FacMethods[j], ]
		tmp3 = tmp3[order(tmp3$AGE), ]
		tmp3$FREQUENCY2 = ((tmp3$FREQUENCY - tmp1$FREQUENCY[match(tmp3$AGE, tmp1$AGE)]))*100 # Error term
		tmp3$MSE = ((tmp3$FREQUENCY - tmp1$FREQUENCY[match(tmp3$AGE, tmp1$AGE)]))^2 # MSE
		tmp3$MRE = ((tmp3$FREQUENCY - tmp1$FREQUENCY[match(tmp3$AGE, tmp1$AGE)])/tmp1$FREQUENCY[match(tmp3$AGE, tmp1$AGE)]) # MRE
		allMethodsPropYear2 = rbind(allMethodsPropYear2, tmp3)
	
	}

}

if(!simulation){

	png(paste0('compareProportions2_', scenarioName, '.png'), height = 700, width = 900, units = 'px', res = 130)
	print(ggplot(allMethodsPropYear2, aes(AGE, FREQUENCY2)) +
			geom_line(aes(color = factor(METHOD))) +
			facet_wrap( ~ factor(YEAR), ncol = 5) +
			theme_bw() +
			scale_x_discrete(limits = 0:agePlus)) 
	dev.off()

	bitmap(paste0('compareProportions2_', scenarioName, '.tiff'), height = 190, width = 190, units = 'mm', res = 900)
	print(ggplot(allMethodsPropYear2, aes(AGE, FREQUENCY2)) +
			geom_line(aes(color = factor(METHOD))) +
			xlab('Age') +
			ylab('Error for estimates of proportion-at-age (%)') +
			facet_wrap( ~ factor(YEAR), ncol = 5) +
			theme_bw() +
			theme(legend.position = c(0.75, 0.08), legend.title = element_blank()) +
			scale_x_discrete(limits = 0:agePlus))
	dev.off() 

}

#  ------------------------------------------------------------------------
# Get table: indicators (WARNING!!!!!)
# THIS SHOULD BE DONE AFTER SIMULATIONS!!!!!!
#  ------------------------------------------------------------------------

# ADD HERE PREVIOUS TREATMENT FOR SIMULATION:
#savePerfInd = rbind(savePerfInd, allMethodsPropYear2)
###
write.csv(allMethodsPropYear, paste0('simPerInd1/replicate_', ix, '_', scenarioName, '.csv'), row.names = FALSE)
write.csv(allMethodsPropYear2, paste0('simPerInd2/replicate_', ix, '_', scenarioName, '.csv'), row.names = FALSE)
#allMethodsPropYear2$RMSE_2 = allMethodsPropYear2$RMSE_1^0.5 # THIS IS THE RMSE_a,y. 
