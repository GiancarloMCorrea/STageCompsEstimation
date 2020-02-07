require(ALKr)
require(reshape2)
require(mgcv)
require(ggplot2)
# --------------------------------------------------------------------------
# Parameters

iniYear = iniYearSam
allYears = iniYear:endYear

#  ------------------------------------------------------------------------
# Read real data:

data2 = read.csv(paste0("simData/paccod_catch_Sim_", scenarioName, '.csv'))
data3 = read.csv(paste0("simData/paccod_len_Sim_", scenarioName, '.csv'))
data4 = read.csv(paste0("simData/paccod_age_Sim_", scenarioName, '.csv'))

#  ------------------------------------------------------------------------
# Compare simulated total abundance and estimated (from sampling) total abundance:

data4 = data4[data4$AGE > 0, ] # delete age 0

data2$STATIONID = as.character(data2$STATIONID)
data3$STATIONID = as.character(data3$STATIONID)
data4$STATIONID = as.character(data4$STATIONID)

data2$ID_HAUL = paste0(data2$YEAR, data2$STATIONID)
data3$ID_HAUL = paste0(data3$YEAR, data3$STATIONID)
data4$ID_HAUL = paste0(data4$YEAR, data4$STATIONID)


# get number of individuals in total catch and subsample:
data5 = aggregate(data3$FREQUENCY, list(ID_HAUL = data3$ID_HAUL), sum)
data3$N_CATCH = data2$NUMBER_FISH[match(data3$ID_HAUL, data2$ID_HAUL)]
data3$N_SUBSAMPLE = data5$x[match(data3$ID_HAUL, data5$ID_HAUL)]

# get lambda:
data3$LAMBDA = data3$N_SUBSAMPLE/data3$N_CATCH

# get c^_l_i
data3$C_L_I = data3$FREQUENCY/data3$LAMBDA # this is c_l_i

# Estimate total abundance:
Factor2 = GridArea/areaSwept
data2$NUMBER_FISH_GRID = data2$NUMBER_FISH * Factor2
write.csv(aggregate(data2$NUMBER_FISH_GRID, list(YEAR = data2$YEAR), sum), 
			paste0('simData/AbundanceEstSim_', scenarioName, '.csv'), row.names = FALSE)

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------
allMethodsPropYear = NULL

#  ------------------------------------------------------------------------
#  ------------------------------------------------------------------------
# METHOD 1: UNIQUE ALK FOR ALL YEARS

allEstAges = minEstAge:maxAge
nages = (maxAge - minEstAge) + 1 # fix number
nlen = (maxLen - minLen) + 1 # use data3 because length data has more len bins
fakeLen = seq(minLen, maxLen, by = 1)

# select just the 50% of the sampling stations per year
data2$INDEX = data2$ID_HAUL
data4$INDEX = data4$ID_HAUL
data6 = data4

# begin analysis
tmp = data6
tmpMat = matrix(NA, ncol = nages, nrow = nlen)

agefac = sort(unique(tmp$AGE)) # THE FIRST AGE MUST BE 1 !!!!!(OK FOR NOW)
for(j in seq_along(agefac)){
    
  tmp2 = tmp[tmp$AGE == agefac[j], ]
  tmp3 = table(tmp2$LENGTH)
  agelen = as.numeric(names(tmp3))
  tmpMat[match(agelen, fakeLen), agefac[j]] = as.numeric(tmp3)
    
}

colnames(tmpMat) = allEstAges 
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

nSamLen = aggregate(data3$STATIONID, list(data3$YEAR), function(x) length(unique(x))) # should it be number of tows with len sample?
names(nSamLen) = c('YEAR', 'NSAM')

#write.csv(nSamLen, 'simData/nSamLen_pcod.csv')


# ------------------------------------------------------------------------------------

# continue:
yearsfactor = data2$YEAR[match(rownames(tempo2), data2$ID_HAUL)]
tempo3 = rowsum(x = tempo2, group = yearsfactor)
#catchtot = aggregate(data2$NUMBER_FISH, list(YEAR = data2$YEAR), sum)
catchtot = data.frame(YEAR = names(rowSums(tempo3)), x = as.vector(rowSums(tempo3)))
finalmat = sweep(tempo3, MARGIN=1, catchtot$x, `/`)
finalmat = as.data.frame(finalmat)
finalmat$YEAR = rownames(finalmat)
met1df = melt(finalmat, id = c('YEAR'))
names(met1df) = c('YEAR', 'AGE', 'FREQUENCY')

met1df$YEAR = as.numeric(as.character(met1df$YEAR))
met1df$AGE = as.numeric(as.character(met1df$AGE))

save1 = met1df
save1$METHOD = 'Method1'
# Save data frame:
allMethodsPropYear = rbind(allMethodsPropYear, save1)

out1 = matrix(NA, ncol = length(minEstAge:agePlus), nrow = length(allYears))
for(ij in seq_along(allYears)){

	tmpx = met1df[met1df$YEAR == allYears[ij], ]
	out1[ij,] = tmpx$FREQUENCY

}

colnames(out1) = minEstAge:agePlus
rownames(out1) = allYears
write.csv(out1, paste0('simData/Method1_pcod_', scenarioName, '.csv'))

#  ------------------------------------------------------------------------
#  ------------------------------------------------------------------------
# METHOD 2: ALK PER YEAR

allEstAges = minEstAge:maxAge
nages = (maxAge - minEstAge) + 1 # fix number
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
	  tmpMat[match(agelen, fakeLen), agefac[j]] = as.numeric(tmp3)
		
	}

	colnames(tmpMat) = allEstAges 
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
	#catchtot = sum(data2tmp$NUMBER_FISH)
	catchtot = sum(tempo3)
	finalmat = tempo3/catchtot
	finalmat = data.frame(YEAR = allYears[k], AGE = names(finalmat), FREQUENCY = as.vector(finalmat))
	met2df = rbind(met2df, finalmat)
	
}

met2df$YEAR = as.numeric(as.character(met2df$YEAR))
met2df$AGE = as.numeric(as.character(met2df$AGE))

save2 = met2df
save2$METHOD = 'Method2'
# Save data frame:
allMethodsPropYear = rbind(allMethodsPropYear, save2)


out2 = matrix(NA, ncol = length(minEstAge:agePlus), nrow = length(allYears))
for(ij in seq_along(allYears)){

	tmpx = met2df[met2df$YEAR == allYears[ij], ]
	out2[ij,] = tmpx$FREQUENCY

}

colnames(out2) = minEstAge:agePlus
rownames(out2) = allYears
write.csv(out2, paste0('simData/Method2_pcod_', scenarioName, '.csv'))


#  ------------------------------------------------------------------------
#  ------------------------------------------------------------------------
# 3) GAM (Lorenzo's approach) per year: AGE ~ s(LENGTH) + s(LON, LAT)

# continue script:
mydatagam = data4

# start loop over years
yearsfac = allYears
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
  data3tmp$AGEROUND = ifelse(test = data3tmp$AGEROUND < minEstAge, minEstAge, data3tmp$AGEROUND)
	
  dwriteAll = rbind(dwriteAll, data3tmp)
  
}

dwriteAll2 = aggregate(dwriteAll$C_L_I, list(YEAR = dwriteAll$YEAR, AGE = dwriteAll$AGEROUND), sum)
catchtot = aggregate(dwriteAll2$x, list(YEAR = dwriteAll2$YEAR), sum)
catchinorder = catchtot$x[match(dwriteAll2$YEAR, catchtot$YEAR)]
dwriteAll2$x = dwriteAll2$x/catchinorder
met3df = dwriteAll2 

names(met3df) = c('YEAR', 'AGE', 'FREQUENCY')

met3df$YEAR = as.numeric(as.character(met3df$YEAR))
met3df$AGE = as.numeric(as.character(met3df$AGE))

save3 = met3df
save3$METHOD = 'Method3'
# Save data frame:
allMethodsPropYear = rbind(allMethodsPropYear, save3)


out3 = matrix(NA, ncol = length(minEstAge:agePlus), nrow = length(allYears))
for(ij in seq_along(allYears)){

	tmpx = met3df[met3df$YEAR == allYears[ij], ]
	out3[ij,] = tmpx$FREQUENCY

}

colnames(out3) = minEstAge:agePlus
rownames(out3) = allYears
write.csv(out3, paste0('simData/Method3_pcod_', scenarioName, '.csv'))



#  ------------------------------------------------------------------------
#  ------------------------------------------------------------------------
# 4) GAM (CRL approach) per year: PROPAGE ~ LENGTH + s(LON, LAT)


mydatagam = data4
#mydatagam = mydatagam[mydatagam$AGE >= minEstAge, ] # CUT IN firstAgeCRL ALL DATA!
# Its better if I set the initial age to model and estimate the final age for each year.
# Its important to discuss this in the manuscript as a drawback of this method.

# define age plus and len vector:
lenvec = fakeLen

AgeYearTab = table(mydatagam$AGE, mydatagam$YEAR)

# start loop over years
yearsfac = allYears
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

  # plusAgeCRL = agesCRL(ages = tmpages, freq = tmpfreq, thr = 5) # 5 ind as a thr?
  plusAgeCRL = agePlus
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
      modtmp = gam(AGEFAC ~ LENGTH + s(LON, LAT, k = 10),family = binomial,
                   data = subdata, method = "ML")
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

met4df$YEAR = as.numeric(as.character(met4df$YEAR))
met4df$AGE = as.numeric(as.character(met4df$AGE))

save4 = met4df
save4$METHOD = 'Method4'
# Save data frame:
allMethodsPropYear = rbind(allMethodsPropYear, save4)


out4 = matrix(NA, ncol = length(minEstAge:agePlus), nrow = length(allYears))
for(ij in seq_along(allYears)){

	tmpx = met4df[met4df$YEAR == allYears[ij], ]
	out4[ij,] = tmpx$FREQUENCY

}

colnames(out4) = minEstAge:agePlus
rownames(out4) = allYears
write.csv(out4, paste0('simData/Method4_pcod_', scenarioName, '.csv'))


# -----------------------------------------------------------------------
# Plot all methods:
allMethodsPropYear$AGE = as.numeric(as.character(allMethodsPropYear$AGE))

ggplot(allMethodsPropYear, aes(x = AGE, y = FREQUENCY)) +
	geom_line(aes(color = METHOD)) +
  facet_wrap( ~ factor(YEAR), nrow = 5) +
  theme_bw()



#aggregate(allMethodsPropYear$FREQUENCY, list(allMethodsPropYear$YEAR, allMethodsPropYear$METHOD), sum)