require(ALKr)
require(reshape2)
require(mgcv)
require(ggplot2)
# --------------------------------------------------------------------------
# Parameters

minEstAge = 1
maxAge = 20
agePlus = 8
minLen = 1
maxLen = 120
iniYear = 1994
endYear = 2016
allYears = iniYear:endYear

#  ------------------------------------------------------------------------
# Read real data:
stratainfo = read.csv('area_stratum.csv')

data2 = read.csv("realData/paccod_catch_1994_2016.csv")
data3 = read.csv("realData/paccod_length_1994_2016.csv")
data4 = read.csv("realData/paccod_age_1994_2016.csv")

#  ------------------------------------------------------------------------
# Compare simulated total abundance and estimated (from sampling) total abundance:

data3$LENGTH = data3$LENGTH/10
data4$LENGTH = data4$LENGTH/10
data4 = data4[data4$AGE > 0, ] # delete age 0

data2$STATIONID = as.character(data2$STATIONID)
data3$STATIONID = as.character(data3$STATIONID)
data4$STATIONID = as.character(data4$STATIONID)

data2$ID_HAUL = paste0(data2$YEAR, data2$STATIONID)
data3$ID_HAUL = paste0(data3$YEAR, data3$STATIONID)
data4$ID_HAUL = paste0(data4$YEAR, data4$STATIONID)


# seeking for lon lat 
data3$LON = data2$START_LONGITUDE[match(data3$ID_HAUL,
                                        data2$ID_HAUL)]
data3$LAT = data2$START_LATITUDE[match(data3$ID_HAUL,
                                       data2$ID_HAUL)]

data4$LON = data2$START_LONGITUDE[match(data4$ID_HAUL,
                                        data2$ID_HAUL)]
data4$LAT = data2$START_LATITUDE[match(data4$ID_HAUL,
                                       data2$ID_HAUL)]


# seeking for stratum
data2$STRATUM_ALT = stratainfo$STRATUM_ALT[match(data2$STRATUM,
												 stratainfo$STRATUM)]
										
data3$STRATUM_ALT = data2$STRATUM_ALT[match(data3$ID_HAUL,
                                        data2$ID_HAUL)]

data4$STRATUM_ALT = data2$STRATUM_ALT[match(data4$ID_HAUL,
                                        data2$ID_HAUL)]


# JUST FOR THE EBS!!!!!!!!!!!!!!!!!!!!!!!!!!!!
data2 = data2[data2$STRATUM_ALT %in% c(10,20,30,40,50,60,82,90), ]
data3 = data3[data3$STRATUM_ALT %in% c(10,20,30,40,50,60,82,90), ]
data4 = data4[data4$STRATUM_ALT %in% c(10,20,30,40,50,60,82,90), ]




# seeking for swept area
data2$SWEPTAREA = ((data2$DISTANCE_FISHED*1000)*data2$NET_WIDTH)*1E-06 # in KM^2
										
data3$SWEPTAREA = data2$SWEPTAREA[match(data3$ID_HAUL,
                                        data2$ID_HAUL)]

data4$SWEPTAREA = data2$SWEPTAREA[match(data4$ID_HAUL,
                                        data2$ID_HAUL)]


# seeking for area stratum:
stratainfo2 = aggregate(stratainfo$AREA_KM, list(STRATUM_ALT = stratainfo$STRATUM_ALT), sum)
data2$AREA_STRATUM_ALT = stratainfo2$x[match(data2$STRATUM_ALT,
												 stratainfo2$STRATUM_ALT)]


# get number of individuals in total catch and subsample:
data5 = aggregate(data3$FREQUENCY, list(ID_HAUL = data3$ID_HAUL), sum)
data3$N_CATCH = data2$NUMBER_FISH[match(data3$ID_HAUL, data2$ID_HAUL)]
data3$N_SUBSAMPLE = data5$x[match(data3$ID_HAUL, data5$ID_HAUL)]

# get lambda:
data3$LAMBDA = data3$N_SUBSAMPLE/data3$N_CATCH

# get c^_l_i
data3$C_L_I = data3$FREQUENCY/data3$LAMBDA # this is c_l_i

# Eq 4 in Stwart and Hommel:
data2$DEN1 = data2$NUMBER_FISH/data2$SWEPTAREA
data2$STRATUMYEAR = paste0(data2$STRATUM_ALT, '_', data2$YEAR)
tempo1 = aggregate(data2$DEN1, list(STRATUMYEAR = data2$STRATUMYEAR), sum)
tempo2 = aggregate(data2$STATIONID, list(STRATUMYEAR = data2$STRATUMYEAR), function(x) length(unique(x)))
tempo1$SECONDPART = tempo1$x/tempo2$x # second part of the equation
tempo3 = aggregate(data2$AREA_STRATUM_ALT, list(STRATUMYEAR = data2$STRATUMYEAR), unique)
tempo1$N_HAT = tempo3$x * tempo1$SECONDPART # This is the equation.
tempo1$YEAR = data2$YEAR[match(tempo1$STRATUMYEAR, data2$STRATUMYEAR)]
tempo1$STRATUM = data2$STRATUM[match(tempo1$STRATUMYEAR, data2$STRATUMYEAR)]


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

# calculate p^_a: (consider stratum here)
tempo = data3$C_L_I * outAL[match(data3$LENGTH, rownames(outAL)), ]
tempo2 = rowsum(x = tempo, group = data3$ID_HAUL)

# continue:
yearsfactor = data2$YEAR[match(rownames(tempo2), data2$ID_HAUL)]
stratumfactor = data2$STRATUM_ALT[match(rownames(tempo2), data2$ID_HAUL)]
uniqfactor = paste0(stratumfactor, '_', yearsfactor)

tempo3 = rowsum(x = tempo2, group = uniqfactor)
catchtotstratum = rowSums(tempo3)
tempo4 = sweep(tempo3, MARGIN=1, catchtotstratum, `/`)

tempo5 = sweep(tempo4, MARGIN=1, tempo1$N_HAT, `*`) # numerator in eq 5 Stwart and Hommel / part 1
tempo6 = aggregate(tempo1$N_HAT, list(YEAR = tempo1$YEAR), sum) # denominator in eq 5 Stwart and Hommel
yearsfactor2 = as.numeric(lapply(strsplit(x = rownames(tempo5), split = '_'), `[[`, 2))
tempo7 = rowsum(x = tempo5, group = yearsfactor2) # numerator in eq 5 Stwart and Hommel / completed
tempo8 = sweep(tempo7, MARGIN=1, tempo6$x, `/`)
tempo9 = as.data.frame(tempo8)
write.csv(tempo9, 'realData/Method1_pcod.csv')
tempo9$YEAR = as.numeric(rownames(tempo9))


met1df = melt(tempo9, id = c('YEAR'))
names(met1df) = c('YEAR', 'AGE', 'FREQUENCY')

met1df$YEAR = as.numeric(as.character(met1df$YEAR))
met1df$AGE = as.numeric(as.character(met1df$AGE))

save1 = met1df
save1$METHOD = 'Method1'
# Save data frame:
allMethodsPropYear = rbind(allMethodsPropYear, save1)



#  ------------------------------------------------------------------------
#  ------------------------------------------------------------------------
# METHOD 2: ALK PER YEAR

allEstAges = minEstAge:maxAge
nages = (maxAge - minEstAge) + 1 # fix number
nlen = (maxLen - minLen) + 1 # use data3 because length data has more len bins
fakeLen = seq(minLen, maxLen, by = 1)

met2df = NULL
tempo9 = NULL
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
	
	
	stratumfactor = data2tmp$STRATUM_ALT[match(rownames(tempo2), data2tmp$ID_HAUL)]
	uniqfactor = stratumfactor

	tempo3 = rowsum(x = tempo2, group = uniqfactor)
	catchtotstratum = rowSums(tempo3)
	tempo4 = sweep(tempo3, MARGIN=1, catchtotstratum, `/`)

	tempo1tmp = tempo1[tempo1$YEAR == allYears[k], ]
	tempo5 = sweep(tempo4, MARGIN=1, tempo1tmp$N_HAT, `*`) # numerator in eq 5 Stwart and Hommel / part 1
	tempo6 = aggregate(tempo1tmp$N_HAT, list(YEAR = tempo1tmp$YEAR), sum) 

	tempo7 = colSums(tempo5) # numerator in eq 5 Stwart and Hommel / completed
	tempo8 = tempo7/tempo6$x
	tempo9tmp = as.data.frame(t(tempo8))
	tempo9 = rbind(tempo9, tempo9tmp)
}

rownames(tempo9) = allYears
write.csv(tempo9, 'realData/Method2_pcod.csv')
tempo9$YEAR = as.numeric(rownames(tempo9))

met2df = melt(tempo9, id = c('YEAR'))
names(met2df) = c('YEAR', 'AGE', 'FREQUENCY')

met2df$YEAR = as.numeric(as.character(met2df$YEAR))
met2df$AGE = as.numeric(as.character(met2df$AGE))

save2 = met2df
save2$METHOD = 'Method2'
# Save data frame:
allMethodsPropYear = rbind(allMethodsPropYear, save2)
	

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

	bitmap(paste0('FinalFigures/Method3_Gau_', scenarioName, '.tiff'), height = 120, width = 120, units = 'mm', res = 900)
	 	gam.check(age_gam)
    dev.off()

  # predict data
  data3tmp$AGE = as.vector(predict(age_gam,newdata=data3tmp,type='response'))

  #Round ages (makes sense?): 
  data3tmp$AGEROUND = round(data3tmp$AGE,0)
  data3tmp$AGEROUND = ifelse(test = data3tmp$AGEROUND > agePlus, agePlus, data3tmp$AGEROUND)
	
  dwriteAll = rbind(dwriteAll, data3tmp)
  
}

dwriteAll2 = aggregate(dwriteAll$C_L_I, list(YEAR = dwriteAll$YEAR, STRATUM_ALT = dwriteAll$STRATUM_ALT, AGE = dwriteAll$AGEROUND), sum)
catchtot = aggregate(dwriteAll2$x, list(YEAR = dwriteAll2$YEAR, STRATUM_ALT = dwriteAll2$STRATUM_ALT), sum)
dwriteAll3 = merge(dwriteAll2, catchtot, by = c('YEAR', 'STRATUM_ALT'))
dwriteAll3$PROP_R_B = dwriteAll3$x.x/dwriteAll3$x.y
dwriteAll3$STRATUMYEAR = paste0(dwriteAll3$STRATUM_ALT, '_', dwriteAll3$YEAR)
dwriteAll3$N_HAT = tempo1$N_HAT[match(dwriteAll3$STRATUMYEAR, tempo1$STRATUMYEAR)]
dwriteAll3$NUM1 = dwriteAll3$PROP_R_B*dwriteAll3$N_HAT
dwriteAll4 = aggregate(dwriteAll3$NUM1, list(YEAR = dwriteAll3$YEAR, AGE = dwriteAll3$AGE), sum)
dwriteAll5 = aggregate(tempo1$N_HAT, list(YEAR = tempo1$YEAR), sum)
dwriteAll4$DEN = dwriteAll5$x[match(dwriteAll4$YEAR, dwriteAll5$YEAR)]
dwriteAll4$OUT = dwriteAll4$x/dwriteAll4$DEN

met3df = dwriteAll4[,c('YEAR', 'AGE', 'OUT')]
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
write.csv(out3, 'realData/Method3_pcod.csv')



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
tempo9 = NULL
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
    
	
	# calculate p^_a: (no consider stratum here)
	data2tmp = data2[data2$YEAR == allYears[j], ]
	data3tmp = data3[data3$YEAR == allYears[j], ]
	
	tempo = sweep(matPreds2, MARGIN=1, data3tmp$C_L_I, `*`)
	tempo2 = rowsum(x = tempo, group = data3tmp$ID_HAUL)
		
	stratumfactor = data2tmp$STRATUM_ALT[match(rownames(tempo2), data2tmp$ID_HAUL)]
	uniqfactor = stratumfactor

	tempo3 = rowsum(x = tempo2, group = uniqfactor)
	catchtotstratum = rowSums(tempo3)
	tempo4 = sweep(tempo3, MARGIN=1, catchtotstratum, `/`)

	tempo1tmp = tempo1[tempo1$YEAR == allYears[j], ]
	tempo5 = sweep(tempo4, MARGIN=1, tempo1tmp$N_HAT, `*`) # numerator in eq 5 Stwart and Hommel / part 1
	tempo6 = aggregate(tempo1tmp$N_HAT, list(YEAR = tempo1tmp$YEAR), sum) 

	tempo7 = colSums(tempo5) # numerator in eq 5 Stwart and Hommel / completed
	tempo8 = tempo7/tempo6$x
	tempo9tmp = as.data.frame(t(tempo8))
	tempo9 = rbind(tempo9, tempo9tmp)
    
}


rownames(tempo9) = allYears
colnames(tempo9) = minEstAge:agePlus
write.csv(tempo9, 'realData/Method4_pcod.csv')
tempo9$YEAR = as.numeric(rownames(tempo9))

met4df = melt(tempo9, id = c('YEAR'))
names(met4df) = c('YEAR', 'AGE', 'FREQUENCY')

met4df$YEAR = as.numeric(as.character(met4df$YEAR))
met4df$AGE = as.numeric(as.character(met4df$AGE))

save4 = met4df
save4$METHOD = 'Method4'
# Save data frame:
allMethodsPropYear = rbind(allMethodsPropYear, save4)
		

# -----------------------------------------------------------------------
# Plot all methods:
allMethodsPropYear$AGE = as.numeric(as.character(allMethodsPropYear$AGE))

ggplot(allMethodsPropYear, aes(x = AGE, y = FREQUENCY)) +
	geom_line(aes(color = METHOD)) +
  facet_wrap( ~ factor(YEAR), nrow = 5, dir = 'v') +
  theme_bw()



#aggregate(allMethodsPropYear$FREQUENCY, list(allMethodsPropYear$YEAR, allMethodsPropYear$METHOD), sum)