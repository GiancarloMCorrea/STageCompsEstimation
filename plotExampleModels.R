
# Read simulated data:
require(ALKr)
require(reshape2)
require(ggplot2)
require(mgcv)

rm(list = ls())

setwd(dir = 'C:/Users/moroncog/Documents/GitHub/STageCompsEstimation/')

data2 = read.csv("simData/paccod_catch_Sim_HighS_HighT.csv")
data3 = read.csv("simData/paccod_len_Sim_HighS_HighT.csv")
data4 = read.csv("simData/paccod_age_Sim_HighS_HighT.csv")

# ------------------------------------------------------------------------
# Parameters:
maxAge = 20
minEstAge = 1
maxLen = 120
minLen = 1
lenBin = 1
nSamLoc = round(349*0.95)
agePlus = 8
allYears = 1994:2016
allLens = seq(from = minLen, to = maxLen, by = lenBin)

# lengths to plot
plotLens = seq(from = 20, to = 35, by = 1)


#  ------------------------------------------------------------------------
# Create matrix to save all prop per year by method
allMethodsProps = NULL

#  ------------------------------------------------------------------------
# METHOD 1: UNIQUE ALK FOR ALL YEARS. (TO DO: select 0.5stations to simulate a 'data poor' case)

allEstAges = minEstAge:maxAge
nages = (maxAge - minEstAge) + 1 # fix number
nlen = (maxLen - minLen) + 1 # use data3 because length data has more len bins
fakeLen = seq(minLen, maxLen, by = 1)

# select just the 50% of the sampling stations per year
data2$INDEX = as.numeric(as.character(paste0(data2$YEAR, data2$STATIONID)))
data4$INDEX = as.numeric(as.character(paste0(data4$YEAR, data4$STATIONID)))
indx = as.vector(unlist(by(data = data2$INDEX, INDICES = data2$YEAR, FUN = sample, size = round(nSamLoc*0.5))))
data6 = data4[data4$INDEX %in% indx, ]

# begin analysis
tmp = data6
tmpMat = matrix(NA, ncol = nages, nrow = nlen)
  
agefac = sort(unique(tmp$AGE)) # THE FIRST AGE MUST BE 1 !!!!!(OK FOR NOW)
for(j in seq_along(agefac)){
    
  tmp2 = tmp[tmp$AGE == agefac[j], ]
  tmp3 = table(tmp2$LENGTH)
  agelen = as.numeric(names(tmp3))
  tmpMat[match(agelen, fakeLen), j] = as.numeric(tmp3)
    
}

colnames(tmpMat) = allEstAges 
rownames(tmpMat) = fakeLen
tmpMat[which(is.na(tmpMat))] = 0

# Apply age plus group:
tmpMat[ ,which(as.numeric(colnames(tmpMat)) == agePlus)] = rowSums(tmpMat[ ,which(as.numeric(colnames(tmpMat)) == agePlus):ncol(tmpMat)])
tmpMat2 = tmpMat[,1:which(as.numeric(colnames(tmpMat)) == agePlus)]
# Calculating ALK:
outAL = calc_ALK(x = tmpMat2)

plotoutAL1 = outAL[rownames(outAL) %in% plotLens,]
plotoutAL1 = melt(data = plotoutAL1)
names(plotoutAL1) = c('LENGTH', 'AGE','PROP')
plotoutAL1$METHOD = 'Method1'

# Save data frame:
allMethodsProps = rbind(allMethodsProps, plotoutAL1)

#  ------------------------------------------------------------------------
# METHOD 2: UNIQUE ALK FOR EACH YEAR.

allEstAges = minEstAge:maxAge
nages = (maxAge - minEstAge) + 1 # fix number
nlen = (maxLen - minLen) + 1 # use data3 because length data has more len bins
fakeLen = seq(minLen, maxLen, by = 1)

met2df = NULL
k = 1

  tmp = data4[data4$YEAR == allYears[k], ]
  tmpMat = matrix(NA, ncol = nages, nrow = nlen)
    
  agefac = sort(unique(tmp$AGE)) # THE FIRST AGE MUST BE 1 !!!!!(OK FOR NOW)
  for(j in seq_along(agefac)){
    
    tmp2 = tmp[tmp$AGE == agefac[j], ]
    tmp3 = table(tmp2$LENGTH)
    agelen = as.numeric(names(tmp3))
    tmpMat[match(agelen, fakeLen), j] = as.numeric(tmp3)
    
  }

  colnames(tmpMat) = allEstAges 
  rownames(tmpMat) = fakeLen
  tmpMat[which(is.na(tmpMat))] = 0

  # Apply age plus group:
  tmpMat[ ,which(as.numeric(colnames(tmpMat)) == agePlus)] = rowSums(tmpMat[ ,which(as.numeric(colnames(tmpMat)) == agePlus):ncol(tmpMat)])
  tmpMat2 = tmpMat[,1:which(as.numeric(colnames(tmpMat)) == agePlus)]
  # Calculating ALK:
  outAL = calc_ALK(x = tmpMat2)

plotoutAL2 = outAL[rownames(outAL) %in% plotLens,]
plotoutAL2 = melt(data = plotoutAL2)
names(plotoutAL2) = c('LENGTH', 'AGE','PROP')
plotoutAL2$METHOD = 'Method2'

# Save data frame:
allMethodsProps = rbind(allMethodsProps, plotoutAL2)


#  ------------------------------------------------------------------------
# 3) GAM (Lorenzo's approach) per year: AGE ~ s(LENGTH) + s(LON, LAT)

# continue script:
mydatagam = data4

# start loop over years
yearsfac = sort(unique(mydatagam$YEAR))
#saveModInd = NULL
dwriteAll = NULL
j = 1

  subdata = mydatagam[mydatagam$YEAR == yearsfac[j], ]
  data3tmp = data.frame(LON = -169.1, LAT = 59, LENGTH = plotLens) # REPEAT THIS LOCATION FOR NEXT METHOD
  
  # run the model GAM:
  age_gam = gam(AGE~s(LENGTH)+s(LON,LAT,k=10), data=subdata, family = tw, 
      method = 'REML')


  # predict data
  data3tmp$AGE = as.vector(predict(age_gam,newdata=data3tmp,type='response'))

  #Round ages (makes sense?): 
  data3tmp$AGEROUND = round(data3tmp$AGE,0)
  data3tmp$AGEROUND = ifelse(test = data3tmp$AGEROUND > agePlus, agePlus, data3tmp$AGEROUND)
  
plotoutAL3 = data3tmp[,c('LENGTH', 'AGEROUND')]
names(plotoutAL3) = c('LENGTH', 'AGE')
plotoutAL3$PROP = 1
plotoutAL3$METHOD = 'Method3'

# Save data frame:
allMethodsProps = rbind(allMethodsProps, plotoutAL3)


#  ------------------------------------------------------------------------
# 4) GAM (CRL approach) per year: PROPAGE ~ LENGTH + s(LON, LAT)

mydatagam = data4
#mydatagam = mydatagam[mydatagam$AGE >= minEstAge, ] # CUT IN firstAgeCRL ALL DATA!
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
j = 1

  tmptab = AgeYearTab[ ,colnames(AgeYearTab) == yearsfac[j]]
  tmpages = as.numeric(names(tmptab))
  tmpfreq = as.vector(tmptab)

  # plusAgeCRL = agesCRL(ages = tmpages, freq = tmpfreq, thr = 5) # 5 ind as a thr?
  plusAgeCRL = agePlus
  #print(plusAgeCRL)
  
  subdata = mydatagam[mydatagam$YEAR == yearsfac[j], ]
  data3sub =  data.frame(LON = -169.1, LAT = 59, LENGTH = plotLens) # REPEAT THIS LOCATION FOR NEXT METHOD
  
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
    
rownames(matPreds2) = plotLens
colnames(matPreds2) = 1:8
plotoutAL4 = melt(data = matPreds2)
names(plotoutAL4) = c('LENGTH', 'AGE','PROP')
plotoutAL4$METHOD = 'Method4'

# Save data frame:
allMethodsProps = rbind(allMethodsProps, plotoutAL4)

# PLOT SECTION
# ----------------------------------------------------------------------

#allMethodsProps2 = allMethodsProps[allMethodsProps$PROP > 1e-09, ]
allMethodsProps2 = allMethodsProps[allMethodsProps$AGE > 0 & allMethodsProps$AGE < 4, ]

bitmap('FinalFigures/CompareMethodsExample.tiff', height = 100, width = 100, units = 'mm', res = 800)
ggplot(allMethodsProps2, aes(LENGTH, PROP)) +
  geom_line(aes(colour = as.factor(AGE))) +
  # scale_colour_gradientn(colours = gradColors) +
  facet_wrap( ~ as.factor(METHOD), ncol = 2) +
  xlab(label = 'Length (cm)') +
  ylab(label = 'Proportion') +
  labs(colour = 'Age') +
  theme_bw()+
  theme(legend.position = c(0.1, 0.15)) 
dev.off()