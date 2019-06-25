#  ------------------------------------------------------------------------
# Compare Methods from simulated data:
# IMPORTANT: DONT APPLY rm FUNCTION BECAUSE THIS SCRIPT USES INFORMATION FROM runSimulation.R
#data2 = read.csv("simData/paccod_catch_Sim.csv")
#data3 = read.csv("simData/paccod_len_Sim.csv")
#data4 = read.csv("simData/paccod_age_Sim.csv")

data2 = allcatchData
data3 = alllenData
data4 = allageData

# fake sex:
data3$SEX = 1
data4$SEX = 1

data2$STATIONID2 = as.character(data2$STATIONID)
data2$ID_HAUL = paste0(data2$YEAR, "_", data2$STATIONID2)

data3$STATIONID2 = as.character(data3$STATIONID)
data3$SEX2 = as.numeric(data3$SEX)
data3$ID_HAUL = paste0(data3$YEAR, "_", data3$STATIONID2)

indSt = unique(data3$ID_HAUL)
nLenHaul = table(data3$ID_HAUL)

data2$effort = areaSwept  # en km2
data2$effort = data2$effort * 100 # en ha
data2$cpueN = data2$NUMBER_FISH/data2$effort

cpueNHaul = data2$cpueN[match(indSt, data2$ID_HAUL)] # cpue in numbers
repsCpue = as.numeric(nLenHaul[match(indSt, names(nLenHaul))])
data3$cpueN = rep(cpueNHaul, times = repsCpue)

nFishLenHaul = tapply(X = data3$FREQUENCY, INDEX = data3$ID_HAUL, 
                      FUN = sum) # number fish sampled: no difference sex
nFishLenHaulOrd = as.numeric(nFishLenHaul[match(indSt, names(nFishLenHaul))])
data3$nSamples = rep(nFishLenHaulOrd, times = repsCpue) # number fish sampled

data3$PROPLEN = data3$FREQUENCY/data3$nSamples # no dimension
data3$CPUELEN = data3$PROPLEN*data3$cpueN # n/ha

areainfo = read.csv("area_stratum.csv")
areainfo$AREA_HC = areainfo$AREA_KM * 100 # in hec

data3$STRATUM_AREA = areainfo$AREA_HC[match(data3$STRATUM, areainfo$STRATUM)]

#data3$STRATUM_ALT = areainfo$STRATUM_ALT[match(data3$STRATUM, areainfo$STRATUM)]
#getAreaALT = tapply(X = areainfo$AREA_HC, INDEX = areainfo$STRATUM_ALT, FUN = sum)  
#data3$STRATUM_AREA_ALT = as.numeric(getAreaALT[match(data3$STRATUM_ALT, as.numeric(names(getAreaALT)))])

# expand len prop to stratum, sex, year
lenbase = seq(0, max(allLens), by = 1)
fakdat = data.frame(len = lenbase, cpue = NA)

outdat = NULL

uniYr = unique(data3$YEAR)

for(i in seq_along(uniYr)){
  
  tmp1 = data3[data3$YEAR == uniYr[i], ]
  
  #create base matrix, comes from data2: it is complete! (counting zeros)
  tempdat_sx1 = data.frame(YEAR = uniYr[i],
                       STRATUM = rep(x = data2[data2$YEAR == uniYr[i], "STRATUM"],
                                     each = length(lenbase)),
                       STATIONID2 = rep(x = data2[data2$YEAR == uniYr[i], "STATIONID2"],
                                  each = length(lenbase)),
                       LENGTH = lenbase,
                       CPUELEN = 0,
                       SEX2 = 1)
  tempdat = tempdat_sx1
  tempdat$STATIONID2 = as.character(tempdat$STATIONID2)
  
  uniStr = sort(unique(tmp1$STRATUM)) # modify if it is necessary
  
  for(j in seq_along(uniStr)){
    
    tmp2 = tmp1[tmp1$STRATUM == uniStr[j], ] # modify if it is necessary
    tempdat2 = tempdat[tempdat$STRATUM == uniStr[j], ]
    
    stratumArea = unique(tmp2$STRATUM_AREA) # must be only 1 number
    # length(stratumArea)
    
    uniSx = unique(tmp2$SEX2)
    
#    for(k in seq_along(uniSx)){
      
      tmp3 = tmp2[tmp2$SEX2 == uniSx[1], ]
      tmp3 = tmp3[, c("YEAR", "STRATUM", "STATIONID2", "LENGTH",
                      "CPUELEN", "SEX2")]
      tempdat3 = tempdat2[tempdat2$SEX2 == uniSx[1], ]
      
      joindata = merge(tmp3, tempdat3, by = c("YEAR", "STRATUM", "STATIONID2", "LENGTH",
                                              "SEX2"), all = T)
      joindata$CPUELEN.x[which(is.na(joindata$CPUELEN.x))] = 0 #se usa cpuelen.x 
      joindata$CPUELEN.y = NULL
      
      cpueLENTAB = tapply(X = joindata$CPUELEN.x, INDEX = joindata$LENGTH, 
                          FUN = mean) # change function if it is required
      fishlengths = as.numeric(names(cpueLENTAB))
      cpuelengths = as.numeric(cpueLENTAB)
      tmp4 = data.frame(len = fishlengths, cpue = cpuelengths)
      
      #merge stand base
      # fakdat2 = merge(fakdat, tmp4, by = "len", all = T)
      # fakdat2$cpue.y[which(is.na(fakdat2$cpue.y))] = 0 
      tmp4$lenabun = tmp4$cpue * stratumArea
      
      lenabund2 = tmp4$lenabun
      tmp5 = data.frame(year = uniYr[i],
                        stratum = uniStr[j],
                        len = tmp4$len,
                        lencpue = tmp4$cpue,
                        lenabun = lenabund2,
                        sex = uniSx[1])
      outdat = rbind(outdat, tmp5)
      
#    }
    
  }
  #print(i)
}

outdat$lenabun = outdat$lenabun/1e+6 # in millons ind

# YEAHHHH!!!
# STRATUM 2: (grouping)
outdat$stratum2 = areainfo$STRATUM_ALT[match(outdat$stratum, areainfo$STRATUM)]

outdat2 = NULL

uniYr = unique(outdat$year)

for(i in seq_along(uniYr)){
  
  tmp1 = outdat[outdat$year == uniYr[i], ]
  
  uniStr = sort(unique(outdat$stratum2)) # modify if it is necessary
  
  for(j in seq_along(uniStr)){
    
    tmp2 = tmp1[tmp1$stratum2 == uniStr[j], ] # modify if it is necessary

    uniSx = unique(tmp2$sex)
    
#    for(k in seq_along(uniSx)){
      
      tmp3 = tmp2[tmp2$sex == uniSx[1], ]
      
      tmp4 = tapply(X = tmp3$lenabun, INDEX = tmp3$len, FUN = sum)
      
      tmp5 = data.frame(year = uniYr[i],
                        stratum = uniStr[j],
                        len = as.numeric(names(tmp4)),
                        lenabun = as.numeric(tmp4),
                        sex = uniSx[1])
      outdat2 = rbind(outdat2, tmp5)
      
#    }
    
  }
  #print(i)
}

if(!simulation){
	write.csv(outdat, "simData/abunlen_stratum1Sim.csv", row.names = F)
	write.csv(outdat2, "simData/abunlen_stratum2Sim.csv", row.names = F)
}