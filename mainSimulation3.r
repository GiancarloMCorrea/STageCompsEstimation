
#L'min in SS equations:
Lminp = minLen - lenBin*0.5 # 0.5 because I am working with cm
bpar = (L1_par-Lminp)/A1_par

Z_par = M_par + (F_par*SelecFish)
#ageStrucAll = matrix(NA, nrow = length(allYears), ncol = length(allAges))#matrix to save the TOTAL age structure at the end of the first half.
NageStrucGrid = array(NA, dim = c(nrow(predictGrid2), length(allAges), length(allYears)))#matrix to save the age structure at the end of the sampling (t1 time)
LageStrucGrid = array(NA, dim = c(nrow(predictGrid2), length(allAges), length(allYears)))#matrix to save the length at age values at the end of the sampling (t1 time)
NageStrucGridSam = array(NA, dim = c(nrow(predictGrid2), length(allAges), length(allYears)))#matrix to save the age structure at the end of the sampling (t1 time)
LageStrucGridSam = array(NA, dim = c(nrow(predictGrid2), length(allAges), length(allYears)))#matrix to save the length at age values at the end of the sampling (t1 time)


#currentDate = format(Sys.time(), "%b %d %Y %X")
#currentDate = gsub(pattern = ' ', replacement = '_', x = currentDate)
#currentDate = gsub(pattern = ':', replacement = '', x = currentDate)

# Simulation --------------------------------------------------------------
allcatchData = NULL # to save catch data
alllenData = NULL # to save len data
allageData = NULL # to save age data
for(k in seq_along(allYears)){

	# Rec for this year:
	R0year = iniR0*exp(rRecTemp[k])
	
	# Rec in density terms
	R0inidengrid = log(R0year/StudyArea) # Nfish/km2: mean density in log scale

	# define age sample locations. ALL RANDOM.
	ageLocations =  sample(x = sampleStations$sampledGrids, size = nSamLoc, replace = FALSE)

	# length counting in samples: max should be 25
	# tmpLenCount = 0

	# create vector of sample for age
	# nFishAgeSampled = numeric(length(allLens))

  for(j in seq_along(yy2@grid.index)){

	if(k == 1){
		# Initial conditions:    
		R0grid = exp(R0inidengrid + Omega1[j] + Epsilon1[j,k]) # add spatial and spatiotemporal structure
		iniNs   = R0grid*exp(-Z_par*allAges) # take care: Z mortality
		#iniLens = ifelse(allAges <= A1_par, Lminp + (bpar*allAges), Linf+(L1_par-Linf)*exp(-(K_par+KparT[k]+yy2$sim1[j])*(allAges-A1_par))) # SS growth, Age vs Len relationship
		iniLens = Linf*(1-exp(-(K_par+KparT[k]+yy2$sim1[j])*(allAges-t0))) # VB growth
		lenatage0 = iniLens[1]
	} else {
		iniNs = toNewYear(vec = NageStrucGrid[j,,(k-1)], firstVal = exp(R0inidengrid + Omega1[j] + Epsilon1[j,k])) # define what is R0x
		iniLens = toNewYear(vec = LageStrucGrid[j,,(k-1)], firstVal = lenatage0)
	}

    # Save state of the population at time t
    nTmp = iniNs
    lTmp = iniLens
    
    # First half: Natural mortality: time step t1*dT
    nTmp2	= nTmp*exp(-dT*t1*Z_par)
    # Individual growth: time step t1*dT
    lTmp2 = lTmp + (lTmp-Linf)*(exp(-dT*t1*(K_par+KparT[k]+yy2$sim1[j])) - 1)
    
    # SD len after first half: 
    sdTmp2 = numeric(length(lTmp2))
    # for(l in seq_along(sdTmp2)){
      
      # if(allAges[l] <= A1_par) sdTmp2[l] = CV1
      # if(allAges[l] > A1_par & allAges[l] < A2_par) sdTmp2[l] = (CV1 + ((lTmp2[l] - lTmp[2])/(lTmp[(maxAge - 1)] - lTmp[2]))*(CV2-CV1))
      # if(allAges[l] >= A2_par) sdTmp2[l] = CV2
      
    # }
    
#	sdTmp2 = ifelse(allAges <= A1_par, CV1, (CV1 + ((lTmp2 - L1_par)/(Linf - L1_par))*(CV2-CV1))) # SS growth, Age vs Len relationship
#	sdTmp2 = ifelse(allAges >= A2_par, CV2, sdTmp2) # SS growth, Age vs Len relationship

	sdTmp2 = CV1 + ((lTmp2 - L1_par)/(Linf - L1_par))*(CV2-CV1) # SS growth, Age vs Len relationship


		nTmp3 = nTmp2
		lTmp3 = lTmp2
#    	sdTmp3 = sdTmp2

    
    # Create Age Length matrix after first half: (just for sampled grids)
	
	if(j %in% sampleStations$sampledGrids & allYears[k] %in% allYearsSam){
	
		AgeLenMatrixProp = matrix(NA, ncol = length(allAges), nrow = length(allLens))
		
		for(i in 1:nrow(AgeLenMatrixProp)){
		
		  if(i == 1){
		  
			#Lminp = minLen - lenBin*0.5 # 0.5 because I am working with 1 cm bin
			Fac1 = (Lminp - lTmp2)/sdTmp2
			AgeLenMatrixProp[i, ] = pnorm(Fac1)
		  
		  }
		  if(i == length(allLens)){
			
			Lmaxp = maxLen - lenBin*0.5
			Fac1 = (Lmaxp - lTmp2)/sdTmp2
			AgeLenMatrixProp[i, ] = 1 - pnorm(Fac1)
		
		  } else {
			
			Ll1p = allLens[i] + lenBin*0.5
			Llp = allLens[i] - lenBin*0.5
			Fac1 = (Ll1p - lTmp2)/sdTmp2
			Fac2 = (Llp - lTmp2)/sdTmp2
			AgeLenMatrixProp[i, ] = pnorm(Fac1) - pnorm(Fac2)
			
		  }
			
		}
		
		AgeLenMatrix = sweep(AgeLenMatrixProp, MARGIN=2, nTmp2, `*`)
		
# Here: sampling strategy:
# --------------------
		# follow Thorson and Haltuch 2019 method: delta model (see equation 16). and Thorson (2018): select length and age at the same time
		SelAbun = sweep(AgeLenMatrix, MARGIN=2, SelecSurv, `*`)
		SelAbun2 = rowSums(SelAbun)
		pL = 1 - exp(-areaSwept*SelAbun)
		pLsim = structure(vapply(pL, rbinom, numeric(1), n = 1, size = 1), dim=dim(pL))
		rL = (areaSwept*SelAbun)/pL # as poisson delta link model. add over pL if it is necessary.
		findNAN = which(is.nan(rL)|is.infinite(rL))
		rL[findNAN] = 1 # to avoid warnings
		#rLsim = structure(vapply(log(rL)-(sigmaM^2)/2, rlnorm, numeric(1), n = 1, sdlog = sigmaM), dim=dim(rL)) # lognormal
		randomNumbers = structure(rnorm(nrow(rL)*ncol(rL), mean = 0, sd = sigmaM), dim=dim(rL))
		rLsim = structure(vapply(rL*exp(randomNumbers), rpois, numeric(1), n = 1), dim=dim(rL)) # poisson
		roundAgeLenSampled = pLsim*rLsim # get round sampled matrix
		nFishLenSampled = rowSums(roundAgeLenSampled) # here the 1s are deleted.
		catchStation = sum(nFishLenSampled)
# --------------------

		# check sampling: (catch data)
		# nFishSampled = sum(roundAgeLenSampled)
		nFishSampled = catchStation
		catchData = data.frame(YEAR = allYears[k], STATIONID = j, START_LONGITUDE = sampleStations$lon[sampleStations$sampledGrids == j], 
							   START_LATITUDE = sampleStations$lat[sampleStations$sampledGrids == j], 
								STRATUM_ALT = sampleStations$stratum[sampleStations$sampledGrids == j],
								STRATUM = sampleStations$stratum2[sampleStations$sampledGrids == j],
								TYPEGRID = sampleStations$typegrid[sampleStations$sampledGrids == j],
								NUMBER_FISH = nFishSampled)
		allcatchData = rbind(allcatchData, catchData)

		# check sampling: (len data): apply function. max th
		
		if(nFishSampled > 0){ # just for positive catches
			# For catch sample less than th:
			if(nFishSampled <= lenCatchTh){
				nFishLenSampled2 = nFishLenSampled
				posLenSam = which(nFishLenSampled2 > 0)
				lenData = data.frame(YEAR = allYears[k], STATIONID = j, LON = sampleStations$lon[sampleStations$sampledGrids == j], 
									   LAT = sampleStations$lat[sampleStations$sampledGrids == j], 
										STRATUM_ALT = sampleStations$stratum[sampleStations$sampledGrids == j],
										STRATUM = sampleStations$stratum2[sampleStations$sampledGrids == j],
										TYPEGRID = sampleStations$typegrid[sampleStations$sampledGrids == j],
										LENGTH = allLens[posLenSam], FREQUENCY = nFishLenSampled2[posLenSam])
				alllenData = rbind(alllenData, lenData)
			}
			# For catch sample greater than th:
			if(nFishSampled > lenCatchTh){
				nFishLenSampled3 = sample(x = rep(allLens, times = nFishLenSampled), size = lenCatchTh, replace = FALSE)
				prev2 = table(nFishLenSampled3)
				nFishLenSampled2 = numeric(length(allLens))
				nFishLenSampled2[allLens %in% as.numeric(names(prev2))] = as.vector(prev2)
			
				posLenSam = which(nFishLenSampled2 > 0)
				lenData = data.frame(YEAR = allYears[k], STATIONID = j, LON = sampleStations$lon[sampleStations$sampledGrids == j], 
									   LAT = sampleStations$lat[sampleStations$sampledGrids == j], 
										STRATUM_ALT = sampleStations$stratum[sampleStations$sampledGrids == j],
										STRATUM = sampleStations$stratum2[sampleStations$sampledGrids == j],
										TYPEGRID = sampleStations$typegrid[sampleStations$sampledGrids == j],
										LENGTH = allLens[posLenSam], FREQUENCY = nFishLenSampled2[posLenSam])
				alllenData = rbind(alllenData, lenData)
			}
			
		}

		# check sampling: (age data): simple strategy: max 25 ind per length IN THE SURVEY. max 50 ind sampled in a station. max 3 ind per length in A STATION.
		
		# ALL AREA ALL STATIONS: RANDOM SAMPLING. 
		if(j %in% ageLocations & nFishSampled > 0){ # no run this part if the n fish sampled for length  = 0
		
			# choose length to be age sampled:
			if(maxNSamAge > sum(nFishLenSampled2)){
				lFishSampled = rep(allLens, times = nFishLenSampled2) # max ind per set for age sample = 40. FIRST CONDITIONAL				
			} else {
				lFishSampled = sample(x = rep(allLens, times = nFishLenSampled2), size = maxNSamAge, replace = FALSE) # max ind per set for age sample = 40. FIRST CONDITIONAL
			}
			
			prev = table(lFishSampled)
			nFishToAge = numeric(length(allLens))
			nFishToAge[allLens %in% as.numeric(names(prev))] = as.vector(prev)
			
					# now, choose age at length to be sampled:
					if(max(nFishToAge) > 0){ # just if there are more lengths to be aged
						agesSam = ageSample(mat = roundAgeLenSampled, vec = nFishToAge)
						lensSam = rep(as.numeric(names(prev)), times = as.numeric(prev)) # length in order

						ageData = data.frame(YEAR = allYears[k], STATIONID = j, LON = sampleStations$lon[sampleStations$sampledGrids == j], 
									   LAT = sampleStations$lat[sampleStations$sampledGrids == j], 
										STRATUM_ALT = sampleStations$stratum[sampleStations$sampledGrids == j],
										STRATUM = sampleStations$stratum2[sampleStations$sampledGrids == j],
										TYPEGRID = sampleStations$typegrid[sampleStations$sampledGrids == j],
										LENGTH = lensSam, AGE = agesSam)
						allageData = rbind(allageData, ageData)
					}

		}
    
		nTmp3 = nTmp2
		lTmp3 = lTmp2

	}
    
	# save the age structure per grid
	#ageStrucMidYear = ageStrucMidYear + nTmp3 # sum the age structure at each grid. at the end of the grid loop I will have the total population structure
	
    # After the survey sample simulation: run the second half of dT:
    nTmp4 = nTmp3*exp(-dT*(1-t1)*Z_par)
    lTmp4 = lTmp3 + (lTmp3-Linf)*(exp(-dT*(1-t1)*(K_par+KparT[k]+yy2$sim1[j])) - 1) # Individual growth: time step 0.5*dT
    
	NageStrucGrid[j,,k] = nTmp4
	LageStrucGrid[j,,k] = lTmp4
	
	NageStrucGridSam[j,,k] = nTmp3
	LageStrucGridSam[j,,k] = lTmp3
	    
  }
  
  #ageStrucAll[k, ] = ageStrucMidYear
  # print(k)
}

if(ix == 1){
	write.csv(allcatchData, paste0('simData/paccod_catch_Sim_', scenarioName, '.csv'), row.names = FALSE)
	write.csv(alllenData, paste0('simData/paccod_len_Sim_', scenarioName, '.csv'), row.names = FALSE)
	write.csv(allageData, paste0('simData/paccod_age_Sim_', scenarioName, '.csv'), row.names = FALSE)
}

# needed for compareMethods.R:
NAgeYearMatrix = t(apply(X = NageStrucGridSam, MARGIN = c(2,3), FUN = sum)*GridArea)
if(ix == 1){
	write.csv(NAgeYearMatrix, paste0('simData/NAgeYearMat_', scenarioName, '.csv'), row.names = FALSE)
}

#save(NageStrucGridSam, file = 'simData/simAgeStructure.RData')
#save(LageStrucGridSam, file = 'simData/simLengthAtAge.RData')

# Plot survey map (for the last year = 2016):

timeSurvey = data.frame(lon = yy2@coords[,1], lat = yy2@coords[,2], time = 1:nrow(yy2@coords))
lenStations = data.frame(lon = yy2@coords[sampleStations$sampledGrids, 1], lat = yy2@coords[sampleStations$sampledGrids, 2])
ageStations = data.frame(lon = yy2@coords[ageLocations, 1], lat = yy2@coords[ageLocations, 2])
#ageStations2 = data.frame(lon = yy2@coords[ageLocations2, 1], lat = yy2@coords[ageLocations2, 2])

# png('surveyDescription.png', height = 700, width = 700, units = 'px', res = 110)
# print(ggplot(timeSurvey, aes(lon, lat)) +
        # geom_point(aes(color = time), size = 1.5) +
		# geom_point(data = lenStations, aes(lon, lat), col = 1) +
		# geom_point(data = ageStations1, aes(lon, lat), col = 4, shape = 2) +
		# geom_point(data = ageStations2, aes(lon, lat), col = 5, shape = 2) +
        # scale_colour_gradientn(colours = gradColors2) +
		# theme_bw())
# dev.off()

if(ix == 1){

	ax1 = map.heatmap(lat = yy2@coords[,2], lon = yy2@coords[,1], yy2@data,
              color_low = "blue", color_high = "red", zeroiswhite = TRUE, xlim = c(-179,-158), ylim = c(54,62.5)) +
			  geom_polygon(data = ak, aes(long, lat, group = group), 
			  fill = 8, color="black") +
			  xlab('longitude') +
			  ylab('latitude') +
			  #xlim(-180,-156) +
			  theme(legend.position = 'none') +
			  theme(plot.margin = unit(c(0,0,0,0),"cm"))
	
	
	ax2 = ggplot(lenStations, aes(lon, lat)) +
			geom_point(size = 1) +
			geom_point(data = ageStations, aes(lon, lat), col = 'red', shape = 2) +
			#geom_point(data = ageStations2, aes(lon, lat), col = 'blue', shape = 2) +
			geom_polygon(data = ak, aes(long, lat, group = group), 
				  fill = 8, color="black") +
			theme_bw() +
			xlab('longitude') +
			ylab('latitude') +
			coord_fixed(ratio=1, xlim=c(-179,-158), ylim=c(54,62.5)) +
			theme(plot.margin = unit(c(0,0,0,0),"cm"))


	bitmap(paste0('surveyDescription3R_', scenarioName, '.tiff'), height = 130, width = 120, units = 'mm', res = 600)
	
	grid.arrange(ax2, ax1, nrow = 2)
			
	dev.off()  
}