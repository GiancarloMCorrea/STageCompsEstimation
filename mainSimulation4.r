# Calculate vector:
Rover_y = rnorm(n = length(allYears), mean = 0, sd = sigmaR)


Z_par = M_par + (F_par*SelecFish)
lambdaP = 1 - exp(-Z_par)
#ageStrucAll = matrix(NA, nrow = length(allYears), ncol = length(allAges))#matrix to save the TOTAL age structure at the end of the first half.
NageStrucGrid = array(NA, dim = c(nrow(predictGrid2), length(minAge:maxAge)-1, length(allYears)))#matrix to save the age structure at the end of the sampling (t1 time)
LageStrucGrid = array(NA, dim = c(nrow(predictGrid2), length(minAge:maxAge)-1, length(allYears)))#matrix to save the length at age values at the end of the sampling (t1 time)
NageStrucGridSam = array(NA, dim = c(nrow(predictGrid2), length(minAge:maxAge), length(allYears)))#matrix to save the age structure at the end of the sampling (t1 time)
LageStrucGridSam = array(NA, dim = c(nrow(predictGrid2), length(minAge:maxAge), length(allYears)))#matrix to save the length at age values at the end of the sampling (t1 time)
cat_len_tot = matrix(NA, ncol = length(allLens), nrow = length(allYears))


# Simulation --------------------------------------------------------------
allcatchData = NULL # to save catch data
alllenData = NULL # to save len data
allageData = NULL # to save age data
catch_year = rep(0, times = length(allYears))
catch_year_biom = rep(0, times = length(allYears))
ssb_year = rep(0, times = length(allYears))
ssb_year_1 = rep(0, times = length(allYears))
r_year = rep(0, times = length(allYears))
Biom_year = rep(0, times = length(allYears))
for(k in seq_along(allYears)){

	# Calculate recruitment for this year
	if(k == 1){
		R_y = R0
	} else {
		iniNs = c(Rgrid[j], nTmp4)
		iniLens = c(lenatage0, lTmp4)
		R_y = ((4*steep*R0*SB_y_1)/(SB_0*(1-steep)+SB_y_1*(5*steep - 1)))*exp(0.5*b_y[k]*sigmaR^2+Rover_y[k])
	}

	# Calculate indicator per grid 
	Epsilon2 = exp(Epsilon1[,k])/sum(exp(Epsilon1[,k]))
	Rgrid = R_y*Epsilon2

	# define age sample locations. ALL RANDOM.
	ageLocations =  sample(x = sampleStations$sampledGrids, size = nSamLoc, replace = FALSE)

	# length counting in samples: max should be 25
	# tmpLenCount = 0

	# create vector of sample for age
	# nFishAgeSampled = numeric(length(allLens))
	SB_y_j = numeric(length(yy2@grid.index))
	SB_y_1_j = numeric(length(yy2@grid.index))
	Cat_j = numeric(length(yy2@grid.index))
	Cat_biom_j = numeric(length(yy2@grid.index))
	Biom_j = numeric(length(yy2@grid.index))
	cat_len_tot_y = numeric(length(allLens))
  for(j in seq_along(yy2@grid.index)){

	if(k == 1){
		# Initial conditions:    
		R0grid = Rgrid[j]/GridArea # density (N/km2)
		iniNstmp   = R0grid*exp(-M_par*(minAge:(maxAgePopMod-1))) # take care: Z mortality
		iniNsPlus = sum(iniNstmp[(maxAge+1-minAge):length(iniNstmp)]) + (iniNstmp[length(iniNstmp)]*exp(-M_par))/(1-exp(-M_par))
		iniNs = c(iniNstmp[1:(maxAge-minAge)], iniNsPlus) # Inital abundance from minAge to maxAge (plus group)

		#Initial lengths
		iniLens = ifelse(allAges <= A1_par, Lminp + (bpar*allAges), Linf+(L1_par-Linf)*exp(-(K_par+KparT[k]+yy2$sim1[j])*(allAges-A1_par))) # SS growth, Age vs Len relationship
		a_new = maxAge:(2*maxAge)
		iniLens[length(iniLens)] = sum(exp(-0.2*(a_new - maxAge))*(iniLens[length(iniLens)] + ((a_new - maxAge)/maxAge)*(Linf-iniLens[length(iniLens)])))/sum(exp(-0.2*(a_new-maxAge)))
		lenatage0 = iniLens[1]

	} else {

		iniNs = c(Rgrid[j]/GridArea, NageStrucGrid[j,,(k-1)]) # density
		iniLens = c(lenatage0, LageStrucGrid[j,,(k-1)])

	}

    # Start of the season: calculate matrix age length and sample of lengths and ages, also SSB:
    nTmp = iniNs
    lTmp = iniLens
    
    # Lenght at the middle of the season
    lTmp2 = lTmp + (lTmp-Linf)*(exp(-dT*t1*(K_par+KparT[k]+yy2$sim1[j])) - 1)


    # SD len at the start of season : 
    sdTmp = numeric(length(lTmp))
	sdTmp = ifelse(allAges <= A1_par, CV1, (CV1 + ((lTmp - L1_par)/(Linf - L1_par))*(CV2-CV1))) # SS growth, Age vs Len relationship
	sdTmp = ifelse(allAges >= A2_par, CV2, sdTmp) # SS growth, Age vs Len relationship


    # SD len after first half: 
    sdTmp2 = numeric(length(lTmp))
	sdTmp2 = ifelse(allAges <= A1_par, CV1, (CV1 + ((lTmp2 - L1_par)/(Linf - L1_par))*(CV2-CV1))) # SS growth, Age vs Len relationship
	sdTmp2 = ifelse(allAges >= A2_par, CV2, sdTmp2) # SS growth, Age vs Len relationship

    # Create Age Length matrix at the start of the season:
	
		AgeLenMatrixProp0 = matrix(NA, ncol = length(allAges), nrow = length(allLens))
		
		for(i in 1:nrow(AgeLenMatrixProp0)){
		
		  if(i == 1){
		  
			#Lminp = minLen - lenBin*0.5 # 0.5 because I am working with 1 cm bin
			Fac1 = (Lminp - lTmp)/sdTmp
			AgeLenMatrixProp0[i, ] = pnorm(Fac1)
		  
		  }
		  if(i == length(allLens)){
			
			Lmaxp = maxLen - lenBin*0.5
			Fac1 = (Lmaxp - lTmp)/sdTmp
			AgeLenMatrixProp0[i, ] = 1 - pnorm(Fac1)
		
		  } else {
			
			Ll1p = allLens[i] + lenBin*0.5
			Llp = allLens[i] - lenBin*0.5
			Fac1 = (Ll1p - lTmp)/sdTmp
			Fac2 = (Llp - lTmp)/sdTmp
			AgeLenMatrixProp0[i, ] = pnorm(Fac1) - pnorm(Fac2)
			
		  }
			
		}

		# Create Age Length matrix at the middle of the season:
	
		AgeLenMatrixProp = matrix(NA, ncol = length(allAges), nrow = length(allLens))
		
		for(i in 1:nrow(AgeLenMatrixProp)){
		
		  if(i == 1){
		  
			#Lminp = minLen - lenBin*0.5 # 0.5 because I am working with 1 cm bin
			Fac3 = (Lminp - lTmp2)/sdTmp2
			AgeLenMatrixProp[i, ] = pnorm(Fac3)
		  
		  }
		  if(i == length(allLens)){
			
			Lmaxp = maxLen - lenBin*0.5
			Fac3 = (Lmaxp - lTmp2)/sdTmp2
			AgeLenMatrixProp[i, ] = 1 - pnorm(Fac3)
		
		  } else {
			
			Ll1p = allLens[i] + lenBin*0.5
			Llp = allLens[i] - lenBin*0.5
			Fac3 = (Ll1p - lTmp2)/sdTmp2
			Fac4 = (Llp - lTmp2)/sdTmp2
			AgeLenMatrixProp[i, ] = pnorm(Fac3) - pnorm(Fac4)
			
		  }
			
		}
		
		# with matrix at the start of the season
		f_a0 = colSums(sweep(AgeLenMatrixProp0, MARGIN=1, Part2Fec, `*`)) # fecundity per age A.1.18
		SB_y_j[j] = sum(fracFem*nTmp*GridArea*f_a0*allMats) # only females in numbers


	if(j %in% sampleStations$sampledGrids & allYears[k] %in% allYearsSam){

		AgeLenMatrix = sweep(AgeLenMatrixProp, MARGIN=2, nTmp, `*`) # abundance per length per age

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
								STRATUM3 = sampleStations$stratum3[sampleStations$sampledGrids == j],
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
										STRATUM3 = sampleStations$stratum3[sampleStations$sampledGrids == j],
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
										STRATUM3 = sampleStations$stratum3[sampleStations$sampledGrids == j],
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
										STRATUM3 = sampleStations$stratum3[sampleStations$sampledGrids == j],
										TYPEGRID = sampleStations$typegrid[sampleStations$sampledGrids == j],
										LENGTH = lensSam, AGE = agesSam)
						allageData = rbind(allageData, ageData)
					}

		}

	}
    

    # Continue after sampling:
    nTmp2	= nTmp

	  # calculate total biomass at the start of the season:
		Biom_j[j] = sum(wt_a*nTmp2*GridArea) # catch in biomass


    if(allYears[k] == iniYear0){
	    # Calculating catch: (CHECK THIS LATER, IT IS NOT THE SAME AS MANUAL)
	    cat_age = rep(0, times = length(allAges))
	    Cat_j[j] = sum(cat_age) # in numbers

	    cat_len = rep(0, times = length(allLens))
	    cat_len_tot_y = cat_len_tot_y + cat_len

	    Cat_biom_j[j] = 0 # catch in biomass


		# save the age structure per grid
		#ageStrucMidYear = ageStrucMidYear + nTmp3 # sum the age structure at each grid. at the end of the grid loop I will have the total population structure
		
	    # THESE TWO CALCULATIONS WILL BE PASSED TO THE NEXT YEAR: (y+1)-----------
	    nTmp4 = nTmp2*exp(-M_par) # THIS SHOULD BE N after fishing

    } else {

    		    # Calculating catch: (CHECK THIS LATER, IT IS NOT THE SAME AS MANUAL)
	    cat_age = (F_par/Z_par)*(nTmp2*GridArea*SelecFish)*lambdaP
	    Cat_j[j] = sum(cat_age) # in numbers

	    cat_len = rowSums(sweep(AgeLenMatrixProp, MARGIN=2, cat_age, `*`))
	    cat_len_tot_y = cat_len_tot_y + cat_len

	    wt_a = colSums(sweep(AgeLenMatrixProp, MARGIN=1, allWts, `*`)) # weight per age A.1.18
	    Cat_biom_j[j] = sum(wt_a*cat_age) # catch in biomass

	    # calculate total biomass at the start of the season:
		Biom_j[j] = sum(wt_a*nTmp2*GridArea) # catch in biomass

		# save the age structure per grid
		#ageStrucMidYear = ageStrucMidYear + nTmp3 # sum the age structure at each grid. at the end of the grid loop I will have the total population structure
		
	    # THESE TWO CALCULATIONS WILL BE PASSED TO THE NEXT YEAR: (y+1)-----------
	    nTmp4 = nTmp2*exp(-Z_par) # THIS SHOULD BE N after fishing
	
    }


    nTmp5 = nTmp4[-length(nTmp4)]
    nTmp5[length(nTmp5)] = nTmp5[length(nTmp5)-1]*exp(-Z_par[length(nTmp5)-1]) +  nTmp5[length(nTmp5)]*exp(-Z_par[length(nTmp5)])# N_y+1 equation A.1.20

    # Length for time y+1 (Equation A.1.10)
    lTmp4 = lTmp + (lTmp-Linf)*(exp(-(K_par+KparT[k]+yy2$sim1[j])) - 1)
    lTmp5 = lTmp4[-length(lTmp4)]
    lTmp5[length(lTmp5)] = (nTmp5[length(nTmp5)-1]*lTmp2[length(lTmp2)] + nTmp5[length(nTmp5)]*(lTmp4[length(lTmp5)] + (lTmp5[length(lTmp5)] - Linf)*(exp(K_par+KparT[k]+yy2$sim1[j])-1)))/(nTmp5[length(nTmp5)-1] + nTmp5[length(nTmp5)])
    # ---------------------------



   # SD len at the end: 
    sdTmp3 = numeric(length(lTmp4))
	sdTmp3 = ifelse(allAges <= A1_par, CV1, (CV1 + ((lTmp4 - L1_par)/(Linf - L1_par))*(CV2-CV1))) # SS growth, Age vs Len relationship
	sdTmp3 = ifelse(allAges >= A2_par, CV2, sdTmp3) # SS growth, Age vs Len relationship

    # Create Age Length matrix at the end of the season:
	
		AgeLenMatrixProp2 = matrix(NA, ncol = length(allAges), nrow = length(allLens))
		
		for(i in 1:nrow(AgeLenMatrixProp2)){
		
		  if(i == 1){
		  
			#Lminp = minLen - lenBin*0.5 # 0.5 because I am working with 1 cm bin
			Fac1 = (Lminp - lTmp4)/sdTmp3
			AgeLenMatrixProp2[i, ] = pnorm(Fac1)
		  
		  }
		  if(i == length(allLens)){
			
			Lmaxp = maxLen - lenBin*0.5
			Fac1 = (Lmaxp - lTmp4)/sdTmp3
			AgeLenMatrixProp2[i, ] = 1 - pnorm(Fac1)
		
		  } else {
			
			Ll1p = allLens[i] + lenBin*0.5
			Llp = allLens[i] - lenBin*0.5
			Fac1 = (Ll1p - lTmp4)/sdTmp3
			Fac2 = (Llp - lTmp4)/sdTmp3
			AgeLenMatrixProp2[i, ] = pnorm(Fac1) - pnorm(Fac2)
			
		  }
			
		}

    # Calculate SSB at the end of the season (or start of the next year)
    f_a2 = colSums(sweep(AgeLenMatrixProp2, MARGIN=1, Part2Fec, `*`)) # fecundity per age A.1.18
	SB_y_1_j[j] = sum(fracFem*nTmp4*GridArea*f_a2*allMats) # only females in numbers

	# Save data:
	NageStrucGrid[j,,k] = nTmp5
	LageStrucGrid[j,,k] = lTmp5
	
	NageStrucGridSam[j,,k] = nTmp # Structure at the start of the year
	LageStrucGridSam[j,,k] = lTmp # Length at the start of the year
	    
  }
  
  catch_year[k] = sum(Cat_j)
  catch_year_biom[k] = sum(Cat_biom_j)

  SB_y = sum(SB_y_j)
  ssb_year[k] = SB_y

  SB_y_1 = sum(SB_y_1_j)
  ssb_year_1[k] = SB_y_1

  r_year[k] = R_y

  Biom_year[k] = sum(Biom_j)

  cat_len_tot[k,] = cat_len_tot_y/sum(cat_len_tot_y)

}


# END OF LOOP 


dir.create('simData', showWarnings = FALSE)

if(ix == 1){
	write.csv(allcatchData, paste0('simData/paccod_catch_Sim_', scenarioName, '.csv'), row.names = FALSE)
	write.csv(alllenData, paste0('simData/paccod_len_Sim_', scenarioName, '.csv'), row.names = FALSE)
	write.csv(allageData, paste0('simData/paccod_age_Sim_', scenarioName, '.csv'), row.names = FALSE)
	write.csv(data.frame(years = allYears, catch = catch_year, catch_biom = catch_year_biom, ssb = ssb_year, ssb_1 = ssb_year_1, 
		R = r_year, biomass = Biom_year), paste0('simData/SAmodel_data_', scenarioName, '.csv'), row.names = FALSE)
	write.csv(cat_len_tot, paste0('simData/SAmodel_catlen_', scenarioName, '.csv'), row.names = FALSE)
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
              color_low = "blue", color_high = "red", zeroiswhite = TRUE, xlim = c(-179,-158), ylim = c(53.5,63)) +
			  geom_polygon(data = ak, aes(long, lat, group = group), 
			  fill = 8, color="black") +
			  xlab('longitude') +
			  ylab('latitude') +
			  #xlim(-180,-156) +
			  theme(legend.position = 'none') +
			  theme(plot.margin = unit(c(0,0,0,0),"cm"))
	
	
	ax2 = ggplot(lenStations, aes(lon, lat)) +
			geom_point(size = 0.9) +
			#geom_point(data = ageStations, aes(lon, lat), col = 'red', shape = 2) +
			#geom_point(data = ageStations2, aes(lon, lat), col = 'blue', shape = 2) +
			geom_polygon(data = ak, aes(long, lat, group = group), 
				  fill = 8, color="black") +
			theme_bw() +
			xlab(' ') +
			ylab('latitude') +
			coord_fixed(ratio=1, xlim=c(-179,-158), ylim=c(53.5,63)) +
			theme(plot.margin = unit(c(0,0,0,0),"cm"))


	bitmap(paste0('surveyDescription3R_', scenarioName, '.tiff'), height = 110, width = 90, units = 'mm', res = 500)
	
	grid.arrange(ax2, ax1, nrow = 2)
			
	dev.off()  
}


if(ix == 1){

	ax3 = map.heatmap2(lat = yy2@coords[,2], lon = yy2@coords[,1], data = yy2@data, data2 = lenStations,
              color_low = "blue", color_high = "red", zeroiswhite = TRUE, xlim = c(-179,-158), ylim = c(53.5,63), pSize = 0.65) +
			  geom_polygon(data = ak, aes(long, lat, group = group), 
			  fill = 8, color="black") +
			  xlab('longitude') +
			  ylab('latitude') +
			  theme(legend.position = c(0.5, 1.18), plot.margin = unit(c(0,0,0,0),"cm"), legend.key.width = unit(0.75, "cm"), 
			  	legend.text=element_text(size=7.5))

	bitmap(paste0('surveyDescription4R_', scenarioName, '.tiff'), height = 55, width = 90, units = 'mm', res = 500)
	
	print(ax3)
			
	dev.off()  
}
