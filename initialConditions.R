
#L'min in SS equations:
Lminp = minLen - lenBin*0.5 # 0.5 because I am working with cm
bpar = (L1_par-Lminp)/A1_par

	# bias adjustement for SR:
	b_y = rep(0, times = length(allYears))
	b_y[which(allYears <= y_3 & allYears >= y_2)] = bmax
	selYear1 = which(allYears < y_2 & allYears > y_1)
	b_y[selYear1] = bmax*(1-((allYears[selYear1] - y_1)/(y_2-y_1)))
	selYear2 = which(allYears < y_4 & allYears > y_3)
	b_y[selYear2] = bmax*(1-((y_3 - allYears[selYear2])/(y_4-y_3)))

		# Initial conditions:    
		iniNstmp   = R0*exp(-M_par*(minAge:(maxAgePopMod-1))) # take care: Z mortality
		iniNsPlus = sum(iniNstmp[(maxAge+1-minAge):length(iniNstmp)]) + (iniNstmp[length(iniNstmp)]*exp(-M_par))/(1-exp(-M_par))
		iniNs = c(iniNstmp[1:(maxAge-minAge)], iniNsPlus) # Inital abundance from minAge to maxAge (plus group)

		#Initial lengths
		iniLens = ifelse(allAges <= A1_par, Lminp + (bpar*allAges), Linf+(L1_par-Linf)*exp(-(K_par)*(allAges-A1_par))) # SS growth, Age vs Len relationship
		a_new = maxAge:(2*maxAge)
		iniLens[length(iniLens)] = sum(exp(-0.2*(a_new - maxAge))*(iniLens[length(iniLens)] + ((a_new - maxAge)/maxAge)*(Linf-iniLens[length(iniLens)])))/sum(exp(-0.2*(a_new-maxAge)))

    # Save state of the population at time t
    nTmp = iniNs
    lTmp = iniLens
    
    # Continue init:
    nTmp2	= nTmp

    # Lenght at the middle of the season
    lTmp2 = lTmp + (lTmp-Linf)*(exp(-dT*t1*(K_par)) - 1)

    # SD len after first half: 
    sdTmp2 = numeric(length(lTmp2))
    
	sdTmp2 = ifelse(allAges <= A1_par, CV1, (CV1 + ((lTmp2 - L1_par)/(Linf - L1_par))*(CV2-CV1))) # SS growth, Age vs Len relationship
	sdTmp2 = ifelse(allAges >= A2_par, CV2, sdTmp2) # SS growth, Age vs Len relationship

		nTmp3 = nTmp2
		lTmp3 = lTmp2

    
    # Create Age Length matrix after first half:
	
		AgeLenMatrixProp = matrix(NA, ncol = length(allAges), nrow = length(allLens))
		
		for(i in 1:nrow(AgeLenMatrixProp)){
		
		  if(i == 1){
		  
			#Lminp = minLen - lenBin*0.5 # 0.5 because I am working with 1 cm bin
			Fac1 = (Lminp - lTmp2)/sdTmp2
			AgeLenMatrixProp[i, ] = pnorm(Fac1)
		  
		  }
		  if(i == length(allLens)){
			
			Lmaxp = maxLen - lenBin*lenBin*0.5
			Fac1 = (Lmaxp - lTmp2)/sdTmp2
			AgeLenMatrixProp[i, ] = 1 - pnorm(Fac1)
		
		  } else {
			
			Ll1p = allLens[i] + lenBin*lenBin*0.5
			Llp = allLens[i] - lenBin*lenBin*0.5
			Fac1 = (Ll1p - lTmp2)/sdTmp2
			Fac2 = (Llp - lTmp2)/sdTmp2
			AgeLenMatrixProp[i, ] = pnorm(Fac1) - pnorm(Fac2)
			
		  }
			
		}
		
		f_a = colSums(sweep(AgeLenMatrixProp, MARGIN=1, Part2Fec, `*`)) # fecundity per age A.1.18
		SB_0 = sum(fracFem*nTmp*f_a) # only females
