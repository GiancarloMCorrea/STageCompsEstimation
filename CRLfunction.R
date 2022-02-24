# This is a function to estimate age proportions based on a set of explanatory variables. 

# It performs a CRL model using the 'gam' function (mgcv package) for estimation. 
# It is designed to incorporate different minimum and maximum ages for different years. 
# Moreover, it is run by year or other temporal variable (e.g. quarter)

# Advice: Check your data before running this function. 
# Make sure that you have the right values for each variable in AgeSubsample (age subsample) and LengthSubsample (length subsample).
# See the manuscript for details: Correa et al. (2020)

# Arguments:

# AgeSubsample = dataset (data.frame), which includes all variables used in the FormulaGAM specification and an 'Age' column.
# LengthSubsample = dataset (data.frame), used for prediction based on FormulaGAM. 
#                  Its column names should match all variables used in FormulaGAM (explanatory variables).      
# FormulaGAM = character object. Formula of ONLY explanatory variables to be passed to the 'gam' function. (e.g. 'LENGTH + s(LON, LAT)')
# AgeMin = an integer or vector (numeric). Minimum age to be estimated. If a integer is specified, it will be used for all TimeVariable classes ('unique'). 
#          If a vector is specified, it should be the same length of all classes in TimeVariable. (i.e. length(AgeMin) = 5 if length(unique(TimeVariable)) = 5). 
#          Moreover, if a vector is specified, it will follow the same order as sort(unique(TimeVariable)).
# AgeMax = an integer or vector (numeric). Maximum age to be estimated. It follows the same logic of AgeMin. length(AgeMin) and length(AgeMax) do not need to have the same length.
# AgeVariable = character. Column name of 'AgeSubsample' to be used as the age variable (e.g. 'AGE')
# TimeVariable = character. Column name of 'AgeSubsample' to be used as the temporal variable (i.e. the model will be run for each time) (e.g. 'YEAR')
# ... = further arguments for passing to 'gam' function.


estimateAgeCRL = function(AgeSubsample, LengthSubsample, FormulaGAM, AgeMin = 1, AgeMax = 8, AgeVariable, TimeVariable, ...){

  require(mgcv)
  require(plyr)
  # Cheking variables and arguments:
  if(missing(AgeVariable)) stop('AgeVariable must be provided')
  if(missing(TimeVariable)) print('No TimeVariable provided. Model is run using the whole dataset')
  
  if(!missing(TimeVariable)){
    if(!(TimeVariable %in% names(AgeSubsample))) stop('TimeVariable is not a column name in AgeSubsample')
    if(!(TimeVariable %in% names(LengthSubsample))) stop('TimeVariable is not a column name in LengthSubsample')
  }

  nAgeMin = length(AgeMin)
  nAgeMax = length(AgeMax)
  nTimeVariable = length(unique(AgeSubsample[ ,TimeVariable]))

  if(nAgeMin > 1) {
    if(nTimeVariable != nAgeMin) stop('AgeMin does not have the same length as the TimeVariable classes')
  }

  if(nAgeMax > 1) {
    if(nTimeVariable != nAgeMax) stop('AgeMax does not have the same length as the TimeVariable classes')
  }


  # Ordering the data:
  #Data2 = Data
  if(missing(TimeVariable)) {
    TimeVariable = 'Time'
    AgeSubsample[ ,TimeVariable] = 1
    LengthSubsample[ ,TimeVariable] = 1
  }

  TimeFac = sort(unique(AgeSubsample[ ,TimeVariable]))
  FormulaGAM2 = as.formula(paste0('Pi_a ~ ', FormulaGAM))

  nAgesAll = max(AgeMax) - min(AgeMin) + 1
  finalReport = data.frame(matrix(ncol = ncol(LengthSubsample) + nAgesAll, nrow = 0))
  colNames2 = c(colnames(LengthSubsample), min(AgeMin):max(AgeMax))
  colnames(finalReport) = colNames2
  ageMatrix_final = NULL

  # Begin loop over time variable:
    for(j in seq_along(TimeFac)){

      if(nAgeMin == 1) AgeMin = rep(AgeMin, times = nTimeVariable)
      if(nAgeMax == 1) AgeMax = rep(AgeMax, times = nTimeVariable)

      subdata = AgeSubsample[AgeSubsample[ ,TimeVariable] == TimeFac[j], ]
      data3sub = LengthSubsample[LengthSubsample[ ,TimeVariable] == TimeFac[j], ]
      
      subdata[,AgeVariable] = ifelse(test = subdata[,AgeVariable] > AgeMax[j], AgeMax[j], subdata[,AgeVariable]) 

      # all ages for that year
      allages = AgeMin[j]:AgeMax[j]

      # for each year
    	PropAgeMat = matrix(NA, ncol = length(allages), nrow = 1) 
    	colnames(PropAgeMat) = allages

      # run the model GAM:
      matPreds = matrix(NA, ncol = length(allages), nrow = nrow(data3sub))
      for(ii in seq_along(allages)){
        
        if(ii == length(allages)){
          predvals = rep(1, times = nrow(data3sub))
        } else {
          
          subdata$Pi_a = ifelse(test = subdata[,AgeVariable] > allages[ii], 0, 1)
          modtmp = gam(FormulaGAM2, family = binomial, data = subdata, method = "ML", ...)
          predtmp = predict(modtmp, newdata = data3sub, type = 'response')
          predvals = as.vector(predtmp)
    	    elimina = which(subdata$Pi_a == 1)
    	  if(length(elimina) > 0) {
    		  subdata = subdata[-which(subdata$Pi_a == 1), ]
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
        
      
    matPreds3 = as.data.frame(matPreds2)
    colnames(matPreds3) = allages
    tmpReport = cbind(data3sub, matPreds3)
    finalReport = rbind.fill(finalReport, tmpReport)
    ageMatrix_final = rbind(ageMatrix_final, matPreds3)
        
    }

  return(list(merged_data = finalReport, prop_age_matrix = ageMatrix_final))

}


