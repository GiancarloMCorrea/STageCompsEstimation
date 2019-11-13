rm(list = ls())

library(foreach)
library(doParallel)
library(doSNOW)

#setwd('C:/Users/moroncog/Documents/GitHub/STageCompsEstimation')
# Number of replicates
nSim = 350

cores = detectCores()
cl = makeCluster(cores[1] - 2)
registerDoSNOW(cl)

foreach(ix = 1:nSim) %dopar% {

	# Required libraries:
	require(sp)
	library(gstat)
	require(BBmisc)
	require(ggplot2)
	require(ALKr)
	require(reshape2)
	library(mgcv)
	library(mapdata)
	library(grid)
	library(RColorBrewer)
	require(geoR)
	require(RandomFields)
	require(gridExtra)
	require(reshape)

	# call aux functions needed for the simulation:
	source('auxFunctionsSimulation.R', local = TRUE)

	# parameters for the simulation and estimation step:
	source('parametersSimulation.R', local = TRUE)

	# simulate Random Fields for recruitment
	source('simulateRandomFields.R', local = TRUE)

	# main code for simulation (population and sampling):
	source('mainSimulation3.R', local = TRUE) # simulation1 is length stratified. simulation2 is random sampling

	# estimates from the sampling output (e.g. total abundance, len abundance):
	source('estimatesSimulation2.R', local = TRUE)

	# Final Step (?): compare age props between different methods
	source('compareMethods2.R', local = TRUE)
	
}

stopCluster(cl)