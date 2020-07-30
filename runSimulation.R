setwd('C:/Users/moroncog/Documents/GitHub/STageCompsEstimation')

# Number of replicates
ix = 1

# parameters and mainsimulation3 is the classic method (in paper)
# parameters2 and mainsimulation4 is the new method (SR relation) 

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
	require(reshape)
	require(gridExtra)
	require(statmod)
	require(fields)


	# call aux functions needed for the simulation:
	source('auxFunctionsSimulation.R')

	# parameters for the simulation and estimation step:
	#source('parametersSimulation.R')
	source('parametersSimulation.R')

	# Calculate SB_0 (only for SS3 like)
	#source('initialConditions.R')

	# simulate Random Fields for recruitment
	source('simulateRandomFields.R')

	# main code for simulation (population and sampling):
	source('mainSimulation3.R') # mainSimulation3 is normal method, mainSimulation4 is SS3 like

	# estimates from the sampling output (e.g. total abundance, len abundance):
	source('estimatesSimulation2.R')

	# Final Step (?): compare age props between different methods
	source('compareMethods2.R')
