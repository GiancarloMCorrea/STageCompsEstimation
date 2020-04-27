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


	# call aux functions needed for the simulation:
	source('auxFunctionsSimulation.R')

	# parameters for the simulation and estimation step:
	#source('parametersSimulation.R')
	source('parametersSimulation2.R')

	# Calculate SB_0
	source('initialConditions.R')

	# simulate Random Fields for recruitment
	source('simulateRandomFields.R')

	# main code for simulation (population and sampling):
	source('mainSimulation4.R') # simulation1 is length stratified. simulation2 is random sampling

	# estimates from the sampling output (e.g. total abundance, len abundance):
	source('estimatesSimulation2.R')

	# Final Step (?): compare age props between different methods
	#source('compareMethods2.R')
	source('studyCase_NoStrata_Simulated.R')

	# Get data base for SA models
	source('order_Data_to_SA.R')

	# Compare SA results with real population
	source('SA_simulated.R')
