#  ------------------------------------------------------------------------
# Compare Methods from simulated data:
# IMPORTANT: DONT APPLY rm FUNCTION BECAUSE THIS SCRIPT USES INFORMATION FROM runSimulation.R
#data2 = read.csv("simData/paccod_catch_Sim.csv")
#data3 = read.csv("simData/paccod_len_Sim.csv")
#data4 = read.csv("simData/paccod_age_Sim.csv")

# read the data:
data2 = allcatchData
data3 = alllenData
data4 = allageData

# create factor for data2:
data2$STATIONID2 = as.character(data2$STATIONID)
data2$ID_HAUL = paste0(data2$YEAR, "_", data2$STATIONID2)

# create factor for data3:
data3$STATIONID2 = as.character(data3$STATIONID)
data3$ID_HAUL = paste0(data3$YEAR, "_", data3$STATIONID2)

# get number of individuals in total catch and subsample:
data5 = aggregate(data3$FREQUENCY, list(ID_HAUL = data3$ID_HAUL), sum)
data3$N_CATCH = data2$NUMBER_FISH[match(data3$ID_HAUL, data2$ID_HAUL)]
data3$N_SUBSAMPLE = data5$x[match(data3$ID_HAUL, data5$ID_HAUL)]

# get lambda:
data3$LAMBDA = data3$N_SUBSAMPLE/data3$N_CATCH

# get c^_l_i
data3$C_L_I = data3$FREQUENCY/data3$LAMBDA

