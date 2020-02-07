# create files
setwd('C:/Users/moroncog/Documents/GitHub/STageCompsEstimation')
scenarioName = 'NoS_NoT'


require(r4ss)
dir.create(path = 'simSA/M1', showWarnings = FALSE)
dir.create(path = 'simSA/M2', showWarnings = FALSE)
dir.create(path = 'simSA/M3', showWarnings = FALSE)
dir.create(path = 'simSA/M4', showWarnings = FALSE)

# delete files if they exist
unlink(file.path('simSA/M1', list.files(file.path('simSA/M1'))))
unlink(file.path('simSA/M2', list.files(file.path('simSA/M2'))))
unlink(file.path('simSA/M3', list.files(file.path('simSA/M3'))))
unlink(file.path('simSA/M4', list.files(file.path('simSA/M4'))))

# read base files:
rawPath = getwd()
starter = list.files('simSA/base')[1]
control = list.files('simSA/base')[2]
datax = list.files('simSA/base')[3]
forecast = list.files('simSA/base')[4]

file.copy(file.path(paste0(rawPath, '/simSA/base/', c(starter, control, forecast))), paste0(rawPath, '/simSA/M1'))
file.copy(file.path(paste0(rawPath, '/simSA/base/', c(starter, control, forecast))), paste0(rawPath, '/simSA/M2'))
file.copy(file.path(paste0(rawPath, '/simSA/base/', c(starter, control, forecast))), paste0(rawPath, '/simSA/M3'))
file.copy(file.path(paste0(rawPath, '/simSA/base/', c(starter, control, forecast))), paste0(rawPath, '/simSA/M4'))


# Read created data:
SAsimdata = read.csv(paste0('simData/SAmodel_data_', scenarioName, '.csv'))
SAsimdata0 = read.csv(paste0('simData/SAmodel_catlen_', scenarioName, '.csv'))
SAsimdata2 = read.csv(paste0('simData/AbundanceEstSim_', scenarioName, '.csv'))

# Read data file from base folder
datainfo = SS_readdat(file.path('simSA/base', datax), version = '3.30')

# Create catch data:
catchdata = data.frame(year = SAsimdata$years, seas = 1, fleet = 1, catch = SAsimdata$catch_biom*1e-03, catch_se = 0.01)
catchdata$years[1] = -999
# create index data:
indexdata = data.frame(year = SAsimdata2$YEAR, seas = 1, index = 2, obs = SAsimdata2$x*1e-2, se_log = 0.1)
# create len composition data fishery:
lendata0 = data.frame(Yr = SAsimdata$years[-1], Seas = 7, FltSvy = 1, FltSvy = 0, Part = 0,  
					Nsamp = 150)
SAsimdata0 = SAsimdata0[-1, 4:120]
lendata = cbind(lendata0, SAsimdata0)

#replace data:
datainfo$catch = catchdata
datainfo$CPUE = indexdata
datainfo$lencomp = lendata


# create age data:
agedata_1 = data.frame(Yr = SAsimdata2$YEAR, Seas = 1, FltSvy = 2, FltSvy = 0, Part = 0, Ageerr = 1, Lbin_lo = 1, lbin_hi = 120, 
					Nsamp = 1400)

	# For method 1:
	agematrix1 = read.csv(paste0('simData/Method1_pcod_', scenarioName, '.csv'))
	agedata1 = cbind(agedata_1, agematrix1[, c(2:ncol(agematrix1))])

	# For method 2:
	agematrix2 = read.csv(paste0('simData/Method2_pcod_', scenarioName, '.csv'))
	agedata2 = cbind(agedata_1, agematrix2[, c(2:ncol(agematrix2))])
	
	# For method 3:
	agematrix3 = read.csv(paste0('simData/Method3_pcod_', scenarioName, '.csv'))
	agedata3 = cbind(agedata_1, agematrix3[, c(2:ncol(agematrix3))])

	# For method 4:
	agematrix4 = read.csv(paste0('simData/Method4_pcod_', scenarioName, '.csv'))
	agedata4 = cbind(agedata_1, agematrix4[, c(2:ncol(agematrix4))])

# replace age comp data and write file
datainfo$agecomp = agedata1
SS_writedat(datlist = datainfo, outfile = file.path(rawPath, 'simSA/M1/simDAT.dat'), version = '3.30', overwrite = TRUE)

datainfo$agecomp = agedata2
SS_writedat(datlist = datainfo, outfile = file.path(rawPath, 'simSA/M2/simDAT.dat'), version = '3.30', overwrite = TRUE)

datainfo$agecomp = agedata3
SS_writedat(datlist = datainfo, outfile = file.path(rawPath, 'simSA/M3/simDAT.dat'), version = '3.30', overwrite = TRUE)

datainfo$agecomp = agedata4
SS_writedat(datlist = datainfo, outfile = file.path(rawPath, 'simSA/M4/simDAT.dat'), version = '3.30', overwrite = TRUE)


# Run models

setwd('C:/Users/moroncog/Documents/GitHub/STageCompsEstimation/simSA/M1')
system(command = "ss", show.output.on.console = FALSE)
#setwd('..')
setwd('C:/Users/moroncog/Documents/GitHub/STageCompsEstimation/simSA/M2')
system(command = "ss", show.output.on.console = FALSE)

setwd('C:/Users/moroncog/Documents/GitHub/STageCompsEstimation/simSA/M3')
system(command = "ss", show.output.on.console = FALSE)

setwd('C:/Users/moroncog/Documents/GitHub/STageCompsEstimation/simSA/M4')
system(command = "ss", show.output.on.console = FALSE)
