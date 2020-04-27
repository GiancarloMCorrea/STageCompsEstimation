# Read SA models
require(r4ss)

setwd('C:/Users/moroncog/Documents/GitHub/STageCompsEstimation')

iniYear0 = 1976
endYear = 2018

mod1 = SS_output(dir = 'simSA/M1', covar = FALSE)
mod2 = SS_output(dir = 'simSA/M2', covar = FALSE)
mod3 = SS_output(dir = 'simSA/M3', covar = FALSE)
mod4 = SS_output(dir = 'simSA/M4', covar = FALSE)

SS_plots(mod1)
SS_plots(mod2)
SS_plots(mod3)
SS_plots(mod4)

# Read real population:

realdatats = read.csv(paste0('simData/SAmodel_data_', scenarioName, '.csv'))
realdataages = read.csv(paste0('simData/NAgeYearMat_', scenarioName, '.csv'))
realdataages = sweep(realdataages, MARGIN=1, rowSums(realdataages), `/`)
realdataages = as.matrix(realdataages)


plot(realdatats$years, realdatats$biomass*1e-03, type = 'l', col = 1, lty = 2, lwd = 2)
lines(mod1$timeseries$Yr, mod1$timeseries$Bio_all, type = 'l', col = 1)
lines(mod2$timeseries$Yr, mod2$timeseries$Bio_all, type = 'l', col = 2)
lines(mod3$timeseries$Yr, mod3$timeseries$Bio_all, type = 'l', col = 3)
lines(mod4$timeseries$Yr, mod4$timeseries$Bio_all, type = 'l', col = 4)
legend('topright', legend = c('SA1', 'SA2','SA3','SA4', 'Real'), col = c(1,2,3,4,1), lty = c(1,1,1,1,2))

plot(realdatats$years, realdatats$R*1e-03, type = 'l', col = 1, lty = 2, lwd = 2)
lines(mod1$timeseries$Yr, mod1$timeseries$Recruit_0, type = 'l', col = 1)
lines(mod2$timeseries$Yr, mod2$timeseries$Recruit_0, type = 'l', col = 2)
lines(mod3$timeseries$Yr, mod3$timeseries$Recruit_0, type = 'l', col = 3)
lines(mod4$timeseries$Yr, mod4$timeseries$Recruit_0, type = 'l', col = 4)
legend('topright', legend = c('SA1', 'SA2','SA3','SA4', 'Real'), col = c(1,2,3,4,1), lty = c(1,1,1,1,2))


sum(abs(mod1$timeseries$Recruit_0[mod1$timeseries$Yr %in% iniYear0:endYear] - realdatats$R*1e-03))
sum(abs(mod2$timeseries$Recruit_0[mod2$timeseries$Yr %in% iniYear0:endYear] - realdatats$R*1e-03))
sum(abs(mod3$timeseries$Recruit_0[mod3$timeseries$Yr %in% iniYear0:endYear] - realdatats$R*1e-03))
sum(abs(mod4$timeseries$Recruit_0[mod4$timeseries$Yr %in% iniYear0:endYear] - realdatats$R*1e-03))


# Compare ages
selMod = mod1
SAages = selMod$natage[selMod$natage[,11] == 'B', ]
SAages = SAages[SAages$Yr %in% iniYear0:endYear, ]
SAages = SAages[ ,13:ncol(SAages)]
SAages = sweep(SAages, MARGIN=1, rowSums(SAages), `/`)
SAages = as.matrix(SAages)
diff1 = abs(SAages - realdataages)

selMod = mod2
SAages = selMod$natage[selMod$natage[,11] == 'B', ]
SAages = SAages[SAages$Yr %in% iniYear0:endYear, ]
SAages = SAages[ ,13:ncol(SAages)]
SAages = sweep(SAages, MARGIN=1, rowSums(SAages), `/`)
SAages = as.matrix(SAages)
diff2 = abs(SAages - realdataages)

selMod = mod3
SAages = selMod$natage[selMod$natage[,11] == 'B', ]
SAages = SAages[SAages$Yr %in% iniYear0:endYear, ]
SAages = SAages[ ,13:ncol(SAages)]
SAages = sweep(SAages, MARGIN=1, rowSums(SAages), `/`)
SAages = as.matrix(SAages)
diff3 = abs(SAages - realdataages)

selMod = mod4
SAages = selMod$natage[selMod$natage[,11] == 'B', ]
SAages = SAages[SAages$Yr %in% iniYear0:endYear, ]
SAages = SAages[ ,13:ncol(SAages)]
SAages = sweep(SAages, MARGIN=1, rowSums(SAages), `/`)
SAages = as.matrix(SAages)
diff4 = abs(SAages - realdataages)


mod1$estimated_non_dev_parameters['VonBert_K_Fem_GP_1', 'Value']
mod2$estimated_non_dev_parameters['VonBert_K_Fem_GP_1', 'Value']
mod3$estimated_non_dev_parameters['VonBert_K_Fem_GP_1', 'Value']
mod4$estimated_non_dev_parameters['VonBert_K_Fem_GP_1', 'Value']

mod1$estimated_non_dev_parameters['L_at_Amax_Fem_GP_1', 'Value']
mod2$estimated_non_dev_parameters['L_at_Amax_Fem_GP_1', 'Value']
mod3$estimated_non_dev_parameters['L_at_Amax_Fem_GP_1', 'Value']
mod4$estimated_non_dev_parameters['L_at_Amax_Fem_GP_1', 'Value']

mod1$estimated_non_dev_parameters['NatM_p_1_Fem_GP_1', 'Value']
mod2$estimated_non_dev_parameters['NatM_p_1_Fem_GP_1', 'Value']
mod3$estimated_non_dev_parameters['NatM_p_1_Fem_GP_1', 'Value']
mod4$estimated_non_dev_parameters['NatM_p_1_Fem_GP_1', 'Value']


mod1$estimated_non_dev_parameters['SR_LN(R0)', 'Value']
mod2$estimated_non_dev_parameters['SR_LN(R0)', 'Value']
mod3$estimated_non_dev_parameters['SR_LN(R0)', 'Value']
mod4$estimated_non_dev_parameters['SR_LN(R0)', 'Value']


