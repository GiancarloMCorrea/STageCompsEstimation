# -------------------------------------------------------------------------
# 1) This is the first part of the study case section
# Run a model with length data and compare numbers at age with estimated age comps by method

setwd(dir = 'C:/Users/moroncog/Documents/GitHub/STageCompsEstimation')
require(r4ss)
require(reshape2)
require(ggplot2)

#  ------------------------------------------------------------------------
rm(list = ls())

modelName1 = 'StockAssessmentModels/Method0_6'

# read the model and check results
tmp1 = SS_output(dir = modelName1)


sel1 = tmp1$ageselex[tmp1$ageselex$Fleet == 2, ]

nest1 = tmp1$natage[tmp1$natage$`Beg/Mid` == 'M', ]
rownames(nest1) = nest1$Yr
nest1sub = nest1[20:42, 13:33]
nest1sub2 = sweep(nest1sub, MARGIN = 2, as.vector(as.matrix(sel1[1,8:28])), `*`) # est comps
nest1sub3 = rowSums(nest1sub2[ , 1:2])
nest1sub4 = rowSums(nest1sub2[ , 9:21])
nest1sub5 = cbind(nest1sub3, nest1sub2[ , 3:8], nest1sub4)
nest1sub6 = sweep(nest1sub5, MARGIN = 1, rowSums(nest1sub5), `/`) # est comps
names(nest1sub6) = c('1', '2', '3', '4', '5', '6', '7', '8')

est1 = read.csv('realData/Method1_pcod.csv', row.names = 1)
est2 = read.csv('realData/Method2_pcod.csv', row.names = 1)
est3 = read.csv('realData/Method3_pcod.csv', row.names = 1)
est4 = read.csv('realData/Method4_pcod.csv', row.names = 1)

names(est1) = c('1', '2', '3', '4', '5', '6', '7', '8')
names(est2) = c('1', '2', '3', '4', '5', '6', '7', '8')
names(est3) = c('1', '2', '3', '4', '5', '6', '7', '8')
names(est4) = c('1', '2', '3', '4', '5', '6', '7', '8')

# MSEtotal:
sum((nest1sub6 - est1)^2)
sum((nest1sub6 - est2)^2)
sum((nest1sub6 - est3)^2)
sum((nest1sub6 - est4)^2)

# Plot:
nest1sub6$YEAR = rownames(nest1sub6)
df0 = melt(nest1sub6, id = c('YEAR'))
df0$Method = 'SSLen'
est1$YEAR = rownames(est1)
df1 = melt(est1, id = c('YEAR'))
df1$Method = 'pooled ALK'
est2$YEAR = rownames(est2)
df2 = melt(est2, id = c('YEAR'))
df2$Method = 'annual ALK'
est3$YEAR = rownames(est3)
df3 = melt(est3, id = c('YEAR'))
df3$Method = 'GAM'
est4$YEAR = rownames(est4)
df4 = melt(est4, id = c('YEAR'))
df4$Method = 'CRL'

mergedf = rbind(df0, df1, df2, df3, df4)
colnames(mergedf) = c('YEAR', 'AGE', 'FREQUENCY', 'Method')
mergedf$AGE = as.numeric(as.character(mergedf$AGE))

bitmap('FinalFigures/compare_estimatedSS3_allMethods.tiff', height = 190, width = 190, units = 'mm', res = 900)
	print(ggplot(mergedf, aes(x = AGE, y = FREQUENCY)) +
	  geom_line(aes(linetype = factor(Method), color = factor(Method))) +
	  facet_wrap( ~ factor(YEAR), nrow = 5) +
	  xlab('Age') +
	  ylab('Proportion of abundance') +
	  scale_x_discrete(limits = 1:8) +
	  scale_linetype_manual(values=c("solid", "solid", "solid", "solid", "twodash")) + 
	  scale_color_manual(values=c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF", "black")) +
	  theme_bw() +
	  theme(legend.position = c(0.8, 0.08), legend.title = element_blank(), legend.text = element_text(size = 10)))
 dev.off()

mergedf = rbind(df1, df2, df3, df4)
colnames(mergedf) = c('YEAR', 'AGE', 'FREQUENCY', 'Method')
mergedf$AGE = as.numeric(as.character(mergedf$AGE))
mergedf$Method = factor(mergedf$Method, levels = c('pooled ALK', 'annual ALK', 'GAM', 'CRL')) # THESE ARE TRUE VALUES

bitmap('FinalFigures/compare_estimatedSS3_allMethods2.tiff', height = 190, width = 190, units = 'mm', res = 900)
	print(ggplot(mergedf, aes(x = AGE, y = FREQUENCY)) +
	  geom_line(aes(linetype = factor(Method), color = factor(Method))) +
	  facet_wrap( ~ factor(YEAR), nrow = 5) +
	  xlab('Age') +
	  ylab('Proportion of abundance') +
	  scale_x_discrete(limits = 1:8) +
	  scale_linetype_manual(values=c("solid", "solid", "solid", "solid")) + 
	  scale_color_manual(values=c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")) +
	  theme_bw() +
	  theme(legend.position = c(0.8, 0.08), legend.title = element_blank(), legend.text = element_text(size = 10)))
 dev.off()


# -------------------------------------------------------------------------
# 2) This is the second part of the study case section
# Run a model with length and age data and compare numbers at age among length model and (length and age) model

setwd(dir = 'C:/Users/moroncog/Documents/Year1_OSUPhD/Simulations_SAmodels/StockAssessmentModels2')
require(r4ss)
require(reshape2)
require(ggplot2)

#  ------------------------------------------------------------------------
rm(list = ls())

modelName1 = 'StockAssessmentModels/Method0_6'

# read the model and check results
tmp1 = SS_output(dir = modelName1)
sel1 = tmp1$ageselex[tmp1$ageselex$Fleet == 2, ]
nest1 = tmp1$natage[tmp1$natage$`Beg/Mid` == 'M', ]
rownames(nest1) = nest1$Yr
nest1sub = nest1[20:42, 13:33]
nest1sub2 = sweep(nest1sub, MARGIN = 2, as.vector(as.matrix(sel1[1,8:28])), `*`) # est comps
nest1sub3 = rowSums(nest1sub2[ , 1:2])
nest1sub4 = rowSums(nest1sub2[ , 9:21])
nest1sub5 = cbind(nest1sub3, nest1sub2[ , 3:8], nest1sub4)
nest1sub6 = sweep(nest1sub5, MARGIN = 1, rowSums(nest1sub5), `/`) # est comps
names(nest1sub6) = c('1', '2', '3', '4', '5', '6', '7', '8')

# read the model 1 and check results
modelName1 = 'Method1_6'
tmp1 = SS_output(dir = modelName1)
sel1 = tmp1$ageselex[tmp1$ageselex$Fleet == 2, ]
nest1 = tmp1$natage[tmp1$natage$`Beg/Mid` == 'M', ]
rownames(nest1) = nest1$Yr
nest1sub = nest1[20:42, 13:33]
nest1sub2 = sweep(nest1sub, MARGIN = 2, as.vector(as.matrix(sel1[1,8:28])), `*`) # est comps
nest1sub3 = rowSums(nest1sub2[ , 1:2])
nest1sub4 = rowSums(nest1sub2[ , 9:21])
nest1sub5 = cbind(nest1sub3, nest1sub2[ , 3:8], nest1sub4)
est1 = sweep(nest1sub5, MARGIN = 1, rowSums(nest1sub5), `/`) # est comps
names(est1) = c('1', '2', '3', '4', '5', '6', '7', '8')

# read the model 2 and check results
modelName1 = 'Method2_6'
tmp1 = SS_output(dir = modelName1)
sel1 = tmp1$ageselex[tmp1$ageselex$Fleet == 2, ]
nest1 = tmp1$natage[tmp1$natage$`Beg/Mid` == 'M', ]
rownames(nest1) = nest1$Yr
nest1sub = nest1[20:42, 13:33]
nest1sub2 = sweep(nest1sub, MARGIN = 2, as.vector(as.matrix(sel1[1,8:28])), `*`) # est comps
nest1sub3 = rowSums(nest1sub2[ , 1:2])
nest1sub4 = rowSums(nest1sub2[ , 9:21])
nest1sub5 = cbind(nest1sub3, nest1sub2[ , 3:8], nest1sub4)
est2 = sweep(nest1sub5, MARGIN = 1, rowSums(nest1sub5), `/`) # est comps
names(est2) = c('1', '2', '3', '4', '5', '6', '7', '8')

# read the model 3 and check results
modelName1 = 'Method3_6'
tmp1 = SS_output(dir = modelName1)
sel1 = tmp1$ageselex[tmp1$ageselex$Fleet == 2, ]
nest1 = tmp1$natage[tmp1$natage$`Beg/Mid` == 'M', ]
rownames(nest1) = nest1$Yr
nest1sub = nest1[20:42, 13:33]
nest1sub2 = sweep(nest1sub, MARGIN = 2, as.vector(as.matrix(sel1[1,8:28])), `*`) # est comps
nest1sub3 = rowSums(nest1sub2[ , 1:2])
nest1sub4 = rowSums(nest1sub2[ , 9:21])
nest1sub5 = cbind(nest1sub3, nest1sub2[ , 3:8], nest1sub4)
est3 = sweep(nest1sub5, MARGIN = 1, rowSums(nest1sub5), `/`) # est comps
names(est3) = c('1', '2', '3', '4', '5', '6', '7', '8')

# read the model 4 and check results
modelName1 = 'Method4_6'
tmp1 = SS_output(dir = modelName1)
sel1 = tmp1$ageselex[tmp1$ageselex$Fleet == 2, ]
nest1 = tmp1$natage[tmp1$natage$`Beg/Mid` == 'M', ]
rownames(nest1) = nest1$Yr
nest1sub = nest1[20:42, 13:33]
nest1sub2 = sweep(nest1sub, MARGIN = 2, as.vector(as.matrix(sel1[1,8:28])), `*`) # est comps
nest1sub3 = rowSums(nest1sub2[ , 1:2])
nest1sub4 = rowSums(nest1sub2[ , 9:21])
nest1sub5 = cbind(nest1sub3, nest1sub2[ , 3:8], nest1sub4)
est4 = sweep(nest1sub5, MARGIN = 1, rowSums(nest1sub5), `/`) # est comps
names(est4) = c('1', '2', '3', '4', '5', '6', '7', '8')

# Plot:
nest1sub6$YEAR = rownames(nest1sub6)
df0 = melt(nest1sub6, id = c('YEAR'))
df0$Method = 'SS3Len'

est1$YEAR = rownames(est1)
df1 = melt(est1, id = c('YEAR'))
df1$Method = 'SSAge1'
est2$YEAR = rownames(est2)
df2 = melt(est2, id = c('YEAR'))
df2$Method = 'SSAge2'
est3$YEAR = rownames(est3)
df3 = melt(est3, id = c('YEAR'))
df3$Method = 'SSAge3'
est4$YEAR = rownames(est4)
df4 = melt(est4, id = c('YEAR'))
df4$Method = 'SSAge4'

mergedf = rbind(df1, df2, df3, df4)
colnames(mergedf) = c('YEAR', 'AGE', 'FREQUENCY', 'Method')
mergedf$AGE = as.numeric(as.character(mergedf$AGE))

bitmap('compare_estimatedSS3_allMethods3.tiff', height = 190, width = 190, units = 'mm', res = 900)
	print(ggplot(mergedf, aes(x = AGE, y = FREQUENCY)) +
	  geom_line(aes(linetype = factor(Method), color = factor(Method))) +
	  facet_wrap( ~ factor(YEAR), nrow = 5) +
	  xlab('Age') +
	  ylab('Proportion of abundance') +
	  scale_x_discrete(limits = 1:8) +
	  scale_linetype_manual(values=c("solid", "solid", "solid", "solid")) + 
	  scale_color_manual(values=c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")) +
	  theme_bw() +
	  theme(legend.position = c(0.8, 0.08), legend.title = element_blank(), legend.text = element_text(size = 10)))
 dev.off()


# -------------------------------------------------------------------------
# 3) This is the second part of the study case section
# Compare the model with the current matrix and 

 ###################### WATCH OUTTTTTTTTTTTTTTTTTT!!!!!!!!!!!!!!!!!
 # FINAL RESULTS ARE IN OTHER FOLDER 

setwd(dir = 'C:/Users/moroncog/Documents/Year1_OSUPhD/Simulations_SAmodels')
require(r4ss)
require(reshape2)
require(ggplot2)
require(gridExtra)

#  ------------------------------------------------------------------------
rm(list = ls())

modelName0 = 'StockAssessmentModels/Method0_6'
modelName1 = 'StockAssessmentModels/Method1_6'
modelName2 = 'StockAssessmentModels/Method2_6'
modelName3 = 'StockAssessmentModels/Method3_6'
modelName4 = 'StockAssessmentModels/Method4_6'

# read the model and check results
tmp0 = SS_output(dir = modelName0)
tmp1 = SS_output(dir = modelName1)
tmp2 = SS_output(dir = modelName2)
tmp3 = SS_output(dir = modelName3)
tmp4 = SS_output(dir = modelName4)

SS_plots(replist = tmp0)
SS_plots(replist = tmp1)
SS_plots(replist = tmp2)
SS_plots(replist = tmp3)
SS_plots(replist = tmp4)



tmp0$likelihoods_used[1,1]
tmp1$likelihoods_used[1,1]
tmp2$likelihoods_used[1,1]
tmp3$likelihoods_used[1,1]
tmp4$likelihoods_used[1,1]


tmp0$parameters['VonBert_K_Fem_GP_1', 'Value']
tmp1$parameters['VonBert_K_Fem_GP_1', 'Value']
tmp2$parameters['VonBert_K_Fem_GP_1', 'Value']
tmp3$parameters['VonBert_K_Fem_GP_1', 'Value']
tmp4$parameters['VonBert_K_Fem_GP_1', 'Value']

tmp0$parameters['VonBert_K_Fem_GP_1', 'Parm_StDev']
tmp1$parameters['VonBert_K_Fem_GP_1', 'Parm_StDev']
tmp2$parameters['VonBert_K_Fem_GP_1', 'Parm_StDev']
tmp3$parameters['VonBert_K_Fem_GP_1', 'Parm_StDev']
tmp4$parameters['VonBert_K_Fem_GP_1', 'Parm_StDev']


tmp0$parameters['L_at_Amax_Fem_GP_1', 'Value']
tmp1$parameters['L_at_Amax_Fem_GP_1', 'Value']
tmp2$parameters['L_at_Amax_Fem_GP_1', 'Value']
tmp3$parameters['L_at_Amax_Fem_GP_1', 'Value']
tmp4$parameters['L_at_Amax_Fem_GP_1', 'Value']

tmp0$parameters['L_at_Amax_Fem_GP_1', 'Parm_StDev']
tmp1$parameters['L_at_Amax_Fem_GP_1', 'Parm_StDev']
tmp2$parameters['L_at_Amax_Fem_GP_1', 'Parm_StDev']
tmp3$parameters['L_at_Amax_Fem_GP_1', 'Parm_StDev']
tmp4$parameters['L_at_Amax_Fem_GP_1', 'Parm_StDev']



tmp0$parameters['SD_young_Fem_GP_1', 'Value']
tmp1$parameters['SD_young_Fem_GP_1', 'Value']
tmp2$parameters['SD_young_Fem_GP_1', 'Value']
tmp3$parameters['SD_young_Fem_GP_1', 'Value']
tmp4$parameters['SD_young_Fem_GP_1', 'Value']

tmp0$parameters['SD_young_Fem_GP_1', 'Parm_StDev']
tmp1$parameters['SD_young_Fem_GP_1', 'Parm_StDev']
tmp2$parameters['SD_young_Fem_GP_1', 'Parm_StDev']
tmp3$parameters['SD_young_Fem_GP_1', 'Parm_StDev']
tmp4$parameters['SD_young_Fem_GP_1', 'Parm_StDev']


tmp0$parameters['SD_old_Fem_GP_1', 'Value']
tmp1$parameters['SD_old_Fem_GP_1', 'Value']
tmp2$parameters['SD_old_Fem_GP_1', 'Value']
tmp3$parameters['SD_old_Fem_GP_1', 'Value']
tmp4$parameters['SD_old_Fem_GP_1', 'Value']

tmp0$parameters['SD_old_Fem_GP_1', 'Parm_StDev']
tmp1$parameters['SD_old_Fem_GP_1', 'Parm_StDev']
tmp2$parameters['SD_old_Fem_GP_1', 'Parm_StDev']
tmp3$parameters['SD_old_Fem_GP_1', 'Parm_StDev']
tmp4$parameters['SD_old_Fem_GP_1', 'Parm_StDev']



tmp0$parameters['SR_LN(R0)', 'Value']
tmp1$parameters['SR_LN(R0)', 'Value']
tmp2$parameters['SR_LN(R0)', 'Value']
tmp3$parameters['SR_LN(R0)', 'Value']
tmp4$parameters['SR_LN(R0)', 'Value']


tmp0$parameters['SR_LN(R0)', 'Parm_StDev']
tmp1$parameters['SR_LN(R0)', 'Parm_StDev']
tmp2$parameters['SR_LN(R0)', 'Parm_StDev']
tmp3$parameters['SR_LN(R0)', 'Parm_StDev']
tmp4$parameters['SR_LN(R0)', 'Parm_StDev']



tmp0$parameters['NatM_p_1_Fem_GP_1', 'Value']
tmp1$parameters['NatM_p_1_Fem_GP_1', 'Value']
tmp2$parameters['NatM_p_1_Fem_GP_1', 'Value']
tmp3$parameters['NatM_p_1_Fem_GP_1', 'Value']
tmp4$parameters['NatM_p_1_Fem_GP_1', 'Value']



tmp0$parameters['NatM_p_1_Fem_GP_1', 'Parm_StDev']
tmp1$parameters['NatM_p_1_Fem_GP_1', 'Parm_StDev']
tmp2$parameters['NatM_p_1_Fem_GP_1', 'Parm_StDev']
tmp3$parameters['NatM_p_1_Fem_GP_1', 'Parm_StDev']
tmp4$parameters['NatM_p_1_Fem_GP_1', 'Parm_StDev']




tmp1$likelihoods_by_fleet[8, 4]
tmp2$likelihoods_by_fleet[8, 4]
tmp3$likelihoods_by_fleet[8, 4]
tmp4$likelihoods_by_fleet[8, 4]



mean(tmp0$stdtable[80:123, 'std'])
mean(tmp1$stdtable[80:123, 'std'])
mean(tmp2$stdtable[80:123, 'std'])
mean(tmp3$stdtable[80:123, 'std'])
mean(tmp4$stdtable[80:123, 'std'])


mean(tmp0$stdtable[129:172, 'std'])
mean(tmp1$stdtable[129:172, 'std'])
mean(tmp2$stdtable[129:172, 'std'])
mean(tmp3$stdtable[129:172, 'std'])
mean(tmp4$stdtable[129:172, 'std'])


mean(tmp0$stdtable[80:123, 'std']/tmp0$timeseries$SpawnBio[1:44])
mean(tmp1$stdtable[80:123, 'std']/tmp1$timeseries$SpawnBio[1:44])
mean(tmp2$stdtable[80:123, 'std']/tmp2$timeseries$SpawnBio[1:44])
mean(tmp3$stdtable[80:123, 'std']/tmp3$timeseries$SpawnBio[1:44])
mean(tmp4$stdtable[80:123, 'std']/tmp4$timeseries$SpawnBio[1:44])


mean(tmp0$stdtable[129:172, 'std']/tmp0$timeseries$Recruit_0[1:44])
mean(tmp1$stdtable[129:172, 'std']/tmp1$timeseries$Recruit_0[1:44])
mean(tmp2$stdtable[129:172, 'std']/tmp2$timeseries$Recruit_0[1:44])
mean(tmp3$stdtable[129:172, 'std']/tmp3$timeseries$Recruit_0[1:44])
mean(tmp4$stdtable[129:172, 'std']/tmp4$timeseries$Recruit_0[1:44])


# PLOT FOR SPAWNING DEPLETION

tdat0 = tmp0$timeseries[,c('Yr', 'SpawnBio')]
tdat0$SD = tmp0$stdtable[80:128, 'std']
tdat0$SpawnBio = tdat0$SpawnBio/tdat0$SpawnBio[1]
tdat0$SD = tdat0$SD/1000000
tdat0$upper = tdat0$SpawnBio + tdat0$SD 
tdat0$lower = tdat0$SpawnBio - tdat0$SD 
tdat0$Method = 'SSLen' 

tdat1 = tmp1$timeseries[,c('Yr', 'SpawnBio')]
tdat1$SD = tmp1$stdtable[80:128, 'std']
tdat1$SpawnBio = tdat1$SpawnBio/tdat1$SpawnBio[1]
tdat1$SD = tdat1$SD/1000000
tdat1$upper = tdat1$SpawnBio + tdat1$SD 
tdat1$lower = tdat1$SpawnBio - tdat1$SD 
tdat1$Method = 'SSAge1' 

tdat2 = tmp2$timeseries[,c('Yr', 'SpawnBio')]
tdat2$SD = tmp2$stdtable[80:128, 'std']
tdat2$SpawnBio = tdat2$SpawnBio/tdat2$SpawnBio[1]
tdat2$SD = tdat2$SD/1000000
tdat2$upper = tdat2$SpawnBio + tdat2$SD 
tdat2$lower = tdat2$SpawnBio - tdat2$SD 
tdat2$Method = 'SSAge2' 

tdat3 = tmp3$timeseries[,c('Yr', 'SpawnBio')]
tdat3$SD = tmp3$stdtable[80:128, 'std']
tdat3$SpawnBio = tdat3$SpawnBio/tdat3$SpawnBio[1]
tdat3$SD = tdat3$SD/1000000
tdat3$upper = tdat3$SpawnBio + tdat3$SD 
tdat3$lower = tdat3$SpawnBio - tdat3$SD 
tdat3$Method = 'SSAge3' 

tdat4 = tmp4$timeseries[,c('Yr', 'SpawnBio')]
tdat4$SD = tmp4$stdtable[80:128, 'std']
tdat4$SpawnBio = tdat4$SpawnBio/tdat4$SpawnBio[1]
tdat4$SD = tdat4$SD/1000000
tdat4$upper = tdat4$SpawnBio + tdat4$SD 
tdat4$lower = tdat4$SpawnBio - tdat4$SD 
tdat4$Method = 'SSAge4' 


tdat = rbind(tdat0[tdat0$Yr %in% 1976:2018, ], 
			tdat1[tdat1$Yr %in% 1976:2018, ], 
			 tdat2[tdat2$Yr %in% 1976:2018, ], 
			 tdat3[tdat3$Yr %in% 1976:2018, ], 
			 tdat4[tdat4$Yr %in% 1976:2018, ])

# tdat2 = rbind(tdat0[tdat0$Yr %in% c(1975, 2019:2023), ], 
# 				tdat1[tdat1$Yr %in% c(1975, 2019:2023), ], 
# 			 tdat2[tdat2$Yr %in% c(1975, 2019:2023), ], 
# 			 tdat3[tdat3$Yr %in% c(1975, 2019:2023), ], 
# 			 tdat4[tdat4$Yr %in% c(1975, 2019:2023), ])

tdat2 = rbind(tdat0[tdat0$Yr %in% c(1975), ], 
				tdat1[tdat1$Yr %in% c(1975), ], 
			 tdat2[tdat2$Yr %in% c(1975), ], 
			 tdat3[tdat3$Yr %in% c(1975), ], 
			 tdat4[tdat4$Yr %in% c(1975), ])

xa1 = ggplot(tdat, aes(x = Yr, y = SpawnBio, color = Method)) +
    geom_line() +
    geom_line(mapping = aes(y = upper), lty = "dashed") +
    geom_line(mapping = aes(y = lower), lty = "dashed") +
    geom_point(data = tdat2, aes(x = Yr, y = SpawnBio, color = Method)) +
    geom_errorbar(data = tdat2, aes(ymin = lower, ymax = upper), width=.35) +
    scale_color_manual(values=c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF",'black')) +
	xlab('') +
	ylab(bquote('Spawning depletion')) +
	xlim(1974, 2019) +
    theme_bw() +
    theme(legend.position = c(0.7, 0.72), legend.title = element_blank(), legend.text = element_text(size = 7))

xa2 = ggplot(tdat, aes(x = Yr)) +
    geom_ribbon(aes(ymin = lower, ymax=upper, fill = Method), alpha=0.2) +
    geom_line(aes(y = SpawnBio, color = Method)) +
    #geom_line(mapping = aes(y = upper), lty = "dashed") +
    #geom_line(mapping = aes(y = lower), lty = "dashed") +
    geom_point(data = tdat2, aes(x = Yr, y = SpawnBio, color = Method)) +
    geom_errorbar(data = tdat2, aes(ymin = lower, ymax = upper, color = Method), width=.35) +
    scale_fill_manual(values=c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF", "grey70"), name = 'fill') +
    scale_color_manual(values=c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF", 'black')) +
	xlab('') +
	ylab(bquote('Spawning depletion')) +
	xlim(1974, 2019) +
    theme_bw() +
    theme(legend.position = c(0.7, 0.72), legend.title = element_blank(), legend.text = element_text(size = 7), 
    	legend.key.size = unit(0.4, "cm")) +
    guides(fill = FALSE)


# PLOT FOR RECRUITMENT

tdat0r = tmp0$timeseries[1:44,c('Yr', 'Recruit_0')]
tdat0r$SD = tmp0$stdtable[129:172, 'std']
tdat0r$Recruit_0 = tdat0r$Recruit_0/1000000
tdat0r$SD = tdat0r$SD/1000000
tdat0r$upper = tdat0r$Recruit_0 + tdat0r$SD 
tdat0r$lower = tdat0r$Recruit_0 - tdat0r$SD 
tdat0r$Method = 'SSLen' 

tdat1r = tmp1$timeseries[1:44,c('Yr', 'Recruit_0')]
tdat1r$SD = tmp1$stdtable[129:172, 'std']
tdat1r$Recruit_0 = tdat1r$Recruit_0/1000000
tdat1r$SD = tdat1r$SD/1000000
tdat1r$upper = tdat1r$Recruit_0 + tdat1r$SD 
tdat1r$lower = tdat1r$Recruit_0 - tdat1r$SD 
tdat1r$Method = 'SSAge1' 

tdat2r = tmp2$timeseries[1:44,c('Yr', 'Recruit_0')]
tdat2r$SD = tmp2$stdtable[129:172, 'std']
tdat2r$Recruit_0 = tdat2r$Recruit_0/1000000
tdat2r$SD = tdat2r$SD/1000000
tdat2r$upper = tdat2r$Recruit_0 + tdat2r$SD 
tdat2r$lower = tdat2r$Recruit_0 - tdat2r$SD 
tdat2r$Method = 'SSAge2' 

tdat3r = tmp3$timeseries[1:44,c('Yr', 'Recruit_0')]
tdat3r$SD = tmp3$stdtable[129:172, 'std']
tdat3r$Recruit_0 = tdat3r$Recruit_0/1000000
tdat3r$SD = tdat3r$SD/1000000
tdat3r$upper = tdat3r$Recruit_0 + tdat3r$SD 
tdat3r$lower = tdat3r$Recruit_0 - tdat3r$SD 
tdat3r$Method = 'SSAge3' 

tdat4r = tmp4$timeseries[1:44,c('Yr', 'Recruit_0')]
tdat4r$SD = tmp4$stdtable[129:172, 'std']
tdat4r$Recruit_0 = tdat4r$Recruit_0/1000000
tdat4r$SD = tdat4r$SD/1000000
tdat4r$upper = tdat4r$Recruit_0 + tdat4r$SD 
tdat4r$lower = tdat4r$Recruit_0 - tdat4r$SD 
tdat4r$Method = 'SSAge4' 


tdat3 = rbind(tdat0r[tdat0r$Yr %in% 1976:2018, ], 
			tdat1r[tdat1r$Yr %in% 1976:2018, ], 
			 tdat2r[tdat2r$Yr %in% 1976:2018, ], 
			 tdat3r[tdat3r$Yr %in% 1976:2018, ], 
			 tdat4r[tdat4r$Yr %in% 1976:2018, ])


tdat4 = rbind(tdat0r[tdat0r$Yr %in% c(1975), ], 
			 tdat1r[tdat1r$Yr %in% c(1975), ], 
			 tdat2r[tdat2r$Yr %in% c(1975), ], 
			 tdat3r[tdat3r$Yr %in% c(1975), ], 
			 tdat4r[tdat4r$Yr %in% c(1975), ])


xa3 = ggplot(tdat3, aes(x = Yr, y = Recruit_0, color = Method)) +
    geom_line() +
    geom_line(mapping = aes(y = upper), lty = "dashed") +
    geom_line(mapping = aes(y = lower), lty = "dashed") +
    geom_point(data = tdat4, aes(x = Yr, y = Recruit_0, color = Method)) +
    geom_errorbar(data = tdat4, aes(ymin = lower, ymax = upper), width=.35) +
    scale_color_manual(values=c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF", 'black')) +
	xlab('Year') +
	ylab(bquote('Age-0 recruits ('*10^9~ind*')')) +
	xlim(1974, 2019) +
    theme_bw() +
    theme(legend.position = 'none')

xa4 = ggplot(tdat3, aes(x = Yr)) +
    geom_ribbon(aes(ymin = lower, ymax=upper, fill = Method), alpha=0.2) +
    geom_line(aes(y = Recruit_0, color = Method)) +
    #geom_line(mapping = aes(y = upper), lty = "dashed") +
    #geom_line(mapping = aes(y = lower), lty = "dashed") +
    geom_point(data = tdat4, aes(x = Yr, y = Recruit_0, color = Method)) +
    geom_errorbar(data = tdat4, aes(ymin = lower, ymax = upper, color = Method), width=.35) +
    scale_fill_manual(values=c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF", "grey70"), name = 'fill') +
    scale_color_manual(values=c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF", 'black')) +
	xlab('Year') +
	ylab(bquote('Age-0 recruits ('*10^9~ind*')')) +
	xlim(1974, 2019) +
    theme_bw() +
    theme(legend.position = 'none')




bitmap('FinalFigures/compare_SSBRec_allMethods2.tiff', height = 120, width = 110, units = 'mm', res = 500)

	grid.arrange(xa1, xa3, nrow = 2)

dev.off()

png('FinalFigures/compare_SSBRec_allMethods3.png', height = 120, width = 110, units = 'mm', res = 500)

	grid.arrange(xa2, xa4, nrow = 2)

dev.off()


# ----------------------------------
# FIND FORECAST CATCH

forecatch0 = SS_ForeCatch(replist = tmp0)
forecatch1 = SS_ForeCatch(replist = tmp1)
forecatch2 = SS_ForeCatch(replist = tmp2)
forecatch3 = SS_ForeCatch(replist = tmp3)
forecatch4 = SS_ForeCatch(replist = tmp4)

mean(forecatch0[,4])
mean(forecatch1[,4])
mean(forecatch2[,4])
mean(forecatch3[,4])
mean(forecatch4[,4])