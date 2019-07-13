require(ggplot2)
require(reshape)
require(gridExtra)

# These two values should be the same as in paramtersSimulation:
scenarioName = 'NoS_NoT'
agePlus = 8

# New parameters
selFolder2 = 'simPerInd2'

# For second figure
allFiles2 = list.files(path = selFolder2)
allFiles22 = allFiles2[grep(pattern = scenarioName, x = allFiles2)]

#  ------------------------------------------------------------------------
# Plot all methods: Figure 3
allMethodsPropYear3 = NULL
for(k in seq_along(allFiles22)){
	tmp = read.csv(file.path(selFolder2, allFiles22[k]))
	allMethodsPropYear3 = rbind(allMethodsPropYear3, tmp)
}

allMethodsPropYear4 = aggregate(allMethodsPropYear3$MRE*100, list(YEAR = allMethodsPropYear3$YEAR, 
								AGE = allMethodsPropYear3$AGE, METHOD = allMethodsPropYear3$METHOD), mean)


	# bitmap(paste0('FinalFigures/compareProportions3_', scenarioName, '.tiff'), height = 190, width = 190, units = 'mm', res = 900)
	# print(ggplot(allMethodsPropYear4, aes(AGE, x)) +
			# geom_line(aes(color = factor(METHOD))) +
			# xlab('Age') +
			# ylab('MRE') +
			# facet_wrap( ~ factor(YEAR), ncol = 5) +
			# theme_bw() +
			# theme(legend.position = c(0.75, 0.085), legend.title = element_blank()) +
			# scale_x_discrete(limits = 1:agePlus))
	# dev.off()


#  ------------------------------------------------------------------------
# Plot all methods: Figure 4: MRE per age

allMethodsPropYear5 = aggregate(allMethodsPropYear3$MRE*100, list(AGE = allMethodsPropYear3$AGE, 
								METHOD = allMethodsPropYear3$METHOD), mean)
allMethodsPropYear55 = aggregate(allMethodsPropYear3$MRE*100, list(AGE = allMethodsPropYear3$AGE, 
								METHOD = allMethodsPropYear3$METHOD), sd)
allMethodsPropYear5$sdx = allMethodsPropYear55$x
# allMethodsPropYear5$menlelb = allMethodsPropYear5$x - allMethodsPropYear5$sdx
# allMethodsPropYear5$menleub = allMethodsPropYear5$x + allMethodsPropYear5$sdx



	ca1 = ggplot(allMethodsPropYear5, aes(AGE, x, color = factor(METHOD))) +
			geom_errorbar(aes(ymin = x - sdx, ymax = x + sdx), width=.35) +
			geom_line(aes(color = factor(METHOD))) +
			geom_point(aes(color = factor(METHOD)), size = 1) +
			xlab('Age') +
			ylab(bquote(MRE[age])) +
			facet_wrap( ~ factor(METHOD), ncol = 5, scales = 'free_y') +
			theme_bw() +
			theme(legend.position="none") +
			scale_x_discrete(limits = 1:agePlus)

#  ------------------------------------------------------------------------
# Plot all methods: Figure 5: MRE per year

allMethodsPropYear6 = aggregate(allMethodsPropYear3$MRE*100, list(YEAR = allMethodsPropYear3$YEAR, 
								METHOD = allMethodsPropYear3$METHOD), mean)
allMethodsPropYear66 = aggregate(allMethodsPropYear3$MRE*100, list(YEAR = allMethodsPropYear3$YEAR, 
								METHOD = allMethodsPropYear3$METHOD), sd)
allMethodsPropYear6$sdx = allMethodsPropYear66$x

	ca2 = ggplot(allMethodsPropYear6, aes(YEAR, x, color = factor(METHOD))) +
			geom_errorbar(aes(ymin = x - sdx, ymax = x + sdx), width=.35) +
			geom_line(aes(color = factor(METHOD))) +
			geom_point(aes(color = factor(METHOD)), size = 1) +
			xlab('Year') +
			ylab(bquote(MRE[year])) +
			facet_wrap( ~ factor(METHOD), ncol = 5, scales = 'free_y') +
			theme_bw() +
			theme(legend.position="none")


# Plot both:
bitmap(paste0('FinalFigures/compareProportions45_', scenarioName, '.tiff'), height = 120, width = 190, units = 'mm', res = 900)

	grid.arrange(ca1, ca2, nrow = 2)

dev.off()

#  ------------------------------------------------------------------------
# Plot all methods: Figure 6: MRE per year per age. One figure per method

allMethodsPropYear44 = aggregate(allMethodsPropYear3$MRE*100, list(YEAR = allMethodsPropYear3$YEAR, 
								AGE = allMethodsPropYear3$AGE, METHOD = allMethodsPropYear3$METHOD), sd)
allMethodsPropYear4$sdx = allMethodsPropYear44$x

tmp1 = allMethodsPropYear4[allMethodsPropYear4$METHOD == 'Method1', ]
	bitmap(paste0('FinalFigures/compareProportions6_1_', scenarioName, '.tiff'), height = 190, width = 190, units = 'mm', res = 900)
	print(ggplot(tmp1, aes(AGE, x)) +
			geom_errorbar(aes(ymin = x - sdx, ymax = x + sdx), width=.35, color = "#F8766D") +
			geom_line(color = "#F8766D") +
			geom_point(color = "#F8766D", size = 1) +
			xlab('Age') +
			ylab(bquote(MRE['a,y'])) +
			facet_wrap( ~ factor(YEAR), ncol = 5) +
			theme_bw() +
			theme(legend.position="none")+
			scale_x_discrete(limits = 1:agePlus))
	dev.off()


tmp2 = allMethodsPropYear4[allMethodsPropYear4$METHOD == 'Method2', ]
	bitmap(paste0('FinalFigures/compareProportions6_2_', scenarioName, '.tiff'), height = 190, width = 190, units = 'mm', res = 900)
	print(ggplot(tmp2, aes(AGE, x)) +
			geom_errorbar(aes(ymin = x - sdx, ymax = x + sdx), width=.35, color = "#7CAE00") +
			geom_line(color = "#7CAE00") +
			geom_point(color = "#7CAE00", size = 1) +
			xlab('Age') +
			ylab(bquote(MRE['a,y'])) +
			facet_wrap( ~ factor(YEAR), ncol = 5) +
			theme_bw() +
			theme(legend.position="none")+
			scale_x_discrete(limits = 1:agePlus))
	dev.off()


tmp3 = allMethodsPropYear4[allMethodsPropYear4$METHOD == 'Method3', ]
	bitmap(paste0('FinalFigures/compareProportions6_3_', scenarioName, '.tiff'), height = 190, width = 190, units = 'mm', res = 900)
	print(ggplot(tmp3, aes(AGE, x)) +
			geom_errorbar(aes(ymin = x - sdx, ymax = x + sdx), width=.35, color = "#00BFC4") +
			geom_line(color = "#00BFC4") +
			geom_point(color = "#00BFC4", size = 1) +
			xlab('Age') +
			ylab(bquote(MRE['a,y'])) +
			facet_wrap( ~ factor(YEAR), ncol = 5) +
			theme_bw() +
			theme(legend.position="none")+
			scale_x_discrete(limits = 1:agePlus))
	dev.off()


tmp4 = allMethodsPropYear4[allMethodsPropYear4$METHOD == 'Method4', ]
	bitmap(paste0('FinalFigures/compareProportions6_4_', scenarioName, '.tiff'), height = 190, width = 190, units = 'mm', res = 900)
	print(ggplot(tmp4, aes(AGE, x)) +
			geom_errorbar(aes(ymin = x - sdx, ymax = x + sdx), width=.35, color = "#C77CFF") +
			geom_line(color = "#C77CFF") +
			geom_point(color = "#C77CFF", size = 1) +
			xlab('Age') +
			ylab(bquote(MRE['a,y'])) +
			facet_wrap( ~ factor(YEAR), ncol = 5) +
			theme_bw() +
			theme(legend.position="none")+
			scale_x_discrete(limits = 1:agePlus))
	dev.off()


#  ------------------------------------------------------------------------
# TABLES MRE
allMethodsPropYear3 = NULL
for(k in seq_along(allFiles22)){
	tmp = read.csv(file.path(selFolder2, allFiles22[k]))
	allMethodsPropYear3 = rbind(allMethodsPropYear3, tmp)
}

# Get table 1: MRE per age
allMethodsPropYear41 = aggregate(allMethodsPropYear3$MRE*100, list(AGE = allMethodsPropYear3$AGE, 
								METHOD = allMethodsPropYear3$METHOD), mean)
allMethodsPropYear41sd = aggregate(allMethodsPropYear3$MRE*100, list(AGE = allMethodsPropYear3$AGE, 
								METHOD = allMethodsPropYear3$METHOD), sd)
tab1 = cast(allMethodsPropYear41, METHOD ~ AGE)
tab1sd = cast(allMethodsPropYear41sd, METHOD ~ AGE)
write.csv(tab1, paste0('FinalTables/MRE_age_', scenarioName, '_mean.csv'), row.names = FALSE)
write.csv(tab1sd, paste0('FinalTables/MRE_age_', scenarioName, '_sd.csv'), row.names = FALSE)

# Get table 2: MRE per year
allMethodsPropYear41 = aggregate(allMethodsPropYear3$MRE*100, list(YEAR = allMethodsPropYear3$YEAR, 
								METHOD = allMethodsPropYear3$METHOD), mean)
allMethodsPropYear41sd = aggregate(allMethodsPropYear3$MRE*100, list(YEAR = allMethodsPropYear3$YEAR, 
								METHOD = allMethodsPropYear3$METHOD), sd)
tab1 = cast(allMethodsPropYear41, METHOD ~ YEAR)
tab1sd = cast(allMethodsPropYear41sd, METHOD ~ YEAR)
write.csv(tab1, paste0('FinalTables/MRE_year_', scenarioName, '_mean.csv'), row.names = FALSE)
write.csv(tab1sd, paste0('FinalTables/MRE_year_', scenarioName, '_sd.csv'), row.names = FALSE)


# Get table 3: MRE total
allMethodsPropYear43 = aggregate(allMethodsPropYear3$MRE*100, list(METHOD = allMethodsPropYear3$METHOD), mean)
allMethodsPropYear43sd = aggregate(allMethodsPropYear3$MRE*100, list(METHOD = allMethodsPropYear3$METHOD), sd)

write.csv(allMethodsPropYear43, paste0('FinalTables/MRE_all_', scenarioName, '_mean.csv'), row.names = FALSE)
write.csv(allMethodsPropYear43sd, paste0('FinalTables/MRE_all_', scenarioName, '_sd.csv'), row.names = FALSE)


#  ------------------------------------------------------------------------
# TABLES MSE

# Table 1: MSE per age
allMethodsPropYear41 = aggregate(allMethodsPropYear3$MSE, list(AGE = allMethodsPropYear3$AGE, 
								METHOD = allMethodsPropYear3$METHOD), mean)
allMethodsPropYear41sd = aggregate(allMethodsPropYear3$MSE, list(AGE = allMethodsPropYear3$AGE, 
								METHOD = allMethodsPropYear3$METHOD), sd)
tab1 = cast(allMethodsPropYear41, METHOD ~ AGE)
tab1sd = cast(allMethodsPropYear41sd, METHOD ~ AGE)
write.csv(tab1, paste0('FinalTables/MSE_age_', scenarioName, '_mean.csv'), row.names = FALSE)
write.csv(tab1sd, paste0('FinalTables/MSE_age_', scenarioName, '_sd.csv'), row.names = FALSE)

# Get table 2: MSE per year
allMethodsPropYear41 = aggregate(allMethodsPropYear3$MSE, list(YEAR = allMethodsPropYear3$YEAR, 
								METHOD = allMethodsPropYear3$METHOD), mean)
allMethodsPropYear41sd = aggregate(allMethodsPropYear3$MSE, list(YEAR = allMethodsPropYear3$YEAR, 
								METHOD = allMethodsPropYear3$METHOD), sd)
tab1 = cast(allMethodsPropYear41, METHOD ~ YEAR)
tab1sd = cast(allMethodsPropYear41sd, METHOD ~ YEAR)
write.csv(tab1, paste0('FinalTables/MSE_year_', scenarioName, '_mean.csv'), row.names = FALSE)
write.csv(tab1sd, paste0('FinalTables/MSE_year_', scenarioName, '_sd.csv'), row.names = FALSE)

# Get table 3: MSE total
allMethodsPropYear43 = aggregate(allMethodsPropYear3$MSE, list(METHOD = allMethodsPropYear3$METHOD), mean)
allMethodsPropYear43sd = aggregate(allMethodsPropYear3$MSE, list(METHOD = allMethodsPropYear3$METHOD), sd)

write.csv(allMethodsPropYear43, paste0('FinalTables/MSE_all_', scenarioName, '_mean.csv'), row.names = FALSE)
write.csv(allMethodsPropYear43sd, paste0('FinalTables/MSE_all_', scenarioName, '_sd.csv'), row.names = FALSE)
