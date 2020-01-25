require(ggplot2)
require(reshape)
require(gridExtra)
require(scales)
require(grid)

#setwd('C:/Users/moroncog/Documents/GitHub/STageCompsEstimation')
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


#  ------------------------------------------------------------------------
# MRE plots:

allMethodsPropYear5 = aggregate(allMethodsPropYear3$MRE*100, list(AGE = allMethodsPropYear3$AGE, 
								METHOD = allMethodsPropYear3$METHOD), mean)
allMethodsPropYear55 = aggregate(allMethodsPropYear3$MRE*100, list(AGE = allMethodsPropYear3$AGE, 
								METHOD = allMethodsPropYear3$METHOD), sd)
allMethodsPropYear5$sdx = allMethodsPropYear55$x
# allMethodsPropYear5$menlelb = allMethodsPropYear5$x - allMethodsPropYear5$sdx
# allMethodsPropYear5$menleub = allMethodsPropYear5$x + allMethodsPropYear5$sdx

# MRE for age:
mreage = 	ggplot(data=allMethodsPropYear5, aes(x=as.factor(AGE), y=x, fill=as.factor(METHOD))) +
				geom_bar(stat="identity", color="black", position=position_dodge(), width = 0.7)+
			#	geom_errorbar(aes(ymin = x-sd, ymax = x+sd), width=.2,
			#                 position=position_dodge(.9)) +
				scale_fill_manual(values = c("#F8766D", "#7CAE00", "#00BFC4", '#C77CFF'), name = '') +
			  	ylab(bquote(MRE['age']~ '(%)')) +
				xlab('') +
				scale_y_continuous(limits=c(max(min(allMethodsPropYear5$x), -10), min(max(allMethodsPropYear5$x), 10)),oob = rescale_none) +
				theme_bw() +
				theme(legend.position='none')

if(max(allMethodsPropYear5$x) > 10 | min(allMethodsPropYear5$x) < -10) {

	subTemp = allMethodsPropYear5[allMethodsPropYear5$x > 10 | allMethodsPropYear5$x < -10, ]
	tmps2 = as.numeric(gsub(pattern = 'Method', replacement = '', x = as.character(subTemp$METHOD)))
	tmps3 = tmps2 - 1
	tmps4 = 0.15*tmps3
	Labs = round(subTemp$x, 0)
	xPos = subTemp$AGE + tmps4 # minAge
	yPos = numeric(nrow(subTemp))
	yPos[subTemp$x > 10] = 10
	yPos[subTemp$x < -10] = -10

	mreage = mreage + annotate(geom="text", x= xPos, y = yPos, label = Labs, size = 2.5, angle = 90)

}


#  ------------------------------------------------------------------------
# Plot all methods: Figure 5: MRE per year

allMethodsPropYear6 = aggregate(allMethodsPropYear3$MRE*100, list(YEAR = allMethodsPropYear3$YEAR, 
								METHOD = allMethodsPropYear3$METHOD), mean)
allMethodsPropYear66 = aggregate(allMethodsPropYear3$MRE*100, list(YEAR = allMethodsPropYear3$YEAR, 
								METHOD = allMethodsPropYear3$METHOD), sd)
allMethodsPropYear6$sdx = allMethodsPropYear66$x

allMethodsPropYear6$CAT = NA
allMethodsPropYear6$CAT[allMethodsPropYear6$YEAR < 2002] = '1994 - 2001'
allMethodsPropYear6$CAT[allMethodsPropYear6$YEAR > 2001 & allMethodsPropYear6$YEAR < 2010] = '2002 - 2009'
allMethodsPropYear6$CAT[allMethodsPropYear6$YEAR > 2009] = '2010 - 2016'
allMethodsPropYear6x = aggregate(allMethodsPropYear6$x, list(METHOD = allMethodsPropYear6$METHOD,
									PERIOD = allMethodsPropYear6$CAT), FUN = mean)


mreyear = 	ggplot(data = allMethodsPropYear6x, aes(y = x, x = factor(PERIOD), fill = factor(METHOD))) +
		geom_bar(stat="identity", color="black", position='dodge', width = 0.7)+
		scale_fill_manual(values = c("#F8766D", "#7CAE00", "#00BFC4", '#C77CFF'), name = '') +
	  	ylab(bquote(MRE['year']~ '(%)')) +
		xlab('') +
		scale_y_continuous(limits=c(max(min(allMethodsPropYear6$x), -10), min(max(allMethodsPropYear6$x), 10)),oob = rescale_none) +
		theme_bw() +
		theme(legend.position='none') 

if(max(allMethodsPropYear6$x) > 10 | min(allMethodsPropYear6$x) < -10) {

	subTemp = allMethodsPropYear6x[allMethodsPropYear6x$x > 10 | allMethodsPropYear6x$x < -10, ]
	
	tmps2 = as.numeric(gsub(pattern = 'Method', replacement = '', x = as.character(subTemp$METHOD)))
	tmps3 = tmps2 - 1
	tmps4 = 0.15*tmps3
	Labs = round(subTemp$x, 0)

	xPos = subTemp$YEAR - 1994 + tmps4 # minYear
	
	yPos = numeric(nrow(subTemp))
	yPos[subTemp$x > 10] = 9.6
	yPos[subTemp$x < -10] = -9.6

	mreyear = mreyear + annotate(geom="text", x= xPos, y = yPos, label = Labs, size = 2)

}


#  ------------------------------------------------------------------------
# TABLES MRE
allMethodsPropYear3 = NULL
for(k in seq_along(allFiles22)){
	tmp = read.csv(file.path(selFolder2, allFiles22[k]))
	allMethodsPropYear3 = rbind(allMethodsPropYear3, tmp)
}

# Get table 3: MRE total
allMethodsPropYear43 = aggregate(allMethodsPropYear3$MRE*100, list(METHOD = allMethodsPropYear3$METHOD), mean)
allMethodsPropYear43sd = aggregate(allMethodsPropYear3$MRE*100, list(METHOD = allMethodsPropYear3$METHOD), sd)

write.csv(allMethodsPropYear43, paste0('FinalTables/MRE_all_', scenarioName, '_mean.csv'), row.names = FALSE)
write.csv(allMethodsPropYear43sd, paste0('FinalTables/MRE_all_', scenarioName, '_sd.csv'), row.names = FALSE)


#  ------------------------------------------------------------------------
# TABLES and FIGURES MSE

# Table 1: MSE per age
allMethodsPropYear41 = aggregate(allMethodsPropYear3$MSE, list(AGE = allMethodsPropYear3$AGE, 
								METHOD = allMethodsPropYear3$METHOD), mean)
allMethodsPropYear41sd = aggregate(allMethodsPropYear3$MSE, list(AGE = allMethodsPropYear3$AGE, 
								METHOD = allMethodsPropYear3$METHOD), sd)
tab1 = cast(allMethodsPropYear41, METHOD ~ AGE)
tab1sd = cast(allMethodsPropYear41sd, METHOD ~ AGE)
write.csv(tab1, paste0('FinalTables/MSE_age_', scenarioName, '_mean.csv'), row.names = FALSE)
write.csv(tab1sd, paste0('FinalTables/MSE_age_', scenarioName, '_sd.csv'), row.names = FALSE)

# MSE per age
MSEdf = allMethodsPropYear41
MSEdf$x = MSEdf$x*10^5

mseage = 	ggplot(data=MSEdf, aes(x=as.factor(AGE), y=x, fill=as.factor(METHOD))) +
				geom_bar(stat="identity", color="black", position=position_dodge(), width = 0.7)+
				scale_fill_manual(values = c("#F8766D", "#7CAE00", "#00BFC4", '#C77CFF'), name = '') +
			  	ylab(bquote(MSE['age'] ~ ' (10'^'-5'~')')) +
				xlab('') +
				scale_y_continuous(limits=c(0, min(max(MSEdf$x), 10)),oob = rescale_none) +
				theme_bw() +
				ggtitle("No S / No T") +
				theme(legend.position = 'none', plot.title = element_text(hjust = 0.5), legend.text=element_text(size=7.5))

if(scenarioName == 'NoS_NoT') mseage = mseage + theme(legend.position=c(0.8, 0.73))


if(max(MSEdf$x) > 10) {

	subTemp = MSEdf[MSEdf$x > 10, ]
	tmps2 = as.numeric(gsub(pattern = 'Method', replacement = '', x = as.character(subTemp$METHOD)))
	tmps3 = tmps2 - 1
	tmps4 = 0.15*tmps3
	Labs = round(subTemp$x, 0)
	xPos = subTemp$AGE + tmps4 # minAge
	yPos = numeric(nrow(subTemp))
	yPos[subTemp$x > 10] = 10

	mseage = mseage + annotate(geom="text", x= xPos, y = yPos, label = Labs, size = 2.5, angle = 90)

}




# Get table 2: MSE per year
allMethodsPropYear41 = aggregate(allMethodsPropYear3$MSE, list(YEAR = allMethodsPropYear3$YEAR, 
								METHOD = allMethodsPropYear3$METHOD), mean)
allMethodsPropYear41sd = aggregate(allMethodsPropYear3$MSE, list(YEAR = allMethodsPropYear3$YEAR, 
								METHOD = allMethodsPropYear3$METHOD), sd)
tab1 = cast(allMethodsPropYear41, METHOD ~ YEAR)
tab1sd = cast(allMethodsPropYear41sd, METHOD ~ YEAR)
write.csv(tab1, paste0('FinalTables/MSE_year_', scenarioName, '_mean.csv'), row.names = FALSE)
write.csv(tab1sd, paste0('FinalTables/MSE_year_', scenarioName, '_sd.csv'), row.names = FALSE)


# MSE per year
MSEdf = allMethodsPropYear41
MSEdf$x = MSEdf$x*10^5

MSEdf$CAT = NA
MSEdf$CAT[MSEdf$YEAR < 2002] = '1994 - 2001'
MSEdf$CAT[MSEdf$YEAR > 2001 & MSEdf$YEAR < 2010] = '2002 - 2009'
MSEdf$CAT[MSEdf$YEAR > 2009] = '2010 - 2016'
MSEdfx = aggregate(MSEdf$x, list(METHOD = MSEdf$METHOD,
									PERIOD = MSEdf$CAT), FUN = mean)



mseyear = 	ggplot(data = MSEdfx, aes(y = x, x = factor(PERIOD), fill = factor(METHOD))) +
		geom_bar(stat="identity", color="black", position='dodge', width = 0.7)+
		scale_fill_manual(values = c("#F8766D", "#7CAE00", "#00BFC4", '#C77CFF'), name = '') +
	  	ylab(bquote(MSE['year'] ~ ' (10'^'-5'~')')) +
		xlab('') +
		scale_y_continuous(limits=c(0, min(max(MSEdf$x), 10)),oob = rescale_none) +
		theme_bw() +
		ggtitle("No S / No T") +
		theme(legend.position='none', plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) 

if(max(MSEdf$x) > 10) {

	subTemp = MSEdf[MSEdf$x > 10, ]
	tmps2 = as.numeric(gsub(pattern = 'Method', replacement = '', x = as.character(subTemp$METHOD)))
	tmps3 = tmps2 - 1
	tmps4 = 0.15*tmps3
	Labs = round(subTemp$x, 0)

	xPos = subTemp$YEAR - 1994 + tmps4 # minYear
	
	yPos = numeric(nrow(subTemp))
	yPos[subTemp$x > 10] = 9.6

	mseyear = mseyear + annotate(geom="text", x= xPos, y = yPos, label = Labs, size = 2.5)

}



# Get table 3: MSE total
allMethodsPropYear43 = aggregate(allMethodsPropYear3$MSE, list(METHOD = allMethodsPropYear3$METHOD), mean)
allMethodsPropYear43sd = aggregate(allMethodsPropYear3$MSE, list(METHOD = allMethodsPropYear3$METHOD), sd)

write.csv(allMethodsPropYear43, paste0('FinalTables/MSE_all_', scenarioName, '_mean.csv'), row.names = FALSE)
write.csv(allMethodsPropYear43sd, paste0('FinalTables/MSE_all_', scenarioName, '_sd.csv'), row.names = FALSE)

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

#setwd('C:/Users/moroncog/Documents/GitHub/STageCompsEstimation')
# These two values should be the same as in paramtersSimulation:
scenarioName = 'HighS_HighT'
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


#  ------------------------------------------------------------------------
# MRE plots:

allMethodsPropYear5 = aggregate(allMethodsPropYear3$MRE*100, list(AGE = allMethodsPropYear3$AGE, 
								METHOD = allMethodsPropYear3$METHOD), mean)
allMethodsPropYear55 = aggregate(allMethodsPropYear3$MRE*100, list(AGE = allMethodsPropYear3$AGE, 
								METHOD = allMethodsPropYear3$METHOD), sd)
allMethodsPropYear5$sdx = allMethodsPropYear55$x
# allMethodsPropYear5$menlelb = allMethodsPropYear5$x - allMethodsPropYear5$sdx
# allMethodsPropYear5$menleub = allMethodsPropYear5$x + allMethodsPropYear5$sdx

# MRE for age:
mreage2 = 	ggplot(data=allMethodsPropYear5, aes(x=as.factor(AGE), y=x, fill=as.factor(METHOD))) +
				geom_bar(stat="identity", color="black", position=position_dodge(), width = 0.7)+
			#	geom_errorbar(aes(ymin = x-sd, ymax = x+sd), width=.2,
			#                 position=position_dodge(.9)) +
				scale_fill_manual(values = c("#F8766D", "#7CAE00", "#00BFC4", '#C77CFF'), name = '') +
			  	ylab('') +
				xlab('') +
				scale_y_continuous(limits=c(max(min(allMethodsPropYear5$x), -10), min(max(allMethodsPropYear5$x), 10)),oob = rescale_none) +
				theme_bw() +
				theme(legend.position='none')

if(max(allMethodsPropYear5$x) > 10 | min(allMethodsPropYear5$x) < -10) {

	subTemp = allMethodsPropYear5[allMethodsPropYear5$x > 10 | allMethodsPropYear5$x < -10, ]
	tmps2 = as.numeric(gsub(pattern = 'Method', replacement = '', x = as.character(subTemp$METHOD)))
	tmps3 = tmps2 - 1
	tmps4 = 0.15*tmps3
	Labs = round(subTemp$x, 0)
	xPos = subTemp$AGE + tmps4 # minAge
	yPos = numeric(nrow(subTemp))
	yPos[subTemp$x > 10] = 10
	yPos[subTemp$x < -10] = -10

	mreage2 = mreage2 + annotate(geom="text", x= xPos, y = yPos, label = Labs, size = 2.5, angle = 90)

}


#  ------------------------------------------------------------------------
# Plot all methods: Figure 5: MRE per year

allMethodsPropYear6 = aggregate(allMethodsPropYear3$MRE*100, list(YEAR = allMethodsPropYear3$YEAR, 
								METHOD = allMethodsPropYear3$METHOD), mean)
allMethodsPropYear66 = aggregate(allMethodsPropYear3$MRE*100, list(YEAR = allMethodsPropYear3$YEAR, 
								METHOD = allMethodsPropYear3$METHOD), sd)
allMethodsPropYear6$sdx = allMethodsPropYear66$x

allMethodsPropYear6$CAT = NA
allMethodsPropYear6$CAT[allMethodsPropYear6$YEAR < 2002] = '1994 - 2001'
allMethodsPropYear6$CAT[allMethodsPropYear6$YEAR > 2001 & allMethodsPropYear6$YEAR < 2010] = '2002 - 2009'
allMethodsPropYear6$CAT[allMethodsPropYear6$YEAR > 2009] = '2010 - 2016'
allMethodsPropYear6x = aggregate(allMethodsPropYear6$x, list(METHOD = allMethodsPropYear6$METHOD,
									PERIOD = allMethodsPropYear6$CAT), FUN = mean)


mreyear2 = 	ggplot(data = allMethodsPropYear6x, aes(y = x, x = factor(PERIOD), fill = factor(METHOD))) +
		geom_bar(stat="identity", color="black", position='dodge', width = 0.7)+
		scale_fill_manual(values = c("#F8766D", "#7CAE00", "#00BFC4", '#C77CFF'), name = '') +
	  	ylab('') +
		xlab('') +
		scale_y_continuous(limits=c(max(min(allMethodsPropYear6$x), -10), min(max(allMethodsPropYear6$x), 10)),oob = rescale_none) +
		theme_bw() +
		theme(legend.position='none') 

if(max(allMethodsPropYear6x$x) > 10 | min(allMethodsPropYear6x$x) < -10) {

	subTemp = allMethodsPropYear6x[allMethodsPropYear6x$x > 10 | allMethodsPropYear6x$x < -10, ]
	tmps2 = as.numeric(gsub(pattern = 'Method', replacement = '', x = as.character(subTemp$METHOD)))
	tmps3 = tmps2 - 1
	tmps4 = 0.15*tmps3
	Labs = round(subTemp$x, 0)

	xPos = c(2,3) - 0.1 # minYear
	
	yPos = numeric(nrow(subTemp))
	yPos[subTemp$x > 10] = 9.6
	yPos[subTemp$x < -10] = -9.6

	mreyear2 = mreyear2 + annotate(geom="text", x= xPos, y = yPos, label = Labs, size = 3, angle = 90)

}


#  ------------------------------------------------------------------------
# TABLES MRE
allMethodsPropYear3 = NULL
for(k in seq_along(allFiles22)){
	tmp = read.csv(file.path(selFolder2, allFiles22[k]))
	allMethodsPropYear3 = rbind(allMethodsPropYear3, tmp)
}

# Get table 3: MRE total
allMethodsPropYear43 = aggregate(allMethodsPropYear3$MRE*100, list(METHOD = allMethodsPropYear3$METHOD), mean)
allMethodsPropYear43sd = aggregate(allMethodsPropYear3$MRE*100, list(METHOD = allMethodsPropYear3$METHOD), sd)

write.csv(allMethodsPropYear43, paste0('FinalTables/MRE_all_', scenarioName, '_mean.csv'), row.names = FALSE)
write.csv(allMethodsPropYear43sd, paste0('FinalTables/MRE_all_', scenarioName, '_sd.csv'), row.names = FALSE)


#  ------------------------------------------------------------------------
# TABLES and FIGURES MSE

# Table 1: MSE per age
allMethodsPropYear41 = aggregate(allMethodsPropYear3$MSE, list(AGE = allMethodsPropYear3$AGE, 
								METHOD = allMethodsPropYear3$METHOD), mean)
allMethodsPropYear41sd = aggregate(allMethodsPropYear3$MSE, list(AGE = allMethodsPropYear3$AGE, 
								METHOD = allMethodsPropYear3$METHOD), sd)
tab1 = cast(allMethodsPropYear41, METHOD ~ AGE)
tab1sd = cast(allMethodsPropYear41sd, METHOD ~ AGE)
write.csv(tab1, paste0('FinalTables/MSE_age_', scenarioName, '_mean.csv'), row.names = FALSE)
write.csv(tab1sd, paste0('FinalTables/MSE_age_', scenarioName, '_sd.csv'), row.names = FALSE)

# MSE per age
MSEdf = allMethodsPropYear41
MSEdf$x = MSEdf$x*10^5

mseage2 = 	ggplot(data=MSEdf, aes(x=as.factor(AGE), y=x, fill=as.factor(METHOD))) +
				geom_bar(stat="identity", color="black", position=position_dodge(), width = 0.7)+
				scale_fill_manual(values = c("#F8766D", "#7CAE00", "#00BFC4", '#C77CFF'), name = '') +
			  	ylab('') +
				xlab('') +
				scale_y_continuous(limits=c(0, min(max(MSEdf$x), 10)),oob = rescale_none) +
				theme_bw() +
				ggtitle("S / T") +
				theme(legend.position = 'none', plot.title = element_text(hjust = 0.5))


if(scenarioName == 'NoS_NoT') mseage2 = mseage2 + theme(legend.position=c(0.85, 0.7))


if(max(MSEdf$x) > 10) {

	subTemp = MSEdf[MSEdf$x > 10, ]
	tmps2 = as.numeric(gsub(pattern = 'Method', replacement = '', x = as.character(subTemp$METHOD)))
	tmps3 = tmps2 - 1
	tmps4 = 0.15*tmps3
	Labs = round(subTemp$x, 0)
	xPos = subTemp$AGE + tmps4 # minAge
	yPos = numeric(nrow(subTemp))
	yPos[subTemp$x > 10] = 9.7

	mseage2 = mseage2 + annotate(geom="text", x= xPos, y = yPos, label = Labs, size = 2.5, angle = 90)

}




# Get table 2: MSE per year
allMethodsPropYear41 = aggregate(allMethodsPropYear3$MSE, list(YEAR = allMethodsPropYear3$YEAR, 
								METHOD = allMethodsPropYear3$METHOD), mean)
allMethodsPropYear41sd = aggregate(allMethodsPropYear3$MSE, list(YEAR = allMethodsPropYear3$YEAR, 
								METHOD = allMethodsPropYear3$METHOD), sd)
tab1 = cast(allMethodsPropYear41, METHOD ~ YEAR)
tab1sd = cast(allMethodsPropYear41sd, METHOD ~ YEAR)
write.csv(tab1, paste0('FinalTables/MSE_year_', scenarioName, '_mean.csv'), row.names = FALSE)
write.csv(tab1sd, paste0('FinalTables/MSE_year_', scenarioName, '_sd.csv'), row.names = FALSE)


# MSE per year
MSEdf = allMethodsPropYear41
MSEdf$x = MSEdf$x*10^5

MSEdf$CAT = NA
MSEdf$CAT[MSEdf$YEAR < 2002] = '1994 - 2001'
MSEdf$CAT[MSEdf$YEAR > 2001 & MSEdf$YEAR < 2010] = '2002 - 2009'
MSEdf$CAT[MSEdf$YEAR > 2009] = '2010 - 2016'
MSEdfx = aggregate(MSEdf$x, list(METHOD = MSEdf$METHOD,
									PERIOD = MSEdf$CAT), FUN = mean)



mseyear2 = 	ggplot(data = MSEdfx, aes(y = x, x = factor(PERIOD), fill = factor(METHOD))) +
		geom_bar(stat="identity", color="black", position='dodge', width = 0.7)+
		scale_fill_manual(values = c("#F8766D", "#7CAE00", "#00BFC4", '#C77CFF'), name = '') +
	  	ylab('') +
		xlab('') +
		scale_y_continuous(limits=c(0, min(max(MSEdf$x), 10)),oob = rescale_none) +
		theme_bw() +
		ggtitle("S / T") +
		theme(legend.position='none', plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) 

if(max(MSEdf$x) > 10) {

	subTemp = MSEdfx[MSEdfx$x > 10, ]
	tmps2 = as.numeric(gsub(pattern = 'Method', replacement = '', x = as.character(subTemp$METHOD)))
	tmps3 = tmps2 - 1
	tmps4 = 0.15*tmps3
	Labs = round(subTemp$x, 0)

	xPos = c(1,2,3) - 0.1 # minYear
	
	yPos = numeric(nrow(subTemp))
	yPos[subTemp$x > 10] = 9.6

	mseyear2 = mseyear2 + annotate(geom="text", x= xPos, y = yPos, label = Labs, size = 3, angle = 90)

}



# Get table 3: MSE total
allMethodsPropYear43 = aggregate(allMethodsPropYear3$MSE, list(METHOD = allMethodsPropYear3$METHOD), mean)
allMethodsPropYear43sd = aggregate(allMethodsPropYear3$MSE, list(METHOD = allMethodsPropYear3$METHOD), sd)

write.csv(allMethodsPropYear43, paste0('FinalTables/MSE_all_', scenarioName, '_mean.csv'), row.names = FALSE)
write.csv(allMethodsPropYear43sd, paste0('FinalTables/MSE_all_', scenarioName, '_sd.csv'), row.names = FALSE)



# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

# Plot MSE and MRE per age:
bitmap(paste0('FinalFigures/MSE_MRE_age_', scenarioName, '_total.tiff'), height = 120, width = 180, units = 'mm', res = 700)

	grid.arrange(mseage, mseage2, mreage, mreage2, nrow = 2, bottom=textGrob("Age", gp=gpar(fontsize=10)))

dev.off()

# Plot MSE and MRE per year:
bitmap(paste0('FinalFigures/MSE_MRE_year_', scenarioName, '_total.tiff'), height = 120, width = 180, units = 'mm', res = 700)

	grid.arrange(mseyear, mseyear2, mreyear, mreyear2, nrow = 2)

dev.off()




# -----------------------------------------------------------------------------

