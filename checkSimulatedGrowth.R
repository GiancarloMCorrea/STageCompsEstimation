require(RColorBrewer)
require(ggplot2)
library(mapdata)
require(BBmisc)
library(gridExtra)

setwd(dir = 'C:/Users/moroncog/Documents/GitHub/STageCompsEstimation')

source('auxFunctionsSimulation.R')
minFigAge = 1
maxFigAge = 6

#  ------------------------------------------------------------------------
#  plot spatial variability in growth:
scenarioName = 'HighS_HighT'

data2 = read.csv(paste0('simData/paccod_catch_Sim_', scenarioName, '.csv'))
data3 = read.csv(paste0('simData/paccod_len_Sim_', scenarioName, '.csv'))
data4 = read.csv(paste0('simData/paccod_age_Sim_', scenarioName, '.csv'))

#  ------------------------------------------------------------------------
# Temporal variability

# Classic palette BuPu, with 4 colors
#coul = rev(brewer.pal(n = 9, name = "Spectral"))
coul = c(rev(brewer.pal(n = 9, name = "Blues")), brewer.pal(n = 9, name = "Reds"))

 
# I can add more tones to this palette :
coul2 = colorRampPalette(coul)(length(unique(data4$YEAR)))

# select ages:
data4xx = data4[data4$AGE < 7, ]

# bitmap(paste0('checkTemporalVariability2_', scenarioName, '.tiff'), height = 80, width = 100, units = 'mm', res = 800)
# print(ggplot(data4xx, aes(x = factor(AGE), y = LENGTH, fill = factor(YEAR), color = factor(YEAR))) +
       # geom_boxplot(outlier.size = 0.2, width = 0.8, position=position_dodge(width=0.6)) +
	   # scale_fill_manual(values = coul2) +
	   # scale_color_manual(values = coul2) +
	   # theme_bw() +
	   # theme(legend.position = c(0.13, 0.69), legend.text = element_text(size=6), legend.title = element_blank()) +
	   # xlab('Age') +
	   # ylab('Length (cm)'))
# dev.off()


bitmap(paste0('FinalFigures/checkTemporalVariability2_', scenarioName, '.tiff'), height = 100, width = 190, units = 'mm', res = 800)
print(ggplot(data4xx, aes(x = factor(YEAR), y = LENGTH)) +
       geom_boxplot(outlier.size = 0.2, width = 0.8, position=position_dodge(width=0.6),  fill = '#A4A4A4') +
	   #scale_fill_manual(values = coul2) +
	   #scale_color_manual(values = coul2) +
	   facet_wrap( ~ AGE, ncol = 3, scales = 'free_y') +
	   theme_bw() +
	   theme(legend.position = 'none') +
	   xlab('Year') +
	   scale_x_discrete(name = "Year", 
                    labels=c(rep(' ', 1), 1995, rep(' ', 4), 2000, rep(' ', 4), 2005, rep(' ', 4), 2010, rep(' ', 4), 2015, rep(' ', 1))) +
	   #theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
	   ylab('Length (cm)'))
dev.off()

#  ------------------------------------------------------------------------
# Spatial Variability

#gradColors = rev(brewer.pal(n = 9, name = "Spectral"))
gradColors = c(rev(brewer.pal(n = 9, name = "Blues")), brewer.pal(n = 9, name = "Reds"))

#data42 = normDF(dat = data4)
data42 = data4
# data42xx = data42[data42$AGE > 0 & data42$AGE < 7, ]
# data42xx$AGE = as.factor(data42xx$AGE)
ak = map_data('worldHires','USA:Alaska')
ak = ak[ak$long < 0, ]
# #temp2 = aggregate(data4xx$LENGTH, list(LON = data4xx$LON, LAT = data4xx$LAT, AGE = data4xx$AGE), mean)

# bitmap(paste0('FinalFigures/checkSpatialVariability2_', scenarioName, '.tiff'), height = 100, width = 190, units = 'mm', res = 800)
# print(ggplot() +
# 		geom_polygon(data = ak, aes(long, lat, group = group), 
# 			   fill = 8, color="black") +
# 		geom_point(data = data42xx, mapping = aes(LON, LAT, color = LENGTH), size = 0.4) +
# 		scale_colour_gradientn(colours = gradColors) +
# 		facet_wrap( ~ AGE, ncol = 3) +
# 		coord_fixed(ratio = 1.1, xlim = c(-179,-158), ylim = c(54,62.5)) +				   
# 		theme_bw() +
# 		xlab('longitude') +
# 		ylab('latitude') +
# 		theme(legend.position = 'none') +
# 		theme(plot.margin = unit(c(0,0,0,0),"cm")))
# dev.off() 

figList = list()
for(k in minFigAge:maxFigAge){

data42xx1 = data42[data42$AGE == k, ]
data42xx1 = aggregate(data42xx1$LENGTH, list(LON = data42xx1$LON, LAT = data42xx1$LAT), mean)
names(data42xx1) = c('LON', 'LAT', 'LENGTH')

if(k %in% c(4,5,6)) { 
	xLab = 'longitude' 
} else {
	xLab = '' 
}

if(k %in% c(1,4)) { 
	yLab = 'latitude' 
} else {
	yLab = '' 
}

	figList[[k]] = ggplot() +
					geom_polygon(data = ak, aes(long, lat, group = group), 
						   fill = 8, color="black") +
					geom_point(data = data42xx1, mapping = aes(LON, LAT, color = LENGTH), size = 0.4) +
					scale_colour_gradientn(colours = gradColors) +
					coord_fixed(ratio = 1.1, xlim = c(-180,-158), ylim = c(52,62.5)) +				   
					theme_bw() +
					xlab(xLab) + 
					ylab(yLab) +
					annotate(geom="text", x=-160, y=62, label=paste0('Age ', k)) +
					theme(legend.position = c(0.11, 0.32), legend.title = element_blank(), legend.key.size = unit(0.2, "cm"), legend.text = element_text(size = 6)) +
					theme(plot.margin = unit(c(0,0,0,0),"cm"))

}

bitmap(paste0('FinalFigures/checkSpatialVariability3_', scenarioName, '.tiff'), height = 100, width = 190, units = 'mm', res = 800)
	do.call("grid.arrange", c(figList, ncol = ceiling(maxFigAge/2)))
dev.off()


#  ------------------------------------------------------------------------
# Make plot for growth variability in a cohort:

len1 = 118.6*(1-exp(-(0.1376)*(0+0.168)))
len2 = len1 + (len1-118.6)*(exp(-(4/12)*(0.1376)) - 1)

len3 = 118.6*(1-exp(-(0.1076)*(2+0.168)))
len4 = 118.6*(1-exp(-(0.1676)*(2+0.168)))


# No S scenario:
simu1 = rnorm(100000, mean = 40, sd = 1.2)
dat1 = data.frame(cohort1 = simu1)

ap1 = ggplot(dat1) + 
	geom_density(aes(x = cohort1), fill = "gray") +
	theme_bw() +
	xlab('') +
	ylab('No S scenario') +
	xlim(30, 47.5) +
	theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) 

ap2 = ggplot(dat1) + 
	geom_density(aes(x = cohort1), linetype = "dashed") +
	geom_text(x = 33, y = 0.3, label="Normal growth", color = 'black', size = 2.6) +
	theme_bw() +
	xlab('') +
	ylab('') +
	xlim(30, 47.5) +
	theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

# S scenario:
simu1 = rnorm(100000, mean = 40, sd = 1.2)
simu2 = rnorm(100000, mean = 37.5, sd = 1.2)
simu3 = rnorm(100000, mean = 42.5, sd = 1.2)
dat1 = data.frame(cohort1 = simu1)
dat2 = data.frame(cohort1 = simu2)
dat3 = data.frame(cohort1 = simu3)
dat4 = data.frame(cohort1 = c(simu1, simu2, simu3))

ap3 = ggplot(dat4) + 
	geom_density(aes(x = cohort1), fill = "gray") +
	theme_bw() +
	xlab('Length (cm)') +
	ylab('S scenario') +
	xlim(30, 47.5) +
	theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

ap4 = ggplot(dat1) + 
	geom_density(aes(x = cohort1), linetype = "dashed") +
	geom_density(data = dat2, aes(x = cohort1), color = 'blue', linetype = "dashed") +
	geom_density(data = dat3, aes(x = cohort1), color = 'red', linetype = "dashed") +
		geom_text(x = 33, y = 0.3, label="Normal growth", color = 'black', size = 2.6) +
		geom_text(x = 33, y = 0.33, label="Slow growth", color = 'blue', size = 2.6) +
		geom_text(x = 33, y = 0.27, label="Fast growth", color = 'red', size = 2.6) +
	theme_bw() +
	xlab('Length (cm)') +
	ylab('') +
	xlim(30, 47.5) +
	theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())


bitmap('FinalFigures/auxPlotGrowthSpatial.tiff', height = 120, width = 120, units = 'mm', res = 600)
	
	grid.arrange(ap1, ap2, ap3, ap4, nrow = 2)
			
dev.off()
