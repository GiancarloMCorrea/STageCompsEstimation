require(RColorBrewer)
require(ggplot2)
library(mapdata)
require(BBmisc)
source('auxFunctionsSimulation.R')

#  ------------------------------------------------------------------------
#  plot spatial variability in growth (include sex):
data2 = read.csv("simData/paccod_catch_Sim.csv")
data3 = read.csv("simData/paccod_len_Sim.csv")
data4 = read.csv("simData/paccod_age_Sim.csv")
scenarioName = 'HighS_HighT'

#  ------------------------------------------------------------------------
# Temporal variability

# Classic palette BuPu, with 4 colors
coul = rev(brewer.pal(n = 9, name = "Spectral"))
 
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


bitmap(paste0('checkTemporalVariability2_', scenarioName, '.tiff'), height = 100, width = 190, units = 'mm', res = 800)
print(ggplot(data4xx, aes(x = factor(YEAR), y = LENGTH, fill = factor(YEAR))) +
       geom_boxplot(outlier.size = 0.2, width = 0.8, position=position_dodge(width=0.6)) +
	   scale_fill_manual(values = coul2) +
	   scale_color_manual(values = coul2) +
	   facet_wrap( ~ AGE, ncol = 3, scales = 'free_y') +
	   theme_bw() +
	   theme(legend.position = 'none') +
	   xlab('Year') +
	   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
	   ylab('Length (cm)'))
dev.off()

#  ------------------------------------------------------------------------
# Spatial Variability

gradColors = rev(brewer.pal(n = 9, name = "Spectral"))

data42 = normDF(dat = data4)
data42xx = data42[data42$AGE > 0 & data42$AGE < 7, ]
data42xx$AGE = as.factor(data42xx$AGE)
ak = map_data('worldHires','USA:Alaska')
ak = ak[ak$long < 0, ]
#temp2 = aggregate(data4xx$LENGTH, list(LON = data4xx$LON, LAT = data4xx$LAT, AGE = data4xx$AGE), mean)

bitmap(paste0('checkSpatialVariability2_', scenarioName, '.tiff'), height = 100, width = 190, units = 'mm', res = 800)
print(ggplot() +
		geom_polygon(data = ak, aes(long, lat, group = group), 
			   fill = 8, color="black") +
		geom_point(data = data42xx, mapping = aes(LON, LAT, color = LENGTH), size = 0.4) +
		scale_colour_gradientn(colours = gradColors) +
		facet_wrap( ~ AGE, ncol = 3) +
		coord_fixed(ratio = 1.1, xlim = c(-179,-158), ylim = c(54,62.5)) +				   
		theme_bw() +
		xlab('longitude') +
		ylab('latitude') +
		theme(legend.position = 'none') +
		theme(plot.margin = unit(c(0,0,0,0),"cm")))
dev.off() 
