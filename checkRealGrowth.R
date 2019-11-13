require(RColorBrewer)
require(ggplot2)
library(mapdata)
require(BBmisc)
library(gridExtra)

#setwd(dir = 'C:/Users/moroncog/Documents/GitHub/STageCompsEstimation')

source('auxFunctionsSimulation.R')
minFigAge = 1
maxFigAge = 6

#  ------------------------------------------------------------------------
#  plot spatial variability in growth:
scenarioName = 'Real'

data2 = read.csv('realData/paccod_catch_1994_2016.csv')
data3 = read.csv('realData/paccod_length_1994_2016.csv')
data4 = read.csv('realData/paccod_age_1994_2016.csv')

# fix some mistakes with real data:
data4$LENGTH = data4$LENGTH/10
data4 = data4[!is.na(data4$AGE), ]
#  ------------------------------------------------------------------------
# Temporal variability

# Classic palette BuPu, with 4 colors
#coul = rev(brewer.pal(n = 9, name = "Spectral"))
coul = c(rev(brewer.pal(n = 9, name = "Blues")), brewer.pal(n = 9, name = "Reds"))

 
# I can add more tones to this palette :
coul2 = colorRampPalette(coul)(length(unique(data4$YEAR)))

# select ages:
data4xx = data4[data4$AGE < 7 & data4$AGE > 0, ]

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


bitmap('FinalFigures/checkRealTemporalVariability2.tiff', height = 100, width = 190, units = 'mm', res = 800)
print(ggplot(data4xx, aes(x = factor(YEAR), y = LENGTH, fill = factor(YEAR))) +
       geom_boxplot(outlier.size = 0.2, width = 0.8, position=position_dodge(width=0.6)) +
	   scale_fill_manual(values = coul2) +
	   scale_color_manual(values = coul2) +
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
# gradColors = c(rev(brewer.pal(n = 9, name = "Blues")), brewer.pal(n = 9, name = "Reds"))

# #data42 = normDF(dat = data4)
# data42 = data4
# ak = map_data('worldHires','USA:Alaska')
# ak = ak[ak$long < 0, ]

# data2$STATIONID = as.character(data2$STATIONID)
# data42$STATIONID = as.character(data42$STATIONID)

# data2$ID_HAUL = paste0(data2$YEAR, data2$STATIONID)
# data42$ID_HAUL = paste0(data42$YEAR, data42$STATIONID)


# # seeking for lon lat 
# data42$LON = data2$START_LONGITUDE[match(data42$ID_HAUL,
#                                         data2$ID_HAUL)]
# data42$LAT = data2$START_LATITUDE[match(data42$ID_HAUL,
#                                        data2$ID_HAUL)]


# figList = list()
# for(k in minFigAge:maxFigAge){

# data42xx1 = data42[data42$AGE == k, ]
# data42xx1 = aggregate(data42xx1$LENGTH, list(LON = data42xx1$LON, LAT = data42xx1$LAT), mean)
# names(data42xx1) = c('LON', 'LAT', 'LENGTH')

# 	figList[[k]] = ggplot() +
# 					geom_polygon(data = ak, aes(long, lat, group = group), 
# 						   fill = 8, color="black") +
# 					geom_point(data = data42xx1, mapping = aes(LON, LAT, color = LENGTH), size = 0.4) +
# 					scale_colour_gradientn(colours = gradColors) +
# 					coord_fixed(ratio = 1.1, xlim = c(-180,-158), ylim = c(52,62.5)) +				   
# 					theme_bw() +
# 				    xlab('') +
# 					ylab('latitude') +
# 					annotate(geom="text", x=-160, y=62, label=paste0('Age ', k)) +
# 					theme(legend.position = c(0.11, 0.32), legend.title = element_blank(), legend.key.size = unit(0.2, "cm"), legend.text = element_text(size = 6)) +
# 					theme(plot.margin = unit(c(0,0,0,0),"cm"))

# }

# bitmap('FinalFigures/checkRealSpatialVariability3.tiff', height = 80, width = 190, units = 'mm', res = 800)
# 	do.call("grid.arrange", c(figList, ncol = ceiling(maxFigAge/2)))
# dev.off()

data4tmp = data4[data4$AGE > 0 & data4$AGE < 13, ]

# For first plot:
require(mgcv)

mod1 = gam(formula = LENGTH ~ factor(AGE), data = data4tmp, family = 'tw')
pred1 = predict(object = mod1, se.fit = TRUE, type = 'response')
uniqpred = aggregate(pred1, list(AGE = data4tmp$AGE), unique)
data4tmp$RESIDUALS = data4tmp$LENGTH - pred1$fit
data4tmp$COLOR = NA 
data4tmp$COLOR[data4tmp$RESIDUALS > 0] = 2
data4tmp$COLOR[data4tmp$RESIDUALS < 0] = 4

# For second plot:
data4tmp2 = aggregate(data4tmp$RESIDUALS, list(STATIONID  = data4tmp$STATIONID), mean)
tmpxlon = aggregate(data4tmp$LON, list(STATIONID  = data4tmp$STATIONID), max)
tmpxlat = aggregate(data4tmp$LAT, list(STATIONID  = data4tmp$STATIONID), max)
data4tmp2$LON = tmpxlon$x
data4tmp2$LAT = tmpxlat$x

data4tmp2$COLOR = NA 
data4tmp2$COLOR[data4tmp2$x > 0] = 2
data4tmp2$COLOR[data4tmp2$x < 0] = 4
data4tmp2$x2 = abs(data4tmp2$x)


# make plot:
bitmap('FinalFigures/checkRealSpatialVariability4.tiff', height = 80, width = 190, units = 'mm', res = 900)
par(mfrow = c(1,2))

par(mar = c(5,6,0.5,0.5))
plot(NA, NA, xlim = c(1,12), ylim = c(5, 120), xlab = '', ylab = '', axes = FALSE)
points(data4tmp$AGE, data4tmp$LENGTH, col = data4tmp$COLOR)
lines(uniqpred$AGE, uniqpred$fit)
lines(uniqpred$AGE, uniqpred$fit + uniqpred$se, lty = 2)
lines(uniqpred$AGE, uniqpred$fit - uniqpred$se, lty = 2)
axis(1, at = 1:12, labels = 1:12, cex.axis = 1.5)
axis(2, las = 2, cex.axis = 1.7)
mtext(text = 'Age', side = 1, line = 2.8, cex = 2)
mtext(text = 'Length (cm)', side = 2, line = 3.6, cex = 2)
box()

par(mar = c(5,6,0.5,0.5))
plot(NA, NA, xlim = c(-179,-158), ylim = c(54.5, 62.5), xlab = '', ylab = '', axes = FALSE)
map("worldHires", add=T, fill=T, col=8)
points(data4tmp2$LON, data4tmp2$LAT, col = data4tmp2$COLOR, cex = (data4tmp2$x2/max(data4tmp2$x2))*2.5, pch = 19)
points(c(-178, -178, -178), c(56,55.5,55), cex = c(1.729, 1.153, 0.576), pch = 19)
text(c(-176, -176, -176), c(56,55.5,55), labels = paste0(c(6,4,2), ' cm'), cex = 1.8)
axis(1, cex.axis = 1.7)
axis(2, las = 2, cex.axis = 1.7)
mtext(text = 'longitude', side = 1, line = 2.8, cex = 2)
mtext(text = 'latitude', side = 2, line = 3.6, cex = 2)
box()

dev.off()
