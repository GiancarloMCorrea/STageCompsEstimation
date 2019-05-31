#  ------------------------------------------------------------------------
#  plot spatial variability in growth (include sex):
data2 = read.csv("simData/paccod_catch_Sim.csv")
data3 = read.csv("simData/paccod_len_Sim.csv")
data4 = read.csv("simData/paccod_age_Sim.csv")

#  ------------------------------------------------------------------------
# Temporal variability

# Classic palette BuPu, with 4 colors
coul = rev(brewer.pal(n = 9, name = "Spectral"))
 
# I can add more tones to this palette :
coul2 = colorRampPalette(coul)(length(unique(data4$YEAR)))
#coul2 = addalpha(coul2, alpha = 0.3)
#coul2 = colorRampPaletteAlpha(coul, length(unique(data4$YEAR)))

data4xx = data4[data4$AGE < 7, ]
png('checkTemporalVariability.png', height = 650, width = 1000, units = 'px', res = 110)
print(ggplot(data4xx, aes(x = factor(AGE), y = LENGTH, fill = factor(YEAR))) +
       geom_boxplot(outlier.size = 0.3) +
	   ylim(0,110))
dev.off()

bitmap('checkTemporalVariability2.tiff', height = 80, width = 100, units = 'mm', res = 800)
print(ggplot(data4xx, aes(x = factor(AGE), y = LENGTH, fill = factor(YEAR), color = factor(YEAR))) +
       geom_boxplot(outlier.size = 0.2, width = 0.8, position=position_dodge(width=0.6)) +
	   scale_fill_manual(values = coul2) +
	   scale_color_manual(values = coul2) +
	   theme_bw() +
	   #scale_colour_gradientn(colours = coul2, values = seq(1994,2016,1)) +
	   theme(legend.position = c(0.13, 0.69), legend.text = element_text(size=6), legend.title = element_blank()) +
	   #guides(fill = FALSE, color = FALSE) +
	   #ylim(0,110) +
	   #labs(color = "") +
	   xlab('Age') +
	   ylab('Length (cm)')
	   )
dev.off()

#  ------------------------------------------------------------------------
# Spatial Variability

gradColors = rev(brewer.pal(n = 9, name = "Spectral"))

data42 = normDF(dat = data4)
data42xx = data42[data42$AGE > 0 & data42$AGE < 7, ]
data42xx$AGE = as.factor(data42xx$AGE)

png('checkSpatialVariability.png', height = 650, width = 800, units = 'px', res = 110)
print(ggplot(data42, aes(LON, LAT)) +
        geom_point(aes(color = LENGTH), size = 1.5) +
        scale_colour_gradientn(colours = gradColors) +
        facet_wrap( ~ factor(AGE), ncol = 3))
dev.off()

bitmap('checkSpatialVariability2.tiff', height = 120, width = 130, units = 'mm', res = 800)
print(ggplot() +
		geom_polygon(data = ak, aes(long, lat, group = group), 
			   fill = 8, color="black") +
		geom_point(data = data42xx, mapping = aes(LON, LAT, color = LENGTH), size = 0.4) +
		scale_colour_gradientn(colours = gradColors) +
		facet_wrap( ~ AGE, ncol = 2) +
		coord_fixed(ratio = 1.1, xlim = c(-179,-158), ylim = c(54,62.5)) +				   
		theme_bw() +
		xlab('longitude') +
		ylab('latitude') +
		theme(legend.title = element_blank()) +
		theme(plot.margin = unit(c(0,0,0,0),"cm")))
dev.off() 
