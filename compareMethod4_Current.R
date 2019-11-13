# Read data:
require(reshape2)
require(ggplot2)
setwd(dir = 'C:/Users/moroncog/Documents/GitHub/STageCompsEstimation')

#  ------------------------------------------------------------------------
rm(list = ls())

method4 = read.csv(file = 'realData/Method4_pcod.csv')
current = read.csv(file = 'realData/Method_Current_pcod.csv')

colnames(method4) = c('YEAR', '1', '2', '3', '4', '5', '6', '7', '8')
colnames(current) = c('YEAR', '1', '2', '3', '4', '5', '6', '7', '8')

met4df = melt(method4, id = c('YEAR'))
met4df$Method = 'Method4'

metcdf = melt(current, id = c('YEAR'))
metcdf$Method = 'Current'

mergedf = rbind(met4df, metcdf)
colnames(mergedf) = c('YEAR', 'AGE', 'FREQUENCY', 'Method')

# -----------------------------------------------------------------------
# Plot all methods:
mergedf$AGE = as.numeric(as.character(mergedf$AGE))

bitmap('FinalFigures/compare_Met4_Current.tiff', height = 190, width = 190, units = 'mm', res = 900)
	print(ggplot(mergedf, aes(x = AGE, y = FREQUENCY)) +
	  geom_line(aes(color = factor(Method))) +
	  facet_wrap( ~ factor(YEAR), nrow = 5) +
	  xlab('Age') +
	  ylab('Proportion of abundance') +
	  scale_x_discrete(limits = 1:8) +
	  scale_color_manual(values=c("black", "#C77CFF")) +
	  theme_bw() +
	  theme(legend.position = c(0.8, 0.1), legend.title = element_blank(), legend.text = element_text(size = 12)))
 dev.off()