#  ------------------------------------------------------------------------
# Read again data:

data2 = read.csv("simData/paccod_catch_Sim.csv")
data3 = read.csv("simData/paccod_len_Sim.csv")
data4 = read.csv("simData/paccod_age_Sim.csv")
#outdat2 = read.csv('simData/abunlen_stratum2Sim.csv')

#  ------------------------------------------------------------------------
# Compare simulated total abundance and estimated (from sampling) total abundance:
tmp1 = by(data = outdat2$lenabun, INDICES = outdat2$year, FUN = sum)
yrs = as.numeric(names(tmp1))
abun = as.vector(tmp1)

NAgeYearMatrix = read.csv('simData/NAgeYearMat.csv')
abunreal = rowSums(NAgeYearMatrix)/1e6

TabunCompare = rbind(data.frame(year = yrs, abun = abun, type1 = 'estimated'),
					 data.frame(year = allYears, abun = abunreal, type1 = 'real'))

png('compareTempTotAbundance.png', height = 550, width = 850, units = 'px', res = 120)
print(ggplot(TabunCompare, aes(year, abun)) +
        geom_line(aes(color = factor(type1))) +
		xlab('Year') +
		ylab('Total Abundance (1e06)') +
		theme_bw())
dev.off()

abunreal2 = rowSums(NAgeYearMatrix[,-1])/1e6

TabunCompare2 = rbind(data.frame(year = yrs, abun = abun, type1 = 'estimated'),
					 data.frame(year = allYears, abun = abunreal2, type1 = 'real'))

png('compareTempTotAbundance_no0.png', height = 550, width = 850, units = 'px', res = 120)
print(ggplot(TabunCompare2, aes(year, abun)) +
        geom_line(aes(color = factor(type1))) +
		xlab('Year') +
		ylab('Total Abundance (1e06)') +
		theme_bw())
dev.off()

recdf = data.frame(years = allYears, rec = NAgeYearMatrix[,1])

png('temporalRecruitmentVar.png', height = 550, width = 850, units = 'px', res = 120)
print(ggplot(recdf, aes(years, rec)) +
        geom_line() +
		xlab('Year') +
		ylab('Recruitment') +
		theme_bw())
dev.off()

