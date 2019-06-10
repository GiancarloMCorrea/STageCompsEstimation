
# New Parameters:

# Read some data:
SelecSurvLen = read.csv('SurveySelecLen.csv') # read survey length selectivity
sampleStations = read.csv('sampledGrids.csv') # read sample stations. WARNING!!!! THIS SHOULD CHANGE IF THE GRID SIZE CHANGE!!!!!!

# Grid size:
dsGrid = 6 # in nautical miles # maybe the best value is 1.5, but for time reasons do it = 6

# Create the grid:

xgrid = seq(from = -179, to = -158, by = dsGrid/60)
ygrid = seq(from = 54, to = 62, by = dsGrid/60)
predictGrid = expand.grid(xgrid,ygrid)
names(predictGrid) = c("x","y")

surveyPolygon = read.csv('polygonSurvey.csv') # this is the survey standard area
predictGrid2 = predictGrid[which(point.in.polygon(point.x = predictGrid$x, point.y = predictGrid$y, 
                                                  pol.x = surveyPolygon$x, pol.y = surveyPolygon$y) == 1), ]
# 9376 rows for 6 mn grid. so, total area = 9376 * 36 = 337536 mn2 = 1159196.5 km2
# 36 nm2 = 123.6 km2. Factor = 3.43

# Rec density:
iniR0 = 451155e03 # This is a unique value for the first year. Average of R0 1994-2016 from SS3 results.
StudyArea = nrow(predictGrid2)*dsGrid*dsGrid*3.43 # in km2 = 1157748 

# R0inidengrid = R0inidenkm2*dsGrid*dsGrid*3.43 # abundance per grid.

# Scenario name:
scenarioName = 'HighS_HighT'

# These parameters will be used later:
iniYear = 1994
endYear = 2016
maxAge = 20
minAge = 0
maxLen = 120
minLen = 1
lenBin = 1
dT = 1
wS = 0.05
wT = 0.04
SpatialScale = 0.3
SD_O = 0.5
SD_E = 0.2
NuMat = 0.5
#thPr = 1e-03
maxNSamAge = 4 # number of ind sampled in station j for age (random sampling)
#maxNSamPerAge = 120 # max number of ind sampled for age in length bin l in the survey. THIS SHOULD NOT BE A RESTRICTION ? (random sampling)
#maxNSamAgeStation = 3 # max number of ind of length bin l sampled in station j for age (random sampling)
t1 = 0.3 # first dt for natural and fishing mortality
Linf = 100.6 # L2 in SS
K_par = 0.175 # 0.195 is the value in SS. seems to be very high.
M_par = 0.34 
CV1 = 3.45 # this is a sd. 3.45 is the value in SS
CV2 = 9.586 # this is a sd. 9.586 is the value in SS
L1_par = 10 # same as SS. L1
A1_par = 0.5 # same as SS. a3
A2_par = 18.5 # a4 in SS
F_par = 0.46 # SS results. Apical F
areaSwept = 0.05 # km2. average of all data: 1994 - 2016. a
nSamLoc = round(nrow(sampleStations)*0.95) # number of age sample locations. 95% all locations = 332 so far
sigmaR = 0.66# for recruitment
sigmaM = 1.3 # for simulated catches
SelecFish = c(0, 0.000231637, 0.00279779, 0.0328583, 0.29149, 0.832831, 0.983694, 0.998633, 0.999887, 
             1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1) # FISHERY. from age 0 to 20. be sure it has same length as allAges
SelecSurv = c(0, 1, 1, 1, 1, 1, 1, 1, 1, 
              1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)# SURVEY. from age 0 to 20. be sure it has same length as allAges

# Derived quantities:
allAges = seq(from = minAge, to = maxAge, by = dT)
allYears = seq(from = iniYear, to = endYear, by = dT)
allLens = seq(from = minLen, to = maxLen, by = lenBin)

# --------------------------------------------------------------
# Parameters for the estimation part

agePlus = 12
firstAgeCRL = 1

# --------------------------------------------------------------
# Create randon numbers for recruitment
# Space and spatiotemporal components

# run colors:
redA = alpha("red", 0.5)
blueA = alpha("blue", 0.5)

n_years = length(allYears)

# Simualte RF
model_O = RandomFields::RMmatern(nu = NuMat, var=SD_O^2, scale=SpatialScale) # spatial model component
model_E = RandomFields::RMmatern(nu = NuMat, var=SD_E^2, scale=SpatialScale) # spatiotemporal model component

# Simulate Omega
Omega = RFsimulate(model = model_O, x=as.matrix(predictGrid2), grid = FALSE)
Omega1 = Omega@data[,1]

png('RandomField_Recs/RandomField_Rec_Omega.png')
print(ggplot() + 
	geom_point(aes(predictGrid2$x, predictGrid2$y, color = Omega1)) +
	scale_colour_gradient2(high="red",mid = 'white', low='blue') +
	theme_bw())
dev.off()

# Simulate Epsilon
n_stations = nrow(predictGrid2)
Epsilon1 = array(NA, dim=c(n_stations,n_years))
dimFig = ceiling(sqrt(n_years))
pdf('RandomField_Recs/RandomField_Rec_Epsilon.pdf')
for(t in 1:n_years){
  Epsilon = RFsimulate(model=model_E, x=as.matrix(predictGrid2), grid = FALSE)
  Epsilon1[,t] = Epsilon@data[,1]
		print(ggplot() + 
			geom_point(aes(predictGrid2$x, predictGrid2$y, color = Epsilon1[,t])) +
			scale_colour_gradient2(high="red",mid = 'white', low='blue') +
			labs(color = allYears[t]) +
			theme_bw())
}
dev.off()
	

# Old version
# in Time:
#rRecTemp = rnorm(n = length(allYears), mean = -(sigmaR^2)/2, sd = sigmaR)
# in Space
#sim12 = grf(n = nrow(predictGrid2), grid = as.matrix(predictGrid2), cov.pars=c(0.25, 1), nug = 0, nsim = length(allYears)) # this is the grf

# --------------------------------------------------------------
# Random Fields for growth parameters 
gDummy2 = gstat(formula=z~1+x+y, locations=~x+y, dummy=T, beta=c(0,0.1,0.3), 
                 model=vgm(psill=0.05, range=5, model='Mat'), nmax = 15)
yy2 = predict(gDummy2, newdata=predictGrid2, nsim=1)
gridded(yy2) = ~x+y
yy2$sim1 = normalize(x = yy2$sim1, method = 'range', range = c(wS,-1*wS)) # this value (0.05) give us a difference of 5 cm between strata.

ak = map_data('worldHires','USA:Alaska')
ak = ak[ak$long < 0, ]
bitmap('RandomField_K.tiff', height = 65, width = 130, units = 'mm', res = 600)
print(map.heatmap(lat = yy2@coords[,2], lon = yy2@coords[,1], yy2@data,
              color_low = "blue", color_high = "red", zeroiswhite = TRUE, xlim = c(-179,-158), ylim = c(54,62.5)) +
			  geom_polygon(data = ak, aes(long, lat, group = group), 
			  fill = 8, color="black") +
			  xlab('longitude') +
			  ylab('latitude') +
			  #xlim(-180,-156) +
			  theme(plot.margin = unit(c(0,0,0,0),"cm")))
dev.off()  
  
# Temporal trend in K parameter:
KparT = normalize(x = allYears, method = 'range', range = c(-1*wT,wT))
