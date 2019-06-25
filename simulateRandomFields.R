# --------------------------------------------------------------
# Create randon numbers for recruitment
# Space and spatiotemporal components

n_years = length(allYears)

# Simualte RF
model_O = RandomFields::RMmatern(nu = NuMat, var=SD_O^2, scale=SpatialScale) # spatial model component
model_E = RandomFields::RMmatern(nu = NuMat, var=SD_E^2, scale=SpatialScale) # spatiotemporal model component

# Simulate Omega
Omega = RFsimulate(model = model_O, x=as.matrix(predictGrid2), grid = FALSE)
Omega1 = Omega@data[,1]

# Simulate Epsilon
n_stations = nrow(predictGrid2)
Epsilon1 = array(NA, dim=c(n_stations,n_years))
dimFig = ceiling(sqrt(n_years))

if(simulation){
	for(t in 1:n_years){
	  Epsilon = RFsimulate(model=model_E, x=as.matrix(predictGrid2), grid = FALSE)
	  Epsilon1[,t] = Epsilon@data[,1]
	}
}

# Plot both fields
	png(paste0('RandomField_Recs/RandomField_Rec_Omega', ix,'.png'), width = 800, height = 700, units = 'px')
	print(ggplot() + 
		geom_point(aes(predictGrid2$x, predictGrid2$y, color = Omega1)) +
		scale_colour_gradient2(high="red",mid = 'white', low='blue') +
		theme_bw())
	dev.off()

if(!simulation){



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
}

# Add temporal variation based on sigmaR
rRecTemp = rnorm(n = length(allYears), mean = -(sigmaR^2)/2, sd = sigmaR)
