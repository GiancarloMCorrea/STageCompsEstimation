# SCRIPT TO COMPARE AGE ABUNDANCE ESTIMATES BY ALK AND CRL
rm(list = ls())

# Libraries
require(mgcv)
require(fishmethods)
require(reshape2)
require(ggplot2)


# Some parameters 
sel_year = 2017
min_age = 1
age_plus = 12
sel_station = c('A-05', 'G-15', 'H-15', 'J-15', 'M-30', 'J-10', 'E-11', 'K-08', 'H-07')

# --------------------------------------------------------------
# Read data:
cdata = read.csv('POLLOCK_CATCH_2015_2019.csv')
ldata = read.csv('POLLOCK_LENGTH_2015_2019.csv')
adata = read.csv('POLLOCK_AGE_2015_2019.csv')

# subset data:
cdata = cdata[cdata$YEAR == sel_year, ]
ldata = ldata[ldata$YEAR == sel_year, ]
adata = adata[adata$YEAR == sel_year, ]

# Analyze age subsamples: IMPORTANT STEP TO DECIDE THE RIGHT AGE PLUS GROUP:
# This might vary across years
#write.csv(table(adata$LENGTH, adata$AGE), 'Examine.csv')

# Estimate C_L per station in length subsample data:
data5 = aggregate(ldata$FREQUENCY, list(STATIONID = ldata$STATIONID), sum)
ldata$N_CATCH = cdata$NUMBER_FISH[match(ldata$STATIONID, cdata$STATIONID)]
ldata$N_SUBSAMPLE = data5$x[match(ldata$STATIONID, data5$STATIONID)]
ldata$LAMBDA = ldata$N_SUBSAMPLE/ldata$N_CATCH
ldata$C_l_i = ldata$FREQUENCY/ldata$LAMBDA

# Apply age plus group:
adata$AGE = ifelse(test = adata$AGE > age_plus, age_plus, adata$AGE) 

# character STATIONID:
cdata$STATIONID = as.character(cdata$STATIONID)
ldata$STATIONID = as.character(ldata$STATIONID)
adata$STATIONID = as.character(adata$STATIONID)


# ALK approach --------------------------------------------------

# Create length vs age freq (entire age subsample data):
ALK_tmp = table(adata$LENGTH, adata$AGE)
len = as.numeric(rownames(ALK_tmp))
nl = numeric(length(len))

# Select a station
tmp_ldata = ldata[ldata$STATIONID %in% sel_station, ]
tmp_station = aggregate(tmp_ldata$FREQUENCY, list(STATIONID = tmp_ldata$STATIONID, LENGTH = tmp_ldata$LENGTH), sum)

# Make a function to calculate Proportions by station
alk_calc = function(station_data){
  
  find_pos = match(station_data$LENGTH, len)
  find_pos_2 = find_pos[!is.na(find_pos)]
  nl[find_pos_2] = station_data$x[!is.na(find_pos)]
  
  alkData = cbind(len, nl, ALK_tmp)
  alkData = as.data.frame(alkData)
  prop_station = alkprop(alkData)
  return(data.frame(AGE = min_age:age_plus, PROP = prop_station$results$prop))
  
}  

# Find propotions by station
prop_by_station = by(data = tmp_station, INDICES = tmp_station$STATIONID, FUN = alk_calc)

# Make data.frame from list
prop_data = do.call(rbind, prop_by_station)
prop_data$STATIONID = gsub("\\..*","",rownames(prop_data))

# Calculate C_a_i:
prop_data$CATCH = cdata$NUMBER_FISH[ match(prop_data$STATIONID, cdata$STATIONID)  ]
prop_data$C_a_i = prop_data$CATCH * prop_data$PROP

# Order data frame:
alk_data = prop_data[,c('STATIONID', 'AGE', 'C_a_i')]
alk_data$METHOD = 'ALK'




# CRL approach ------------------------------------------------------------

# Create new data frame:
tmp_adata = adata

# Create matrix to save results
all_ages = sort(unique(tmp_adata$AGE))
PropAgeMat = matrix(NA, ncol = length(all_ages), nrow = 1) 
colnames(PropAgeMat) = all_ages

# run the model GAM:
matPreds = matrix(NA, ncol = length(all_ages), nrow = nrow(ldata))
for(ii in seq_along(all_ages)){
  
  if(ii == length(all_ages)){
    predvals = rep(1, times = nrow(ldata))
  } else {
    
    tmp_adata$AGEFAC = ifelse(test = tmp_adata$AGE > all_ages[ii], 0, 1)
    modtmp = gam(AGEFAC ~ s(LENGTH) + s(LON, LAT, k = 10), family = binomial,
                 data = tmp_adata)
    predtmp = predict(modtmp, newdata = ldata, type = 'response')
    predvals = as.vector(predtmp)
    elimina = which(tmp_adata$AGEFAC == 1)
    if(length(elimina) > 0) {
      tmp_adata = tmp_adata[-which(tmp_adata$AGEFAC == 1), ]
    } else {
      tmp_adata = tmp_adata
    }
    
  }
  
  matPreds[,ii] = predvals
  
}

matPreds2 = matrix(NA, ncol = length(all_ages), nrow = nrow(ldata))
for(kk in seq_along(all_ages)){
  
  if(kk == 1){
    matPreds2[,kk] = matPreds[,kk]
  } else {
    mattmp = 1 - as.matrix(matPreds[,(kk-1):1])
    matPreds2[,kk] =  matPreds[,kk]*apply(X = mattmp, MARGIN = 1, FUN = prod) 
  }
  
}


tempo2 = sweep(matPreds2, MARGIN = 1, ldata$C_l_i, `*`)
rownames(tempo2) = ldata$STATIONID
colnames(tempo2) = min_age:age_plus

tmp_station = tempo2[rownames(tempo2) %in% sel_station, ]
tmp_station = reshape2::melt(tmp_station)
crl_data = aggregate(tmp_station$value, list(STATIONID = tmp_station$Var1, AGE = tmp_station$Var2), sum)
colnames(crl_data)[3] = 'C_a_i'
crl_data$METHOD = 'CRL'


# MERGE ALK AND CRL DATASETS: ------------------------------------------
plot_data = rbind(alk_data, crl_data)


# Plot: ----------------------------------------------------------------

ggplot(data = plot_data, aes(AGE, C_a_i, color = METHOD)) +
  geom_line() +
  theme_bw() +
  xlab('Age') +
  labs(color = "Method") +
  scale_x_continuous(breaks = min_age:age_plus) +
  facet_wrap(~ STATIONID, scales = 'free_y')

