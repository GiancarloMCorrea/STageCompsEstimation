# SET YOUR WORKING DIRECTORY !!

rm(list = ls())

# Call function
source('CRLfunction.R')

# --------------------------------------------------------------
# Read data:
cdata = read.csv('POLLOCK_CATCH_2015_2019.csv')
ldata = read.csv('POLLOCK_LENGTH_2015_2019.csv')
adata = read.csv('POLLOCK_AGE_2015_2019.csv')

# ---------------------------------------------------------------------
# Example 1: 

newData = estimateAgeCRL(DataModel = adata, DataEstimation = ldata, FormulaGAM = 's(LENGTH) + s(LON, LAT, k = 10)', 
                         AgeMin = 1, AgeMax = 9, AgeVariable = 'AGE', TimeVariable = 'YEAR')


