###################################################################################################################
#### Run the aggregate and ecological models on a test data set with constant and varying pollution-disease effects
###################################################################################################################
#### Specifically the following 4 models are fitted to the sample data:
## Ecological model (1) with a constant air pollution and healh effect - using the R package CARBayes.
## Ecological model (1) with a varying air pollution and healh effect - using the R package CARBayes.
## Aggregate model (8) with a constant air pollution and health effect - using poisson.lerouxCARaggregate.R.
## Aggregate model (8) with a varying air pollution and health effect - using poisson.lerouxCARaggregatevcm.R.


################################################################
#### Load the libraries and functions required to run the models
################################################################
#### Libraries
library(sp)
library(spdep)
library(rgdal)
library(rgeos)
library(CARBayes)
library(coda)
library(spam)
library(truncdist)
library(truncnorm)
library(splines)


#### Functions
## The functions to run the aggregate models
source('poisson.lerouxCARaggregatevcm.R')
source('poisson.lerouxCARaggregate.R')

## Helper functions that the above model running functions use
Rcpp::sourceCpp('CARBayesagg.cpp')
source('common.functions.R')
source('print.CARBayes.R')




###########################################
#### Read in and format the areal unit data
###########################################
#### Load the data
load(file="IZ level data.Rdata")
K <- nrow(sp.dat@data)
## Note, the disease incidence data Y have been jittered as the real data are not 
## available for confidentiality reasons. The variables in this data set include:
## IZ - unique intermediate zone code.
## Y - jittered disease incidence data.
## e - expected disesae counts to adjust for population size and demographics.
## no2 / pm10 / pm25 - Areal unit averaged pollution concentrations obtained from the PCM model.
## income / crime / housing / access - SIMD based indicators of socio-economic deprivation.
## propertyprice - average property price used as the effect modifier.


#### Compute the spatial neighbourhood matrix W
W.nb <- poly2nb(sp.dat, row.names = rownames(sp.dat@data))
W <- nb2mat(W.nb, style = "B")



#################################################
#### Read in and format the grid level quantities
#################################################
#### Load the grid level population and pollution data containing the following columns
## bothcoords - unique grid square identifier.
## pop - Estimated population total.
## easting  / northing - grid square coordinates.
## no2 / pm10 / pm25 - grid square pollution concentrations obtained from the PCM model.
load(file="pollution and population grid data.Rdata")


#### Change the grid level data into a spatialpolygons data frame
dat.grid$easting <- as.numeric(dat.grid$easting)
dat.grid$northing <- as.numeric(dat.grid$northing)
sp.grid.pop <- SpatialPixelsDataFrame(points=dat.grid[ ,3:4], data=dat.grid[ ,c(1,2,5:8)])
proj4string(sp.grid.pop) <- proj4string(sp.dat)
sp.grid <- as(sp.grid.pop, "SpatialPolygonsDataFrame")
M <- nrow(sp.grid@data)


#### Read in the area of intersection matrix between each grid square and each areal unit
area <- read.csv(file="Area of intersection.csv")
rownames(area) <- area$X
area$X <- NULL
area <- as.matrix(area)
area <- area[-c(760:762), ]


#### Create the population weighted areas of intersection for use in the aggregate model (8)
area <- area / 1000000
population.weights <- area * matrix(rep(sp.grid@data$pop, K), nrow=K, ncol=M, byrow=TRUE)
weights <- population.weights / matrix(rep(rowSums(population.weights), M), nrow=K, ncol=M, byrow=FALSE)



###############################
#### Run all 4 models with pm10
###############################
#### MCMC quantities
burnin <- 100000
n.sample <- 500000
thin <- 200


#### Ecological model with a constant pollution effect 
form1 <- Y ~ offset(log(e)) + income + crime + housing + access +  pm10
mod1 <- S.CARleroux(formula=form1, data=sp.dat@data, family="poisson", W=W, burnin=burnin, n.sample=n.sample, thin=thin, MALA=FALSE)
print(mod1)


#### Ecological model with a varying pollution effect 
EM <- log(sp.dat@data$propertyprice) - mean(log(sp.dat@data$propertyprice))
vcm <-  EM * sp.dat@data$pm10
form2 <- Y ~ offset(log(e)) + income + crime + housing + access +  pm10 + vcm
mod2 <- S.CARleroux(formula=form2, data=sp.dat@data, family="poisson", W=W, burnin=burnin, n.sample=n.sample, thin=thin, MALA=FALSE)
print(mod2)


#### Aggregate model with a constant pollution effect 
form3 <-  Y ~ offset(log(e)) + income + crime + housing + access
pollutant <- sp.grid@data$pm10
mod3 <- poisson.lerouxCARaggregate(formula=form3, data=sp.dat@data, x.grid=pollutant, weights=weights, W=W, burnin=burnin, n.sample=n.sample, thin=thin)
print(mod3)


#### Aggregate model with a varying pollution effect 
EM <- log(sp.dat@data$propertyprice) - mean(log(sp.dat@data$propertyprice))
form4 <-  Y ~ offset(log(e)) + income + crime + housing + access
pollutant <- sp.grid@data$pm10
mod4 <- poisson.lerouxCARaggregatevcm(formula=form4, data=sp.dat@data,  x.grid=pollutant, EM=EM, weights=weights, W=W, burnin=burnin, n.sample=n.sample, thin=thin,  MALA=FALSE)
print(mod4)


