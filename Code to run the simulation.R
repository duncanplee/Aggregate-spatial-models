############################################################################
#### Simulation study to see what effect MAUP has on regression coefficients
############################################################################

###################
#### Load libraries
###################
library(spdep)
library(sp)
library(MASS)
library(spam)
library(CARBayes)
library(coda)
library(truncnorm)

source('poisson.lerouxCARaggregate.R')
source('common.functions.R')
Rcpp::sourceCpp('CARBayesagg.cpp')
source('print.CARBayes.R')

set.seed(1)

####################
#### Set up the grid
####################
#### Set up a square lattice region
x.easting <- 1:71
x.northing <- 1:71
Grid <- as.matrix(expand.grid(x.easting, x.northing))
K.grid <- nrow(Grid)

#### Create the W matrix
W.grid.nb.temp <- knn2nb(knearneigh(Grid, k=4), row.names=1:K.grid)
W.grid.nb <- make.sym.nb(W.grid.nb.temp)
W.grid <- nb2mat(W.grid.nb, style = "B")
W.grid.list <- nb2listw(W.grid.nb)


#### Create the distance matrix for all grid squares
D <- as.matrix(dist(Grid))



##########################
#### Simulation quantities
##########################
#### Create the quantities of the simulation
## Amount of spatial correlation
nu1 <- 0.5 # Covariate x
nu2 <- 0.1 # Random effects phi

## Population size P()
pop.sd <- 10
pop.min <- 30

## Disease prevalence
e.prop <- 0.05 # 0.01 or 0.05

## Covariate standard deviation
x.sd <- 1 # 1 or 5

## Pollution - disease effect size
beta <- 0.1

## Random effects standard deviation (phi)
phi.sd <- 0.2 # 0.04 or 0.2


## Number of areal units K
K <- 100 # 100 or 500 or 1000


#### Compute the spatial variance matrices and their cholesky factors
Sigma1 <- exp(-nu1 * D)
Sigma2 <- exp(-nu2 * D)
chol.Sigma1 <-t(chol(Sigma1))
chol.Sigma2 <-t(chol(Sigma2))


#### Create matrices to save the results in
n.dat <- 50
n.replicates <- 20
n.sim <- n.dat * n.replicates
beta.results <- array(NA, c(n.sim,9))
colnames(beta.results) <- c("x_m1_est", "x_m2_est", "x_m1_LCI", "x_m2_LCI", "x_m1_UCI", "x_m2_UCI", "dataset", "within_var_x1","between_var_x1")
beta.results[ ,7] <- kronecker(1:n.dat, rep(1, n.replicates))

#### The two models are:
## M1 - Naive ecological model
## M2 - Proposed aggregate model


#### MCMC quantities
burnin <- 100000
n.sample <- 300000
thin <- 100


#############################
#### Run the simulation study
#############################
counter <- 1

for(z in 1:n.dat)
{
    ####################################
    #### Generate grid level "true" data
    ####################################
    #### Population and offset
    pop.temp <- pop.sd * chol.Sigma2 %*% rnorm(K.grid) + pop.min
    pop.grid <- pop.temp^2 
    e.grid <- pop.grid * e.prop

    #### Covariates
    ## Use the first two lines for an independent x covariate and the 
    ## last two for a spatially correlated x covariate
    #x.grid <- rnorm(K.grid)
    #x.grid <- x.sd * (x.grid - mean(x.grid)) / sd(x.grid)
    x.grid <- chol.Sigma1 %*% rnorm(K.grid)
    x.grid <- x.sd * (x.grid - mean(x.grid)) / sd(x.grid)
    
    #### Random effects
    phi.grid <- chol.Sigma2 %*% rnorm(K.grid)
    phi.grid <- phi.sd * (phi.grid - mean(phi.grid)) / sd(phi.grid)

    #### Response
    lp.grid <- beta * x.grid  + phi.grid
    lp.grid <- lp.grid - mean(lp.grid)
    risk.grid <- exp(lp.grid)
    fitted.grid <- e.grid * risk.grid
    Y.grid <- rpois(n=K.grid, lambda=fitted.grid)
    sp.grid <- SpatialPixelsDataFrame(points=Grid, data=data.frame(region=1:K.grid, x.grid, phi.grid, pop.grid, e.grid, risk.grid, fitted.grid, Y.grid))
    #spplot(sp.grid, "pop.grid")

    
    for(v in 1:n.replicates)
    {
        ###########################################
        #### Aggregate the grid data to areal units
        ###########################################
        #### Choose the areal unit centroids
        area.centroids <- Grid[sample(x=1:K.grid, size=K, replace=FALSE), ]
        
        #### Create the distance matrix between each centroid and each grid square
        Dist.mat <- array(NA, c(K, K.grid))
            for(i in 1:K)
            {
            centroid.temp <- area.centroids[i, ]
            Dist.mat[i, ] <- as.numeric(sqrt((centroid.temp[1] - Grid[ ,1])^2 + (centroid.temp[2] - Grid[ ,2])^2))
            }
        
        #### Assign each grid square to its nearest centroid
        ordering <- apply(Dist.mat, 2, order)
        area <- ordering[1, ]
        sp.grid@data$area <- area
        spplot(sp.grid, "area")
        
        #### Create the area of intersection matrix and hence the weights
        area.mat <- array(0, c(K, K.grid))
            for(r in 1:K)
            {
            temp <- which(sp.grid@data$area==r)
            area.mat[r, temp] <- 10000
            }
        
        
        #### Create the aggregted data set
        Y.area <- as.numeric(tapply(sp.grid@data$Y.grid, area, sum))
        e.area <- as.numeric(tapply(sp.grid@data$e.grid, area, sum))
        x.area <- as.numeric(tapply(sp.grid@data$x.grid, area, mean)) 
        dat.area <- data.frame(region=1:K, Y.area, e.area, x.area)
        
        #### Create the W matrix
        W.area.nb.temp <- knn2nb(knearneigh(area.centroids, k=4), row.names=1:K)
        W.area.nb <- make.sym.nb(W.area.nb.temp)
        W.area <- nb2mat(W.area.nb, style = "B")
        W.area.list <- nb2listw(W.area.nb)

        
        #### Summarise the within and between variation in both x variables
        beta.results[counter, 8] <- sd(x.area)
        beta.results[counter, 9] <- mean(tapply(sp.grid@data$x.grid, area, sd), na.rm=TRUE)

        
                
        #################
        #### Fit model M2
        #################
        formula.M1 <- Y.area ~ offset(log(e.area)) + x.area
        M1 <- S.CARleroux(formula=formula.M1, family="poisson", W=W.area, burnin=burnin, n.sample=n.sample, thin=thin)
        beta.results[counter, 1] <- M1$summary.results[2, 1]
        beta.results[counter, 3] <- M1$summary.results[2, 2]
        beta.results[counter, 5] <- M1$summary.results[2, 3]
        

        
        
        #################    
        #### Fit model M3
        #################
        area.mat <- area.mat / 1000000
        population.weights <- area.mat * matrix(rep(sp.grid@data$pop.grid, K), nrow=K, ncol=K.grid, byrow=TRUE)
        weights <- population.weights / matrix(rep(rowSums(population.weights), K.grid), nrow=K, ncol=K.grid, byrow=FALSE)

        formula.M2 <- Y.area ~ offset(log(e.area))
        x.grid <- x.grid
        M2 <- poisson.lerouxCARaggregate(formula=formula.M2, x.grid=x.grid, weights=weights, W=W.area, burnin=burnin, n.sample=n.sample, thin=thin, verbose=TRUE)
        beta.results[counter, 2] <- M2$summary.results[2, 1]
        beta.results[counter, 4] <- M2$summary.results[2, 2]
        beta.results[counter, 6] <- M2$summary.results[2, 3]
        
        
        
        ###############################
        #### Save the results to a file
        ###############################
        counter <- counter + 1
        write.csv(beta.results, file="simulation results.csv")
    }
}

