#######################################################
# This is a workflow to compare models that use 
#  various functional transformation of a metric (P) 
#  used to approximate patch specific emmigration.
#
# P must be a matrix with n.site rows and n.years 
#  columns. It cant have missing values or NAs. 
#
# n.site: the number of sites
# n.years: the number of years

############################################################################
# 1. SETUP and SOURCE
############################################################################

# first you need to load these packages:
library(devtools)
library(mapview)
library(raster)
library(jagsUI)
library(nimble)
library(RCurl)
library(coda)

#need this bit of code to source function from git
source_https <- function(url, ...) {
  sapply(c(url, ...), function(u) {
    eval(parse(text = getURL(u, followlocation = TRUE, 
                             cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl"))), 
         envir = .GlobalEnv)
})}
#the code for the model
source_https("https://raw.githubusercontent.com/chrissuthy/emmigration_transform/master/code/spatial_colext.R")




############################################################################
# 2. Loading the data 
############################################################################

# now load the data:
occ <- getURL("https://raw.githubusercontent.com/chrissuthy/emmigration_transform/master/data/occ.csv")
occ <- read.csv(text = occ)

metrics <- getURL("https://raw.githubusercontent.com/chrissuthy/emmigration_transform/master/data/patch_metrics.csv")
metrics <- read.csv(text = metrics)

occ.dat <- merge(occ, metrics, by.x = "Mpatch", by.y = "Mpatch")

bng <- '+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +datum=OSGB36 +units=m +no_defs'
sp.all <- SpatialPointsDataFrame(coords = occ.dat[,c("XCoord","YCoord")], data = occ.dat, proj4string = crs(bng))
mapview(sp.all, cex=sqrt(apply(sp.all@data[,paste0("y",1999:2012)],1,mean,na.rm=T)*100), col.regions=as.numeric(sp.all@data$Area), 
        map.types = "OpenTopoMap", label = sp.all@data$Mpatch)




############################################################################
# 3. Creating data objects
#    z in an npatch x nyear matrix
#    P is an npatch x nyear matrix  
############################################################################

# z: the occupancy states (NAs will be 'estimated')
z <- as.matrix(sp.led@data[,paste0("y",1999:2012)])

# P: what ever the 'surrogate' for effective dispersers:

#for example: area
P1 <- matrix(sp.led@data$Len_m/1000,nrow=nrow(sp.led@data),ncol=ncol(z),byrow=FALSE)

#for example: area^2
P2 <- matrix(sp.led@data$Len_m/1000,nrow=nrow(sp.led@data),ncol=ncol(z),byrow=FALSE)^2

#for example: sqrt(area)
P3 <- matrix(sp.led@data$Len_m/1000,nrow=nrow(sp.led@data),ncol=ncol(z),byrow=FALSE)^0.5

#choose the one I will use:
P <- P2


#need to define the dimensions of the data:
nyear <- ncol(z)
nsite <- nrow(z)

#need to provide the distance matrix
dmat <- as.matrix(dist(sp.led@data[,c("XCoord","YCoord")]))


############################################################################
# 4. package up the data objects and model settings 
############################################################################


data <- list(z = z, P = P, dmat = dmat)           #data object
consts <- list(nyear = nyear, nsite = nsite)      #constants
params <- c("delta", "beta0", "beta1")   #parameters we want to keep track of

#a function to generate initial values
inits <- function(){ 
  list(delta = runif(1,-4,-2), 
       beta0 = runif(1,1,3), 
       beta1 = runif(1,-3,-2))}





############################################################################
# 5. fit the model 
############################################################################

sample <- nimbleMCMC(code = sp.colext,
                     constants = consts,
                     data = data,
                     inits = inits,
                     monitors = params,
                     niter = 15000,
                     nburnin = 5000,
                     thin = 5, 
                     summary = TRUE, 
                     WAIC = TRUE)

coda.samples <- as.mcmc(sample$samples)
par(mfrow=c(2,2))
traceplot(coda.samples[,'delta'] )
traceplot(coda.samples[,'gamma'])
traceplot(coda.samples[,'beta0'])
traceplot(coda.samples[,'beta1'])
sample$WAIC

