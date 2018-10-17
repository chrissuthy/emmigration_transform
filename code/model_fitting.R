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

# first you need to load these packages:
library(devtools)
library(mapview)
library(raster)
library(jagsUI)
library(nimble)
library(RCurl)
library(coda)



# now load the data:
led <- getURL("https://raw.githubusercontent.com/chrissuthy/emmigration_transform/master/data/led.csv")
led <- read.csv(text = led)
metrics <- getURL("https://raw.githubusercontent.com/chrissuthy/emmigration_transform/master/data/patch_metrics.csv")
metrics <- read.csv(text = metrics)
occ.dat <- merge(led, metrics, by.x = "Patch", by.y = "Mpatch")


# we can plot the data too (look in the 'viewer' tab to see the map):
bng <- '+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +datum=OSGB36 +units=m +no_defs'
sp.all <- SpatialPointsDataFrame(coords = metrics[,c("XCoord","YCoord")], data = metrics, proj4string = crs(bng))
mapview(sp.all, cex=sqrt(sp.all@data$Len_m/10), col.regions=ifelse(sp.all@data$Area == "led", "green","red"), 
        map.types = "OpenTopoMap", label = sp.all@data$Mpatch)



#now we need to load the R scripts:
bugs <- getURL("https://raw.githubusercontent.com/chrissuthy/emmigration_transform/master/code/spatial_colext.R")
eval(parse(text = bugs), envir=.GlobalEnv)


# create the occupancy data:
z <- as.matrix(sp.led@data[,paste0("X",2010:2016)])


#CHANGE THIS TO WHATEVER YOU WANT:
P <- matrix(sp.led@data$Len_m/1000,nrow=nrow(sp.led@data),ncol=ncol(z),byrow=FALSE)


nyear <- ncol(z)
nsite <- nrow(z)
dmat <- as.matrix(dist(sp.led@data[,c("XCoord","YCoord")]))
data <- list(z = z, P = P, dmat = dmat)
consts <- list(nyear = nyear, nsite = nsite)
inits <- function() list(delta = runif(1,-4,-2), gamma = runif(1,2,4),
                         beta0 = runif(1,1,3), beta1 = runif(1,-3,-2))

params <- c("delta", "gamma", "beta0", "beta1")

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

