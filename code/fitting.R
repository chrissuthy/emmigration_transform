#analyze the Assynt data:

library(mapview)
library(raster)
library(jagsUI)
library(nimble)
library(coda)


occ <- read.csv("led.csv")
metrics <- read.table("patchmetrics.txt", header=T)
led <- merge(occ, metrics, by.x = "Patch", by.y = "Mpatch")

bng <- '+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +datum=OSGB36 +units=m +no_defs'
sp.led <- SpatialPointsDataFrame(coords = led[,c("XCoord","YCoord")], data = led, proj4string = crs(bng))
sp.all <- SpatialPointsDataFrame(coords = metrics[,c("XCoord","YCoord")], data = metrics, proj4string = crs(bng))
#mapview(sp.led, 
#        cex=apply(sp.led@data[,paste0("X",2010:2016)],1,sum)*2, 
#        col.regions="green")
mapview(sp.all, 
        cex=sqrt(sp.all@data$Len_m/10), 
        col.regions=ifelse(sp.all@data$Area == "led", "green","red"),
        map.types = "OpenTopoMap", label = sp.all@data$Mpatch)

source("spatial_colext.R")



z <- as.matrix(sp.led@data[,paste0("X",2010:2016)])

#CHANGE THIS TO WHATEVER YOU WANT:
P <- matrix(sp.led@data$Len_m/1000,nrow=nrow(sp.led@data),ncol=ncol(z),byrow=FALSE)
#

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

par(mfrow=c(1,2))
aseq <- seq(min(data$P),max(data$P),length=100)
plot(aseq,aseq,ylim=c(0,1), axes=F, xlab="",ylab="",type="n")
box(bty="o");axis(1);axis(2)
for(i in 1:nrow(coda.samples)){
  lines(aseq, plogis(coda.samples[i,"beta0"] + aseq*coda.samples[i,"beta1"]),
        col=adjustcolor(4,0.01))
}
lines(aseq, plogis(mean(coda.samples[,"beta0"]) + aseq*mean(coda.samples[,"beta1"])),
      col=4, lwd=3)

plot(aseq,aseq, ylim=c(0,8), xlim=c(0,2), axes=F, xlab="",ylab="",type="n")
box(bty="o");axis(1);axis(2)
for(i in 1:nrow(coda.samples)){
  lines(aseq, aseq^coda.samples[i,"gamma"],col=adjustcolor(2,0.005))
}
lines(aseq, aseq^mean(coda.samples[i,"gamma"]),
      col=2, lwd=3)
abline(0,1,lwd=2,lty=2)
