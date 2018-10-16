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
library(mapview)
library(raster)
library(jagsUI)
library(nimble)
library(coda)



# now load the data:
read.table(url("https://www.dropbox.com/sh/rt630kpps8mc7fb/AABwYNLHp1Z2MzfOIblvdArOa?dl=0"))
read.table(url("https://www.dropbox.com/s/fb5pdqfosgh7xgj/patch_metrics.csv?dl=0

