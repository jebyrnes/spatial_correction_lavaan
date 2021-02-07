######################################################################
##### Code for Calculating a Corrected SE using Moran's I in SEM
##### covering Composite Variables & Advanced Topics
#####
##### Jarrett E.K. Byrnes
#####
##### Last Modified 1/13/15
######################################################################



## @knitr load-data
# Boreality data from http://www.highstat.com/book2.htm
# Mixed Effects Models and Extensions in Ecology with R (2009). 
# Zuur, Ieno, Walker, Saveliev and Smith. Springer
boreal <- read.table("./Boreality.txt", header=T)

#For later
source("./lavSpatialCorrect.R")
source("./predict_lavaan.R") #prediction by regression functions

## @knitr visualize-data
#Let's look at the spatial structure
library(ggplot2)

qplot(x, y, data=boreal, size=Wet, color=NDVI) +
  theme_bw(base_size=18) + 
  scale_size_continuous("Index of Wetness", range=c(0,10)) + 
  scale_color_gradient("NDVI", low="lightgreen", high="darkgreen")


## @knitr sem-model
library(lavaan)

# A simple model where NDVI is determined
# by nTot, temperature, and Wetness
# and nTot is related to temperature
borModel <- '
  NDVI ~ nTot + T61 + Wet 
  nTot ~ T61
'

#note meanstructure=T to obtain intercepts
borFit <- sem(borModel, data=boreal, meanstructure=T)

## @knitr residuals
# residuals are key for the analysis
borRes <- as.data.frame(residuals_lavaan(borFit))
#raw visualization of NDVI residuals
qplot(x, y, data=boreal, color=borRes$NDVI, size=I(5)) +
  theme_bw(base_size=17) + 
  scale_color_gradient("NDVI Residual", low="blue", high="yellow")

## @knitr residual-analysis-sign
#raw visualization of sign of residuals
qplot(x, y, data=boreal, color=borRes$NDVI>0, size=I(5)) +
  theme_bw(base_size=17) + 
  scale_color_manual("NDVI Residual >0", values=c("blue", "red"))



## @knitr generate-spatial-weight-matrix
#Evaluate Spatial Residuals
#First create a distance matrix
library(ape)
distMat <- as.matrix(dist(cbind(boreal$x, boreal$y)))

#invert this matrix for weights
distsInv <- 1/distMat
diag(distsInv) <- 0

## @knitr moran_i_ndvi
#calculate Moran's I just for NDVI
mi.ndvi <- Moran.I(c(borRes$NDVI), as.matrix(distsInv))
mi.ndvi

## @knitr spatial_corrected_sample
#What is our corrected sample size?
n.ndvi <- nrow(boreal)*(1-mi.ndvi$observed)/(1+mi.ndvi$observed)


## @knitr spatial_corrected_se
#Where did we get the SE from?
sqrt(diag(vcov(borFit)))

#New SE
ndvi.var <- diag(vcov(borFit))[1:3]

ndvi.se <- sqrt(ndvi.var*nrow(boreal)/n.ndvi)

ndvi.se


#compare to old SE
sqrt(diag(vcov(borFit)))[1:3]


## @knitr spatial_corrected_z
#new z values
z <- coef(borFit)[1:3]/ndvi.se

2*pnorm(abs(z), lower.tail=F)

#compare with lavSpatialCorrect function
spat_cor<-lavSpatialCorrect(borFit, boreal$x, boreal$y)
spat_cor$parameters$NDVI[,c(5,6)]
