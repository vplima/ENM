##ENM for native potential agroforestry species
#Valdeir Lima
#May 27, 2021

##My Library
library(CoordinateCleaner)
library(PresenceAbsence)
library(maptools)
library(spatstat)
library(stats)
library(dismo)
library(raster)
library(terra)
library(rgdal)
library(usdm)
library(rgeos)
library(spThin)
library(kernlab)
library(mgcv)
library(rJava)
library(speciesgeocodeR)
library(ConR)
library(readr)
library(dplyr)
library(rredlist)
library(jsonlite)
library(rnaturalearthdata)

####Download GBIF and Cleaning pipeline
{sp_occ <- gbif(genus="genus", species="epithet", extent(-85, -35, -55, 12))
sp_occ01 <- subset(sp_occ, !is.na(lon) & !is.na(lat))
write.csv(sp_occ01,"~/Desktop/spp.csv")
  
c <- read.csv("~/Desktop/spp.csv", header = TRUE)
c01 <- clean_coordinates(c, lon = "lon", lat = "lat", species = "species",
                           tests = c("capitals","centroids", "equal", "gbif", 
                                 "zeros","institutions", "outliers", "seas"))
write.csv(c01,"~/Desktop/c01.csv")
d <- read.csv("~/Desktop/c01.csv")
d01 <- cc_sea(d, lon = "lon", lat = "lat",
                ref = NULL, scale = 110, value = "clean", speedup = TRUE,
                verbose = TRUE)
write.csv(d01,"~/Desktop/step01.csv")
  
step02 <- read.csv("~/Desktop/step01.csv")
Crs      = CRS("+proj=longlat +datum=WGS84")
remove.spatial.outliers = function (step02, quant = 0.075, plot.map = T){
    
    force(quant)
    
    X = SpatialPoints(step02, Crs)
    Y = ppp(x=step02[,1], y=step02[,2], c(-85, -35),c(-55, 12) )
    
    if(quant == -1){
      q_test <- intersect(X, regions)
      n.regions = length(unique(q_test@data[1][,1]))
      if(n.regions == 5) quant = 0.025
      if(n.regions == 6) quant = 0.010
      #cat("q=", quant, "\n")
    }
    
    D = density(Y, kernel = "gaussian")
    
    d = extract(raster(D), X)
    d[is.na(d)] = 0
    qp = quantile(d, quant)
    
    if(plot.map == T){
      plot(D, main= "density kernel estimation")
      points(X, pch=16, cex=0.5, col="green")
      contour(D, axes = FALSE, add = T)
      #add.geography(draw.forestborder = T)
      points(X[d<qp], pch=16, cex=0.5, col="red")
    }
    
    return(step02[d>=qp,])
  }
  dt = step02[,c("lon","lat")]
  dt_kernel <- remove.spatial.outliers(dt)
  write.csv(dt_kernel, "~/Desktop/spp.csv")}

##Sampling bias and spatial autocorrelation
x <- list.files("~/Desktop/ENM_R/Records")

for (i in x) {
  
  SDM_records <- read.csv(paste0("~/Desktop/ENM_R/Records/",i,""))
  thin(loc.data = SDM_records, 
       lat.col = "lat", long.col = "lon", 
       spec.col = "X", 
       thin.par = 20, reps = 100, 
       locs.thinned.list.return = TRUE, 
       write.files = TRUE, 
       max.files = 1, 
       out.dir = "~/Desktop/Results", out.base = paste0("",i,""), 
       write.log.file = TRUE,
       log.file = paste0("",i,"_thinned_full_log_file.txt" ))
  
}

##Predictors
wc <- stack("~/Desktop/ENM_R/wc_current/bio1.bil",
            "~/Desktop/ENM_R/wc_current/bio2.bil",
            "~/Desktop/ENM_R/wc_current/bio3.bil",
            "~/Desktop/ENM_R/wc_current/bio4.bil", 
            "~/Desktop/ENM_R/wc_current/bio6.bil",
            "~/Desktop/ENM_R/wc_current/bio7.bil",
            "~/Desktop/ENM_R/wc_current/bio12.bil",
            "~/Desktop/ENM_R/wc_current/bio14.bil",
            "~/Desktop/ENM_R/wc_current/bio15.bil",
            "~/Desktop/ENM_R/wc_current/bio19.bil")

wp <- extent(90, -25, -70, 20)
plot(worldclim <- stack(crop(wc, wp)))

##VIF
v1 <- vifstep(worldclim)
v1
str(v1)
class(v1)
v1@results

###Most relevant predictors to different IUCN growth forms
spp <- read.csv("~/Desktop/ENM_R/Records/x.csv", header = TRUE)
spp.sp <- SpatialPoints(cbind(spp$lon, spp$lat))
proj4string(spp.sp) <- "+proj=longlat +datum=WGS84 +no_defs"
DataSpecies <- spp.sp
hull <- SpatialPoints(cbind(spp$lon, spp$lat))
spp_hull <- gConvexHull(hull)
hull_plus_buffer <- gBuffer(spp_hull,width=2.7)
wc <- stack("~/Desktop/ENM_R/wc_current/bio1.bil",
            "~/Desktop/ENM_R/wc_current/bio2.bil",
            "~/Desktop/ENM_R/wc_current/bio3.bil",
            "~/Desktop/ENM_R/wc_current/bio4.bil", 
            "~/Desktop/ENM_R/wc_current/bio6.bil",
            "~/Desktop/ENM_R/wc_current/bio7.bil",
            "~/Desktop/ENM_R/wc_current/bio12.bil",
            "~/Desktop/ENM_R/wc_current/bio14.bil",
            "~/Desktop/ENM_R/wc_current/bio15.bil",
            "~/Desktop/ENM_R/wc_current/bio19.bil")
wp <- mask(wc,hull_plus_buffer)
myExpl <- stack(crop(wp, extent(-90, -25, -70, 20) ))
plot(myExpl)

##VIF
v1 <- vifstep(myExpl)
v1
str(v1)
class(v1)
v1@results

##Recreating the data
predictors <- stack(worldclim[[c("bio1", "bio4", "bio7", "bio13", "bio14", "bio15")]])
plot(predictors)
Araucaria <- read.table("~/Desktop/ENM_R/Records_1/x.csv",  header=TRUE,  sep=',')
Araucaria <- Araucaria[,-1]
head(Araucaria)

pred_nf <- dropLayer(predictors,'')

##Training and testing set
set.seed(0)
group <- kfold(Araucaria, 5)
pres_train <- Araucaria[group != 1, ]
pres_test <- Araucaria[group == 1, ]

##Restricting predictions
ext <- extent(-90, -25, -70, 20) 

#Background data for training and a testing set
set.seed(10)
backg <- randomPoints(pred_nf, n=10000, ext=ext, extf = 1.25)
colnames(backg) = c('lon', 'lat')
group <- kfold(backg, 5)
backg_train <- backg[group != 1, ]
backg_test <- backg[group == 1, ]

#Adding predictors
r <- raster(predictors, 1)
plot(!is.na(r), col=c('white', 'light grey'), legend=FALSE)
plot(ext, add=TRUE, col='red', lwd=2)
points(backg_train, pch='-', cex=0.5, col='yellow')
points(backg_test, pch='-',  cex=0.5, col='black')
points(pres_train, pch= '+', col='green')
points(pres_test, pch='+', col='blue')

#Models
maxent()
xm <- maxent(predictors, pres_train)
plot(xm)
response(xm)

e <- evaluate(pres_test, backg_test, xm, predictors)
e

px <- predict(predictors, xm, ext=ext, progress='')
par(mfrow=c(1,2))
plot(px, main='Maxent')
plot(wrld_simpl, add=TRUE, border='dark grey')
tr <- threshold(e, 'spec_sens')
plot(px > tr, main='Binary')
plot(wrld_simpl, add=TRUE, border='dark grey')
points(pres_train, pch='+')

x <- extract(predictors, Araucaria)
nr <- nullRandom(x, bioclim, n=25, rep=25, pa=FALSE)
mean(sapply(nr, function(x)x@auc))

##Conservation assessment

dat <- read_csv("spp_data.csv")%>%
  dplyr::select(X, #species by X
                decimallongitude = lon,
                decimallatitude = lat)

inp <- dat%>%
  dplyr::select(ddlat = decimallatitude,
                ddlon = decimallongitude, 
                tax = X) #species by X
ev <- IUCN.eval(inp)
ev


