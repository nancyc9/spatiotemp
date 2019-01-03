
#=====================================#
#	Eliminates the observations in 'input.covars' that
#	fall outside the domain specified by
#	'regionPoly'
#=====================================#
screen.outOfRegion <- function(input.covars, regionPoly, checkStates = c("CA", "IL", "MD", "MN" ,"NY", "NC"), addRegionTag=TRUE){
	sp.installed <- require("sp")
	if(!sp.installed) stop("Package 'sp' is required")
	if (any(!checkStates %in% names(regionPoly))) stop("Checking state that has no polygon region")

	input.covars$outOfMapFlag = TRUE
	input.covars$region = "None"
	for (r in 1:length(checkStates)){
		point.x <- input.covars$longitude
		point.y <- input.covars$latitude
		mp <- regionPoly[[checkStates[r]]] # mp = modelpoly
		cutoff <- which(is.na(mp$x) & is.na(mp$y))
		if (length(cutoff)>0){
			poly.list <-vector(mode="list", length=length(cutoff) + 1)
			ind1 <- 1
			ind2 <- cutoff[1] - 1
			for (i in 1:(length(cutoff)+1)){
				if( i ==3) cat("Three polys not supported\n")
				poly.list[[i]]$x <- mp$x[ind1:ind2]
				poly.list[[i]]$y <- mp$y[ind1:ind2]
				ind1 <- cutoff[i]+1
				ind2 <- min(cutoff[i+1]-1, length(mp$x), na.rm=T)
			}
		} else {
			poly.list <- list(mp)
		}
		in.region <- vector(mode="list",length=length(poly.list))
		for (j in 1:length(poly.list)){
			in.region[[j]] <- point.in.polygon(point.x, point.y,
									 poly.list[[j]]$x, poly.list[[j]]$y)
		}
		in.region <- apply(matrix(unlist(in.region), ncol=length(in.region)), 1, any)
        input.covars$region[as.logical(in.region)] <- checkStates[r]
        input.covars$outOfMapFlag <- input.covars$outOfMapFlag & !in.region
	}
	output.covars <- input.covars[input.covars$outOfMapFlag==0,]
	output.covars$outOfMapFlag <- NULL
	return(output.covars)
}




# Function for computing a lambert projection of locations,
# given latitude and longitude.  Author is believed to be P. Sampson,
# but this documentation is written J. Keller.
#
# geoconfig should be a matrix/dataframe of two columns:
# The first is Latitude, the second is -Longitude
#
######
##		NOTE!!!  The second column of 'geoconfig' should be -Longtiude,
##               which is -1*Longitude
#######
#
#	If computing a new lambert project, leave the other arguments blank.  If you
# want the projections for new locations given a previous projection, insert
# all arguments using the appropriate parameters.
#
Flamb2<-function(geoconfig, latrf1 = NA, latrf2 = NA, latref = NA, lngref = NA)
{
      geo <- as.matrix(geoconfig)
      if(dim(geo)[[2.]] != 2.)
              stop("the input should be an nx2 matrix")
      xy.coord <- array(0., dim(geo))
      if(is.na(latref))
              latref <- sum(range(geo[, 1.], na.rm = T))/2.
      if(is.na(lngref))
              lngref <- sum(range(geo[, 2.], na.rm = T))/2.
      if(is.na(latrf1))
              latrf1 <- latref - (3./10.) * diff(range(geo[, 1.], na.rm = T))
      if(is.na(latrf2))
              latrf2 <- latref + (3./10.) * diff(range(geo[, 1.], na.rm = T))
      lat <- geo[, 1.]
      long <- geo[, 2.]
      pi <- 3.14159265
      a <- 6378137.
      b <- 6356752.
      radlf1 <- (pi/180.) * latrf1
      radlf2 <- (pi/180.) * latrf2
      radlgf <-  - (pi/180.) * lngref
      radltf <- (pi/180.) * latref
      eccen <- sqrt((a^2. - b^2.)/a^2.)
      capr <- (a * (1. - eccen^2.))/((1. - eccen^2. * sin(radltf)^2.)^1.5)
      n <- log(cos(radlf1)/cos(radlf2))/(log(tan(pi/4. + radlf2/2.)/tan(pi/4. + radlf1/2.)))
      capf <- (cos(radlf1) * ((tan(pi/4. + radlf1/2.))^n))/n
      #
      rho0 <- (capr * capf)/#
      ((tan(pi/4. + radltf/2.))^n)
      radlat <- (pi/180.) * lat
      radlng <-  - (pi/180.) * long
      theta <- n * (radlng - radlgf)
      rho <- (capr * capf)/((tan(pi/4. + radlat/2.))^n)
      x <- (0.001) * rho * sin(theta)
      y <- (0.001) * (rho0 - rho * cos(theta))
      return(list(xy = cbind(x, y),
                  latrf1 = latrf1,
                  latrf2 = latrf2,
                  latref = latref,
                  lngref = lngref))
}



### Standard function for calculating canyon depth from height and at-floor data:
canyonHeight2depth=function(wallheight,atfloor,width,maxval=10,freeheight=6,freewidth=100)
{
depth=wallheight-atfloor
depth[depth>0 & width>=freewidth]=0 # If there's no canyon to speak of (too wide), you're not inside one
depth[depth>0 & atfloor>=freeheight]=0 # If you're too high above ground, you're not in a canyon trafficwise
depth[depth>maxval]=maxval # upper cutoff...
depth[depth<(-maxval)]=-maxval #... and lower

return(depth)
}


canyon_process <- function(covars, includeCanyonDepth=TRUE, maxwidth=120, minwidth=20, roadmin=25){

x <- covars

canyon_cols <- grep("^canyon",names(x))


# Lower truncate all values at 0.5
# This includes height, width, # of floors
for (col in canyon_cols){
	x[,col] <- ifelse(x[,col]<0.5, 0.5, x[,col])
}

# Put in 0.5 for locations missing building height
x$canyon_bldg_hgt[is.na(x$canyon_bldg_hgt)] <- 0.5

# For sites missing canyon width, impute as maximum observed canyon width
x$canyon_width[is.na(x$canyon_width)] <- maxwidth
# Lower truncate canyon width at 20
x$canyon_width[x$canyon_width<minwidth] <- minwidth
x$canyon_width[x$canyon_width> maxwidth] <- maxwidth

x$canyon_bldg_meanflrs[is.na(x$canyon_bldg_meanflrs)] <- 0.5
x$canyon_opp_meanflrs[is.na(x$canyon_opp_meanflrs)] <- 0.5

buffers <- x[,grep("a00",names(x))]
buffers[is.na(buffers)] <- 0.5
x[,grep("a00",names(x))] <- buffers


# Add in minimum road length variable
x$canyon_road_len[x$canyon_road_len<roadmin] <- roadmin
x$canyon_road_len[is.na(x$canyon_road_len)]<- roadmin

if (includeCanyonDepth){
	x$canyon_depth_a00050 <- canyonHeight2depth(x$canyon_height_a00050, covars$living_floor, width=x$canyon_width)
	x$canyon_depth_a00100 <- canyonHeight2depth(x$canyon_height_a00100, covars$living_floor, width=x$canyon_width)
	x$canyon_depth_a00150 <- canyonHeight2depth(x$canyon_height_a00150, covars$living_floor, width=x$canyon_width)
	x$canyon_depth_a00300 <- canyonHeight2depth(x$canyon_height_a00300, covars$living_floor, width=x$canyon_width)

}

return(x)
}
