require('sp')

###########################################################
# Approximate endpoints of line                           #
# Wells&Coppersmith 1994 FLR/FWR(Surface Rupture)         # 
###########################################################
getEnds <- function(lat,lon,mag,strike){
  # Estimating surface rupture length and width
  rl <- -2.42 + 0.58 * mag
  rw <- -1.61 + 0.41 * mag
  
  flt <- c(10^rl,10^rw)
  
  lat1 <- lat
  lon1 <- lon
  len <- flt[1]/2
  theta <- strike*pi/180
  for (d in seq(0.05,4,0.05)){
    
    lat2 <- lat1+d*sin(theta)
    lon2 <- lon1+d*cos(theta)
    if (exists("est")) est <- rbind(est,d) else est <- d
    if (exists("dst")){
      tmp <- spDistsN1(as.matrix(cbind(lon1,lat1)),c(lon2,lat2),longlat=TRUE)
      dst <- rbind(dst,tmp)
    }else{
      dst <- spDistsN1(as.matrix(cbind(lon1,lat1)),c(lon2,lat2),longlat=TRUE)
    }
  }
  dsel <- which(dst >= len)[1]
  p1 <- c(lon1+est[dsel]*cos(theta) , lat1+est[dsel]*sin(theta)) 
  p2 <- c(lon1-est[dsel]*cos(theta) , lat1-est[dsel]*sin(theta))
  
  return(rbind(p1,p2))
}

###########################################################
# Gives back endpoints of rupture according to WC1994     #
###########################################################
mag <- 6.38
lon <- 35.5
lat <- 31.8
strike <- 246.9

endpoints <- getEnds(lat,lon,mag,strike)

