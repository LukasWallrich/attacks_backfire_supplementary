###### Survey analysis helper functions ######
Int2Factor <- function(x)
{
    if(!is.null(attr(x, "value.labels"))){
        vlab <- attr(x, "value.labels")
        if(sum(duplicated(vlab)) > 0) 
            cat("Duplicated levels:", vlab, "\n")
        else if(sum(duplicated(names(vlab))) > 0) {
            #cat("Duplicated labels:", names(vlab)[duplicated(names(vlab))], "\n")  
          #browser()
            names(vlab) <- 
              paste0(vlab, names(vlab))
          #  cat("Separated labels:", names(vlab), "\n")            
            x <- factor(x, levels = as.numeric(vlab),
                        labels = names(vlab))
           }
        else {
            y <- factor(x, levels = as.numeric(vlab),
                        labels = names(vlab))
            #browser()
            if(sum(is.na(y))==sum(is.na(x))) x <- y
        }
    }
    x
}

#Code missing values
Neg2NA <- function(x) {
  x[as.numeric(x)<0] <- NA
  x
}


## Helper functions for spatial analysis

# Calculates the geodesic distance between two points specified by radian latitude/longitude 
# using the Haversine formula (hf)
gcd.hf <- function(long1, lat1, long2, lat2) {
  long1 %<>% deg2rad
  long2 %<>% deg2rad
  lat1 %<>% deg2rad
  lat2 %<>% deg2rad
  R <- 6371 # Earth mean radius [km]
  delta.long <- (long2 - long1)
  delta.lat <- (lat2 - lat1)
  a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.long/2)^2
  c <- 2 * asin(min(1,sqrt(a)))
  d = R * c
  return(d) # Distance in km
}

deg2rad <- function(deg) return(deg*pi/180)

#Returns logical vector of rows within specific distance of start point
within_km <- function(lon, lat, km, data = arvig) {
  data %>% select(longitude, latitude) %>% {map2(.$longitude, .$latitude, gcd.hf, lon, lat)}  %>% unlist() %>% {.<km}
}

#Summarises geo-dataset based on various parameters
arvig_geo_summary <- function(lon_start, lat_start, km, startD = NULL, endD = NULL, data = arvig, counts = T) {
  #browser()
  if(!(is.null(startD)&is.null(endD))) {
    if(is.null(startD)) startD <- 20140101
    if(is.null(endD)) endD <- 20500101
    startD <- lubridate::ymd(startD)
    endD <- lubridate::ymd(endD)
    data %<>% filter(date >= startD & date <= endD)
  }
  x <- within_km(lon_start, lat_start, km, data)
  if (sum(x)==0) {
    NULL
  } else {
    if (counts) {
      data[within_km(lon_start, lat_start, km, data),] %>% dplyr::count(category_en)
    } else {
      data[within_km(lon_start, lat_start, km, data),]
    }  
  }  
}

#Specific summary functions for assault data
arvig_ags_summary <- function(ags, geoD = geoData, ...) {
  x <- arvig_geo_summary(geoD$lon[geoD$ags == ags], geoD$lat[geoD$ags == ags], ...)
  x$ags <- ags
  x
}
arvig_ags_df <- function(ags_list, ...) {
  x <- map_df(ags_list, arvig_ags_summary, ...)
  spread(x, category_en, n, fill=0)
}

#Specific summary functions for newspaper data
newspaper_ags_df <- function(ags_list, ...) {
  map_df(ags_list, arvig_ags_summary, ..., data = newspapers, counts = FALSE)
}

