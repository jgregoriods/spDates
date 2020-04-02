library(automap)
library(data.table)
library(dplyr)
library(gdistance)
library(ggplot2)
library(gstat)
library(maptools)
library(parallel)
library(raster)
library(rcarbon)
library(smatr)
library(viridisLite)


#' Filter archaeological site coordinates and dates, retaining only the
#' earliest radiocarbon date per site.
#'
#' @param sites A SpatialPointsDataFrame object with archaeological sites and
#' associated radiocarbon ages.
#' @param c14bp A string. Name of the field with the radiocarbon ages in C14 BP
#' format.
#' @return A SpatialPointsDataFrame object with the earliest C14 date for every
#' site.
#' @export
filterDates <- function(sites, c14bp) {
    x <- c(colnames(coordinates(sites)))[1]
    y <- c(colnames(coordinates(sites)))[2]

    clusters <- zerodist(sites, zero = 0, unique.ID = TRUE)

    sites$clusterID <- clusters
    sites.df <- as.data.frame(sites)
    sites.max <- as.data.frame(sites.df %>% group_by(clusterID) %>%
                               top_n(1, get(c14bp)))

    xy <- cbind(sites.max[[x]], sites.max[[y]])

    sites.max.spdf <- SpatialPointsDataFrame(xy, sites.max)
    proj4string(sites.max.spdf) <- proj4string(sites)

    return(sites.max.spdf)
}


#' Perform regression of dates versus distances from multiple potential origins
#' in order to find the best model. It is also possible to specify multiple
#' distances for the spatial binning of the dates.
#'
#' @param ftrSites A SpatialPointsDataFrame object with associated earliest
#' C14 dates per site and respective calibrated distributions (CalDates
#' objects) in a field named "cal".
#' @param c14bp A string. Name of the field with the radiocarbon ages in C14
#' BP format.
#' @param origins A SpatialPointsDataFrame object. The sites to be tested for
#' the most likely origin of expansion.
#' @param siteNames A string. Name of the field with the site names or ids for
#' the potential origins.
#' @param binWidths A number or vector of numbers. Width(s) of the spatial bins
#' in km.
#' @param nsim A number. Number of simulations to be run during the
#' bootstrapping procedure. Default is 999.
#' @param cost Optional. A RasterLayer with friction values to calculate
#' least cost path distances instead of great circle distances.
#' @param method A string. Method to be used in the regression. One of "rma"
#' or "ols". Default is "rma".
#' @return a list with two elements, the result of the iteration over all
#' potential origins and the best model selected among those.
#' @export
iterateSites <- function(ftrSites, c14bp, origins, siteNames, binWidths = 0,
                         nsim = 999, cost = NULL, method = "rma") {

    datalen <- length(binWidths) * length(origins)

    res <- data.table("r" = numeric(datalen), "p" = numeric(datalen),
                      "bin" = numeric(datalen), "n" = numeric(datalen),
		              "site" = character(datalen))

    origins$r <- 0

    # Perform reduced major axis regression for each bin width and for each
    # hypothetical origin, if more than one is specified
    pb <- txtProgressBar(min = 0, max = datalen, style = 3)
    counter <- 0
    for (i in (1:length(binWidths))) {
        for (j in (1:length(origins))) {
            counter <- counter + 1

            dateModel <- modelDates(ftrSites, c14bp = c14bp,
                                    origin = origins[j,],
                                    binWidth = binWidths[i], nsim = nsim,
                                    cost = cost, method = method)

            if (method == "rma") {
                slope <- mean(dateModel$model[,"slo"])
            } else if (method == "ols") {
                # For practical reasons, consider only time-versus-distance
                # for ols models
                slope <- mean(dateModel$td.model[,"slo"])
            }

            # Only proceed if the slope is negative, i.e. dates become more
            # recent with increasing distance
            if (slope < 0) {

                if (method == "rma") {
                    meanr <- mean(dateModel$model[,"r"])
                    meanp <- mean(dateModel$model[,"p"])
                } else if (method == "ols") {
                    # For practical reasons, consider only time-versus-distance
                    # for ols models
                    meanr <- mean(dateModel$td.model[,"r"])
                    meanp <- mean(dateModel$td.model[,"p"])
                }

                # Save the best model
                if (!exists("bestModel")) {
                    bestModel <- dateModel
                    bestModel$meanr <- meanr
                } else if (meanp < 0.05 & meanr > bestModel$meanr) {
                    bestModel <- dateModel
                    bestModel$meanr <-meanr
                }

                if (meanp < 0.01) {
                    ast <- "**"
                } else if (meanp < 0.05) {
                    ast <- "* "
                } else {
                    ast <- "  "
                }

                if (binWidths[i] > 0) {
                    nsites <- nrow(dateModel$binSites)
                } else {
                    nsites <- nrow(dateModel$allSites)
                }

                # Extract the r and p of each model and store it in the result,
                # together with the corresponding bin width
                rval <- round(as.double(sqrt(meanr)), 4)
                res[counter, "r"] <- rval
                res[counter, "p"] <- round(as.double(meanp), 2)
                res[counter, "bin"] <- binWidths[i]
                res[counter, "n"] <- nsites
                res[counter, "site"] <- paste(toString(origins[[siteNames]][j]),
                                              ast, sep="")
                
                if (rval > origins$r[j]) {
                    origins$r[j] <- rval
                }
            }
            setTxtProgressBar(pb, counter)
        }
    }
    close(pb)

    # Reproject to Robinson for Kriging
    #origins.m <- spTransform(origins, CRS("+proj=robin +lon_0=0 +x_0=0 +y_0=0
    #                                       +ellps=WGS84 +datum=WGS84 +units=m
    #                                       +no_defs"))
    
    # Perform autoKriging
    #krg.m <- autoKrige(r~1, origins.m)
    
    # Reproject back to lonlat
    #newr <- raster(extent(origins))
    #proj4string(newr) <- proj4string(origins)
    #res(newr) <- 0.1
    #krg <- projectRaster(raster(krg.m$krige), newr)

    origins.idw <- interpolateIDW(origins, "r")

    # Create map object
    map <- list("sites" = ftrSites, "origins" = origins, "idw" = origins.idw)
    class(map) <- "dateMap"

    res <- res[order(-r)]
    return(list("results" = res, "model" = bestModel, "map" = map))
}


#' Interpolate using inverse distance weighting.
#' @param points A SpatialPointsDataFrame object.
#' @param attr A string. Name of the field with values to interpolate.
#' @return A RasterLayer.
interpolateIDW <- function(points, attr) {
    grd <- as.data.frame(spsample(points, "regular", n = 10000))
    names(grd) <- c("x", "y")
    coordinates(grd) <- ~x+y
    proj4string(grd) <- proj4string(points)
    gridded(grd) <- TRUE
    fullgrid(grd) <- TRUE

    points.idw <- idw(get(attr) ~ 1, points, newdata = grd, idp = 2.0)
    return(raster(points.idw))
}


#' Perform regression of archaeological dates on great circle distances from
#' a hypothetical origin. Dates can be filtered to retain only the earliest
#' dates per distance bins (Hamilton and Buchanan 2007). Bootstrap is executed
#' to account for uncertainty in calibrated dates. If a cost surface is
#' provided, distances are calculated using least cost paths instead of great
#' circle distances. Regression can be either reduced major axis or ordinary
#' least squares. If using ordinary least squares, regression is performed
#' both on time-versus-distance and on distance-versus-time.
#' 
#' @param ftrSites A SpatialPointsDataFrame object with associated earliest
#' C14 dates per site and respective calibrated distributions (CalDates
#' objects) in a field named "cal". Result of applying filterDates() and
#' calAll().
#' @param c14bp A string. Name of the field with the radiocarbon ages in C14
#' BP format.
#' @param origin A SpatialPointsDataFrame object. The site considered as
#' hypothetical origin of expansion.
#' @param binWidth A number. Width of the spatial bins in km, calculated as
#' distance intervals from the hypothetical origin. Default is 0 (no bins).
#' @param nsim A number. Number of simulations to be run during the
#' bootstrapping procedure. Default is 999.
#' @param cost Optional. A RasterLayer with friction values to calculate
#' least cost path distances instead of great circle distances.
#' @param method A string. Method to be used in the regression. One of "rma"
#' or "ols". Default is "rma".
#' @return a dateModel object.
#' @export
modelDates <- function(ftrSites, c14bp, origin, binWidth = 0, nsim = 999,
                       cost = NULL, method = "rma") {

    no_cores <- detectCores()
    cl <- makeCluster(no_cores - 1)
    clusterEvalQ(cl, library("rcarbon"))
    clusterEvalQ(cl, library("smatr"))
    clusterExport(cl, "sampleCalDates", envir = .GlobalEnv)

    # If a cost surface is provided, calculate cost distances
    if (!missing(cost) & class(cost) == "RasterLayer") {
        tr <- transition(cost, function(x) 1/mean(x), 16)
        tr <- geoCorrection(tr)

        # Create cost paths and calculate their distances
        paths <- shortestPath(tr, origin, ftrSites, output = "SpatialLines")
        ftrSites$dists <- SpatialLinesLengths(paths)
    } else {
        # If no cost surface is provided, caculate great circle distances
        ftrSites$dists <- spDistsN1(ftrSites, origin, longlat = TRUE)
    }

    # Perform spatial binning, retain only earliest date per bin
    if (binWidth > 0)  {
        ftrSites$bin <- floor(ftrSites$dists / binWidth)
        ftrSites.df <- as.data.frame(ftrSites)
        ftrSites.max <- ftrSites.df %>% group_by(bin) %>% top_n(1, get(c14bp))
        binSites <- as.data.frame(ftrSites.max)
        dists <- binSites$dists
        calDates <- binSites$cal
    } else {
        binSites <- NULL
        dists <- ftrSites$dists
        calDates <- ftrSites$cal
    }

    if (method == "rma") {
        
        models <- vector("list", length = nsim)
    
        clusterExport(cl, list("models", "dists", "calDates"),
                      envir = environment())
    
        models <- parLapply(cl, 1:nsim, function(k) {
            dates <- sampleCalDates(calDates)
            model <- sma(dates ~ dists, robust = TRUE)
        })

        stopCluster(cl)
        gc()

        res <- matrix(nrow = nsim, ncol = 4)

        # Save the coefficient, probability, slope and intercept of the
        # regressions
        colnames(res) <- c("r", "p", "slo", "int")
        lapply(1:nsim, function(k) {
            res[k, "r"] <<- models[[k]]$r[[1]]
            res[k, "p"] <<- models[[k]]$p[[1]]
            res[k, "slo"] <<- models[[k]]$coef[[1]][2, 1]
            res[k, "int"] <<- models[[k]]$coef[[1]][1, 1]
        })

        dateModel <- list("model" = res, "binSites" = binSites,
                          "allSites" = ftrSites, "binWidth" = binWidth,
                          "method" = method)

    } else if (method == "ols") {

        td.models <- vector("list", length = nsim)
        dt.models <- vector("list", length = nsim)
    
        clusterExport(cl, list("td.models", "dt.models", "dists", "calDates"),
                      envir = environment())

        # Time-versus-distance
        td.models <- parLapply(cl, 1:nsim, function(k) {
            dates <- sampleCalDates(calDates)
            model <- lm(dates ~ dists)
        })

        # Distance-versus-time
        dt.models <- parLapply(cl, 1:nsim, function(k) {
            dates <- sampleCalDates(calDates)
            model <- lm(dists ~ dates)
        })

        stopCluster(cl)
        gc()

        td.res <- matrix(nrow = nsim, ncol = 4)
        dt.res <- matrix(nrow = nsim, ncol = 4)

        # Save the coefficient, probability, slope and intercept of the
        # regressions
        colnames(td.res) <- c("r", "p", "slo", "int")
        colnames(dt.res) <- c("r", "p", "slo", "int")

        lapply(1:nsim, function(k) {
            # Time-versus-distance
            td.res[k, "r"] <<- sqrt(summary(td.models[[k]])$r.squared)
            td.res[k, "p"] <<- summary(td.models[[k]])$coefficients[2, 4]
            td.res[k, "slo"] <<- td.models[[k]]$coefficients[[2]]
            td.res[k, "int"] <<- td.models[[k]]$coefficients[[1]]

            # Distance-versus-time
            int <- dt.models[[k]]$coefficients[[1]]
            slo <- dt.models[[k]]$coefficients[[2]]
            dt.res[k, "r"] <<- sqrt(summary(dt.models[[k]])$r.squared)
            dt.res[k, "p"] <<- summary(dt.models[[k]])$coefficients[2, 4]
            
            # Invert the intercept and slope for the plot
            dt.res[k, "slo"] <<- 1 / slo
            dt.res[k, "int"] <<- -int / slo
        })
        
        dateModel <- list("td.model" = td.res, "dt.model" = dt.res,
                          "binSites" = binSites, "allSites" = ftrSites,
                          "binWidth" = binWidth, "method" = method)
    }

    class(dateModel) <- "dateModel"

    return(dateModel)
}


#' Plot a map with the results of the iteration over multiple sites, showing
#' an interpolated surface of r values for the sites considered as potential
#' origins.
#'
#' @param dateMap a dateMap object. One of the list elements returned from
#' the iterateSites() function.
#' @return a spplot object.
#' @export
plot.dateMap <- function(dateMap) {
    data(land)
    e <- extent(dateMap$sites)
    sites <- list("sp.points", dateMap$sites, cex = 0.25, col = "black",
                  alpha = 0.5)
    continent <- list("sp.polygons", land, fill = "white")
    plt <- spplot(mask(dateMap$idw, land), cuts = 8,
                  col.regions = magma,
                  xlim = c(e[1], e[2]), ylim = c(e[3], e[4]),
                  xlab = list("R", cex = 0.75),
                  par.settings = list(panel.background = list(col = "lightcyan2"),
                                      layout.heights = list(xlab.key.padding = 1)),
                  sp.layout = list(sites, continent),
                  scales = list(draw = TRUE, cex = 0.5, font = 1,
                                tck = c(1, 0)),
                  colorkey = list(height = 0.75, space = "bottom"),
                  main = list(label = "Interpolation map of potential origins' R values",
                              cex = 1))
    return(plt)
}


#' Plot the results of the regression bootstrap on archaeological dates versus
#' distances from a hypothetical origin, with fitted lines.
#'
#' @param dateModel a dateModel object created with modelDates().
#' @return a ggplot object.
#' @export
plot.dateModel <- function(dateModel) {

    if (!is.null(dateModel$binSites)) {
        binSites <- data.frame("dists" = dateModel$binSites$dists,
                               "dates" = dateModel$binSites$med,
                               "binned" = 1)
    } else {
        binSites <- NULL
    }

    allSites <- data.frame("dists" = dateModel$allSites$dists,
                           "dates" = dateModel$allSites$med,
                           "binned" = 0)

    if (!is.null(binSites)) {
        points <- rbind(allSites, binSites)
        binLabel <- paste(dateModel$binWidth, " km bins")
    } else {
        points <- allSites
        binLabel <- NULL
    }

    if (dateModel$method == "rma") {

        model <- as.data.frame(dateModel$model)

        # Mean intercept and slope of bootstrap for display
        int.m <- mean(model$int)
        slo.m <- mean(model$slo)

        # Display significant p value as either < 0.01 or < 0.05. If not
        # significant, display actual value up to two decimals.
        if (mean(model$p) < 0.01) {
            pval <- "p < 0.01"
        } else if (mean(model$p) < 0.05) {
            pval <- "p < 0.05"
        } else {
            pval <- paste("p = ", format(mean(model$p), digits = 2,
                                        scientific = FALSE))
        }

        rval <- paste("r = -", format(sqrt(mean(mean(model$r))), dig = 2),
                      sep="")

        plt <- ggplot() + xlim(0, max(points$dists)) +
                          ylim(min(points$dates), max(points$dates)) +
                          geom_abline(data = model,
                                      aes(intercept = int, slope = slo),
                                      alpha = 0.05, lwd = 0.5,
                                      colour = "#f8766d") +
                          geom_abline(aes(intercept = int.m, slope = slo.m)) +
                          geom_point(data = points, aes(x = dists, y = dates,
                                                        fill = factor(binned)),
                                     shape = 21, size = 2) +
                          labs(x = "Distance from origin (km)", y = "Cal yr BP",
                               fill = binLabel) +
                          scale_fill_manual(labels = c("all dates",
                                                       "earliest per bin"),
                                            values = c("white", "black")) +
                          annotate("text", x = min(points$dists),
                                           y = min(points$dates),
                                           label = paste(rval, pval),
                                           hjust = 0) +
                          theme(legend.position = c(1, 1),
                                legend.justification = c(1, 1))

    } else if (dateModel$method == "ols") {

        # Time-versus-distance
        td.model <- as.data.frame(dateModel$td.model)
        td.int.m <- mean(td.model$int)
        td.slo.m <- mean(td.model$slo)

        if (mean(td.model$p) < 0.01) {
            td.pval <- "p < 0.01"
        } else if (mean(td.model$p) < 0.05) {
            td.pval <- "p < 0.05"
        } else {
            td.pval <- paste("p = ", format(mean(td.model$p), digits = 2,
                                            scientific = FALSE))
        }

        td.rval <- paste("r = -", format(sqrt(mean(mean(td.model$r))), dig = 2),
                         sep="")

        # Distance-versus-time
        dt.model <- as.data.frame(dateModel$dt.model)
        dt.int.m <- mean(dt.model$int)
        dt.slo.m <- mean(dt.model$slo)

        if (mean(dt.model$p) < 0.01) {
            dt.pval <- "p < 0.01"
        } else if (mean(dt.model$p) < 0.05) {
            dt.pval <- "p < 0.05"
        } else {
            dt.pval <- paste("p = ", format(mean(dt.model$p), digits = 2,
                                            scientific = FALSE))
        }

        dt.rval <- paste("r = -", format(sqrt(mean(mean(dt.model$r))), dig = 2),
                         sep="")

        plt <- ggplot() + xlim(0, max(points$dists)) +
                          ylim(min(points$dates), max(points$dates)) +
                          geom_abline(data = dt.model,
                                      aes(intercept = int, slope = slo),
                                      alpha = 0.05, lwd = 0.5,
                                      colour = "#00bfc4") +
                          geom_abline(data = td.model,
                                      aes(intercept = int, slope = slo),
                                      alpha = 0.05, lwd = 0.5,
                                      colour = "#f8766d") +
                          geom_abline(aes(intercept = dt.int.m,
                                          slope = dt.slo.m), lty = 2) +
                          geom_abline(aes(intercept = td.int.m,
                                          slope = td.slo.m)) +
                          geom_point(data = points, aes(x = dists, y = dates,
                                                        fill = factor(binned)),
                                     shape = 21, size = 2) +
                          labs(x = "Distance from origin (km)", y = "Cal yr BP",
                               fill = binLabel) +
                          scale_fill_manual(labels = c("all dates",
                                                       "earliest per bin"),
                                            values = c("white", "black")) +
                          annotate("text", x = min(points$dists),
                                           y = min(points$dates),
                                           label = paste(td.rval, td.pval),
                                           hjust = 0) +
                          theme(legend.position = c(1, 1),
                                legend.justification = c(1, 1))

    }

    return(plt)
}


#' Extract a single year estimate for ranges of calibrated dates with a
#' probability given by the calibrated probability distribution.
#'
#' @param calDates A CalDates object or a vector of CalDates.
#' @return A vector of cal BP single year estimates.
#' @export
sampleCalDates <- function(calDates) {

    years <- numeric(length(calDates))
    
	for (i in 1:length(calDates)) {
        years[i] <- sample(calDates[i]$grids[[1]]$calBP, size = 1,
                           prob = calDates[i]$grids[[1]]$PrDens)
    }

    return(years)
}


#' Return summary statistics for a dateModel object.
#'
#' @param dateModel A dateModel object created with modelDates().
#' @return A dataframe with estimated start date and speed of advance.
#' @export
summary.dateModel <- function(dateModel) {

    if (dateModel$method == "rma") {

        start <- mean(dateModel$model[, "int"])
        start.SD <- sd(dateModel$model[, "int"])
        speed <- 1 / mean(dateModel$model[, "slo"])
        speed.SD <- (1 / ((mean(dateModel$model[, "slo"])) +
                    (1.96 * sd(dateModel$model[, "slo"])))) - speed

        df <- data.frame(c(paste(format(start, digits = 0, scientific = FALSE),
                                 "+/-", format(start.SD, digits = 0,
                                               scientific = FALSE), "cal BP"),
                           paste(format(abs(speed), digits = 2,
                                        scientific = FALSE),
                                 "+/-", format(abs(speed.SD), digits = 2,
                                               scientific = FALSE),
                                 "km/yr")))

        rownames(df) <- c("Start date:", "Speed of advance:")
        colnames(df) <- NULL

    } else if (dateModel$method == "ols") {

        # Time-versus-distance
        td.start <- mean(dateModel$td.model[, "int"])
        td.start.SD <- sd(dateModel$td.model[, "int"])
        td.speed <- 1 / mean(dateModel$td.model[, "slo"])
        td.speed.SD <- (1 / ((mean(dateModel$td.model[, "slo"])) +
                       (1.96 * sd(dateModel$td.model[, "slo"])))) - td.speed
        
        # Distance-versus-time
        dt.start <- mean(dateModel$dt.model[, "int"])
        dt.start.SD <- sd(dateModel$dt.model[, "int"])
        dt.speed <- 1 / mean(dateModel$dt.model[, "slo"])
        dt.speed.SD <- (1 / ((mean(dateModel$dt.model[, "slo"])) +
                       (1.96 * sd(dateModel$dt.model[, "slo"])))) - dt.speed

        df <- data.frame(c(paste(format(td.start, digits = 0,
                                        scientific = FALSE), "+/-",
                                 format(td.start.SD, digits = 0,
                                        scientific = FALSE), "cal BP"),
                           paste(format(abs(td.speed), digits = 2,
                                        scientific = FALSE),
                                 "+/-", format(abs(td.speed.SD), digits = 2,
                                               scientific = FALSE),
                                 "km/yr")),
                         c(paste(format(dt.start, digits = 0,
                                        scientific = FALSE), "+/-",
                                 format(dt.start.SD, digits = 0,
                                        scientific = FALSE), "cal BP"),
                           paste(format(abs(dt.speed), digits = 2,
                                        scientific = FALSE),
                                 "+/-", format(abs(dt.speed.SD), digits = 2,
                                               scientific = FALSE),
                                 "km/yr")))

        rownames(df) <- c("Start date:", "Speed of advance:")
        colnames(df) <- c("Time-versus-distance", "Distance-versus-time")

    }

    return(df)
}


#' Radiocarbon dates and coordinates of 717 Neolithic sites in the Near East
#' and Europe. Modified from Pinhasi et al. (2005). Only the earliest dates
#' per site are included.
#' 
#' @format A data frame with 717 rows and 13 variables.
#' \itemize{
#'   \item Latitude. Site latitude in decimal degrees.
#'   \item Longitude. Site longitude in decimal degrees.
#'   \item Site. Site name.
#'   \item Location. Region where the site is located (Near East, Europe etc).
#'   \item Country. Country where the site is located.
#'   \item Period. Site period or culture (PPNA, PPNB, LBK etc.).
#'   \item LabNumber. Laboratory number of the C14 date.
#'   \item C14Age. Date in C14 years BP.
#'   \item C14SD. Standard error of the radiocarbon date.
#'   \item Material. Material dated (Charcoal, shell etc.).
#'   \item Curve. Curve to be used in the calibration of each date (intcal13,
#'                marine13).
#'   \item cal. Calibrated dates as CalDates objects.
#'   \item med. Median of the calibrated date in cal yr BP.
#' }
"neof"


#' Coordinates of 9 sites considered as potential centers of origin of the
#' Neolithic expansion. Modified from Pinhasi et al. (2005).
#' 
#' @format A data frame with 9 rows and 3 variables.
#' \itemize{
#'   \item Latitude. Site latitude in decimal degrees.
#'   \item Longitude. Site longitude in decimal degrees.
#'   \item Site. Site name.
#' }
"centers"


#' Cost surface to calculate shortest paths (least-cost paths) of Neolithic
#' expansion from the Near East to Europe. Values are 1 (easiest, coast), 2
#' (land below 1750 m) and 3 (hardest, ocean and land above 1750 m).
#' 
#' @format A RasterLayer object.
"cost"


#' Land polygons.
#' 
#' @format A SpatialPolygonsDataFrame object.
"land"