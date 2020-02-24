library(base)
library(data.table)
library(dplyr)
library(gdistance)
library(ggplot2)
library(parallel)
library(raster)
library(rcarbon)
library(smatr)


#' Calibrate all radiocarbon dates in a dataset, differentiating between
#' terrestrial and marine (shell) samples.
#'
#' @param sites A dataframe with, minimally, the C14 ages, standard errors and
#' materials dated.
#' @param c14bp A string. Name of the field with the radiocarbon ages in C14
#' BP format.
#' @param sd A string. Name of the field with the standard errors of the C14
#' dates.
#' @param material A string. Name of the field with the materials dated.
#' @param curve A string. Calibration curve as per the rcarbon package.
#' Default is "intcal13".
#' @return A dataframe with an added column "cal" with CalDates objects from
#' the rcarbon package.
#' @export
calAll <- function(sites, c14bp, sd, material, curve = "intcal13") {
    terrestrial <- sites[!grepl("shell", sites[[material]],
                         ignore.case = TRUE),]
    marine <- sites[grepl("shell", sites[[material]], ignore.case = TRUE),]
    terrestrial$cal <- calibrate(terrestrial[[c14bp]], terrestrial[[sd]],
                                 calCurves = curve)
    terrestrial$med <- medCal(terrestrial$cal)
    if (nrow(marine) > 0) {
        marine$cal <- calibrate(marine[[c14bp]], marine[[sd]],
                                calCurves = "marine13")
        marine$med <- medCal(marine$cal)
        return(rbind(terrestrial, marine))
    } else {
        return(terrestrial)
    }
}


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
#' objects) in a field named "cal". Result of applying filterDates() and
#' calAll().
#' @param c14bp A string. Name of the field with the radiocarbon ages in C14
#' BP format.
#' @param siteNames A string. Name of the field with the site names or ids.
#' @param origins A SpatialPointsDataFrame object. The sites to be tested for
#' the most likely origin of expansion.
#' @param binWidths A number or vector of numbers. Width(s) of the spatial bins
#' in km.
#' @param nsim A number. Number of simulations to be run during the
#' bootstrapping procedure. Default is 999.
#' @param cost A RasterLayer with friction values to calculate least cost
#' path distances instead of great circle distances.
#' @return a list with two elements, the result of the iteration over all
#' potential origins and the best model selected among those.
#' @export
iterateSites <- function(ftrSites, c14bp, siteNames, origins, binWidths,
                         nsim = 999, cost = NULL) {

    datalen <- length(binWidths) * length(origins)

    res <- data.table("r" = numeric(datalen), "p" = numeric(datalen),
                      "bin" = numeric(datalen), "n" = numeric(datalen),
		              "site" = character(datalen))

    # Perform reduced major axis regression for each bin width and for each
    # hypothetical origin, if more than one is specified
    pb <- txtProgressBar(min = 0, max = datalen, style = 3)
    counter <- 0
    for (i in (1:length(binWidths))) {
        for (j in (1:length(origins))) {
            counter <- counter + 1

            dateModel <- rmaDates(ftrSites, c14bp = c14bp, origin = origins[j,], 
                                  binWidth = binWidths[i], nsim = nsim,
                                  cost = cost)

            # Only proceed if the slope is negative, i.e. dates become more
            # recent with increasing distance
            if (mean(dateModel$model[,"slo"]) < 0) {
                
                # Extract the r and p of each model and store it in the result,
                # together with the corresponding bin width 
                meanr <- mean(dateModel$model[,"r"])
                meanp <- mean(dateModel$model[,"p"])

                # Save the best model
                if (!exists("bestModel")) {
                    bestModel <- dateModel
                } else if (meanp < 0.05 & meanr > mean(bestModel$model[,"r"])) {
                    bestModel <- dateModel
                }

                if (meanp < 0.01) {
                    ast <- "**"
                } else if (meanp < 0.05) {
                    ast <- "* "
                } else {
                    ast <- "  "
                }

                res[counter, "r"] <- round(as.double(sqrt(meanr)), 4)
                res[counter, "p"] <- round(as.double(meanp), 4)
                res[counter, "bin"] <- binWidths[i]
                res[counter, "n"] <- nrow(dateModel$binSites)
                res[counter, "site"] <- paste(toString(origins[[siteNames]][j]),
                                              ast, sep="")
            }
            setTxtProgressBar(pb, counter)
        }
    }
    close(pb)

    res <- res[order(-r)]
    return(list("results" = res, "model" = bestModel))
}


#' Plot the results of the RMA bootstrap on archaeological dates versus
#' distances from a hypothetical origin, with fitted lines.
#'
#' @param dateModel a dateModel object created with rmaDates().
#' @return a ggplot object.
#' @export
plot.dateModel <- function(dateModel) {

    binSites <- data.frame("dists" = dateModel$binSites$dists,
                           "dates" = dateModel$binSites$med,
                           "binned" = 1)
    allSites <- data.frame("dists" = dateModel$allSites$dists,
                           "dates" = dateModel$allSites$med,
                           "binned" = 0)

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

    rval <- paste("r = -", format(sqrt(mean(mean(model$r))), dig = 2), sep="")

    points <- rbind(allSites, binSites)

    plt <- ggplot() + xlim(0, max(allSites$dists)) +
                      ylim(min(allSites$dates), max(model$int)) +
                      geom_abline(data = model,
                                  aes(intercept = int, slope = slo),
                                  alpha = 0.05, lwd=0.5,
                                  colour = "lightcoral") +
                      geom_abline(aes(intercept = int.m, slope = slo.m)) +
                      geom_point(data = points, aes(x = dists, y = dates,
                                 fill = factor(binned)), shape = 21, size = 2) +
                      labs(x = "Distance from origin (km)", y = "Cal yr BP",
                           fill = paste(dateModel$binWidth, " km bins")) +
                      scale_fill_manual(labels = c("all dates",
                                                   "earliest per bin"),
                                        values = c("white", "black")) +
                      annotate("text", x = min(allSites$dists),
                                       y = min(allSites$dates),
                                       label = paste(rval, pval), hjust = 0) +
                      theme(legend.position = c(1, 1),
                            legend.justification = c(1, 1))
    return(plt)
}


#' Perform reduced major axis regression of archaeological dates versus
#' great circle distances from a hypothetical origin. Dates can be filtered
#' to retain only the earliest dates per distance bins (Hamilton and
#' Buchanan 2007). Bootstrap is executed to account for uncertainty in
#' calibrated dates. If a cost surface is provided, distances are calculated
#' using least cost paths instead of great circle distances.
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
#' @param cost A RasterLayer with friction values to calculate least cost
#' path distances instead of great circle distances.
#' @return a dateModel object.
#' @export
rmaDates <- function(ftrSites, c14bp, origin, binWidth = 0, nsim = 999,
                     cost = NULL) {

    no_cores <- detectCores()
    cl <- makeCluster(no_cores - 1)
    clusterEvalQ(cl, library("rcarbon"))
    clusterEvalQ(cl, library("smatr"))
    clusterExport(cl, "sampleDates", envir = .GlobalEnv)

    models <- vector("list", length = nsim)

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
    }

    dists <- binSites$dists
    calDates <- binSites$cal

    clusterExport(cl, list("models", "dists", "calDates"),
                  envir = environment())

    models <- parLapply(cl, 1:nsim, function(k) {
        dates <- sampleDates(calDates)
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
                      "allSites" = ftrSites, "binWidth" = binWidth)
    class(dateModel) <- "dateModel"

    return(dateModel)
}


#' Extract a single year estimate for ranges of calibrated dates with a
#' probability given by the calibrated probability distribution.
#'
#' @param calDates A CalDates object or a vector of CalDates.
#' @return A vector of cal BP single year estimates.
#' @export
sampleDates <- function(calDates) {

    years <- numeric(length(calDates))
    
	for (i in 1:length(calDates)) {
        years[i] <- sample(calDates[i]$grids[[1]]$calBP, size = 1,
                           prob = calDates[i]$grids[[1]]$PrDens)
    }

    return(years)
}


#' Return summary statistics for a dateModel object.
#'
#' @param dateModel A dateModel object created with rmaDates().
#' @return A dataframe with estimated start date and speed of advance.
#' @export
summary.dateModel <- function(dateModel) {
    start <- mean(dateModel$model[, "int"])
    start.SD <- sd(dateModel$model[, "int"])
    speed <- 1 / mean(dateModel$model[, "slo"])
    speed.SD <- (1 / ((mean(dateModel$model[, "slo"])) +
                (1.96 * sd(dateModel$model[, "slo"])))) - speed
    df <- data.frame(c(paste(format(start, digits = 0, scientific = FALSE),
                             "+/-", format(start.SD, digits = 0,
                             scientific = FALSE), "cal BP"),
                       paste(format(abs(speed), digits = 2, scientific = FALSE),
                             "+/-", format(abs(speed.SD), digits = 2,
                             scientific = FALSE), "km per year")))
    rownames(df) <- c("Start date:", "Speed of advance:")
    colnames(df) <- NULL
    return(df)
}


# Load sample Neolithic dates
neo <- read.csv("./data/neolithic.csv")
coordinates(neo) <- ~Longitude+Latitude
projection(neo) <- CRS("+init=epsg:4326")

neo.f <- filterDates(neo, "C14Age")

# Load sample cost surface
cost <- raster("./data/cost.tif")