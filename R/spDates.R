library(base)
library(dplyr)
library(ggplot2)
library(parallel)
library(raster)
library(rcarbon)
library(smatr)


#' Filter archaeological site coordinates and dates, retaining only the earliest
#' radiocarbon date per site.
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


#' Calibrate all radiocarbon dates in a dataset, differentiating between
#' terrestrial and marine (shell) samples.
#'
#' @param sites A dataframe with, minimally, the C14 ages, standard errors and
#' materials dated.
#' @param c14bp A string. Name of the field with the radiocarbon ages in C14 BP
#' format.
#' @param sd A string. Name of the field with the standard errors of the C14
#' dates.
#' @param material A string. Name of the field with the materials dated.
#' @param curve A string. Calibration curve as per the rcarbon package. Default
#' is "intcal13".
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


#' Perform reduced major axis regression of archaeological dates versus
#' great circle distances from a hypothetical origin. Dates can be filtered
#' to retain only the earliest dates per distance bins (Hamilton and
#' Buchanan 2007). Bootstrap is executed to account for uncertainty in
#' calibrated dates.
#' 
#' @param ftrSites A SpatialPointsDataFrame object with associated earliest
#' C14 dates per site and respective calibrated distributions (CalDates objects)
#' in a field named "cal". Result of applying filterDates() and calAll().
#' @param c14bp A string. Name of the field with the radiocarbon ages in C14 BP
#' format.
#' @param origin A SpatialPointsDataFrame object. The site considered as
#' hypothetical origin of expansion.
#' @param binWidth A number. Width of the spatial bins in km, calculated as
#' distance intervals from the hypothetical origin. Default is 0 (no bins).
#' @param nsim A number. Number of simulations to be run during the
#' bootstrapping procedure. Default is 999.
#' @return a dateModel object.
#' @export
rmaDates <- function(ftrSites, c14bp, origin, binWidth = 0, nsim = 999) {

    no_cores <- detectCores()
    cl <- makeCluster(no_cores - 1)
    clusterEvalQ(cl, library("rcarbon"))
    clusterEvalQ(cl, library("smatr"))
    clusterExport(cl, "sampleDates", envir = .GlobalEnv)

    models <- vector("list", length = nsim)

    # Calculate distances from origin
    ftrSites$dists <- spDistsN1(ftrSites, origin, longlat = TRUE)

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

    # Save the coefficient, probability, slope and intercept of the regressions
    colnames(res) <- c("r", "p", "slo", "int")
    lapply(1:nsim, function(k) {
        res[k, "r"] <<- models[[k]]$r[[1]]
        res[k, "p"] <<- models[[k]]$p[[1]]
        res[k, "slo"] <<- models[[k]]$coef[[1]][2, 1]
        res[k, "int"] <<- models[[k]]$coef[[1]][1, 1]
    })

    dateModel <- list("model" = res, "binSites" = binSites,
                      "allSites" = ftrSites)
    class(dateModel) <- "dateModel"

    return(dateModel)
}


plot.dateModel <- function(model1, model2 = NULL, title = NULL,
						   origin = NULL, bin = NULL) {

    sts <- data.frame("dists" = model1$binSites$dists,
                      "dates" = model1$binSites$med)
    asts <- data.frame("dists" = model1$allSites$dists,
                       "dates" = model1$allSites$med)
    mdl1 <- as.data.frame(model1$model)
    plt <- ggplot() + xlim(0, max(sts$dists)) + ylim(min(model1$allSites$med),
                                                     max(mdl1$int)) +
                      xlab("Distance from origin (km)") + ylab("Cal yr BP")
    
    if (!missing(model2)) {
        mdl <- as.data.frame(model2$model)
        sts2 <- data.frame("dists" = model2$sites$dists,
                           "dates" = model2$sites$med)

        mItc <- mean(mdl$int)
        mSlo <- mean(mdl$slo)        

        plt <- plt + geom_abline(data = mdl,
                                 aes(intercept = int, slope = slo),
                                 alpha = 0.05, lwd=0.5,
                                 colour = "#E64B35B2") +
                     geom_abline(aes(intercept = mItc, slope = mSlo))
    }

    mdl <- as.data.frame(model1$model)
    mItc <- mean(mdl$int)
    mSlo <- mean(mdl$slo)     

    plt <- plt +                  geom_abline(data = mdl,
                             aes(intercept = int, slope = slo),
                             alpha = 0.05, lwd=0.5,
                             colour = "#E64B35B2") +
                 geom_abline(aes(intercept = mItc, slope = mSlo)) +
                geom_point(data = asts, aes(x = dists, y = dates),
                            shape = 21, size = 3, fill = "white") +
                 geom_point(data = sts, aes(x = dists, y = dates),
                            shape = 21, size = 3, fill = "black")
  
    if (!missing(model2)) {

        rval <- paste("r = -", format(sqrt(mean(mean(mdl$r))),
                                               dig = 2), sep="")

        # Display significant p value as either < 0.01 or < 0.05. If not
        # significant, display actual value up to two decimals.
        if (mean(mdl$p) < 0.01) {
            pval <- "p < 0.01"
        } else if (mean(mdl$p) < 0.05) {
            pval <- "p < 0.05"
        } else {
            pval <- paste("p = ", format(mean(mdl$p), digits = 2,
                                        scientific = FALSE))
        }
        
        plt <- plt + geom_point(data = sts2, aes(x = dists, y = dates),
                                shape = 21, size = 3, fill = "black") +
                     annotate("text", x = max(sts2$dists) / 1.1, y = max(sts2$dates),
                              label = paste(rval, pval))
	}
    plt + labs(title = title, subtitle = paste("Origin: ", origin,
											   ", Bins: ", bin, " km",  sep=""))
}


neo <- read.csv("euroevol.csv")
coordinates(neo) <- ~Longitude+Latitude
projection(neo) <- CRS("+init=epsg:4326")

neo.f <- filterDates(neo, "C14Age")