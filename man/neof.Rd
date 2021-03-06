% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spDates.R
\docType{data}
\name{neof}
\alias{neof}
\title{Radiocarbon dates and coordinates of 717 Neolithic sites in the Near East
and Europe. Modified from Pinhasi et al. (2005). Only the earliest dates
per site are included.}
\format{
A data frame with 717 rows and 13 variables.
\itemize{
  \item Latitude. Site latitude in decimal degrees.
  \item Longitude. Site longitude in decimal degrees.
  \item Site. Site name.
  \item Location. Region where the site is located (Near East, Europe etc).
  \item Country. Country where the site is located.
  \item Period. Site period or culture (PPNA, PPNB, LBK etc.).
  \item LabNumber. Laboratory number of the C14 date.
  \item C14Age. Date in C14 years BP.
  \item C14SD. Standard error of the radiocarbon date.
  \item Material. Material dated (Charcoal, shell etc.).
  \item Curve. Curve to be used in the calibration of each date (intcal13,
               marine13).
  \item cal. Calibrated dates as CalDates objects.
  \item med. Median of the calibrated date in cal yr BP.
}
}
\usage{
neof
}
\description{
Radiocarbon dates and coordinates of 717 Neolithic sites in the Near East
and Europe. Modified from Pinhasi et al. (2005). Only the earliest dates
per site are included.
}
\keyword{datasets}
