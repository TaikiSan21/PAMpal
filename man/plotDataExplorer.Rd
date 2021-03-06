% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotDataExplorer.R
\name{plotDataExplorer}
\alias{plotDataExplorer}
\title{Explore Data in an Interactive Plot}
\usage{
plotDataExplorer(x, callType = NULL, maxCategories = 15)
}
\arguments{
\item{x}{data to plot, can be an \code{AcousticStudy}, \code{AcousticEvent},
data.frame or a list of \code{AcousticEvent} objects}

\item{callType}{the specific type of call to plot. If \code{NULL} (default),
will prompt user to choose which type if more than one is present.}

\item{maxCategories}{maximum number of categories to color and facet by. Only
character and factor data with a number of unique values less than or equal
to this number will be shown as options for selecting colors and facets. Not
recommended to increase this value much beyond 20, trying to plot a large number
of colors will cause R to be sad.}
}
\value{
nothing, just plots
}
\description{
Creates an interactive plot of detector data. Allows user to
  choose which numeric data to plot, and will allow user to both color and
  facet the plot by any columns that are characters or factors
}
\examples{

data(exStudy)

if(interactive()) plotDataExplorer(exStudy)
if(interactive()) plotDataExplorer(exStudy, callType='click')

}
\author{
Taiki Sakai \email{taiki.sakai@noaa.gov}
}
