% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClass.R
\docType{class}
\name{AcousticEvent-class}
\alias{AcousticEvent-class}
\title{\code{AcousticEvent} Class}
\description{
An S4 class storing acoustic detections from an Acoustic Event
  as well as other related metadata
}
\section{Slots}{

\describe{
\item{\code{id}}{unique id or name for this event}

\item{\code{detectors}}{a list of data frames that have acoustic detections and
any measurements calculated on those detections. Each data frame is named
by the detector that made the detection}

\item{\code{localizations}}{a named list storing localizations, named by method}

\item{\code{settings}}{a list of recorder settings}

\item{\code{species}}{a list of species classifications for this event, named by
classification method (ie. BANTER model, visual ID)}

\item{\code{files}}{a list of files used to create this object, named by the type of
file (ie. binaries, database)}

\item{\code{ancillary}}{a list of miscellaneous extra stuff. Store whatever you want here}
}}

\author{
Taiki Sakai \email{taiki.sakai@noaa.gov}
}
