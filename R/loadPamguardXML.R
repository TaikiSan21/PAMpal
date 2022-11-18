#' @title Load Pamguard XML Settings
#'
#' @description Loads in relevant settings and formats for use in PAMpal 
#'
#' @param x an XML file created by Pamguard's "Export XML Configuration"
#' 
#' @return A list with settings for audio \code{sources} (sound acuisition, decimators,
#'   FFT, and cepstrum) and \code{detectors} (click detector and whistle and moan 
#'   detector). Also stores the entire XML file as \code{raw} and the file name as
#'   \code{file}
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @examples
#'
#' xmlFile <- system.file('extdata', 'Example.xml', package='PAMpal')
#' xmlList <- loadPamguardXML(xmlFile)
#' str(xmlList)
#'
#' @importFrom xml2 read_xml xml_find_all xml_attr xml_name xml_attrs xml_children
#' @export
#'
loadPamguardXML <- function(x) {
  if(!file.exists(x)) {
    warning('File ', x, ' does not exist')
    return(NULL)
  }
  x <- normalizePath(x, winslash = '/')
  pgxml <- read_xml(x)
  psf <- xml_find_all(pgxml, 'INFO') %>% 
    xml_attr('CONFIGURATION')
  pxSources <- xml_getSources(pgxml)
  pxDetectors <- xml_getDetectors(pgxml, pxSources)
  names(pxDetectors) <- gsub(' ', '_', names(pxDetectors))
  list(file = x,
       sources=pxSources,
       detectors=pxDetectors,
       raw = as.character(pgxml))
}

xml_getFFT <- function(pgxml) {
  fftNodes <- xml_find_all(pgxml, '//MODULE[@UnitType="FFT Engine"]')
  if(length(fftNodes) == 0) {
    return(NULL)
  }
  fftList <- vector('list', length = length(fftNodes))
  names(fftList) <- xml_find_all(fftNodes, 'PROCESS[starts-with(@Name, "FFT -")]/Output') %>% 
    xml_attr('Name')
  fftHops <- xml_find_all(fftNodes,'//CONFIGURATION/SETTINGS/fftHop') %>% 
    xml_attr('Value') %>% 
    as.numeric()
  fftInputs <- xml_find_all(fftNodes, '//CONFIGURATION/SETTINGS/dataSourceName') %>% 
    xml_attr('Value')
  fftLengths <- xml_find_all(fftNodes, '//CONFIGURATION/SETTINGS/fftLength') %>% 
    xml_attr('Value') %>% 
    as.numeric()
  for(i in seq_along(fftList)) {
    fftList[[i]] <- list(source=fftInputs[i],
                         length=fftLengths[i],
                         hop=fftHops[i])
  }
  fftList
}

xml_getSa <- function(pgxml) {
  saNode <- xml_find_all(pgxml, '//MODULE[@UnitName="Sound Acquisition"]/PROCESS/Output[starts-with(@Type, "PamDetection")]') %>% 
    xml_attrs()
  if(length(saNode) == 0) {
    return(NULL)
  }
  saList <- vector('list', length=1)
  saList[[1]] <- list(sr = as.numeric(saNode[[1]][['SampleRate']]))
  names(saList) <- saNode[[1]][['Name']]
  saList
}

xml_getCepstrum <- function(pgxml, fftList) {
  cepsNodes <- xml_find_all(pgxml, '//MODULE[@UnitType="Cepstrum"]')
  if(length(cepsNodes) == 0) {
    return(NULL)
  }
  cepsList <- vector('list', length = length(cepsNodes))
  names(cepsList) <- xml_find_all(cepsNodes, 'PROCESS/Output') %>% 
    xml_attr('Name')
  cepsInput <- xml_find_all(cepsNodes, 'PROCESS/Input') %>% 
    xml_attr('Name')
  for(i in seq_along(cepsList)) {
    thisFft <- fftList[[cepsInput[i]]]
    cepsList[[i]] <- list(source = cepsInput[i],
                          length = thisFft$length,
                          hop = thisFft$hop,
                          sr = thisFft$sr)
  }
  cepsList
}

xml_getDecimator <- function(pgxml) {
  deciNodes <- xml_find_all(pgxml, '//MODULE[@UnitType="Decimator"]')
  if(length(deciNodes) == 0) {
    return(NULL)
  }
  deciList <- vector('list', length = length(deciNodes))
  names(deciList) <- xml_find_all(deciNodes, 'PROCESS/Output') %>% 
    xml_attr('Name')
  deciInput <- xml_find_all(deciNodes, 'PROCESS/Input') %>% 
    xml_attr('Name')
  deciSr <- xml_find_all(deciNodes, 'PROCESS/Output') %>% 
    xml_attr('SampleRate') %>% 
    as.numeric()
  for(i in seq_along(deciList)) {
    deciList[[i]] <- list(source = deciInput[i],
                          sr = deciSr[i])
  }
  deciList
}

xml_getWMD <- function(pgxml) {
  wmdNodes <- xml_find_all(pgxml, '//MODULE[@UnitType="WhistlesMoans"]')
  if(length(wmdNodes) == 0) {
    return(NULL)
  }
  wmdList <- vector('list', length = length(wmdNodes))
  names(wmdList) <- xml_attr(wmdNodes, 'UnitName')
  wmdInputs <- xml_find_all(wmdNodes, 'CONFIGURATION/SETTINGS/dataSource') %>% 
    xml_attr('Value')
  wmdMinFreq <- xml_find_all(wmdNodes, 'CONFIGURATION/SETTINGS/minFrequency') %>% 
    xml_attr('Value') %>% 
    as.numeric()
  wmdMinLength <- xml_find_all(wmdNodes, 'CONFIGURATION/SETTINGS/minLength') %>% 
    xml_attr('Value') %>% 
    as.numeric()
  wmdMinPixels <- xml_find_all(wmdNodes, 'CONFIGURATION/SETTINGS/minPixels') %>% 
    xml_attr('Value') %>% 
    as.numeric()
  wmdMaxCross <- xml_find_all(wmdNodes, 'CONFIGURATION/SETTINGS/maxCrossLength') %>% 
    xml_attr('Value') %>% 
    as.numeric()
  wmdShortLength <- xml_find_all(wmdNodes, 'CONFIGURATION/SETTINGS/shortLength') %>% 
    xml_attr('Value') %>% 
    as.numeric()
  wmdMedianFilter <- xml_find_all(wmdNodes, 'CONFIGURATION/SETTINGS/specNoiseSettings//filterLength') %>% 
    xml_attr('Value') %>% 
    as.numeric()
  wmdAvgSub <- xml_find_all(wmdNodes, 'CONFIGURATION/SETTINGS/specNoiseSettings//updateConstant') %>% 
    xml_attr('Value') %>% 
    as.numeric()
  wmdThreshold <- xml_find_all(wmdNodes, 'CONFIGURATION/SETTINGS/specNoiseSettings//thresholdDB') %>% 
    xml_attr('Value') %>% 
    as.numeric()
  for(i in seq_along(wmdNodes)) {
    wmdList[[i]] <- list(
      source = wmdInputs[i],
      type = ifelse(grepl('Cepstrum', wmdInputs[i]), 'cepstrum', 'whistle'),
      minFrequency = wmdMinFreq[i],
      minLength = wmdMinLength[i],
      minPixels = wmdMinPixels[i],
      maxCross = wmdMaxCross[i],
      shortLength = wmdShortLength[i],
      medianFilter = wmdMedianFilter[i],
      avgSubtraction = wmdAvgSub[i],
      threshold = wmdThreshold[i]
    )
  }
  wmdList
}

xml_getClick <- function(pgxml) {
  clickNodes <- xml_find_all(pgxml, '//MODULE[@UnitType="Click Detector"]')
  if(length(clickNodes) == 0) {
    return(NULL)
  }
  clickList <- vector('list', length = length(clickNodes))
  names(clickList) <- xml_attr(clickNodes, 'UnitName')
  longFilter <- xml_find_all(clickNodes, 'CONFIGURATION/SETTINGS[@Type="Click Detector"]/longFilter') %>% 
    xml_attr('Value') %>% 
    as.numeric()
  longFilter2 <- xml_find_all(clickNodes, 'CONFIGURATION/SETTINGS[@Type="Click Detector"]/longFilter2') %>% 
    xml_attr('Value') %>% 
    as.numeric()
  shortFilter <- xml_find_all(clickNodes, 'CONFIGURATION/SETTINGS[@Type="Click Detector"]/shortFilter') %>% 
    xml_attr('Value') %>% 
    as.numeric()
  threshold <- xml_find_all(clickNodes, 'CONFIGURATION/SETTINGS[@Type="Click Detector"]/dbThreshold') %>% 
    xml_attr('Value') %>% 
    as.numeric()
  preFilterNodes <- xml_find_all(clickNodes, 'CONFIGURATION/SETTINGS[@Type="Click Detector"]/preFilter')
  preFilter <- lapply(preFilterNodes, function(x) {
    vals <- xml_attr(xml_children(x), 'Value')
    names(vals) <- xml_name(xml_children(x))
    vals
  })
  
  triggerFilterNodes <- xml_find_all(clickNodes, 'CONFIGURATION/SETTINGS[@Type="Click Detector"]/triggerFilter')
  triggerFilter <- lapply(triggerFilterNodes, function(x) {
    vals <- xml_attr(xml_children(x), 'Value')
    names(vals) <- xml_name(xml_children(x))
    vals
  })
  basicIdNodes <- xml_find_all(clickNodes, 'CONFIGURATION/SETTINGS[@Type="BasicClickIdParams"]/clickTypeParams')
  basicId <- sapply(basicIdNodes, function(x) {
    tmp <- lapply(xml_children(x), function(y) {
      nds <- xml_find_all(y, '*[@Value]')
      lst <- xml_attr(nds, 'Value')
      names(lst) <- xml_name(nds)
      lst
    })
    names(tmp) <- sapply(tmp, function(y) y[['speciesCode']])
    tmp
  })

  for(i in seq_along(clickList)) {
    input <- xml_find_all(clickNodes, paste0('PROCESS[@Name="',
                                             names(clickList)[i],
                                             '"]/Input')) %>% 
      xml_attr('Name')
    result <- list(longFilter=longFilter[i],
                   longFilter2=longFilter2[i],
                   shortFilter=shortFilter[i],
                   threshold=threshold[i],
                   preFilter=preFilter[[i]],
                   triggerFilter=triggerFilter[[i]])
    result$source <- input
    result$type <- 'click'
    result$basicClass <- basicId[[i]]
    # sweeps may not be present
    sweepIdNodes <- xml_find_all(clickNodes[[i]], 'CONFIGURATION/SETTINGS[@Type="ClickSweepClassifier"]/classifierSets')
    if(length(sweepIdNodes) != 0) {
      sweepId <- lapply(sweepIdNodes, function(x) {
        tmp <- lapply(xml_children(x), function(y) {
          # browser()
          nds <- xml_find_all(y, '*[@Value]')
          lst <- xml_attr(nds, 'Value')
          names(lst) <- xml_name(nds)
          lst
        })
        names(tmp) <- sapply(tmp, function(y) y[['speciesCode']])
        tmp
      })[[1]]
    } else {
      sweepId <- list()
    }
    result$sweepClass <- sweepId
    clickList[[i]] <- result
  }
  clickList
}


xml_getSources <- function(pgxml) {
  sa <- xml_getSa(pgxml)
  deci <- xml_getDecimator(pgxml)
  wavSource <- c(sa, deci)
  fft <- xml_getFFT(pgxml)
  for(i in seq_along(fft)) {
    fft[[i]]$sr <- wavSource[[fft[[i]]$source]]$sr
    fft[[i]]$decimated <- grepl('Decimator', fft[[i]]$source)
  }
  ceps <- xml_getCepstrum(pgxml, fft)
  for(i in seq_along(ceps)) {
    ceps[[i]]$sr <- fft[[ceps[[i]]$source]]$sr
    ceps[[i]]$decimated <- fft[[ceps[[i]]$source]]$decimated
  }
  c(wavSource, fft, ceps)
}

xml_getDetectors <- function(pgxml, sources) {
  wmd <- xml_getWMD(pgxml)
  for(i in seq_along(wmd)) {
    thisFft <- sources[[wmd[[i]]$source]]
    wmd[[i]]$sr <- thisFft$sr
    wmd[[i]]$length <- thisFft$length
    wmd[[i]]$hop <- thisFft$hop
    wmd[[i]]$decimated <- thisFft$decimated
  }
  click <- xml_getClick(pgxml)
  for(i in seq_along(click)) {
    thisWav <- sources[[click[[i]]$source]]
    click[[i]]$sr <- thisWav$sr
    click[[i]]$decimated <- grepl('Decimator', click[[i]]$source)
  }
  c(wmd, click)
}
