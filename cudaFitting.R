library(QuasarFitCuda)


#TODO: separate reglin kernels
continuumFitting <- function(wavelengthsMatrix, spectrumsMatrix, errorsMatrix, sizesVector, fitParams) {
  
  filteredMatrices <- cppFilterWithWavelengthWindows(
    wavelengthsMatrix,
    spectrumsMatrix,
    errorsMatrix,
    sizesVector,
    continuumWindowsVectorR = fitParams$windows
  )
  # filter again for non-logarithmic values
  filteredMatrices <- cppFilterWithValues(
    filteredMatrices$wavelengthsMatrix,
    filteredMatrices$spectrumsMatrix,
    filteredMatrices$errorsMatrix,
    sizesVector
  )
  sizesModified <- cppCountNInfSize(filteredMatrices$spectrumsMatrix)
  maxSize <- max(sizesModified)
  
  filteredWavelenghts <- cppCopyNInf(filteredMatrices$wavelengthsMatrix, maxSize)
  filteredSpectrums <- cppCopyNInf(filteredMatrices$spectrumsMatrix, maxSize)
  filteredErrors <- cppCopyNInf(filteredMatrices$errorsMatrix, maxSize)
  
  wavelenghtsLog10 <- cppMatrixLog10(filteredWavelenghts);
  spectrumsLog10 <- cppMatrixLog10(filteredSpectrums);
  cReglinResults <- cppReglin(wavelenghtsLog10, spectrumsLog10, sizesModified)
  
  lambdaAmp <- fitParams$lambdaAmplitude
  if(lambdaAmp > 0) {
    lambdaAmpLog <- log10(lambdaAmp)
    wavelenghtsLog10 <- cppMatrixAddScalar(wavelenghtsLog10, lambdaAmpLog)
  }
  
  reglin <- cppReglin(wavelenghtsLog10, spectrumsLog10, sizesModified)
  fixedReglin <- cppFixReglin(cReglinResults, reglin)
  filteredContinuum <- cppCalculateContinuumMatrix(filteredWavelenghts, fixedReglin$cReglinFixed)
  chisqFiltered <- cppChisq(filteredSpectrums, filteredContinuum, filteredErrors, sizesModified)
  chisqFilteredReduced <- cppReduceContinuumChisq(chisqFiltered, sizesModified)
  
  continuumMatrixDcMatrix <- cppCalculateCfunDcfun(
    wavelengthsMatrix, 
    fixedReglin$cReglinFixed, 
    fixedReglin$reglinFixed
  )
  
  result <- list(
    cfun=continuumMatrixDcMatrix$cfun,
    dcfun=continuumMatrixDcMatrix$dcfun,
    fitParams=list(
      cReglin=fixedReglin$cReglinFixed,
      chisq=chisqFilteredReduced
    )
  )
  if(lambdaAmp > 0) {
    result$fitParams$reglin <- fixedReglin$reglinFixed
  }
  result
}

feFitting <- function(spectrum, continuum, feTemplate, feWindows, feFitParams) {
  #expand template
  spectrumsMatrix <- spectrum$spectrums
  wavelengthsMatrix <- spectrum$wavelengths
  sizes <- spectrum$sizes
  expandedTemplate <- cppExpandTemplate(
    wavelengthsMatrix,
    sizes,
    feTemplate$lambda,
    feTemplate$flux,
    feFitParams
  )
  if(feFitParams$fitType == 'fwin' || feFitParams$fitType == 'win') {
    filteredMatrices <- cppFilterWithWavelengthWindows(
      wavelengthsMatrix = wavelengthsMatrix,
      spectrumsMatrix = spectrumsMatrix,
      errorsMatrix = errorsMatrix,
      sizesVectorR = sizes,
      continuumWindowsVectorR = feWindows
    )
    
    spectrumsMatrix <- filteredMatrices$spectrumsMatrix
    wavelengthsMatrix <- filteredMatrices$wavelengthsMatrix
    errorsMatrix <- filteredMatrices$errorsMatrix
    
  }
  # filter template with windows
  
  # scale template 
  
  # calculate params
  
  print('here')
  #scale template
}

elementFitting <- function() {
  
}

loadMatrices <- function() {
  
}

testFitting <- function() {
  # load n spectra
  quas_list <- 'quasar_list'
  folder <- 'spectra'
  quasars <- loadQuasarsFromFolder(folder, quas_list)
  continuumParams <- list(
    windows = loadDefaultContinuumWindows(),
    lambdaAmplitude = 900.0
  )
  
  feFitParams <- loadDefaultFeFitParams()
  feTemplate <- loadDefaultFeTemplate()
  feWindows <- loadDefaultFeWindows()
  
  fluxMatrix <- getFluxMatrix(quasars)
  errorMatrix <- getErrorMatrix(quasars)
  sizesVector <- getSizesVector(quasars)
  
  params <- getLambdaParams(quasars)
  wavelengthsMatrix <- cppGenerateWavelenghtMatrix(params)
  #TODO: correct!
  wavelengthsMatrix <- t(wavelengthsMatrix)
  #testFittingQ()
  fluxMatrix <- cppMovingAverage(fluxMatrix, sizesVector, 10)
  
  # filtering non-negative values
  filteredMatrices <- cppFilterWithValues(wavelengthsMatrix, fluxMatrix, errorMatrix, sizesVector)
  sizesVector <- cppCountNInfSize(filteredMatrices$spectrumsMatrix)
  spectrums <- rFilterInfs(filteredMatrices$spectrumsMatrix)
  wavelengths <- rFilterInfs(filteredMatrices$wavelengthsMatrix)
  errors <- rFilterInfs(filteredMatrices$errorsMatrix)
  
  
  for(i in seq(1:3)) {
    continuumResult <- continuumFitting(
      wavelengths,
      spectrums,
      errors,
      sizesVector,
      continuumParams
    )
    
    feFitResults <- feFitting(
      list(
        wavelengths=wavelengths,
        spectrums=cppMatrixMinusMatrix(spectrums, continuumResult$cfun),
        sizes=sizesVector
        ),
      continuum = continuumResult$cfun,
      feTemplate,
      feWindows,
      feFitParams = feFitParams
    )
  }
  
  
  print('finish')
  
}

testFittingQ <- function() {
  # load n spectra
  quas_list <- 'quasar_list'
  folder <- 'spectra'
  quasars <- loadQuasarsFromFolder(folder, quas_list)
  # fill till blockSize with 0 -s
  #q
  params <- getLambdaParams(quasars)
  wm <- cppGenerateWavelenghtMatrix(params)
  wm <- t(wm)
  
  list(
    flux = getFluxMatrix(quasars),
    error = getErrorMatrix(quasars),
    lambda = wm,
    sizes = getSizesVector(quasars)
  )
}

testFQ <- function() {
  library(QuasarFitCuda)
  m <- testFittingQ()
  cppFilterWithValues(m$lambda, m$flux, m$error, m$sizes)
}

qP <- function() {
  d1 <- list(c(3.57840, 0.00010, 1.80677, 0.00000))
  d2 <- list(c(3.57900, 0.00010, 1.51807, 0.00000))
  dm <- list(rep(0, 4))
  lmx <- rep(dm, 30)
  c(d1, d2, lmx)
}

testFitting()
