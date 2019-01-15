library(QuasarFitCuda)

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
  
  
  print('here')
}

feFitting <- function() {
  
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
  filteredSpectrums <- rFilterInfs(filteredMatrices$spectrumsMatrix)
  filteredWavelengths <- rFilterInfs(filteredMatrices$wavelengthsMatrix)
  filteredErrors <- rFilterInfs(filteredMatrices$errorsMatrix)
  
  for(i in seq(1:3)) {
    continuumResult <- continuumFitting(
      filteredWavelengths,
      filteredSpectrums,
      filteredErrors,
      sizesVector,
      continuumParams
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
