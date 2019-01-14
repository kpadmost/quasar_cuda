library(QuasarFitCuda)

continuumFitting <- function() {
  
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
  # fill till blockSize with 0 -s
  #q
  #
  #
  fluxMatrix <- getFluxMatrix(quasars)
  errorMatrix <- getErrorMatrix(quasars)
  sizesVector <- getSizesVector(quasars)
  
  params <- getLambdaParams(quasars)
  wm <- cppGenerateWavelenghtMatrix(params)
  wm <- t(wm)
  #testFittingQ()
  fluxMatrix <- cppMovingAverage(fluxMatrix, sizesVector, 10)
  matrices <- cppFilterWithValues(wm, fluxMatrix, errorMatrix, sizesVector)
  flux <- matrices$spectrumsMatrix
  wavelength <- matrices$wavelengthMatrix
  f1 <- flux[1,]
  l1 <- wm[1,]
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
