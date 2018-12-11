library(QuasarFitCuda)


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
  params <- getLambdaParams(quasars)
  print('finish')
  # $wavelengthMatrix
  # $spectrumMatrix
  # $errorMatrix
}

testFitting()
