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
  cppGenerateWavelenghtMatrix(params)
  print('finish')

}

qP <- function() {
  dm <- list(c(1:4))
  rep(dm, 4)
}

testFitting()
