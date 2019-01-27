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
    wavelenghtsLog10 <- cppMatrixAddScalar(wavelenghtsLog10, -lambdaAmpLog)
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
  errorsMatrix <- spectrum$errors
  sizes <- spectrum$sizes
  continuumC <- continuum
  expandedTemplate <- cppExpandTemplate(
    wavelengthsMatrix,
    sizes,
    feTemplate$lambda,
    feTemplate$flux,
    feFitParams
  )
  expandedTemplateC <- expandedTemplate
  if(feFitParams$fitType == 'fwin' || feFitParams$fitType == 'win') {
    filteredMatrices <- cppFilterWithWavelengthWindows(
      wavelengthsMatrix = wavelengthsMatrix,
      spectrumsMatrix = spectrumsMatrix,
      errorsMatrix = errorsMatrix,
      sizesVectorR = sizes,
      continuumWindowsVectorR = feWindows
    )
    wavelengthsMatrix <- filteredMatrices$wavelengthsMatrix
    spectrumsMatrix <- filteredMatrices$spectrumsMatrix
    errorsMatrix <- filteredMatrices$errorsMatrix
    filteredMatrices <- cppFilterWithParagon(
      spectrumsMatrix,
      continuum,
      expandedTemplate,
      sizes
    )
    continuum <- filteredMatrices$aMatrixFiltered
    expandedTemplate <- filteredMatrices$bMatrixFiltered
  }
  
  if(feFitParams$fitType == 'fwin') {
    filteredMatrices <- cppFilterWithValues(
      wavelengthMatrix = expandedTemplate,
      spectrumMatrix = wavelengthsMatrix,
      errorMatrix = errorsMatrix,
      sizes = sizes
    )
    expandedTemplate <- filteredMatrices$wavelengthsMatrix
    wavelengthsMatrix <- filteredMatrices$spectrumsMatrix
    errorsMatrix <- filteredMatrices$errorsMatrix
    
    filteredMatrices <- cppFilterWithParagon(
      paragonMatrix = expandedTemplate,
      aMatrix = continuum,
      bMatrix = spectrumsMatrix,
      sizes
    )
    continuum <- filteredMatrices$aMatrixFiltered
    spectrumsMatrix <- filteredMatrices$bMatrixFiltered
    
  }
  
  spectrumsMatrixFiltered <- spectrumsMatrix;
  wavelengthsMatrixFiltered <- wavelengthsMatrix;
  errorsMatrixFiltered <- errorsMatrix;
  continuumsMatrixFiltered <- continuum;
  expandedFeMatrixFiltered <- expandedTemplate;
  sizesFe <- sizes
  maxSize <- ASTRO_OBJ_SIZE
  if(feFitParams$fitType == 'fwin') {
    sizesFe <- cppCountNInfSize(spectrumsMatrixFiltered)
    maxSize <- max(sizesFe)
    
    spectrumsMatrixFiltered <- cppCopyNInf(spectrumsMatrix, maxSize)
    wavelengthsMatrixFiltered <- cppCopyNInf(wavelengthsMatrix, maxSize)
    errorsMatrixFiltered <- cppCopyNInf(errorsMatrix, maxSize)
    continuumsMatrixFiltered <- cppCopyNInf(continuum, maxSize)
    expandedFeMatrixFiltered <- cppCopyNInf(expandedTemplate, maxSize)
  }
  
  # scale template 
  feScaleRates <- cppCalculateFeScaleRates(
    spectrumsMatrix = spectrumsMatrixFiltered,
    templateMatrix = expandedFeMatrixFiltered,
    sizes = sizesFe,
    feFitParams = feFitParams
  )
  
  expandedTemplateC<- cppMatrixMultiplyCol(expandedTemplateC, feScaleRates)
  expandedFeMatrixFiltered <- cppMatrixMultiplyCol(expandedFeMatrixFiltered, feScaleRates)
  
  # calculate fit params
  reducedChisqFiltered <- cppCalculateFeReducesChisq(
    spectrumsMatrix = spectrumsMatrixFiltered,
    templateMatrix = expandedFeMatrixFiltered,
    errorsMatrix = errorsMatrixFiltered,
    sizes = sizesFe  
  )
  
  reducedChisqFull <- cppCalculateFeReducesChisq(
    spectrumsMatrix = spectrum$spectrums,
    templateMatrix = expandedTemplateC,
    errorsMatrix = spectrum$errors,
    sizes = sizes
  )
  equiwalentWidth <- cppCalculateTrapz(
    x = spectrum$wavelengths,
    y = cppMatrixDivideMatrix(expandedTemplateC, continuumC),
    sizes = sizes
  )
  
  #TODO: add fitting within range
  list(
    feFitted = expandedTemplateC,
    fitParams = list(
      scaleRates = feScaleRates,
      sizesFeWin = sizesFe,
      reducedChisqWindows = reducedChisqFiltered,
      reducedChisqFull = reducedChisqFull,
      equiwalentWidthFull = equiwalentWidth
    )
    #within range
  )
}

elementFitting <- function(
  wavelengthsMatrix, 
  spectrumsMatrix, 
  continuumMatrix,
  errorsMatrix,
  sizesVector,
  element
) {
  filteredMatrices <- cppFilterWithWavelengthWindows(
    wavelengthsMatrix = wavelengthsMatrix,
    spectrumsMatrix = spectrumsMatrix,
    errorsMatrix = continuumMatrix,
    sizesVectorR = sizesVector,
    continuumWindowsVectorR = list(c(element[['left']], element[['right']]))
  )
  sizes <- cppCountNInfSize(filteredMatrices$spectrumsMatrix)
  maxSize <- max(sizes)
  #TODO: return fit repeated
  if(maxSize == 0) 
    return(list(fit=F))
  spectrumsFiltered <- cppCopyNInf(filteredMatrices$spectrumsMatrix, maxSize)
  wavelengthsFiltered <- cppCopyNInf(filteredMatrices$wavelengthsMatrix, maxSize)
  continuumFiltered <- cppCopyNInf(filteredMatrices$errorsMatrix, maxSize)
  
  fitInitialResults <- rep(list(c(element[['a']], element[['b']], element[['c']], 0)), nrow(wavelengthsMatrix))
  fitResults <- cppFitGaussian(
    x = wavelengthsFiltered,
    y = spectrumsFiltered,
    sizes = sizes,
    results = fitInitialResults
  )
  gaussianMatrix <- cppCalculateGaussian(
    wavelengthsFiltered,
    fitResults,
    sizes
  )
  gaussianMatrix <- cppMatrixDivideMatrix(gaussianMatrix, continuumFiltered)
  gaussianChisq <- cppCalculateGaussianChisq(
    wavelengthsMatrix = wavelengthsMatrix,
    spectrumsMatrix = spectrumsMatrix,
    errorsMatrix = errorsMatrix,
    fitGaussianParams = fitResults,
    sizes = sizesVector
  )
  gaussianEquiwalentWidth <- cppCalculateTrapz(
    x = wavelengthsFiltered,
    y = gaussianMatrix,
    sizes = sizes
  )
  gaussianFWHM <- cppCalculateFWHM(fitResults)
  list(
    fitResults = fitResults,
    equiwalentWidth = gaussianEquiwalentWidth,
    chisq = gaussianChisq,
    fwhm  = gaussianFWHM,
    fit=T
  )
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
  
  elements <- loadDefaultSpectralLines()
  
  spectrumsC <- spectrums
  iter <- 3
  for(i in seq(1:iter)) {
    continuumResult <- continuumFitting(
      wavelengths,
      spectrums,
      errors,
      sizesVector,
      continuumParams
    )
    spectrums <- spectrumsC
    feFitResults <- feFitting(
      list(
        wavelengths=wavelengths,
        spectrums=cppMatrixMinusMatrix(spectrumsC, continuumResult$cfun),
        errors=errors,
        sizes=sizesVector
        ),
      continuum = continuumResult$cfun,
      feTemplate,
      feWindows,
      feFitParams = feFitParams
    )
    if(i < iter) {
      spectrums <- cppMatrixMinusMatrix(spectrums, feFitResults$feFitted)
    }
  }
  allFitParams <- append(list(continuum = continuumResult$fitParams), list(feFit=feFitResults$fitParams))
  spectrums <- cppMatrixMinusMatrix(spectrumsC, feFitResults$feFitted)
  spectrums <- cppMatrixMinusMatrix(spectrums, continuumResult$cfun)
  elementsFitParams <- list()
  for(i in 1:nrow(elements)) {
    element <- elements[i,]
    fitParams <- elementFitting(
      wavelengthsMatrix = wavelengths,
      spectrumsMatrix = spectrums,
      continuumMatrix = continuumResult$cfun,
      errorsMatrix = errors,
      sizesVector = sizesVector,
      element = element
    )
    if(fitParams$fit) {
        fitParams$fit <- NULL
        elementsFitParams <- append(elementsFitParams, list(list(element=element$name, params=fitParams)))
      }
    }
    
  allFitParams <- append(allFitParams, list(elements = elementsFitParams))
  
  filenames <- readLines(quas_list)
  saveResults(allFitParams, filenames, 'results')
}




t <- system.time({testFitting()})
print(t)
