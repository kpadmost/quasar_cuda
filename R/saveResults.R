#' @importFrom RJSONIO toJSON


#' @export
saveResults <- function(fitResults, objectList, savename) {
  results <- formatResults(fitResults, objectList)
  write(toJSON(results), paste(savename, '.json', sep=''))
}

formatResults <- function(fitResults, objectsList) {
  result <- list()
  for(i in 1:length(objectsList)) {
    result <- append(result, list(stobjectName = getResultsGlobal(fitResults, i)))
  }
  names(result) <- objectsList
  result
}

getResultsGlobal <- function(allFitParams, objectIndex) {
  result <- list(
    continuum = list(
      cReglin = allFitParams$continuum$cReglin[objectIndex],
      cReglin = allFitParams$continuum$reglin[objectIndex],
      chisq = allFitParams$continuum$chisq[objectIndex]
    ),
    feFit = list(
      scale = allFitParams$feFit$scaleRates[objectIndex],
      sizesWindows = allFitParams$feFit$sizesFeWin[objectIndex],
      reducedChisqWindows = allFitParams$feFit$reducedChisqWindows[objectIndex],
      reducedChisqFull = allFitParams$feFit$reducedChisqFull[objectIndex],
      equiwalentWidthFull = allFitParams$feFit$equiwalentWidthFull[objectIndex]
    )
  )
  elements <- list()
  elements <- lapply(allFitParams$elements, function(elm) {
    list(
      name = elm$element,
      fitResults = list(
        result = elm$params$fitResults[objectIndex],
        equiwalentWidth = elm$params$equiwalentWidth[objectIndex],
        chisq = elm$params$chisq[objectIndex],
        fwhm = elm$params$fwhm[objectIndex]
      )
    )
  })
  append(result, list(elements=elements))
}
