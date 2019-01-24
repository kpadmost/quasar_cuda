#' @importFrom parseFITS readQuasarInfo

BLOCK_SIZE <- 512
ASTRO_OBJ_SIZE <- 4096

#' @export
loadQuasarsFromFolder <- function(folder, datafile) {
  quasarFiles <- readLines(datafile)
  files <- file.path(folder, quasarFiles)
  lapply(files, readQuasarInfo)
}

getMatrixFromVector <- function(quasars, values, complemented = T) {
  data <- lapply(lapply(quasars,  `[[`, values), 
                      function(vec) { c(vec, rep(Inf, ASTRO_OBJ_SIZE - length(vec))) })
  dataM <- do.call(rbind, data)
  if(complemented)
    dataM <- complementMatrix(dataM)
  dataM
}

#' @export
getFluxMatrix <- function(quasars, complemented=T) 
{
  getMatrixFromVector(quasars, 'flux', complemented)
}

#' @export
getSizesVector <- function(quasars, complemented=T) {
  data <- sapply(lapply(quasars,  `[[`, 'flux'), length)
  if(complemented)
    data <- complementVector(data)
  data
}


#' @export
getErrorMatrix <- function(quasars, complemented=T) 
{
  getMatrixFromVector(quasars, 'error', complemented)
}


#' @export
getLambdaParams <- function(quasars, complemented=T) {
  params <- lapply(
    lapply(
      quasars,
      `[[`,
      'specs'
    ), function(quasarParams) c(quasarParams$a, quasarParams$b, quasarParams$z, 0)
  )
  if(complemented) {
    comSize <- BLOCK_SIZE - (length(params) %% BLOCK_SIZE)
    if(comSize == BLOCK_SIZE) return(params)
    dm <- list(rep(0, 4))
    params <- c(params, rep(dm, comSize))
  }
  params
}

complementVector <- function(sVector) {
  size <- length(sVector)
  complementedSize <- BLOCK_SIZE  - (size %% BLOCK_SIZE)
  c(sVector, rep(ASTRO_OBJ_SIZE, complementedSize))
}

complementMatrix <- function(sMatrix) {
  rows <- nrow(sMatrix)
  complementedRSize <- BLOCK_SIZE  - (rows %% BLOCK_SIZE)
  if(complementedRSize == BLOCK_SIZE) {
    return(sMatrix)
  }
  complementedMatrix <- matrix(0, nrow = complementedRSize, ncol = ASTRO_OBJ_SIZE)
  rbind(sMatrix, complementedMatrix)
}

complementQuasarMatrix <- function(quasarMatrix) {
  lapply(quasarMatrix, complementMatrix)
}