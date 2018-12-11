
require(FITSio)

UNIT <- 1e-17

extractFITSHeaderValue <- function(header, key) {
  valueIndex <- which(header == key)
  return(header[valueIndex + 1])
}

getHeaderInfo <- function(fitFile) {
  header <- fitFile$hdr
  a <- as.numeric(extractFITSHeaderValue(header, 'COEFF0'))
  b <- as.numeric(extractFITSHeaderValue(header, 'COEFF1'))
  mjd <- as.numeric(extractFITSHeaderValue(header, "MJD"))
  fiber <- as.numeric(extractFITSHeaderValue(header, "FIBERID"))
  plate <- as.numeric(extractFITSHeaderValue(header, "PLATEID"))
  z <- as.numeric(extractFITSHeaderValue(header, "Z"))
  ra <- as.numeric(extractFITSHeaderValue(header, "RAOBJ"))
  dec <- as.numeric(extractFITSHeaderValue(header, "DECOBJ"))
  type <- as.numeric(extractFITSHeaderValue(header, "OBJTYPE"))
  list(mjd=mjd, a=a, b=b, fiber=fiber, plate=plate, z=z, ra=ra, dec=dec, type=type)
}

getImageInfo <- function(fitFile) {
  noise <- fitFile$imDat[, 3] * UNIT
  flux <- fitFile$imDat[, 1] * UNIT
  bitmask <- fitFile$imDat[, 4]
  list(flux=flux, noise=noise, bitmask=bitmask)
}

# function to read single fits file
#' @export
readQuasarInfo <- function(filepath) {
  # since we need first header, is the best choise
  fitFile <- readFITS(filepath, hdu = 1)
  headerInfo <- getHeaderInfo(fitFile)
  imageInfo <- getImageInfo(fitFile)
  
  # restoring lambda from params
  #lambda <- restoreLambda(length(imageInfo$flux), as.double(headerInfo$a), as.double(headerInfo$b))
  #imageInfo$lambda <- lambda
  list(specs=headerInfo, flux=imageInfo$flux, error=imageInfo$error, bitmask=imageInfo$bitmask)
}