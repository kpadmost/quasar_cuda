

loadWindows <- function(filename) {
  raw <- readLines(filename)
  raw <- raw[sapply(raw, nchar) > 0]
  lapply(raw, function(x) {
    slice <- strsplit(x, ' ')[[1]]
    c(as.numeric(slice[1]), as.numeric(slice[2]))
  })
}

#TODO: remove!
#' @export
qP <- function() {
  d1 <- list(c(3.57840, 0.00010, 1.80677, 0.00000))
  d2 <- list(c(3.57900, 0.00010, 1.51807, 0.00000))
  dm <- list(rep(0, 4))
  lmx <- rep(dm, 30)
  c(d1, d2, lmx)
}

#' @export
loadFeTemplate <- function(filename='iron_emission_temp.bin') {
  template <- read.delim(filename, header = F, sep = " ")
  lambda <- as.double(template[, 1])
  flux <- as.double(template[, 3])
  # transposing template
  return(list(lambda=lambda, flux=flux))
}

#' @export
loadContinuumWindows <- function(file='contwind') {
  loadWindows(file)
}

#' @export
loadFeWindows <- function(file='iron_emission_windows') {
  loadWindows(file)
}

#' @export
loadElements <- function(file='spectral_lines') {
  read.table(file,
             header= F,
             col.names = c('name', 'left', 'right', 'a', 'b', 'c'),
             colClasses = c("character", "numeric", "numeric","numeric","numeric","numeric"))
  
}

#' @export
loadFeFitParams <- function() {
  list(
    fwhmn = 1600.0,
    fwhmt = 900.0,
    feFitRange = c(2200.0, 2650.0),
    fitType="fwin"
  )
}