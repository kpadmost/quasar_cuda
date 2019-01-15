loadWindows <- function(filename) {
  raw <- readLines(filename)
  raw <- raw[sapply(raw, nchar) > 0]
  lapply(raw, function(x) {
    slice <- strsplit(x, ' ')[[1]]
    c(as.numeric(slice[1]), as.numeric(slice[2]))
  })
}

#' @export
loadDefaultFeTemplate <- function() {
  template <- system.file("data", "iron_emission_temp.bin", package = "QuasarRfit")
  loadFeTemplate(template)
}

loadFeTemplate <- function(filename='iron_emission_temp.bin') {
  template <- read.delim(filename, header = F, sep = " ")
  lambda <- as.double(template[, 1])
  flux <- as.double(template[, 3])
  # transposing template
  list(lambda=lambda, flux=flux)
}

#' @export
loadDefaultContinuumWindows <- function() {
  windows <- system.file("data", 'contwind', package = "QuasarRfit")
  loadContinuumWindows(windows)
}

loadContinuumWindows <- function(file) {
  loadWindows(file)
}

#' @export
loadDefaultFeWindows <- function() {
  windows <- system.file("data", 'iron_emission_windows', package = "QuasarRfit")
  loadFeWindows(windows)
}

loadFeWindows <- function(file) {
  loadWindows(file)
}

#' @export
loadDefaultSpectralLines <- function() {
  lines <- system.file("data", 'spectral_lines', package = "QuasarRfit")
  loadSpectralLines(lines)
}

loadSpectralLines <- function(file) {
  read.table(file,
             header= F,
             col.names = c('name', 'left', 'right', 'a', 'b', 'c'),
             colClasses = c("character", "numeric", "numeric","numeric","numeric","numeric"))
  
}

loadDefaultFeFitParams <- function() {
  list(
    fwhmn = 1600.0,
    fwhmt = 900.0,
    feFitRange = c(2200.0, 2650.0),
    isSubC = FALSE,
    fitType="fwin"
  )
}
