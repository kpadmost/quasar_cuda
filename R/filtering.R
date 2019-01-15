rFilterInfs <- function(inputMatrix, newSize=ASTRO_OBJ_SIZE) {
  cppCopyNInf(inputMatrix, newSize)
}