#' @title read spreadsheets from Excel Workbook (xlsx) formatted for spmm
#' 
#' @description Read in spmm worksheet data
#' 
#' @param path path to xlsx workbook
#' @param filename name of xlsx workbook
#' 
#' @note Workbook should follow template found in extdata. (**ADD TEMPLATE**)
#' 
#' @references
#' Schauberger P, Walker A (2022). openxlsx: Read, Write and Edit xlsx Files. 
#' https://ycphs.github.io/openxlsx/index.html, https://github.com/ycphs/openxlsx. 
#' 
#' @export
spmm.readxl <- function(path, filename) {
  source(system.file("R/spmm.readxl.internal.R", package = "metapopbio"))
}

