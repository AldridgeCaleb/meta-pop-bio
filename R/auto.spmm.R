#' @title automatically construct and analyze an spmm 
#' 
#' @description This function automates construction and analysis of a spatial 
#' population matrix model, from construction of the vec-permutation matrix to
#' plotting. THe user provides a `path` and `filename` for a spmm formatted
#' Excel Workbook (xlsx). 
#' 
#' Future plans include a "report" option that would include an RMarkdown (.Rmd)
#' and Word Document (.docx) as outputs. Additionally, I am exploring the ability
#' to "update" a spmm as a convenience function. 
#' 
#' @param path path to xlsx workbook
#' @param filename name of xlsx workbook
#' @param plot optional; output plots from `auto.spmm`?
#' @param ylabs optional; labels for y-axis labels
#' @param xlabs optional; labels for x-axis labels
#' 
#' @note Workbook should follow template found in extdata.
#' 
#' @references
#' Schauberger P, Walker A (2022). openxlsx: Read, Write and Edit xlsx Files. 
#' https://ycphs.github.io/openxlsx/index.html, https://github.com/ycphs/openxlsx. 
#' 
#' @export
auto.spmm <- function(path, filename, plot = FALSE, ylabs = NA, xlabs = NA) {
  source(system.file("R/auto.spmm.internal.R", package = "metapopbio"))
}

