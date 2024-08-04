#' @title automatically construct and analyze an spmm 
#' 
#' @description This function automates construction and analysis of a spatial 
#' population matrix model, from construction of the vec-permutation matrix to
#' plotting. The user provides a `path` and `filename` for a spmm formatted
#' Excel Workbook (.xlsx). 
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
#' @param ... optional; additional arguments passed to `spmm.project`
#' 
#' @note Workbook should follow template found in extdata. Objects can be 
#' assigned from the list using the unpacking / destructuring assignment "%<-%" 
#' from the zeallot package. (metapopbio currently lists zeallot as a dependency).
#' 
#' @references
#' Schauberger P, Walker A (2022). openxlsx: Read, Write and Edit xlsx Files. 
#' https://ycphs.github.io/openxlsx/index.html, https://github.com/ycphs/openxlsx. 
#' 
#' @export
auto.spmm <- function(path, filename, plot = FALSE, ylabs = NA, xlabs = NA, ... ) {
  # source("C:/Users/caldridge/Documents/R/meta-pop-bio/R/auto.spmm.internal.R", local = TRUE)
  
  sheetNames <- openxlsx::getSheetNames(paste0(path, filename))
  patch_idx <- grep("patch", sheetNames)
  stage_idx <- grep("stage", sheetNames)
  
  tmp <- spmm.readxl(path, filename)
  
  c(n_stages, n_patches, group_by, lh_order, n_timesteps, stage_names, patch_names, n, matrices) %<-% tmp[1:9]
  
  P <-
    metapopbio::vec.perm(n_stages = n_stages,
                         n_patches = n_patches,
                         group_by = group_by)
  
  for (i in patch_idx) {
    assign(sheetNames[i], matrices[sheetNames[i]])
    assign(sheetNames[i], as.matrix(get(sheetNames[i])[[1]]))
  }
  
  for (i in stage_idx) {
    assign(sheetNames[i], matrices[sheetNames[i]])
    assign(sheetNames[i], as.matrix(get(sheetNames[i])[[1]]))
  }
  
  # if ()
  BB <- metapopbio::blk.diag(mget(sheetNames[patch_idx]))
  
  MM <- metapopbio::blk.diag(mget(sheetNames[stage_idx]))
  
  A <-
    metapopbio::spmm.project.matrix(
      P = P,
      BB = BB,
      MM = MM,
      group_by = group_by,
      lh_order = lh_order
    )
  
  if (!exists("ddf")) {
    ddf <- NA
  } 
  if (!exists("H")) {
    H <- NA 
  }
  if (!exists("D")) {
    D <- NA
  }
  
  projs <-
    metapopbio::spmm.project(
      n = n,
      A = A,
      n_timesteps = n_timesteps,
      n_stages = n_stages,
      n_patches = n_patches,
      ddf = ddf,
      H = H,
      D = D,
      P = P,
      BB = BB,
      MM= MM
    )
  
  if (plot == TRUE) {
    p <- metapopbio::spmm.plot(
      projections = projs,
      ylabs = ifelse(is.na(ylabs), "Value", ylabs),
      xlabs = ifelse(is.na(xlabs), "Time", xlabs),
      stage_names = stage_names,
      patch_names = patch_names
    )
    l <- list(tmp, P, BB, MM, A, projs, p)
  } else {
    l <- list(tmp, P, BB, MM, A, projs)
  }
  return(l)
}

