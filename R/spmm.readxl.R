#' @title read spreadsheets from Excel Workbook (xlsx) formatted for spmm
#'
#' @description Read in spmm worksheet data
#'
#' @param path path to xlsx workbook
#' @param filename name of xlsx workbook
#'
#' @note Workbook should follow template found in extdata.
#' 
#' @returns Returns a list. Objects can be assigned from the list using the 
#' unpacking / destructuring assignment "%<-%" from the zeallot package.
#' (metapopbio currently lists zeallot as a dependency).
#'
#' @references
#' Schauberger P, Walker A (2022). openxlsx: Read, Write and Edit xlsx Files.
#' https://ycphs.github.io/openxlsx/index.html, https://github.com/ycphs/openxlsx.
#'
#' @export
spmm.readxl <- function(path, filename) {
  sheetNames <- openxlsx::getSheetNames(paste0(path, filename))
  # metadata
  # metadata part 1
  metadata <-
    openxlsx::read.xlsx(
      xlsxFile = paste0(path, filename),
      sheet = "metadata",
      rows = c(1:5),
      cols = c(1:2),
      colNames = FALSE
    )
  for (i in 1:nrow(metadata)) {
    assign(metadata[i, 1], metadata[i, 2])
  }
  n_stages <- as.numeric(n_stages)
  n_patches <- as.numeric(n_patches)
  n_timesteps <- as.numeric(n_timesteps)
  # metadata part 2
  stage_names <-
    as.vector(
      openxlsx::read.xlsx(
        xlsxFile = paste0(path, filename),
        sheet = "metadata",
        rows = c(7:250),
        cols = c(1),
        colNames = TRUE
      )[, 1])
  patch_names <-
    as.vector(
      openxlsx::read.xlsx(
        xlsxFile = paste0(path, filename),
        sheet = "metadata",
        rows = c(7:10000),
        cols = c(2),
        colNames = TRUE
      )[, 1])
  n <- as.vector(t(
    rev(openxlsx::read.xlsx(
      xlsxFile = paste0(path, filename),
      sheet = "metadata",
      rows = c(7:10000),
      cols = c(3:252),
      colNames = TRUE
    ))))
  comment(n) <- group_by
  
  # movement matrices
  stage_idx <- grep("stage", sheetNames)
  for (i in stage_idx) {
    as.matrix(assign(
      sheetNames[i],
      openxlsx::read.xlsx(
        xlsxFile = paste0(path, filename),
        sheet = i,
        colNames = FALSE
      )
    ))
  }
  # mget(sheetNames[stage_idx])

  # demographic matrices
  patch_idx <- grep("patch", sheetNames)
  for (i in patch_idx) {
    assign(sheetNames[i],
           as.matrix(
             openxlsx::read.xlsx(
               xlsxFile = paste0(path, filename),
               sheet = i,
               colNames = FALSE
             )
           ))
  }
  
  B_list <- lapply(mget(sheetNames[patch_idx]), function(x) as.matrix(x))
  BB <- blk.diag(B_list)
  M_list <- lapply(mget(sheetNames[stage_idx]), function(x) as.matrix(x))
  MM <- blk.diag(M_list)
  
  return(
    list(
      n_stages = n_stages,
      n_patches = n_patches,
      group_by = group_by,
      lh_order = lh_order,
      n_timesteps = n_timesteps,
      stage_names = stage_names,
      patch_names = patch_names,
      n = n,
      mget(sheetNames[-1]),
      MM = MM,
      BB = BB
    )
  )
}
