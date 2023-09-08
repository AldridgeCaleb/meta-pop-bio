#' internal / backend script to read spreadsheets formatted for spmm
#' @param path path to xlsx workbook
#' @param filename name of xlsx workbook
#' @noRd
sheetNames <- getSheetNames(paste0(path, filename))
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
      rows = c(6:250),
      cols = c(1),
      colNames = TRUE
    )[, 1]
  )
patch_names <-
  as.vector(
    openxlsx::read.xlsx(
      xlsxFile = paste0(path, filename),
      sheet = "metadata",
      rows = c(6:250),
      cols = c(2),
      colNames = TRUE
    )[, 1]
  )
n <-
  as.vector(
    openxlsx::read.xlsx(
      xlsxFile = paste0(path, filename),
      sheet = "metadata",
      rows = c(6:250),
      cols = c(3),
      colNames = TRUE
    )[, 1]
  )

# demographics matrices
stage_idx <- grep("stage", sheetNames)
for (i in stage_idx) {
  as.matrix(assign(sheetNames[i],
                   openxlsx::read.xlsx(
                     xlsxFile = paste0(path, filename),
                     sheet = i,
                     colNames = FALSE
                   )))
}

# movement matrices
patch_idx <- grep("patch", sheetNames)
for (i in stage_idx) {
  assign(sheetNames[i],
         as.matrix(openxlsx::read.xlsx(
           xlsxFile = paste0(path, filename),
           sheet = i,
           colNames = FALSE
         )))
}