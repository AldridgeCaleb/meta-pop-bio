#' internal / backend script to automatically construct and analyze an spmm
#' @param path path to xlsx workbook
#' @param filename name of xlsx workbook
#' @noRd
spmm.readxl(path, filename)

P <-
    metapopbio::vec.perm(n_stages = n_stages,
                         n_patches = n_patches,
                         group_by = group_by)

BB <- metapopbio::blk.diag(get(mget(sheetNames[patch_idx])))

MM <- metapopbio::blk.diag(get(mget(sheetNames[stage_idx])))

A <-
  metapopbio::spmm.project.matrix(
    P = P,
    BB = BB,
    MM = MM,
    group_by = group_by,
    lh_order = lh_order
  )

projs <-
  metapopbio::spmm.project(
    n = n,
    A = A,
    n_timesteps = n_timesteps,
    n_stages = n_stages,
    n_patches = n_patches
  )

if (plot == TRUE) {
  metapopbio::spmm.plot(
    projections = projs,
    ylabs = if(!is.na(ylabs)){"Value"},
    xlabs = if(!is.na(xlabs)){"Time"},
    stage_names = stage_names,
    patch_names = patch_names
  )
}
