# Build package website
# Install released version from CRAN
install.packages("pkgdown")
# Run once to configure your package to use pkgdown
# usethis::use_pkgdown()
usethis::use_logo("metapopbio.png")

pkgdown::build_site(examples = TRUE,
                    lazy = FALSE,
                    preview = TRUE,
                    devel = FALSE)

usethis::use_pkgdown_github_pages()
