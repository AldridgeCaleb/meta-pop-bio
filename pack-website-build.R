# Build package website
# Install released version from CRAN
install.packages("pkgdown")
# Run once to configure your package to use pkgdown
# usethis::use_pkgdown()

pkgdown::build_site(examples = FALSE,
                    preview = FALSE)

usethis::use_pkgdown_github_pages()
