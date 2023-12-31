# Helpers building package
# 1. make changes and save; commit and push to GitHub
# 2. run roxygenise
roxygen2::roxygenise()  # builds NAMESPACE and man files
# 3. commit and push to GitHub (again)
# 4. download from GitHub (may have to remove from library prior)
# 5. run check
devtools::check()
# 6. build site locally
pkgdown::build_site(examples = TRUE,
                    lazy = FALSE,
                    preview = TRUE,
                    devel = FALSE)
# 7. publish site to gh-pages branch
usethis::use_pkgdown_github_pages()




# Preamble/one-run --------------------------------------------------------
# Build package website
# Install released version from CRAN
# install.packages("pkgdown")
# Run once to configure your package to use pkgdown
# usethis::use_pkgdown()
# usethis::use_logo("metapopbio.png")