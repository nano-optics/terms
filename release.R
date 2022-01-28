library(here)
setwd(here('.'))
getwd()
library(zip)
library(fs)

fs::dir_ls(recurse = 0)

release_name <- 'terms_code'
library(desc)
(release_version <- desc::desc_get('Version', 'DESCRIPTION'))

rmarkdown::render('vignettes/Keywords.Rmd', output_format = 'html_document')
rmarkdown::render('vignettes/Keywords.Rmd', output_format = 'pdf_document')

build <- fs::dir_ls('build', recurse = 0)
test <- fs::dir_ls('test', recurse = 0)
src <- setdiff(fs::dir_ls('src', recurse = 0), 'src/Makefile') # exclude dummy Makefile

files <- c(build, 'CMakeLists.txt', src, test, 'README.md', 'vignettes/Keywords.html', 'vignettes/Keywords.pdf', 'publications.bib')

zip::zip(glue::glue('{release_name}_{release_version}.zip'), files)
utils::tar(glue::glue('{release_name}_{release_version}.tar'), files, compression = 'xz')
