#BiocManager::install("GenomicFeatures")
#browseVignettes("GenomicFeatures")

# install.packages("remotes")
# remotes::install_github("asrinivasan-oa/gread")
# BiocManager::install("GenomicRanges")
# 
# library(gread)
# path <- system.file("tests", package="gread")
# gtf_file <- file.path(path, "sample.gtf")
# gtf <- read_format(gtf_file)
# 
# # update gtf with intron coordinats from the exon gaps 
# # and return the updated object
# ans <- construct_introns(gtf, update=TRUE)[] # default
# # same as above, but returns only constructed introns
# introns <- construct_introns(gtf, update=FALSE)
