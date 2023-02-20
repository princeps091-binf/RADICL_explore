library(renv)

renv::init()
renv::install("tidyverse")
renv::install("dbplyr")
renv::install("duckdb")
renv::install("DBI")
renv::install("fs")
renv::install("vroom")
renv::install("igraph")
renv::install("dbscan")
renv::install("bioc::GenomicRanges")
renv::install("bioc::AnnotationHub")
renv::install("furrr")
renv::install("valr@0.6.3")
renv::install("remotes")


renv::snapshot()
