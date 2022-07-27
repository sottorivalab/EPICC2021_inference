on_alma = TRUE

# load purity data
dfiles = list.files("input_data/purity_data/", "[.]rds$", full.names = TRUE)
for (i in seq_along(dfiles)) {
  assign(gsub("[.]rds$", "", basename(dfiles[i])), readRDS(dfiles[i]))
}

# load functions
null = list.files("functions", full.names = TRUE) %>% lapply(source)

# find input files:
result_dirs = "inference_results"
result_dirs = list.files(result_dirs, full.names = TRUE, include.dirs = TRUE)
result_dirs = result_dirs[grepl("C[0-9]{3}", basename(result_dirs))]
result_dirs = list.files(result_dirs, full.names = TRUE, include.dirs = TRUE)
result_dirs = result_dirs[file.exists(file.path(result_dirs, "target.rds"))]

