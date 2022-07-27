on_alma = TRUE

# load purity data
dfiles = list.files("input_data/purity_data/", "[.]rds$", full.names = TRUE)
for (i in seq_along(dfiles)) {
  assign(gsub("[.]rds$", "", basename(dfiles[i])), readRDS(dfiles[i]))
}

# load functions
null = list.files("functions", full.names = TRUE) %>% lapply(source)

# find input files:
result_source = "inference_results"
result_dirs = file.path(result_source, dirname(list.files(result_source, "target.rds", recursive = TRUE)))
result_dirs = result_dirs[grepl("C[0-9]{3}", basename(dirname(result_dirs)))]

