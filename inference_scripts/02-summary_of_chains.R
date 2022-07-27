library(dplyr)
library(cowplot)
library(ggplot2)
library(CHESS)
theme_set(theme_cowplot())
source("setup_environment.R")


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Options ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

result_dir = c("inference_results")

set_order = 
  paste0(
    rep(c("neutral", "selection", "selection_2"), each = 2), 
    "_",
    rep(c("no_death_push", "death_push"), 3)
  )

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Main script ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

sim_dirs = 
  result_dir %>% 
  list.files("target[.]rds", recursive = TRUE, full.names = TRUE) %>% 
  dirname()

smc_abc_summary = 
  data.frame(
    result_dir=dirname(dirname(sim_dirs)),
    case=factor(basename(dirname(sim_dirs))),
    set=basename(sim_dirs),
    set_formated=factor(format_model_string(basename(sim_dirs), short = TRUE), model_string_levels),
    option_set=basename(dirname(dirname(sim_dirs))),
    length=sapply(sim_dirs, get_chain_length), 
    finished=sapply(sim_dirs, smc_abc_finished),
    row.names = NULL
  ) %>% dplyr::filter(length>0) %>% 
  dplyr::mutate(set=factor(set, unique(c(set_order, set)), ordered = TRUE))


# save csv file to result dirs
for (c_result_dir in unique(smc_abc_summary$result_dir)) {
  smc_abc_summary %>% 
    dplyr::filter(result_dir == c_result_dir) %>%
    dplyr::select(-result_dir) %>%
    readr::write_csv(file.path(c_result_dir, "chain_summary.csv"))
}


plot_chain = function(d) {
  d %>% ggplot(aes(x=case, y=length, color=finished)) + 
    geom_point() + 
    facet_wrap(~set, ncol=1, scales = "free_y") + 
    theme(axis.text.x = element_text(angle=90, vjust = 0.5)) + 
    background_grid() + 
    labs(color="Finished") + 
    scale_color_brewer(palette = "Set1", direction = 1) + 
    xlab("Case") + 
    ylab("Chain length") + 
    ylim(0, NA) + 
    theme(legend.position = "bottom") + 
    scale_x_discrete(drop=FALSE)
}
 
# save csv file to result dirs
for (c_result_dir in unique(smc_abc_summary$result_dir)) {
  
  plot_smc_abc_summary = 
    smc_abc_summary %>% 
    dplyr::filter(result_dir == c_result_dir) %>%
    plot_chain()
  
  out_file = file.path(as.character(c_result_dir), "chain_summary.pdf")
  height = 10/9 * length(unique(plot_smc_abc_summary$data$set)) + 1.2
  ggsave(out_file, plot_smc_abc_summary, width=5, height=height)
  
}

plot_smc_abc_summary = smc_abc_summary %>% plot_chain()
height = 10/9 * length(unique(plot_smc_abc_summary$data$set)) + 1.2
out_file = file.path(result_dir, "chain_summary.pdf")
ggsave(out_file, plot_smc_abc_summary, width=5, height=height)


# save csv file to result dirs
for (c_result_dir in unique(dirname(smc_abc_summary$result_dir))) {
  
  plot_smc_abc_summary = 
    smc_abc_summary %>% 
    dplyr::filter(dirname(result_dir) == c_result_dir) %>%
    plot_chain()
  
  out_file = file.path(as.character(c_result_dir), "chain_summary.pdf")
  height = 10/9 * length(unique(plot_smc_abc_summary$data$set)) + 1.2
  ggsave(out_file, plot_smc_abc_summary, width=5, height=height)
  
}
