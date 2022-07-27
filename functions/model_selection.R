format_model_string = function(x, short=TRUE) {
  x = as.character(x)
  
  if (length(x) > 1) 
    return(sapply(x, format_model_string))
  
  x = gsub("no_death_", "", x) %>% 
    gsub(pattern = "selection_(?=2)", replacement = "selection x ", perl = TRUE) %>% 
    gsub(pattern = "dispersion", replacement = "D") %>% 
    gsub(pattern = "regionsize", replacement = "RS") %>% 
    strsplit("_") %>%
    unlist() %>% Hmisc::capitalize() %>%
    paste0(collapse = " + ") 
  
  if (short)  {
    x = gsub(" [+] Push", "", x)
    x = gsub(" [+] D( |$)", " ", x)
    x = gsub(" [+] RS( |$)", " ", x)
    x = gsub(" $", "", x)
  }
  
  x
}
