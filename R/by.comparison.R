by.comparison <- function(data,sm){
  #Binding variables to function
  treatments <- comparisons <- V1 <-V2 <-
    trt.e <-trt.c <- NULL
  if (sm %in%c("MD", "SMD" )){
    data %<>% select(study, trt, outcome, n)
    data.st <- data %>%select( study, trt)
    data.st %<>% nest(treatments=c(trt))
    data.st %<>%
      mutate(comparisons = map(treatments, function(x) combn(x[[1]], 2) %>% t() %>% `colnames<-`(c("V1", "V2")) %>% as_tibble)) %>%
      select(-treatments) %>%
      unnest(cols = c(comparisons)) %>%
      rename(trt.e=V1, trt.c=V2) %>%
      left_join(data %>% select(study, outcome, n, trt), by = c("study", "trt.e" = "trt")) %>%
      left_join(data %>% select(study, outcome, n, trt), by = c("study", "trt.c" = "trt"), suffix=c(".e",".c"))%>%
      mutate(comparison = ifelse(trt.e < trt.c,
                                 paste(trt.e, trt.c, sep = " vs. "),
                                 paste(trt.c, trt.e, sep = " vs. ")))

  } else if (sm %in%c("OR", "RR" )){
    data %<>% select(study, trt, outcome, n)
    data.st <- data %>%select( study, trt)
    data.st %<>% nest(treatments=c(trt))
    data.st %<>%
      mutate(comparisons = map(treatments, function(x) combn(x[[1]], 2) %>% t() %>% `colnames<-`(c("V1", "V2")) %>% as_tibble)) %>%
      select(-treatments) %>%
      unnest(cols = c(comparisons)) %>%
      rename(trt.e=V1, trt.c=V2) %>%
      left_join(data %>% select(study, outcome, n, trt), by = c("study", "trt.e" = "trt")) %>%
      left_join(data %>% select(study, outcome, n, trt), by = c("study", "trt.c" = "trt"), suffix=c(".e",".c"))%>%
      mutate(comparison = ifelse(trt.e < trt.c,
                                 paste(trt.e, trt.c, sep = " vs. "),
                                 paste(trt.c, trt.e, sep = " vs. ")))

  }
  return(data.st)
}
