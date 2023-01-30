#This function quantifies the number of significantly enriched GO terms
#by AD biological domain

bd.tally <- function( enrVct, biodomDefTbl){
  bdt <- bind_cols(
    domain = unique(biodomDefTbl$Biodomain),
    n_term = map_dbl( unique(biodomDefTbl$Biodomain),
                      ~ biodomDefTbl %>% filter(Biodomain == .x) %>% 
                        dplyr::select(GOterm_Name) %>% distinct() %>% nrow()),
    n_sig_term = map_dbl( unique(biodomDefTbl$Biodomain),
                          ~ enrVct[ enrVct %in% 
                                      biodomDefTbl$GOterm_Name[
                                        biodomDefTbl$Biodomain == .x]] %>% 
                            length()) ) %>% 
    bind_rows(
      ., 
      tibble(domain = 'none', 
             n_term = setdiff(enrVct, biodomDefTbl$GOterm_Name) %>% length(), 
             n_sig_term = setdiff(enrVct, biodomDefTbl$GOterm_Name) %>% length())
      ) %>% 
    mutate(domain = fct_reorder(domain, n_sig_term, .desc = F)) %>% 
    arrange(domain)
  bdt$prop <- bdt$n_sig_term / bdt$n_term
  return(bdt)
}
