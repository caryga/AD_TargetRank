# this function 
#	(1) calculates GSEA for all GO terms,
# (2) annotates results with biodomain-specific information

bd.enr <- function(gene_list, key_type = 'SYMBOL', score_type = 'pos'){
    
    # run GSEA analysis
    enr <- clusterProfiler::gseGO(
      geneList = gene_list, 
      ont = 'all',
      OrgDb = org.Hs.eg.db::org.Hs.eg.db, 
      keyType = key_type,
      # minGSSize = 10,
      # maxGSSize = 500,
      eps = 0,
      pvalueCutoff = 1,
      scoreType = score_type
      )
    
    enr <- enr@result %>% 
      rename(pathway = Description, 
             pval = pvalue, 
             padj = p.adjust, 
             ES = enrichmentScore, 
             size = setSize, 
             leadingEdge = core_enrichment) %>% 
      left_join(., biodom %>% select(pathway=GOterm_Name, Biodomain), by='pathway') %>% 
      full_join(., dom.cols, by = c('Biodomain'='domain')) %>% 
      relocate(Biodomain, .after = pathway) %>% 
      mutate(
        Biodomain = case_when(is.na(Biodomain)~'none', T~ Biodomain),
        color = case_when(is.na(color)~'#7f7f7f', T~color)
        # , Biodomain = fct_relevel(Biodomain, as.character(bdt$domain) )
        ) %>% 
      arrange(Biodomain) %>% 
      rowwise() %>% mutate(leadingEdge = str_split(leadingEdge, '/')) 
    
  return(enr)
  
}

