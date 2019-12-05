#' Takes a topGO object and a summary table, and goes down from the top to find the N most signicant non-go terms that have a less than specified overlap
topNonRedundantTerms = function(GO_object, 
                                res_table, 
                                topN = 6, 
                                avoid_terms = NULL, 
                                maxN = 200,
                                minN = 5,
                                pval_thresh = 1e-2,
                                overlap = 0.3){
 
  
  gene_terms_all = genesInTerm(GO_object)
  GO_terms_all = c()
  countn = 0
  
  for (idx in c(1:nrow(res_table))){
    GO_tmp = res_table[idx,"GO.ID"]
    Nannot = res_table[idx,"Annotated"]
    pval = as.numeric(res_table[idx,"pval"])
    if (is.na(pval)) pval = 1e-30
    if (pval > pval_thresh) break
    
    if ((Nannot >= maxN) | (Nannot <= minN) ) next
    if ((GO_tmp %in% avoid_terms)) next
    
    if (length(GO_terms_all) == 0){
      GO_terms_all = c(GO_terms_all, GO_tmp)
    } else {
      
      genes_in_term0 = gene_terms_all[[GO_tmp]]
      
      flag_tmp = 0
      for (gterm in GO_terms_all){
        genes_in_term1 = gene_terms_all[[gterm]]
        
        if ((length(intersect(genes_in_term0, genes_in_term1)) / min(length(genes_in_term0),length(genes_in_term1))) > overlap ){
          flag_tmp=1
        }
          
      }
      if (flag_tmp == 0) GO_terms_all = c(GO_terms_all, GO_tmp)
      
    }
    
    if (length(GO_terms_all) >= topN) break
    
  }
  
  return(GO_terms_all)
  
}
