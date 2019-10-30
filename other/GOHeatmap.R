dirpath = "~/Dropbox/CompRepos/sc_R_sparse/Oct2017/"
source(paste0(dirpath,"sc.R"))
source(paste0(dirpath, "scFunctions.R"))
source(paste0(dirpath,"scR_DE.R"))
source(paste0(dirpath,"scR_plotting.R"))
#source(paste0(dirpath,"load_SIMLR.R"))
library(dplyr)

# Create a Heatmap of GO terms

time_points = c("ctrl","12h","1d","2d","4d","1w","2w")

genes_temporal = list()
for (tt in time_points){
  genes_temporal[[tt]] = read.csv(paste0(tt,"_microglia.csv"), header = FALSE, stringsAsFactors = FALSE)$V1
  
}

bg_genes = unique(unlist(genes_temporal))

# GO analysis
topGO_res = list()
for (tt in time_points){
  debugonce("topGOterms")
  topGO_res[[tt]] = topGOterms(genes_temporal[[tt]], bg.genes = bg_genes, topnodes.print = 1620)
}

GO_terms_to_plot = c(); GO_term_names = c()
for (tt in time_points){
  topGO_res[[tt]]$res.table$fdr = p.adjust(as.numeric(topGO_res[["ctrl"]]$res.table$pval), method="fdr")
  GO_terms_temp  = subset(topGO_res[[tt]]$res.table, pval < 0.01)$GO.ID
  GO_terms_to_plot = unique(c(GO_terms_to_plot, GO_terms_temp))
  GO_terms_temp  = subset(topGO_res[[tt]]$res.table, pval < 0.01)$Term
  GO_term_names = unique(c(GO_term_names, GO_terms_temp))
}

PvalMat = matrix(NA, ncol = 7, nrow = length(GO_terms_to_plot) )
colnames(PvalMat) = time_points; rownames(PvalMat) = GO_terms_to_plot
for (tt in time_points){
  GOtemp = topGO_res[[tt]]$res.table; rownames(GOtemp) = GOtemp$GO.ID
  PvalMat[,tt] = as.numeric(GOtemp[GO_terms_to_plot, "pval"]) # change to fdr
}

PvalMat = -log(PvalMat)
rownames(PvalMat) = paste0(GO_terms_to_plot,":",GO_term_names)

pdf("GOHeatmap_microglia.pdf",w=10,h=10, useDingbats = FALSE)
p=DataHeatmap(data=t(PvalMat), return.plot = TRUE, cex.row = 8, col.low = "blue", col.mid = "white",col.high = "red", quantile.thresh = c(0.01,0.99) )
print(p)
dev.off()

rows.use = apply(PvalMat,1, max) > 0.1
p=DataHeatmap(data=t(PvalMat[rows.use,]), return.plot = TRUE, cex.row = 8, col.low = "blue", col.mid = "white",col.high = "red", quantile.thresh = c(0.01,0.99) )
