# Save as h5
library(DropletUtils)

Counts_E13toP0 = readRDS("../Objects/Counts_E13toP0.rds")
Counts_P5 = readRDS("../Objects/Counts_P5.rds")
Counts_P5 = Counts_P5$All

# Filter Count matrix
genes.common = intersect(rownames(Counts_E13toP0), rownames(Counts_P5))
Count.mat = cbind(Counts_E13toP0[genes.common,], Counts_P5[genes.common,])
ngenes = colSums(Count.mat > 0)
Count.mat = Count.mat[,ngenes > 700]


write10xCounts("RGC_dev_full.h5", Count.mat, barcodes = colnames(Count.mat), gene.id = rownames(Count.mat),
               gene.symbol = rownames(Count.mat), type = "HDF5", genome="mm10")
