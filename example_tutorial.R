source("~/Dropbox/CompRepos/sc_R_sparse/Oct2017/sc.R")
source("~/Dropbox/CompRepos/sc_R_sparse/Oct2017/scFunctions.R")
source("~/Dropbox/CompRepos/sc_R_sparse/Oct2017/scR_plotting.R")
source("~/Dropbox/CompRepos/sc_R_sparse/Oct2017/FastPCA.R")
source("~/Dropbox/CompRepos/sc_R_sparse/Oct2017/scR_DE.R")
source("~/Dropbox/CompRepos/sc_R_sparse/Oct2017/load_SIMLR.R")

# Load Old object
load("~/Dropbox/projects/MouseRGC/analysis/AdultRGC_10X/May2017/adultRGC_rgcs_only_cluster_170729_removeclust40.Rdata")

# Update to new version
new.object = UpdateObjectOct17(rgcs_only_recluster)
rgcs_only_recluster = new.object

rm(new.object)
rgcs_only_recluster = set.all.ident(rgcs_only_recluster, id="m") # Revert to cluster identities


# Let's find DE genes using the wilcoxon rank sum test
p=proc.time()
markers =find.markers(rgcs_only_recluster, ident.1 = 6, ident.2 = 13, min.perc = 0.3, test.use="wilcox")
print(proc.time() - p)

# Let's find DE genes using the wilcoxon rank sum test by restricting the max number of cells per cluster to 200
p=proc.time()
markers_downsamp =find.markers(rgcs_only_recluster, ident.1 = 6, ident.2 = 13, min.perc = 0.3, test.use="wilcox", max.cells.per.ident = 200)
print(proc.time() - p)

# We get almost identical results while saving 1/3rd of the time!!

# Find markers for two clusters and attach a table
markers_17_18 =find.all.markers(rgcs_only_recluster, idents.use = c(17,18), test.use = "t", add.to.slot = FALSE)

# Find markers for two clusters and add it to the new de slot in the object (convenient way to store DE info)
# Note that we perform DE only for two clusters idents.use = c(17,18) but not specifying this parameter, will perform the calculation for all clusters
rgcs_only_recluster =find.all.markers(rgcs_only_recluster, idents.use = c(17,18), test.use = "wilcox", add.to.slot = TRUE)

# Find markers and add to a different slot based on the edgeR QLF test
rgcs_only_recluster =find.all.markers(rgcs_only_recluster, idents.use = c(17,18), test.use = "edgeRQLF", add.to.slot = TRUE, max.cells.per.ident = 500)

# We now have two DE lists (I plan to write a function, where we could take the "Union" of all DE lists)
names(rgcs_only_recluster@de.list)
rgcs_only_recluster_merge = MergeClustersByDendro(rgcs_only_recluster, var.genes = var.genes_rgcs_only)

# Get silhouette coefficients for clustering
sourceCpp("~/Dropbox/CompRepos/sc_R_sparse/Oct2017/Rcpp_matmult.cpp")
sk1 = GetSilhouette(rgcs_only_recluster, reduction.use="pca", dims.use = c(1:40),max.cells.per.ident = 300, return.avg = TRUE, do.plot = TRUE)
sk2 = GetSilhouette(rgcs_only_recluster, reduction.use="tsne", dims.use = c(1:2),max.cells.per.ident = 300, return.avg = TRUE, do.plot = TRUE)
sk3 = GetSilhouette(rgcs_only_recluster, reduction.use="raw", genes.use = var.genes_rgcs_only,max.cells.per.ident = 300, return.avg = TRUE, do.plot = TRUE)






