#' Update old R object to accomodate new features
#' @param object scR object
#' @return Returns a scR object compatible with latest changes
#' @examples
#' \dontrun{
#' updated_object = UpdateObject2017(object = old_object)
#' }
#'
ConvertSeuratToscR <- function(object, ident.fxn=getStat1) {
  if (.hasSlot(object, "version")) {
    print(paste0("Input object is of Seurat version ", object@version))
  }
  Count.mat = object@raw.data
  new.object <- scR(count.data = Count.mat, ident.fxn = ident.fxn)
  new.object@count.data = object@raw.data
  new.object@data = object@data
  new.object@data.info = object@meta.data
  colnames(new.object@data.info)[2] = "nTranscripts"
  new.object@data.info$orig = new.object@data.info$orig.ident
  new.object@is.expr = object@is.expr
  new.object@var.genes = object@var.genes
  new.object@cell.names = object@cell.names
  new.object@cluster.tree = object@cluster.tree
  new.object@ident = object@ident 

  
  if (!is.null(object@dr$pca)){
    new.object@pca.x = as.data.frame(object@dr$pca@gene.loadings)
    new.object@pca.rot = as.data.frame(object@dr$pca@cell.embeddings)
  }
  
  if (!is.null(object@dr$ica)){
    new.object@ica.x = as.data.frame(object@dr$ica@gene.loadings)
    new.object@ica.rot = as.data.frame(object@dr$ica@cell.embeddings)
  }
  
  if (!is.null(object@dr$tsne)){
    new.object@tsne.rot = as.data.frame(object@dr$tsne@cell.embeddings)
  }

  
  if (!is.null(object@dr$FItSNE)){
    new.object@tsne.rot = as.data.frame(object@dr$FItSNE@cell.embeddings)
    colnames(new.object@tsne.rot) = c("tsne_1","tsne_2")
  }
  detach("package:Seurat", unload=TRUE)
  return(new.object)
}

FinishObject <- function(object) {

  new.object <-object
  
  # if(length(object@mix.probs) > 0 | length(object@mix.param) > 0 |
  #    length(object@final.prob) > 0 | length(object@insitu.matrix) > 0) {
  #   new.object@spatial <- new(
  #     "spatial.info",
  #     mix.probs = object@mix.probs,
  #     mix.param = object@mix.param,
  #     final.prob = object@final.prob,
  #     insitu.matrix = object@insitu.matrix
  #   )
  # }
  # Conversion from development versions prior to 2.0.0
  if (nrow(object@pca.rot)> 0){
    pca.obj <- new(
      Class = "dim.reduction",
      gene.loadings = as.matrix(object@pca.x),
      cell.embeddings = as.matrix(object@pca.rot),
      key = "PC"
    )
    new.object@dr$pca <- pca.obj
  }
  
  if (nrow(object@ica.rot)> 0){
    ica.obj <- new(
      Class = "dim.reduction",
      gene.loadings = as.matrix(object@ica.x),
      cell.embeddings = as.matrix(object@ica.rot),
      key = "IC"
    )
    new.object@dr$ica <- ica.obj
  }

  if (nrow(object@tsne.rot)> 0){
    tsne.obj <- new(
      Class = "dim.reduction",
      cell.embeddings = as.matrix(object@tsne.rot),
      key = "tsne_"
    )
    new.object@dr$tsne <- tsne.obj
  }
  
  
  return(new.object)
}

# Fills slot in new object with equivalent slot in old object if it still exists
#
# @param slot.name   slot to fill
# @param old.object  object to get slot value from
# @param new.slot    object to set slot value in
#
# @return            returns new object with slot filled
#
FillSlot <- function(slot.name, old.object, new.object){
  new.slot <- tryCatch(
    {
      slot(object = old.object, name = slot.name)
    },
    error = function(err){
      return(NULL)
    }
  )
  if(!is.null(x = new.slot)) {
    slot(new.object, slot.name) <- new.slot
  }
  return(new.object)
}
