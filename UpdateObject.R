#' Update old R object to accomodate new features
#' @param object scR object
#' @return Returns a scR object compatible with latest changes
#' @examples
#' \dontrun{
#' updated_object = UpdateObject2017(object = old_object)
#' }
#'
UpdateObjectOct17 <- function(object) {
  if (.hasSlot(object, "version")) {
    if (object@version == "Oct2017"){
    cat("Object representation is consistent with the most current version.\n")
    return(object)
    }
  }
  Count.mat = object@count.data[rownames(object@data), colnames(object@data)]
  new.object <- scR(count.data = Count.mat, ident.fxn = object@ident.fxn)
  new.object <- setup(new.object, min.cells = 0, min.genes = 0, verbose = FALSE)
  
  # Fill in new slots
  new.slots <- slotNames(new.object)
  for(s in new.slots){
    new.object <- FillSlot(
      slot.name = s,
      old.object = object,
      new.object = new.object
    )
  }
  
  print(1)
  new.object@ident = new.object@data.info$orig
  names(new.object@ident) = rownames(new.object@data.info)
  
  
  # Copy over old slots if they have info stored
  if(length(object@kmeans.obj) > 0){
    new.object@kmeans$gene.kmeans.obj <- object@kmeans.obj
  }
  if(length(object@kmeans.col) >0 ){
    new.object@kmeans$cell.kmeans.obj <- object@kmeans.col
  }

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
  if ((.hasSlot(object, "dr"))) {
    for (i in 1:length(object@dr)) {
      new.object@dr[[i]]@cell.embeddings <- object@dr[[i]]@rotation
      new.object@dr[[i]]@gene.loadings <- object@dr[[i]]@x
      new.object@dr[[i]]@gene.loadings.full <- object@dr[[i]]@x.full
      new.object@dr[[i]]@sdev <- object@dr[[i]]@sdev
      new.object@dr[[i]]@key <- object@dr[[i]]@key
      new.object@dr[[i]]@misc <- object@dr[[i]]@misc
    }
  }
  # Conversion from release versions prior to 2.0.0
  # Slots to replace: pca.x, pca.rot, pca.x.full, tsne.rot, ica.rot, ica.x,
  #                   tsne.rot
  else{
    pca.sdev <- object@pca.obj[[1]]$sdev
    if (is.null(x = pca.sdev)) {
      pca.sdev <- object@pca.obj[[1]]$d
    }
    pca.obj <- new(
      Class = "dim.reduction",
      gene.loadings = as.matrix(object@pca.x),
      gene.loadings.full = as.matrix(object@pca.x.full),
      cell.embeddings = as.matrix(object@pca.rot),
      sdev = pca.sdev,
      key = "PC"
    )
    new.object@dr$pca <- pca.obj
    ica.obj <- new(
      Class = "dim.reduction",
      gene.loadings = as.matrix(object@ica.x),
      cell.embeddings = as.matrix(object@ica.rot),
      key = "IC"
    )
    new.object@dr$ica <- ica.obj
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
