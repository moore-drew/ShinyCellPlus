#' Generate data files required for shiny app
#'
#' Generate data files required for shiny app. Five files will be generated, 
#' namely (i) the shinycell config \code{prefix_conf.rds}, (ii) the gene 
#' mapping object config \code{prefix_gene.rds}, (iii) the single-cell gene 
#' expression \code{prefix_gexpr.h5}, (iv) the single-cell metadata 
#' \code{prefix_meta.rds} and (v) the defaults for the Shiny app 
#' \code{prefix_def.rds}. A prefix is specified in each file to allow for the 
#' use of multiple single-cell datasets in one Shiny app. Note that both 
#' \code{makeShinyFiles} and \code{makeShinyCodes} functions are ran when 
#' running the wrapper function \code{makeShinyApp}.
#'
#' @param obj input single-cell object for Seurat (v3+) / SingleCellExperiment 
#'   data or input file path for h5ad / loom files
#' @param scConf shinycell config data.table
#' @param gex.assay assay in single-cell data object to use for plotting 
#'   gene expression, which must match one of the following:
#'   \itemize{
#'     \item{Seurat objects}: "RNA" or "integrated" assay, 
#'       default is "RNA"
#'     \item{SCE objects}: "logcounts" or "normcounts" or "counts", 
#'       default is "logcounts"
#'     \item{h5ad files}: "X" or any assay in "layers",
#'       default is "X"
#'     \item{loom files}: "matrix" or any assay in "layers",
#'       default is "matrix"
#'   }
#' @param gex.slot slot in single-cell assay to plot. This is only used 
#'   for Seurat objects (v3+). Default is to use the "data" slot
#' @param gene.mapping specifies whether to convert human / mouse Ensembl gene 
#'   IDs (e.g. ENSG000xxx / ENSMUSG000xxx) into "user-friendly" gene symbols. 
#'   Set this to \code{TRUE} if you are using Ensembl gene IDs. Default is 
#'   \code{FALSE} which is not to perform any conversion. Alternatively, users 
#'   can supply a named vector where \code{names(gene.mapping)} correspond 
#'   to the actual gene identifiers in the gene expression matrix and 
#'   \code{gene.mapping} correspond to new identifiers to map to
#' @param shiny.prefix specify file prefix 
#' @param shiny.dir specify directory to create the shiny app in
#' @param default.gene1 specify primary default gene to show
#' @param default.gene2 specify secondary default gene to show
#' @param default.multigene character vector specifying default genes to 
#'   show in bubbleplot / heatmap
#' @param default.dimred character vector specifying the two default dimension 
#'   reductions. Default is to use UMAP if not TSNE embeddings
#' @param chunk.size number of genes written to h5file at any one time. Lower 
#'   this number to reduce memory consumption. Should not be less than 10
#' @param markers.all boolean flag as to whether to create the 
#'   "Cluster Markers, All" tab and prepare associated Seurat data
#' @param markers.top20 boolean flag as to whether to create the
#'   "Cluster Markers, Top 20" tab and prepare associated Seurat data
#' @param de.genes boolean flag as to whether to create the "Diff. Exp. Genes"
#'   tab and prepre associated Seurat data
#' @param gene.ranks boolean flag as to whether to create the "Gene Signature"
#'   tab and prepare the associated Seurat data
#' @param volc.plot boolean flag as to whether to create the 
#'   "Diff. Gene Exp., Volcano" tab and prepare associated Seurat data
#' @param gene.ont boolean flag as to whether to create the "ToppGene Ontology"
#'   tab and prepate the associated Seurat data
#' @param pval.cutoff upper limit of pval to filter cluster gene expression by (pvals
#'    greater than this are filtered out)
#' @param num.genes max number of most expressed genes to include in ToppGene query
#' @return data files required for shiny app
#'
#' @author John F. Ouyang, Drew Moore
#'
#' @import data.table hdf5r reticulate hdf5r scToppR
#' @importFrom dplyr rename
#'
#' @examples
#' makeShinyFiles(seu, scConf, gex.assay = "RNA", gex.slot = "data",
#'                shiny.prefix = "sc1", shiny.dir = "shinyApp/",
#'                default.gene1 = "GATA3", default.gene2 = "DNMT3L",
#'                default.multigene = c("ANPEP","NANOG","ZIC2","NLGN4X","DNMT3L",
#'                                      "DPPA5","SLC7A2","GATA3","KRT19"),
#'                default.dimred = c("UMAP_1", "UMAP_2"))
#'
#' @export
makeShinyFiles <- function(
  obj, scConf, gex.assay = NA, gex.slot = c("data", "scale.data", "counts"), 
  gene.mapping = FALSE, shiny.prefix = "sc1", shiny.dir = "shinyApp/",
  default.gene1 = NA, default.gene2 = NA, default.multigene = NA, 
  default.dimred = NA, chunk.size = 500, markers.all=FALSE, markers.top20=FALSE, de.genes=FALSE, gene.ranks=FALSE, volc.plot=FALSE, gene.ont=FALSE, pval.cutoff=0.5, num.genes=400){
  ### Preprocessing and checks
  # Generate defaults for gex.assay / gex.slot
  if(class(obj)[1] == "Seurat"){
    # Seurat Object
    if(is.na(gex.assay[1])){gex.assay = "RNA"} # check for existence of this assay, then "Spatial?"  is it possible for a Seurat to have both RNA and Spatial assays in the same object?
    # alternative to these conditionals is to use GetAssayData(obj, assay, slot), 
    # but the slot parameter would still need to be conditionally chosen in some way
    # to either "data" or "counts"
    if(obj@version >= "5.0.0") {
      gex.matdim = dim(LayerData(obj, assay=gex.assay[1], layer=gex.slot[1]))
      gex.rownm = rownames(LayerData(obj, assay=gex.assay[1], layer=gex.slot[1]))
      gex.colnm = colnames(LayerData(obj, assay=gex.assay[1], layer=gex.slot[1]))
    }
    else {
      gex.matdim = dim(slot(obj@assays[[gex.assay[1]]], gex.slot[1]))
      gex.rownm = rownames(slot(obj@assays[[gex.assay[1]]], gex.slot[1]))
      gex.colnm = colnames(slot(obj@assays[[gex.assay[1]]], gex.slot[1]))
    }
    # defGenes = obj@assays[[gex.assay[1]]]@var.features[1:10]
    defGenes = Seurat::VariableFeatures(obj)[1:10]
    if(is.na(defGenes[1])){
      warning(paste0("Variable genes for seurat object not found! Have you ",
                     "ran `FindVariableFeatures` or `SCTransform`?"))
        defGenes = gex.rownm[1:10]
    }
    sc1meta = data.table(sampleID = rownames(obj@meta.data), obj@meta.data)
    
  } else if (class(obj)[1] == "SingleCellExperiment"){
    # SCE Object
    if(is.null(colnames(obj)[1])){
      colnames(obj) = paste0("cell_", seq(ncol(obj)))
    }    # Populate cell IDs if they are not present
    if(is.na(gex.assay[1])){gex.assay = "logcounts"}
    gex.matdim = dim(SummarizedExperiment::assay(obj, gex.assay[1]))
    gex.rownm = rownames(SummarizedExperiment::assay(obj, gex.assay[1]))
    gex.colnm = colnames(SummarizedExperiment::assay(obj, gex.assay[1]))
    defGenes = gex.rownm[1:10]
    sc1meta = SingleCellExperiment::colData(obj)
    sc1meta = data.table(sampleID = rownames(sc1meta), 
                         as.data.frame(sc1meta))
    
  } else if (tolower(tools::file_ext(obj)) == "h5ad"){
    # h5ad file
    if(is.na(gex.assay[1])){gex.assay = "X"}
    # Can just check X since inpH5$layers should have same dimensions
    ad <- import("anndata", convert = FALSE)
    sp <- import('scipy.sparse', convert = FALSE)
    inpH5 = ad$read_h5ad(obj)
    gex.matdim = rev(unlist(py_to_r(inpH5$X$shape)))  
    gex.rownm = py_to_r(inpH5$var_names$values)
    gex.colnm = py_to_r(inpH5$obs_names$values)
    defGenes = gex.rownm[1:10]
    sc1meta = data.table(sampleID = gex.colnm)
    sc1meta = cbind(sc1meta, data.table(py_to_r(inpH5$obs$values)))
    colnames(sc1meta) = c("sampleID", py_to_r(inpH5$obs$columns$values))
    for(i in colnames(sc1meta)[-1]){
      sc1meta[[i]] = unlist(sc1meta[[i]])   # unlist and refactor
      if(as.character(inpH5$obs[i]$dtype) == "category"){
        sc1meta[[i]] = factor(sc1meta[[i]], levels = 
                                py_to_r(inpH5$obs[i]$cat$categories$values))
      }
    } 

  } else if (tolower(tools::file_ext(obj)) == "loom"){
    # loom file
    if(is.na(gex.assay[1])){gex.assay = "matrix"}
    # Can just check matrix since inpLM[["layers"]] should have same dimensions
    inpLM = H5File$new(obj, mode = "r+")
    gex.matdim = rev(inpLM[["matrix"]]$dims)
    gex.rownm = inpLM[["row_attrs"]][["Gene"]]$read()
    for(i in unique(gex.rownm[duplicated(gex.rownm)])){
      gex.rownm[gex.rownm == i] = paste0(i, "-", seq(sum(gex.rownm == i)))
    } # make unique gene names
    gex.colnm = inpLM[["col_attrs"]][["CellID"]]$read()
    defGenes = gex.rownm[1:10]
    cellIdx = which(inpLM[["col_attrs"]]$names == "CellID")
    sc1meta = data.table(sampleID = gex.colnm)
    for(i in inpLM[["col_attrs"]]$names[-cellIdx]){
      tmp = inpLM[["col_attrs"]][[i]]$read()
      if(length(tmp) == nrow(sc1meta)){sc1meta[[i]] = tmp}
    }
     
  } else {
    stop("Only Seurat/SCE objects or h5ad/loom file paths are accepted!")
  }
  
  # Perform gene.mapping if specified (also map defGenes)
  if(gene.mapping[1] == TRUE){
    if(sum(grepl("^ENSG000", gex.rownm)) >= sum(grepl("^ENMUSG000", gex.rownm))){
      tmp1 = fread(system.file("extdata", "geneMapHS.txt.gz", 
                               package = "ShinyCellPlus"))
    } else {
      tmp1 = fread(system.file("extdata", "geneMapMM.txt.gz", 
                               package = "ShinyCellPlus"))
    }
    gene.mapping = tmp1$geneName
    names(gene.mapping) = tmp1$geneID
  }
  # Check if gene.mapping is partial or not
  if(gene.mapping[1] == FALSE){
    gene.mapping = gex.rownm      
    names(gene.mapping) = gex.rownm    # Basically no mapping
  } else {
    if(!all(gex.rownm %in% names(gene.mapping))){
      # warning("Mapping for some gene identifiers are not provided!")
      tmp1 = gex.rownm[gex.rownm %in% names(gene.mapping)]
      tmp1 = gene.mapping[tmp1]
      tmp2 = gex.rownm[!gex.rownm %in% names(gene.mapping)]
      names(tmp2) = tmp2
      gene.mapping = c(tmp1, tmp2)
    } 
    gene.mapping = gene.mapping[gex.rownm]
  }
  defGenes = gene.mapping[defGenes]
  
  # Check default.gene1 / default.gene2 / default.multigene
  default.gene1 = default.gene1[1]
  default.gene2 = default.gene2[1]
  if(is.na(default.gene1)){default.gene1 = defGenes[1]}
  if(is.na(default.gene2)){default.gene2 = defGenes[2]}
  if(is.na(default.multigene[1])){default.multigene = defGenes}
  if(default.gene1 %in% gene.mapping){
    default.gene1 = default.gene1
  } else if(default.gene1 %in% names(gene.mapping)){
    default.gene1 = gene.mapping[default.gene1]
  } else {
    warning("default.gene1 doesn't exist in gene expression, using defaults...")
    default.gene1 = defGenes[1]
  }
  if(default.gene2 %in% gene.mapping){
    default.gene2 = default.gene2
  } else if(default.gene2 %in% names(gene.mapping)){
    default.gene2 = gene.mapping[default.gene2]
  } else {
    warning("default.gene2 doesn't exist in gene expression, using defaults...")
    default.gene2 = defGenes[2]
  }
  if(all(default.multigene %in% gene.mapping)){
    default.multigene = default.multigene
  } else if(all(default.multigene %in% names(gene.mapping))){
    default.multigene = gene.mapping[default.multigene]
  } else {
    warning(paste0("default.multigene doesn't exist in gene expression, ", 
                   "using defaults..."))
    default.multigene = defGenes
  }
  

  ### Actual object generation
  # Make XXXmeta.rds and XXXconf.rds (updated with dimred info)
  sc1conf = scConf
  sc1conf$dimred = FALSE
  sc1meta = sc1meta[, c("sampleID", as.character(sc1conf$ID)), with = FALSE]
  # Factor metadata again
  for(i in as.character(sc1conf[!is.na(fID)]$ID)){
    sc1meta[[i]] = factor(sc1meta[[i]],
                          levels = strsplit(sc1conf[ID == i]$fID, "\\|")[[1]])
    levels(sc1meta[[i]]) = strsplit(sc1conf[ID == i]$fUI, "\\|")[[1]]
    sc1conf[ID == i]$fID = sc1conf[ID == i]$fUI
  }
  # Extract dimred and append to both XXXmeta.rds and XXXconf.rds...
  if(class(obj)[1] == "Seurat"){
    # Seurat Object
    for(iDR in names(obj@reductions)){
      drMat = obj@reductions[[iDR]]@cell.embeddings
      if(ncol(drMat) > 5){drMat = drMat[, 1:5]}  # Take first 5 components only
      drMat = drMat[sc1meta$sampleID, ]          # Ensure ordering
      drMat = as.data.table(drMat)
      sc1meta = cbind(sc1meta, drMat)            
      
      # Update sc1conf accordingly
      tmp = data.table(ID = colnames(drMat), UI = colnames(drMat),
                       fID = NA, fUI = NA, fCL = NA, fRow = NA, 
                       default = 0, grp = FALSE, split = FALSE, dimred = TRUE)
      tmp$UI = gsub("_", "", tmp$UI)
      sc1conf = rbindlist(list(sc1conf, tmp))
    }

  } else if (class(obj)[1] == "SingleCellExperiment"){
    # SCE Object
    for(iDR in SingleCellExperiment::reducedDimNames(obj)){
      drMat = SingleCellExperiment::reducedDim(obj, iDR)
      if(ncol(drMat) > 5){drMat = drMat[, 1:5]}  # Take first 5 components only
      if(is.null(colnames(drMat))){
        colnames(drMat) = paste0(iDR, seq(ncol(drMat)))
      }
      drMat = drMat[sc1meta$sampleID, ]          # Ensure ordering
      drMat = as.data.table(drMat)
      sc1meta = cbind(sc1meta, drMat)            
      
      # Update sc1conf accordingly
      tmp = data.table(ID = colnames(drMat), UI = colnames(drMat),
                       fID = NA, fUI = NA, fCL = NA, fRow = NA, 
                       default = 0, grp = FALSE, dimred = TRUE, split = FALSE)
      tmp$UI = gsub("_", "", tmp$UI)
      sc1conf = rbindlist(list(sc1conf, tmp))
    }
    
  } else if (tolower(tools::file_ext(obj)) == "h5ad"){
    # h5ad file
    for(iDR in py_to_r(inpH5$obsm_keys())){
      drMat = py_to_r(inpH5$obsm[iDR])
      tmpName = gsub("pca", "pc", gsub("X_", "", iDR))
      tmpName = paste0(tmpName, "_", 1:ncol(drMat))
      colnames(drMat) = tmpName
      if(ncol(drMat) > 5){drMat = drMat[, 1:5]}  # Take first 5 components only
      drMat = as.data.table(drMat)
      sc1meta = cbind(sc1meta, drMat)
      
      # Update sc1conf accordingly
      tmp = data.table(ID = colnames(drMat), UI = colnames(drMat),
                       fID = NA, fUI = NA, fCL = NA, fRow = NA, 
                       default = 0, grp = FALSE, dimred = TRUE, split = FALSE)
      tmp$UI = gsub("_", "", tmp$UI)
      sc1conf = rbindlist(list(sc1conf, tmp))
    }
    
  } else if (tolower(tools::file_ext(obj)) == "loom"){
    # loom file
    nDR = inpLM[["col_attrs"]]$names[
      grep("pca|tsne|umap", inpLM[["col_attrs"]]$names, ignore.case = TRUE)]
    for(iDR in nDR){
      drMat = t(inpLM[["col_attrs"]][[iDR]]$read())
      colnames(drMat) = paste0(iDR, "_", 1:ncol(drMat))
      if(ncol(drMat) > 5){drMat = drMat[, 1:5]}  # Take first 5 components only
      drMat = as.data.table(drMat)
      sc1meta = cbind(sc1meta, drMat)
      
      # Update sc1conf accordingly
      tmp = data.table(ID = colnames(drMat), UI = colnames(drMat),
                       fID = NA, fUI = NA, fCL = NA, fRow = NA, 
                       default = 0, grp = FALSE, dimred = TRUE, split = FALSE)
      tmp$UI = gsub("_", "", tmp$UI)
      sc1conf = rbindlist(list(sc1conf, tmp))
    }

  }
  sc1conf$ID = as.character(sc1conf$ID)     # Remove levels
  
  # Make XXXgexpr.h5
  if(!dir.exists(shiny.dir)){dir.create(shiny.dir)}
  filename = paste0(shiny.dir, "/", shiny.prefix, "gexpr.h5")
  sc1gexpr <- H5File$new(filename, mode = "w")
  sc1gexpr.grp <- sc1gexpr$create_group("grp")
  sc1gexpr.grp.data <- sc1gexpr.grp$create_dataset(
    "data",  dtype = h5types$H5T_NATIVE_FLOAT,
    space = H5S$new("simple", dims = gex.matdim, maxdims = gex.matdim),
    chunk_dims = c(1,gex.matdim[2]))
  chk = chunk.size
  while(chk > (gex.matdim[1]-8)){
    chk = floor(chk / 2)     # Account for cases where nGene < chunk.size
  } 
  if(class(obj)[1] == "Seurat"){
    # Seurat Object
    if(obj@version >= "5.0.0") {
      for(i in 1:floor((gex.matdim[1]-8)/chk)){
        sc1gexpr.grp.data[((i-1)*chk+1):(i*chk), ] <- as.matrix(
          LayerData(obj, assay=gex.assay[1], layer=gex.slot[3])[((i-1)*chk+1):(i*chk),])
      }
      sc1gexpr.grp.data[(i*chk+1):gex.matdim[1], ] <- as.matrix(
        LayerData(obj, assay=gex.assay[1], layer=gex.slot[3])[(i*chk+1):gex.matdim[1],])
    }
    else {
      for(i in 1:floor((gex.matdim[1]-8)/chk)){
        sc1gexpr.grp.data[((i-1)*chk+1):(i*chk), ] <- as.matrix(
          slot(obj@assays[[gex.assay[1]]], gex.slot[1])[((i-1)*chk+1):(i*chk),])
      }
      sc1gexpr.grp.data[(i*chk+1):gex.matdim[1], ] <- as.matrix(
        slot(obj@assays[[gex.assay[1]]], gex.slot[1])[(i*chk+1):gex.matdim[1],])
    }
  } else if (class(obj)[1] == "SingleCellExperiment"){
    # SCE Object
    for(i in 1:floor((gex.matdim[1]-8)/chk)){
      sc1gexpr.grp.data[((i-1)*chk+1):(i*chk), ] <- as.matrix(
        SummarizedExperiment::assay(obj, gex.assay[1])[((i-1)*chk+1):(i*chk),])
    }
    sc1gexpr.grp.data[(i*chk+1):gex.matdim[1], ] <- as.matrix(
      SummarizedExperiment::assay(obj, gex.assay[1])[(i*chk+1):gex.matdim[1],])
    
  } else if (tolower(tools::file_ext(obj)) == "h5ad"){
    # h5ad file
    if(gex.assay == "X"){
      scGEX = Matrix::t(py_to_r(sp$csc_matrix(inpH5$X)))
    } else {
      scGEX = Matrix::t(py_to_r(sp$csc_matrix(inpH5$layers[[gex.assay]])))
    }
    for(i in 1:floor((gex.matdim[1]-8)/chk)){
      sc1gexpr.grp.data[((i-1)*chk+1):(i*chk), ] <- as.matrix(
        scGEX[((i-1)*chk+1):(i*chk),])
    }
    sc1gexpr.grp.data[(i*chk+1):gex.matdim[1], ] <- as.matrix(
      scGEX[(i*chk+1):gex.matdim[1],])

  } else if (tolower(tools::file_ext(obj)) == "loom"){
    # loom file
    if(gex.assay == "matrix"){
      for(i in 1:floor((gex.matdim[1]-8)/chk)){
        sc1gexpr.grp.data[((i-1)*chk+1):(i*chk), ] <- t(
          inpLM[["matrix"]][, ((i-1)*chk+1):(i*chk)])
      }
      sc1gexpr.grp.data[(i*chk+1):gex.matdim[1], ] <- t(
        inpLM[["matrix"]][, (i*chk+1):gex.matdim[1]])
    } else {
      for(i in 1:floor((gex.matdim[1]-8)/chk)){
        sc1gexpr.grp.data[((i-1)*chk+1):(i*chk), ] <- t(
          inpLM[["layers"]][[gex.assay]][, ((i-1)*chk+1):(i*chk)])
      }
      sc1gexpr.grp.data[(i*chk+1):gex.matdim[1], ] <- t(
        inpLM[["layers"]][[gex.assay]][, (i*chk+1):gex.matdim[1]])
    }
  }

  # sc1gexpr.grp.data[, ] <- as.matrix(gex.matrix[,])
  sc1gexpr$close_all()
  if(!isTRUE(all.equal(sc1meta$sampleID, gex.colnm))){
    sc1meta$sampleID = factor(sc1meta$sampleID, levels = gex.colnm)
    sc1meta = sc1meta[order(sampleID)]
    sc1meta$sampleID = as.character(sc1meta$sampleID)
  }
  
  # Make XXXgenes.rds
  sc1gene = seq(gex.matdim[1])
  names(gene.mapping) = NULL
  names(sc1gene) = gene.mapping
  sc1gene = sc1gene[order(names(sc1gene))]
  sc1gene = sc1gene[order(nchar(names(sc1gene)))]
  
  # Make XXXdef.rds (list of defaults)
  if(all(default.dimred %in% sc1conf[dimred == TRUE]$ID)){
    default.dimred[1] = sc1conf[ID == default.dimred[1]]$UI
    default.dimred[2] = sc1conf[ID == default.dimred[2]]$UI
  } else if(all(default.dimred %in% sc1conf[dimred == TRUE]$UI)) {
    default.dimred = default.dimred    # Nothing happens
  } else {
    warn = TRUE
    if(is.na(default.dimred[1])){
      default.dimred = "umap"
      warn = FALSE
    }
    # Try to guess... and give a warning
    guess = gsub("[0-9]", "", default.dimred[1])
    if(length(grep(guess, sc1conf[dimred == TRUE]$UI, ignore.case=TRUE)) >= 2){
      default.dimred = sc1conf[dimred == TRUE]$UI[
        grep(guess, sc1conf[dimred == TRUE]$UI, ignore.case = TRUE)[1:2]]
    } else {
      nDR = length(sc1conf[dimred == TRUE]$UI)
      default.dimred = sc1conf[dimred == TRUE]$UI[(nDR-1):nDR]
    }
    if(warn){
      warning(paste0("default.dimred not found, switching to ", 
                     default.dimred[1], " and ", default.dimred[1]))
    } # Warn if user-supplied default.dimred is not found
  }
  # Note that we stored the display name here
  sc1def = list()
  sc1def$meta1 = sc1conf[default == 1]$UI   # Use display name
  sc1def$meta2 = sc1conf[default == 2]$UI   # Use display name 
  sc1def$gene1 = default.gene1              # Actual == Display name
  sc1def$gene2 = default.gene2              # Actual == Display name
  sc1def$genes = default.multigene          # Actual == Display name
  sc1def$dimred = default.dimred            # Use display name 
  tmp = nrow(sc1conf[default != 0 & grp == TRUE])
  if(tmp == 2){
    sc1def$grp1 = sc1def$meta1
    sc1def$grp2 = sc1def$meta2
  } else if(tmp == 1){
    sc1def$grp1 = sc1conf[default != 0 & grp == TRUE]$UI
    if(nrow(sc1conf[default == 0 & grp == TRUE]) == 0){
      sc1def$grp2 = sc1def$grp1
    } else {
      sc1def$grp2 = sc1conf[default == 0 & grp == TRUE]$UI[1]
    }
  } else {
    sc1def$grp1 = sc1conf[default == 0 & grp == TRUE]$UI[1]
    if(nrow(sc1conf[default == 0 & grp == TRUE]) < 2){
      sc1def$grp2 = sc1def$grp1
    } else {
      sc1def$grp2 = sc1conf[default == 0 & grp == TRUE]$UI[2]
    }
  }
  sc1conf = sc1conf[, -c("fUI", "default"), with = FALSE]
  
  

  ### Saving objects
  #saveRDS(sc1conf, file = paste0(shiny.dir, "/", shiny.prefix, "conf.rds"))
  saveRDS(sc1meta, file = paste0(shiny.dir, "/", shiny.prefix, "meta.rds"))
  saveRDS(sc1gene, file = paste0(shiny.dir, "/", shiny.prefix, "gene.rds"))
  saveRDS(sc1def,  file = paste0(shiny.dir, "/", shiny.prefix, "def.rds"))

  sc1conf$extra_tabs <- FALSE
  sc1conf$DEs <- FALSE
  sc1conf$gene_ranks <- FALSE
  sc1conf$markers <- FALSE

  
  if(class(obj)[1]=='Seurat') {
    if(markers.all==TRUE) { 
      #sc1conf$extra_tabs[1] = TRUE
      if(!is.null(obj@misc$markers$presto$all)) {
        cat("creating .rds for all presto markers...\n\n")
        m_all <- obj@misc$markers$presto$all
        names(m_all)[names(m_all) == 'group'] <- 'cluster'
        #sc1conf$extra_tabs[1] = TRUE
        saveRDS(m_all, file=paste0(shiny.dir, "/", shiny.prefix, "m_all.rds"))
        sc1conf$markers[1] <- TRUE
        sc1conf$extra_tabs[1] <- TRUE
      }
      else if(!is.null(obj@misc$markers$seurat$all)) {
        cat("creating .rds for Seurat::FindAllMarkers() results...\n\n")
        m_all <- obj@misc$markers$seurat$all
        saveRDS(m_all, file=paste0(shiny.dir, "/", shiny.prefix, "m_all.rds"))
        sc1conf$markers[1] <- TRUE
        sc1conf$extra_tabs[1] <- TRUE
      }
      else {
        sc1conf$markers[1] <- "ERROR"
        cat("Warning: 'markers.all' data not found in Seurat (structure expected: seurat@misc$markers$presto$all);\ncorresponding Shiny tab ('Cluster Markers, All') will be created but with an error message instead of what is expected...\n\n")
        #sc1conf$extra_tabs[1] = FALSE
      }
    }

    if(markers.top20==TRUE) {
      #sc1conf$extra_tabs[2] = TRUE
      if(!is.null(obj@misc$markers$presto$top_20)) {
        cat("creating .rds for top 20 presto markers...\n\n")
        m_t20 <- obj@misc$markers$presto$top_20
        names(m_t20)[names(m_t20) == 'group'] <- 'cluster'

        if("10" %in% colnames(m_t20)) {
          rank_i <- which(colnames(m_t20) == "rank")
          sorted_colnames <- c("rank", sort(as.numeric(colnames(m_t20[,-c(rank_i)]))))
        } else {
          sorted_colnames <- colnames(m_t20)
        }
        m_t20 <- m_t20[,sorted_colnames]

        #sc1conf$extra_tabs[2] = TRUE
        saveRDS(m_t20, file=paste0(shiny.dir, "/", shiny.prefix, "m_t20.rds"))
        sc1conf$markers[2] <- TRUE
        sc1conf$extra_tabs[2] = TRUE
      }
      else if(!is.null(obj@misc$markers$seurat$top_20)) {
        cat("creating .rds for Seurat::FindAllMarkers() top 20 results...\n\n")
        m_t20 <- obj@misc$markers$seurat$top_20
        saveRDS(m_t20, file=paste0(shiny.dir, "/", shiny.prefix, "m_t20.rds"))
        sc1conf$markers[2] <- TRUE
        sc1conf$extra_tabs[2] = TRUE
      }
      else {
        sc1conf$markers[2] <- "ERROR"
        cat("Warning: 'markers.top20' data not found in Seurat (structure expected: seurat@misc$markers$presto$top_20);\ncorresponding Shiny tab ('Cluster Markers, Top 20') will be created but with an error message instead of what is expected...\n\n")
        #sc1conf$extra_tabs[2] = FALSE
      }
    }

    if(de.genes==TRUE) { 
      #sc1conf$extra_tabs[3] = TRUE
      if(!is.null(obj@misc$DE_genes$libra$overall)) {
        cat("creating .rds for differentially expressed genes...\n\n")
        de_genes <- obj@misc$DE_genes$libra$overall
        #sc1conf$extra_tabs[3] = TRUE
        de_genes$de_family <- NULL
        de_genes$de_method <- NULL
        de_genes$de_type <- NULL
        de_genes$de_name <- as.factor(de_genes$de_name)
        sc1conf$DEs[1] <- paste0(levels(de_genes$de_name), collapse="|")
        saveRDS(de_genes, file=paste0(shiny.dir, "/", shiny.prefix, "de_genes.rds"))
        sc1conf$extra_tabs[3] <- TRUE
      }
      else if(!is.null(obj@misc$DE_genes$seurat$overall)) {
        cat("creating .rds for differentially expressed genes...\n\n")
        de_genes <- obj@misc$DE_genes$seurat$overall
        #sc1conf$extra_tabs[3] = TRUE
        de_genes$de_family <- NULL
        de_genes$de_method <- NULL
        de_genes$de_type <- NULL
        de_genes$de_name <- as.factor(de_genes$de_name)
        sc1conf$DEs[1] <- paste0(levels(de_genes$de_name), collapse="|")
        saveRDS(de_genes, file=paste0(shiny.dir, "/", shiny.prefix, "de_genes.rds"))
        sc1conf$extra_tabs[3] <- TRUE
      }
      else {
        sc1conf$DEs[1] <- "ERROR" # |
        cat("Warning: 'de.genes' data not found in Seurat (structure expected: seurat@misc$DE_genes$libra$overall);\ncorresponding Shiny tab ('Diff. Exp. Genes') will be created but with an error message instead of what is expected...\n\n")
        #sc1conf$extra_tabs[3] = FALSE
      }
    }

    if(gene.ranks==TRUE) {
      #sc1conf$extra_tabs[4] = TRUE
      if(!is.null(obj@misc$gene_ranks$aucell$all)) {
        cat("creating .rds for gene ranks (may take a while)...\n\n")
        gene_ranks <- obj@misc$gene_ranks$aucell$all
        #sc1conf$extra_tabs[4] = TRUE
        saveRDS(gene_ranks, file=paste0(shiny.dir, "/", shiny.prefix, "gene_ranks.rds"))
        sc1conf$gene_ranks[1] <- TRUE
        sc1conf$extra_tabs[4] <- TRUE
      }
      else {
        sc1conf$gene_ranks[1] <- "ERROR"
        cat("Warning: 'gene.ranks' data not found in Seurat (structure expected: seurat@misc$gene_ranks$aucell$all);\ncorresponding Shiny tab ('Gene Signature') will be created but with an error message instead of what is expected...\n\n")
        #sc1conf$extra_tabs[4] = FALSE
      }
    }

    if(volc.plot==TRUE) {
      #sc1conf$extra_tabs[5] = TRUE
      #TODO: tighten these two conditionals up, there's only one line of difference
      if(!is.null(obj@misc$DE_genes$libra$overall)) {
        cat("creating .rds of differentially expressed genes formatted for ggvolc...\n\n")
        # check to see if has correct column names?
        # or maybe wrap this in a try-catch for rename
        # errors?
        de_genes <- obj@misc$DE_genes$libra$overall
        de_genes <- dplyr::rename(de_genes, genes = gene) 
        de_genes <- dplyr::rename(de_genes, log2FoldChange = avg_logFC) 
        de_genes <- dplyr::rename(de_genes, pvalue = p_val_adj)
        de_genes$de_name <- as.factor(de_genes$de_name)
        sc1conf$DEs[2] <- paste0(levels(de_genes$de_name), collapse="|")
        columns <- c("genes", "log2FoldChange", "p_val", "pvalue", "de_family", "de_method", "de_type", "de_name")
        uniq_cols <- names(de_genes)[ !(names(de_genes) %in% columns) ]
        sc1conf$DEs[3] <- paste0(uniq_cols, collapse="|")
        saveRDS(de_genes, file=paste0(shiny.dir, "/", shiny.prefix, "de_genes_ggvolc.rds"))
        sc1conf$extra_tabs[5] <- TRUE
      }
      else if(!is.null(obj@misc$DE_genes$seurat$overall)) {
        cat("creating .rds of differentially expressed genes formatted for ggvolc...\n\n")
        de_genes <- obj@misc$DE_genes$seurat$overall
        de_genes <- dplyr::rename(de_genes, genes = gene)
        de_genes <- dplyr::rename(de_genes, log2FoldChange = avg_log2FC)
        de_genes <- dplyr::rename(de_genes, pvalue = p_val_adj)
        de_genes$de_name <- as.factor(de_genes$de_name)
        sc1conf$DEs[2] <- paste0(levels(de_genes$de_name), collapse="|")
        columns <- c("genes", "log2FoldChange", "p_val", "pvalue", "de_family", "de_method", "de_type", "de_name")
        uniq_cols <- names(de_genes)[ !(names(de_genes) %in% columns) ]
        sc1conf$DEs[3] <- paste0(uniq_cols, collapse="|")
        saveRDS(de_genes, file=paste0(shiny.dir, "/", shiny.prefix, "de_genes_ggvolc.rds"))
        sc1conf$extra_tabs[5] <- TRUE
      }
      else { # if want to make plot but data is not in seurat
        # if(file.exists(paste0(shiny.dir, "/", shiny.prefix, "de_genes_ggvolc.rds"))) { # but a correctly named file was previously created
        #   # ask them whether to use the file, delete the file (do not auto delete for safety / liability reasons), or abort the session
        #   cat("\nData for 'volc.plot' currently not found in Seurat object (structure expected: seurat@misc$DE_genes$libra$overall)...\nThere is a previously created data file (\"de_genes_ggvolc.rds\"), would you like to use this file?\n(WARNING: Make sure this file belongs to your Seurat dataset or it will lead to incorrect data presentation!) \n\t1. use this file (make sure data corresponds to your Seurat object's state!)\n\t2. delete this file and continue (an error will be shown on the tab for missing file)\n\t3. do not use this tab (file will still exist)\n\t4. to abort")
        #   while(TRUE) {
        #     resp <- readline("\t>")
        #     if(resp=="1") {
        #       cat("\n")
        #       de_genes <- readRDS(paste0(shiny.dir, "/", shiny.prefix, "de_genes_ggvolc.rds"))
        #       sc1conf$DEs[2] <- paste0(levels(de_genes$de_name), collapse="|")
        #       columns <- c("genes", "log2FoldChange", "p_val", "pvalue", "de_family", "de_method", "de_type", "de_name")
        #       uniq_cols <- names(de_genes)[ !(names(de_genes) %in% columns) ]
        #       sc1conf$DEs[3] <- paste0(uniq_cols, collapse="|")
        #       break
        #     }
        #     else if(resp=="2") {
        #       cat("\n")
        #       file.remove(paste0(shiny.dir, "/", shiny.prefix, "de_genes_ggvolc.rds"))
        #       # needed if you are already deleting the file?
        #       # should only be needed if data isn't in seurat,
        #       # but the associated file exists, which this 
        #       # section should be taking care of...
        #       # could use this to catch some sort of unknown
        #       # state that shouldn't be reached?
        #       sc1conf$DEs[2] <- "ERROR" # |
        #       sc1conf$DEs[3] <- "ERROR" # |
        #       break
        #     }
        #     else if(resp=="3") {
        #       # problem with this is that the volc.plot flag
        #       # stays TRUE and gets carried over to the 
        #       # makeShinyCodes.R script and still makes the
        #       # tab.  hotfixed with the double sc1conf$DE
        #       # ERRORS, but that's obviously not what is
        #       # expected by user... need to find a way to
        #       # transfer this info to writer.R without
        #       # changing function signatures for ease?
        #       # could make it so that the page flags
        #       # get condensed into a list returned by this
        #       # script and used in the makeShinyCodes.R
        #       # script?  
        #       cat("\n")
        #       sc1conf$extra_tabs[5] <- FALSE
        #       sc1conf$DEs[2] <- "ERROR" # |
        #       sc1conf$DEs[3] <- "ERROR" # |
        #       break
        #     }
        #     else if(resp=="4") {
        #       abort()
        #     }
        #     else {
        #       cat("\nIncorrect response, please choose from the following:\n\t1. use this file (make sure data corresponds to your Seurat object's state!)\n\t2. delete this file and continue (an error will be shown on the tab for missing file)\n\t3. do not use this tab (file will still exist)\n\t4. abort")
        #     }
        #   }
        # }
        # else {
          # instead of all of that above, what about just setting extra_tabs[5]
          # to FALSE, that way it doesn't load a file, sc1conf$DEs[2] and [3] 
          # to "ERROR", and then running that to show the error?
          # how would you get this to work with the spreadsheet tabs?

          #sc1conf$extra_tabs[5] = FALSE
          sc1conf$DEs[2] <- "ERROR" # |
          sc1conf$DEs[3] <- "ERROR" # |

          # "or pre-existing file (name expected: de_genes_ggvolc.rds)"
          cat("Warning: 'volc.plot' data not found in Seurat (structure expected: seurat@misc$DE_genes$libra$overall);\ncorresponding Shiny tab ('Diff. Gene Exp., Volcano') will be created but with an error message instead of what is expected...\n\n")
        #}
        #sc1conf$extra_tabs[5] = FALSE
      }
    }

    if(gene.ont==TRUE) {
      #sc1conf$extra_tabs[6] = TRUE
      if(!is.null(obj@misc$DE_genes$libra$overall)) {
        cat("creating .rds for all ToppGene ontology...\n\n")
        de_genes <- obj@misc$DE_genes$libra$overall
        de_names <- levels(as.factor(de_genes$de_name))
        gene_ont <- data.frame()
        for(name in de_names) {
          print(name)
          de_genes_subset <- filter(de_genes, de_name == name)
          de_genes_subset$cell_type <- as.character(de_genes_subset$cell_type)
          gene_ont_subset <- toppFun(de_genes_subset, cluster_col="cell_type", p_val_col="p_val_adj", num_genes=num.genes, pval_cutoff=pval.cutoff, min_genes=1)
          gene_ont_subset$de_name <- name
          gene_ont <- rbind(gene_ont, gene_ont_subset)
        }
        sc1conf$DEs[4] <- paste0(de_names, collapse="|")
        saveRDS(gene_ont, file=paste0(shiny.dir, "/", shiny.prefix, "gene_ont.rds"))
        sc1conf$extra_tabs[6] <- TRUE
      }
      else {
        sc1conf$DEs[4] <- "ERROR"
        cat("Warning: 'de.genes' data not found in Seurat (structure expected: seurat@misc$DE_genes$libra$overall);\ncorresponding Shiny tab ('ToppGene Ontology') will be created but with an error message instead of what is expected...\n\n")
      }
    }
  }

  saveRDS(sc1conf, file = paste0(shiny.dir, "/", shiny.prefix, "conf.rds"))
  return(sc1conf)
}


