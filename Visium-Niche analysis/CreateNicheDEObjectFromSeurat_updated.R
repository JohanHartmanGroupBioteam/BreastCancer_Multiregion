
#' CreateNicheDEObjectFromSeurat
#'
#' This function creates a niche-DE object from a seurat object
#'
#' @param seurat_object A spatial seurat object.Coordinate matrix will be extracted via the
#' seurat function 'GetTissueCoordinates'. The coordiantes will be scaled such that the median nearest neighbor distance is 100.
#' @param assay The assay from which to extract the counts matrix from. The counts matrix
#' will be extracted from the counts slot.
#' @param library_mat Matrix indicating average expression profile for each cell type in the sample
#' @param deconv_mat Deconvolution or cell type assignment matrix of data
#' @param sigma List of kernel bandwidths to use in calculating the effective niche
#' @param Int Boolean of if counts data supplied is integer. Default is true. When performing niche-DE,
#'  Negative binomial regression is performed if True. Linear regression with a gene specific variance is performed if False.
#' @return A niche-DE object
#' @export
CreateNicheDEObjectFromSeurat = function(seurat_object,assay,coordinates, library_mat,deconv_mat,sigma, Int = T){
  print("Creating Niche-DE object")
  #make sure that counts matrix is provided
  if (missing(x = seurat_object)) {
    stop("Must provide seurat object matrix")
  }
  #extract raw counts matrix from seurat object
  sobj_assay = Seurat::GetAssay(seurat_object,assay)
  counts_mat = Matrix::t(sobj_assay@counts)
  #make sure that counts_mat is integers
  if (Int ==T & sum(counts_mat%%1)!=0){
    stop('counts matrix must contain only integers')
  }
  
  #extract coordinate matrix from seurat object
  #coordinate_mat = Seurat::GetTissueCoordinates(seurat_object,image = names(seurat_object@images)[1])
  #slice = names(seurat_object@images)[1]
  coordinates = coordinates
  x = as.numeric(coordinates$x)
  y = as.numeric(coordinates$y)
  coordinate_mat = cbind(x,y)
  rownames(coordinate_mat) = rownames(coordinates)
  #make sure that cell names (rownames) are not null
  if (is.null(x = rownames(x = counts_mat))){
    stop('cell/spot names (rownames) of counts matrix must be non-null')
  }
  if (is.null(x = colnames(x = counts_mat))){
    stop('gene names (colnames) of counts matrix must be non-null')
  }
  #make sure that cell names are unique
  if (anyDuplicated(x = rownames(x = counts_mat))){
    stop('cell/spot names (rownames) of counts matrix must be unique')
  }
  #make sure that gene names are unique
  if (anyDuplicated(x = colnames(x = counts_mat))){
    stop('gene names (colnames) of counts matrix must be unique')
  }
  
  
  #make sure that counts_mat and coordinate_mat have the same cell names in the same order
  if(mean(rownames(counts_mat)==rownames(coordinate_mat))!=1){
    stop('cell/spot names (rownames) of counts matrix and coordinate matrix do not match')
  }
  
  #make sure that counts_mat and deconv_mat have the same cell names in the same order
  if(mean(rownames(counts_mat)==rownames(deconv_mat))!=1){
    stop('cell/spot names (rownames) of counts matrix and deconvolution matrix do not match')
  }
  
  #make sure that deconv_mat and library_mat have the same cell types in the same order
  if(mean(colnames(deconv_mat)==rownames(library_mat))!=1){
    stop('celltypes of deconvolution matrix and reference expression matrix do not match')
  }
  
  
  
  #make sure that sigma is a vector
  if(length(sigma)>0){
    if((is.vector(sigma) && is.atomic(sigma))==F){
      stop('Sigma must be a vector')
    }
    #make sure that sigma is numeric
    if(is.numeric(sigma)==F){
      stop('sigma must be numeric')
    }
  }
  #make sure that sigma has positive length
  if(length(sigma)==0){
    warning('No sigma(kernel bandwidth) values selected. Default values will be used.
            These default values are only appropriate for data that has a similar resolution to 10X VISIUM
            (55 micrometers in diameter).')
  }
  
  
  #calculate number of cells per spot and expected expression
  #get genes that are shared between data and reference expression
  sim_gene = which(colnames(library_mat) %in% colnames(counts_mat))
  #get the gene list
  gene_list = colnames(library_mat)[sim_gene]
  #filter library_mat by removing genes that are not in the counts_mat
  library_mat = library_mat[,sim_gene]
  #get rowsums of the library_matrix (i.e expected library size of a cell type)
  LM = rowSums(library_mat)
  #Get library size of spots
  #get genes that are in the library_mat
  sim_gene = which(colnames(counts_mat)%in% gene_list)
  #filter counts_mat by removing genes that are not in the library_mat
  countsM = as.matrix(counts_mat[,sim_gene])
  colnames(countsM) = colnames(counts_mat)[sim_gene]
  rownames(countsM) = rownames(counts_mat)
  #get library size of each spot
  Lib_spot = rowSums(countsM)
  #make matrix sparse
  countsM = Matrix::Matrix(counts_mat[,sim_gene], sparse=TRUE)
  
  #Get expected library size given pi(deconvolution estimate for each spot)
  EL = deconv_mat%*%as.matrix(LM)
  #get expected number of total cells in a spot
  num_cell = Lib_spot/EL
  #get effective niche
  nst = diag(num_cell[,1])%*%as.matrix(deconv_mat)
  rownames(nst) = rownames(countsM)
  #get expected gene expression given pi
  #EEX = as.matrix(nst)%*%as.matrix(library_mat)
  #rownames(EEX) = rownames(countsM)
  #reorder columns of count data to match that of library reference
  col.order = colnames(library_mat)
  countsM = countsM[,col.order]
  
  #get min spot distance
  D = as.matrix(dist(coordinate_mat),diag = T)
  min_dist = mean(apply(D,2,function(x){sort(x,decreasing = F)[3]}))
  
  #scale coordiante matrix so that min_dist = 100
  scale = 100/min_dist
  coordinate_mat = coordinate_mat*scale
  min_dist = 100
  
  #make sure that counts_mat and library_mat have the same gene names in the same order
  if(mean(colnames(countsM)==colnames(library_mat))!=1){
    stop('gene names (colnames) of counts matrix and library expression matrix do not match')
  }
  
  object = new(Class = 'Niche_DE',counts = countsM,coord = coordinate_mat,
               sigma = sigma,num_cells = nst,ref_expr = library_mat,#null_expected_expression = EEX,
               cell_names = rownames(countsM),cell_types = colnames(deconv_mat),
               gene_names = colnames(countsM),batch_ID = rep(1,nrow(countsM)),
               spot_distance = min_dist,scale = scale,Int = Int)
  #make sure that counts_mat and
  
  A = paste0('Niche-DE object created with ',nrow(object@counts),' observations, ', ncol(object@counts),' genes, ',
             length(unique(object@batch_ID)), ' batch(es), and ', length(object@cell_types), ' cell types.')
  print(A)
  return(object)
}




CreateLibraryMatrixFromSeurat = function(seurat_object,assay){
  #get desired assay from seurat object
  sobj_assay = Seurat::GetAssay(seurat_object,assay)
  #get counts matrix
  data = Matrix::t(GetAssayData(sobj_assay, layer = "counts"))
  #get cell type vector
  cell_type = Seurat::Idents(seurat_object)
  if(mean(rownames(data)== names(cell_type))!=1){
    stop('Data rownames and Cell type matrix names do not match')
  }
  print('Computing average expression profile matrix')
  #get unique cell types
  CT = unique(cell_type)
  n_CT = length(CT)
  L = matrix(NA,n_CT,ncol(data))
  rownames(L) = CT
  colnames(L) = colnames(data)
  #iterate over cell types
  for (j in c(1:n_CT)){
    #get cells that belong to this cell type
    cells = which(cell_type==CT[j])
    #if there are too many cells, downsample
    if(length(cells)>1000){
      print(paste0("Too many cell of type ",CT[j],", downsampling to 1000."))
      cells = sample(cells,1000,replace = F)
    }
    cells = data[cells,]
    L[j,] = apply(cells,2,function(x){mean(x)})
  }
  print('Average expression profile matrix computed.')
  return(L)
}