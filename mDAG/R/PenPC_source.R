connectedComp <-
  function(G,x) {
    # purpose : search for connected components of the vertices set x
    # Input   :
    #     - G : undirected adjacency matrix (edges =TRUE, gaps = FALSE)
    #     - x : set of vertices
    # Output  : the union of connected components of all elements of x (x included)
    if (!identical(dim(G)[1],dim(G)[2])) stop('G must be square matrix')
    if (!all(G==t(G))) stop('G must be symmetric')
    if (sum(diag(G))!=0) stop('diag(G) must be all zero')
    if (!is.vector(x)) stop('x must be a vector')
    p = nrow(G)
    id = seq_len(p)
    
    num_addcomp = 1L
    Ccomp = upx = x # connected components to be updated
    while (num_addcomp != 0) {
      adjmat = as.matrix(G[,upx])
      candx = id[rowSums(adjmat)!=0]
      upx = candx[!candx%in%Ccomp]
      #print(upx)
      num_addcomp = length(upx)
      Ccomp = unique(c(Ccomp,upx))
    }
    return((id %in% Ccomp))
  }
