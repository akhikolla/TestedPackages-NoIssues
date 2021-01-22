#' Root Mean Squared Deviation (RMSD) of Two Protein Conformations
#' @description 
#' RMSD calculation between the atoms of two PDB objects.
#' @usage  
#' pdbrmsd(pdb1, pdb2, start, end, type='all', optimal=FALSE)
#' @param pdb1 PDB object containing reference coordinates of atoms in protein conformation.
#' @param pdb2 PDB object containing coordinates of atoms in protein conformation to compare with pdb1. 
#' @param start The starting residue position for the RMSD calculation. If not supplied, defaults to first residue of chain.
#' @param end The ending residue position for the RMSD calculation. If not supplied, defaults to final residue of chain.
#' @param type Specifies atoms to be included in the calculation.  Can be 'all', 'CA' (CA atoms only), or 'backbone' (CA, N, C, O).
#' @param optimal Apply optimal rotation and superposition? As described in 
#' \url{https://cnx.org/contents/HV-RsdwL@23/Molecular-Distance-Measures}
#
#' @return
#' Returns a list with calculated RMSD value and the optimal rotation matrix.
#'
#' @details
#' Similar to \link[bio3d]{rmsd}, but with implementation in C++.
#'
#' @examples
#' pdbrmsd(nat879, pred879, start=10, end=20, 'all', optimal=TRUE)
#'
#' @export


pdbrmsd <- function(pdb1, pdb2, start=NULL, end=NULL, type='all', optimal=FALSE){
  stopifnot(bio3d::is.pdb(pdb1))
  stopifnot(bio3d::is.pdb(pdb2))
  
  ####Combine into one df ----------------------------------------------------------####
  pdb1new <- pdb1$atom[c("elety", "resid", "resno", "x", "y", "z")]
  pdb2new <- pdb2$atom[c("elety", "resid", "resno", "x", "y", "z")]
  total <- merge(pdb1new,pdb2new,by=c("elety", "resid", "resno")) # , all.x=TRUE)
  total <- as.data.frame.list(total)
  
  #####-------------------Select a range-----------------------------####  
  #Specific Resno/ Resno Range
  if(!is.null(start)){
    if(!is.null(end)){
      df<-subset(total, total$resno >= start & total$resno <= end)
    }
    else{
      df<-subset(total, total$resno >= start)
    }
  }
  else{
    if(!is.null(end)){
      df<-subset(total, total$resno <= end)
    }
    else{
      df<- total
    }
  }
  
  #####----------------------atom type selection and calculation of rmsd----------------------####
  #CA only 
  if('CA' %in% type){
    df<- subset(df, total$elety=="CA")
  }
  
  #CA, C, N, O
  if('backbone' %in% type){
    df<- subset(df, total$elety=="CA" | total$elety=="C" | total$elety=="N" | total$elety=="O")
  }
  ####----------------------------Chose the df (optimal/not optimal)----------------------####
  #If the user types the word optimal, then calculate the optimal matrix, else do nothing
  if(optimal){
    result <- LRMSD(t(as.matrix(df[c("x.x", "y.x", "z.x")])), t(as.matrix(df[c("x.y", "y.y", "z.y")])))
    rmsd <- result$RMSD
    U <- t(result$rotation)
  } 
  else{
    rmsd <- RMSD(t(as.matrix(df[c("x.x", "y.x", "z.x")])), t(as.matrix(df[c("x.y", "y.y", "z.y")])))
    U<-NULL
  } 
  resno_out<-list("RMSD"=rmsd, "optmatrix"= U)
  return(resno_out)
}

