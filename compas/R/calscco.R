#' Calculate Cartesian coordinates of side chains
#' @description 
#' Rotates the free side chain dihedral angles of an amino acid to the specified values. Calculates the updated Cartesian coordinates of all the atoms of that amino acid side chain.
#
#' @param pdb A PDB object 
#' @param resno The residue number of the amino acid side chain to rotate
#' @param chi A vector of dihedral angles (in degrees), with length matching the number of free side chain dihedral angles for that amino acid type.  See \link{atomdeps} for definitions.
#
#' @return
#' Returns a PDB object with updated coordinates of side chain atoms in "resno".
#'
#' @details
#' Calls \link{calCo} successively for each atom in the amino acid side chain, using the bond parameters defined in \link{atomdeps}.
#'
#' @examples
#' ## Position 10 of nat879 is ASP with 2 side chain dihedrals chi1 and chi2
#' nat879$atom[nat879$atom$resno==10,]
#' pdbn <- calscco(nat879,10,c(60.0,-80.0))
#' pdbn$atom[pdbn$atom$resno==10,]
#'
#' @export

calscco <-function(pdb, resno, chi){
  
  atyperows<-which(pdb$atom$resno == resno, arr.ind = TRUE)
  if (length(atyperows) == 0)
    stop(paste("Resno", resno, "does not exist in pdb"))
  
  rtype <- pdb$atom$resid[atyperows[1]]

  numchi <- sum(is.na(compas::atomdeps[[rtype]]$tangle))
  if (length(chi) != numchi)
    stop(paste("Expecting a vector of length", numchi, "chi angles for", rtype))

  
  #####go to atomdep and find which atoms are needed (names)
  xnew<-compas::atomdeps[[rtype]]$names
  
  chicount <- 0
  ####Run calco on each of the atoms----
  for(j in 1:length(xnew)){
    
    
    elements<-compas::atomdeps[[rtype]]$matx[j,]
    
    
    xpdbca1<-as.double(pdb$atom[pdb$atom$resno==resno & pdb$atom$elety==elements[1], c('x', 'y', 'z')])
    xpdbca2<-as.double(pdb$atom[pdb$atom$resno==resno & pdb$atom$elety==elements[2], c('x', 'y', 'z')])
    xpdbca3<-as.double(pdb$atom[pdb$atom$resno==resno & pdb$atom$elety==elements[3], c('x', 'y', 'z')])
    
    
    xa_matrix <- rbind(xpdbca1,xpdbca2, xpdbca3)
    
    
    xbangle<-compas::atomdeps[[rtype]]$bangle[j]        
    xblength<-compas::atomdeps[[rtype]]$blength[j]
    if (is.na(compas::atomdeps[[rtype]]$tangle[j])) {
      chicount <- chicount + 1
      xbchi <- chi[chicount]  
    } else {
      xbchi <- compas::atomdeps[[rtype]]$tangle[j]
    }
    
    
    result<-calCo(xa_matrix, xblength, xbangle, xbchi)
    pdb$atom[pdb$atom$resno==resno & pdb$atom$elety==xnew[j], c('x', 'y', 'z')]<-result
    
    
    #print(result)
  }
  
  return(pdb)
}


