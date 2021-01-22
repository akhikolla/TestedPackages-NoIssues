#' DFIRE Energy Evaluation for Protein Conformations
#' @description 
#' Calculates the energy of a protein conformation using the DFIRE potential.
#' @usage  
#' dfireE(pdb)
#' @param pdb A PDB object, read using \link[bio3d]{read.pdb}.
#'
#' @return
#' Returns the DFIRE energy.
#'
#'
#' @examples
#' dfireE(nat879)
#' 
#' @references 
#' Zhou, Hongyi, and Yaoqi Zhou. "Distance-scaled, finite ideal-gas reference state improves structure-derived potentials of mean force for structure selection and stability prediction." \emph{Protein science} 11.11 (2002): 2714-2726.
#' 
#' @export

dfireE <- function(pdb) {
  
  atypevec <- c()
  for (i in 1:nrow(pdb$atom)) {
    atypevec[i] <- which(compas::atomtype[,1] == pdb$atom$resid[i] & compas::atomtype[,2] ==  pdb$atom$elety[i], arr.ind=T)-1
  }

  dfeval(pdb$atom,atypevec,dfiretable)
  
}
