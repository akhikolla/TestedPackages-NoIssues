#' Sample Protein Conformation 1
#' 
#' A sample protein conformation in PDB format.  
#' Excerpted from the native structure of CASP12 target 879, residues 5-24
#' 
#' @format A PDB object, read using \link[bio3d]{read.pdb}
#' @source CASP12 data archive on Prediction Center: \url{http://predictioncenter.org/download_area/CASP12/targets/casp12.targets_T0.releaseDec022016.tgz}

"nat879"

#' Sample Protein Conformation 2
#' 
#' A sample protein conformation in PDB format.  
#' Excerpted from a structure prediction of CASP12 target 879, residues 5-24
#' 
#' @format A PDB object, read using \link[bio3d]{read.pdb}
#' @source CASP12 data archive on Prediction Center: \url{http://predictioncenter.org/download_area/CASP12/targets/casp12.targets_TR.releaseDec022016.tgz}

"pred879"


#' Atom Type Table
#' 
#' Table listing the 167 standard atom types in protein structures (amino acid type and atom identifier pairs).
#' Hydrogen atoms are not considered.
#' 
#' @format A two-column data frame, with "resid" providing the 3-letter amino acid abbreviation and "atomid" 
#' providing the component atoms of each amino acid.
#' @source PDB ATOM entry: \url{http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM}

"atomtype"

#' Atom parameters and definitions for side chains
#' 
#' List of the atom dependencies and typical bond parameters in protein structures for side chains of the
#' 20 standard amino acid types.
#' 
#' @format A list for the 20 amino acid types, each having the attributes
#' \itemize{
#'  \item{"tangle"}: {Dihedral angle defined by atoms A-B-C-D as described below. \code{NA}'s represent the free side chain dihedral angles \emph{chi} for that amino acid.}
#'  \item{"names"}: {The component atoms of the amino acid side chain.  Represents atom D in the dihedral.}
#'  \item{"matx"}: {The names of atoms A,B,C in the dihedral, with rows corresponding to the atoms in "names".}
#'  \item{"bangle"}: {Planar bond angle formed between B-C-D}
#'  \item{"blength"}: {Bond length between C-D}
#' }
#'
#' @references Engh, Richard A., and Robert Huber. "Accurate bond and angle parameters for X-ray protein structure refinement."
#' \emph{Acta Crystallographica} Section A 47.4 (1991): 392-400.

"atomdeps"
