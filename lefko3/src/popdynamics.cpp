#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//' Vectorize Matrix for Historical Mean Matrix Estimation
//' 
//' Function \code{flagrantcrap()} vectorizes core indices of matrices
//' input as list elements.
//' 
//' @param Xmat A matrix originally a part of a list object.
//' @param allindices A vector of indices to remove from the matrix
//' 
//' @return A column vector of certain elements from the input matrix.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::vec flagrantcrap(arma::mat Xmat, arma::uvec allindices) {
  
  arma::vec newcol = Xmat.elem(allindices);
  
  return newcol;
}

//' Vectorize Matrix for Ahistorical Mean Matrix Estimation
//' 
//' Function \code{moreflagrantcrap()} vectorizes matrices input as list elements.
//' 
//' @param Xmat A matrix originally a part of a list object.
//' 
//' @return A column vector of the input matrix.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::vec moreflagrantcrap(arma::mat Xmat) {
  
  arma::vec newcol = arma::vectorise(Xmat);
  
  return newcol;
}

//' Estimates Mean LefkoMat Object for Historical MPM
//' 
//' Function \code{turbogeodiesel()} estimates mean historical population projection matrices,
//' treating the mean as element-wise arithmetic.
//' 
//' @param loy A data frame denoting the population, patch, and time step designation
//' of each matrix. Includes a total of 9 variables.
//' @param Umats A matrix with all U matrices turned into columns.
//' @param Fmats A matrix with all F matrices turned into columns.
//' @param stages This is the core stageframe held by \code{mats}, equivalent
//' to \code{ahstages}.
//' @param hstages This is the \code{hstages} object held by \code{mats}.
//' @param modelqc This is the \code{modelqc} or \code{dataqc} portion of \code{mats}.
//' @param patchmats A logical value stating whether to estimate patch-level means.
//' @param popmats A logical value stating whether to estimate population-level means.
//' 
//' @return A list using ther basic blueprint of a lefkoMat object.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List turbogeodiesel(DataFrame loy, List Umats, List Fmats, DataFrame stages, DataFrame hstages,
                    DataFrame modelqc, bool patchmats, bool popmats) {
  
  arma::uvec pops = loy["pop"];
  arma::uvec patches = loy["patch"];
  arma::uvec year2 = loy["year2"];
  arma::uvec poppatchc = loy["poppatchc"];
  arma::uvec patchesinpop = loy["patchesinpop"];
  arma::uvec yearsinpatch = loy["yearsinpatch"];
  arma::uvec uniquepops = unique(pops);
  arma::uvec uniquepoppatches = unique(poppatchc);
  int loydim = pops.n_elem;
  int numofpops = uniquepops.n_elem;
  int numofpatches = uniquepoppatches.n_elem;
  
  arma::uvec poporderlong(loydim);
  arma::uvec patchorderlong(loydim);
  arma::uvec annmatriceslong(loydim);
  arma::uvec meanassign(loydim);
  poporderlong.zeros();
  patchorderlong.zeros();
  annmatriceslong.zeros();
  meanassign.zeros();
  
  poporderlong(0) = pops(0);
  patchorderlong(0) = patches(0);
  annmatriceslong(0) = 1;
  meanassign(0) = 1;
  
  int counter {0};
  
  // In this beginning bit, we assess how many mean matrices we will need, and the overall order of means
  if (loydim > 1) {
    for (int i = 1; i < loydim; i++) {
      
      if (poppatchc(i) != poppatchc(i-1)) {
        
        counter++;
        poporderlong(counter) = pops(i);
        patchorderlong(counter) = patches(i);
        annmatriceslong(counter) = 1;
        meanassign(i) = meanassign(i-1) + 1;
        
      } else {
        
        annmatriceslong(counter) = annmatriceslong(counter) + 1;
        meanassign(i) = meanassign(i-1);
        
      }
    }
  }
  
  arma::uvec toestimate = find(poporderlong);
  
  int totalmatrices = toestimate.n_elem + numofpops;
  
  if (patchmats == 1 && popmats == 0) {
    totalmatrices = toestimate.n_elem;
  } else if (patchmats == 0 && popmats == 1) {
    totalmatrices = numofpops;
  }
  
  arma::uvec poporder = poporderlong.elem(toestimate);
  arma::uvec patchorder = patchorderlong.elem(toestimate);
  arma::uvec annmatrices = annmatriceslong.elem(toestimate);
  
  // This next chunk predicts which elements will be targeted for arithmetic mean estimation
  arma::uvec astages = stages["stage_id"];
  int numstages = astages.n_elem;
  
  arma::uvec hstage3in = hstages["stage_id_2"];
  arma::uvec hstage2nin = hstages["stage_id_1"];
  int numhstages = hstage3in.n_elem;
  
  int predictedsize = numstages * numstages * numstages;
  
  arma::uvec hsindexl(predictedsize);
  hsindexl.zeros();
  
  counter = 0;
  
  for (int i1 = 0; i1 < numhstages; i1++) {
    for (int i2 = 0; i2 < numhstages; i2++) {
      if (hstage3in(i1) == hstage2nin(i2)) {
        
        hsindexl(counter) = (i1 * numhstages) + i2;
        
        counter++;
      }
    }
  }
  
  arma::uvec hsgood = find(hsindexl);
  arma::uvec hsindex = hsindexl.elem(hsgood);
  arma::uvec zerovec(1);
  zerovec.zeros();
  arma::uvec allindices = join_cols(zerovec, hsindex);
  
  // Now we build a U and F matrices of element-wise arithmetic means, where
  // each column corresponds to the predicted non-zero elements of each mean
  // matrix, and each matrix is presented as a column vector within the 
  // overall matrix. The A matrix is the sum of U and F.
  
  int core_elem = counter;
  
  arma::mat umatvec(core_elem, totalmatrices);
  arma::mat fmatvec(core_elem, totalmatrices);
  umatvec.zeros();
  fmatvec.zeros();
  
  int patchchoice {0};
  int popchoice {0};
  
  for (int i = 0; i < loydim; i ++) {
    
    if (patchmats == 1) {
      
      patchchoice = poppatchc(i);
      
      umatvec.col(patchchoice) = umatvec.col(patchchoice) + (flagrantcrap(Umats[i], allindices) / yearsinpatch(i));
      fmatvec.col(patchchoice) = fmatvec.col(patchchoice) + (flagrantcrap(Fmats[i], allindices) / yearsinpatch(i));
      
    }
    
    if (popmats == 1) {
      if (patchmats == 1) {
        
        popchoice = numofpatches + pops(i) - 1;
        
      } else {
        
        popchoice = pops(i) - 1;
        
      }
      
      umatvec.col(popchoice) = umatvec.col(popchoice) + (flagrantcrap(Umats[i], allindices) / (yearsinpatch(i) * patchesinpop(i)));
      fmatvec.col(popchoice) = fmatvec.col(popchoice) + (flagrantcrap(Fmats[i], allindices) / (yearsinpatch(i) * patchesinpop(i)));
      
    }
  }
  
  arma::mat amatvec = umatvec + fmatvec;
  
  // Here we create the cheat sheet algorithm
  arma::uvec poporder_redone = join_cols(poporder, uniquepops);
  
  arma::uvec newpatchones(numofpops);
  newpatchones.zeros();
  arma::uvec patchorder_redone = join_cols(patchorder, newpatchones);
  
  if (patchmats == 1 && popmats == 0) {
    poporder_redone = poporder;
    patchorder_redone = patchorder;
  } else if (patchmats == 0 && popmats == 1) {
    poporder_redone = uniquepops;
    patchorder_redone = newpatchones;
  }
  
  DataFrame cheatsheet = DataFrame::create(Named("pop") = poporder_redone, _["patch"] = patchorder_redone);
  
  // Now we will create the main list objects holding the matrices
  arma::mat umat_base(numhstages, numhstages);
  arma::mat fmat_base(numhstages, numhstages);
  arma::mat amat_base(numhstages, numhstages);
  
  umat_base.zeros();
  fmat_base.zeros();
  amat_base.zeros();
  
  umat_base.elem(allindices) = umatvec.col(0);
  fmat_base.elem(allindices) = fmatvec.col(0);
  amat_base.elem(allindices) = amatvec.col(0);
  
  List U = List::create(umat_base);
  List F = List::create(fmat_base);
  List A = List::create(amat_base);
  
  if (totalmatrices > 1) {
    for (int i = 1; i < totalmatrices; i++) {
      umat_base.zeros();
      fmat_base.zeros();
      amat_base.zeros();
      
      umat_base.elem(allindices) = umatvec.col(i);
      fmat_base.elem(allindices) = fmatvec.col(i);
      amat_base.elem(allindices) = amatvec.col(i);
      
      U.push_back(umat_base);
      F.push_back(fmat_base);
      A.push_back(amat_base);
    }
  }
  
  // Matrix QC output
  
  arma::uvec utrans = find(umatvec);
  arma::uvec ftrans = find(fmatvec);
  int totalutrans = utrans.n_elem;
  int totalftrans = ftrans.n_elem;
  
  arma::vec matrixqc(3);
  matrixqc(0) = totalutrans; // summed number of u transitions
  matrixqc(1) = totalftrans; // summed number of f transitions
  matrixqc(2) = totalmatrices;
  
  
  // Final output
  
  List output = List::create(Named("A") = A, _["U"] = U, _["F"] = F, _["hstages"] = hstages,
                             _["ahstages"] = stages, _["labels"] = cheatsheet, _["matrixqc"] = matrixqc,
                             _["modelqc"] = modelqc);
  
  return output;
}

//' Estimates Mean LefkoMat Object for Ahistorical MPM
//' 
//' Function \code{geodiesel()} estimates mean ahistorical population projection matrices,
//' treating the mean as element-wise arithmetic.
//' 
//' @param loy A data frame denoting the population, patch, and time step designation
//' of each matrix. Includes a total of 9 variables.
//' @param Umats A matrix with all U matrices turned into columns.
//' @param Fmats A matrix with all F matrices turned into columns.
//' @param stages This is the core stageframe held by \code{mats}, equivalent
//' to \code{ahstages}.
//' @param modelqc This is the \code{modelqc} or \code{dataqc} portion of \code{mats}.
//' @param patchmats A logical value stating whether to estimate patch-level means.
//' @param popmats A logical value stating whether to estimate population-level means.
//' 
//' @return A list using th4 basic blueprint of a LefkoMat object.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List geodiesel(DataFrame loy, List Umats, List Fmats, DataFrame stages, DataFrame modelqc, bool patchmats, bool popmats) {
  
  arma::uvec pops = loy["pop"];
  arma::uvec patches = loy["patch"];
  arma::uvec year2 = loy["year2"];
  arma::uvec poppatchc = loy["poppatchc"];
  arma::uvec patchesinpop = loy["patchesinpop"];
  arma::uvec yearsinpatch = loy["yearsinpatch"];
  arma::uvec uniquepops = unique(pops);
  arma::uvec uniquepoppatches = unique(poppatchc);
  int loydim = pops.n_elem;
  int numofpops = uniquepops.n_elem;
  int numofpatches = uniquepoppatches.n_elem;
  
  arma::uvec poporderlong(loydim);
  arma::uvec patchorderlong(loydim);
  arma::uvec annmatriceslong(loydim);
  arma::uvec meanassign(loydim);
  poporderlong.zeros();
  patchorderlong.zeros();
  annmatriceslong.zeros();
  meanassign.zeros();
  
  poporderlong(0) = pops(0);
  patchorderlong(0) = patches(0);
  annmatriceslong(0) = 1;
  meanassign(0) = 1;
  
  int counter {0};
  
  // In this beginning bit, we assess how many mean matrices we will need, and the overall order of means
  if (loydim > 1) {
    for (int i = 1; i < loydim; i++) {
      
      if (poppatchc(i) != poppatchc(i-1)) {
        
        counter++;
        poporderlong(counter) = pops(i);
        patchorderlong(counter) = patches(i);
        annmatriceslong(counter) = 1;
        meanassign(i) = meanassign(i-1) + 1;
        
      } else {
        
        annmatriceslong(counter) = annmatriceslong(counter) + 1;
        meanassign(i) = meanassign(i-1);
        
      }
    }
  }
  
  arma::uvec toestimate = find(poporderlong);
  
  int totalmatrices = toestimate.n_elem + numofpops;
  
  if (patchmats == 1 && popmats == 0) {
    totalmatrices = toestimate.n_elem;
  } else if (patchmats == 0 && popmats == 1) {
    totalmatrices = numofpops;
  }
  
  arma::uvec poporder = poporderlong.elem(toestimate);
  arma::uvec patchorder = patchorderlong.elem(toestimate);
  arma::uvec annmatrices = annmatriceslong.elem(toestimate);
  
  // This next chunk predicts which elements will be targeted for arithmetic mean estimation
  arma::uvec astages = stages["stage_id"];
  int numstages = astages.n_elem;
  
  // Now we build a U and F matrices of element-wise arithmetic means, where
  // each column corresponds to the predicted non-zero elements of each mean
  // matrix, and each matrix is presented as a column vector within the 
  // overall matrix. The A matrix is the sum of U and F.
  
  int core_elem = numstages * numstages;
  
  arma::mat umatvec(core_elem, totalmatrices);
  arma::mat fmatvec(core_elem, totalmatrices);
  umatvec.zeros();
  fmatvec.zeros();
  
  int patchchoice {0};
  int popchoice {0};
  
  for (int i = 0; i < loydim; i ++) {
    
    if (patchmats == 1) {
      
      patchchoice = poppatchc(i);
      
      umatvec.col(patchchoice) = umatvec.col(patchchoice) + (moreflagrantcrap(Umats[i]) / yearsinpatch(i));
      fmatvec.col(patchchoice) = fmatvec.col(patchchoice) + (moreflagrantcrap(Fmats[i]) / yearsinpatch(i));
      
    }
    
    if (popmats == 1) {
      if (patchmats == 1) {
        
        popchoice = numofpatches + pops(i) - 1;
        
      } else {
        
        popchoice = pops(i) - 1;
        
      }
      
      umatvec.col(popchoice) = umatvec.col(popchoice) + (moreflagrantcrap(Umats[i]) / (yearsinpatch(i) * patchesinpop(i)));
      fmatvec.col(popchoice) = fmatvec.col(popchoice) + (moreflagrantcrap(Fmats[i]) / (yearsinpatch(i) * patchesinpop(i)));
      
    }
  }
  
  arma::mat amatvec = umatvec + fmatvec;
  
  // Here we create the cheat sheet algorithm
  arma::uvec poporder_redone = join_cols(poporder, uniquepops);
  
  arma::uvec newpatchones(numofpops);
  newpatchones.zeros();
  arma::uvec patchorder_redone = join_cols(patchorder, newpatchones);
  
  if (patchmats == 1 && popmats == 0) {
    poporder_redone = poporder;
    patchorder_redone = patchorder;
  } else if (patchmats == 0 && popmats == 1) {
    poporder_redone = uniquepops;
    patchorder_redone = newpatchones;
  }
  
  DataFrame cheatsheet = DataFrame::create(Named("pop") = poporder_redone, _["patch"] = patchorder_redone);
  
  // Now we will create the main list objects holding the matrices
  
  arma::mat umat_base = umatvec.col(0);
  umat_base.reshape(numstages, numstages);
  
  arma::mat fmat_base = fmatvec.col(0);
  fmat_base.reshape(numstages, numstages);
  
  arma::mat amat_base = amatvec.col(0);
  amat_base.reshape(numstages, numstages);
  
  List U = List::create(umat_base);
  List F = List::create(fmat_base);
  List A = List::create(amat_base);
  
  if (totalmatrices > 1) {
    for (int i = 1; i < totalmatrices; i++) {
      umat_base.zeros();
      fmat_base.zeros();
      amat_base.zeros();
      
      umat_base = umatvec.col(i);
      fmat_base = fmatvec.col(i);
      amat_base = amatvec.col(i);
      
      umat_base.reshape(numstages, numstages);
      fmat_base.reshape(numstages, numstages);
      amat_base.reshape(numstages, numstages);
      
      U.push_back(umat_base);
      F.push_back(fmat_base);
      A.push_back(amat_base);
    }
  }
  
  // Matrix QC output
  
  arma::uvec utrans = find(umatvec);
  arma::uvec ftrans = find(fmatvec);
  int totalutrans = utrans.n_elem;
  int totalftrans = ftrans.n_elem;
  
  arma::vec matrixqc(3);
  matrixqc(0) = totalutrans; // summed number of u transitions
  matrixqc(1) = totalftrans; // summed number of f transitions
  matrixqc(2) = totalmatrices;
  
  
  // Final output
  
  List output = List::create(Named("A") = A, _["U"] = U, _["F"] = F, _["ahstages"] = stages, 
                             _["labels"] = cheatsheet, _["matrixqc"] = matrixqc, _["modelqc"] = modelqc);
  
  return output;
}

//' Complete Full Eigen Analysis of a Single Dense Matrix
//' 
//' \code{decomp3()} returns all eigenvalues, right eigenvectors, and left
//' eigenvectors estimated for a matrix by the \code{eig_gen}() function
//' in the C++ Armadillo library. Works with dense matrices.
//' 
//' @param Amat A population projection matrix of class \code{matrix}.
//'
//' @return This function returns all estimated eigenvalues, right
//' eigenvectors, and left eigenvectors of a single matrix. This output is
//' provided as a list with three parts, named appropriately.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List decomp3(arma::mat Amat) {
  arma::cx_vec Aeigval;
  arma::cx_mat Aeigvecl;
  arma::cx_mat Aeigvecr;
  
  eig_gen(Aeigval, Aeigvecl, Aeigvecr, Amat);
  
  List output = List::create(Named("eigenvalues") = Aeigval, _["left_eigenvectors"] = Aeigvecl,
                             _["right_eigenvectors"] = Aeigvecr);
  
  return output;
}

//' Complete Full Eigen Analysis of a Single Sparse Matrix
//' 
//' \code{decomp3sp()} returns all eigenvalues, right eigenvectors, and left
//' eigenvectors estimated for a matrix by the \code{eigs_gen}() function
//' in the C++ Armadillo library. Works with sparse matrices.
//' 
//' @param Amat A population projection matrix of class \code{matrix}.
//'
//' @return This function returns all estimated eigenvalues, right
//' eigenvectors, and left eigenvectors of a single matrix. This output is
//' provided as a list with three parts, named appropriately.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List decomp3sp(arma::mat Amat) {
  arma::sp_mat spAmat(Amat);
  arma::sp_mat t_spAmat = spAmat.t();
  
  arma::cx_vec Aeigval;
  arma::cx_vec Aeigvall;
  arma::cx_mat Aeigvecl;
  arma::cx_mat Aeigvecr;
  
  eigs_gen(Aeigval, Aeigvecr, spAmat, 1);
  
  eigs_gen(Aeigvall, Aeigvecl, t_spAmat, 1);
  
  List output = List::create(Named("eigenvalues") = Aeigval, _["left_eigenvectors"] = Aeigvecl,
                             _["right_eigenvectors"] = Aeigvecr);
  
  return output;
}

//' Estimate Deterministic Population Growth Rate of a Dense Matrix
//' 
//' \code{lambda3matrix()} returns the dominant eigenvalue of a single
//' dense projection matrix.
//' 
//' @param Amat A population projection matrix of class \code{matrix}.
//'
//' @return This function returns the dominant eigenvalue of the matrix. This
//' is given as the largest real part of all eigenvalues estimated via the 
//' \code{eig_gen}() function in the C++ Armadillo library.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
double lambda3matrix(arma::mat Amat) {
  List eigenstuff = decomp3(Amat);
  
  cx_vec Eigenvals = eigenstuff["eigenvalues"];
  arma::vec realeigenvals = real(Eigenvals);
  
  double lambda = max(realeigenvals);
  
  return lambda;
}

//' Estimate Deterministic Population Growth Rate of a Sparse Matrix
//' 
//' \code{lambda3matrixsp()} returns the dominant eigenvalue of a single
//' sparse projection matrix. This function can handle large and sparse 
//' matrices, and so can be used with large historical matrices, IPMs, 
//' age x stage matrices, as well as smaller ahistorical matrices.
//' 
//' @param Amat A population projection matrix of class \code{matrix}.
//'
//' @return This function returns the dominant eigenvalue of the matrix. This
//' is given as the largest real part of all eigenvalues estimated via the 
//' \code{eigs_gen}() function in the C++ Armadillo library.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
double lambda3matrixsp(arma::mat Amat) {
  List eigenstuff = decomp3sp(Amat);
  
  cx_vec Eigenvals = eigenstuff["eigenvalues"];
  arma::vec realeigenvals = real(Eigenvals);
  
  double lambda = max(realeigenvals);
  
  return lambda;
}

//' Estimate Stable Stage Distribution for a Dense Population Matrix
//' 
//' \code{ss3matrix()} returns the stable stage distribution for a 
//' dense population matrix.
//' 
//' @param Amat A population projection matrix of class \code{matrix}.
//' 
//' @return This function returns the stable stage distribution corresponding to
//' the input matrix. The stable stage distribution is given as the right 
//' eigenvector associated with largest real part of the eigenvalues estimated 
//' for the matrix via the \code{eig_gen}() function in the C++ Armadillo 
//' library, divided by the sum of the associated right eigenvector. 
//' 
//' @seealso \code{\link{stablestage3}()}
//' @seealso \code{\link{stablestage3.lefkoMat}()}
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::vec ss3matrix(arma::mat Amat) {
  List eigenstuff = decomp3(Amat);
  
  cx_vec Eigenvals = eigenstuff["eigenvalues"];
  arma::vec realeigenvals = real(Eigenvals);
  
  int lambda1 = realeigenvals.index_max();
  
  cx_mat rightvec = eigenstuff["right_eigenvectors"];
  arma::vec realrightvec = real(rightvec.col(lambda1));
  
  double rvsum = sum(realrightvec);
  int rvel = realrightvec.n_elem;
  
  arma::vec wcorr (rvel);
  
  for (int i = 0; i < rvel; i++) {
    wcorr(i) = realrightvec(i) / rvsum;
  }
  
  return wcorr;
}

//' Estimate Stable Stage Distribution for a Sparse Population Matrix
//' 
//' \code{ss3matrixsp()} returns the stable stage distribution for a 
//' sparse population matrix. This function can handle large and sparse 
//' matrices, and so can be used with large historical matrices, IPMs, 
//' age x stage matrices, as well as smaller ahistorical matrices.
//' 
//' @param Amat A population projection matrix of class \code{matrix}.
//' 
//' @return This function returns the stable stage distribution corresponding to
//' the input matrix. The stable stage distribution is given as the right 
//' eigenvector associated with largest real part of the eigenvalues estimated 
//' for the matrix via the \code{eigs_gen}() function in the C++ Armadillo 
//' library, divided by the sum of the associated right eigenvector. 
//' 
//' @seealso \code{\link{stablestage3}()}
//' @seealso \code{\link{stablestage3.lefkoMat}()}
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::vec ss3matrixsp(arma::mat Amat) {
  List eigenstuff = decomp3sp(Amat);
  
  cx_vec Eigenvals = eigenstuff["eigenvalues"];
  arma::vec realeigenvals = real(Eigenvals);
  
  int lambda1 = realeigenvals.index_max();
  
  cx_mat rightvec = eigenstuff["right_eigenvectors"];
  arma::vec realrightvec = real(rightvec.col(lambda1));
  realrightvec.clean(0.00000000000001); // This line replaces all numbers lower than 1 x 10-14 with 0
  
  double rvsum = sum(realrightvec);
  int rvel = realrightvec.n_elem;
  
  arma::vec wcorr (rvel);
  
  for (int i = 0; i < rvel; i++) {
    wcorr(i) = realrightvec(i) / rvsum;
  }
  
  return wcorr;
}

//' Estimate Reproductive Value for a Dense Population Matrix
//' 
//' \code{rv3matrix()} returns the reproductive values for stages in a 
//' dense population matrix. The function provides standard reproductive 
//' values, meaning that the overall reproductive values of basic life 
//' history stages in a historical matrix are not provided (the 
//' \code{\link{repvalue3.lefkoMat}()} function estimates these on the basis 
//' of stage description information provided in the \code{lefkoMat} object 
//' used as input in that function).
//' 
//' @param Amat A population projection matrix.
//' 
//' @return This function returns a vector characterizing the 
//' reproductive values for stages of a population projection matrix. This is 
//' given as the left eigenvector associated with largest real part of the
//' dominant eigenvalue estimated via the \code{eig_gen}() function in the C++ 
//' Armadillo library, divided by the first non-zero element of the left 
//' eigenvector. 
//' 
//' @seealso \code{\link{repvalue3}()}
//' @seealso \code{\link{repvalue3.lefkoMat}()}
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::vec rv3matrix(arma::mat Amat) {
  List eigenstuff = decomp3(Amat);
  
  cx_vec Eigenvals = eigenstuff["eigenvalues"];
  arma::vec realeigenvals = real(Eigenvals);
  
  int lambda1 = realeigenvals.index_max();
  
  cx_mat leftvec = eigenstuff["left_eigenvectors"];
  arma::vec realleftvec = real(leftvec.col(lambda1));
  arma::vec rlvabs = abs(realleftvec);
  
  // This section identifies the first non-zero element of the reproductive value vector
  arma::uvec rlvabsalt = find(rlvabs);
  int rlvminelem = rlvabsalt(0);
  double rlvmin = realleftvec(rlvminelem);
  
  int lvel = realleftvec.n_elem;
  
  arma::vec vcorr (lvel);
  
  for (int i = 0; i < lvel; i++) {
    vcorr(i) = realleftvec(i) / rlvmin;
  }
  
  return vcorr;
}

//' Estimate Reproductive Value for a Sparse Population Matrix
//' 
//' \code{rv3matrixsp()} returns the reproductive values for stages in a 
//' sparse population matrix. The function provides standard reproductive 
//' values, meaning that the overall reproductive values of basic life 
//' history stages in a historical matrix are not provided (the 
//' \code{\link{repvalue3.lefkoMat}()} function estimates these on the basis 
//' of stage description information provided in the \code{lefkoMat} object 
//' used as input in that function). This function can handle large and 
//' sparse matrices, and so can be used with large historical matrices, IPMs, 
//' age x stage matrices, as well as smaller ahistorical 
//' matrices.
//' 
//' @param Amat A population projection matrix.
//' 
//' @return This function returns a vector characterizing the 
//' reproductive values for stages of a population projection matrix. This is 
//' given as the left eigenvector associated with largest real part of the
//' dominant eigenvalue estimated via the \code{eigs_gen}() function in the C++ 
//' Armadillo library, divided by the first non-zero element of the left 
//' eigenvector. 
//' 
//' @seealso \code{\link{repvalue3}()}
//' @seealso \code{\link{repvalue3.lefkoMat}()}
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::vec rv3matrixsp(arma::mat Amat) {
  List eigenstuff = decomp3sp(Amat);
  
  cx_vec Eigenvals = eigenstuff["eigenvalues"];
  arma::vec realeigenvals = real(Eigenvals);
  
  int lambda1 = realeigenvals.index_max();
  
  cx_mat leftvec = eigenstuff["left_eigenvectors"];
  arma::vec realleftvec = real(leftvec.col(lambda1));
  realleftvec.clean(0.00000000000001); // This line replaces all numbers lower than 1 x 10-14 with 0
  arma::vec rlvabs = abs(realleftvec);
  
  // This section identifies the first non-zero element of the reproductive value vector
  arma::uvec rlvabsalt = find(rlvabs);
  int rlvminelem = rlvabsalt(0);
  double rlvmin = realleftvec(rlvminelem);
  
  int lvel = realleftvec.n_elem;
  
  arma::vec vcorr (lvel);
  
  for (int i = 0; i < lvel; i++) {
    vcorr(i) = realleftvec(i) / rlvmin;
  }
  
  return vcorr;
}

//' Estimate Sensitivities for a Dense Population Matrix
//' 
//' \code{sens3matrix()} returns the sensitivity of lambda with respect
//' to each element in a dense matrix. This is accomplished via the
//' \code{eig_gen}() function in the C++ Armadillo library.
//' 
//' @param Amat A population projection matrix.
//' 
//' @return This function returns a matrix of sensitivities. 
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::mat sens3matrix(arma::mat Amat) {
  List eigenstuff = decomp3(Amat);
  
  cx_vec Eigenvals = eigenstuff["eigenvalues"];
  arma::vec realeigenvals = real(Eigenvals);
  
  int lambda1 = realeigenvals.index_max();
  
  // This is the w vector
  cx_mat rightvec = eigenstuff["right_eigenvectors"];
  arma::vec realrightvec = real(rightvec.col(lambda1));
  
  double rvsum = sum(realrightvec);
  int rvel = realrightvec.n_elem;
  
  // This is the v vector
  cx_mat leftvec = eigenstuff["left_eigenvectors"];
  arma::vec realleftvec = real(leftvec.col(lambda1));
  arma::vec rlvabs = abs(realleftvec);
  
  arma::uvec rlvabsalt = find(rlvabs);
  int rlvminelem = rlvabsalt(0);
  double rlvmin = realleftvec(rlvminelem);
  
  arma::vec wcorr (rvel);
  arma::vec vcorr (rvel);
  arma::vec vwprod (rvel);
  arma::mat smat (rvel, rvel);
  smat.zeros();
  
  // This loop and the following line create the scalar product vw
  for (int i = 0; i < rvel; i++) {
    wcorr(i) = realrightvec(i) / rvsum;
    vcorr(i) = realleftvec(i) / rlvmin;
    
    vwprod(i) = wcorr(i) * vcorr(i);
  }
  
  double vwscalar = sum(vwprod);
  
  // This loop populates the sensitivity matrix
  for (int i = 0; i < rvel; i++) {
    for (int j = 0; j < rvel; j++) {
      
      smat(i, j) = vcorr(i) * wcorr(j) / vwscalar;
    }
  }
  
  return smat;
}

//' Estimate Sensitivities for a Sparse Population Matrix
//' 
//' \code{sens3matrixsp()} returns the sensitivity of lambda with respect
//' to each element in a sparse matrix. This is accomplished via the
//' \code{eig_gen}() function in the C++ Armadillo library.
//' 
//' @param Amat A population projection matrix.
//' 
//' @return This function returns a matrix of sensitivities. 
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::mat sens3matrixsp(arma::mat Amat) {
  List eigenstuff = decomp3sp(Amat);
  
  cx_vec Eigenvals = eigenstuff["eigenvalues"];
  arma::vec realeigenvals = real(Eigenvals);
  
  int lambda1 = realeigenvals.index_max();
  
  // This is the w vector
  cx_mat rightvec = eigenstuff["right_eigenvectors"];
  arma::vec realrightvec = real(rightvec.col(lambda1));
  realrightvec.clean(0.00000000000001); // This line replaces all numbers lower than 1 x 10-14 with 0
  
  double rvsum = sum(realrightvec);
  int rvel = realrightvec.n_elem;
  
  // This is the v vector
  cx_mat leftvec = eigenstuff["left_eigenvectors"];
  arma::vec realleftvec = real(leftvec.col(lambda1));
  realleftvec.clean(0.00000000000001); // This line replaces all numbers lower than 1 x 10-14 with 0
  arma::vec rlvabs = abs(realleftvec);
  
  arma::uvec rlvabsalt = find(rlvabs);
  int rlvminelem = rlvabsalt(0);
  double rlvmin = realleftvec(rlvminelem);
  
  arma::vec wcorr (rvel);
  arma::vec vcorr (rvel);
  arma::vec vwprod (rvel);
  arma::mat smat (rvel, rvel);
  smat.zeros();
  
  // This loop and the following line create the scalar product vw
  for (int i = 0; i < rvel; i++) {
    wcorr(i) = realrightvec(i) / rvsum;
    vcorr(i) = realleftvec(i) / rlvmin;
    
    vwprod(i) = wcorr(i) * vcorr(i);
  }
  
  double vwscalar = sum(vwprod);
  
  // This loop populates the sensitivity matrix
  for (int i = 0; i < rvel; i++) {
    for (int j = 0; j < rvel; j++) {
      
      smat(i, j) = vcorr(i) * wcorr(j) / vwscalar;
    }
  }
  
  return smat;
}

//' Estimate Sensitivities for a Historical LefkoMat Object
//' 
//' \code{sens3hlefko()} returns the sensitivity of lambda with respect
//' to each historical stage-pair in the matrix, and the associated
//' sensitivity for each life history stage. This is accomplished via the 
//' \code{eigs_gen}() function in the C++ Armadillo library.
//' 
//' @param Amat A population projection matrix.
//' @param ahstages An integar vector of unique ahistorical stages.
//' @param hstages An integar vector of unique historical stage pairs.
//' 
//' @return This function returns a list with two sensitivity matrices:
//' \item{h_smat}{Matrix of sensitivities corresponding to the historical matrix.}
//' \item{ah_smat}{Matrix of sensitivities corresponding to the ahistorical
//' matrix.}
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List sens3hlefko(arma::mat Amat, DataFrame ahstages, DataFrame hstages) {
  
  arma::uvec stage_id = ahstages["stage_id"];
  arma::uvec h_stage_2 = hstages["stage_id_2"];
  arma::uvec h_stage_1 = hstages["stage_id_1"];
  
  List eigenstuff = decomp3sp(Amat);
  
  cx_vec Eigenvals = eigenstuff["eigenvalues"];
  arma::vec realeigenvals = real(Eigenvals);
  
  //double lambda = max(realeigenvals);
  int lambda1 = realeigenvals.index_max();
  
  // This is the w vector
  cx_mat rightvec = eigenstuff["right_eigenvectors"];
  arma::vec realrightvec = real(rightvec.col(lambda1));
  realrightvec.clean(0.00000000000001); // This line replaces all numbers lower than 1 x 10-14 with 0
  
  double rvsum = sum(realrightvec);
  int rvel = realrightvec.n_elem;
  
  // This is the v vector
  cx_mat leftvec = eigenstuff["left_eigenvectors"];
  arma::vec realleftvec = real(leftvec.col(lambda1));
  realleftvec.clean(0.00000000000001); // This line replaces all numbers lower than 1 x 10-14 with 0
  arma::vec rlvabs = abs(realleftvec);
  
  arma::uvec rlvabsalt = find(rlvabs);
  int rlvminelem = rlvabsalt(0);
  double rlvmin = realleftvec(rlvminelem);
  
  int ahstagelength = stage_id.n_elem;
  
  arma::vec wcorr (rvel);
  arma::vec vcorr (rvel);
  arma::vec vwprod (rvel);
  arma::mat smat (rvel, rvel);
  smat.zeros();
  
  arma::vec wcorrah (ahstagelength);
  wcorrah.zeros();
  arma::vec vcorrah (ahstagelength);
  vcorrah.zeros();
  arma::vec vwprodah (ahstagelength);
  vwprodah.zeros();
  arma::mat ahsens(ahstagelength, ahstagelength);
  ahsens.zeros();
  
  int ahrows {0};
  //int ahcols {0};
  
  // This loop and the following line create the scalar product vw and the ahistorical stable stage distribution w
  for (int i = 0; i < rvel; i++) {
    wcorr(i) = realrightvec(i) / rvsum;
    vcorr(i) = realleftvec(i) / rlvmin;
    
    vwprod(i) = wcorr(i) * vcorr(i);
    
    ahrows = h_stage_2(i) - 1;
    
    wcorrah(ahrows) = wcorrah(ahrows) + wcorr(i);
  }
  
  double vwscalar = sum(vwprod);
  
  // This loop and the following line create the scalar product vw and the ahistorical stable stage distribution w
  for (int i = 0; i < rvel; i++) {
    ahrows = h_stage_2(i) - 1;
    
    vcorrah(ahrows) = vcorr(i) * wcorr(i) / wcorrah(ahrows);
    
  }
  
  // This loop populates the historical and ahistorical elasticity matrices
  for (int i = 0; i < rvel; i++) {
    for (int j = 0; j < rvel; j++) {
      
      smat(i, j) = vcorr(i) * wcorr(j) / vwscalar;
      
    }
  }
  
  for (int i = 0; i < ahstagelength; i++) {
    vwprodah(i) = wcorrah(i) * vcorrah(i);
  }
  
  double vwscalarah = sum(vwprodah);
  
  for (int i = 0; i < ahstagelength; i++) {
    for (int j = 0; j < ahstagelength; j++) {
      ahsens(i, j) = vcorrah(i) * wcorrah(j) / vwscalarah;
    }
  }
  
  List output = List::create(Named("h_smat") = smat, _["ah_smat"] = ahsens);
  
  return output;
}

//' Estimate Elasticities for a Dense Population Matrix
//' 
//' \code{elas3matrix()} returns the elasticity of lambda with respect
//' to each element in a dense matrix. This is accomplished via the
//' \code{eig_gen}() function in the C++ Armadillo library.
//' 
//' @param Amat A population projection matrix.
//' 
//' @return This function returns a matrix of elasticities. 
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::mat elas3matrix(arma::mat Amat) {
  List eigenstuff = decomp3(Amat);
  
  cx_vec Eigenvals = eigenstuff["eigenvalues"];
  arma::vec realeigenvals = real(Eigenvals);
  
  double lambda = max(realeigenvals);
  int lambda1 = realeigenvals.index_max();
  
  // This is the w vector
  cx_mat rightvec = eigenstuff["right_eigenvectors"];
  arma::vec realrightvec = real(rightvec.col(lambda1));
  
  double rvsum = sum(realrightvec);
  int rvel = realrightvec.n_elem;
  
  // This is the v vector
  cx_mat leftvec = eigenstuff["left_eigenvectors"];
  arma::vec realleftvec = real(leftvec.col(lambda1));
  arma::vec rlvabs = abs(realleftvec);
  
  arma::uvec rlvabsalt = find(rlvabs);
  int rlvminelem = rlvabsalt(0);
  double rlvmin = realleftvec(rlvminelem);
  
  arma::vec wcorr (rvel);
  arma::vec vcorr (rvel);
  arma::vec vwprod (rvel);
  arma::mat emat (rvel, rvel);
  emat.zeros();
  
  // This loop and the following line create the scalar product vw
  for (int i = 0; i < rvel; i++) {
    wcorr(i) = realrightvec(i) / rvsum;
    vcorr(i) = realleftvec(i) / rlvmin;
    
    vwprod(i) = wcorr(i) * vcorr(i);
  }
  
  double vwscalar = sum(vwprod);
  
  // This loop populates the elasticity matrix
  for (int i = 0; i < rvel; i++) {
    for (int j = 0; j < rvel; j++) {
      
      emat(i, j) = (vcorr(i) * wcorr(j) * Amat(i, j)) / (vwscalar * lambda);
    }
  }
  
  return emat;
}

//' Estimate Elasticities for a Sparse Population Matrix
//' 
//' \code{elas3matrixsp()} returns the elasticity of lambda with respect
//' to each element in a sparse matrix. This is accomplished via the
//' \code{eigs_gen}() function in the C++ Armadillo library.
//' 
//' @param Amat A population projection matrix.
//' 
//' @return This function returns a matrix of elasticities. 
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::mat elas3matrixsp(arma::mat Amat) {
  List eigenstuff = decomp3sp(Amat);
  
  cx_vec Eigenvals = eigenstuff["eigenvalues"];
  arma::vec realeigenvals = real(Eigenvals);
  
  double lambda = max(realeigenvals);
  int lambda1 = realeigenvals.index_max();
  
  // This is the w vector
  cx_mat rightvec = eigenstuff["right_eigenvectors"];
  arma::vec realrightvec = real(rightvec.col(lambda1));
  realrightvec.clean(0.00000000000001); // This line replaces all numbers lower than 1 x 10-14 with 0
  
  double rvsum = sum(realrightvec);
  int rvel = realrightvec.n_elem;
  
  // This is the v vector
  cx_mat leftvec = eigenstuff["left_eigenvectors"];
  arma::vec realleftvec = real(leftvec.col(lambda1));
  realleftvec.clean(0.00000000000001); // This line replaces all numbers lower than 1 x 10-14 with 0
  arma::vec rlvabs = abs(realleftvec);
  
  arma::uvec rlvabsalt = find(rlvabs);
  int rlvminelem = rlvabsalt(0);
  double rlvmin = realleftvec(rlvminelem);
  
  arma::vec wcorr (rvel);
  arma::vec vcorr (rvel);
  arma::vec vwprod (rvel);
  arma::mat emat (rvel, rvel);
  emat.zeros();
  
  // This loop and the following line create the scalar product vw
  for (int i = 0; i < rvel; i++) {
    wcorr(i) = realrightvec(i) / rvsum;
    vcorr(i) = realleftvec(i) / rlvmin;
    
    vwprod(i) = wcorr(i) * vcorr(i);
  }
  
  double vwscalar = sum(vwprod);
  
  // This loop populates the elasticity matrix
  for (int i = 0; i < rvel; i++) {
    for (int j = 0; j < rvel; j++) {
      
      emat(i, j) = (vcorr(i) * wcorr(j) * Amat(i, j)) / (vwscalar * lambda);
    }
  }
  
  return emat;
}

//' Estimate Elasticities for a Historical LefkoMod Object
//' 
//' \code{elas3hlefko()} returns the elasticity of lambda with respect
//' to each historical stage-pair in the matrix, and the summed elasticities
//' for each life history stage. This is accomplished via the \code{eigs_gen}()
//' function in the C++ Armadillo library.
//' 
//' @param Amat A population projection matrix.
//' @param ahstages An integar vector of unique ahistorical stages.
//' @param hstages An integar vector of unique historical stage pairs.
//' 
//' @return This function returns a list with two elasticity matrices:
//' \item{h_emat}{Matrix of elasticities corresponding to the historical matrix.}
//' \item{ah_emat}{Matrix of elasticities corresponding to the ahistorical
//' matrix, but using summed historical elasticities as the basis of estimation.}
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List elas3hlefko(arma::mat Amat, DataFrame ahstages, DataFrame hstages) {
  
  arma::uvec stage_id = ahstages["stage_id"];
  arma::uvec h_stage_2 = hstages["stage_id_2"];
  arma::uvec h_stage_1 = hstages["stage_id_1"];
  
  List eigenstuff = decomp3sp(Amat);
  
  cx_vec Eigenvals = eigenstuff["eigenvalues"];
  arma::vec realeigenvals = real(Eigenvals);
  
  double lambda = max(realeigenvals);
  int lambda1 = realeigenvals.index_max();
  
  // This is the w vector
  cx_mat rightvec = eigenstuff["right_eigenvectors"];
  arma::vec realrightvec = real(rightvec.col(lambda1));
  realrightvec.clean(0.00000000000001); // This line replaces all numbers lower than 1 x 10-14 with 0
  
  double rvsum = sum(realrightvec);
  int rvel = realrightvec.n_elem;
  
  // This is the v vector
  cx_mat leftvec = eigenstuff["left_eigenvectors"];
  arma::vec realleftvec = real(leftvec.col(lambda1));
  realleftvec.clean(0.00000000000001); // This line replaces all numbers lower than 1 x 10-14 with 0
  arma::vec rlvabs = abs(realleftvec);
  
  arma::uvec rlvabsalt = find(rlvabs);
  int rlvminelem = rlvabsalt(0);
  double rlvmin = realleftvec(rlvminelem);
  
  arma::vec wcorr (rvel);
  arma::vec vcorr (rvel);
  arma::vec vwprod (rvel);
  arma::mat emat (rvel, rvel);
  emat.zeros();
  
  // This loop and the following line create the scalar product vw
  for (int i = 0; i < rvel; i++) {
    wcorr(i) = realrightvec(i) / rvsum;
    vcorr(i) = realleftvec(i) / rlvmin;
    
    vwprod(i) = wcorr(i) * vcorr(i);
  }
  
  double vwscalar = sum(vwprod);
  
  // The next few lines set up the empty ahistorical matrix
  int ahstagelength = stage_id.n_elem;
  
  int ahrows {0};
  int ahcols {0};
  
  arma::mat ahelas(ahstagelength, ahstagelength);
  ahelas.zeros();
  
  // This loop populates the historical and ahistorical elasticity matrices
  for (int i = 0; i < rvel; i++) {
    for (int j = 0; j < rvel; j++) {
      
      emat(i, j) = (vcorr(i) * wcorr(j) * Amat(i, j)) / (vwscalar * lambda);
      
      ahrows = h_stage_2(i) - 1;
      ahcols = h_stage_1(i) - 1;
      
      ahelas(ahrows, ahcols) = ahelas(ahrows, ahcols) + emat(i, j);
    }
  }
  
  List output = List::create(Named("h_emat") = emat, _["ah_emat"] = ahelas);
  
  return output;
}

