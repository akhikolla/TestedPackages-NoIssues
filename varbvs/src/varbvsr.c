// Part of the varbvs package, https://github.com/pcarbo/varbvs
//
// Copyright (C) 2012-2018, Peter Carbonetto
//
// This program is free software: you can redistribute it under the
// terms of the GNU General Public License; either version 3 of the
// License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANY; without even the implied warranty of
// MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Public License for more details.
//
#include <R.h>
#include <Rinternals.h>
#include "misc.h"
#include "varbvs.h"

// FUNCTION DEFINITIONS
// ---------------------------------------------------------------------
// This function is used to implement the R function varbvsnormupdate.
// It is called in R using the .Call interface.
SEXP varbvsnormupdate_Call (SEXP Xp, SEXP sigmap, SEXP sap, SEXP logoddsp,
			    SEXP xyp, SEXP dp, SEXP alphap, SEXP mup,
			    SEXP Xrp, SEXP ip) {

  // (1) GET INPUTS AND OUTPUTS
  // --------------------------
  double* X       = REAL(Xp);        // Input X.
  double  sigma   = *REAL(sigmap);   // Input sigma.
  double  sa      = *REAL(sap);      // Input sa.
  double* logodds = REAL(logoddsp);  // Input logodds.
  double* xy      = REAL(xyp);       // Input xy.
  double* d       = REAL(dp);        // Input d.
  double* alpha   = REAL(alphap);    // Input and output alpha.
  double* mu      = REAL(mup);       // Input and output mu.
  double* Xr      = REAL(Xrp);       // Input and output Xr.
  int*    i       = INTEGER(ip);     // Input i.

  // Get the number of samples (n) and the number of coordinate ascent
  // updates (m).
  R_xlen_t n = length(Xrp);
  R_xlen_t m = length(ip);

  // (2) CYCLE THROUGH COORDINATE ASCENT UPDATES
  // -------------------------------------------
  for (R_xlen_t j = 0; j < m; j++) {
    R_xlen_t k = (R_xlen_t) i[j];

    // Get the kth column of matrix X.
    const double* x = getConstColumn(X,k,n);

    // Perform the update.
    varbvsnormupdate(x,xy[k],d[k],sigma,sa,logodds[k],alpha + k,mu + k,Xr,n);
  }

  return R_NilValue;
}

// ---------------------------------------------------------------------
// This function is used to implement the R function varbvsbinupdate.
// It is called in R using the .Call interface.
SEXP varbvsbinupdate_Call (SEXP Xp, SEXP sap, SEXP logoddsp, SEXP dp,
			   SEXP xdxp, SEXP xyp, SEXP xdp, SEXP alphap,
			   SEXP mup, SEXP Xrp, SEXP ip) {

  // (1) GET INPUTS AND OUTPUTS
  // --------------------------
  double* X       = REAL(Xp);        // Input X.
  double  sa      = *REAL(sap);      // Input sa.
  double* logodds = REAL(logoddsp);  // Input logodds.
  double* d       = REAL(dp);        // Input d.
  double* xdx     = REAL(xdxp);      // Input xdx.
  double* xy      = REAL(xyp);       // Input xy.
  double* xd      = REAL(xdp);       // Input xd.
  double* alpha   = REAL(alphap);    // Input and output alpha.
  double* mu      = REAL(mup);       // Input and output mu.
  double* Xr      = REAL(Xrp);       // Input and output Xr.
  int*    i       = INTEGER(ip);     // Input i.

  // Get the number of samples (n) and the number of coordinate ascent
  // updates (m).
  R_xlen_t n = length(Xrp);
  R_xlen_t m = length(ip);

  // (2) CYCLE THROUGH COORDINATE ASCENT UPDATES
  // -------------------------------------------
  for (R_xlen_t j = 0; j < m; j++) {
    R_xlen_t k = (R_xlen_t) i[j];

    // Get the kth column of matrix X.
    const double* x = getConstColumn(X,k,n);

    // Perform the update.
    varbvsbinupdate(x,xy[k],xd[k],xdx[k],d,sa,logodds[k],alpha+k,mu+k,Xr,n);
  }

  return R_NilValue;
}

// ---------------------------------------------------------------------
// This function is used to implement the R function varbvsbinzupdate.
// It is called in R using the .Call interface.
SEXP varbvsbinzupdate_Call (SEXP Xp, SEXP sap, SEXP logoddsp, SEXP dp,
			    SEXP xdxp, SEXP xyp, SEXP dzrp, SEXP alphap,
			    SEXP mup, SEXP Xrp, SEXP ip) {

  // (1) GET INPUTS AND OUTPUTS
  // --------------------------
  double* X       = REAL(Xp);        // Input X.
  double  sa      = *REAL(sap);      // Input sa.
  double* logodds = REAL(logoddsp);  // Input logodds.
  double* d       = REAL(dp);        // Input d.
  double* xdx     = REAL(xdxp);      // Input xdx.
  double* xy      = REAL(xyp);       // Input xy.
  double* dzr     = REAL(dzrp);      // Input dzr.
  double* alpha   = REAL(alphap);    // Input and output alpha.
  double* mu      = REAL(mup);       // Input and output mu.
  double* Xr      = REAL(Xrp);       // Input and output Xr.
  int*    i       = INTEGER(ip);     // Input i.

  // Get the number of samples (n), the number of covariates (m), and
  // the number of coordinate ascent updates (numiter).
  SEXP     r       = getAttrib(dzrp,R_DimSymbol);
  R_xlen_t m       = INTEGER(r)[1];
  R_xlen_t n       = length(Xrp);
  R_xlen_t numiter = length(ip);

  // These arrays are used to store some intermediate calculations.
  double a[m];
  double b[m];

  // (2) CYCLE THROUGH COORDINATE ASCENT UPDATES
  // -------------------------------------------
  for (R_xlen_t j = 0; j < numiter; j++) {
    R_xlen_t k = (R_xlen_t) i[j];

    // Get the kth column of matrix X.
    const double* x = getConstColumn(X,k,n);

    // Perform the update.
    varbvsbinzupdate(x,xy[k],xdx[k],d,dzr,sa,logodds[k],alpha + k,mu + k,Xr,
		     a,b,n,m);
  }

  return R_NilValue;
}

// ---------------------------------------------------------------------
// This function is used to implement the R function varbvsmixupdate.
// It is called in R using the .Call interface.
SEXP varbvsmixupdate_Call (SEXP Xp, SEXP sigmap, SEXP sap, SEXP wp, SEXP xyp,
			   SEXP dp, SEXP alphap, SEXP mup, SEXP Xrp, SEXP ip, 
			   SEXP epsp) {

  // (1) GET INPUTS AND OUTPUTS
  // --------------------------
  double* X       = REAL(Xp);       // Input X.
  double  sigma   = *REAL(sigmap);  // Input sigma.
  double* sa      = REAL(sap);      // Input sa.
  double* w       = REAL(wp);       // Input w.
  double* xy      = REAL(xyp);      // Input xy.
  double* d       = REAL(dp);       // Input d.
  double* alpha   = REAL(alphap);   // Input and output alpha.
  double* mu      = REAL(mup);      // Input and output mu.
  double* Xr      = REAL(Xrp);      // Input and output Xr.
  int*    I       = INTEGER(ip);    // Input i.
  double  eps     = *REAL(epsp);    // Input eps.

  // Get the number of samples (n), the number of mixture components
  // (k), and the number of coordinate ascent updates (numiter).
  R_xlen_t n       = length(Xrp);
  R_xlen_t k       = length(wp);
  R_xlen_t numiter = length(ip);

  // These arrays are used to store some intermediate calculations.
  double s[k];
  double logw[k];

  // (2) CYCLE THROUGH COORDINATE ASCENT UPDATES
  // -------------------------------------------
  for (R_xlen_t j = 0; j < numiter; j++) {
    R_xlen_t i = (R_xlen_t) I[j];

    // Get the kth column of matrix X.
    const double* x = getConstColumn(X,i,n);

    // Perform the update.
    varbvsmixupdate(x,xy[i],d[i],sigma,sa,w,getColumn(alpha,i,k),
		    getColumn(mu,i,k),Xr,s,logw,n,k,eps);
  }

  return R_NilValue;
}

