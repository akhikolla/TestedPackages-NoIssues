#include "rpf.h"

// the code in this file is copied from OpenMx

static int elementEqualsDataframe(SEXP column, int offset1, int offset2) {
	switch (TYPEOF(column)) {
	case REALSXP:
		if(ISNA(REAL(column)[offset1])) return ISNA(REAL(column)[offset2]);
		if(ISNA(REAL(column)[offset2])) return ISNA(REAL(column)[offset1]);
		return(REAL(column)[offset1] == REAL(column)[offset2]);
	case LGLSXP:
	case INTSXP:
		return(INTEGER(column)[offset1] == INTEGER(column)[offset2]);
	}
	stop("elementEqualsDataframe: unknown type %d", TYPEOF(column)); //nocov
}

static int testRowDataframe(SEXP data, int numrow, int numcol, int i, int *row, int base) {
	SEXP column;
	int j, equal = TRUE;

	if (i == numrow) {
		equal = FALSE;
	} else {
		for(j = 0; j < numcol && equal; j++) {
			column = VECTOR_ELT(data, j);
			equal = elementEqualsDataframe(column, base, i);
		}
	}

	if (!equal) {
		int gap = i - base;
		for(j = 0; j < gap; j++) {
			row[base + j] = gap - j;
		}
		base = i;
	}
	return(base);
}

static SEXP findIdenticalDataFrame(SEXP data, SEXP missing, SEXP defvars,
			    SEXP skipMissingExp, SEXP skipDefvarsExp) {

	SEXP retval, identicalRows, identicalMissing, identicalDefvars;
	int i, numrow, numcol;
	int *irows, *imissing, *idefvars;
	int baserows;
	int skipMissing, skipDefvars;

	skipMissing = LOGICAL(skipMissingExp)[0];
	skipDefvars = LOGICAL(skipDefvarsExp)[0];
	numrow = Rf_length(VECTOR_ELT(data, 0));
	numcol = Rf_length(data);
	Rf_protect(retval = Rf_allocVector(VECSXP, 3));
	Rf_protect(identicalRows = Rf_allocVector(INTSXP, numrow));
	Rf_protect(identicalMissing = Rf_allocVector(INTSXP, numrow));
	Rf_protect(identicalDefvars = Rf_allocVector(INTSXP, numrow));
	irows = INTEGER(identicalRows);
	imissing = INTEGER(identicalMissing);
	idefvars = INTEGER(identicalDefvars);
	if (skipMissing) {
		for(i = 0; i < numrow; i++) {
			imissing[i] = numrow - i;
		}
	}
	if (skipDefvars) {
		for(i = 0; i < numrow; i++) {
			idefvars[i] = numrow - i;
		}
	}
	baserows = 0;
	for(i = 1; i <= numrow; i++) {
		baserows = testRowDataframe(data, numrow, numcol, i, irows, baserows); 
	}
	SET_VECTOR_ELT(retval, 0, identicalRows);
	SET_VECTOR_ELT(retval, 1, identicalMissing);
	SET_VECTOR_ELT(retval, 2, identicalDefvars);
	Rf_unprotect(4); // retval, identicalRows, identicalMissing, identicalDefvars
	return retval;
}

// [[Rcpp::export]]
SEXP findIdenticalRowsData(SEXP data, SEXP missing, SEXP defvars,
			   SEXP skipMissingness, SEXP skipDefvars)
{
	if (Rf_isMatrix(data)) {
		stop("Only data.frame is implemented");
	} else {
		return(findIdenticalDataFrame(data, missing, defvars, skipMissingness, skipDefvars));
	}
}
