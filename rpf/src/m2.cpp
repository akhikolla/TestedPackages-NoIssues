#include "rpf.h"

// M2 is not implemented yet

struct ch2012 {
	ifaGroup grp;
	bool pearson;
	double stat;
	double weightSum;
	std::vector<bool> rowMask;

	ch2012(bool twotier, SEXP Rgrp);
	void run(const char *method);
	void accumulate(double observed, double expected);
};

ch2012::ch2012(bool twotier, SEXP Rgrp)
	: grp(twotier)
{
	grp.quad.setNumThreads(1);
	grp.import(Rgrp);
	rowMask.reserve(grp.getNumUnique());
	for (int rx=0; rx < grp.getNumUnique(); ++rx) {
		bool missing = false;
		for (int cx=0; cx < (int) grp.dataColumns.size(); ++cx) {
			if (grp.dataColumns[cx][rx] == NA_INTEGER) {
				missing = true;
				break;
			}
		}
		rowMask.push_back(!missing);
	}
}

void ch2012::accumulate(double observed, double expected)
{
	if (pearson) {
		double diff = observed-expected;
		stat += (diff*diff) / expected;
	} else {
		stat += 2 * observed * (log(observed) - log(expected));
	}
	checkUserInterrupt();  // could loop for a long time
}

void ch2012::run(const char *method)
{
	/*
	std::vector<int> &itemOutcomes = grp.itemOutcomes;

	int numFirstOrder = 0;
	for (int ix=0; ix < grp.numItems(); ++ix) {
		numFirstOrder += itemOutcomes[ix] - 1;
	}

	int numSecondOrder = 0;
	for (int i1=1; i1 < grp.numItems(); ++i1) {
		for (int i2=0; i2 < i1; ++i2) {
			numSecondOrder += (itemOutcomes[i1] - 1) * (itemOutcomes[i2] - 1);
		}
	}
	*/

	if (strEQ(method, "pearson")) {
		pearson = true;
	} else if (strEQ(method, "lr")) {
		pearson = false;
	} else {
		stop("Unknown method '%s'", method);
	}

	// Data need to be compressed TODO
	//if (!grp.rowWeight) stop("weightColumn required");  ?? TODO
	ba81NormalQuad &quad = grp.quad;

	weightSum = 0;
	for (int rx=0; rx < grp.getNumUnique(); ++rx) {
		if (!rowMask[rx]) continue;
		weightSum += grp.getRowWeight(rx);
	}

	stat = 0;
	quad.cacheOutcomeProb(grp.param, false);
	quad.allocBuffers();
	for (int px=0; px < grp.getNumUnique(); ++px) {
		if (!rowMask[px]) continue;
		double patternLik1 = quad.computePatternLik(0, px);
		accumulate(grp.getRowWeight(px), patternLik1 * weightSum);
	}
}

// [[Rcpp::export]]
List CaiHansen2012_cpp(SEXP Rgrp, const CharacterVector &Rmethod, bool twotier)
{
	ch2012 engine(twotier, Rgrp);
	engine.run(Rmethod[0]);
	
	//obMargin1(col);
	//obMargin2(col1, col2);

	return List::create(_["stat"] = wrap(engine.stat),
											_["n"] = wrap(engine.weightSum));
}
