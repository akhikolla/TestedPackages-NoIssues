#include "rpf.h"

struct eap {
	std::vector<double *> scoresOut;
	int numLatents;
	Eigen::ArrayXXd thrScore;

	void setup(ba81NormalQuad &quad, int numThr)
	{
		numLatents = quad.abilities() + triangleLoc1(quad.abilities());
		thrScore.resize(numLatents, numThr);
	}
};

template <typename T>
struct BA81LatentScores {
	void normalizeWeights(class ifaGroup *state, T extraData, int px, double weight, int thrId);
	void end(class ifaGroup *state, T extraData);
	bool hasEnd() { return true; }
	bool wantSummary() { return true; }
};

template <typename T>
void BA81LatentScores<T>::normalizeWeights(class ifaGroup *state, T extraData,
					   int px, double patternLik1, int thrId)
{
	ba81NormalQuad &quad = state->quad;
	ba81NormalQuad::layer &layer = quad.getLayer();
	const int maxAbilities = quad.abilities();

	Eigen::Map< Eigen::ArrayXd > wvec(&layer.Qweight.coeffRef(0, thrId), layer.Qweight.rows());
	double *scorePad = &extraData.thrScore.coeffRef(0, thrId);
	Eigen::Map< Eigen::ArrayXd > EscorePad(scorePad, extraData.numLatents);

	quad.EAP(wvec, patternLik1, EscorePad);

	std::vector<double*> &out = extraData.scoresOut;

	for (int ax=0; ax < maxAbilities; ++ax) {
		out[ax][px] = scorePad[ax];
	}
	for (int ax=0; ax < maxAbilities; ++ax) {
		out[maxAbilities + ax][px] = sqrt(scorePad[maxAbilities + triangleLoc0(ax)]);
	}
	for (int ax=0; ax < triangleLoc1(maxAbilities); ++ax) {
		out[2*maxAbilities + ax][px] = scorePad[maxAbilities + ax];
	}
}

template <typename T>
void BA81LatentScores<T>::end(class ifaGroup *state, T extraData)
{
	std::vector<int> &rowMap = state->rowMap;
	const int numUnique = (int) rowMap.size();
	std::vector<double*> &out = extraData.scoresOut;

	for (int px=0; px < numUnique; px++) {
		if (state->patternLik[px]) continue;
		for (int ax=0; ax < int(out.size()); ++ax) {
			out[ax][px] = NA_REAL;
		}
	}
}

// [[Rcpp::export]]
DataFrame eap_wrapper(SEXP Rgrp)
{
	eap eapContext;

	ifaGroup grp(true);
	grp.quad.setNumThreads(GlobalNumberOfCores);
	grp.import(Rgrp);
	grp.buildRowSkip();
	if (grp.getNumUnique() == 0) {
		stop("EAP requested but there are no data rows");
	}
	ba81NormalQuad &quad = grp.quad;
	quad.cacheOutcomeProb(grp.param, false);

	// TODO Wainer & Thissen. (1987). Estimating ability with the wrong
	// model. Journal of Educational Statistics, 12, 339-368.

	/*
	int numQpoints = state->targetQpoints * 2;  // make configurable TODO

	if (numQpoints < 1 + 2.0 * sqrt(state->itemSpec->cols)) {
		// Thissen & Orlando (2001, p. 136)
		Rf_warning("EAP requires at least 2*sqrt(items) quadrature points");
	}

	ba81SetupQuadrature(oo, numQpoints, 0);
	ba81Estep1(oo);
	*/

	int maxAbilities = quad.abilities();
	if (maxAbilities == 0) stop("At least 1 factor is required");
	int rows = grp.getNumUnique();  // allow indexvector for compressed tables TODO
	int cols = 2 * maxAbilities + triangleLoc1(maxAbilities);

	List Rscores(cols);
	for (int cx=0; cx < cols; ++cx) {
		NumericVector vec(rows);
		Rscores[cx] = vec;
		eapContext.scoresOut.push_back(vec.begin());
	}

	const int SMALLBUF = 20;
	char buf[SMALLBUF];
	CharacterVector names(cols);
	for (int nx=0; nx < maxAbilities; ++nx) {
		names[nx] = grp.factorNames[nx].c_str();
		snprintf(buf, SMALLBUF, "se%d", nx+1);
		names[maxAbilities + nx] = buf;
	}
	for (int nx=0; nx < triangleLoc1(maxAbilities); ++nx) {
		snprintf(buf, SMALLBUF, "cov%d", nx+1);
		names[maxAbilities*2 + nx] = buf;
	}
	Rscores.attr("names") = names;

	if (grp.dataRowNames) {
		Rscores.attr("row.names") = grp.dataRowNames;
	}

	eapContext.setup(quad, GlobalNumberOfCores);
	BA81Engine<eap&, BA81LatentScores, BA81OmitEstep> engine;
	engine.ba81Estep1(&grp, eapContext);

	return DataFrame(Rscores);
}
