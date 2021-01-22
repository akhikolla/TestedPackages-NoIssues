#include "rpf.h"

class ssEAP {
	int lastItem;
	void tpbw1995Prep();

	template <typename T1, typename T2>
	void aggregateSpecific(Eigen::ArrayBase<T1> &inMat, Eigen::ArrayBase<T2> &Eis);

	template <typename T1, typename T2, typename T3>
	void tt2prod(Eigen::ArrayBase<T1> &slCur, Eigen::ArrayBase<T2> &buffer, Eigen::ArrayBase<T3> &perSpecific);

	template <typename T1, typename T2, typename T3>
	void prod2ss(Eigen::ArrayBase<T1> &buffer, Eigen::ArrayBase<T2> &ssCur, Eigen::ArrayBase<T3> &perSpecific);

public:
	ifaGroup grp;
	const int *mask;
	int maxScore;
	std::vector<int> items;

	Eigen::ArrayXXd ttCur;
	Eigen::ArrayXi ttCurMax;
	Eigen::ArrayXXd slCur;
	Eigen::ArrayXd ssProbCur;

	Eigen::ArrayXXd ttPrev;
	Eigen::ArrayXi ttPrevCurMax;
	Eigen::ArrayXXd slPrev;
	Eigen::ArrayXd ssProbPrev;

	ssEAP(bool twotier) : grp(twotier) { grp.quad.setNumThreads(1); }
	void setup(SEXP grp, const int *_mask);
	void setLastItem(int which);
	void tpbw1995Vanilla();
	void tpbw1995TwoTier();
	void tpbw1995();

	// tt2ss == two-tier to sum-score
	template <typename T1, typename T2, typename T3>
	void tt2ss(Eigen::ArrayBase<T1> &curMax, Eigen::ArrayBase<T2> &slCur, Eigen::ArrayBase<T3> &ssProbCur);
};

void ssEAP::setup(SEXP robj, const int *_mask)
{
	lastItem = -1;
	mask = _mask;

	grp.import(robj);
}

void ssEAP::setLastItem(int which)
{
	lastItem = which;
}

void ssEAP::tpbw1995Prep()
{
	maxScore = 0;
	for (int cx = 0; cx < grp.numItems(); cx++) {
		const double *spec = grp.spec[cx];
		int no = spec[RPF_ISpecOutcomes];
		if ((lastItem != -1 && cx == lastItem) || mask[cx]) {
			maxScore += no - 1;
			if (cx != lastItem) items.push_back(cx);
		}
	}

	if (lastItem >= 0) items.push_back(lastItem);
}

void ssEAP::tpbw1995()
{
	tpbw1995Prep();

	if (!grp.quad.hasBifactorStructure) {
		tpbw1995Vanilla();
	} else {
		tpbw1995TwoTier();
	}
}

void ssEAP::tpbw1995Vanilla()
{
	ba81NormalQuad::layer &layer = grp.quad.getLayer();
	slCur.resize(layer.totalQuadPoints, 1+maxScore);
	Eigen::Map<Eigen::ArrayXd> areaCol(layer.priQarea.data(), layer.priQarea.size());
	slCur.colwise() = areaCol;

	Eigen::VectorXi abx(layer.maxDims);
	Eigen::VectorXd where(layer.maxDims);
	int curMax = 0;
	int item0 = items[0];
	{
		const double *spec = grp.spec[item0];
		int id = spec[RPF_ISpecID];
		const int dims = spec[RPF_ISpecDims];
		Eigen::VectorXd ptheta(dims);
		double *iparam = grp.getItemParam(item0);
		int outcomes = grp.itemOutcomes[item0];
		Eigen::VectorXd oprob(outcomes);
		for (int qx=0; qx < layer.totalQuadPoints; ++qx) {
			layer.pointToGlobalAbscissa(qx, abx, where);
			for (int dx=0; dx < dims; dx++) {
				ptheta[dx] = where[std::min(dx, layer.maxDims-1)];
			}
			Glibrpf_model[id].prob(spec, iparam, ptheta.data(), oprob.data());
			for (int ox=0; ox < outcomes; ++ox) {
				slCur(qx, ox) *= oprob[ox];
			}
		}
		curMax += outcomes - 1;
	}

	slPrev.resize(layer.totalQuadPoints, 1+maxScore);

	for (int curItem=1; curItem < int(items.size()); ++curItem) {
		int ix = items[curItem];
		slCur.swap(slPrev);
		const double *spec = grp.spec[ix];
		int id = spec[RPF_ISpecID];
		const int dims = spec[RPF_ISpecDims];
		Eigen::VectorXd ptheta(dims);
		double *iparam = grp.getItemParam(ix);
		int outcomes = grp.itemOutcomes[ix];
		Eigen::VectorXd oprob(outcomes);
		slCur.topLeftCorner(slCur.rows(), curMax + outcomes).setZero();
		for (int qx=0; qx < layer.totalQuadPoints; ++qx) {
			layer.pointToGlobalAbscissa(qx, abx, where);
			for (int dx=0; dx < dims; dx++) {
				ptheta[dx] = where[std::min(dx, layer.maxDims-1)];
			}
			Glibrpf_model[id].prob(spec, iparam, ptheta.data(), oprob.data());
			for (int cx=0; cx <= curMax; cx++) {
				for (int ox=0; ox < outcomes; ox++) {
					slCur(qx, cx + ox) += slPrev(qx, cx) * oprob[ox];
				}
			}
		}
		curMax += outcomes - 1;
		checkUserInterrupt();
	}

	ssProbCur.resize(1+maxScore);
	ssProbCur = slCur.colwise().sum();
	ssProbPrev.resize(1+maxScore);     // make smaller TODO
	ssProbPrev = slPrev.colwise().sum();
}

template <typename T1, typename T2>
void ssEAP::aggregateSpecific(Eigen::ArrayBase<T1> &inMat, Eigen::ArrayBase<T2> &Eis)
{
	ba81NormalQuad &quad = grp.quad;
	ba81NormalQuad::layer &layer = grp.quad.getLayer();

	Eis.setZero();
	for (int qx=0, qloc = 0, eisloc = 0; qx < layer.totalPrimaryPoints; qx++) {
		for (int sx=0; sx < quad.gridSize; sx++) {
			for (int sgroup=0; sgroup < layer.numSpecific; ++sgroup) {
				Eis.row(eisloc + sgroup) += inMat.row(qloc);
				++qloc;
			}
		}
		eisloc += layer.numSpecific;
	}
}

template <typename T1, typename T2, typename T3>
void ssEAP::tt2prod(Eigen::ArrayBase<T1> &tt, Eigen::ArrayBase<T2> &buffer, Eigen::ArrayBase<T3> &perSpecific)
{
	ba81NormalQuad::layer &layer = grp.quad.getLayer();

	int combinations = perSpecific.prod();
	int destRows = tt.rows() / perSpecific.rows();
	buffer.setOnes();

	for (int qx=0; qx < destRows; qx++) {
		for (int cx=0; cx < combinations; ++cx) {
			int chip = cx;
			for (int sgroup=0; sgroup < layer.numSpecific; ++sgroup) {
				int col = chip % perSpecific[sgroup];
				chip /= perSpecific[sgroup];
				buffer(qx, cx) *= tt(qx * perSpecific.rows() + sgroup, col);
			}
		}
	}
}

template <typename T1, typename T2, typename T3>
void ssEAP::prod2ss(Eigen::ArrayBase<T1> &buffer, Eigen::ArrayBase<T2> &ssMat, Eigen::ArrayBase<T3> &perSpecific)
{
	ba81NormalQuad::layer &layer = grp.quad.getLayer();
	int combinations = perSpecific.prod();

	ssMat.setZero();
	for (int cx=0; cx < combinations; ++cx) {
		int chip = cx;
		int ss = 0;
		for (int sgroup=0; sgroup < layer.numSpecific; ++sgroup) {
			ss += chip % perSpecific[sgroup];
			chip /= perSpecific[sgroup];
		}
		ssMat.col(ss) += buffer.col(cx);
	}
}

template <typename T1, typename T2, typename T3>
void ssEAP::tt2ss(Eigen::ArrayBase<T1> &curMax1, Eigen::ArrayBase<T2> &curTbl,
		  Eigen::ArrayBase<T3> &outTbl)
{
	int numScores = 1+maxScore;
	ba81NormalQuad::layer &layer = grp.quad.getLayer();

	Eigen::ArrayXXd Eis(layer.totalPrimaryPoints * layer.numSpecific, numScores);
	aggregateSpecific(curTbl, Eis);

	Eigen::ArrayXi perSpecific = curMax1 + 1;
	int combinations = perSpecific.prod();

	Eigen::ArrayXXd prodEis(layer.totalPrimaryPoints, combinations);
	tt2prod(Eis, prodEis, perSpecific);

	outTbl.derived().resize(layer.totalPrimaryPoints, numScores);
	prod2ss(prodEis, outTbl, perSpecific);

	Eigen::Map<Eigen::ArrayXd> areaCol(layer.priQarea.data(), layer.priQarea.size());
	outTbl.colwise() *= areaCol;
}

void ssEAP::tpbw1995TwoTier()
{
	int numScores = 1+maxScore;

	ba81NormalQuad &quad = grp.quad;
	ba81NormalQuad::layer &layer = grp.quad.getLayer();

	Eigen::ArrayXd speAreaCol(layer.totalQuadPoints * layer.numSpecific);
	for (int qx=0, qloc = 0; qx < layer.totalPrimaryPoints; qx++) {
		for (int sx=0; sx < quad.gridSize * layer.numSpecific; sx++) {
			speAreaCol[qloc] = layer.speQarea[sx];
			++qloc;
		}
	}

	ttCur.resize(layer.totalQuadPoints * layer.numSpecific, numScores);
	ttPrev.resize(layer.totalQuadPoints * layer.numSpecific, numScores);
	ttCur.colwise() = speAreaCol;
	ttPrev.colwise() = speAreaCol;

	ttCurMax.resize(layer.numSpecific);
	ttCurMax.setZero();
	ttPrevCurMax = ttCurMax;

	Eigen::VectorXi abx(layer.maxDims);
	Eigen::VectorXd where(layer.numAbil());
	int item0 = items[0];
	{
		const double *spec = grp.spec[item0];
		int id = spec[RPF_ISpecID];
		const int dims = spec[RPF_ISpecDims];
		Eigen::VectorXd ptheta(dims);
		double *iparam = grp.getItemParam(item0);
		int outcomes = grp.itemOutcomes[item0];
		int Sgroup = layer.Sgroup[item0];
		Eigen::VectorXd oprob(outcomes);
		for (int qx=0; qx < layer.totalQuadPoints; ++qx) {
			layer.pointToGlobalAbscissa(qx, abx, where);
			for (int dx=0; dx < dims; dx++) {
				ptheta[dx] = where[std::min(dx, layer.maxDims-1)];
			}
			Glibrpf_model[id].prob(spec, iparam, ptheta.data(), oprob.data());
			for (int ox=0; ox < outcomes; ++ox) {
				ttCur(qx * layer.numSpecific + Sgroup, ox) *= oprob[ox];
			}
		}
		ttCurMax(Sgroup) += outcomes - 1;
	}

	for (int curItem=1; curItem < int(items.size()); ++curItem) {
		int ix = items[curItem];
		ttPrev = ttCur; // can't swap because we only update the item's Sgroup
		ttPrevCurMax = ttCurMax;
		const double *spec = grp.spec[ix];
		int id = spec[RPF_ISpecID];
		const int dims = spec[RPF_ISpecDims];
		Eigen::VectorXd ptheta(dims);
		double *iparam = grp.getItemParam(ix);
		int outcomes = grp.itemOutcomes[ix];
		int Sgroup = layer.Sgroup[ix];
		Eigen::VectorXd oprob(outcomes);
		for (int qx=0; qx < layer.totalQuadPoints; ++qx) {
			int row = qx * layer.numSpecific + Sgroup;
			ttCur.row(row).setZero();
			layer.pointToGlobalAbscissa(qx, abx, where);
			for (int dx=0; dx < dims; dx++) {
				ptheta[dx] = where[std::min(dx, layer.maxDims-1)];
			}
			Glibrpf_model[id].prob(spec, iparam, ptheta.data(), oprob.data());
			for (int cx=0; cx <= ttCurMax(Sgroup); cx++) {
				for (int ox=0; ox < outcomes; ox++) {
					ttCur(row, cx + ox) += ttPrev(row, cx) * oprob[ox];
				}
			}
		}
		ttCurMax(Sgroup) += outcomes - 1;
		checkUserInterrupt();
	}

	tt2ss(ttCurMax, ttCur, slCur);
	ssProbCur.resize(numScores);
	ssProbCur = slCur.colwise().sum();

	tt2ss(ttPrevCurMax, ttPrev, slPrev);
	ssProbPrev.resize(numScores);
	ssProbPrev = slPrev.colwise().sum();
}

template <typename T1, typename T2>
void otMix(ssEAP &myeap, int Sgroup, int ox, Eigen::ArrayBase<T1> &iProb, Eigen::ArrayBase<T2> &numer)
{
	ba81NormalQuad &quad = myeap.grp.quad;
	ba81NormalQuad::layer &layer = quad.getLayer();

	if (layer.numSpecific == 0) {
		numer = (myeap.slPrev.colwise() * iProb.col(ox)).colwise().sum();
	} else {
		Eigen::ArrayXXd ttPrev = myeap.ttPrev;
		for (int qx = 0; qx < layer.totalQuadPoints; ++qx) {
			ttPrev.row(qx*layer.numSpecific + Sgroup) *= iProb(qx, ox);
		}
		Eigen::ArrayXXd ssPrev;
		myeap.tt2ss(myeap.ttPrevCurMax, ttPrev, ssPrev);
		numer = ssPrev.colwise().sum();
	}
}

// [[Rcpp::export]]
SEXP ot2000(SEXP robj, int iPlusOne, bool alter, const LogicalVector &mask, bool twoTier)
{
	int interest = iPlusOne - 1;

	ssEAP myeap(twoTier);
	myeap.setup(robj, mask.begin());
	myeap.setLastItem(interest);
	myeap.tpbw1995();

	int outcomes = myeap.grp.itemOutcomes[interest];

	ifaGroup &grp = myeap.grp;
	ba81NormalQuad::layer &layer = grp.quad.getLayer();
	Eigen::ArrayXXd iProb(layer.totalQuadPoints, outcomes);
	Eigen::VectorXi abx(layer.maxDims);
	Eigen::VectorXd where(layer.numAbil());

	{
		const double *spec = grp.spec[interest];
		int id = spec[RPF_ISpecID];
		const int dims = spec[RPF_ISpecDims];
		Eigen::VectorXd ptheta(dims);
		double *iparam = grp.getItemParam(interest);
		Eigen::ArrayXd oprob(outcomes);
		for (int qx=0; qx < layer.totalQuadPoints; ++qx) {
			layer.pointToGlobalAbscissa(qx, abx, where);
			for (int dx=0; dx < dims; dx++) {
				ptheta[dx] = where[std::min(dx, layer.maxDims-1)];
			}
			Glibrpf_model[id].prob(spec, iparam, ptheta.data(), oprob.data());
			for (int ox=0; ox < outcomes; ox++) {
				iProb(qx, ox) = oprob[ox];
			}
		}
	}

	int Sgroup = 0;
	if (layer.numSpecific) Sgroup = layer.Sgroup[interest];

	if (alter) {
		// as documented in various publications
		Eigen::ArrayXd &ssProb = myeap.ssProbCur;

		NumericMatrix Rexpected(1+myeap.maxScore, outcomes);
		Eigen::Map< Eigen::ArrayXXd > out(Rexpected.begin(), 1+myeap.maxScore, outcomes);
		out.setZero();

		for (int ox=0; ox < outcomes; ++ox) {
			Eigen::ArrayXd numer;
			otMix(myeap, Sgroup, ox, iProb, numer);
			for (int startScore=0; startScore+outcomes <= 1+myeap.maxScore; ++startScore) {
				out(startScore+ox, ox) = numer(startScore) / ssProb(startScore + ox);
			}
		}

		return Rexpected;
	} else {
		// slightly more powerful for small number of items (IRTPRO/flexMIRT)
		Eigen::ArrayXd &ssProb = myeap.ssProbPrev;

		int prevMaxScore = myeap.maxScore - (outcomes-1);
		NumericMatrix Rexpected(prevMaxScore + 1, outcomes);
		Eigen::Map< Eigen::ArrayXXd > out(REAL(Rexpected), prevMaxScore + 1, outcomes);

		for (int ox=0; ox < outcomes; ++ox) {
			Eigen::ArrayXd numer;
			otMix(myeap, Sgroup, ox, iProb, numer);
			for (int startScore=0; startScore <= prevMaxScore; ++startScore) {
				out(startScore, ox) = numer(startScore) / ssProb(startScore);
			}
		}

		return Rexpected;
	}
}

// [[Rcpp::export]]
NumericMatrix ssEAP_wrapper(SEXP robj, SEXP Rmask, bool twoTier)
{
	ssEAP myeap(twoTier);
	myeap.setup(robj, LOGICAL(Rmask));
	myeap.tpbw1995();

	ba81NormalQuad &quad = myeap.grp.quad;
	ba81NormalQuad::layer &layer = quad.getLayer();
	ifaGroup &grp = myeap.grp;

	Eigen::ArrayXXd &slCur = myeap.slCur;

	int curMax = myeap.maxScore;
	int outRows = 1 + curMax;
	int outCols = 1 + 2 * layer.primaryDims + triangleLoc1(layer.primaryDims);

	List dimnames(2);
	CharacterVector names(outCols);
	names[0] = "p";
	for (int ax=0; ax < layer.primaryDims; ++ax) {
		const int SMALLBUF = 20;
		char buf[SMALLBUF];
		names[1+ax] = grp.factorNames[ax].c_str();
		snprintf(buf, SMALLBUF, "se%d", 1+ax);
		names[1+layer.primaryDims+ax] = buf;
	}
	for (int cx=0; cx < triangleLoc1(layer.primaryDims); ++cx) {
		const int SMALLBUF = 20;
		char buf[SMALLBUF];
		snprintf(buf, SMALLBUF, "cov%d", 1+cx);
		names[1+2*layer.primaryDims+cx] = buf;
	}
	dimnames[1] = names;

	Eigen::ArrayXd &ssProb = myeap.ssProbCur;

	ba81NormalQuad pquad;  //primary only
	pquad.setNumThreads(1);
	Eigen::Map< Eigen::VectorXd > meanVec(grp.mean, layer.primaryDims);
	Eigen::Map<Eigen::MatrixXd> fullCov(grp.cov, grp.quad.abilities(), grp.quad.abilities());
	Eigen::MatrixXd priCov = fullCov.block(0, 0, layer.primaryDims, layer.primaryDims);
	Eigen::Map< Eigen::ArrayXXd > Eparam(grp.param, grp.paramRows, grp.numItems());
	pquad.setStructure(grp.qwidth, grp.qpoints, Eparam, meanVec, priCov, false);
	pquad.refresh(meanVec, priCov);

	NumericMatrix Rout(outRows, outCols);
	Rout.attr("dimnames") = dimnames;
	double *out = Rout.begin();
	memcpy(out, ssProb.data(), sizeof(double) * outRows);
	for (int cx=0; cx <= curMax; cx++) {
		Eigen::ArrayXd pad(layer.primaryDims + triangleLoc1(layer.primaryDims));
		Eigen::Map< Eigen::ArrayXd > slCurCol(&slCur.coeffRef(0, cx), slCur.rows());
		pquad.EAP(slCurCol, ssProb[cx], pad);
		for (int sx=0; sx < layer.primaryDims; ++sx) {
			out[(1+sx) * outRows + cx] = pad[sx];
			out[(1+layer.primaryDims+sx) * outRows + cx] = sqrt(pad[layer.primaryDims + triangleLoc0(sx)]);
		}
		for (int sx=0; sx < triangleLoc1(layer.primaryDims); ++sx) {
			out[(1+2*layer.primaryDims+sx) * outRows + cx] = pad[layer.primaryDims + sx];
		}
		checkUserInterrupt();
	}

	return Rout;
}

// [[Rcpp::export]]
NumericMatrix pairwiseExpected_cpp(SEXP robj, IntegerVector items, bool twoTier)
{
	if (items.size() != 2) stop("A pair of items must be specified");

	ifaGroup grp(twoTier);
	grp.quad.setNumThreads(1);
	grp.import(robj);

	ba81NormalQuad &quad = grp.quad;
	ba81NormalQuad::layer &layer = quad.getLayer();

	int i1 = items[0];
	int i2 = items[1];
	if (i1 < 0 || i1 >= (int) grp.spec.size()) stop("Item %d out of range", i1);
	if (i2 < 0 || i2 >= (int) grp.spec.size()) stop("Item %d out of range", i2);
	if (i1 == i2) Rf_warning("Request to create bivariate distribution of %d with itself", i1);

	double *i1par = &grp.param[i1 * grp.paramRows];
	double *i2par = &grp.param[i2 * grp.paramRows];

	int specific1 = -1;
	int specific2 = -1;
	if (layer.numSpecific) {
		int priDims = layer.maxDims-1;
		for (int ax=priDims; ax < quad.abilities(); ax++) {
			if (i1par[ax] != 0) {
				specific1 = ax - priDims;
			}
			if (i2par[ax] != 0) {
				specific2 = ax - priDims;
			}
		}
	}

	const double *spec1 = grp.spec[i1];
	int id1 = spec1[RPF_ISpecID];
	int outcomes1 = spec1[RPF_ISpecOutcomes];

	const double *spec2 = grp.spec[i2];
	int id2 = spec2[RPF_ISpecID];
	int outcomes2 = spec2[RPF_ISpecOutcomes];

	NumericMatrix Rexpected(outcomes1, outcomes2);
	Eigen::Map<Eigen::MatrixXd> out(Rexpected.begin(), outcomes1, outcomes2);
	out.setZero();

	// See Cai & Hansen (2012) Eqn 25, 26

	Eigen::VectorXi abx(layer.maxDims);
	Eigen::VectorXd where(layer.numAbil());
	Eigen::VectorXd o1(outcomes1);
	Eigen::VectorXd o2(outcomes2);

	if (specific1 == -1 && specific2 == -1) {
		int specificIncr = layer.numSpecific? quad.gridSize : 1;
		for (int qx=0; qx < layer.totalPrimaryPoints; ++qx) {
			layer.pointToGlobalAbscissa(qx * specificIncr, abx, where);
			(*Glibrpf_model[id1].prob)(spec1, i1par, where.data(), o1.data());
			(*Glibrpf_model[id2].prob)(spec2, i2par, where.data(), o2.data());
			out += (o1 * o2.transpose()) * layer.priQarea[qx];
		}
	} else if (specific1 == specific2) {
		Eigen::VectorXd ptheta(quad.abilities());
		for (int qloc=0, qx=0; qx < layer.totalPrimaryPoints; ++qx) {
			for (int sx=0; sx < quad.gridSize; ++sx) {
				layer.pointToGlobalAbscissa(qloc, abx, where);
				for (int dx=0; dx < quad.abilities(); dx++) {
					ptheta[dx] = where[std::min(dx, layer.maxDims-1)];
				}
				(*Glibrpf_model[id1].prob)(spec1, i1par, ptheta.data(), o1.data());
				(*Glibrpf_model[id2].prob)(spec2, i2par, ptheta.data(), o2.data());
				double area = layer.priQarea[qx] * layer.speQarea[sx * layer.numSpecific + specific1];
				out += (o1 * o2.transpose()) * area;
				++qloc;
			}
		}
	} else if (specific1 != specific2) {
		Eigen::VectorXd spo1(outcomes1);
		Eigen::VectorXd spo2(outcomes2);
		Eigen::VectorXd ptheta(quad.abilities());
		for (int qloc=0, qx=0; qx < layer.totalPrimaryPoints; ++qx) {
			o1.setZero();
			o2.setZero();
			for (int sx=0; sx < quad.gridSize; ++sx) {
				layer.pointToGlobalAbscissa(qloc, abx, where);
				for (int dx=0; dx < quad.abilities(); dx++) {
					ptheta[dx] = where[std::min(dx, layer.maxDims-1)];
				}
				(*Glibrpf_model[id1].prob)(spec1, i1par, ptheta.data(), spo1.data());
				(*Glibrpf_model[id2].prob)(spec2, i2par, ptheta.data(), spo2.data());
				if (specific1 == -1) {
					if (sx==0) o1 = spo1;
				} else {
					o1 += spo1 * layer.speQarea[sx * layer.numSpecific + specific1];
				}
				if (specific2 == -1) {
					if (sx==0) o2 = spo2;
				} else {
					o2 += spo2 * layer.speQarea[sx * layer.numSpecific + specific2];
				}
				++qloc;
			}
			out += (o1 * o2.transpose()) * layer.priQarea[qx];
		}
	}

	return Rexpected;
}

static const double KANG_CHEN_MIN_EXPECTED = 1.0;  // customizable parameter?

struct ManhattenCollapse {
	Eigen::Map<Eigen::ArrayXXd> obs;
	Eigen::Map<Eigen::ArrayXXd> expected;

	Eigen::DenseIndex smr, smc;
	double bestFit;
	Eigen::DenseIndex bestR, bestC;
	double minExpected;

	ManhattenCollapse(int rows, int cols, double *oMem, double *eMem)
		: obs(oMem, rows, cols), expected(eMem, rows, cols),
		  minExpected(KANG_CHEN_MIN_EXPECTED) {};
	void setMinExpected(double th) { minExpected = th; };
	void probe(Eigen::DenseIndex br, Eigen::DenseIndex bc);
	double findMinCoeff(Eigen::DenseIndex *br, Eigen::DenseIndex *bc);
	int run();
};

void ManhattenCollapse::probe(Eigen::DenseIndex br, Eigen::DenseIndex bc)
{
	bool outside = (br < 0 || bc < 0 ||
			br >= expected.rows() || bc >= expected.cols());
	if (outside) return;
	if (expected(br, bc) < bestFit) {
		bestFit = expected(br, bc);
		bestR = br;
		bestC = bc;
	}
}

double ManhattenCollapse::findMinCoeff(Eigen::DenseIndex *br, Eigen::DenseIndex *bc)
{
	// Eigen's minCoeff is not defined when some cells contain NaN
	double mc = 1e100;
	for (int cx=0; cx < expected.cols(); ++cx) {
		for (int rx=0; rx < expected.rows(); ++rx) {
			if (expected(rx,cx) < mc) {
				mc = expected(rx,cx);
				*br = rx;
				*bc = cx;
			}
		}
	}
	return mc;
}

int ManhattenCollapse::run()
{
	const int maxDist = obs.rows() + obs.cols();
	int collapsed = 0;

	while (findMinCoeff(&smr, &smc) < minExpected) {
		bestFit = 1e100;
		for (int dist=1; dist < maxDist; ++dist) {
			for (int updown=0; updown <= dist; ++updown) {
				int leftright = dist - updown;
				probe(smr + updown, smc + leftright);
				probe(smr + updown, smc - leftright);
				probe(smr - updown, smc + leftright);
				probe(smr - updown, smc - leftright);
			}
			if (bestFit < minExpected) break;
		}

		expected(bestR, bestC) += expected(smr, smc);
		obs(bestR, bestC) += obs(smr, smc);
		expected(smr, smc) = NA_REAL;
		obs(smr, smc) = NA_REAL;
		++collapsed;
	}

	return collapsed;
}

// [[Rcpp::export]]
List collapse(const NumericMatrix &r_observed_orig,
							const NumericMatrix &r_expected_orig, const NumericVector &r_min)
{
  int rows = r_expected_orig.nrow();
	int cols = r_expected_orig.ncol();

  {
    int orows = r_observed_orig.nrow();
		int ocols = r_observed_orig.ncol();
    if (rows != orows || cols != ocols)
	    stop("Observed %dx%d and expected %dx%d matrices must have same dimensions",
		     orows, ocols, rows, cols);
  }

	NumericMatrix r_observed = clone(r_observed_orig);
	NumericMatrix r_expected = clone(r_expected_orig);

  ManhattenCollapse mcollapse(rows, cols, r_observed.begin(), r_expected.begin());
  if (r_min.size()) mcollapse.setMinExpected(r_min[0]);
  int collapsed = mcollapse.run();

	return List::create(_["O"] = r_observed,
											_["E"] = r_expected,
											_["collapsed"] = wrap(collapsed));
}

static int maxObservedSumScore(ifaGroup &grp, const int *itemMask)
{
	int curMax = 0;
	for (int ix=0; ix < int(grp.spec.size()); ++ix) {
		if (!itemMask[ix]) continue;
		const double *spec = grp.spec[ix];
		int no = spec[RPF_ISpecOutcomes];
		curMax += no - 1;
	}
	return curMax;
}

static bool computeObservedSumScore(ifaGroup &grp, const int *itemMask, int row, int *sumOut)
{
	int sum = 0;
	for (int ix=0; ix < int(grp.spec.size()); ++ix) {
		if (!itemMask[ix]) continue;
		const int *resp = grp.dataColumn(ix);
		if (resp[row] == NA_INTEGER) return true;
		sum += resp[row] - 1;
	}
	*sumOut = sum;
	return false;
}

// [[Rcpp::export]]
NumericMatrix fast_tableWithWeights(IntegerVector Ritem1, IntegerVector Ritem2,
													 RObject Rweight)
{
	int rows = Ritem1.size();
	if (rows != Ritem2.size()) stop("Data are of different lengths");

	Eigen::Map<Eigen::ArrayXi> item1(Ritem1.begin(), rows);
	Eigen::Map<Eigen::ArrayXi> item2(Ritem2.begin(), rows);
	double *wvec = 0;
	if (!Rweight.isNULL()) {
		NumericVector weight = as<NumericVector>(Rweight);
		if (weight.size() != rows) stop("Weight vector must be length %d", rows);
		wvec = weight.begin();
	}

	CharacterVector lev1 = Ritem1.attr("levels");
	CharacterVector lev2 = Ritem2.attr("levels");
	int nlev1 = lev1.size();
	int nlev2 = lev2.size();

	NumericMatrix Rdist(nlev1, nlev2);
	Eigen::Map<Eigen::ArrayXXd> result(Rdist.begin(), nlev1, nlev2);
	result.setZero();

	for (int rx=0; rx < rows; ++rx) {
		if (item1[rx] == NA_INTEGER || item2[rx] == NA_INTEGER) continue;
		int i1 = item1[rx] - 1;
		int i2 = item2[rx] - 1;
		double weight = 1;
		if (wvec) weight = wvec[rx];
		result(i1,i2) += weight;
	}

	return Rdist;
}

// [[Rcpp::export]]
List observedSumScore_cpp(SEXP Rgrp, const LogicalVector &Rmask)
{
	ifaGroup grp(false);
	grp.quad.setNumThreads(1);
	grp.import(Rgrp);
	if (grp.getNumUnique() == 0) stop("observedSumScore requires data");
  grp.buildRowMult();

	if (Rmask.size() != int(grp.spec.size())) {
		stop("Mask must be of length %d not %d", int(grp.spec.size()), Rmask.size());
	}
	const int *itemMask = Rmask.begin();

	int numScores = 1+maxObservedSumScore(grp, itemMask);

	NumericVector Rdist(numScores);
	Eigen::Map<Eigen::ArrayXd> distOut(Rdist.begin(), numScores);
	distOut.setZero();

	double rowsIncluded = 0;
	for (int rx=0; rx < grp.getNumUnique(); ++rx) {
		int ss;
		if (computeObservedSumScore(grp, itemMask, rx, &ss)) continue;
		double weight = grp.rowMult(rx);
		distOut[ss] += weight;
		rowsIncluded += weight;
	}

	return List::create(_["dist"] = Rdist,
											_["n"] = wrap(rowsIncluded));
}

// [[Rcpp::export]]
List itemOutcomeBySumScore_cpp(SEXP Rgrp, const LogicalVector &Rmask, int interestPlusOne)
{
	ifaGroup grp(false);
	grp.quad.setNumThreads(1);
	grp.import(Rgrp);
	if (grp.getNumUnique() == 0) stop("itemOutcomeBySumScore requires data");
  grp.buildRowMult();

	if (Rmask.size() != int(grp.spec.size())) {
		stop("Mask must be of length %d not %d", int(grp.spec.size()), Rmask.size());
	}
	const int *itemMask = Rmask.begin();

	int numScores = 1+maxObservedSumScore(grp, itemMask);

	int interest = interestPlusOne - 1;
	if (interest < 0 || interest >= int(grp.spec.size())) {
		stop("Item of interest %d must be between 1 and %d", 1+interest, int(grp.spec.size()));
	}

	const double *spec = grp.spec[interest];
	int outcomes = spec[RPF_ISpecOutcomes];

	NumericMatrix r_ans(numScores, outcomes);
	Eigen::Map<Eigen::ArrayXXd> out(r_ans.begin(), numScores, outcomes);
	out.setZero();

	const int *iresp = grp.dataColumn(interest);

	double rowsIncluded = 0;
	for (int rx=0; rx < grp.getNumUnique(); ++rx) {
		int pick = iresp[rx];
		if (pick == NA_INTEGER) continue;
		int ss;
		if (computeObservedSumScore(grp, itemMask, rx, &ss)) continue;
		double weight = grp.rowMult(rx);
		out(ss, pick-1) += weight;
		rowsIncluded += weight;
	}

	return List::create(_["table"] = r_ans,
											_["n"] = wrap(rowsIncluded));
}

static double table_concordance(const NumericMatrix &mat,
																int rows, int cols, int ii, int jj)
{
  double sum=0;
  for (int hh=ii+1; hh < rows; ++hh) {
    for (int kk=jj+1; kk < cols; ++kk) {
      sum += mat[kk * rows + hh];
    }
  }
  return sum;
}

static double table_discordance(const NumericMatrix &mat,
																int rows, int cols, int ii, int jj)
{
  double sum=0;
  for (int hh=ii+1; hh < rows; ++hh) {
    for (int kk=0; kk < jj; ++kk) {
      sum += mat[kk * rows + hh];
    }
  }
  return sum;
}

// See Agresti (1990, p. 22)
// [[Rcpp::export]]
double gamma_cor(const NumericMatrix &mat)
{
	int rows = mat.nrow();
	int cols = mat.ncol();

  double concord = 0;
  for (int ii=0; ii < rows-1; ++ii) {
    for (int jj=0; jj < cols-1; ++jj) {
      concord += mat[jj * rows + ii] * table_concordance(mat, rows, cols, ii, jj);
    }
  }

  double discord = 0;
  for (int ii=0; ii < rows-1; ++ii) {
    for (int jj=1; jj < cols; ++jj) {
      discord += mat[jj * rows + ii] * table_discordance(mat, rows, cols, ii, jj);
    }
  }

  double gamma = (concord - discord) / (concord + discord);
	return gamma;
}

template <typename T1, typename T2, typename T3>
static inline double crosstabMS(Eigen::ArrayBase<T1> &observed,
				Eigen::ArrayBase<T2> &expected,
				Eigen::ArrayBase<T3> &rowSum)
{
	Eigen::ArrayXXd diff(observed.rows(), observed.cols());
	observed.colwise() /= rowSum;
	diff = observed - expected;
	if (observed.rows() == 1) {
		return ((diff * diff).rowwise().sum() * rowSum).sum();
	} else {
		// not sure if this is correct TODO
		diff.colwise() *= rowSum;
		return (diff * diff).sum();
	}
}

// [[Rcpp::export]]
double crosstabTest_cpp(const NumericMatrix &Robserved,
												const NumericMatrix &Rexpected, int trials)
{
	int rows = Robserved.nrow();
	int cols = Robserved.ncol();
	{
		int erows = Rexpected.nrow();
		int ecols = Rexpected.ncol();
		if (rows != erows || cols != ecols) {
			stop("observed and expected matrices must be the same dimension");
		}
	}

	Eigen::ArrayXXd observed(rows, cols);
	memcpy(observed.data(), Robserved.begin(), sizeof(double) * rows * cols);

	Eigen::ArrayXXd expected(rows, cols);
	memcpy(expected.data(), Rexpected.begin(), sizeof(double) * rows * cols);

	Eigen::ArrayXd rowSum(rows);
	rowSum = observed.rowwise().sum();
	if (((expected.rowwise().sum() - rowSum).abs() > 1e-6).any()) {
		stop("observed and expected row sums must match");
	}

	expected.colwise() /= rowSum;

	Eigen::VectorXi simSize(rows);
	simSize = rowSum.cast<int>();

	if (rows == 1) {
		// Perkins, Tygert, Ward (2011, p. 12)
		for (int sx=0; sx < rows; ++sx) simSize(sx) = std::min(simSize(sx), 185);
	}
	Eigen::ArrayXd simSizeD(rows);
	simSizeD = simSize.cast<double>();

	double refMS = crosstabMS(observed, expected, rowSum);

	Eigen::ArrayXXi Eprob(rows, cols);
	Eprob = (expected * RAND_MAX).cast<int>();

	// SE = 1/(.05 * .95 * sqrt(trials))
	Eigen::ArrayXd mcMS(trials);
	for (int tx=0; tx < trials; ++tx) {
		Eigen::ArrayXXd draw(rows, cols);
		draw.setZero();
		for (int rx=0; rx < rows; ++rx) {
			for (int sx=0; sx < simSize(rx); ++sx) {
				int r1 = RAND_MAX * unif_rand();
				for (int cx=0; cx < cols; ++cx) {
					int threshold = Eprob(rx, cx);
					if (r1 <= threshold) {
						draw(rx, cx) += 1;
						break;
					} else {
						r1 -= threshold;
					}
				}
			}
		}
		mcMS(tx) = crosstabMS(draw, expected, simSizeD);
	}
	Eigen::DenseIndex cnt = (mcMS >= refMS).count();
	return double(cnt) / trials;
}
