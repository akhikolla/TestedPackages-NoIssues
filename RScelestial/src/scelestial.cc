#include <Rcpp.h>
#define logger Rcpp::Rcerr
#include "scelestial.h"
#include "synthesis.h"
using namespace Rcpp;

#ifdef assert
#undef assert
#endif
# define assert(EX) (void)((EX) || (__assert (#EX, __FILE__, __LINE__),0))
void __assert (const char *msg, const char *file, int line) {
    char buffer [100];
    snprintf( buffer, 100, "Assert Failure: %s at %s line #%d", msg, file, line );
    ::Rf_error( buffer );
	throw std::invalid_argument(buffer);
}

typedef std::mt19937 MyRNG;  // the Mersenne Twister with a popular choice of parameters
uint32_t seed_val;           // populate somehow

MyRNG rng;                   // e.g. keep one global instance (per thread)
std::uniform_int_distribution<int> charDistribution('A','Z'); // range [0,10]

char randomCharacter() {
	return charDistribution(rng);
}

vector<string> getColumnNames(DataFrame data) {
	vector<string> ret;
	SEXP colNames = data.names();
	if (Rf_isNull(colNames)) stop("input data does not have column names is null");
	R_xlen_t n = Rf_xlength(colNames);
	for (R_xlen_t i=0; i < n; ++i) {
		const char* name = CHAR(STRING_ELT(colNames, i));
		ret.push_back(name);
		// cerr << "   colName " << i << " " << name << endl;
	}
	return ret;
}


void load(UniverseVertexSet& universeVertexSet, DataFrame data) {
	universeVertexSet.cells.clear();
	// cerr << "  cells cleared" << endl;
	// for (auto i = data.begin(); i != data.end(); i++) {
	for (int i=0; i<data.cols(); i++) {
		// cerr << "  iterate col=" << i << endl;
		IntegerVector col = data[i];
		std::string colClass = col.attr("class");
		if (colClass != std::string("factor")) {
			stop("columns of dataframe should have factor type");
		}
		// cerr << "  iterate col=" << i << " " << col << endl;
		CharacterVector levels = col.attr("levels");
		// cerr << "levels: " << levels << endl;
		universeVertexSet.cells.push_back(Cell());
		for (int j=0; j<col.size(); j++) {
			// std::string s = levels(col(j));
			int e = col(j);
			String s = levels(e-1);
			std::string ss;
			ss = s.get_cstring();
			universeVertexSet.cells.back().append(allelCode(ss));
		}
	}

}

string generateName(int v, set<string>& usedNames) {
	string s;
	for (s = "I" + to_string(v); usedNames.find(s) != usedNames.end() && s.length() < 100; s = "I" + s)
		;
	if (s.length() >= 100) {
		do {
			s = "";
			for (int i=0; i<10; i++) {
				s += randomCharacter();
			}
		} while (usedNames.find(s) != usedNames.end());
	}
	usedNames.insert(s);
	return s;
}

List mapToKeyValueList(const map<string, string>& m, string keyName, string valueName) {
	vector<string> keys;
	vector<string> values;
	for (auto seqIt : m) {
		keys.push_back(seqIt.first);
		values.push_back(seqIt.second);
	}
	return List::create(_[keyName] = keys, _[valueName] = values);
}

List getResultAsGraph(UniverseVertexSet& universeVertexSet, const vector<EdgeWeight>& edges, double cost, 
		const vector<int>& cells, const map<int, Cell>& imputation,
		const vector<string>& inputCellNames) {
	tuple<vector<EdgeWeight>, vector<int>> tt = compressGraph(universeVertexSet, edges, cells);

	if (logLevel > 0)
		logger << "Tree compressed" << endl;

	vector<EdgeWeight> e = get<0>(tt);
	vector<int> v = get<1>(tt);

	set<int> inputCells;
	for (auto c: cells) {
		inputCells.insert(c);
	}

	set<string> usedNames;
	map<int, string> nodeNames;
	for (int j=0; j<(int)v.size(); j++) {
		int i = v[j];
		if (inputCells.find(i) != inputCells.end()) {
			nodeNames[i] = inputCellNames[i];
			usedNames.insert(nodeNames[i]);
		} else
			nodeNames[i] = generateName(i, usedNames);
	}

	if (logLevel > 0)
		logger << "Node names created" << endl;


	// std::cout << v.size() << endl;
	map<string, string> seq, input;
	for (int j=0; j<(int)v.size(); j++) {
		int i = v[j];
		string iName = nodeNames[i];
		// seq.push_back(universeVertexSet.getVertex(i).toString());
		if (inputCells.find(i) != inputCells.end()) {
			input[iName] = universeVertexSet.getVertex(i).toString();
			//following line of code is based on the fact that input vertices
			// are stored with the same order of their input and as the first vertices of universeVertexSet
		}

		seq[iName] = imputation.find(i) != imputation.end() ? imputation.find(i)->second.toString() : universeVertexSet.getVertex(i).toString();

		// std::cout << i << " " << 
		// 	(inputCells.find(i) == inputCells.end() ? 0 : 1) << " " << 
		// 	universeVertexSet.getVertex(i).toString() << " " << 
		// 	(imputation.find(i) != imputation.end() ? imputation.find(i)->second.toString() : "-") << endl;
	}

	// LOGGER( logger << " seq and input created" << endl; )

	List seqX = mapToKeyValueList(seq, "node", "seq");
	List inputX = mapToKeyValueList(input, "node", "seq");

	// LOGGER( logger << " seqX and inputX created" << endl; )

	// std::cout << e.size() << endl;
	vector<string> Esrc, Edst;
	vector<double> Ew;
	for (auto ee: e) {
		vector<int> _ee;
		Esrc.push_back(nodeNames[ee.v]);
		Edst.push_back(nodeNames[ee.u]);
		Ew.push_back(ee.w);
		// E.push_back(_ee);
		// std::cout << ee.v << " " << ee.u << " " << ee.w << endl;
	}

	if (logLevel > 0)
		logger << " E created" << endl;

	CharacterVector EsrcX(Esrc.begin(), Esrc.end());
	// IntegerVector EsrcX(10);
	//  = IntegerVector::create(1, 2, 3, 4, 5);
	CharacterVector EdstX(Edst.begin(), Edst.end());
	NumericVector EwX(Ew.begin(), Ew.end());

	// LOGGER( logger << "data frame created" << endl;)

	List l = List::create(Named("input") = inputX , _["sequence"] = seqX, _["Esrc"] = EsrcX, _["Edst"] = EdstX, _["Ew"] = EwX);
	return l;
}

//' Internal function for running scelestial algorithm.
//'
//' @param data The data
//' @param minK,maxK Minimum and maximum number of vertices to be considered 
//'   for k-restricted Steiner tree.
//' @return The tree as well as missing value imputation
//' @export
// [[Rcpp::export(.scelestial)]]
List _scelestial(DataFrame data, int minK=3, int maxK=4) {
	// LOGGER(
	// 	logger << "data.colnames.size() = " << colnames(data).size() << endl;
	// )
// 
	try {
		init();
		UniverseVertexSet universeVertexSet;
		load(universeVertexSet, data);

		kRestrictionSteinerTreeMin = std::max(3, minK);
		kRestrictionSteinerTreeMax = std::max(kRestrictionSteinerTreeMin, maxK);

		if (data.cols() < kRestrictionSteinerTreeMin) {
			Rcerr << "Error: Number of columns should be at least minK" << endl;
			return R_NilValue;
		}

		assert(universeVertexSet.size() < MAX_SEQUENCE);
		assert(kRestrictionSteinerTreeMax < MAXTREELEAFS);
		if (logLevel > 0)
			Rcerr << "Loaded " << universeVertexSet.cells << endl;

		vector<int> cells;
		for (int i=0; i<universeVertexSet.size(); i++)
			cells.push_back(i);

		tuple<vector<EdgeWeight>, double> t = optimizeTree(universeVertexSet, cells, kRestrictionSteinerTreeMin, kRestrictionSteinerTreeMax );

		if (logLevel > 0)
			logger << "Tree optimized" << " cost=" << get<1>(t) << endl;

		map<int,Cell> imputation = calculateImputation(universeVertexSet, get<0>(t), cells);

		//following will change the graph!
		List l = getResultAsGraph(universeVertexSet, get<0>(t), get<1>(t), cells, imputation, getColumnNames(data));

		
		//TODO: ADD ROOT OPTION
		if (logLevel > 0)
			logger << "Done" << endl;
		return l;
	} catch(const ExitException& e)  {  
        stop(e.what());
    }
}

//' Internal function for generating synthetic single-cell
//' data through simulation of tumor growth and evolution.
//'
//' @param sample Number of samples
//' @param site Number of sites
//' @param evolutionSteps Number of non-root nodes in the evolutionary tree 
//'   to be generated.
//' @param mutationRate The rate of mutation on each evolutionary step in evolutionary tree synthesis.
//' @param advantageIncreaseRatio,advantageDecreaseRatio,advantageKeepRatio A child node
//'   in the evolutionary tree is chosen for increase/decrease/keep its parent advantage with
//'   probabilities proportional to \code{advantage.increase.ratio}/\code{advantage.decrease.ratio}/\code{advantage.keep.ratio}.
//' @param advantageIncreaseStep,advantageDecreaseStep The amount of 
//'   increasing or decreasing the advantage of a cell relative to its parent.
//' @param mvRate Rate of missing value to be added to the resulting sequences.
//' @param fpRate,fnRate Rate of false positive (0 -> 1) and false negative (1 -> 0)
//'   in the sequences.
//' @param seed The seed for randomization.
//' @return The function returns a list. The list consists of 
//'   \itemize{
//'     \item \code{sequence}: A data frame representing
//'       result of sequencing. The data frame has a row for each locus and a column for each sample.
//'     \item \code{true.sequence}: The actual sequence for the sample before adding errors and missing values.
//'     \item \code{true.clone}: A list that stores index of sampled cells for each node in the evolutionary tree.
//'     \item \code{true.tree}: The evolutionary tree that the samples are sampled from. It is a data frame
//'       with \code{src}, \code{dest}, and \code{len} columns representing source, destination and weight of edges of the tree,
//'       respectively.
//'  }
//' 
//' @export
// [[Rcpp::export(.synthesis)]]
List _synthesis(int sample, int site, int evolutionSteps, 
	double mutationRate = 0.01, 
	double advantageIncreaseRatio = 1, double advantageDecreaseRatio = 10, double advantageKeepRatio = 100, 
	double advantageIncreaseStep = 0.01, double advantageDecreaseStep = 0.01, 
	double mvRate = 0.5, 
	double fpRate = 0.2, double fnRate = 0.1,
	int seed = -1) {

	if (logLevel > 0)
		Rcerr << "synthesis C++ start..." << endl;

	if (seed == -1) {
		seed = time(0);
	}
	int seedError = seed;
	synth::generator.seed(seedError);
	// srand(seed);
	//Rcerr << "synthesis C++ seed setted" << endl;

	int sampleCount = sample;
	synth::locusCount = site;
	synth::step = evolutionSteps;

	synth::advIncStep = advantageIncreaseStep;
	synth::advDecStep = advantageDecreaseStep;

	synth::stepMutationRate = mutationRate;
	synth::missingValueRate = mvRate;
	synth::zeroToOneRate = fpRate;
	synth::oneToZeroRate = fnRate;

	double probSum = advantageIncreaseRatio + advantageDecreaseRatio + advantageKeepRatio;
	synth::incAdvProb = advantageIncreaseRatio / probSum;
	synth::decAdvProb = advantageDecreaseRatio / probSum;
	synth::keepAdvProb = advantageKeepRatio / probSum;

	if (logLevel > 0)
		Rcerr << "synthesis C++ before simulate" << endl;
	synth::simulate();
	if (logLevel > 0)
		Rcerr << "synthesis C++ after simulate" << endl;
	synth::Output o = synth::sample(sampleCount);

	if (logLevel > 0)
		Rcerr << "synthesis C++ sampled" << endl;
	//Named("parameters") =  ,
	//Named("input") =  ,
	// Clone
	// Tree
	// Seq
	// TrueSeq

	CharacterVector locusNames;
	for (int j=0; j<synth::locusCount; j++) {
		locusNames.push_back(string("L") + to_string(j+1));
	}

	if (logLevel > 0)
		Rcerr << "synthesis C++ locus names created" << endl;

	List sequences;
	for (int i=0; i<o.sampleCount; i++) {
		NumericVector nv;
		for (int j=0; j<synth::locusCount; j++) {
			nv.push_back(o.output[j][i]);
		}
		sequences.push_back(nv, string("C") + to_string(i+1));
	}
	// rownames(sequences) = locusNames;
	
	if (logLevel > 0)
		Rcerr << "synthesis C++ sequences created" << endl;

	List trueSequences;
	for (int i=0; i<o.sampleCount; i++) {
		NumericVector nv;
		for (int j=0; j<synth::locusCount; j++) {
			nv.push_back(synth::sequence[o.sampleClone[i]][j]);
		}
		trueSequences.push_back(nv, string("C") + to_string(i+1));
	}
	// rownames(trueSequences) = locusNames;

	if (logLevel > 0)
		Rcerr << "synthesis C++ true sequences created" << endl;

	List clones;
	for (int i=0; i<synth::n; i++) {
		vector<string> oneClone;
		for (auto j: o.cloneSamples[i]) {
			oneClone.push_back(string("C") + to_string(j+1));
		}
		if (oneClone.size() > 0)
			clones.push_back(oneClone, "N" + to_string(i+1));
		if (logLevel > 0)
			Rcerr << " oneClone " << oneClone << endl;
	}
	
	if (logLevel > 0)
		Rcerr << "synthesis C++ clones create" << endl;

	List tree;
	{
		// NumericVector treeNodes;
		CharacterVector treeEdgeSrc, treeEdgeDst;
		NumericVector treeEdgeWeight;
		// for (int i=0; i<synth::n; i++) {
		// 	if (synth::parent[i] == -1 || o.cloneSamples[i].size() > 0)
		// 		treeNodes.push_back(i+1);
		// }

		for (int i=0; i<synth::n; i++) {
			if (synth::parent[i] != -1 && o.cloneSamples[i].size() > 0) {
				treeEdgeDst.push_back("N" + to_string(i+1));
				treeEdgeSrc.push_back("N" + to_string(o.parentCompressed[i]+1));
				treeEdgeWeight.push_back(o.parentCompressedDistance[i]);
			}
		}
		tree = List::create(_["s"] = treeEdgeSrc, _["d"] =  treeEdgeDst, _["w"] = treeEdgeWeight);
	}

	if (logLevel > 0)
		Rcerr << "synthesis C++ tree create" << endl;


	List l = List::create(_["sequence"] = sequences, 
		_["true.sequence"] = trueSequences, 
		_["true.clone"] = clones,
		_["true.tree"] = tree,
		_["locus.names"] = locusNames
		);
	
	if (logLevel > 0)
		Rcerr << "synthesis C++ result create" << endl;

	return l;
}

/*** R
*/
