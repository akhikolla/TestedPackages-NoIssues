#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <map>
#include <set>
#include <iterator>
#include <algorithm>
//#include <boost/range/adaptor/map.hpp>
//#include <boost/range/algorithm/copy.hpp>
//#include <boost/assign.hpp>
#include <cassert>
#include <cmath>
#include <chrono>

using std::string;
using std::ostream;
using std::vector;
using std::set;
using std::map;
using std::tuple;
using std::ostringstream;
using std::to_string;
using std::istream;
using std::ostream;
using std::cin;
using std::get;
using std::istringstream;
using std::stoi;
using std::stod;
using std::min;

class ExitException: public std::exception {
	public:
	ExitException(string _reason = "") : reason(_reason) {}
	virtual const char* what() const throw() {
		return reason.c_str();
	}
	private:
	string reason;
};

// #define USE_HAMILTONIAN_PATH_LOWER_BOUND_FOR_SMALL_STEINER_TREES_PROUNING
// #define DP_DOUBLE_TO_INT_APPROX // It has some BUGS!! Do not use it!

template<typename T, typename C>
ostream& operator<<(ostream& os, const set<T, C>& c) {
	for (auto const &i : c)
		os << i << " ";
	return os;
}

template<typename T>
ostream& operator<<(ostream& os, const vector<T>& c) {
	for (auto const &i : c)
		os << i << " ";
	return os;
}

template<typename T, typename K>
ostream& operator<<(ostream& os, const map<T, K>& c) {
	for (auto &i : c)
		os << i.first << ":" << i.second << " ";
	return os;
}

template<typename T>
string array2string(const T* a, const T* b) {
	ostringstream os;
	for (const T* i=a; i<b; i++)
		os << *i << " ";
	return os.str();
}

enum debug_option
{
    DEBUG_DISABLE,
    DEBUG_ENABLE,
    DEBUG_ENABLE_V
};
debug_option logLevel = DEBUG_DISABLE;
using std::endl;

const int MAXTREELEAFS = 10;
const int MAXNODE = 3 * MAXTREELEAFS;
const double EPSILON = 1e-7;

char nucleotideAcids[] = {'A', 'T', 'C', 'G'};
const int nucleotideAcidCount = 4;
const int allelCodingSize = 11;
string allelCoding[][2] = {
	{"A","A/A"},		
						{"T","T/T"},
						{"C","C/C"},
						{"G","G/G"},
						{"K","A/C"},
						{"L","A/G"},
						{"M","C/T"},
						{"N","C/G"},
						{"O","T/G"},
						{"P","T/A"},
						{"X","./."}};
map<string, char> allelCodingMap;
map<char,int> allelCoding2Int;
// map<int, char> int2AllelCoding;
char int2AllelCoding[allelCodingSize];
int xCode;
void init() {
	for (int i=0; i<allelCodingSize; i++) {
		allelCodingMap[allelCoding[i][1]] = allelCoding[i][0][0];
		string r = string("") + allelCoding[i][1][2] + allelCoding[i][1][1] + allelCoding[i][1][0];
		allelCodingMap[r] = allelCoding[i][0][0];
	}

	// int j = 0;
	for (int i=0; i<allelCodingSize; i++) {
		allelCoding2Int[allelCoding[i][0][0]] = i;
		int2AllelCoding[i] = allelCoding[i][0][0];
	}

	xCode = allelCoding2Int['X'];

	// logger << "Allel: " << allelCoding2Int << "\n" << int2AllelCoding << endl
	//     << allelCoding2Int << endl
	//     << int2AllelCoding << endl;
}

char allelCode(string allel) {
	if (allelCodingMap.find(allel) == allelCodingMap.end()) {
		if (logLevel > 0)
			logger << ("Invalid Allel " + allel) << endl;
		throw ExitException("Invalid Allel " + allel);
	}
	//logger << "C " << allel << "->" << allelCodingMap[allel];
	return allelCodingMap[allel];
}

long long power(int a, int b, int mod) {
	if (b == 0)
		return 1;
	long long p = power(a, b/2, mod);
	p = (p * p) % mod;
	if (b % 2 == 1)
		p = p * a;
	return p % mod;
}

const double IMPUTATION_COST = 0.50001;
const double IMPUTATION_X_X_FACTOR = 2 * IMPUTATION_COST;
const int MAX_SEQUENCE = 1000;
int kRestrictionSteinerTreeMax = 4, kRestrictionSteinerTreeMin = 3;
// We know that MST/2 < Stiener-tree, thus by setting this to 0.5, this is only a valid bound,
// If we let this value to be more than 0.5, we may miss some subsets, but the algorithm would finish faster. 
double steinerTreeMSTLowerBoundRate = 0.5;

#if !defined(IMPUT_COST_2) && !defined(IMPUT_COST_3)
inline double charDistance(char a, char b) {
	if (a == 'X' && b == 'X')
		return IMPUTATION_X_X_FACTOR * IMPUTATION_COST;
	else
		return (a=='X' || b=='X') ? IMPUTATION_COST : a == b ? 0 : 1;
}
#endif
#ifdef IMPUT_COST_2
inline double charDistance(char a, char b) {
	if (a == 'X' && b == 'X')
		return 0;
	else
		return (a=='X' || b=='X') ? IMPUTATION_COST : a == b ? 0 : 1;
}
#endif
#ifdef IMPUT_COST_3
inline double charDistance(char a, char b) {
	if (a == 'X' && b == 'X')
		return IMPUTATION_COST + 0.00001;
	else
		return (a=='X' || b=='X') ? IMPUTATION_COST : a == b ? 0 : 1;
}
#endif

struct Cell {
	vector<char> s;
	// vector<int> sIndex;

	Cell() {
		// hashValue = 0;
	}

	int size() const {
		return (int) s.size();
	}

	void load(string str) {
		// hashValue = 0;
		s.clear();
		for (auto& c: str)
			append(c);
	}

	void append(char a) {
		// hashValue = (hashValue + allelCoding2Int[a] * power(hashP, s.size(), hashM)) % hashM;
		// if (hashValue < 0) {
		//     logger << "OH!!!" << " " << hashValue << " " << power(hashP, s.size(), hashM) << " " 
		//         << allelCoding2Int[a] * power(hashP, s.size(), hashM) << endl;
		// }
		// assert(hashValue >= 0);
		s.push_back(a);
		// sIndex.push_back(allelCoding2Int[a]);
	}

	void set(char a, int pos) {
		// hashValue = (hashValue + 
		//             power(hashP, pos, hashM) * 
		//                 (allelCoding2Int[a] - allelCoding2Int[s[pos]]) 
		//             + 100 * hashM) % hashM;
		s[pos] = a;
		// sIndex[pos] = allelCoding2Int[a];
	}

	char get(int pos) const {
		return s[pos];
	}

	double distance(const Cell& c) const {
		double cnt = 0;
		for (int i=0; i<(int)s.size(); i++) {
			cnt += charDistance(s[i], c.s[i]);
			// if (s[i] == 'X' && c.s[i] == 'X')
			//     cnt += 100 * IMPUTATION_COST; //DEBUG
			// else
			//     cnt += (s[i] == c.s[i]) ? 0 : ((s[i] == 'X' || c.s[i] == 'X') ? IMPUTATION_COST : s[i] != c.s[i]);
		}
		return cnt;
	}

	string toString() const {
		string r(s.begin(), s.end());
		// return r + ":" + to_string(hashValue);
		return r;
	}

	// static const int hashP;
	// static const int hashM;
	
	// int hashValue;

	// int hash() const {
	//     return hashValue;
	// }

};

// const int Cell::hashP = 1610612741L; 
// const int Cell::hashM = 1000000000 + 9;

ostream& operator<<(ostream& os, const Cell& s) {
	// return os << s.toString() << " " << s.s.size() << endl;
	return os << s.toString();
}

// template<typename T>
// struct DisjointSet {
//     mutable map<T, T> parent;

//     DisjointSet() {}

//     template<typename It>
//     DisjointSet(It begin, It end) {
//         for (It i = begin; i != end; i++) {
//             parent[*i] = *i;
//         }
//     }

//     void add(const T& t) {
//         parent[t] = t;
//     }

//     T parentStar(const T& t) const {
//         if (parent[t] == t)
//             return t;
//         return parent[t] = parentStar(parent[t]);
//     }

//     void join(const T& a, const T& b) {
//         parent[parentStar(a)] = parent[parentStar(b)];
//     }

//     bool isJoint(const T& a, const T& b) const {
//         return parentStar(a) == parentStar(b);
//     }
// };
// template<typename T>
// ostream& operator<<(ostream& os, const DisjointSet<T>& ds) {
//     os << "{";
//     for (auto i: ds.parent) {
//         os << i.first << ":" << i.second << " ";
//     }
//     return os << "}";
// }


struct DisjointSetArray {
	// mutable int parent[MAX_SEQUENCE];
	mutable vector<int> parent;

	DisjointSetArray() : parent(MAX_SEQUENCE) {
		for(int i=0; i<MAX_SEQUENCE; i++)
			parent[i] = -1;
	}

	template<typename It>
	DisjointSetArray(It begin, It end) : parent(MAX_SEQUENCE)  {
		for (It i = begin; i != end; i++) {
			parent[*i] = *i;
		}
	}

	void add(int t) {
		parent[t] = t;
	}

	int parentStar(int t) const {
		if (parent[t] == t)
			return t;
		return parent[t] = parentStar(parent[t]);
	}

	void join(int a, int b) {
		parent[parentStar(a)] = parent[parentStar(b)];
	}

	bool isJoint(int a, int b) const {
		return parentStar(a) == parentStar(b);
	}
};

ostream& operator<<(ostream& os, const DisjointSetArray& ds) {
	os << "{";
	for (int i=0; i<MAX_SEQUENCE; i++) {
		if (ds.parent[i] != -1) {
			os << i << ":" << ds.parent[i] << " ";
		}
	}
	return os << "}";
}


struct GenerateAllTrees;

ostream& operator<<(ostream& os, const GenerateAllTrees& a);

struct GenerateAllTrees {
	//parameters
	int n;
	bool unRooted;

	//inner variables
	vector<int> nodeIndex;
	int nodeCount;
	vector<set<int>> availableLeafs;
	/** children of node v in tree is in childs[v]. If childs[v].size() == 1, v is a leaf and childs[v][0] is the index of the leaf, 
	 * on which this tree is built. 
	 */ 
	vector<vector<int>> childs;

	//output
	vector<string> trees;
	char treeRepresentation[100 * 8 * 3]; // MAXN * 3 * MAXLGN
	int treeRepresentationEnd;
	vector< vector<vector<int>> > outChilds;

	//debug
	int depth;

	GenerateAllTrees(int _n, bool _unRooted = true) : n(_n), unRooted(_unRooted) {
	}

	vector<vector<vector<int>>> run() {
		trees.clear();
		depth = 0;
		treeRepresentationEnd = 0;
		availableLeafs.push_back(set<int>());
		for (int i=0; i<n; i++) {
			availableLeafs.back().insert(i);
		}
		treeRepresentationAppend("(");
		childs.clear();
		childs.push_back(vector<int>());
		nodeIndex.push_back(0);
		nodeCount = 1;
		outChilds.clear();
		rec();
		return outChilds;
	}

	int treeRepresentationAppend(string s) {
		int l = treeRepresentationEnd;
		// strcpy(treeRepresentation + treeRepresentationEnd, s.c_str());
		for (int i=0, j=treeRepresentationEnd; i< (int)s.length(); i++, j++) {
			treeRepresentation[j] = s[i];
		}
		treeRepresentationEnd += (int) s.length();
		return l;
	}

	void treeRepresentationMoveBack(int lt) {
		treeRepresentation[lt] = 0;
		treeRepresentationEnd = lt;
	}

	void forAllSubsets(const set<int>& refSet, set<int>::const_iterator i, set<int>& s1, set<int>& s2) {
		//logger << " 2^S: " << *i << " [" << s1 << "] [" << s2 << "]" << endl;
		if (i == refSet.end()) {
			if (s2.size() > 0 || refSet.size() == 1) {
				//Set is ready.
				availableLeafs.push_back(s2); // replaces the current available leafs for the current node with the remaining nodes, s2. 
				//Node that, current set, refSet, is already removed from the list.
				availableLeafs.push_back(s1); // leafs for the newly generated node.

				childs[nodeIndex.back()].push_back(nodeCount);
				childs.push_back(vector<int>());
				nodeIndex.push_back(nodeCount);
				nodeCount++;
				int lt = treeRepresentationAppend("(");
				rec();
				treeRepresentationMoveBack(lt);
				nodeCount--;
				nodeIndex.pop_back();
				childs.pop_back();
				childs[nodeIndex.back()].pop_back();

				availableLeafs.pop_back();
				availableLeafs.pop_back();
			}
			return;
		}
		set<int>::const_iterator j = i;
		j++;
		s1.insert(*i);
		forAllSubsets(refSet, j, s1, s2);
		s1.erase(s1.find(*i));
		if (i != refSet.begin()) {
			//always put first to the first set
			s2.insert(*i);
			forAllSubsets(refSet, j, s1, s2);
			s2.erase(s2.find(*i));
		}
	}

	void rec() {
		depth++;
		//logger << "Rec: \n" << *this << endl;

		if (availableLeafs.size() == 0) {
			//Tree is ready
			if (!unRooted || childs[0].size() >= 3) {
				// std::cout << treeRepresentation << endl;
				trees.push_back(treeRepresentation);
				outChilds.push_back(childs);
				// logger << "  ";
				// for (auto v: childs) {
				//     logger << v << " /";
				// }
				// logger << endl;
			}
		} else if (availableLeafs.back().size() == 0) {
			//we return to a node, all its available childs are already handled.
			int p = nodeIndex.back();
			nodeIndex.pop_back();
			availableLeafs.pop_back();
			int lt = treeRepresentationAppend(")");
			rec();
			treeRepresentationMoveBack(lt);
			availableLeafs.push_back(set<int>());
			nodeIndex.push_back(p);
		} else if (availableLeafs.back().size() == 1 && childs[nodeIndex.back()].size() == 0) {
			//this is a new child.
			int n = *availableLeafs.back().begin();
			availableLeafs.pop_back();

			int p = nodeIndex.back();
			nodeIndex.pop_back();
			childs[p].push_back(n); // no need to be removed after rec, it will be removed from forAllSubsets imediately after
			int lt = treeRepresentationAppend(to_string(n) + ")");
			rec();
			treeRepresentationMoveBack(lt);
			nodeIndex.push_back(p);

			set<int> r;
			r.insert(n);
			availableLeafs.push_back(r);
		} else {
			set<int> r = availableLeafs.back();
			availableLeafs.pop_back();

			set<int> s1, s2;
			forAllSubsets(r, r.begin(), s1, s2);

			availableLeafs.push_back(r);
		}

		depth--;
	}
};

ostream& operator<<(ostream& os, const GenerateAllTrees& a) {
	char tab[100];
	for (int i=0; i<a.depth * 2; i++)
		tab[i] = ' ';
	tab[a.depth * 2] = 0;

	int j=0;
	for (auto i: a.availableLeafs) {
		os << tab << j << ": " << i << endl;
		j++;
	}

	os << tab << "  T:   " << a.treeRepresentation << endl;
	return os;
}

//TODO: make a tree struct, instead of vector<vector<int>>

struct UniverseVertexSet {

	vector<Cell> cells;

	int length() const {
		if (cells.size() == 0)
			return 0;
		return cells[0].size();
	}

	int size() const {
		return (int) cells.size();
	}

	int add(const Cell& c) {
		cells.push_back(c);
		return (int)(cells.size()) - 1;
	}

	Cell& getVertex(int index) {
		return cells[index];
	}

	const Cell& getVertex(int index) const {
		return cells[index];
	}

	double distance(int v, int u) const {
		return cells[v].distance(cells[u]);
	}

	int precalculatedPairwiseDistanceSize;
	double precalculatedPairwiseDistance[MAX_SEQUENCE][MAX_SEQUENCE];

	void precalculatePairwiseDistances() {
		assert(cells.size() < MAX_SEQUENCE);
		precalculatedPairwiseDistanceSize = (int) cells.size();
		for (int i=0; i<(int)cells.size(); i++) {
			precalculatedPairwiseDistance[i][i] = 0;
			for (int j=i+1; j<(int)cells.size(); j++)
				precalculatedPairwiseDistance[i][j] = precalculatedPairwiseDistance[j][i] = distance(i, j);
		}
	}

	double inputSequencesDistance(int v, int u) const {
		if (v >= precalculatedPairwiseDistanceSize || u >= precalculatedPairwiseDistanceSize) {
			if (logLevel > 0)
				logger << "!! inputSeqDistance requested but not pre calculated " << v << " " << u << " " << precalculatedPairwiseDistanceSize << endl;
			return distance(v, u);
		}
		return precalculatedPairwiseDistance[v][u];
	}

};

struct EdgeWeight;
ostream& operator<<(ostream& os, const EdgeWeight& e);

struct EdgeWeight {
	int v, u;
	double w;

	EdgeWeight(int _v=0, int _u=0, double _w=0) : v(_v), u(_u), w(_w) {}

	bool operator<(const EdgeWeight& b) const {
		const EdgeWeight& a = *this;
		long long w1 = (long long)(a.w/EPSILON),
			w2 = (long long)(b.w/EPSILON);
		if (w1 != w2) {
			return w1 < w2;
		}
		if (a.v != b.v) {
			return a.v < b.v;
		}
		return a.u < b.u;
	}
};

ostream& operator<<(ostream& os, const EdgeWeight& e) {
	return os << "(" << e.v << " " << e.u << " " << e.w << ")";
}

#ifndef DP_DOUBLE_TO_INT_APPROX
#else
const int DOUBLE_COST_TO_COST_APPROX_FACTOR = 1000;
#endif

/**
 */ 
tuple<vector<EdgeWeight>, double> imputeTree(
		const vector<int>& leaves, 
		const vector<vector<int>>& tree, 
		UniverseVertexSet& universeVertexSet, 
		bool onlyCost) {
	int len = universeVertexSet.length();

	double finalCost = 0;
	vector<vector<int>> imputedTree(tree.size(), vector<int>(len, -1));

	for (int l=0; l<len; l++) {
		#ifndef DP_DOUBLE_TO_INT_APPROX
		// HADI APPROX: 
		double cost[MAXNODE][allelCodingSize];
		#else
		long long cost[MAXNODE][allelCodingSize];
		#endif
		int path[MAXNODE][allelCodingSize][MAXNODE];
		for (int v=(int)(tree.size())-1; v>=0; v--) {
			for (int a = 0; a<allelCodingSize; a++) {
				if (tree[v].size() == 1) {
					// a leaf!
					int c = tree[v][0];
					const Cell& leaf = universeVertexSet.getVertex(leaves[c]);

					#ifndef DP_DOUBLE_TO_INT_APPROX
					// HADI APPROX: 
					cost[v][a] = charDistance(leaf.s[l], int2AllelCoding[a]);
					#else
					cost[v][a] = (long long)(DOUBLE_COST_TO_COST_APPROX_FACTOR * charDistance(leaf.s[l], int2AllelCoding[a]));
					#endif
					// (leaf.sIndex[l] == a) ? 0 : ((leaf.s[l] == 'X') ? IMPUTATION_COST : 100000);
				} else {
					#ifndef DP_DOUBLE_TO_INT_APPROX
					// HADI APPROX: 
					double cstSum = 0;
					#else
					long long cstSum = 0;
					#endif
					for (auto c: tree[v]) {
						#ifndef DP_DOUBLE_TO_INT_APPROX
						// HADI APPROX: 
						double cst = cost[c][a];
						#else
						long long cst = cost[c][a];
						#endif
						path[v][a][c] = a;
						for (int aa=0; aa<allelCodingSize; aa++)
							if (aa != xCode && cst > cost[c][aa] + 1) {
								cst = cost[c][aa]+1;
								path[v][a][c] = aa;
							}
						cstSum += cst;
					}
					cost[v][a] = cstSum;
				}
			}
		}

		// logger << "IMP: " << l << endl;
		// for (int v=0; v<tree.size(); v++) {
		//     logger << " " << v << ": ";
		//     for (int a = 0; a<allelCodingSize; a++) {
		//         logger << cost[v][a] << " ";
		//     }
		//     logger << endl;
		// }


		int minA = -1;
		for (int a=0; a<allelCodingSize; a++) {
			if (a != xCode && (minA == -1 || cost[0][a] < cost[0][minA]))
				minA = a;
		}

		// HADI APPROX: 
		#ifndef DP_DOUBLE_TO_INT_APPROX
		finalCost += cost[0][minA];
		#else
		finalCost += cost[0][minA] * 1.0 / DOUBLE_COST_TO_COST_APPROX_FACTOR;
		#endif

		// logger << "  Path: ";

		if (!onlyCost) {
			vector<int> minAs(tree.size(), -1);
			minAs[0] = minA;
			for (int v=0; v<(int)tree.size(); v++) {
				// logger << v << ":" << minAs[v] << "/ ";
				imputedTree[v][l] = minAs[v];
				if (tree[v].size() == 1) {
					// imputed[tree[v][0]][l] = minAs[v];
				} else {
					for (auto c: tree[v]) {
						// int cst = cost[c][minAs[v]],
						//     cstA = minAs[v];
						// for (int aa=0; aa<allelCodingSize; aa++)
						//     if (aa != xCode && cst > cost[c][aa]+1) {
						//         cst = cost[c][aa]+1;
						//         cstA = aa;
						//     }
						// minAs[c] = cstA;
						minAs[c] = path[v][minAs[v]][c];
					}
				}
			}
			// logger << endl;
		}

	}
	// return imputedTree;


	if (!onlyCost) {
		// logger << "Generating graph ..." << endl;
		vector<EdgeWeight> gEdges;
		vector<int> treeNode2UniverseIndex(tree.size(), -1);
		for (int v=(int)(imputedTree.size())-1; v>=0; v--) {
			// logger << "  " << v << ": ";
			Cell c;
			for (auto loc : imputedTree[v]) {
				c.append(int2AllelCoding[loc]);
				// logger << loc << " " << int2AllelCoding[loc] << " ";
			}
			// logger << endl;
			int h = universeVertexSet.add(c);
			treeNode2UniverseIndex[v] = h;

			if (tree[v].size() == 1) {
				int leafIndex = leaves[tree[v][0]];
				gEdges.push_back(EdgeWeight(
						leafIndex, h, 
						universeVertexSet.distance(leafIndex, h)));
					// logger << "      +ge " << h << " " << leafIndex << endl;
			} else {
				for (auto u : tree[v]) {
					int uIndex = treeNode2UniverseIndex[u];
					gEdges.push_back(EdgeWeight(
						h, uIndex,
						universeVertexSet.distance(h, uIndex)
					));
					// logger << "      +ge " << h << " " << uIndex << endl;
				}
			}
		}
		// logger << "TREE: " << endl;
		// for (int i=0; i<tree.size(); i++) {
		//     logger  << " " << i << ":" << tree[i];
		// //     if (tree[i].size() == 1)
		// //         logger << "[" << cells[tree[i][0]] << "]" ;
		// //     logger
		// //         << " -> ";
		// //     for (auto a: imputedTree[i])
		// //         logger << int2AllelCoding[a] << " ";
		//     // logger << endl;
		// logger << g << endl;
		// }

		// return make_tuple(g, gNodeCells, finalCost);
		return make_tuple(gEdges, finalCost);
	}
	return make_tuple(vector<EdgeWeight>(), finalCost);
	// return make_tuple(Graph(), map<int, Cell>(), finalCost);
}

tuple<vector<EdgeWeight>, double> imputeTree(
		const vector<int>& leaves, 
		const vector<vector<int>>& tree, 
		UniverseVertexSet& universeVertexSet) {
	return imputeTree(leaves, tree, universeVertexSet, false);
}

double imputeTreeCost(
		const vector<int>& leaves, 
		const vector<vector<int>>& tree, 
		UniverseVertexSet& universeVertexSet) {
	return get<1>(imputeTree(leaves, tree, universeVertexSet, true));
}

// vector<int> pickSomeNumbers(int n, int k) {
// 	set<int> r;
// 	for (int i=0; i<k; i++) {
// 		int x = rand() % n;
// 		for (; r.find(x) != r.end(); ) 
// 			x = rand() % n;
// 		r.insert(x);
// 	}
// 	// logger << "picked " << n << " " << k << " // " << r << endl;
// 	return vector<int>(r.begin(), r.end());
// }

template<typename T>
struct SubsetIterator {
	vector<int> l; 
	int n, k;

	const vector<T>& referenceSet;

	SubsetIterator(int _n, int _k, const vector<T>& _referenceSet) : n(_n), k(_k), referenceSet(_referenceSet) {
		init();
	}

	void init() {
		l.clear();
		for (int i=0; i<k; i++)
			l.push_back(i);
	}

	bool isValid() const {
		return l.size() > 0;
	}

	void next() {
		for (int nn=n; l.size() > 0 && l.back() == nn-1; ) {
			l.pop_back();
			nn--;
		}
		if (l.size() == 0)
			return;
		l.back()++;
		while ((int)l.size() < k)
			l.push_back(l.back()+1);
		return;
	}

	vector<T> get() {
		vector<T> li;
		for (auto i: l){
			li.push_back(referenceSet[i]);
		}
		return li;
	}

	static int size(int n, int k1, int k2) {
		int sum = 0;
		for (int k=k1; k<=k2; k++) {
			int soorat = 1, makhraj = 1;
			for (int i=0; i<k; i++) {
				soorat *= (n-i);
				makhraj *= (i+1);
			}
			sum += soorat / makhraj;
		}
		return sum;
	}
};

// bool nextChoose(vector<int>& l, int n, int k) {
//     for (int nn=n; l.size() > 0 && l.back() == nn-1; ) {
//         l.pop_back();
//         nn--;
//     }
//     if (l.size() == 0)
//         return false;
//     l.back()++;
//     while (l.size() < k)
//         l.push_back(l.back()+1);
//     return true;
// }

// bool isValidChoose(vector<int>& l, int n, int k) {
//     return l.size() > 0;
// }

double treeCost(const vector<EdgeWeight>& edges, UniverseVertexSet& universeVertexSet) {
	double cst = 0;
	for (auto e: edges) {
		cst += e.w;
	}
	return cst;
}

template<typename T, typename C>
bool includesInSet(const set<T>& A, const set<T, C>& B) {
	for (auto i: A) {
		if (B.find(i) == B.end())
			return false;
	}
	return true;
}

template<typename T, typename C>
bool includesInSet(const vector<T>& A, const set<T, C>& B) {
	for (auto i: A) {
		if (B.find(i) == B.end())
			return false;
	}
	return true;
}

template<typename T, typename C>
void subtractEq(set<T, C>& A, const vector<T>& B) {
	for (auto i: B) {
		// HADI: I commented following lines from previous version, is it correct?!
		// if (A.find(i) == A.end()) {
		//     // logger << "  " << A << " in " << i << endl;
		// }
		// assert(A.find(i) != A.end());
		if (A.find(i) != A.end())
			A.erase(A.find(i));
	}
}

// template<typename V, typename T>
// void add_to_container(V& container, const T& element) {
// }

template<typename T>
void add_to_container(vector<T>& container, const T& element) {
	container.push_back(element);
}

template<typename T, typename C>
void add_to_container(set<T, C>& container, const T& element) {
	container.insert(element);
}

template<typename T>
void sort_container(vector<T>& container) {
	sort(container.begin(), container.end());
}


template<typename T, typename C>
void sort_container(set<T, C>& container) {
}


// V == vector<EdgeWeigt> or set<EdgeWeight, ???>
template<typename V>
double mstEq(V& M) {
	DisjointSetArray ds;
	set<int> vertices;
	for (auto e : M) {
		ds.add(e.v);
		ds.add(e.u);
		vertices.insert(e.v);
		vertices.insert(e.u);
	}
	int n = (int) vertices.size();
	sort_container(M);

	double cost = 0;
	V R;
	int addedEdge = 0;
	for (auto e: M) {
		if (addedEdge >= n - 1)
			break;
		if (!ds.isJoint(e.v, e.u)) {
			ds.join(e.v, e.u);
			add_to_container(R, e);
			cost += e.w;
			// R.insert(e);
		}
	}
	M = R;
	return cost;
}

tuple<vector<EdgeWeight>, double> bermenApplyCandidateTrees(UniverseVertexSet& universeVertexSet,
			const vector<int>& input, int n, const vector<EdgeWeight>& T, 
			const vector<tuple<vector<vector<int>>, vector<EdgeWeight>, vector<EdgeWeight>, vector<int>>>& candidateStack,
			bool applyOnUniverse);

double bermenCandidateTreeCost(UniverseVertexSet& universeVertexSet,
			const vector<int>& input, int n, const vector<EdgeWeight>& T, 
			const vector<tuple<vector<vector<int>>, vector<EdgeWeight>, vector<EdgeWeight>, vector<int>>>& candidateStack) {
	auto t = bermenApplyCandidateTrees(universeVertexSet,
		input, n, T, 
		candidateStack,
		false);
	return get<1>(t);
}

vector<tuple<vector<vector<int>>, vector<EdgeWeight>, vector<EdgeWeight>, vector<int>>> bermenGenerateCandidateTrees(
	UniverseVertexSet& universeVertexSet,
			const vector<int>& input, int n, int minkk, int kk, vector<EdgeWeight>& T
) {

	vector<tuple<vector<vector<int>>, vector<EdgeWeight>, vector<EdgeWeight>, vector<int>>> candidateStack;

	if (logLevel > 0)
		logger << "Init Tree: " << T << " for k=" << kk << " cost=" << bermenCandidateTreeCost(
					universeVertexSet,
						input, n, T, 
						candidateStack
				) << endl; 


	int stepsPassed = 0, stepsLBDd = 0;
	int stepsTotal = SubsetIterator<int>::size(n, minkk, kk);
	auto startTime = std::chrono::steady_clock::now();

	for (int k = minkk; k<=kk; k++) {

		GenerateAllTrees alg(k);
		vector<vector<vector<int>>> trees = alg.run();

		//fill candidateStack
		SubsetIterator<int> choose(n, k, input);


		// vector<int> l;
		// for (int i=0; i<k; i++)
		//     l.push_back(i);
		// for (int tt=0; isValidChoose(l, n, k); tt++, nextChoose(l, n, k)) {
		for (int tt=0; choose.isValid(); tt++, choose.next()) {


			if (logLevel > 0) {
				stepsPassed++;
				if (stepsPassed % 1000 == 0) {
					auto end = std::chrono::steady_clock::now();

					//Following is too time consuming, 
					// but is good to have a quote about the performance of the algorithm during steps
					double cost = bermenCandidateTreeCost(
						universeVertexSet,
							input, n, T, 
							candidateStack
					);

					logger << "step " << stepsPassed << "/" << stepsTotal << " lbd=" << stepsLBDd << " time=" << std::chrono::duration_cast<std::chrono::seconds>(end - startTime).count() << " k=" << k << " cost=" << cost << endl; 
					// return candidateStack; //BAD BAD
				}
			}

			// vector<int> li; // index of vertices
			// for (auto i: l){
			//     li.push_back(input[i]);
			// }
			vector<int> li = choose.get();
			// // logger << /*"l=" << l <<*/ " li=" << li << endl;

			sort(T.begin(), T.end());
			/**
			 * terminalTreeDiscardingEdges = bridges
			 */ 
			vector<EdgeWeight> terminalTreeRemainingEdges, terminalTreeDiscardingEdges;
			double bridgeCost = 0; // the cost of separating l from each other in T
			{
				// logger << "calculating remaining edges ... lc=" << li << endl;
				DisjointSetArray ds;
				for (auto n : input) {
					ds.add(n);
				}
				for (int i=0; i+1<(int)li.size(); i++) {
					ds.join(li[i], li[i+1]);
				}

				// logger << "  T: " << T << endl;
				for (auto &e: T) {
					if (!ds.isJoint(e.v, e.u)) {
						ds.join(e.v, e.u);
						terminalTreeRemainingEdges.push_back(e);
						// logger << "    keep " << e << endl;
					} else {
						//this should be heaviest edge between two vertices of l
						// logger << "    bridge: " << e << endl;
						terminalTreeDiscardingEdges.push_back(e);
						bridgeCost += e.w;
					}
				}
				// logger << "  " << "cost: " << bridgeCost << " " << minC << " rem: " << terminalTreeRemainingEdges <<  " discard: " << terminalTreeDiscardingEdges << endl;
			}

#ifdef USE_HAMILTONIAN_PATH_LOWER_BOUND_FOR_SMALL_STEINER_TREES_PROUNING
			double aLowerBoundForMinC = 0;
			{
				vector<EdgeWeight> clique;
				for (auto v: li) {
					for (auto u: li) 
						if (v < u) {
							clique.push_back(EdgeWeight(v, u, universeVertexSet.inputSequencesDistance(v, u)));
						}
				}
				mstEq(clique);
				double mstCost = 0;
				for (auto e: clique) {
					mstCost += e.w;
				}
				// aLowerBoundForMinC = mstCost/2;
				aLowerBoundForMinC = mstCost * steinerTreeMSTLowerBoundRate;
			}
			// in the following case, we know bridgeCost <= aLowerBoundForMinC <= minC,
			// thus, gain - bridgeCost - minC <= 0, which cause the l to not be selected,
			// thus, with this prouning we ignore calculation of the dynamic programming.
			if (bridgeCost <= aLowerBoundForMinC) {

				// vector<vector<int>> minT;
				// double minC = 1000000000;
				// for (auto tree : trees) {
				// 	double c = imputeTreeCost(li, tree, universeVertexSet);
				// 	if (c < minC) {
				// 		minC = c;
				// 		minT = tree;
				// 	}
				// }
				// if (aLowerBoundForMinC > minC) {
				// 	logger << "!!lb bc=" << bridgeCost << " lb=" << aLowerBoundForMinC << " " << minC << endl;
				// 	logger << "     list ";
				// 	int lastV = li.back();
				// 	for (auto v: li) {
				// 		logger <<  v << " " << universeVertexSet.inputSequencesDistance(v, lastV) << ":" << universeVertexSet.getVertex(v) << " ";
				// 		lastV = v;
				// 	}
				// 	logger << endl;

				// 	tuple<vector<EdgeWeight>, double> t = imputeTree(li, minT, universeVertexSet);

				// 	logger << "  tree = " << minT << " imp = " << get<0>(t) << " " << get<1>(t) << endl;
				// }
				// assert(aLowerBoundForMinC <= minC);
				// // logger << "lbd ";
				stepsLBDd++;

				continue;
			}
#endif

			double minC = 1000000000;
			vector<vector<int>> minT;

			for (auto tree : trees) {
				double c = imputeTreeCost(li, tree, universeVertexSet);
				if (c < minC) {
					minC = c;
					minT = tree;
				}
			}

			// logger << "  imputed." << endl;

			// logger << "    minT: " << minT << " cost: " << minC << endl;
			// assert(minC == get<2>(t));
			{

				double gain = bridgeCost - minC;
				// logger << "  gain = " << gain << endl;
				if (gain > EPSILON) {
					// logger << "lc=" << li << " bridge: " << bridgeCost << " c: " << minC << endl;
					// logger << "  rem: " << terminalTreeRemainingEdges << " disc: " << terminalTreeDiscardingEdges << endl;
					// logger << "  g: " << g << endl;
					// logger << "    nodeCells: " << nodeCells << endl;
					// logger << "    lc: " << lc << endl;
					// logger << "    tau edges : " << edges << endl;
					DisjointSetArray ds, ds2;
					for (auto n : input) {
						ds.add(n);
						ds2.add(n);
					}
					for (int i=0; i+1<(int)li.size(); i++) {
						ds.join(li[i], li[i+1]);
						// logger << "    join tau: " << li[i] << " " << li[i+1] << endl;
					}

					for (auto &e: T) {
						if (!ds.isJoint(e.v, e.u)) {
							ds.join(e.v, e.u);
							ds2.join(e.v, e.u);
							// logger << "    join " << e.v << " " << e.u << endl;
						} else {
							// logger << "    join! " << e.v << " " << e.u << endl;
						}
					}

					// logger << "    ds: " << ds << " ds2: " << ds2 << endl;

					//returns terminal of a component
					map<int, int> dsRoot2TerminalMap;
					for (auto t : li) {
						// logger << "      root: " << *h2VertexMap[t.hash()] << " " << *h2VertexMap[ds2.parentStar(t.hash())] << " " << endl;
						dsRoot2TerminalMap[ds2.parentStar(t)] = t;
					}

					// logger << "    dsRoot2TerminalMap: " << dsRoot2TerminalMap << endl;
					vector<EdgeWeight> terminalNewEdges;
					for (auto e: terminalTreeDiscardingEdges) {
						if (ds2.isJoint(e.v, e.u)) {
							if (logLevel > 0)
								logger << "Discarding edges are connected " << e << endl;
							ostringstream os;
							os << "Discarding edges are connected " << e;
							throw ExitException(os.str());
						}

						int t1 = dsRoot2TerminalMap[ds2.parentStar(e.v)],
							t2 = dsRoot2TerminalMap[ds2.parentStar(e.u)];
						// logger << "    terminals: " << e.c1 << "," << e.c2 << " " << t1 << "," << t2 << " / " << ds2.parentStar(e.c1) << " " << ds2.parentStar(e.c2) << endl;
						// logger << "    terminals: " << *h2VertexMap[e.c1] << "," << *h2VertexMap[e.c2] << " " << *h2VertexMap[t1] << "," << *h2VertexMap[t2] << endl;

						//We can keep the track of the tree here too!
						EdgeWeight ne = EdgeWeight(t1, t2, e.w - gain);
						terminalNewEdges.push_back(ne);

						// logger << "    ne: " << ne << endl;
					}

					// A good thing to see on logger
					
					candidateStack.push_back(make_tuple(minT, terminalTreeDiscardingEdges, terminalNewEdges, li));
					

					T = terminalTreeRemainingEdges;
					T.insert(T.end(), terminalNewEdges.begin(), terminalNewEdges.end());

					if (logLevel > 0) {
						double cost = bermenCandidateTreeCost(
							universeVertexSet,
								input, n, T, 
								candidateStack
						);

						logger << " candid: " << li << " gain=" << gain;
						logger << "  s=" << stepsPassed << " cost=" << cost << endl; 
					}

					assert((int)T.size() == n-1);

					// logger << "  T: " << T << " disc: " << terminalTreeDiscardingEdges << " new: " << terminalNewEdges << endl;

					// auto actualTree = actualTreeOfTerminalCollapsedTree(input, T);

					// logger << "  " << "cost: " << bridgeCost << " " << minC << 
					//     // " rem: " << terminalTreeRemainingEdges <<  
					//     " newCost: " << get<2>(actualTree) << 
					//     endl;

				} 
			}

			// logger << "Updated " << t << " " << minC << " " << t1.w() ;
			// logger << " : [" << l << "] G:" << g << " T1:" << t1;
			// logger << endl;
		}
			
	}


	if (logLevel > 0) {
		auto end = std::chrono::steady_clock::now();
		logger << "step " << stepsPassed << "/" << stepsTotal << " lbd=" << stepsLBDd << " time=" << std::chrono::duration_cast<std::chrono::seconds>(end - startTime).count() << endl; 
	}
	if (logLevel > 0) {
		logger << "Candidate trees calculated" << endl;
	}

	return candidateStack;
}

tuple<vector<EdgeWeight>, double> bermenApplyCandidateTrees(UniverseVertexSet& universeVertexSet,
			const vector<int>& input, int n, const vector<EdgeWeight>& T, 
			const vector<tuple<vector<vector<int>>, vector<EdgeWeight>, vector<EdgeWeight>, vector<int>>>& candidateStack,
			bool applyOnUniverse) {
	
	vector<EdgeWeight> gEdges;
	double cost = 0;

//check candidateStack
	//and apply candidateStack
	{
		// logger << "Applying candidate Stack" << endl;
		// logger << "  M: " << T << endl;

		// map<int, Cell> gMap;
		DisjointSetArray ds;

		set<EdgeWeight> E = set<EdgeWeight>(T.begin(), T.end()),
			M = E;

		for (auto n : input) {
			ds.add(n);
		}


		// for (; !candidateStack.empty(); candidateStack.pop_back()) {
		// 	auto t = candidateStack.back();
		for (int i=(int)(candidateStack.size())-1; i >=0; i--) {
			auto t = candidateStack[i];

			vector<vector<int>> tree = get<0>(t);
			vector<EdgeWeight> terminalTreeDiscardingEdges = get<1>(t), 
				terminalNewEdges = get<2>(t);
			vector<int> leavesIndices = get<3>(t);

			assert(includesInSet(M, E));
			
			// logger << "  ? " << "new: " << terminalNewEdges << " dics: " << terminalTreeDiscardingEdges << endl;
			// logger << "    M: " << M << " E: " << E << endl;
			
			// E += terminalTreeDiscardingEdges;
			//TODO: check for duplicates.
			E.insert(terminalTreeDiscardingEdges.begin(), terminalTreeDiscardingEdges.end());
			// logger << "    E += " << terminalTreeDiscardingEdges << endl;
			if (includesInSet(terminalNewEdges, M)) {
				// logger << "    A: " << tree << endl;
				//apply
				// logger << "    imputeTree ... " << leavesIndices << " " << tree << endl;
				if (applyOnUniverse) {
					tuple<vector<EdgeWeight>, double> t = imputeTree(
						leavesIndices,
						tree,
						universeVertexSet,
						false
					);
					// logger << "    imputeTree ! " << endl;
					vector<EdgeWeight> ee = get<0>(t);
					gEdges.insert(gEdges.end(), ee.begin(), ee.end());
					
					cost += get<1>(t);
				} else {
					double t = imputeTreeCost(
						leavesIndices,
						tree,
						universeVertexSet
					);
					cost += t;
				}
				// logger << "  ! new: " << terminalNewEdges << " dics: " << terminalTreeDiscardingEdges << " ee: " << ee << endl;
				// logger << "   universeVertexSet.cells.size()= " << universeVertexSet.cells.size() << endl;
				for (auto e: terminalNewEdges) {
					// logger << "    +e: " << e << endl;

					ds.join(e.v, e.u);
				}
				// logger << "  G: " << g << endl;
				// logger << "  M(f): " << M << " E: " << E << endl;

			} else {
				subtractEq(E, terminalNewEdges);
				M = E; // HADI: I changed it from previous version, is it correct?!
				// logger << "  MST ... " << endl;
				mstEq(M);
				// logger << "      done " << endl;
			}

			assert((int)M.size() == n-1);
		}


		// logger << "  applying M: " << M << endl;
		// logger << "    ds: " << ds << endl;

		for (auto e: M) {
			// logger << "  +?e:M " << e.name1 << " " << e.name2 << endl;
			if (!ds.isJoint(e.v, e.u)) {
				// logger << "  +e:M " << e.name1 << " " << e.name2 << endl;

				double ew = universeVertexSet.distance(e.v, e.u);
				gEdges.push_back(EdgeWeight(
					e.v, e.u,
					ew
				));
				cost += ew;


				ds.join(e.v, e.u);
			}
		}


	}

	return make_tuple(gEdges, cost);
}

/**
 * Implementation of: mproved Approximations for the Steiner Tree Problem
 * Authors: Piotr Berman and Viswanathan Ramaiyer
 */ 
tuple<vector<EdgeWeight>, double> optimizeTree(
		UniverseVertexSet& universeVertexSet,
		const vector<int>& input, int minkk, int kk) {
	/**
	 * T always contain a tree between terminal nodes.
	 */
	// vector<EdgeExt> T = initTree(input);

	assert(kk < MAXTREELEAFS);

	double distanceMatrix[MAX_SEQUENCE][MAX_SEQUENCE];

	vector<EdgeWeight> T;
	for (auto v: input) {
		for (auto u: input) {
			double w = universeVertexSet.distance(v, u);
			T.push_back(EdgeWeight(v, u, w));
			distanceMatrix[v][u] = distanceMatrix[u][v] = w;
		}
	}
	//double initTreeCost = 
	mstEq(T);

	int n = (int)input.size();

	//Tree, distartedEdges, newEdges, treeLeaveIndices(tau)
	vector<tuple<vector<vector<int>>, vector<EdgeWeight>, vector<EdgeWeight>, vector<int>>> candidateStack = 
		bermenGenerateCandidateTrees(universeVertexSet,
			input, n, minkk, kk, T);

	tuple<vector<EdgeWeight>, double> t = bermenApplyCandidateTrees(
		universeVertexSet,
			input, n, T, 
			candidateStack,
			true
	);
	if (logLevel > 0) {
		logger << "Trees applied" << endl; 
	}


	// double cost = treeCost(gEdges, universeVertexSet);
	// logger << "result G: " << gEdges << " " << " cost: " << cost << endl;
	return t;

	// return actualTreeOfTerminalCollapsedTree(input, T);
}

// Removes degree two edges which does not provide any other information regarding input vertices
tuple<vector<EdgeWeight>, vector<int>> compressGraph(UniverseVertexSet& universeVertexSet, const vector<EdgeWeight>& edges, const vector<int>& inputCells) {
	DisjointSetArray ds;
	for (int i=0; i<universeVertexSet.size(); i++) {
		ds.add(i);
	}

	for (auto e: edges) {
		if (e.w == 0) {
			// assert(universeVertexSet.distance(e.v, e.u) == 0);
			// logger << "  0-edge: " << e.v << " " << e.u << ": " << universeVertexSet.getVertex(e.v) << " " << universeVertexSet.getVertex(e.u) << endl;
			ds.join(e.v, e.u);
		}
	}

	map<int, int> minOfSet;
	for (int i=0; i<universeVertexSet.size(); i++) {
		int p = ds.parentStar(i);
		if (minOfSet.find(p) == minOfSet.end())
			minOfSet[p] = p;
		minOfSet[p] = min(minOfSet[p], i);
	}

	// logger << "minOfSet: " << minOfSet << endl;

	vector<EdgeWeight> E;
	for (auto e: edges) {
		if (e.w != 0) {
			E.push_back(EdgeWeight(
				minOfSet[ds.parentStar(e.v)],
				minOfSet[ds.parentStar(e.u)],
				e.w
			));
		}
	}

	set<int> V;
	for (auto m: minOfSet) {
	   V.insert(m.second); 
	}

	//We should add input vertices, even there are more than one with equal sequences.
	// logger << "  minOfSet: input cells " << inputCells << endl;
	for (auto v: inputCells) {
		if (V.find(v) == V.end()) {
			// logger << "  minOfSet: Add input " << v << endl;
			V.insert(v);
			E.push_back(EdgeWeight(
				v,
				minOfSet[ds.parentStar(v)],
				0
			));
		}
	}

	return make_tuple(E, vector<int>(V.begin(), V.end()));
}

void printResultAsGraph(ostream& os, UniverseVertexSet& universeVertexSet, const vector<EdgeWeight>& edges, double cost, const vector<int>& cells, const map<int, Cell>& imputation) {

	tuple<vector<EdgeWeight>, vector<int>> tt = compressGraph(universeVertexSet, edges, cells);

	if (logLevel > 0)
		logger << "Tree compressed" << endl;

	vector<EdgeWeight> e = get<0>(tt);
	vector<int> v = get<1>(tt);

	set<int> inputCells;
	for (auto c: cells) {
		inputCells.insert(c);
	}
	
	os << v.size() << endl;
	for (int j=0; j<(int)v.size(); j++) {
		int i = v[j];
		os << i << " " << 
			(inputCells.find(i) == inputCells.end() ? 0 : 1) << " " << 
			universeVertexSet.getVertex(i).toString() << " " << 
			(imputation.find(i) != imputation.end() ? imputation.find(i)->second.toString() : "-") << endl;
	}
	os << e.size() << endl;
	for (auto ee: e) {
		os << ee.v << " " << ee.u << " " << ee.w << endl;
	}
}

void rewriteCellXValuesTo(Cell& to, const Cell& from) {
	for (int i=0; i<to.size(); i++) {
		if (to.s[i] == 'X' && from.s[i] != 'X')
			to.s[i] = from.s[i];
	}
}


map<int,Cell> calculateImputation(UniverseVertexSet& universeVertexSet, const vector<EdgeWeight>& edges, const vector<int>& inputCells) {
	set<int> input(inputCells.begin(), inputCells.end());
	map<int, Cell> imputation;
	for (auto v: inputCells) {
		imputation[v] = universeVertexSet.getVertex(v);
	}
	for (auto e: edges) {
		if (input.find(e.v) != input.end()) {
			rewriteCellXValuesTo(imputation[e.v], universeVertexSet.getVertex(e.u));
			// logger << "imp " << e.v << " -> " << e.u << endl;
		} else if (input.find(e.u) != input.end()) {
			rewriteCellXValuesTo(imputation[e.u], universeVertexSet.getVertex(e.v));
			// logger << "imp " << e.u << " -> " << e.v << endl;
		}
	}
	Cell aCell;
	for (int i=0; i<universeVertexSet.length(); i++) {
		aCell.append('A');
	}
	for (auto v: inputCells) {
		rewriteCellXValuesTo(imputation[v], aCell);
	}
	return imputation;
}
