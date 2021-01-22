#pragma once

#include <RcppArmadillo.h>
#include <limits>

using namespace std;
using namespace arma;

class TabProbsResults
{
public:
	TabProbsResults(int N, int kr, int J, int kc);
	~TabProbsResults();
	mat _tabprobaV;
	mat _tabprobaW;
};

