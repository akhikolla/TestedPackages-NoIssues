#include "TabProbsResults.h"



TabProbsResults::TabProbsResults(int N, int kr, int J, int kc)
{
	this->_tabprobaV.zeros(N, kr);
	this->_tabprobaW.zeros(J, kc);
}


TabProbsResults::~TabProbsResults()
{
}
