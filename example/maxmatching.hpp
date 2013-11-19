#ifndef __MAXMATCHING_HPP
#include <algorithm>
#include <lpsolve/lp_lib.h>
#include <map>
#include <vector>
#include <utility>
#include <set>
#include <iostream>

using namespace std;

map<int, int> getMaxMatching(vector<int> labels1, vector<int> labels2);
map<int, int> getWeightedMaxMatching(vector<int> labels1, vector<int> labels2, vector<double> weights);
map<int, int> getMaxMatchingConsistentWithOldMatching(vector<int> labels1, vector<int> labels2, map<int, int> oldmatchings);


//map<int, int> getMinWeightLeftSidePerfectMatching(vector< pair<int, int> > nodePairs, vector<double> edgeWeights);

#include "maxmatching_impl.hpp"
#define __MAXMATCHING_HPP
#endif /* __MAXMATCHING_HPP */
