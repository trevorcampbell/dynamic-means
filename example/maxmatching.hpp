#ifndef __MAXMATCHING_HPP
#include <algorithm>
#include <lpsolve/lp_lib.h>
#include <map>
#include <vector>
#include <utility>
#include <set>
#include <iostream>

using namespace std;
class MaxMatching{
	public:
		MaxMatching();
		map<int, int> getMaxMatching(vector<int> labels1, vector<int> labels2);
		map<int, int> getWeightedMaxMatching(vector<int> labels1, vector<int> labels2);
		map<int, int> getConsistentMaxMatching(vector<int> labels1, vector<int> labels2);
	private:
		map<int, int> oldmatchings;
		Simplex s;
};

#include "maxmatching_impl.hpp"
#define __MAXMATCHING_HPP
#endif /* __MAXMATCHING_HPP */
