#ifndef __MAXMATCHING_HPP
#include <algorithm>
#include <map>
#include <vector>
#include <utility>
#include <set>
#include <iostream>

using namespace std;
class MaxMatching{
	public:
		MaxMatching();
		map<int, int> getLabelMatching(vector<int> labels1, vector<int> labels2);
		map<int, int> getWeightedMaxMatching(vector<int> labels1, vector<int> labels2);
		map<int, int> getConsistentLabelMatching(vector<int> labels1, vector<int> labels2);
	private:
		map<int, int> oldmatchings;
		Simplex s;
		class InvalidLabelsSizeException{
			public:
				InvalidLabelsSizeException(int s1, int s2){
					std::cout << "The two label vectors have an invalid size for max matching!" << std::endl;
					std::cout << "size of labels1: " << s1 << " size of labels2: " << s2 << std::endl;
				}
				InvalidLabelsSizeException(int s1, int s2, int w){
					std::cout << "The two label vectors/weight vector have an invalid size for max matching!" << std::endl;
					std::cout << "size of labels1: " << s1 << " size of labels2: " << s2 << " size of weights: " w <<  std::endl;
				}
		};
};



/*
 *
 * This is testing code for the getMaxMatching function below
 * getMaxMatching has been verified to work.
int main(int argc, char** argv){
	vector<int> l1, l2;

	//l1 labels
	for (int i = 0; i < 3; i++){
		l1.push_back(0);
	}
	for (int i = 0; i < 7; i++){
		l1.push_back(5);
	}
	for (int i = 0; i < 4; i++){
		l1.push_back(3);
	}
	for (int i = 0; i < 5; i++){
		l1.push_back(2);
	}


	//l2 labels
	for (int i = 0; i < 3; i++){
		l2.push_back(8);
	}
	for (int i = 0; i < 7; i++){
		l2.push_back(4);
	}
	for (int i = 0; i < 4; i++){
		l2.push_back(0);
	}
	for (int i = 0; i < 5; i++){
		l2.push_back(1);
	}

	map<int, int> matching = getMaxMatching(l1, l2);
	cout << "Matchings:" << endl;
	for (map<int, int>::iterator it = matching.begin(); it != matching.end(); it++){
		cout << it->first << ":" << it->second << endl;
	}
}*/

//this function finds the best correspondance between labels1 and labels2
//not *all* labels in either labels1 or labels2 is guaranteed to be used
//the map will not contain a key for any unused label in labels1
//the map will not map any key to an unused label in labels2

//further, if any label in labels1 or labels2 is < 0, those elements
//are removed from the proceedings at the start (it is assumed that negative labels
//stand for *unknown* / *unlabelled* data that is "wrong" by default)



#include "maxmatching_impl.hpp"
#define __MAXMATCHING_HPP
#endif /* __MAXMATCHING_HPP */
