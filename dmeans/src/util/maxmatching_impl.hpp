#ifndef __MAXMATCHING_IMPL_HPP

void MaxMatching::resetOldMatchings(){
	oldmatchings.clear();
}

double MaxMatching::getObjective(){
	return objective;
}

void MaxMatching::pruneInvalidLabelPairs(vector<int>& labels1, vector<int>& labels2, vector<double>& weights){
	//prune data with negative labels
	for (uint64_t i = 0; i < labels1.size(); i++){
		if(labels1[i] < 0 || labels2[i] < 0){
			labels1.erase(labels1.begin()+i);
			labels2.erase(labels2.begin()+i);
			weights.erase(weights.begin()+i);
			i--;
		}
	}
}

void MaxMatching::pruneInconsistentLabelPairs(vector<int>& labels1, vector<int>& labels2, vector<double>& weights){
	//prune labels already in the old matching
	for (map<int, int>::iterator it = oldmatchings.begin(); it != oldmatchings.end(); ++it){
		int l1 = it->first, l2 = it->second;
		for (uint64_t i = 0; i < labels1.size(); i++){
			if(labels1[i] == l1 || labels2[i] == l2){
				labels1.erase(labels1.begin()+i);
				labels2.erase(labels2.begin()+i);
				weights.erase(weights.begin()+i);
				i--;
			}
		}
	}
}

set<int> MaxMatching::getUniqueLabels(const vector<int>& labels){
	set<int> lset;
	for (uint64_t i = 0; i < labels.size(); i++){
		lset.insert(labels[i]);
	}
	return lset;
}

void MaxMatching::getMaps(const vector<int>& labels1, const set<int>& l1set, const vector<int>& labels2, const set<int>& l2set, const vector<double>& weights,
							map< pair<int, int>, int>& varMap, map<int, pair<int, int> >& invvarMap, map<int, double>& weightMap){
	int k = 0;
	for (set<int>::iterator it1 = l1set.begin(); it1 != l1set.end(); it1++){
	for (set<int>::iterator it2 = l2set.begin(); it2 != l2set.end(); it2++){
			//variable maps
			varMap[pair<int, int>(*it1, *it2)] =  k;
			invvarMap[k] = pair<int,int>(*it1, *it2);
			//weight map
			//increment the weight of this variable (graph edge) based on the number
			//of matching observations between the two
			double weight = 0;
			for (uint64_t i = 0; i < labels1.size(); i++){
				if (labels1[i] == *it1 && labels2[i] == *it2){
					weight+= weights[i];
				}
			}
			weightMap[k] = weight;
			//cout << "Var map (" << *it1 << "," << *it2 << ") = " << k << ", weight map (" << k << ") = " << weight << endl;
			k++;
	}
	}
}

map<int, int> 
MaxMatching::getMaxMatching(vector<int> labels1, vector<int> labels2, vector<double> weights){
	if (labels1.size() != labels2.size() || labels1.size() == 0 || (weights.size() > 0 && weights.size() != labels1.size())){
		throw InvalidLabelsSizeException(labels1.size(), labels2.size(), weights.size());
	}
	if (weights.size() == 0){
		weights.resize(labels1.size(), 1);
	}
	this->pruneInvalidLabelPairs(labels1, labels2, weights);

	/*cout << "Labels 1:"<< endl;
	for (uint64_t i =0 ; i < labels1.size(); i++){
		cout << labels1[i] << " ";
	}
	cout << endl;
	cout << "Labels 2:"<< endl;
	for (uint64_t i =0 ; i < labels2.size(); i++){
		cout << labels2[i] << " ";
	}
	cout << endl;*/
	//first create ordered sets of labels1 and labels2
	//this pulls out the unique labels
	set<int> l1set = this->getUniqueLabels(labels1);
	set<int> l2set = this->getUniqueLabels(labels2);
	/*cout << "Unique Labels 1" << endl;
	for (set<int>::iterator it = l1set.begin(); it != l1set.end(); it++){
		cout << *it << " ";
	}
	cout << endl;
	cout << "Unique Labels 2" << endl;
	for (set<int>::iterator it = l2set.begin(); it != l2set.end(); it++){
		cout << *it << " ";
	}
	cout << endl;*/
	//now create maps from the label indices in l1 and l2 to variable indices in the LP 
	//and simultaneously create the weight map
	map<pair<int, int>, int> varMap;
	map<int, pair<int, int> > invvarMap;
	map<int, double> weightMap;
	this->getMaps(labels1, l1set, labels2, l2set, weights, varMap, invvarMap, weightMap);
	

	//construct the linear program
	optimization::Simplex splx("MaxMatching");
	//add variables
	for(uint64_t i = 0; i < varMap.size(); i++){
		splx.add_variable(new optimization::Variable(&splx, std::to_string(i).c_str()));
	}

	//add constraints
	pilal::Matrix coeffs(1, varMap.size(), 0);
	//constrant type 1: the sum of outgoing edges from each of the A vertices = 1
	for (set<int>::iterator it1 = l1set.begin(); it1 != l1set.end(); it1++){
		for (set<int>::iterator it2 = l2set.begin(); it2 != l2set.end(); it2++){
			coeffs(0, varMap[pair<int,int>(*it1, *it2)]) = 1;
		}
		splx.add_constraint(optimization::Constraint(coeffs, optimization::ConstraintType::CT_LESS_EQUAL, 1));
	}
	//constraint type 2: the sum of incoming edges to each fo the B vertices = 1
	coeffs.empty();
	for (set<int>::iterator it2 = l2set.begin(); it2 != l2set.end(); it2++){
		for (set<int>::iterator it1 = l1set.begin(); it1 != l1set.end(); it1++){
			coeffs(0, varMap[pair<int,int>(*it1, *it2)]) = 1;
		}
		splx.add_constraint(optimization::Constraint(coeffs, optimization::ConstraintType::CT_LESS_EQUAL, 1));
	}
	//add objective and set to maximization mode
	coeffs.empty();
	for (map<pair<int, int>, int>::iterator it = varMap.begin(); it != varMap.end(); it++){
		coeffs(0, it->second) = weightMap[it->second];
	}
    splx.set_objective_function(optimization::ObjectiveFunction(optimization::ObjectiveFunctionType::OFT_MAXIMIZE, coeffs));

	//solve
	splx.solve();
	if (splx.must_be_fixed() || !splx.has_solutions() || (splx.has_solutions() && splx.is_unlimited())){
		throw LinearProgrammingException(splx.must_be_fixed(), splx.has_solutions(), splx.is_unlimited());
    }
	coeffs = splx.get_solution();
	this->objective = splx.get_objective();

	/*cout << "raw output: " << endl;
	for (uint64_t i = 0; i <invvarMap.size(); i++){
		cout << varweights[i] << endl;
	}*/
	//create the output
	map<int, int> retMap;
	set<int> usedL1Labels, usedL2Labels;
	for (uint64_t i = 0; i < invvarMap.size(); i++){
		if (fabs(coeffs(i)-1.0) < 1e-6){
			//ensure no duplicate L1 labels
			if(usedL1Labels.find(invvarMap[i].first) != usedL1Labels.end()){
				throw InvalidMatchingException(invvarMap[i].first, retMap[invvarMap[i].first], invvarMap[i].second, true);
			}
			//ensure no duplicate L2 labels
			if(usedL2Labels.find(invvarMap[i].second) != usedL2Labels.end()){
				int l1label1 = -1;
				for (map<int, int>::iterator it = retMap.begin(); it != retMap.end(); it++){
					if (it->second == invvarMap[i].second){
						l1label1 = it->first;
						break;
					}
				}
				throw InvalidMatchingException(invvarMap[i].second, l1label1, invvarMap[i].first, false);
			}
			retMap[invvarMap[i].first] = invvarMap[i].second;
			usedL1Labels.insert(invvarMap[i].first);
			usedL2Labels.insert(invvarMap[i].second);
		}
	}
	//output
	return retMap;
}


map<int, int> MaxMatching::getMaxConsistentMatching(vector<int> labels1, vector<int> labels2, vector<double> weights){
	if (labels1.size() != labels2.size() || labels1.size() == 0 || (weights.size() > 0 && weights.size() != labels1.size())){
		throw InvalidLabelsSizeException(labels1.size(), labels2.size(), weights.size());
	}
	if (weights.size() == 0){
		weights.resize(labels1.size(), 1);
	}
	this->pruneInconsistentLabelPairs(labels1, labels2, weights);

	map<int, int> retMap = this->getMaxMatching(labels1, labels2, weights);

	//update the oldMatchings
	for (map<int, int>::iterator it = retMap.begin(); it != retMap.end(); ++it){
		oldmatchings[it->first] = it->second;
	}

	//output the updated oldmatchings
	return oldmatchings;
}

#define __MAXMATCHING_IMPL_HPP
#endif
