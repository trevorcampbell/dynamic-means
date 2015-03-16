#ifndef __DYNMEANS_IMPL_HPP
template<class Vec>
DynMeans<Vec>::DynMeans(double lambda, double Q, double tau, bool verbose){
	this->verbose = verbose;
	this->lambda = lambda;
	this->Q = Q;
	this->tau = tau;
	this->ages.clear();
	this->oldprms.clear();
	this->oldprmlbls.clear();
	this->observations.clear();
	this->weights.clear();
	this->nextLbl = 0;
}

template<class Vec>
DynMeans<Vec>::~DynMeans(){
}

template<class Vec>
void DynMeans<Vec>::reset(){
	this->ages.clear();
	this->oldprms.clear();
	this->oldprmlbls.clear();
	this->observations.clear();
	this->weights.clear();
	this->nextLbl = 0;
}

//This function is used when sampling parameters - it returns a vector of the observations in the next cluster,
//along with the index of that cluster.
//If this is the last parameter to be sampled, the function returns true; otherwise, false.
template<class Vec>
std::vector<Vec> DynMeans<Vec>::getObsInCluster(int idx, std::vector<int> lbls){
	//std::cout << "Getting obs set in next cluster" << std::endl;
	std::vector<Vec> obsInCluster;
	obsInCluster.reserve(lbls.size());
	//std::cout << "Getting obs for param with idx " << idx << std::endl;
	for (int i = 0; i < lbls.size(); i++){
		if (lbls[i] == idx){
			obsInCluster.push_back(this->observations[i]);
		}
	}
	return obsInCluster;
}

template<class Vec>
std::vector<int> DynMeans<Vec>::updateState(std::vector<int> lbls, std::vector<int> cnts, std::vector<Vec> prms){
	this->oldprms = prms;
	std::vector<int> outLbls; //stores the label output using oldprmlbls
	//update the weights/ages
	for (int i = 0; i < prms.size(); i++){
		if (i < this->weights.size() && cnts[i] > 0){
			//this is an instantiated cluster from a previous time; set age to 0 and update weights
			this->weights[i] = 1.0/(1.0/this->weights[i] + this->ages[i]*this->tau) + cnts[i];
			this->ages[i] = 0;
		} else if (i >= this->weights.size() ) {
			//new cluster
			//push back a 0 for the age, and set the weight to the number of observations
			this->ages.push_back(0);
			this->weights.push_back(cnts[i]);
			this->oldprmlbls.push_back(this->nextLbl++);
		}
		//increment the age for all clusters (including old clusters with i < weights.size() && cnts = 0)
		this->ages[i]++;
	}
	for(int i = 0; i < lbls.size(); i++){
		outLbls.push_back(this->oldprmlbls[lbls[i]]);
	}
	//now that all is updated, check to see which clusters are permanently dead
	for(int i = 0; i < this->oldprms.size(); i++){
		if (this->ages[i]*this->Q > this->lambda){
			this->oldprms.erase(this->oldprms.begin()+i);
			this->oldprmlbls.erase(this->oldprmlbls.begin()+i);
			this->weights.erase(this->weights.begin()+i);
			this->ages.erase(this->ages.begin()+i);
			i--;
		}
	}
	return outLbls;
}

template<class Vec>
void DynMeans<Vec>::cluster(std::vector<Vec>& newobservations, int nRestarts, 
		std::vector<int>& finalLabels, std::vector<Vec>& finalParams, double& finalObj, double& tTaken){
	timeval tStart;
	gettimeofday(&tStart, NULL);


	//set the new obs
	this->observations = newobservations;

	if (newobservations.size() == 0){
		std::cout << "libdynmeans: ERROR: newobservations is empty" << std::endl;
		return;
	}
	if (nRestarts <= 0){
		std::cout << "libdynmeans: ERROR: Cannot have nRestarts <= 0" << std::endl;
		return;
	}

	//generate the orderings
	std::vector< std::vector<int> > randOrderings;
	std::vector<int> tmpindices;
	for (int i = 0; i < newobservations.size(); i++){
		tmpindices.push_back(i);
	}
	std::srand( unsigned( std::time(0) ) );
	for (int i = 0; i < nRestarts; i++){
		//cout << "Generating random order " << i+1 << "/" << nOrderings << "               \r" << flush;
		random_shuffle(tmpindices.begin(), tmpindices.end());
		randOrderings.push_back(tmpindices);
	}

	//stuff for storing best clustering
	finalObj =  std::numeric_limits<double>::max();
	std::vector<int> finalCnts;

	if (verbose){
		std::cout << "libdynmeans: Clustering " << newobservations.size() << " datapoints with " << nRestarts << " restarts." << std::endl;
	}
	for (int i = 0; i < nRestarts; i++){

		//create working variables
		std::vector<Vec> prms;
		std::vector<int> cnts;
		std::vector<int> lbls;

		for (int j = 0; j < this->oldprms.size(); j++){
			prms.push_back(this->oldprms[j]); //just placeholders for updated parameters if the old ones get instantiated
			cnts.push_back(0); //start with count 0
		}

		//Initialization: no label on anything
		for (int j = 0; j < this->observations.size(); j++){
			lbls.push_back(-1); // start with no labels on anything
		}

		double obj, prevobj;
		obj = prevobj = std::numeric_limits<double>::max();

		do {
			//do the kmeans iteration
			prevobj = obj;
			this->assignObservations(randOrderings[i], lbls, cnts, prms);
			obj = this->setParameters(lbls, cnts, prms);
			if (obj > prevobj){
				std::cout << "Error: obj > prevobj - monotonicity violated! Check your distance/set parameter functions..." << std::endl;
				std::cout << "obj: " << obj << " prevobj: " << prevobj << std::endl;
			}
			//std::cout << "KMeans -- Restart: " << i+1 << "/" << nRestarts << " Iteration: " << iter << " Objective: " << obj << std::endl;//"\r" << std::flush;
			if (verbose){
				int numinst = 0;
				for (int ii = 0; ii < cnts.size(); ii++){
					if (cnts[ii] > 0){
						numinst++;
					}
				}
				int numnew = prms.size() - this->oldprms.size();
				int numoldinst = numinst - numnew;
				int numolduninst = cnts.size() - numinst;
			std::cout << "libdynmeans: Trial: " << i+1 << "/" << nRestarts << " Objective: " << obj << " Old Uninst: " << numolduninst  << " Old Inst: " << numoldinst  << " New: " << numnew <<  "                   \r" << std::flush; 
			}
		} while(prevobj > obj);

		if (obj < finalObj){
			finalObj = obj;
			finalParams = prms;
			finalLabels = lbls;
			finalCnts = cnts;
		}
	}
	if (verbose){
		int numinst = 0;
		for (int ii = 0; ii < finalCnts.size(); ii++){
			if (finalCnts[ii] > 0){
				numinst++;
			}
		}
		int numnew = finalParams.size() - this->oldprms.size();
		int numoldinst = numinst - numnew;
		int numolduninst = finalCnts.size() - numinst;
		std::cout << "libdynmeans: Done clustering. Min Objective: " << finalObj << " Old Uninst: " << numolduninst  << " Old Inst: " << numoldinst  << " New: " << numnew <<  std::endl;
	}
	//update the stored results to the one with minimum cost
	finalLabels = this->updateState(finalLabels, finalCnts, finalParams);
	timeval tCur;
	gettimeofday(&tCur, NULL);
	tTaken = (double)(tCur.tv_sec - tStart.tv_sec) + (double)(tCur.tv_usec - tStart.tv_usec)/1.0e6;

	return;
}


template<class Vec>
void DynMeans<Vec>::assignObservations(std::vector<int> assgnOrdering, std::vector<int>& lbls, std::vector<int>& cnts, std::vector<Vec>& prms){
	for (int i = 0; i < assgnOrdering.size(); i++){
		//get the observation idx from the random ordering
		int idx = assgnOrdering[i];

		//store the old lbl for possibly deleting clusters later
		int oldlbl = lbls[idx];

		//calculate the distances to all the parameters
		int minind = 0;
		double mindistsq = std::numeric_limits<double>::max();
		for (int j = 0; j < prms.size(); j++){
			double tmpdistsq = (prms[j] - this->observations[idx]).squaredNorm();
			if (cnts[j] == 0){//the only way cnts can get to 0 is if it's an old parameter
				double gamma = 1.0/(1.0/this->weights[j] + this->ages[j]*this->tau);
				tmpdistsq = gamma/(1.0+gamma)*tmpdistsq + this->ages[j]*this->Q;
			}
			if(tmpdistsq < mindistsq){
				minind = j;
				mindistsq = tmpdistsq;
			}
		}

		//if the minimum distance is stil greater than lambda + startup cost, start a new cluster
		if (mindistsq > this->lambda){
			prms.push_back(this->observations[idx]);
			lbls[idx] = prms.size()-1;
			cnts.push_back(1);
		} else {
			if (cnts[minind] == 0){ //if we just instantiated an old cluster
						//update its parameter to the current timestep
						//so that upcoming assignments are valid
				double gamma = 1.0/(1.0/this->weights[minind] + this->ages[minind]*this->tau);
				prms[minind] = (this->oldprms[minind]*gamma + this->observations[idx])/(gamma + 1);
			}
			lbls[idx] = minind;
			cnts[minind]++;
		}


		//if obs was previously assigned to something, decrease the count the clus it was assigned to
		//we do cluster deletion *after* assignment to prevent corner cases with monotonicity
		if (oldlbl != -1){
			cnts[oldlbl]--;
			//if this cluster now has no observations, but was a new one (no age recording for it yet)
			//remove it and shift labels downwards
			if (cnts[oldlbl] == 0 && oldlbl >= this->oldprms.size()){
				prms.erase(prms.begin() + oldlbl);
				cnts.erase(cnts.begin() + oldlbl);
				for (int j = 0; j < lbls.size(); j++){
					if (lbls[j] > oldlbl){
						lbls[j]--;
					}
				}
			} else if (cnts[oldlbl] == 0){//it was an old parameter, reset it to the oldprm
				prms[oldlbl] = this->oldprms[oldlbl];
			}
		}
	}
	return;	
}

template<class Vec>
double DynMeans<Vec>::setParameters(std::vector<int>& lbls, std::vector<int>& cnts, std::vector<Vec>& prms){
	double objective = 0;
	for (int i = 0; i < prms.size(); i++){
		if (cnts[i] > 0){
			//add cost for new clusters - lambda
			// or add cost for old clusters - Q
			if (i < this->oldprms.size()){
				objective += this->Q*this->ages[i];
			} else {
				objective += this->lambda;
			}
			std::vector<Vec> obsInClus = this->getObsInCluster(i, lbls);
			Vec tmpvec = obsInClus[0];
			for (int j = 1; j < obsInClus.size(); j++){
				tmpvec = tmpvec + obsInClus[j];
			}
			tmpvec = tmpvec / obsInClus.size();
			if (i < this->oldprms.size()){ //updating an old param
				double gamma = 1.0/(1.0/this->weights[i] + this->ages[i]*this->tau);
				prms[i] = (this->oldprms[i]*gamma + tmpvec*cnts[i])/(gamma + cnts[i]);
				//add parameter lag cost
				double tmpsqdist = (prms[i] - this->oldprms[i]).squaredNorm();
				objective += gamma*tmpsqdist;
			} else { //just setting a new param
				prms[i] = tmpvec;
				//no lag cost for new params
			}
			//get cost for prms[i]
			for (int j = 0; j < obsInClus.size(); j++){
				objective += (prms[i] - obsInClus[j]).squaredNorm();
			}
		}
	}
	return objective;
}


#define __DYNMEANS_IMPL_HPP
#endif /* __DYNMEANS_IMPL_HPP */
