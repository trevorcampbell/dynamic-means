//#include <ros/ros.h>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <list>
#include <Eigen/Dense>
#include <sys/time.h>
#include <cmath>
#include <random>

using namespace std;

typedef Eigen::Vector2d V2d;

//function declarations
double computeAccuracy(vector<int> labels1, vector<int> labels2, map<int, int> matchings);
void birthDeathMotionProcesses(vector<V2d>& clusterCenters, vector<bool>& aliveClusters, double birthProbability, double deathProbability, double motionStdDev);
void generateData(vector<V2d> clusterCenters, vector<bool> aliveClusters, int nDataPerClusterPerStep, double likelihoodstd, vector<V2d>& clusterData, vector<int>& trueLabels);

//random number generator
mt19937 rng;//uses the same seed every time, 5489u


class VectorGraph{
	public:
		std::vector<V2d> data, oldprms;
		std::vector<int> oldprmlbls;
		VectorGraph(){
			data.clear();
			oldprms.clear();
			oldprmlbls.clear();
		}
		void updateData(std::vector<V2d> data){
			this->data = data;
		}
		void updateOldParameters(std::vector<V2d> data, std::vector<int> lbls, std::vector<double> gammas, std::vector<int> prmlbls){
			std::vector<V2d> updatedoldprms;
			for (int i = 0; i < prmlbls.size(); i++){
				//if there is no data assigned to this cluster, must be old/uninstantiated
				if (find(lbls.begin(), lbls.end(), prmlbls[i]) == lbls.end()){
					int oldidx = distance(oldprmlbls.begin(), find(oldprmlbls.begin(), oldprmlbls.end(), prmlbls[i]));
					updatedoldprms.push_back(oldprms[oldidx]);
				//if the label is not in oldprmlbls, must be a new cluster
				//furthermore, must have at least one label in lbls = prmlbls[i]
				} else if (find(oldprmlbls.begin(), oldprmlbls.end(), prmlbls[i]) == oldprmlbls.end()) {
					V2d tmpprm = V2d::Zero();
					int tmpcnt = 0;
					for (int j =0 ; j < lbls.size(); j++){
						if (lbls[j] == prmlbls[i]){
							tmpprm += data[j];
							tmpcnt++;
						}
					}
					tmpprm /= tmpcnt;
					updatedoldprms.push_back(tmpprm);
				//old instantiated cluster
				} else {
					int oldidx = distance(oldprmlbls.begin(), find(oldprmlbls.begin(), oldprmlbls.end(), prmlbls[i]));
					V2d tmpprm = gammas[i]*oldprms[oldidx];
					int tmpcnt = 0;
					for (int j =0 ; j < lbls.size(); j++){
						if (lbls[j] == prmlbls[i]){
							tmpprm += data[j];
							tmpcnt++;
						}
					}
					tmpprm /= (gammas[i]+tmpcnt);
					updatedoldprms.push_back(tmpprm);
				}
			}
			this->oldprms = updatedoldprms;
			this->oldprmlbls = prmlbls;
		}
		double diagSelfSimDD(const int i) const{
			return data[i].transpose()*data[i];
		}
		double offDiagSelfSimDD(const int i) const{
			return 0;
		}
		double selfSimPP(const int i) const{
			return oldprms[i].transpose()*oldprms[i];
		}
		double simDD(const int i, const int j) const{
			return data[i].transpose()*data[j];
		}
		double simDP(const int i, const int j) const{
			return data[i].transpose()*oldprms[j];
		}
		int getNodeCt(const int i) const{
			return 1;
		}
		int getNNodes() const {
			return data.size();
		}
		int getNOldPrms() const {
			return oldprms.size();
		}
};

#include <dynmeans/kerndynmeans.hpp>
#include "maxmatching.hpp"



int main(int argc, char** argv){

	//generates clusters that jump around on the domain R^2
	//they move stochastically with a normal distribution w/ std deviation 0.05
	//they die stochastically with probability 0.05 at each step
	//a cluster is created stochastically with probability 0.10 at each step in the area [0,1]x[0,1]
	
	//note: when computing accuracy, ***proper tracking is taken into account***
	//i.e. if true label 1 is matched to learned label 3, then that matching is fixed from that point forward
	//later pairings that attempt to match true label 1 to something else
	//or learned label 3 to something else are discarded
	

	//constants
	//play with these to change the data generation process
	double birthProbability = 0.10;
	double deathProbability = 0.05;
	double motionStdDev = 0.05;
	double clusterStdDev = 0.05;
	int nDataPerClusterPerStep = 10; // datapoints generated for each cluster per timestep
	int nSteps = 100;//run the experiment for nSteps steps
	int initialClusters = 4; //the number of clusters to start with

	
	//data structures for the cluster centers
	vector<V2d> clusterCenters;
	vector<bool> aliveClusters;

	//start with some number of initial clusters
	uniform_real_distribution<double> uniformDist01(0, 1);
	for (int i = 0; i < initialClusters; i++){
			V2d newCenter;
			newCenter(0) = uniformDist01(rng);
			newCenter(1) = uniformDist01(rng);
			clusterCenters.push_back(newCenter);
			aliveClusters.push_back(true);
	}

	//the Dynamic Means object
	//play with lambda/Q/tau to change Dynamic Means' performance
	double lambda = 0.05;
	double T_Q = 6.8;
	double K_tau = 1.01;
	double Q = lambda/T_Q;
	double tau = (T_Q*(K_tau-1.0)+1.0)/(T_Q-1.0);
	int nRestarts = 10;
	int nCoarsest = 10;
	KernDynMeans<VectorGraph> kdm(lambda, Q, tau, true);
	VectorGraph vgr;

	//run the experiment
	double cumulativeAccuracy = 0;//stores the accuracy accumulated for each step
	map<int, int> matchings;//stores the saved matchings from previous timesteps
							//enables proper label tracking (see note at line 27)
	for (int i = 0; i < nSteps; i++){
		//----------------------------
		//birth/death/motion processes
		//----------------------------
		cout << "Step " << i << ": Clusters undergoing birth/death/motion..." << endl;
		birthDeathMotionProcesses(clusterCenters, aliveClusters, birthProbability, deathProbability, motionStdDev);
		//------------------------------------------
		//generate the data for the current timestep
		//------------------------------------------
		cout << "Step " << i << ": Generating data from the clusters..." << endl;
		vector<V2d> clusterData;
		vector<int> trueLabels;
		generateData(clusterCenters, aliveClusters, nDataPerClusterPerStep, clusterStdDev, clusterData, trueLabels);

		//------------------------------------------
		//Take the vectors that we just created, and package them in the KD class defined above
		//------------------------------------------
		vgr.updateData(clusterData);

		//----------------------------
		//cluster using Dynamic Means
		//---------------------------
		vector<int> learnedLabels;
		vector<double> gammas;
		vector<int> prmlbls;
		double tTaken, obj;
		cout << "Step " << i << ": Clustering..." << endl;
		kdm.cluster(vgr, nRestarts, nCoarsest, learnedLabels, obj, gammas, prmlbls, tTaken);
		vgr.updateOldParameters(clusterData, learnedLabels, gammas, prmlbls);
		//----------------------------------------------------
		//calculate the accuracy via linear programming
		//including proper cluster label tracking (see above)
		//----------------------------------------------------
		matchings = getMaxMatchingConsistentWithOldMatching(learnedLabels, trueLabels, matchings);
		double acc = computeAccuracy(learnedLabels, trueLabels, matchings);
		cout << "Step " << i << ": Accuracy = " << acc <<  "\%" << endl;
		cumulativeAccuracy += acc;
	}
	cout << "Average Accuracy: " << cumulativeAccuracy/(double)nSteps << "\%" << endl;
	cout << "Done!" << endl;

	return 0;
}


//this function takes two label sets and a matching and outputs 
//the accuracy of the labelling (assuming labels1 = learned, labels2 = true)
double computeAccuracy(vector<int> labels1, vector<int> labels2, map<int, int> matchings){
	if (labels1.size() != labels2.size()){
		cout << "Error: computeAccuracy requires labels1/labels2 to have the same size" << endl;
		return -1;
	}

	double acc = 0;
	for (int i = 0; i < labels1.size(); i++){
		if (matchings[labels1[i]] == labels2[i]){
			acc += 1.0;
		}
	}
	//compute the accuracy
	return 100.0*acc/(double)labels1.size();
}


//this function takes a set of cluster centers and runs them through a birth/death/motion model
void birthDeathMotionProcesses(vector<V2d>& clusterCenters, vector<bool>& aliveClusters, double birthProbability, double deathProbability, double motionStdDev){

	//distributions
	uniform_real_distribution<double> uniformDistAng(0, 2*M_PI);
	uniform_real_distribution<double> uniformDist01(0, 1);
	normal_distribution<double> transitionDistRadial(0, motionStdDev);


	for (int j = 0; j < clusterCenters.size(); j++){
		//for each cluster center, decide whether it dies
		if (aliveClusters[j] && uniformDist01(rng) < deathProbability){
			cout << "Cluster " << j << " died." << endl;
			aliveClusters[j] = false;
		} else if (aliveClusters[j]) {
			//if it survived, move it stochastically
			//radius sampled from normal
			double steplen = transitionDistRadial(rng);
			//angle sampled from uniform
			double stepang = uniformDistAng(rng);
			clusterCenters[j](0) += steplen*cos(stepang);
			clusterCenters[j](1) += steplen*sin(stepang);
			cout << "Cluster " << j << " moved to " << clusterCenters[j].transpose() << endl;
		}
	}
	//decide whether to create a new cluster
	if (uniformDist01(rng) < birthProbability || all_of(aliveClusters.begin(), aliveClusters.end(), [](bool b){return !b;}) ) {
		V2d newCenter;
		newCenter(0) = uniformDist01(rng);
		newCenter(1) = uniformDist01(rng);
		clusterCenters.push_back(newCenter);
		aliveClusters.push_back(true);
		cout << "Cluster " << clusterCenters.size()-1 << " was created at " << clusterCenters.back().transpose() << endl;
	}
}

//this function takes a set of cluster centers and generates data from them
void generateData(vector<V2d> clusterCenters, vector<bool> aliveClusters, int nDataPerClusterPerStep, double clusterStdDev, vector<V2d>& clusterData, vector<int>& trueLabels){

	//distributions
	uniform_real_distribution<double> uniformDistAng(0, 2*M_PI);
	normal_distribution<double> likelihoodDistRadial(0, clusterStdDev);
	//loop through alive centers, generate nDataPerClusterPerStep datapoints for each
	for (int j = 0; j < clusterCenters.size(); j++){
		if (aliveClusters[j]){
			for (int k = 0; k < nDataPerClusterPerStep; k++){
				V2d newData = clusterCenters[j];
				double len = likelihoodDistRadial(rng);
				double ang = uniformDistAng(rng);
				newData(0) += len*cos(ang);
				newData(1) += len*sin(ang);
				clusterData.push_back(newData);
				trueLabels.push_back(j);
			}
		}
	}
}

