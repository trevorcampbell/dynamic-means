#include <fstream>
#include<algorithm>
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

#include <dmeans/core>
#include <dmeans/matching>

using namespace std;

typedef Eigen::Vector2d V2d;
typedef dmeans::VectorSpaceModel<2> VSModel;

//function declarations
double computeAccuracy(vector<int> labels1, vector<int> labels2, map<int, int> matchings);
void birthDeathMotionProcesses(vector<V2d>& clusterCenters, vector<bool>& aliveClusters, double birthProbability, double deathProbability, double motionStdDev);
void generateData(vector<V2d> clusterCenters, vector<bool> aliveClusters, int nDataPerClusterPerStep, double likelihoodstd, vector<VSModel::Data>& data, vector<uint64_t>& trueLabels);

//random number generator
mt19937 rng;//uses the same seed every time, 5489u

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
	int nDataPerClusterPerStep = 15; // datapoints generated for each cluster per timestep
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
	Config cfg;
	const double T_Q = 6.8;
	const double K_tau = 1.01;
	const double lambda = 0.05;
	const double Q = lambda/T_Q;
	const double tau = (T_Q*(K_tau-1.0)+1.0)/(T_Q-1.0);

	cfg.set("lambda", lambda);
	cfg.set("Q", Q);
	cfg.set("tau", tau);
	cfg.set("nRestarts", 10);
	cfg.set("verbose", true);
	dmeans::DMeans<VSModel, dmeans::IterativeWithMonotonicityChecks> dynm(cfg);

	//run the experiment
	double cumulativeAccuracy = 0;//stores the accuracy accumulated for each step
	map<int, int> matchings;//stores the saved matchings from previous timesteps
							//enables proper label tracking (see note at line 27)
	dmeans::MaxMatching maxm;
	for (int i = 0; i < nSteps; i++){
		//****************************
		//birth/death/motion processes
		//****************************
		cout << "Step " << i << ": Clusters undergoing birth/death/motion..." << endl;
		birthDeathMotionProcesses(clusterCenters, aliveClusters, birthProbability, deathProbability, motionStdDev);
		//******************************************
		//generate the data for the current timestep
		//******************************************
		cout << "Step " << i << ": Generating data from the clusters..." << endl;
		vector<VSModel::Data> data;
		vector<uint64_t> trueLabels;
		generateData(clusterCenters, aliveClusters, nDataPerClusterPerStep, clusterStdDev, data, trueLabels);

		//***************************
		//cluster using Dynamic Means
		//***************************
		cout << "Step " << i << ": Clustering " << data.size() << " datapoints..." << endl;
		dmeans::Results<VSModel> res = dynm.cluster(data);

		//***************************************************
		//calculate the accuracy via linear programming
		//including proper cluster label tracking (see above)
		//***************************************************
		vector<int> rlblint, tlblint;
		for (uint64_t j = 0; j < trueLabels.size(); j++){
			rlblint.push_back(res.lbls[j]);
			tlblint.push_back(trueLabels[j]);
		}
		//std::transform(res.lbls.begin(), res.lbls.end(), rlblint.begin(), [](uint64_t u) -> int { int ui = u; return ui;});
		//std::transform(trueLabels.begin(), trueLabels.end(), tlblint.begin(), [](uint64_t u) -> int { int ui = u; return ui;});
		matchings = maxm.getMaxConsistentMatching(rlblint, tlblint, std::vector<double>());
		double acc = 100.0*(double)maxm.getObjective()/ (double)data.size();
		//double acc = computeAccuracy(learnedLabels, trueLabels, matchings);
		cout << "Step " << i << ": Accuracy = " << acc <<  "\%" << endl;
		cumulativeAccuracy += acc;
	}
	cout << "Average Accuracy: " << cumulativeAccuracy/(double)nSteps << "\%" << endl;
	cout << "Done!" << endl;

	return 0;
}


//this function takes a set of cluster centers and runs them through a birth/death/motion model
void birthDeathMotionProcesses(vector<V2d>& clusterCenters, vector<bool>& aliveClusters, double birthProbability, double deathProbability, double motionStdDev){

	//distributions
	uniform_real_distribution<double> uniformDistAng(0, 2*M_PI);
	uniform_real_distribution<double> uniformDist01(0, 1);
	normal_distribution<double> transitionDistRadial(0, motionStdDev);


	for (uint64_t j = 0; j < clusterCenters.size(); j++){
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
void generateData(vector<V2d> clusterCenters, vector<bool> aliveClusters, int nDataPerClusterPerStep, double clusterStdDev, vector<VSModel::Data>& data, vector<uint64_t>& trueLabels){

	//distributions
	uniform_real_distribution<double> uniformDistAng(0, 2*M_PI);
	normal_distribution<double> likelihoodDistRadial(0, clusterStdDev);
	//loop through alive centers, generate nDataPerClusterPerStep datapoints for each
	for (uint64_t j = 0; j < clusterCenters.size(); j++){
		if (aliveClusters[j]){
			for (int k = 0; k < nDataPerClusterPerStep; k++){
				V2d newData = clusterCenters[j];
				double len = likelihoodDistRadial(rng);
				double ang = uniformDistAng(rng);
				newData(0) += len*cos(ang);
				newData(1) += len*sin(ang);
				VSModel::Data d;
				d.v = newData;
				data.push_back(d);
				trueLabels.push_back(j);
			}
		}
	}
}







