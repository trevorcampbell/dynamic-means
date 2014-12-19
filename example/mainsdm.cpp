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
#include <dmeans/spectral>
#include <dmeans/model>
#include <dmeans/utils>

#include "movinggaussdatagen.hpp"
#include "movingringdatagen.hpp"
using namespace std;

typedef Eigen::Vector2d V2d;
typedef dmeans::ExponentialKernelModel<2> ESModel;

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
	int nSteps = 100;//run the experiment for nSteps steps
	//play with the below constants to change the data generation process
	dmeans::Config data_cfg;
	data_cfg.set("birthProbability", 0.10);
	data_cfg.set("deathProbability", 0.05);
	data_cfg.set("motionStdDev", 0.05);
	data_cfg.set("clusterStdDev", 0.05);
	data_cfg.set("nDataPerClusterPerStep", 15);
	data_cfg.set("initialClusters", 4);
	MovingDataGenerator* datagen = new MovingGaussianDataGenerator(data_cfg);
	//MovingDataGenerator* datagen = new MovingRingGenerator(data_cfg);

	//the Dynamic Means object
	//play with lambda/Q/tau to change Dynamic Means' performance
	dmeans::Config dynm_cfg;
	double kernelWidth = 0.07;
	double lambda = 10;
	double T_Q = 5;
	double K_tau = 1.05;
	double Q = lambda/T_Q;
	double tau = (T_Q*(K_tau-1.0)+1.0)/(T_Q-1.0);
	dynm_cfg.set("lambda", lambda);
	dynm_cfg.set("Q", Q);
	dynm_cfg.set("tau", tau);
	dynm_cfg.set("nRestarts", 1);
	dynm_cfg.set("nProjectionRestarts", 10);
	dynm_cfg.set("verbose", true);
	dynm_cfg.set("kernelWidth", kernelWidth);
	dynm_cfg.set("sparseApproximationSize", 15);
	dynm_cfg.set("eigenSolverType", dmeans::EigenSolver::Type::REDSVD);
	dynm_cfg.set("eigenSolverDimension", 40);
	dmeans::DMeans<ESModel, dmeans::SpectralWithMonotonicityChecks> dynm(dynm_cfg);

	//run the experiment
	double cumulativeAccuracy = 0.0;//stores the accuracy accumulated for each step
	double cumulativeTime = 0.0; //store the cpu time accumulated at each step
	map<int, int> matchings;//stores the saved matchings from previous timesteps
							//enables proper label tracking (see note at line 27)
	dmeans::MaxMatching maxm;
	for (int i = 0; i < nSteps; i++){
		//****************************
		//birth/death/motion processes
		//****************************
		cout << "Step " << i << ": Clusters undergoing birth/death/motion..." << endl;
		datagen->step();
		//******************************************
		//generate the data for the current timestep
		//******************************************
		cout << "Step " << i << ": Generating data from the clusters..." << endl;
		vector<V2d> vdata;
		vector<uint64_t> trueLabels;
		vector<ESModel::Data> data;
		datagen->get(vdata, trueLabels);
		for(uint64_t i = 0; i < vdata.size(); i++){
			ESModel::Data d;
			d.v = vdata[i];
			data.push_back(d);
		}

		//***************************
		//cluster using Dynamic Means
		//***************************
		cout << "Step " << i << ": Clustering " << data.size() << " datapoints..." << endl;
		dmeans::Results<ESModel> res = dynm.cluster(data);

		//***************************************************
		//calculate the accuracy via linear programming
		//including proper cluster label tracking (see above)
		//***************************************************
		vector<int> rlblint, tlblint; //just a not-so-clever way to convert vector<uint64_t> to vector<int>
		for (uint64_t j = 0; j < trueLabels.size(); j++){
			rlblint.push_back(res.lbls[j]);
			tlblint.push_back(trueLabels[j]);
		}
		matchings = maxm.getMaxConsistentMatching(rlblint, tlblint, std::vector<double>());
		double acc = 100.0*(double)maxm.getObjective()/ (double)data.size();
		//double acc = computeAccuracy(learnedLabels, trueLabels, matchings);
		cout << "Step " << i << ": Accuracy = " << acc <<  "\%" << " CPU Time = " << res.tTaken << "s" << endl;
		cumulativeAccuracy += acc;
		cumulativeTime += res.tTaken;
	}
	cout << "Average Accuracy: " << cumulativeAccuracy/(double)nSteps << "\% Total CPU Time = " << cumulativeTime << "s" << endl;
	cout << "Done!" << endl;
	delete datagen;

	return 0;
}








