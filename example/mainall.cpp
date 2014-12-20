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
#include <dmeans/iterative>
#include <dmeans/spectral>
#include <dmeans/kernelized>
#include <dmeans/model>
#include <dmeans/utils>

#include "movinggaussdatagen.hpp"
#include "movingringdatagen.hpp"
#include "movingshapedatagen.hpp"

using namespace std;

typedef Eigen::Vector2d V2d;
typedef dmeans::VectorSpaceModel<2> VSModel;
typedef dmeans::ExponentialKernelModel<2> EKModel;

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
	data_cfg.set("clusterStdDev", 0.0);
	data_cfg.set("nDataPerClusterPerStep", 15);
	data_cfg.set("initialClusters", 4);
	//MovingDataGenerator* datagen = new MovingGaussianDataGenerator(data_cfg);
	//MovingDataGenerator* datagen = new MovingRingDataGenerator(data_cfg);
	data_cfg.set("radius", 0.05);
	MovingDataGenerator* datagen = new MovingShapeDataGenerator(data_cfg);

	//the Dynamic Means object
	//play with lambda/Q/tau to change Dynamic Means' performance
	dmeans::Config dynm_cfg;
	double T_Q = 6.8;
	double K_tau = 1.05;
	double lambda = 0.05;
	double Q = lambda/T_Q;
	double tau = (T_Q*(K_tau-1.0)+1.0)/(T_Q-1.0);
	dynm_cfg.set("lambda", lambda);
	dynm_cfg.set("Q", Q);
	dynm_cfg.set("tau", tau);
	dynm_cfg.set("nRestarts", 10);
	dynm_cfg.set("verbose", true);
	dmeans::DMeans<VSModel, dmeans::IterativeWithMonotonicityChecks> dynm(dynm_cfg);

	dmeans::Config kdynm_cfg;
	lambda = 10;
	T_Q = 5;
	K_tau = 1.05;
	Q = lambda/T_Q;
	tau = (T_Q*(K_tau-1.0)+1.0)/(T_Q-1.0);
	kdynm_cfg.set("lambda", lambda);
	kdynm_cfg.set("Q", Q);
	kdynm_cfg.set("tau", tau);
	kdynm_cfg.set("nRestarts", 10);
	kdynm_cfg.set("kernelWidth", 0.07);
	kdynm_cfg.set("sparseApproximationSize", 15);
	kdynm_cfg.set("verbose", true);
	dmeans::DMeans<EKModel, dmeans::KernelizedWithMonotonicityChecks> kdynm(kdynm_cfg);

	dmeans::Config sdynm_cfg;
	lambda = 10;
	T_Q = 5;
	K_tau = 1.05;
	Q = lambda/T_Q;
	tau = (T_Q*(K_tau-1.0)+1.0)/(T_Q-1.0);
	sdynm_cfg.set("lambda", lambda);
	sdynm_cfg.set("Q", Q);
	sdynm_cfg.set("tau", tau);
	sdynm_cfg.set("nRestarts", 1);
	sdynm_cfg.set("nProjectionRestarts", 10);
	sdynm_cfg.set("verbose", true);
	sdynm_cfg.set("kernelWidth", 0.07);
	sdynm_cfg.set("sparseApproximationSize", 15);
	sdynm_cfg.set("eigenSolverType", dmeans::EigenSolver::Type::REDSVD);
	sdynm_cfg.set("eigenSolverDimension", 40);
	dmeans::DMeans<EKModel, dmeans::SpectralWithMonotonicityChecks> sdynm(sdynm_cfg);

	//run the experiment
	double cumulativeAccuracyD = 0.0;//stores the accuracy accumulated for each step
	double cumulativeAccuracyS = 0.0;//stores the accuracy accumulated for each step
	double cumulativeAccuracyK = 0.0;//stores the accuracy accumulated for each step
	double cumulativeTimeD = 0.0; //store the cpu time accumulated at each step
	double cumulativeTimeS = 0.0; //store the cpu time accumulated at each step
	double cumulativeTimeK = 0.0; //store the cpu time accumulated at each step
	map<int, int> matchingsD, matchingsS, matchingsK;//stores the saved matchings from previous timesteps
							//enables proper label tracking (see note at line 27)
	dmeans::MaxMatching maxmD, maxmS, maxmK;
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
		vector<VSModel::Data> dataVS;
		vector<EKModel::Data> dataEK;
		datagen->get(vdata, trueLabels);
		ofstream dataout("data.log", ios_base::app);
		dataout << vdata.size() << endl;
		for(uint64_t i = 0; i < vdata.size(); i++){
			VSModel::Data d;
			EKModel::Data f;
			d.v = vdata[i];
			f.v = vdata[i];
			dataVS.push_back(d);
			dataEK.push_back(f);
			dataout << vdata[i].transpose() << endl;
		}
		dataout.close();

		//***************************
		//cluster using Dynamic Means
		//***************************
		cout << "Step " << i << ": Clustering " << dataVS.size() << " datapoints..." << endl;
		dmeans::Results<VSModel> resD = dynm.cluster(dataVS);
		dmeans::Results<EKModel> resS = sdynm.cluster(dataEK);
		dmeans::Results<EKModel> resK = kdynm.cluster(dataEK);

		//***************************************************
		//calculate the accuracy via linear programming
		//including proper cluster label tracking (see above)
		//***************************************************
		vector<int> rlblintD, tlblintD; //just a not-so-clever way to convert vector<uint64_t> to vector<int>
		for (uint64_t j = 0; j < trueLabels.size(); j++){
			rlblintD.push_back(resD.lbls[j]);
			tlblintD.push_back(trueLabels[j]);
		}
		matchingsD = maxmD.getMaxConsistentMatching(rlblintD, tlblintD, std::vector<double>());
		double accD = 100.0*(double)maxmD.getObjective()/ (double)dataVS.size();

		vector<int> rlblintS, tlblintS; //just a not-so-clever way to convert vector<uint64_t> to vector<int>
		for (uint64_t j = 0; j < trueLabels.size(); j++){
			rlblintS.push_back(resS.lbls[j]);
			tlblintS.push_back(trueLabels[j]);
		}
		matchingsS = maxmS.getMaxConsistentMatching(rlblintS, tlblintS, std::vector<double>());
		double accS = 100.0*(double)maxmS.getObjective()/ (double)dataEK.size();

		vector<int> rlblintK, tlblintK; //just a not-so-clever way to convert vector<uint64_t> to vector<int>
		for (uint64_t j = 0; j < trueLabels.size(); j++){
			rlblintK.push_back(resK.lbls[j]);
			tlblintK.push_back(trueLabels[j]);
		}
		matchingsK = maxmK.getMaxConsistentMatching(rlblintK, tlblintK, std::vector<double>());
		double accK = 100.0*(double)maxmK.getObjective()/ (double)dataEK.size();

		//double acc = computeAccuracy(learnedLabels, trueLabels, matchings);
		cout << "Step " << i << ": Accuracy(D-Means) = " << accD <<  "\%" << " CPU Time = " << resD.tTaken << "s" << endl;
		cout << "Step " << i << ": Accuracy(SD-Means) = " << accS <<  "\%" << " CPU Time = " << resS.tTaken << "s" << endl;
		cout << "Step " << i << ": Accuracy(KD-Means) = " << accK <<  "\%" << " CPU Time = " << resK.tTaken << "s" << endl;
		cumulativeAccuracyD += accD;
		cumulativeAccuracyS += accS;
		cumulativeAccuracyK += accK;
		cumulativeTimeD += resD.tTaken;
		cumulativeTimeS += resS.tTaken;
		cumulativeTimeK += resK.tTaken;
	}
	cout << "Average Accuracy(D-Means): " << cumulativeAccuracyD/(double)nSteps << "\% Total CPU Time = " << cumulativeTimeD << "s" << endl;
	cout << "Average Accuracy(SD-Means): " << cumulativeAccuracyS/(double)nSteps << "\% Total CPU Time = " << cumulativeTimeS << "s" << endl;
	cout << "Average Accuracy(KD-Means): " << cumulativeAccuracyK/(double)nSteps << "\% Total CPU Time = " << cumulativeTimeK << "s" << endl;
	cout << "Done!" << endl;
	delete datagen;

	return 0;
}

