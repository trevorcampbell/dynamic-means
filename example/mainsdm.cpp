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
#include <omp.h>

#include <dmeans/core>
#include <dmeans/spectral>
#include <dmeans/model>
#include <dmeans/utils>

#include "movinggaussdatagen.hpp"
#include "movingringdatagen.hpp"
using namespace std;

typedef Eigen::Vector2d V2d;
typedef dmeans::MSTKernelModel<2> MSTModel;

int main(int argc, char** argv){
	//paramsweep();
	//generates clusters that jump around on the domain R^2
	//they move stochastically with a normal distribution w/ std deviation 0.05
	//they die stochastically with probability 0.05 at each step
	//a cluster is created stochastically with probability 0.10 at each step in the area [0,1]x[0,1]
	
	//note: when computing accuracy, ***proper tracking is taken into account***
	//i.e. if true label 1 is matched to learned label 3, then that matching is fixed from that point forward
	//later pairings that attempt to match true label 1 to something else
	//or learned label 3 to something else are discarded
	

	//constants
	int nSteps = 10;//run the experiment for nSteps steps
	//play with the below constants to change the data generation process
	dmeans::Config data_cfg;
	data_cfg.set("birthProbability", 0.0);
	data_cfg.set("deathProbability", 0.0);
	data_cfg.set("radius", 0.4);
	data_cfg.set("motionStdDev", 0.05);
	data_cfg.set("clusterStdDev", 0.01);
	data_cfg.set("nDataPerClusterPerStep", 200);
	data_cfg.set("initialClusters", 3);
	MovingDataGenerator* datagen = new MovingRingDataGenerator(data_cfg);

	//the Dynamic Means object
	//play with lambda/Q/tau to change Dynamic Means' performance
	dmeans::Config dynm_cfg;
	double lambda = 55;
	double T_Q = 13;
	double K_tau = 4.5;
	double Q = lambda/T_Q;
	double tau = (T_Q*(K_tau-1.0)+1.0)/(T_Q-1.0);
	double jumpThresh = 0.07;
	//double omega = 0.07;
	dynm_cfg.set("lambda", lambda);
	dynm_cfg.set("Q", Q);
	dynm_cfg.set("tau", tau);
	dynm_cfg.set("nRestarts", 1);
	dynm_cfg.set("nProjectionRestarts", 10);
	dynm_cfg.set("verbose", true);
	dynm_cfg.set("jumpThreshold", jumpThresh);
//	dynm_cfg.set("kernelWidth", omega);
	dynm_cfg.set("sparseApproximationSize", 20);
	dynm_cfg.set("eigenSolverType", dmeans::EigenSolver::Type::EIGEN_SELF_ADJOINT);
	dynm_cfg.set("eigenSolverDimension", 40);
	dmeans::DMeans<MSTModel, dmeans::MatchingSpectralWithMonotonicityChecks> dynm(dynm_cfg);

	//run the experiment
	double cumulativeAccuracy = 0.0;//stores the accuracy accumulated for each step
	double cumulativeTime = 0.0; //store the cpu time accumulated at each step
	map<int, int> matchings;//stores the saved matchings from previous timesteps
							//enables proper label tracking (see note at line 27)
	dmeans::MaxMatching maxm;
	ofstream dataout("data.log", ios_base::trunc);
	ofstream lblout("lbls.log", ios_base::trunc);
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
		vector<MSTModel::Data> data;
		datagen->get(vdata, trueLabels);
		dataout << vdata.size() << endl;
		for(uint64_t i = 0; i < vdata.size(); i++){
			MSTModel::Data d;
			d.v = vdata[i];
			data.push_back(d);
			dataout << vdata[i].transpose() << endl;
		}


		//***************************
		//cluster using Dynamic Means
		//***************************
		cout << "Step " << i << ": Clustering " << data.size() << " datapoints..." << endl;
		dmeans::Results<MSTModel> res = dynm.cluster(data);
		lblout << vdata.size() << endl;
		for (uint64_t i = 0; i < vdata.size(); i++){
			lblout << res.lbls[i] << endl;
		}

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

		//MSTModel::MST treetest;
		//treetest.construct(vdata);
		//treetest.write("tree.log", jumpThresh);
		//return 0;

	}
	cout << "Average Accuracy: " << cumulativeAccuracy/(double)nSteps << "\% Total CPU Time = " << cumulativeTime << "s" << endl;
	cout << "Done!" << endl;
	dataout.close();
	lblout.close();
	delete datagen;

	return 0;
}

void paramsweep(){
	std::vector<double> lambdas, T_Qs, K_taus;
	lambdas.push_back(15.0);
	lambdas.push_back(20.0);
	lambdas.push_back(25.0);
	lambdas.push_back(30.0);
	lambdas.push_back(35.0);
	lambdas.push_back(40.0);
	lambdas.push_back(45.0);
	lambdas.push_back(50.0);
	lambdas.push_back(55.0);

	T_Qs.push_back(3.0);
	T_Qs.push_back(5.0);
	T_Qs.push_back(7.0);
	T_Qs.push_back(9.0);
	T_Qs.push_back(11.0);
	T_Qs.push_back(13.0);
	T_Qs.push_back(15.0);
	T_Qs.push_back(17.0);

	K_taus.push_back(1.0);
	K_taus.push_back(1.5);
	K_taus.push_back(2.0);
	K_taus.push_back(2.5);
	K_taus.push_back(3.0);
	K_taus.push_back(3.5);
	K_taus.push_back(4.0);
	K_taus.push_back(4.5);


	std::ofstream paramsweepout("paramsweep.log", std::ios_base::app);
	#pragma omp parallel for schedule(dynamic, 1) collapse(4)
	for (int iii = 0; iii < lambdas.size(); iii++){
	for (int jjj = 0; jjj < T_Qs.size(); jjj++){
	for (int kkk = 0; kkk < K_taus.size(); kkk++){
	for (int nnn = 0; nnn < 20; nnn++){
		//constants
		int nSteps = 10;//run the experiment for nSteps steps
		//play with the below constants to change the data generation process
		dmeans::Config data_cfg;
		data_cfg.set("birthProbability", 0.0);
		data_cfg.set("deathProbability", 0.0);
		data_cfg.set("radius", 0.4);
		data_cfg.set("motionStdDev", 0.02);
		data_cfg.set("clusterStdDev", 0.01);
		data_cfg.set("nDataPerClusterPerStep", 200);
		data_cfg.set("initialClusters", 3);
		MovingDataGenerator* datagen = new MovingRingDataGenerator(data_cfg);

		//the Dynamic Means object
		//play with lambda/Q/tau to change Dynamic Means' performance
		dmeans::Config dynm_cfg;
		double lambda = lambdas[iii];
		double T_Q = T_Qs[jjj];
		double K_tau = K_taus[kkk];
		double Q = lambda/T_Q;
		double tau = (T_Q*(K_tau-1.0)+1.0)/(T_Q-1.0);
		double jumpThresh = 0.07;
		//double omega = 0.07;
		dynm_cfg.set("lambda", lambda);
		dynm_cfg.set("Q", Q);
		dynm_cfg.set("tau", tau);
		dynm_cfg.set("nRestarts", 1);
		dynm_cfg.set("nProjectionRestarts", 10);
		dynm_cfg.set("verbose", true);
		dynm_cfg.set("jumpThreshold", jumpThresh);
//		dynm_cfg.set("kernelWidth", omega);
		dynm_cfg.set("sparseApproximationSize", 100);
		dynm_cfg.set("eigenSolverType", dmeans::EigenSolver::Type::EIGEN_SELF_ADJOINT);
		dynm_cfg.set("eigenSolverDimension", 100);
		dmeans::DMeans<MSTModel, dmeans::MatchingSpectralWithMonotonicityChecks> dynm(dynm_cfg);

		//run the experiment
		double cumulativeAccuracy = 0.0;//stores the accuracy accumulated for each step
		double cumulativeTime = 0.0; //store the cpu time accumulated at each step
		map<int, int> matchings;//stores the saved matchings from previous timesteps
								//enables proper label tracking (see note at line 27)
		dmeans::MaxMatching maxm;
		ofstream dataout("data.log", ios_base::trunc);
		ofstream lblout("lbls.log", ios_base::trunc);
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
			vector<MSTModel::Data> data;
			datagen->get(vdata, trueLabels);
			dataout << vdata.size() << endl;
			for(uint64_t i = 0; i < vdata.size(); i++){
				MSTModel::Data d;
				d.v = vdata[i];
				data.push_back(d);
				dataout << vdata[i].transpose() << endl;
			}


			//***************************
			//cluster using Dynamic Means
			//***************************
			cout << "Step " << i << ": Clustering " << data.size() << " datapoints..." << endl;
			dmeans::Results<MSTModel> res = dynm.cluster(data);
			lblout << vdata.size() << endl;
			for (uint64_t i = 0; i < vdata.size(); i++){
				lblout << res.lbls[i] << endl;
			}

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

			//MSTModel::MST treetest;
			//treetest.construct(vdata);
			//treetest.write("tree.log", jumpThresh);
			//return 0;

		}
		cout << "Average Accuracy: " << cumulativeAccuracy/(double)nSteps << "\% Total CPU Time = " << cumulativeTime << "s" << endl;
		cout << "Done!" << endl;
		#pragma omp critical
		{ //first three outputs are a unique ID for this parameter combination
		paramsweepout << iii << jjj << kkk << " " << lambdas[iii] << " " << T_Qs[jjj] << " " << K_taus[kkk] << " " << cumulativeAccuracy/(double)nSteps << " " << cumulativeTime/(double)nSteps << endl;
		}

		dataout.close();
		lblout.close();
		delete datagen;
	}
	}
	}
	}


	paramsweepout.close();
	return;
}







