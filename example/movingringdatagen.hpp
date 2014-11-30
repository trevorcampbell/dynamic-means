#ifndef __MOVINGRINGDATAGEN_HPP
#include<Eigen/Dense>
#include<dmeans/utils>
#include<vector>

class MovingRingDataGenerator{
	typedef Eigen::Vector2d V2d;
	private:
	std::vector<V2d> centers;
	std::vector<bool> alive;
	double birthProbability, deathProbability, motionStdDev, clusterStdDev, radius;
	int initialClusters, nDataPerClusterPerStep;
	bool verbose;

	public:
		MovingRingDataGenerator(dmeans::Config cfg){
			radius = cfg.get("radius", dmeans::Config::OPTIONAL, 0.1);
			birthProbability = cfg.get("birthProbability", dmeans::Config::OPTIONAL, 0.10);
			deathProbability = cfg.get("deathProbability", dmeans::Config::OPTIONAL, 0.05);
			motionStdDev = cfg.get("motionStdDev", dmeans::Config::OPTIONAL, 0.05);
			clusterStdDev = cfg.get("clusterStdDev", dmeans::Config::OPTIONAL, 0.05);
			initialClusters = cfg.get("initialClusters", dmeans::Config::OPTIONAL, 4);
			nDataPerClusterPerStep = cfg.get("nDataPerClusterPerStep", dmeans::Config::OPTIONAL, 15);
			verbose = cfg.get("verbose", dmeans::Config::OPTIONAL, false);

			std::uniform_real_distribution<double> uniformDist01(0, 1);
			for (int i = 0; i < initialClusters; i++){
					V2d newCenter;
					newCenter(0) = uniformDist01(dmeans::RNG::get());
					newCenter(1) = uniformDist01(dmeans::RNG::get());
					centers.push_back(newCenter);
					alive.push_back(true);
			}
		}

		void step(){
			//distributions
			std::uniform_real_distribution<double> uniformDistAng(0, 2*M_PI);
			std::uniform_real_distribution<double> uniformDist01(0, 1);
			std::normal_distribution<double> transitionDistRadial(0, motionStdDev);


			for (uint64_t j = 0; j < centers.size(); j++){
				//for each cluster center, decide whether it dies
				if (alive[j] && uniformDist01(dmeans::RNG::get()) < deathProbability){
					if(verbose){cout << "Cluster " << j << " died." << endl;}
					alive[j] = false;
				} else if (alive[j]) {
					//if it survived, move it stochastically
					//radius sampled from normal
					double steplen = transitionDistRadial(dmeans::RNG::get());
					//angle sampled from uniform
					double stepang = uniformDistAng(dmeans::RNG::get());
					centers[j](0) += steplen*cos(stepang);
					centers[j](1) += steplen*sin(stepang);
					cout << "Cluster " << j << " moved to " << centers[j].transpose() << endl;
				}
			}
			//decide whether to create a new cluster
			if (uniformDist01(dmeans::RNG::get()) < birthProbability || all_of(alive.begin(), alive.end(), [](bool b){return !b;}) ) {
				V2d newCenter;
				newCenter(0) = uniformDist01(dmeans::RNG::get());
				newCenter(1) = uniformDist01(dmeans::RNG::get());
				centers.push_back(newCenter);
				alive.push_back(true);
				cout << "Cluster " << centers.size()-1 << " was created at " << centers.back().transpose() << endl;
			}
		}
		void get(std::vector<V2d >& vdata, std::vector<uint64_t>& trueLabels) const{
			vdata.clear(); trueLabels.clear();
			//distributions
			uniform_real_distribution<double> uniformDistAng(0, 2*M_PI);
			normal_distribution<double> likelihoodDistRadial(0, clusterStdDev);
			//loop through alive centers, generate nDataPerClusterPerStep datapoints for each
			for (uint64_t j = 0; j < centers.size(); j++){
				if (alive[j]){
					for (int k = 0; k < nDataPerClusterPerStep; k++){
						V2d newData = centers[j];
						double len = radius+likelihoodDistRadial(dmeans::RNG::get());
						double ang = uniformDistAng(dmeans::RNG::get());
						newData(0) += len*cos(ang);
						newData(1) += len*sin(ang);
						vdata.push_back(newData);
						trueLabels.push_back(j);
					}
				}
			}
		}
};
#define __MOVINGRINGDATAGEN_HPP
#endif /* __MOVINGRINGDATAGEN_HPP */
