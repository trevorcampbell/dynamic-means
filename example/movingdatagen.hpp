#ifndef __MOVINGDATAGEN_HPP
class MovingDataGenerator{
	typedef Eigen::Vector2d V2d;
	protected:
	std::vector<bool> alive;
	double birthProbability, deathProbability;
	int initialClusters, nDataPerClusterPerStep;
	bool verbose;

	public:
		MovingDataGenerator(dmeans::Config cfg){
			birthProbability = cfg.get("birthProbability", dmeans::Config::OPTIONAL, 0.10);
			deathProbability = cfg.get("deathProbability", dmeans::Config::OPTIONAL, 0.05);
			initialClusters = cfg.get("initialClusters", dmeans::Config::OPTIONAL, 4);
			nDataPerClusterPerStep = cfg.get("nDataPerClusterPerStep", dmeans::Config::OPTIONAL, 15);
			verbose = cfg.get("verbose", dmeans::Config::OPTIONAL, false);

			for (int i = 0; i < initialClusters; i++) alive.push_back(true);
		}

		void step(){
			std::uniform_real_distribution<double> uniformDist01(0, 1);
			for (uint64_t j = 0; j < alive.size(); j++){
				//for each cluster center, decide whether it dies
				if (alive[j] && uniformDist01(dmeans::RNG::get()) < deathProbability){
					if(verbose){cout << "Cluster " << j << " died." << endl;}
					alive[j] = false;
				} else if (alive[j]) {
					this->transitionCluster(j);
				}
			}
			//decide whether to create a new cluster
			if (uniformDist01(dmeans::RNG::get()) < birthProbability || all_of(alive.begin(), alive.end(), [](bool b){return !b;}) ) {
				this->createCluster();
				alive.push_back(true);
			}
		}
		void get(std::vector<V2d>& vdata, std::vector<uint64_t>& trueLabels) const{
			vdata.clear(); trueLabels.clear();
			//distributions
			//loop through alive centers, generate nDataPerClusterPerStep datapoints for each
			for (uint64_t j = 0; j < alive.size(); j++){
				if (alive[j]){
					std::vector<V2d> ndata = this->genData(j, nDataPerClusterPerStep);
					vdata.insert(vdata.end(), ndata.begin(), ndata.end());
					trueLabels.insert(trueLabels.end(), ndata.size(), j);
				}
			}
		}

		virtual void transitionCluster(int j) = 0;
		virtual void createCluster() = 0;
		virtual std::vector<V2d> genData(int j, int n) const = 0;
};
#define __MOVINGDATAGEN_HPP
#endif /* __MOVINGDATAGEN_HPP */
