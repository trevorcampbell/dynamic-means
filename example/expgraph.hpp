#ifndef __EXPGRAPH_HPP
#include <Eigen/Dense>
#include<vector>
using namespace std;
class ExpGraph{
	public:
		//None of data, oldprms, oldprmlbls, updateData(), and updateOldParameters are mandatory 
		//I implemented them for my own application
		//data holds the current timestep's vectors, oldprms/oldprmlbls keeps track of old parameters and their labels
		//updateData() does what it advertises
		//updateOldParameters makes sure that the KernDynMeans and VectorGraph objects have the same ordering of oldprmlbls,
		//and also updates oldprms (weighted by gammas) to reflect the clustering that was just completed
		//
		//basically, the pattern for clustering batch-sequential data is:
		//1) VectorGraph::updateData( [the data for this timestep] );
		//2) KernDynMeans::cluster( [the vector graph] );
		//3) VectorGraph::updateOldPrms( [output of KernDynMeans::cluster] );
		//Go back to 1 for the next time step
		std::vector<V2d> data, oldprms;
		std::vector<int> oldprmlbls;
		ExpGraph(){
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

		//---ALL FUNCTIONS BELOW ARE MANDATORY---
		//Any affinity class must implement all of the below functions, as KernDynMeans calls them explicitly

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

#define __EXPGRAPH_HPP
#endif /* __EXPGRAPH_HPP */
