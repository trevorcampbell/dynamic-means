#ifndef __MATCHINGSPECTRAL_IMPL_HPP
template<class Model, bool monoCheck>
_MatchingSpectral<Model, monoCheck>::_MatchingSpectral(const Config& cfg){
	this->solverType = static_cast<EigenSolver::Type>(this->cfg.get("eigenSolverType", Config::Type::OPTIONAL, static_cast<int>(EigenSolver::Type::EIGEN_SELF_ADJOINT)));
	if (this->solverType == EigenSolver::Type::EIGEN_SELF_ADJOINT){
		this->nEigs = 0;
		std::cout << "USING EIGEN SELF ADJOINT" << std::endl;
	} else {
		this->nEigs = this->cfg.get("eigenSolverDimension", Config::Type::REQUIRED, 0);
		std::cout << "USING REDSVD" << std::endl;
	}
	this->verbose = cfg.get("verbose", Config::Type::OPTIONAL, false);
	this->nProjectionRestarts = this->cfg.get("nProjectionRestarts", Config::Type::OPTIONAL, 1);
}

template<class Model, bool monoCheck>
double _MatchingSpectral<Model, monoCheck>::cluster(const std::vector<typename Model::Data>& obs, std::vector<Clus>& clus, const Model& model) const{
	if(obs.size() == 0){return 0.0;}

	//compute the kernel matrix
	MXd KUp = this->getKernelMatUpper(obs, clus, model);
	//solve the eigensystem
	EigenSolver eigsol(KUp, this->solverType, this->nEigs, model.getEigenvalueLowerThreshold());

	MXd Z; VXd eigvals;
	eigsol.getResults(eigvals, Z);

	if (eigvals.size() == 0){
		//all eigenvalues/vectors were thresholded out -- just use a constant single column
		Z = MXd::Ones(KUp.rows(), 1);
		eigvals = VXd::Ones(1);
	}

	uint64_t nA = obs.size();

	const uint64_t nZCols = Z.cols(); //number of clusters currently instantiated

	//conservative resize Z
	Z.conservativeResize(nA, Z.cols());

	//normalize the new rows of Z
	for (uint64_t j = 0; j < nA; j++){
		double rownorm = sqrt(Z.row(j).squaredNorm());
		if (rownorm == 0){ //if rownorm is a hard zero, just set the row to ones -- it's the only thing we can do, lambda was set too high
			Z.row(j) = MXd::Ones(1, nZCols);
			rownorm = sqrt(Z.row(j).squaredNorm());
		}
		Z.row(j) *= (1.0/rownorm);
	}

	//initialize X (constrained version of Z) 
	MXd X(nA, nZCols);
	//initialize V (rotation matrix on Z to make Z*V close to X)
	MXd V(nZCols, nZCols);

	//stuff for storing best clustering (best clustering is determined by normalized cuts objective)
	double minObj = numeric_limits<double>::infinity();
	vector<uint64_t> minLbls;

	//propose nProjectionRestarts V trials
	for (uint64_t i = 0; i < nProjectionRestarts; i++){
		V.setZero();
		//initialize unitary V via ``most orthogonal rows'' method
		int rndRow = (RNG::get())()%nA;
		V.col(0) = Z.row(rndRow).transpose();
		MXd c(nA, 1);
		c.setZero();
		for (uint64_t j = 1; j < nZCols; j++){
			c += (Z*V.col(j-1)).cwiseAbs();
			int unused, nxtRow;
			c.minCoeff(&nxtRow, &unused);
			V.col(j) = Z.row(nxtRow).transpose();
		}
		this->orthonormalize(V);

		//initialize X
		X.setZero(); //need this because after each projection restart X has some value

		//solve the alternating minimization for X
		double obj, prevobj;
		obj = prevobj = numeric_limits<double>::infinity();
		if (!monoCheck){
			do{
				prevobj = obj;

				this->findClosestConstrained(Z*V, X);

				this->findClosestRelaxed(Z, X, V); 

				obj = (X-Z*V).squaredNorm();
			} while( fabs(obj - prevobj)/obj > 1e-6);
		} else {
			double tmpobj;
			do{
				tmpobj = prevobj = obj;
				this->findClosestConstrained(Z*V, X);
				obj = (X-Z*V).squaredNorm();
				if (obj > (1.0+1e-6)*tmpobj){
					throw MonotonicityViolationException(tmpobj, obj, "findClosestConstrained()");
				}

				tmpobj = obj;
				this->findClosestRelaxed(Z, X, V); 
				obj = (X-Z*V).squaredNorm();
				if (obj > (1.0+1e-6)*tmpobj){
					throw MonotonicityViolationException(tmpobj, obj, "findClosestRelaxed()");
				}
			} while( fabs(obj - prevobj)/obj > 1e-6);
		}
		//compute the labels for the data
		if (obj < minObj){
			minObj = obj;
			minLbls = this->getLblsFromIndicatorMat(X);
			//there may be missing labels now (if a column was empty)
			//but getOldNewMatching below handles that case
		}
	}

	std::map<uint64_t, uint64_t> lblmap;
	double clusterCost = this->getOldNewMatching(minLbls, obs, clus, model, lblmap);

	//assign the data
	for (uint64_t i = 0; i < obs.size(); i++){
		while ( lblmap[minLbls[i]] >= clus.size()){ //use a while loop because getOldNewMatching guarantees all labels between 0 and max are used
													//therefore pushing a few empty clusters is not problematic, they will be filled after
			Clus newclus;
			clus.push_back(newclus);
		}
		clus[lblmap[minLbls[i]]].assignData(i, obs[i]);
	}
	for(uint64_t i = 0; i < clus.size(); i++){
		assert( !(clus[i].isEmpty() && clus[i].isNew()) );
	}

	//update the parameters
	for (uint64_t i = 0; i < clus.size(); i++){
		model.updatePrm(clus[i]);
	}
	return clusterCost;
}

template<class Model, bool monoCheck>
MXd _MatchingSpectral<Model, monoCheck>::getKernelMatUpper(const std::vector<typename Model::Data>& obs, const std::vector<Clus>& clus, const Model& model) const{
	const int nA = obs.size();
	const int nB = clus.size();
	MXd KUp = MXd::Zero(nA+nB, nA+nB);

	//insert ATA matrix
	for (int i = 0; i < nA; i++){
		for (int j = i; j < nA; j++){
			double sim = model.kernelDD(obs[i], obs[j]);
			if (sim > 1e-16){
				KUp(i, j) = sim;
			}
		}
	}
	//insert ATB
	for (int i = 0; i < nA; i++){
		for (int j = 0; j < nB; j++){
			double sim = model.kernelDOldP(obs[i], clus[j]);
			if (sim > 1e-16){
				KUp(i, j+nA) =  sqrt(model.oldWeight(clus[j]))*sim;
			}
		}
	}
	//insert BTB
	for (int i = 0; i < nB; i++){
		double sim = model.kernelOldPOldP(clus[i]);
		KUp(i+nA, i+nA) = model.oldWeight(clus[i])*sim+model.getOldPenalty(clus[i]);
	}
	return KUp;
}

template<class Model, bool monoCheck>
void _MatchingSpectral<Model, monoCheck>::findClosestConstrained(const MXd& ZV, MXd& X) const{
	const int nRows = ZV.rows();

	//initialize X to zero; it will be the indicator mat output
	X.setZero();
	//nonmaximum suppression
	for (int i = 0; i < nRows; i++){
		int unused, maxCol;
		ZV.row(i).maxCoeff(&unused, &maxCol);
		X(i, maxCol) = 1;
	}
	return;
}

template<class Model, bool monoCheck>
void _MatchingSpectral<Model, monoCheck>::findClosestRelaxed(const MXd& Z, const MXd& X, MXd& V) const{
	V = X.transpose()*Z;
	this->orthonormalize(V);
}

template<class Model, bool monoCheck>
void _MatchingSpectral<Model, monoCheck>::orthonormalize(MXd& V) const{
	//do a safe svd, with guaranteed orthonormal columns even in the event of numerically small singular values
	Eigen::JacobiSVD<Eigen::MatrixXd, Eigen::HouseholderQRPreconditioner> svd(V, Eigen::ComputeFullU | Eigen::ComputeFullV); //slowest/safest preconditioning
	MXd U = svd.matrixU();
	MXd W = svd.matrixV();
	MXd id = MXd::Identity(V.rows(), V.cols());
	if ( (U.transpose()*U - id).squaredNorm() > 1e-6 || 
		 (W.transpose()*W - id).squaredNorm() > 1e-6  ){
		//buggy eigen sometimes spits out non-unitary U/W matrices
		//however, this only happens when the singular value is numerical-precision small, like 1e-34ish
		//so just pull those columns cause they don't matter anyway

		//find good columns
		vector<int> goodCols;
		for (int i = 0; i < U.cols(); i++){
			if ( fabs(U.col(i).squaredNorm()-1.0) < 1e-6 && fabs(W.col(i).squaredNorm()-1.0) < 1e-6){
				goodCols.push_back(i);
			}
		}
		//shift them up -- aliasing doesn't occur since goodcols increases by at least 1 each time
		for (uint64_t i = 0; i < goodCols.size(); i++){
			U.col(i) = U.col(goodCols[i]);
			W.col(i) = W.col(goodCols[i]);
		}
		//gram schmidt them just to be sure
		for (uint64_t i = 0; i < goodCols.size(); i++){
			for (uint64_t j = 0; j < i; j++){
				double dd = U.col(i).transpose()*U.col(j);
				U.col(i) -= U.col(j)*dd;
			}
			U.col(i) /= U.col(i).norm();
			for (uint64_t j = 0; j < i; j++){
				double dd = W.col(i).transpose()*W.col(j);
				W.col(i) -= W.col(j)*dd;
			}
			W.col(i) /= W.col(i).norm();
		}
		
		//graham schmidt to create the last few columns
		//U first
		int nGoodCols = goodCols.size();
		for (int i = 0; i < U.cols(); i++){
			if (nGoodCols == U.cols()){
				break;
			}
			VXd gscol = VXd::Zero(U.rows());
			gscol(i) = 1.0;
			//subtract the projections
			for (int j = 0; j < nGoodCols; j++){
				double dd = gscol.transpose()*U.col(j);
				gscol -= dd*U.col(j);
			}
			if (gscol.squaredNorm() > 1e-6){
				U.col(nGoodCols) = gscol/gscol.norm();
				nGoodCols++;
			}
		}
		assert(nGoodCols == U.cols());
		//now W
		nGoodCols = goodCols.size();
		for (int i = 0; i < W.cols(); i++){
			if (nGoodCols == W.cols()){
				break;
			}
			VXd gscol = VXd::Zero(U.rows());
			gscol(i) = 1.0;
			//subtract the projections
			for (int j = 0; j < nGoodCols; j++){
				double dd = gscol.transpose()*W.col(j);
				gscol -= dd*W.col(j);
			}
			if (gscol.squaredNorm() > 1e-6){
				W.col(nGoodCols) = gscol/gscol.norm();
				nGoodCols++;
			}
		}
		assert(nGoodCols == W.cols());
	}
	V  = W*U.transpose();

	return; 
}

template<class Model, bool monoCheck>
vector<uint64_t> _MatchingSpectral<Model, monoCheck>::getLblsFromIndicatorMat(const MXd& X) const{
	vector<uint64_t> lbls;
	const uint64_t nR = X.rows();
	const uint64_t nC = X.cols();
	for (uint64_t i = 0; i < nR; i++){
		for (uint64_t j = 0; j < nC; j++){
			if (fabs(X(i, j)-1.0) < 1.0e-6){
				lbls.push_back(j);
				break;
			}
		}
	} 
	return lbls;
}


template<class Model, bool monoCheck>
double _MatchingSpectral<Model, monoCheck>::getOldNewMatching(const std::vector<uint64_t>& lbls, const std::vector<typename Model::Data>& obs, 
		const std::vector<Clus>& clus, const Model& model, std::map<uint64_t, uint64_t>& lblmap) const{

	//get a mapping from lbl to assigned data idcs
	std::map<uint64_t, std::vector<uint64_t> > inClus;
	for(uint64_t j = 0; j < lbls.size(); j++){
		inClus[lbls[j]].push_back(j);
	}

	//cache the old cluster kernels
	vector<double> oldKs;
	for(uint64_t i = 0; i < clus.size(); i++){
		oldKs.push_back(model.kernelOldPOldP(clus[i]));
	}

	double obsSumWt = 0.0;
	vector<int> newlbls, oldlbls;
	vector<double> wts;
	double minWt = std::numeric_limits<double>::infinity();
	for(auto it = inClus.begin(); it != inClus.end(); ++it){
		//calculate the within-cluster kernel
		double clusK = 0.0;
		for(uint64_t j = 0; j < it->second.size(); j++){
			double owt = model.kernelDD(obs[it->second[j]], obs[it->second[j]]);
			clusK += owt;
			obsSumWt += owt;
			for(uint64_t k = j+1; k < it->second.size(); k++){
				clusK += 2.0*model.kernelDD(obs[it->second[j]], obs[it->second[k]]);
			}
		}
		//push back the new cost
		newlbls.push_back(it->first);
		oldlbls.push_back(clus.size()+it->first);//unique label for this new cluster
		wts.push_back(clusK/it->second.size() - model.getNewPenalty());
		if (wts.back() < minWt){
			minWt = wts.back();
		}

		//calculate the cluster-old prm kernels
		for(uint64_t j = 0; j < clus.size(); j++){
			double ctok = 0.0;
			for(uint64_t k = 0; k < it->second.size(); k++){
				ctok += 2.0*model.kernelDOldP(obs[it->second[k]], clus[j]);
			}

			double gamma = model.oldWeight(clus[j]);
			newlbls.push_back(it->first);
			oldlbls.push_back(j);
			wts.push_back( (clusK+ gamma*(ctok - it->second.size()*oldKs[j]))/(gamma+it->second.size()) - model.getOldPenalty(clus[j]));
			if (wts.back() < minWt){
				minWt = wts.back();
			}
		}
	}

	minWt = minWt >= 0 ? minWt*0.9 : minWt*1.1;
	for (uint64_t i = 0; i < wts.size(); i++){
		wts[i] -= minWt;
	}

	MaxMatching maxm;
	map<int, int> opt = maxm.getMaxMatching(newlbls, oldlbls, wts);
	lblmap.clear();
	uint64_t nextNewLbl = clus.size();
	for(auto it = opt.begin(); it != opt.end(); ++it){
		if (it->second >= (int)clus.size()){
			lblmap[it->first] = nextNewLbl++;
		} else {
			lblmap[it->first] = it->second;
		}
	}

	return obsSumWt- (maxm.getObjective() + minWt*opt.size());
}


#define __MATCHINGSPECTRAL_IMPL_HPP
#endif /* __MATCHINGSPECTRAL_IMPL_HPP */
