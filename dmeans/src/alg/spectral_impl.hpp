#ifndef __SPECTRAL_IMPL_HPP
template<class Model, bool monoCheck>
_Spectral<Model, monoCheck>::_Spectral(const Config& cfg){
	this->solverType = this->cfg.get("eigenSolverType", Config::Type::OPTIONAL, EigenSolver::Type::EIGEN_SELF_ADJOINT);
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
double _Spectral<Model, monoCheck>::cluster(const std::vector<typename Model::Data>& obs, std::vector<Clus>& clus, const Model& model) const{
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
	uint64_t nB = clus.size();

	//premultiply Z with \hat{Gamma}^{-1/2}
	for (uint64_t j = nA; j < nA+nB; j++){
		Z.row(j) *= 1.0/sqrt(model.oldWeight(clus[j-nA])); 
	}
	const int nZCols = Z.cols(); //number of clusters currently instantiated

	//normalize the new rows of Z
	for (uint64_t j = 0; j < nA; j++){
		double rownorm = sqrt(Z.row(j).squaredNorm());
		if (rownorm == 0){ //if rownorm is a hard zero, just set the row to ones -- it's the only thing we can do, lambda was set too high
			Z.row(j) = MXd::Ones(1, nZCols);
			rownorm = sqrt(Z.row(j).squaredNorm());
		}
		Z.row(j) *= (1.0/rownorm);
	}

	//for the old rows, decide whether the sq norm is too small to be considered instantiated
	for (uint64_t j = nA; j < nA+nB; j++){
		double rowsqnorm = Z.row(j).squaredNorm();
		if (0.5/(nA+model.oldWeight(clus[j-nA])) <= rowsqnorm){
			//it's a non-degenerate row, normalize it
			Z.row(j) *= (1.0/sqrt(rowsqnorm));
		}
	}

	//initialize X (constrained version of Z) 
	MXd X(nA+nB, nZCols);
	//initialize V (rotation matrix on Z to make Z*V close to X)
	MXd V(nZCols, nZCols);

	//stuff for storing best clustering (best clustering is determined by normalized cuts objective)
	double minNCutsObj = numeric_limits<double>::infinity();
	vector<int> minLbls;

	//propose nProjectionRestarts V trials
	for (uint64_t i = 0; i < nProjectionRestarts; i++){
		V.setZero();
		//initialize unitary V via ``most orthogonal rows'' method
		int rndRow = (RNG::get())()%(nA+nB);
		V.col(0) = Z.row(rndRow).transpose();
		MXd c(nA+nB, 1);
		c.setZero();
		for (int j = 1; j < nZCols; j++){
			c += (Z*V.col(j-1)).cwiseAbs();
			int unused, nxtRow;
			c.minCoeff(&nxtRow, &unused);
			V.col(j) = Z.row(nxtRow).transpose();
		}
		this->orthonormalize(V);

		//initialize X
		X.setZero();

		//solve the alternating minimization for X
		double obj, prevobj;
		obj = prevobj = numeric_limits<double>::infinity();
		if (!monoCheck){
			do{
				prevobj = obj;

				this->findClosestConstrained(Z*V, X, nB);

				this->findClosestRelaxed(Z, X, V); 

				obj = (X-Z*V).squaredNorm();
			} while( fabs(obj - prevobj)/obj > 1e-6);
		} else {
			double tmpobj;
			do{
				tmpobj = prevobj = obj;
				this->findClosestConstrained(Z*V, X, nB);
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
		//compute the normalized cuts objective
		vector<int> tmplbls = this->getLblsFromIndicatorMat(X, nB);
		double nCutsObj = this->getNormalizedCutsObj(KUp, tmplbls);
		if (nCutsObj < minNCutsObj){
			minNCutsObj = nCutsObj;
			minLbls = tmplbls;
		}
	}
	//assign the data
	int nToAdd = (*std::max_element(minLbls.begin(), minLbls.end()))+1-clus.size();
	if (nToAdd > 0){
		for (int i = 0; i < nToAdd; i++){
			Clus newclus;
			clus.push_back(newclus);
		}
	}
	for (uint64_t i = 0; i < obs.size(); i++){
		clus[minLbls[i]].assignData(i, obs[i]);
	}
	//update the parameters
	for (uint64_t i = 0; i < clus.size(); i++){
		model.updatePrm(clus[i]);
	}
	return minNCutsObj;
}

template<class Model, bool monoCheck>
MXd _Spectral<Model, monoCheck>::getKernelMatUpper(const std::vector<typename Model::Data>& obs, const std::vector<Clus>& clus, const Model& model) const{
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
void _Spectral<Model, monoCheck>::findClosestConstrained(const MXd& ZV, MXd& X, const int nOld) const{
	const int nRows = ZV.rows();
	const int nCols = ZV.cols();
	const int nB = nOld;
	const int nA = nRows-nB;

	//initialize X to zero; it will be the indicator mat output
	X.setZero();

	//nonmaximum suppression
	for (int i = 0; i < nA; i++){
		int unused, maxCol;
		ZV.row(i).maxCoeff(&unused, &maxCol);
		X(i, maxCol) = 1;
	}

	//if there are no old rows, we're done
	if (nB == 0){
		return;
	}

	vector<int> rows, cols;
	vector<double> wts;
	for (int kk = 0; kk < nB; kk++){
		for (int jj = 0; jj < nCols; jj++){
			rows.push_back(kk+nA);
			cols.push_back(jj);
			wts.push_back(-(1.0-2.0*ZV(kk+nA, jj)));
		}
	}
	MaxMatching maxm;
	map<int, int> constrainedSoln = maxm.getMaxMatching(rows, cols, wts);

	for (auto it = constrainedSoln.begin(); it != constrainedSoln.end(); it++){
			X(it->first, it->second) = 1;
	}
	return;
}

template<class Model, bool monoCheck>
void _Spectral<Model, monoCheck>::findClosestRelaxed(const MXd& Z, const MXd& X, MXd& V) const{
	V = X.transpose()*Z;
	this->orthonormalize(V);
}

template<class Model, bool monoCheck>
void _Spectral<Model, monoCheck>::orthonormalize(MXd& V) const{
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
double _Spectral<Model, monoCheck>::getNormalizedCutsObj(const MXd& KUp, const vector<int>& lbls) const{
	const int nA = lbls.size();
	map<int, double> nums, denoms;
	for (int i = 0; i < KUp.rows(); ++i){
		for (int j = 0; j < KUp.cols(); j++){
			int l1, l2;
			//decide whether the row corresponds to an old cluster or a new datapt
			if (i >= nA){
				//the label is determined by the row
				l1 = i-nA;
			} else {
				l1 = lbls[i];
			}
			//decide whether the col corresponds to an old cluster or a new datapt
			if (j >= nA){
				//the label is determined by the col
				l2 = j-nA;
			} else {
				l2 = lbls[j];
			}
			//accumulate the numerator/denominator
			double value = i < j ? KUp(i, j) : KUp(j, i); //account for KUp being uppertri
			if (l1 != l2) {
				nums[l1] += value;//this relies on nums[] default value constructor setting to 0
			}
			denoms[l1] += value;//this relies on denoms[] default value constructor setting to 0
		}
	}
	double obj = 0;
	for (auto it = denoms.begin(); it != denoms.end(); ++it){
		if (it->second > 0){
			obj += nums[it->first]/it->second;//this relies on nums[] default value constructor setting to 0
		}
	}
	return obj;
}

template<class Model, bool monoCheck>
vector<int> _Spectral<Model, monoCheck>::getLblsFromIndicatorMat(const MXd& X, const int nOld) const{
	const int nA = X.rows()-nOld;
	const int nR = X.rows();
	const int nC = X.cols();
	vector<int> lbls;
	int nextNewLbl = nOld;
	//first find the correspondance between columns and labels
	map<int, int> colToLbl;
	for (int i = 0; i < nC; i++){
		//search for a "1" in the old entries for this column
		int oneRow = -1;
		for(int j = nA; j < nR; j++){
			if (fabs(X(j, i)-1.0) < 1.0e-6){
				oneRow = j-nA;//old cluster index
				break;
			}
		}
		if (oneRow == -1){ //if no one was found, search to make sure this column is actually instantiated at all
			bool foundOne = false;
			for (int j = 0; j <	nA; j++){
				if (fabs(X(j, i)-1.0) < 1.0e-6){
					foundOne = true;
					break;
				}
			}
			if (foundOne){
				//this is a new cluster, since no 1 was found in the old set and the new set instantiated it
				colToLbl[i] = nextNewLbl;
				nextNewLbl++;
			}
		} else { //old cluster
			colToLbl[i] = oneRow;
		}
	}
	//now push back the labels
	lbls.clear();
	for (int i = 0; i < nA; i++){
		for (int j = 0; j < nC; j++){
			if (fabs(X(i, j)-1.0) < 1e-6){
				lbls.push_back(colToLbl[j]);
				break;
			}
		}
	}
	return lbls;
}


#define __SPECTRAL_IMPL_HPP
#endif /* __SPECTRAL_IMPL_HPP */
