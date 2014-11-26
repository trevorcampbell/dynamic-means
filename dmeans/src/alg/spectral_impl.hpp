#ifndef __SPECTRAL_IMPL_HPP
Spectral::Spectral(const Config& cfg){
	this->solverType = this->cfg.get("eigenSolverType", Config::Type::OPTIONAL, EigenSolver::Type::EIGEN_SELF_ADJOINT);
	if (this->solverType == EigenSolver::Type::EIGEN_SELF_ADJOINT){
		this->nEigs = 0;
	} else {
		this->nEigs = this->cfg.get("eigenDimension", Config::Type::REQUIRED, 0);
	}
	this->verbose = cfg.get("verbose", Config::Type::OPTIONAL, false);
	this->nProjectionRestarts = this->cfg.get("nProjectionRestarts", Config::Type::OPTIONAL, 1);
}

double Spectral::cluster(const std::vector<typename Model::Data>& obs, std::vector<Clus>& clus, const Model& model) const{
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

	//premultiply Z with \hat{Gamma}^{-1/2}
	for (int j = nA; j < nA+nB; j++){
		Z.row(j) *= 1.0/sqrt(this->gammas[j-nA]);
	}
	const int nZCols = Z.cols(); //number of clusters currently instantiated

	//normalize the new rows of Z
	for (int j = 0; j < nA; j++){
		double rownorm = sqrt(Z.row(j).squaredNorm());
		if (rownorm == 0){ //if rownorm is a hard zero, just set the row to ones -- it's the only thing we can do, lambda was set too high
			Z.row(j) = MXd::Ones(1, nZCols);
			rownorm = sqrt(Z.row(j).squaredNorm());
		}
		Z.row(j) *= (1.0/rownorm);
	}

	//for the old rows, decide whether the sq norm is too small to be considered instantiated
	for (int j = nA; j < nA+nB; j++){
		double rowsqnorm = Z.row(j).squaredNorm();
		if (0.5/(nA+this->gammas[j-nA]) <= rowsqnorm){
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

	//propose nRestarts V trials
	if (verbose){
		cout << "libspecdynmeans: Finding discretized solution with " << nRestarts << " restarts" << endl;
	}
	for (int i = 0; i < nProjectionRestarts; i++){
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
		do{
			prevobj = obj;

			this->findClosestConstrained(Z*V, X);

			this->findClosestRelaxed(Z, X, V); 

			obj = (X-Z*V).squaredNorm();
		} while( fabs(obj - prevobj)/obj > 1e-6);
		//compute the normalized cuts objective
		vector<int> tmplbls = this->getLblsFromIndicatorMat(X);
		double nCutsObj = this->getNormalizedCutsObj(kUpper, tmplbls);
		if (nCutsObj < minNCutsObj){
			minNCutsObj = nCutsObj;
			minLbls = tmplbls;
		}
	}

	return;
}

MXd& Spectral::getKernelMatUpper(const std::vector<typename Model::Data>& obs, const std::vector<Clus>& clus, const Model& model) const{
	const int nA = obs.size();
	const int nB = clus.size();
	MXd KUp = MXd::Zeros(nA+nB, nA+nB);

	//insert ATA matrix
	for (int i = 0; i < nA; i++){
		for (int j = i; j < nA; j++){
			double sim = aff.simDD(i, j);
			if (sim > 1e-16){
				KUp(i, j) = sim;
			}
		}
	}
	//insert ATB
	for (int i = 0; i < nA; i++){
		for (int j = 0; j < nB; j++){
			double sim = aff.simDP(i, j);
			if (sim > 1e-16){
				KUp(i, j+nA) =  sqrt(this->gammas[j])*sim;
			}
		}
	}
	//insert BTB
	for (int i = 0; i < nB; i++){
		double sim = aff.simPP(i);
		KUp(i+nA, i+nA) = this->gammas[i]*sim+this->agecosts[i];
	}
	return KUp;
}

void Spectral::findClosestConstrained(const MXd& ZV, MXd& X) const{
	const int nRows = ZV.rows();
	const int nCols = ZV.cols();
	const int nB = this->ages.size();
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
	//otherwise, create the LP and solve it
	vector< pair<int, int> > nodePairs;
	vector<double> edgeWeights;
	for (int kk = 0; kk < nB; kk++){
		//push back the null cost
		double rsqnorm = ZV.row(kk+nA).squaredNorm();
		nodePairs.push_back( pair<int, int>(kk+nA, nCols+kk) ); //nCols+kk is the unique ID for this row's null choice
		edgeWeights.push_back(rsqnorm); //nCols+kk is the unique ID for this row's null choice
		rsqnorm += 1.0;
		for (int jj = 0; jj < nCols; jj++){
			nodePairs.push_back(pair<int,int>(kk+nA, jj));
			edgeWeights.push_back(rsqnorm-2.0*ZV(kk+nA,jj));
		}
	}
	map<int, int> constrainedSoln = getOldNewMatching(nodePairs, edgeWeights);
	for (auto it = constrainedSoln.begin(); it != constrainedSoln.end(); it++){
		if (it->second < nCols){ //if it's >= then it's a null row so do nothing since X was already set to zero earlier
			X(it->first, it->second) = 1;
		}
	}
	return;


}

void Spectral::findClosestRelaxed(const MXd& Z, const MXd& X, MXd& V) const{
	V = X.transpose()*Z;
	this->orthonormalize(V);
}

void Spectral::orthonormalize(MXd& V) const{
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
		for (int i = 0; i < goodCols.size(); i++){
			U.col(i) = U.col(goodCols[i]);
			W.col(i) = W.col(goodCols[i]);
		}
		//gram schmidt them just to be sure
		for (int i = 0; i < goodCols.size(); i++){
			for (int j = 0; j < i; j++){
				double dd = U.col(i).transpose()*U.col(j);
				U.col(i) -= U.col(j)*dd;
			}
			U.col(i) /= U.col(i).norm();
			for (int j = 0; j < i; j++){
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
		if(nGoodCols != U.cols()){
			cout << "libspecdynmeans: ERROR: GRAHAM SCHMIDT DID NOT ORTHONORMALIZE U!" << endl;
			int ddd;
			cin >> ddd;
		}
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
		if(nGoodCols != W.cols()){
			cout << "libspecdynmeans: ERROR: GRAHAM SCHMIDT DID NOT ORTHONORMALIZE W!" << endl;
			int ddd;
			cin >> ddd;
		}
	}
	V  = W*U.transpose();

	return; 


}

double Spectral::getNormalizedCutsObj(const MXd& mUp, const vector<int>& lbls) const{
const int nA = lbls.size();
	const int nB = this->ages.size();
	map<int, double> nums, denoms;
	for (int i = 0; i < spmatUpper.outerSize(); ++i){
		for (SMXd::InnerIterator it(spmatUpper, i); it; ++it){
			int l1, l2;
			//decide whether the row corresponds to an old cluster or a new datapt
			if (it.row() >= nA){
				//the label is determined by the row
				l1 = this->oldprmlbls[it.row() - nA];
			} else {
				l1 = lbls[it.row()]; 
			}
			//decide whether the col corresponds to an old cluster or a new datapt
			if (it.col() >= nA){
				//the label is determined by the col
				l2 = this->oldprmlbls[it.col() - nA];
			} else {
				l2 = lbls[it.col()];
			}
			//accumulate the numerator/denominator
			if (l1 != l2) {
				nums[l1] += it.value();//this relies on nums[] default value constructor setting to 0
			}
			denoms[l1] += it.value();//this relies on denoms[] default value constructor setting to 0
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

vector<int> Spectral::getLblsFromIndicatorMat(const MXd& X) const{
const int nB = this->ages.size();
	const int nA = X.rows()-nB;
	const int nR = X.rows();
	const int nC = X.cols();
	vector<int> lbls;
	int nextNewLbl = this->maxLblPrevUsed+1;
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
			colToLbl[i] = this->oldprmlbls[oneRow];
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
