#ifndef __SPECDYNMEANS_IMPL_HPP


template <typename G>
SpecDynMeans<G>::SpecDynMeans(double lamb, double Q, double tau, bool verbose /* = false*/, int seed /*= -1*/) {
	this->verbose = verbose;
	if (lamb < 0 || Q < 0 || tau < 0){
		cout << "libspecdynmeans: ERROR: Parameters of Spectral Dynamic Means cannot be < 0." << endl;
		cout << "libspecdynmeans: Lambda: " << lamb << " Q: " << Q << " tau: " << tau << endl;
	}
	this->lamb = lamb;
	this->Q = Q;
	this->tau = tau;
	this->ages.clear();
	this->weights.clear();
	this->gammas.clear();
	this->agecosts.clear();
	this->oldprmlbls.clear();
	this->maxLblPrevUsed = -1;
	//seed the random number generator with time now (seed < 0) or seed (seed > 0)
	if(seed < 0){
		this->rng.seed(unsigned( time(0) ) );
	} else{
		this->rng.seed(seed);
	}
	grbenv = new GRBEnv();
	grbenv->set(GRB_IntParam_OutputFlag, 0);//controls the output of Gurobi - 0 means no output, 1 means normal output
	//grbenv->set(GRB_IntParam_Method, 1);//controls which method Gurobi uses - default -1 (auto), 0=primal simplex, 1=dual simplex, 2=barrier, 3=concurrent, 4=deterministic concurrent
	grbenv->set(GRB_IntParam_Threads, 1);//controls the number of threads Gurobi uses - I force it to use 1 since the optimization in this algorithm is fairly small/simple
									    // 												and it just ends up wasting time constantly creating/deleting threads otherwise
}

//Just the destructor for the class

template <typename G>
SpecDynMeans<G>::~SpecDynMeans(){
	delete grbenv;//make sure the gurobi environment is destroyed
}

//Reset returns the class object to its initial state, ready to start a new DDP chain

template <typename G>
void SpecDynMeans<G>::reset(){
	this->ages.clear();
	this->weights.clear();
	this->gammas.clear();
	this->agecosts.clear();
	this->oldprmlbls.clear();
	this->maxLblPrevUsed = -1;
}

//This function updates the weights/ages of all the clusters after each clustering step is complete
//This function updates the weights/ages of all the clusters after each clustering step is complete
template <typename G>
void SpecDynMeans<G>::finalizeStep(const G& aff, const vector<int>& lbls, vector<double>& prevgammas_out, vector<int>& prmlbls_out){
	//first increment the age of everything
	//this will be undone below for any current clusters
	for (int i = 0; i < this->ages.size(); i++){
		this->ages[i]++;
	}
	//get the set of current cluster labels
	vector<int> unqlbls = lbls;
	sort(unqlbls.begin(), unqlbls.end());
	unqlbls.erase(unique(unqlbls.begin(), unqlbls.end()), unqlbls.end());
	//for every current cluster
	for (int i = 0; i < unqlbls.size(); i++){
		const int& lbl = unqlbls[i];
		//get the number of nodes assigned to the cluster
		int nclus = 0;
		for (int j = 0; j < lbls.size(); j++){
			nclus += (lbls[j] == lbl ? aff.getNodeCt(j) : 0);
		}
		//check if it's an old label
		auto it2 = find(this->oldprmlbls.begin(), this->oldprmlbls.end(), lbl);
		if (it2 == this->oldprmlbls.end()){
			//it's a new cluster, so create stuff for it
			this->ages.push_back(1);
			this->weights.push_back(nclus);
			this->oldprmlbls.push_back(lbl);
		} else {
			//it's an old cluster, so update old stuff
			int oldidx = distance(this->oldprmlbls.begin(), it2);
			this->ages[oldidx] = 1;
			this->weights[oldidx] = this->gammas[oldidx] + nclus;
		}
	}
	//save the old gammas for output, update this->gammas
	//prevgammas_out still needs to be cleaned up and modified a bit (occurs in the rest of this function)
	prevgammas_out = this->gammas;
	//pad prevgammas with zeros for all new clusters
	prevgammas_out.insert(prevgammas_out.end(), this->ages.size()-prevgammas_out.size(), 0);
	//update the internal gammas
	this->gammas.clear();
	this->gammas.reserve(this->ages.size());
	for (int i = 0; i < this->ages.size(); i++){
			this->gammas.push_back(1.0/(1.0/this->weights[i] + this->ages[i]*this->tau)); 
	}
	//update the Age Costs
	this->agecosts.clear();
	this->agecosts.reserve(this->ages.size());
	for (int i = 0; i < this->ages.size(); i++){
			this->agecosts.push_back(this->Q*this->ages[i]); 
	}

	//delete any old cluster whose age cost exceeds lambda
	for (int i = 0; i < this->ages.size(); i++){
		if (this->agecosts[i] > this->lamb){
			this->weights.erase(this->weights.begin()+i);
			this->ages.erase(this->ages.begin()+i);
			this->gammas.erase(this->gammas.begin()+i);
			prevgammas_out.erase(prevgammas_out.begin()+i);
			this->agecosts.erase(this->agecosts.begin()+i);
			this->oldprmlbls.erase(this->oldprmlbls.begin()+i);
			i--;
		}
	}
	prmlbls_out = this->oldprmlbls;//save the parameter labels output 

	//update this->maxLblPrevUsed if new clusters were created
	int maxlbl = *max_element(lbls.begin(), lbls.end());
	if (maxlbl >= this->maxLblPrevUsed){
		this->maxLblPrevUsed = maxlbl;
	}

	return;
}

template <typename G>
void SpecDynMeans<G>::cluster(const G& aff, const int nRestarts, const int nClusMax, EigenSolverType type, vector<int>& finalLabels, double& finalObj, 
				std::vector<double>& finalGammas, std::vector<int>& finalPrmLbls, double& tTaken){

	timeval tStart;
	gettimeofday(&tStart, NULL);

	const int nB = this->oldprmlbls.size();
	const int nA = aff.getNNodes();

	if (nA <= 0){
		cout << "libspecdynmeans: WARNING: data size <=0 (= " << nA << "); Returning empty labels."<<  endl;
		finalLabels = vector<int>();
		tTaken = 0;
		finalObj = 0;
		return;
	}
	if (nRestarts <= 0){
		cout << "libspecdynmeans: ERROR: nRestarts <=0 (= " << nRestarts << ")"<<  endl;
		return;
	}
	if (verbose){
		cout << "libspecdynmeans: Clustering " << nA << " datapoints with " << nRestarts << " restarts." << endl;
	}


	//compute the kernel matrix
	SMXd kUpper;
	this->getKernelMat(aff, kUpper);
	//solve the eigensystem
	MXd Z;
	VXd eigvals;
	if (verbose){
		cout << "libspecdynmeans: Solving the eigensystem..." << flush;
	}
	this->solveEigensystem(kUpper, nClusMax+nB, type, Z); //number of eigvecs = # of clusters to track + 
													//number of old prms (for possible old clus indicator vecs)
	if (verbose){
		cout << "Done!" << endl;
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
	for (int i = 0; i < nRestarts; i++){
		V.setZero();
		//initialize unitary V via ``most orthogonal rows'' method
		int rndRow = this->rng()%(nA+nB);
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

		if (verbose){
			cout << "libspecdynmeans: Attempt " << i+1 << "/" << nRestarts << ", Obj = " << nCutsObj << " MinObj = " << minNCutsObj << "                                       \r" <<  flush;
		}
	}

	if (verbose){
		vector<int> unqlbls = minLbls;
		sort(unqlbls.begin(), unqlbls.end());
		unqlbls.erase(unique(unqlbls.begin(), unqlbls.end()), unqlbls.end());
		int numnew = 0;
		for (int i = 0; i < unqlbls.size(); i++){
			if (find(this->oldprmlbls.begin(), this->oldprmlbls.end(), unqlbls[i]) == this->oldprmlbls.end()){
				numnew++;
			}
		}
		int numoldinst = unqlbls.size() - numnew;
		int numolduninst = this->ages.size() - numoldinst;
		cout << endl << "libspecdynmeans: Done clustering. Min Objective: " << minNCutsObj << " Old Uninst: " << numolduninst  << " Old Inst: " << numoldinst  << " New: " << numnew <<  endl;
	}

	//update the state of the ddp chain
	this->finalizeStep(aff, minLbls, finalGammas, finalPrmLbls);

	//output results
	finalObj =  minNCutsObj;
	finalLabels = minLbls;
	//get final time taken
	timeval tCur;
	gettimeofday(&tCur, NULL);
	tTaken = (double)(tCur.tv_sec - tStart.tv_sec) + (double)(tCur.tv_usec - tStart.tv_usec)/1.0e6;
	return;
}

template <typename G>
void SpecDynMeans<G>::getKernelMat(const G& aff, SMXd& kUpper){
	const int nA = aff.getNNodes();
	const int nB = this->oldprmlbls.size();
	kUpper.resize(nA+nB, nA+nB);

	vector<TD> kUpperTrips;
	kUpperTrips.reserve(nA+nB);//just start with some heuristic size
	//insert ATA matrix
	for (int i = 0; i < nA; i++){
		for (int j = i; j < nA; j++){
			double sim = aff.simDD(i, j);
			if (sim > 1e-16){
				kUpperTrips.push_back(TD(i, j, sim));
			}
		}
	}
	//insert ATB
	for (int i = 0; i < nA; i++){
		for (int j = 0; j < nB; j++){
			double sim = aff.simDP(i, j);
			if (sim > 1e-16){
				kUpperTrips.push_back(TD(i, j+nA, sqrt(this->gammas[j])*sim));
			}
		}
	}
	//insert BTB
	for (int i = 0; i < nB; i++){
		kUpperTrips.push_back(TD(i+nA, i+nA, this->gammas[i]+this->agecosts[i]));
	}

	kUpper.setFromTriplets(kUpperTrips.begin(), kUpperTrips.end());
	return;
}


template <typename G>
void SpecDynMeans<G>::solveEigensystem(SMXd& kUpper, const int nEigs, EigenSolverType type, MXd& eigvecs){
	const int nB = this->ages.size();
	const int nA = kUpper.rows()-nB;

	VXd eigvals;
	if (type == EIGEN_SELF_ADJOINT){
			Eigen::SelfAdjointEigenSolver<MXd> eigB;
			eigB.compute(MXd(kUpper));
			//since the eigenvalues are sorted in increasing order, chop off the ones at the front
			eigvals = eigB.eigenvalues();
			eigvecs = eigB.eigenvectors();
			int chopIdx = 0;
			while (chopIdx < eigvals.size() && eigvals(chopIdx) < this->lamb) chopIdx++; 
			if (chopIdx == eigvals.size()){
				eigvals = eigvals.tail(1).eval();
				eigvecs = eigvecs.col(eigvecs.cols()-1).eval();
			} else {
				int nLeftOver = eigvals.size()-chopIdx;
				eigvals = eigvals.tail(nLeftOver).eval();
				eigvecs = eigvecs.topRightCorner(eigvecs.rows(), nLeftOver).eval();
			}
	} else {
			std::tie(eigvecs, eigvals) = std::move(redsvdEigenSolver(kUpper, nEigs));
	}

    if (verbose){
    	cout << "Eigenvalues: ";
    	for (int i = 0; i < eigvals.rows(); i++){
    		cout << eigvals(i) << " ";
		}
		cout << endl;
	}
	return;
}

//This code was adapted from https://code.google.com/p/redsvd/wiki/English Copyright (c) 2010 Daisuke Okanohara
//I made a number of modifications to make it more efficient (sparse matrix computation, faster gram schmidt, faster gaussian sampling)
//To learn about it, please see "Finding structure with randomness: Stochastic algorithms for constructing approximate matrix
//decompositions", N. Halko, P.G. Martinsson, J. Tropp, arXiv 0909.4061
template <typename G>
tuple<MXd, VXd> SpecDynMeans<G>::redsvdEigenSolver(SMXd& AUp, int r){
	r = (r < AUp.cols()) ? r : AUp.cols();
	//compute gaussian matrix
	normal_distribution<> nrm(0, 1);
	MXd M(AUp.rows(), r);
	for (int i = 0; i < AUp.rows(); i++){
		for (int j =0 ; j < r; j++){
			M(i, j) = nrm(this->rng);
		}
	}
	//compute Y 
	MXd Y = (AUp.selfadjointView<Eigen::Upper>())*M;
	//orthonormalize Y -- Gram Schmidt
	this->gramschmidt(Y);
	//run it again if necessary (numerical precision increase)
	if (  (Y.transpose()*Y-MXd::Identity(Y.cols(), Y.cols())).cwiseAbs().maxCoeff() > 1e-6 ){
		this->gramschmidt(Y);
	}
	MXd B = Y.transpose()*((AUp.selfadjointView<Eigen::Upper>())*Y);
	Eigen::SelfAdjointEigenSolver<MXd> eigB(B);
	//since the eigenvalues are sorted in increasing order, chop off the ones at the front
	VXd eigvals = eigB.eigenvalues();
	MXd eigvecs = Y*(eigB.eigenvectors());
	int chopIdx = 0;
	while (chopIdx < eigvals.size() && eigvals(chopIdx) < this->lamb) chopIdx++; 
	if (chopIdx == eigvals.size()){
		return tuple<MXd, VXd>(eigvecs.col(eigvecs.cols()-1), eigvals.tail(1));
	} else {
		int nLeftOver = eigvals.size()-chopIdx;
		return tuple<MXd, VXd>(eigvecs.topRightCorner(eigvecs.rows(), nLeftOver), eigvals.tail(nLeftOver));
	}
}


template <typename G>
void SpecDynMeans<G>::findClosestConstrained(const MXd& ZV, MXd& X) const{
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


template <typename G>
map<int, int> SpecDynMeans<G>::getOldNewMatching(vector< pair<int, int> > nodePairs, vector<double> edgeWeights ) const{
	//get params
	int nVars = edgeWeights.size();

	//start up GRB
	try{
	GRBModel grbmodel(*grbenv);
	//add variables/objective
	double* obj = new double[nVars];
	for (int i = 0; i < nVars; i++){
		obj[i] = edgeWeights[i];
	}
	GRBVar* grbvars = grbmodel.addVars(NULL, NULL,obj, NULL, NULL, nVars);

	grbmodel.update();

	//one constraint for each A/B node, plus one constraint for each edge
	vector<int> A, B;
	for (int i = 0; i < nodePairs.size(); i++){
		if (find(A.begin(), A.end(), nodePairs[i].first) == A.end()){
			A.push_back(nodePairs[i].first);
		}
		if (find(B.begin(), B.end(), nodePairs[i].second) == B.end()){
			B.push_back(nodePairs[i].second);
		}
	}
	//add constraints
	//constraint type 1: sum of outgoing edges from A nodes = 1
	for (int i = 0; i < A.size(); i++){
		GRBLinExpr constrlhs;
		for (int j = 0; j < nVars; j++){
			if (nodePairs[j].first == A[i]){
				constrlhs += 1.0*grbvars[j];
			}
		}
		grbmodel.addConstr(constrlhs, GRB_EQUAL, 1);
	}
	
	//constraint type 2: sum of incoming edges to B nodes <= 1
	for (int i = 0; i < B.size(); i++){
		GRBLinExpr constrlhs;
		for (int j = 0; j < nVars; j++){
			if (nodePairs[j].second == B[i]){
				constrlhs += 1.0*grbvars[j];
			}
		}
		grbmodel.addConstr(constrlhs, GRB_LESS_EQUAL, 1);
	}
	//constraint type 3: all edge variables >= 0
	//simplex has an implicit bound of >=0 on all variables, don't need this

	grbmodel.optimize();

	map<int, int> retmap;
	for (int j = 0; j < edgeWeights.size(); j++){
		double val = grbvars[j].get(GRB_DoubleAttr_X);
		if (fabs(val - 1.0) < 1e-10){
			retmap[nodePairs[j].first] = nodePairs[j].second;
			A.erase(remove(A.begin(), A.end(), nodePairs[j].first), A.end());
		}
	}
	delete[] grbvars;
	delete[] obj;
	return retmap;
	} catch (GRBException e){
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	} catch (...){
		cout << "Exception during optimization" << endl;
	}
}


template <typename G>
void SpecDynMeans<G>::findClosestRelaxed(const MXd& Z, const MXd& X, MXd& V) const{ 
	V = X.transpose()*Z;
	this->orthonormalize(V);

	return; 
}


template <typename G>
void SpecDynMeans<G>::orthonormalize(MXd& V) const{ 
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


template <typename G>
vector<int> SpecDynMeans<G>::getLblsFromIndicatorMat(const MXd& X) const {
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


template <typename G>
double SpecDynMeans<G>::getNormalizedCutsObj(const SMXd& spmatUpper, const vector<int>& lbls) const{
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

template<typename G>
void SpecDynMeans<G>::gramschmidt(MXd& m){
	for (int i = 0; i < m.cols(); i++){
		if (m.col(i).norm() < 1e-8){
			//remove it if it's super small (pretty much dependent with previous vectors)
			if (i < m.cols()-1){
				m.block(0, i, m.rows(), m.cols()-i-1) = m.block(0, i+1, m.rows(), m.cols()-i-1).eval();
			}
			m.conservativeResize(m.rows(), m.cols()-1);
			i--;
			continue;
		}
		m.col(i) /= m.col(i).norm();
		for (int j = i+1; j < m.cols(); j++){
			m.col(j) -= m.col(j).dot(m.col(i))*m.col(i);
		}
	}
}

#define __SPECDYNMEANS_IMPL_HPP
#endif /* __SPECDYNMEANS_IMPL_HPP */
