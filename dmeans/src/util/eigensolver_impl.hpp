#ifndef __EIGENSOLVER_IMPL_HPP
EigenSolver::EigenSolver(MXd& A_UpperTriangle, Type t, uint64_t nEigs, double lowerThresh){
	if (A_UpperTriangle.rows() != A_UpperTriangle.cols()){
		throw MatrixNotSquareException(A_UpperTriangle.rows() , A_UpperTriangle.cols());
	}
	switch(t){
		case EIGEN_SELF_ADJOINT:
			this->selfadjointSolver(A_UpperTriangle);
		case REDSVD:
			nEigs = (nEigs==0 ? A_UpperTriangle.rows() : nEigs);
			nEigs = (nEigs > A_UpperTriangle.rows() ? A_UpperTriangle.rows() : nEigs);
			this->redsvdSolver(A_UpperTriangle, nEigs);
	}
	if (lowerThresh > 0){
		this->pruneSmallEigs(lowerThresh);
	}
}

//This code was adapted from https://code.google.com/p/redsvd/wiki/English Copyright (c) 2010 Daisuke Okanohara
//I made a number of modifications to make it more efficient (sparse matrix computation, faster gram schmidt, faster gaussian sampling)
//To learn about it, please see "Finding structure with randomness: Stochastic algorithms for constructing approximate matrix
//decompositions", N. Halko, P.G. Martinsson, J. Tropp, arXiv 0909.4061
void EigenSolver::redsvdSolver(MXd& AUp, uint64_t r){
	r = (r < AUp.cols()) ? r : AUp.cols();
	//compute gaussian matrix
	normal_distribution<> nrm(0, 1);
	MXd M(AUp.rows(), r);
	for (int i = 0; i < AUp.rows(); i++){
		for (int j =0 ; j < r; j++){
			M(i, j) = nrm(RNG::get());
		}
	}
	//compute Y 
	MXd Y = (AUp.selfadjointView<Eigen::Upper>())*M;
	//orthonormalize Y -- Gram Schmidt
	this->columnGramSchmidt(Y);
	//run it again if necessary (numerical precision increase)
	if (  (Y.transpose()*Y-MXd::Identity(Y.cols(), Y.cols())).cwiseAbs().maxCoeff() > 1e-6 ){
		this->columnGramSchmidt(Y);
	}
	MXd B = Y.transpose()*((AUp.selfadjointView<Eigen::Upper>())*Y);
	Eigen::SelfAdjointEigenSolver<MXd> eigB(B);
	//since the eigenvalues are sorted in increasing order, chop off the ones at the front
	this->eigvals = eigB.eigenvalues();
	this->eigvecs = Y*(eigB.eigenvectors());
}

void EigenSolver::selfadjointSolver(MXd& AUp){
	Eigen::SelfAdjointEigenSolver<MXd> eigB;
	eigB.compute(AUp.selfadjointView<Eigen::Upper>());
	this->eigvals = eigB.eigenvalues();
	this->eigvecs = eigB.eigenvectors();
}

void EigenSolver::columnGramSchmidt(MXd& m){
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

void EigenSolver::pruneSmallEigs(double thresh){
	//assumes the eigenvalues are sorted in increasing order
	if (eigvals(eigvals.size()-1) < thresh){
		this->eigvecs = VXd();
		this->eigvals = MXd();
		return;
	}
	int chopIdx = 0;
	while (chopIdx < eigvals.size() && eigvals(chopIdx) < thresh) chopIdx++; 

	int nLeftOver = eigvals.size()-chopIdx;
	this->eigvecs = eigvecs.topRightCorner(eigvecs.rows(), nLeftOver);
	this->eigvals = eigvals.tail(nLeftOver);
}

void EigenSolver::getResults(VXd& eigvals, MXd& eigvecs){
	eigvals = this->eigvals;
	eigvecs = this->eigvecs;
}

#define __EIGENSOLVER_IMPL_HPP
#endif /* __EIGENSOLVER_IMPL_HPP */
