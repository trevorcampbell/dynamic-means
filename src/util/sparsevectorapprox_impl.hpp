#ifndef __SPARSEVECTORAPPROX_IMPL_HPP

#include <iostream>
#include <numeric>

SparseVectorApproximation::SparseVectorApproximation(int K, double eps){
	this->K = K;
	this->eps = eps;
}

void SparseVectorApproximation::fromKernelMatrix(MXd m, VXd acoeffs){
	std::cout << "SpVecApprox: Computing sparse linear combination approximation via Kernel Gram-Schmidt..." << std::endl;
	std::cout << "SpVecApprox: Max basis size = " << this->K << std::endl;
	std::cout << "SpVecApprox: Error threshold = " << this->eps << std::endl; 
	MXd msave = m;
	//clear old results
	this->bvecs.clear();
	this->coeffs.clear();

	//find the alpha-weighted column sum for computing the cost reduction for each column
	VXd acolsums = m*acoeffs;
	//find the initial cost (norm of the linear combination)
	double cost = acoeffs.transpose()*acolsums;
	std::cout << "SpVecApprox: Initial cost = " << cost << std::endl;

	//store the current set of indices
	std::vector<int> idcs(acoeffs.size()); //idcs stores the mapping from current matrix row/col to original matrix row/col
	std::iota(idcs.begin(), idcs.end(), 0); 
	std::vector<int> srtidcs = idcs; //srtidcs stores a list of current matrix row/cols sorted by cost reduction (best is .back())

	//gram schmidt loop 
	while (cost > this->eps && this->bvecs.size() < this->K){
		//sort to find the index with maximum cost reduction
		VXd costreds = acolsums.array()*acolsums.array()/(m.diagonal().array());
		std::sort(srtidcs.begin(), srtidcs.end(), 
				[&costreds](const int & a, const int & b) -> bool
				{ return costreds(a) < costreds(b); });
		//add the best vector to the basis; project everything onto orthogonal subspace and remove the vector from further consideration
		int kMax = srtidcs.back();
		cost -= costreds(kMax);
		this->bvecs.push_back(idcs[kMax]);
		this->projectOntoOrthogonalSubspace(acolsums, m, idcs, srtidcs, kMax);
		std::cout << "SpVecApprox: Basis size = " << this->bvecs.size() << " Cost = " << cost << std::endl;
	}

	//final computation of coefficients
	MXd W = MXd::Zero(this->bvecs.size(), this->bvecs.size());
	MXd B = MXd::Zero(this->bvecs.size(), msave.cols());
	for(int i = 0; i < this->bvecs.size(); i++){
		for (int j = 0; j < this->bvecs.size(); j++){
			W(i, j) = msave(this->bvecs[i], this->bvecs[j]);
		}
		for (int j = 0; j < msave.cols(); j++){
			B(i, j) = msave(this->bvecs[i], j);
		}
	}
	//for (int i =0 ; i < bvecs.size(); i++){
	//	std::cout << bvecs[i] << " ";
	//}
	//std::cout << std::endl;
	//std::cout << "W: " << std::endl << W << std::endl;
	//std::cout << "B: " << std::endl << B << std::endl;

	std::cout << "SpVecApprox: Computing coefficients using LDLT decomposition" << std::endl;
	VXd finalcoeffs = W.ldlt().solve(B*acoeffs);
	for(int i = 0; i < finalcoeffs.size(); i++){
		this->coeffs.push_back(finalcoeffs(i));
	}
	std::cout << "SpVecApprox: Done!"<< std::endl;
}


void SparseVectorApproximation::projectOntoOrthogonalSubspace(VXd& acolsums, MXd& m, std::vector<int>& idcs, std::vector<int>& srtidcs, const int kMax){
	//rank one update to the column sums
	double ack = acolsums(kMax)/m(kMax, kMax);
	for(int i =0 ; i < acolsums.size(); i++){
		acolsums(i) -= ack*m(kMax, i);
	}
	//rank one update to the gram matrix
	m -= (m.col(kMax)*m.row(kMax)/m(kMax, kMax)).eval(); //eval to prevent aliasing

	//remove the row/column corresponding to kMax
	removeRowCol(m, kMax);
	removeElem(acolsums, kMax);
	idcs.erase(idcs.begin()+kMax);
	for(int i =0 ; i < srtidcs.size(); i++){
		if (srtidcs[i] >= kMax){
			srtidcs[i]--;
		}
	}
	srtidcs.pop_back();
	return;
}

void SparseVectorApproximation::removeRowCol(MXd& m, int k){
	uint32_t nRows = m.rows();
	uint32_t nCols = m.cols();
	if (k < nRows-1){
		m.block(k, 0, nRows-1-k, nCols) = m.block(k+1, 0, nRows-1-k, nCols).eval();
		m.block(0, k, nRows, nCols-1-k) = m.block(0, k+1, nRows, nCols-1-k).eval();
	}
	m.conservativeResize(nRows-1, nCols-1);
	return;
}

void SparseVectorApproximation::removeElem(VXd& v, int k){
	uint32_t nRows = v.size();
	if (k < nRows-1){
		v.segment(k, nRows-1-k) = v.tail(nRows-1-k).eval();
	}
	v.conservativeResize(nRows-1);
	return;
}


void SparseVectorApproximation::getApprox(std::vector<int>& vecs, std::vector<double>& acoeffs){
	vecs = this->bvecs;
	acoeffs = this->coeffs;
	return;
}



#define __SPARSEVECTORAPPROX_IMPL_HPP
#endif /* __SPARSEAPPROX_IMPL_HPP */
