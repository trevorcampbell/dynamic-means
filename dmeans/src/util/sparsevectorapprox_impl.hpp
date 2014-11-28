#ifndef __SPARSEVECTORAPPROX_IMPL_HPP

SparseVectorApproximation::SparseVectorApproximation(uint64_t K, double eps){
	this->K = K;
	this->eps = eps;
}

void SparseVectorApproximation::fromKernelMatrix(MXd m, VXd acoeffs){
	//std::cout << "SpVecApprox: Computing sparse linear combination approximation via Kernel Gram-Schmidt..." << std::endl;
	//std::cout << "SpVecApprox: Max basis size = " << this->K << std::endl;
	//std::cout << "SpVecApprox: Relative error threshold = " << this->eps << std::endl; 
	MXd msave = m;
	//clear old results
	this->bvecs.clear();
	this->coeffs.clear();

	//jump out early if K is larger than the dimension of the matrix
	if ((uint64_t)m.rows() < this->K){
		//std::cout << "SpVecApprox: Input vector size < K, returning the full set."<< std::endl;
		for(int i = 0; i < m.rows(); i++){
			this->bvecs.push_back(i);
			this->coeffs.push_back(acoeffs(i));
		}
		return;
	}

	//find the alpha-weighted column sum for computing the cost reduction for each column
	VXd acolsums = m*acoeffs;
	//find the initial cost (sqnorm of the linear combination we're trying to approximate)
	double vecSqNorm = acoeffs.transpose()*acolsums;
	//std::cout << "SpVecApprox: Vector Square Norm = " << vecSqNorm << std::endl;

	//store the current set of indices
	std::vector<uint64_t> idcs(acoeffs.size()); //idcs stores the mapping from current matrix row/col to original matrix row/col
	std::iota(idcs.begin(), idcs.end(), 0); 
	std::vector<uint64_t> srtidcs = idcs; //srtidcs stores a list of current matrix row/cols sorted by cost reduction (best is .back())

	//gram schmidt loop 
	double errSqNorm = vecSqNorm;
	while (sqrt(errSqNorm/vecSqNorm) > this->eps && this->bvecs.size() < this->K){
		//sort to find the index with maximum cost reduction
		VXd costreds = acolsums.array()*acolsums.array()/(m.diagonal().array());
		std::sort(srtidcs.begin(), srtidcs.end(), 
				[&costreds](const uint64_t & a, const uint64_t & b) -> bool
				{ return costreds(a) < costreds(b); });
		//add the best vector to the basis; project everything onto orthogonal subspace and remove the vector from further consideration
		uint64_t kMax = srtidcs.back();
		errSqNorm -= costreds(kMax);
		this->bvecs.push_back(idcs[kMax]);
		this->projectOntoOrthogonalSubspace(acolsums, m, idcs, srtidcs, kMax);
		//std::cout << "SpVecApprox: Basis size = " << this->bvecs.size() << " Relative Error = " << sqrt(errSqNorm/vecSqNorm) << std::endl;
	}

	//final computation of coefficients
	MXd W = MXd::Zero(this->bvecs.size(), this->bvecs.size());
	MXd B = MXd::Zero(this->bvecs.size(), msave.cols());
	for(uint64_t i = 0; i < this->bvecs.size(); i++){
		for (uint64_t j = 0; j < this->bvecs.size(); j++){
			W(i, j) = msave(this->bvecs[i], this->bvecs[j]);
		}
		for (uint64_t j = 0; j < (uint64_t)msave.cols(); j++){
			B(i, j) = msave(this->bvecs[i], j);
		}
	}
	//for (uint64_t i =0 ; i < bvecs.size(); i++){
	//	std::cout << bvecs[i] << " ";
	//}
	//std::cout << std::endl;
	//std::cout << "W: " << std::endl << W << std::endl;
	//std::cout << "B: " << std::endl << B << std::endl;

	//std::cout << "SpVecApprox: Computing coefficients using LDLT decomposition" << std::endl;
	VXd finalcoeffs = W.ldlt().solve(B*acoeffs);
	for(uint64_t i = 0; i < (uint64_t)finalcoeffs.size(); i++){
		this->coeffs.push_back(finalcoeffs(i));
	}
	//std::cout << "SpVecApprox: Done!"<< std::endl;
}


void SparseVectorApproximation::projectOntoOrthogonalSubspace(VXd& acolsums, MXd& m, std::vector<uint64_t>& idcs, std::vector<uint64_t>& srtidcs, const uint64_t kMax){
	//rank one update to the column sums
	double ack = acolsums(kMax)/m(kMax, kMax);
	for(uint64_t i =0 ; i < (uint64_t)acolsums.size(); i++){
		acolsums(i) -= ack*m(kMax, i);
	}
	//rank one update to the gram matrix
	m -= (m.col(kMax)*m.row(kMax)/m(kMax, kMax)).eval(); //eval to prevent aliasing

	//remove the row/column corresponding to kMax
	removeRowCol(m, kMax);
	removeElem(acolsums, kMax);
	idcs.erase(idcs.begin()+kMax);
	for(uint64_t i =0 ; i < srtidcs.size(); i++){
		if (srtidcs[i] >= kMax){
			srtidcs[i]--;
		}
	}
	srtidcs.pop_back();
	return;
}

void SparseVectorApproximation::removeRowCol(MXd& m, uint64_t k){
	uint64_t nRows = m.rows();
	uint64_t nCols = m.cols();
	if (k < nRows-1){
		m.block(k, 0, nRows-1-k, nCols) = m.block(k+1, 0, nRows-1-k, nCols).eval();
		m.block(0, k, nRows, nCols-1-k) = m.block(0, k+1, nRows, nCols-1-k).eval();
	}
	m.conservativeResize(nRows-1, nCols-1);
	return;
}

void SparseVectorApproximation::removeElem(VXd& v, uint64_t k){
	uint64_t nRows = v.size();
	if (k < nRows-1){
		v.segment(k, nRows-1-k) = v.tail(nRows-1-k).eval();
	}
	v.conservativeResize(nRows-1);
	return;
}


void SparseVectorApproximation::getApprox(std::vector<uint64_t>& vecs, std::vector<double>& acoeffs){
	vecs = this->bvecs;
	acoeffs = this->coeffs;
	return;
}



#define __SPARSEVECTORAPPROX_IMPL_HPP
#endif /* __SPARSEAPPROX_IMPL_HPP */
