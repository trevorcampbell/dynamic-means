dynamic-means
=============

###Clustering for Temporally Evolving Datasets

####Introduction

The [Dynamic Means algorithm](http://arxiv.org/abs/1305.6659) is a k-means-like algorithm for clustering large or temporally evolving datasets.
It operates batch-sequentially, i.e. it processes "windows" of data rather than processing an entire dataset
at once. This allows it to capture clusters that change over time (via (a) motion, (b) creation, and (c) deletion), or to
cluster large datasets by considering small chunks of data at a time. Spectral Dynamic Means is an extension 
to Dynamic Means that captures cluster motion, creation and deletion in general similarity graphs -- this allows
it to capture more general data types and cluster shapes, at the expense of increased computational cost.
Efficient C++ implementations of both algorithms are provided in this repository.

<p align="center">
<img src="https://github.com/tc4mp/dynamic-means/blob/master/imgs/clustermotion.png?raw=true"/>
</p>


####Usage
1. Clone this repository:
	<pre>
    git clone http://github.com/tc4mp/dynamic-means
    </pre>

2. Navigate to the directory and run the install script:
	<pre>
	sudo ./install
	</pre>
3. In your code, include the header:
	<pre>
	#include &lt;dynmeans/dynmeans.hpp>
	</pre>
   for Dynamic Means, and:
	<pre>
	#include &lt;specdynmeans/specdynmeans.hpp>
	</pre>
   for Spectral Dynamic Means. Make sure [Gurobi](www.gurobi.com) is installed (free for academic use) if 
   you want to use Spectral Dynamic Means. If you can't get access to Gurobi, feel free to modify
   the `SpecDynMeans::getOldNewMatching` function in `src/specdynmeans_impl.hpp` to use a different
   LP solver.
4. Create a DynMeans and/or SpecDynMeans object:
	<pre>
	double lambda = .05;
	double T_Q = 6.8;
	double K_tau = 1.01;
	double Q = lambda/T_Q;
	double tau = (T_Q*(K_tau - 1.0)+1.0)/(T_Q-1.0);
	DynMeans&lt;Eigen::Vector2d> dynm(lambda, Q, tau);
	SpecDynMeans&lt;DatType, PrmType> sdynm(lambda, Q, tau);
	</pre>
	where Eigen::Vector2d is the vector data type that Dynamic Means is going to cluster.
	Note that other types can be used in place of Eigen::Vector2d, but they must
	implement vector addition, and scalar multiplcation/division. For Spectral Dynamic Means,
	DatType and PrmType are the data and parameter types, respectively. To see which functions
	DatType and PrmType must implement, see the example in examples/mainsdm.cpp.
	See [the Dynamic Means paper](http://arxiv.org/abs/1305.6659) for a description
	of the values `lambda, T_Q, K_tau, Q, tau`.

5. To cluster the first window of data with Dynamic Means, just call the `DynMeans::cluster` function
	<pre>
	vector&lt;Eigen::Vector2d> dataWindow1;
	...
	int nRestarts = 10;
	vector&lt;Eigen::Vector2d> learnedParams1;
	vector&lt;int> learnedLabels1;
	double obj1, tTaken1;
	dynm.cluster(dataWindow1, nRestarts, learnedLabels1, learnedParams1, obj1, tTaken1);
	</pre>
	where `obj1` is the clustering cost output, `tTaken1` is the clustering time output, 
	`learnedLabels1` is the data labels output, and `learnedParams1` is the cluster parameters output.
	`nRestarts` is the number of random label assignment orders Dynamic Means will try.
	
	To cluster the first window of data with Spectral Dynamic Means, just call the `SpecDynMeans::cluster` function
	<pre>
	vector&lt;Eigen::Vector2d> dataWindow1;
	...
	int nRestarts = 10;
	vector&lt;Eigen::Vector2d> learnedParams1;
	vector&lt;int> learnedLabels1;
	double obj1, tTaken1;
	sdynm.cluster(dataWindow1, nRestarts, nClusMax, SpecDynMeans<DatType,PrmType>::EigenSolverType::REDSVD, learnedLabels1, obj1, tTaken1);
	</pre>
	where `obj1` is the clustering cost output, `tTaken1` is the clustering time output, 
	`learnedLabels1` is the data labels output, and `learnedParams1` is the cluster parameters output.
	`nRestarts` is the number of random orthogonal matrix initializations Spectral Dynamic Means will try,
	and `nClusMax` is (intuitively) the maximum number of new clusters expected in each timestep (mathematically,
	it is the rank approximation to use when doing eigendecompositions). `EigenSolverType::REDSVD` tells
	the algorithm to use an approximate eigendecomposition (adapted from [redsvd](https://code.google.com/p/redsvd/)).

6. To cluster another window of data, just call `DynMeans::cluster` and/or `SpecDynMeans::cluster` again
	<pre>
	vector&lt;Eigen::Vector2d> dataWindow2;
	...
	vector&lt;Eigen::Vector2d> learnedParams2;
	vector&lt;int> learnedLabels2;
	double obj2, tTaken2;
	dynm.cluster(dataWindow2, nRestarts, learnedLabels2, learnedParams2, obj2, tTaken2);
	</pre>

7. Repeat step 6 as many times as required (e.g., split a dataset of 1,000,000 datapoints into chunks of 1,000 and call `DynMeans::cluster` and/or `SpecDynMeans::cluster` on each)

####Example Code
To run the example, first make sure liblpsolve is installed (required for label accuracy computations):
    
    sudo apt-get install liblpsolve55-dev
   


Navigate to the `examples` folder to compile and run the example for Dynamic Means:
    
    make config=release DynMeansExample
    ./DynMeansExample

Or for Spectral Dynamic Means, run

    make config=release SpecDynMeansExample
    ./SpecDynMeansExample
 

If you want to change how the example compiles, a [premake](http://industriousone.com/premake) 
Makefile generation script is included.

####Citation

If you use Dynamic Means or Spectral Dynamic Means for a paper or project, please use the following BibTeX entry for citation:
	<pre>
    @inproceedings{Campbell13_NIPS,
    	Author = {Trevor Campbell and Miao Liu and Brian Kulis and Jonathan P.~How and Lawrence Carin},
    	Title = {Dynamic Clustering via Asymptotics of the Dependent Dirichlet Process Mixture},
    	Year = {2013},
    	Booktitle = {Advances in Neural Information Processing Systems 26}}
   	</pre>


