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
Kernelized Dynamic Means is another extension which also operates on similarity graphs, but does not require
the computationally expensive step of computing eigenvectors that the spectral code performs.
Efficient C++ implementations of all algorithms are provided in this repository.

<p align="center">
<img src="https://github.com/trevorcampbell/dynamic-means/blob/master/imgs/clustermotion.png?raw=true"/>
</p>


####Usage
1. Clone this repository:
	<pre>
    git clone https://github.com/trevorcampbell/dynamic-means.git
    </pre>

2. Navigate to the directory and run the install script:
	<pre>
	sudo ./install
	</pre>
3. In your code, include one or more of the following headers,
	<pre>
	#include &lt;dynmeans/dynmeans.hpp>
	#include &lt;dynmeans/specdynmeans.hpp>
	#include &lt;dynmeans/kerndynmeans.hpp>
	</pre>
    depending on whether you want to use regular Dynamic Means,
    Spectral Dynamic Means, or 
    Kernelized Dynamic Means. Make sure [Gurobi](www.gurobi.com) is installed (free for academic use) if 
   you want to use Spectral/Kernel Dynamic Means. If you can't get access to Gurobi, feel free to modify
   the `SpecDynMeans::getOldNewMatching` function in `src/specdynmeans_impl.hpp` and
   the `KernDynMeans::getMinWtMatching` function in `src/kerndynmeans_impl.hpp` to use a different
   LP solver.
4. Create a DynMeans, KernDynMeans, and/or SpecDynMeans object:
	<pre>
	double lambda = .05;
	double T_Q = 6.8;
	double K_tau = 1.01;
	double Q = lambda/T_Q;
	double tau = (T_Q*(K_tau - 1.0)+1.0)/(T_Q-1.0);
	DynMeans&lt;Eigen::Vector2d> dynm(lambda, Q, tau);
	SpecDynMeans&lt;YourAffType> sdynm(lambda, Q, tau);
	KernDynMeans&lt;YourAffType> kdynm(lambda, Q, tau);
	</pre>
	where Eigen::Vector2d is the vector data type that Dynamic Means is going to cluster.
	Note that other types can be used in place of Eigen::Vector2d, but they must
	implement vector addition, and scalar multiplcation/division. For Spectral/Kernel Dynamic Means,
	YourAffType is a wrapper you must write to abstract the computation of node->node and node->cluster affinities. To 
	find out which functions YourAffType must implement, see the example in examples/mainkdm.cpp. 
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
	
	To cluster the first window of data with Spectral/Kernel Dynamic Means, just call the `SpecDynMeans::cluster`/`KernDynMeans::cluster` function
	<pre>
	vector&lt;Eigen::Vector2d> dataWindow1;
	int nRestarts = 10;
	int nCoarsest = 20;
	int nClusMax = 5;
	...
	YourAffinityWrapper aff;
	aff.yourDataUpdateFunction(dataWindow1);
	vector&lt;int> learnedLabels1, prmlbls1;
	vector&lt;double> gammas;
	double obj1, tTaken1;
	//for spectral dynamic means, call sdynm.cluster
	sdynm.cluster(affinities, nRestarts, nClusMax, SpecDynMeans<YourAffinityWrapper>::EigenSolverType::REDSVD learnedLabels1, obj1, gammas1, prmlbls1, tTaken1);
	//for kernel dynamic means, call kdynm.cluster
	kdynm.cluster(affinities, nRestarts, nCoarsest, learnedLabels1, obj1, gammas1, prmlbls1, tTaken1);
	//when done, make sure to update the affinities wrapper with the new clustering
	affinities.yourUpdateFunction(learnedLabels, gammas, prmlbls);
	</pre>
	where `obj1` is the clustering cost output, `tTaken1` is the clustering time output, 
	`learnedLabels1` is the data labels output, `gammas1` is the old cluster weights output, and `prmlbls1` 
	is the old parameter labels output. 
	`nRestarts` is the number of random orthogonal matrix initializations Spectral Dynamic Means will try,
	and `nClusMax` is (intuitively) the maximum number of new clusters expected in each timestep (mathematically,
	it is the rank approximation to use when doing eigendecompositions). `EigenSolverType::REDSVD` tells
	the algorithm to use an approximate eigendecomposition (adapted from [redsvd](https://code.google.com/p/redsvd/)).
	To reduce sensitivity to initialization, 
	Kernel Dynamic Means iteratively coarsifies the graph until it has `nCoarsest` nodes, applies a spectral clustering to that coarse graph,
	and then iteratively refines the graph and applies k-means-like labelling updates. Pick `nCoarsest` to be a number that is small enough
	such that computing the eigenvectors of an `nCoarsest` by `nCoarsest` matrix is fast enough for your application.

6. To cluster another window of data, just call `DynMeans::cluster`/`SpecDynMeans::cluster`/`KernDynMeans::Cluster` again
	<pre>
	vector&lt;Eigen::Vector2d> dataWindow2;
	...
	vector&lt;Eigen::Vector2d> learnedParams2;
	vector&lt;int> learnedLabels2;
	double obj2, tTaken2;
	dynm.cluster(dataWindow2, nRestarts, learnedLabels2, learnedParams2, obj2, tTaken2);
	</pre>
	and make sure to use `YourAffinityWrapper::yourUpdateFunction` afterwards to compute the new set of "old parameter nodes"
	to prepare Spectral/Kernel Dynamic Means for the next clustering step

7. Repeat step 6 as many times as required (e.g., split a dataset of 1,000,000 datapoints into chunks of 1,000 and cluster each sequentially) 

####Example Code
To run the example, first make sure liblpsolve is installed (required for label accuracy computations):
    
    sudo apt-get install liblpsolve55-dev

make sure liblpsolve55.so (not liblpsolve55.a) is located in your /usr/lib/ folder when the installation is complete, and if it is not, move it to that location. If you compile against the static library (.a), you will
get linker errors about undefined references to dlclose, dlopen, colamd_*, etc.
   

Navigate to the `examples` folder to compile and run the example for Dynamic Means:
    
    make config=release DynMeansExample
    ./DynMeansExample

For Spectral Dynamic Means, run

    make config=release SpecDynMeansExample
    ./SpecDynMeansExample
 
For Kernel Dynamic Means, run

    make config=release KernDynMeansExample
    ./KernDynMeansExample

If you want to change how the example compiles, a [premake](http://industriousone.com/premake) 
Makefile generation script is included.

####Citation

If you use Dynamic Means/Spectral Dynamic Means/Kernel Dynamic Means for a paper or project, please use the following BibTeX entry for citation:
	<pre>
    @inproceedings{Campbell13_NIPS,
    	Author = {Trevor Campbell and Miao Liu and Brian Kulis and Jonathan P.~How and Lawrence Carin},
    	Title = {Dynamic Clustering via Asymptotics of the Dependent Dirichlet Process Mixture},
    	Year = {2013},
    	Booktitle = {Advances in Neural Information Processing Systems 26}}
   	</pre>


