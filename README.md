dynamic-means
=============

###Clustering for Temporally Evolving Datasets

####Introduction

The [Dynamic Means algorithm](http://arxiv.org/abs/1305.6659) is a k-means-like algorithm for clustering large or temporally evolving datasets.
It operates batch-sequentially, i.e. it processes "windows" of data rather than processing an entire dataset
at once. This allows it to capture clusters that change over time (via (a) motion, (b) creation, and (c) deletion), or to
cluster large datasets by considering small chunks of data at a time.

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
4. Create a DynMeans object:
	<pre>
	double lambda = .05;
	double T_Q = 6.8;
	double K_tau = 1.01;
	double Q = lambda/T_Q;
	double tau = (T_Q*(K_tau - 1.0)+1.0)/(T_Q-1.0);
	DynMeans&lt;Eigen::Vector2d> dynm(lambda, Q, tau);
	</pre>
	where Eigen::Vector2d is the vector data type that Dynamic Means is going to cluster.
	Note that other types can be used in place of Eigen::Vector2d, but they must
	implement vector addition, and scalar multiplcation/division.
	See [the Dynamic Means paper](http://arxiv.org/abs/1305.6659) for a description
	of the values `lambda, T_Q, K_tau, Q, tau`.

5. To cluster the first window of data, just call the `DynMeans::cluster` function
	<pre>
	vector&lt;Eigen::Vector2d> dataWindow1;
	...
	int nRestarts = 10;
	dynm.cluster(dataWindow1, nRestarts);
	</pre>

6. To get results from the latest clustering step, call `DynMeans::getClustering`
	<pre>
	vector&lt;int> labels1;
	vector&lt;Eigen::Vector2d> parameters1;
	dynm.getClustering(labels1, parameters1);
	</pre>

7. To cluster another window of data, just call `DynMeans::cluster` and `DynMeans::getClustering` again
	<pre>
	vector&lt;Eigen::Vector2d> dataWindow2;
	...
	dynm.cluster(dataWindow2, nRestarts);
	vector&lt;int> labels2;
	vector&lt;Eigen::Vector2d> parameters2;
	dynm.getClustering(labels2, parameters2);
	</pre>

8. Repeat step 7 as many times as required (e.g., split a dataset of 1,000,000 datapoints into chunks of 1,000 and run `DynMeans::cluster` on each)

7. (optional) To run the example, first make sure liblpsolve is installed (required for label accuracy computations):
	<pre>
	sudo apt-get install liblpsolve55-dev
	</pre>
	and navigate to the `examples` folder to compile and run the example:
	<pre>
	make config=release
	./DynMeansExample
	</pre>

####Citation

If you use Dynamic Means for a paper or project, please use the following BibTeX entry for citation:
	<pre>
    @inproceedings{Campbell13_NIPS,
    	Author = {Trevor Campbell and Miao Liu and Brian Kulis and Jonathan P.~How and Lawrence Carin},
    	Title = {Dynamic Clustering via Asymptotics of the Dependent Dirichlet Process Mixture},
    	Year = {2013},
    	Booktitle = {Neural Information Processing Systems (NIPS)}}
   	</pre>


