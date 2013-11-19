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

2. Install the liblpsolve dependency:
	<pre>
	sudo apt-get install liblpsolve55-dev
	</pre>

3. Navigate to the directory and run the install script:
	<pre>
	sudo ./install
	</pre>

4. To test out the code, check out the `examples` folder:
	<pre>
	make config=release
	./DynMeansExample
	</pre>

	Also included is a premake4.lua script, in case you have premake and want to generate the makefiles for yourself

####Citation

If you use Dynamic Means for a paper or project, please use the following BibTeX entry for citation:
	<pre>
    @inproceedings{Campbell13_NIPS,
    	Author = {Trevor Campbell and Miao Liu and Brian Kulis and Jonathan P.~How and Lawrence Carin},
    	Title = {Dynamic Clustering via Asymptotics of the Dependent Dirichlet Process Mixture},
    	Year = {2013},
    	Booktitle = {Neural Information Processing Systems (NIPS)}}
   	</pre>


