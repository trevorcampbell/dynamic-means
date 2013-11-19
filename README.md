dynamic-means
=============

###Clustering for Temporally Evolving Datasets

####Introduction

The Dynamic Means algorithm is a k-means-like algorithm for clustering large or temporally evolving datasets.
It operates batch-sequentially, i.e. it processes "windows" of data rather than processing an entire dataset
at once. This allows it to capture clusters that change over time (via (a) motion, (b) creation, and (c) deletion), or to
cluster large datasets by considering small chunks of it at a time.

![](https://github.com/tc4mp/dynamic-means/blob/master/imgs/clustermotion.png?raw=true)

####Citation

If you use Dynamic Means for a paper or project, please use the following bibtex entry for citation:

    @inproceedings{Campbell13_NIPS,
    	Author = {Trevor Campbell and Miao Liu and Brian Kulis and Jonathan P.~How and Lawrence Carin},
    	Title = {Dynamic Clustering via Asymptotics of the Dependent Dirichlet Process Mixture},
    	Year = {2013},
    	Booktitle = {Neural Information Processing Systems (NIPS)}}


