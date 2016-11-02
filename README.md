# GOS-SynCut
Code Accompanying the paper "The Geometry of Synchronization Problems and Learning Group Actions"

#### Brief Description of Each Script
+ `batchSynCut.m` Batch script for repeated experiments and statistics charts.
+ `demoMultiWaySynCutLemurTeeth.m` Demo script for SynCut on the 50 lemur teeth dataset. Default number of clusters is 3, but the code runs with general number of clusters.
+ `demoTwoWaySynCut.m` Demo script for SynCut with two clusters, to illustrate the ideas spectral clustering with edge-wise frustrations as weights.
+ `demoTwoWaySynCut_SBM.m` Demo script for SynCut with two clusters, same functionality as `demoTwoWaySynCut.m` except for using [Stochastic Block Model](https://en.wikipedia.org/wiki/Stochastic_block_model) to generate the random graph.


