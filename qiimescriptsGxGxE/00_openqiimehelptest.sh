#!/bin/bash

#start qiime and open help

#the following must be run each time to start qiime on the symbiont cluster.
#other computing scenarios will vary!
source /symbiont/apps/miniconda3/etc/profile.d/conda.sh
conda activate qiime2-2022.2

#qiime commands, these will be the same across clusters, except for the filepaths
qiime --help 
