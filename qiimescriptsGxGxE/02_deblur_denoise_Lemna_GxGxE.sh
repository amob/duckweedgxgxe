#!/bin/bash

#GOAL: run deblur on trimmed, merged, and quality filtered sequences

#NOTE: deblur takes more time for more samples. For this dataset, on the [cluster] machine it ran on, only 1 hour, maybe less.
#if re-running analysis, double check results of previous script in qiime2view  before setting truncate length in the below

##run script with i.e.
#nohup sh qiimescriptsGxGxE/02_deblur_denoise_Lemna_GxGxE.sh > deb_den.out 2> deb_den.err &

source /symbiont/apps/miniconda3/etc/profile.d/conda.sh
conda activate qiime2-2022.2


#run deblur

#run slower with 1 job (ex. this one), or run faster with 8 jobs (example next)
qiime deblur denoise-16S --i-demultiplexed-seqs /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_filter.qza \
        --p-trim-length 402 --p-sample-stats --p-min-reads 10 --p-jobs-to-start 1 \
        --o-table /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_deblurtable.qza \
        --o-representative-sequences /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_deblurrep-seqs.qza \
        --o-stats /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_deblur-stats.qza

qiime deblur denoise-16S --i-demultiplexed-seqs /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_filter.qza \
        --p-trim-length 402 --p-sample-stats --p-min-reads 2 --p-jobs-to-start 8 \
        --o-table /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_deblurtable_2.qza \
        --o-representative-sequences /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_deblurrep-seqs_2.qza \
        --o-stats /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_deblur-stats_2.qza



#summarize and visualize each output deblur file, first for deblur with 10 min then 2 min.
qiime feature-table summarize --i-table /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_deblurtable.qza \
        --o-visualization /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_deblurtable.qzv

qiime feature-table tabulate-seqs --i-data /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_deblurrep-seqs.qza \
        --o-visualization /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_deblurrep-seqs.qzv

qiime deblur visualize-stats --i-deblur-stats /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_deblur-stats.qza \
        --o-visualization /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_deblur-stats.qzv

qiime feature-table summarize --i-table /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_deblurtable_2.qza \
        --o-visualization /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_deblurtable_2.qzv

qiime feature-table tabulate-seqs --i-data /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_deblurrep-seqs_2.qza \
        --o-visualization /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_deblurrep-seqs_2.qzv

qiime deblur visualize-stats --i-deblur-stats /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_deblur-stats_2.qza \
        --o-visualization /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_deblur-stats_2.qzv
        
#trunc length 201 would be most similar to the deblur option? reads are ~ 250, so this is essentially no truncation
#run denoise
qiime dada2 denoise-paired --i-demultiplexed-seqs /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_imported_demux_trim.qza \
		--p-trunc-len-f 0 --p-trunc-len-r 0 --p-pooling-method 'pseudo' --p-n-threads 8 \
        --o-table /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_denoisetable_2.qza \
        --o-representative-sequences /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_denoiserep-seqs_2.qza \
        --o-denoising-stats /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_denoise-stats_2.qza


qiime feature-table summarize --i-table /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_denoisetable_2.qza \
        --o-visualization /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_denoisetable_2.qzv

qiime feature-table tabulate-seqs --i-data /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_denoiserep-seqs_2.qza \
        --o-visualization /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_denoiserep-seqs_2.qzv

