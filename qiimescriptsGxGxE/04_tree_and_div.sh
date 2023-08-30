#!/bin/bash


##GOAL of this script is to build a phylogenetic tree of ASVs in the dataset and calculate standard diversity metrics within and between samples
## the first command uses sepp-refs-gg-13-8.qza, which is a database downloaded from https://docs.qiime2.org/2021.4/data-resources/ November, 2020 

##run script with i.e.
#nohup sh qiimescriptsGxGxE/04_tree_and_div.sh > treeanddiv.out 2> treeanddiv.err &

source /symbiont/apps/miniconda3/etc/profile.d/conda.sh
conda activate qiime2-2022.2

###the first command gets phylogeny info; note this depends on the header to set threads, a user may need to alter this.
###the second command cleans features from table that can't be placed in phylogeny
###the third and fourth commands calculate basic set of diversity metrics
###note that alpha diversity metrics may not be ideal for 16s sequence data, see e.g. https://doi.org/10.1016/j.csbj.2021.12.036
###the fifth and sixth commands subset to inocula-only data, followed by a seventh command to re-compute community distance diversity metrics

#deblurred data, with 10 sequences across samples as a minimum to consider as ASV
qiime fragment-insertion sepp --p-threads 8 \
	 --i-reference-database /symbiont/annamob/qiimeinputsGxGxE/sepp-refs-gg-13-8.qza \
     --i-representative-sequences /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_filteredtaxa-seqs_deb.qza \
     --o-tree /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_sepptree_deb.qza \
     --o-placements /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_seppplace_deb.qza
qiime fragment-insertion filter-features --verbose \
	--i-table /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_filteredtaxa-table_deb.qza \
	--i-tree /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_sepptree_deb.qza \
	--o-filtered-table /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_sepp_keptfeatures_table_deb.qza \
	--o-removed-table /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_sepp_rmdfeatures_table_deb.qza
qiime diversity core-metrics-phylogenetic --i-table /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_sepp_keptfeatures_table_deb.qza \
	--i-phylogeny /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_sepptree_deb.qza --p-sampling-depth 1000 \
	--m-metadata-file /symbiont/annamob/qiimeinputsGxGxE/Lemna_microbiome_barebones_meta_wkmcC_ladptFM.tsv \
	--output-dir /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_coremetrics_deb
qiime diversity alpha-rarefaction  --i-table /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_sepp_keptfeatures_table_deb.qza \
	--i-phylogeny /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_sepptree_deb.qza  --p-max-depth 5000 \
	--m-metadata-file /symbiont/annamob/qiimeinputsGxGxE/Lemna_microbiome_barebones_meta_wkmcC_ladptFM.tsv \
	--o-visualization /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_alpha-rarefaction_deb.qzv
qiime feature-table filter-samples --i-table /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_sepp_keptfeatures_table_deb.qza \
	--m-metadata-file /symbiont/annamob/qiimeinputsGxGxE/Lemna_microbiome_barebones_meta_wkmcC_ladptFM.tsv --p-where "fieldormaster='M'" \
	--o-filtered-table /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_sepp_keptfeatures_table_debMa.qza
qiime feature-table filter-features --i-table /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_sepp_keptfeatures_table_debMa.qza \
	--p-min-frequency 1 --o-filtered-table /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_sepp_keptfeatures_table_debM.qza
qiime diversity core-metrics-phylogenetic --i-table /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_sepp_keptfeatures_table_debM.qza \
	--i-phylogeny /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_sepptree_deb.qza --p-sampling-depth 1000 \
	--m-metadata-file /symbiont/annamob/qiimeinputsGxGxE/Lemna_microbiome_barebones_meta_wkmcC_ladptM.tsv \
	--output-dir /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_coremetrics_deb_M



#deblurred data, with 2 sequences across samples as a minimum to consider as ASV
qiime fragment-insertion sepp --p-threads 8 \
	--i-reference-database /symbiont/annamob/qiimeinputsGxGxE/sepp-refs-gg-13-8.qza \
     --i-representative-sequences /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_filteredtaxa-seqs_deb2.qza \
     --o-tree /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_sepptree_deb2.qza \
     --o-placements /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_seppplace_deb2.qza
qiime fragment-insertion filter-features --verbose \
	--i-table /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_filteredtaxa-table_deb2.qza \
	--i-tree /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_sepptree_deb2.qza \
	--o-filtered-table /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_sepp_keptfeatures_table_deb2.qza \
	--o-removed-table /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_sepp_rmdfeatures_table_deb2.qza
qiime diversity core-metrics-phylogenetic --i-table /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_sepp_keptfeatures_table_deb2.qza \
	--i-phylogeny /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_sepptree_deb2.qza --p-sampling-depth 1000 \
	--m-metadata-file /symbiont/annamob/qiimeinputsGxGxE/Lemna_microbiome_barebones_meta_wkmcC_ladptFM.tsv \
	--output-dir /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_coremetrics_deb2
qiime diversity alpha-rarefaction  --i-table /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_sepp_keptfeatures_table_deb2.qza \
	--i-phylogeny /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_sepptree_deb2.qza  --p-max-depth 5000 \
	--m-metadata-file /symbiont/annamob/qiimeinputsGxGxE/Lemna_microbiome_barebones_meta_wkmcC_ladptFM.tsv \
	--o-visualization /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_alpha-rarefaction_deb2.qzv
qiime feature-table filter-samples --i-table /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_sepp_keptfeatures_table_deb2.qza \
	--m-metadata-file /symbiont/annamob/qiimeinputsGxGxE/Lemna_microbiome_barebones_meta_wkmcC_ladptFM.tsv --p-where "fieldormaster='M'" \
	--o-filtered-table /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_sepp_keptfeatures_table_deb2Ma.qza
qiime feature-table filter-features --i-table /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_sepp_keptfeatures_table_deb2Ma.qza \
	--p-min-frequency 1 --o-filtered-table /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_sepp_keptfeatures_table_deb2M.qza
qiime diversity core-metrics-phylogenetic --i-table /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_sepp_keptfeatures_table_deb2M.qza \
	--i-phylogeny /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_sepptree_deb2.qza --p-sampling-depth 1000 \
	--m-metadata-file /symbiont/annamob/qiimeinputsGxGxE/Lemna_microbiome_barebones_meta_wkmcC_ladptM.tsv \
	--output-dir /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_coremetrics_deb2_M


#denoised data
qiime fragment-insertion sepp --p-threads 8 \
	--i-reference-database /symbiont/annamob/qiimeinputsGxGxE/sepp-refs-gg-13-8.qza \
     --i-representative-sequences /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_filteredtaxa-seqs_den.qza \
     --o-tree /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_sepptree_den.qza \
     --o-placements /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_seppplace_den.qza
qiime fragment-insertion filter-features --verbose \
	--i-table /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_filteredtaxa-table_den.qza \
	--i-tree /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_sepptree_den.qza \
	--o-filtered-table /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_sepp_keptfeatures_table_den.qza \
	--o-removed-table /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_sepp_rmdfeatures_table_den.qza 
qiime diversity core-metrics-phylogenetic --i-table /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_sepp_keptfeatures_table_den.qza \
	--i-phylogeny /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_sepptree_den.qza --p-sampling-depth 1000 \
	--m-metadata-file /symbiont/annamob/qiimeinputsGxGxE/Lemna_microbiome_barebones_meta_wkmcC_ladptFM.tsv \
	--output-dir /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_coremetrics_den
qiime diversity alpha-rarefaction  --i-table /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_sepp_keptfeatures_table_den.qza \
	--i-phylogeny /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_sepptree_den.qza  --p-max-depth 5000 \
	--m-metadata-file /symbiont/annamob/qiimeinputsGxGxE/Lemna_microbiome_barebones_meta_wkmcC_ladptFM.tsv \
	--o-visualization /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_alpha-rarefaction_den.qzv
qiime feature-table filter-samples --i-table /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_sepp_keptfeatures_table_den.qza \
	--m-metadata-file /symbiont/annamob/qiimeinputsGxGxE/Lemna_microbiome_barebones_meta_wkmcC_ladptFM.tsv --p-where "fieldormaster='M'" \
	--o-filtered-table /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_sepp_keptfeatures_table_denMa.qza
qiime feature-table filter-features --i-table /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_sepp_keptfeatures_table_denMa.qza \
	--p-min-frequency 1 --o-filtered-table /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_sepp_keptfeatures_table_denM.qza
qiime diversity core-metrics-phylogenetic --i-table /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_sepp_keptfeatures_table_denM.qza \
	--i-phylogeny /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_sepptree_den.qza --p-sampling-depth 1000 \
	--m-metadata-file /symbiont/annamob/qiimeinputsGxGxE/Lemna_microbiome_barebones_meta_wkmcC_ladptM.tsv \
	--output-dir /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_coremetrics_den_M

	