#!/bin/bash

##the whole script should take only a few min to run
#GOAL: get and compare taxonomy for deblurred, denoised data

##run script with i.e.
#nohup sh qiimescriptsGxGxE/03_taxonomy_Lemna_GxGxE.sh > taxonomy.out 2> taxonomy.err &

source /symbiont/apps/miniconda3/etc/profile.d/conda.sh
conda activate qiime2-2022.2

####GET CLASSIFIERS for deblurred & denoised data
##import greengenes data
qiime tools import --type 'FeatureData[Sequence]' \
    --input-path /symbiont/annamob/qiimeinputsGxGxE/gg_13_8_otus/rep_set/99_otus.fasta \
    --output-path /symbiont/annamob/qiimeoutputsGxGxE/otus99.qza
# input greengenes are from https://docs.qiime2.org/2021.4/data-resources/  downloaded november 2020 for analysis in manuscript.
# the zipped greengenes refence database include both of the inputs used in the above and below command:  99_otus.fasta &  99_otu_taxonomy.txt 
qiime tools import --type 'FeatureData[Taxonomy]' \
	--input-path /symbiont/annamob/qiimeinputsGxGxE/gg_13_8_otus/taxonomy/99_otu_taxonomy_header.txt \
	--output-path /symbiont/annamob/qiimeoutputsGxGxE/ref-taxonomy99.qza
## the 99_otu_taxonomy file from greengenes was modified to have a header row 'Feature ID', 'Taxon', tab separated, as this was required by qiime2


###trim it based on our reads, then ##fit the classifiers
#with the truncation performed in deblur, and on both types of deblur 
qiime feature-classifier extract-reads --i-sequences /symbiont/annamob/qiimeoutputsGxGxE/otus99.qza  \
	--p-f-primer  CCTACGGGNGGCWGCAG --p-r-primer  GACTACHVGGGTATCTAATCC \
	--p-trunc-len 402  --p-min-length 402 --o-reads /symbiont/annamob/qiimeoutputsGxGxE/ref-seqs99_402.qza
qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads /symbiont/annamob/qiimeoutputsGxGxE/ref-seqs99_402.qza \
	--i-reference-taxonomy /symbiont/annamob/qiimeoutputsGxGxE/ref-taxonomy99.qza --o-classifier /symbiont/annamob/qiimeoutputsGxGxE/classifier_402.qza
#with the truncation performed in denoising
qiime feature-classifier extract-reads --i-sequences /symbiont/annamob/qiimeoutputsGxGxE/otus99.qza  \
	--p-f-primer  CCTACGGGNGGCWGCAG --p-r-primer  GACTACHVGGGTATCTAATCC \
	--p-trunc-len 429  --p-min-length 260 --o-reads /symbiont/annamob/qiimeoutputsGxGxE/ref-seqs99_den.qza
qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads /symbiont/annamob/qiimeoutputsGxGxE/ref-seqs99_den.qza \
	--i-reference-taxonomy /symbiont/annamob/qiimeoutputsGxGxE/ref-taxonomy99.qza --o-classifier /symbiont/annamob/qiimeoutputsGxGxE/classifier_den.qza
# length range of 260 429 to reflect range for denoised sequences

####GET TAXONOMY CALLS IN DATASET
##identify the taxonomy in the dataset
#for deblurred, (min 10 seq across samples vs min 2)
qiime feature-classifier classify-sklearn --i-classifier /symbiont/annamob/qiimeoutputsGxGxE/classifier_402.qza \
	--i-reads /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_deblurrep-seqs_2.qza --p-confidence 0.7 \
	--o-classification /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_taxonomy_deb2.qza
qiime feature-classifier classify-sklearn --i-classifier /symbiont/annamob/qiimeoutputsGxGxE/classifier_402.qza \
	--i-reads /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_deblurrep-seqs.qza --p-confidence 0.7 \
	--o-classification /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_taxonomy_deb.qza
#for denoised
qiime feature-classifier classify-sklearn --i-classifier /symbiont/annamob/qiimeoutputsGxGxE/classifier_402.qza \
	--i-reads /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_denoiserep-seqs_2.qza --p-confidence 0.7 \
	--o-classification /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_taxonomy_den.qza

##get taxonomy visual, all three methods
qiime metadata tabulate --m-input-file /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_taxonomy_deb.qza \
	--o-visualization /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_taxonomy_deb.qzv
qiime metadata tabulate --m-input-file /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_taxonomy_deb2.qza \
	--o-visualization /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_taxonomy_deb2.qzv
qiime metadata tabulate --m-input-file /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_taxonomy_den.qza \
	--o-visualization /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_taxonomy_den.qzv


###FILTER SEQUENCES from groups, VISUALIZE taxonomy calling results.
##get barplots before filtering, just one example
qiime taxa barplot --i-table /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_deblurtable.qza \
	--i-taxonomy /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_taxonomy_deb.qza --m-metadata-file /symbiont/annamob/qiimeinputsGxGxE/MANIFEST.csv \
	--o-visualization /symbiont/annamob/qiimeoutputsGxGxE//Lemna_GxGxE_taxonomy_deb_taxabarplot_before_filteraxa.qzv

##remove land plant chloroplast and mitochondria -- some mitochondria may be insect or lab technician contamination; most mitochondria will likely be duckweed, some small fraction may be interesting in another analysis
##get barplots after filtering out chloroplast and mitochondria

#deblurred "normal" with 10 sequences across samples informing retention of ASV
qiime taxa filter-seqs --i-sequences /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_deblurrep-seqs.qza --i-taxonomy /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_taxonomy_deb.qza \
	--p-exclude Streptophyta,mitochondria --o-filtered-sequences /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_filteredtaxa-seqs_deb.qza
qiime taxa filter-table --i-table /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_deblurtable.qza --i-taxonomy /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_taxonomy_deb.qza \
	--p-exclude Streptophyta,mitochrondria --o-filtered-table /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_filteredtaxa-table_deb.qza
qiime taxa barplot --i-table /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_filteredtaxa-table_deb.qza --i-taxonomy /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_taxonomy_deb.qza \
	--m-metadata-file /symbiont/annamob/qiimeinputsGxGxE/MANIFEST.csv \
	--o-visualization /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_filteredtaxa-barplot_deb.qzv
##same for deblur with two sequences across samples informing retention of ASV
qiime taxa filter-seqs --i-sequences /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_deblurrep-seqs_2.qza --i-taxonomy /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_taxonomy_deb2.qza \
	--p-exclude Streptophyta,mitochondria --o-filtered-sequences /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_filteredtaxa-seqs_deb2.qza
qiime taxa filter-table --i-table /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_deblurtable_2.qza --i-taxonomy /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_taxonomy_deb2.qza \
	--p-exclude Streptophyta,mitochrondria --o-filtered-table /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_filteredtaxa-table_deb2.qza
qiime taxa barplot --i-table /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_filteredtaxa-table_deb2.qza --i-taxonomy /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_taxonomy_deb2.qza \
	--m-metadata-file /symbiont/annamob/qiimeinputsGxGxE/MANIFEST.csv \
	--o-visualization /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_filteredtaxa-barplot_deb2.qzv
#same for denoised data
qiime taxa filter-seqs --i-sequences /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_denoiserep-seqs_2.qza --i-taxonomy /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_taxonomy_den.qza \
	--p-exclude Streptophyta,mitochondria --o-filtered-sequences /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_filteredtaxa-seqs_den.qza
qiime taxa filter-table --i-table /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_denoisetable_2.qza --i-taxonomy /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_taxonomy_den.qza \
	--p-exclude Streptophyta,mitochrondria --o-filtered-table /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_filteredtaxa-table_den.qza
qiime taxa barplot --i-table /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_filteredtaxa-table_den.qza --i-taxonomy /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_taxonomy_den.qza \
	--m-metadata-file /symbiont/annamob/qiimeinputsGxGxE/MANIFEST.csv \
	--o-visualization /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_filteredtaxa-barplot_den.qzv


