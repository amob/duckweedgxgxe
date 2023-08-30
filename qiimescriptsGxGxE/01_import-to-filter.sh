#!/bin/bash


##GOAL: summarize demultiplexed reads in the GxGxE project, trim adapters, join pairs, quality filter and re-summarize
#quality filtering should be run before deblur, which is in the next step.


##run script with i.e.
#nohup sh qiimescriptsGxGxE/01_import-to-filter.sh > import-to-filter.out 2> import-to-filter.err &

source /symbiont/apps/miniconda3/etc/profile.d/conda.sh
conda activate qiime2-2022.2

qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path /symbiont/annamob/qiimeinputsGxGxE/MANIFEST.csv \
  --input-format PairedEndFastqManifestPhred33 \
  --output-path /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_imported_demux.qza
#import data using sequencing files and a manifest file

qiime demux summarize  --i-data /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_imported_demux.qza \
   --o-visualization /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_demux.qzv
#outputs file Lemna_microbiome_view_demux.qzv, a qiime readable object, all .qzv can be visualizd on the https://view.qiime2.org website
#.qza objects can be uploaded to the site as well to view information on provenance.

qiime cutadapt trim-paired --i-demultiplexed-sequences /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_imported_demux.qza \
   --p-cores 40 --p-front-f CCTACGGGNGGCWGCAG --p-front-r GACTACHVGGGTATCTAATCC \
   --o-trimmed-sequences /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_imported_demux_trim.qza 
#this is common region of the primer used -- the variable region comes "before" (in the reading direction) these start

#merging paired ends first is required, quality filtering is reccommended by deblur authors
##merge
qiime vsearch join-pairs --i-demultiplexed-seqs /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_imported_demux_trim.qza \
          --p-minovlen 15 --p-maxdiffs 10 --o-joined-sequences /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_join.qza

##quality filter
qiime quality-filter q-score  --i-demux /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_join.qza \
        --o-filtered-sequences /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_filter.qza \
        --o-filter-stats /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_filter-stats.qza

##demultiplex summmarize to see what truncation should be
qiime demux summarize --i-data /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_filter.qza \
           --o-visualization /symbiont/annamob/qiimeoutputsGxGxE/Lemna_GxGxE_filter.qzv
###based on this visualization, truncation was set at 402 (see summary table, can retain 91% of sequences with this length)
##read quality is good for longer, but deblur and taxonomy assignment require a fixed length of sequence across all reads, and many reads are lost after this length
###TRUNCATION HAPPENS DURING DEBLUR, this number is determined here, but used in the next script.