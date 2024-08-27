#!/bin/bash
geneA="CD9"
geneB="NRG1"
breakpointA=6334700
breakpointB=32453346
TranscriptA="ENST00000382518"
TranscriptB="ENST00000356819"
strandA="+"
strandB="+"
BAM="/data/admsenzha/BAM/example.bam"
split=499
span=0

# plot read coverage using reads with markduplicates
sudo docker run --rm -v /data:/data senzhao/fusviz_shiny_app:1.7.0 \
  R -e "library(FuSViz); \
	options(uscsChromosomeName=FALSE); \
	pdf(file=paste('/data/admsenzha/', '$geneA', '-', '$geneB', '_', $split, '_', $span, '-markduplicates.pdf', sep=''), height=7, width=14); \
	plot_separate_individual_bam(first_name='$geneA', second_name='$geneB', breakpoint_A=$breakpointA, breakpoint_B=$breakpointB, transcriptA='$TranscriptA', transcriptB='$TranscriptB', coverage_plot_trans=F, version='hg19', properpair=NA, duplicate=F, unmappedmate=F, notpassQC=F, rna_bam_path='$BAM', split=$split, span=$span, fusion_strandA='$strandA', fusion_strandB='$strandB'); \
	dev.off();"

# plot read coverage using reads with unmarkduplicates
sudo docker run --rm -v /data:/data senzhao/fusviz_shiny_app:1.7.0 \
  R -e "library(FuSViz); \
	options(uscsChromosomeName=FALSE); \
	pdf(file=paste('/data/admsenzha/', '$geneA', '-', '$geneB', '_', $split, '_', $span, '-unmarkduplicates.pdf', sep=''), height=7, width=14); \
	plot_separate_individual_bam(first_name='$geneA', second_name='$geneB', breakpoint_A=$breakpointA, breakpoint_B=$breakpointB, transcriptA='$TranscriptA', transcriptB='$TranscriptB', coverage_plot_trans=F, version='hg19', properpair=NA, duplicate=T, unmappedmate=F, notpassQC=F, rna_bam_path='$BAM', split=$split, span=$span, fusion_strandA='$strandA', fusion_strandB='$strandB'); \
	dev.off();"

