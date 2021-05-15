#' Plot straight lines to link breakpoints
#'
#' @description Plot straight lines to link breakpoints of geneA and geneB transcripts (reducing introns) in 'domain' tab-panel
#'
#' @param upstream_xy A data.frame (i.e. a returned object \code{'A1_xy'} in FUNCTION "plot_separate_domain_geneA" - exon coordinates(x1, y1, x2, y2) of intron-reduced transcripts).
#' @param upstream_flag A list (i.e. the same object \code{'first'} used in FUNCTION "plot_separate_domain_geneA")
#'        with two elements (e.g. \code{'first\$pos\$transcript'} and \code{'first\$pos\$select_region'}, in which introns were reduced in the transcript).
#' @param downstream_xy A data.frame (i.e. a returned object \code{'B1_xy'} in FUNCTION "plot_separate_domain_geneB" - exon coordinates(x1, y1, x2, y2) of intron-reduced transcripts).
#' @param downstream_flag A list (i.e. the same object \code{'second'} used in FUNCTION "plot_separate_domain_geneB")
#'        with two elements (e.g. \code{'second\$pos\$transcript'} and \code{'second\$pos\$select_region'}, in which introns were reduced in the transcript).
#' @param breakpoint_xy A data.frame with one row (e.g. 'breakpoint_xy[i,1]' - breakpoint pos in geneA transcript; 'breakpoint_xy[i,2]' - breakpoint pos in geneB transcript).
#'        NOTE: breakpoint coordinates are adjusted by reducing introns.
#'
#' @return A data.frame with five columns (i.e. \code{'x_pos_gene_upstream'}, \code{'x_pos_gene_downstream'}, \code{'y_pos_gene_upstream'}, \code{'y_pos_gene_downstream'}, \code{'tag'}).
#'         NOTE: \code{'tag'} suggests different biological consequences of fusion with four categories (e.g. 1:Unknown, 2:Inframe, 3:Outframe, 4:Truncate-loss)
#'
.plot_curve_reduce <- function(upstream_xy, upstream_flag, downstream_xy, downstream_flag, breakpoint_xy) {
	# For testing: upstream_xy=A1_xy; upstream_flag=first; downstream_xy=B1_xy; downstream_flag=second
	break_A=breakpoint_xy[, 1];
	break_B=breakpoint_xy[, 2];
	tag = 0;

	exon_xy_A1 = Gviz::coords(upstream_xy$GeneRegionTrack); #// coordinates of all exons in geneA plotting
	exon_xy_B1 = Gviz::coords(downstream_xy$GeneRegionTrack);  #// coordinates of all exons in geneB plotting 		
	#// extract exon_id harboring/close to breakpoint of geneA
	exon_match_A1 = grep(paste(upstream_flag[[as.character(break_A)]]$transcript$exon, collapse='|'), rownames(exon_xy_A1), value=T); 
	exon_match_pos_A1 = exon_xy_A1[exon_match_A1, ,drop=FALSE]; #// get coordinates of extracted exon_id for geneA
	#// extract exon_id harboring/close to breakpoint of geneB
	exon_match_B1 = grep(paste(downstream_flag[[as.character(break_B)]]$transcript$exon, collapse='|'), rownames(exon_xy_B1), value=T); 
	exon_match_pos_B1 = exon_xy_B1[exon_match_B1, ,drop=FALSE]; #// get coordinates of extracted exon_id for geneB

	#// process exon coordinates for geneA
	y_pos_gene_upstream = max(exon_match_pos_A1[, 4]); #// get coordinates of extracted exon_id at the bottom of plotting for geneA
	target_exon_A1 = rownames(subset(exon_match_pos_A1, exon_match_pos_A1[, 4] %in% y_pos_gene_upstream)); #// get name of the exon_id at the bottom of plotting for geneA
	target_exon_rename_A1 = sub("\\.[0-9]+", "", target_exon_A1, ignore.case=T); #// select core regular expression part (remove .d) of exon_id
	exon_tx_A1 = names(upstream_xy$GeneRegionTrack@imageMap@tags$origExonId); #// 'exon_tx_A1' is a vector of all exon_id
	names(exon_tx_A1) = upstream_xy$GeneRegionTrack@imageMap@tags$title; #// assign transcript_id to the name of each element in 'exon_tx_A1'
			
	tx_A1 = exon_tx_A1[exon_tx_A1 %in% target_exon_A1[1]]; #// get transcript_id at the bottom of plotting (if more than one transcripts on the bottom, please select first one)
	tx_exon_A1 = exon_tx_A1[names(exon_tx_A1) %in% names(tx_A1)]; #// select all exon_id of the transcript at the bottom of plotting
	#// select the exon_id (one exon may be selected twice due to split of utr-cds) related to the breakpoint in the vector 'tx_exon_A1'
	tx_exon_select_A1 = grep(target_exon_rename_A1[1], tx_exon_A1, value=T); 
	#// 0 <= fract_A1 <= 1
	fract_A1 = upstream_flag[[as.character(break_A)]]$transcript[upstream_flag[[as.character(break_A)]]$transcript$TXNAME==names(tx_A1),]$Prop; 
	exon_coord_A1 = subset(exon_match_pos_A1, rownames(exon_match_pos_A1) %in% tx_exon_select_A1); #// get coordinates of selected exon_id in 'tx_exon_select_A1'
	#// upstream_type = 1, 2 or 3
	upstream_type = upstream_flag[[as.character(break_A)]]$transcript[upstream_flag[[as.character(break_A)]]$transcript$TXNAME==names(tx_A1),]$Judge;
		
	#// calculate the x-coordinates of breakpoint for geneA
	if ( upstream_type == 1 || upstream_type == 2 || upstream_type == 3) { #// breakpoint at exon boundary (1), falls within exon (2) or intron (3)
		x2 = max(exon_coord_A1[,"x2"]); x1 = min(exon_coord_A1[,"x1"])
		x_pos_gene_upstream = x1 + (x2 - x1)*fract_A1;
	} else {
		stop("Judgement of breakpoint at upstream gene error in FUNCTION plot_curve_reduce!");
	}
	#// judge which transcript interval (utr5, CDS, utr3 and other) breakpoints fall in (for geneA)
	coding_A = upstream_flag[[as.character(break_A)]]$select_region;
	domain_breakpoint_A = upstream_flag[[as.character(break_A)]]$breakpoint_collect;
	cds_start_A = coding_A[coding_A$feature == "CDS", ][1,]$start; # or NA
	type_A = apply(coding_A, 1, function(x){ #// utr5, CDS, utr3 and other
		feature_A = ifelse(test = domain_breakpoint_A >= as.numeric(x[2]) && domain_breakpoint_A <= as.numeric(x[3]), yes = as.character(x[6]), no = "");
		return(feature_A);
	})
	type_A = type_A[type_A != ""];

	#// process exon coordinates for geneB
	y_pos_gene_downstream = min(exon_match_pos_B1[, 2]); #// get coordinates of extracted exon_id at the top of plotting for geneB
	target_exon_B1 = rownames(subset(exon_match_pos_B1, exon_match_pos_B1[,2] %in% y_pos_gene_downstream)); #// get name of the exon_id at the top of plotting for geneB
	target_exon_rename_B1 = sub("\\.[0-9]+", "", target_exon_B1, ignore.case=T); #// select core regular expression part (remove .d) of exon name			
	exon_tx_B1 = names(downstream_xy$GeneRegionTrack@imageMap@tags$origExonId); #// 'exon_tx_B1' is a vector of all exon_id
	names(exon_tx_B1) = downstream_xy$GeneRegionTrack@imageMap@tags$title; #// assign transcript_id to the name of each element in 'exon_tx_B1'
			
	tx_B1 = exon_tx_B1[exon_tx_B1 %in% target_exon_B1[1]]; #// get transcript_id at the top of plotting (if more than one transcripts on the top, please select first one)
	tx_exon_B1 = exon_tx_B1[names(exon_tx_B1) %in% names(tx_B1)]; #// select all exon_id of the transcript at the top of plotting
	#// select the exon_id (one exon may be selected twice due to split of utr-cds) related to the breakpoint in the vector 'tx_exon_B1' 
	tx_exon_select_B1 = grep(target_exon_rename_B1[1], tx_exon_B1, value=T); 
	#// 0 <= fract_B1 <= 1
	fract_B1 = downstream_flag[[as.character(break_B)]]$transcript[downstream_flag[[as.character(break_B)]]$transcript$TXNAME==names(tx_B1),]$Prop; 
	exon_coord_B1 = subset(exon_match_pos_B1, rownames(exon_match_pos_B1) %in% tx_exon_select_B1); #// get coordinates of selected exon_id in 'tx_exon_select_B1'
	#// upstream_type = 1, 2 or 3
	downstream_type = downstream_flag[[as.character(break_B)]]$transcript[downstream_flag[[as.character(break_B)]]$transcript$TXNAME==names(tx_B1),]$Judge;

	#// calculate the x-coordinates of breakpoint for geneB
	if ( downstream_type == 1 || downstream_type == 2 || downstream_type == 3) { #// breakpoint at exon boundary (1), falls within exon (2) or intron (3)
		x2 = max(exon_coord_B1[,"x2"]); x1 = min(exon_coord_B1[,"x1"]);
		x_pos_gene_downstream = x1 + (x2 - x1)*fract_B1;
	} else {
		stop("Judgement of breakpoint at downstream gene error in FUNCTION plot_curve_reduce!");		
	}
	#// judge which transcript interval (utr5, CDS, utr3 and other) breakpoints fall in (for geneB)	
	coding_B = downstream_flag[[as.character(break_B)]]$select_region;
	domain_breakpoint_B = downstream_flag[[as.character(break_B)]]$breakpoint_collect;
	cds_start_B = coding_B[coding_B$feature == "CDS", ][1,]$start; # or NA
	type_B = apply(coding_B, 1, function(x){ #// utr5, CDS, utr3 and other
		feature_B = ifelse(test = domain_breakpoint_B >= as.numeric(x[2]) && domain_breakpoint_B <= as.numeric(x[3]), yes = as.character(x[6]), no = "");
		return(feature_B);
	})
	type_B = type_B[type_B != ""];

	#// assign biological consequence of breakpoint
	#// tag = 0 (initial - "black"), 1 (Unknown - "black"), 2 (Inframe - "blue"), 3 (Outframe - "red"), 4 (Truncate-loss - "#008080")
	if ( length(type_A) > 0 && length(type_B) > 0 ) {
		if ( type_A == "utr5" ) {
			if ( type_B == "utr5" ) {
				tag = 2; # "Inframe";
			} else if ( type_B == "CDS" ) {
				tag = 4; # "Truncate-loss";
			} else if ( type_B == "utr3" ) {
				tag = 4; # "Truncate-loss";
			} else {
				tag = 4; # "Truncate-loss";
			}
		} else if ( type_A == "CDS" ) {
			if ( type_B == "utr5" ) {
				tag = 1; # "Unknown";
			} else if ( type_B == "CDS" ) {
				if ( upstream_type == 3 || downstream_type == 3 ) {
					tag = 1; # "Unknown";
				} else {
					if (! is.na(cds_start_A) && ! is.na(cds_start_B) ) {
						Yushu_A = (domain_breakpoint_A - 1 - cds_start_A + 1) %% 3;
						Yushu_B = (domain_breakpoint_B - cds_start_B + 1) %% 3;
						if ( Yushu_A == 0 && Yushu_B == 1 ) {
							tag = 2; # "Inframe";
						} else if ( Yushu_A == 1 && Yushu_B == 2 ) {
							tag = 2; # "Inframe";
						} else if ( Yushu_A == 2 && Yushu_B == 0 ) {
							tag = 2; # "Inframe";
						} else {
							tag = 3; # "Outframe";
						}
					} else {
						stop("CDS annotation is wrong!");
					}
				}
			} else if ( type_B == "utr3" ) {
				tag = 4; # "Truncate-loss";
			} else {
				tag = 4; # "Truncate-loss";
			}
		} else if ( type_A == "utr3" ) {
			if ( type_B == "utr5" ) {
				tag = 2; # "Inframe";
			} else if ( type_B == "CDS" ) {
				tag = 4; # "Truncate-loss";
			} else if ( type_B == "utr3" ) {
				tag = 2; # "Inframe";
			} else {
				tag = 2; # "Inframe";
			}
		} else {
			if ( type_B == "utr5" ) {
				tag = 2; # "Inframe";
			} else if ( type_B == "CDS" ) {
				tag = 4; # "Truncate-loss";
			} else if ( type_B == "utr3" ) {
				tag = 4; # "Truncate-loss";
			} else {
				tag = 1; # "Unknown";
			}
		}
	}

	#// Choose y position closest to the top of the plot. That is the lowest y value of the two, since the y axis is flipped
	y_pos_gene_upstream = max(exon_xy_A1[, 4])
	y_pos_gene_downstream = min(exon_xy_B1[, 2]);
	coordinate = c(x_pos_gene_upstream, x_pos_gene_downstream, y_pos_gene_upstream, y_pos_gene_downstream, tag); 
	names(coordinate) = c("x_pos_gene_upstream", "x_pos_gene_downstream", "y_pos_gene_upstream", "y_pos_gene_downstream", "tag");
	return(coordinate);
}	

#' Plot transcript-model(reducing introns), domain and motif
#'
#' @description Plot transcript-model(reducing introns), domain and motif of geneA in 'domain' tab-panel
#'
#' @param first A list: \strong{'key'} is a string \code{'\$pos'} (i.e. breakpoint pos in geneA); \code{'first\$pos\$transcript'} is a data.frame object
#'        that collects evaluation of breakpoint pos in geneA transcripts after reducing intron; \code{'first\$pos\$select_region'} is a data.frame object
#'        (5utr-cds-utr3 annotation after reducing intron) for plotting and it is a constant variable for different \code{'\$pos'} values.
#' @param first_name A string - symbol name of geneA (e.g. 'TMPRSS2').
#' @param first_domain A data.frame for collecting domain annotation (e.g. \code{'first_domain\$Start'}: domain start; \code{'first_domain\$End'}: domain end; 
#'        \code{'first_domain\$Domain_name'}: domain name abbreviation).
#' @param first_motif A data.frame for collecting motif annotation (e.g. \code{'first_motif\$Start'}: motif start; \code{'first_motif\$End'}: motif end; 
#'        \code{'first_motif\$Domain_name'}: motif name abbreviation).
#'
#' @return A data.frame (i.e. \code{'A1_xy'}) that contains exon coordinates in the plot layout.
#'
#' @export
plot_separate_domain_geneA <- function(first, first_name, first_domain, first_motif) {
	# For testing: first=geneA; first_name=symbol_A; first_domain=domain_geneA; first_motif=motif_geneA
	chrom_f = first[[1]]$transcript$Chrom[1] #// chromosome name of first: Ensembl
		
	#// Set plotting tracks for geneA after reducing intron
	grTrack_f = Gviz::GeneRegionTrack(first[[1]]$select_region, chromosome=chrom_f, showId=T, stacking="squish", geneSymbols=F, transcriptAnnotation="transcript", 
		just.group="above", fontsize.group = 11, fontsize=7, fontcolor="black", fontcolor.group="black", showTitle=F, col="black", col.line="black", fill="green", col.symbol="black", 
		lwd=0.1, lwd.border=1, lex=1, cex=0.6, cex.axis=0.6, background.title="transparent", background.panel="transparent", stackHeight=0.8, min.width=0.1, min.height=3, min.distance=0);	
	first_vis_s = as.numeric(first[[1]]$transcript$TXSTART) #// start position of geneA for visualization
	first_vis_e = as.numeric(first[[1]]$transcript$TXEND) #// end position of geneA for visualization	
	#// Set plotting tracks of domain and motif for geneA
	domain_f = NULL; motif_f = NULL;
	if (! is.na(first_domain[1,1]) ) {
		domain_f = Gviz::AnnotationTrack(start=first_domain$Start, end=first_domain$End, id=first_domain$Domain_name, name="Domain", col.title="black", 
			chromosome=chrom_f, fill="green", featureAnnotation = "id", shape="box", fontsize=7, fontcolor.feature = "black", cex.title=1.4)
	} else { #// if domain annotation not present
		domain_f = Gviz::AnnotationTrack(start=1, end=2, name="Domain", col.title="black", col="white", chromosome=chrom_f, fill="white", fontsize=7, cex.title=1.4)
	}
	if (! is.na(first_motif[1,1]) ) {
		motif_f = Gviz::AnnotationTrack(start=first_motif$Start, end=first_motif$End, id=first_motif$Domain_name, name="Motif", col.title="black", 
			chromosome=chrom_f, fill="green", featureAnnotation = "id", shape="box", fontsize=7, fontcolor.feature = "black", cex.title=1.4)
	} else { #// if motif annotation not present
		motif_f = Gviz::AnnotationTrack(start=1, end=2, name="Motif", col.title="black", col="white", chromosome=chrom_f, fill="white", fontsize=7, cex.title=1.4)
	}
		
	#// set the axis coordinate
	axis = Gviz::GenomeAxisTrack(name="Axis", labelPos="beside", fontcolor="black", littleTicks=F, cex=0.5, cex.id=0.5, cex.axis=0.5, 
		fontsize=16, col="black", distFromAxis=0.1, lwd=1, lwd.border=1, min.width=0.1, min.height=2, min.distance=0); 

	#////////////////////////////////////////////////////////////////////
	#// perform the plot ('motif', 'domain', 'axis' and 'gene' tracks) //
	#////////////////////////////////////////////////////////////////////
	offset = 0; #// define plot region [start-offset, end+offset]
	nf = grid::grid.layout(nrow=1, ncol = 1)
	grid::grid.newpage(); #// start a new page
	grid::pushViewport(grid::viewport(layout = nf)); #// divide the page to different section (nf specifies the section composite)

	#// plot 'motif', 'domain', 'axis' and 'gene' tracks in proportion as 0.6 : 0.6 : 0.25 : 1
	grid::pushViewport(grid::viewport(layout.pos.row = 1)) #// load on the section (1, 1)
	A1_xy = Gviz::plotTracks(c(motif_f, domain_f, axis, grTrack_f), chromosome=chrom_f, from=first_vis_s-offset, to=first_vis_e+offset, margin=0, 
			innerMargin=1, sizes=c(0.6, 0.6, 0.3, 1), add=T);
	grid::popViewport(1); #// close the section (1, 1)
		
	#// 'A1_xy$GeneRegionTrack': coordinates(x1, y1, x2, y2) of all exons of transcripts for geneA after reducing intron
	return(A1_xy)
}

#' Plot transcript-model(reducing introns), domain and motif
#'
#' @description Plot transcript-model(reducing introns), domain and motif of geneB in 'domain' tab-panel
#'
#' @param second A list: \strong{'key'} is a string \code{'\$pos'} (i.e. breakpoint pos in geneB); \code{'second\$pos\$transcript'} is a data.frame object
#'        that collects evaluation of breakpoint pos in geneB transcript after reducing intron; \code{'second\$pos\$select_region'} is a data.frame object
#'        (5utr-cds-utr3 annotation after reducing intron) for plotting and it is a constant variable for different \code{'\$pos'} values.
#' @param second_name A string - symbol name of geneB (e.g. 'ERG').
#' @param second_domain A data.frame for collecting domain annotation (e.g. \code{'second_domain\$Start'}: domain start; \code{'second_domain\$End'}: domain end; 
#'        \code{'second_domain\$Domain_name'}: domain name abbreviation).
#' @param second_motif A data.frame for collecting motif annotation (e.g. \code{'second_motif\$Start'}: motif start; \code{'secondt_motif\$End'}: motif end;
#'        \code{'second_motif\$Domain_name'}: motif name abbreviation).
#'
#' @return A data.frame (i.e. \code{'B1_xy'}) that contains exon coordinates in the plot layout.
#'
#' @export
plot_separate_domain_geneB <- function(second, second_name, second_domain, second_motif) {
	# For testing: second=geneB; second_name=symbol_B; second_domain=domain_geneB; second_motif=motif_geneB; cytoband=chrTrack
	chrom_s = second[[1]]$transcript$Chrom[1] #// chromosome name of second: Ensembl

	#// Set plotting tracks for geneB after reducing intron
	grTrack_s = Gviz::GeneRegionTrack(second[[1]]$select_region, chromosome=chrom_s, showId=T, stacking="squish", geneSymbols=F, transcriptAnnotation="transcript", just.group="below", 
		fontsize.group = 11, fontsize=7, fontcolor="black", fontcolor.group="black", showTitle=F, col="black", col.line="black", fill="orange", col.symbol="black", lwd=0.1, lwd.border=1, 
		lex=1, cex=0.6, cex.axis=0.6, background.title="transparent", background.panel="transparent", stackHeight=0.8, min.width=0.1, min.height=3, min.distance=0);
	second_vis_s = as.numeric(second[[1]]$transcript$TXSTART)	#// start position of geneB for visualization
	second_vis_e = as.numeric(second[[1]]$transcript$TXEND) #// start position of geneB for visualization
	#// Set plotting tracks of domain and motif for geneB
	domain_s = NULL; motif_s = NULL;
	if (! is.na(second_domain[1,1]) ) {
		domain_s = Gviz::AnnotationTrack(start=second_domain$Start, end=second_domain$End, id=second_domain$Domain_name, name="Domain", col.title="black", 
			chromosome=chrom_s, fill="orange", featureAnnotation = "id", shape="box", fontsize=7, fontcolor.feature = "black", cex.title=1.4)
	} else { #// if domain annotation not present
		domain_s = Gviz::AnnotationTrack(start=1, end=2, name="Domain", col.title="black", col="white", chromosome=chrom_s, fill="white", fontsize=7, cex.title=1.4)
	}
	if (! is.na(second_motif[1,1]) ) {
		motif_s = Gviz::AnnotationTrack(start=second_motif$Start, end=second_motif$End, id=second_motif$Domain_name, name="Motif", col.title="black", 
			chromosome=chrom_s, fill="orange", featureAnnotation = "id", shape="box", fontsize=7, fontcolor.feature = "black", cex.title=1.4)
	} else { #// if motif annotation not present
		motif_s = Gviz::AnnotationTrack(start=1, end=2, name="Motif", col.title="black", col="white", chromosome=chrom_s, fill="white", fontsize=7, cex.title=1.4)
	}

	#// set the axis coordinate
	axis = Gviz::GenomeAxisTrack(name="Axis", labelPos="beside", fontcolor="black", littleTicks=F, cex=0.5, cex.id=0.5, cex.axis=0.5, 
		fontsize=16, col="black", distFromAxis=0.1, lwd=1, lwd.border=1, min.width=0.1, min.height=2, min.distance=0); 

	#////////////////////////////////////////////////////////////////////
	#// perform the plot ('gene', 'axis', 'domain' and 'motif' tracks) //
	#////////////////////////////////////////////////////////////////////
	offset = 0; #// define axis plot region [start-offset, end+offset]
	nf = grid::grid.layout(nrow=1, ncol = 1)
	grid::grid.newpage(); #// start a new page
	grid::pushViewport(grid::viewport(layout = nf)); #// divide the page to different section (nf specifies the section composite)

	#// plot 'gene', 'axis', 'domain' and 'motif' in proportion as 1 : 0.25 : 0.6 : 0.6
	grid::pushViewport(grid::viewport(layout.pos.row = 1)) #// load on the section (1, 2)
	B1_xy = Gviz::plotTracks(c(grTrack_s, axis, domain_s, motif_s), chromosome=chrom_s, from=second_vis_s-offset, to=second_vis_e+offset, 
			margin=0, innerMargin=1, sizes=c(1, 0.3, 0.6, 0.6), add=T); # c(1, 0.25, 0.6, 0.6)
	grid::popViewport(1); #// close the section (1, 1)
	
	#// 'B1_xy$GeneRegionTrack': coordinates(x1, y1, x2, y2) of all exons of transcripts for geneB after reducing intron
	return(B1_xy)
}

#' Plot multiple straight lines that link breakpoints of geneA and geneB transcripts
#'
#' @description Plot multiple straight lines that link breakpoints of geneA and geneB transcripts (reducing introns) in 'domain' tab-panel
#'
#' @param A1_xy A data.frame: \code{'A1_xy\$GeneRegionTrack'} - exon coordinates(x1, y1, x2, y2) of intron-reduced transcripts,
#'        \code{'A1_xy'} is the return object of FUNCTION "plot_separate_domain_geneA".
#' @param B1_xy A data.frame: \code{'B1_xy\$GeneRegionTrack'} - exon coordinates(x1, y1, x2, y2) of intron-reduced transcripts,
#'        \code{'B1_xy'} is the return object of FUNCTION "plot_separate_domain_geneB".
#' @param first A list: \strong{'key'} is a string \code{'\$pos'} (i.e. breakpoint pos in geneA); \code{'first\$pos\$transcript'} is a data.frame object
#'        that collects evaluation of breakpoint pos in geneA transcript after reducing intron; \code{'first\$pos\$select_region'} is a data.frame object
#'        (5utr-cds-utr3 annotation after reducing intron) for plotting and it is a constant variable for different \code{'\$pos'} values.
#' @param second A list: \strong{'key'} is a string \code{'\$pos'} (i.e. breakpoint pos in geneB); \code{'second\$pos\$transcript'} is a data.frame object
#'        that collects evaluation of breakpoint pos in geneB transcript after reducing intron; \code{'second\$pos\$select_region'} is a data.frame object
#'        (5utr-cds-utr3 annotation after reducing intron) for plotting and it is a constant variable for different \code{'\$pos'} values.
#' @param breakpoint A data.frame with two columns (e.g. 'breakpoint[i,1]' - breakpoint pos in geneA transcript; 'breakpoint[i,2]' - breakpoint pos in geneB transcript).
#'
#' @return NULL
#'         NOTE: In terms of biological consequences of fusion line, there are four corresponding colors (i.e. black - Unknown, blue - Inframe, red - Outframe and '#008080' - Truncate-loss).
#'
#' @export
plot_separate_domain_arrow <- function(A1_xy, B1_xy, first, second, breakpoint) {
	# For testing: print("GeneA_exon_coord: "); print(Gviz::coords(A1_xy$GeneRegionTrack)); print("GeneB_exon_coord: "); print(Gviz::coords(B1_xy$GeneRegionTrack))
	if ( nrow(breakpoint) > 0 ) { #// make sure 'breakpoint' data.frame is not empty
		for (i in 1:length(breakpoint[,1])) {
			curve_coordinate = .plot_curve_reduce(A1_xy, first, B1_xy, second, breakpoint[i,]);
			#// curve_coordinate - a vector class that is composed of four elements (i.e. "x_pos_gene_upstream": x-coordinate of geneA breakpoint, 
			#// "x_pos_gene_downstream": x-coordinate of geneA breakpoint, 
			#// "y_pos_gene_upstream": y-coordinate of geneA breakpoint, "y_pos_gene_downstream": y-coordinate of geneB breakpoint)
				
			if (! is.null(curve_coordinate) ) { 
				#// NOTE: y-pos in 'curve_coordinate' data.frame are coordinates of output$domainA and output$domainB, which should be converted to that of output$fuseline
				curve_coordinate["y_pos_gene_upstream"] = 0
				#// NOTE: this value should be identical to the setting "div(class="fuseline", plotOutput("linking", height=60))" in ui.R
				curve_coordinate["y_pos_gene_downstream"] = 60;
				#// Draw the straight lines between the transcripts
				grid::pushViewport(grid::viewport(xscale = c(0, grDevices::dev.size(units = "px")[1]), yscale = c(grDevices::dev.size(units = "px")[2], 0)));
				if ( curve_coordinate[5] == 1 ) { #// 'Unknown' type
					grid::grid.lines(c(curve_coordinate["x_pos_gene_upstream"], curve_coordinate["x_pos_gene_downstream"]), #// x positions
						c(curve_coordinate["y_pos_gene_upstream"], curve_coordinate["y_pos_gene_downstream"]), #// y positions
						default.units = "native", arrow=grid::arrow(angle=50, length=grid::unit(0.05, "inches"), ends="last", type="closed"), 
						gp = grid::gpar(col="black", fill="black", lwd=0.8));
				} else if ( curve_coordinate[5] == 2 ) { #// 'Inframe' type
					grid::grid.lines(c(curve_coordinate["x_pos_gene_upstream"], curve_coordinate["x_pos_gene_downstream"]), #// x positions
						c(curve_coordinate["y_pos_gene_upstream"], curve_coordinate["y_pos_gene_downstream"]), #// y positions
						default.units = "native", arrow=grid::arrow(angle=50, length=grid::unit(0.05, "inches"), ends="last", type="closed"),
						gp = grid::gpar(col="blue", fill="blue", lwd=0.8));
				} else if ( curve_coordinate[5] == 3 ) { #// 'Outframe' type
					grid::grid.lines(c(curve_coordinate["x_pos_gene_upstream"], curve_coordinate["x_pos_gene_downstream"]), #// x positions
						c(curve_coordinate["y_pos_gene_upstream"], curve_coordinate["y_pos_gene_downstream"]), #// y positions
						default.units = "native", arrow=grid::arrow(angle=50, length=grid::unit(0.05, "inches"), ends="last", type="closed"),
						gp = grid::gpar(col="red", fill="red", lwd=0.8));
				} else if ( curve_coordinate[5] == 4 ) { #// 'Truncate-loss' type
					grid::grid.lines(c(curve_coordinate["x_pos_gene_upstream"], curve_coordinate["x_pos_gene_downstream"]), #// x positions
						c(curve_coordinate["y_pos_gene_upstream"], curve_coordinate["y_pos_gene_downstream"]), #// y positions
						default.units = "native", arrow=grid::arrow(angle=50, length=grid::unit(0.05, "inches"), ends="last", type="closed"),
						gp = grid::gpar(col="#008080", fill="#008080", lwd=0.8));
				} else {
					grid::grid.lines(c(curve_coordinate["x_pos_gene_upstream"], curve_coordinate["x_pos_gene_downstream"]), #// x positions
						c(curve_coordinate["y_pos_gene_upstream"], curve_coordinate["y_pos_gene_downstream"]), #// y positions
						default.units = "native", arrow=grid::arrow(angle=50, length=grid::unit(0.05, "inches"), ends="last", type="closed"),
						gp = grid::gpar(col="black", fill="black", lwd=0.8));
				}
			}
		}
	}
}

#' Plot transcript-model(reducing introns), domain, motif and multiple straight lines that link breakpoints of geneA and geneB transcripts (plot downloading version)
#'
#' @description Plot transcript-model(reducing introns), domain, motif and multiple straight lines that link breakpoints of geneA and geneB transcripts together.
#'
#' @export
plot_separate_domain_download <- function(first, first_name, first_domain, first_motif, second, second_name, second_domain, second_motif, breakpoint) {

	chrom_f = first[[1]]$transcript$Chrom[1] #// chromosome name of first: Ensembl
	#// Set plotting tracks for geneA after reducing intron
	grTrack_f = Gviz::GeneRegionTrack(first[[1]]$select_region, chromosome=chrom_f, showId=T, stacking="squish", geneSymbols=F, transcriptAnnotation="transcript", just.group="above", 
		fontsize.group = 11, fontsize=7, fontcolor="black", fontcolor.group="black", showTitle=F, col="black", col.line="black", fill="green", col.symbol="black", lwd=0.1, 
		lwd.border=1, lex=1, cex=0.6, cex.axis=0.6, background.title="transparent", background.panel="transparent", stackHeight=0.8, min.width=0.1, min.height=3, min.distance=0);     
	first_vis_s = as.numeric(first[[1]]$transcript$TXSTART) #// start position of geneA for visualization
	first_vis_e = as.numeric(first[[1]]$transcript$TXEND) #// end position of geneA for visualization       
	#// Set plotting tracks of domain and motif for geneA
	domain_f = NULL; motif_f = NULL;
	if (! is.na(first_domain[1,1]) ) {
		domain_f = Gviz::AnnotationTrack(start=first_domain$Start, end=first_domain$End, id=first_domain$Domain_name, name="Domain", col.title="black", 
			chromosome=chrom_f, fill="green", featureAnnotation = "id", shape="box", fontsize=7, fontcolor.feature = "black", cex.title=1.4)
	} else { #// if domain annotation not present
		domain_f = Gviz::AnnotationTrack(start=1, end=2, name="Domain", col.title="black", col="white", chromosome=chrom_f, fill="white", fontsize=7, cex.title=1.4)
	}
	if (! is.na(first_motif[1,1]) ) {
		motif_f = Gviz::AnnotationTrack(start=first_motif$Start, end=first_motif$End, id=first_motif$Domain_name, name="Motif", col.title="black", 
			chromosome=chrom_f, fill="green", featureAnnotation = "id", shape="box", fontsize=7, fontcolor.feature = "black", cex.title=1.4)
	} else { #// if motif annotation not present
		motif_f = Gviz::AnnotationTrack(start=1, end=2, name="Motif", col.title="black", col="white", chromosome=chrom_f, fill="white", fontsize=7, cex.title=1.4)
	}

	chrom_s = second[[1]]$transcript$Chrom[1] #// chromosome name of second: Ensembl
	#// Set plotting tracks for geneB after reducing intron
	grTrack_s = Gviz::GeneRegionTrack(second[[1]]$select_region, chromosome=chrom_s, showId=T, stacking="squish", geneSymbols=F, transcriptAnnotation="transcript", just.group="below", 
		fontsize.group = 11, fontsize=7, fontcolor="black", fontcolor.group="black", showTitle=F, col="black", col.line="black", fill="orange", col.symbol="black", lwd=0.1, 
		lwd.border=1, lex=1, cex=0.6, cex.axis=0.6, background.title="transparent", background.panel="transparent", stackHeight=0.8, min.width=0.1, min.height=3, min.distance=0);
	second_vis_s = as.numeric(second[[1]]$transcript$TXSTART)       #// start position of geneB for visualization
	second_vis_e = as.numeric(second[[1]]$transcript$TXEND) #// start position of geneB for visualization
	#// Set plotting tracks of domain and motif for geneB
	domain_s = NULL; motif_s = NULL;
	if (! is.na(second_domain[1,1]) ) {
		domain_s = Gviz::AnnotationTrack(start=second_domain$Start, end=second_domain$End, id=second_domain$Domain_name, name="Domain", col.title="black", 
			chromosome=chrom_s, fill="orange", featureAnnotation = "id", shape="box", fontsize=7, fontcolor.feature = "black", cex.title=1.4)
	} else { #// if domain annotation not present
		domain_s = Gviz::AnnotationTrack(start=1, end=2, name="Domain", col.title="black", col="white", chromosome=chrom_s, fill="white", fontsize=7, cex.title=1.4)
	}
	if (! is.na(second_motif[1,1]) ) {
		motif_s = Gviz::AnnotationTrack(start=second_motif$Start, end=second_motif$End, id=second_motif$Domain_name, name="Motif", col.title="black", 
			chromosome=chrom_s, fill="orange", featureAnnotation = "id", shape="box", fontsize=7, fontcolor.feature = "black", cex.title=1.4)
	} else { #// if motif annotation not present
		motif_s = Gviz::AnnotationTrack(start=1, end=2, name="Motif", col.title="black", col="white", chromosome=chrom_s, fill="white", fontsize=7, cex.title=1.4)
	}

	#// set the axis coordinate
	axis = Gviz::GenomeAxisTrack(name="Axis", labelPos="beside", fontcolor="black", littleTicks=F, cex=0.5, cex.id=0.5, cex.axis=0.5,
		fontsize=14, col="black", distFromAxis=0.1, lwd=1, lwd.border=1, min.width=0.1, min.height=2, min.distance=0);

	#// plot domain and fusion events
	offset = 0; #// define plot region [start-offset, end+offset]
	nf = grid::grid.layout(nrow=3, ncol=2, widths=grid::unit(c(1, 1), "null"), heights = grid::unit(c(36, 1, 36), "null"))
	grid::grid.newpage(); #// start a new page
	grid::pushViewport(grid::viewport(layout = nf)); #// divide the page to different sections (nf specifies the section composite)

	grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = NULL)); #// load on the section (3, null)
	A1_xy = Gviz::plotTracks(c(motif_f, domain_f, axis, grTrack_f), chromosome=chrom_f, from=first_vis_s-offset, to=first_vis_e+offset, margin=0, 
			innerMargin=1, sizes=c(0.6, 0.6, 0.25, 1), add=T);
	grid::popViewport(1); #// close the section (1)

	grid::pushViewport(grid::viewport(layout.pos.row = 3, layout.pos.col = NULL)) #// load on the section (1, 2)
	B1_xy = Gviz::plotTracks(c(grTrack_s, axis, domain_s, motif_s), chromosome=chrom_s, from=second_vis_s-offset, to=second_vis_e+offset, margin=0,
			innerMargin=1, sizes=c(1, 0.25, 0.6, 0.6), add=T);
	grid::popViewport(1); #// close the section (3)

	if ( nrow(breakpoint) > 0 ) { #// make sure 'breakpoint' data.frame is not empty
		for (i in 1:length(breakpoint[,1])) {

			curve_coordinate = .plot_curve_reduce(A1_xy, first, B1_xy, second, breakpoint[i,]);
			if (! is.null(curve_coordinate) ) { 
				#// NOTE: y-pos in 'curve_coordinate' data.frame are coordinates of output$domainA and output$domainB, which should be converted to that of output$fuseline
				curve_coordinate["y_pos_gene_upstream"] = curve_coordinate["y_pos_gene_upstream"] + 1;
				#// Draw the straight lines between the transcripts
				grid::pushViewport(grid::viewport(xscale = c(0, grDevices::dev.size(units = "px")[1]), yscale = c(grDevices::dev.size(units = "px")[2], 0)));
				if ( curve_coordinate[5] == 1 ) { #// 'Unknown' type
					grid::grid.lines(c(curve_coordinate["x_pos_gene_upstream"], curve_coordinate["x_pos_gene_downstream"]), #// x positions
								c(curve_coordinate["y_pos_gene_upstream"], curve_coordinate["y_pos_gene_downstream"]), #// y positions
								default.units = "native", arrow=grid::arrow(angle=50, length=grid::unit(0.05, "inches"), ends="last", type="closed"),
								gp = grid::gpar(col="black", fill="black", lwd=0.6));
				} else if ( curve_coordinate[5] == 2 ) { #// 'Inframe' type
					grid::grid.lines(c(curve_coordinate["x_pos_gene_upstream"], curve_coordinate["x_pos_gene_downstream"]), #// x positions
								c(curve_coordinate["y_pos_gene_upstream"], curve_coordinate["y_pos_gene_downstream"]), #// y positions
								default.units = "native", arrow=grid::arrow(angle=50, length=grid::unit(0.05, "inches"), ends="last", type="closed"),
								gp = grid::gpar(col="blue", fill="blue", lwd=0.6));
				} else if ( curve_coordinate[5] == 3 ) { #// 'Outframe' type
					grid::grid.lines(c(curve_coordinate["x_pos_gene_upstream"], curve_coordinate["x_pos_gene_downstream"]), #// x positions
								c(curve_coordinate["y_pos_gene_upstream"], curve_coordinate["y_pos_gene_downstream"]), #// y positions
								default.units = "native", arrow=grid::arrow(angle=50, length=grid::unit(0.05, "inches"), ends="last", type="closed"),
								gp = grid::gpar(col="red", fill="red", lwd=0.6));
				} else if ( curve_coordinate[5] == 4 ) { #// 'Truncate-loss' type
					grid::grid.lines(c(curve_coordinate["x_pos_gene_upstream"], curve_coordinate["x_pos_gene_downstream"]), #// x positions
								c(curve_coordinate["y_pos_gene_upstream"], curve_coordinate["y_pos_gene_downstream"]), #// y positions
								default.units = "native", arrow=grid::arrow(angle=50, length=grid::unit(0.05, "inches"), ends="last", type="closed"), 
								gp = grid::gpar(col="#008080", fill="#008080", lwd=0.6));
				} else {
					grid::grid.lines(c(curve_coordinate["x_pos_gene_upstream"], curve_coordinate["x_pos_gene_downstream"]), #// x positions
								c(curve_coordinate["y_pos_gene_upstream"], curve_coordinate["y_pos_gene_downstream"]), #// y positions
								default.units = "native", arrow=grid::arrow(angle=50, length=grid::unit(0.05, "inches"), ends="last", type="closed"),
								gp = grid::gpar(col="black", fill="black", lwd=0.6));
				}
			}
		}
	}
}

