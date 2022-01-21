#' Plot transcript-model and chromosome ideogram of partner genes
#'
#' @description Plot transcript-model and chromosome ideogram of partner genes in 'overview' tab-panel
#'
#' @param first A list - \strong{'key'} is \code{'\$pos'} (i.e. breakpoint pos of geneA); \code{'first\$pos\$transcript'} is a data.frame object
#'        (only one row) that collects evaluation of breakpoint pos for geneA transcripts; \code{'first\$pos\$select_region'} is a data.frame object
#'        (5utr-cds-utr3 annotation) for plotting and it is a constant variable for different \code{'\$pos'} values.
#' @param first_name A string - symbol name of geneA (e.g. 'TMPRSS2').
#' @param second A list - \strong{'key'} is \code{'\$pos'} (i.e. breakpoint pos of geneB); \code{'second\$pos\$transcript'} is a data.frame object
#'        (only one row) that collects evaluation of breakpoint pos for geneB transcripts; \code{'second\$pos\$select_region'} is a data.frame object
#'        (5utr-cds-utr3 annotation) for plotting and it is a constant variable for different \code{'\$pos'} values.
#' @param second_name A string - symbol name of geneB (e.g. 'ERG').
#' @param cytoband A GRange object for plotting cytoband of chromosome ideogram.
#' @param offset A numeric value - extend the partner gene interval (default: 2000, e.g. [start-offset, end+offset]).
#'
#' @return A list with two data.frame elements (i.e. \code{'A1_xy'} and \code{'B1_xy'}) that contain exon coordinates in the plot layout.
#'
#' @export
plot_separate_overview <- function(first, first_name, second, second_name, cytoband, offset=2000) {
	chrom_f = first[[1]]$transcript$Chrom[1] #// chromosome name of first: Ensembl
	chrom_s = second[[1]]$transcript$Chrom[1] #// chromosome name of second: Ensembl

	#// Set plotting tracks for geneA
	grTrack_f = Gviz::GeneRegionTrack(first[[1]]$select_region, chromosome=chrom_f, showId=T, stacking="squish", geneSymbols=F, transcriptAnnotation="transcript", just.group="below", 
		fontsize=5, fontcolor="black", fontcolor.group="black", showTitle=F, col="black", col.line="black", fill="green", col.symbol="black", lwd=0.1, lwd.border=1, 
		lex=1, cex=0.6, cex.axis=0.6, background.title="transparent", background.panel="transparent", stackHeight=0.8, min.width=0.1, min.height=3, min.distance=0);	
	first_vis_s = min(first[[1]]$transcript$GStart) #// start position of geneA for visualization 
	first_vis_e = max(first[[1]]$transcript$GEnd) #// end position of geneA for visualization 

	#// Set plotting tracks for geneB
	grTrack_s = Gviz::GeneRegionTrack(second[[1]]$select_region, chromosome=chrom_s, showId=T, stacking="squish", geneSymbols=F, transcriptAnnotation="transcript", 
		just.group="below", fontsize=5, fontcolor="black", fontcolor.group="black", showTitle=F, col="black", col.line="black", fill="orange", col.symbol="black", lwd=0.1, 
		lwd.border=1, lex=1, cex=0.6, cex.axis=0.6, background.title="transparent", background.panel="transparent", stackHeight=0.8, min.width=0.1, min.height=3, min.distance=0);
	second_vis_s = min(second[[1]]$transcript$GStart) #// start position of geneB for visualization 
	second_vis_e = max(second[[1]]$transcript$GEnd) #// end position of geneB for visualization

	#// adjust font size of transcript_id (it will probably be adjusted by customers in UI in further version)
	if ( length(first[[1]]$transcript$TXNAME) > 0 &&  length(first[[1]]$transcript$TXNAME) <= 5 ) { 
		grTrack_f@dp@pars$fontsize.group = 12; 
	} else if ( length(first[[1]]$transcript$TXNAME) > 5 &&	 length(first[[1]]$transcript$TXNAME) <= 15 ) { 
		grTrack_f@dp@pars$fontsize.group = 10; 
	} else if ( length(first[[1]]$transcript$TXNAME) > 15 &&  length(first[[1]]$transcript$TXNAME) <= 20 ) {
		grTrack_f@dp@pars$fontsize.group = 9; 
	} else if ( length(first[[1]]$transcript$TXNAME) > 20 &&  length(first[[1]]$transcript$TXNAME) <= 25 ) {
		grTrack_f@dp@pars$fontsize.group = 8;
	} else if ( length(first[[1]]$transcript$TXNAME) > 25 &&  length(first[[1]]$transcript$TXNAME) <= 30 ) {
		grTrack_f@dp@pars$fontsize.group = 6;
	} else {
		grTrack_f@dp@pars$fontsize.group = 4;
	}
	if ( length(second[[1]]$transcript$TXNAME) > 0 &&  length(second[[1]]$transcript$TXNAME) <= 5 ) { 
		grTrack_s@dp@pars$fontsize.group = 12; 
	} else if ( length(second[[1]]$transcript$TXNAME) > 5 &&  length(second[[1]]$transcript$TXNAME) <= 15 ) { 
		grTrack_s@dp@pars$fontsize.group = 10; 
	} else if ( length(second[[1]]$transcript$TXNAME) > 15 &&  length(second[[1]]$transcript$TXNAME) <= 20 ) {
		grTrack_s@dp@pars$fontsize.group = 9; 
	} else if ( length(second[[1]]$transcript$TXNAME) > 20 &&  length(second[[1]]$transcript$TXNAME) <= 25 ) {
		grTrack_s@dp@pars$fontsize.group = 8;
	} else if ( length(second[[1]]$transcript$TXNAME) > 25 &&  length(second[[1]]$transcript$TXNAME) <= 30 ) {
		grTrack_s@dp@pars$fontsize.group = 6;
	} else {
		grTrack_s@dp@pars$fontsize.group = 4;
	}

	#// set the chromosome axis coordinate
	axis = Gviz::GenomeAxisTrack(name="Axis", labelPos="alternating", fontcolor="black", littleTicks=F, cex=0.5, cex.id=0.5, cex.axis=0.5, 
		fontsize=16, col="black", distFromAxis=0.4, lwd=1.4, lwd.border=1, min.width=0.1, min.height=3, min.distance=0); 

	#//////////////////////////////////////////////////////////////////////
	#// perform the plot ('gene+axis', 'title' and 'chromosome' tracks)	 //
	#//////////////////////////////////////////////////////////////////////
	#// define plot region [start-offset, end+offset]
	#// plot 'gene+axis', 'title' and 'chromosome' tracks in proportion as 10 : 0.2 : 0.8 (can be adjusted)
	nf = grid::grid.layout(nrow=3, ncol=2, widths=grid::unit(c(1, 1), "null"), heights = grid::unit(c(10, 0.2, 0.8), "null"))
	grid::grid.newpage(); #// start a new page
	grid::pushViewport(grid::viewport(layout = nf)); #// divide the page to different sections (nf specifies the section composite)

	#// plot 'geneA' and 'axis' tracks in proportion as 1 : 0.2 (can be adjusted)
	grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1)) #// load on the section (1, 1)
	if ( first[[1]]$transcript$Strand[1] == "+" ) {
		A1_xy = Gviz::plotTracks(c(grTrack_f, axis), chromosome=chrom_f, from=first_vis_s-offset, to=first_vis_e+offset, 
			margin=0, innerMargin=1, sizes=c(1, 0.2), add=T);
	} else {
		A1_xy = Gviz::plotTracks(c(grTrack_f, axis), chromosome=chrom_f, from=first_vis_s-offset, to=first_vis_e+offset, 
			margin=0, innerMargin=1, sizes=c(1, 0.2), reverseStrand=TRUE, add=T);
	}
	grid::popViewport(1); #// close the section (1, 1)
	#// plot 'geneB' and 'axis' tracks in proportion as 1 : 0.2 (can be adjusted)
	grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 2)) #// load on the section (1, 2)
	if ( second[[1]]$transcript$Strand[1] == "+" ) {
		B1_xy = Gviz::plotTracks(c(grTrack_s, axis), chromosome=chrom_s, from=second_vis_s-offset, to=second_vis_e+offset, 
			margin=0, innerMargin=1, sizes=c(1, 0.2), add=T);
	} else {
		B1_xy = Gviz::plotTracks(c(grTrack_s, axis), chromosome=chrom_s, from=second_vis_s-offset, to=second_vis_e+offset, 
			margin=0, innerMargin=1, sizes=c(1, 0.2), reverseStrand=TRUE, add=T);
	}
	grid::popViewport(1); #// close the section (1, 2)
		
	#// add sub-title for geneA
	grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = 1)); #// load on the section (2, 1)
	grid::grid.text(label=paste("Upstream gene: ", first_name, " (", first[[1]]$transcript$Strand[1], ")", sep=""), 
			x=grid::unit(0.5, "npc"), y=grid::unit(0.85, "npc"), gp=grid::gpar(fontsize=10, fontface="bold", col="black"));
	grid::popViewport(1); #// close the section (2, 1)
	#// add sub-title for geneB
	grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = 2)); #// load on the section (2, 2)
	grid::grid.text(label=paste("Downstream gene: ", second_name, " (", second[[1]]$transcript$Strand[1], ")", sep=""), 
			x=grid::unit(0.5, "npc"), y=grid::unit(0.85, "npc"), gp=grid::gpar(fontsize=10, fontface="bold", col="black"));		
	grid::popViewport(1); #// close the section (2, 2)

	if ( chrom_f == chrom_s ) { #// if both geneA and geneB at the same chromosome
		#// plot cytoband of the chromosome
		names(cytoband) = chrom_f;
		grid::pushViewport(grid::viewport(layout.pos.row = 3, layout.pos.col = NULL)); #// load on the section (3, null)
		if ( first[[1]]$transcript$GStart[1] < second[[1]]$transcript$GStart[1] ) {
			Gviz::plotTracks(cytoband, chromosome=chrom_f, from=first_vis_s, to=second_vis_e, add=T, main="", cex.main=0.6, fontsize=15, fontface.main=2, margin=1, innerMargin=1)
		} else {
			Gviz::plotTracks(cytoband, chromosome=chrom_s, from=second_vis_s, to=first_vis_e, add=T, main="", cex.main=0.6, fontsize=15, fontface.main=2, margin=1, innerMargin=1)
		}
		grid::popViewport(1); 
		#// plot fusion gene title
		grid::pushViewport(grid::viewport(layout.pos.row = 3, layout.pos.col = NULL));
		grid::grid.text(label=paste("Fusion: ", first_name, "-", second_name, sep=""), x=grid::unit(0.5, "npc"), y=grid::unit(0.85, "npc"), 
				gp=grid::gpar(fontsize=10, fontface="bold", col="black"));
		grid::popViewport(1); #// close the section (3, null)
	} else { #// if geneA and geneB at the different chromosomes
		cytoband_A = cytoband;	names(cytoband_A) = chrom_f;
		cytoband_B = cytoband;	names(cytoband_B) = chrom_s;
		#// plot cytoband of geneA chromosome 
		grid::pushViewport(grid::viewport(layout.pos.row = 3, layout.pos.col = 1)); #// load on the section (3, 1)
		Gviz::plotTracks(cytoband_A, chromosome=chrom_f, from=first_vis_s, to=first_vis_e, add=T, cex.main=0.6, fontsize=15, fontface.main=2, margin=1, innerMargin=1);
		grid::popViewport(1); #// close the section (3, 1)
		#// plot cytoband of geneB chromosome
		grid::pushViewport(grid::viewport(layout.pos.row = 3, layout.pos.col = 2)); #// load on the section (3, 2)
		Gviz::plotTracks(cytoband_B, chromosome=chrom_s, from=second_vis_s, to=second_vis_e, add=T, cex.main=0.6, fontsize=15, fontface.main=2, margin=1, innerMargin=1);
		grid::popViewport(1); #// close the section (3, 2)
		#// plot fusion gene title
		grid::pushViewport(grid::viewport(layout.pos.row = 3, layout.pos.col = NULL)); #// load on the section (3, null)
		grid::grid.text(label=paste("Fusion: ", first_name, "-", second_name, sep=""), x=grid::unit(0.5, "npc"), y=grid::unit(0.85, "npc"), 
				gp=grid::gpar(fontsize=10, fontface="bold", col="black"))
		grid::popViewport(1); #// close the section (3, null)
	}
	#// 'A1_xy$GeneRegionTrack' is a data.frame that collects coordinates(x1, y1, x2, y2) of all exons of transcripts for plotting geneA 
	#// 'B1_xy$GeneRegionTrack' is a data.frame that collects coordinates(x1, y1, x2, y2) of all exons of transcripts for plotting geneB
	plotbody = list(A1_xy=A1_xy, B1_xy=B1_xy)
	return(plotbody)
}

#' Plot multiple curved lines that link breakpoints of geneA and geneB transcripts
#'
#' @description Plot multiple curved lines that link breakpoints of geneA and geneB transcripts in 'overview' tab-panel
#'
#' @param collect A list that is a return object in FUNCTION "plot_separate_overview" with two elements (i.e. \code{'collect\$A1_xy\$GeneRegionTrack'} and 
#'        \code{'collect\$B1_xy\$GeneRegionTrack'}) - exon coordinates(x1, y1, x2, y2) of geneA and geneB transcripts.
#' @param first A list - \strong{'key'} is \code{'\$pos'} (i.e. breakpoint pos of geneA); \code{'first\$pos\$transcript'} is a data.frame object
#'        (only one row) that collects evaluation of breakpoint pos for geneA transcripts; \code{'first\$pos\$select_region'} is a data.frame object
#'        (5utr-cds-utr3 annotation) for plotting and it is a constant variable for different \code{'\$pos'} values.
#' @param second A list - \strong{'key'} is \code{'\$pos'} (i.e. breakpoint pos of geneB); \code{'second\$pos\$transcript'} is a data.frame object
#'        (only one row) that collects evaluation of breakpoint pos for geneB transcript; \code{'second\$pos\$select_region'} is a data.frame object
#'        (5utr-cds-utr3 annotation) for plotting and it is a constant variable for different \code{'\$pos'} values.
#' @param breakpoint A data.frame with three columns (e.g. \code{'breakpoint[i,1]'} - breakpoint pos of geneA; \code{'breakpoint[i,2]'} - breakpoint pos of geneB; 
#'        \code{'breakpoint[i,3]'} is the number of samples with such a fusion).
#' @param zoom A numeric value indicates zoom degree of plotting.
#'
#' @return NULL
#'
#' @export
plot_arrow_overview <- function(collect, first, second, breakpoint, zoom) {
	A1_xy = collect$A1_xy
	B1_xy = collect$B1_xy
	zoom = zoom/50;
	bezier_control_point_offset = 40*zoom; #// Control the degree of bezier curve (can be adjusted)
			
	#// Process curve coordinate (x1, y, x2, y2)
	if ( nrow(breakpoint) > 0 ) { #// make sure 'breakpoint' data.frame is not empty
		for (i in 1:length(breakpoint[,1])) {
			curve_coordinate = .plot_curve(A1_xy, first, B1_xy, second, breakpoint[i,]);
			#// 'curve_coordinate' is a vector class with three elements (i.e. x_pos_gene_upstream: x-coordinate of geneA breakpoint, 
			#// x_pos_gene_downstream: x-coordinate of geneB breakpoint, y_pos: y-coordinate of geneA and geneB breakpoint)	
			if (! is.null(curve_coordinate) ) {
				#// NOTE: this value should be identical to the setting 'plotOutput("chimerics_up", height=50)' in ui.R
				curve_coordinate["y_pos"] = 50*zoom;
				#// Draw the bezier curve between transcripts						
				grid::pushViewport(grid::viewport(xscale = c(0, grDevices::dev.size(units = "px")[1]), yscale = c(grDevices::dev.size(units = "px")[2], 0)));
				grid::grid.bezier(rep(c(curve_coordinate["x_pos_gene_upstream"], curve_coordinate["x_pos_gene_downstream"]), each = 2), #// x coordinates
					c(curve_coordinate["y_pos"], curve_coordinate["y_pos"] - bezier_control_point_offset, curve_coordinate["y_pos"] - bezier_control_point_offset, 
					curve_coordinate["y_pos"] ), default.units = "native", arrow=grid::arrow(angle=50, length=grid::unit(0.05, "inches"), ends="last", type="closed"), 
					gp = grid::gpar(col="red", fill="red", lwd=0.8));

				number_position_x = curve_coordinate["x_pos_gene_upstream"] + (curve_coordinate["x_pos_gene_downstream"] - curve_coordinate["x_pos_gene_upstream"]) / 2
				number_position_y = curve_coordinate["y_pos"] - bezier_control_point_offset + 4*zoom;
				#// 'supporting_reads_text' is the number of samples supporting a given breakpoint combination
				supporting_reads_text = paste("(", breakpoint[i,3], ")", sep="")
				grid::grid.text(supporting_reads_text, x = number_position_x / grDevices::dev.size(units = "px")[1], y = 1 - number_position_y / grDevices::dev.size(units = "px")[2],
					vp = grid::viewport(xscale = c(0, grDevices::dev.size(units = "px")[1]), yscale = c(grDevices::dev.size(units = "px")[2], 0)),
					gp = grid::gpar(fontsize = 7))
			}
		}
	}
}

#' Plot transcript-model, chromosome ideogram and multiple curved lines that link breakpoints of geneA and geneB transcripts (plot downloading version)
#'
#' @description Plot transcript-model, chromosome ideogram and multiple curved lines that link breakpoints of geneA and geneB transcripts together in 'overview' tab-panel
#'
#' @param first A list - \strong{'key'} is \code{'\$pos'} (i.e. breakpoint pos of geneA); \code{'first\$pos\$transcript'} is a data.frame object
#'        (only one row) that collects evaluation of breakpoint pos for geneA transcripts; \code{'first\$pos\$select_region'} is a data.frame object
#'        (5utr-cds-utr3 annotation) for plotting and it is a constant variable for different \code{'\$pos'} values.
#' @param first_name A string - symbol name of geneA (e.g. 'TMPRSS2').
#' @param second A list - \strong{'key'} is \code{'\$pos'} (i.e. breakpoint pos of geneB); \code{'second\$pos\$transcript'} is a data.frame object
#'        (only one row) that collects evaluation of breakpoint pos for geneB transcripts; \code{'second\$pos\$select_region'} is a data.frame object
#'        (5utr-cds-utr3 annotation) for plotting and it is a constant variable for different \code{'\$pos'} values.
#' @param second_name A string - symbol name of geneB (e.g. 'ERG').
#' @param cytoband A GRange object for plotting cytoband of chromosome ideogram.
#' @param breakpoint A data.frame with three columns (e.g. \code{'breakpoint[i,1]'} - breakpoint pos of geneA; \code{'breakpoint[i,2]'} - breakpoint pos of geneB;
#' @param offset A numeric value - extend the partner gene interval (default: 2000, e.g. [start-offset, end+offset]).
#'
#' @export
plot_separate_overview_download <- function(first, first_name, second, second_name, cytoband, breakpoint, offset=2000) {
	chrom_f = first[[1]]$transcript$Chrom[1] #// chromosome name of first: Ensembl
	chrom_s = second[[1]]$transcript$Chrom[1] #// chromosome name of second: Ensembl

	#// Set plotting tracks for geneA
	grTrack_f = Gviz::GeneRegionTrack(first[[1]]$select_region, chromosome=chrom_f, showId=T, stacking="squish", geneSymbols=F, transcriptAnnotation="transcript", just.group="below",
			fontsize=5, fontcolor="black", fontcolor.group="black", showTitle=F, col="black", col.line="black", fill="green", col.symbol="black", lwd=0.1, lwd.border=1,
			lex=1, cex=0.6, cex.axis=0.6, background.title="transparent", background.panel="transparent", stackHeight=0.8, min.width=0.1, min.height=3, min.distance=0);
	first_vis_s = min(first[[1]]$transcript$GStart) #// start position of geneA for visualization
	first_vis_e = max(first[[1]]$transcript$GEnd) #// end position of geneA for visualization 

	#// Set plotting tracks for geneB
	grTrack_s = Gviz::GeneRegionTrack(second[[1]]$select_region, chromosome=chrom_s, showId=T, stacking="squish", geneSymbols=F, transcriptAnnotation="transcript",
			just.group="below", fontsize=5, fontcolor="black", fontcolor.group="black", showTitle=F, col="black", col.line="black", fill="orange", col.symbol="black", lwd=0.1,
			lwd.border=1, lex=1, cex=0.6, cex.axis=0.6, background.title="transparent", background.panel="transparent", stackHeight=0.8, min.width=0.1, min.height=3, min.distance=0);
	second_vis_s = min(second[[1]]$transcript$GStart) #// start position of geneB for visualization
	second_vis_e = max(second[[1]]$transcript$GEnd) #// end position of geneB for visualization

	#// adjust font size of transcript_id (it will probably be adjusted by customers in UI in further version)
	if ( length(first[[1]]$transcript$TXNAME) > 0 &&  length(first[[1]]$transcript$TXNAME) <= 5 ) { 
		grTrack_f@dp@pars$fontsize.group = 11; 
	} else if ( length(first[[1]]$transcript$TXNAME) > 5 &&  length(first[[1]]$transcript$TXNAME) <= 15 ) { 
		grTrack_f@dp@pars$fontsize.group = 9; 
	} else if ( length(first[[1]]$transcript$TXNAME) > 15 &&  length(first[[1]]$transcript$TXNAME) <= 20 ) {
		grTrack_f@dp@pars$fontsize.group = 8; 
	} else if ( length(first[[1]]$transcript$TXNAME) > 20 &&  length(first[[1]]$transcript$TXNAME) <= 25 ) {
		grTrack_f@dp@pars$fontsize.group = 7;
	} else if ( length(first[[1]]$transcript$TXNAME) > 25 &&  length(first[[1]]$transcript$TXNAME) <= 30 ) {
		grTrack_f@dp@pars$fontsize.group = 5;
	} else {
		grTrack_f@dp@pars$fontsize.group = 3;
	}
	if ( length(second[[1]]$transcript$TXNAME) > 0 &&  length(second[[1]]$transcript$TXNAME) <= 5 ) {
		grTrack_s@dp@pars$fontsize.group = 11;
	} else if ( length(second[[1]]$transcript$TXNAME) > 5 &&  length(second[[1]]$transcript$TXNAME) <= 15 ) {
		grTrack_s@dp@pars$fontsize.group = 9;
	} else if ( length(second[[1]]$transcript$TXNAME) > 15 &&  length(second[[1]]$transcript$TXNAME) <= 20 ) {
		grTrack_s@dp@pars$fontsize.group = 8;
	} else if ( length(second[[1]]$transcript$TXNAME) > 20 &&  length(second[[1]]$transcript$TXNAME) <= 25 ) {
		grTrack_s@dp@pars$fontsize.group = 7;
	} else if ( length(second[[1]]$transcript$TXNAME) > 25 &&  length(second[[1]]$transcript$TXNAME) <= 30 ) {
		grTrack_s@dp@pars$fontsize.group = 5;
	} else {
		grTrack_s@dp@pars$fontsize.group = 3;
	}

	#// set the chromosome axis coordinate
	axis = Gviz::GenomeAxisTrack(name="Axis", labelPos="alternating", fontcolor="black", littleTicks=F, cex=0.5, cex.id=0.5, cex.axis=0.5,
			fontsize=14, col="black", distFromAxis=0.4, lwd=1.4, lwd.border=1, min.width=0.1, min.height=3, min.distance=0);

	#// define plot region [start-offset, end+offset]
	#// plot 'blank' 'gene+axis', 'title' and 'chromosome' tracks in proportion as 0.95 : 9 : 0.2 : 0.7 (can be adjusted)
	nf = grid::grid.layout(nrow=4, ncol=2, widths=grid::unit(c(1, 1), "null"), heights = grid::unit(c(0.90, 9, 0.2, 0.7), "null"))
	grid::grid.newpage(); #// start a new page
	grid::pushViewport(grid::viewport(layout = nf)); #// divide the page to different sections (nf specifies the section composite)

	#// plot 'geneA' and 'axis' tracks in proportion as 1 : 0.2 (can be adjusted)
	grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = 1)) #// load on the section (1, 1)
	if ( first[[1]]$transcript$Strand[1] == "+" ) {
		A1_xy = Gviz::plotTracks(c(grTrack_f, axis), chromosome=chrom_f, from=first_vis_s-offset, to=first_vis_e+offset,
				margin=1, innerMargin=1, sizes=c(1, 0.2), add=T);
	} else {
		A1_xy = Gviz::plotTracks(c(grTrack_f, axis), chromosome=chrom_f, from=first_vis_s-offset, to=first_vis_e+offset,
				margin=1, innerMargin=1, sizes=c(1, 0.2), reverseStrand=TRUE, add=T);
	}
	grid::popViewport(1); #// close the section (1, 1)
	#// plot 'geneB' and 'axis' tracks in proportion as 1 : 0.2 (can be adjusted)
	grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = 2)) #// load on the section (1, 2)
	if ( second[[1]]$transcript$Strand[1] == "+" ) {
		B1_xy = Gviz::plotTracks(c(grTrack_s, axis), chromosome=chrom_s, from=second_vis_s-offset, to=second_vis_e+offset, 
				margin=1, innerMargin=1, sizes=c(1, 0.2), add=T);
	} else {
		B1_xy = Gviz::plotTracks(c(grTrack_s, axis), chromosome=chrom_s, from=second_vis_s-offset, to=second_vis_e+offset, 
				margin=1, innerMargin=1, sizes=c(1, 0.2), reverseStrand=TRUE, add=T);
	}
	grid::popViewport(1); #// close the section (1, 2)

	#// add sub-title for geneA
	grid::pushViewport(grid::viewport(layout.pos.row = 3, layout.pos.col = 1)); #// load on the section (2, 1)
	grid::grid.text(label=paste("Upstream gene: ", first_name, " (", first[[1]]$transcript$Strand[1], ")", sep=""),
						x=grid::unit(0.5, "npc"), y=grid::unit(0.85, "npc"), gp=grid::gpar(fontsize=10, fontface="bold", col="black"));
	grid::popViewport(1); #// close the section (2, 1)
	#// add sub-title for geneB
	grid::pushViewport(grid::viewport(layout.pos.row = 3, layout.pos.col = 2)); #// load on the section (2, 2)
	grid::grid.text(label=paste("Downstream gene: ", second_name, " (", second[[1]]$transcript$Strand[1], ")", sep=""),
						x=grid::unit(0.5, "npc"), y=grid::unit(0.85, "npc"), gp=grid::gpar(fontsize=10, fontface="bold", col="black"));
	grid::popViewport(1); #// close the section (2, 2)

	if ( chrom_f == chrom_s ) { #// if both geneA and geneB at the same chromosome
		#// plot cytoband of the chromosome
		names(cytoband) = chrom_f;
		grid::pushViewport(grid::viewport(layout.pos.row = 4, layout.pos.col = NULL)); #// load on the section (3, null)
		if ( first[[1]]$transcript$GStart[1] < second[[1]]$transcript$GStart[1] ) {
			Gviz::plotTracks(cytoband, chromosome=chrom_f, from=first_vis_s, to=second_vis_e, add=T, main="", cex.main=0.6, fontsize=15, fontface.main=2, margin=1, innerMargin=1)
		} else {
			Gviz::plotTracks(cytoband, chromosome=chrom_s, from=second_vis_s, to=first_vis_e, add=T, main="", cex.main=0.6, fontsize=15, fontface.main=2, margin=1, innerMargin=1)
		}
		grid::popViewport(1);
		#// plot fusion gene title
		grid::pushViewport(grid::viewport(layout.pos.row = 4, layout.pos.col = NULL));
		grid::grid.text(label=paste("Fusion: ", first_name, "-", second_name, sep=""), x=grid::unit(0.5, "npc"), y=grid::unit(0.85, "npc"),
								gp=grid::gpar(fontsize=10, fontface="bold", col="black"));
		grid::popViewport(1); #// close the section (3, null)
	} else { #// if geneA and geneB at the different chromosomes
		cytoband_A = cytoband;	names(cytoband_A) = chrom_f;
		cytoband_B = cytoband;	names(cytoband_B) = chrom_s;
		#// plot cytoband of geneA chromosome
		grid::pushViewport(grid::viewport(layout.pos.row = 4, layout.pos.col = 1)); #// load on the section (3, 1)
		Gviz::plotTracks(cytoband_A, chromosome=chrom_f, from=first_vis_s, to=first_vis_e, add=T, cex.main=0.6, fontsize=15, fontface.main=2, margin=1, innerMargin=1);
		grid::popViewport(1); #// close the section (3, 1)
		#// plot cytoband of geneB chromosome
		grid::pushViewport(grid::viewport(layout.pos.row = 4, layout.pos.col = 2)); #// load on the section (3, 2)
		Gviz::plotTracks(cytoband_B, chromosome=chrom_s, from=second_vis_s, to=second_vis_e, add=T, cex.main=0.6, fontsize=15, fontface.main=2, margin=1, innerMargin=1);
		grid::popViewport(1); #// close the section (3, 2)
		#// plot fusion gene title
		grid::pushViewport(grid::viewport(layout.pos.row = 4, layout.pos.col = NULL)); #// load on the section (3, null)
		grid::grid.text(label=paste("Fusion: ", first_name, "-", second_name, sep=""), x=grid::unit(0.5, "npc"), y=grid::unit(0.85, "npc"),
								gp=grid::gpar(fontsize=10, fontface="bold", col="black"))
		grid::popViewport(1); #// close the section (3, null)
	}

	bezier_control_point_offset = 40; #// Control the degree of bezier curve (can be adjusted)
	#// Process curve coordinate (x1, y, x2, y2)
	if ( nrow(breakpoint) > 0 ) { #// make sure 'breakpoint' data.frame is not empty
		for (i in 1:length(breakpoint[,1])) {
			curve_coordinate = .plot_curve(A1_xy, first, B1_xy, second, breakpoint[i,]);
			#// 'curve_coordinate' is a vector class with three elements (i.e. x_pos_gene_upstream: x-coordinate of geneA breakpoint,
			#// x_pos_gene_downstream: x-coordinate of geneB breakpoint, y_pos: y-coordinate of geneA and geneB breakpoint)
			if (! is.null(curve_coordinate) ) {
				#// NOTE: this value should be identical to the setting 'plotOutput("chimerics_up", height=50)' in ui.R
				curve_coordinate["y_pos"] = 50;
				#// Draw the bezier curve between transcripts
				grid::pushViewport(grid::viewport(xscale = c(0, grDevices::dev.size(units = "px")[1]), yscale = c(grDevices::dev.size(units = "px")[2], 0)));
				grid::grid.bezier(rep(c(curve_coordinate["x_pos_gene_upstream"], curve_coordinate["x_pos_gene_downstream"]), each = 2), #// x coordinates
					c(curve_coordinate["y_pos"], curve_coordinate["y_pos"] - bezier_control_point_offset, curve_coordinate["y_pos"] - bezier_control_point_offset,
					curve_coordinate["y_pos"] ), default.units = "native", arrow=grid::arrow(angle=50, length=grid::unit(0.05, "inches"), ends="last", type="closed"),
					gp = grid::gpar(col="red", fill="red", lwd=0.4));

				number_position_x = curve_coordinate["x_pos_gene_upstream"] + (curve_coordinate["x_pos_gene_downstream"] - curve_coordinate["x_pos_gene_upstream"]) / 2
				number_position_y = curve_coordinate["y_pos"] - bezier_control_point_offset + 4;
				#// 'supporting_reads_text' is the number of samples supporting a given breakpoint combination
				supporting_reads_text = paste("(", breakpoint[i,3], ")", sep="")
				grid::grid.text(supporting_reads_text, x = number_position_x / grDevices::dev.size(units = "px")[1], y = 1 - number_position_y / grDevices::dev.size(units = "px")[2],
					vp = grid::viewport(xscale = c(0, grDevices::dev.size(units = "px")[1]), yscale = c(grDevices::dev.size(units = "px")[2], 0)),
					gp = grid::gpar(fontsize = 7))
				#// Add dash lines showing the breakpoint of upstream partner
				grid::grid.lines(c(curve_coordinate["x_pos_gene_upstream"], curve_coordinate["x_pos_gene_upstream"]), #// x positions
					c(curve_coordinate["y_pos"] - 5, 450), #// y positions
					default.units = "native", arrow=NULL, gp = grid::gpar(col="black", lty=3, lwd=0.3));
				#// Add dash lines showing the breakpoint of downstream partner
				grid::grid.lines(c(curve_coordinate["x_pos_gene_downstream"], curve_coordinate["x_pos_gene_downstream"]), #// x positions
					c(curve_coordinate["y_pos"] - 5, 450), #// y positions
					default.units = "native", arrow=NULL, gp = grid::gpar(col="black", lty=3, lwd=0.3));
			}
		}
	}
}

