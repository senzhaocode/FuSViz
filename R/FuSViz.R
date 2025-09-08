#' FuSViz
#'
#' Make a HTML widget for an embeddable interactive genome visualization javascript library (igv.js)
#'

configure <- NULL

#' FuSViz widget
#'
#' @description Defined a main class to build widget
#' 
#' @param genomeName A string shows human genome reference version, either \code{hg19} or \code{hg38}.
#' @param trackHeight A number defines initial height of a track (e.g. default value: \code{300}).
#' @param displayMode Annotation track display mode (e.g. \code{COLLAPSED}, \code{EXPANDED}, \code{SQUISHED}).
#' @param initialLocus A string defines initial genomic location when widget loading is launched (e.g. default value: \code{TP53}).
#' @param width,height Must be a valid CSS unit (like \code{100\%},
#'        \code{400px}, \code{auto}) or a number, which will be coerced to a string and have \code{px} appended.
#' @param elementId A string represents the name of a HTML element id used for DOM manipulation.
#'
#' @export
FuSViz <- function(genomeName, trackHeight=300, displayMode="EXPANDED", initialLocus=NULL, width=NULL, height=NULL, elementId=NULL) {
	# parameters "genomeName, trackHeight, displayMode and initialLocus" are related to setting in igv.js
	print("Loading parameters from FuSViz class.");
	if ( is.null(genomeName) ) { genomeName = ""; }
	# create a list 
	x <- list(genomeName=genomeName, displayMode=displayMode, trackHeight=trackHeight);

	# binding to FuSViz.js using htmlwidgets::createWidget
	htmlwidgets::createWidget(name = 'FuSViz', x, width = width, height = height, elementId = elementId, package = 'FuSViz');
	
}

#' Shiny bindings for FuSViz
#'
#' @description Output element embedded in ui { }
#'
#' @param outputId Manipulate as a DOM element id 
#' @param width,height Must be a valid CSS unit (like \code{100\%},
#'   \code{400px}, \code{auto}) or a number, which will be coerced to a string and have \code{px} appended.
#'
#' @export
FuSVizOutput <- function(outputId, width = '100%', height = '400px') {
	htmlwidgets::shinyWidgetOutput(outputId, 'FuSViz', width, height, package = 'FuSViz')
}

#' Shiny bindings for FuSViz
#'
#' @description Render igv browser element in server { }
#'
#' @param expr An expression that generates a FuSViz.
#' @param env The environment in which to evaluate \code{expr}.
#' @param quoted Is \code{expr} a quoted expression (with \code{FALSE} or \code{TRUE})? This
#'        is useful if you want to save an expression in a variable.
#'
#' @export
renderFuSViz <- function(expr, env = parent.frame(), quoted = FALSE) {
	if (!quoted){ expr <- substitute(expr) }
	htmlwidgets::shinyRenderWidget(expr, FuSVizOutput, env, quoted = TRUE)
}

#' Remove track with replicate names
#'
#' @description Remove track with replicate names
#'
#' @param session An environment object that is used to access information and functionality by shiny.
#' @param name A character vector represents track name.
#'
#' @export
DuplicateTrackRemove <- function(session, name) {
	connect.igvjs <- list(addtracks=name)
	session$sendCustomMessage("DuplicateTrackRemove", connect.igvjs)
}

#' Remove track added by users
#'
#' @description Remove track added by users
#'
#' @param session An environment object that is used to access information and functionality by shiny.
#'
#' @export
UserTracksRemove <- function(session) {
	DuplicateTrackRemove(session, configure)
	configure <- NULL
}

#' Load a bed format track
#'
#' @description Load a bed format track representing SV breakpoint interval
#'
#' @param session An environment object that is used to access information and functionality by shiny.
#' @param datastr A data.frame with at least three columns (e.g. \code{chr}, \code{start}, \code{end}, from the fourth to the last column is not essential).
#' @param name A string - track name (e.g. \code{RNA breakpoint track}).
#' @param color CSS color value for track feature (e.g. \code{#ff0000} or \code{rgb(100,0,100)}).
#' @param trackHeight Initial height of track viewport in pixels (e.g. \code{100}).
#' @param displayMode Annotation track display mode (e.g. \code{COLLAPSED}, \code{EXPANDED}, \code{SQUISHED}).
#' @param SameNameDel A logical value controls whether to add a track if its name has been existed (e.g. \code{FALSE} or \code{TRUE}).
#'
#' @export
TrackinBed <- function(session, datastr, name, color="gray", trackHeight=100, displayMode="EXPANDED", SameNameDel=TRUE) {
	# check contents of datastr
	print("load user-defined track in bed format")
	if ( SameNameDel ) { #// if userdefined track has been loaded, please remove the previous one
		DuplicateTrackRemove(session, name)
	}
	configure <- unique(c(configure, name))

	if (! is(datastr$chr, "character") ) { stop("datastr$chr not character!"); }
	if (! is(datastr$start, "numeric") ) { stop("datastr$start not numeric!"); }
	if (! is(datastr$end, "numeric") ) { stop("datastr$end not numeric!"); }

	datastr = datastr[order(datastr$chr, datastr$start), ] # sort by chromosome and position
	connect.igvjs <- list(name=name, datastr=jsonlite::toJSON(datastr), color=color, trackHeight=trackHeight, displayMode=displayMode) # convert datastr to JSON format
	session$sendCustomMessage("TrackinBed", connect.igvjs)
}

#' Load a bedgraph format track
#'
#' @description Load a bedgraph format track - count the freq of SV breakpoints
#'
#' @param session An environment object that is used to access information and functionality by shiny.
#' @param datastr A data.frame with four columns (e.g. \code{chr}, \code{start} and \code{value}).
#' @param name A string - track name (e.g. \code{RNA breakpoint freq track}).
#' @param color CSS color value for track feature (e.g. \code{#ff0000} or \code{rgb(100,0,100)}).
#' @param trackHeight Initial height of track viewport in pixels (e.g. \code{100}).
#' @param autoscale A logical value of y-axis scale (e.g. \code{FALSE} or \code{TRUE}).
#' @param min A number represents minimum value of y-axis.
#' @param max A number represents maximum value of y-axis.
#' @param displayMode Display mode of annotation track (e.g. \code{COLLAPSED}, \code{EXPANDED}, \code{SQUISHED}).
#' @param SameNameDel A logical value controls whether to add track if its name has been existed (e.g. \code{FALSE} or \code{TRUE}).
#'
#' @export
TrackinBedGraph <- function(session, datastr, name, color="gray", trackHeight=50, autoscale=TRUE, min=0, max=50, displayMode="EXPANDED", SameNameDel=TRUE) {
	print("load user-defined track in bedgraph format");
	if ( SameNameDel ) { #// if userdefined track has been loaded, please remove the previous one
		DuplicateTrackRemove(session, name)
	}
	configure <- unique(c(configure, name))

	if ( ncol(datastr) > 4 ) { stop("column number for bedgraph format <= 4!"); } # make sure the col number <= 4
	if (! is(datastr$chr, "character") ) { stop("datastr$chr not character!"); }
	if (! is(datastr$start, "numeric") ) { stop("datastr$start not numeric!"); }
	if (! is(datastr$end, "numeric") ) { stop("datastr$end not numeric!"); }
	if (! is(datastr$value, "numeric") ) { stop("datastr$value not numeric!"); }
	if ( autoscale == FALSE ) { stopifnot(!is.na(max)) }

	datastr = datastr[order(datastr$chr, datastr$start), ] # sort by chromosome and position
	connect.igvjs <- list(name=name, datastr=jsonlite::toJSON(datastr), color=color, trackHeight=trackHeight, 
						  autoscale=autoscale, min=min, max=max, displayMode=displayMode) # convert datastr to JSON format
	session$sendCustomMessage("TrackinBedGraph", connect.igvjs)
}

#' Load a bedpe format track
#'
#' @description Load a bedpe format track - count the distribution of SVs (type including: DEL, DUP, INS and INV; BND is not shown due to translocation)
#'
#' @param session An environment object that is used to access information and functionality by shiny.
#' @param datastr A data.frame with nine columns, e.g. \code{chr1}, \code{start1}, \code{end1}, \code{chr2}, \code{start2}, \code{end2}, 
#'        \code{name}(multiple samples separated by ";"), \code{score}(number of samples), \code{type}.
#' @param name A string - track name (e.g. \code{DNA bedpe track}).
#' @param color CSS color value for track feature (e.g. \code{#ff0000} or \code{rgb(100,0,100)}).
#' @param thickness A number represents thickness of curve (e.g. \code{0.5}).
#' @param trackHeight Initial height of track viewport in pixels (e.g. \code{100}).
#' @param logScale A logical value (e.g. \code{FALSE} or \code{TRUE}).
#' @param alpha A number (0-1) that controls arc of curve (e.g. \code{0}).
#' @param min A number represents minimum value of y-axis.
#' @param displayMode Display mode of annotation track (e.g. \code{COLLAPSED}, \code{EXPANDED}, \code{SQUISHED}).
#' @param SameNameDel A logical value controls whether to add track if its name has been existed (e.g. \code{FALSE} or \code{TRUE}).
#'
#' @export
TrackinBedPe <- function(session, datastr, name, color="blue", thickness=0.5, trackHeight=200, logScale=TRUE, alpha=0, min=0, displayMode="EXPANDED", SameNameDel=TRUE) {
	print("load user-defined track in bedpe format")
	if ( SameNameDel ) { #// if userdefined track has been loaded, please remove the previous one first
		DuplicateTrackRemove(session, name)
	}
	configure <- unique(c(configure, name))

	# check contents of datastr
	if (! is(datastr$chr1, "character") ) { stop("datastr$chr1 not character!"); }
	if (! is(datastr$start1, "numeric") ) { stop("datastr$start1 not numeric!"); }
	if (! is(datastr$end1, "numeric") ) { stop("datastr$end1 not numeric!"); }
	if (! is(datastr$chr2, "character") ) { stop("datastr$chr2 not character!"); }
	if (! is(datastr$start2, "numeric") ) { stop("datastr$start2 not numeric!"); }
	if (! is(datastr$end2, "numeric") ) { stop("datastr$end2 not numeric!"); }
	if (! is(datastr$name, "character") ) { stop("datastr$name not character!"); }
	if (! is(datastr$score, "numeric") ) { stop("datastr$score not numeric!"); }
	if (! is(datastr$type, "character") ) { stop("datastr$type not character!"); }

	datastr = datastr[order(datastr$chr1, datastr$start1), ] # sort by chromosome and position
	connect.igvjs <- list(name=name, datastr=jsonlite::toJSON(datastr), color=color, trackHeight=trackHeight,
			logScale=logScale, min=min, alpha=alpha, thickness=thickness, displayMode=displayMode); # convert datastr to JSON format
	session$sendCustomMessage("TrackinBedPe", connect.igvjs)
}

#' Load a BAF plot track
#'
#' @description Load a compressed VCF (*.vcf.gz) file - import the BAF plot track
#'
#' @param session An environment object that is used to access information and functionality by shiny.
#' @param filedata A data.frame of imported compressed vcf files with four columns (e.g. \code{name}, \code{size}, \code{type} and \code{datapath}).
#' @param local A string - the tmp path of shiny app (e.g. \code{/var/folders/tt/qx4jdszs5sj533vn92lrcx4m0000gn/T//Rtmpr6bFPt}).
#' @param BinSize Bin size value for BAF plot, available values: 10000, 100000, 1000000, 10000000 or 100000000 (e.g. \code{default: 1000000}).
#' @param trackHeight Initial height of track viewport in pixels (e.g. \code{100}).
#' @param displayMode Display mode of annotation track (e.g. \code{COLLAPSED}, \code{EXPANDED}, \code{SQUISHED}).
#'
#' @export
FileBAF <- function(session, filedata, local, BinSize, trackHeight=200, displayMode="EXPANDED", SameNameDel=TRUE) {
	print("load a BAF-plot file in vcf format")
	for ( name in filedata$name ) {
		if ( SameNameDel ) { #// if userdefined track has been loaded, please remove the previous one first
			DuplicateTrackRemove(session, name)
		}
		configure <- unique(c(configure, name))
	}

	connect.igvjs <- list(filedata=jsonlite::toJSON(filedata), local=local, BinSize=BinSize, trackHeight=trackHeight, displayMode=displayMode); # convert datastr to JSON format
	session$sendCustomMessage("FileBAF", connect.igvjs)
}

#' Update a BAF plot track
#'
#' @description Update BAF plot track with more customized settings
#'
#' @param session An environment object that is used to access information and functionality by shiny.
#' @param BinSize Bin size value for BAF plot, available values: 10000, 100000, 1000000, 10000000 or 100000000 (e.g. \code{default: 1000000}).
#' @param trackHeight Initial height of track viewport in pixels (e.g. \code{100}).
#' @param displayMode Display mode of annotation track (e.g. \code{COLLAPSED}, \code{EXPANDED}, \code{SQUISHED}).
#'
#' @export
UpdateBAF <- function(session, BinSize, trackHeight=200, displayMode="EXPANDED") {
	print("Update the BAF plot track")

	connect.igvjs <- list(BinSize=BinSize, trackHeight=trackHeight, displayMode=displayMode);
	session$sendCustomMessage("UpdateBAF", connect.igvjs)
}

#' Load a splice junction track
#'
#' @description Load junction bed (*.SJ.out.bed.gz and *.SJ.out.bed.gz.tbi) files - import splice junction tracks
#'
#' @param session An environment object that is used to access information and functionality by shiny.
#' @param filedata A data.frame of imported splice junction bed files with four columns (e.g. \code{name}, \code{size}, \code{type} and \code{datapath}).
#' @param local A string - the tmp path of shiny app (e.g. \code{/var/folders/tt/qx4jdszs5sj533vn92lrcx4m0000gn/T//Rtmpr6bFPt}).
#' @param unique Junction must be supported by at least this number of uniquely-mapped reads (e.g. \code{10}).
#' @param total Junction must be supported by at least this number of uniquely-mapped + multi-mapped reads (e.g. \code{10}).
#' @param percet (uniquely-mapped reads)/(total reads) must be <= this threshold (e.g. \code{0-1}).
#' @param overhang Mininum spliced alignment overhang in base pairs (e.g. \code{20}).
#' @param Selectcolor Splice junction color, available values: 'numUniqueReads', 'numReads', 'isAnnotatedJunction', 'strand', 'motif' (e.g. \code{default: numUniqueReads}).
#' @param Selectthick Splice junction line thickness, available values: 'numUniqueReads', 'numReads', 'isAnnotatedJunction' (e.g. \code{default: numUniqueReads}).
#' @param Selectcurve Splice junction curve height, available values: 'random', 'distance', 'thickness' (e.g. \code{default: random}).
#' @param trackHeight Initial height of track viewport in pixels (e.g. \code{100}).
#' @param displayMode Display mode of annotation track (e.g. \code{COLLAPSED}, \code{EXPANDED}, \code{SQUISHED}).
#'
#' @export
FileinJunction <- function(session, filedata, local, unique, total, percet, overhang, Selectcolor, Selectthick, Selectcurve, trackHeight=200, displayMode="EXPANDED", SameNameDel=TRUE) {
	print("load splice junction file in bed format")
	for ( name in filedata$name ) {
		if ( SameNameDel ) { #// if userdefined track has been loaded, please remove the previous one first
			DuplicateTrackRemove(session, name)
		}
		configure <- unique(c(configure, name))
	}

	connect.igvjs <- list(filedata=jsonlite::toJSON(filedata), local=local, unique=unique, total=total, percet=percet, overhang=overhang, Selectcolor=Selectcolor, Selectthick=Selectthick, 
			Selectcurve=Selectcurve, trackHeight=trackHeight, displayMode=displayMode); # convert datastr to JSON format
	session$sendCustomMessage("FileinJunction", connect.igvjs)
}

#' Update a splice junction track
#'
#' @description Update splice junction track with more customized settings
#'
#' @param session An environment object that is used to access information and functionality by shiny.
#' @param unique Junction must be supported by at least this number of uniquely-mapped reads (e.g. \code{10}).
#' @param total Junction must be supported by at least this number of uniquely-mapped + multi-mapped reads (e.g. \code{10}).
#' @param percet (uniquely-mapped reads)/(total reads) must be <= this threshold (e.g. \code{0-1}).
#' @param overhang Mininum spliced alignment overhang in base pairs (e.g. \code{20}).
#' @param Selectcolor Splice junction color, available values: 'numUniqueReads', 'numReads', 'isAnnotatedJunction', 'strand', 'motif' (e.g. \code{default: numUniqueReads}).
#' @param Selectthick Splice junction line thickness, available values: 'numUniqueReads', 'numReads', 'isAnnotatedJunction' (e.g. \code{default: numUniqueReads}).
#' @param Selectcurve Splice junction curve height, available values: 'random', 'distance', 'thickness' (e.g. \code{default: random}).
#' @param trackHeight Initial height of track viewport in pixels (e.g. \code{100}).
#' @param displayMode Display mode of annotation track (e.g. \code{COLLAPSED}, \code{EXPANDED}, \code{SQUISHED}).
#'
#' @export
UpdateJunction <- function(session, unique, total, percet, overhang, Selectcolor, Selectthick, Selectcurve, trackHeight=200, displayMode="EXPANDED") {
	print("Update splice junction track")

	connect.igvjs <- list(unique=unique, total=total, percet=percet, overhang=overhang, Selectcolor=Selectcolor, Selectthick=Selectthick, 
			Selectcurve=Selectcurve, trackHeight=trackHeight, displayMode=displayMode);
	session$sendCustomMessage("UpdateJunction", connect.igvjs)
}

#' Load a seg format track
#'
#' @description Load a seg format track - count the DEL and DUP events of DNA SVs (small copy number variations)
#'
#' @param session An environment object that is used to access information and functionality by shiny.
#' @param datastr A data.frame with five columns (e.g. \code{chr}, \code{start}, \code{end}, \code{value} and \code{sample}).
#' @param name A string - track name (e.g. \code{DNA seg track}).
#' @param trackHeight Initial height of track viewport in pixels (e.g. \code{100}).
#' @param isLog A logical value (e.g. \code{FALSE} or \code{TRUE} - if segment values are 2*log2(copyNumber/2)).
#' @param displayMode Display mode of annotation track (e.g. \code{COLLAPSED}, \code{EXPANDED}, \code{SQUISHED}).
#' @param SameNameDel A logical value controls whether to add track if its name has been existed (e.g. \code{FALSE} or \code{TRUE}).
#'
#' @export
TrackinSeg <- function(session, datastr, name, trackHeight=50, isLog=TRUE, displayMode="EXPANDED", SameNameDel=TRUE) {
	print("load user-defined track in seg format")
	if ( SameNameDel ) { #// if userdefined track has been loaded, please remove the previous one first
		DuplicateTrackRemove(session, name)
	}
	configure <- unique(c(configure, name))

	# check contents of datastr
	if (! is(datastr$chr, "character") ) { stop("datastr$chr not character!"); }
	if (! is(datastr$start, "numeric") ) { stop("datastr$start not numeric!"); }
	if (! is(datastr$end, "numeric") ) { stop("datastr$end not numeric!"); }
	if (! is(datastr$value, "numeric") ) { stop("datastr$value not numeric!"); }
	if (! is(datastr$sample, "character") ) { stop("datastr$sample not character!"); }

	datastr = datastr[order(datastr$chr, datastr$start), ] # sort by chromosome and position
	connect.igvjs <- list(name=name, datastr=jsonlite::toJSON(datastr), trackHeight=trackHeight, isLog=isLog, displayMode=displayMode) # convert datastr to JSON format
	session$sendCustomMessage("TrackinSeg", connect.igvjs)
}

#' Load a BAM format from url
#'
#' @description Load a BAM format from url
#'
#' @param session An environment object that is used to access information and functionality by shiny.
#' @param bam An url link of BAM file.
#' @param bamindex An url link of BAM index file.
#' @param trackHeight Initial height of track viewport in pixels (e.g. \code{200}).
#'
#' @export
TrackinBAM <- function(session, bam, bamindex, trackHeight=200) {
	print("Load http:// or s3:// BAM and CRAM url");

	connect.igvjs <- list(bam=bam, bamindex=bamindex, trackHeight=trackHeight)
	session$sendCustomMessage("TrackinBAM", connect.igvjs)
}

#' Load gene track in offline mode
#'
#' @description Load gene track in offline mode (only for hg19_offline or hg38_offline)
#'
#' @param session An environment object that is used to access information and functionality by shiny.
#' @param version A string shows genome reference version (either 'hg19' or 'hg38')
#' @param name Initial height of track viewport in pixels (e.g. \code{100}).
#' @param trackHeight Initial height of track viewport in pixels (e.g. \code{100}).
#' @param displayMode Display mode of annotation track (e.g. \code{COLLAPSED}, \code{EXPANDED}, \code{SQUISHED}).
#' @param SameNameDel A logical value controls whether to add a track if its name has been existed (e.g. \code{FALSE} or \code{TRUE}).
#'
#' @export
Trackoffline <- function(session, version, name, trackHeight=200, displayMode="EXPANDED", SameNameDel=TRUE) {
	print("load gene track in offline mode")
	if ( SameNameDel ) { #// if userdefined track has been loaded, please remove the previous one
		DuplicateTrackRemove(session, name)
	}
	configure <- unique(c(configure, name))

	connect.igvjs <- list(version=version, name=name, trackHeight=trackHeight, displayMode=displayMode) # 
	session$sendCustomMessage("Trackoffline", connect.igvjs)
}

#' Genomic coordinate retrieve
#'
#' @description Obtain genomic coordinate of a zoom region
#'
#' @param session An environment object that is used to access information and functionality by shiny.
#' @param inputid A string represents the name of DOM element (inputId) defined in ui.R (e.g. \code{textInput(inputId, label, value = "")}).
#'
#' @export
Coordinate <- function(session, inputid) {
	# check contents of inputid
	if (! is(inputid, "character") ) { stop("inputid not character!"); }

	connect.igvjs <- list(inputid=inputid)
	session$sendCustomMessage("Coordinate", connect.igvjs)
}

