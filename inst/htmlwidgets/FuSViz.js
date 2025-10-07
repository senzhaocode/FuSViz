HTMLWidgets.widget({

	name: 'FuSViz',
	type: 'output',

	factory: function(el, width, height) {
		/*
		 * igv.min.js v2.10.5 (modified version)
		 * (the raw code were modified to meet some requirements of SV visualization, then compiled by Nodejs)
		 * igv.css v2.10.5
		 */

		var widgets = null; // initilize widgets

		return {
			renderValue: function(x) {
				console.log("el:" + el + '!'); // contents of el object
				var div = el;
				var objectid = el.id;
				console.log("el.id:" + el.id + '!'); // contents of el object
				var options = selectIGVoptions(x.genomeName, x.initialLocus, x.displayMode, x.trackHeight);
				
				igv.removeAllBrowsers();
				if ( x.genomeName == 'hg38' || x.genomeName == 'hg19' || x.genomeName == 'GRCm39') {
					igv.createBrowser(div, options).
					then(function (browser) {
						igv.browser = browser;
						// assign igv browser to html container id
						document.getElementById(objectid).igvFuSViz = browser;

						// activate trackclick function 
						igv.browser.on('trackclick', function (track, popoverData) {
							console.log("popoverData contents: "  + popoverData);
							var markup = "<table class=\"igv-popover-table\">";
							// Do not show a pop-over when there's no data (click region without object)
							if (!popoverData || !popoverData.length) { return false; }
							popoverData.forEach(function (element) {
								if (element.name) {
									var value = element.name.toLowerCase() === 'name' ? '<a href="https://www.genecards.org/Search/Keyword?queryString=' + element.value + '" target="_blank">' + element.value + '</a>' : element.value;
									markup += "<tr><td class=\"igv-popover-td\">"
										+ "<div class=\"igv-popover-name-value\">"
										+ "<span class=\"igv-popover-name\">" + element.name + ": " + "</span>"
										+ "<span class=\"igv-popover-value\">" + value + "</span>"
										+ "</div>" + "</td></tr>";
								} else {
									markup += "<tr><td>" + element.toString() + "</td></tr>";
								}
							})
							markup += "</table>";
							// By returning a string from the trackclick handler we're asking IGV to use our custom HTML in its pop-over.
							return markup;
						})
						console.log("Create IGV browser embedded in Shiny app!");
					})
				}
			},

			resize: function(width, height) {
				// The only reason not to implement this method is if your widget naturally scales
			}
		};
	}
});

//----------------------------------------------
// function for loading track in bed format
//----------------------------------------------
Shiny.addCustomMessageHandler("TrackinBed",
	function(message){
		console.log("TrackinBed in js");
		var name = message.name;
		var color = message.color;
		var trackHeight = message.trackHeight;
		var displayMode = message.displayMode;
		var datastr = message.datastr;

		var setting = { type: "annotation", format: "bed", features: datastr, indexed: false, name: name, visibilityWindow: -1,
						order: Number.MAX_VALUE, displayMode: displayMode, color: color, height: trackHeight, removable: true };
		igv.browser.loadTrack(setting);
	}
);

//----------------------------------------------
// function for loading track in bedgrapheformat
//----------------------------------------------
Shiny.addCustomMessageHandler("TrackinBedGraph",
	function(message){
		console.log("TrackinBedGraph in js");
		var name = message.name;
		var color = message.color;
		var trackHeight = message.trackHeight;
		var autoscale = message.autoscale; // R TRUE/FALSE || js true/false
		var displayMode = message.displayMode;
		var min = message.min;
		var max = message.max;
		var datastr = message.datastr;

		var setting = { type: "wig", format: "bedgraph", features: datastr, indexed: false, name: name, autoscale: autoscale, visibilityWindow: -1,
						order: Number.MAX_VALUE, displayMode: displayMode, color: color, height: trackHeight, min: min, max: max, removable: true };
		igv.browser.loadTrack(setting);
	}
);

//-----------------------------------------
// function for load track in bedpe format
//-----------------------------------------
Shiny.addCustomMessageHandler("TrackinBedPe",
	function(message){
		console.log("TrackinBedPe in js");
		var name = message.name;
		var color = message.color;
		var trackHeight = message.trackHeight;
		var logScale = message.logScale; // R TRUE/FALSE || js true/false
		var min = message.min;
		var displayMode = message.displayMode;
		var alpha = message.alpha;
		var thickness = message.thickness;
		var datastr = message.datastr;
		console.log(datastr)

		var setting = { type: "interact", format: "bedpe", features: datastr, indexed: false, name: name, logScale: logScale, order: Number.MAX_VALUE, visibilityWindow: -1,
						color: color, height: trackHeight, showBlocks: true, arcType: "nested", alpha: alpha, displayMode: displayMode, thickness: thickness, removable: true };
		igv.browser.loadTrack(setting);
	}
);

//-----------------------------------------
// function for load track in seg format
//-----------------------------------------
Shiny.addCustomMessageHandler("TrackinSeg",
	function(message){
		console.log("TrackinSeg in js");
		var name = message.name;
		var trackHeight = message.trackHeight;
		var displayMode = message.displayMode;
		var isLog = message.isLog;
		var datastr = message.datastr;

		var setting = {type: "seg", format: "seg", name: name, isLog: isLog, features: datastr, indexed: false, 
					displayMode: displayMode, height: trackHeight, order: Number.MAX_VALUE, removable: true};
		igv.browser.loadTrack(setting);
	}
);

//-------------------------------------------------------
// function for load BAF plot track in vcf format
//-------------------------------------------------------
Shiny.addCustomMessageHandler("FileBAF",
	function(message){
		console.log("FileBAF in js");
		var filedata = message.filedata;
		var BinSize = message.BinSize;
		var trackHeight = message.trackHeight;
		var displayMode = message.displayMode;

		var baf = {}; // link *.vcf.gz per sample
		for (let one of filedata) {
			if ( one.name.endsWith(".baf.vcf.gz") ) {
				if ( one.name in baf ) {
					baf[one.name]["file"] = one;
				} else {
					baf[one.name] = {};
					baf[one.name]["file"] = one;
				}
			} else {
				alert("File format is not accetped: " + one.name);
			}
		}

		var filebaf = []; // set splice junction track configuration
		for (let one in baf) {
			var file = baf[one]["file"];
			if ( file !== undefined ) {
				var host_path = window.location.href + "tmp";
				var local_path = message.local;
				file.datapath=file.datapath.replace(local_path, host_path);
				filebaf.push({ name: one, type: "cnvpytor", 
					url: file.datapath, 
					bin_size: BinSize,
					height: trackHeight,
					visibilityWindow: -1
				});
			} else {
				alert("BAF VCF file does not meet requirement: " + one);
			}
		}

		if ( filebaf.length > 0 ) {
			console.log(filebaf);
			igv.browser.loadTrackList(filebaf);
		}
	}
);

//-------------------------------------------------------
// function for load splice junction track in bed format
//-------------------------------------------------------
Shiny.addCustomMessageHandler("FileinJunction",
	function(message){
		console.log("FileinJunction in js");
		var filedata = message.filedata;
		var unique = message.unique;
		var total = message.total;
		var percet = message.percet;
		var overhang = message.overhang;
		var Selectcolor = message.Selectcolor;
		var Selectthick = message.Selectthick;
		var Selectcurve = message.Selectcurve;
		var trackHeight = message.trackHeight;
		var displayMode = message.displayMode;

		var junction = {}; // link *.SJ.out.bed.gz and *.SJ.out.bed.gz.tbi together per sample
		for (let one of filedata) {
			if ( one.name.endsWith(".SJ.out.bed.gz") ) {
				if ( one.name in junction ) {
					junction[one.name]["file"] = one;
				} else {
					junction[one.name] = {};
					junction[one.name]["file"] = one;
				}
			} else if ( one.name.endsWith(".SJ.out.bed.gz.tbi") ) {
				var tmp_name = one.name.replace(/\.tbi$/, "");
				 if ( tmp_name in junction ) {
					junction[tmp_name]["index"] = one;
				} else {
					junction[tmp_name] = {};
					junction[tmp_name]["index"] = one;
				}
			} else {
				alert("File format is not accetped: " + one.name);
			}				
		}

		var filejunction = []; // set splice junction track configuration
		for (let one in junction) {
			var file = junction[one]["file"];
			var index = junction[one]["index"];
			if ( file !== undefined && index !== undefined ) {
				var host_path = window.location.href + "tmp";
				var local_path = message.local;
				file.datapath=file.datapath.replace(local_path, host_path);
				index.datapath=index.datapath.replace(local_path, host_path);
				filejunction.push({ name: one, type: "junction", format: "bed", 
					url: file.datapath, 
					indexURL: index.datapath,
					minUniquelyMappedReads: unique, 
					minTotalReads: total, 
					maxFractionMultiMappedReads: percet, 
					minSplicedAlignmentOverhang: overhang, 
					thicknessBasedOn: Selectthick, 
					bounceHeightBasedOn: Selectcurve,
					colorBy: Selectcolor,
					labelUniqueReadCount: true,
					labelMultiMappedReadCount: false,
					labelTotalReadCount: false,
					labelMotif: false,
					labelIsAnnotatedJunction: " [A]",
					hideAnnotatedJunctions: false,
					hideUnannotatedJunctions: false,
					hideMotifs: ['GT/AT', 'non-canonical'],
					height: trackHeight,
					visibilityWindow: -1
				});
			} else {
				alert("Junction bed and index files do not match: " + one);
			}
		}

		if ( filejunction.length > 0 ) {
			console.log(filejunction);
			igv.browser.loadTrackList(filejunction);
		}
	}
);

//------------------------------------------------
// function for load gene track in offline mode 
//------------------------------------------------
Shiny.addCustomMessageHandler("Trackoffline",
	function(message){
		console.log("Trackoffline in js");
		var version = message.version;
		var name = message.name;
		var trackHeight = message.trackHeight;
		var displayMode = message.displayMode;
		var setting = undefined;
		
		if ( version === 'hg19' ) {
			setting = {format: "refgene", name: name, url: window.location.href + "Reference/ncbiRefSeq_hg19.txt", indexed: false,
					visibilityWindow: -1, height: trackHeight, searchable: true, displayMode: displayMode};
		} else if ( version === 'hg38' ) {
			setting = {format: "refgene", name: name, url: window.location.href + "Reference/ncbiRefSeq_hg38.txt", indexed: false,
					visibilityWindow: -1, height: trackHeight, searchable: true, displayMode: displayMode};
		} else if ( version === 'GRCm39' ) {
			setting = {format: "refgene", name: name, url: window.location.href + "Reference/ncbiRefSeq_mm39.txt", indexed: false,
					visibilityWindow: -1, height: trackHeight, searchable: true, displayMode: displayMode};
		}
		igv.browser.loadTrack(setting);
	}
);

//-------------------------------------------------------
// function for update BAF plot track in vcf format
//-------------------------------------------------------
Shiny.addCustomMessageHandler("UpdateBAF",
	function(message){
		console.log("UpdateBAF in js");
		var BinSize = message.BinSize;
		var trackHeight = message.trackHeight;
		var displayMode = message.displayMode;

		var i = 0; // count loop step
		var baf = [];
		while ( i <= (igv.browser.trackViews.length - 1) ) {
			var trackName = igv.browser.trackViews[i].track.name;
			if ( trackName ) {
				if ( trackName.match(/.baf.vcf.gz$/) ) { // match to splice junction track
					baf.push({name: trackName, url: igv.browser.trackViews[i].track.config.url});
					igv.browser.removeTrack(igv.browser.trackViews[i].track);
					console.log("Remove replicated track: "+ trackName);
				} else {
					i++;
				}
			} else {
				i++;
			}
		}

		var filebaf = []; // set splice junction track configuration
		if ( baf ) {
			for (let one of baf) {
				filebaf.push({ name: one.name, type: "cnvpytor", 
					url: one.url,
					bin_size: BinSize,
					height: trackHeight,
					visibilityWindow: -1
				});
			}
		} else {
			alert("No BAF plot tracks are available!");
		}

		if ( filebaf.length > 0 ) {
			console.log(filebaf);
			igv.browser.loadTrackList(filebaf);
		}
	}
);
			
//-------------------------------------------------------
// function for update splice junction track in bed format
//-------------------------------------------------------
Shiny.addCustomMessageHandler("UpdateJunction",
	function(message){
		console.log("UpdateJunction in js");
		var unique = message.unique;
		var total = message.total;
		var percet = message.percet;
		var overhang = message.overhang;
		var Selectcolor = message.Selectcolor;
		var Selectthick = message.Selectthick;
		var Selectcurve = message.Selectcurve;
		var trackHeight = message.trackHeight;
		var displayMode = message.displayMode;

		var i = 0; // count loop step
		var junction = [];
		while ( i <= (igv.browser.trackViews.length - 1) ) {
			var trackName = igv.browser.trackViews[i].track.name;
			if ( trackName ) {
				if ( trackName.match(/SJ.out.bed.gz$/) ) { // match to splice junction track
					junction.push({name: trackName, url: igv.browser.trackViews[i].track.config.url, indexURL: igv.browser.trackViews[i].track.config.indexURL});
					igv.browser.removeTrack(igv.browser.trackViews[i].track);
					console.log("Remove replicated track: "+ trackName);
				} else {
					i++;
				}
			} else {
				i++;
			}
		}

		var filejunction = []; // set splice junction track configuration
		if ( junction ) {
			for (let one of junction) {
				filejunction.push({ name: one.name, type: "junction", format: "bed", 
					url: one.url, 
					indexURL: one.indexURL,
					minUniquelyMappedReads: unique, 
					minTotalReads: total, 
					maxFractionMultiMappedReads: percet, 
					minSplicedAlignmentOverhang: overhang, 
					thicknessBasedOn: Selectthick, 
					bounceHeightBasedOn: Selectcurve,
					colorBy: Selectcolor,
					labelUniqueReadCount: true,
					labelMultiMappedReadCount: false,
					labelTotalReadCount: false,
					labelMotif: false,
					labelIsAnnotatedJunction: " [A]",
					hideAnnotatedJunctions: false,
					hideUnannotatedJunctions: false,
					hideMotifs: ['GT/AT', 'non-canonical'],
					height: trackHeight,
					visibilityWindow: -1
				});
			}
		} else {
			alert("No splice junction tracks are available!");
		}

		if ( filejunction.length > 0 ) {
			console.log(filejunction);
			igv.browser.loadTrackList(filejunction);
		}
	}
);

//------------------------------------------------------------------------
// remove tracks with duplicate names, and only keep most upgraded track
//------------------------------------------------------------------------
Shiny.addCustomMessageHandler("DuplicateTrackRemove",
	function(message){
		console.log("DuplicateTrackRemove in js")
		var addtracks = message.addtracks;
		if (! Array.isArray(addtracks) ) { addtracks = [addtracks]; }

		var n = 0; // count loop step
		while ( n <= (igv.browser.trackViews.length - 1) ) {
			var trackName = igv.browser.trackViews[n].track.name;
			if ( trackName ) { // if not undefined
				var replicate = addtracks.indexOf(trackName) >= 0;
				if ( replicate ) { // matched
					igv.browser.removeTrack(igv.browser.trackViews[n].track);
					console.log("Remove replicated track: "+ trackName);
				}
			}
			n++;
		}
	}
);

//--------------------------------
// Genomic coordinate retrieve
//--------------------------------
Shiny.addCustomMessageHandler("Coordinate",
	function(message){
		console.log("Retrieve genomic coordinate in js")

		var inputid = message.inputid;
		var position = igv.browser.currentLoci();

		if ( typeof(position) == "string" ) {
			var posnew = position.replace(/\.[\d]+/g, "");
			Shiny.setInputValue(inputid, posnew, {priority: "event"});
			console.log("current position: " + position);
		} else {
			Shiny.setInputValue(inputid, "", {priority: "event"});
		}
	}
);

//------------------------------------
// function for load track of bam url
//------------------------------------
Shiny.addCustomMessageHandler("TrackinBAM",
	function(message){
		var bam = message.bam;
		var bamindex = message.bamindex;
		var trackHeight = message.trackHeight;
		var bam_name = null;
		var bam_type = null;

		if ( bam.endsWith(".bam") && (bamindex.endsWith(".bam.bai") || bamindex.endsWith(".bam.csi")) ) { // BAM format
			bam_name = bam.split(/\/|\\/).pop();
			bam_type = "bam";
		} else if ( bam.endsWith(".cram") && bamindex.endsWith(".cram.crai")  ) { // CRAM format
			bam_name = bam.split(/\/|\\/).pop();
			bam_type = "cram";
		} else {
			alert("BAM / CRAM format not correct: " + bam);
		}

		var setting = {type: "alignment", format: bam_type, url: bam, indexURL: bamindex, name: bam_name, order: Number.MAX_VALUE, height: trackHeight};
		igv.browser.loadTrack(setting);
	}
);

//----------------------------------------------------
//  function for select genome version (hg19 or hg38)
//---------------------------------------------------- 
function selectIGVoptions(genomeName, initialLocus, displayMode, trackHeight) {
	// setting for hg19 version
	var genome_hg19 = {
		minimumBases: 5,
		flanking: 1000,
		showRuler: true,
		reference: {
			id: "hg19",
			name: "Human (CRCh37/hg19)",
			fastaURL: "https://s3.amazonaws.com/igv.broadinstitute.org/genomes/seq/hg19/hg19.fasta",
			indexURL: "https://s3.amazonaws.com/igv.broadinstitute.org/genomes/seq/hg19/hg19.fasta.fai",
			cytobandURL: "https://s3.amazonaws.com/igv.broadinstitute.org/genomes/seq/hg19/cytoBand.txt"
		},
		tracks: [
			{
				name: 'RefSeq Genes (hg19)',
				format: "refgene",
				url: "https://s3.amazonaws.com/igv.org.genomes/hg19/ncbiRefSeq.sorted.txt.gz",
				indexURL: "https://s3.amazonaws.com/igv.org.genomes/hg19/ncbiRefSeq.sorted.txt.gz.tbi",
				visibilityWindow: -1,
				removable: true,
				height: trackHeight,
				displayMode: displayMode
			}
		]
	}
	// setting for hg38 version
	var genome_hg38 = {
		minimumBases: 5,
		flanking: 1000,
		showRuler: true,
		reference: {
			id: "hg38",
			name: "Human (GRCh38/hg38)",
			fastaURL: "https://s3.amazonaws.com/igv.broadinstitute.org/genomes/seq/hg38/hg38.fa",
			indexURL: "https://s3.amazonaws.com/igv.broadinstitute.org/genomes/seq/hg38/hg38.fa.fai",
			cytobandURL: "https://s3.amazonaws.com/igv.org.genomes/hg38/annotations/cytoBandIdeo.txt.gz"
		},
		tracks: [
			{
				name: 'RefSeq Genes (hg38)',
				format: "refgene",
				url: "https://s3.amazonaws.com/igv.org.genomes/hg38/ncbiRefSeq.sorted.txt.gz",
				indexURL: "https://s3.amazonaws.com/igv.org.genomes/hg38/ncbiRefSeq.sorted.txt.gz.tbi",
				visibilityWindow: -1,
				removable: true,
				height: trackHeight,
				displayMode: displayMode
			}
		]
	}
	// setting for mm39 version
	var genome_mm39 = {
		minimumBases: 5,
		flanking: 1000,
		showRuler: true,
		reference: {
			id: "mm39",
			name: "Mouse (GRChm39 / mm39)",
			fastaURL: "https://s3.amazonaws.com/igv.org.genomes/mm39/mm39.fa",
			indexURL: "https://s3.amazonaws.com/igv.org.genomes/mm39/mm39.fa.fai",
			cytobandURL: "https://s3.amazonaws.com/igv.org.genomes/mm39/cytoBandIdeo.txt.gz"
		},
		tracks: [
			{
				name: 'RefSeq Genes (mm39)',
				format: "refgene",
				url: "https://s3.amazonaws.com/igv.org.genomes/mm39/ncbiRefSeq.txt.gz",
				indexURL: "https://s3.amazonaws.com/igv.org.genomes/mm39/ncbiRefSeq.txt.gz.tbi",
				visibilityWindow: -1,
				removable: true,
				height: trackHeight,
				displayMode: displayMode
			}
		]
	}

	if ( genomeName === "hg19" ) {
		return(genome_hg19);
	} else if ( genomeName === "hg38" ) {
		return(genome_hg38);
	} else if ( genomeName === "GRCm39" ) {
		return(genome_mm39);
	} else {
		console.log("Genome version " + genomeName + " is not present!");
		return(undefined);
	}
}

