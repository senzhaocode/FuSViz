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
				
				igv.removeAllBrowsers()
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

		var setting = { type: "annotation", format: "bed", features: datastr, indexed: false, name: name, 
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

		var setting = { type: "wig", format: "bedgraph", features: datastr, indexed: false, name: name, autoscale: autoscale,
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

		var setting = { type: "interaction", format: "bedpe", features: datastr, indexed: false, name: name, logScale: logScale, order: Number.MAX_VALUE,
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
			setting = {format: "bed", name: name, url: window.location.href + "Reference/refGene.hg19.bed", indexed: false,
					visibilityWindow: -1, height: trackHeight, searchable: true, displayMode: displayMode};
		} else if ( version === 'hg38' ) {
			setting = {format: "refgene", name: name, url: window.location.href + "Reference/refGene.sorted.txt", indexed: false,
					visibilityWindow: -1, height: trackHeight, searchable: true, displayMode: displayMode};
		} else if ( version === 'GRCm39' ) {
			setting = {format: "refgene", name: name, url: window.location.href + "Reference/ncbiRefSeq.txt", indexed: false,
					visibilityWindow: -1, height: trackHeight, searchable: true, displayMode: displayMode};
		}
		igv.browser.loadTrack(setting);
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

		if ( typeof(position[0]) == "string" ) {
			Shiny.setInputValue(inputid, position[0], {priority: "event"});
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

		if ( bam.endsWith(".bam") && bamindex.endsWith(".bam.bai") ) { // BAM format
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
	// setting for hg38 version
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

