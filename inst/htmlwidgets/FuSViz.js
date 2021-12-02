HTMLWidgets.widget({

	name: 'FuSViz',
	type: 'output',

	factory: function(el, width, height) {
		/*
		 * igv.min.js v2.10.4 (modified version)
		 * (the raw code were modified to meet some requirements of SV visualization, then compiled by Nodejs)
		 * igv.css v2.10.4
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
											+ "<span class=\"igv-popover-name\">" + element.name + "</span>"
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
					console.log("Remove replicated track: "+ igv.browser.trackViews[n].track);
				}
			}
			n++;
		}
	}
)

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
)

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
		} else if ( bam.endsWith(".bam.cram") && bamindex.endsWith(".bam.cram.crai")  ) { // CRAM format
			bam_name = bam.split(/\/|\\/).pop();
			bam_type = "cram";
		} else {
			alert("BAM / CRAM format not correct: " + bam);
		}

		var setting = {type: "alignment", format: bam_type, url: bam, indexURL: bamindex, name: bam_name, order: Number.MAX_VALUE, height: trackHeight};
		igv.browser.loadTrack(setting);
	}
)

//-------------------------------------------------------------------
// function for loading and processing userdefined bed files locally
//-------------------------------------------------------------------
Shiny.addCustomMessageHandler("TrackinFile",
	function(message){
		var fileobj = message.fileobj;
		var trackHeight = message.trackHeight;

		// list upload file names
		var filename = "<ul>";
		for (let one of fileobj) {
			filename += "<li>" + one.name + "</li>";
		}
		filename += "</ul>";
		document.getElementById("filename").innerHTML = filename;

		var vcf = {}; // annotation in VCF format
		var bed = {}; // annotation in BED format
		var gtf = {}; // annotation in GTF format
		for (let one of fileobj) {
			if ( one.name.endsWith(".bed.gz") ) {
				if ( one.name in bed ) {
					bed[one.name]["file"] = one.newpath;
				} else {
					bed[one.name] = {};
					bed[one.name]["file"] = one.newpath;
				}
			} else if ( one.name.endsWith(".bed.gz.tbi") ) {
				var tmp_name = one.name.replace(/\.tbi$/, "");
				if ( tmp_name in bed ) {
					bed[tmp_name]["index"] = one.newpath;
				} else {
					bed[tmp_name] = {};
					bed[tmp_name]["index"] = one.newpath;
				}
			} else if ( one.name.endsWith(".vcf.gz") ) {
				if ( one.name in vcf ) {
					vcf[one.name]["file"] = one.newpath;
				} else {
					vcf[one.name] = {};
					vcf[one.name]["file"] = one.newpath;
				}
			} else if ( one.name.endsWith(".vcf.gz.tbi") ) {
				var tmp_name = one.name.replace(/\.tbi$/, "");
				if ( tmp_name in vcf ) {
					vcf[tmp_name]["index"] = one.newpath;
				} else {
					vcf[tmp_name] = {};
					vcf[tmp_name]["index"] = one.newpath;
				}
			} else if ( one.name.endsWith(".gtf.gz") ) {
				if ( one.name in gtf ) {
					gtf[one.name]["file"] = one.newpath;
				} else {
					gtf[one.name] = {};
					gtf[one.name]["file"] = one.newpath;
				}
			} else if ( one.name.endsWith(".gtf.gz.tbi") ) {
				var tmp_name = one.name.replace(/\.tbi$/, "");
				if ( tmp_name in gtf ) {
					gtf[tmp_name]["index"] = one.newpath;
				} else {
					gtf[tmp_name] = {};
					gtf[tmp_name]["index"] = one.newpath;
				}
			} else {
				alert("File type not accetped: " + one.name);
			}
		}

		var uploadtrack = [];
		// set BED format configuration
		for (let one of bed) {
			var file = bed[one]["file"];
			var index = bed[one]["index"];

			if ( file !== undefined && index !== undefined ) {
				uploadtrack.push({
					type: "annotation",
					format: "bed",
					sourceType: "file",
					height: trackHeight,
					name: one,
					url: window.location.href + file,
					indexURL: window.location.href + index,
					order: Number.MAX_VALUE
				})
			} else {
				alert("BED and index files not match: " + one);
			}
		}
		// set GTF format configuration
		for (let one of gtf) {
			var file = gtf[one]["file"];
			var index = gtf[one]["index"];

			if ( file !== undefined && index !== undefined ) {
				uploadtrack.push({
					type: "annotation",
					format: "gtf",
					sourceType: "file",
					height: trackHeight,
					name: one,
					url: window.location.href + file,
					indexURL: window.location.href + index,
					order: Number.MAX_VALUE
				})
			} else {
				alert("GTF and index files not match: " + one);
			}
		}
		// set VCF format configuration
		for (let one of vcf) {
			var file = vcf[one]["file"];
			var index = vcf[one]["index"];

			if ( file !== undefined && index !== undefined ) {
				uploadtrack.push({
					type: "variant",
					format: "vcf",
					sourceType: "file",
					height: trackHeight,
					name: one,
					url: window.location.href + file,
					indexURL: window.location.href + index,
					order: Number.MAX_VALUE
				})
			} else {
				alert("VCF and index files not match: " + one);
			}
		}

		//console.log("uploadtrack: " + uploadtrack + "!");
		if (uploadtrack.length > 0) {
			igv.browser.loadTrackList(uploadtrack);
		}
	}
)

//----------------------------------------------------
//  function for select genome version (hg19 or hg38)
//---------------------------------------------------- 
function selectIGVoptions(genomeName, initialLocus, displayMode, trackHeight) {
	// setting for hg19 version
	var genome_hg19 = {
		locus: initialLocus,
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
				url: "https://s3.amazonaws.com/igv.broadinstitute.org/annotations/hg19/genes/refGene.hg19.bed.gz",
				indexURL: "https://s3.amazonaws.com/igv.broadinstitute.org/annotations/hg19/genes/refGene.hg19.bed.gz.tbi",
				visibilityWindow: -1,
				removable: false,
				height: trackHeight,
				displayMode: displayMode
			}
		]
	}	
	// setting for hg19 version (offline)
	var genome_hg19_offline = {
		locus: initialLocus,
		minimumBases: 5,
		flanking: 1000,
		showRuler: true,
		reference: {
			id: "hg19",
			name: "Human [CRCh37/hg19]",
			fastaURL: window.location.href + "Reference/hg19.fasta",
			indexURL: window.location.href + "Reference/hg19.fasta.fai",
			cytobandURL: window.location.href + "Reference/cytoBand.txt"
		},
		tracks: [
			{
				name: 'RefSeq Genes [hg19]',
				url: window.location.href + "Reference/refGene.hg19.bed.gz",
				indexURL: window.location.href + "Reference/refGene.hg19.bed.gz.tbi",
				visibilityWindow: -1,
				removable: false,
				height: trackHeight,
				displayMode: displayMode
			}
		]
	}
	// setting for hg38 version
	var genome_hg38 = {
		locus: initialLocus,
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
				url: "https://s3.amazonaws.com/igv.org.genomes/hg38/refGene.sorted.txt.gz",
				indexURL: "https://s3.amazonaws.com/igv.org.genomes/hg38/refGene.sorted.txt.gz.tbi",
				visibilityWindow: -1,
				removable: false,
				height: trackHeight,
				displayMode: displayMode
			}
		]
	}
	// setting for hg38 version (offline)
	var genome_hg38_offline = {
		locus: initialLocus,
		minimumBases: 5,
		flanking: 1000,
		showRuler: true,
		reference: {
			id: "hg38",
			name: "Human [GRCh38/hg38]",
			fastaURL: window.location.href + "Reference/hg38.fa",
			indexURL: window.location.href + "Reference/hg38.fa.fai",
			cytobandURL: window.location.href + "Reference/cytoBandIdeo.txt.gz"
		},
		tracks: [
			{
				name: 'RefSeq Genes [hg38]',
				url: window.location.href + "Reference/refGene.sorted.txt.gz",
				indexURL: window.location.href + "Reference/refGene.sorted.txt.gz.tbi",
				visibilityWindow: -1,
				removable: false,
				height: trackHeight,
				displayMode: displayMode
			}
		]
	}

	if ( genomeName === "hg19" ) {
		return(genome_hg19);
	} else if ( genomeName === "hg38" ) {
		return(genome_hg38);
	} else if ( genomeName === "hg19_offline" ) {
		return(genome_hg19_offline);
	} else if ( genomeName === "hg38_offline" ) {
		return(genome_hg38_offline);
	} else {
		console.log("Genome version " + genomeName + " is not present!");
		return(undefined);
	}
}

