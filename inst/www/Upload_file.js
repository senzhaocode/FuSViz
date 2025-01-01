// All alignments and annotation files are loaded directly from a user’s local file system, or via URL from web servers and cloud providers. 
// They are run by the web browser and no data is ever uploaded to the host site.
function myFunction(element) {
	$('a[href="#shiny-tab-igv"]').on('shown.bs.tab', function() {
		$('.sidebar-menu a').filter('a[href="#' + 'shiny-tab-igv' + '"]').tab('show');
		if ( igv.browser != undefined && element != undefined ) {
			igv.browser.search(element);
			element = undefined;
		}
	});
}

function load() {
	var fileobj = document.getElementById("uploadfile");
	var elements = fileobj.files;

	var filename = "<ul>";
	for (let one of elements) {
		filename += "<li>" + one.name + "</li>";
	}
	filename += "</ul>";
	document.getElementById("filename").innerHTML = filename;

// e.g. var elements = [{"name":"senz.bam"}, {"name":"senz.vcf.gz"}, {"name":"senz.bed.gz"}, {"name":"senz.gtf.gz"}, 
//		{"name":"senz.bam.bai"}, {"name":"senz.vcf.gz.tbi"}, {"name":"senz.bed.gz.tbi"}, {"name":"senz.gtf.gz.tbi"}];

	// assign upload files to different formats
	var bam = {}; // alignment in BAM format
	var cram = {}; // alignment in CRAM format
	var vcf = {}; // annotation in VCF format
	var bed = {}; // annotation in BED format
	var gtf = {}; // annotation in GTF format
	var localref = {}; // reference sequence in FASTA format

	for (let one of elements) {
		if ( one.name.endsWith(".fasta") ) {
			localref["file"] = one;
		} else if ( one.name.endsWith(".fai") ) {
			localref["index"] = one;
		} else if ( one.name.startsWith("cytoBand") ) {
			localref["cytoband"] = one;
		} else if ( one.name.endsWith(".cram") ) {
			if ( one.name in cram ) {
				cram[one.name]["file"] = one;
			} else {
				cram[one.name] = {};
				cram[one.name]["file"] = one;
			}
		} else if ( one.name.endsWith(".crai") ) {
			var tmp_name = one.name.replace(/\.crai$/, "");
			if ( tmp_name in cram ) {
				cram[tmp_name]["index"] = one;
			} else {
				cram[tmp_name] = {};
				cram[tmp_name]["index"] = one;
			}
		} else if ( one.name.endsWith(".bam") ) {
			if ( one.name in bam ) {
				bam[one.name]["file"] = one;
			} else {
				bam[one.name] = {};
				bam[one.name]["file"] = one;
			}			
		} else if ( one.name.endsWith(".bai") ) {
			var tmp_name = one.name.replace(/\.bai$/, "");
			if ( tmp_name in bam ) {
				bam[tmp_name]["index"] = one;
			} else {
				bam[tmp_name] = {};
				bam[tmp_name]["index"] = one;
			}
		} else if ( one.name.endsWith(".vcf.gz") ) {
			if ( one.name in vcf ) {
				vcf[one.name]["file"] = one;
			} else {
				vcf[one.name] = {};
				vcf[one.name]["file"] = one;
			}
		} else if ( one.name.endsWith(".vcf.gz.tbi") ) {
			var tmp_name = one.name.replace(/\.tbi$/, "");
			if ( tmp_name in vcf ) {
				vcf[tmp_name]["index"] = one;
			} else {
				vcf[tmp_name] = {};
				vcf[tmp_name]["index"] = one;
			}
		} else if ( one.name.endsWith(".bed.gz") ) {
			if ( one.name in bed ) {
				bed[one.name]["file"] = one;
			} else {
				bed[one.name] = {};
				bed[one.name]["file"] = one;
			}
		} else if ( one.name.endsWith(".bed.gz.tbi") ) {
			var tmp_name = one.name.replace(/\.tbi$/, "");
			if ( tmp_name in bed ) {
				bed[tmp_name]["index"] = one;
			} else {
				bed[tmp_name] = {};
				bed[tmp_name]["index"] = one;
			}
		} else if ( one.name.endsWith(".gtf.gz") ) {
			if ( one.name in gtf ) {
				gtf[one.name]["file"] = one;
			} else {
				gtf[one.name] = {};
				gtf[one.name]["file"] = one;
			}
		} else if ( one.name.endsWith(".gtf.gz.tbi") ) {
			var tmp_name = one.name.replace(/\.tbi$/, "");
			if ( tmp_name in gtf ) {
				gtf[tmp_name]["index"] = one;
			} else {
				gtf[tmp_name] = {};
				gtf[tmp_name]["index"] = one;
			}
		} else {
			alert("File type not accetped: " + one.name);
		}
	}

	var options = undefined;
	// set reference genome sequence configuration
	if ( Object.keys(localref).length !== 0 ) {
		if ( 'file' in localref && 'index' in localref && 'cytoband' in localref ) {
			if ( localref['file'].name.startsWith('hg19') ) {
				options = {
					loadDefaultGenomes: false,
					flanking: 1000,
					minimumBases: 5,
					showRuler: true,
					reference: { id: "hg19", name: "local", fastaURL: localref['file'], indexURL: localref['index'], cytobandURL: localref['cytoband'] }}
			} else if ( localref['file'].name.startsWith('hg38') ) {
				options = { 
					loadDefaultGenomes: false,
					flanking: 1000,
					minimumBases: 5,
					showRuler: true,
					reference: { id: "hg38", name: "local", fastaURL: localref['file'], indexURL: localref['index'], cytobandURL: localref['cytoband'] }}
			} else if ( localref['file'].name.startsWith('mm39') ) { 
				options = { 
					loadDefaultGenomes: false,
					flanking: 1000,
					minimumBases: 5,
					showRuler: true,
					reference: { id: "mm39", name: "local", fastaURL: localref['file'], indexURL: localref['index'], cytobandURL: localref['cytoband'] }}
			} else {
				alert("Genome reference version not available!!!");
			}
		} else {
			alert("FASTA, fai index and cytoband files is absent!!!");
		}
	}
	var uploadtrack = [];
	// set BAM format configuration
	for (let one in bam) {
		var file = bam[one]["file"];
		var index = bam[one]["index"];

		if ( file !== undefined && index !== undefined ) {
			uploadtrack.push({ name: one, type: "alignment", format: "bam", url: file, indexURL: index, alleleFreqThreshold: 0.05, samplingDepth: 2000, alignmentRowHeight: 10, height: 800, visibilityWindow: 10000, order: Number.MAX_VALUE });
		} else {
			alert("BAM and index files not match: " + one);
		}
	}
	// set CRAM format configuration
	for (let one in cram) {
		var file = cram[one]["file"];
		var index = cram[one]["index"];

		if ( file !== undefined && index !== undefined ) {
			uploadtrack.push({ name: one, type: "alignment", format: "cram", url: file, indexURL: index, alleleFreqThreshold: 0.05, samplingDepth: 2000, alignmentRowHeight: 10, height: 800, visibilityWindow: 10000, order: Number.MAX_VALUE });
		} else {
			alert("CRAM and index files not match: " + one);
		}
	}
	// set VCF format configuration
	for (let one in vcf) {
		var file = vcf[one]["file"];
		var index = vcf[one]["index"];

		if ( file !== undefined && index !== undefined ) {
			uploadtrack.push({ name: one, type: "variant", format: "vcf", url: file, indexURL: index, order: Number.MAX_VALUE });
		} else {
			alert("VCF and index files not match: " + one);
		}
	}
	// set BED format configuration
	for (let one in bed) {
		var file = bed[one]["file"];
		var index = bed[one]["index"];

		if ( file !== undefined && index !== undefined ) {
			uploadtrack.push({ name: one, type: "annotation", format: "bed", url: file, indexURL: index, order: Number.MAX_VALUE });
		} else {
			alert("BED and index files not match: " + one);
		}
	}
	// set GTF format configuration
	for (let one in gtf) {
		var file = gtf[one]["file"];
		var index = gtf[one]["index"];

		if ( file !== undefined && index !== undefined ) {
			uploadtrack.push({ name: one, type: "annotation", format: "gtf", url: file, indexURL: index, order: Number.MAX_VALUE });
		} else {
			alert("GTF and index files not match: " + one);
		}
	}

	// Load tracks to igv widget
	if ( options !== undefined ) {
		const div = document.getElementById("FuSViz");
		igv.removeAllBrowsers()
		igv.createBrowser(div, options).
			then(function (browser) {
				igv.browser = browser;
				document.getElementById("FuSViz").igvFuSViz = browser;

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
			})
	} else {
		var igvFuSViz = document.getElementById('FuSViz').igvFuSViz;	
		if ( uploadtrack.length > 0 ) {
			igvFuSViz.loadTrackList(uploadtrack);
		}
	}
}
