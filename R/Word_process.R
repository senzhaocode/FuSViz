#' Create a data.frame for wordcloud visualization
#'
#' @description Create a data.frame for wordcloud visualization
#'
#' @param word A data.frame object with two columns (i.e. \code{'name'} - sample; \code{'gene'} - symbol).
#' @param gene_freq A numeric value for data filtering using recurrent freq of partner genes.
#' @param type A character string (e.g. 'RNA', 'DNA' or 'Mut').
#' @param cancergenes A list of cancer genes (e.g. \code{'oncogene'}, \code{'tumorsuppress'} or \code{'related'}).
#' @param onco_color A character string represents color code for oncogenes (default: '#ff8566').
#' @param supp_color A character string represents color code for tumor suppressed genes (default: '#00ccff').
#' @param rela_color A character string represents color code for cancer-related genes (default: '#ffcc33').
#'
#' @return A list with two elements (i.e. \code{'freq'} - a data.frame with two columns; \code{'colorlist'} - a vector)
#'
#' @export
wordcloud_processs <- function(word, gene_freq, type, cancergenes, onco_color, supp_color, rela_color) {
	stopifnot(is.character(type));
	stopifnot(is.numeric(gene_freq));

	if ( type == 'RNA' || type == 'DNA' ) {
		# assemble new data.frame with 'name' and 'gene'
		word = data.frame(name=c(as.character(word$name), as.character(word$name)), gene=c(as.character(word$gene1), as.character(word$gene2)), stringsAsFactors = F);
	}
	#// remove intergenic annotation (*) / gene with ensembl_id (ENSGXXXXX)
	word = word[! word$gene == '*', ];
	word = word[! grepl("ENSG00", word$gene), ];
	word = unique(word);	word = word$gene;	word = table(word);
	#// control the frequency
	word = word[word >= gene_freq];
	word_onco = word[names(word) %in% names(cancergenes$oncogene)];
	if ( length(word_onco) == 0 ) { word_onco = NULL; }
	word_supp = word[names(word) %in% names(cancergenes$tumorsuppress)];
	if ( length(word_supp) == 0 ) { word_supp = NULL; }
	word_rela = word[names(word) %in% names(cancergenes$related)];
	if ( length(word_rela) == 0 ) { word_rela = NULL; }
	#// if 'word_onco', 'word_supp' and 'word_rela' are NULL
	if ( is.null(word_onco) && is.null(word_supp) && is.null(word_rela) ) {
		return(list(freq=NULL, colorlist=NULL));
	} else {
		freq = rbind(data.frame(word = names(word_onco), Freq = as.vector(word_onco), stringsAsFactors = F),
				data.frame(word = names(word_supp), Freq = as.vector(word_supp), stringsAsFactors = F),
				data.frame(word = names(word_rela), Freq = as.vector(word_rela), stringsAsFactors = F));
		#// sqrt(raw freq number) -> reduce scale for better visualization
		freq[,2] = round(sqrt(freq[,2]), 3);
		colorlist = c(rep(onco_color, length(word_onco)), rep(supp_color, length(word_supp)), rep(rela_color, length(word_rela)));
		return(list(freq=freq, colorlist=colorlist));
	}
}

