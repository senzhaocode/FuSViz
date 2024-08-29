#' Create input format for visNetwork package 
#'
#' @description Create input format of DNA and RNA SVs for visNetwork visualization
#'
#' @param tmp A data.frame with three columns (i.e. \code{'gene1'}, \code{'gene2'}, \code{'name'}).
#' @param type A string (e.g. \code{'DNA'} or \code{'RNA'}).
#' @param cancergenes A list of cancer genes (e.g. \code{'oncogene'}, \code{'tumorsuppress'} or \code{'related'}).
#' @param color_onco A string of color value (e.g. \code{'#ff0000'} or \code{'rgb(100,0,100)'}); if partner gene is oncogene (defalut: '#ff8566').
#' @param color_supp A string of color value (e.g. \code{'#ff0000'} or \code{'rgb(100,0,100)'}); if partner gene is tumor suppress gene (defalut: '#00ccff').
#' @param color_rela A string of color value (e.g. \code{'#ff0000'} or \code{'rgb(100,0,100)'}); if partner gene is cancer-related gene (defalut: '#ffcc33').
#' @param color_inter A string of color value (e.g. \code{'#ff0000'} or \code{'rgb(100,0,100)'}); if breakpoint is in intergenic region (defalut: '#262626').
#' @param color_other A string of color value (e.g. \code{'#ff0000'} or \code{'rgb(100,0,100)'}); if not belonging to any category above (defalut: '#f2f2f2').
#'
#' @return A list with three data.frame elements (i.e. \code{'nodes'}, \code{'edges'} and \code{'degree_score'}).
#'
#' @export
network_process <- function(tmp, type, cancergenes, color_onco, color_supp, color_rela, color_inter, color_other) {
	tmp = tmp[! (tmp$gene1 == '*' & tmp$gene2 == '*'), ]; #// remove the entry with both breakpoints at integenic region
	tmp = tmp[! (grepl("ENSG00", tmp$gene1) & grepl("ENSG00", tmp$gene2)), ]; #// remove the entry with both partners as ensembl_id
	geneset = c(names(cancergenes$oncogene), names(cancergenes$tumorsuppress), names(cancergenes$related)); #// merge pre-defined geneset

	#// only keep gene pairs if either gene1 or gene2 (at least one of them) is in pre-defined geneset
	# tmp = tmp[! (tmp$gene1 == '*' & !(tmp$gene2 %in% geneset)), ]; #// if gene1 == *; gene2 shall be in pre-defined geneset
	# tmp = tmp[! (tmp$gene2 == '*' & !(tmp$gene1 %in% geneset)), ]; #// if gene2 == *; gene1 shall be in pre-defined geneset
	tmp = tmp[! (!(tmp$gene1 %in% geneset) & !(tmp$gene2 %in% geneset)), ];

	#// unique pair of gene1-gene2 / gene2-gene1 (ignoring gene1 and gene2 order) -- only for DNA SV
	if ( type == "DNA") {
		tmp_gene = tmp[, c("gene1", "gene2")];
		#// create unique pair index
		two_no_order <- function(x) unique(cbind(pmin(x[,1], x[,2]), pmax(x[,1], x[,2])))
		tmp_gene_unq = two_no_order(tmp_gene);
		tmp_gene_unq = unique(as.data.frame(tmp_gene_unq)); 
		names(tmp_gene_unq) = c("gene1", "gene2");

		#// TRUE if gene1-gene2 matches unique pair index; otherwise as FALSE (will swap the order of gene1 and gene2)
		tag_bilateral = apply(tmp, 1, function(x) {
			tmp_gene1 = as.character(x[1]);
			tmp_gene2 = as.character(x[2]);
			forward = nrow(tmp_gene_unq[tmp_gene_unq$gene1==tmp_gene1 & tmp_gene_unq$gene2==tmp_gene2, ]);
			if ( forward > 0 ) { return(TRUE); } else { return(FALSE); }
		})
		tmp = cbind(tmp, tag_bilateral); 
		tmp[tmp$tag_bilateral == FALSE, c("gene1", "gene2")] = tmp[tmp$tag_bilateral == FALSE, c("gene2", "gene1")] 
		tmp$tag_bilateral = NULL;
	}

	#// only keep the unique fusion partner pairs per sample (if multiple breakpoint combination present)
	tmp = unique(tmp); 
	tmp$gene1 = ifelse(tmp$gene1 == '*', yes=paste("**", tmp$gene2, "**", sep=""), no=tmp$gene1); #// if gene1 == '*'; change as '**gene2**'
	tmp$gene2 = ifelse(tmp$gene2 == '*', yes=paste("**", tmp$gene1, "**", sep=""), no=tmp$gene2); #// if gene2 == '*'; change as '**gene1**'

	# initilise nodes, edges and nodes_degree_score
	nodes = NULL;	edges = NULL;
	nodes_degree_score = data.frame(nodes=NULL, degree=NULL, score=NULL, stringsAsFactors=F);
	if ( nrow(tmp) > 0 ) {
		#// build data.frame for nodes
		nodes = c(tmp$gene1, tmp$gene2);    nodes = table(nodes);
		nodes = as.data.frame(nodes, stringsAsFactors = F);
		#// build data.frame for edges
		# edges = as.data.frame(tmp %>% dplyr::group_by(gene1, gene2) %>% dplyr::summarise(value = n())); #// value = num of sample harboring fusions
		edges = aggregate(.~gene1+gene2, tmp, FUN=length);	names(edges)[3] = "value";

		#// tag node with its degree = 1 as FALSE and assign color to partner based on gene category
		tag_node = apply(nodes, 1, function(x) {
			color = NULL; #// tag the color category of a node (e.g. red, blue, orange, black, grey)
			hidden = NULL; #// set as T if color = grey, otherwise set as F
			genename = as.character(x[1]);
			tmp_edge = edges[edges$gene1 == genename | edges$gene2 == genename, ]; #// select edge contain given node
			num_sample = sum(as.numeric(tmp_edge$value));
			judge = as.logical(ifelse(num_sample > 0, yes = TRUE, no = FALSE));

			if ( genename %in% names(cancergenes$oncogene) ) {
				color = color_onco;  hidden = FALSE;
			} else if ( genename %in% names(cancergenes$tumorsuppress) ) {
				color = color_supp; hidden = FALSE;
			} else if ( genename %in% names(cancergenes$related) ) {
				color = color_rela;   hidden = FALSE;
			} else {
				if ( grepl('^\\*', genename) ) {
					color = color_inter; hidden = TRUE; #// one of SV partner in intergenic region
				} else { 
					color = color_other; hidden = TRUE;
				}
			}
			return(c(judge, color, hidden));
		})
		tag_node = t(tag_node); # transposition 
		nodes = cbind(nodes, tag_node);
		names(nodes)[4] = "tag_color";  
		names(nodes)[3] = "tag_node"; 
		nodes$tag_node = as.logical(nodes$tag_node);
		names(nodes)[5] = "tag_hidden";  
		nodes$tag_hidden = as.logical(nodes$tag_hidden);

		#// tag both nodes with degree = 1 as FALSE in edges
		tag_edge = apply(edges, 1, function(x) {
			judge = NULL; #// if both nodes with tag_node 'F', set as F; otherwise set as T
			hidden = NULL; #// if either of two nodes is gray, set as T; otherwise set as T
			judge1 = nodes[nodes$nodes==as.character(x[1]), ]$tag_node;
			judge2 = nodes[nodes$nodes==as.character(x[2]), ]$tag_node;
			hidden1 = nodes[nodes$nodes==as.character(x[1]), ]$tag_hidden;
			hidden2 = nodes[nodes$nodes==as.character(x[2]), ]$tag_hidden;
			if ( judge1 == FALSE && judge2 == FALSE ) { 
				if ( hidden1 == TRUE || hidden2 == TRUE ) {
					judge = FALSE;  hidden = TRUE;
				} else {
					judge = FALSE;  hidden = FALSE;
				}
			} else { 
				if ( hidden1 == TRUE || hidden2 == TRUE ) {
					judge = TRUE;  hidden = TRUE;
				} else {
					judge = TRUE;  hidden = FALSE;
				}
			}
			return(c(judge, hidden));
		})
		tag_edge = t(tag_edge); # convert 
		edges = cbind(edges, tag_edge);
		names(edges)[4] = "tag_edge";   
		edges$tag_edge = as.logical(edges$tag_edge);
		names(edges)[5] = "tag_hidden";   
		edges$tag_hidden = as.logical(edges$tag_hidden);
		edges = edges[edges$tag_edge == T, ];
		nodes = nodes[nodes$nodes %in% unique(c(edges$gene1, edges$gene2)), ];
   
		#// calculate degree and score (recurrence) per each node
		tag_degree_score = apply(nodes, 1, function(x) {
			gene_name = as.character(x[1]);
			select_edge = edges[edges$gene1 == gene_name | edges$gene2 == gene_name, ];
			degree = nrow(select_edge);
			score = sum(select_edge$value);
			return(c(degree, score));
		})
		if ( length(tag_degree_score) > 0 ) {
			tag_degree_score = t(tag_degree_score); # transposition
			nodes_degree_score = data.frame(nodes=as.character(nodes$nodes), degree=as.numeric(tag_degree_score[,1]), score=as.numeric(tag_degree_score[,2]), stringsAsFactors=F);
		}
	} 
	tmp_new <- list(nodes=nodes, edges=edges, degree_score=nodes_degree_score);
	return(tmp_new);
}

