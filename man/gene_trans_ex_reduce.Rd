% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Breakpoint_domain.R
\name{gene_trans_ex_reduce}
\alias{gene_trans_ex_reduce}
\title{Tag breakpoint position in each transcript (concatenated exons by reducing introns) of partner genes}
\usage{
gene_trans_ex_reduce(breakpoint, object, whole_tx, AB_order, offset = 5)
}
\arguments{
\item{breakpoint}{A numeric vector (only one element) represents genomic coordinate of geneA / geneB breakpoint.}

\item{object}{A list collects geneA / geneB breakpoint annotation, two key elements (i.e. a GRange object: \code{object\$txTr_f}, a data.frame object: \code{object\$dataset}).}

\item{whole_tx}{A GRange object - exon interval grouped by transcript per gene, which was generated by processing \code{txdb_ref} object.}

\item{AB_order}{A string value - either 'upstream' or 'downstream'.}

\item{offset}{A numeric value - set a offset value to adjust breakpoint deviation (default: 5).}
}
\value{
A list has three data.frame elements (i.e. \code{transcript}, \code{select_region} and \code{breakpoint_collect}) with transcript and breakpoint annotations of partner genes
}
\description{
Tag breakpoint position in each transcript (e.g. at exon boundary/within exon or within intron, exon rank) of partner genes
for two-way plotting in domain tab-panel
}
