% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Word_process.R
\name{wordcloud_processs}
\alias{wordcloud_processs}
\title{Create a data.frame for wordcloud visualization}
\usage{
wordcloud_processs(word, gene_freq, type, cancergenes, onco_color, supp_color, rela_color)
}
\arguments{
\item{word}{A data.frame object with two columns (i.e. \code{'name'} - sample; \code{'gene'} - symbol).}

\item{gene_freq}{A numeric value for data filtering using recurrent freq of partner genes.}

\item{type}{A character string (e.g. 'RNA', 'DNA' or 'Mut').}

\item{cancergenes}{A list of cancer genes (e.g. \code{'oncogene'}, \code{'tumorsuppress'} or \code{'related'}).}

\item{onco_color}{A character string represents color code for oncogenes (default: '#ff8566').}

\item{supp_color}{A character string represents color code for tumor suppressed genes (default: '#00ccff').}

\item{rela_color}{A character string represents color code for cancer-related genes (default: '#ffcc33').}
}
\value{
A list with two elements (i.e. \code{'freq'} - a data.frame with two columns; \code{'colorlist'} - a vector)
}
\description{
Create a data.frame for wordcloud visualization
}
