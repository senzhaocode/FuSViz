suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(shiny)))
suppressWarnings(suppressPackageStartupMessages(library(shinydashboard)))
suppressWarnings(suppressPackageStartupMessages(library(shinyWidgets)))
suppressWarnings(suppressPackageStartupMessages(library(shinycssloaders)))
suppressWarnings(suppressPackageStartupMessages(library(jsonlite)))
suppressWarnings(suppressPackageStartupMessages(library(FuSViz)))
suppressWarnings(suppressPackageStartupMessages(library(DT)))
suppressWarnings(suppressPackageStartupMessages(library(BioCircos)))
suppressWarnings(suppressPackageStartupMessages(library(wordcloud2)))
suppressWarnings(suppressPackageStartupMessages(library(visNetwork)))
suppressWarnings(suppressPackageStartupMessages(library(datamods)))

#-------------------------------------
# set and intialize global variables
#-------------------------------------
#// System path of installed FuSViz package
system_path = path.package("FuSViz");
extdata = system.file("extdata", package = "FuSViz");

#// load canonical transcript id
load(file=file.path(extdata, "canonical.Rd"));

#// Load gene name data: gene_id a data.frame class with two columns: 'ensembl_gene_id' and 'gene_symbol'
load(file=file.path(extdata, "ensembl_symbol.Rd"));
#// assign the ensembl_id to gene symbol
symbol_ensem = gene_id$Gene_name;
names(symbol_ensem) = gene_id$Gene_ID;

#// chrom control for circular plot
chrom_cir = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY");

#// Load predefined gene dataset for oncogenes, tumor-suppress genes and cancer-related genes
load(file=file.path(extdata, "cancergenes.Rd"));
oncogenes = cancergenes$oncogene;
tumorsupress = cancergenes$tumorsuppress;
related = cancergenes$related;
#// define the default color for 'oncogene', 'tumorsuppress gene' and 'cancer-related gene'
onco_color = "#ff8566";
supp_color = "#00ccff";
rela_color = "#ffcc33";
not_onco_supp_rela = "black";
intergenic_color = "#262626";
other_color = "#f2f2f2";

#// laod drug target information
drug_target = readRDS(file=file.path(extdata, "opentargets_cancer_drugs_202102_simple.rds"));
