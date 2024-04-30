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

#// define the default color for 'oncogene', 'tumorsuppress gene' and 'cancer-related gene'
onco_color = "#ff8566";
supp_color = "#00ccff";
rela_color = "#ffcc33";
not_onco_supp_rela = "black";
intergenic_color = "#262626";
other_color = "#f2f2f2";

#// laod drug target information
drug_target = readRDS(file=file.path(extdata, "opentargets_cancer_drugs_202102_simple.rds"));
