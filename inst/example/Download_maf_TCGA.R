library(TCGAbiolinks)
query <- GDCquery(
    project = "TCGA-PRAD", 
    data.category = "Simple Nucleotide Variation", 
    access = "open",
    data.type = "Masked Somatic Mutation", 
    workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)
GDCdownload(query)
maf <- GDCprepare(query)

maf_dataframe <- as.data.frame(maf)
write.table(maf_dataframe, file="output.maf", sep="\t", col.names=TRUE, row.names=FALSE, quote=F)