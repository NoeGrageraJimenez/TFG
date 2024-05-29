if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")
library(biomaRt)


# 1. generate gene to transcript mapping
#
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "hsapiens_gene_ensembl",
                         host = "https://www.ensembl.org",
                         verbose = TRUE)
# attributes = listAttributes(mart)
# hgnc_symbol gives less than external_gene_name
working_attributes <- c("ensembl_gene_id",
                        "external_gene_name",
                        "gene_biotype",
                        "description",
                        "start_position",
                        "end_position")
t2g <- biomaRt::getBM(attributes = working_attributes,
                      mart = mart,
                      verbose = TRUE)
t2g$geneLength <- t2g$end_position - t2g$start_position
drop <- c("start_position", "end_position")
t2g <- t2g[, !(names(t2g) %in% drop)]
dim(t2g)
View(t2g)
write.csv(t2g, file = "C:/Users/NOE/Desktop/TFG-Example/Islet data/t2g.csv", row.names = FALSE)
