## Load packages ----
library(biomaRt)
library(tidyverse)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(msigdbr)
library(fgsea)

## Load data ----

### make metadata table
sra_table <- read.delim(file = "data/SraRunTable.txt", header = TRUE, sep = ",")
sra_table <- sra_table[, c("Sample.Name", "publication_id", "Age", "sex", "disease", "disease_type", "tissue_collection")]

### load expression data
files = list.files(path = "data/GSE235236_RAW/") # get the list of file names

mtx_list <- list()
for (i in 1:length(files)){
  mtx <- read.delim(file = paste0("data/GSE235236_RAW/", files[i]), header = TRUE, sep = "\t")
  mtx <- mtx[, c("gene_id", "TPM")] 
  mtx$gene_id <- sapply(strsplit(mtx$gene_id, split = "\\."), "[", 1) # take only ENSG gene_id
  mtx <- distinct(mtx, gene_id, .keep_all = TRUE)
  sample_id <- sapply(strsplit(files[i], split = "_"), "[", 1) # GSM number as sample_id
  colnames(mtx) <- c("gene_id", sample_id)
  mtx_list[[i]] <- mtx
}
cts <- purrr::reduce(mtx_list, inner_join, by = "gene_id") #combine all the mtx

### convert ENSG to gene symbol, remove genes without gene symbol, remove duplicated genes
mart <- useMart('ENSEMBL_MART_ENSEMBL')
mart <- useDataset('hsapiens_gene_ensembl', mart)
lookup <- getBM(mart = mart, attributes = c('hgnc_symbol','ensembl_gene_id'),uniqueRows = TRUE) # lookup table to map ids to names
gene_names <- lookup$hgnc_symbol[match(cts$gene_id, lookup$ensembl_gene_id)]
gene_names[gene_names == ""] <- NA
cts$gene_name <- gene_names
cts <- filter(cts, !is.na(gene_name)) %>% distinct(gene_name, .keep_all = TRUE)
rownames(cts) <- cts$gene_name
cts <- dplyr::select(cts, -c("gene_id", "gene_name"))

### make gct file

#### sample information from sra_table
coldata <- sra_table[, c("Sample.Name", "disease")]
coldata$disease <- factor(coldata$disease, levels = c("HC","UC", "CD"))    # make disease column factor with order of HC<UC<CD
coldata <- coldata[order(coldata$disease),]    #sort coldata rows by disease

#### format cts as gct
nsample <- nrow(coldata)
ngene <- nrow(cts)
cts <- cts[,coldata$Sample.Name]    # sort cts columns (i.e. samples) as the same order as coldata
colnames(cts) <- paste(coldata$disease, coldata$Sample.Name, sep = "_")    # prefix sample name with disease type
gct <- cbind(rownames(cts), rep("NA", ngene), cts)
header1 <- c("#1.2", rep("", nsample + 1))
header2 <- c(ngene, nsample, rep("", nsample))
header3 <- c(c("gene_name", "description"), colnames(cts))
gct <- rbind(header1, header2, header3, gct)
write.table(gct,"results/gene_expression/bulk_TPM.gct", sep ="\t", col.names = FALSE, row.names = FALSE)