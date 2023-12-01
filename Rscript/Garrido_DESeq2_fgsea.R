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
sra_table <- dplyr::select(sra_table, Sample.Name, publication_id, Age, sex, disease, disease_type, tissue_collection)

### load expression data
files = list.files(path = "data/GSE235236_RAW/") # get the list of file names

mtx_list <- list()
for (i in 1:length(files)){
  mtx <- read.delim(file = paste0("data/GSE235236_RAW/", files[i]), header = TRUE, sep = "\t")
  mtx <- dplyr::select(mtx, gene_id, expected_count) 
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

## DESeq2 ----

### make coldata from sra_table
coldata <- dplyr::select(sra_table, Sample.Name, disease)
rownames(coldata) <- coldata$Sample.Name
coldata <- dplyr::select(coldata, disease)
colnames(coldata) <- "condition"
coldata$condition <- factor(coldata$condition)

### make sure cts and coldata are in the same order
all(rownames(coldata) %in% colnames(cts))
all(rownames(coldata) == colnames(cts))
cts <- cts[, rownames(coldata)]
all(rownames(coldata) == colnames(cts))

### make DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = round(cts), colData = coldata, design = ~condition)

#### pre-fileter rows with very few counts
smallestGroupSize <- 8 # there are 8 HC samples
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

#### set factor levels of condition
dds$condition <- factor(dds$condition, levels = c("HC","UC", "CD"))

### DE analysis
dds <- DESeq(dds)

### PCA plot
vst = varianceStabilizingTransformation(dds)
pcaData <- plotPCA(vst, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_bw()

## UC vs HC ----

### get result, volcano plot
res <- results(dds, contrast = c("condition", "UC", "HC")) %>% data.frame()
res <- rownames_to_column(res, var = "gene_name")
write.table(res,"results/DEG/bulk_DESeq2_UC_vs_HC.txt", sep ="\t", col.names = T, row.names = F)
ggplot(res, aes(log2FoldChange, -log10(pvalue))) + geom_point(size = 0.5, alpha = 0.5) + theme_bw() +
  ylab("-log10(unadjusted p-value)") + geom_text_repel(aes(label = ifelse(padj < 0.05, gene_name,"")), colour = "red", size = 3)

### GSEA by fgsea

#### make gene rank
rank <- res %>% 
  dplyr::dplyr::select(gene_name, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(gene_name) %>% 
  summarize(stat=mean(stat))
rank <- deframe(rank)

#### prepare gene sets
collections <- list()
collections$BIOCARTA <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "BIOCARTA")
collections$CGP <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CGP")
collections$HALLMARKS <- msigdbr(species = "Homo sapiens", category = "H")
collections$C6 <- msigdbr(species = "Homo sapiens", category = "C6")
collections <- lapply(collections, function(x) {
  out <- split(x = x$gene_symbol, f = x$gs_name)
})

#### run fgsea
fgseaRes <- fgsea(pathways = collections$C6, stats = rank, eps=0.0, minSize = 10, maxSize = 500)
fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES))
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.25)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()
fgseaRes$leadingEdge <- fgseaRes$leadingEdge %>% lapply(paste, collapse = ",") %>% unlist(recursive = FALSE)
write.table(fgseaRes,"results/DEG/bulk_fgsea_C6_UC_vs_HC.txt", sep ="\t", col.names = TRUE, row.names = FALSE)
# plotEnrichment(collections$C6[["EGFR_UP.V1_UP"]], ranks) + labs(title="C6_EGFR_UP.V1_UP")

## UC vs HC with log fold change shrinkage, whose usage is not clear ----
# resultsNames(dds)
# resLFC <- lfcShrink(dds, coef="condition_UC_vs_HC", type="apeglm") %>% data.frame()
# ggplot(resLFC, aes(log2FoldChange, -log10(pvalue))) + geom_point(size = 0.5, alpha = 0.5) + theme_bw() +
#   ylab("-log10(unadjusted p-value)") + geom_text_repel(aes(label = ifelse(padj < 0.05, rownames(resLFC),"")), colour = "red", size = 3)

## CD vs HC ----

### get result, volcano plot
res <- results(dds, contrast = c("condition", "CD", "HC")) %>% data.frame()
res <- rownames_to_column(res, var = "gene_name")
write.table(res,"results/DEG/bulk_DESeq2_CD_vs_HC2.txt", sep ="\t", col.names = T, row.names = F)
ggplot(res, aes(log2FoldChange, -log10(pvalue))) + geom_point(size = 0.5, alpha = 0.5) + theme_bw() +
  ylab("-log10(unadjusted p-value)") + geom_text_repel(aes(label = ifelse(padj < 0.05, gene_name,"")), colour = "red", size = 3)

### GSEA by fgsea

#### make gene rank
rank <- res %>% 
  dplyr::select(gene_name, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(gene_name) %>% 
  summarize(stat=mean(stat))
rank <- deframe(rank)

#### prepare gene sets
collections <- list()
collections$BIOCARTA <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "BIOCARTA")
collections$CGP <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CGP")
collections$HALLMARKS <- msigdbr(species = "Homo sapiens", category = "H")
collections$C6 <- msigdbr(species = "Homo sapiens", category = "C6")
collections <- lapply(collections, function(x) {
  out <- split(x = x$gene_symbol, f = x$gs_name)
})

#### run fgsea
fgseaRes <- fgsea(pathways = collections$HALLMARKS, stats = rank, eps=0.0, minSize = 10, maxSize = 500)
fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES))
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.25)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()
fgseaRes$leadingEdge <- fgseaRes$leadingEdge %>% lapply(paste, collapse = ",") %>% unlist(recursive = FALSE)
write.table(fgseaRes,"results/DEG/bulk_fgsea_H_CD_vs_HC_2.txt", sep ="\t", col.names = TRUE, row.names = FALSE)
# plotEnrichment(collections$C6[["EGFR_UP.V1_UP"]], ranks) + labs(title="C6_EGFR_UP.V1_UP")

## plotCounts ----
plotCounts(dds, gene = "OLFM4", intgroup = "condition")
