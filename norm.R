#### Load libraries ####

library(SummarizedExperiment)
library(TCGAbiolinks)
require(EDASeq)
require(dplyr)
require(NOISeq)
library(DESeq2)
library(biomaRt)

#### Prepare functions ####

filter_TCGA <- function(x) {dataFilt <- TCGAanalyze_Filtering(tabDF = x,
                                                              method = "quantile",
                                                              qnt.cut = 0.25)
threshold <- round(dim(x)[2]/2)
ridx <- rowSums(dataFilt == 0) <= threshold
dataFilt <- dataFilt[ridx, ]
ridx <- rowMeans(dataFilt) >= 10
dataFilt <- dataFilt[ridx, ]
x <- x[rownames(x) %in% rownames(dataFilt), ]
print(dim(x))
return(x)
}

get_annot <- function(x,y) {
  inter <- intersect(rownames(x), y$HGNC_symbol)
  length(inter)
  annot1 <- y[y$HGNC_symbol  %in% inter,]
  print(dim(annot1))
  annot1 <- annot1[!duplicated(annot1$HGNC_symbol),]
  print(dim(annot1))
  x <- x[rownames(x) %in% annot1$HGNC_symbol,]
  print(dim(x))
  x <- x[!duplicated(rownames(x)),] #
  print(dim(x))
  print(head(rownames(x)))
  print(head(annot1$HGNC_symbol))
  annot1 <- annot1[match(rownames(x), annot1$HGNC_symbol), ]
  print(dim(annot1))
  print(dim(x))
  print(head(rownames(x)))
  print(head(annot1$HGNC_symbol))
  return(list(x,annot1))
}

norm <- function(x, y, z) {
  ln.data <- withinLaneNormalization(x, y$Length, which = "full")
  gcn.data <- withinLaneNormalization(ln.data , y$GC, which = "full")
  Btwn.Norm <- betweenLaneNormalization(gcn.data, which = "full") 
  norm.counts <- tmm(Btwn.Norm, long = 1000, lc = 0, k = 0)
  noiseqData <- NOISeq::readData(norm.counts, factors = as.data.frame(z$Group))
  mydata2corr1 = NOISeq::ARSyNseq(noiseqData, norm = "n",  logtransf = FALSE)
  rnas2 <- exprs(mydata2corr1)
  return(rnas2)
}

Get_raw_matrix <- function(x,y) {
  z <- cbind(assay(x), assay(y))
  print(dim(z))
  print(head(rownames(z)))
  rownames(z) <- rowData(x)$gene_name
  print(head(rownames(z)))
  z <- z[!duplicated(rownames(z)),]
  print(dim(z))
  print(head(rownames(z)))
  return(z)
}

Get_factors_objects <- function(x,y){
  factors <- rbind(x, y)
  rownames(factors) <- factors$Sample
  return(factors)
}

get_norm_matrices <- function(x, a, y, z) { 
  rownames(x) <- a$ensembl_gene_id[match(rownames(x), a$HGNC_symbol)]
  x <- x[, y$Group==z]
  print(dim(x))
  x <- cbind(rownames(x), x)
  colnames(x)[1] <- "gene"
  print(dim(x))
  return(x)
}

#### Get annotation file ####
httr::set_config(httr::config(ssl_verifypeer = FALSE))
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "www")

features <- c("ensembl_gene_id", "chromosome_name", 
              "start_position", "end_position", "hgnc_symbol",	
              "percentage_gene_gc_content", "gene_biotype", "ensembl_gene_id_version", "hgnc_id")
chrs <- c(1:22, "X", "Y")

annot <- getBM(attributes = features,
               filters = "chromosome_name",
               values = chrs, 
               mart = ensembl)

colnames(annot)<-c("ensembl_gene_id", "Chr", "Start", "End", "HGNC_symbol", "GC", "Type", "Ensembl_ID_Version", "HGNC_ID")
annot$Length <- abs(annot$End - annot$Start)
annot <- annot[!duplicated(annot$Ensembl_ID_Version),]
dim(annot)
annot <- annot[!duplicated(annot$HGNC_symbol),]
dim(annot)

#### Get TCGA data ####

#In server, avoid download the data again
#setwd("/Users/kenzuke/Documents/R_projects/B_cell_malignancies/")

# Normal Bone Marrow

cases_doublezero <- readRDS("cases_doublezero.RDS")
qry.rna_Normal_BM <- GDCquery(project = "TARGET-AML",
                              data.category= "Transcriptome Profiling",
                              data.type = "Gene Expression Quantification",
                              workflow.type = "STAR - Counts",
                              sample.type = "Bone Marrow Normal",
                              barcode = cases_doublezero)
GDCdownload(qry.rna_Normal_BM)
Normal_BoneMarrow <- GDCprepare(qry.rna_Normal_BM, summarizedExperiment = TRUE) 
rm(qry.rna_Normal_BM)

# MM
#clinical_info <- readRDS("clinical_info_primary.RDS")
qry.rna_MM <- GDCquery(project = "MMRF-COMMPASS",
                           data.category= "Transcriptome Profiling",
                           data.type = "Gene Expression Quantification",
                           workflow.type = "STAR - Counts",
                           sample.type = c("Primary Blood Derived Cancer - Bone Marrow",
                                           "Recurrent Blood Derived Cancer - Bone Marrow")
)
GDCdownload(qry.rna_MM)
MM_BM <- GDCprepare(qry.rna_MM, summarizedExperiment = TRUE)
ALL_BM_noNA <- ALL_BM[,!(is.na(ALL_BM$sample_id))]
ALL_BM_noNA_1 <- ALL_BM_noNA[ , (which(!is.na(ALL_BM_noNA$primary_diagnosis)))]
MM_BM <- ALL_BM_noNA_1[ , ALL_BM_noNA_1$primary_diagnosis == "Precursor B-cell lymphoblastic leukemia"]
dim(MM_BM)
#rm(qry.rna_ALL_BM, ALL_BM_noNA, ALL_BM, ALL_BM_noNA_1)
\
MM_BM <- MM_BM[ , MM_BM$sample_type == "Primary Blood Derived Cancer - Bone Marrow"]
dim(MM_BM)

clinical <- as.data.frame(colData(MM_BM))
dim(clinical)
clinical <- clinical[!duplicated(clinical$patient),]
dim(clinical)

dim(MM_BM)
MM_BM <- MM_BM[,colnames(MM_BM) %in% clinical$barcode]
dim(MM_BM)

#Set the rigth directory
#setwd("/home/anakamura/SSN_pseudogenes")
# Error in setwd("/home/anakamura/SSN_pseudogenes") : 
#   cannot change working directory

##### Get raw matrices ####

MM_NBM <- Get_raw_matrix(MM_BM, Normal_BoneMarrow)

#### Get factors objects ####

factors_MM <- data.frame(Group = "MM", Sample = colnames(MM_BM))
factors_NBM <- data.frame(Group = "NormalBM", Sample = colnames(Normal_BoneMarrow))

factors_MM_NBM <- Get_factors_objects(factors_MM, factors_NBM)

#### Filter ####
MM_NBM <- filter_TCGA(MM_NBM)

#### Get annotation objects ####

MM_NBM <- get_annot(MM_NBM, annot)
saveRDS(MM_NBM[[1]], "MM_NBM_raw_counts.RDS")
saveRDS(MM_NBM[[2]], "MM_NBM_annot.RDS")

#### Normalize ####

MM_NBM_norm <- norm(MM_NBM[[1]], MM_NBM[[2]], factors_MM_NBM)

#### Get expressio matrices ####

MM_norm <- get_norm_matrices(MM_NBM_norm, MM_NBM[[2]], factors_MM_NBM, "MM")
NBM_MM_norm <- get_norm_matrices(MM_NBM_norm, MM_NBM[[2]], factors_MM_NBM, "NormalBM")

write.table(MM_norm, file = "rnas_norm_MM_only_primary.tsv", row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE)
write.table(NBM_MM_norm, file = "rnas_norm_NBM_MM_only_primary.tsv", row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE)
#saveRDS(MM_NBM[[2]], "/Users/kenzuke/Documents/R_projects/SSN_PSEUDOGENES/TEST_RUN/annot_MM.RDS")

MM_norm <- read.delim("/Users/kenzuke/Documents/R_projects/SSN_PSEUDOGENES/TEST_RUN/rnas_norm_MM.tsv")
rownames(MM_norm) <- MM_norm$X
MM_norm <- MM_norm[,-c(1,2)]
rownames(MM_norm) <- MM_annot$HGNC_symbol[match(rownames(MM_norm), MM_annot$ensembl_gene_id)]
MM_annot <- readRDS("/Users/kenzuke/Documents/R_projects/SSN_PSEUDOGENES/TEST_RUN/annot_MM.RDS")

#### DEG ####
dds <- DESeqDataSetFromMatrix(countData = round(MM_NBM[[1]]),
                              colData = factors_MM_NBM,
                              design = ~ Group)
dds <- DESeq(dds)
res <-  results(dds)
res
summary(res)
resLFC <- lfcShrink(dds, coef="Group_NormalBM_vs_MM", type="apeglm")
#El control va a la derecha, o sea el grupo referencia

write.table(resLFC, file = "/Users/kenzuke/Documents/R_projects/SSN_PSEUDOGENES/TEST_RUN/resLFC_NormalBM_vs_MM.tsv", row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE)

##### Filter pseudogenes #####

MM_annot <- MM_NBM[[2]]

pseudogene_in_data <- MM_annot[grep(".*pseudogene", MM_annot$Type),]


MM_norm <- MM_NBM_norm[, colnames(MM_NBM_norm) %in% factors_MM$Sample]
pseudogene_in_expr_matrix <- MM_norm[rownames(MM_norm) %in% pseudogene_in_data$HGNC_symbol,]

rownames(pseudogene_in_expr_matrix) <- MM_annot$HGNC_symbol[match(rownames(pseudogene_in_expr_matrix), MM_annot$ensembl_gene_id)]

saveRDS(pseudogene_in_expr_matrix, "pseudogene_in_expr_matrix.RDS")
write.table(pseudogene_in_expr_matrix, file = "rnas_pseudogenes_norm_MM.tsv", row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE)



ledges <- function(x){
  net <- rcorr(t(x), type = "spearman")
  
  coexpr_val <- net$r
  coexpr_val[lower.tri(coexpr_val, diag = TRUE)] <- NA
  
  coexpr_P <- net$P
  coexpr_P[lower.tri(coexpr_P, diag = TRUE)] <- NA
  
  #Converting the square matrix into a list of unique values
  val_melt <- melt(coexpr_val)
  val_melt <- val_melt[!is.na(val_melt$value), ]
  
  P_melt <- melt(coexpr_P)
  P_melt <- P_melt[!is.na(P_melt$value), ]
  
  colnames(val_melt) <- c("Source", "Target", "Sp_Coeff")
  val_melt <- val_melt %>% relocate(Sp_Coeff, .after=Source)
  Sp_matrix <- cbind(val_melt, P_melt$value)
  colnames(Sp_matrix) <- c("Source","Sp_Coeff", "Target", "P_value")
  
  Sp_matrix$abs <- abs(Sp_matrix$Sp_Coeff)
  sorted <- Sp_matrix[order(Sp_matrix$abs, decreasing = TRUE),]
  
  sorted <- sorted[sorted$P_value <= 0.00000001,]
  
  return(sorted)
  
}

pseudogene_in_expr_matrix <- readRDS("../pseudogene_in_expr_matrix.RDS")

pseudogene_network <- ledges(pseudogene_in_expr_matrix)
dim(pseudogene_network)

netFun_1 <- function(x, ...) {
  stats::cor(t(x), method="spearman") 
}

pseudogene_network_lioness <- lioness(pseudogene_in_expr_matrix, netFun_1)

pseudogene_in_expr_matrix_2 <- matrix(as.numeric(pseudogene_in_expr_matrix), nrow = nrow(pseudogene_in_expr_matrix), ncol = ncol(pseudogene_in_expr_matrix))
colnames(pseudogene_in_expr_matrix_2) <- colnames(pseudogene_in_expr_matrix)
rownames(pseudogene_in_expr_matrix_2) <- rownames(pseudogene_in_expr_matrix)

genes_in_network <- unique(c(pseudogene_network$Source, pseudogene_network$Target))
pseudogene_in_expr_matrix_2 <- pseudogene_in_expr_matrix_2[rownames(pseudogene_in_expr_matrix_2) %in% genes_in_network,]

pseudogene_network_lioness <- lioness(pseudogene_in_expr_matrix_2, netFun_1)
pseudogene_network_lioness <- readRDS("../../B_cell_malignancies/pseudogene_network_lioness.RDS")

##### Filter significant interactions #####

pseudogene_network$Source <- as.character(pseudogene_network$Source)
pseudogene_network$Target <- as.character(pseudogene_network$Target)

pseudogene_network$S_T <- paste( pmin(pseudogene_network$Source, pseudogene_network$Target), 
                                 pmax(pseudogene_network$Source, pseudogene_network$Target) 
                                ,sep="_")

edges <- intersect(rownames(pseudogene_network_lioness), pseudogene_network$S_T)

pseudogene_network_lioness_filtered <- pseudogene_network_lioness[rownames(pseudogene_network_lioness) %in% pseudogene_network$S_T  ,]

write.table(assay(pseudogene_network_lioness_filtered), file = "/Users/kenzuke/Documents/R_projects/SSN_PSEUDOGENES/TEST_RUN/pseudogene_network_lioness_filtered.tsv", row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)

pseudogene_network_lioness_filtered <- read.delim("pseudogene_network_lioness_filtered.tsv")
