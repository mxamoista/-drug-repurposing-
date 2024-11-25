library(SingleCellExperiment)
library(SC3)
library(grid)
library(DESeq2)

#==========================
#We cited the code in  "https://github.com/Ellen1101/OA-subtype".
#Using "cartilage.count.201905113.txt" , OA patients were classified into four subtypes with distinct molecular signatures.
raw_matrix <- read.table("./Bone Research/cartilage.count.201905113.txt",header = TRUE , row.names =1)
length(colnames(raw_matrix))
raw_matrix <- raw_matrix[,6:ncol(raw_matrix)]
dim(raw_matrix)
colSums(raw_matrix)[order(colSums(raw_matrix))]
geneNumber <- colSums(raw_matrix>0)
geneNumber[order(geneNumber)]
sample_matrix <- read.table("./Bone Research/cart_ann_20190513.txt",header = TRUE , row.names =1)
sample_matrix <- sample_matrix[colnames(raw_matrix),]
dim(sample_matrix)
sampleInfor <- sample_matrix[order(rownames(sample_matrix)),]
treat <- rep("OA",nrow(sampleInfor))
sampleInfor <- cbind(sampleInfor,treat )
sampleInfor$treat <- as.character(sampleInfor$treat)
sampleInfor[c("P2_050","P2_197","P2_201","P2_205","P2_209"),]$treat<-rep("Normal",5)
# filter samples according totoal reads 
matrix <- raw_matrix[,colSums(raw_matrix) > 500000 ]
matrix <- matrix[,colSums(matrix) < 2000000]
dim(matrix)  #[1] 58684   135
matrix.nor <- t(t(matrix)/colSums(matrix)*1000000)    # calculate CPM
matrix  <- matrix.nor  

sce <- SingleCellExperiment(assays = list(counts = as.matrix(matrix),logcounts = log2(as.matrix(matrix) + 1)), colData = sampleInfor[colnames(matrix),] )
rowData(sce)$feature_symbol <- rownames(sce)

# filter genes var lage in normal
sampleInfor.normal <- sampleInfor[c("P2_197","P2_201","P2_205","P2_209"),]
sce.normal <- sce[,rownames(sampleInfor.normal)]
summary (apply(assay(sce.normal,"logcounts"), 1, var) )
sum(apply(assay(sce.normal,"logcounts"), 1, var) > 3)  #[1] 4423
drop_genes.01 <- apply(assay(sce.normal,"logcounts"), 1, function(x) {var(x) > 3 }) 
matrix.filter01 <- matrix[!drop_genes.01, ]       # filter genes var lage in normal
col_data <- colData(sce)

# filter samples 
matrix.filter02 <- matrix.filter01[,!(colnames(matrix.filter01) %in% c("P2_197","P2_201","P2_205","P2_209")) ]  # delet normal samples
dim(matrix.filter02)  #[1] 54261   131
new_sampleInfor <- col_data[colnames(matrix.filter02),]

matrix.filter03 <- matrix.filter02[rowSums(matrix.filter02>5)>15,]  # filter genes with low count
dim(matrix.filter03)    # [1] 9028  131

#==========================
varGene <- apply(matrix.filter03, 1, var)[order(apply(matrix.filter03, 1, var), decreasing = T)]
matrix.filter.final <- matrix.filter03[names(varGene[1:4000]),]

sce.final <- SingleCellExperiment(assays = list(counts = matrix.filter.final,logcounts = log2(as.matrix(matrix.filter.final) + 1)), colData = new_sampleInfor)
rowData(sce.final)$feature_symbol <- rownames(sce.final)
sce.final <- sc3_estimate_k(sce.final)
str(metadata(sce.final)$sc3)       # $ k_estimation: num 3
sce.final <- sc3(sce.final, ks = 2:8, biology = TRUE,n_cores=3,gene_filter = F) #
col_data <- colData(sce.final)
table(col_data[,"sc3_4_clusters"])   # 81 24 10 16

k<-4
dataset <- get_processed_dataset(sce.final)
hc <- metadata(sce.final)$sc3$consensus[[as.character(k)]]$hc
ht_ann <- col_data
ht_ann$sc3_4_clusters  <- sub("^1","C1",ht_ann$sc3_4_clusters)
ht_ann$sc3_4_clusters  <- sub("^2","C2",ht_ann$sc3_4_clusters)
ht_ann$sc3_4_clusters  <- sub("^3","C3",ht_ann$sc3_4_clusters)
ht_ann$sc3_4_clusters  <- sub("^4","C4",ht_ann$sc3_4_clusters)

annotation_col = data.frame(row.names = rownames(ht_ann), cluster=as.factor(ht_ann$sc3_4_clusters)) 

C1 <- annotation_col[annotation_col$cluster == "C1", , drop = FALSE]
C2 <- annotation_col[annotation_col$cluster == "C2", , drop = FALSE]
C3 <- annotation_col[annotation_col$cluster == "C3", , drop = FALSE]
C4 <- annotation_col[annotation_col$cluster == "C4", , drop = FALSE]

C1.matrix <- matrix.filter.final[,rownames(C1)]
C2.matrix <- matrix.filter.final[,rownames(C2)]
C3.matrix <- matrix.filter.final[,rownames(C3)]
C4.matrix <- matrix.filter.final[,rownames(C4)]

C1VSOthers.matrix <- cbind(C2.matrix,C3.matrix,C4.matrix)
C1VSOthers<-cbind(C1VSOthers.matrix ,C1.matrix)
matrix.roundC1VSOthers <- round(C1VSOthers)
conditionC1VSOthers <- factor(c(rep("Others",51),rep("C1",80)))
colDataC1VSOthers <- data.frame(
  sample_id = colnames(C1VSOthers),
  conditionC1VSOther = conditionC1VSOthers,
  stringsAsFactors = FALSE
)

# Ensure the condition column is a factor with "Others" as the reference level
colDataC1VSOthers$conditionC1VSOthers <- factor(colDataC1VSOthers$conditionC1VSOther, levels = c("Others", "C1"))

ddsC1VSOthers <- DESeqDataSetFromMatrix(countData = matrix.roundC1VSOthers,
                                        colData = colDataC1VSOthers,
                                        design = ~ conditionC1VSOthers)
# Run the DESeq analysis
ddsC1VSOthers <- DESeq(ddsC1VSOthers)
#Get the results, with "C1" compared to "Others" (the default behavior)
resC1VSOthers <- results(ddsC1VSOthers)
resC1VSOthers <- data.frame(resC1VSOthers, stringsAsFactors = FALSE, check.names = FALSE)# res format conversion: use data.frame to convert to table form
resC1VSOthers$gene <- rownames(resC1VSOthers)
resC1VSOthers <- subset(resC1VSOthers, padj < 1E-2 & abs(log2FoldChange) > 1 & abs(log2FoldChange) != Inf)

#Transfer gene into GeneID
library(biomaRt)
# Connect to the human genes dataset
mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
# Get the mapping
gene_ids_C1 <- getBM(attributes=c("hgnc_symbol", "entrezgene_id"),
                     filters="hgnc_symbol",
                     values=rownames(resC1VSOthers),
                     mart=mart)

# Merge the "gene_ids_C1" with DESeq2 result
C1_dz_signature <- merge(resC1VSOthers, gene_ids_C1, by.x="gene", by.y="hgnc_symbol", all.x=TRUE)
C1_dz_signature$entrezgene_id
C1_dz_signature <- subset(C1_dz_signature, !is.na(entrezgene_id) & entrezgene_id !='?')
C1_dz_signature <- subset(C1_dz_signature, select=c("entrezgene_id","gene", "log2FoldChange", "padj"))
names(C1_dz_signature) <- c("GeneID", "Symbol", "value", "padj")
C1_dz_signature <- subset(C1_dz_signature, !is.na(value))
C1_dz_signature <- C1_dz_signature[order(C1_dz_signature$value),]
C1_dz_signature$up_down <- "up"
C1_dz_signature$up_down[C1_dz_signature$value < 0] <- "down"
save(C1_dz_signature, file='./lincs_myself/github/C1_dz_signature.RData')

#gene.list
load(file = './lincs_myself/github/C1_lincs_experiments_PC3_6h_allconcentration.RData')
C1_gene.list <- C1_lincs_experiments$pr_gene_id

#load gist genes
C1_dz_signature <- subset(C1_dz_signature, GeneID %in% C1_gene.list)
C1_dz_genes_up <- subset(C1_dz_signature,up_down=="up",select="GeneID")
C1_dz_genes_down <- subset(C1_dz_signature,up_down=="down",select="GeneID")
save(C1_dz_signature, file='./lincs_myself/github/C1_dz_signature.RData')
