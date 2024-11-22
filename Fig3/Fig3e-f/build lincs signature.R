library(cmapR)
library(tidyverse)
library(dplyr)

#Load differential gene expression signatures
col_meta <- read_gctx_meta("GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx", dim="col")

#Load Metadata for Level 5 data 
sig_info <- read.table("GSE92742_Broad_LINCS_sig_info.txt", 
                       header = TRUE, 
                       sep = "\t", 
                       stringsAsFactors = FALSE,
                       quote = "", 
                       fill = TRUE, 
                       na.strings = c("NA", ""))

# Create a match vector
match_vector <- match(col_meta$id, sig_info$sig_id)

# Reorder sig_info
sig_info_ordered <- sig_info[match_vector, ]

#Load Metadata for rows / genes of matrices
lincs_experiments <- read.table("GSE92742_Broad_LINCS_gene_info.txt", 
                            header = TRUE,  # If the first row is a header
                            sep = "\t",     # If it's tab-separated
                            quote = "",     # If there are no quoted fields
                            fill = TRUE,    # To handle uneven column numbers
                            comment.char = "") # If # is used in the data

library(signatureSearch)
gctx <- "./GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"

#==========================
#Filter the dataframe with specific cell lines and treatment
#Test1 in C1
C1_sig_info_selected <- sig_info_ordered[
  sig_info_ordered$cell_id %in% c("PC3") &
    sig_info_ordered$pert_itime  == "6 h",]

gctx2h5(gctx, cid=C1_sig_info_selected$sig_id,
        h5file="./data/C1_lincs.h5", overwrite=TRUE)
C1_h5_file <- "./data/C1_lincs.h5"
h5ls(C1_h5_file)
C1_full_assay <- h5read(C1_h5_file, "assay")
C1_full_assay_colnames <- h5read(C1_h5_file, "colnames")
C1_full_assay_rownames <- h5read(C1_h5_file, "rownames")

# To access C1_lincs_experiments according to C1_sig_info_selected
C1_lincs_signature_selected <- C1_full_assay
# Reorder and select lincs_experiments and sig_info_selected as the same as C1_lincs_signature_selected
C1_lincs_experiments <- lincs_experiments[match(C1_full_assay_rownames,lincs_experiments$pr_gene_id), ]#Select gene_id
save(C1_lincs_experiments, file='./lincs_myself/github/C1_lincs_experiments_PC3_6h_allconcentration.RData')
C1_sig_info_selected <- C1_sig_info_selected[match(C1_full_assay_colnames, C1_sig_info_selected$sig_id),]#Select sig_id
save(C1_sig_info_selected, file='./lincs_myself/github/C1_sig_info_selected.RData')

# Modify and add the rownames and colnames of C1_lincs_signature_selected
identical(as.integer(as.vector(C1_full_assay_rownames)),C1_lincs_experiments$pr_gene_id) #[1] TRUE
identical(as.character(as.vector(C1_full_assay_colnames)),C1_sig_info_selected$sig_id) #[1] TRUE
rownames(C1_lincs_signature_selected) <- C1_lincs_experiments$pr_gene_id
colnames(C1_lincs_signature_selected) <- rownames(C1_sig_info_selected)

save(C1_lincs_signature_selected, file='./lincs_myself/github/C1_lincs_signature_PC3_6h_allconcentration.RData')

#==========================
#Test1 in C4
C4_sig_info_selected <- sig_info_ordered[
  sig_info_ordered$cell_id %in% c("MCF7") &
    sig_info_ordered$pert_itime  == "6 h",]

gctx2h5(gctx, cid=C4_sig_info_selected$sig_id,
        h5file="./data/C4_lincs.h5", overwrite=TRUE)
C4_h5_file <- "./data/C4_lincs.h5"
h5ls(C4_h5_file)
C4_full_assay <- h5read(C4_h5_file, "assay")
C4_full_assay_colnames <- h5read(C4_h5_file, "colnames")
C4_full_assay_rownames <- h5read(C4_h5_file, "rownames")

# To access C4_lincs_experiments according to C4_sig_info_selected
C4_lincs_signature_selected <- C4_full_assay

#Reorder lincs_experiments and sig_info_selected as the same as lincs_signature_selected
C4_lincs_experiments <- lincs_experiments[match(C4_full_assay_rownames, lincs_experiments$pr_gene_id), ]
save(C4_lincs_experiments, file='./lincs_myself/github/C4_lincs_experiments_MCF7_6h_allconcentration.RData')
C4_sig_info_selected <- C4_sig_info_selected[match(C4_full_assay_colnames, C4_sig_info_selected$sig_id),]
save(C4_sig_info_selected, file='./lincs_myself/github/C4_sig_info_selected.RData')

# Modify and add the rownames and colnames of C1_lincs_signature_selected
identical(as.integer(as.vector(C4_full_assay_rownames)),C4_lincs_experiments$pr_gene_id) #[1] TRUE
identical(as.character(as.vector(C4_full_assay_colnames)),C4_sig_info_selected$sig_id) #[1] TRUE
rownames(C4_lincs_signature_selected) <- C4_lincs_experiments$pr_gene_id
colnames(C4_lincs_signature_selected) <- rownames(C4_sig_info_selected)

save(C4_lincs_signature_selected, file='./lincs_myself/github/C4_lincs_signature_MCF7_6h_allconcentration.RData')

