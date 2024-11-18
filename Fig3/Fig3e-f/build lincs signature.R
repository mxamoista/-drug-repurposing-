library(cmapR)
library(tidyverse)
library(cmapR)
library(dplyr)

col_meta <- read_gctx_meta("GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx", dim="col")
gctx_demo = parse_gctx("GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx", 
                       cid=1:4, rid=1:4)
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

#==========================
#Filter the dataframe with specific cell lines and treatment
#Test1 in C1
sig_info_selected <- sig_info_ordered[
  sig_info_ordered$cell_id %in% c("PC3") &
    sig_info_ordered$pert_itime  == "6 h",]

landmark_data <- read.table("GSE92742_Broad_LINCS_gene_info.txt", 
                            header = TRUE,  # If the first row is a header
                            sep = "\t",     # If it's tab-separated
                            quote = "",     # If there are no quoted fields
                            fill = TRUE,    # To handle uneven column numbers
                            comment.char = "") # If # is used in the data

library(signatureSearch)
gctx <- "./GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"
gctx2h5(gctx, cid=sig_info_selected$sig_id,
        h5file="./data/lincs.h5", overwrite=TRUE)
h5_file <- "./data/lincs.h5"
h5ls(h5_file)
full_assay <- h5read(h5_file, "assay")
full_assay_colnames <- h5read(h5_file, "colnames")
full_assay_rownames <- h5read(h5_file, "rownames")

# To access a lanmark gene subset of the lincs_signature_selected 
#lincs_signature_selected <- full_assay[1:978, ]
C1_lincs_signature_selected <- full_assay
#landmark_data <- landmark_data[landmark_data$pr_gene_id %in% full_assay_rownames[1:978], ]
C1_landmark_data <- landmark_data[landmark_data$pr_gene_id %in% full_assay_rownames, ]

# reorder landmark_data and sig_info_selected as the same as lincs_signature_selected
C1_landmark_data <- C1_landmark_data[match(full_assay_rownames, C1_landmark_data$pr_gene_id), ]
save(C1_landmark_data, file='./lincs_myself/C1_landmark_data_PC3_6h_allconcentration.RData')
C1_sig_info_selected <- sig_info_selected[match(full_assay_colnames, sig_info_selected$sig_id),]
save(C1_sig_info_selected, file='./lincs_myself/C1_sig_info_selected.RData')

# make the rownames and colnames of lincs_signature_selected
rownames(C1_lincs_signature_selected) <- C1_landmark_data$pr_gene_id
colnames(C1_lincs_signature_selected) <- rownames(C1_sig_info_selected)

save(C1_lincs_signature_selected, file='./lincs_myself/C1_lincs_signature_PC3_6h_allconcentration.RData')

#==========================
#Test1 in C4
sig_info_selected <- sig_info_ordered[
  sig_info_ordered$cell_id %in% c("MCF7") &
    sig_info_ordered$pert_itime  == "6 h",]

landmark_data <- read.table("GSE92742_Broad_LINCS_gene_info.txt", 
                            header = TRUE,  # If the first row is a header
                            sep = "\t",     # If it's tab-separated
                            quote = "",     # If there are no quoted fields
                            fill = TRUE,    # To handle uneven column numbers
                            comment.char = "") # If # is used in the data

library(signatureSearch)
gctx <- "./GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"
gctx2h5(gctx, cid=sig_info_selected$sig_id,
        h5file="./data/lincs.h5", overwrite=TRUE)
h5_file <- "./data/lincs.h5"
h5ls(h5_file)
full_assay <- h5read(h5_file, "assay")
full_assay_colnames <- h5read(h5_file, "colnames")
full_assay_rownames <- h5read(h5_file, "rownames")

# To access a lanmark gene subset of the lincs_signature_selected 
#lincs_signature_selected <- full_assay[1:978, ]
C4_lincs_signature_selected <- full_assay
#landmark_data <- landmark_data[landmark_data$pr_gene_id %in% full_assay_rownames[1:978], ]
C4_landmark_data <- landmark_data[landmark_data$pr_gene_id %in% full_assay_rownames, ]

# reorder landmark_data and sig_info_selected as the same as lincs_signature_selected
C4_landmark_data <- C4_landmark_data[match(full_assay_rownames, C4_landmark_data$pr_gene_id), ]
save(C4_landmark_data, file='./lincs_myself/C4_landmark_data_PC3_6h_allconcentration.RData')
C4_sig_info_selected <- sig_info_selected[match(full_assay_colnames, sig_info_selected$sig_id),]
save(C4_sig_info_selected, file='./lincs_myself/C4_sig_info_selected.RData')

# make the rownames and colnames of lincs_signature_selected
rownames(C4_lincs_signature_selected) <- C4_landmark_data$pr_gene_id
colnames(C4_lincs_signature_selected) <- rownames(C4_sig_info_selected)

save(C4_lincs_signature_selected, file='./lincs_myself/C4_lincs_signature_MCF7_6h_allconcentration.RData')

#identical(rownames(lincs_signatures), rownames(lincs_signature_selected)[1:4000])
#[1] TRUE
