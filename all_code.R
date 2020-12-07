library(readr)
library(dplyr)
library(ggplot2)
library(WGCNA)

########## READING IN P-VALUES AND FOLD CHANGE OF TOP 10000 DIFERENTIALLY EXPRESSED PROBES ########## 

# SNpc target, VTA contrast
files_1 <- list.files(path = "/Users/rezaadibnia/Desktop/PD NM Project/mDBR - VTA vs SNpc/SNpc target, VTA contrast/", pattern = "Probes_.*csv", full.names = TRUE, ignore.case = FALSE)
SNt_VTAc_probe_list <- lapply(files_1, read_csv)
SNt_VTAc_probes <- data.frame(bind_rows(SNt_VTAc_probe_list))

# VTA target, SNpc contrast
files_2 <- list.files(path = "/Users/rezaadibnia/Desktop/PD NM Project/mDBR - VTA vs SNpc/VTA target, SNpc contrast/", pattern = "Probes_.*csv", full.names = TRUE, ignore.case = FALSE)
VTAt_SNc_probe_list <- lapply(files_2, read_csv)
VTAt_SNc_probes <- data.frame(bind_rows(VTAt_SNc_probe_list))

########## CORRECTING P-VALUES AND TRANSFORMING FOLD CHANGES ########## 

# BH-correcting p-values (Benjamini-Hochberg Procedure)
SNt_VTAc_probes$p <- p.adjust(SNt_VTAc_probes$p, method="BH")
VTAt_SNc_probes$p <- p.adjust(VTAt_SNc_probes$p, method="BH")

# Converting VTA target, SNpc contrast fold change to SNpc target, VTA contrast fold change
VTAt_SNc_probes$fold.change <- 1/(VTAt_SNc_probes$fold.change)

# Combining the data frames containing p values and fold change
p_and_FC <- data.frame(bind_rows(SNt_VTAc_probes, VTAt_SNc_probes))
p_and_FC <- p_and_FC[,c(4,10:11)]

########## VISUALIZING DIFFERENTIALLY EXPRESSED PROBES - VOLCANO PLOT ########## 

# Vertical lines for log2(fold change) thresholds, horizontal line for the p-value threshold
# Top 10% of upregulated probes and bottom 10% of downregulated probes are chosen
horz_min = quantile(log2(VTAt_SNc_probes$fold.change), 0.10)
horz_max = quantile(log2(SNt_VTAc_probes$fold.change), 0.90)

# Add a column to p_and_FC data frame indicating whether probe is upregulated, downregulated, or not significant
p_and_FC$differentially_expressed <- "Not Significant"
# If log2Foldchange > horz_max and pvalue < 0.05, set as "UP" 
p_and_FC$differentially_expressed[log2(p_and_FC$fold.change) > horz_max & p_and_FC$p < 0.05] <- "Upregulated"
# If log2Foldchange < horz_min and pvalue < 0.05, set as "DOWN"
p_and_FC$differentially_expressed[log2(p_and_FC$fold.change) < horz_min & p_and_FC$p < 0.05] <- "Downregulated"

# Construct volcano plot
volcano_plot <- ggplot(data=p_and_FC, aes(x=log2(fold.change), y=-log10(p), col=differentially_expressed)) + 
  geom_point() + 
  theme_minimal() +
  geom_vline(xintercept=c(horz_min, horz_max), col="black") +
  geom_hline(yintercept=-log10(0.05), col="black")

# The significantly differentially expressed probes are  found in the upper-left and upper-right corners of the volcano plot
volcano_plot


########## FILTERING FOR P < 0.05 and FC > 90th PERCENTILE (TOP 10%) ########## 

# Filters for top 10% of genes with highest FC compared to all probes, not just significantly differentially expressed prboes VTA target, SNpc contrast is filtered for fold changes less than the 10th percentile because fold change are all less than 1. The lowest 10% of fold change values indicate the top 10% of probes that are the most downregulated in the SNpc compared to the VTA.

# SNpc target, VTA contrast
SNt_VTAc_probes <- SNt_VTAc_probes %>% filter(SNt_VTAc_probes$p < 0.05 & SNt_VTAc_probes$fold.change > quantile(SNt_VTAc_probes$fold.change, 0.90))

# VTA target, SNpc contrast
VTAt_SNc_probes <- VTAt_SNc_probes %>% filter(VTAt_SNc_probes$p < 0.05 & VTAt_SNc_probes$fold.change < quantile(VTAt_SNc_probes$fold.change, 0.10))

########## READING IN EXPRESSION DATA FOR TOP 10000 DIFFERENTIALLY EXPRESSED PROBES ########## 

read_csv_labeled_header <- function(file) {
  labels = c("Probe_IDs", "SN_1", "SN_2", "SN_3", "SN_4", "SN_5", "SN_6", "VTA_1", "VTA_2", "VTA_3", "VTA_4", "VTA_5", "VTA_6")
  read_csv(file, col_names = labels)
}

# SNpc target, VTA contrast
files_expr_1 <- list.files(path = "/Users/rezaadibnia/Desktop/PD NM Project/mDBR - VTA vs SNpc/SNpc target, VTA contrast/", pattern = "Expression_.*csv", full.names = TRUE, ignore.case = FALSE)
SNt_VTAc_expr_list <- lapply(files_expr_1, read_csv_labeled_header)
SNt_VTAc_expr <- data.frame(bind_rows(SNt_VTAc_expr_list))

# VTA target, SNpc contrast
files_expr_2 <- list.files(path = "/Users/rezaadibnia/Desktop/PD NM Project/mDBR - VTA vs SNpc/VTA target, SNpc contrast/", pattern = "Expression_.*csv", full.names = TRUE, ignore.case = FALSE)
VTAt_SNc_expr_list <- lapply(files_expr_2, read_csv_labeled_header)
VTAt_SNc_expr <- data.frame(bind_rows(VTAt_SNc_expr_list))

########## JOINING EXPRESSION AND PROBE DATA FRAMES INTO ONE ########## 

SNt_VTAc_expr <- SNt_VTAc_expr %>% inner_join(SNt_VTAc_probes, by = c("Probe_IDs" = "id"))
VTAt_SNc_expr <- VTAt_SNc_expr %>% inner_join(VTAt_SNc_probes, by = c("Probe_IDs" = "id"))

expr <- data.frame(bind_rows(SNt_VTAc_expr, VTAt_SNc_expr))


########## ALIGNING ALL PROBES TO GENES ########## 

# Map probes to genes
# First download 'Complete normalized microarray datasets' for all six donors from http://human.brain-map.org/static/download
# Store the downloaded files in a AHBA download directory

# Set donor names
donorNames <- c("donor9861", "donor10021", "donor12876", "donor14380", "donor15496", "donor15697")
names(donorNames) <- donorNames

# Location of AHBA directories
ahba_download <- "/Users/rezaadibnia/Desktop/AHBA Full Data/Microarray" # Where downloaded files are stored

# Read in probe info (files are the same for each donor)
probeInfo <- read.csv(paste0(ahba_download, "/normalized_microarray_donor9861/Probes.csv")) # Same for each donor

# Read expression data for each donor (probes x samples)
brainExpr <- lapply(donorNames, function(d){
  file1 <- paste0(ahba_download, "/normalized_microarray_", d, "/MicroarrayExpression.csv")
  e <- read.csv(file1, header = FALSE)
  rownames(e) <- e[,1]
  e <- e[,-1]
  file2 <- paste0(ahba_download, "/normalized_microarray_", d, "/SampleAnnot.csv")
  sample_annotation <- read.csv(file2)
  colnames(e) <- sample_annotation$structure_id
  e
})

# Probes with missing Entrez IDs
probes_missing_entrezID <- is.na(probeInfo$entrez_id) # same in all donors
print(paste(sum(probes_missing_entrezID), "probes with missing entrez IDs"))

# Probes with expression well above background in at least 1% of samples in all donors
pa_call <- lapply(donorNames, function(d){
  file <- paste0(ahba_download, "/normalized_microarray_", d, "/PACall.csv")  
  pa <- read.csv(file, header = FALSE)
  rownames(pa) <- pa[,1]
  pa[,-1]
})
pa_concat <- Reduce(cbind, pa_call)
sum <- rowSums(pa_concat)
low_presence <- sum < 0.01*ncol(pa_concat)
print(paste(sum(low_presence), "probes present in <1% of samples"))

# Concatenate data across all donors and filter probes by ID
expr_concat <- Reduce(cbind, brainExpr)
expr_concat <- expr_concat[!(probes_missing_entrezID | low_presence),]
probeInfo <- probeInfo[!(probes_missing_entrezID | low_presence), ]

probe2gene <- collapseRows(expr_concat, 
                           rowGroup = probeInfo$entrez_id, 
                           rowID = probeInfo$probe_id,
                           method = "maxRowVariance",
                           connectivityBasedCollapsing = TRUE
)
selected_probes <- probe2gene$selectedRow
entrez_id <- probeInfo$entrez_id[selected_probes]
entrez_order <- order(entrez_id)

probeInfo <- probeInfo[selected_probes, ]
probeInfo <- probeInfo[entrez_order,]

########## ALIGNING SIGNIFICANTLY DIFFERENTIALLY EXPRESSED PROBES TO GENES ########## 

# Combine data frames by probe ID to filter out all probes not corresponding to a gene
final_expr <- expr %>% inner_join(probeInfo, by = c("Probe_IDs" = "probe_id"))

# Remove unnecessary and repeated columns
final_expr <- final_expr[,-c(14, 18:21, 24:29)]
final_expr$fold.change <- sort(final_expr$fold.change, decreasing = TRUE)

final_expr[,c(15, 17:18)]



