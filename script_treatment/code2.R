library(SNPRelate)


# --- 1. SETUP & DATA LOADING ---
showfile.gds(closeall=TRUE) # Ensure handles are closed
vcf_fn <- "/courses/master_b2s/raw_data/s02_vc/treatment/Supp File 1 - Sniffles SV.vcf"
gds_fn <- "~/vcf_treamt/treatment01/multisample_12G.gds"

# Convert and open (using biallelic only for clean clustering signal)
# snpgdsVCF2GDS(vcf_fn, gds_fn, method="biallelic.only") 
genofile <- snpgdsOpen(gds_fn)

# --- 2. CALCULATE ALLELE FREQUENCY (AF) FOR ALL GROUPS ---
# results$cluster contains your 8 clusters from the Haplotype Clustering step
af_list <- list()
for(i in 1:8) {
  s_ids <- results$sample.id[results$cluster == i]
  # Calculate AF while allowing non-autosomes (Chr 5, 12, 17)
  af_data <- snpgdsSNPRateFreq(genofile, sample.id = s_ids, with.id=TRUE)
  af_list[[paste0("G", i)]] <- af_data$AlleleFreq
}

# Combine into a Master Comparison Matrix
af_matrix <- as.data.frame(do.call(cbind, af_list))
snps <- snpgdsSNPList(genofile)
af_matrix <- cbind(snps[,c("chromosome", "position")], af_matrix)

# --- 3. DEFINE "GOOD OUTCOME" BASELINE ---
# We combine Healthy (G1) and Success (G8) to see what 'normal' treatment looks like
af_matrix$good_baseline <- (af_matrix$G1 + af_matrix$G8) / 2

# --- 4. EXTRACT UNIQUE MARKERS FOR EACH GROUP ---
# We find SNPs where a group is DIFFERENT from the Good Baseline by > 40%
# And DIFFERENT from all other negative groups combined

unique_report <- list()

for(i in 2:7) {
  target_group <- paste0("G", i)
  other_neg_groups <- setdiff(paste0("G", 2:7), target_group)
  
  # Calculate the max frequency found in any OTHER negative group
  af_matrix$other_neg_max <- apply(af_matrix[, other_neg_groups], 1, max, na.rm=TRUE)
  
  # Logic for "Unique Marker":
  # 1. Significantly different from G1 & G8 (the healthy/success baseline)
  # 2. Significantly different from the other side-effect groups
  markers <- af_matrix[abs(af_matrix[[target_group]] - af_matrix$good_baseline) > 0.40 & 
                         abs(af_matrix[[target_group]] - af_matrix$other_neg_max) > 0.30, ]
  
  # Cleanup and save
  markers <- na.omit(markers)
  unique_report[[target_group]] <- markers
  write.csv(markers, paste0("Unique_Markers_", target_group, ".csv"), row.names=FALSE)
}

# View top markers for Group 8 (Success) vs Group 1 (Healthy)
efficacy_markers <- af_matrix[abs(af_matrix$G8 - af_matrix$G1) > 0.50, ]
print(head(efficacy_markers))