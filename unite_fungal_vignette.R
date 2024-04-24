# Sample filtering and trimming with DADA2 ----

# Import reference FASTAs
unite.its <- "./vignettes/sh_general_release_dynamic_29.11.2022.fasta" # ITS fungal reference FASTA

# Import sample metadata
meta <- read.csv("vignettes/Mapping file for ITS sequencing.csv")

# Adjust metadata formatting
colnames(meta) <- NULL
colnames(meta) <- meta[2,]
meta <- meta[-c(1,2),]

# Define filter path for the forward and reverse reads of each sample fastq file
seq_path <- "./vignettes/CC_Seq"  # replace with the directory for fastq file after unzipping the folder

# Sort forward and reverse reads to be in the same order
fnFs <- sort(list.files(seq_path, pattern = "_R1_001.fastq.gz"))
fnRs <- sort(list.files(seq_path, pattern = "_R2_001.fastq.gz"))

# Specify the full path to the fnFs and fnRs
fnFs <- file.path(seq_path, fnFs)
fnRs <- file.path(seq_path, fnRs)

# Extract sample names, assuming file paths have format: FOL_DER/SAMPLENAME_XXX.fastq.gz OR FOL_DER/SAMPLENAME_XXX.fastq
sNames <- sapply(strsplit(fnFs, "_"), `[`, 2)
sNames <- sapply(strsplit(sNames, "/"), `[`, 2)

## DON'T RUN ##
# Inspect the first 2 forward reads:
dada2::plotQualityProfile(fnFs[1:2])  ## forward reads maintain high throughput quality, trimmed at position 245
# Inspect the first 2 reverse reads:
dada2::plotQualityProfile(fnRs[1:2])  ## reverse read quality drops at ~position 160, trimmed at position 160
## first 10 bps also removed due to potential pathological issues

# Define file names for filtered fastq.gz files and assign file paths
filt_path <- file.path(seq_path, "filtered") # Place filtered files in new subdirectory if one does not already exist
if (!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sNames, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sNames, "_R_filt.fastq.gz"))

# Filter and trim files based on the previous plots
out <- dada2::filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                            truncLen = c(240,160),
                            maxN = 0,
                            maxEE = c(2,2),
                            truncQ = 2,
                            rm.phix = TRUE,
                            compress = TRUE,
                            multithread = FALSE) # On Windows set multithread = FALSE

# UNITE ITS general release fungal reference library ----

# Empty list for storing the results
UNITE.fungi <- list()

# Loop the analysis process for every file in the directory
for (i in 1:length(fnFs)){

  # Assign metadata to be input as a list for ps.net argument
  meta_data <- list()
  for (j in 1:length(fnFs)){

    ID <- sapply(strsplit(fnFs[i], "_"), `[`, 2)
    ID <- sapply(strsplit(ID, "/"), `[`, 2)
    meta_data[[j]] <- as.list(meta[meta$`sample-Id` %in% ID,])
  }

  # Create a new row of results for each output from the `ps.net`
  UNITE.fungi[[i]] <- mga::mga(fastq.Fs = fnFs[i],
                             fastq.Rs = fnRs[i], # file paths for forward and reverse raw fastq files
                             filtFs = filtFs[i],
                             filtRs = filtRs[i], # file paths for filtered fastq files
                             refFasta = unite.its,
                             metadata = meta_data[[j]],
                             make.unique = TRUE,
                             tree.args = list(k = 4,
                                              inv = 0.2,
                                              model = "GTR",
                                              rearrangement = "stochastic"),
                             network = TRUE,
                             network.args = list(type = "taxa",
                                                 distance = "jaccard",
                                                 max.dist = 0.35,
                                                 keep.isolates = TRUE))

  # Save the results as they are produced, in case of crashes/errors
  save(UNITE.fungi, file = "vignettes/RDS/unite.its_mga_results.rda")
}

# Extract sample results from the saved list
RES <- list()
for (i in 1:length(UNITE.fungi)){

  RES[[i]] <- UNITE.fungi[[i]]$results.samples
}

# Build combined data frame of metrics
results.df <- do.call(rbind, RES)
results.df <- dplyr::relocate(results.df, "sample.ID", .before = "Shannon") # move sample.ID column to leftmost

# Save data frame of results
saveRDS(results.df, file = "vignettes/RDS/results.df_unite.its_mga.rds")


## Plotting ----

# Import data
# load(file = "vignettes/RDS/silva132_mga_results.rda")

# Plot phylogenetic tree
png(filename = "vignettes/figures/unite.its_tree1.png", width = 6.86, height = 6.86, units = "in", res = 600)
plot(UNITE.fungi[[2]], type = "tree", cex = 0.25) # smaller tip label font
dev.off()

png(filename = "vignettes/figures/unite.its_tree2.png", width = 6.86, height = 6.86, units = "in", res = 600)
plot(UNITE.fungi[[2]], type = "tree", show.tip.label = FALSE) # remove tip labels
dev.off()

# Plot co-occurrence network
png(filename = "vignettes/figures/unite.its_net.png", width = 6.86, height = 6.86, units = "in", res = 600)
plot(UNITE.fungi[[2]], type = "network", point_size = 2, label = NULL) # remove labels and shrink points
dev.off()

## Taxonomic Filtering ----

# Import .csv file for taxa to keep
pathogens <- read.csv(file = "vignettes/Plant_Pathogens.csv")

# Filter for pathogens
test_filter <- list()
for (i in 1:length(UNITE.fungi)) {
  test_filter[[i]] <- mga::drop_taxa(UNITE.fungi[[i]], taxa = pathogens)
}

plot(test_filter[[2]], type = "tree", cex = 0.5)

save(test_filter, file = "vignettes/RDS/pathogens_unite.its.rda")


# Load fungal data and results
load(file = "vignettes/RDS/unite.its_mga_results.rda")
results.df4 <- readRDS("vignettes/RDS/results.df_unite.its_mga.rds")


# Empty list to store subset data
UNITE <- list()
# Subset list for only buckwheat and buffalo grass samples
for (i in 1:length(UNITE.fungi)) {

  if ("Buckwheat" %in% UNITE.fungi[[i]]$sampledata$cover_crop |
      "Buffalo Grass" %in% UNITE.fungi[[i]]$sampledata$cover_crop) {
    UNITE[[i]] <- UNITE.fungi[[i]]
  }
}

UNITE <- Filter(Negate(is.null), UNITE) # remove null elements in list

# Subset data frame
sub.df <- results.df4[results.df4$cover_crop == "Buckwheat" |
                        results.df4$cover_crop == "Buffalo Grass",]

# Filter for pathogens
test_filter <- list()
for (i in 1:length(UNITE)) {
  test_filter[[i]] <- mga::drop_taxa(UNITE[[i]], taxa = pathogens)
}

save(test_filter, file = "vignettes/RDS/pathogens_unite.its_fix.rda")
