# Sample filtering and trimming with DADA2 ----

# Import reference FASTAs
silva.138.1 <- "./vignettes/silva_nr99_v138.1_train_set.fa.gz" # 16S silva reference FASTA
silva.138.1S <- "./vignettes/silva_nr99_v138.1_wSpecies_train_set.fa.gz" # 16S silva reference FASTA with species
silva.132 <- "./vignettes/silva_nr_v132_train_set.fa.gz" # old 16S silva reference FASTA
GG.13.8 <- "./vignettes/gg_13_8_train_set_97.fa.gz" # 16S greengenes reference FASTA

# Import sample metadata
meta <- read.csv("https://raw.githubusercontent.com/spholmes/F1000_workflow/master/data/MIMARKS_Data_combined.csv")

# Reformat sample names in metadata
meta$sample.ID <- paste0(gsub("00", "", meta$host_subject_id), "D", meta$age-21)
meta <- meta[!duplicated(meta$sample.ID),]  # Remove duplicate entries for reverse reads
meta <- meta[,-c(1:30,45,46)]  # Remove empty columns

# Define filter path for the forward and reverse reads of each sample fastq file
seq_path <- "./vignettes/MiSeq_SOP"  # replace with the directory for fastq file after unzipping the folder

# Sort forward and reverse reads to be in the same order
fnFs <- sort(list.files(seq_path, pattern = "_R1_001.fastq"))
fnRs <- sort(list.files(seq_path, pattern = "_R2_001.fastq"))

# Remove "Mock" fastq files
fnFs <- fnFs[-20]
fnRs <- fnRs[-20]

# Specify the full path to the fnFs and fnRs
fnFs <- file.path(seq_path, fnFs)
fnRs <- file.path(seq_path, fnRs)

# Extract sample names, assuming file paths have format: FOL_DER/SAMPLENAME_XXX.fastq.gz OR FOL_DER/SAMPLENAME_XXX.fastq
sNames <- sapply(strsplit(fnFs, "_"), `[`, 2)
sNames <- sapply(strsplit(sNames, "/"), `[`, 2)

## DON'T RUN ##
# Inspect the first 2 forward reads:
# dada2::plotQualityProfile(fnFs[1:2])  ## forward reads maintain high throughput quality, trimmed at position 245
# Inspect the first 2 reverse reads:
# dada2::plotQualityProfile(fnRs[1:2])  ## reverse read quality drops at ~position 160, trimmed at position 160
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

# 16S Silva v132 reference library ----

# Empty list for storing the results
SILVA.OLD <- list()

# Loop the analysis process for every file in the directory
for (i in 1:length(fnFs)){

  # Assign metadata to be input as a list for ps.net argument
  meta_data <- list()
  for (j in 1:length(fnFs)){

    ID <- sapply(strsplit(fnFs[i], "_"), `[`, 2)
    ID <- sapply(strsplit(ID, "/"), `[`, 2)
    meta_data[[j]] <- as.list(meta[meta$sample.ID %in% ID,])
  }

  # Create a new row of results for each output from the `ps.net`
  SILVA.OLD[[i]] <- mga::mga(fastq.Fs = fnFs[i],
                            fastq.Rs = fnRs[i], # file paths for forward and reverse raw fastq files
                            filtFs = filtFs[i],
                            filtRs = filtRs[i], # file paths for filtered fastq files
                            refFasta = silva.132,
                            metadata = meta_data[[j]])

  # Save the results as they are produced, in case of crashes/errors
  save(SILVA.OLD, file = "vignettes/RDS/silva132_mga_results.rda")
}

# Extract sample results from the saved list
RES <- list()
for (i in 1:length(SILVA.OLD)){

  RES[[i]] <- SILVA.OLD[[i]]$results.samples
}

# Build combined data frame of metrics
results.df <- do.call(rbind, RES)
results.df <- dplyr::relocate(results.df, "sample.ID", .before = "Shannon") # move sample.ID column to leftmost

# Add additional metadata
meta2 <- read.delim(file = "vignettes/MiSeq_SOP/mouse.dpw.metadata")
results.df$dpw <- meta2$dpw
results.df$dpwgroup <- NA
results.df[results.df$dpw <= 9,]$dpwgroup <- "0-9"
results.df[results.df$dpw >= 141,]$dpwgroup <- "141-150"

# Save data frame of results
saveRDS(results.df, file = "vignettes/RDS/results.df_silva132_mga.rds")


## Plotting ----

# Import data
# load(file = "vignettes/RDS/silva132_mga_results.rda")

# Plot phylogenetic tree
png(filename = "vignettes/figures/silva132_tree1.png", width = 6.86, height = 6.86, units = "in", res = 600)
plot(SILVA.OLD[[2]], type = "tree", cex = 0.25) # smaller tip label font
dev.off()

png(filename = "vignettes/figures/silva132_tree2.png", width = 6.86, height = 6.86, units = "in", res = 600)
plot(SILVA.OLD[[2]], type = "tree", show.tip.label = FALSE) # remove tip labels
dev.off()

# Plot co-occurrence network
png(filename = "vignettes/figures/silva132_net.png", width = 6.86, height = 6.86, units = "in", res = 600)
plot(SILVA.OLD[[2]], type = "network", point_size = 2, label = NULL) # remove labels and shrink points
dev.off()


# 16S Silva NR99 v138.1 reference library with species ----

# Empty list for storing the results
SILVA.NEWS <- list()

# Loop the analysis process for every file in the directory
for (i in 1:length(fnFs)){

  # Assign metadata to be input as a list for ps.net argument
  meta_data <- list()
  for (j in 1:length(fnFs)){

    ID <- sapply(strsplit(fnFs[i], "_"), `[`, 2)
    ID <- sapply(strsplit(ID, "/"), `[`, 2)
    meta_data[[j]] <- as.list(meta[meta$sample.ID %in% ID,])
  }

  # Create a new row of results for each output from the `ps.net`
  SILVA.NEWS[[i]] <- mga::mga(fastq.Fs = fnFs[i],
                             fastq.Rs = fnRs[i], # file paths for forward and reverse raw fastq files
                             filtFs = filtFs[i],
                             filtRs = filtRs[i], # file paths for filtered fastq files
                             refFasta = silva.138.1S, # silva ref FASTA with species
                             metadata = meta_data[[j]])

  # Save the results as they are produced, in case of crashes/errors
  save(SILVA.NEWS, file = "vignettes/RDS/silva138.1S_mga_results.rda")
}

# Extract sample results from the saved list
RES2 <- list()
for (i in 1:length(SILVA.NEWS)){

  RES2[[i]] <- SILVA.NEWS[[i]]$results.samples
}

# Build combined data frame of metrics
results.df2 <- do.call(rbind, RES2)
results.df2 <- dplyr::relocate(results.df2, "sample.ID", .before = "Shannon") # move sample.ID column to leftmost

# Add additional metadata
results.df2$dpw <- meta2$dpw
results.df2$dpwgroup <- NA
results.df2[results.df2$dpw <= 9,]$dpwgroup <- "0-9"
results.df2[results.df2$dpw >= 141,]$dpwgroup <- "141-150"

# Save data frame of results
saveRDS(results.df2, file = "vignettes/RDS/results.df_silva138.1S_mga.rds")


## Plotting ----

# Import data
# load(file = "vignettes/RDS/silva138.1S_mga_results.rda")

# Plot phylogenetic tree
png(filename = "vignettes/figures/silva138.1S_tree1.png", width = 6.86, height = 6.86, units = "in", res = 600)
plot(SILVA.NEWS[[2]], type = "tree", cex = 0.25) # smaller tip label font
dev.off()

png(filename = "vignettes/figures/silva138.1S_tree2.png", width = 6.86, height = 6.86, units = "in", res = 600)
plot(SILVA.NEWS[[2]], type = "tree", show.tip.label = FALSE) # remove tip labels
dev.off()

# Plot co-occurrence network
png(filename = "vignettes/figures/silva138.1S_net.png", width = 6.86, height = 6.86, units = "in", res = 600)
plot(SILVA.NEWS[[2]], type = "network", point_size = 2, label = NULL) # remove labels and shrink points
dev.off()


# 16S Silva NR99 v138.1 reference library (no species, no network) ----

# Empty list for storing the results
SILVA.NEW <- list()

# Loop the analysis process for every file in the directory
for (i in 1:length(fnFs)){

  # Assign metadata to be input as a list for ps.net argument
  meta_data <- list()
  for (j in 1:length(fnFs)){

    ID <- sapply(strsplit(fnFs[i], "_"), `[`, 2)
    ID <- sapply(strsplit(ID, "/"), `[`, 2)
    meta_data[[j]] <- as.list(meta[meta$sample.ID %in% ID,])
  }

  # Create a new row of results for each output from the `ps.net`
  SILVA.NEW[[i]] <- mga::mga(fastq.Fs = fnFs[i],
                             fastq.Rs = fnRs[i], # file paths for forward and reverse raw fastq files
                             filtFs = filtFs[i],
                             filtRs = filtRs[i], # file paths for filtered fastq files
                             refFasta = silva.138.1, # silva ref FASTA no species
                             metadata = meta_data[[j]],
                             network = FALSE)

  # Save the results as they are produced, in case of crashes/errors
  save(SILVA.NEW, file = "vignettes/RDS/silva138.1_mga_noNet.rda")
}

# Extract sample results from the saved list
RES3 <- list()
for (i in 1:length(SILVA.NEW)){

  RES3[[i]] <- SILVA.NEW[[i]]$results.samples
}

# Build combined data frame of metrics
results.df3 <- do.call(rbind, RES3)
results.df3 <- dplyr::relocate(results.df3, "sample.ID", .before = "Shannon") # move sample.ID column to leftmost

# Add additional metadata
results.df3$dpw <- meta2$dpw
results.df3$dpwgroup <- NA
results.df3[results.df3$dpw <= 9,]$dpwgroup <- "0-9"
results.df3[results.df3$dpw >= 141,]$dpwgroup <- "141-150"

# Save data frame of results
saveRDS(results.df3, file = "vignettes/RDS/results.df_silva138.1_mga_noNet.rds")


## Plotting ----

# Import data
# load(file = "vignettes/RDS/silva138.1_mga_noNet.rda")

# Plot phylogenetic tree
png(filename = "vignettes/figures/silva138.1_tree1.png", width = 6.86, height = 6.86, units = "in", res = 600)
plot(SILVA.NEW[[2]], type = "tree", cex = 0.25) # smaller tip label font
dev.off()

png(filename = "vignettes/figures/silva138.1_tree2.png", width = 6.86, height = 6.86, units = "in", res = 600)
plot(SILVA.NEW[[2]], type = "tree", show.tip.label = FALSE) # remove tip labels
dev.off()

# Plot co-occurrence network -- WILL GIVE WARNING (no network to plot)
plot(SILVA.NEW[[2]], type = "network", point_size = 2, label = NULL) # remove labels and shrink points


# 16S Greengenes 12.10 reference library ----

# Empty list for storing the results
GREEN <- list()

# Loop the analysis process for every file in the directory
for (i in 1:length(fnFs)){

  # Assign metadata to be input as a list for ps.net argument
  meta_data <- list()
  for (j in 1:length(fnFs)){

    ID <- sapply(strsplit(fnFs[i], "_"), `[`, 2)
    ID <- sapply(strsplit(ID, "/"), `[`, 2)
    meta_data[[j]] <- as.list(meta[meta$sample.ID %in% ID,])
  }

  # Create a new row of results for each output from the `ps.net`
  GREEN[[i]] <- mga::mga(fastq.Fs = fnFs[i],
                         fastq.Rs = fnRs[i], # file paths for forward and reverse raw fastq files
                         filtFs = filtFs[i],
                         filtRs = filtRs[i], # file paths for filtered fastq files
                         refFasta = GG.13.8, # greengenes ref FASTA
                         metadata = meta_data[[j]])

  # Save the results as they are produced, in case of crashes/errors
  save(GREEN, file = "vignettes/RDS/gg13.8_mga_results.rda")
}

# Extract sample results from the saved list
RES4 <- list()
for (i in 1:length(GREEN)){

  RES4[[i]] <- GREEN[[i]]$results.samples
}

# Build combined data frame of metrics
results.df4 <- do.call(rbind, RES2)
results.df4 <- dplyr::relocate(results.df4, "sample.ID", .before = "Shannon") # move sample.ID column to leftmost

# Add additional metadata
results.df4$dpw <- meta2$dpw
results.df4$dpwgroup <- NA
results.df4[results.df4$dpw <= 9,]$dpwgroup <- "0-9"
results.df4[results.df4$dpw >= 141,]$dpwgroup <- "141-150"

# Save data frame of results
saveRDS(results.df4, file = "vignettes/RDS/results.df_gg13.8_mga.rds")


# mga(fastq.Fs = fnFs[i],
#     fastq.Rs = fnRs[i],
#     filtFs = filtFs[i],
#     filtRs = filtRs[i],
#     refFasta = silva.138.1S,
#     metadata = meta_data[[j]])


# Single sample ----

unite.its <- "./vignettes/sh_general_release_dynamic_29.11.2022.fasta" # ITS fungal

# Define filter path for the forward and reverse reads of each sample fastq file
seq_path <- "./vignettes/MiSeq_SOP1"  # replace with the directory for fastq file after unzipping the folder

# Sort forward and reverse reads to be in the same order
fnFs <- sort(list.files(seq_path, pattern = "_R1_001.fastq"))
fnRs <- sort(list.files(seq_path, pattern = "_R2_001.fastq"))

# Specify the full path to the fnFs and fnRs
fnFs <- file.path(seq_path, fnFs)
fnRs <- file.path(seq_path, fnRs)

# Extract sample names, assuming file paths have format: FOL_DER/SAMPLENAME_XXX.fastq.gz OR FOL_DER/SAMPLENAME_XXX.fastq
sNames <- sapply(strsplit(fnFs, "_"), `[`, 2)
sNames <- sapply(strsplit(sNames, "/"), `[`, 2)

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

# Empty list for storing the results
TEST2 <- list()


# Loop the analysis process for every file in the directory
for (i in 1:length(fnFs)){

  # Assign metadata to be input as a list for ps.net argument
  meta_data <- list()
  for (j in 1:length(fnFs)){

    ID <- sapply(strsplit(fnFs[i], "_"), `[`, 2)
    ID <- sapply(strsplit(ID, "/"), `[`, 2)
    meta_data[[j]] <- as.list(meta[meta$sample.ID %in% ID,])
  }

  # Create a new row of results for each output from the `ps.net`
  TEST2[[i]] <- mga::mga(fastq.Fs = fnFs[i],
                             fastq.Rs = fnRs[i], # file paths for forward and reverse raw fastq files
                             filtFs = filtFs[i],
                             filtRs = filtRs[i], # file paths for filtered fastq files
                             refFasta = unite.its,
                             metadata = meta_data[[j]])

  # Save the results as they are produced, in case of crashes/errors
  save(TEST2, file = "vignettes/RDS/TEST2_mga_results.rda")
}

# Extract sample results from the saved list
RES <- list()
for (i in 1:length(TEST2)){

  RES[[i]] <- TEST2[[i]]$results.samples
}

# Build combined data frame of metrics
results.df6 <- do.call(rbind, RES)
results.df6 <- dplyr::relocate(results.df6, "sample.ID", .before = "Shannon") # move sample.ID column to leftmost

# Add additional metadata
meta2 <- read.delim(file = "vignettes/MiSeq_SOP/mouse.dpw.metadata")
results.df6$dpw <- meta2$dpw[1]
results.df6$dpwgroup <- NA
results.df6[results.df6$dpw <= 9,]$dpwgroup <- "0-9"


# Save data frame of results
saveRDS(results.df6, file = "vignettes/RDS/results.df_TEST2_mga.rds")


# Computing speed ----

# Define filter path for the forward and reverse reads of each sample fastq file
seq_path <- "./vignettes/MiSeq_SOP"  # replace with the directory for fastq file after unzipping the folder

# Sort forward and reverse reads to be in the same order
fnFs <- sort(list.files(seq_path, pattern = "_R1_001.fastq"))
fnRs <- sort(list.files(seq_path, pattern = "_R2_001.fastq"))

# Remove "Mock" fastq files
fnFs <- fnFs[-20]
fnRs <- fnRs[-20]

# Specify the full path to the fnFs and fnRs
fnFs <- file.path(seq_path, fnFs)
fnRs <- file.path(seq_path, fnRs)

# Extract sample names, assuming file paths have format: FOL_DER/SAMPLENAME_XXX.fastq.gz OR FOL_DER/SAMPLENAME_XXX.fastq
sNames <- sapply(strsplit(fnFs, "_"), `[`, 2)
sNames <- sapply(strsplit(sNames, "/"), `[`, 2)

# Define file names for filtered fastq.gz files and assign file paths
filt_path <- file.path(seq_path, "filtered") # Place filtered files in new subdirectory if one does not already exist
if (!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sNames, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sNames, "_R_filt.fastq.gz"))

## Per sample ----
# Empty list for storing the results
TICTOC <- list()

comp_time <- list()

# Loop the analysis process 10 times
for (i in 1:10) {

  # Assign metadata to be input as a list for ps.net argument
  meta_data <- list()
  for (j in 1:length(fnFs)) {

    ID <- sapply(strsplit(fnFs[i], "_"), `[`, 2)
    ID <- sapply(strsplit(ID, "/"), `[`, 2)
    meta_data[[j]] <- as.list(meta[meta$sample.ID %in% ID,])
  }

  # Add a sample to every new i loop
  Fs <- fnFs[1:i]
  Rs <- fnRs[1:i]
  FFs <- filtFs[1:i]
  FRs <- filtRs[1:i]

  tictoc::tic()

  # Loop for every file in Fs
  for (n in 1:length(Fs)) {

  # Create a new row of results for each output from the `ps.net`
  TICTOC[[n]] <- mga::mga(fastq.Fs = fnFs[n],
                              fastq.Rs = fnRs[n], # file paths for forward and reverse raw fastq files
                              filtFs = filtFs[n],
                              filtRs = filtRs[n], # file paths for filtered fastq files
                              refFasta = silva.138.1S, # silva ref FASTA with species
                              metadata = meta_data[[j]])

  # Save the results as they are produced, in case of crashes/errors
  save(TICTOC, file = "vignettes/RDS/TICTOC_mga_results.rda")
  }

  comp_time[[i]] <- tictoc::toc()

  saveRDS(comp_time, file = "vignettes/RDS/comp_time.rds")
}

comp_time <- readRDS(file = "vignettes/RDS/comp_time.rds")

time.df <- list()

for (i in 1:length(comp_time)) {

  time.df[[i]] <- c(paste(i), comp_time[[i]]$callback_msg)

}

time.df2 <- t(as.data.frame(time.df))
time.df2 <- as.data.frame(time.df2)
rownames(time.df2) <- NULL
colnames(time.df2) <- c("Samples", "Run_time")
time.df2$Run_time <- gsub(" sec elapsed", "", time.df2$Run_time)
time.df2$Samples <- as.factor(time.df2$Samples)
time.df2$Run_time <- as.numeric(time.df2$Run_time)

plot(x = time.df2$Samples, y = time.df2$Run_time)

time_plot <- ggplot(data = time.df2, aes(x = Samples, y = Run_time)) +
  geom_point(shape = 19, size = 3) +
  labs(x = "No. Samples", y = "Run Time (s)") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 15, family = "sans", face = "bold"),
        axis.text.y = element_text(size = 12, family = "sans"),
        axis.text.x = element_text(size = 12, family = "sans"),
        axis.title.x = element_text(size = 15, family = "sans", face = "bold"),
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = 0.2, fill = NA),
        plot.background = element_blank(),
        plot.margin = unit(c(0.2,0.1,0,0.2), "cm")) # top, right, bottom, left)

ggsave(file = "vignettes/figures/comp_time.png", plot = time_plot, width = 6.86, height = 6.86)


# No metadata ----

# Empty list for storing the results
NO.META <- list()

# Loop the analysis process for every file in the directory
for (i in 1:3){

  # Create a new row of results for each output from the `ps.net`
  NO.META[[i]] <- mga::mga(fastq.Fs = fnFs[i],
                              fastq.Rs = fnRs[i], # file paths for forward and reverse raw fastq files
                              filtFs = filtFs[i],
                              filtRs = filtRs[i], # file paths for filtered fastq files
                              refFasta = silva.138.1S) # silva ref FASTA with species

  # Save the results as they are produced, in case of crashes/errors
  save(NO.META, file = "vignettes/RDS/no.meta_mga_results.rda")
}

# Wrong metadata format ----

pathogens <- read.csv(file = "vignettes/Plant_Pathogens.csv")

# Empty list for storing the results
WRONG.META <- list()

# Loop the analysis process for every file in the directory
for (i in 1:3){

  # Assign metadata to be input as a list for ps.net argument
  # meta_data <- list()
  # for (j in 1:length(fnFs)){
  #
  #   ID <- sapply(strsplit(fnFs[i], "_"), `[`, 2)
  #   ID <- sapply(strsplit(ID, "/"), `[`, 2)
  #   meta_data[[j]] <- as.list(pathogens[pathogens$sample.ID %in% ID,])
  # }

  # Create a new row of results for each output from the `ps.net`
  WRONG.META[[i]] <- mga::mga(fastq.Fs = fnFs[i],
                           fastq.Rs = fnRs[i], # file paths for forward and reverse raw fastq files
                           filtFs = filtFs[i],
                           filtRs = filtRs[i], # file paths for filtered fastq files
                           refFasta = silva.138.1S, # silva ref FASTA with species
                           metadata = pathogens[i])

  # Save the results as they are produced, in case of crashes/errors
  save(WRONG.META, file = "vignettes/RDS/wrong.meta_mga_results.rda")
}

