# ---------------------------------------------------------------------------------------------------------------#
# Integrated FST Pipeline - J.R. Phaneuf 2025
#
# Full pipeline for testing four sequencing depth normalization methods (CSS with metagenomeSeq, TMM with
# edgeR, and MED/VST with DESeq2), calculating Bray-Curtis distances, investigating beta diversity significance
# with ANOSIM, NMDS visualization, running FEAST, calculating correlation coefficients between bioinformatic
# tools and dPCR markers, completing differential abundance analyses with ANCOM-BC2, and creating figures.
#
# For Step 1, only Bracken outputs and a grouping metadata file are required to begin.
# For Step 2, tool metadata and "Source" OTUs in the read table are required.
# For Step 3, data file must contain all tool outputs and dPCR marker concentrations.
# For Step 4, a metadata file continaing state, site, and land use is required.
# For Step 5, finalized NMDS plots and heatmaps must be saved and called from a directory.
# For Step 6, SourceApp results must be matched to land use for each sample for plotting.
# ---------------------------------------------------------------------------------------------------------------#

# ---------------------------------------------------------------------------------------------------------------#
############## Step 0: Load Required Packages and Format Data ##############
# ---------------------------------------------------------------------------------------------------------------#
library(metagenomeSeq)
library(edgeR)
library(DESeq2)
library(vegan)
library(dplyr)
library(tidyr)
library(tidyverse)
library(tools)
library(FEAST)
library(ANCOMBC)
library(shadowtext)
library(pheatmap)
library(ggsci)
library(cowplot)
library(patchwork)
library(magick)
library(ellipse)

# Load and merge Bracken files from a directory into a wide count matrix.
# taxa_col: Name of the column containing taxa in the bracken output
# count_col: column carrying read counts 
# suffix_pattern: Strip from the filename to get a clean sample ID
load_bracken_files <- function(dir, file_pattern, taxa_col, count_col, suffix_pattern) {
  files <- list.files(path = dir, pattern = file_pattern, full.names = TRUE)
  parsed <- lapply(files, function(f) {
    df       <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    samp_id  <- gsub(suffix_pattern, "", basename(f))
    df %>%
      dplyr::select(all_of(c(taxa_col, count_col))) %>%
      dplyr::rename(taxon = all_of(taxa_col), !!samp_id := all_of(count_col))
  })
  merged <- Reduce(function(x, y) full_join(x, y, by = "taxon"), parsed)
  merged[is.na(merged)] <- 0
  merged %>% column_to_rownames("taxon") %>% as.matrix()
}

# Compute 95% confidence ellipses for NMDS coordinates grouped by a factor.
# The "coords_df" dataframe must contain columns named denoting NMDS coordinates
compute_ellipses <- function(coords_df, group_col, x_col = "NMDS1", y_col = "NMDS2") {
  coords_df %>%
    group_by(.data[[group_col]]) %>%
    filter(n() >= 3) %>%
    group_map(~ {
      pts <- .x[, c(x_col, y_col)]
      ell <- as.data.frame(ellipse(cov(pts), centre = colMeans(pts),
                                   level = 0.95, npoints = 100))
      colnames(ell) <- c("x", "y")
      ell[[group_col]] <- .y[[group_col]]
      ell
    }) %>%
    bind_rows()
}

# Shared ggplot2 theme for the NMDS plots.
theme_nmds <- function(base_size = 18, legend_position = "right",
                       title_size = 16, legend_text_size = 16,
                       axis_text_size = 12, axis_title_size = 14) {
  theme(
    axis.text.y      = element_text(colour = "black", size = axis_text_size,   face = "bold"),
    axis.text.x      = element_text(colour = "black", size = axis_text_size,   face = "bold"),
    axis.title.y     = element_text(face = "bold",    size = axis_title_size),
    axis.title.x     = element_text(face = "bold",    size = axis_title_size,  colour = "black"),
    legend.text      = element_text(size = legend_text_size, face = "bold", colour = "black"),
    legend.title     = element_text(size = legend_text_size, colour = "black", face = "bold"),
    legend.position  = legend_position,
    legend.key       = element_blank(),
    panel.background = element_blank(),
    panel.border     = element_rect(colour = "black", fill = NA, linewidth = 1.2),
    plot.title       = element_text(size = title_size, face = "bold", hjust = 0.5)
  )
}

# Save a ggplot to TIFF with project-standard settings.
save_tiff <- function(plot, filename, width = 8.5, height = 6, dpi = 1200) {
  ggsave(filename = filename, plot = plot, device = "tiff",
         dpi = dpi, width = width, height = height, units = "in", compression = "lzw")
}

# ---------------------------------------------------------------------------------------------------------------#
############## Step 1: Testing different normalization methods: CSS, TMM, MED, and VST ##############
# ---------------------------------------------------------------------------------------------------------------#

# --- 1a: Load metadata and merge bracken files ---
bracken_dir <- "br_genera/"
file_list   <- list.files(path = bracken_dir, pattern = "_genus\\.br$", full.names = TRUE)

read_bracken <- function(file) {
  df          <- read.delim(file)
  sample_name <- file_path_sans_ext(basename(file)) %>% str_replace(".genus$", "")
  df %>%
    filter(taxonomy_lvl == "G") %>%
    dplyr::select(name, new_est_reads) %>%
    dplyr::rename(!!sample_name := new_est_reads)
}

merged       <- reduce(lapply(file_list, read_bracken), full_join, by = "name")
merged[is.na(merged)] <- 0
count_matrix <- merged %>% column_to_rownames("name") %>% as.matrix()
otu_table    <- count_matrix
colnames(otu_table) <- sub("_genus", "", colnames(otu_table))

grouping_file <- "norm_metadata.csv"
groupings     <- read.csv(grouping_file, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)

# --- 1b: Match and order samples ---
sample_num   <- ncol(otu_table)
grouping_num <- nrow(groupings)

if (sample_num != grouping_num) {
  message("Number of samples differs between OTU table and metadata. Attempting intersection.")
}

rows_to_keep <- intersect(colnames(otu_table), rownames(groupings))
if (length(rows_to_keep) == 0) stop("No matching sample IDs found.")

otu_table <- otu_table[, rows_to_keep, drop = FALSE]
groupings <- groupings[rows_to_keep, , drop = FALSE]

if (!identical(colnames(otu_table), rownames(groupings))) {
  stop("Failed to match sample ordering after intersection.")
}
message("Samples successfully matched: ", length(rows_to_keep))

group_col    <- "State"
group_col_ga <- "GASites"
group_col_or <- "ORSites"

# --- 1c: Pairwise ANOSIM function ---
pairwise_anosim <- function(dist_matrix, groups, permutations = 999) {
  groups <- as.factor(groups)
  combs  <- combn(levels(groups), 2)
  results <- data.frame(Group1 = character(), Group2 = character(),
                        R = numeric(), p = numeric(), stringsAsFactors = FALSE)
  for (i in seq_len(ncol(combs))) {
    g1 <- combs[1, i]
    g2 <- combs[2, i]
    subset_idx <- which(groups %in% c(g1, g2))
    if (length(subset_idx) < 3) {
      res_p <- NA; res_r <- NA
    } else {
      sub_dist <- as.dist(as.matrix(dist_matrix)[subset_idx, subset_idx])
      anos     <- anosim(sub_dist, groups[subset_idx], permutations = permutations)
      res_p    <- anos$signif
      res_r    <- anos$statistic
    }
    results <- rbind(results, data.frame(
      Group1 = g1, Group2 = g2,
      R = ifelse(is.na(res_r), NA, round(res_r, 3)),
      p = res_p
    ))
  }
  results$FDR <- p.adjust(results$p, method = "fdr")
  results$p   <- signif(results$p,   3)
  results$FDR <- signif(results$FDR, 3)
  return(results)
}

# --- 1d: Normalization, Bray-Curtis, and ANOSIM analysis functions ---
run_analysis <- function(otu_table, groupings,
                         group_col    = "State",
                         group_col_ga = "GASites",
                         group_col_or = "ORSites",
                         method = c("CSS", "TMM", "MED", "VST")) {

  method <- match.arg(method)
  cat("Running method:", method, "\n")

  if (method == "CSS") {
    pheno        <- AnnotatedDataFrame(groupings)
    feature_data <- AnnotatedDataFrame(data.frame(OTU = rownames(otu_table), stringsAsFactors = FALSE))
    rownames(feature_data) <- rownames(otu_table)
    obj          <- newMRexperiment(counts = otu_table, phenoData = pheno, featureData = feature_data)
    p            <- cumNormStat(obj, pFlag = FALSE)
    obj_norm     <- cumNorm(obj, p = p)
    mat          <- log2(MRcounts(obj_norm, norm = TRUE) + 1)

  } else if (method == "TMM") {
    dge                  <- DGEList(counts = otu_table)
    dge$samples$group    <- groupings[[group_col]]
    dge                  <- calcNormFactors(dge, method = "TMM")
    mat                  <- cpm(dge, normalized.lib.sizes = TRUE)

  } else if (method == "MED") {
    dds <- DESeqDataSetFromMatrix(countData = otu_table,
                                  colData   = data.frame(row.names = rownames(groupings), groupings),
                                  design    = ~1)
    dds <- estimateSizeFactors(dds, type = "ratio")
    mat <- counts(dds, normalized = TRUE)

  } else if (method == "VST") {
    dds <- DESeqDataSetFromMatrix(countData = otu_table,
                                  colData   = data.frame(row.names = rownames(groupings), groupings),
                                  design    = ~1)
    vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
    mat <- assay(vsd)
  }

  write.csv(mat, paste0(method, "_otus.csv"), quote = FALSE)

  subsets <- list(
    "All" = rep(TRUE, ncol(mat)),
    "GA"  = !is.na(groupings[[group_col_ga]]) & groupings[[group_col_ga]] != "",
    "OR"  = !is.na(groupings[[group_col_or]]) & groupings[[group_col_or]] != ""
  )

  results_list <- list()
  for (subset_name in names(subsets)) {
    idx     <- subsets[[subset_name]]
    mat_sub <- mat[, idx, drop = FALSE]
    if (ncol(mat_sub) == 0) {
      cat("Skipping subset", subset_name, "- no samples.\n"); next
    }
    dist_mat <- vegdist(t(mat_sub), method = "bray")

    group_vec <- switch(subset_name,
      "All" = factor(groupings[[group_col]][idx]),
      "GA"  = factor(groupings[[group_col_ga]][idx]),
      "OR"  = factor(groupings[[group_col_or]][idx])
    )
    site_vec <- factor(groupings$Site[idx])

    anosim_res   <- anosim(dist_mat, grouping = group_vec, permutations = 999)
    pairwise_res <- pairwise_anosim(dist_mat, group_vec)
    nmds_res     <- metaMDS(dist_mat, k = 2, trymax = 100)
    coords       <- scores(nmds_res, display = "sites")

    results_list[[subset_name]] <- list(
      dist           = dist_mat,
      anosim         = anosim_res,
      pairwise_anosim = pairwise_res,
      coords         = coords,
      stress         = nmds_res$stress,
      group_vector   = group_vec,
      site_vector    = site_vec
    )
    cat("Completed", method, "->", subset_name, "\n")
  }
  return(results_list)
}

# --- 1e: Run all four normalization methods ---
set.seed(123)
css_results <- run_analysis(otu_table, groupings, method = "CSS")
tmm_results <- run_analysis(otu_table, groupings, method = "TMM")
med_results <- run_analysis(otu_table, groupings, method = "MED")
vst_results <- run_analysis(otu_table, groupings, method = "VST")

# --- 1f: Collect and save NMDS stress values ---
collect_nmds_stress <- function(results_list, method_name) {
  stress_df <- data.frame(Method = character(), Subset = character(),
                          Stress = numeric(), stringsAsFactors = FALSE)
  for (subset_name in names(results_list)) {
    nmds <- metaMDS(results_list[[subset_name]]$dist, k = 2, trymax = 100, trace = FALSE)
    stress_df <- rbind(stress_df, data.frame(
      Method = method_name, Subset = subset_name, Stress = nmds$stress
    ))
  }
  return(stress_df)
}

stress_all <- rbind(
  collect_nmds_stress(css_results, "CSS"),
  collect_nmds_stress(tmm_results, "TMM"),
  collect_nmds_stress(med_results, "MED"),
  collect_nmds_stress(vst_results, "VST")
)
write.csv(stress_all, "nmds_stress_values.csv", row.names = FALSE)
# The MED normalization method produced the lowest stress values, used here onward. 

# --- 1g: Collect and save MED ANOSIM results (UPDATE BASED ON CHOSEN NORMALIZATION METHOD) ---
collect_anosim_results <- function(results_list, method_name = "MED") {
  all_res <- data.frame()
  for (subset_name in names(results_list)) {
    anosim_res <- results_list[[subset_name]]$anosim
    overall    <- data.frame(Method = method_name, Subset = subset_name,
                             Comparison = "Overall", R = anosim_res$statistic,
                             p = anosim_res$signif, FDR = NA)
    pairwise   <- results_list[[subset_name]]$pairwise_anosim
    if (!is.null(pairwise) && nrow(pairwise) > 0) {
      pairwise_df <- data.frame(
        Method     = method_name, Subset = subset_name,
        Comparison = paste(pairwise$Group1, "vs", pairwise$Group2),
        R = pairwise$R, p = pairwise$p, FDR = pairwise$FDR
      )
      all_res <- rbind(all_res, rbind(overall, pairwise_df))
    } else {
      all_res <- rbind(all_res, overall)
    }
  }
  return(all_res)
}

med_anosim <- collect_anosim_results(med_results)
write.csv(med_anosim, "anosim_med_results.csv", row.names = FALSE)

# --- 1h: NMDS visualizations ---

# Color / shape palettes
all_shapes <- c(1,0,2,3,6,7,9,10,11,5,12,13,16,14,8,4)
# Bea, Cap, Cha, Chi, Cow, Dai, Daw, Fan, Gal, Joh, McK, Roc, Sco, Tua, Uto, Wil
all_colors <- c("Georgia" = "#D62728", "Oregon" = "#1F77B4")

# MED: Georgia vs Oregon
df_nm_all      <- cbind(as.data.frame(med_results$All$coords),
                         group_vector = med_results$All$group_vector,
                         Site         = med_results$All$site_vector)
ellipse_df_all <- compute_ellipses(df_nm_all, "group_vector")

nmds_state <- ggplot(df_nm_all, aes(x = NMDS1, y = NMDS2, color = group_vector)) +
  geom_point(aes(shape = Site), size = 2) +
  geom_path(data = ellipse_df_all, aes(x = x, y = y, color = group_vector),
            linetype = "dashed", linewidth = 1) +
  scale_color_manual(values = all_colors) +
  scale_shape_manual(values = all_shapes) +
  xlim(-1.3, 1.3) + ylim(-1.3, 1.3) +
  theme_light(base_size = 18) +
  ggtitle(paste0("a) Georgia vs Oregon NMDS (Stress = ", round(med_results$All$stress, 3), ")")) +
  labs(x = "NMDS1", y = "NMDS2", colour = "State", shape = "Site") +
  theme_nmds(legend_text_size = 16)
nmds_state
save_tiff(nmds_state, "tiff_plots/nmds_state.tiff")

# MED: Georgia sites
df_nm_ga      <- cbind(as.data.frame(med_results$GA$coords),
                        group_vector = med_results$GA$group_vector,
                        Site         = med_results$GA$site_vector)
ellipse_df_ga <- compute_ellipses(df_nm_ga, "group_vector")

nmds_ga <- ggplot(df_nm_ga, aes(x = NMDS1, y = NMDS2, color = group_vector)) +
  geom_point(size = 2) +
  geom_path(data = ellipse_df_ga, aes(x = x, y = y, color = group_vector),
            linetype = "dashed", linewidth = 1) +
  scale_color_d3(palette = "category20") +
  xlim(-1, 1) + ylim(-1, 1) +
  theme_light(base_size = 18) +
  ggtitle(paste0("Georgia NMDS (Stress = ", round(med_results$GA$stress, 4), ")")) +
  labs(x = "NMDS1", y = "NMDS2", colour = "Site") +
  theme_nmds(legend_text_size = 12, title_size = 12)
nmds_ga

# MED: Oregon sites
df_nm_or <- cbind(as.data.frame(med_results$OR$coords),
                   group_vector = med_results$OR$group_vector,
                   Site         = med_results$OR$site_vector)

nmds_or <- ggplot(df_nm_or, aes(x = NMDS1, y = NMDS2, color = group_vector)) +
  geom_point(size = 3) +
  scale_color_d3(palette = "category20") +
  xlim(-0.75, 0.75) + ylim(-0.75, 0.75) +
  theme_light(base_size = 18) +
  ggtitle(paste0("Oregon NMDS (Stress = ", round(med_results$OR$stress, 4), ")")) +
  labs(x = "NMDS1", y = "NMDS2", colour = "Site") +
  theme_nmds(legend_text_size = 12, title_size = 12)
nmds_or

# --- 1i: Land use significance on diversity ---
med_mat <- read.csv("MED_otus.csv", row.names = 1, check.names = FALSE)

analyze_landuse <- function(med_mat, groupings, state_col, buffer_col, permutations = 999) {
  valid    <- !is.na(groupings[[buffer_col]]) & groupings[[buffer_col]] != "" &
              !is.na(groupings[[state_col]])  & groupings[[state_col]]  != ""
  dist_mat <- vegdist(t(med_mat[, valid]), method = "bray")
  nmds_res <- metaMDS(dist_mat, k = 2, trymax = 100, trace = FALSE)
  coords   <- as.data.frame(nmds_res$points)
  coords$Sample <- rownames(coords)
  coords$Group  <- factor(groupings[[buffer_col]][valid],
                           levels = unique(groupings[[buffer_col]][valid]))
  coords$Site   <- factor(groupings[[state_col]][valid],
                           levels = unique(groupings[[state_col]][valid]))
  group_vec    <- factor(groupings[[buffer_col]][valid])
  anosim_res   <- anosim(dist_mat, grouping = group_vec, permutations = permutations)
  pairwise_res <- pairwise_anosim(dist_mat, group_vec, permutations = permutations)
  list(dist     = dist_mat, nmds   = nmds_res, coords = coords,
       stress   = round(nmds_res$stress, 3),
       anosim   = anosim_res, pairwise = pairwise_res, valid = valid)
}

lu_results <- list(
  GA_lu = analyze_landuse(med_mat, groupings, "GASites",  "QuartMile"),
  OR_lu = analyze_landuse(med_mat, groupings, "ORSites",  "QuartMile")
)

collect_lu_anosim <- function(results_list) {
  all_res <- data.frame()
  for (subset_name in names(results_list)) {
    anosim_res <- results_list[[subset_name]]$anosim
    overall    <- data.frame(Subset = subset_name, Comparison = "Overall",
                             R = anosim_res$statistic, p = anosim_res$signif, FDR = NA)
    pairwise   <- results_list[[subset_name]]$pairwise
    if (!is.null(pairwise) && nrow(pairwise) > 0) {
      pairwise_df <- data.frame(
        Subset     = subset_name,
        Comparison = paste(pairwise$Group1, "vs", pairwise$Group2),
        R = pairwise$R, p = pairwise$p, FDR = pairwise$FDR
      )
      all_res <- rbind(all_res, rbind(overall, pairwise_df))
    } else {
      all_res <- rbind(all_res, overall)
    }
  }
  return(all_res)
}

landuse_anosim <- collect_lu_anosim(lu_results)
write.csv(landuse_anosim, "anosim_landuse_results.csv", row.names = FALSE)

# --- 1j: NMDS land use visualizations ---
lu_colors <- c("Agriculture" = "#FF7F0E", "Industrial" = "#D62728",
                "Forest/Parks" = "#2CA02C", "Residential" = "#1F77B4")

# Georgia land use
df_nm_luga      <- lu_results$GA_lu$coords
ellipse_df_luga <- compute_ellipses(df_nm_luga, "Group", x_col = "MDS1", y_col = "MDS2")

nmds_luga <- ggplot(df_nm_luga, aes(x = MDS1, y = MDS2, color = Group)) +
  geom_point(aes(shape = Site), size = 2.2) +
  geom_path(data = ellipse_df_luga, aes(x = x, y = y, color = Group),
            linetype = "dashed", linewidth = 1) +
  scale_color_manual(values = lu_colors) +
  scale_shape_manual(values = c(0,2,6,5,8,4)) +
  # Cap, Cha, Cow, Joh, Uto, Wil
  xlim(-1.1, 1.1) + ylim(-1.1, 1.1) +
  theme_light(base_size = 18) +
  ggtitle(paste0("Georgia Land Use NMDS (Stress = ", lu_results$GA_lu$stress, ")")) +
  labs(x = "NMDS1", y = "NMDS2", colour = "Land Use", shape = "Site") +
  theme_nmds(legend_text_size = 12, title_size = 12)
nmds_luga
save_tiff(nmds_luga, "tiff_plots/nmds_luga.tiff")

# Oregon land use
df_nm_luor      <- lu_results$OR_lu$coords
ellipse_df_luor <- compute_ellipses(df_nm_luor, "Group", x_col = "MDS1", y_col = "MDS2")

nmds_luor <- ggplot(df_nm_luor, aes(x = MDS1, y = MDS2, color = Group)) +
  geom_point(aes(shape = Site), size = 2.6) +
  geom_path(data = ellipse_df_luor, aes(x = x, y = y, color = Group),
            linetype = "dashed", linewidth = 1) +
  scale_color_manual(values = lu_colors) +
  scale_shape_manual(values = c(1,3,7,9,10,11,12,13,16,14)) +
  # Bea, Chi, Dai, Daw, Fan, Gal, McK, Roc, Sco, Tua
  xlim(-0.75, 0.75) + ylim(-0.75, 0.75) +
  theme_light(base_size = 18) +
  ggtitle(paste0("Oregon Land Use NMDS (Stress = ", lu_results$OR_lu$stress, ")")) +
  labs(x = "NMDS1", y = "NMDS2", colour = "Land Use", shape = "Site") +
  theme_nmds(legend_text_size = 12, title_size = 12)
nmds_luor
save_tiff(nmds_luor, "tiff_plots/nmds_luor.tiff")

# ---------------------------------------------------------------------------------------------------------------#
############## Step 2: Running FEAST and Calculating Relative Source Contributions ##############
# ---------------------------------------------------------------------------------------------------------------#

# --- 2a: Build OTU table from bracken source files ---
bracken_directory <- "br_sources/"
feast_files       <- list.files(path = bracken_directory, pattern = "*_genus.br", full.names = TRUE)

feast_bracken <- function(file) {
  sample_id <- tools::file_path_sans_ext(basename(file))
  read.delim(file, sep = "\t", header = TRUE) %>%
    dplyr::select(name, fraction_total_reads) %>%
    rename(Genus = name, !!sample_id := fraction_total_reads)
}

otu_table <- reduce(lapply(feast_files, feast_bracken), full_join, by = "Genus") %>%
  replace_na(as.list(setNames(rep(0, length(.)), names(.)))) %>%
  column_to_rownames("Genus")

write.table(otu_table, file = "usgs_otus.tsv", sep = "\t", quote = FALSE, col.names = NA)

# Remove "_genus" from columns and rename rows to "OTU_x" before proceeding to Step 2c.

# --- 2b: Load metadata and OTU table ---
feast_metadata <- Load_metadata(metadata_path = "feast_metadata_usgs.txt")
otus_dec       <- Load_CountMatrix(CountMatrix_path = "usgs_otus.tsv")
otus           <- ceiling(otus_dec * 1000000)

write.table(t(otus), file = "feast_results/usgs_otus_mst2.tsv", # Used for mST2 input in HPC environment
            sep = "\t", quote = FALSE, col.names = NA)

# --- 2c: Run FEAST (x3) ---
FEAST_output <- FEAST(C = otus, metadata = feast_metadata, different_sources_flag = 0,
                      dir_path = "feast_results/",
                      outfile  = paste0("usgs_FEAST", run_i))

# ---------------------------------------------------------------------------------------------------------------#
############## Step 3: Spearman's Rank Correlation for Significance Testing ##############
# ---------------------------------------------------------------------------------------------------------------#

# --- 3a: Load and pivot apportionments ---
markers <- read.csv("markers_apportions.csv", stringsAsFactors = FALSE)
markers_clean <- markers %>%
  mutate(across(-Sample, ~ as.numeric(str_remove(., "%"))))

apportions_long <- markers_clean %>%
  pivot_longer(cols = -Sample, names_to = "SourceTool", values_to = "Apportion") %>%
  mutate(
    Tool   = ifelse(str_ends(SourceTool, "_f"), "FEAST", "mST2"),
    Source = sub("_[fm]$", "", SourceTool)
  ) %>%
  dplyr::select(-SourceTool)

# --- 3b: Define sources per marker ---
marker_sources <- list(
  EC23S  = c("cat", "chick", "cow", "dog", "goat", "pig", "sep", "ww"),
  HF183  = c("sep", "ww"),
  Gull4  = c("chick"),
  CowM2  = c("cow"),
  BacCan = c("dog")
)

# --- 3c: Summed apportionments per Sample × Marker × Tool ---
apportions_sum <- apportions_long %>%
  crossing(Marker = names(marker_sources)) %>%
  rowwise() %>%
  mutate(keep = Source %in% marker_sources[[Marker]]) %>%
  filter(keep) %>%
  group_by(Sample, Marker, Tool) %>%
  summarise(Apportion = sum(Apportion, na.rm = TRUE), .groups = "drop")

# --- 3d: Join marker concentrations ---
markers_long <- markers_clean %>%
  pivot_longer(cols = c("EC23S", "HF183", "Gull4", "CowM2", "BacCan"),
               names_to = "Marker", values_to = "dPCR_conc")

ma_long <- markers_long %>%
  left_join(apportions_sum, by = c("Sample", "Marker"))

ma_long_extended <- ma_long %>%
  bind_rows(
    ma_long %>%
      filter(Marker == "HF183") %>%
      mutate(Tool = "SourceApp", Apportion = SourceApp)
  )

# --- 3e: Compute Spearman correlations with FDR correction ---
safe_cor <- function(x, y) {
  tryCatch(cor.test(x, y, method = "spearman"),
           error = function(e) list(estimate = NA, p.value = NA))
}

ma_results <- ma_long_extended %>%
  group_by(Marker, Tool) %>%
  summarise(cor_test = list(safe_cor(dPCR_conc, Apportion)), .groups = "drop") %>%
  mutate(
    rho     = sapply(cor_test, function(x) if (is.list(x)) x$estimate  else NA),
    p_value = sapply(cor_test, function(x) if (is.list(x)) x$p.value   else NA)
  ) %>%
  group_by(Marker) %>%
  mutate(p_adj_FDR = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  dplyr::select(Marker, Tool, rho, p_value, p_adj_FDR) %>%
  bind_rows(
    markers_clean %>%
      dplyr::select(Sample, imp, EC23S, HF183, Gull4, CowM2, BacCan, SourceApp) %>%
      pivot_longer(cols = c("EC23S", "HF183", "Gull4", "CowM2", "BacCan", "SourceApp"),
                   names_to = "Marker", values_to = "Value") %>%
      group_by(Marker) %>%
      summarise(cor_test = list(tryCatch(
        cor.test(imp, Value, method = "spearman", exact = FALSE),
        error = function(e) list(estimate = NA, p.value = NA)
      )), .groups = "drop") %>%
      mutate(
        Tool      = "%Imp",
        rho       = sapply(cor_test, function(x) if (is.list(x)) x$estimate else NA),
        p_value   = sapply(cor_test, function(x) if (is.list(x)) x$p.value  else NA),
        p_adj_FDR = p.adjust(p_value, method = "BH")
      ) %>%
      dplyr::select(Marker, Tool, rho, p_value, p_adj_FDR)
  ) %>%
  bind_rows(
    combn(c("EC23S", "HF183", "Gull4", "CowM2", "BacCan"), 2, simplify = FALSE) %>%
      map_dfr(function(pair) {
        tryCatch({
          test <- cor.test(markers_clean[[pair[1]]], markers_clean[[pair[2]]],
                           method = "spearman", exact = FALSE)
          tibble(Marker = pair[1], Tool = pair[2], rho = test$estimate, p_value = test$p.value)
        }, error = function(e) {
          tibble(Marker = pair[1], Tool = pair[2], rho = NA, p_value = NA)
        })
      }) %>%
      mutate(p_adj_FDR = p.adjust(p_value, method = "BH"))
  )

print(ma_results)
write.csv(ma_results, "correlation_results.csv", row.names = FALSE)

# --- 3f: Plot LOESS smooth ribbons ---
source_map <- c("Cat" = "cat", "Chicken" = "chick", "Cow" = "cow", "Dog" = "dog",
                "Goat" = "goat", "Pig" = "pig", "Septage" = "sep", "Wastewater" = "ww")

ma_long_tidy <- ma_long %>%
  dplyr::select(-Apportion) %>%
  tidyr::pivot_longer(cols = matches("(_f|_m)$"), names_to = "Source_Tool", values_to = "Apportion") %>%
  mutate(
    Tool     = ifelse(grepl("_f$", Source_Tool), "FEAST", "mST2"),
    Source   = gsub("(_f|_m)$", "", Source_Tool),
    Source   = dplyr::recode(Source, !!!source_map),
    Apportion = as.numeric(str_remove_all(Apportion, "%"))
  ) %>%
  distinct(Sample, Marker, Source, Tool, .keep_all = TRUE)

custom_colors  <- c("cat" = "#8C564B", "chick" = "#9467BD", "cow" = "#2CA02C",
                    "dog" = "#E377C2", "goat" = "#1F77B4", "pig" = "#FF7F0E",
                    "sep" = "#17BECF", "ww"   = "#D62728")
custom_linetype <- c("FEAST" = "solid", "mST2" = "dashed")
source_labels   <- c("cat" = "Cat", "chick" = "Chicken", "cow" = "Cow", "dog" = "Dog",
                     "goat" = "Goat", "pig" = "Pig", "sep" = "Septage", "ww" = "Wastewater")

axis_limits <- list(
  EC23S  = list(x = c(0, 96.8), y = c(0, 50),  breaks = NULL, label_accuracy = 1),
  HF183  = list(x = c(0, 29.0), y = c(0, 100), breaks = NULL, label_accuracy = 1),
  BacCan = list(x = c(0, 75.4), y = c(0, 0.5), breaks = NULL, label_accuracy = 1),
  Gull4  = list(x = c(0, 1),    y = c(0, 6),   breaks = c(0, 0.25, 0.5, 0.75, 1), label_accuracy = 0.01),
  CowM2  = list(x = c(0, 4),    y = c(0, 50),  breaks = NULL, label_accuracy = 1)
)

markers_to_plot <- c("EC23S", "HF183", "CowM2", "BacCan", "Gull4")
marker_labels   <- c(EC23S = "a) ", HF183 = "b) ", BacCan = "c) ", Gull4 = "d) ", CowM2 = "e) ")

make_marker_plot <- function(mk) {
  marker_df <- ma_long_tidy %>%
    filter(Marker == mk, Source %in% marker_sources[[mk]], !is.na(Apportion))
  has_conc <- sum(!is.na(marker_df$dPCR_conc))

  if (has_conc > 0) {
    marker_df %>%
      filter(!is.na(dPCR_conc)) %>%
      ggplot(aes(x = dPCR_conc, y = Apportion, color = Source, fill = Source, linetype = Tool)) +
      geom_point(aes(shape = Tool), size = 1.5, alpha = 0.3, stroke = 0.3, color = "black") +
      geom_smooth(method = "loess", se = TRUE, alpha = 0.15, linewidth = 0.8, fullrange = TRUE) +
      scale_color_manual(values = custom_colors, labels = source_labels) +
      scale_fill_manual(values  = custom_colors, labels = source_labels) +
      scale_linetype_manual(values = custom_linetype, labels = c("FEAST", "mST2")) +
      scale_shape_manual(values = c("FEAST" = 22, "mST2" = 21)) +
      scale_x_continuous(
        limits   = axis_limits[[mk]]$x,
        breaks   = if (!is.null(axis_limits[[mk]]$breaks)) axis_limits[[mk]]$breaks else waiver(),
        expand   = c(0, 0),
        labels   = scales::label_number(accuracy = axis_limits[[mk]]$label_accuracy)
      ) +
      scale_y_continuous(limits = axis_limits[[mk]]$y, expand = c(0, 0), name = "RSC (%)") +
      labs(title = paste0(marker_labels[[mk]], mk), x = NULL) +
      theme_minimal(base_size = 14) +
      theme(legend.position = "none")
  } else {
    ggplot(marker_df, aes(x = Tool, y = Apportion, fill = Source)) +
      geom_boxplot(alpha = 0.2, outlier.shape = NA, color = "grey40") +
      geom_jitter(aes(shape = Tool), width = 0.2, height = 0, alpha = 0.6,
                  size = 2, stroke = 0.5, color = "black") +
      scale_fill_manual(values  = custom_colors, labels = source_labels) +
      scale_color_manual(values = custom_colors, labels = source_labels) +
      scale_shape_manual(values = c("FEAST" = 22, "mST2" = 21)) +
      scale_y_continuous(limits = axis_limits[[mk]]$y, expand = c(0, 0), name = "RSC (%)") +
      labs(title = if (mk == "CowM2") "e) CowM2 (All Samples BLoD)" else paste0(marker_labels[[mk]], mk),
           x = "Tool") +
      theme_minimal(base_size = 14) +
      theme(legend.position = "none")
  }
}

smooth_plots <- lapply(markers_to_plot, make_marker_plot)

smooth_plots_labeled <- lapply(seq_along(smooth_plots), function(i) {
  p <- smooth_plots[[i]]
  if (i %in% c(4, 5)) p <- p + labs(x = "Concentration (cp/mL)")
  if (i %in% c(2, 3, 5)) p <- p + labs(y = NULL)
  p
})

# Build shared legend from EC23S
make_legend_plot <- function(mk) {
  marker_df <- ma_long_tidy %>%
    filter(Marker == mk, Source %in% marker_sources[[mk]],
           !is.na(dPCR_conc), !is.na(Apportion))
  ggplot(marker_df, aes(x = dPCR_conc, y = Apportion,
                        color = Source, fill = Source, linetype = Tool)) +
    geom_point(aes(shape = Tool), size = 1.5, alpha = 0.3, stroke = 0.3, color = "black") +
    geom_smooth(method = "loess", se = TRUE, alpha = 0.15, linewidth = 0.8, fullrange = TRUE) +
    scale_color_manual(values = custom_colors, labels = source_labels) +
    scale_fill_manual(values  = custom_colors, labels = source_labels) +
    scale_linetype_manual(values = custom_linetype, labels = c("FEAST", "mST2")) +
    scale_shape_manual(values = c("FEAST" = 22, "mST2" = 21)) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "right")
}

source_legend   <- get_legend(make_legend_plot("EC23S") + guides(linetype = "none", shape = "none"))
tool_legend     <- get_legend(make_legend_plot("EC23S") + guides(color = "none", fill = "none"))
combined_legend <- plot_grid(source_legend, tool_legend, ncol = 1,
                             rel_heights = c(2, 1), align = "v", axis = "lr")
legend_element  <- wrap_elements(full = combined_legend)

combined_plot <- (smooth_plots_labeled[[1]] +  # EC23S
                  smooth_plots_labeled[[2]] +  # HF183
                  smooth_plots_labeled[[3]] +  # CowM2 (boxplot)
                  smooth_plots_labeled[[4]] +  # BacCan
                  smooth_plots_labeled[[5]] +  # Gull4
                  legend_element) +
  plot_layout(design = "ABC\nDEF") +
  plot_annotation(theme = theme(plot.margin = margin(5, 5, 5, 5)))

print(combined_plot)
ggsave("tiff_plots/markers.tiff", combined_plot, width = 13, height = 7.5, dpi = 1200)

# ---------------------------------------------------------------------------------------------------------------#
############## Step 4: ANCOM-BC2 Differential Abundance Analyses ##############
# ---------------------------------------------------------------------------------------------------------------#

# --- 4a: Shared metadata ---
ancom_metadata <- read.table("ancom_metadata.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
  mutate(Sample = trimws(Sample)) %>%
  column_to_rownames("Sample")

# --- 4b: Shared filtering function ---
prevalence_thresholds  <- c(1.0)
rel_abund_thresholds   <- c(1e-4)
min_reads_thresholds   <- c(10)

generate_filtered_list <- function(otus_mat, samples, prefix) {
  otus_sub <- otus_mat[, samples, drop = FALSE]
  out_list <- list()
  for (p in prevalence_thresholds) {
    for (min_reads in min_reads_thresholds) {
      taxa_prev <- rowMeans(otus_sub >= min_reads)
      keep_prev <- which(taxa_prev >= p)
      if (length(keep_prev) == 0) next
      otus_prev <- otus_sub[keep_prev, , drop = FALSE]
      for (ra in rel_abund_thresholds) {
        rel_ab   <- sweep(otus_prev, 2, colSums(otus_prev), FUN = "/")
        rel_ab[is.na(rel_ab)] <- 0
        keep_ra  <- which(apply(rel_ab, 1, function(x) any(x >= ra)))
        if (length(keep_ra) == 0) next
        otus_filtered <- otus_prev[keep_ra, , drop = FALSE]
        otus_name     <- paste0(prefix, "_p", round(p * 100),
                                "_reads", min_reads,
                                "_ra", format(ra, scientific = FALSE))
        out_list[[otus_name]] <- otus_filtered
      }
    }
  }
  return(out_list)
}

drop_samples  <- c("CHICKEN_SCH_SHER_9_7_21", "SCOGGINS_DAMDIS_8_9_22") # ANCOM-BC2 requires >= samples
run_locations <- c("State", "GASite", "ORSite", "QuartMile")
corrections   <- c("bonferroni")

get_samples_for_run <- function(loc, metadata, drop_samples) {
  all_samples <- rownames(metadata)
  switch(loc,
    "State"     = all_samples,
    "GASite"    = rownames(subset(metadata, State == "Georgia" & !rownames(metadata) %in% drop_samples)),
    "ORSite"    = rownames(subset(metadata, State == "Oregon"  & !rownames(metadata) %in% drop_samples)),
    "QuartMile" = rownames(subset(metadata, !is.na(QuartMile))),
    stop(paste("Unknown loc:", loc))
  )
}

get_meta_col <- function(loc) {
  switch(loc,
    "State"     = "State",
    "GASite"    = "GASite",
    "ORSite"    = "ORSite",
    "QuartMile" = "QuartMile",
    stop(paste("Unknown loc:", loc))
  )
}

# --- 4c: Single ANCOM-BC2 run loop, parameterized by taxonomic level ---
run_ancombc2_level <- function(ancom_otus, level_label, out_suffix) {
  summary_sig <- data.frame()
  run_counter <- 1
  total_runs  <- length(run_locations) * length(prevalence_thresholds) *
                 length(min_reads_thresholds) * length(rel_abund_thresholds) * length(corrections)
  for (loc in run_locations) {
    samples_use <- intersect(get_samples_for_run(loc, ancom_metadata, drop_samples),
                             colnames(ancom_otus))
    if (length(samples_use) < 2) next
    otus_filtered_list <- generate_filtered_list(ancom_otus, samples_use, prefix = loc)
    meta_col <- get_meta_col(loc)
    for (otus_name in names(otus_filtered_list)) {
      otus_sub     <- otus_filtered_list[[otus_name]][, samples_use, drop = FALSE]
      meta_sub_run <- ancom_metadata[samples_use, , drop = FALSE]
      prevalence_val <- as.numeric(gsub(".*_p(\\d+)_.*",       "\\1", otus_name))
      min_reads_val  <- as.numeric(gsub(".*_reads(\\d+)_.*",   "\\1", otus_name))
      rel_abund_val  <- as.numeric(gsub(".*_ra([0-9.e-]+)",    "\\1", otus_name))
      for (corr in corrections) {
        cat(sprintf("Running ANCOM-BC2 (%s): run %d of %d (%s, %s, %s)\n",
                    level_label, run_counter, total_runs, loc, otus_name, corr))
        ancom_res <- ancombc2(
          data          = otus_sub,
          meta_data     = meta_sub_run,
          fix_formula   = meta_col,
          group         = meta_col,
          taxa_are_rows = TRUE,
          p_adj_method  = corr,
          struc_zero    = FALSE,
          neg_lb        = FALSE,
          alpha         = 0.01,
          verbose       = FALSE,
          pairwise      = FALSE
        )
        res_df         <- ancom_res$res
        site_diff_cols <- res_df %>% select(starts_with("diff_")) %>% select(-contains("Intercept"))
        res_df$Significant <- apply(site_diff_cols, 1, any)
        n_sig <- sum(res_df$Significant, na.rm = TRUE)
        summary_sig <- rbind(summary_sig, data.frame(
          Location      = loc, Prevalence = prevalence_val,
          Min_Reads     = min_reads_val, Rel_Abundance = rel_abund_val,
          Correction    = corr, Taxa = n_sig
        ))
        sig_taxa <- res_df[res_df$Significant == TRUE, , drop = FALSE]
        if (nrow(sig_taxa) > 0) {
          write.csv(sig_taxa,
                    sprintf("ancombc2_results/ancombc2_%s%s_taxa.csv", loc, out_suffix),
                    row.names = FALSE)
        }
        run_counter <- run_counter + 1
      }
    }
  }
  return(summary_sig)
}

set.seed(123)

ancom_otus_family <- load_bracken_files("br_family",  "_family.br",  "taxon", "new_est_reads", "_family\\.br$")
summary_sig_family <- run_ancombc2_level(ancom_otus_family, "family",  "")

# --- 4d: Bar plot function ---
make_lfc_plot <- function(csv_path, y_label, x_limits = c(-1.5, 1.5)) {
  taxa_df  <- read_csv(csv_path, show_col_types = FALSE)
  top_taxa <- taxa_df %>%
    mutate(abs_lfc = abs(lfc_StateOregon)) %>%
    slice_max(order_by = abs_lfc, n = 12) %>%
    arrange(lfc_StateOregon) %>%
    mutate(Taxon = factor(taxon, levels = taxon),
           State = ifelse(lfc_StateOregon > 0, "Oregon", "Georgia"))

  ggplot(top_taxa, aes(x = lfc_StateOregon, y = Taxon, fill = State)) +
    geom_col(color = "black", linewidth = 0.4) +
    scale_fill_manual(values = c("Georgia" = "#D62728", "Oregon" = "#1F77B4")) +
    geom_vline(xintercept = 0, color = "black") +
    shadowtext::geom_shadowtext(
      aes(x = lfc_StateOregon / 2, label = round(lfc_StateOregon, 2)),
      color = "white", bg.color = "black", bg.r = 0.1,
      vjust = 0.5, size = 5, fontface = "bold"
    ) +
    labs(x = "Log-Fold Change", y = y_label, fill = "Higher in") +
    scale_x_continuous(limits = x_limits, breaks = seq(x_limits[1], x_limits[2], 0.3)) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.y  = element_text(face = "italic", size = 10),
      axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 15)),
      axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 15))
    )
}

# Family
state_taxa_plot <- make_lfc_plot("ancombc2_results/family/state_family.csv", "Family")
state_taxa_plot
save_tiff(state_taxa_plot, "tiff_plots/state_tax.tiff")

# Order
state_taxa_order_plot <- make_lfc_plot("ancombc2_results/order/state_order.csv", "Order", 
                                       x_limits = c(-1.3, 1.3))
state_taxa_order_plot
save_tiff(state_taxa_order_plot, "tiff_plots/state_tax_order.tiff")

# Class
state_taxa_class_plot <- make_lfc_plot("ancombc2_results/class/state_class.csv", "Class",
                                        x_limits = c(-1.2, 1.2))
state_taxa_class_plot
save_tiff(state_taxa_class_plot, "tiff_plots/state_tax_class.tiff")

# --- 4e: Heatmap function ---
plot_ancom_heatmap <- function(file_path, top_n = 12, save_tiff = NULL, custom_title = NULL) {
  heat_df  <- read_csv(file_path, show_col_types = FALSE)
  lfc_cols <- grep("^lfc_", colnames(heat_df), value = TRUE)
  if (length(lfc_cols) == 0) stop("No lfc_ columns found in ", file_path)

  heat_df <- heat_df %>%
    rowwise() %>%
    mutate(max_abs_lfc = max(abs(c_across(all_of(lfc_cols))), na.rm = TRUE)) %>%
    ungroup()

  top_taxa <- heat_df %>%
    slice_max(order_by = max_abs_lfc, n = top_n) %>%
    arrange(max_abs_lfc) %>%
    select(taxon, all_of(lfc_cols))

  colnames(top_taxa)[-1] <- str_remove(colnames(top_taxa)[-1], "^lfc_")
  top_taxa <- top_taxa %>%
    mutate(taxon = factor(taxon, levels = taxon)) %>%
    column_to_rownames("taxon")

  mat        <- as.matrix(top_taxa)
  min_val    <- min(mat, na.rm = TRUE)
  max_val    <- max(mat, na.rm = TRUE)
  breaks     <- c(seq(min_val, 0, length.out = 51), seq(0, max_val, length.out = 51)[-1])
  title_text <- if (!is.null(custom_title)) custom_title else tools::file_path_sans_ext(basename(file_path))

  if (!is.null(save_tiff)) {
    tiff(save_tiff, width = 10, height = 8, units = "in", res = 1200, compression = "lzw")
    pheatmap(mat,
             cluster_rows  = FALSE, cluster_cols = TRUE,
             color         = colorRampPalette(c("#006CD1", "#FFFFFF", "#D62728"))(100),
             breaks        = breaks, border_color = "black",
             labels_row    = sapply(rownames(mat), function(x)
                               parse(text = paste0("italic('", x, "')"))),
             main          = title_text,
             fontsize      = 16, fontsize_row = 14, fontsize_col = 14)
    dev.off()
  }
}

# --- 4f: Run heatmaps for all three levels ---
heat_level_configs <- list(
  family = list(
    files      = c("ancombc2_results/family/Land Use.csv",
                   "ancombc2_results/family/a) Georgia Sites.csv",
                   "ancombc2_results/family/b) Oregon Sites.csv"),
    tiff_prefix = "tiff_plots/heat_"
  ),
  order  = list(
    files      = c("ancombc2_results/order/Land Use.csv",
                   "ancombc2_results/order/a) Georgia Sites.csv",
                   "ancombc2_results/order/b) Oregon Sites.csv"),
    tiff_prefix = "tiff_plots/heat_order_"
  ),
  class  = list(
    files      = c("ancombc2_results/class/Land Use.csv",
                   "ancombc2_results/class/a) Georgia Sites.csv",
                   "ancombc2_results/class/b) Oregon Sites.csv"),
    tiff_prefix = "tiff_plots/heat_class_"
  )
)

for (cfg in heat_level_configs) {
  for (file in cfg$files) {
    name      <- tools::file_path_sans_ext(basename(file))
    tiff_path <- paste0(cfg$tiff_prefix, tolower(gsub(" ", "_", name)), ".tiff")
    plot_ancom_heatmap(file, top_n = 12, save_tiff = tiff_path, custom_title = name)
  }
}

# ---------------------------------------------------------------------------------------------------------------#
############## Step 5: Facet NMDS Plots and Heatmaps ##############
# ---------------------------------------------------------------------------------------------------------------#

# --- 5a: Edit NMDS plots for faceting (legend stripped) ---
luga_edit <- ggplot(lu_results$GA_lu$coords, aes(x = MDS1, y = MDS2, color = Group)) +
  geom_point(aes(shape = Site), size = 2.2) +
  geom_path(data = ellipse_df_luga, aes(x = x, y = y, color = Group),
            linetype = "dashed", linewidth = 1) +
  scale_color_manual(values = lu_colors, guide = "none") +
  scale_shape_manual(values = c(0,2,6,5,8,4)) +
  # Cap, Cha, Cow, Joh, Uto, Wil
  xlim(-1, 1) + ylim(-1, 1) +
  theme_light(base_size = 18) +
  ggtitle(paste0("b) Georgia Land Use NMDS (Stress = ", lu_results$GA_lu$stress, ")")) +
  labs(x = "NMDS1", y = "NMDS2", colour = "Land Use", shape = "Site") +
  theme_nmds(legend_position = "none", title_size = 16) +
  theme(plot.title = element_text(hjust = -0.05))
print(luga_edit)
save_tiff(luga_edit, "tiff_plots/luga_edit.tiff", width = 6, height = 6)

luor_edit <- ggplot(lu_results$OR_lu$coords, aes(x = MDS1, y = MDS2, color = Group)) +
  geom_point(aes(shape = Site), size = 2.6) +
  geom_path(data = ellipse_df_luor, aes(x = x, y = y, color = Group),
            linetype = "dashed", linewidth = 1) +
  scale_color_manual(values = lu_colors, guide = "none") +
  scale_shape_manual(values = c(7,8,9,10,11,12,13,14,16,3)) +
  xlim(-0.8, 0.8) + ylim(-0.8, 0.8) +
  theme_light(base_size = 18) +
  ggtitle(paste0("c) Oregon Land Use NMDS (Stress = ", lu_results$OR_lu$stress, ")")) +
  labs(x = "NMDS1", y = "NMDS2", colour = "Land Use", shape = "Site") +
  theme_nmds(legend_position = "none", title_size = 16) +
  theme(plot.title = element_text(hjust = -0.05))
print(luor_edit)
save_tiff(luor_edit, "tiff_plots/luor_edit.tiff", width = 6, height = 6)

# --- 5b: Assemble NMDS facet with magick ---
top     <- image_read("tiff_plots/nmds_state.tiff")
bottom1 <- image_read("tiff_plots/luga_edit.tiff")
bottom2 <- image_read("tiff_plots/luor_edit.tiff")

bottom     <- image_trim(image_append(c(bottom1, bottom2)))
top_padded <- image_extent(top,
                            geometry = geometry_area(width  = image_info(bottom)$width,
                                                     height = image_info(top)$height),
                            gravity  = "Center")
combined_nmds <- image_trim(image_append(c(top_padded, bottom), stack = TRUE))

# --- 5c: Land use legend ---
lu_levels <- names(lu_colors)
lu_legend <- get_legend(
  ggplot(data.frame(Group = factor(lu_levels, levels = lu_levels), x = 1, y = 1),
         aes(x = x, y = y, color = Group)) +
    geom_point(size = 3) +
    scale_color_manual(values = lu_colors) +
    guides(color = guide_legend(title = "Land Use", title.position = "top",
                                nrow = 1, byrow = TRUE)) +
    theme_void() +
    theme(legend.position = "bottom",
          legend.text      = element_text(size = 14, face = "bold"),
          legend.title     = element_text(size = 16, face = "bold"))
)
plot(lu_legend)
save_tiff(lu_legend, "tiff_plots/lu_legend.tiff", width = 6, height = 1)

lu_legend_img <- image_trim(image_read("tiff_plots/lu_legend.tiff"))

# --- 5d: Combine NMDS and legend ---
nm_width     <- image_info(combined_nmds)$width
legend_width <- image_info(lu_legend_img)$width
legend_height <- image_info(lu_legend_img)$height
nudge_x      <- 150
vertical_pad <- 150
pad_bottom   <- 150

blank_canvas   <- image_blank(width = nm_width,
                               height = legend_height + vertical_pad + pad_bottom, color = "white")
center_offset  <- round((nm_width - legend_width) / 2)
lu_legend_nudged <- image_composite(blank_canvas, lu_legend_img,
                                     offset = paste0("+", center_offset + nudge_x, "+", vertical_pad))
nmds_facet <- image_append(c(combined_nmds, lu_legend_nudged), stack = TRUE)
magick::image_write(nmds_facet, path = "tiff_plots/nmds_combined.tiff",
                    format = "tiff", compression = "lzw", density = "1200x1200")

# --- 5e: Heatmap facets ---
heat_lu <- image_read("tiff_plots/family/heat_lu.tiff")
heat_ga <- image_read("tiff_plots/family/heat_ga.tiff")
heat_or <- image_read("tiff_plots/family/heat_or.tiff")

# Equalize widths
max_width <- max(image_info(heat_lu)$width, image_info(heat_ga)$width, image_info(heat_or)$width)
heat_lu   <- image_extent(heat_lu, geometry = geometry_size_pixels(width = max_width), gravity = "Center")
heat_ga   <- image_extent(heat_ga, geometry = geometry_size_pixels(width = max_width), gravity = "Center")
heat_or   <- image_extent(heat_or, geometry = geometry_size_pixels(width = max_width), gravity = "Center")

# Full facet (land use and states)
bottom_row      <- image_trim(image_append(c(heat_ga, heat_or)))
heat_lu_centered <- image_extent(heat_lu,
                                  geometry = geometry_size_pixels(width = image_info(bottom_row)$width),
                                  gravity  = "Center")
facet_combined  <- image_append(c(heat_lu_centered, bottom_row), stack = TRUE)
facet_padded    <- image_extent(facet_combined,
                                 geometry = geometry_size_pixels(
                                   width  = image_info(facet_combined)$width,
                                   height = image_info(facet_combined)$height + 100),
                                 gravity = "North")
image_write(facet_padded, path = "tiff_plots/heat_combined.tiff",
            format = "tiff", compression = "lzw", density = "1200x1200")

# States-only facet
max_w_state <- max(image_info(heat_ga)$width, image_info(heat_or)$width)
heat_ga_s   <- image_extent(heat_ga, geometry = geometry_size_pixels(width = max_w_state), gravity = "Center")
heat_or_s   <- image_extent(heat_or, geometry = geometry_size_pixels(width = max_w_state), gravity = "Center")
states_combined <- image_append(c(heat_ga_s, heat_or_s), stack = TRUE)
states_padded   <- image_extent(states_combined,
                                 geometry = geometry_size_pixels(
                                   width  = image_info(states_combined)$width,
                                   height = image_info(states_combined)$height + 100),
                                 gravity = "North")
image_write(states_padded, path = "tiff_plots/heat_states.tiff",
            format = "tiff", compression = "lzw", density = "1200x1200")

# Maps facet
map_ga <- image_read("tiff_maps/or_full.tiff")
map_or <- image_read("tiff_maps/or_facet.tiff")
max_width_map <- max(image_info(map_ga)$width, image_info(map_or)$width)
map_ga <- image_extent(map_ga, geometry = geometry_size_pixels(width = max_width_map), gravity = "Center")
map_or <- image_extent(map_or, geometry = geometry_size_pixels(width = max_width_map), gravity = "Center")
spacer <- image_blank(width = max_width_map, height = 100, color = "white")
facet_map_padded <- image_extent(
  image_append(c(map_ga, spacer, map_or), stack = TRUE),
  geometry = geometry_size_pixels(
    width  = image_info(image_append(c(map_ga, spacer, map_or), stack = TRUE))$width,
    height = image_info(image_append(c(map_ga, spacer, map_or), stack = TRUE))$height + 100),
  gravity = "North")
image_write(facet_map_padded, path = "tiff_maps/facet_maps.tiff",
            format = "tiff", compression = "lzw", density = "1200x1200")

# ---------------------------------------------------------------------------------------------------------------#
############## Step 6: SourceApp MAGs Plot Between States ##############
# ---------------------------------------------------------------------------------------------------------------#
sourceapp_results <- read_csv("sourceapp_results.csv")
sourceapp_results$site <- factor(sourceapp_results$site, levels = unique(sourceapp_results$site))

landuse_colors <- c("Agriculture"  = "#FF7F0E", "Industrial"   = "#D62728",
                    "Forest/Parks" = "#2CA02C", "Residential"  = "#1F77B4")

sample_counts <- sourceapp_results %>%
  group_by(state, site) %>%
  summarise(n_samples = n(), max_mags = max(mags), .groups = "drop") %>%
  mutate(label_y = max_mags + 1.6)

transition_lines <- tribble(
  ~state,    ~xintercept, ~landuse,
  "Georgia",  0.55,       "Agriculture",
  "Georgia",  1.5,        "Forest/Parks",
  "Georgia",  3.5,        "Industrial",
  "Georgia",  4.5,        "Residential",
  "Oregon",   0.65,       "Agriculture",
  "Oregon",   4.5,        "Forest/Parks",
  "Oregon",   5.5,        "Industrial",
  "Oregon",   6.5,        "Residential"
) %>%
  mutate(line_color = landuse_colors[landuse])

hf183_vals  <- sort(unique(na.omit(sourceapp_results$hf183)))
all_breaks  <- sort(unique(c(pretty(hf183_vals, n = 4), min(hf183_vals), max(hf183_vals))))

new_mags_plot <- ggplot(sourceapp_results, aes(x = site, y = mags)) +
  geom_vline(data = transition_lines, aes(xintercept = xintercept),
             color = transition_lines$line_color, linetype = "dashed", linewidth = 1) +
  geom_text(data = transition_lines, x = transition_lines$xintercept,
            label = transition_lines$landuse, y = 16, angle = 90,
            vjust = -0.3, hjust = 1, size = 3.5, color = "black", inherit.aes = FALSE) +
  geom_jitter(width = 0.2, height = 0, shape = 21, size = 6,
              aes(fill = hf183), color = "black", alpha = 0.8) +
  geom_hline(yintercept = 1, color = "black", linetype = "dashed", linewidth = 1) +
  geom_text(data = sample_counts, aes(x = site, y = label_y, label = paste0("n=", n_samples)),
            inherit.aes = FALSE, size = 5, fontface = "plain") +
  facet_wrap(~state, ncol = 1, scales = "free_x") +
  scale_fill_viridis_c(name      = "HF183 [cp/mL]", option = "viridis", direction = -1,
                       na.value  = "grey60", breaks = all_breaks,
                       labels    = scales::label_comma(accuracy = 0.01)(all_breaks)) +
  scale_y_continuous(limits = c(0, 16)) +
  labs(x = "Site", y = "Detected Wastewater MAGs") +
  theme_bw() +
  theme(axis.text.x      = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y      = element_text(size = 16),
        axis.title.y     = element_text(angle = 90, vjust = 0.5, size = 18,
                                        face = "bold", margin = margin(r = 15)),
        axis.title.x     = element_text(size = 18, face = "bold", margin = margin(t = 15)),
        strip.text       = element_text(size = 18, face = "bold", color = "white"),
        strip.background = element_rect(fill = "#9467BD"),
        legend.position  = "right",
        legend.title     = element_text(size = 16, face = "bold", margin = margin(b = 20)),
        legend.text      = element_text(size = 14),
        legend.key.size  = unit(1.2, "cm"),
        panel.spacing    = unit(0.5, "lines"))
print(new_mags_plot)
ggsave("tiff_plots/mags_plot.tiff", new_mags_plot, width = 10, height = 12, dpi = 1200)