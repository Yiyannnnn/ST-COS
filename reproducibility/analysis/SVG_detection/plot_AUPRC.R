library(data.table)
library(ggplot2)
library(tidyverse)
library(PRROC)

get_reps <- function(S) {
  sim_dir <- file.path(
    "./Data",
    paste0(S, "_gene/simulated")
  )
  pattern <- paste0("^", S, "_gene_expression_([0-9]+)\\.csv$")
  files <- list.files(sim_dir, pattern = pattern)
  if (length(files) > 0) {
    reps <- as.integer(sub(pattern, "\\1", files))
    reps <- sort(unique(reps[!is.na(reps)]))
    if (length(reps) > 0) {
      return(reps)
    }
  }
  if (S %in% c("S7", "S8")) {
    return(1:3)
  }
  1:5
}

get_truth_genes <- function(S, rep_id) {
  truth_path <- file.path(
    "./Data",
    paste0(S, "_gene/simulated"),
    paste0(S, "_marker_gene_truth.csv")
  )
  if (file.exists(truth_path)) {
    truth_dt <- fread(truth_path)
    genes <- truth_dt[scenario == rep_id, gene_id]
    if (length(genes) > 0) {
      return(genes)
    }
  }
  if (S %in% c("S2", "S5")) {
    return(character(0))
  }
  if (S %in% c("S1", "S3", "S4", "S6")) {
    return(paste0("gene_", 1:15))
  }
  stop("No SVG truth definition found for ", S, " replicate ", rep_id, ".")
}

args <- commandArgs(trailingOnly = TRUE)
scenarios <- if (length(args) > 0) args else c("S7", "S8")

method <- c("moranI","HEARTSVG","SPARK","SPARKX","SpatialDE","SpatialDE2")
S <- scenarios[1]
if (S %in% c("S2", "S5")) {
  stop("AUPRC is undefined for ", S, " because the scenario has no positive SVG truth.")
}
moranI_AUPRC <- c()
HEARTSVG_AUPRC <- c()
SPARK_AUPRC <- c()
SPARKX_AUPRC <- c()
SpatialDE_AUPRC <- c()
SpatialDE2_AUPRC <- c()

for (i in get_reps(S)) {
  truth_genes <- get_truth_genes(S, i)
  res_moranI <- fread(paste0("./Data/",S,"_gene/","moranI","_",i,".csv"))
  res_moranI$SVG <- ifelse(res_moranI$pval_norm_fdr_bh <= 0.05, 1, 0)
  # moranI_AUPRC <- c(moranI_AUPRC,pr.curve(scores.class0 = res_moranI[res_moranI$V1 %in% paste0("gene_",1:15),]$SVG, 
  # scores.class1 = res_moranI[!res_moranI$V1 %in% paste0("gene_",1:15),]$SVG,curve=T)$auc.integral)
  moranI_AUPRC <- c(moranI_AUPRC,pr.curve(scores.class0 = res_moranI[res_moranI$V1 %in% truth_genes,]$I, 
  scores.class1 = res_moranI[!res_moranI$V1 %in% truth_genes,]$I,curve=T)$auc.integral)

  res_HEARTSVG <- fread(paste0("./Data/",S,"_gene/","HEARTSVG","_",i,".csv"))
  res_HEARTSVG$SVG <- ifelse(res_HEARTSVG$p_adj <= 0.05, 1, 0)
  HEARTSVG_AUPRC <- c(HEARTSVG_AUPRC,pr.curve(scores.class0 = res_HEARTSVG[res_HEARTSVG$gene %in% truth_genes,]$SVG, 
  scores.class1 = res_HEARTSVG[!res_HEARTSVG$gene %in% truth_genes,]$SVG,curve=T)$auc.integral)

  res_SPARK <- fread(paste0("./Data/",S,"_gene/","SPARK","_",i,".csv"))
  res_SPARK$SVG <- ifelse(res_SPARK$adjusted_pvalue <= 0.05, 1, 0)
  SPARK_AUPRC <- c(SPARK_AUPRC,pr.curve(scores.class0 = res_SPARK[res_SPARK$gene %in% truth_genes,]$SVG, 
  scores.class1 = res_SPARK[!res_SPARK$gene %in% truth_genes,]$SVG,curve=T)$auc.integral)

  res_SPARKX <- fread(paste0("./Data/",S,"_gene/","SPARKX","_",i,".csv"))
  res_SPARKX$SVG <- ifelse(res_SPARKX$adjustedPval <= 0.05, 1, 0)
  SPARKX_AUPRC <- c(SPARKX_AUPRC,pr.curve(scores.class0 = res_SPARKX[res_SPARKX$gene %in% truth_genes,]$SVG, 
  scores.class1 = res_SPARKX[!res_SPARKX$gene %in% truth_genes,]$SVG,curve=T)$auc.integral)

  res_SpatialDE <- fread(paste0("./Data/",S,"_gene/","spatialDE","_",i,".csv"))
  SpatialDE_AUPRC <- c(SpatialDE_AUPRC,pr.curve(scores.class0 = res_SpatialDE[res_SpatialDE$g %in% truth_genes,]$FSV, 
  scores.class1 = res_SpatialDE[!res_SpatialDE$g %in% truth_genes,]$FSV,curve=T)$auc.integral)

  res_spatialDE2 <- fread(paste0("./Data/",S,"_gene/","spatialDE2","_",i,".csv"))
  SpatialDE2_AUPRC <- c(SpatialDE2_AUPRC,pr.curve(scores.class0 = res_spatialDE2[res_spatialDE2$V1 %in% truth_genes,]$FSV, 
  scores.class1 = res_spatialDE2[!res_spatialDE2$V1 %in% truth_genes,]$FSV,curve=T)$auc.integral)
  
}


AUPRC <- as.data.frame(t(cbind(moranI_AUPRC,HEARTSVG_AUPRC,SPARK_AUPRC,SPARKX_AUPRC,SpatialDE_AUPRC,SpatialDE2_AUPRC)))
AUPRC$method <- method
AUPRC <- as.data.table(AUPRC)
AUPRC <- data.table::melt(AUPRC, id.vars = "method", variable.name = "layer", value.name = "AUPRC")

pdf(paste0("./SVG_detection/", S, "_AUPRC.pdf"), width = 6, height = 6)
ggplot(AUPRC, aes(x = method, y = AUPRC, fill = method)) +
  geom_boxplot() +
    labs(title = "AUPRC of SVG Detection Methods", x = "Method", y = "AUPRC") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()



library(data.table)
library(ggplot2)
library(tidyverse)
library(PRROC)

method <- c("moranI","HEARTSVG","SPARK","SPARKX","SpatialDE","SpatialDE2")
S <- scenarios[1]

# Initialize empty lists to store F1 scores for each method
moranI_F1 <- c()
HEARTSVG_F1 <- c()
SPARK_F1 <- c()
SPARKX_F1 <- c()
SpatialDE_F1 <- c()
SpatialDE2_F1 <- c()

for (i in get_reps(S)) {
  truth_genes <- get_truth_genes(S, i)
  # Function to calculate F1 score
  calc_f1 <- function(tp, fp, fn) {
    precision <- ifelse((tp + fp) == 0, 0, tp / (tp + fp))
    recall <- ifelse((tp + fn) == 0, 0, tp / (tp + fn))
    ifelse((precision + recall) == 0, 0, 2 * (precision * recall) / (precision + recall))
  }
  
  # Load results and calculate F1 score for each method
  
  # MoranI
  res_moranI <- fread(paste0("./Data/",S,"_gene/","moranI","_",i,".csv"))
  res_moranI$SVG <- ifelse(res_moranI$pval_norm_fdr_bh <= 0.05, 1, 0)
  tp <- sum(res_moranI$SVG[res_moranI$V1 %in% truth_genes] == 1)
  fp <- sum(res_moranI$SVG[!res_moranI$V1 %in% truth_genes] == 1)
  fn <- sum(res_moranI$SVG[res_moranI$V1 %in% truth_genes] == 0)
  moranI_F1 <- c(moranI_F1, calc_f1(tp, fp, fn))

  # HEARTSVG
  res_HEARTSVG <- fread(paste0("./Data/",S,"_gene/","HEARTSVG","_",i,".csv"))
  res_HEARTSVG$SVG <- ifelse(res_HEARTSVG$p_adj <= 0.05, 1, 0)
  tp <- sum(res_HEARTSVG$SVG[res_HEARTSVG$gene %in% truth_genes] == 1)
  fp <- sum(res_HEARTSVG$SVG[!res_HEARTSVG$gene %in% truth_genes] == 1)
  fn <- sum(res_HEARTSVG$SVG[res_HEARTSVG$gene %in% truth_genes] == 0)
  HEARTSVG_F1 <- c(HEARTSVG_F1, calc_f1(tp, fp, fn))

  # SPARK
  res_SPARK <- fread(paste0("./Data/",S,"_gene/","SPARK","_",i,".csv"))
  res_SPARK$SVG <- ifelse(res_SPARK$adjusted_pvalue <= 0.05, 1, 0)
  tp <- sum(res_SPARK$SVG[res_SPARK$gene %in% truth_genes] == 1)
  fp <- sum(res_SPARK$SVG[!res_SPARK$gene %in% truth_genes] == 1)
  fn <- sum(res_SPARK$SVG[res_SPARK$gene %in% truth_genes] == 0)
  SPARK_F1 <- c(SPARK_F1, calc_f1(tp, fp, fn))

  # SPARKX
  res_SPARKX <- fread(paste0("./Data/",S,"_gene/","SPARKX","_",i,".csv"))
  res_SPARKX$SVG <- ifelse(res_SPARKX$adjustedPval <= 0.05, 1, 0)
  tp <- sum(res_SPARKX$SVG[res_SPARKX$gene %in% truth_genes] == 1)
  fp <- sum(res_SPARKX$SVG[!res_SPARKX$gene %in% truth_genes] == 1)
  fn <- sum(res_SPARKX$SVG[res_SPARKX$gene %in% truth_genes] == 0)
  SPARKX_F1 <- c(SPARKX_F1, calc_f1(tp, fp, fn))

  # SpatialDE
  res_SpatialDE <- fread(paste0("./Data/",S,"_gene/","spatialDE","_",i,".csv"))
  res_SpatialDE$SVG <- ifelse(res_SpatialDE$qval <= 0.05, 1, 0)
  tp <- sum(res_SpatialDE$SVG[res_SpatialDE$g %in% truth_genes] == 1)
  fp <- sum(res_SpatialDE$SVG[!res_SpatialDE$g %in% truth_genes] == 1)
  fn <- sum(res_SpatialDE$SVG[res_SpatialDE$g %in% truth_genes] == 0)
  SpatialDE_F1 <- c(SpatialDE_F1, calc_f1(tp, fp, fn))

  # SpatialDE2
  res_spatialDE2 <- fread(paste0("./Data/",S,"_gene/","spatialDE2","_",i,".csv"))
  res_spatialDE2$SVG <- ifelse(res_spatialDE2$FSV <= 0.05, 1, 0)
  tp <- sum(res_spatialDE2$SVG[res_spatialDE2$V1 %in% truth_genes] == 1)
  fp <- sum(res_spatialDE2$SVG[!res_spatialDE2$V1 %in% truth_genes] == 1)
  fn <- sum(res_spatialDE2$SVG[res_spatialDE2$V1 %in% truth_genes] == 0)
  SpatialDE2_F1 <- c(SpatialDE2_F1, calc_f1(tp, fp, fn))
}

F1_scores <- as.data.frame(t(cbind(moranI_F1, HEARTSVG_F1, SPARK_F1, SPARKX_F1, SpatialDE_F1, SpatialDE2_F1)))
F1_scores$method <- method
F1_scores <- as.data.table(F1_scores)
F1_scores <- data.table::melt(F1_scores, id.vars = "method", variable.name = "layer", value.name = "F1_score")

pdf(paste0("./SVG_detection/",S,"_F1_score.pdf"), width = 6, height = 6)
ggplot(F1_scores, aes(x = method, y = F1_score, fill = method)) +
  geom_boxplot() +
  labs(title = "F1 Score of SVG Detection Methods", x = "Method", y = "F1 Score") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

library(data.table)
library(ggplot2)
library(tidyverse)
library(PRROC)

# Define methods and scenarios
method <- c("moranI", "HEARTSVG", "SPARK", "SPARKX", "SpatialDE", "SpatialDE2")

# Initialize an empty data frame to hold F1 scores across all scenarios
all_F1_scores <- data.frame()

# Function to calculate F1 score
calc_f1 <- function(tp, fp, fn) {
  precision <- ifelse((tp + fp) == 0, 0, tp / (tp + fp))
  recall <- ifelse((tp + fn) == 0, 0, tp / (tp + fn))
  ifelse((precision + recall) == 0, 0, 2 * (precision * recall) / (precision + recall))
}

# Loop over each scenario
for (S in scenarios) {
  moranI_F1 <- c()
  HEARTSVG_F1 <- c()
  SPARK_F1 <- c()
  SPARKX_F1 <- c()
  SpatialDE_F1 <- c()
  SpatialDE2_F1 <- c()
  
  # Calculate F1 for each method within each scenario
  for (i in get_reps(S)) {
    truth_genes <- get_truth_genes(S, i)
    # MoranI
    res_moranI <- fread(paste0("./Data/", S, "_gene/moranI_", i, ".csv"))
    res_moranI$SVG <- ifelse(res_moranI$pval_norm_fdr_bh <= 0.05, 1, 0)
    tp <- sum(res_moranI$SVG[res_moranI$V1 %in% truth_genes] == 1)
    fp <- sum(res_moranI$SVG[!res_moranI$V1 %in% truth_genes] == 1)
    fn <- sum(res_moranI$SVG[res_moranI$V1 %in% truth_genes] == 0)
    moranI_F1 <- c(moranI_F1, calc_f1(tp, fp, fn))

    # HEARTSVG
    res_HEARTSVG <- fread(paste0("./Data/", S, "_gene/HEARTSVG_", i, ".csv"))
    res_HEARTSVG$SVG <- ifelse(res_HEARTSVG$p_adj <= 0.05, 1, 0)
    tp <- sum(res_HEARTSVG$SVG[res_HEARTSVG$gene %in% truth_genes] == 1)
    fp <- sum(res_HEARTSVG$SVG[!res_HEARTSVG$gene %in% truth_genes] == 1)
    fn <- sum(res_HEARTSVG$SVG[res_HEARTSVG$gene %in% truth_genes] == 0)
    HEARTSVG_F1 <- c(HEARTSVG_F1, calc_f1(tp, fp, fn))

    # SPARK
    res_SPARK <- fread(paste0("./Data/", S, "_gene/SPARK_", i, ".csv"))
    res_SPARK$SVG <- ifelse(res_SPARK$adjusted_pvalue <= 0.05, 1, 0)
    tp <- sum(res_SPARK$SVG[res_SPARK$gene %in% truth_genes] == 1)
    fp <- sum(res_SPARK$SVG[!res_SPARK$gene %in% truth_genes] == 1)
    fn <- sum(res_SPARK$SVG[res_SPARK$gene %in% truth_genes] == 0)
    SPARK_F1 <- c(SPARK_F1, calc_f1(tp, fp, fn))

    # SPARKX
    res_SPARKX <- fread(paste0("./Data/", S, "_gene/SPARKX_", i, ".csv"))
    res_SPARKX$SVG <- ifelse(res_SPARKX$adjustedPval <= 0.05, 1, 0)
    tp <- sum(res_SPARKX$SVG[res_SPARKX$gene %in% truth_genes] == 1)
    fp <- sum(res_SPARKX$SVG[!res_SPARKX$gene %in% truth_genes] == 1)
    fn <- sum(res_SPARKX$SVG[res_SPARKX$gene %in% truth_genes] == 0)
    SPARKX_F1 <- c(SPARKX_F1, calc_f1(tp, fp, fn))

    # SpatialDE
    res_SpatialDE <- fread(paste0("./Data/", S, "_gene/spatialDE_", i, ".csv"))
    res_SpatialDE$SVG <- ifelse(res_SpatialDE$qval <= 0.05, 1, 0)
    tp <- sum(res_SpatialDE$SVG[res_SpatialDE$g %in% truth_genes] == 1)
    fp <- sum(res_SpatialDE$SVG[!res_SpatialDE$g %in% truth_genes] == 1)
    fn <- sum(res_SpatialDE$SVG[res_SpatialDE$g %in% truth_genes] == 0)
    SpatialDE_F1 <- c(SpatialDE_F1, calc_f1(tp, fp, fn))

    # SpatialDE2
    res_spatialDE2 <- fread(paste0("./Data/", S, "_gene/spatialDE2_", i, ".csv"))
    res_spatialDE2$SVG <- ifelse(res_spatialDE2$FSV >= quantile(res_spatialDE2$FSV, probs = 0.9), 1, 0)
    tp <- sum(res_spatialDE2$SVG[res_spatialDE2$V1 %in% truth_genes] == 1)
    fp <- sum(res_spatialDE2$SVG[!res_spatialDE2$V1 %in% truth_genes] == 1)
    fn <- sum(res_spatialDE2$SVG[res_spatialDE2$V1 %in% truth_genes] == 0)
    SpatialDE2_F1 <- c(SpatialDE2_F1, calc_f1(tp, fp, fn))
  }
  
  # Combine F1 scores for the current scenario
  rep_count <- length(get_reps(S))
  F1_scores <- data.frame(
    method = rep(method, each = rep_count),
    layer = rep(seq_len(rep_count), times = length(method)),
    F1_score = c(moranI_F1, HEARTSVG_F1, SPARK_F1, SPARKX_F1, SpatialDE_F1, SpatialDE2_F1),
    scenario = S
  )
  
  # Append to the overall data frame
  all_F1_scores <- rbind(all_F1_scores, F1_scores)
}

custom_colors <- c(
  "HEARTSVG" = "#e19d92", 
  "moranI" = "#eac47c",             
  "SPARK" = "#dbc0af",      
  "SPARKX" = "#7e8cad",            
  "SpatialDE" = "#a5c2cd", 
  "SpatialDE2" = "#dadae4"
)

# Plot combined F1 score plot for all scenarios
pdf("./SVG_detection/all_scenarios_F1_score.pdf", width = 24, height = 6)
ggplot(all_F1_scores, aes(x = method, y = F1_score, fill = method)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.7) +  # Add jitter for data points
  facet_wrap(~ scenario, nrow = 1, scales = "free_x") +
  scale_fill_manual(values = custom_colors) + 
  scale_y_continuous(position = "left") +
  labs(x = "Method", 
       y = "F1 Score",
       fill = "Methods") +
  theme_classic() +
  theme(
    strip.background = element_rect(fill = "#f0ffff"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),
    strip.text = element_text(size = 25, face = "bold"),
    plot.title = element_text(size = 25, face = "bold"),
    legend.text = element_text(size = 20),        # Change legend text size
    legend.title = element_text(size = 25)        # Change legend title size
  )
dev.off()
