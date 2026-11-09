


library(data.table)
library(ggplot2)
library(tidyverse)
library(PRROC)

get_reps <- function(S) {
  if (S %in% c("S7", "S8")) {
    return(1:3)
  }
  1:5
}

args <- commandArgs(trailingOnly = TRUE)
scenarios <- if (length(args) > 0) args else c("S7", "S8")

# Define methods and scenarios
method <- c("TPR", "FPR")

all_eval <- data.frame()
# Loop over each scenario
for (S in scenarios) {

  FPR <- c()
  TPR <- c()

  correlation_matrix <- matrix(0, nrow = 150, ncol = 150)
  if (S %in% c("S4","S5","S6")){
      correlation_matrix[1:5, 1:5] <- 1
      correlation_matrix[6:10, 6:10] <- 1
      correlation_matrix[11:15, 11:15] <- 1
      }
  for (i in get_reps(S)) {
    
    res_SpaceX <- fread(paste0("./Data/", S, "_gene/SpaceX_", i, ".csv"))
    threshold <- 0.25
    true_binary <- ifelse(abs(correlation_matrix) >= threshold, 1, 0)
    estimated_binary <- ifelse(abs(res_SpaceX) >= threshold, 1, 0)
    false_positives <- sum((true_binary == 0) & (estimated_binary == 1))
    total_negatives <- sum(true_binary == 0)  # Count all true negatives
    true_positives <- sum((true_binary == 1) & (estimated_binary == 1))
    actual_positives <- sum(true_binary == 1) 
    FPR <- c(FPR, ifelse(total_negatives == 0, 0, false_positives / total_negatives))
    TPR <- c(TPR, ifelse(actual_positives == 0, 0, true_positives / actual_positives))
  }
  rep_count <- length(get_reps(S))
  eval <- data.frame(
    method = rep(method, each = rep_count),
    layer = rep(seq_len(rep_count), times = length(method)),
    eval = c(TPR,FPR),
    scenario = S)

   all_eval <- rbind(all_eval, eval)


}


custom_colors <- c(
  "TPR" = "#4c9390", 
  "FPR" = "#d65d48"            
)

pdf("./coexpression/SpaceX_eval.pdf", width = 15, height = 6)
ggplot(all_eval, aes(x = method, y = eval, fill = method)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.7) +  # Add jitter for data points
  facet_wrap(~ scenario, nrow = 1, scales = "free_x") +
  scale_fill_manual(values = custom_colors,labels = c("False Positive Rate (FPR)","True Positive Rate (TPR)")) + 
  scale_y_continuous(position = "left") +
  labs(y = "Rate",
       fill = "Evaluation Metrics") +
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
