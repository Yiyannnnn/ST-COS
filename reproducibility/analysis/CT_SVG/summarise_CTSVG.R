library(tidyverse)
# Create empty list to store results
all_results <- list()

# Loop through scenarios
for (s in 1:6) {
    base_dir <- file.path("./Data", paste0("S", s, "_gene"))
    
    # Loop through replicates
    for (rep in 1:5) {
        # Try CELINA files
        celina_file <- file.path(base_dir, paste0("CELINA_", rep, ".RDS"))
        if (file.exists(celina_file)) {
            tryCatch({
                celina_data <- readRDS(celina_file)
                
                # Process CELINA data
                for (ct in c("ct1", "ct2", "ct3")) {
                    if (!is.null(celina_data[[ct]])) {
                        all_results[[length(all_results) + 1]] <- list(
                            method = "CELINA",
                            Scenario = s,
                            rep = rep,
                            cell_type = ct,
                            gene_list = celina_data[[ct]]
                        )
                    } else {
                        all_results[[length(all_results) + 1]] <- list(
                            method = "CELINA",
                            Scenario = s,
                            rep = rep,
                            cell_type = ct,
                            gene_list = NULL
                        )}
                }
            }, error = function(e) {
                message(sprintf("Could not process CELINA file for S%d, rep %d", s, rep))
            })
        } else {
            for (ct in c("ct1", "ct2", "ct3")) {
                all_results[[length(all_results) + 1]] <- list(
                    method = "CELINA",
                    Scenario = s,
                    rep = rep,
                    cell_type = ct,
                    gene_list = NULL
                )
            }
        }
        
        # Try CSIDE files
        cside_file <- file.path(base_dir, paste0("CSIDE_", rep, ".RDS"))
        if (file.exists(cside_file)) {
            tryCatch({
                cside_data <- readRDS(cside_file)
                
                # Process CSIDE data
                for (ct in c("ct1", "ct2", "ct3")) {
                    if (!is.null(cside_data[[ct]])) {
                        all_results[[length(all_results) + 1]] <- list(
                            method = "CSIDE",
                            Scenario = s,
                            rep = rep,
                            cell_type = ct,
                            gene_list = cside_data[[ct]]
                        )
                    } else {
                        all_results[[length(all_results) + 1]] <- list(
                            method = "CSIDE",
                            Scenario = s,
                            rep = rep,
                            cell_type = ct,
                            gene_list = NULL
                        )
                    }
                }
            }, error = function(e) {
                message(sprintf("Could not process CSIDE file for S%d, rep %d", s, rep))
            })
        } else {
            for (ct in c("ct1", "ct2", "ct3")) {
                all_results[[length(all_results) + 1]] <- list(
                    method = "CSIDE",
                    Scenario = s,
                    rep = rep,
                    cell_type = ct,
                    gene_list = NULL
                )
            }
        }
    }
}

# Convert to data frame
results_df <- do.call(rbind, lapply(all_results, function(x) {
    data.frame(
        method = x$method,
        Scenario = x$Scenario,
        rep = x$rep,
        cell_type = x$cell_type,
        gene_list = ifelse(length(x$gene_list) > 0, paste(x$gene_list, collapse = ","), 'None')
    )
}))

# Create true gene sets
true_genes <- list(
    ct1 = paste0("gene_", 1:5),
    ct2 = paste0("gene_", 6:10),
    ct3 = paste0("gene_", 11:15)
)

# Calculate performance metrics
results_df <- results_df %>%
    mutate(
        predicted_genes = strsplit(as.character(gene_list), ","),
        true_set = case_when(
            Scenario %in% c(2, 5) ~ list(character(0)),
            cell_type == "ct1" ~ list(true_genes$ct1),
            cell_type == "ct2" ~ list(true_genes$ct2),
            cell_type == "ct3" ~ list(true_genes$ct3)
        ),
        has_positive_truth = lengths(true_set) > 0,
        TP = mapply(function(pred, true) {
            if (identical(pred, "None")) return(0)
            length(intersect(pred, true))
        }, predicted_genes, true_set),
        FP = mapply(function(pred, true) {
            if (identical(pred, "None")) return(0)
            length(setdiff(pred, true))
        }, predicted_genes, true_set),
        FN = mapply(function(pred, true) {
            if (identical(pred, "None")) return(length(true))
            length(setdiff(true, pred))
        }, predicted_genes, true_set),
        TN = mapply(function(pred, true) {
            if (identical(pred, "None")) return(150 - length(true))
            150 - length(union(pred, true))
        }, predicted_genes, true_set),
        precision = TP / (TP + FP),
        recall = TP / (TP + FN),
        F1 = 2 * (precision * recall) / (precision + recall)
    ) %>%
    mutate(
        precision = if_else(has_positive_truth, replace_na(precision, 0), NA_real_),
        recall = if_else(has_positive_truth, replace_na(recall, 0), NA_real_),
        F1 = if_else(has_positive_truth, replace_na(F1, 0), NA_real_)
    )

    # Count non-zero F1 scores by method
    f1_counts <- results_df %>%
        group_by(method, Scenario, rep) %>%
        summarise(
            nonzero_f1 = sum(F1 > 0),
            total = n(),
            proportion = nonzero_f1 / total
        )

        # Set custom colors
        custom_colors <- c("#E69F00", "#56B4E9")

        # Plot F1 scores
        pdf("CTSVG_F1_scores.pdf", width = 24, height = 6)
        ggplot(results_df, aes(x = method, y = F1, fill = method)) +
            geom_boxplot(outlier.shape = NA) +
            geom_jitter(width = 0.2, size = 1, alpha = 0.7) +
            facet_wrap(~ Scenario, nrow = 1, scales = "free_x") +
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
                legend.text = element_text(size = 20),
                legend.title = element_text(size = 25)
            )
        dev.off()
# Save results
saveRDS(results_df, "CTSVG_summary.RDS")
