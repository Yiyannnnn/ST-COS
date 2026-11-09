
library(STANCE)
library(PRROC)
S <- "S1"
STANCE_AUPRC <- c()
STANCE_F1 <- c()
for (suffix in paste0("_",1:5)){
mySTANCE <- readRDS(paste0("./Data/", S, "_gene/STANCE", suffix, ".RDS"))
mySTANCE_1 <- CT_topGenes(object = mySTANCE, 
                        CT_of_interest = "ct1", numTopGenes = 15)
mySTANCE_2 <- CT_topGenes(object = mySTANCE, 
                        CT_of_interest = "ct2", numTopGenes = 15)
mySTANCE_3 <- CT_topGenes(object = mySTANCE, 
                        CT_of_interest = "ct3", numTopGenes = 15)
ct_svg <- c(mySTANCE_1@cell_type_top_genes$top_genes,
mySTANCE_2@cell_type_top_genes$top_genes,
mySTANCE_3@cell_type_top_genes$top_genes)

STANCE_AUPRC <- c(STANCE_AUPRC,pr.curve(scores.class0 = c(paste0("gene_",1:15) %in% ct_svg), 
scores.class1 = c(paste0("gene_",16:150) %in% ct_svg),curve=T)$auc.integral)
    
    tp <- sum(ct_svg %in% paste0("gene_",1:15))
    fp <- sum(!ct_svg %in% paste0("gene_",1:15))
    fn <- sum(!paste0("gene_",1:15) %in% ct_svg)
    STANCE_F1 <- c(STANCE_F1, calc_f1(tp, fp, fn))
}

