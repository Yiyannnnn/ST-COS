

source("./Data/simulator_YZ.R")

args <- commandArgs(trailingOnly = TRUE)
rep_id <- if (length(args) >= 1L) suppressWarnings(as.integer(args[1])) else 1L
if (is.na(rep_id) || !rep_id %in% 1:5) stop("replicate must be an integer from 1 to 5")
rep_seed <- 10000L + rep_id
set.seed(rep_seed)
svg_num=15
noise_num=135

x_len <- 100
y_len <- 100
nGenes <- 150


gene_pattern = data.frame(ct1 = c(rep('layer',5),rep('no_pattern',145)),
                          ct2= c(rep('no_pattern',5),rep('layer',5),rep('no_pattern',140)),
                          ct3 = c(rep('no_pattern',10),rep('layer',5),rep('no_pattern',135)))



layer_param = data.frame(ct1 = c(1/3,1/3,1/3),
                          ct2 = c(1/3,1/3,1/3),
                          ct3 = c(1/3,1/3,1/3))

#define layer
xmin <- c(0, 33, 66)
xmax <- c(33, 66, 100)
ymin <- c(0, 0, 0)
ymax <- c(100, 100, 100)

rownames(gene_pattern) = paste0('gene',1:nGenes)

coord = generate_coordinates(x_len,y_len,'hex',spot_distance=5,seed=rep_seed)
mean_param = generate_gene_pattern(coord=coord,gene_pattern=gene_pattern,seed=rep_seed)
true_count <- generate_true_count(coord,mean_param,gene_relationship = NULL)
layer_coord <- generate_layer(coord, xmin, xmax, ymin, ymax)
celltype_proportion <- generate_cell_prop(layer_coord,layer_param)
true_count <- generate_batch_count(
  true_count, p0 = 0.05, celltype_proportion = celltype_proportion,
  seq_depth = 1e6
)
spot_count_layers <- generate_spot_count_abundance_weighted(
  true_count = true_count, nGenes = nGenes, cell_number = 10,
  celltype_proportion = celltype_proportion
)

rownames(coord) <- paste0("spot_",1:nrow(coord))
rownames(spot_count_layers) <- paste0("spot_",1:nrow(spot_count_layers))
colnames(spot_count_layers) <- paste0("gene_",1:ncol(spot_count_layers))
rownames(layer_coord) <- paste0("spot_",1:nrow(layer_coord))


# fwrite(coord,"./Data/S1_gene/simulated/S1_coord_1.csv")
# fwrite(spot_count_layers,"./Data/S1_gene/simulated/S1_gene_expression_1.csv")
# fwrite(layer_coord,"./Data/S1_gene/simulated/S1_layer_1.csv")

# fwrite(coord,"./Data/S1_gene/simulated/S1_coord_2.csv")
# fwrite(spot_count_layers,"./Data/S1_gene/simulated/S1_gene_expression_2.csv")
# fwrite(layer_coord,"./Data/S1_gene/simulated/S1_layer_2.csv")

# fwrite(coord,"./Data/S1_gene/simulated/S1_coord_3.csv")
# fwrite(spot_count_layers,"./Data/S1_gene/simulated/S1_gene_expression_3.csv")
# fwrite(layer_coord,"./Data/S1_gene/simulated/S1_layer_3.csv")

# fwrite(coord,"./Data/S1_gene/simulated/S1_coord_4.csv")
# fwrite(spot_count_layers,"./Data/S1_gene/simulated/S1_gene_expression_4.csv")
# fwrite(layer_coord,"./Data/S1_gene/simulated/S1_layer_4.csv")

# fwrite(coord,"./Data/S1_gene/simulated/S1_coord_5.csv")
# fwrite(spot_count_layers,"./Data/S1_gene/simulated/S1_gene_expression_5.csv")
# fwrite(layer_coord,"./Data/S1_gene/simulated/S1_layer_5.csv")

sc_ref <- rbindlist(lapply(names(true_count), function(x) data.table(Cell_type = x, gene = t(true_count[[x]]))))
rownames(sc_ref) <- paste0("cell_",1:nrow(sc_ref))
colnames(sc_ref)[c(2:ncol(sc_ref))] <- paste0("gene_",1:(ncol(sc_ref)-1))

output_dir <- "./Data/S1_gene/simulated"
realized_prop <- celltype_proportion
region_dt <- cbind(
  layer_coord[, .(x_coord, y_coord)],
  data.table(region = paste0("layer_", layer_coord$layer)),
  as.data.table(realized_prop)
)
fwrite(coord, file.path(output_dir, sprintf("S1_coord_%d.csv", rep_id)))
fwrite(spot_count_layers, file.path(output_dir, sprintf("S1_gene_expression_%d.csv", rep_id)))
fwrite(layer_coord, file.path(output_dir, sprintf("S1_layer_%d.csv", rep_id)))
fwrite(region_dt, file.path(output_dir, sprintf("S1_region_%d.csv", rep_id)))
fwrite(sc_ref, file.path(output_dir, sprintf("S1_scref_%d.csv", rep_id)))

#fwrite(sc_ref,"./Data/S1_gene/simulated/S1_scref_1.csv")
#fwrite(sc_ref,"./Data/S1_gene/simulated/S1_scref_2.csv")
#fwrite(sc_ref,"./Data/S1_gene/simulated/S1_scref_3.csv")
#fwrite(sc_ref,"./Data/S1_gene/simulated/S1_scref_4.csv")
#fwrite(sc_ref,"./Data/S1_gene/simulated/S1_scref_5.csv")

pdf(file.path(output_dir, sprintf("S1_generation_diagnostics_%d.pdf", rep_id)), width = 20)
param_all = c()
for(ct in colnames(gene_pattern)){
  temp = mean_param[[ct]] %>% as.data.table()
  temp$cell_type = ct
  param_all = rbind(param_all,temp)
}
coord_param = cbind(coord,param_all)
coord_param = melt(coord_param,id.vars = c('x_coord','y_coord','cell_type'))
colnames(coord_param) <- c("x_coord","y_coord","cell_type","variable","mean")

coord_param |> filter(variable %in% paste0("gene",1:15)) |> ggplot()+geom_point(aes(x=x_coord,y=y_coord,fill=mean),shape=21,size=1)+
  facet_wrap(cell_type~variable,ncol=15)+
  coord_fixed()+xlim(-1,x_len+1)+ylim(-1,y_len+1)+theme_article()+scale_fill_distiller(palette = "Spectral")

prop_all = c()
for(ct in colnames(celltype_proportion)){
  temp = celltype_proportion[,ct] %>% as.data.table()
  prop_all = cbind(prop_all,temp)
  colnames(prop_all)[ncol(prop_all)] <- ct
}
coord_prop = cbind(coord,prop_all)
coord_prop = melt(coord_prop,id.vars = c('x_coord','y_coord'))
colnames(coord_prop) = c('x_coord','y_coord','cell_type','proportion')
ggplot(coord_prop)+geom_point(aes(x=x_coord,y=y_coord,fill=proportion),shape=21,size=2)+
  facet_wrap(as.factor(coord_prop$cell_type),ncol=3)+
  coord_fixed()+xlim(-1,x_len+1)+ylim(-1,y_len+1)+theme_article()+scale_fill_distiller(palette = "Spectral")
layer_factor = as.numeric(as.factor(layer_coord$layer))


dev.off()
