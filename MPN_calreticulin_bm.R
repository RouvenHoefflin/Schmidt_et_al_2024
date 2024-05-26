# Load libraries ----------------------------------------------------------
library(Matrix)
library(Seurat)
library(scalop)
library(ggpubr)
library(ComplexHeatmap)
library(tidyverse)

an_cols <- c("#cb5629",
             "#845bcd",
             "#72be46",
             "#cc55bb",
             "#4a8a35",
             "#d3477d",
             "#59c187",
             "#d1444b",
             "#3fc1bf",
             "#d59032",
             "#6873c5",
             "#bcae37",
             "#984f8a",
             "#a1ae62",
             "#d590ce",
             "#39855f",
             "#b66068",
             "#5f9cd5",
             "#756d28",
             "#c9855b")

an_cols <- c(
  "CLP" = "#e8e79a",
  "late_pro_BC" = "#c2d89a",
  "pre_BC_S_phase" = "#8cbf9a",
  "pre_BC_G2M_phase" = "#5fa2a4",
  "pro_BC" =  "#477b95",
  "pre_BC_small" = "#315b88",
  "im_BC" =  "#24396b",
  "Pcell" = "#191f40",
  "myeloblast" = "#da9fb8", 
  "myelocyte" = "#ae2565",
  "neutrophil" = "#c1447e",
  "monoblast" = "#d9d2cc",
  "monocyte" =  "#adbe7c",
  "macrophage" =  "#6e8537",
  "basophil" = "#cc55bb",
  "T_NK" = "#c3a016",
  "pro_erythroblast" = "#f4a464",
  "erythroblast" =   "#bf2729",
  "lq" = "grey50")

an_cols_condition <- c(
  "mpn" = "#d8537d",
  "normal" = "#6DC5B2"
)

# Load custom functions ---------------------------------------------------
Julie_enricher_mm <- function (test_gene_sets,
                               ref_gene_sets, universe = NULL,
                               minGSSize = 10, 
                               maxGSSize = 500,
                               pAdjustMethod = "BH",
                               pvalueCutoff = 0.05, 
                               qvalueCutoff = 0.2) 
{
  
  
  
  if (is.character(test_gene_sets)) {
    test_gene_sets = list(test_gene_sets)
  }
  
  
  test_gene_sets = sapply(test_gene_sets, function(s) limma::alias2Symbol(s, species = "Mm"), 
                          simplify = F)
  
  universe = limma::alias2Symbol(universe, species = "Mm")
  
  ref_universe = unique(unlist(ref_gene_sets))
  
  universe = intersect(universe, ref_universe)
  
  test_gene_sets = sapply(test_gene_sets, function(x) x[x %in% universe], simplify = F)
  
  ref_gene_sets = sapply(ref_gene_sets, function(x) x[x %in% universe], simplify = F)
  
  term2gene = data.frame(term = names(Unlist(ref_gene_sets)), 
                         gene = Unlist(ref_gene_sets), stringsAsFactors = F)
  
  result = sapply(test_gene_sets, function(gene_set) {
    clusterProfiler::enricher(gene = gene_set, universe = universe, 
                              TERM2GENE = term2gene, minGSSize = minGSSize, maxGSSize = maxGSSize, 
                              pAdjustMethod = pAdjustMethod, pvalueCutoff = pvalueCutoff, 
                              qvalueCutoff = pvalueCutoff)
  }, simplify = F)
  
  sapply(1:length(result), function(i) {
    as.data.frame(result[[i]]) %>% dplyr::mutate(name = rep(names(result)[i], 
                                                            nrow(result[[i]])))
  }, simplify = F)
  
  
}

#to generate the ref in the global env
msigdb_sigs_mm = scalop::msigdb(species = "Mus musculus", category = c("H", "C2", "C3", "C5", "C6") ) %>% 
  unlist(.,  recursive = F) %>% .[(sapply(., function(i) length(i) < 250))]

enricher_func_mm <- function(genes_to_test = NULL,
                             all_genes_within_matrix = rownames(centered_data_mm),
                             ref_sigs_to_add = sigs,
                             title = "delete me",
                             title_for_added_ref_sigs = "Custom",
                             ref_sigs = msigdb_sigs_mm,
                             print_plot = T,
                             save_plot = F){
  
  
  
  if(!base::exists("msigdb_sigs_mm")){
    msigdb_sigs_mm = scalop::msigdb(species = "Mus musculus", category = c("H", "C1", "C2", "C3", "C5", "C6", "C8") ) %>% 
      unlist(.,  recursive = F) %>% .[(sapply(., function(i) length(i) < 250))]
    # %>% sapply(., toupper, simplify = F)
  }
  
  
  if(!is.null(ref_sigs_to_add)){
    names(ref_sigs_to_add) <- paste0(title_for_added_ref_sigs,".",toupper(names(ref_sigs_to_add) ) )
    ref_sigs_to_add <- ref_sigs_to_add %>% sapply(., firstup, simplify = F)
    ref_sigs <- append(msigdb_sigs_mm, ref_sigs_to_add)
  }
  
  
  enrichment_df <- Julie_enricher_mm(test_gene_sets = genes_to_test,
                                     universe =all_genes_within_matrix, 
                                     ref_gene_sets = ref_sigs, )[[1]]
  
  enrichment_df <- enrichment_df[order(enrichment_df$qvalue),] # order the programs by the q-val
  #break the geneID column to string
  enrichment_df <- enrichment_df %>% mutate(geneID = str_split(pattern = "/", simplify = F, string = geneID) )
  enrichment_df_for_plot <- enrichment_df[1:50, c("Description", "qvalue", "GeneRatio", "BgRatio")] # The only the top 50 programs and the relevant data
  enrichment_df_for_plot[,2] <- -log10(enrichment_df_for_plot$qvalue) # -log10 the q val
  enrichment_df_for_plot <- enrichment_df_for_plot[order(enrichment_df_for_plot$qvalue),] # order again - maybe no need
  colnames(enrichment_df_for_plot) <- c("Enricher program", "-log10(qvalue)", "GeneRatio", "BgRatio") # change the colnames
  enrichment_df_for_plot$Collection <- scalop::substri(x = enrichment_df_for_plot$`Enricher program`, pos = 1) # add a collection column
  enrichment_df_for_plot$`Enricher program` <- scalop::substri(x = enrichment_df_for_plot$`Enricher program`, pos = -(1)) # remove the collection from the name
  
  enrichment_df_for_plot$Collection <-gsub(pattern = "C2", replacement = "Curated", x = enrichment_df_for_plot$Collection)
  enrichment_df_for_plot$Collection <-gsub(pattern = "C5", replacement = "Ontology", x = enrichment_df_for_plot$Collection)
  enrichment_df_for_plot$Collection <-gsub(pattern = "H", replacement = "Hallmark", x = enrichment_df_for_plot$Collection)
  enrichment_df_for_plot$Collection <-gsub(pattern = "C8", replacement = "Cell type", x = enrichment_df_for_plot$Collection)
  enrichment_df_for_plot$Collection <-gsub(pattern = "C1", replacement = "Chro. Pos.", x = enrichment_df_for_plot$Collection)
  enrichment_df_for_plot$Collection <-gsub(pattern = "C3", replacement = "Regulatory", x = enrichment_df_for_plot$Collection)
  enrichment_df_for_plot$Collection <-gsub(pattern = "C6", replacement = "Oncogenic", x = enrichment_df_for_plot$Collection)
  
  
  
  enricher_plot <- ggbarplot(enrichment_df_for_plot,
                             x = 'Enricher program', y = '-log10(qvalue)',
                             orientation = 'horiz', fill = "Collection",
                             label = paste(substri(enrichment_df_for_plot$GeneRatio, "/", 1), "/", substri(enrichment_df_for_plot$BgRatio, "/", 1)),
                             lab.hjust = -(0.3), lab.vjust = 0.6) +
    labs(title = paste("Enrichment plot for", title), subtitle = paste0(length(genes_to_test),"genes tested")) +
    scale_x_discrete(breaks=enrichment_df_for_plot$`Enricher program`, labels = scalop::substri(x = enrichment_df_for_plot$`Enricher program`, pos = 1:6)) 
  
  
  if(save_plot == T){
    write.csv(x = enrichment_df, file = paste0("/home/labs/tirosh/rotemta/Golan_PDXs/enrichments/", format(Sys.Date(),'%y%m%d'),"_", title, "_enrichment", ".csv"))
    ggsave(enricher_plot, file = paste0("/home/labs/tirosh/rotemta/Golan_PDXs/plots/", format(Sys.Date(),'%y%m%d'), "_emrichment_", title, ".png"), dpi = 300, type = "cairo", scale = 3)
  }
  
  if(print_plot){
    plot(enricher_plot)
  }
  
  
  return(enrichment_df)
}

firstup <- function(x) {
  x <- tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

##Differential Expression Analysis
dea <- function(x, g1, g2, name) {
  de_res <- tibble(Gene = rownames(x),
                   log2FC = rowMeans(x[, g1]) - rowMeans(x[, g2]),
                   p.val = sapply(1:nrow(x), function(i) wilcox.test(x = x[i, g1], y = x[i, g2])$p.value),
                   p.adj = p.adjust(p.val, "fdr"),
                   Name = name)
  de_res
}

outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}


# Import Data -------------------------------------------------------------

genes_ev2 <- read_delim("/home/EV2_features.tsv.gz", col_names = c("ENSEMBL_ID", "symbol", "type")) 
rmat_ev2 <- readMM("/home/EV2_matrix.mtx.gz") 
cells_ev2 <- read_delim("/home/EV2_barcodes.tsv.gz", delim = "\t", col_names = "cell_name")
cells_ev2$cell_name <- str_replace_all(cells_ev2$cell_name, "-", "_")
cells_ev2$cell_name <- paste0("ev2_", cells_ev2$cell_name)

colnames(rmat_ev2) <- cells_ev2$cell_name
rownames(rmat_ev2) <- genes_ev2$symbol

genes_ev3 <- read_delim("/home/EV3_features.tsv.gz", col_names = c("ENSEMBL_ID", "symbol", "type")) 
rmat_ev3 <- readMM("/home/EV3_matrix.mtx.gz") 
cells_ev3 <- read_delim("/home/EV3_barcodes.tsv.gz", delim = "\t", col_names = "cell_name")
cells_ev3$cell_name <- str_replace_all(cells_ev3$cell_name, "-", "_")
cells_ev3$cell_name <- paste0("ev3_", cells_ev3$cell_name)

colnames(rmat_ev3) <- cells_ev3$cell_name
rownames(rmat_ev3) <- genes_ev3$symbol

genes_mpn3 <- read_delim("/home/MPN3_features.tsv.gz", col_names = c("ENSEMBL_ID", "symbol", "type")) 
rmat_mpn3 <- readMM("/home/MPN3_matrix.mtx.gz") 
cells_mpn3 <- read_delim("/home/MPN3_barcodes.tsv.gz", delim = "\t", col_names = "cell_name")
cells_mpn3$cell_name <- str_replace_all(cells_mpn3$cell_name, "-", "_")
cells_mpn3$cell_name <- paste0("mpn3_", cells_mpn3$cell_name)

colnames(rmat_mpn3) <- cells_mpn3$cell_name
rownames(rmat_mpn3) <- genes_mpn3$symbol

genes_mpn5 <- read_delim("/home/MPN5_features.tsv.gz", col_names = c("ENSEMBL_ID", "symbol", "type")) 
rmat_mpn5 <- readMM("/home/MPN5_matrix.mtx.gz") 
cells_mpn5 <- read_delim("/home/MPN5_barcodes.tsv.gz", delim = "\t", col_names = "cell_name")
cells_mpn5$cell_name <- str_replace_all(cells_mpn5$cell_name, "-", "_")
cells_mpn5$cell_name <- paste0("mpn5_", cells_mpn5$cell_name)

colnames(rmat_mpn5) <- cells_mpn5$cell_name
rownames(rmat_mpn5) <- genes_mpn5$symbol

## Intersect genes and cbind mats
genes <- reduce(list(genes_ev2$symbol, genes_ev3$symbol, genes_mpn3$symbol, genes_mpn5$symbol), intersect)
umi_data <- cbind(rmat_ev2[genes, ], rmat_ev3[genes, ], rmat_mpn3[genes, ], rmat_mpn5[genes, ])
umi_data <- as.matrix(umi_data)
rmat <- (sweep(umi_data, 2, colSums(umi_data), "/"))*1e6  ### UMI count to CPM --> The data is TPM normalized, e.g. the sum of each column equals 10^6.
rmat_log2 <- log2(rmat/10+1) # log transform
rmat_centered <- rmat_log2 - rowMeans(rmat_log2) # centering

cells <- tibble(cell_name = colnames(rmat), mouse = if_else(str_detect(colnames(rmat), "^ev2_.."), true = "ev2", false = "NA"))
cells <- cells %>% 
  mutate(mouse = if_else(str_detect(colnames(rmat), "^ev3_.."), true = "ev3", false = mouse)) %>% 
  mutate(mouse = if_else(str_detect(colnames(rmat), "^mpn5_.."), true = "mpn5", false = mouse)) %>% 
  mutate(mouse = if_else(str_detect(colnames(rmat), "^mpn3_.."), true = "mpn3", false = mouse))
cells <- cells %>% 
  mutate(sample = if_else(mouse %in% c("ev2", "ev3"), true = "normal", false = "mpn"))
cells$complexity <- apply(rmat, 2, function (x) length(which(x !=0))) # add complexity to metadata

rm(cells_ev2)
rm(cells_ev3)
rm(cells_mpn3)
rm(cells_mpn5)
rm(genes_ev2)
rm(genes_ev3)
rm(genes_mpn3)
rm(genes_mpn5)
rm(rmat_ev2)
rm(rmat_ev3)
rm(rmat_mpn3)
rm(rmat_mpn5)


# QC ----------------------------------------------------------------------

## Plot complexity distribution and threshold
ggplot(cells, aes(x=complexity)) + 
  geom_histogram(binwidth = 25) + 
  geom_vline(xintercept = 600, color = "red") +
  xlim(0,9000) +
  ggtitle("Complexity")

cells <- cells[cells$complexity >= 600, ]
top_genes <- rownames(rmat)[log2(rowMeans(rmat[ ,cells$cell_name])+1) > 4]

fmat <- rmat_centered[top_genes, cells$cell_name]
dim(fmat)
range(fmat)

rm(rmat_log2)
rm(rmat)
gc()


# Seurat pipeline with Tirosh QC ------------------------------------------
# We create a Seurat object with the cells and genes that passed our QC procedure (thereby bypassing Seurat's QC).
# The following steps (until the cell type inference) are pretty standard and straight forward

sds <- CreateSeuratObject(counts = umi_data[top_genes, cells$cell_name],
                          project = "Zeiser",
                          min.cells = 3, min.features = 200)

rm(umi_data)
# This is Seurat's standard normalization method - returns a log-transformed (natural logarithm) counts-per-10K matrix
sds <-NormalizeData(sds,
                    normalization.method = "LogNormalize",
                    scale.factor = 10000)

# Computes the most highly variable genes which will be used for dimensionality redcution
sds <- FindVariableFeatures(sds, 
                            selection.method = "vst", 
                            nfeatures = 3000)

# Centers and scales the expression matrix
all.genes <- rownames(sds)
sds <- ScaleData(sds, 
                 features = all.genes)

npcs <- 100

# Computes PCA (basis for dimensionality reduction)
set.seed(100)
sds <- RunPCA(sds, 
              features = VariableFeatures(object = sds), 
              npcs = npcs)

# Cluster the data (Louvain clustering)
set.seed(100)
sds <- FindNeighbors(sds, 
                     dims = 1:npcs, 
                     graph.name = "snn")
set.seed(100)
sds <- FindClusters(sds, 
                    resolution = 0.95, 
                    graph.name = "snn")

# Compute UMAP
set.seed(100)
sds <- RunUMAP(sds, 
               dims = 1:npcs)

# Plot the UMAP result, by default project the clustering result
DimPlot(sds, 
        reduction = "umap", 
        label = T)

## Subclustering of cluster 17
sds <- FindSubCluster(sds, 
               cluster = c(15), 
               graph.name = "snn",
               subcluster.name = "snn_sub_15")


sds <- FindSubCluster(sds, 
               cluster = c(10), 
               graph.name = "snn",
               subcluster.name = "snn_sub_10", 
               resolution = 0.2)

sds <- FindSubCluster(sds, 
                      cluster = c(18), 
                      graph.name = "snn",
                      subcluster.name = "snn_sub_18", 
                      resolution = 0.2)


# Plotting cell characteristics -------------------------------------------
backup_cells <- cells

cells <- Embeddings(sds, reduction = "umap")
cells <- cbind(cells, FetchData(sds, vars = c("ident", "snn_sub_15", "snn_sub_10"))) %>% 
  as.data.frame() %>% 
  mutate(snn_sub = case_when(
    snn_sub_15 %in% c("10") ~ snn_sub_10,
    T ~ snn_sub_15
  ))

cells <- tibble(cell_name = rownames(cells), 
                       UMAP_1 = cells$UMAP_1, 
                       UMAP_2 = cells$UMAP_2, 
                       seurat_cluster = cells$snn_sub)
cells <- cells %>% left_join(backup_cells, by = "cell_name")

## PLotting
ggplot(cells, aes(x=UMAP_1, y=UMAP_2, color = cells$seurat_cluster)) + 
  geom_point(size = 0.3, alpha = 0.5) + 
  directlabels::geom_dl(aes(label = cells$seurat_cluster), method = "smart.grid") +
  guides(colour = guide_legend(override.aes = list(size=4, alpha = 1))) +
  labs(color = "cluster") +
  ggtitle("UMAP")

ggplot(cells, aes(x=UMAP_1, y=UMAP_2, color = cells$complexity)) + 
  geom_point(size = 0.3, alpha = 0.5) + 
  scale_color_viridis_c(option = "A", limits = c(0, 5000), oob = scales::squish) +
  labs(color = "cluster") +
  ggtitle("UMAP - Complexity")

ggplot(cells, aes(x=UMAP_1, y=UMAP_2, color = cells$mouse)) + 
  geom_point(size = 0.3, alpha = 0.5) + 
  guides(colour = guide_legend(override.aes = list(size=4, alpha = 1))) +
  labs(color = "Mouse") +
  ggtitle("UMAP - Sample")

ggplot(cells, aes(x=UMAP_1, y=UMAP_2, color = cells$sample)) + 
  geom_point(size = 0.3, alpha = 0.5) + 
  guides(colour = guide_legend(override.aes = list(size=15, alpha = 1))) +
  labs(color = "Condition") +
  ggtitle("UMAP - Condition") +
  theme_classic()

# DEGs --------------------------------------------------------------------

##Match fmat to cells
fmat <- fmat[ , cells$cell_name]

## DEGs
cells$seurat_cluster <- as.factor(cells$seurat_cluster)

DEGs <- tibble(1:100) # Create df in which the DEG for each cluster will be inserted

for (i in levels(cells$seurat_cluster)) {    # i iterates all clusters
  top_DE_ix <- sort( (rowMeans(fmat[, which(i == cells$seurat_cluster) ]) - rowMeans(fmat[,which(i != cells$seurat_cluster)]) ), decreasing = T, index.return = T )$ix[1:100]   # sort gene ix
  DEGs <- cbind(DEGs, rownames(fmat)[top_DE_ix[1:100]]  )   # select top gene names
}

DEGs <- DEGs[,-1] # Delete redundant first column
colnames(DEGs) <- levels(cells$seurat_cluster)

## Enricher using msigdb
DEGs_upper <- as_tibble(apply(DEGs, 2, str_to_upper))

# Most important categories from msigdb: H : hallmark, C2: curated, C5: GO annotations, C1: positional gene sets - by chromosome locus C8: cell type - specific signatures
msigdb_sigs = scalop::msigdb(category = c("C8"), species = "human")
# flatten list
msigdb_sigs = unlist(msigdb_sigs, recursive = F)
#### make universe
# Here universe is all possible genes, but you might want to define a smaller universe. Defining the universe is important for the hypergeometric test 
infercna::useGenome('hg19')
gene_hg19 <- infercna::retrieveGenome(name = 'hg19')
universe = gene_hg19$symbol
enrichments = scalop::enricher(test_gene_sets = DEGs_upper, 
                               ref_gene_sets = msigdb_sigs, 
                               universe = universe
                               )

names(enrichments) <- colnames(DEGs_upper)

## Load mouse bone marrow specific sigs
sig_stumpf <- read_delim("Stumpf_2020_sigs.csv")
sig_stumpf <- sig_stumpf %>% select(-Rank) %>% as.list(.)
names(sig_stumpf) <- paste0("stumpf_", names(sig_stumpf))

sig_harris <- read_delim("Harris_2021_sigs.csv")
sig_harris <- sig_harris %>% 
  filter(rank <=50) %>% 
  select(cell_type, gene) %>% 
  pivot_wider(names_from = cell_type, values_from = gene, values_fn = list) %>% 
  flatten(.)
names(sig_harris) <- paste0("harris_", names(sig_harris))
  
sig_lee <- map(readxl::excel_sheets("B_cell_dev_stages.xlsx"), ~readxl::read_excel("B_cell_dev_stages.xlsx", sheet = .))
names(sig_lee) <- c("lee_pre_pro_Bcell",
                    "lee_cycling_pro_Bcell",
                    "lee_pro_B_VDJ",
                    "lee_preBCRd",
                    "lee_preBCRi",
                    "lee_preBCRi_II_S_phase",
                    "lee_preBCRi_II_G2_M_phase",
                    "lee_kappa_pre_Bcell",
                    "lee_lambda_pre_Bcell",
                    "lee_immature_Bcell",
                    "lee_cycling_immature_Bcell",
                    "lee_mature_Bcell",
                    "lee_Pcell",
                    "lee_high_mito_Bcell")

sig_lee <- lapply(sig_lee, function(tbl) as.character(head(tbl[[1]], 100)))


sig_custom <- c(sig_harris, sig_stumpf, sig_lee)






## Enrichments for custom sigs. Cave: When running against mouse universe, Rotems modified enricher function needs to be used
## because Julie's is only for human

msigdb_sigs_mm = scalop::msigdb(species = "Mus musculus", category = c("C8") ) %>% 
  unlist(.,  recursive = F) %>% .[(sapply(., function(i) length(i) < 250))]

## Generate list of ref_sigs including custom sigs + HAY:
sig_HAY <- msigdb_sigs_mm %>% keep(str_detect(names(msigdb_sigs_mm), pattern = "HAY"))



enricher_func_mm(
  genes_to_test = DEGs$`18`,
  all_genes_within_matrix = rownames(fmat),
  ref_sigs_to_add = sig_custom,
  ref_sigs = sig_HAY,
  print_plot = T
)

## Generate heatmap for enrichments
tmp <- Julie_enricher_mm(test_gene_sets = DEGs, 
                         ref_gene_sets = sig_custom, 
                         universe = rownames(fmat), pvalueCutoff = 1, qvalueCutoff = 1
)
names(tmp) <- colnames(DEGs)

tmp <- do.call(rbind, tmp) %>% 
  as_tibble()

tmp2 <- tmp %>% select(sig = ID, qvalue, cluster = name) %>% 
  mutate(qvalue = -log10(qvalue)) %>%  
  pivot_wider(names_from = sig, values_from = qvalue)
  

tmp2 <- tmp2 %>% 
  column_to_rownames(., var = "cluster") %>% 
  mutate_all(~ replace_na(., 0))

Heatmap(tmp2,
        col = circlize::colorRamp2(breaks = seq(0,20,length.out = 9), 
                                   RColorBrewer::brewer.pal(9, "Reds")),
        name = "-log10 qvalue"
        )

## For Be cells only
Heatmap(tmp2[rownames(tmp2) %in% c("1", "0", "5", "15_2","13","14","11", "16","15_4"),],
        col = circlize::colorRamp2(breaks = seq(0,25,length.out = 9), 
                                   RColorBrewer::brewer.pal(9, "Reds")),
        name = "-log10 qvalue"
)

# Cluster Annotation ------------------------------------------------------
DimPlot(sds, 
        reduction = "umap", 
        label = T, 
        group.by = "snn_sub_15")

cells <- cells %>% 
  mutate(seurat_cluster = as.character(seurat_cluster)) %>% 
  mutate(cell_type = as.character(seurat_cluster)) %>% 
  mutate(cell_type = if_else(seurat_cluster %in% c(5, "15_2"), true = "CLP", false = seurat_cluster)) %>% 
  mutate(cell_type = if_else(seurat_cluster %in% c(13), true = "late_pro_BC", false = cell_type)) %>% 
  mutate(cell_type = if_else(seurat_cluster %in% c(14), true = "pre_BC_S_phase", false = cell_type)) %>% 
  mutate(cell_type = if_else(seurat_cluster %in% c(11), true = "pre_BC_G2M_phase", false = cell_type)) %>% 
  mutate(cell_type = if_else(seurat_cluster %in% (0), true = "pre_BC_small", false = cell_type)) %>% 
  mutate(cell_type = if_else(seurat_cluster %in% c(1), true = "im_BC", false = cell_type)) %>% 
  mutate(cell_type = if_else(seurat_cluster %in% c(16,"15_4"), true = "pro_BC", false = cell_type)) %>% #17 Mast?
  mutate(cell_type = if_else(seurat_cluster %in% c(18), true = "T_NK", false = cell_type)) %>% 
  mutate(cell_type = if_else(seurat_cluster %in% c(6), true = "pro_erythroblast", false = cell_type)) %>% 
  mutate(cell_type = if_else(seurat_cluster %in% c(8), true = "erythroblast", false = cell_type)) %>% 
  mutate(cell_type = if_else(seurat_cluster %in% c(2,3,17,4), true = "neutrophil", false = cell_type)) %>% 
  mutate(cell_type = if_else(seurat_cluster %in% c(7), true = "myelocyte", false = cell_type)) %>% 
  mutate(cell_type = if_else(seurat_cluster %in% c(12), true = "myeloblast", false = cell_type)) %>% 
  mutate(cell_type = if_else(seurat_cluster %in% c(10, "10_0"), true = "monoblast", false = cell_type)) %>% 
  mutate(cell_type = if_else(seurat_cluster %in% c("10_2"), true = "macrophage", false = cell_type)) %>% 
  mutate(cell_type = if_else(seurat_cluster %in% c(19), true = "basophil", false = cell_type)) %>% 
  mutate(cell_type = if_else(seurat_cluster %in% c(9,21,20, "10_1"), true = "monocyte", false = cell_type)) %>% 
  mutate(cell_type = if_else(seurat_cluster %in% c("15_0", "15_1"), true = "lq", false = cell_type)) %>% 
  #mutate(cell_type = if_else(seurat_cluster %in% c("15_2"), true = "excluded", false = cell_type)) %>% 
  mutate(cell_type = if_else(seurat_cluster %in% c("15_3"), true = "Pcell", false = cell_type))
  
sds$cell_type <- cells$cell_type
sds$condition <- cells$sample

## Plot
ggarrange(
  DimPlot(
    sds,
    reduction = "umap",
    label = T,
    group.by = "cell_type",
    cols = an_cols,
    repel = T,
    label.box = T,
    raster = T,
    pt.size = 0.7,
    raster.dpi = c(1000, 1000)
  )
  ,
  DimPlot(
    sds,
    reduction = "umap",
    label = F,
    group.by = "condition",
    cols = an_cols_condition,
    pt.size = 0.7,
    raster = T,
    raster.dpi = c(1000, 1000)
  ),
  ncol = 2
)

# DEGs cell types ---------------------------------------------------------

DEGs_ct <- tibble(1:50) # Create df in which the DEG for each cluster will be inserted
for (i in unique(cells$cell_type)) {    # i iterates all clusters
  top_DE_ix <- sort( (rowMeans(fmat[, which(i == cells$cell_type) ]) - rowMeans(fmat[,which(i != cells$cell_type)]) ), decreasing = T, index.return = T )$ix[1:50]   # sort gene ix
  DEGs_ct <- cbind(DEGs_ct, rownames(fmat)[top_DE_ix[1:50]]  )   # select top gene names
}
DEGs_ct <- DEGs_ct[,-1] # Delete redundant first column
colnames(DEGs_ct) <- unique(cells$cell_type)

## Generate heatmap for enrichments
tmp <- Julie_enricher_mm(test_gene_sets = DEGs_ct, 
                         ref_gene_sets = sig_custom, 
                         universe = rownames(fmat), pvalueCutoff = 1, qvalueCutoff = 1
)
names(tmp) <- colnames(DEGs_ct)

tmp <- do.call(rbind, tmp) %>% 
  as_tibble()

tmp2 <- tmp %>% select(sig = ID, qvalue, cluster = name) %>% 
  mutate(qvalue = -log10(qvalue)) %>%  
  pivot_wider(names_from = sig, values_from = qvalue)


tmp2 <- tmp2 %>% 
  column_to_rownames(., var = "cluster") %>% 
  mutate_all(~ replace_na(., 0))

Heatmap(tmp2[c("CLP",
               "pre_BC_G2M_phase",
               "late_pro_BC",
               "pre_BC_S_phase",
               "pro_BC",
               "pre_BC_small",
               "im_BC",
               "Pcell",
               "macrophage",
               "monocyte",
               "monoblast",
               "myeloblast",
               "basophil",
               "myelocyte",
               "neutrophil",
               "T_NK",
               "pro_erythroblast",
               "erythroblast",
               "lq"
               ),
             c("lee_preBCRi",
                "lee_preBCRi_II_G2_M_phase",
                "lee_preBCRi_II_S_phase",
                "lee_pro_B_VDJ",
                "lee_kappa_pre_Bcell", 
                "harris_immature B cell",
                "lee_Pcell",
                "harris_macrophage", 
                "harris_monocyte",
                "stumpf_Monoblasts",
                "stumpf_Myeloblasts",
                "harris_basophil",
                "stumpf_Myelocytes",
                "harris_granulocyte", 
                "harris_T cell",
                "harris_proerythroblast",
                "harris_erythroblast"
                )],
        col = circlize::colorRamp2(breaks = seq(0,20,length.out = 9), 
                                   RColorBrewer::brewer.pal(9, "Reds")),
        row_names_side = "left",
        row_dend_side = "right",
        name = "-log10 qvalue", cluster_columns = F, cluster_rows = F
)


universe = rownames(fmat)

# Calculate DEGs between conditions with statistics -----------------------
## run DE for each cluster against all other cells
de_genes <- lapply(levels(as.factor(cells$cell_type)), 
                   function(dbc) dea(fmat,
                                     cells$cell_name[cells$cell_type == dbc & cells$sample == "mpn"],
                                     cells$cell_name[cells$cell_type == dbc & cells$sample == "normal"],
                                     dbc))

## bind all tables into one table
de_genes <- do.call(rbind, de_genes)

## Significant genes in B cell cluster
de_genes %>% 
  filter(Name %in% c("im_BC", "CLP", "Pcell", "pre_BC_small", "pre_BC_G2M_phase", "pre_BC_S_phase", "pro_BC", "late_pro_BC")) %>% 
  filter(!(Gene %in% c("Hba-a2","Hba-a1", "Hbb-bt", "Hbb-bs"))) %>% 
  filter(p.adj < 0.1 | log2FC > 2 | log2FC < -2) %>% 
  filter(log2FC > 1 | log2FC < -1)



# Volcano plots -----------------------------------------------------------
list_volcano_tgfb <- vector("list")

for(i in levels(as.factor(cells$cell_type))) {
  
  list_volcano_tgfb[[i]] <-  ggplot(de_genes %>% 
                                      filter(!(Gene %in% c("Hba-a2","Hba-a1", "Hbb-bt", "Hbb-bs")) ) %>% 
                                      filter(Name == i), 
                                    aes(x = log2FC, y = -log10(p.adj), label = Gene)) +
    geom_point(alpha = 0.6, size = 2) +
    theme_classic() +
    labs(x="log2(FC)", y="-log10(FDR)") + 
    theme(axis.text = element_text(size=20), axis.title = element_text(size=20),
          plot.title = element_text(size=24), legend.text = element_text(size=16),
          legend.key=element_blank(), panel.background = element_rect(fill = "white"),
          axis.line = element_line(colour = "black"), aspect.ratio = 1) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", size = 1) + 
    geom_vline(xintercept = 1, linetype = "dashed", size = 1) +
    geom_vline(xintercept = -1, linetype = "dashed", size = 1) +
    ggrepel::geom_label_repel(data = subset(de_genes, Gene == "Tgfb1" & Name == i),
                              box.padding = 0.03, 
                              point.padding = 0.5,
                              segment.color = "red",
                              color = "red",
                              force_pull = 0.1,
                              nudge_y = 10) +
    ggtitle(i)
}

list_volcano_tgfb$pro_erythroblast
ggarrange(plotlist = list_volcano_tgfb[c(2:12)])

## Volcano plots B cell cluster only
list_volcano_b_cell <- vector("list")

for(i in c("im_BC", "CLP", "Pcell", "pre_BC_small", "pre_BC_G2M_phase", "pre_BC_S_phase", "pro_BC", "late_pro_BC")) {
  
  list_volcano_b_cell[[i]] <-  ggplot(de_genes %>% 
                                      filter(!(Gene %in% c("Hba-a2","Hba-a1", "Hbb-bt", "Hbb-bs")) ) %>% 
                                      filter(Name == i), 
                                    aes(x = log2FC, y = -log10(p.adj), label = Gene)) +
    geom_point(alpha = 0.6, size = 2) +
    theme_classic() +
    labs(x="log2(FC)", y="-log10(FDR)") + 
    theme(axis.text = element_text(size=20), axis.title = element_text(size=20),
          plot.title = element_text(size=24), legend.text = element_text(size=16),
          legend.key=element_blank(), panel.background = element_rect(fill = "white"),
          axis.line = element_line(colour = "black"), aspect.ratio = 1) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", size = 1) + 
    geom_vline(xintercept = 1, linetype = "dashed", size = 1) +
    geom_vline(xintercept = -1, linetype = "dashed", size = 1) +
    ggrepel::geom_label_repel(data = subset(de_genes %>% 
                                              filter(!(Gene %in% c("Hba-a2","Hba-a1", "Hbb-bt", "Hbb-bs")) ) %>% 
                                              filter(Name %in% c(i)), 
                                            p.adj < 0.05 & log2FC > 1 | p.adj < 0.05 & log2FC < -1),
                              box.padding = 0.03, 
                              point.padding = 0.5,
                              segment.color = "red",
                              color = "red",
                              force_pull = 0.1,
                              nudge_y = 10, max.overlaps = Inf, direction = "y") +
    ggtitle(i)
}

ggarrange(plotlist = list_volcano_b_cell[c(5:8)], ncol = 4)

subset(de_genes %>% 
         filter(!(Gene %in% c("Hba-a2","Hba-a1", "Hbb-bt", "Hbb-bs")) ) %>% 
         filter(Name %in% c("Pcell")), 
       log2FC > 2 | log2FC < -2)


# Cell type fractions -----------------------------------------------------
tmp <- cells %>% 
  group_by(sample, cell_type) %>% 
  count() %>% 
  ungroup() %>% 
  group_by(sample) %>% 
  mutate(frac = n / sum(n))

hotmap <- readRDS('/home/labs/tirosh/alissag/st_glioma/scripts_workspaces/hotmap.rds')

ggplot(tmp, aes(x=factor(sample, levels = c("normal", "mpn")), y=fct_reorder(cell_type, frac, .desc = F), fill=frac)) +
  geom_tile() +
  scale_fill_gradientn(colors = hotmap, 
                       limits = c(0,0.15),
                       oob = scales::squish,
                       guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
  theme_classic() +
  labs(x="condition", y="cell type", fill="fraction of cells")



