#### Jinxin Meng, 20240801, 20250917 ####
setwd('F:/project/20240731_BC_IBD_ILA_ZhangYN/Code-Data-available/git/Figure2/')
pacman::p_load(tidyverse, ggpubr, openxlsx)
source('/Code/R_func/profile_process.R')

#### Figure 2c ####
fpkm <- read.delim('RNASeq_fpkm.tsv', row.names = 1, check.names = F) %>% 
  dplyr::select(contains('7d'))

group <- read.delim('RNASeq_group.tsv') %>% 
  filter(sample %in% colnames(fpkm))
group_level <- c('Veh_7d', 'BC_7d')

markers <- list(
  aISC = c('Lgr5','Sox9','Ascl2'),
  rISC = c('Hopx','Clu','Bmi1'),
  PC = c('Ang4','Lyz1'),
  GC = c('Muc2','Ccl9'),
  TC = c('Dclk1','Cd24a','Trpm5'),
  EEC = c('Chga', 'Neurog3')
)

marker_level <- c('aISC', 'rISC', 'PC', 'GC', 'TC', 'EEC')

fpkm <- filter(fpkm, rownames(fpkm) %in% unlist(markers, use.names = F)) 

row_data <- map2_dfr(markers, names(markers), ~ data.frame(name = .x, class = .y)) %>% 
  mutate(name = factor(name, rownames(fpkm))) %>% 
  arrange(name) %>% 
  column_to_rownames('name') %>% 
  mutate(class = factor(class, marker_level))

col_data <- dplyr::select(group, sample, group) %>% 
  mutate(sample = factor(sample, colnames(fpkm))) %>% 
  arrange(sample) %>% 
  column_to_rownames('sample') %>% 
  mutate(group = factor(group, group_level))

row_split <- row_data$class
col_split <- col_data$group

pdf('Figure 2c.pdf', width = 4, height = 4)
ComplexHeatmap::pheatmap(
  fpkm, scale = 'row',
  color = colorRampPalette(c('#4196aa', 'white', '#d18b0f'))(100),
  split = row_split, column_split = col_split,
  cluster_row_slices = F, cluster_column_slices = F, 
  row_title_rot = 0,
  row_title_gp = grid::gpar(fontsize = 8), column_title_gp = grid::gpar(fontsize = 8),
  row_gap = unit(1, 'mm'), column_gap = unit(1, 'mm'),
  treeheight_col = 15, treeheight_row = 15,
  cellheight = 12, cellwidth = 12,
  border_color = 'white',
  border_gp = grid::gpar(col = 'black'),
  heatmap_legend_param = list(title = 'Scale FPKM'))
dev.off()

#### Figure 2b ####
group <- read.delim('RNASeq_group.tsv') %>% 
  filter(group %in% c('Veh_7d', 'BC_7d')) %>% 
  mutate(group = factor(group, c('Veh_7d', 'BC_7d')))

fpkm <- read.delim('RNASeq_fpkm.tsv', row.names = 1)

gene_info <- AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, rownames(fpkm), keytype = 'SYMBOL',
                                   columns = c('GENETYPE', 'ENTREZID')) %>% 
  filter(GENETYPE == 'protein-coding')

fpkm <- filter(fpkm, rownames(fpkm) %in% gene_info$SYMBOL) %>% 
  select(all_of(group$sample))
fpkm <- log2(fpkm[rowSums(fpkm > 0) == ncol(fpkm), ])

group_level <- group$group
design <- model.matrix(~ group_level)
colnames(design) <- gsub('group_level', '', colnames(design))

fit <- limma::lmFit(fpkm, design)
fit <- limma::eBayes(fit)

difference <- limma::topTable(fit, coef = 'BC_7d', number = Inf, adjust.method = 'BH') %>% 
  rownames_to_column('name') %>% 
  mutate(enriched = ifelse(logFC > 1 & P.Value < 0.05, 'BC_7d', 
                           ifelse(logFC < -1 & P.Value < 0.05, 'Veh_7d', 'none')))
write.xlsx(diff, 'D7_difference.xlsx')

eKEGG <- filter(difference, enriched == 'BC_7d') %>% 
  mutate(ENTREZID = gene_info$ENTREZID[match(name, gene_info$SYMBOL)]) %>% 
  pull(ENTREZID) %>% 
  clusterProfiler::enrichKEGG(organism = 'mmu', pvalueCutoff = 1, qvalueCutoff = 1) %>% 
  data.frame(row.names = NULL) %>% 
  mutate(geneName = map_vec(strsplit(geneID, '/'), ~ 
                              paste(gene_info$SYMBOL[match(.x, gene_info$ENTREZID)], 
                                    collapse = '/')), .after = geneID)
write.xlsx(eKEGG, 'D7.eKEGG.up_gene.xlsx')

eGO <- filter(difference, enriched == 'BC_7d') %>% 
  pull(name) %>% 
  clusterProfiler::enrichGO(org.Mm.eg.db::org.Mm.eg.db, keyType = 'SYMBOL', 
                            pvalueCutoff = 1, qvalueCutoff = 1, ont = 'ALL') %>% 
  data.frame(row.names = NULL)
write.xlsx(eGO, 'D7.eGO.up_gene.xlsx')

plot_data <- read.xlsx('D7.enrichment.xlsx') %>% 
  dplyr::select(name = Description, FoldEnrichment, qvalue) %>% 
  arrange(desc(FoldEnrichment)) %>% 
  mutate(name = factor(name, .$name))

ggscatter(plot_data, 'FoldEnrichment', 'name', fill = 'qvalue', size = 'FoldEnrichment',
          shape = 21, legend = 'right', ylab = '', 
          title = 'GO and KEGG pathway enrichment analyses for up-regulated\ngenes in D7 post-DSS administration') +
  scale_size_continuous(range = c(4, 6)) +
  scale_fill_viridis_c(direction = 1, begin = .5) +
  theme(axis.ticks.length = unit(2, 'mm'),
        plot.title = element_text(face = 'bold', vjust = .5),
        panel.grid.major = element_line(linewidth = .4, color = 'grey90'))
ggsave('Figure 1a.pdf', width = 8, height = 4)

#### Figure 2d ####
fpkm <- read.delim('RNASeq_fpkm.tsv', row.names = 1, check.names = F) %>% 
  dplyr::select(contains('0d'))

group <- read.delim('RNASeq_group.tsv') %>% 
  filter(sample %in% colnames(fpkm))
group_level <- c('Veh_0d','BC_0d')

markers <- list(
  aISC = c('Lgr5','Sox9','Ascl2'),
  rISC = c('Hopx','Bmi1','Clu'),
  PC = c('Ang4','Lyz1'),
  GC = c('Muc2','Ccl9'),
  TC = c('Dclk1','Cd24a','Trpm5'),
  EEC = c('Chga', 'Neurog3')
)

marker_level <- c('aISC', 'rISC', 'PC', 'GC', 'TC', 'EEC')

fpkm <- filter(fpkm, rownames(fpkm) %in% unlist(markers, use.names = F)) 

row_data <- map2_dfr(markers, names(markers), ~ data.frame(name = .x, class = .y)) %>% 
  mutate(name = factor(name, rownames(fpkm))) %>% 
  arrange(name) %>% 
  column_to_rownames('name') %>% 
  mutate(class = factor(class, marker_level))

col_data <- dplyr::select(group, sample, group) %>% 
  mutate(sample = factor(sample, colnames(fpkm))) %>% 
  arrange(sample) %>% 
  column_to_rownames('sample') %>% 
  mutate(group = factor(group, group_level))

row_split <- row_data$class
col_split <- col_data$group

pdf('Figure 2d.pdf', width = 4.5, height = 4)
ComplexHeatmap::pheatmap(
  fpkm, scale = 'row',
  color = colorRampPalette(c('#4196aa', 'white', '#d18b0f'))(100),
  split = row_split, column_split = col_split,
  cluster_row_slices = F, cluster_column_slices = F, 
  row_title_rot = 0,
  row_title_gp = grid::gpar(fontsize = 8), column_title_gp = grid::gpar(fontsize = 8),
  row_gap = unit(1, 'mm'), column_gap = unit(1, 'mm'),
  treeheight_col = 15, treeheight_row = 15,
  cellheight = 12, cellwidth = 12,
  border_color = 'white',
  border_gp = grid::gpar(col = 'black'),
  heatmap_legend_param = list(title = 'Scale FPKM'))
dev.off()

#### Figure 2e ####
fpkm <- read.delim('RNASeq_fpkm.tsv', row.names = 1, check.names = F) %>% 
  dplyr::select(contains('3d'))

group <- read.delim('RNASeq_group.tsv') %>% 
  filter(sample %in% colnames(fpkm))
group_level <- c('Veh_3d', 'BC_3d')

markers <- list(
  aISC = c('Lgr5','Sox9','Ascl2'),
  rISC = c('Hopx','Clu','Bmi1'),
  PC = c('Ang4','Lyz1'),
  GC = c('Muc2','Ccl9'),
  TC = c('Dclk1','Cd24a','Trpm5'),
  EEC = c('Chga', 'Neurog3')
)

marker_level <- c('aISC', 'rISC', 'PC', 'GC', 'TC', 'EEC')

fpkm <- filter(fpkm, rownames(fpkm) %in% unlist(markers, use.names = F)) 

row_data <- map2_dfr(markers, names(markers), ~ data.frame(name = .x, class = .y)) %>% 
  mutate(name = factor(name, rownames(fpkm))) %>% 
  arrange(name) %>% 
  column_to_rownames('name') %>% 
  mutate(class = factor(class, marker_level))

col_data <- dplyr::select(group, sample, group) %>% 
  mutate(sample = factor(sample, colnames(fpkm))) %>% 
  arrange(sample) %>% 
  column_to_rownames('sample') %>% 
  mutate(group = factor(group, group_level))

row_split <- row_data$class
col_split <- col_data$group

pdf('Figure 2e.pdf', width = 5, height = 4)
ComplexHeatmap::pheatmap(
  fpkm, scale = 'row',
  color = colorRampPalette(c('#4196aa', 'white', '#d18b0f'))(100),
  split = row_split, column_split = col_split,
  cluster_row_slices = F, cluster_column_slices = F, 
  row_title_rot = 0,
  row_title_gp = grid::gpar(fontsize = 8), column_title_gp = grid::gpar(fontsize = 8),
  row_gap = unit(1, 'mm'), column_gap = unit(1, 'mm'),
  treeheight_col = 15, treeheight_row = 15,
  cellheight = 12, cellwidth = 12,
  border_color = 'white',
  border_gp = grid::gpar(col = 'black'),
  heatmap_legend_param = list(title = 'Scale FPKM'))
dev.off()
