#### Jinxin Meng, 20240801, 20241113 ####
setwd('F:/project/20240731_BC_IBD_ILA_ZhangYN/Code-Data-available/git/Figure4/')
pacman::p_load(tidyverse, ggpubr)
source('/code/R_func/calcu_difference.R')
source('/code/R_func/profile_process.R')

#### Figure 4e ####
metabolite_info <- read.delim('MBX_Mus_metabolite_info.tsv')
profile <- read.delim('MBX_Mus_profile.tsv', row.names = 1)
group <- read.delim('MBX_Mus_group.tsv')

plsda <- ropls::opls(data.frame(t(profile), check.names = F), pull(group, group),
                     log10L = T, orthoI = 0, predI = NA)

mean_data <- profile_smp2grp(profile, group) %>% 
  rownames_to_column('cpd_id') %>% 
  mutate(FC = BC_3d / Veh_3d,
         log2FC = log2(FC))

results <- data.frame(VIP = plsda@vipVn) %>% 
  rownames_to_column('cpd_id') %>% 
  left_join(mean_data, ., by = 'cpd_id') %>% 
  mutate(enriched = ifelse(VIP > 1 & log2FC > 1, 'BC_3d', 
                           ifelse(VIP > 1 & log2FC < -1, 'Veh_3d', 'none')))
write.xlsx(results, 'difference.xlsx')

plot_data <- select(results, cpd_id, log2FC, VIP, enriched) %>% 
  mutate(enriched = factor(enriched, c('BC_3d', 'none', 'Veh_3d')),
         log2FC = ifelse(log2FC > 3, 3, log2FC),
         log2FC = ifelse(log2FC < -3, -3, log2FC)) %>% 
  left_join(metabolite_info %>% select(1:2, cpd_class), by = 'cpd_id') %>% 
  mutate(label = ifelse(enriched != 'none', cpd_name, ''))

ggscatter(plot_data, 'log2FC', 'VIP', color = 'enriched', legend = 'right',
          palette = c('#eeacec', 'grey', '#21c1dc'), size = 2,
          xlab = 'log2FoldChange', ylab = 'Variable Importance in Projection') +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  geom_vline(xintercept = c(-1, 1), linetype = 'dashed') +
  geom_text(aes(label = label), size = 1, fontface = 'italic') +
  theme(aspect.ratio = 1,
        axis.ticks.length = unit(2, 'mm'))
ggsave('Figure 4e.pdf', width = 5.5, height = 4)

#### Figure 4d ####
bg_data <- read.delim('../data/cpd2path_enrichment.tsv')

cpds <- filter(results, enriched == 'BC_3d') %>% 
  select(cpd_id) %>% 
  left_join(metabolite_info, by = 'cpd_id') %>% 
  filter(cpd_KEGG != '') %>% 
  pull(cpd_KEGG)

eKEGG <- clusterProfiler::enricher(gene = cpds, TERM2GENE = bg_data,
                                   minGSSize = 1, pvalueCutoff = 1, qvalueCutoff = 1) %>% 
  data.frame()

plot_data <- read.xlsx('eKEGG_plot.xlsx') %>% 
  dplyr::select(all_of(c('Description','qvalue'))) %>% 
  mutate(Description = reorder(Description, qvalue))

ggbarplot(plot_data, 'Description', 'qvalue', fill = 'qvalue', rotate = T, width = .6,
          xlab = '', legend = 'right', x.text.angle = 30) +
  scale_fill_gradient(low = '#abd9e9', high = '#fdae61') +
  theme(aspect.ratio = 2,
        axis.ticks.length = unit(2, 'mm'))
ggsave('Figure 4d.pdf', width = 8, height = 6)
