#### Jinxin Meng, 20250925, 20250925 ####
setwd('F:/project/20240731_BC_IBD_ILA_ZhangYN/Code-Data-available/git/Figure5/')
pacman::p_load(tidyverse, ggpubr, openxlsx)
source('/code/R_func/calcu_difference.R')

#### Figure 5b ####
profile <- read.xlsx('MBX_Trp_data.xlsx', sheet = 'profile', rowNames = T)
group <- read.xlsx('MBX_Trp_data.xlsx', sheet = 'group')

group_level <- c('BC', 'CMC')
group_color <- c('BC' = '#EFD496', 'CMC' = '#9BB89C')

plsda <- ropls::opls(t(log10(profile)), y = pull(group, group), crossvalI = 6, predI = 1)
difference <- difference_analysis(profile, group, comparison = group_level, method = 't')

result <- data.frame(vip = plsda@vipVn) %>% 
  rownames_to_column('name') %>% 
  left_join(difference, ., by = 'name') %>% 
  mutate(enriched = ifelse(vip > 1 & padj < 0.05 & log2FC > 0, 'BC', 
                           ifelse(vip > 1 & padj < 0.05 & log2FC < 0, 'CMC', 'none')))

difference <- filter(difference, grepl('^Indole', name)) %>% 
  arrange(desc(log2FC)) %>% 
  mutate(label = paste0('Fold=', round(FC, 3), ', FDR=', round(padj, 3)))

plot_data <- data.frame(t(profile[difference$name,]), check.names = F) %>% 
  rownames_to_column('sample') %>% 
  gather('name', 'value', -sample) %>% 
  left_join(group, by = 'sample') %>% 
  mutate(name = factor(name, difference$name))

map(difference$name, ~ filter(plot_data, name == .x) %>% 
      ggplot(aes(x = value, y = name, fill = group)) +
      ggridges::geom_density_ridges(rel_min_height = 0.001, scale = 100, alpha = .8) +
      labs(x = 'Metabolite Content (nmol/mg)', y = '', fill = '',
           title = 'Indole-related metabolites') +
      xlim(0, NA) +
      scale_fill_manual(values = group_color) +
      scale_y_discrete(expand = c(.01, 0), 
                       labels = paste0(.x, '\n', difference[difference$name == .x, 'label'])) +
      theme_pubr() +
      theme(aspect.ratio = 1/7,
            axis.ticks.length = unit(2, 'mm'),
            legend.position = 'right') ) %>% 
  cowplot::plot_grid(plotlist = ., ncol = 1, align = 'v')
ggsave('Figure 5b', width = 7, height = 7)

#### Figure 5k and S4g ####
data <- read.delim('minimap2.rc.tsv') %>% 
  filter(reads_count > 5000000) %>% 
  mutate(Bco_rela_ab = Bcoccoides / reads_count * 100,
         Csp_rela_ab = Csporogenes / reads_count * 100,
         Pru_rela_ab = Prussellii / reads_count * 100) %>% 
  dplyr::select(name, ends_with('rela_ab')) %>% 
  left_join(read.delim('project_info.tsv'), by = c('name' = 'sample')) %>% 
  filter(!name %in% c('SRR27217299','ERR1620278','ERR1620262','SRR27029799','SRR27029565',
                      'SRR27029396','SRR5650177','SRR5650125','SRR6468521','SRR6468501',
                      'SRR6468624','SRR6504857','SRR6468572','SRR10029190','SRR27029841',
                      'SRR6468596'))

project <- split(data$name, data$project_ID)

plots <- map2(project, names(project), ~ {
  .data <- filter(data, name %in% .x)
  ggscatter(.data, 'Pru_rela_ab', 'Bco_rela_ab', cor.coef = T, cor.method = 'spearman', 
            color = '#66c2a5', shape = 16, xlab = 'Relative abundance of P. russellii (%)', 
            ylab = 'Relative abundance of B. coccoides (%)', title = .y,
            add = 'reg.line') +
    theme(aspect.ratio = 1,
          plot.title = element_text(face = 'bold', hjust = .5),
          axis.ticks.length = unit(2, 'mm')) })
cowplot::plot_grid(plotlist = plots, align = 'v', nrow = 1)
ggsave('Figure 5k.pdf', width = 20, height = 4)

plots <- map2(project, names(project), ~ {
  .data <- filter(data, name %in% .x)
  ggscatter(.data, 'Bco_rela_ab', 'Csp_rela_ab', cor.coef = T, cor.method = 'spearman', 
            color = '#66c2a5', shape = 16, xlab = 'Relative abundance of B. coccoides (%)', 
            ylab = 'Relative abundance of C. sporogenes (%)', title = .y,
            add = 'reg.line') +
    theme(aspect.ratio = 1,
          plot.title = element_text(face = 'bold', hjust = .5),
          axis.ticks.length = unit(2, 'mm')) })
cowplot::plot_grid(plotlist = plots, align = 'v', nrow = 1)
ggsave('Figure S4g.pdf', width = 20, height = 4)
