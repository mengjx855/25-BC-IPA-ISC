#### Jinxin Meng, 20240125, 20260121 ####
pacman::p_load(tidyverse, ggpubr)
setwd('F:/project/20240731_BC_IBD_ILA_ZhangYN/Code-Data-available/git/Figure1/')
source('/code/R_func/calcu_difference.R')
source('/code/R_func/profile_process.R')

#### Figure 1a ####
group_level <- c('CD', 'Control')
group_color <- structure(c('#E69F00','#69a7bc'), names = group_level)
  
profile <- read.delim('16S_profile_g.txt', row.names = 1, check.names = F) %>% 
  profile_transRA()
group <- read.delim('16S_group.txt')

plot_data <- profile %>% 
  filter(rownames(.) %in% 'g__Blautia') %>% 
  t() %>% 
  data.frame() %>% 
  dplyr::rename(value = 1) %>% 
  rownames_to_column('sample') %>% 
  left_join(group, by = 'sample') %>%
  mutate(group = factor(group, c('Control', 'CD')))

ggviolin(plot_data, 'group', 'value', color = 'group', legend = 'none', xlab = '', 
         size = 1, width = .6, ylab = 'Relative abundance of Blautia (%)', 
         palette = group_color, outlier.shape = NA) +
  geom_boxplot(aes(color = group), width = .1, size = 1, outlier.shape = NA) + 
  lims(y = c(NA, 30)) +
  geom_jitter(aes(color = group), width = .2) +
  stat_compare_means(comparisons = list(c('CD', 'Control')), label.y = 27, 
                     tip.length = .01, vjust = .2, size = 3) +
  theme(aspect.ratio = 1.5)
ggsave('Figure 1a.pdf', width = 4, height = 4)

#### Figure S1a ####
source('/code/R_func/model_randomforest.R')
source('/code/R_func/plot_roc.R')

profile <- read.delim('16S_profile_g.txt', row.names = 1, check.names = F) %>% 
  profile_transRA()
group <- read.delim('16S_group.txt')

pred <- rf_loom(profile, group, seed = 2024, ntree = 1000)
pred <- left_join(pred, group, by = 'sample')
roc <- roc(pred$group, pred$CD)
plot_roc(roc, plot_se = T)
ggsave('Figure S1a-1.pdf', width = 4, height = 4)

difference <- difference_analysis(profile, group, comparison = c('CD', 'Control')) %>% 
  mutate(enriched = ifelse(avg_ab1 > avg_ab2, 'CD', 'Control'))

data <- rf_importance(profile, group, seed = 2024, ntree = 1000)

left_join(data, difference, by = 'name') %>% 
  head(10) %>% 
  ggbarplot('name', 'MeanDecreaseAccuracy', fill = 'enriched', rotate = T, sort.val = 'asc',
            sort.by.groups = F, width = .6, palette = group_color, xlab = '') +
  theme(aspect.ratio = 1.4)
ggsave('Figure S1a-2.pdf', width = 5.5, height = 5)

#### Figure 1b ####
# microeco是基于R6class开发的
library(ggtree)
library(microeco)

set.seed(2024)

profile <- read.delim('16S_table_rarefied.tsv', row.names = 1, check.names = F) %>%
  profile_transRA()
group <- read.delim('16S_group.txt', row.names = 1)
tax <- read.delim('16S_taxonomy.tsv', row.names = 1) %>% microeco::tidy_taxonomy()
tr <- treeio::read.tree('16S_tree.nwk')

dataset <- microeco::microtable$new(
  sample_table = group, otu_table = profile, 
  tax_table = tax, phylo_tree = tr)

lefse <- microeco::trans_diff$new(
  dataset = dataset, method = 'lefse', 
  group = 'group', alpha = 0.05, lefse_subgroup = NULL)

lefse$plot_diff_bar(use_number = 1:30, width = 0.6, group_order = c('Control', 'CD'), 
                    color_values = c('#69a7bc','#E69F00'))
ggsave('Figure 1b-1', width = 8, height = 10)

lefse$plot_diff_cladogram(use_taxa_num = 100, use_feature_num = 30, clade_label_level = 4, 
                          group_order = c('Control', 'CD'), color = c('#69a7bc','#E69F00'))
ggsave('Figure 1b-2.pdf', width = 10, height = 8)

#### Figure 1c LloydPriceJ_2019 ####
group <- read.xlsx('LloydPriceJ_2019.xlsx', sheet = 'group') %>% 
  dplyr::select(sample = external, group) %>% 
  mutate(sample = as.character(sample))

profile <- read.xlsx('LloydPriceJ_2019.xlsx', sheet = 'profile', rowNames = T)

taxonomy <- data.frame(name = rownames(profile), taxonomy = profile$taxonomy) %>% 
  mutate(genus = str_split_i(taxonomy, ';', 6),
         genus = ifelse(grepl('__\\w$', genus), 'Unknown', gsub('^\\s+[_]+', '', genus)))

profile <- profile[, -ncol(profile)]

data <- profile_filter(profile, dplyr::select(group, sample, group), 
                       by_group = T, n_group = 1, min_n = 3) %>% 
  profile_transRA()

features <- taxonomy %>% 
  filter(genus == 'Blautia') %>% 
  pull(name)

plot_data <- data %>% 
  filter(rownames(.) %in% features) %>% 
  colSums() %>% 
  data.frame(value = .) %>% 
  rownames_to_column('sample') %>% 
  left_join(group, by = 'sample') %>% 
  mutate(group = factor(group, c('nonIBD', 'CD', 'UC')))

ggviolin(plot_data, 'group', 'value', color = 'group', legend = 'none', xlab = '', size = 1, width = .5,
         ylab = 'Relative abundance of Blautia (%)', palette = c('#69a7bc','#E69F00','#f08178'),
         outlier.shape = NA, title = 'LloydPriceJ_2019') +
  geom_boxplot(aes(color = group), width = .18, size = 1, outlier.shape = NA) + 
  geom_jitter(aes(color = group), width = .2) +
  lims(y = c(NA, 10.5)) +
  stat_compare_means(comparisons = list(c('nonIBD', 'CD'), c('nonIBD', 'UC'),c('CD', 'UC')), 
                     label.y = 8.5, tip.length = .01, step.increase = .022, vjust = .9, size = 3) +
  theme(aspect.ratio = 1.5,
        axis.ticks.length = unit(2, 'mm'),
        plot.title = element_text(color = 'black', face = 'bold', hjust = .5))
ggsave('Figure 1c.pdf', width = 3.5, height = 4.5)

rstatix::wilcox_test(plot_data, value ~ group)

#### Figure 1d SchirmerM_2018 ####
group <- read.xlsx('SchirmerM_2018.PRJNA389280.xlsx', sheet = 'group') %>% 
  select(sample = run, group = group2)
data <- read.xlsx('SchirmerM_2018.PRJNA389280.xlsx', sheet = 'profile') %>% 
  mutate(rela_ab = rc / libs * 100) %>% 
  left_join(group, by = 'sample') %>% 
  mutate(group = factor(group, c('CT', 'CD', 'UC')))

ggviolin(data, 'group', 'rela_ab', color = 'group', legend = 'none', xlab = '', size = 1, width = .5,
         ylab = 'Relative abundance of B.coccoides (%)',  palette = c('#69a7bc','#E69F00','#f08178'),
         outlier.shape = NA, title = 'SchirmerM_2018') +
  geom_boxplot(aes(color = group), width = .2, size = 1, outlier.shape = NA) + 
  geom_jitter(aes(color = group), width = .2) +
  lims(y = c(NA, .88)) +
  stat_compare_means(comparisons = combn(c('CT', 'CD', 'UC'), 2, simplify = F), 
                     label.y = .75, tip.length = .02,  step.increase = .07,
                     vjust = .1, size = 3) +
  theme(aspect.ratio = 1.5,
        axis.ticks.length = unit(2, 'mm'),
        plot.title = element_text(color = 'black', face = 'bold', hjust = .5))
ggsave('Figure 1d.pdf', width = 3.5, height = 4.5)

rstatix::wilcox_test(data, rela_ab ~ group)

#### Figure 1e HeQ_2017 ####
group <- read.xlsx('HeQ_2017.PRJEB15371.xlsx', sheet = 'group')
data <- read.xlsx('HeQ_2017.PRJEB15371.xlsx', sheet = 'profile') %>% 
  mutate(rela_ab = rc / libs * 100) %>% 
  left_join(group, by = 'sample') %>% 
  mutate(group = factor(group, c('Control', 'CD')))

ggviolin(data, 'group', 'rela_ab', color = 'group', legend = 'none', xlab = '', size = 1, width = .5,
         ylab = 'Relative abundance of B.coccoides (%)', palette = c('#69a7bc','#E69F00','#f08178'),
         outlier.shape = NA, title = 'HeQ_2017') +
  geom_boxplot(aes(color = group), width = .2, size = 1, outlier.shape = NA) + 
  geom_jitter(aes(color = group), width = .2) +
  lims(y = c(NA, .8)) +
  stat_compare_means(comparisons = list(c('CD','Control')), label.y = .65,
                     tip.length = .004, step.increase = .042, vjust = .9, size = 3) +
  theme(aspect.ratio = 1.8,
        axis.ticks.length = unit(2, 'mm'),
        plot.title = element_text(color = 'black', face = 'bold', hjust = .5))
ggsave('Figure 1e.pdf', width = 3, height = 4.5)

rstatix::wilcox_test(data, rela_ab ~ group)

#### Figure 1f WengY_2019 ####
group <- read.xlsx('WengY_2019.PRJNA429990.xlsx', sheet = 'group')
data <- read.xlsx('WengY_2019.PRJNA429990.xlsx', sheet = 'profile') %>% 
  mutate(rela_ab = rc / libs * 100) %>% 
  left_join(group, by = 'sample') %>% 
  mutate(group = factor(group, c('Control', 'CD', 'UC')))

ggviolin(data, 'group', 'rela_ab', color = 'group', legend = 'none', xlab = '', size = 1, width = .5,
         ylab = 'Relative abundance of B.coccoides (%)',  palette = c('#69a7bc','#E69F00','#f08178'),
         outlier.shape = NA, title = 'WengY_2019') +
  geom_boxplot(aes(color = group), width = .2, size = 1, outlier.shape = NA) + 
  geom_jitter(aes(color = group), width = .2) +
  lims(y = c(NA, .58)) +
  stat_compare_means(comparisons = list(c('CD', 'Control'), c('UC', 'Control'), c('CD', 'UC')),
                     label.y = .38, tip.length = .003, step.increase = .01, vjust = .9, size = 3) +
  theme(aspect.ratio = 1.5,
        axis.ticks.length = unit(2, 'mm'),
        plot.title = element_text(color = 'black', face = 'bold', hjust = .5))
ggsave('Figure 1f.pdf', width = 3.5, height = 4.5)

rstatix::wilcox_test(data, rela_ab ~ group)

#### Figure 1g YanQ_2023c ####
group <- read.xlsx('YanQ_2023c.PRJEB67456.xlsx', sheet = 'group') %>% 
  select(sample = run, group)
data <- read.xlsx('YanQ_2023c.PRJEB67456.xlsx', sheet = 'profile') %>% 
  mutate(rela_ab = rc / libs * 100) %>% 
  left_join(group, by = 'sample') %>% 
  mutate(group = factor(group, c('healthy', 'CD', 'UC')))

ggviolin(data, 'group', 'rela_ab', color = 'group', legend = 'none', xlab = '', size = 1, width = .5,
         ylab = 'Relative abundance of B.coccoides (%)', palette = c('#69a7bc','#E69F00','#f08178'),
         outlier.shape = NA, title = 'YanQ_2023') +
  geom_boxplot(aes(color = group), width = .2, size = 1, outlier.shape = NA) + 
  geom_jitter(aes(color = group), width = .2) +
  lims(y = c(NA, 1.61)) +
  stat_compare_means(comparisons = list(c('healthy', 'CD'), c('healthy', 'UC'), c('UC', 'CD')), 
                     label.y = 1.35, tip.length = .01, step.increase = .042, vjust = .9, size = 3) +
  theme(aspect.ratio = 1.5,
        axis.ticks.length = unit(2, 'mm'),
        plot.title = element_text(color = 'black', face = 'bold', hjust = .5))
ggsave('Figure 1g.pdf', width = 3.5, height = 4.5)

rstatix::wilcox_test(data, rela_ab ~ group)
