#### Jinxin Meng, 20241027, 20250903 ####
setwd('F:/project/20240731_BC_IBD_ILA_ZhangYN/Code-Data-available/git/Figure6/')
pacman::p_load(tidyverse, ggpubr)
library(Biostrings)
library(rlang)

#### IPA code ####
msa <- Biostrings::readAAMultipleAlignment('key_fldH.msa.trim')%>% 
  as.character() %>% 
  as.list() %>% 
  map2_dfc(names(.), ~ data.frame(aa = unlist(str_split(.x, ''))) %>%
             dplyr::rename(!!.y := aa)) %>% 
  rownames_to_column('loc') %>% 
  select(1:3)

data <- map2(seq(0, 320, 80), seq(80, 400, 80), \(x, y)
     filter(msa, loc %in% seq(x + 1, y, 1)) %>% 
       gather(key = 'seq', value = 'aa', -loc) %>% 
       mutate(loc = factor(loc, seq(x + 1, y, 1)))) %>% 
  setNames(paste0('seq', 1:5))

aa_col <- c('E' = '#ff6d6d', 'D' = '#ff6d6d', 'P' = '#f2be3c', 'A' = '#f2be3c', 'V' = '#f2be3c',
            'M' = '#f2be3c', 'L' = '#f2be3c', 'I' = '#f2be3c', 'G' = '#f2be3c', 'K' = '#769dcc',
            'R' = '#769dcc', 'H' = '#769dcc', 'N' = '#74ce98', 'T' = '#74ce98', 'C' = '#74ce98',
            'Q' = '#74ce98', 'S' = '#74ce98', 'F' = '#ffff66', 'Y' = '#ffff66', 'W' = '#ffff66',
            " " = 'grey95')

map(data, \(x) 
    ggplot(x, aes(loc, seq, fill = aa)) +
  geom_tile(color = 'white', show.legend = F) +
  geom_text(aes(label = aa), size = 1.8) +
  labs(x = '', y = '') +
  scale_x_discrete(breaks = c(as.character(min(as.numeric(as.character(x$loc)))),
                              as.character(max(as.numeric(as.character(x$loc))))),
                   labels = c(as.character(min(as.numeric(as.character(x$loc)))),
                              as.character(max(as.numeric(as.character(x$loc)))))) +
  scale_fill_manual(values = aa_col) +
  coord_fixed() +
  theme_pubr() +
  theme(axis.line = element_blank(),
        axis.text = element_text(size = 5),
        axis.ticks = element_blank())
    ) %>% cowplot::plot_grid(plotlist = ., ncol = 1)

ggsave('Figure 6a', width = 10, height = 6)
