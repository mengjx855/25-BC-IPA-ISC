#### Jinxin Meng, 20240308, 20250729 v0.3.3 ####

# 20231101: add all_group parameter in profile_filter.
# 20231219: add function profile_replace.
# 20250227: add function profile_top_n/frac.
# 20250605: add function profile_aggregate.
# 20250729: add function CLR(), profile_transCLR().
# 20250823: add digits parameter in profile_trans* functions

library(tidyverse)

#### LOG ####
# LOG transformation method in MaAsLin2;
# The default log transformation incorporated into MaAsLin does add a pseudo-count;
# As is best-known practice currently, the pseudo-count is half the minimum feature; 
# x [numeric]: a numeric vector.
LOG2 <- function(x, pseudocount = 1e-6, use_half_minimum = TRUE) {
  if (isFALSE(use_half_minimum)) x = replace(x, x == 0, pseudocount)
  if (isTRUE(use_half_minimum)) x = replace(x, x == 0, min(x[x > 0]) / 2)
  log2(x)
}

LOG10 <- function(x, pseudocount = 1e-6, use_half_minimum = TRUE) {
  if (isFALSE(use_half_minimum)) x = replace(x, x == 0, pseudocount)
  if (isTRUE(use_half_minimum)) x = replace(x, x == 0, min(x[x > 0]) / 2)
  log10(x)
}

CLR <- function(x, pseudocount = 1e-6, use_half_minimum = FALSE) {
  if (isFALSE(use_half_minimum)) x = replace(x, x == 0, pseudocount)
  if (isTRUE(use_half_minimum)) x = replace(x, x == 0, min(x[x > 0]) / 2)
  gm <- exp(mean(log(x))) # 计算每行的几何均值
  log(x) - log(gm) # 执行CLR转换：log(x) - log(几何均值)
}

#### profile_transXXX ####
# conduct LOG transformation for a otu_table, rows represent features, 
# and columns represent samples.
profile_transLOG2 <- function(profile, pseudocount = 1e-6, use_half_minimum = FALSE) {
  apply(profile, 1, \(x) LOG2(x, pseudocount = pseudocount, use_half_minimum = use_half_minimum)) %>% 
    t() %>% 
    data.frame(check.names = F)
  }
  
profile_transLOG10 <- function(profile, pseudocount = 1e-6, use_half_minimum = FALSE) {
  apply(profile, 1, \(x) LOG10(x, pseudocount = pseudocount, use_half_minimum = use_half_minimum)) %>% 
    t() %>% 
    data.frame(check.names = F)
  }

profile_transCLR <- function(profile, pseudocount = 1e-6, use_half_minimum = FALSE) {
  apply(profile, 2, \(x) CLR(x, pseudocount = pseudocount, use_half_minimum = use_half_minimum)) %>% 
    data.frame(check.names = F)
 }

profile_transSqrt <- function(profile) {
  sqrt(profile)
}

#### profile_transRA ####
# Jinxin Meng, 20231119
# replace some value meeting parameter to specified value.
profile_transRA <- function(profile, base = 100, digits = 8) {
  apply(data.frame(profile, check.names = F), 2, \(x) 
        round(x/sum(na.omit(x))*base, digits = digits)) %>%
    data.frame(check.names = F)
  }
  
#### profile_filter ####
# profile: input a data.frame of relative abundance profile.
# group: mapping (sample|group), also specify by map_names parameter.
# by_group: filter feature in each group.
# all_group prevalence only meet all group are outputted.
# min_prevalence: threshold of prevalence of features in all sample.
# min_abundance: threshold of abundance in a sample is considered a feature presenting in the sample.
profile_filter <- function(profile, group, group_rename = NULL, 
                           by_group = F, all_group = F, n_group = 1, 
                           min_prevalence = 0.1, min_n = NULL,
                           min_abundance = 0.0) {
  profile <- data.frame(profile, check.names = F)
  
  if (!is.null(min_n) & is.numeric(min_n)) 
    min_prevalence <- NULL
  
  if (isTRUE(all_group)) 
    n_group <- NULL
 
  # 按照分组过滤，每个组内计算流行率或者计数，然后进行过滤
  if (isTRUE(by_group) & !missing(group)) { 
    if (!all(c('sample', 'group') %in% colnames(group)) & is.null(group_rename)) 
      stop('group field (sample|group)')
    
    if (!is.null(group_rename)) 
      group <- data.frame(group, check.names = F) %>% 
        dplyr::rename(all_of(group_rename))
    
    if (!is.null(min_prevalence) & is.numeric(min_prevalence)) {
      prevalence <- profile %>% 
        t() %>% 
        data.frame(check.names = F) %>% 
        mutate(group = group$group[match(rownames(.), group$sample)]) %>% 
        group_by(group) %>% 
        group_modify(~ purrr::map_df(.x, \(x) sum(x > min_abundance)/length(x))) %>%
        ungroup() %>% 
        dplyr::select(-group)
      
      if (isTRUE(all_group))
        .vec <- map_vec(prevalence, \(x) all(x > min_prevalence))
      
      if (isFALSE(all_group))
        .vec <- map_vec(prevalence, \(x) sum(x > min_prevalence) >= n_group)
      
    } else if (!is.null(min_n) & is.numeric(min_n)) {
      
      count <- profile %>% 
        t() %>% 
        data.frame(check.names = F) %>% 
        mutate(group = group$group[match(rownames(.), group$sample)]) %>% 
        group_by(group) %>% 
        group_modify(~ purrr::map_df(.x, \(x) sum(x > min_abundance))) %>%
        ungroup() %>% 
        dplyr::select(-group)
      
      if (isTRUE(all_group))
        .vec <- map_vec(count, \(x) all(x >= min_n))
        
      if (isFALSE(all_group))
        .vec <- map_vec(count, \(x) sum(x >= min_n) >= n_group)
    }
    
    .vec <- names(.vec[.vec])
    profile <- profile[.vec, ]
  } 
  
  # 不按照分组。计算整体的流行率
  if (!isTRUE(by_group)) {
    if (!is.null(min_prevalence) & is.numeric(min_prevalence))
      .vec <- data.frame(t(profile), check.names = F) %>% 
        purrr::map_vec(\(x) sum(x > min_abundance)/length(x) > min_prevalence)
    
    if (!is.null(min_n) & is.numeric(min_n))
      .vec <- data.frame(t(profile), check.names = F) %>% 
        purrr::map_vec(\(x) sum(x > min_abundance) >= min_n)
    
    .vec <- names(.vec[.vec])
    profile <- profile[.vec,]
    
  }
  return(profile)
}
#### profile_smp2grp ####
# Jinxin Meng, 20231118
# method: merge column method. please get help in summarise_all() function.
profile_smp2grp <- function(profile, group, group_rename = NULL, method = 'mean'){
  if (!all(c('sample', 'group') %in% colnames(group)) & is.null(group_rename)) 
    stop('group field (sample|group)')
  
  if (!is.null(group_rename)) 
    group <- data.frame(group, check.names = F) %>% 
      dplyr::rename(all_of(group_rename))
  
  profile <- data.frame(t(profile), check.names = F) %>% 
    mutate(group = group$group[match(rownames(.), group$sample)]) %>% 
    aggregate(. ~ group, ., method) %>% 
    ungroup() %>% 
    column_to_rownames('group') %>%
    t() %>% 
    data.frame(check.names = F)
  return(profile)
}


#### profile_smp2grp_1 ####
# Jinxin Meng, 20240814, for profile containing NA or NAN .....
profile_smp2grp_1 <- function(profile, group, group_rename = NULL, method = "mean"){
  if (!all(c("sample", "group") %in% colnames(group)) & is.null(group_rename)) 
    stop("group field (sample|group)")
  if (!is.null(group_rename)) group <- data.frame(group, check.names = F) %>% 
      dplyr::rename(all_of(group_rename))
  data <- data.frame(t(profile), check.names = F) %>% 
    mutate(group = group$group[match(rownames(.), group$sample)])
  
  if (method == "mean") out <- data %>% 
    group_by(group) %>% 
    group_modify(~purrr::map_df(.x, \(x) ifelse(sum(!is.na(x)) == 0, 0, mean(na.omit(x))))) %>% 
    ungroup()

  if (method == "median") out <- data %>% 
      group_by(group) %>% 
      group_modify(~purrr::map_df(.x, \(x) ifelse(sum(!is.na(x)) == 0, 0, median(na.omit(x))))) %>% 
      ungroup()
  out <- out %>% 
    column_to_rownames("group") %>% 
    t() %>% 
    data.frame(check.names = F) %>% 
    rownames_to_column("name")
  return(out)
}

#### profile_top_n ####
profile_top_n <- function(profile, n = 12, out_other = F, other_name = 'Other', 
                          sort_method = 'mean') {
  .feats <- apply(profile, 1, sort_method) %>% 
    sort(decreasing = T) %>% 
    head(n = n) %>% 
    names()

  if (isFALSE(out_other)) profile <- filter(profile, rownames(profile) %in% .feats)
  
  if(isTRUE(out_other)) {
    profile <- profile %>% 
      mutate(name = rownames(.),
             name = replace(name, !name %in% .feats, other_name)) %>% 
      aggregate(. ~ name, ., sum) %>% 
      column_to_rownames('name')
    }
  
  return(profile)
  }

#### profile_top_frac ####
profile_top_frac <- function(profile, frac = .1, out_other = F,
                             other_name = 'Other', sort_method = 'mean') {
  .feats <- apply(profile, 1, sort_method) %>% 
    sort(decreasing = T) %>% 
    head(n = floor(nrow(profile) * frac)) %>% 
    names()
  
  if (isFALSE(out_other)) profile <- filter(profile, rownames(profile) %in% .feats)
  
  if(isTRUE(out_other)) {
    profile <- profile %>% 
      mutate(name = rownames(.),
             name = replace(name, !name %in% .feats, other_name)) %>% 
      aggregate(. ~ name, ., sum) %>% 
      column_to_rownames('name')
    }
  
  return(profile)
  }

#### profile_replace ####
# Jinxin Meng, 20231119
# replace some value meeting parameter to specified value.
profile_replace <- function(profile, min_value = 1, fill_value = 0, transRA = F){
  profile <- data.frame(profile, check.names = F)
  
  if (isTRUE(transRA)) 
    profile <- apply(profile, 2, x/sum(x) * 100) %>% 
      data.frame(check.names = F)
  
  profile <- apply(profile, 2, \(x) ifelse(x > min_value, x, fill_value)) %>% 
    data.frame(check.names = F)
  
  profile <- profile[rowSums(profile) != 0, colSums(profile) != 0]
  return(profile)
}

#### profile_adjacency ####
# Jinxin Meng, 20231119
# profile to adjacent matrix
profile_adjacency <- function(profile, logical = FALSE) {
  if (isTRUE(logical)) {
    profile <- apply(profile, 2, \(x) ifelse(x > 0, TRUE, FALSE)) %>% 
      data.frame(check.names = F)
  } else {
    profile <- apply(profile, 2, \(x) ifelse(x > 0, 1, 0)) %>% 
      data.frame(check.names = F)
  }
  return(profile)
}
  
#### profile_prevalence ####
# prevalence of each feature in all samples or samples belonged to each group
profile_prevalence <- function(profile, group, by_group = T, min_abundance = 0, count = F) {
  if (isTRUE(by_group) & missing(group)) stop("if by_group = TRUE, group (field: sample|group) should be provided ..")
  if (isFALSE(count)) {
    if (isTRUE(by_group)) {
      prevalence <- data.frame(t(profile), check.names = F) %>% 
        mutate(group = group$group[match(rownames(.), group$sample)]) %>% 
        group_by(group) %>% group_modify(~ purrr::map_df(.x, \(x) sum(x > min_abundance)/length(x))*100) %>%
        column_to_rownames(var = "group") %>% t %>% data.frame(check.names = F)
    } else if (isFALSE(by_group)) {
      prevalence <- data.frame(t(profile), check.names = F) %>% 
        map_vec(\(x) sum(x > min_abundance)/length(x)*100) %>% 
        data.frame(prevalence = .) %>% rownames_to_column(var = "name")
    }
  } else {
    if (isTRUE(by_group)) {
      prevalence <- data.frame(t(profile), check.names = F) %>% 
        mutate(group = group$group[match(rownames(.), group$sample)]) %>% 
        group_by(group) %>% group_modify(~ purrr::map_df(.x, \(x) sum(x > min_abundance))) %>%
        column_to_rownames(var = "group") %>% t %>% data.frame(check.names = F)
    } else if (isFALSE(by_group)) {
      prevalence <- data.frame(t(profile), check.names = F) %>% 
        map_vec(\(x) sum(x > min_abundance)) %>% 
        data.frame(prevalence = .) %>% rownames_to_column(var = "name")
    }
  }
  return(prevalence)
}

#### profile_statistics ####
profile_statistics <- function(profile, group, group_rename = NULL, 
                               by_group = TRUE) {
  profile <- data.frame(profile, check.names = F)
  
  if (isTRUE(by_group) & !missing(group)) { 
    if (!all(c("sample", "group") %in% colnames(group)) & is.null(group_rename)) stop("group field (sample|group)")
    if (!is.null(group_rename)) group <- data.frame(group, check.names = F) %>% dplyr::rename(all_of(group_rename))
    data <- profile %>% 
      t %>% 
      data.frame(check.names = F) %>% 
      mutate(group = group$group[match(rownames(.), group$sample)]) %>% 
    var = unique(group$group)[1]
    
    data <- cbind(
      data %>% 
        group_by(group) %>% 
        group_modify(~map_df(.x, \(x) mean(x))) %>% 
        column_to_rownames(var = "group") %>% 
        t %>% 
        data.frame(check.names = F) %>% 
        dplyr::rename_with(~paste0(.x, "_mean")),
      data %>% 
        group_by(group) %>% 
        group_modify(~map_df(.x, \(x) sd(x))) %>% 
        column_to_rownames(var = "group") %>% 
        t %>% 
        data.frame(check.names = F) %>% 
        dplyr::rename_with(~paste0(.x, "_sd")),
      data %>% 
        group_by(group) %>% 
        group_modify(~map_df(.x, \(x) median(x))) %>% 
        column_to_rownames(var = "group") %>% 
        t %>% 
        data.frame(check.names = F) %>% 
        dplyr::rename_with(~paste0(.x, "_median")),
      data %>% 
        group_by(group) %>% 
        group_modify(~map_df(.x, \(x) sum(x>0)/length(x)*100)) %>% 
        column_to_rownames(var = "group") %>% 
        t %>% 
        data.frame(check.names = F) %>% 
        dplyr::rename_with(~paste0(.x, "_prev"))
    ) %>% 
      relocate(starts_with(var))
  } else if (!isTRUE(by_group)) {
    
  } else {
    stop("if by_group = TRUE, group (field: sample|group) should be provided ..")
  }
  
}

#### profile_aggregate ####
# 按照行特征进行分类汇总
profile_aggregate <- function(profile, metadata, by = NULL, method = 'sum') {
  if (is.null(by))
    stop('what field used to aggregate ? ? ?')
  
  if (!all(c('name', by) %in% colnames(metadata)))
    stop('metadata field not found.')
  
  data <- data.frame(profile, check.names = F) %>%
    mutate(class = metadata[match(rownames(.), metadata$name), by],
           class = ifelse(is.na(class), 'Unknown', class)) %>% 
    aggregate(. ~ class, ., method) %>% 
    column_to_rownames('class')
  
  return(data)
}