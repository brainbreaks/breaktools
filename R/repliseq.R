repliseq_read = function(path) {
  repliseq_df = as.data.frame(t(readr::read_tsv(path, col_names=F))) %>%
    dplyr::rename(repliseq_chrom="V1", repliseq_start="V2", repliseq_end="V3") %>%
    dplyr::mutate(repliseq_start=as.numeric(repliseq_start), repliseq_end=as.numeric(repliseq_end)) %>%
    reshape2::melt(id.vars=c("repliseq_chrom", "repliseq_start", "repliseq_end"), variable.name="repliseq_fraction", value.name="repliseq_value") %>%
    dplyr::mutate(repliseq_fraction=as.numeric(gsub("V", "", repliseq_fraction))-3, repliseq_value=as.numeric(repliseq_value))
  repliseq_df.keep = repliseq_df %>%
    dplyr::arrange(repliseq_start) %>%
    dplyr::group_by(repliseq_chrom, repliseq_start, repliseq_end) %>%
    dplyr::summarize(repliseq_isna=sum(is.na(repliseq_value))>10) %>%
    dplyr::group_by(repliseq_chrom) %>%
    dplyr::mutate(first_value_position=min(which(!repliseq_isna)), last_value_position=max(which(!repliseq_isna)), keep_value=dplyr::between(1:dplyr::n(), first_value_position[1], last_value_position[1])) %>%
    dplyr::filter(keep_value) %>%
    dplyr::ungroup() %>%
    dplyr::select(repliseq_chrom, repliseq_start, repliseq_end)
  repliseq_df.f = repliseq_df %>%
    dplyr::inner_join(repliseq_df.keep, by=c("repliseq_chrom", "repliseq_start", "repliseq_end"))

  repliseq_df.f
}

repliseq_summarize = function(repliseq_df, window=5) {
  th.repliseq_value_norm = 0.1
  repliseq_df = repliseq_df %>%
    dplyr::group_by(repliseq_chrom, repliseq_start, repliseq_end) %>%
    dplyr::mutate(repliseq_value_norm=((repliseq_value)/max(repliseq_value))) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(repliseq_chrom, repliseq_start, repliseq_end) %>%
    do((function(z) {
      zz<<-z

      i.max = z$repliseq_fraction[which.max(z$repliseq_value_norm)]
      lb = which(z$repliseq_value_norm < th.repliseq_value_norm & z$repliseq_fraction < i.max)
      lb = z$repliseq_fraction[ifelse(length(lb), lb[length(lb)]+1, 1)]
      ub = which(z$repliseq_value_norm < th.repliseq_value_norm & z$repliseq_fraction > i.max)
      ub = z$repliseq_fraction[ifelse(length(ub), ub[1]-1, nrow(z))]

      z$repliseq_value_in_scope = dplyr::between(z$repliseq_fraction, lb, ub)
      z
    })(.)) %>%
    dplyr::ungroup()

  repliseq_time_df = repliseq_df %>%
    dplyr::mutate(repliseqTime_chrom=repliseq_chrom, repliseqTime_start=repliseq_start, repliseqTime_end=repliseq_end) %>%
    dplyr::group_by(repliseqTime_chrom, repliseqTime_start, repliseqTime_end) %>%
    dplyr::mutate(repliseqTime_avg=weighted.mean(repliseq_fraction[repliseq_value_in_scope], repliseq_value_norm[repliseq_value_in_scope], na.rm=T)) %>%
    dplyr::mutate(lb=dplyr::between(repliseq_fraction, floor(repliseqTime_avg[1]-4), ceiling(repliseqTime_avg[1]-2)), repliseqTime_min=ifelse(any(lb), weighted.mean(repliseq_fraction[lb], repliseq_value_norm[lb]), NA_real_)) %>%
    dplyr::mutate(ub=dplyr::between(repliseq_fraction, floor(repliseqTime_avg[1]+2), ceiling(repliseqTime_avg[1]+4)), repliseqTime_max=ifelse(any(ub), weighted.mean(repliseq_fraction[ub], repliseq_value_norm[ub]), NA_real_)) %>%
    dplyr::summarize(repliseqTime_avg=repliseqTime_avg[1], repliseqTime_min=repliseqTime_min[1], repliseqTime_max=repliseqTime_max[1], repliseqTime_lb=min(which(repliseq_value_in_scope)), repliseqTime_ub=max(which(repliseq_value_in_scope))) %>%
    dplyr::group_by(repliseqTime_chrom) %>%
    dplyr::mutate(repliseqTime_avg=smoother::smth.gaussian(repliseqTime_avg, window=window), repliseqTime_min=smoother::smth.gaussian(repliseqTime_min, window=window), repliseqTime_max=smoother::smth.gaussian(repliseqTime_max, window=window)) %>%
    dplyr::mutate(repliseqTime_avg=zoo::na.fill(repliseqTime_avg, "extend"), repliseqTime_min=zoo::na.fill(repliseqTime_min, "extend"), repliseqTime_max=zoo::na.fill(repliseqTime_max, "extend")) %>%
    dplyr::ungroup()

  repliseq_time_df
}

repliseq_preprocess = function() {
  # Load repliseq data
  repliseq_ESC_df = repliseq_read("~/Workspace/Datasets/zhao_bmc_repliseq_2020/GSE137764_mESC_Gaussiansmooth_scaled_autosome.mat") %>% dplyr::mutate(repliseq_celltype="esc")
  repliseq_NPC_df = repliseq_read("~/Workspace/Datasets/zhao_bmc_repliseq_2020/GSE137764_mNPC_Gaussiansmooth_scaled_autosome.mat") %>% dplyr::mutate(repliseq_celltype="npc")
  repliseq_ESC2NPC_df = dplyr::bind_rows(repliseq_ESC_df, repliseq_NPC_df)
  readr::write_tsv(repliseq_ESC2NPC_df,"data/zhao_bmc_repliseq_2020/preprocessed/repliseq.tsv")

  repliseqTime_ESC_df = repliseq_summarize(repliseq_ESC_df) %>% dplyr::mutate(repliseqTime_celltype="esc")
  repliseqTime_NPC_df = repliseq_summarize(repliseq_NPC_df) %>% dplyr::mutate(repliseqTime_celltype="npc")
  repliseqTime_ESC2NPC_df = dplyr::bind_rows(repliseqTime_ESC_df, repliseqTime_NPC_df)
  readr::write_tsv(repliseqTime_ESC2NPC_df,"data/zhao_bmc_repliseq_2020/preprocessed/repliseqTime.tsv")


  repliseq_merged_df = repliseq_ESC2NPC_df %>%
    dplyr::group_by(repliseq_chrom, repliseq_start, repliseq_end, repliseq_fraction) %>%
    dplyr::summarize(repliseq_value=mean(repliseq_value, na.rm=T)) %>%
    dplyr::arrange(repliseq_chrom, repliseq_start, repliseq_end, repliseq_fraction)
  repliseqTime_merged_df = repliseq_summarize(repliseq_merged_df) %>% dplyr::mutate(repliseqTime_celltype="esc/npc")
  readr::write_tsv(repliseqTime_merged_df,"data/zhao_bmc_repliseq_2020/preprocessed/repliseqTime_merged.tsv")
}