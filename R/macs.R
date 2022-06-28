#' @export
macs_cols = function() {
  readr::cols(
    macs_chrom=readr::col_character(), macs_start=readr::col_double(), macs_end=readr::col_double(), macs_length=readr::col_character(), macs_summit_abs=readr::col_double(),
    macs_pileup=readr::col_double(), macs_pvalue=readr::col_double(), macs_fc=readr::col_double(), macs_qvalue=readr::col_double(), macs_name=readr::col_character(), macs_comment=readr::col_character()
  )
}

macs2_params = function(extsize=1e5, exttype="symetrical", llocal=1e7, minqvalue=NA_real_, minpvalue=NA_real_, maxgap=5e5, minlen=200, effective_size=1.87e9, baseline=0) {
  as.data.frame(list(exttype=exttype, extsize=extsize, llocal=llocal, minqvalue=minqvalue, minpvalue=minpvalue, maxgap=maxgap, minlen=minlen, effective_size=effective_size, baseline=baseline))
}


#' @export
macs_blank = function() {
  blank_tibble(macs_cols()) %>% dplyr::mutate(macs_sample=NA_character_, macs_group=NA_character_)
}

#' @export
macs2 = function(name, sample, effective_size, control=NULL, maxgap=NULL, qvalue=0.01, extsize=2000, slocal=50000, llocal=10000000, output_dir="data/macs2") {
  bed_sample = paste("-t", sample)
  bed_control = ifelse(is.null(control), "", paste("-c", control))
  maxgap = ifelse(is.null(maxgap), "", paste("--max-gap", sprintf("%0.0f", maxgap)))
  effective_size = sprintf("%0.0f", effective_size)

  cmd = paste("macs2 callpeak ", bed_sample, bed_control, "--seed 123 -f BED --keep-dup all --nomodel --bdg --mfold 5 100 --trackline", maxgap, "-g", effective_size, "-n", name, "--outdir", output_dir, "--slocal", sprintf("%0.0f", slocal), "--extsize", sprintf("%0.0f", extsize), "-q", qvalue, "--llocal", sprintf("%0.0f", llocal))
  # cmd = paste0("macs2 callpeak {bed_sample} {bed_control} --seed 123 {maxgap} -f BED -g {effsize} --keep-dup all -n {name} --outdir {output_dir} --nomodel --slocal {slocal} --extsize {extsize} -q {qvalue} --llocal {llocal} --bdg --trackline", bed_sample=bed_sample, bed_control=bed_control, name=name, output_dir=output_dir, extsize=extsize, qvalue=qvalue, maxgap=maxgap, llocal=sprintf("%0.0f", llocal), slocal=sprintf("%0.0f", slocal), effsize=effective_size)
  log(cmd)
  output = system(paste(cmd, " 2>&1"), intern = T)
  output = paste0(output, collapse="\n")
  log(output)

  output_df = readr::read_tsv(paste0(output_dir, "/", name, "_peaks.xls"), comment="#", col_names=names(macs_cols()$cols), col_types=macs_cols())
  output_df %>%
    dplyr::slice(-1) %>%
    dplyr::select(-dplyr::matches("macs_comment"))
}

macs2_coverage = function(sample_ranges, control_ranges=NULL, params, tmp_prefix=NULL, plot=F)
{
  if(is.null(tmp_prefix)) {
    tmp_prefix = file.path("tmp", basename(tempfile()))
  }

  cutoff = ifelse(!is.na(params$minqvalue), params$minqvalue, params$minpvalue)

  sample_path = paste0(tmp_prefix, "-sample.bdg")
  sample_df = as.data.frame(sample_ranges) %>% dplyr::mutate(sample_chrom=seqnames, sample_start=start, sample_end=end, sample_score=score)
  sample_ranges = sample_df %>% GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=T)
  readr::write_tsv(sample_df %>% dplyr::select(sample_chrom, sample_start, sample_end, sample_score), file=sample_path, col_names=F)

  bgmodel_data_df = suppressWarnings(sample_df %>%
    dplyr::mutate(sample_chrom=droplevels(sample_chrom)) %>%
    df2ranges(sample_chrom, sample_start, sample_end) %>%
    ranges2pure_tiles(column="sample_score", width=200, step=200, diff=0.1) %>%
    as.data.frame() %>%
    dplyr::filter(sample_score>0) %>%
    dplyr::rename(sample_chrom="seqnames", sample_start="start", sample_end="end") %>%
    dplyr::select(dplyr::matches("^sample_")))
  bgmodel_df = suppressWarnings(bgmodel_data_df %>%
    dplyr::group_by(sample_chrom, .drop=F) %>%
    dplyr::do((function(z) {
      yy<<-z

      fit_gamma = MASS::fitdistr(z$sample_score, densfun="gamma")
      cbind(z[1,"sample_chrom", drop=F], as.data.frame(t(fit_gamma$estimate)) %>% setNames(paste0("bgmodel_", colnames(.))))
    })(.)) %>%
    dplyr::ungroup())

  if(plot==T) {
    bgmodel_rand_df = bgmodel_df %>%
      dplyr::group_by(sample_chrom) %>%
      dplyr::summarize(sample_score=rgamma(10000, rate=bgmodel_rate, shape=bgmodel_shape))
    bgmodel_cutoff_df = bgmodel_df %>%
      dplyr::group_by(sample_chrom) %>%
      dplyr::summarize(cutoff=qgamma(cutoff, rate=bgmodel_rate, shape=bgmodel_shape, lower.tail=F))
    ggplot() +
      geom_density(aes(sample_score), bw=4, data=bgmodel_data_df %>% sample_n(100000)) +
      geom_density(aes(sample_score, color="fit"), data=bgmodel_rand_df) +
      geom_vline(aes(xintercept=cutoff), data=bgmodel_cutoff_df) +
      facet_wrap(~sample_chrom, scales="free")
  }

  writeLines("Detected bgmodel is \n==============================================\n")
  writeLines(knitr::kable(bgmodel_df %>% dplyr::mutate(dplyr::across(dplyr::matches("bgmodel_"), round, 4))))

  if(0 + !is.na(params$minqvalue) + !is.na(params$minpvalue) != 1) {
    stop("Please provide either minimal q-value or p-value")
  }

  control_path = paste0(tmp_prefix, "-control.bdg")
  qvalue_path = gsub(".bdg$", "_qvalue.bdg", sample_path)
  peaks_path = gsub(".bdg$", ".peaks", sample_path)
  bedpeaks_path = gsub(".bdg$", "_peaks.bed", sample_path)

  if(is.null(control_ranges) | length(control_ranges)==0) {
    control_df = sample_df %>%
      dplyr::left_join(bgmodel_df, by="sample_chrom") %>%
      # dplyr::mutate(score=tidyr::replace_na(sample_bgmodel, 0), control_score=score) %>%
      dplyr::select(control_chrom=sample_chrom, control_start=sample_start, control_end=sample_end, dplyr::matches("bgmodel_"))
    control_ranges = control_df %>% df2ranges(control_chrom, control_start, control_end)
  }
  # control_df = as.data.frame(control_ranges) %>% dplyr::select(seqnames, start, end, score)
  # readr::write_tsv(control_df, file=control_path, col_names=F)

  if(F) {
    # control_df = as.data.frame(control_ranges) %>% dplyr::select(seqnames, start, end, score)
    # readr::write_tsv(control_df, file=control_path, col_names=F)
    # This doesn't work anymore !!!! (because bgmodel is not always a poisson with just one argument
    # cmd_bdgcmp = stringr::str_glue("macs3 bdgcmp -t {sample} -c {control} -m qpois -o {output} -p 0.00000000000000001",
    #    sample=sample_path, control=control_path, output=qvalue_path)
    # writeLines(cmd_bdgcmp)
    # system(cmd_bdgcmp)
  } else {
    sample_df %>%
      dplyr::select(sample_chrom, sample_start, sample_end, sample_score) %>%
      readr::write_tsv("sample.bedgraph", col_names = F)
    qvalues_df = sample_df %>%
      dplyr::inner_join(control_df, by=c("sample_chrom"="control_chrom", "sample_start"="control_start", "sample_end"="control_end")) %>%
      dplyr::group_by(bgmodel_shape, bgmodel_rate) %>%
      dplyr::mutate(pvalue=pgamma(sample_score, shape=bgmodel_shape, rate=bgmodel_rate, lower.tail=F), bgmodel_signal=qgamma(cutoff,  bgmodel_shape, bgmodel_rate, lower.tail=F)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(qvalue_pvalue=pmin(pmax(pvalue, 0), 1)) %>%
      dplyr::mutate(qvalue_pvalue=ifelse(qvalue_pvalue<=0, 315, -log10(qvalue_pvalue))) %>%
      dplyr::mutate(qvalue=qvalue::qvalue(pvalue)$qvalues) %>%
      dplyr::mutate(qvalue_qvalue=ifelse(qvalue<=0, 315, -log10(qvalue))) %>%
      dplyr::ungroup() %>%
      dplyr::select(qvalue_chrom=seqnames, qvalue_start=start, qvalue_end=end, qvalue_qvalue, qvalue_pvalue, dplyr::starts_with("bgmodel_")) %>%
      dplyr::mutate(qvalue_score=dplyr::case_when(!is.na(params$minqvalue)~qvalue_qvalue, T~qvalue_pvalue))

    qvalues_df %>%
      dplyr::select(qvalue_chrom, qvalue_start, qvalue_end, qvalue_score) %>%
      readr::write_tsv(file=qvalue_path, col_names=F)
  }

  cmd_bdgpeakcall = stringr::str_glue("macs3 bdgpeakcall -i {qvalue} -c {cutoff} --min-length {format(minlen, scientific=F)} --max-gap {format(maxgap, scientific=F)} -o {output}",
     qvalue=qvalue_path, output=peaks_path, cutoff=-log10(cutoff), maxgap=params$maxgap, minlen=params$minlen)
  writeLines(cmd_bdgpeakcall)
  system(cmd_bdgpeakcall)

  # qvalue_cols = cols(qvalue_chrom=col_character(), qvalue_start=col_double(), qvalue_end=col_double(), qvalue_score=col_double())
  # qvalues_df = readr::read_tsv(qvalue_path, col_names=names(qvalue_cols$cols), col_types=qvalue_cols)

  peaks_cols = cols(
    island_chrom=col_character(), island_start=col_double(), island_end=col_double(), island_peak=col_character(), island_summit_abs=col_double(),
    island_score=col_character(), island_fc=col_double(), island_pvalue_log10=col_double(), island_qvalue_log10=col_double(), island_sammit_offset=col_double()
  )
  islands_df = readr::read_tsv(peaks_path, skip=1, col_names=names(peaks_cols$cols), col_types=peaks_cols) %>%
    dplyr::mutate(island_length=island_end-island_start, island_summit_pos=island_start + island_sammit_offset) %>%
    dplyr::select(island_chrom, island_start, island_end, island_summit_abs, island_sammit_offset, island_length, island_summit_pos)

  if(nrow(islands_df)>0) {
    islands_df$island_name = paste0("RDC_", stringr::str_pad(1:nrow(islands_df), 3, pad="0"))
    qvalues_ranges = qvalues_df %>% df2ranges(qvalue_chrom, qvalue_start, qvalue_end)
    islands_results_df = islands_df %>%
      df2ranges(island_chrom, island_start, island_end) %>%
      innerJoinByOverlaps(qvalues_ranges) %>%
      dplyr::group_by_at(colnames(islands_df)) %>%
      dplyr::mutate_at(dplyr::vars(dplyr::starts_with("bgmodel_")), mean) %>%
      dplyr::mutate(island_qvalue=max(qvalue_qvalue), island_pvalue=max(qvalue_pvalue)) %>%
      dplyr::ungroup() %>%
      dplyr::distinct_at(colnames(islands_df), .keep_all=T) %>%
      df2ranges(island_chrom, island_start, island_end) %>%
      innerJoinByOverlaps(sample_ranges) %>%
      dplyr::group_by_at(colnames(islands_df)) %>%
      dplyr::mutate(island_sammit_pos=round(sample_start[which.max(sample_score)]/2+sample_end[which.max(sample_score)]/2), island_summit_abs=max(sample_score)) %>%
      dplyr::ungroup() %>%
      dplyr::distinct_at(colnames(islands_df), .keep_all=T) %>%
      dplyr::select(island_name, dplyr::matches("^(island_|bgmodel_)"))

    #
    # Extended island coordinates
    #
    # sample_ranges1 = as.data.frame(sample_ranges) %>%
    #   dplyr::select(sample_chrom=seqnames, sample_start=start, sample_end=end, sample_score=score) %>%
    #   df2ranges(sample_chrom, sample_start, sample_end)
    # control_ranges1 = as.data.frame(control_ranges) %>%
    #   dplyr::select(control_chrom=seqnames, control_start=start, control_end=end, control_score=score) %>%
    #   df2ranges(control_chrom, control_start, control_end)
    # bgmodel_ranges = innerJoinByOverlaps(sample_ranges1, control_ranges1) %>%
    #   dplyr::filter(sample_score>=control_score) %>%
    #   df2ranges(sample_chrom, sample_start, sample_end) %>%
    #   GenomicRanges::reduce() %>%
    #   as.data.frame() %>%
    #   dplyr::select(island_extended_chrom=seqnames, island_extended_start=start, island_extended_end=end) %>%
    #   df2ranges(island_extended_chrom, island_extended_start, island_extended_end)
    # islands_df = islands_df %>%
    #   df2ranges(island_chrom, island_start, island_end) %>%
    #   leftJoinByOverlaps(bgmodel_ranges) %>%
    #   dplyr::group_by(island_chrom, island_start, island_end) %>%
    #   dplyr::mutate(island_extended_start=min(island_extended_start), island_extended_end=max(island_extended_end)) %>%
    #   dplyr::slice(1) %>%
    #   dplyr::select(-island_extended_chrom) %>%
    #   dplyr::ungroup()
  } else {
    islands_df$island_name = character()
    islands_df$island_summit_qvalue = double()
    islands_df$island_summit_abs = double()
    islands_df$island_summit_pos = double()
    islands_df$island_extended_start = double()
    islands_df$island_extended_end = double()
  }

  list(qvalues=qvalues_df, islands=islands_results_df)
}
