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

macs2_coverage = function(sample_ranges, control_ranges=NULL, params, tmp_prefix=NULL)
{
  if(is.null(tmp_prefix)) {
    tmp_prefix = file.path("tmp", basename(tempfile()))
  }

  sample_path = paste0(tmp_prefix, "-sample.bdg")
  sample_df = as.data.frame(sample_ranges) %>% dplyr::mutate(sample_chrom=seqnames, sample_start=start, sample_end=end, sample_score=score)
  sample_ranges = sample_df %>% GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=T)
  readr::write_tsv(sample_df %>% dplyr::select(sample_chrom, sample_start, sample_end, sample_score), file=sample_path, col_names=F)

  baseline_df = sample_df %>%
    dplyr::mutate(sample_chrom=droplevels(sample_chrom)) %>%
    df2ranges(sample_chrom, sample_start, sample_end) %>%
    ranges_tile(column="sample_score", width=params$extsize, ceiling(params$extsize/4)) %>%
    as.data.frame() %>%
    dplyr::rename(sample_chrom="seqnames", sample_start="start", sample_end="end") %>%
    dplyr::select(dplyr::matches("^sample_")) %>%
    dplyr::filter(sample_end-sample_start>0 & sample_score>1e-6) %>%
    dplyr::group_by(sample_chrom, .drop=F) %>%
    dplyr::mutate(sample_density=sample_score/(sample_end-sample_start)) %>%
    dplyr::filter(quantile(sample_density, 0.05)<=sample_density & sample_density<=median(sample_density, na.rm=T)*2) %>%
    dplyr::summarise(sample_baseline=sum(sample_score*(sample_end-sample_start))/sum(sample_end-sample_start)) %>%
    dplyr::mutate(sample_baseline=tidyr::replace_na(sample_baseline, 0), sample_baseline=pmax(sample_baseline, params$baseline))

  # baseline_df = sample_df %>%
  #   dplyr::mutate(sample_chrom=droplevels(sample_chrom)) %>%
  #   dplyr::filter(sample_end-sample_start>0 & sample_score>1e-6) %>%
  #   dplyr::group_by(sample_chrom, .drop=F) %>%
  #   dplyr::mutate(sample_density=sample_score/(sample_end-sample_start)) %>%
  #   dplyr::filter(quantile(sample_density, 0.05)<=sample_density & sample_density<=median(sample_density, na.rm=T)*2) %>%
  #   dplyr::summarise(sample_baseline=sum(sample_score*(sample_end-sample_start))/sum(sample_end-sample_start)) %>%
  #   dplyr::mutate(sample_baseline=tidyr::replace_na(sample_baseline, 0), sample_baseline=pmax(sample_baseline, params$baseline))
  writeLines(paste0("Detected baseline is \n", paste(paste0("    ", baseline_df$sample_chrom, "=", format(round(baseline_df$sample_baseline, 5), nsmall=5)), collapse="\n")))

  if(0 + !is.na(params$minqvalue) + !is.na(params$minpvalue) != 1) {
    stop("Please provide either minimal q-value or p-value")
  }

  control_path = paste0(tmp_prefix, "-control.bdg")
  if(is.null(control_ranges) | length(control_ranges)==0) {
    control_ranges = sample_df %>%
      dplyr::left_join(baseline_df, by="sample_chrom") %>%
      dplyr::mutate(score=tidyr::replace_na(sample_baseline, 0), control_score=score) %>%
      dplyr::select(control_chrom=sample_chrom, control_start=sample_start, control_end=sample_end, score, control_score) %>%
      df2ranges(control_chrom, control_start, control_end)
  }
  control_df = as.data.frame(control_ranges) %>%
    dplyr::mutate(score=score) %>%
    dplyr::select(seqnames, start, end, score)
  readr::write_tsv(control_df, file=control_path, col_names=F)

  qvalue_path = gsub(".bdg$", "_qvalue.bdg", sample_path)
  peaks_path = gsub(".bdg$", ".peaks", sample_path)
  bedpeaks_path = gsub(".bdg$", "_peaks.bed", sample_path)

  if(F) {
    cmd_bdgcmp = stringr::str_glue("macs3 bdgcmp -t {sample} -c {control} -m qpois -o {output} -p 0.00000000000000001",
       sample=sample_path, control=control_path, output=qvalue_path)
    writeLines(cmd_bdgcmp)
    system(cmd_bdgcmp)
  } else {
    qvalues_df = control_df %>%
      dplyr::select(seqnames, start, end, control_score=score) %>%
      dplyr::inner_join(sample_df, by=c("seqnames", "start", "end")) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(pvalue=pgamma(sample_score, shape=control_score, rate=1, lower.tail=F)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(qvalue_pvalue=pmin(pmax(pvalue, 0), 1)) %>%
      dplyr::mutate(qvalue_pvalue=ifelse(qvalue_pvalue<=0, 315, -log10(qvalue_pvalue))) %>%
      dplyr::group_by(seqnames) %>%
      dplyr::mutate(qvalue_qvalue=qvalue::qvalue(pvalue)$qvalues) %>%
      dplyr::mutate(qvalue_qvalue=ifelse(qvalue_qvalue<=0, 315, -log10(qvalue_qvalue))) %>%
      dplyr::ungroup() %>%
      # dplyr::mutate(qvalue_score=qvalue::qvalue(pvalue)$qvalues) %>%
      # dplyr::mutate(qvalue_score=ifelse(qvalue_score==0, 315, -log10(qvalue_score))) %>% # TODO: change 315 to something better
      # dplyr::mutate(qvalue_score=tidyr::replace_na(qvalue_score, 0)) %>%
      # dplyr::mutate(pvalue_score=ifelse(pvalue==0, 315, -log10(pvalue))) %>%
      # dplyr::mutate(qvalue_score=0) %>%
      dplyr::select(qvalue_chrom=seqnames, qvalue_start=start, qvalue_end=end, qvalue_qvalue, qvalue_pvalue)
    qvalues_df %>%
      dplyr::mutate(score=dplyr::case_when(!is.na(params$minqvalue)~qvalue_qvalue, T~qvalue_pvalue)) %>%
      dplyr::select(qvalue_chrom, qvalue_start, qvalue_end, score) %>%
      readr::write_tsv(file=qvalue_path, col_names=F)
  }
  # qvalues_df = readr::read_tsv(qvalue_path, col_names=names(qvalue_cols$cols), col_types=qvalue_cols)
  # table(qvalues_df$qvalue_score)
  # x = qvalues_df %>%
  #   df2ranges(qvalue_chrom, qvalue_start, qvalue_end) %>%
  #   innerJoinByOverlaps(sample_ranges)
  # plot(x$qvalue_score, x$sample_score)

  cutoff = ifelse(!is.na(params$minqvalue), params$minqvalue, params$minpvalue)
  cmd_bdgpeakcall = stringr::str_glue("macs3 bdgpeakcall -i {qvalue} -c {cutoff} --min-length {format(minlen, scientific=F)} --max-gap {format(maxgap, scientific=F)} -o {output}",
     qvalue=qvalue_path, output=peaks_path, cutoff=-log10(cutoff), maxgap=params$maxgap, minlen=params$minlen)
  writeLines(cmd_bdgpeakcall)
  system(cmd_bdgpeakcall)

  qvalue_cols = cols(qvalue_chrom=col_character(), qvalue_start=col_double(), qvalue_end=col_double(), qvalue_score=col_double())
  qvalues_df = readr::read_tsv(qvalue_path, col_names=names(qvalue_cols$cols), col_types=qvalue_cols)

  peaks_cols = cols(
    island_chrom=col_character(), island_start=col_double(), island_end=col_double(), island_peak=col_character(), island_summit_abs=col_double(),
    island_score=col_character(), island_fc=col_double(), island_pvalue_log10=col_double(), island_qvalue_log10=col_double(), island_sammit_offset=col_double()
  )
  islands_df = readr::read_tsv(peaks_path, skip=1, col_names=names(peaks_cols$cols), col_types=peaks_cols) %>%
    dplyr::mutate(island_length=island_end-island_start, island_summit_pos=island_start + island_sammit_offset)

  if(nrow(islands_df)>0) {
    islands_df$island_name = paste0("MACS3_", stringr::str_pad(1:nrow(islands_df), 3, pad="0"))
    qvalues_ranges = qvalues_df %>% df2ranges(qvalue_chrom, qvalue_start, qvalue_end)
    summit_ranges = islands_df %>% df2ranges(island_chrom, island_summit_pos, island_summit_pos)
    summit2qvalue_df = as.data.frame(IRanges::mergeByOverlaps(qvalues_ranges, summit_ranges)) %>%
      dplyr::group_by(island_name) %>%
      dplyr::summarise(island_summit_qvalue=max(qvalue_score), .groups="keep")
    summit2sample_df = as.data.frame(IRanges::mergeByOverlaps(sample_ranges, summit_ranges)) %>%
      dplyr::group_by(island_name) %>%
      dplyr::summarise(island_summit_abs=max(score), .groups="keep")

    islands_df = islands_df %>%
      dplyr::select(-dplyr::matches("island_summit_abs|island_summit_qvalue")) %>%
      dplyr::inner_join(summit2qvalue_df, by="island_name") %>%
      dplyr::inner_join(summit2sample_df, by="island_name")

    sample_ranges1 = as.data.frame(sample_ranges) %>%
      dplyr::select(sample_chrom=seqnames, sample_start=start, sample_end=end, sample_score=score) %>%
      df2ranges(sample_chrom, sample_start, sample_end)
    control_ranges1 = as.data.frame(control_ranges) %>%
      dplyr::select(control_chrom=seqnames, control_start=start, control_end=end, control_score=score) %>%
      df2ranges(control_chrom, control_start, control_end)
    baseline_ranges = innerJoinByOverlaps(sample_ranges1, control_ranges1) %>%
      dplyr::filter(sample_score>=control_score) %>%
      df2ranges(sample_chrom, sample_start, sample_end) %>%
      GenomicRanges::reduce() %>%
      as.data.frame() %>%
      dplyr::select(island_extended_chrom=seqnames, island_extended_start=start, island_extended_end=end) %>%
      df2ranges(island_extended_chrom, island_extended_start, island_extended_end)
    islands_df = islands_df %>%
      df2ranges(island_chrom, island_start, island_end) %>%
      leftJoinByOverlaps(baseline_ranges) %>%
      dplyr::group_by(island_chrom, island_start, island_end) %>%
      dplyr::mutate(island_extended_start=min(island_extended_start), island_extended_end=max(island_extended_end)) %>%
      dplyr::slice(1) %>%
      dplyr::select(-island_extended_chrom) %>%
      dplyr::ungroup()
  } else {
    islands_df$island_name = character()
    islands_df$island_summit_qvalue = double()
    islands_df$island_summit_abs = double()
    islands_df$island_summit_pos = double()
    islands_df$island_extended_start = double()
    islands_df$island_extended_end = double()
  }

  # Find coverage areas not overlapping with islands
  notisland_sample_ranges = sample_ranges[!(IRanges::overlapsAny(sample_ranges, islands_df %>% df2ranges(island_chrom, island_start, island_end)))]

  # Calculate baseline and signal-to-noise ratio for each island separately
  island_llocal_ranges = islands_df %>%
    df2ranges(island_chrom, island_start-params$llocal/2, island_end+params$llocal/2)
  island2baseline_df = as.data.frame(IRanges::mergeByOverlaps(island_llocal_ranges, notisland_sample_ranges)) %>%
    dplyr::select(-dplyr::matches("_ranges\\.")) %>%
    dplyr::group_by(island_name) %>%
    dplyr::mutate(island_density=sample_score/(sample_end-sample_start)) %>%
    dplyr::filter(0.01>=island_density & island_density<=median(island_density, na.rm=T)*2) %>%
    dplyr::summarise(island_baseline=sum(sample_score*(sample_end-sample_start))/sum(sample_end-sample_start))
  islands_results_df = islands_df %>%
    dplyr::inner_join(island2baseline_df, by="island_name") %>%
    dplyr::mutate(island_snr=island_summit_abs/island_baseline) %>%
    dplyr::select(island_name, island_chrom, island_start, island_end, island_extended_start, island_extended_end, island_length, island_summit_pos, island_summit_qvalue, island_summit_abs, island_baseline, island_snr)

  islands_results_df %>%
    dplyr::mutate(strand=".") %>%
    dplyr::select(island_chrom, island_start, island_end, island_name, island_summit_qvalue, strand) %>%
    readr::write_tsv(bedpeaks_path, col_names=F)

  list(qvalues=qvalues_df, islands=islands_results_df)
}
