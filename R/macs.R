#' @export
macs_cols = function() {
  readr::cols(
    macs_chrom=readr::col_character(), macs_start=readr::col_double(), macs_end=readr::col_double(), macs_length=readr::col_character(), macs_summit_abs=readr::col_double(),
    macs_pileup=readr::col_double(), macs_pvalue=readr::col_double(), macs_fc=readr::col_double(), macs_qvalue=readr::col_double(), macs_name=readr::col_character(), macs_comment=readr::col_character()
  )
}

macs2_params = function(extsize=1e5, exttype="symetrical", llocal=1e7, minqvalue=NA_real_, minpvalue=NA_real_, maxgap=5e5, minlen=100e3, seedlen=50e3, seedgap=50e3, effective_size=1.87e9, baseline=0) {
  as.data.frame(list(exttype=exttype, extsize=extsize, llocal=llocal, minqvalue=minqvalue, minpvalue=minpvalue, maxgap=maxgap, seedgap=seedgap, minlen=minlen, seedlen=seedlen, effective_size=effective_size, baseline=baseline)) %>%
    dplyr::mutate(minsignif=ifelse(!is.na(minqvalue), minqvalue, minpvalue))
}


#' @export
macs_blank = function() {
  blank_tibble(macs_cols()) %>% dplyr::mutate(macs_sample=NA_character_, macs_group=NA_character_)
}

macs2_coverage_bgmodel = function(coverage_ranges, distr, coverage_column="score", mask_ranges=NULL, debug_plots=F, plot_cutoff=0.01, plot_lower_tail=F)
{
  if(!(coverage_column %in% colnames(GenomicRanges::mcols(coverage_ranges)))) {
    stop(paste0("Column '", coverage_column, "' not found in coverage ranges data columns"))
  }

  coverage_df = as.data.frame(coverage_ranges) %>%
    dplyr::select_at(c(coverage_chrom="seqnames", coverage_start="start", coverage_end="end", bgmodel_strand="strand", coverage_score=coverage_column))
  if(is.null(mask_ranges)) {
    mask_ranges = GenomicRanges::GRanges()
  }

  # as.data.frame(mask_ranges) %>% dplyr::select(mask_chrom, mask_start, mask_end) %>% readr::write_tsv("reports/detect_rdc/mask.bed", col_names=F)

  bgmodel_data_df = suppressWarnings(coverage_ranges %>%
    ranges_sample(column=coverage_column, mask_ranges=mask_ranges, ntile=1e5) %>%
    as.data.frame() %>%
    dplyr::select_at(c(bgmodel_chrom="seqnames", coverage_start="start", coverage_end="end", bgmodel_strand="strand", coverage_score=coverage_column)) %>%
    df2ranges(bgmodel_chrom, coverage_start, coverage_end, bgmodel_strand) %>%
    leftJoinByOverlaps(mask_ranges) %>%
    dplyr::filter(is.na(mask_start)) %>%
    dplyr::select(dplyr::matches("^(bgmodel_|coverage_)")))
  bgmodel_quantiles_df = bgmodel_data_df %>%
    dplyr::group_by(bgmodel_chrom, bgmodel_strand) %>%
    dplyr::summarize(coverage_score.q999=quantile(coverage_score[coverage_score>0], 0.999), coverage_score.q500=quantile(coverage_score[coverage_score>0], 0.5), dataset="Observed data") %>%
    data.frame()

  bgmodel_df = suppressWarnings(bgmodel_data_df %>%
    dplyr::filter(bgmodel_chrom!="chrY") %>%
    dplyr::group_by(bgmodel_chrom, bgmodel_strand) %>%
    dplyr::do((function(z) {
      yy<<-z
      fixargs = NULL
      probs = c(pmin(0.1, mean(z$coverage_score==0)), pmax(1-mean(z$coverage_score>=3), 0.9))
      fit_gamma = fitdistrplus::fitdist(z$coverage_score, distr=distr, fix.arg=fixargs, method="qme", probs=probs)
      if(!is.null(fixargs) & length(fixargs)>0) {
        fit_gamma$estimate = unlist(c(fit_gamma$estimate, fixargs[setdiff(names(fixargs), names(fit_gamma$estimate))]))
      }

      cbind(bgmodel_distr=distr, as.data.frame(t(fit_gamma$estimate)) %>% setNames(paste0("bgmodel_", colnames(.))))
    })(.)) %>%
    dplyr::ungroup())

  writeLines("Detected bgmodel is \n==============================================\n")
  writeLines(knitr::kable(bgmodel_df))

  if(debug_plots==T) {
    bgmodel_quant_df = bgmodel_df %>%
      dplyr::group_by(bgmodel_chrom, bgmodel_strand) %>%
      dplyr::summarize(coverage_quantile=seq(0, 1, length.out=500), coverage_density=distr_call(ver="q", distr=distr, x=coverage_quantile, params=dplyr::cur_data_all()))
    bgmodel_data_quant_df = bgmodel_data_df %>%
      dplyr::group_by(bgmodel_chrom, bgmodel_strand) %>%
      dplyr::summarize(coverage_quantile=seq(0, 1, length.out=500), coverage_density=quantile(coverage_score, coverage_quantile))
    bgmodel_qplot_df = bgmodel_quant_df %>%
      dplyr::inner_join(bgmodel_data_quant_df, by=c("bgmodel_chrom", "bgmodel_strand", "coverage_quantile"))
    bgmodel_data_prepared_df = bgmodel_data_df %>%
      dplyr::mutate(dataset="Observed") %>%
      dplyr::inner_join(bgmodel_quantiles_df %>% dplyr::distinct(bgmodel_chrom, bgmodel_strand, coverage_score.q500), by=c("bgmodel_chrom", "bgmodel_strand"))

    print(ggplot2::ggplot(bgmodel_qplot_df) +
      ggplot2::geom_hex(ggplot2::aes(x=coverage_density.x, y=coverage_density.y), alpha=0.5, shape=1, size=2) +
      ggplot2::facet_wrap(~bgmodel_chrom+bgmodel_strand, scales="free") +
      ggplot2::geom_abline(intercept=0, slope=1) +
      ggplot2::labs(x="Fitted model", y="Observed dataodel", title=paste0("Q-Q plot -- ", coverage_df$tlx_group[1], ", ", distr, " distribution")))

    bgmodel_rand_df = bgmodel_df %>%
      dplyr::inner_join(bgmodel_quantiles_df %>% dplyr::distinct(bgmodel_chrom, bgmodel_strand, coverage_score.q500), by=c("bgmodel_chrom", "bgmodel_strand")) %>%
      dplyr::inner_join(bgmodel_data_df %>% dplyr::distinct(bgmodel_chrom, bgmodel_strand), by=c("bgmodel_chrom", "bgmodel_strand")) %>%
      dplyr::group_by(bgmodel_chrom, bgmodel_strand) %>%
      dplyr::summarize(dataset="Fitted model", coverage_score=round(seq(0, coverage_score.q500*10, length.out=10000)), coverage_density=distr_call(ver="d", distr=distr, x=coverage_score, params=dplyr::cur_data_all()))
    bgmodel_cutoff_df = bgmodel_df %>%
      dplyr::group_by(bgmodel_chrom) %>%
      dplyr::summarize(cutoff_score=distr_call("q", distr, plot_cutoff, params=dplyr::cur_data_all(), lower.tail=plot_lower_tail), dataset=paste0("Sign. cutoff (p<=", round(plot_cutoff, 3), ")"))
    # data=bgmodel_data_df %>% dplyr::filter(coverage_score<coverage_score.q500*10)
    print(ggplot2::ggplot(mapping=ggplot2::aes(x=coverage_score, color=bgmodel_strand, linetype=dataset)) +
      ggplot2::geom_density(bw=1, data=bgmodel_data_prepared_df %>% dplyr::filter(coverage_score<coverage_score.q500*10) %>% dplyr::sample_n(100000, replace=T)) +
      ggplot2::geom_line(ggplot2::aes(y=coverage_density), data=bgmodel_rand_df) +
      ggplot2::geom_vline(ggplot2::aes(xintercept=cutoff_score, color=dataset), data=bgmodel_cutoff_df) +
      ggplot2::facet_wrap(~bgmodel_chrom, scales="free") +
      ggplot2::labs(x="Pileup", y="Density", title=paste0("Pileup density plot ", distr, " distribution")))
  }

  bgmodel_df
}

macs2_coverage = function(sample_ranges, control_ranges=NULL, params, bgmodel_df, extended_islands=F, extended_islands_dist=1e6, extended_islands_significance=0.1, debug_plots=F)
{
  sample_df = as.data.frame(sample_ranges) %>% dplyr::mutate(sample_chrom=seqnames, sample_start=start, sample_end=end, sample_strand=strand, sample_score=score)
  sample_ranges = sample_df %>% GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=T)

  if(0 + !is.na(params$minqvalue) + !is.na(params$minpvalue) != 1) {
    stop("Please provide either minimal q-value or p-value")
  }

  if(is.null(control_ranges) | length(control_ranges)==0) {
    control_df = sample_df %>%
      dplyr::inner_join(bgmodel_df, by=c("sample_chrom"="bgmodel_chrom", "sample_strand"="bgmodel_strand")) %>%
      dplyr::select(control_chrom=sample_chrom, control_start=sample_start, control_end=sample_end, control_strand=sample_strand, dplyr::matches("bgmodel_"))
  }

  if(F) {
    # control_df = as.data.frame(control_ranges) %>% dplyr::select(seqnames, start, end, score)
    # readr::write_tsv(control_df, file=control_path, col_names=F)
    # This doesn't work anymore !!!! (because bgmodel is not always a poisson with just one argument
    # cmd_bdgcmp = stringr::str_glue("macs3 bdgcmp -t {sample} -c {control} -m qpois -o {output} -p 0.00000000000000001",
    #    sample=sample_path, control=control_path, output=qvalue_path)
    # writeLines(cmd_bdgcmp)
    # system(cmd_bdgcmp)
  } else {
    pvalues_df = sample_df %>%
      dplyr::inner_join(control_df, by=c("sample_chrom"="control_chrom", "sample_start"="control_start", "sample_end"="control_end", "sample_strand"="control_strand")) %>%
      dplyr::group_by_at(colnames(control_df)[grepl("bgmodel_", colnames(control_df))]) %>%
      dplyr::mutate(pvalue=distr_call("p", distr=bgmodel_distr[1], x=sample_score, params=dplyr::cur_data_all(), lower.tail=F), bgmodel_signal=distr_call("q", distr=bgmodel_distr[1], x=params$minsignif, params=dplyr::cur_data_all(), lower.tail=F)) %>%
      dplyr::mutate(pvalue=pmax(pvalue, dpois(sample_score, 2))) %>%
      dplyr::ungroup()

    if(debug_plots==T) {
      print(ggplot2::ggplot(pvalues_df) +
        ggplot2::geom_histogram(ggplot2::aes(x=pvalue)) +
        ggplot2::facet_wrap(~sample_chrom+sample_strand, scales="free") +
        ggplot2::labs(title=sample_df$tlx_group[1]))
    }

    qvalues_df = pvalues_df %>%
      dplyr::mutate(qvalue_pvalue=pmin(pmax(pvalue, 0), 1)) %>%
      dplyr::mutate(qvalue_pvalue=ifelse(qvalue_pvalue<=0, 315, -log10(qvalue_pvalue))) %>%
      dplyr::mutate(qvalue=p.adjust(pvalue, "BH")) %>%
      # dplyr::mutate(qvalue=qvalue::qvalue(pvalue)$qvalues) %>%
      dplyr::mutate(qvalue_qvalue=ifelse(qvalue<=0, 315, -log10(qvalue))) %>%
      dplyr::ungroup() %>%
      dplyr::select(qvalue_chrom=seqnames, qvalue_start=start, qvalue_end=end, qvalue_strand=strand, qvalue_qvalue, qvalue_pvalue, dplyr::starts_with("bgmodel_")) %>%
      dplyr::mutate(qvalue_score=tidyr::replace_na(dplyr::case_when(!is.na(params$minqvalue)~qvalue_qvalue, T~qvalue_pvalue), 0))

    # Old MACS3 code
    # x = qvalues_df %>%
    #   dplyr::mutate(qvalue_chrom=paste0(qvalue_chrom, ";", qvalue_strand)) %>%
    #   dplyr::select(qvalue_chrom, qvalue_start, qvalue_end, qvalue_score) %>%
    #   readr::write_tsv(file=qvalue_path, col_names=F)
  }

  #
  # Find islands
  #
  qvalues_sign_ranges = qvalues_df %>%
    dplyr::filter(qvalue_score>=-log10(params$minsignif)) %>%
    df2ranges(qvalue_chrom, qvalue_start, qvalue_end, qvalue_strand)

  seed_ranges = qvalues_sign_ranges %>% GenomicRanges::reduce(min.gapwidth=params$seedgap, ignore.strand=F)
  islands_ranges = seed_ranges[GenomicRanges::width(seed_ranges)>=params$seedlen] %>%
    GenomicRanges::reduce(min.gapwidth=params$maxgap, ignore.strand=F)  %>%
    as.data.frame() %>%
    dplyr::select(island_chrom=seqnames, island_strand=strand, island_start=start, island_end=end, island_length=width) %>%
    dplyr::filter(extended_islands | island_length>=params$minlen) %>%
    df2ranges(island_chrom, island_start, island_end, island_strand) %>%
    innerJoinByOverlaps(qvalues_sign_ranges) %>%
    dplyr::mutate(island_score=qvalue_score, island_pvalue=qvalue_pvalue, island_qvalue=qvalue_qvalue) %>%
    dplyr::arrange(dplyr::desc(island_score)) %>%
    dplyr::distinct(island_chrom, island_start, island_end, island_strand, .keep_all=T) %>%
    dplyr::select(dplyr::starts_with("island_")) %>%
    df2ranges(island_chrom, island_start, island_end, island_strand)


  #
  # Find peaks
  #
  peaks_ranges = suppressWarnings(innerJoinByOverlaps(islands_ranges, sample_ranges) %>%
    dplyr::group_by(island_chrom, island_start, island_end, island_strand, island_length) %>%
    dplyr::mutate(island_summit_abs=max(sample_score, na.rm=T))) %>%
    dplyr::filter(sample_score==island_summit_abs) %>%
    df2ranges(sample_chrom, sample_start, sample_end, sample_strand) %>%
    GenomicRanges::reduce(min.gapwidth=10) %>%
    as.data.frame() %>%
    dplyr::select(peak_chrom=seqnames, peak_strand=strand, peak_start=start, peak_end=end, peak_length=width) %>%
    df2ranges(peak_chrom, peak_start, peak_end, peak_strand) %>%
    innerJoinByOverlaps(sample_ranges) %>%
    dplyr::group_by(peak_chrom, peak_start, peak_end, peak_strand, peak_length) %>%
    dplyr::summarize(island_summit_pos=round(peak_start/2+peak_end/2)[1], island_summit_abs=max(sample_score, na.rm=T)) %>%
    dplyr::ungroup() %>%
    df2ranges(peak_chrom, peak_start, peak_end, peak_strand)

  #
  # Final islands data.frame with peak position and pileup
  #
  islands_df = innerJoinByOverlaps(islands_ranges, peaks_ranges) %>%
    dplyr::arrange(dplyr::desc(peak_length)) %>%
    dplyr::distinct(island_chrom, island_start, island_end, island_strand, island_length, .keep_all=T) %>%
    dplyr::select(dplyr::starts_with("island_"))
  writeLines(paste0("Found ", nrow(islands_df), " initial islands (", paste0(names(table(islands_df$island_strand)), ":", table(islands_df$island_strand), collapse=","), ")"))

  #
  # Calculate extended island coordinates and reduce overlapping RDC to single instance
  #
  if(nrow(islands_df)>0) {
    if(extended_islands) {
      #
      # Extend each island with smoothing function
      #
      writeLines("Extending islands...")
      extended_significance_log10 = -log10(extended_islands_significance)
      significance_log10 = -log10(params$minsignif)
      qvalues_ranges = qvalues_df %>% df2ranges(qvalue_chrom, qvalue_start, qvalue_end, qvalue_strand)
      islands_extended_ranges = islands_df %>% df2ranges(island_chrom, island_start-extended_islands_dist, island_end+extended_islands_dist, island_strand)

      island_qvalues_df = islands_extended_ranges %>%
        innerJoinByOverlaps(qvalues_ranges) %>%
        dplyr::filter(island_strand==qvalue_strand) %>%
        dplyr::group_by_at(colnames(islands_df)) %>%
        dplyr::mutate(
          island_qvalue=max(island_qvalue),
          island_pvalue=max(island_pvalue),
          island_score=max(island_score)) %>%
        dplyr::ungroup()

      islands_extended_df = island_qvalues_df %>%
        # dplyr::filter(island_chrom=="chr10" & island_start>=80110000 & island_end<=80200000 & island_strand=="+") %>%
        # dplyr::filter(island_chrom=="chr6" & island_start==33348072 & island_end==33398072) %>%
        # dplyr::filter(island_chrom=="chr8" & island_start==16283635 & island_end==17255978) %>%
        # dplyr::filter(island_chrom=="chr8" & island_start==48737078 & island_end==48769021) %>%
        # dplyr::filter(island_chrom=="chr8" & island_start==120249366 & island_end==120296331) %>%
        dplyr::group_by_at(colnames(islands_df)) %>%
        dplyr::do((function(y){
          yy<<-y
          y_island = y %>% dplyr::select(island_chrom, island_start, island_end, island_strand) %>% dplyr::slice(1)

          # Smoothen pileup
          y_loess_span = 200e3/(max(y$qvalue_end)-min(y$qvalue_start))
          y_loess_step = 2000
          y_ranges = y %>%
            dplyr::rename(raw_chrom="qvalue_chrom", raw_start="qvalue_start", raw_end="qvalue_end") %>%
            df2ranges(raw_chrom, raw_start, raw_end)
          y_long = data.frame(qvalue_chrom=y_island$island_chrom, qvalue_start=seq(min(y$qvalue_start), max(y$qvalue_end), by=y_loess_step)) %>%
            dplyr::mutate(qvalue_end=qvalue_start+1) %>%
            df2ranges(qvalue_chrom, qvalue_start, qvalue_end) %>%
            innerJoinByOverlaps(y_ranges) %>%
            dplyr::mutate(qvalue_end=qvalue_start+y_loess_step-1) %>%
            dplyr::select(dplyr::matches("island_|qvalue_")) %>%
            dplyr::mutate(qvalue_pos=qvalue_start)

          y_model = loess(qvalue_pvalue~qvalue_pos, data=y_long, span=y_loess_span)
          y_smooth = y_long %>%
            dplyr::mutate(qvalue_pvalue = predict(y_model, .)) %>%
            dplyr::mutate(qvalue_chrom=island_chrom, qvalue_strand=island_strand) %>%
            dplyr::arrange(qvalue_pos) %>%
            dplyr::mutate(qvalue_start=qvalue_pos, qvalue_end=dplyr::lead(qvalue_pos, 1)) %>%
            dplyr::arrange(qvalue_start) %>%
            dplyr::mutate(qvalue_pvalue.diff=(qvalue_pvalue-dplyr::lag(qvalue_pvalue, 1))/(qvalue_end-qvalue_start)) %>%
            dplyr::filter(!is.na(qvalue_pvalue) & !is.na(qvalue_end) & !is.na(qvalue_pvalue.diff))

          # Find extended island boundaries
          y_extended = y_smooth %>%
            dplyr::filter(
              qvalue_pvalue >= extended_significance_log10 & (qvalue_start <= island_start & qvalue_pvalue.diff>=1e-6 | qvalue_end >= island_end & qvalue_pvalue.diff<=-1e-6)
            ) %>%
            df2ranges(island_chrom, qvalue_start, qvalue_end, island_strand) %>%
            GenomicRanges::reduce(min.gapwidth=params$seedgap, ignore.strand=F) %>%
            as.data.frame() %>%
            dplyr::rename(extended_chrom="seqnames", extended_start="start", extended_end="end", extended_strand="strand") %>%
            dplyr::bind_cols(y_island) %>%
            dplyr::filter(extended_strand==island_strand & (extended_start<=island_start & island_start-extended_end<=params$maxgap | extended_end>=island_end & extended_start-island_end<=params$maxgap))

          res_df = data.frame(island_extended_start=min(c(y$island_start, y_extended$extended_start)), island_extended_end=max(c(y$island_end, y_extended$extended_end)))

          cbind(y_island, res_df)
        })(.)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(island_extended_length=island_extended_end-island_extended_start) %>%
        dplyr::distinct_at(colnames(islands_df), .keep_all=T) %>%
        df2ranges(island_chrom, island_start, island_end, island_strand) %>%
        innerJoinByOverlaps(sample_ranges) %>%
        dplyr::group_by_at(colnames(islands_df)) %>%
        dplyr::mutate(island_sammit_pos=round(sample_start[which.max(sample_score)]/2+sample_end[which.max(sample_score)]/2), island_summit_abs=max(sample_score)) %>%
        dplyr::ungroup() %>%
        dplyr::distinct_at(colnames(islands_df), .keep_all=T) %>%
        dplyr::select(dplyr::matches("^(island_|bgmodel_)"))

    islands_extended_ranges = islands_extended_df %>%
      df2ranges(island_chrom, island_extended_start, island_extended_end, island_strand)
    islands_results_df = islands_extended_ranges %>%
      GenomicRanges::reduce(ignore.strand=F) %>%
      as.data.frame() %>%
      dplyr::rename(reduced_chrom="seqnames", reduced_start="start", reduced_end="end", reduced_strand="strand") %>%
      df2ranges(reduced_chrom, reduced_start, reduced_end, reduced_strand) %>%
      innerJoinByOverlaps(islands_extended_ranges) %>%
      dplyr::group_by(island_chrom=reduced_chrom, island_extended_start=reduced_start, island_extended_end=reduced_end, island_strand=reduced_strand) %>%
      dplyr::summarize(
        island_start=min(island_start),
        island_end=max(island_end),
        island_length=island_end-island_start,
        island_extended_length=island_extended_end[1]-island_extended_start[1],
        island_summit_pos=island_summit_pos[which.max(island_summit_abs)],
        island_summit_abs=island_summit_abs[which.max(island_summit_abs)],
        island_qvalue=island_qvalue[which.max(island_score)],
        island_pvalue=island_pvalue[which.max(island_score)],
        island_score=island_score[which.max(island_score)]) %>%
      dplyr::ungroup() %>%
      dplyr::filter(island_extended_length>=params$minlen) %>%
      dplyr::mutate(island_name=paste0("RDC_", stringr::str_pad(1:dplyr::n(), 3, pad="0")))

    writeLines(paste0("Extending islands reduced to ", nrow(islands_extended_ranges), " initial islands (", paste0(names(table(islands_extended_ranges$island_strand)), ":", table(islands_extended_ranges$island_strand), collapse=","), ")"))
    } else {
      islands_results_df = islands_df %>%
        dplyr::filter(island_length>=params$minlen) %>%
        dplyr::mutate(island_name=paste0("RDC_", stringr::str_pad(1:dplyr::n(), 3, pad="0")))
    }
  } else {
    islands_results_df = islands_df
    islands_results_df$island_name = character()
    islands_results_df$island_summit_abs = double()
    islands_results_df$island_summit_pos = double()
    islands_results_df$island_score = double()
    islands_results_df$island_qvalue = double()
    islands_results_df$island_pvalue = double()
    islands_results_df$island_length = double()
    if(extended_islands) {
      islands_results_df$island_extended_start = double()
      islands_results_df$island_extended_end = double()
      islands_results_df$island_extended_length = double()
    }
  }

  islands_results_reordered_df = islands_results_df %>%
    dplyr::select(island_name, island_chrom, island_start, island_end, dplyr::matches("island_summit_abs"), dplyr::matches("island_extended_"), dplyr::matches("island_"))
  list(qvalues=qvalues_df, islands=islands_results_reordered_df)
}
