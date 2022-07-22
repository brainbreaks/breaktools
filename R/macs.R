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

macs2_coverage_bgmodel = function(coverage_ranges, mask_ranges, debug_plots=F)
{
  coverage_df = as.data.frame(coverage_ranges) %>%
    dplyr::select(coverage_chrom=seqnames, coverage_start=start, coverage_end=end, coverage_strand=strand, coverage_score=score)

  mask_ranges = coverage_find_empty_intervals(coverage_ranges, coverage_column="score", minlen=1e3, mincoverage=0)

  # as.data.frame(mask_ranges) %>% dplyr::select(mask_chrom, mask_start, mask_end) %>% readr::write_tsv("reports/detect_rdc/mask.bed", col_names=F)

  bgmodel_data_df = suppressWarnings(coverage_ranges %>%
    ranges_sample(column="score", mask_ranges=mask_ranges, ntile=1e5) %>%
    as.data.frame() %>%
    dplyr::select(coverage_chrom=seqnames, coverage_start=start, coverage_end=end, coverage_strand=strand, coverage_score=score) %>%
    df2ranges(coverage_chrom, coverage_start, coverage_end, coverage_strand) %>%
    leftJoinByOverlaps(mask_ranges) %>%
    dplyr::filter(is.na(mask_start)) %>%
    dplyr::select(dplyr::matches("^coverage_")))
  table(bgmodel_data_df$coverage_chrom, bgmodel_data_df$coverage_score>0)
  bgmodel_quantiles_df = bgmodel_data_df %>%
    dplyr::group_by(coverage_chrom, coverage_strand) %>%
    dplyr::summarize(coverage_score.q999=quantile(coverage_score[coverage_score>0], 0.999), coverage_score.q500=quantile(coverage_score[coverage_score>0], 0.5), dataset="Observed data") %>%
    data.frame()

  bgmodel_data_df %>%
    dplyr::filter(coverage_chrom!="chrY") %>%
    dplyr::group_by(coverage_chrom) %>%
    dplyr::summarize(prop=mean(coverage_score>0))

  bgmodel_df = suppressWarnings(bgmodel_data_df %>%
    dplyr::filter(coverage_chrom!="chrY") %>%
    dplyr::group_by(coverage_chrom, coverage_strand) %>%
    dplyr::do((function(z) {
      yy<<-z
      asdasd()
      # z = bgmodel_data_df %>% dplyr::filter(sample_chrom=="chr8")

      fixargs = NULL
      probs = c(pmin(0.1, mean(z$coverage_score==0)), pmax(1-mean(z$coverage_score>=3), 0.9))
      range = quantile(z$coverage_score, probs)
      # , lower=c(1.99, 0), start=list(mu=pmax(2, mean(z$coverage_score.q500)), size=2)
      fit_gamma = fitdistrplus::fitdist(z$coverage_score, distr=distr, fix.arg=fixargs, method="qme", probs=probs)
      # if(fit_gamma$estimate["mu"]<=2) {
      #   fixargs = list(mu=2)
      #   # fit_gamma = fitdistrplus::fitdist(z$coverage_score, distr=distr, method="mle")
      #   fit_gamma = fitdistrplus::fitdist(z$coverage_score, distr=distr, fix.arg=fixargs, method="mle")
      # }
      dnbinom(3, size=fit_gamma$estimate["size"], mu=fit_gamma$estimate["mu"])
      plot(density(z$coverage_score[z$coverage_score>2]))

      # ggplot2::ggplot(z, aes(x=coverage_score, y=1)) +
      #   ggridges::stat_density_ridges(quantile_lines=T, calc_ecdf=T, geom="density_ridges_gradient", quantiles=probs, bandwidth=1, n=10000) +
      #   scale_x_continuous(limits = c(0, z$coverage_score.q500*10))

      # fit_gamma = MASS::fitdistr(z$coverage_score, densfun="gamma")
      if(!is.null(fixargs) & length(fixargs)>0) {
        fit_gamma$estimate = unlist(c(fit_gamma$estimate, fixargs[setdiff(names(fixargs), names(fit_gamma$estimate))]))
      }

      cbind(z[1,"coverage_chrom", drop=F], as.data.frame(t(fit_gamma$estimate)) %>% setNames(paste0("bgmodel_", colnames(.))))
    })(.)) %>%
    dplyr::ungroup())
  bgmodel_df %>%
    reshape2::melt(measure.vars=c("bgmodel_size", "bgmodel_mu")) %>%
    reshape2::dcast(variable+coverage_chrom ~ coverage_strand, value.var="value") %>%
    ggplot() +
    geom_point(aes(x=`+`, y=`-`)) +
      facet_wrap(~variable, scales="free")

  writeLines("Detected bgmodel is \n==============================================\n")
  writeLines(knitr::kable(bgmodel_df %>% dplyr::mutate(dplyr::across(dplyr::matches("bgmodel_"), round, 4))))

  if(debug_plots==T) {
    bgmodel_quant_df = bgmodel_df %>%
      dplyr::group_by(coverage_chrom, coverage_strand) %>%
      dplyr::summarize(coverage_quantile=seq(0, 1, length.out=500), coverage_density=distr_call(ver="q", distr=distr, x=coverage_quantile, params=dplyr::cur_data_all()))
    bgmodel_data_quant_df = bgmodel_data_df %>%
      dplyr::group_by(coverage_chrom, coverage_strand) %>%
      dplyr::summarize(coverage_quantile=seq(0, 1, length.out=500), coverage_density=quantile(coverage_score, coverage_quantile))
    bgmodel_qplot_df = bgmodel_quant_df %>%
      dplyr::inner_join(bgmodel_data_quant_df, by=c("coverage_chrom", "coverage_strand", "coverage_quantile"))
    print(ggplot2::ggplot(bgmodel_qplot_df) +
      ggplot2::geom_hex(ggplot2::aes(x=coverage_density.x, y=coverage_density.y), alpha=0.5, shape=1, size=2) +
      ggplot2::facet_wrap(~coverage_chrom+coverage_strand, scales="free") +
      ggplot2::geom_abline(intercept=0, slope=1) +
      ggplot2::labs(x="Fitted model", y="Observed dataodel", title=paste0("Q-Q plot -- ", coverage_df$tlx_group[1], ", ", distr, " distribution")))

    bgmodel_rand_df = bgmodel_df %>%
      dplyr::inner_join(bgmodel_quantiles_df %>% dplyr::distinct(coverage_chrom, coverage_strand), by=c("coverage_chrom", "coverage_strand")) %>%
      dplyr::inner_join(bgmodel_data_df %>% dplyr::distinct(coverage_chrom, coverage_strand), by=c("coverage_chrom", "coverage_strand")) %>%
      dplyr::group_by(coverage_chrom, coverage_strand) %>%
      dplyr::summarize(dataset="Fitted model", coverage_score=round(seq(0, coverage_score.q500*10, length.out=10000)), coverage_density=distr_call(ver="d", distr=distr, x=coverage_score, params=dplyr::cur_data_all()))
    bgmodel_cutoff_df = bgmodel_df %>%
      dplyr::group_by(coverage_chrom) %>%
      dplyr::summarize(cutoff_score=distr_call("q", distr, cutoff, params=dplyr::cur_data_all(), lower.tail=distr_lower.tail), dataset=paste0("Sign. cutoff (p<=", round(cutoff, 3), ")"))
    print(ggplot2::ggplot(mapping=ggplot2::aes(x=coverage_score, color=coverage_strand, linetype=dataset)) +
      ggplot2::geom_density(bw=1, data=bgmodel_data_df %>% dplyr::filter(coverage_score<coverage_score.q500*10) %>% dplyr::coverage_n(100000, replace=T)) +
      ggplot2::geom_line(ggplot2::aes(y=coverage_density), data=bgmodel_rand_df) +
      ggplot2::geom_vline(ggplot2::aes(xintercept=cutoff_score, color=dataset), data=bgmodel_cutoff_df) +
      ggplot2::facet_wrap(~coverage_chrom, scales="free") +
      ggplot2::labs(x="Pileup", y="Density", title=paste0("Pileup density plot -- ", coverage_df$tlx_group[1], ", ", distr, " distribution")))
  }

  bgmodel_df
}

macs2_coverage = function(sample_ranges, control_ranges=NULL, params, bgmodel_df, tmp_prefix=NULL, extended_islands=F, debug_plots=F)
{
  if(is.null(tmp_prefix)) {
    tmp_prefix = file.path("tmp", basename(tempfile()))
  }

  cutoff = ifelse(!is.na(params$minqvalue), params$minqvalue, params$minpvalue)
  sample_path = paste0(tmp_prefix, "-sample.bdg")
  sample_df = as.data.frame(sample_ranges) %>% dplyr::mutate(sample_chrom=seqnames, sample_start=start, sample_end=end, sample_strand=strand, sample_score=score)
  sample_ranges = sample_df %>% GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=T)
  # readr::write_tsv(sample_df %>% dplyr::select(sample_chrom, sample_start, sample_end, sample_score), file=sample_path, col_names=F)

  # distr = "nbinom"
  # distr = "gamma"
  # distr = "weibull"
  distr_lower.tail = F

  if(0 + !is.na(params$minqvalue) + !is.na(params$minpvalue) != 1) {
    stop("Please provide either minimal q-value or p-value")
  }

  control_path = paste0(tmp_prefix, "-control.bdg")
  qvalue_path = gsub(".bdg$", "_qvalue.bdg", sample_path)
  peaks_path = gsub(".bdg$", ".peaks", sample_path)
  bedpeaks_path = gsub(".bdg$", "_peaks.bed", sample_path)

  if(is.null(control_ranges) | length(control_ranges)==0) {
    control_df = sample_df %>%
      dplyr::left_join(bgmodel_df, by=c("sample_chrom"="bgmodel_chrom", "sample_strand"="bgmodel_strand")) %>%
      # dplyr::mutate(score=tidyr::replace_na(sample_bgmodel, 0), control_score=score) %>%
      dplyr::select(control_chrom=sample_chrom, control_start=sample_start, control_end=sample_end, control_strand=sample_strand, dplyr::matches("bgmodel_"))
    control_ranges = control_df %>% df2ranges(control_chrom, control_start, control_end, control_strand)
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
    pvalues_df = sample_df %>%
      dplyr::inner_join(control_df, by=c("sample_chrom"="control_chrom", "sample_start"="control_start", "sample_end"="control_end", "sample_strand"="control_strand")) %>%
      dplyr::group_by_at(colnames(control_df)[grepl("bgmodel_", colnames(control_df))]) %>%
      dplyr::mutate(pvalue=distr_call("p", distr=bgmodel_distr[1], x=sample_score, params=dplyr::cur_data_all(), lower.tail=distr_lower.tail), bgmodel_signal=distr_call("q", distr=bgmodel_distr[1], x=cutoff, params=dplyr::cur_data_all(), lower.tail=distr_lower.tail)) %>%
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
      dplyr::select(qvalue_chrom=seqnames, qvalue_start=start, qvalue_end=end, qvalue_strand=strand, qvalue_qvalue, qvalue_pvalue, dplyr::starts_with("bgmodel_"))
    if(!is.na(params$minqvalue)) {
      qvalues_df$qvalue_score = qvalues_df$qvalue_qvalue
    } else {
      qvalues_df$qvalue_score = qvalues_df$qvalue_pvalue
    }
    if(any(is.na(qvalues_df$qvalue_score))) {
      qvalues_df$qvalue_score[is.na(qvalues_df$qvalue_score)] = 0
    }

    x = qvalues_df %>%
      dplyr::mutate(qvalue_chrom=paste0(qvalue_chrom, ";", qvalue_strand)) %>%
      dplyr::select(qvalue_chrom, qvalue_start, qvalue_end, qvalue_score) %>%
      readr::write_tsv(file=qvalue_path, col_names=F)
  }


  # islands_df = qvalues_df %>%
  #   dplyr::filter(qvalue_score>=-log10(cutoff)) %>%
  #   dplyr::group_by(qvalue_chrom, qvalue_strand) %>%
  #   dplyr::do(GenomicRanges::reduce(df2ranges(., qvalue_chrom, qvalue_start, qvalue_end), min.gapwidth=params$seedgap) %>% as.data.frame()) %>%
  #   dplyr::ungroup() %>%
  #   dplyr::select(island_chrom=seqnames, island_strand=strand, island_start=start, island_end=end, island_length=width) %>%
  #   dplyr::mutate(island_start=island_start-1) %>%
  #   df2ranges(island_chrom, island_start, island_end) %>%
  #   innerJoinByOverlaps(qvalues_df %>% dplyr::rename(qvalue_tlx_group=tlx_group) %>% df2ranges(qvalue_chrom, qvalue_start, qvalue_end)) %>%
  #   dplyr::filter(island_strand==qvalue_strand) %>%
  #   dplyr::group_by(island_chrom, island_start, island_end, island_strand, island_length) %>%
  #   dplyr::summarize(island_summit_pos=round(qvalue_end/2+qvalue_start/2)[which.max(qvalue_score)], island_summit_abs=max(qvalue_score))

  #
  # Old MACS3 usage code
  #
  cmd_bdgpeakcall = stringr::str_glue("macs3 bdgpeakcall -i {qvalue} -c {cutoff} --min-length {format(minlen, scientific=F)} --max-gap {format(maxgap, scientific=F)} -o {output}",
     qvalue=qvalue_path, output=peaks_path, cutoff=-log10(cutoff), maxgap=params$seedgap, minlen=params$seedlen)
  writeLines(cmd_bdgpeakcall)
  system(cmd_bdgpeakcall)
  peaks_cols = cols(
    island_chrom=col_character(), island_start=col_double(), island_end=col_double(), island_peak=col_character(), island_summit_abs=col_double(),
    island_score=col_character(), island_fc=col_double(), island_pvalue_log10=col_double(), island_qvalue_log10=col_double(), island_sammit_offset=col_double()
  )
  islands_df = readr::read_tsv(peaks_path, skip=1, col_names=names(peaks_cols$cols), col_types=peaks_cols) %>%
    dplyr::mutate(island_length=island_end-island_start, island_summit_pos=island_start + island_sammit_offset) %>%
    dplyr::mutate(island_strand=gsub(".*;", "", island_chrom), island_chrom=gsub(";.*", "", island_chrom)) %>%
    dplyr::select(island_chrom, island_start, island_end, island_strand, island_summit_abs, island_sammit_offset, island_length, island_summit_pos) %>%
    dplyr::filter(island_length>=params$seedlen)


  #
  # Calculate extended island coordinates and reduce overlapping RDC to single instance
  #
  if(nrow(islands_df)>0) {
    if(extended_islands) {
      qvalues_ranges = qvalues_df %>% df2ranges(qvalue_chrom, qvalue_start, qvalue_end, qvalue_strand)
      island_qvalues_df = islands_results_df = islands_df %>%
        # dplyr::filter(island_chrom=="chr6" & island_start>=71199394 & island_end<=78418469) %>%
        df2ranges(island_chrom, island_start-1e6, island_end+1e6, island_strand) %>%
        innerJoinByOverlaps(qvalues_ranges) %>%
        dplyr::group_by_at(colnames(islands_df)) %>%
        dplyr::mutate_at(dplyr::vars(dplyr::starts_with("bgmodel_")), mean) %>%
        dplyr::mutate(island_qvalue=max(qvalue_qvalue), island_pvalue=max(qvalue_pvalue)) %>%
        dplyr::mutate(island_score_mean=weighted.mean(!is.na(qvalue_score) & qvalue_score>=-log10(cutoff)), qvalue_end-qvalue_start)

      islands_preresults_df = island_qvalues_df %>%
        dplyr::group_by_at(colnames(islands_df)) %>%
        dplyr::do((function(y){
          yy<<-y
          # y = islands_results_df %>% dplyr::filter(island_chrom=="chrX" & island_start==129965403)
          # y = islands_results_df %>% dplyr::filter(island_chrom=="chr9" & island_start==124232867)
          # y = islands_results_df %>% dplyr::filter(island_chrom=="chr9" & island_start==3580254)
          # y = island_qvalues_df %>% dplyr::filter(island_chrom=="chr6" & island_start==36742609)
          y_island = y %>% dplyr::select(island_chrom, island_start, island_end, island_qvalue) %>% dplyr::slice(1)
          y_long = y %>%
            reshape2::melt(measure.vars=c("qvalue_start", "qvalue_end"), value.name="qvalue_pos") %>%
            dplyr::arrange(qvalue_pos)
          y_model = loess(qvalue_pvalue~qvalue_pos, data=y_long, span=0.2)
          y_smooth = y_long %>%
            dplyr::distinct(island_chrom, island_start, island_end) %>%
            tidyr::crossing(qvalue_pos=seq(min(y_long$qvalue_pos), max(y_long$qvalue_pos), length.out=1000)) %>%
            dplyr::mutate(qvalue_pvalue = predict(y_model, .)) %>%
            dplyr::arrange(qvalue_pos) %>%
            dplyr::mutate(qvalue_start=qvalue_pos, qvalue_end=dplyr::lead(qvalue_pos, 1)) %>%
            dplyr::arrange(qvalue_start) %>%
            dplyr::mutate(qvalue_pvalue.diff=(qvalue_pvalue-dplyr::lag(qvalue_pvalue, 1))/(qvalue_end-qvalue_start)) %>%
            dplyr::filter(!is.na(qvalue_pvalue) & !is.na(qvalue_end))

          y_extended = y_smooth %>%
            dplyr::mutate(island_mid=(island_end+island_start)/2) %>%
            dplyr::filter(
              qvalue_pvalue >= -log10(cutoff) |
              qvalue_pvalue >= -log10(0.1) & (qvalue_start <= island_start & qvalue_pvalue.diff>=1e-6 | qvalue_end >= island_end & qvalue_pvalue.diff<=-1e-6)
            ) %>%
            dplyr::rename(extended_chrom="island_chrom", extended_start="qvalue_start", extended_end="qvalue_end") %>%
            df2ranges(extended_chrom, extended_start, extended_end) %>%
            # GenomicRanges::reduce() %>%
            # dplyr::filter(width>=params$seedlen) %>%
            GenomicRanges::reduce(min.gapwidth=params$maxgap) %>%
            as.data.frame() %>%
            dplyr::rename(extended_chrom="seqnames", extended_start="start", extended_end="end") %>%
            # dplyr::mutate(island_chrom=y$island_chrom[1], island_start=y$island_start[1], island_end=y$island_end[1])
            df2ranges(extended_chrom, extended_start, extended_end) %>%
            innerJoinByOverlaps(y_long %>% dplyr::distinct(island_chrom, island_start, island_end) %>% df2ranges(island_chrom, island_start, island_end))

          res_df = data.frame(island_extended_start=min(c(y$island_start, y_extended$extended_start)), island_extended_end=max(c(y$island_end, y_extended$extended_end)))
          #
          # p_min = min(y_long$qvalue_pvalue)
          # p_range = diff(range(y_long$qvalue_pvalue))
          # p = ggplot(y_long) +
          #     # geom_hline(yintercept=-log10(0.5), color="#CCCCCC") +
          #     geom_line(aes(x=qvalue_pos, y=qvalue_pvalue, color=island_strand)) +
          #     geom_line(aes(x=qvalue_pos, y=qvalue_pvalue, color="smooth fit"), data=y_smooth %>% dplyr::mutate(qvalue_pvalue=pmax(qvalue_pvalue, 0))) +
          #     # geom_segment(aes(x=extended_start, xend=extended_end, color="extended2"), y=p_min+p_range*0.01, yend=p_min+p_range*0.01, data=y_extended2) +
          #     geom_segment(aes(x=extended_start, xend=extended_end, color="extended"), y=p_min+p_range*0.02, yend=p_min+p_range*0.02, data=y_extended) +
          #     geom_segment(aes(x=island_start, xend=island_end, color="island"), y=p_min+p_range*0.035, yend=p_min+p_range*0.035, data=y_island) +
          #     geom_segment(aes(x=island_extended_start, xend=island_extended_end), y=p_min+p_range*0.05, yend=p_min+p_range*0.05, data=res_df) +
          #     labs(title=paste0(y$island_chrom[1], ":", y$island_start[1], "-", y$island_end[1]))
          # print(p)

          cbind(y_island, res_df)
        })(.)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(island_extended_length=island_extended_end-island_extended_start) %>%
        dplyr::filter(island_extended_length>=params$minlen) %>%
        dplyr::distinct_at(colnames(islands_df), .keep_all=T) %>%
        df2ranges(island_chrom, island_start, island_end, island_strand) %>%
        innerJoinByOverlaps(sample_ranges) %>%
        dplyr::group_by_at(colnames(islands_df)) %>%
        dplyr::mutate(island_sammit_pos=round(sample_start[which.max(sample_score)]/2+sample_end[which.max(sample_score)]/2), island_summit_abs=max(sample_score)) %>%
        dplyr::ungroup() %>%
        dplyr::distinct_at(colnames(islands_df), .keep_all=T) %>%
        dplyr::select(dplyr::matches("^(island_|bgmodel_)"))
    islands_preresults_ranges = islands_preresults_df %>%
      df2ranges(island_chrom, island_extended_start, island_extended_end, island_strand)
    islands_results_df = islands_preresults_ranges %>%
      GenomicRanges::reduce() %>%
      as.data.frame() %>%
      dplyr::rename(reduced_chrom="seqnames", reduced_start="start", reduced_end="end", reduced_strand="strand") %>%
      df2ranges(reduced_chrom, reduced_start, reduced_end, reduced_strand) %>%
      innerJoinByOverlaps(islands_preresults_ranges) %>%
      dplyr::group_by(island_chrom=reduced_chrom, island_extended_start=reduced_start, island_extended_end=reduced_end, island_strand=reduced_strand) %>%
      dplyr::summarize(
        island_start=min(island_start),
        island_end=max(island_end),
        island_length=island_end-island_start,
        island_extended_length=island_extended_end[1]-island_extended_start[1],
        # island_sammit_offset=island_sammit_offset[which.max(island_qvalue)],
        island_summit_abs=island_summit_abs[which.max(island_qvalue)],
        # island_summit_pos=island_summit_pos[which.max(island_qvalue)],
        island_qvalue=island_qvalue[which.max(island_qvalue)]) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(island_name=paste0("RDC_", stringr::str_pad(1:dplyr::n(), 3, pad="0")))
    } else {
      islands_results_df = islands_df %>% dplyr::mutate(island_name=paste0("RDC_", stringr::str_pad(1:dplyr::n(), 3, pad="0")))
    }
  } else {
    islands_results_df = islands_df
    islands_results_df$island_name = character()
    islands_results_df$island_summit_abs = double()
    islands_results_df$island_summit_pos = double()
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
