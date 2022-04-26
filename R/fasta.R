#' @export
fasta_duplicates_count = function(paths, breakpoints=c(1, 10, 50, 100, 1000, 10000), palette=NULL, plot=T) {
  duplicates_sumdf = data.frame(fasta_path=paths) %>%
    dplyr::mutate(fasta_name=basename(fasta_path)) %>%
    dplyr::rowwise() %>%
    dplyr::do((function(f){
      ff<<-f
      fasta = ShortRead::readFastq(f$fasta_path)
      sequences = as.character(ShortRead::sread(fasta))
      r = data.frame(sequence=sequences, fasta_path=f$fasta_path, fasta_name=f$fasta_name) %>%
        dplyr::group_by(fasta_path, fasta_name) %>%
        dplyr::mutate(fastq_size=dplyr::n()) %>%
        dplyr::group_by(fasta_path, fasta_name, sequence) %>%
        dplyr::summarize(sequence_duplicates=dplyr::n(), sequence_duplicates_prop=sequence_duplicates/fastq_size[1]) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(sequence_duplicates_group_ub=as.numeric(gsub(".*,([^]]+)]", "\\1", cut(sequence_duplicates, breaks=c(0, breakpoints))))) %>%
        dplyr::group_by(fasta_path, fasta_name, sequence_duplicates_group_ub) %>%
        dplyr::summarize(group_size=sum(sequence_duplicates), group_unique_size=length(unique(sequence)))
      return(r)
    })(.))%>%
    dplyr::ungroup() %>%
    dplyr::mutate(sequence_duplicates_group=dplyr::case_when(
      sequence_duplicates_group_ub==1 ~ "=1",
      is.na(sequence_duplicates_group_ub) ~ paste0(">",max(breakpoints)),
      T ~ paste0("<=", sequence_duplicates_group_ub))) %>%
    dplyr::arrange(sequence_duplicates_group_ub) %>%
    dplyr::mutate(sequence_duplicates_group=factor(sequence_duplicates_group, unique(sequence_duplicates_group))) %>%
    dplyr::mutate(fasta_name=gsub("(_R1|_R2)\\.(fasta|fastq|fq)(\\.gz)?$", "", fasta_name))

  if(plot==T) {
    if(is.null(palette)) {
      palette = RColorBrewer::brewer.pal(pmax(3, length(breakpoints)+1), "YlGnBu")[0:length(breakpoints)+1]
      names(palette) = c(ifelse(breakpoints==1, "=1", paste0("<=", breakpoints)), paste0(">",max(breakpoints)))
    }

    duplicates_ggplot = ggplot2::ggplot(duplicates_sumdf) +
      ggplot2::geom_bar(ggplot2::aes(x=fasta_name, y=group_size, fill=sequence_duplicates_group), stat="identity") +
      ggplot2::scale_fill_manual(values=palette) +
      ggplot2::scale_y_continuous(labels=scales::label_number(accuracy=0.1, scale=1e-3, suffix="k")) +
      ggplot2::theme_minimal() +
      ggplot2::labs(y="Reads", x="", fill="Number of duplicates") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1, vjust=1))

    return(duplicates_ggplot)
  } else {
    return(duplicates_sumdf)
  }
}

plot_logos_coordinates = function(fastq_paths, sample_names, widths=list("Beginning"=c(1,30), "End"=c(-25, -1)), max_sequences=100000) {
  plots = list()
  for(i in 1:length(fastq_paths)) {
    sample_path = fastq_paths[i]
    sample_name = sample_names[i]
    writeLines(paste0(i, "/", length(fastq_paths), " : ", sample_name, "    ", sample_path))

    fasta_stream = ShortRead::FastqSampler(sample_path, max_sequences)
    fasta = ShortRead::yield(fasta_stream)
    fasta_reads = as.character(ShortRead::sread(fasta))
    fasta_unique_reads = unique(fasta_reads)

    # fasta = ShortRead::readFastq(sample_path)
    # fasta_reads = ShortRead::sread(fasta)
    plist = list()
    plist[[1]] = cowplot::ggdraw() +
      cowplot::draw_label(sample_names[i], size=10)
    for(wname in names(widths)) {
      writeLines(paste0("  ", wname))

      # If width is longer for all or some sequences then extend the sequences to the needed length
      fasta_reads_long = fasta_unique_reads
      fasta_widths = Biostrings::width(fasta_reads_long)
      fasta_ends = ifelse(fasta_widths>widths[[wname]][2], widths[[wname]][2], fasta_widths)
      missing_nchar = widths[[wname]][2] - fasta_ends
      if(any(missing_nchar>0)) {
        fasta_missing = Biostrings::DNAStringSet(sapply(missing_nchar[missing_nchar>0], function(x) paste(replicate(x, expr="N"), collapse="")))
        fasta_reads_long[missing_nchar>0] = Biostrings::xscat(fasta_reads_long[missing_nchar>0], fasta_missing)
      }

      fasta_w = Biostrings::subseq(fasta_reads_long, start=widths[[wname]][1], end=widths[[wname]][2])
      plist[[length(plist)+1]] = ggplot2::ggplot() +
          ggseqlogo::geom_logo(as.character(fasta_w)) +
          ggplot2::labs(title=wname) +
          ggseqlogo::theme_logo() +
          ggplot2::guides(fill="none") +
          ggplot2::theme(axis.title.x=ggplot2::element_blank(), axis.text.x=ggplot2::element_blank(), axis.ticks.x=ggplot2::element_blank(), axis.title.y=ggplot2::element_blank(), axis.text.y=ggplot2::element_blank(), axis.ticks.y=ggplot2::element_blank())
    }


    defaultW = getOption("warn")
    options(warn = -1)
    p = cowplot::plot_grid(plotlist=plist, ncol=1+length(widths), rel_widths=c(1,3,3))
    options(warn = defaultW)
    plots[[sample_name]] = p
  }
  plots
}