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