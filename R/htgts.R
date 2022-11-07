#' @export
htgts_calculate_positions = function(sequences_df, database_path, primer_sequence_column="Primer", sgRNA_sequence_column="Sequence", PAM_sequence_column="Sequence_PAM") {
  sequences_cols = colnames(sequences_df)
  sequences_df = sequences_df %>% data.frame() %>% dplyr::mutate(sequence_latent_id=1:dplyr::n())
  sequences_long_df = sequences_df %>%
    dplyr::rename(Sequence=sgRNA_sequence_column, Sequence_PAM=PAM_sequence_column, RED_primer=primer_sequence_column) %>%
    dplyr::mutate(Sequence_and_PAM=paste0(Sequence, Sequence_PAM))  %>%
    reshape2::melt(measure.vars=c("RED_primer", "Sequence_and_PAM"), value.name="Seq")

  blat_df = get_blat(unique(sequences_long_df$Seq), fasta=database_path)
  blat_identical_df = blat_df %>% dplyr::filter(blat_mismatch==0 & blat_target_gap_bases==0)
  sequences_long2blat_df = sequences_long_df %>%
    dplyr::inner_join(blat_identical_df %>% dplyr::select(blat_sequence, blat_target_name, blat_target_start, blat_target_end, blat_strand, blat_mismatch, blat_target_gap_bases), by=c("Seq"="blat_sequence")) %>%
    dplyr::inner_join(blat_identical_df %>% dplyr::select(blat_sequence.other=blat_sequence, blat_target_name.other=blat_target_name, blat_target_start.other=blat_target_start), by=c("blat_target_name"="blat_target_name.other")) %>%
    dplyr::filter(blat_sequence.other!=Seq) %>%
    dplyr::arrange(sequence_latent_id, variable, abs(blat_target_start-blat_target_start.other)) %>%
    dplyr::distinct(sequence_latent_id, variable, .keep_all=T)

  sequences_result_df = sequences_df %>%
    dplyr::inner_join(sequences_long2blat_df %>% dplyr::mutate(variable=paste0(variable, ".start")) %>% reshape2::dcast(sequence_latent_id ~ variable, value.var="blat_target_start"), by="sequence_latent_id") %>%
    dplyr::inner_join(sequences_long2blat_df %>% dplyr::mutate(variable=paste0(variable, ".end")) %>% reshape2::dcast(sequence_latent_id ~ variable, value.var="blat_target_end"), by="sequence_latent_id") %>%
    dplyr::inner_join(sequences_long2blat_df %>% dplyr::mutate(variable=paste0(variable, ".strand")) %>% reshape2::dcast(sequence_latent_id ~ variable, value.var="blat_strand"), by="sequence_latent_id") %>%
    dplyr::mutate(Cut_position=ifelse(Sequence_and_PAM.strand=="+", Sequence_and_PAM.end-nchar(Sequence_PAM)-3, Sequence_and_PAM.start+nchar(Sequence_PAM)+2)) %>%
    dplyr::mutate(Start_calculated=ifelse(RED_primer.strand=="+", RED_primer.start, Cut_position+1),
                  End_calculated=ifelse(RED_primer.strand=="+", Cut_position, RED_primer.end+1)) %>%
    data.frame()


  sequences_result_final_df = (sequences_result_df %>% dplyr::mutate(Start=Start_calculated, End=End_calculated, Strand=RED_primer.strand, SequenceStrand=Sequence_and_PAM.strand))[,unique(c(sequences_cols, "Start", "End", "Strand", "SequenceStrand"))]
  sequences_result_final_df
}

#' @export
htgts_barcodes_detect = function(fastq_paths, primers, max_sequences=10000) {
  if(length(fastq_paths) != length(primers)) {
    stop("Vectors of fasta files and primes should be of the same length")
  }
  # writeLines("Extracting number of reads from each sample")
  # fasta_count = ShortRead::countFastq(fastq_paths)$records
  fasta_count = rep(max_sequences, length(fastq_paths))

  barcodes_results_df = data.frame()
  for(i in 1:length(fastq_paths)) {
    writeLines(paste0(i, "/", length(fastq_paths)))

    fasta_stream = ShortRead::FastqSampler(fastq_paths[i], pmin(max_sequences, fasta_count[i]))
    fasta = ShortRead::yield(fasta_stream)
    fasta_reads = as.character(ShortRead::sread(fasta))
    fasta_reads_count = length(fasta_reads)
    fasta_unique_reads = unique(fasta_reads)
    fasta_unique_reads_count = length(fasta_unique_reads)
    close(fasta_stream)
    writeLines(paste0("Read ", fasta_reads_count, "/", fasta_count[i], " lines (unique: ", fasta_unique_reads_count, ") from file '", fastq_paths[i], "'"))
    fasta_reads = fasta_unique_reads

    # Align primer to reads
    primer_alignment = Biostrings::pairwiseAlignment(pattern=fasta_reads, subject=primers[i])
    primer_alignment_width = Biostrings::width(Biostrings::pattern(primer_alignment))
    primer_alignment_score = Biostrings::score(primer_alignment)
    primer_alignment_offset = Biostrings::start(Biostrings::pattern(primer_alignment))
    primer_present = primer_alignment_score <= -200 & primer_alignment_width <= nchar(primers[i])*1.3
    primer_alignment_offset = primer_alignment_offset[primer_present]
    fasta_aligned_reads = fasta_reads[primer_present]

    plot(sort(primer_alignment_width))
    plot(sort(primer_alignment_score))

    primer_offset_table = sort(table(primer_alignment_offset), decreasing=T)
    primer_offset = as.numeric(names(primer_offset_table)[1])
    fasta_offset_reads = fasta_aligned_reads[primer_alignment_offset==primer_offset]
    mid_all = substr(fasta_offset_reads, 0, stop=primer_offset-1)

    mid_table = sort(table(mid_all), decreasing=T)
    mid_percent = mid_table[1]/length(fasta_reads)
    barcodes_results_df.i = data.frame(barcode_fasta=fastq_paths[i], barcode_primer=primers[i], primer_percent=mean(primer_present), barcode_percent=mid_percent, barcode_sequence=names(mid_table)[1])
    barcodes_results_df = rbind(barcodes_results_df, barcodes_results_df.i)
  }
  barcodes_results_df
}