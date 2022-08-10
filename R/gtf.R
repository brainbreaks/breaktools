#' @export
gtf_read = function(path, dtype="gene", tmp_dir="tmp", max_age=NULL) {
  # Load gene annotations

  clear_tmpdir(tmp_dir=tmp_dir, max_age=max_age)
  results = list()
  if("transcript" %in% dtype) {
    transcript_cache = file.path(tmp_dir, "gtf_transcript.rda")
    if(file.exists(transcript_cache)) {
      load(transcript_cache)
    } else {
      genome_txdb = GenomicFeatures::makeTxDbFromGFF(path, format="gtf")
      transcript_df = as.data.frame(GenomicFeatures::transcripts(genome_txdb, columns=c("tx_id", "gene_id")))
      save(transcript_df, file=transcript_cache)
    }

    results[["transcript"]] = transcript_df %>%
      dplyr::mutate(promoter_chrom=seqnames, promoter_start=start, promoter_end=end, promoter_strand=strand) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(tx_id=paste0(unique(unlist(tx_id)), collapse=","), gene_id=paste0(unique(unlist(gene_id)), collapse=","))
  }
  if("gene" %in% dtype) {
    genes_cache = file.path(tmp_dir, "gtf_read.rda")
    if(file.exists(genes_cache)) {
      load(genes_cache)
    } else {
      genome_txdb = GenomicFeatures::makeTxDbFromGFF(path, format="gtf")
      genes_df = as.data.frame(GenomicFeatures::genes(genome_txdb))
      save(genes_df, file=genes_cache)
    }

    results[["gene"]] = genes_df %>%
      dplyr::mutate(gene_chrom=as.character(seqnames), gene_start=start, gene_end=end, gene_strand=as.character(strand), gene_length=gene_end-gene_start)

    #
    # Give a cluster position number to genes that overlap (for visualization)
    #
    genes_ranges = GenomicRanges::makeGRangesFromDataFrame(results[["gene"]], keep.extra.columns=T)
    genes_reduced_ranges = GenomicRanges::reduce(genes_ranges, ignore.strand=T)
    genes_reduced_ranges$gene_cluster = 1:length(genes_reduced_ranges)
    results[["gene"]] = as.data.frame(IRanges::mergeByOverlaps(genes_ranges, genes_reduced_ranges)) %>%
      dplyr::select(-dplyr::matches("_ranges\\.")) %>%
      dplyr::arrange(dplyr::desc(gene_length)) %>%
      dplyr::group_by(gene_cluster) %>%
      dplyr::mutate(gene_cluster_i=1:dplyr::n()) %>%
      dplyr::ungroup()
  }

  if(length(results)==0) {
    return(NULL)
  }

  if(length(results)==1) {
    return(results[[1]])
  }

  results
}