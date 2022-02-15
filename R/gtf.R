#' @export
gtf_read = function(path, dtype="gene") {
  # Load gene annotations
  genome_txdb = GenomicFeatures::makeTxDbFromGFF(path, format="gtf")

  results = list()
  if("transcript" %in% dtype) {
    results[["transcript"]] = as.data.frame(GenomicFeatures::transcripts(genome_txdb, columns=c("tx_id", "gene_id"))) %>%
      dplyr::mutate(promoter_chrom=seqnames, promoter_start=start, promoter_end=end, promoter_strand=strand) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(tx_id=paste0(unique(unlist(tx_id)), collapse=","), gene_id=paste0(unique(unlist(gene_id)), collapse=","))
  }
  if("gene" %in% dtype) {
    results[["gene"]]  = as.data.frame(GenomicFeatures::genes(genome_txdb)) %>%
      dplyr::mutate(gene_chrom=as.character(seqnames), gene_start=start, gene_end=end, gene_strand=strand, gene_length=gene_end-gene_start)

    #
    # Give a cluster position number to genes that overlap (for visualization)
    #
    genes_ranges = GenomicRanges::makeGRangesFromDataFrame(results[["gene"]], keep.extra.columns=T)
    genes_reduced_ranges = GenomicRanges::reduce(genes_ranges, ignore.strand=T) %>% plyranges::mutate(gene_cluster= 1:plyranges::n())
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