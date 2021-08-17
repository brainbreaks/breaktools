library(stringr)
library(dplyr)
library(readr)
library(bedr)
library(digest)


blank_tibble = function(cols) {
  stopifnot(class(tlx_cols)=="col_spec")

  x = tibble::tibble()
  for(n in names(cols$cols)) {
    n_class = class(cols$cols[[n]])
    n_type = match.fun(gsub("collector_", "", n_class[grepl("collector_", n_class)]))
    x[[n]] = n_type()
  }

  x
}

macs_cols = cols(
  macs_chrom=col_character(), macs_start=col_double(), macs_end=col_double(), macs_length=col_character(), macs_summit_abs=col_double(),
  macs_pileup=col_double(), macs_pvalue=col_double(), macs_fc=col_double(), macs_qvalue=col_double(), macs_name=col_character(), macs_comment=col_character()
)

macs_blank = function() {
  blank_tibble(macs_cols) %>% dplyr::mutate(macs_sample=NA_character_, macs_group=NA_character_)
}

bed_read = function(path) {
  bed = rtracklayer::import.bed(path)
  GenomicRanges::start(bed) = GenomicRanges::start(bed)-1
  bed
}

# one-based index
get_seq = function(fasta, ranges) {
  bed_df = data.frame(
    chr=as.character(GenomicRanges::seqnames(ranges)),
    start=as.numeric(GenomicRanges::start(ranges))-1,
    end=as.numeric(GenomicRanges::end(ranges)),
    strand=as.character(GenomicRanges::strand(ranges))) %>%
    dplyr::mutate(name="", score=0) %>%
    dplyr::select(chr, start, end, name, score, strand)
  bed_df.order = with(bed_df, order(chr, start, end, strand))
  res = bedr::get.fasta(bed_df[bed_df.order,], fasta=fasta, strand=T, check.zero.based=F, check.chr=F, check.valid=F, check.sort=T, check.merge=F)
  ranges$sequence = res$sequence[match(1:nrow(bed_df), bed_df.order)]

  ranges
}

join_offtarget2bait = function(offtargets_df, baits_df, genome_path) {
  baits_ranges = GenomicRanges::makeGRangesFromDataFrame(baits_df %>% dplyr::mutate(seqnames=bait_chrom, start=bait_start, end=bait_end, strand=bait_strand))
  offtargets_ranges = GenomicRanges::makeGRangesFromDataFrame(offtargets_df %>% dplyr::mutate(seqnames=offtarget_chrom, start=offtarget_start, end=offtarget_end, strand=offtarget_strand))

  # Combine baits and offtarget ranges so that sequences can be retrieve in a single call to bedtools (performance optimization)
  ranges_seq = get_seq(fasta=genome_path, ranges=BiocGenerics::append(baits_ranges, offtargets_ranges))
  baits_df$bait_sequence = ranges_seq$sequence[1:length(baits_ranges)]
  offtargets_df$offtarget_sequence = ranges_seq$sequence[-(1:length(baits_ranges))]

  offtarget2bait_df = offtargets_df %>%
    tidyr::crossing(baits_df) %>%
    dplyr::mutate(bait2offtarget_alignment=Biostrings::pairwiseAlignment(offtarget_sequence, bait_sequence, type="global", scoreOnly=T)) %>%
    dplyr::arrange(bait2offtarget_alignment) %>%
    dplyr::distinct(bait_sample, offtarget_chrom, offtarget_start, offtarget_end, offtarget_strand, .keep_all=T) %>%
    data.frame()

  offtarget2bait_df
}

offtargets_read = function(path) {
  offtargets_ranges = bed_read(path)
  as.data.frame(offtargets_ranges) %>%
    dplyr::select(offtarget_chrom=seqnames, offtarget_start=start, offtarget_end=end, offtarget_strand=strand)
}