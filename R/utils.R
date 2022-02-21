blank_tibble = function(cols) {
  stopifnot(class(cols)=="col_spec")

  x = tibble::tibble()
  for(n in names(cols$cols)) {
    n_class = class(cols$cols[[n]])
    n_type = match.fun(gsub("collector_", "", n_class[grepl("collector_", n_class)]))
    x[[n]] = n_type()
  }

  x
}

#' @export
blast_columns = function() {
  readr::cols(
    blast_qseqid=readr::col_character(),
    blast_sseqid=readr::col_character(),
    blast_pident=readr::col_integer(),
    blast_length=readr::col_integer(),
    blast_mismatch=readr::col_integer(),
    blast_gapopen=readr::col_integer(),
    blast_qstart=readr::col_integer(),
    blast_qend=readr::col_integer(),
    blast_sstart=readr::col_integer(),
    blast_send=readr::col_integer(),
    blast_evalue=readr::col_double(),
    blast_bitscore=readr::col_double()
  )
}

#' @export
blast_blank = function() {
  blast_columns(blast_columns()) %>% dplyr::mutate(blast_sequence="")
}

#' @export
blat_columns = function() {
  readr::cols(
    blat_match = readr::col_integer(),
    blat_mismatch = readr::col_integer(),
    blat_rep_match = readr::col_integer(),
    blat_Ns = readr::col_integer(),
    blat_query_gap_count = readr::col_integer(),
    blat_query_gap_bases = readr::col_integer(),
    blat_target_gap_count = readr::col_integer(),
    blat_target_gap_bases = readr::col_integer(),
    blat_strand = readr::col_character(),
    blat_query_name = readr::col_character(),
    blat_query_size = readr::col_integer(),
    blat_query_start = readr::col_integer(),
    blat_query_end = readr::col_integer(),
    blat_target_name = readr::col_character(),
    blat_target_size = readr::col_integer(),
    blat_target_start = readr::col_integer(),
    blat_target_end = readr::col_integer(),
    blat_block_count = readr::col_integer(),
    blat_blockSizes = readr::col_character(),
    blat_query_starts = readr::col_character(),
    blat_target_starts = readr::col_character()
  )
}

#' @export
blat_blank = function() {
  blank_tibble(blat_columns()) %>% dplyr::mutate(blat_sequence="")
}


#' @export
leftJoinByOverlaps = function(query, subject, ...) {
  query$query_id = 1:length(query)
  subject$subject_id = 1:length(subject)
  result_df = as.data.frame(IRanges::mergeByOverlaps(query, subject, ...))
  result_df = dplyr::bind_rows(result_df, as.data.frame(query) %>% dplyr::anti_join(result_df %>% dplyr::select(query_id), by="query_id"))
  result_df = result_df %>% dplyr::select(-dplyr::matches("^(query|subject)\\."))

  result_df
}

#' @export
liftOverRanges = function(ranges, chain) {
  chain = rtracklayer::import.chain("genomes/mm9/mm9ToMm10.over.chain")
  ranges$ranges_id = 1:length(ranges)
  as.data.frame(unlist(rtracklayer::liftOver(ranges, chain))) %>%
    dplyr::group_by(ranges_id) %>%
    dplyr::mutate(start=max(start), end=min(end)) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(ranges_id, .keep_all=T) %>%
    dplyr::select(-ranges_id) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=T)
}


#' @export
bed_read = function(path) {
  bed = rtracklayer::import.bed(path)
  GenomicRanges::start(bed) = GenomicRanges::start(bed)-1
  bed
}

#' @export
get_seq = function(fasta, ranges) {
  if(!bedr::check.binary("bedtools")) {
    stop("bedtools executible not found")
  }

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

#' @export
get_blat = function(sequences, database, minscore=20, stepSize=5, tmp_dir="tmp") {
  if(!bedr::check.binary("blat")) {
    stop("blat executible not found")
  }
  if(!bedr::check.binary("bowtie")) {
    stop("bowtie executible not found")
  }
  if(!bedr::check.binary("samtools")) {
    stop("samtools executible not found")
  }

  if(!file.exists(database)) { stop(paste0("Database file '", database, "' doesn't exist")) }

  dir.create(tmp_dir, recursive=T, showWarnings=F)

  fasta_path = file.path(tmp_dir, basename(tempfile()))
  psl_path = paste0(fasta_path, ".psl")
  sequences_df = data.frame(blat_query_name=paste0("seq", 1:length(sequences)), blat_sequence=sequences, stringsAsFactors=F)
  writeLines(paste0(">", sequences_df$blat_query_name, "\n", sequences_df$blat_sequence), con=fasta_path)

  cmd = paste0("blat -stepSize=", stepSize, " -repMatch=2253 -minScore=", minscore, " -minIdentity=0 -out=psl -dots=1 -noHead ", database, " ", fasta_path, " ", psl_path)
  writeLines(cmd)
  system(cmd)

  if(file.info(psl_path)$size > 0) {
    result = readr::read_tsv(psl_path, col_types=blat_columns(), col_names=names(blat_columns()$cols))
  } else {
    result = blat_blank()
  }

  result %>%
    dplyr::mutate(blat_target_start=blat_target_start+1) %>%
    dplyr::inner_join(sequences_df, by="blat_query_name")
}

#' @export
get_bowtie = function(sequences, database, threads=30, tmp_dir="tmp") {
  if(!bedr::check.binary("bowtie")) {
    stop("bowtie executible not found")
  }
  if(!bedr::check.binary("samtools")) {
    stop("samtools executible not found")
  }
  if(!file.exists(database)) { stop(paste0("Database file '", database, "' doesn't exist")) }

  database_name = gsub("\\.[^.]+", "", database)
  if(!file.exists(paste0(database_name, ".1.ebwt"))) {
    cmd_makedb = paste("bowtie-build --threads", threads, database, database_name)
    system(cmd_makedb)
  }

  dir.create(tmp_dir, recursive=T, showWarnings=F)
  fasta_path = file.path(tmp_dir, basename(tempfile()))
  sam_path = paste0(fasta_path, ".sam")
  bam_path = paste0(fasta_path, ".bam")
  sequences_df = data.frame(bowtie_qid=paste0("seq", 1:length(sequences)), bowtie_sequence=sequences, stringsAsFactors=F)
  sequences_reduced_df = sequences_df %>% dplyr::distinct(bowtie_sequence, .keep_all=T)
  writeLines(paste0(">", sequences_reduced_df$bowtie_qid, "\n", sequences_reduced_df$bowtie_sequence), con=fasta_path)

  cmd = paste0("bowtie --threads ", threads, " --tryhard --all -v 3 --seedlen 5 -f --sam --no-unal ", database_name, " ", fasta_path, " > ", sam_path)
  system(cmd)

  sort_cmd = paste0("samtools sort -@ ", threads, " ", sam_path, " > ", bam_path)
  system(sort_cmd)

  # Analyze genome position and give final output
  primers_alignments_param = Rsamtools::ScanBamParam(tag=c("AS", "NM"), what=c("qname", "rname", "strand", "flag", "pos", "qwidth",  "cigar", "mapq", "qual"))
  primers_alignments = lapply(Rsamtools::scanBam(bam_path, param=primers_alignments_param), function(z) {
    d = data.frame(z[names(z)!="tag"])
    d$bowtie_mismatches=z$tag$NM
    cbind(d, do.call(rbind, lapply(d$flag, SamSeq::samFlags))) })[[1]] %>%
    dplyr::mutate(bowtie_qid=qname, bowtie_end=pos+qwidth-1, bowtie_primary=!NOT_PRIMARY_ALIGNMENT) %>%
    dplyr::filter(!READ_UNMAPPED & rname!="chrM") %>%
    dplyr::select(bowtie_qid, bowtie_chrom=rname, bowtie_strand=strand, bowtie_start=pos, bowtie_end, bowtie_flag=flag, bowtie_cigar=cigar, bowtie_length=qwidth, bowtie_mismatches, bowtie_primary) %>%
    dplyr::mutate(bowtie_qid=as.character(bowtie_qid), bowtie_chrom=as.character(bowtie_chrom), bowtie_strand=as.character(bowtie_strand), bowtie_start=as.numeric(bowtie_start), bowtie_end=as.numeric(bowtie_end), bowtie_flag=as.numeric(bowtie_flag), bowtie_cigar=as.character(bowtie_cigar)) %>%
    dplyr::inner_join(sequences_reduced_df, by="bowtie_qid") %>%
    dplyr::select(-bowtie_qid)

  primers_alignments
}

#' @export
get_blast = function(sequences, database, word_size=4, perc_identity=100, tmp_dir="tmp") {
  if(!bedr::check.binary("blastn")) {
    stop("blast executible not found")
  }
  if(!file.exists(database)) { stop(paste0("Database file '", database, "' doesn't exist")) }

  database_name = gsub("\\.[^.]+", "", database)
  if(!file.exists(paste0(database_name, ".nsq"))) {
    cmd_makedb = paste0(paste0("makeblastdb -in ", database, " -out ", database_name, " -dbtype 'nucl' -hash_index"))
    system(cmd_makedb)
  }

  dir.create(tmp_dir, recursive=T, showWarnings=F)
  fasta_path = file.path(tmp_dir, basename(tempfile()))
  sequences_df = data.frame(blast_qseqid=paste0("seq", 1:length(sequences)), blast_sequence=sequences, stringsAsFactors=F)
  writeLines(paste0(">", sequences_df$blast_query_name, "\n", sequences_df$blast_sequence), con=fasta_path)

  out_path = paste0(fasta_path, ".blast")
  cmd = paste0("blastn -db ", database_name, " -query ", fasta_path, " -out ", out_path, " -word_size ", word_size, " -perc_identity ", perc_identity, " -outfmt 6")
  writeLines(cmd)
  system(cmd)

  if(file.info(out_path)$size > 0) {
    result = readr::read_tsv(out_path, col_types=blast_columns(), col_names=names(blast_columns()$cols))
  } else {
    result = blast_blank()
  }

  result %>%
    dplyr::mutate(blast_sstart=blast_sstart+1) %>%
    dplyr::inner_join(sequences_df, by="blast_qseqid")
}

#' @export
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

#' @export
offtargets_read = function(path) {
  offtargets_ranges = bed_read(path)
  as.data.frame(offtargets_ranges) %>%
    dplyr::select(offtarget_chrom=seqnames, offtarget_start=start, offtarget_end=end, offtarget_strand=strand)
}

#' @export
file_count_lines = function(path) {
  f = file(path, open="rb")
  nlines = 0L
  while (length(chunk <- readBin(f, "raw", 65536)) > 0) {
     nlines = nlines + sum(chunk == as.raw(10L))
  }
  close(f)
  nlines
}