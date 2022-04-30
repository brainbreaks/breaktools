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

blast_blank = function() {
  blast_columns(blast_columns()) %>% dplyr::mutate(blast_sequence="")
}

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

blat_blank = function() {
  blank_tibble(blat_columns()) %>% dplyr::mutate(blat_sequence="")
}

separate_target_subject_columns = function(query_ranges, subject_ranges) {
  query_mcols = GenomicRanges::mcols(query_ranges)
  subject_mcols = GenomicRanges::mcols(subject_ranges)
  intersected_columns = intersect(colnames(query_mcols), colnames(subject_mcols))
  if(length(intersected_columns)>0) {
    colnames(query_mcols)[match(intersected_columns, colnames(query_mcols))] = paste0(intersected_columns, ".x")
    colnames(subject_mcols)[match(intersected_columns, colnames(subject_mcols))] = paste0(intersected_columns, ".y")
    GenomicRanges::mcols(query_ranges) = query_mcols
    GenomicRanges::mcols(subject_ranges) = subject_mcols
  }
  list(query_ranges=query_ranges, subject_ranges=subject_ranges)
}


#' @export
#' @description Left-join query and subject ranges based on overlap between them
#'
#' @param query_ranges Query ranges
#' @param subject_ranges Subject ranges
#'
#' @return A data frame with two ranges objects left-joined based on the overlap between them
#' @examples
#' query_ranges = data.frame(query_chrom="chr1", query_start=1:2, query_end=1:2, col="AAA") %>% df2ranges(query_chrom, query_start, query_end)
#' subject_ranges = data.frame(subject_chrom="chr1", subject_start=2:3, subject_end=2:3, col="AAA") %>% df2ranges(subject_chrom, subject_start, subject_end)
#' leftJoinByOverlaps(query_ranges, subject_ranges)
leftJoinByOverlaps = function(query_ranges, subject_ranges) {
  r = separate_target_subject_columns(query_ranges, subject_ranges)
  query_ranges = r$query_ranges
  subject_ranges = r$subject_ranges

  query_ranges$query_id = 1:length(query_ranges)
  subject_ranges$subject_id = 1:length(subject_ranges)
  result_df = as.data.frame(IRanges::mergeByOverlaps(query_ranges, subject_ranges)) %>% dplyr::select(-dplyr::matches("^(query|subject)_ranges\\."))
  result_df = dplyr::bind_rows(result_df, as.data.frame(GenomicRanges::mcols(query_ranges)) %>% dplyr::anti_join(result_df %>% dplyr::select(query_id), by="query_id"))
  result_df = result_df %>% dplyr::select(-dplyr::matches("^(query|subject)_ranges\\."))

  result_df
}

#' @export
#' @description Full-join query and subject ranges based on overlap between them
#'
#' @param query_ranges Query ranges
#' @param subject_ranges Subject ranges
#'
#' @return A data frame with two ranges objects full-joined based on the overlap between them
#' @examples
#' query_ranges = data.frame(query_chrom="chr1", query_start=1:2, query_end=1:2, col="AAA") %>% df2ranges(query_chrom, query_start, query_end)
#' subject_ranges = data.frame(subject_chrom="chr1", subject_start=2:3, subject_end=2:3, col="AAA") %>% df2ranges(subject_chrom, subject_start, subject_end)
#' fullJoinByOverlaps(query_ranges, subject_ranges)
fullJoinByOverlaps = function(query_ranges, subject_ranges) {
  r = separate_target_subject_columns(query_ranges, subject_ranges)
  query_ranges = r$query_ranges
  subject_ranges = r$subject_ranges

  query_ranges$query_id = 1:length(query_ranges)
  subject_ranges$subject_id = 1:length(subject_ranges)
  result_df = as.data.frame(IRanges::mergeByOverlaps(query_ranges, subject_ranges))
  result_df = dplyr::bind_rows(
    result_df,
    as.data.frame(GenomicRanges::mcols(query_ranges)) %>% dplyr::anti_join(result_df %>% dplyr::select(query_id), by="query_id"),
    as.data.frame(GenomicRanges::mcols(subject_ranges)) %>% dplyr::anti_join(result_df %>% dplyr::select(subject_id), by="subject_id"))
  result_df = result_df %>% dplyr::select(-dplyr::matches("^(query|subject)_ranges\\."), -dplyr::matches("^(subject_id|query_id)$"))

  result_df
}

#' @export
#' @description Inner-join query and subject ranges based on overlap between them
#'
#' @param query_ranges Query ranges
#' @param subject_ranges Subject ranges
#'
#' @return A data frame with two ranges objects inner-joined based on the overlap between them
#' @examples
#' query_ranges = data.frame(query_chrom="chr1", query_start=1:2, query_end=1:2, col="AAA") %>% df2ranges(query_chrom, query_start, query_end)
#' subject_ranges = data.frame(subject_chrom="chr1", subject_start=2:3, subject_end=2:3, col="AAA") %>% df2ranges(subject_chrom, subject_start, subject_end)
#' innerJoinByOverlaps(query_ranges, subject_ranges)
innerJoinByOverlaps = function(query_ranges, subject_ranges) {
  r = separate_target_subject_columns(query_ranges, subject_ranges)
  query_ranges = r$query_ranges
  subject_ranges = r$subject_ranges

  result_ranges = IRanges::mergeByOverlaps(query_ranges, subject_ranges)
  as.data.frame(result_ranges) %>% dplyr::select(-dplyr::matches("_ranges\\."))
}

test = function()
{
  query_ranges = GenomicRanges::makeGRangesFromDataFrame(data.frame(seqnames="chr1", start=1, end=1, col="AAA"),  keep.extra.columns = T)
  subject_ranges = GenomicRanges::makeGRangesFromDataFrame(data.frame(seqnames="chr1", start=1, end=1, col="BBB"),  keep.extra.columns = T)
  ranges_list = list(query_ranges, subject_ranges)
  devtools::load_all('~/Workspace/breaktools/')
  innerJoinManyByOverlaps(ranges_list)

  devtools::load_all('~/Workspace/breaktools/')
  query_ranges = data.frame(query_chrom="chr1", query_start=1:2, query_end=1:2, col="AAA") %>% df2ranges(query_chrom, query_start, query_end)
  subject_ranges = data.frame(subject_chrom="chr1", subject_start=2:3, subject_end=2:3, col="AAA") %>% df2ranges(subject_chrom, subject_start, subject_end)
  fullJoinByOverlaps(query_ranges, subject_ranges)
}

#' @export
#' @description Convert data.frame-like object to GRanges
#'
#' @param df Data.frame that will be converted to GRanges
#' @param chrom Name of the column or expression that can be evaluated into chromosome name
#' @param start Name of the column or expression that can be evaluated into start coordinate
#' @param end Name of the column or expression that can be evaluated into end coordinate
#' @param strand Name of the column or expression that can be evaluated into strand
#' @param keep_coordinates_columns Specifies whether Columns used to specify coordinatess should be removed from GRanges mcols
#'
#' @return GRanges object
#' @examples
#' df = data.frame(query_chrom="chr1", query_start=1:2, query_end=1:2, col="AAA")
#' df2ranges(df, query_chrom, query_start, query_end)
df2ranges = function(df, chrom, start, end, strand=NULL, keep_coordinates_columns=T) {
  env = parent.frame()
  for(l in names(df)) env[[l]] = df[[l]]

  seqnames.field = eval(substitute(chrom), envir=env)
  start.field = eval(substitute(start), envir=env)
  end.field = eval(substitute(end), envir=env)
  coordinates.cols = c(substitute(chrom), substitute(start), substitute(end))
  has_strand = deparse(substitute(strand)) != "NULL"
  if(has_strand) {
    coordinates.cols = c(coordinates.cols, substitute(strand))
    strand.field = eval(substitute(strand), envir=env)
  } else {
    strand.field = rep("*", nrow(df))
  }

  ranges = GenomicRanges::GRanges(
    seqnames=seqnames.field,
    ranges=IRanges::IRanges(start=start.field, end=end.field),
    strand=strand.field)

  df.cols = colnames(df)
  if(!keep_coordinates_columns) {
    df = df[,setdiff(df.cols, coordinates.cols), drop=F]
  }

  GenomicRanges::mcols(ranges) = df

  ranges
}

#' @export
#' @description Lift over coordinates from one model to another
#'
#' @param ranges GRanges of original coordinates
#' @param chain_path Path to chain file used to migrate to new coordinates
#' @return Lifted GRanges object
liftOverRanges = function(ranges, chain_path) {
  chain = rtracklayer::import.chain(chain_path)
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
#' @description Read BED file
#' @param path Path to BED file
#' @return GRanges object
bed_read = function(path) {
  bed = rtracklayer::import.bed(path)
  GenomicRanges::start(bed) = GenomicRanges::start(bed)-1
  bed
}

#' @export
#' @description Get sequences from specified fasta file using coordinates
#' @param fasta Path to FASTA file
#' @param ranges GRanges object with coordinates
#' @return GRanges object with additional column `sequence` containing extracted sequences
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
#' @description Use BLAT to search for sequences in FASTA file
#'
#' @param sequences List of sequences to search for
#' @param fasta Path to FASTA file
#' @param minscore Smallest BLAT score to call a matching sequence (See BLAT manual)
#' @param stepSize Spacing between tiles (See BLAT manual)
#' @param tmp_dir Directory where the temporary BLAT-results file will be stored
#'
#' @return Data.frame with BLAT results
get_blat = function(sequences, fasta, minscore=30, stepSize=5, tmp_dir="tmp") {
  if(!bedr::check.binary("blat")) {
    stop("blat executible not found")
  }
  if(!bedr::check.binary("bowtie")) {
    stop("bowtie executible not found")
  }
  if(!bedr::check.binary("samtools")) {
    stop("samtools executible not found")
  }

  if(!file.exists(fasta)) { stop(paste0("Database file '", fasta, "' doesn't exist")) }

  dir.create(tmp_dir, recursive=T, showWarnings=F)

  fasta_path = file.path(tmp_dir, basename(tempfile()))
  psl_path = paste0(fasta_path, ".psl")
  sequences_df = data.frame(blat_query_name=paste0("seq", 1:length(sequences)), blat_sequence=sequences, stringsAsFactors=F)
  writeLines(paste0(">", sequences_df$blat_query_name, "\n", sequences_df$blat_sequence), con=fasta_path)

  cmd = paste0("blat -stepSize=", stepSize, " -repMatch=2253 -minScore=", minscore, " -minIdentity=0 -out=psl -dots=1 -noHead ", fasta, " ", fasta_path, " ", psl_path)
  writeLines(cmd)
  system(cmd)

  if(file.info(psl_path)$size > 0) {
    result = readr::read_tsv(psl_path, col_types=blat_columns(), col_names=names(blat_columns()$cols))
  } else {
    result = blat_blank()
  }

  file.remove(fasta_path)

  result %>%
    dplyr::mutate(blat_target_start=blat_target_start+1) %>%
    dplyr::inner_join(sequences_df, by="blat_query_name")
}

#' @export
#' @description Use BOWTIE to search for sequences in FASTA file
#'
#' @param sequences List of sequences to search for
#' @param fasta Path to FASTA file
#' @param threads Number of threads to use by BOWTIE aligner
#' @param tmp_dir Directory where the temporary BLAT-results file will be stored
#'
#' @return Data.frame with BOWTIE results
get_bowtie = function(sequences, fasta, threads=30, tmp_dir="tmp") {
  if(!bedr::check.binary("bowtie")) {
    stop("bowtie executible not found")
  }
  if(!bedr::check.binary("samtools")) {
    stop("samtools executible not found")
  }
  if(!file.exists(fasta)) { stop(paste0("Database file '", fasta, "' doesn't exist")) }

  fasta_name = gsub("\\.[^.]+", "", fasta)
  if(!file.exists(paste0(fasta_name, ".1.ebwt"))) {
    cmd_makedb = paste("bowtie-build --threads", threads, fasta, fasta_name)
    system(cmd_makedb)
  }

  dir.create(tmp_dir, recursive=T, showWarnings=F)
  fasta_path = file.path(tmp_dir, basename(tempfile()))
  sam_path = paste0(fasta_path, ".sam")
  bam_path = paste0(fasta_path, ".bam")
  sequences_df = data.frame(bowtie_qid=paste0("seq", 1:length(sequences)), bowtie_sequence=sequences, stringsAsFactors=F)
  sequences_reduced_df = sequences_df %>% dplyr::distinct(bowtie_sequence, .keep_all=T)
  writeLines(paste0(">", sequences_reduced_df$bowtie_qid, "\n", sequences_reduced_df$bowtie_sequence), con=input_path)

  cmd = paste0("bowtie --threads ", threads, " --tryhard --all -v 3 --seedlen 5 -f --sam --no-unal ", fasta_name, " ", input_path, " > ", sam_path)
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

  file.remove(bam_path)
  file.remove(input_path)

  primers_alignments
}

#' @export
#' @description Use BOWTIE to search for sequences in FASTA file
#'
#' @param sequences List of sequences to search for
#' @param fasta Path to FASTA file
#' @param word_size Word size in BLAST alignment (see BLASt documentation)
#' @param perc_identity Percent of identity to call matching alignment (see BLASt documentation)
#' @param tmp_dir Directory where the temporary BLAT-results file will be stored
#'
#' @return Data.frame with BLAST results
get_blast = function(sequences, fasta, word_size=4, perc_identity=100, tmp_dir="tmp") {
  if(!bedr::check.binary("blastn")) {
    stop("blast executible not found")
  }
  if(!file.exists(fasta)) { stop(paste0("Database file '", fasta, "' doesn't exist")) }

  fasta_name = gsub("\\.[^.]+", "", fasta)
  if(!file.exists(paste0(fasta_name, ".nsq"))) {
    cmd_makedb = paste0(paste0("makeblastdb -in ", fasta, " -out ", fasta_name, " -dbtype 'nucl' -hash_index"))
    system(cmd_makedb)
  }

  dir.create(tmp_dir, recursive=T, showWarnings=F)
  input_path = file.path(tmp_dir, basename(tempfile()))
  sequences_df = data.frame(blast_qseqid=paste0("seq", 1:length(sequences)), blast_sequence=sequences, stringsAsFactors=F)
  writeLines(paste0(">", sequences_df$blast_query_name, "\n", sequences_df$blast_sequence), con=input_path)

  out_path = paste0(input_path, ".blast")
  cmd = paste0("blastn -db ", fasta_name, " -query ", input_path, " -out ", out_path, " -word_size ", word_size, " -perc_identity ", perc_identity, " -outfmt 6")
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
#' @description Efficiently count number of lines in specified text file
#' @param path Path to file
#' @return Number of lines in text file
file_count_lines = function(path) {
  f = file(path, open="rb")
  nlines = 0L
  while (length(chunk <- readBin(f, "raw", 65536)) > 0) {
     nlines = nlines + sum(chunk == as.raw(10L))
  }
  close(f)
  nlines
}

#' @export
#' @description Trim all values to be in lower- and upper- bound intervals. Values that are outside that interval
#' are set to `lb` and `ub` respectfully
#'
#' @param value Numeric vector
#' @param lb Lower bound
#' @param ub Upper bound
#'
trim = function(value, lb, ub) {
  pmin(pmax(value, lb), ub)
}