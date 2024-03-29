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


#' @title leftJoinByOverlaps
#' @export
#' @description Left-join query and subject ranges based on overlap between them
#'
#' @param query_ranges Query ranges
#' @param subject_ranges Subject ranges
#' @param ... Further parameters are passed to IRanges::mergeByOverlaps
#'
#' @return A data frame with two ranges objects left-joined based on the overlap between them
#' @examples
#' query_ranges = data.frame(query_chrom="chr1", query_start=1:2, query_end=1:2, col="AAA") %>% df2ranges(query_chrom, query_start, query_end)
#' subject_ranges = data.frame(subject_chrom="chr1", subject_start=2:3, subject_end=2:3, col="AAA") %>% df2ranges(subject_chrom, subject_start, subject_end)
#' leftJoinByOverlaps(query_ranges, subject_ranges)
leftJoinByOverlaps = function(query_ranges, subject_ranges, ...) {
  r = separate_target_subject_columns(query_ranges, subject_ranges)
  query_ranges = r$query_ranges
  subject_ranges = r$subject_ranges

  if(length(query_ranges)>0) {
    query_ranges$query_id = 1:length(query_ranges)
  }
  if(length(subject_ranges)>0) {
    subject_ranges$subject_id = 1:length(subject_ranges)
  }
  result_df = as.data.frame(IRanges::mergeByOverlaps(query_ranges, subject_ranges, ...)) %>% dplyr::select(-dplyr::matches("^(query|subject)_ranges\\."))
  result_df = dplyr::bind_rows(result_df, as.data.frame(GenomicRanges::mcols(query_ranges)) %>% dplyr::anti_join(result_df %>% dplyr::select(query_id), by="query_id"))
  result_df = result_df %>% dplyr::select(-dplyr::matches("^(query|subject)_ranges\\.")) %>% dplyr::select(-dplyr::matches("^query_id|subject_id$"))

  result_df
}

#' @title fullJoinByOverlaps
#' @export
#' @description Full-join query and subject ranges based on overlap between them
#'
#' @param query_ranges Query ranges
#' @param subject_ranges Subject ranges
#' @param ... Further parameters are passed to IRanges::mergeByOverlaps
#'
#' @return A data frame with two ranges objects full-joined based on the overlap between them
#' @examples
#' query_ranges = data.frame(query_chrom="chr1", query_start=1:2, query_end=1:2, col="AAA") %>% df2ranges(query_chrom, query_start, query_end)
#' subject_ranges = data.frame(subject_chrom="chr1", subject_start=2:3, subject_end=2:3, col="AAA") %>% df2ranges(subject_chrom, subject_start, subject_end)
#' fullJoinByOverlaps(query_ranges, subject_ranges)
fullJoinByOverlaps = function(query_ranges, subject_ranges, ...) {
  r = separate_target_subject_columns(query_ranges, subject_ranges)
  query_ranges = r$query_ranges
  subject_ranges = r$subject_ranges

  if(length(query_ranges)>0) {
    query_ranges$query_id = 1:length(query_ranges)
  }
  if(length(subject_ranges)>0) {
    subject_ranges$subject_id = 1:length(subject_ranges)
  }

  result_df = as.data.frame(IRanges::mergeByOverlaps(query_ranges, subject_ranges, ...))
  result_df = dplyr::bind_rows(
    result_df,
    as.data.frame(GenomicRanges::mcols(query_ranges)) %>% dplyr::anti_join(result_df %>% dplyr::select(query_id), by="query_id"),
    as.data.frame(GenomicRanges::mcols(subject_ranges)) %>% dplyr::anti_join(result_df %>% dplyr::select(subject_id), by="subject_id"))
  result_df = result_df %>% dplyr::select(-dplyr::matches("^(query|subject)_ranges\\."), -dplyr::matches("^(subject_id|query_id)$"))

  result_df
}

#' @title innerJoinByOverlaps
#' @export
#' @description Inner-join query and subject ranges based on overlap between them
#'
#' @param query_ranges Query ranges
#' @param subject_ranges Subject ranges
#' @param ... Further parameters are passed to IRanges::mergeByOverlaps
#'
#' @return A data frame with two ranges objects inner-joined based on the overlap between them
#' @examples
#' query_ranges = data.frame(query_chrom="chr1", query_start=1:2, query_end=1:2, col="AAA") %>% df2ranges(query_chrom, query_start, query_end)
#' subject_ranges = data.frame(subject_chrom="chr1", subject_start=2:3, subject_end=2:3, col="AAA") %>% df2ranges(subject_chrom, subject_start, subject_end)
#' innerJoinByOverlaps(query_ranges, subject_ranges)
innerJoinByOverlaps = function(query_ranges, subject_ranges, ...) {
  r = separate_target_subject_columns(query_ranges, subject_ranges)
  query_ranges = r$query_ranges
  subject_ranges = r$subject_ranges

  result_ranges = IRanges::mergeByOverlaps(query_ranges, subject_ranges, ...)
  as.data.frame(result_ranges) %>% dplyr::select(-dplyr::matches("_ranges\\."))
}

#' @export
innerJoinManyByOverlaps = function(ranges_list) {
  results_ranges = ranges_list[[1]]
  for(r_ranges in ranges_list[2:length(ranges_list)]) {
    results_ranges = IRanges::mergeByOverlaps(results_ranges, r_ranges)
    results_df = as.data.frame(results_ranges)
    results_ranges = df2ranges(results_df, results_ranges.seqnames, results_ranges.start, results_ranges.end, results_ranges.strand)
    results_cols = colnames(results_df)[!grepl("_ranges\\.", colnames(results_df))]
    results_ranges = results_ranges[,results_cols]
  }

  as.data.frame(results_ranges)
}

#' @title df2ranges
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

#' @title liftOverRanges
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


#' @title bed_read
#' @export
#' @description Read BED file
#' @param path Path to BED file
#' @return GRanges object
bed_read = function(path) {
  bed = rtracklayer::import.bed(path)
  GenomicRanges::start(bed) = GenomicRanges::start(bed)-1
  bed
}

#' @title get_seq
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

#' @title get_pairwise_alignment
#' @export
#' @description Find pairwise identity between two sets of sequences
#'
#' @param seq1 First set of sequences (pattern)
#' @param seq2 Second set of sequences (subject)
#'
#' @return Vector with identity scores between 1 and 100
get_pairwise_alignment = function(seq1, seq2, gapOpening=10, gapExtension=4) {
  y1 = Biostrings::pairwiseAlignment(seq1, seq2, type="global-local", gapOpening=gapOpening, gapExtension=gapExtension)
  y2 = Biostrings::pairwiseAlignment(Biostrings::reverseComplement(Biostrings::DNAStringSet(seq1)), seq2, type="global-local", gapOpening=gapOpening, gapExtension=gapExtension)

  y1_df = cbind(data.frame(score=Biostrings::score(y1), pid=Biostrings::pid(y2)), as.data.frame(IRanges::ranges(Biostrings::subject(y1))))
  y2_df = cbind(data.frame(score=Biostrings::score(y2), pid=Biostrings::pid(y2)), as.data.frame(IRanges::ranges(Biostrings::subject(y2))))
  ret_df = y1_df
  ret_df[y2_df$score>y1_df$score,] = y2_df[y2_df$score>y1_df$score,]

  ret_df %>% dplyr::select(-width)
}

#' @title get_blat
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

#' @title get_bowtie
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

#' @title clear_tmpdir
#' @export
#' @description Clear cache directory from old files
#'
#' @param tmp_dir Path to cache directory
#' @param max_age Max age of files in seconds. Remove everything older than specified age
clear_tmpdir = function(tmp_dir="tmp", max_age=2419200) {
  if(is.null(max_age)) {
    max_age = 604800
  }
  for(f in list.files(tmp_dir, recursive=T, full.names=T)) {
    ctime = file.info(f)$ctime
    cdiff = as.numeric(difftime(Sys.time(), ctime, units = "secs"))
    if(cdiff>max_age) {
      file.remove(f)
    }
  }
}

#' @title get_blast
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

#' @title file_count_lines
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

#' @title trim
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

distr_call = function(ver, distr, x, params, lower.tail=T) {
  # z.ver <<- ver
  # z.distr <<- distr
  # z.x <<- x
  # z.params <<- params
  # z.lower.tail <<- lower.tail
  # ver = z.ver
  # distr=z.distr
  # x=z.x
  # params=z.params
  # lower.tail=z.lower.tail


  params_clean = params %>%
    dplyr::select(dplyr::starts_with("bgmodel_"), -dplyr::matches("bgmodel_distr|bgmodel_chrom|bgmodel_strand")) %>%
    dplyr::slice(1) %>%
    setNames(gsub("bgmodel_", "", colnames(.))) %>%
    as.list()
  if(ver=="d") params_clean[["x"]] = x
  if(ver=="r") params_clean[["n"]] = x
  if(ver=="p") params_clean[["q"]] = x
  if(ver=="q") params_clean[["p"]] = x
  if(ver=="p" || ver=="q") params_clean[["lower.tail"]] = lower.tail

  do.call(paste0(ver, distr), params_clean)
}

#' @title coverage_find_empty_intervals
#' @export
#'
#' @description Finds ranges in pileup data where there is not enough coverage
#'
#' @param coverage_ranges Pileup data that is used to find the intervals
#' @param coverage_column Column name with pileup numbers
#' @param minlen Minimal length of the range with data sattisfying maxcoverage critaria
#' @param maxcoverage A threshold on amount of data that would fullfills coverage criteria. Areas with less than this threshold are considered empty
#'
#' @return GRanges object with ranges sattisfying coverage criteria
coverage_find_empty_intervals = function(coverage_ranges, coverage_column="score", minlen=1e6, maxcoverage=0)
{
  if(!(coverage_column %in% colnames(GenomicRanges::mcols(coverage_ranges)))) {
    stop(paste0("Column '", coverage_column, "' not found in coverage ranges data columns"))
  }
  mask_ranges = GenomicRanges::reduce(coverage_ranges[GenomicRanges::mcols(coverage_ranges)[[coverage_column]]<=maxcoverage]) %>%
    as.data.frame() %>%
    dplyr::group_by(seqnames) %>%
    dplyr::filter(width>=minlen) %>%
    dplyr::ungroup() %>%
    dplyr::select(mask_chrom=seqnames, mask_start=start, mask_end=end) %>%
    df2ranges(mask_chrom, mask_start, mask_end)

  mask_ranges
}

coverage_merge_strands = function(coverage_ranges, aggregate_fun, score_column="score")
{
  coverage_df = as.data.frame(coverage_ranges) %>% dplyr::select_at(c(raw_seqnames="seqnames", raw_start="start", raw_end="end", raw_strand="strand", raw_score=score_column))
  coverage_clean_ranges = coverage_df %>% df2ranges(raw_seqnames, raw_start, raw_end, raw_strand)
  merged_long_df = coverage_df %>%
    reshape2::melt(measure.vars=c("raw_start", "raw_end"), value.name="sum_pos") %>%
    dplyr::distinct_at(c(sum_seqnames="raw_seqnames", "sum_pos")) %>%
    dplyr::arrange(sum_seqnames, sum_pos) %>%
    dplyr::mutate(sum_start=dplyr::lag(sum_pos), sum_end=sum_pos-1, sum_start=ifelse(sum_start>sum_end, 1, sum_start)) %>%
    dplyr::filter(!is.na(sum_start)) %>%
    dplyr::select(sum_seqnames, sum_start, sum_end) %>%
    df2ranges(sum_seqnames, sum_start, sum_end) %>%
    innerJoinByOverlaps(coverage_clean_ranges)

  if(any(duplicated(merged_long_df[,c("sum_seqnames", "sum_start", "sum_end", "raw_strand")]))) {
    stop(paste0("Found overlaping regions in single strand data.\nCheck that strand information is not missing from provided ranges object.\nStrands present in ranges object: ", paste(as.character(unique(coverage_df$raw_strand)), collapse=",")))
  }

  res_strands = as.character(unique(merged_long_df$raw_strand))
  merged_wide_df = merged_long_df %>%
    reshape2::dcast(sum_seqnames+sum_start+sum_end ~ raw_strand, value.var="raw_score") %>%
    replace(is.na(.), 0)
  merged_wide_df$score = apply(merged_wide_df[, res_strands, drop=F], 1, FUN=aggregate_fun)
  merged_wide_ranges = GenomicRanges::makeGRangesFromDataFrame(merged_wide_df, ignore.strand=T, keep.extra.columns=T)

  as(GenomicRanges::coverage(merged_wide_ranges, weight=GenomicRanges::mcols(merged_wide_ranges)$score), "GRanges")
}

ranges_sample = function(ranges, mask_ranges, column, ntile=100000) {
  ranges_df = as.data.frame(ranges)

  #
  # Prepare a strand-specific version of sampled ranges and calculate coverage coverage
  #
  ranges_strand_df = ranges_df %>%
    dplyr::mutate(group_name=paste0(seqnames, ";", strand)) %>%
    dplyr::select_at(c(group_name="group_name", group_start="start", group_end="end", column))
  ranges_strand = ranges_strand_df %>%
    df2ranges(group_name, group_start, group_end)
  df_score = GenomicRanges::coverage(ranges_strand, weight=GenomicRanges::mcols(ranges_strand)[[column]])

  #
  # Prepare a strand-specific version of mask that shall be removed from sampling
  #
  mask_strand_ranges = as.data.frame(mask_ranges)%>%
    dplyr::select(mask_chrom, mask_start, mask_end) %>%
    tidyr::crossing(data.frame(strand=unique(ranges_df$strand)))%>%
    dplyr::mutate(mask_chrom=paste0(mask_chrom, ";", strand))%>%
    df2ranges(mask_chrom, mask_start, mask_end)

  #
  # Prepare a strand-specific genome ranges object to be sampled
  #
  seqlengts_strand_df = ranges_strand_df %>%
    dplyr::group_by(group_name) %>%
    dplyr::summarize(seqlengths=max(group_end))
  genome_strand_ranges = seqlengts_strand_df %>%
    df2ranges(group_name, 1, seqlengths)

  sampled_ranges = GenomicRanges::setdiff(genome_strand_ranges, mask_strand_ranges)
  tile_ranges = unlist(GenomicRanges::tile(sampled_ranges, width=round(sum(GenomicRanges::width(sampled_ranges))/ntile)))
  GenomicRanges::end(tile_ranges) = GenomicRanges::start(tile_ranges)
  GenomeInfoDb::seqlevels(tile_ranges) = names(df_score)
  # tile_ranges = unlist(GenomicRanges::tileGenome(seqlengts, ntile=ntile, cut.last.tile.in.chrom=F))
  # GenomicRanges::end(tile_ranges) = GenomicRanges::start(tile_ranges)

  ret = as.data.frame(GenomicRanges::binnedAverage(tile_ranges, df_score, column)) %>%
    dplyr::mutate(strand=gsub(".*;", "", seqnames), seqnames=gsub(";.*", "", seqnames)) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=T)

  ret
}

seqlengths2tiles = function(seqlengths, width, step) {
  genome_tiles_template = GenomicRanges::tileGenome(seqlengths, tilewidth=width, cut.last.tile.in.chrom=T)
  genome_tiles = GenomicRanges::GRanges()
  for(str in seq(0, width-1, by=step)) {
    genome_tiles = suppressWarnings(IRanges::append(genome_tiles, GenomicRanges::trim(IRanges::shift(genome_tiles_template, str))))
  }
  genome_tiles = as.data.frame(GenomicRanges::sort(genome_tiles)) %>%
    dplyr::rename(tile_chrom="seqnames", tile_start="start", tile_end="end") %>%
    dplyr::select(tile_chrom, tile_start, tile_end) %>%
    df2ranges(tile_chrom, tile_start, tile_end)
}

ranges2pure_tiles = function(ranges, column, width=1000, step=width, diff=0.1) {
  df_score = GenomicRanges::coverage(ranges, weight=GenomicRanges::mcols(ranges)[[column]])

  seqlengts = as.data.frame(ranges) %>%
    dplyr::group_by(seqnames) %>%
    dplyr::summarize(seqlengths=max(end)) %>%
    tibble::deframe()
  df_bins_template = GenomicRanges::tileGenome(seqlengts, tilewidth=width, cut.last.tile.in.chrom=T)
  df_bins = df_bins_template
  for(str in seq(0, width-1, by=step)[-1]) {
    df_bins = IRanges::append(df_bins, GenomicRanges::trim(IRanges::shift(df_bins_template, str)))
  }
  df_bins = GenomicRanges::sort(df_bins)
  # GenomicRanges::binnedAverage(df_bins, df_score, column)

  bins_per_chrom <- split(GenomicRanges::ranges(df_bins), GenomicRanges::seqnames(df_bins))
  means_list <- lapply(names(df_score), function(seqname) {
      ss <<- seqname
      v = IRanges::Views(df_score[[seqname]], bins_per_chrom[[seqname]])
      x = IRanges::viewMaxs(v) - IRanges::viewMeans(v) <= diff & is.finite(IRanges::viewMaxs(v))
      m = data.frame(col=IRanges::viewMaxs(v))
      colnames(m) = column

      GenomicRanges::mcols(bins_per_chrom[[seqname]]) = m
      GenomicRanges::GRanges(seqnames=seqname, bins_per_chrom[[seqname]][x])
  })

  ret = suppressWarnings(do.call(what = c, args = means_list))
  ret
}


ranges2seqlengths = function(ranges) {
  as.data.frame(ranges) %>%
    dplyr::group_by(seqnames) %>%
    dplyr::summarize(seqlengths=max(end)) %>%
    tibble::deframe()
}

bootstrap_data_overlaps = function(evaluated_ranges, data_ranges, excluded_ranges=evaluated_ranges, genome_tiles_step=10000, genome_tiles_width=50000, n_samples=1000) {
  seqlengths = ranges2seqlengths(data_ranges)
  genome_tiles = seqlengths2tiles(seqlengths, genome_tiles_width, genome_tiles_step)
  genome_tiles = genome_tiles[GenomicRanges::width(genome_tiles)==genome_tiles_width]

  evaluated_df = as.data.frame(evaluated_ranges) %>%
    dplyr::select(evaluated_chrom=seqnames, evaluated_start=start, evaluated_end=end) %>%
    dplyr::mutate(evaluated_tile_count=ceiling((evaluated_end-evaluated_start) / genome_tiles_width))
  evaluated_ranges = evaluated_df %>% df2ranges(evaluated_chrom, evaluated_start, evaluated_end)
  data_df = as.data.frame(data_ranges) %>%
    dplyr::select(data_chrom=seqnames, data_start=start, data_end=end)
  data_ranges = data_df %>% df2ranges(data_chrom, data_start, data_end)
  excluded_df = as.data.frame(excluded_ranges) %>%
    dplyr::select(excluded_chrom=seqnames, excluded_start=start, excluded_end=end)
  excluded_ranges = excluded_df %>% df2ranges(excluded_chrom, excluded_start, excluded_end)

  bg_tiles_df = genome_tiles %>%
    leftJoinByOverlaps(excluded_ranges) %>%
    dplyr::filter(is.na(excluded_start)) %>%
    dplyr::select(dplyr::starts_with("tile_")) %>%
    dplyr::right_join(evaluated_df, by=c("tile_chrom"="evaluated_chrom")) %>%
    dplyr::group_by(bootstrap_chrom=tile_chrom, bootstrap_start=evaluated_start, bootstrap_end=evaluated_end) %>%
    dplyr::sample_n(evaluated_tile_count[1]*n_samples, replace=T) %>%
    dplyr::group_by(bootstrap_chrom, bootstrap_start, bootstrap_end) %>%
    dplyr::mutate(bootstrap_sample_num=floor((1:dplyr::n()-1)/evaluated_tile_count)+1) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(has_data=!is.na(tile_chrom) & !is.na(tile_start) & !is.na(tile_end))
  bg_df = bg_tiles_df %>%
    dplyr::mutate(bootstrap_data_count=0) %>%
    dplyr::group_by(has_data) %>%
    dplyr::do((function(x){
      xx <<-x
      if(!x$has_data[1]) return(x)
      x$bootstrap_data_count = x %>%
        df2ranges(tile_chrom, tile_start, tile_end) %>%
        GenomicRanges::countOverlaps(data_ranges)
      x
    })(.)) %>%
    dplyr::group_by(bootstrap_chrom, bootstrap_start, bootstrap_end, bootstrap_sample_num) %>%
    dplyr::summarize(bootstrap_data_count=sum(bootstrap_data_count), bootstrap_type="background") %>%
    dplyr::group_by(bootstrap_chrom, bootstrap_start, bootstrap_end) %>%
    dplyr::mutate(bootstrap_data_mean=mean(bootstrap_data_count)) %>%
    dplyr::ungroup()
  sg_df = tibble::tibble(
    bootstrap_chrom=evaluated_ranges$evaluated_chrom,
    bootstrap_start=evaluated_ranges$evaluated_start,
    bootstrap_end=evaluated_ranges$evaluated_end,
    bootstrap_sample_num=1,
    bootstrap_type="signal",
    bootstrap_data_count=evaluated_ranges %>% GenomicRanges::countOverlaps(data_ranges)) %>%
    dplyr::mutate(bootstrap_data_mean=bootstrap_data_count)
  dplyr::bind_rows(bg_df, sg_df)
}