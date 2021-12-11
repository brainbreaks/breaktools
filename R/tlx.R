validate_group = function(group) {
  valid_group = c("none", "group", "sample")
  if(!(group %in% valid_group)) { stop(simpleError(paste0("group should be one of: ", paste(valid_group, collapse=", "), " (Actual:", group, ")"))) }
}

validate_group_within = function(within) {
  within_group = c("all", "group", "none")
  if(!(within %in% within_group)) { stop(simpleError(paste0("group should be one of: ", paste(within_group, collapse=", "), " (Actual:", within, ")"))) }
}


validate_normalization_target = function(target) {
  valid_targets = c("smallest", "largest", "mean", "median")
  if(!(target %in% valid_targets)) { stop(simpleError(paste0("normalization target should be one of: ", paste(valid_targets, collapse=", "), " (Actual:", target, ")"))) }
}


validate_exttype = function(exttype) {
  valid_exttype = c("symmetrical", "along", "none")
  if(!(exttype %in% valid_exttype)) { stop(simpleError(paste0("exttype should be one of: ", paste(valid_exttype, collapse=", "), " (Actual:", exttype, ")"))) }
}

#' @export
tlx_cols = function() {
  readr::cols(
    Qname=readr::col_character(), JuncID=readr::col_character(), Rname=readr::col_character(), Junction=readr::col_double(),
    Strand=readr::col_character(), Rstart=readr::col_double(), Rend=readr::col_double(),
    B_Rname=readr::col_character(), B_Rstart=readr::col_double(), B_Rend=readr::col_double(), B_Strand=readr::col_double(),
    B_Qstart=readr::col_double(), B_Qend=readr::col_double(), Qstart=readr::col_double(), Qend=readr::col_double(), Qlen=readr::col_double(),
    B_Cigar=readr::col_character(), Cigar=readr::col_character(), Seq=readr::col_character(), J_Seq=readr::col_character(), Barcode=readr::col_logical(),
    unaligned=readr::col_double(), baitonly=readr::col_double(), uncut=readr::col_double(), misprimed=readr::col_double(), freqcut=readr::col_double(),
    largegap=readr::col_double(), mapqual=readr::col_double(), breaksite=readr::col_double(), sequential=readr::col_double(), repeatseq=readr::col_double(), duplicate=readr::col_double()
  )
}

#' @export
tlx_blank = function() {
  blank_tibble(tlx_cols()) %>%
    dplyr::mutate(tlx_sample=NA_character_, tlx_path=NA_character_, tlx_group=NA_character_, tlx_control=NA) %>%
    dplyr::mutate(tlx_is_bait_chromosome=NA, tlx_is_bait_junction=NA, tlx_is_offtarget=NA) %>%
    dplyr::slice(0)
}

#' @export
tlx_read = function(path, sample, group="", group_i=1, control=F) {
  readr::read_tsv(path, comment="#", skip=16, col_names=names(tlx_cols()$cols), col_types=tlx_cols()) %>%
    dplyr::mutate(tlx_strand=ifelse(Strand<0, "-", "+")) %>%
    dplyr::mutate(Seq_length=nchar(Seq), tlx_sample=sample, tlx_path=path, tlx_group=group, tlx_group_i=group_i, tlx_control=control)
}

#' @export
tlx_read_many = function(samples_df) {
  tlx_df.all = data.frame()
  for(f in 1:nrow(samples_df)) {
    log("Reading tlx file ", f, "/", nrow(samples_df), ": ",  samples_df$path[f])
    tlx_df.f = tlx_read(samples_df$path[f], sample=samples_df$sample[f], control=samples_df$control[f], group=samples_df$group[f], group_i=samples_df$group_i[f])
    tlx_df.all = dplyr::bind_rows(tlx_df.all, tlx_df.f)
  }

  tlx_df.all
}

tlx_generate_filename_col = function(df, include_sample=F, include_group=F, include_strand=F, include_treatment=T) {
  if(!include_sample & !include_group & !include_strand & !include_treatment) {
    stop(simpleError("At least one parameter must be used to generate file name (include_sample, include_group, include_strand, include_treatment"))
  }
  if(include_sample & !("tlx_sample" %in% colnames(df))) stop(simpleError("tlx_sample column is not pressent in data.frame"))
  if(include_group & !("tlx_group" %in% colnames(df))) stop(simpleError("tlx_group column is not pressent in data.frame"))
  if(include_strand & !("tlx_strand" %in% colnames(df))) stop(simpleError("tlx_strand column is not pressent in data.frame"))
  if(include_treatment & !("tlx_control" %in% colnames(df))) stop(simpleError("tlx_control column is not pressent in data.frame"))

  map_strand = c("+"="plus", "-"="minus")
  filenames_df = df %>% dplyr::select()
  if(include_group) filenames_df$group = df$tlx_group
  if(include_sample) filenames_df$sample = df$tlx_sample
  if(include_treatment) filenames_df$control = ifelse(df$tlx_control, "ctrl", "trmnt")
  if(include_strand) filenames_df$strand = map_strand[df$tlx_strand]
  table(filenames_df$group)
  table(df$tlx_group)
  filenames_df %>% dplyr::distinct()
  # generated_path = apply(filenames_df, 1, paste, collapse="_")
  # generated_path = unlist(purrr::pmap(filenames_df, paste, sep = '_'))
  generated_path = do.call(paste, c(filenames_df, sep = "_"))
  unique_generated_path = unique(generated_path)
  map_unique_generated_path = gsub("[/.]", "-", unique_generated_path)
  map_unique_generated_path = gsub("[^A-Za-z0-9 -]+", "_", map_unique_generated_path)
  map_unique_generated_path = gsub("(_| )+", "\\1", map_unique_generated_path)
  names(map_unique_generated_path) = unique_generated_path
  unname(map_unique_generated_path[generated_path])
}

#' @export
tlx_write_bedgraph = function(tlx_df, path, group_within, exttype, extsize, normalize_within, normalization_target="smallest", split_strand=F) {
  validate_group_within(group_within)
  validate_group_within(normalize_within)
  validate_normalization_target(normalization_target)
  validate_exttype(exttype)

  writeLines("Calculating coverage...")
  tlxcov_df = tlx_coverage(tlx_df, group_within=group_within, exttype=exttype, extsize=extsize, normalize_within=normalize_within, normalization_target=normalization_target, ignore.strand=!split_strand) %>%
    dplyr::arrange(tlxcov_chrom, tlxcov_start)

  writeLines("calculating filenames(s)...")
  if(group_within=="all") tlxcov_df = tlxcov_df %>% dplyr::mutate(g=tlx_generate_filename_col(., include_group=F, include_sample=F, include_treatment=T, include_strand=split_strand))
  if(group_within=="all") tlxcov_df = tlxcov_df %>% dplyr::mutate(g=tlx_generate_filename_col(., include_group=F, include_sample=F, include_treatment=T, include_strand=split_strand))
  if(group_within=="group") tlxcov_df = tlxcov_df %>% dplyr::mutate(g=tlx_generate_filename_col(., include_group=T, include_sample=F, include_treatment=T, include_strand=split_strand))
  if(group_within=="sample") tlxcov_df = tlxcov_df %>% dplyr::mutate(g=tlx_generate_filename_col(., include_group=F, include_sample=T, include_treatment=T, include_strand=split_strand))
  if(split_strand) tlxcov_df = tlxcov_df %>% dplyr::mutate(tlxcov_pileup=ifelse(tlx_strand=="+", 1, -1)*tlxcov_pileup)


  writeLines("Writing bedgraph file(s)...")
  if(!dir.exists(path)) dir.create(path, recursive=T)
  tlxcov_df %>%
    dplyr::group_by(g) %>%
    dplyr::do((function(z){
      z.out = z %>% dplyr::select(tlxcov_chrom, tlxcov_start, tlxcov_end, tlxcov_pileup)
      z.path = file.path(path, paste0(z$g[1], ".bedgraph"))
      writeLines(paste0("Writing to file '", z.path, "'"))
      readr::write_tsv(z.out, file=z.path, col_names=F)
      data.frame()
    })(.))
  return(0)
}

#' @export
tlx_write_bed = function(tlx_df, path, group_within, split_strand=F) {
  validate_group_within(group_within)

  tlx_bed_df = tlx_df %>%
    dplyr::mutate(start=Junction, end=Junction, name=paste0(Qname, " (", tlx_sample, ")"))
    # dplyr::mutate(start=ifelse(Strand=="-1", Junction-1, Junction), end=ifelse(Strand=="-1", Junction, Junction+1))

  writeLines("calculating filenames(s)...")
  if(group_within=="all") tlx_bed_df = tlx_bed_df %>% dplyr::mutate(g=tlx_generate_filename_col(., include_group=F, include_sample=F, include_treatment=T, include_strand=split_strand))
  if(group_within=="group") tlx_bed_df = tlx_bed_df %>% dplyr::mutate(g=tlx_generate_filename_col(., include_group=T, include_sample=F, include_treatment=T, include_strand=split_strand))
  if(group_within=="sample") tlx_bed_df = tlx_bed_df %>% dplyr::mutate(g=tlx_generate_filename_col(., include_group=F, include_sample=T, include_treatment=T, include_strand=split_strand))

  writeLines("Writing bedgraph file(s)...")
  if(!dir.exists(path)) dir.create(path, recursive=T)
  tlx_bed_df %>%
    dplyr::group_by(g) %>%
    dplyr::do((function(z){
      z.out = z %>% dplyr::select(Rname, start, end, name, mapqual, tlx_strand)
      z.path = file.path(path, paste0(z$g[1], ".bed"))
      writeLines(paste0("Writing to file '", z.path, "'"))
      readr::write_tsv(z.out, file=z.path, col_names=F)
      data.frame()
    })(.))
  return(0)
}

#' @title tlx_coverage
#'
#' @description Calculates coverage from tlx dataframe. The coverage can be pulled from multiple samples based on \code{group} parameter and
#' normalized base on \code{normalize_within} and \code{normalize_between} and \code{normalization_target}
#' on group_within
#'
#' @param tlx_df TLX dataframe produced by HTGTS translocation pipeline
#' @param group Grouping columns. Possible values are:
#'    all         - group together all treatment and control samples
#'    group       - group together all treatment and control samples from each group separately
#'    none - Do not group anything
#' @param extsize Each junction coordinate is extended by \code{extsize} bais pairs to calculate pileup
#' @param exttype Type of junction extention for pileup calculation:
#'    symmetrical - Extention is done around the center of the junction extending by \code{extsize/2} each direction
#'    along       - Extention is done in the direction corresponding to pray strand and extending \code{extsize} that direction
#' @param normalize_within Reads normalization strategy within the group. It signifies how each sample will be normalized
#' before pulling it together with other samples from this group. Possible values are:
#'    all         - Group together all treatment and control samples
#'    group       - Group together each group treatment and control samples separately
#'    none        - Do not group anything
#' @param normalize_between Reads normalization strategy between the groups. It signifies how each group is normalized
#' to allow direct comparison to other groups. Possible values are:
#'    all         - All groups control and treetment libraries are normalized to have same size so that they can be directly compared
#'    group       - Control and treatment libraries are normalized within each group separately so that they can be compared within the group
#'    none        - No normalization is done
#' @param normalization_target Normalization strategy.
#'    smallest    - Normalize to smallest sample/group (make other samples smaller)
#'    largest     - Normalize to largest sample/group (make other samples bigger)
#'    mean        - Normalize to mean
#'    meadian     - Normalize to median
#'
#' @return A data frame with coverages for each sample or group
#' @examples
#' no examples yet
tlx_coverage = function(tlx_df, group, extsize, exttype, normalize_within=NULL, normalize_between=NULL, normalization_target="smallest", ignore.strand=T) {
  if(is.null(normalize_within)) normalize_within = group
  if(is.null(normalize_between)) normalize_between = group
  validate_group_within(group)
  validate_group_within(normalize_within)
  validate_group_within(normalize_between)
  validate_exttype(exttype)
  validate_normalization_target(normalization_target)

  tlx_coverage_ = function(x, extsize, exttype) {
    validate_exttype(exttype)

    if(exttype[1]=="along") {
      x_ranges  = GenomicRanges::makeGRangesFromDataFrame(x %>% dplyr::mutate(seqnames=Rname, sstart=ifelse(Strand=="-1", Junction-extsize, Junction-1), end=ifelse(Strand=="-1", Junction, Junction+extsize-1)), ignore.strand=T, keep.extra.columns=T)
    } else {
      if(exttype[1]=="symmetrical") {
        x_ranges  = GenomicRanges::makeGRangesFromDataFrame(x %>% dplyr::mutate(seqnames=Rname, start=Junction-ceiling(extsize/2), end=Junction+ceiling(extsize/2)), ignore.strand=T, keep.extra.columns=T)
      } else {
        x_ranges  = GenomicRanges::makeGRangesFromDataFrame(x %>% dplyr::mutate(seqnames=Rname, start=Junction, end=Junction+1), ignore.strand=T, keep.extra.columns=T)
      }
    }

    cov_ranges = as(GenomicRanges::coverage(x_ranges), "GRanges")
    ret_df = as.data.frame(cov_ranges) %>%
      dplyr::rename(tlxcov_chrom="seqnames", tlxcov_start="start", tlxcov_end="end", tlxcov_pileup="score") %>%
      dplyr::select(matches("tlxcov_"))
    ret_df
  }

  if(group=="all") group_cols = c("tlx_control")
  if(group=="group") group_cols = c("tlx_group", "tlx_control")
  if(group %in% c("sample", "none")) group_cols = c("tlx_group", "tlx_group_i", "tlx_control", "tlx_sample")
  if(!ignore.strand) group_cols = c(group_cols, "tlx_strand")
  if(normalize_within=="all") normalize_within_cols = c("tlx_control")
  if(normalize_within=="group") normalize_within_cols = c("tlx_group", "tlx_control")
  if(normalize_within %in% "none") normalize_within_cols = c("tlx_sample", "tlx_control")
  if(normalize_between=="all") normalize_between_cols = c()
  if(normalize_between=="group") normalize_between_cols = setdiff(group_cols, "tlx_control")
  if(normalize_between %in% "none") normalize_between_cols = group_cols

  normalization_target_fun = c("smallest"=min, "largest"=max, "mean"=mean, "median"=median)[[normalization_target]]

  # Calculate library sizes for each sample and a normalization factor according to normalize argument
  writeLines("Calculating normalization factor for sample...")
  libsizes_df = tlx_df %>%
    dplyr::ungroup() %>%
    dplyr::mutate(tlx_control=ifelse(tlx_control, "Control", "Treatment")) %>%
    dplyr::group_by(tlx_sample, tlx_group, tlx_control) %>%
    dplyr::summarize(library_size=n(), .groups="keep") %>%
    dplyr::ungroup() %>%
    dplyr::group_by_at(normalize_within_cols) %>%
    dplyr::arrange(tlx_group, tlx_control, tlx_sample) %>%
    dplyr::mutate(library_factor=normalization_target_fun(library_size)/library_size, library_target=ifelse(normalization_target_fun(library_size)==library_size,normalization_target, "")) %>%
    dplyr::ungroup()
  groupsizes_df = libsizes_df %>%
    dplyr::group_by_at(group_cols) %>%
    dplyr::summarize(library_groupsize=sum(library_size*library_factor)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by_at(normalize_between_cols) %>%
    dplyr::mutate(library_groupfactor=normalization_target_fun(library_groupsize)/library_groupsize) %>%
    data.frame()
  libsizes_df = libsizes_df %>%
    dplyr::inner_join(groupsizes_df, by=group_cols) %>%
    dplyr::mutate(library_factor.adj=library_factor*library_groupfactor)
  writeLines(knitr::kable(libsizes_df))


  # Calculate coverage for each sample
  writeLines("Calculating each sample coverage...")
  tlxcov_df = tlx_df %>%
    dplyr::group_by(tlx_group, tlx_group_i, tlx_sample, tlx_control, tlx_path, tlx_strand) %>%
    dplyr::do(tlx_coverage_(., extsize=extsize, exttype=exttype)) %>%
    dplyr::ungroup() %>%
    dplyr::inner_join(libsizes_df %>% dplyr::select(tlx_sample, library_factor), by="tlx_sample") %>%
    dplyr::mutate(tlxcov_pileup.norm=tlxcov_pileup*library_factor)

  # Summarize group coverage by summing all samples in the group with each sample having a weight decided by library size
  writeLines("Adding up coverages from sample(s)...")
  zret = tlxcov_df %>%
    dplyr::group_by_at(group_cols) %>%
    dplyr::do((function(z){
      z
      z_ranges = GenomicRanges::makeGRangesFromDataFrame(z %>% dplyr::mutate(seqnames=tlxcov_chrom, start=tlxcov_start, end=tlxcov_end), ignore.strand=T, keep.extra.columns=T)
      cov_ranges = as(GenomicRanges::coverage(z_ranges, weight=z$tlxcov_pileup.norm), "GRanges")
      ret_df = as.data.frame(cov_ranges) %>%
        dplyr::rename(tlxcov_chrom="seqnames", tlxcov_start="start", tlxcov_end="end", tlxcov_pileup="score") %>%
        dplyr::select(matches("tlxcov_"))
      ret_df
    })(.)) %>%
    dplyr::ungroup()
}

#' @export
tlx_remove_rand_chromosomes = function(tlx_df) {
  tlx_df %>%
    dplyr::filter(Rname %in% paste0("chr", c(1:40, "X", "Y")))
}

#' @export
tlx_mark_rand_chromosomes = function(tlx_df) {
  tlx_df %>%
    dplyr::mutate(tlx_is_rand_chrom = !(Rname %in% paste0("chr", c(1:40, "X", "Y"))))
}

tlx_calc_copynumber = function(tlx_df, bowtie_index, max_hits=500, threads=8, tmp_dir="tmp") {
  tlx_df = tlx_df %>%
    dplyr::mutate(QSeq=substr(Seq, Qstart, Qend))

  dir.create(tmp_dir, recursive=T, showWarnings=F)
  qnames_hash = openssl::md5(paste0(tlx_df$Qname, collapse=""))
  qseq_fasta = file.path(tmp_dir, paste0(qnames_hash, ".fa"))
  qseq_count = file.path(tmp_dir, paste0(qnames_hash, ".count"))
  qseq_cumcount = file.path(tmp_dir, paste0(qnames_hash, ".cumcount"))

  if(!file.exists(qseq_cumcount)) {
    if(!file.exists(qseq_count)) {
      if(!file.exists(qseq_fasta)) {
        writeLines(with(tlx_df, paste0(">", Qname, "\n", QSeq)), con=qseq_fasta)
      }

      cmd = paste0("bowtie2 -f -x ", bowtie_index, " -U ", qseq_fasta ," -k ", max_hits, " --threads ", threads, " -S ", qseq_count)
      system(cmd)
    }

    qseq_count_df = readr::read_tsv(qseq_count, col_names=F, skip=68)
    qseq_cumcount_df = qseq_count_df %>%
      dplyr::group_by(X1) %>%
      dplyr::summarize(n=n())%>%
      setNames(c("Qname", "tlx_copynumber"))
    readr::write_tsv(qseq_cumcount_df, file=qseq_cumcount)
  } else {
    qseq_cumcount_df = readr::read_tsv(qseq_cumcount)
  }

  tlx_df %>%
    dplyr::select(-dplyr::matches("tlx_copynumber")) %>%
    dplyr::left_join(qseq_cumcount_df, by="Qname")
}

#' @export
tlx_mark_dust = function(tlx_df, tmp_dir="tmp") {
  dir.create(tmp_dir, recursive=T, showWarnings=F)
  qnames_hash = openssl::md5(paste0(tlx_df$Qname, collapse=""))
  qnames_fasta = file.path(tmp_dir, paste0(qnames_hash, ".fa"))
  qnames_dust = file.path(tmp_dir, paste0(qnames_hash, ".dust"))

  if(!file.exists(qnames_dust)) {
    if(!file.exists(qnames_fasta)) {
      writeLines("Writing temporary fasta file with sequences...")
      writeLines(with(tlx_df, paste0(">", Qname, "\n", Seq)), con=qnames_fasta)
    }
    writeLines("Using dustmasker to calculate low complexity regions...")
    system(paste0("dustmasker -in ", qnames_fasta, " -out ", qnames_dust, " -outfmt acclist"))
  }
  tlx_dust_df = readr::read_tsv(qnames_dust, col_names=c("Qname", "dust_start", "dust_end")) %>%
    dplyr::mutate(dust_length=dust_end-dust_start+1, Qname=gsub("^>", "", Qname)) %>%
    dplyr::group_by(Qname) %>%
    dplyr::summarize(tlx_dust_dust_regions=n(), tlx_dust_length=sum(dust_length), tlx_dust_coordinates=paste0(dust_start, "-", dust_end, collapse="; ")) %>%
    dplyr::mutate(tlx_has_dust=T) %>%
    dplyr::ungroup()
  tlx_df = tlx_df %>%
    dplyr::select(-dplyr::matches("tlx_dust_dust_regions|tlx_dust_length|tlx_dust_prop|tlx_dust_coordinates|tlx_has_dust")) %>%
    dplyr::left_join(tlx_dust_df, by="Qname") %>%
    dplyr::mutate(tlx_has_dust=tidyr::replace_na(tlx_has_dust, F), tlx_dust_prop=tlx_dust_length/Seq_length)
  tlx_df
}

#' @export
tlx_identify_baits = function(tlx_df, breaksite_size=19, genome_fasta="") {
  if(is.null(tlx_df) || nrow(tlx_df)==0) {
    return(data.frame(bait_sample=NA, bait_chrom=NA, bait_strand=NA, bait_start=NA, bait_end=NA) %>% dplyr::slice(0))
  }

  baits_df = tlx_df %>%
    dplyr::group_by(tlx_group, tlx_sample, B_Rname, B_Strand) %>%
    dplyr::do((function(z){
      misprimed_max = max(z$misprimed-z$uncut)
      if(z$B_Strand[1]<0) bait_start = unique(z$B_Rstart + z$misprimed - misprimed_max-2)
      else bait_start = unique(z$B_Rend - z$misprimed + misprimed_max)
      data.frame(bait_start=bait_start, bait_end=bait_start+breaksite_size - 1)
    })(.)) %>%
    dplyr::mutate(bait_strand=ifelse(B_Strand<0, "-", "+")) %>%
    dplyr::ungroup() %>%
    dplyr::select(bait_group=tlx_group, bait_sample=tlx_sample, bait_chrom=B_Rname, bait_strand, bait_start, bait_end)

  if(genome_fasta!="") {
    if(!file.exists(genome_fasta)) {
      log("Could not find genome file '", genome_fasta, "'")
    }
    baits_ranges = GenomicRanges::makeGRangesFromDataFrame(baits_df %>% dplyr::select(seqnames=bait_chrom, start=bait_start, end=bait_end, strand=bait_strand))
    baits_df$bait_sequence = get_seq(genome_fasta, baits_ranges)$sequence
  }

  baits_df
}

#' @export
tlx_test_hits = function(tlx_df, hits_ranges, paired_samples=T, paired_controls=T, extsize=10000, exttype="along") {
  validate_exttype(exttype)

  if(exttype[1]=="along") {
    tlx_ranges  = GenomicRanges::makeGRangesFromDataFrame(tlx_df %>% dplyr::mutate(seqnames=Rname, sstart=ifelse(Strand=="-1", Junction-extsize, Junction-1), end=ifelse(Strand=="-1", Junction, Junction+extsize-1)), ignore.strand=T, keep.extra.columns=T)
  } else {
    if(exttype[1]=="symmetrical") {
      tlx_ranges  = GenomicRanges::makeGRangesFromDataFrame(tlx_df %>% dplyr::mutate(seqnames=Rname, start=Junction-ceiling(extsize/2), end=Junction+ceiling(extsize/2)), ignore.strand=T, keep.extra.columns=T)
    } else {
      tlx_ranges  = GenomicRanges::makeGRangesFromDataFrame(tlx_df %>% dplyr::mutate(seqnames=Rname, start=Junction, end=Junction+1), ignore.strand=T, keep.extra.columns=T)
    }
  }

  hits_df = as.data.frame(hits_ranges) %>% dplyr::mutate(compare_chrom=seqnames, compare_start=start, compare_end=end)
  hits_ranges = GenomicRanges::makeGRangesFromDataFrame(hits_df, keep.extra.columns=T)

  hits_ranges_reduced = GenomicRanges::makeGRangesFromDataFrame(as.data.frame(GenomicRanges::reduce(hits_ranges)) %>% dplyr::mutate(compare_chrom=seqnames, compare_start=start, compare_end=end), keep.extra.columns=T)
  hits_reduced_df = as.data.frame(hits_ranges_reduced) %>% dplyr::select(compare_chrom, compare_start, compare_end)

  tlxsum_df = tlx_df %>%
    dplyr::group_by(tlx_sample, .drop=F) %>%
    dplyr::summarize(compare_total=sum(!tlx_is_bait_junction)) %>%
    dplyr::ungroup()

  # Prepare overlap counts table (add compare_n=0 for missing entries)
  counts_df_incomplete = as.data.frame(IRanges::mergeByOverlaps(hits_ranges_reduced, tlx_ranges)) %>%
    dplyr::rename(compare_group="tlx_group", compare_group_i="tlx_group_i", compare_sample="tlx_sample") %>%
    dplyr::group_by(compare_chrom, compare_start, compare_end, compare_group, compare_group_i, compare_sample, tlx_control, .drop=F) %>%
    dplyr::summarize(compare_n=n())
  counts_df = dplyr::bind_rows(
    counts_df_incomplete,
    hits_reduced_df %>%
      dplyr::select(compare_chrom, compare_start, compare_end) %>%
      tidyr::crossing(tlx_df %>% dplyr::distinct(compare_group=tlx_group, compare_group_i=tlx_group_i, compare_sample=tlx_sample, tlx_control)) %>%
      dplyr::mutate(compare_n=0)) %>%
    dplyr::distinct(compare_chrom, compare_start, compare_end, compare_group, compare_group_i, compare_sample, tlx_control, .keep_all=T) %>%
    dplyr::inner_join(tlxsum_df, by=c("compare_sample"="tlx_sample")) %>%
    data.frame()

  #
  # Calculate breaks count adjusted with control (by substracting control breaks)
  #
  counts_df.input = counts_df %>% dplyr::filter(!tlx_control)
  counts_df.control = counts_df %>% dplyr::filter(tlx_control)
  if(paired_controls) {
   normcounts_df = counts_df.input %>%
     dplyr::inner_join(counts_df.control %>% select(compare_chrom, compare_start, compare_end, compare_group, compare_group_i, compare_total.control=compare_total, compare_n.control=compare_n), by=c("compare_chrom", "compare_start", "compare_end", "compare_group", "compare_group_i")) %>%
     dplyr::mutate(compare_n.control_adj=compare_n.control*(compare_total/compare_total.control),  compare_n.norm=compare_n-compare_n.control_adj, compare_frac.norm=compare_n.norm/compare_total, compare_frac=compare_n/compare_total) %>%
     dplyr::arrange(compare_group, compare_group_i)
  } else {
   normcounts_df = counts_df.input %>%
     dplyr::left_join(counts_df.control %>% dplyr::group_by(compare_chrom, compare_start, compare_end, compare_group) %>% summarize(compare_total.control=sum(compare_total), compare_n.control=sum(compare_n)), by=c("compare_chrom", "compare_start", "compare_end", "compare_group")) %>%
     dplyr::mutate(compare_total.control=ifelse(!is.na(compare_total.control), compare_total.control, compare_total), compare_n.control=tidyr::replace_na(compare_n.control, 0)) %>%
     dplyr::mutate(compare_n.control_adj=compare_n.control*(compare_total/compare_total.control),  compare_n.norm=compare_n-compare_n.control_adj, compare_frac.norm=compare_n.norm/compare_total, compare_frac=compare_n/compare_total) %>%
     dplyr::arrange(compare_group, compare_group_i)
  }

  if(paired_samples)
  {
    normcounts_sum_df = normcounts_df %>%
      dplyr::select(compare_chrom, compare_start, compare_end, compare_group, compare_group_i, compare_n.norm, compare_n, compare_total, compare_n.control, compare_total.control, compare_n.control_adj)
  } else {
    normcounts_sum_df = normcounts_df %>%
      dplyr::group_by(compare_chrom, compare_start, compare_end, compare_group) %>%
      dplyr::summarize(compare_group_i=1, compare_n=sum(compare_n), compare_n.norm=sum(compare_n.norm), compare_total=sum(compare_total), compare_n.control=sum(compare_n.control), compare_total.control=sum(compare_total.control)) %>%
      dplyr::mutate(compare_n.control_adj=compare_n.control*(compare_total/compare_total.control))
  }

  if(length(unique(tlx_df$tlx_group)))
  {
    z_sum.test = as.data.frame(t(apply(combn(unique(normcounts_sum_df$compare_group), 2), 2, sort))) %>%
      dplyr::rename(compare_group1="V1", compare_group2="V2") %>%
      dplyr::inner_join(normcounts_sum_df %>% dplyr::rename(compare_n.norm1="compare_n.norm", compare_n1="compare_n", compare_total1="compare_total", compare_n.control1="compare_n.control", compare_n.control_adj1="compare_n.control_adj", compare_total.control1="compare_total.control"), by=c("compare_group1"="compare_group")) %>%
      dplyr::inner_join(normcounts_sum_df %>% dplyr::rename(compare_n.norm2="compare_n.norm", compare_n2="compare_n", compare_total2="compare_total", compare_n.control2="compare_n.control", compare_n.control_adj2="compare_n.control_adj", compare_total.control2="compare_total.control"), by=c("compare_group2"="compare_group", "compare_group_i", "compare_chrom", "compare_start", "compare_end")) %>%
      dplyr::group_by(compare_group1, compare_group2, compare_chrom, compare_start, compare_end) %>%
      dplyr::do((function(z){
        zz<<-z

        z.groups = c(z$compare_group1[1], z$compare_group2[1])
        z.fold = mean(z$compare_n1/z$compare_total1) / mean(z$compare_n2/z$compare_total2)

        # Reapeated measures ANOVA
        z.test_data = z %>%
          reshape2::melt(measure.vars=c("compare_n1", "compare_n.control_adj1", "compare_n2", "compare_n.control_adj2")) %>%
          dplyr::select(compare_group_i, treatment=variable, breaks=value) %>%
          dplyr::mutate(group=z.groups[as.numeric(gsub(".*([0-9])$", "\\1", treatment))], treatment=gsub("([0-9])$", "", treatment)) %>%
          dplyr::mutate(treatment=c(compare_n.control="control", compare_n="treatment", compare_n.control_adj="control")[treatment]) %>%
          dplyr::mutate(group=factor(group), treatment=factor(treatment), compare_group_i=factor(compare_group_i)) %>%
          tibble::tibble()
        z.aov = rstatix::anova_test(data=z.test_data, dv=breaks, wid=compare_group_i, within=c(treatment, group))
        z.aov_pval = data.frame(z.aov) %>% dplyr::filter(Effect=="group") %>% .$p

        i.contignency = lapply(split(z, 1:nrow(z)), function(y) matrix(as.numeric(y[c("compare_n1", "compare_n2", "compare_total1", "compare_total2")]), ncol=2))
        if(length(i.contignency)>=2) {
          i.contignency = abind::abind(i.contignency, along=3)
          i.test = mantelhaen.test(i.contignency)
        } else {
          i.test = fisher.test(i.contignency[[1]])
        }

        z %>%
          dplyr::slice(1) %>%
          dplyr::mutate(compare_pvalue=i.test$p.value, compare_odds=i.test$estimate, compare_aov_pvalue=z.aov_pval, compare_fold=z.fold) %>%
          dplyr::select(compare_group1, compare_group2, compare_chrom, compare_start, compare_end, compare_pvalue, compare_odds, compare_aov_pvalue, compare_fold)
      })(.)) %>%
      data.frame()
  } else {
    compare_test.cols = readr::cols(compare_group1=readr::col_character(), compare_group2=readr::col_character(), compare_chrom=readr::col_character(), compare_start=readr::col_double(), compare_end=readr::col_double(), compare_pvalue=readr::col_double(), compare_odds=readr::col_double(), compare_aov_pvalue=readr::col_double(), compare_fold=readr::col_double())
    z_sum.test = blank_tibble(compare_test.cols)
  }

  list(test=z_sum.test, data=normcounts_df)
}

#' @export
tlx_mark_bait_chromosome = function(tlx_df) {
  tlx_df %>%
    dplyr::select(-dplyr::matches("tlx_is_bait_chromosome")) %>%
    dplyr::mutate(tlx_is_bait_chromosome=B_Rname==Rname)
}

#' @export
tlx_mark_bait_junctions = function(tlx_df, bait_region) {
  tlx_df %>%
    dplyr::select(-dplyr::matches("tlx_is_bait_junction")) %>%
    dplyr::mutate(tlx_is_bait_junction=B_Rname==Rname & (abs(B_Rstart-Rstart)<=bait_region/2 | abs(Rend-B_Rend)<=bait_region/2))
}

#' @export
tlx_mark_offtargets = function(tlx_df, offtarget2bait_df) {
  # @todo: Change 100 to something meaningful?
  tlx_df$tlx_id = 1:nrow(tlx_df)
  tlx_bait_ranges = GenomicRanges::makeGRangesFromDataFrame(tlx_df %>% dplyr::mutate(seqnames=B_Rname, start=B_Rstart-50000, end=B_Rend+50000), keep.extra.columns=T, ignore.strand=T)
  tlx_junc_ranges = GenomicRanges::makeGRangesFromDataFrame(tlx_df %>% dplyr::mutate(seqnames=Rname, start=Rstart-50000, end=Rend+50000), keep.extra.columns=T, ignore.strand=T)
  offtarget2bait_bait_ranges = GenomicRanges::makeGRangesFromDataFrame(offtarget2bait_df %>% dplyr::mutate(seqnames=bait_chrom, start=bait_start, end=bait_end), ignore.strand=T)
  offtarget2bait_offt_ranges = GenomicRanges::makeGRangesFromDataFrame(offtarget2bait_df %>% dplyr::mutate(seqnames=offtarget_chrom, start=offtarget_start, end=offtarget_end), ignore.strand=T)

  tlx_offtarget_ids = as.data.frame(IRanges::findOverlaps(tlx_bait_ranges, offtarget2bait_bait_ranges)) %>%
    dplyr::rename(tlx_id="queryHits", o2b_id="subjectHits") %>%
    dplyr::inner_join(as.data.frame(IRanges::findOverlaps(tlx_junc_ranges, offtarget2bait_offt_ranges)), by=c(tlx_id="queryHits", o2b_id="subjectHits")) %>%
    dplyr::distinct(tlx_id) %>%
    .$tlx_id

  tlx_df$tlx_is_offtarget = tlx_df$tlx_id %in% tlx_offtarget_ids

  tlx_df
}

#' @export
tlx_mark_repeats = function(tlx_df, repeatmasker_df) {
  # @todo: make group_by faster using data.table
  repeatmasker_ranges = GenomicRanges::makeGRangesFromDataFrame(repeatmasker_df %>% dplyr::mutate(seqnames=repeatmasker_chrom, start=repeatmasker_start, end=repeatmasker_end), keep.extra.columns=T)
  tlx_df = tlx_df %>%
    dplyr::mutate(tlx_id=1:n()) %>%
    dplyr::select(-dplyr::matches("tlx_repeatmasker_")) %>%
    dplyr::mutate(seqnames=Rname, start=Rstart, end=Rend)
  tlx_ranges = GenomicRanges::makeGRangesFromDataFrame(tlx_df, keep.extra.columns=T, ignore.strand=T)
  r1 = as.data.frame(IRanges::findOverlaps(tlx_ranges, repeatmasker_ranges)) %>%
    dplyr::inner_join(repeatmasker_df, by=c("subjectHits"="repeatmasker_id"))

  data.table::setDT(r1)[,list(tlx_repeatmasker_class=paste0(unique(repeatmasker_class),collapse=", ")), by=list(queryHits)] %>%
    dplyr::right_join(tlx_df, by=c("queryHits"="tlx_id")) %>%
    dplyr::select(-queryHits) %>%
    data.frame()
  # data.table::setDT(r1)[,.(tlx_repeatmasker_class=paste0(unique(repeatmasker_class),collapse=", ")), by = .(queryHits)] %>%
  #   dplyr::right_join(tlx_df, by=c("queryHits"="tlx_id")) %>%
  #   dplyr::select(-queryHits) %>%
  #   data.frame()
}

#' @export
tlx_macs2 = function(tlx_df, effective_size, maxgap=NULL, qvalue=0.01, pileup=1, extsize=2000, slocal=50000, llocal=10000000, exclude_bait_region=F, exclude_repeats=F, exclude_offtargets=F, exttype, grouping) {
  if(exclude_bait_region && !("tlx_is_bait_junction" %in% colnames(tlx_df))) {
    stop("tlx_is_bait_junction is not found in tlx data frame")
  }

  validate_group(grouping)
  validate_exttype(exttype)

  macs2_tlx_df = tlx_df

  if(exclude_offtargets) {
    if(!("tlx_is_offtarget" %in% colnames(macs2_tlx_df))) {
      stop("tlx_is_offtarget is not found in tlx data frame")
    }
    macs2_tlx_df = macs2_tlx_df %>% dplyr::filter(!tlx_is_offtarget)
  }
  if(exclude_repeats) {
    if(!("tlx_repeatmasker_class" %in% colnames(macs2_tlx_df))) {
      stop("tlx_repeatmasker_class is not found in tlx data frame")
    }
    macs2_tlx_df = macs2_tlx_df %>% dplyr::filter(is.na(tlx_repeatmasker_class))
  }

  macs2_tlx_df = macs2_tlx_df %>%
    dplyr::filter(!exclude_bait_region | !tlx_is_bait_junction) %>%
    dplyr::mutate(bed_strand=ifelse(Strand=="1", "-", "+"))

  # @TODO: I think macs does this internally (NO!)
  if(exttype[1]=="along") {
    macs2_tlx_df = macs2_tlx_df %>% dplyr::mutate(bed_start=ifelse(Strand=="-1", Junction-extsize, Junction-1), bed_end=ifelse(Strand=="-1", Junction, Junction+extsize-1))
  } else {
    if(exttype[1]=="symmetrical") {
      macs2_tlx_df = macs2_tlx_df %>% dplyr::mutate(bed_start=Junction-ceiling(extsize/2), bed_end=Junction+ceiling(extsize/2))
    } else {
      macs2_tlx_df = macs2_tlx_df %>% dplyr::mutate(bed_start=Junction, bed_end=Junction+1)
    }
  }

  if(is.null(maxgap) || maxgap==0 || maxgap=="") {
    maxgap = NULL
  }

  if(grouping=="sample") {
    macs2_tlx_df$grouping = paste(macs2_tlx_df$tlx_group, macs2_tlx_df$tlx_group_i)
  }
  if(grouping=="group") {
    macs2_tlx_df$grouping = macs2_tlx_df$tlx_group
  }

  macs_df.all = data.frame()
  for(gr in unique(macs2_tlx_df$grouping)) {
    tlx_df.gr = macs2_tlx_df %>% dplyr::filter(grouping==gr)

    f_input_bed = tempfile()
    f_control_bed = tempfile()
    # f_input_bed = "tmp/input.bed"
    # f_control_bed = "tmp/control.bed"

    tlx_df.gr %>%
      dplyr::filter(!tlx_control) %>%
      dplyr::select(Rname, bed_start, bed_end, Qname, 0, bed_strand) %>%
      readr::write_tsv(file=f_input_bed, na="", col_names=F)

    if(any(tlx_df.gr$tlx_control)) {
      tlx_df.gr %>%
        dplyr::filter(tlx_control) %>%
        dplyr::select(Rname, bed_start, bed_end, Qname, 0, bed_strand) %>%
        readr::write_tsv(file=f_control_bed, na="", col_names=F)

      log("Running MACS with control")
      print("Asdad")
      macs_df = macs2(name=basename(f_input_bed), sample=f_input_bed, control=f_control_bed, maxgap=maxgap, effective_size=length(unique(tlx_df.gr$tlx_path))*effective_size, extsize=extsize, qvalue=qvalue, slocal=slocal, llocal=llocal, output_dir=dirname(f_input_bed))
    } else {
      log("Running MACS without control")
      macs_df = macs2(name=basename(f_input_bed), sample=f_input_bed, maxgap=maxgap, effective_size=effective_size, extsize=extsize, qvalue=qvalue, slocal=slocal, llocal=llocal, output_dir=dirname(f_input_bed))
    }

    macs_df$macs_group = tlx_df.gr$tlx_group[1]
    if(grouping=="sample") {
      macs_df$tlx_group_i = tlx_df.gr$tlx_group_i[1]
    }

    macs_df.all = rbind(macs_df.all, macs_df %>% dplyr::mutate(macs_group=gr))
  }

  macs_df.all = macs_df.all %>% dplyr::filter(macs_pileup>=pileup)

  macs_df.all
}

plot_logos_coordinates = function(fastq_paths, sample_names, widths=list("Beginning"=c(1,30), "End"=c(-25, -1))) {
  plots = list()
  for(i in 1:length(fastq_paths)) {
    sample_path = fastq_paths[i]
    sample_name = sample_names[i]
    writeLines(paste0(i, "/", length(fastq_paths), " : ", sample_name, "    ", sample_path))

    fasta = ShortRead::readFastq(sample_path)
    fasta_reads = ShortRead::sread(fasta)
    plist = list()
    plist[[1]] = cowplot::ggdraw() +
      cowplot::draw_label(sample_names[i], size=10)
    for(wname in names(widths)) {
      writeLines(paste0("  ", wname))

      # If width is longer for all or some sequences then extend the sequences to the needed length
      fasta_reads_long = fasta_reads
      fasta_widths = Biostrings::width(fasta_reads_long)
      fasta_ends = ifelse(fasta_widths>widths[[wname]][2], widths[[wname]][2], fasta_widths)
      missing_nchar = widths[[wname]][2] - fasta_ends
      if(any(missing_nchar>0)) {
        fasta_missing = Biostrings::DNAStringSet(sapply(missing_nchar[missing_nchar>0], function(x) paste(replicate(x, expr="N"), collapse="")))
        fasta_reads_long[missing_nchar>0] = Biostrings::xscat(fasta_reads_long[missing_nchar>0], fasta_missing)
      }

      fasta_w = Biostrings::subseq(fasta_reads_long, start=widths[[wname]][1], end=widths[[wname]][2])
      plist[[length(plist)+1]] = ggplot() +
          ggseqlogo::geom_logo(as.character(fasta_w)) +
          labs(title=wname) +
          ggseqlogo::theme_logo() +
          guides(fill="none") +
          theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
    }


    defaultW = getOption("warn")
    options(warn = -1)
    p = cowplot::plot_grid(plotlist=plist, ncol=1+length(widths), rel_widths=c(1,3,3))
    options(warn = defaultW)
    plots[[sample_name]] = p
  }
  plots
}