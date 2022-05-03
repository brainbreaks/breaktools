validate_group = function(group) {
  valid_group = c("all", "group", "treatment", "none")
  if(!(group %in% valid_group)) { stop(simpleError(paste0("group should be one of: ", paste(valid_group, collapse=", "), " (Actual:", group, ")"))) }
}

validate_group_within = function(within) {
  within_group = c("group", "treatment", "none")
  if(!(within %in% within_group)) { stop(simpleError(paste0("group should be one of: ", paste(within_group, collapse=", "), " (Actual:", within, ")"))) }
}

validate_group_between = function(within) {
  within_group = c("all", "group", "treatment", "none")
  if(!(within %in% within_group)) { stop(simpleError(paste0("group should be one of: ", paste(within_group, collapse=", "), " (Actual:", within, ")"))) }
}

validate_bed_mode = function(mode) {
  valid_mode = c("junction", "alignment")
  if(!(mode %in% valid_mode)) { stop(simpleError(paste0("group should be one of: ", paste(valid_mode, collapse=", "), " (Actual:", mode, ")"))) }
}

validate_normalization_target = function(target) {
  valid_targets = c("min", "max")
  if(!(target %in% valid_targets)) { stop(simpleError(paste0("normalization target should be one of: ", paste(valid_targets, collapse=", "), " (Actual:", target, ")"))) }
}


validate_exttype = function(exttype) {
  valid_exttype = c("symmetrical", "along", "opposite", "none")
  if(!(exttype %in% valid_exttype)) { stop(simpleError(paste0("exttype should be one of: ", paste(valid_exttype, collapse=", "), " (Actual:", exttype, ")"))) }
}

tlx_get_group_cols = function(group, ignore.strand=T, ignore.control=F) {
  validate_group(group)
  if(group=="all") group_cols = c("")[-1]
  if(group=="group") group_cols = c("tlx_group")
  if(group=="treatment") group_cols = c("tlx_group", "tlx_control")
  if(group %in% c("sample", "none")) group_cols = c("tlx_group", "tlx_group_i", "tlx_sample")

  if(!ignore.control) group_cols = unique(c(group_cols, "tlx_control"))
  if(!ignore.strand) group_cols = unique(c(group_cols, "tlx_strand"))

  return(group_cols)
}

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

tlx_blank = function() {
  blank_tibble(tlx_cols()) %>%
    dplyr::mutate(tlx_sample=NA_character_, tlx_path=NA_character_, tlx_group=NA_character_, tlx_control=NA) %>%
    dplyr::mutate(tlx_is_bait_chromosome=NA, tlx_is_bait_junction=NA, tlx_is_offtarget=NA) %>%
    dplyr::slice(0)
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

  generated_path = do.call(paste, c(filenames_df, sep = "_"))
  unique_generated_path = unique(generated_path)
  map_unique_generated_path = gsub("[/.]", "-", unique_generated_path)
  map_unique_generated_path = gsub("[^A-Za-z0-9 -]+", "_", map_unique_generated_path)
  map_unique_generated_path = gsub("(_| )+", "\\1", map_unique_generated_path)
  names(map_unique_generated_path) = unique_generated_path
  unname(map_unique_generated_path[generated_path])
}

#' @title tlx_read_samples
#' @export
#' @description Read tab-separated-file with description and location of TLX files (generated with HTGTS), one sample per row
#'
#' @param annotation_path Path of tab-separated-file with samples information
#' @param samples_path Path to a folder where actual TLX files are located in case the information in samples tab-separated-information doesn't contain a full path
#'
#' @details
#' Original HTGTS pipeline is hosted on \href{https://github.com/robinmeyers/transloc_pipeline}{github}. A docker image is available
#' through \href{https://hub.docker.com/repository/docker/sandrejev/htgts}{DockerHUB} (docker://sandrejev/htgts:latest)
#'
#' The file should contain following columns:
#'   \strong{path} - Path to the TLX file. TLX file path is appended to samples_path to get an actual path to TLX sample
#'   \strong{sample} - Name of the sample
#'   \strong{control} - A boolean column that specifies whether sample is a control sample or not (values TRUE/FALSE)
#'   \strong{group} - Name of a group to which a sample belongs (for example it can be different experiments)
#'
#' @seealso tlx_read tlx_read_many
#'
#' @return A data frame with sample information
tlx_read_samples = function(annotation_path, samples_path=".") {
  samples_df = readr::read_tsv(annotation_path, comment="#") %>%
    dplyr::mutate(path=file.path(samples_path, path)) %>%
    dplyr::mutate(tlx_exists=file.exists(path))
  samples_df
}

#' @title tlx_read
#' @export
#' @description Read TLX file (generated with HTGTS)
#'
#' @param sample Path to the TLX file
#' @param sample Name of the sample
#' @param group Name of a group to which a sample belongs (for example it can be different experiments)
#' @param group_i Index of the sample inside the group
#' @param control A boolean value specifying whether sample is a control sample or not
#'
#' @details
#' Original HTGTS pipeline is hosted on \href{https://github.com/robinmeyers/transloc_pipeline}{github}. A docker image is available
#' through \href{https://hub.docker.com/repository/docker/sandrejev/htgts}{DockerHUB} (docker://sandrejev/htgts:latest)
#'
#' @seealso tlx_read_many tlx_read_samples
#'
#' @return Data frame from TLX file
tlx_read = function(path, sample, group="", group_i=1, control=F) {
  tlx_single_df = readr::read_tsv(path, comment="#", skip=1, col_names=names(tlx_cols()$cols), col_types=tlx_cols()) %>%
    dplyr::mutate(tlx_strand=as.character(ifelse(Strand<0, "-", "+"))) %>%
    dplyr::mutate(Seq_length=as.numeric(nchar(Seq)), tlx_sample=as.character(sample), tlx_path=as.character(path), tlx_group=as.character(group), tlx_group_i=as.numeric(group_i), tlx_control=as.logical(control)) %>%
    dplyr::mutate(tlx_duplicated=duplicated(paste0(Rname, B_Rstart, B_Rend, ifelse(tlx_strand=="+", Rstart, Rend)))) %>%
    dplyr::mutate(QSeq=substr(Seq, Qstart, Qend))
}


#' @title tlx_read_many
#' @export
#' @description Read multiple TLX files (generated with HTGTS) using information in samples data, one sample per row
#'
#' @param samples_df Path of tab-separated-file with samples information
#' @param threads Use multiple threads to read TLX files in paralel
#'
#' @details
#' Original HTGTS pipeline is hosted on \href{https://github.com/robinmeyers/transloc_pipeline}{github}. A docker image is available
#' through \href{https://hub.docker.com/repository/docker/sandrejev/htgts}{DockerHUB} (docker://sandrejev/htgts:latest)
#'
#' Data frame with sample information should contain following columns:
#'   \strong{path} - Path to the TLX file. TLX file path is appended to samples_path to get an actual path to TLX sample
#'   \strong{sample} - Name of the sample
#'   \strong{control} - A boolean column that specifies whether sample is a control sample or not (values TRUE/FALSE)
#'   \strong{group} - Name of a group to which a sample belongs (for example it can be different experiments)
#'
#' @seealso tlx_read tlx_read_samples
#'
#' @return Data frame from multiple TLX file
tlx_read_many = function(samples_df, threads=1) {
  if(!all(samples_df$tlx_exists)) {
    files_str = paste(samples_df %>% dplyr::filter(!file.exists(samples_df$path)) %>% .$path, collapse="\n")
    stop(paste0("Some TLX files do not exist: \n", files_str))
  }

  tlx_df.all = data.frame()
  if(threads > 1) {
    doParallel::registerDoParallel(cores=threads)
    tlx_df.all = foreach(f=1:nrow(samples_df)) %dopar% {

      df = tlx_read(
        path=samples_df$path[f],
        sample=samples_df$sample[f],
        control=samples_df$control[f],
        group=samples_df$group[f],
        group_i=ifelse("group_i" %in% colnames(samples_df), samples_df$group_i[f], 1))
      log("Read tlx file ", f, "/", nrow(samples_df), ": ",  samples_df$path[f])
      df
    }
    log("Merging TLX files...")
    tlx_df.all = do.call(dplyr::bind_rows, tlx_df.all)
  } else {
    for(f in 1:nrow(samples_df)) {
      log("Reading tlx file ", f, "/", nrow(samples_df), ": ",  samples_df$path[f])
      tlx_df.f = tlx_read(
        path=samples_df$path[f],
        sample=samples_df$sample[f],
        control=samples_df$control[f],
        group=samples_df$group[f],
        group_i=ifelse("group_i" %in% colnames(samples_df), samples_df$group_i[f], 1))
      tlx_df.all = dplyr::bind_rows(tlx_df.all, tlx_df.f)
    }
  }

  tlx_df.all %>%
    dplyr::inner_join(samples_df, by=c("tlx_sample"="sample"))
}

#' @export
#'
tlxcov_write_bedgraph = function(tlxcov_df, path, group) {
  ignore.strand = !("tlx_strand" %in% colnames(tlxcov_df))

  writeLines("calculating filenames(s)...")
  if(group=="all") tlxcov_df = tlxcov_df %>% dplyr::mutate(g=paste0(dirname(path), "/", basename(path), "-", tlx_generate_filename_col(., include_group=F, include_sample=F, include_treatment=T, include_strand=!ignore.strand), ".bedgraph"))
  if(group=="group") tlxcov_df = tlxcov_df %>% dplyr::mutate(g=paste0(dirname(path), "/", basename(path), "-", tlx_generate_filename_col(., include_group=T, include_sample=F, include_treatment=T, include_strand=!ignore.strand), ".bedgraph"))
  if(group=="sample") tlxcov_df = tlxcov_df %>% dplyr::mutate(g=paste0(dirname(path), "/", basename(path), "-", tlx_generate_filename_col(., include_group=F, include_sample=T, include_treatment=T, include_strand=!ignore.strand), ".bedgraph"))
  if(!ignore.strand) tlxcov_df = tlxcov_df %>% dplyr::mutate(tlxcov_pileup=ifelse(tlx_strand=="+", 1, -1)*tlxcov_pileup)


  writeLines("Writing bedgraph file(s)...")
  if(!dir.exists(dirname(path))) dir.create(dirname(path), recursive=T)
  tlxcov_df %>%
    dplyr::group_by(g) %>%
    dplyr::do((function(z){
      z.out = z %>% dplyr::select(tlxcov_chrom, tlxcov_start, tlxcov_end, tlxcov_pileup)
      z.path = z$g[1]
      writeLines(paste0("Writing to file '", z.path, "'"))
      readr::write_tsv(z.out, file=z.path, col_names=F)
      data.frame()
    })(.))

  tlxcov_df %>%
    dplyr::select(bedgraph_path=g, dplyr::matches("tlx_group|tlx_sample|tlx_control")) %>%
    dplyr::distinct(bedgraph_path, .keep_all=T)
}

test = function()
{
  devtools::load_all('~/Workspace/breaktools/')
  samples_df = tlx_read_samples("~/Workspace/Datasets/HTGTS/samples/All_samples.tsv", "~/Workspace/Datasets/HTGTS") %>%
    dplyr::filter(grepl("concentration", experiment))

  tlx_df = tlx_read_many(samples_df, threads=30)
  tlx_all_df = tlx_df %>%
    dplyr::mutate(tlx_group=ifelse(tlx_group=="DMSO", "APH 0.2 uM 96h", tlx_group))
  tlx_libfactors = tlx_libfactors(tlx_all_df, group="group", normalize_within="treatment", normalize_between="all")

  samples_df = tlx_read_samples("samples.tsv", "TLX")
  tlx_df = tlx_read_many(samples_df, threads=30)
  tlx_libfactors = tlx_libfactors(tlx_all_df, normalize_within="treatment", normalize_between="treatment")
  tlxcov_df = tlx_coverage(tlx_df, group="treatment", extsize=20e3, exttype="symmertrical", libfactors_df=tlx_libfactors, ignore.strand=T)
}

#' @title tlx_libfactors
#' @export
#' @description Calculate normalization factors for samples
#'
#' @param tlx_df Data frame with information from multiple TLX files
#' @param normalize_within Strategy of data normalization within each group (possible values are:
#'   \strong{none} - do not perform any sample normalization
#'   \strong{group} - normalize all samples within group realative to their size
#'   \strong{treatment} - normalizes samples within each group like group option but analyze control and treatment separately)
#' @param normalize_between Strategy of data normalization between groups or treatments (possible values are:
#'   \strong{none} - do not perform any between-groups normalization
#'   \strong{treatment} - normalize treatments within each group
#'   \strong{group} - normalizes all groups and treatments accross all dataset
#' @param normalization_target Name of the function that should be used for normalization (\strong{group} - normalize to smallest group/sample, \strong{max} - normalize to largest group/sample)
#'
#' @seealso tlx_coverage
#'
#' @return Data frame with normalization factors for each sample
tlx_libfactors = function(tlx_df, normalize_within, normalize_between, normalization_target="min")
{
  validate_group_within(normalize_within)
  validate_group_between(normalize_between)
  validate_normalization_target(normalization_target)

  if(normalize_within=="group") normalize_within_cols = c("tlx_group")
  if(normalize_within=="treatment") normalize_within_cols = c("tlx_group", "tlx_control")
  if(normalize_within=="none") normalize_within_cols = c("tlx_sample")
  if(normalize_between=="group") normalize_between_cols = c()
  if(normalize_between=="treatment") normalize_between_cols = setdiff(normalize_within_cols, "tlx_control")
  if(normalize_between=="none") normalize_between_cols = normalize_within_cols

  normalization_target_fun = match.fun(normalization_target)

  # Calculate library sizes for each sample and a normalization factor according to normalize argument
  libsizes_df = tlx_df %>%
    dplyr::ungroup() %>%
    dplyr::mutate(tlx_control=ifelse(tlx_control, "Control", "Treatment")) %>%
    dplyr::group_by(tlx_sample, tlx_group, tlx_group_i, tlx_control) %>%
    dplyr::summarize(library_size=dplyr::n(), .groups="keep") %>%
    dplyr::ungroup() %>%
    dplyr::group_by_at(normalize_within_cols) %>%
    dplyr::mutate(library_within_factor=normalization_target_fun(library_size)/library_size, library_target=ifelse(normalization_target_fun(library_size)==library_size,normalization_target, "")) %>%
    dplyr::ungroup()
  groupsizes_df = libsizes_df %>%
    dplyr::group_by_at(normalize_within_cols) %>%
    dplyr::summarize(library_groupsize=sum(library_size*library_within_factor)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by_at(normalize_between_cols) %>%
    dplyr::mutate(library_between_factor=normalization_target_fun(library_groupsize)/library_groupsize) %>%
    data.frame()
  libsizes_df = libsizes_df %>%
    dplyr::inner_join(groupsizes_df, by=intersect(normalize_within_cols, normalize_between_cols)) %>%
    dplyr::mutate(library_factor=library_within_factor*library_between_factor)

  libsizes_df
}

#' @title tlx_write_bed
#' @export
#' @description Export TLX data frame as single, multiple groupped or multiple single-sample BED file(s)
#'
#' @param tlx_df Data frame with information from multiple TLX files
#' @param path Path to a folder for exporting BED file(s)
#' @param group Data groupping strategy
#'   \strong{none} - Export each sample separately
#'   \strong{group} - Export each group separately
#'   \strong{all} - Pull data from all samples together
#' @param mode Specifies which coordinates to use in BED file. Possible values are:
#'   \strong{junction} - Single nucleotide position of a junction
#'   \strong{alignment} - Coordinates of aligned prey sequence
#' @param ignore.strand Boolean value indicating whether a separate BED file for each strand junctions should be created
#' @param include_strand Boolean value indicating whether a separate BED file for each control/treatment should be created
tlx_write_bed = function(tlx_df, path, group="all", mode="junction", ignore.strand=T, ignore.treatment=F) {
  validate_group(group)
  validate_bed_mode(mode)

  if(mode=="junction") tlx_bed_df = tlx_df %>% dplyr::mutate(start=Junction, end=Junction, name=paste0(Qname, " (", tlx_sample, ")"))
  if(mode=="alignment") tlx_bed_df = tlx_df %>% dplyr::mutate(start=Rstart, end=Rend, name=paste0(Qname, " (", tlx_sample, ")"))

  writeLines("calculating filenames(s)...")
  if(group=="all") tlx_bed_df = tlx_bed_df %>% dplyr::mutate(g=tlx_generate_filename_col(., include_group=F, include_sample=F, include_treatment=!ignore.treatment, include_strand=!ignore.strand))
  if(group=="group") tlx_bed_df = tlx_bed_df %>% dplyr::mutate(g=tlx_generate_filename_col(., include_group=T, include_sample=F, include_treatment=!ignore.treatment, include_strand=!ignore.strand))
  if(group=="sample") tlx_bed_df = tlx_bed_df %>% dplyr::mutate(g=tlx_generate_filename_col(., include_group=F, include_sample=T, include_treatment=!ignore.treatment, include_strand=!ignore.strand))

  writeLines("Writing bedgraph file(s)...")
  if(!dir.exists(dirname(path))) dir.create(dirname(path), recursive=T)
  tlx_bed_df %>%
    dplyr::group_by(g) %>%
    dplyr::do((function(z){
      z.out = z %>% dplyr::select(Rname, start, end, name, mapqual, tlx_strand)
      z.path = file.path(dirname(path), paste0(basename(path), "-", z$g[1], ".bed"))
      writeLines(paste0("Writing to file '", z.path, "'"))
      readr::write_tsv(z.out, file=z.path, col_names=F)
      data.frame()
    })(.))

  return(NULL)
}

#' @title tlx_coverage
#' @export
#'
#' @description Calculates coverage (pileup) from tlx data frame. The coverage can be pulled from multiple samples based on \code{group} parameter and
#' normalized base factors specified for each sample in \code{libfactors_df}. The pileup is constructed upon summing overlapping
#' extended (extended based on \code{extsize} and \code{exttype} parameters) junctions.
#'
#' @param tlx_df Data frame with information from multiple TLX files
#' @param group Grouping columns. Possible values are:
#'   \strong{all} - group together all groups
#'   \strong{group} - group together all samples within each group
#'   \strong{treatment} - group together all samples within each group, separately for control and treatment
#'   \strong{none} - do not perform any groupping
#' @param extsize To build coverage (pileup) data frame it is often usefull to extend the size of the junction (which is just one nucleotide representing the double stranded break).
#'    Single nucleotide position is extended to a size large enough to overlap with other junctions, so that density could be calculated
#' @param exttype Type of extension
#'    \strong{symmetrical} - Extention is done around the center of the junction extending by \code{extsize/2} each direction
#'    \strong{along} - Extention is done in the direction corresponding to pray strand and extending \code{extsize} that direction
#'    \strong{opposite} - Extention is done opposite of pray strand direction and extending \code{extsize}
#'    \strong{none} - No extention is performed
#' @param libfactors_df A data frame with \code{tlx_sample} and \code{library_factor} columns used for sample normalization
#' @param ignore.strand Ignore strand information or calculate coverage for each strand individually
#'
#' @return A data frame with coverages for each sample or group
#' @examples
#' samples_df = tlx_read_samples("samples.tsv", "TLX")
#' tlx_df = tlx_read_many(samples_df, threads=30)
#' tlx_libfactors = tlx_libfactors(tlx_all_df, normalize_within="treatment", normalize_between="treatment")
#' tlxcov_df = tlx_coverage(tlx_df, group="treatment", extsize=20e3, exttype="symmertrical", libfactors_df=tlx_libfactors, ignore.strand=T)
tlx_coverage = function(tlx_df, group, extsize, exttype, libfactors_df=NULL, ignore.strand=T) {
  validate_group(group)
  validate_exttype(exttype)

  if(is.null(libfactors_df)) {
    libfactors_df = tlx_df %>%
      dplyr::mutate(library_factor=1) %>%
      dplyr::distinct(tlx_sample, library_factor)
  }

  if(!all(tlx_df$tlx_sample %in% libfactors_df$tlx_sample)) {
    stop(paste0("Samples are missing from libfactors_df data.frame: ", paste(setdiff(tlx_df$tlx_sample, libfactors_df$tlx_sample), collapsse=", ")))
  }

  tlx_coverage_ = function(x, extsize, exttype) {
    if(exttype[1]=="along") {
      x_ranges  = GenomicRanges::makeGRangesFromDataFrame(x %>% dplyr::mutate(seqnames=Rname, start=ifelse(tlx_strand=="-", Junction-extsize, Junction-1), end=ifelse(tlx_strand=="-", Junction, Junction+extsize-1)) %>% dplyr::select(seqnames, start, end), ignore.strand=T, keep.extra.columns=T)
    } else {
      if(exttype[1]=="opposite") {
        x_ranges  = GenomicRanges::makeGRangesFromDataFrame(x %>% dplyr::mutate(seqnames=Rname, start=ifelse(tlx_strand=="+", Junction-extsize, Junction-1), end=ifelse(tlx_strand=="+", Junction, Junction+extsize-1)) %>% dplyr::select(seqnames, start, end), ignore.strand=T, keep.extra.columns=T)
      } else {
        if(exttype[1]=="symmetrical") {
          x_ranges  = GenomicRanges::makeGRangesFromDataFrame(x %>% dplyr::mutate(seqnames=Rname, start=Junction-ceiling(extsize/2), end=Junction+ceiling(extsize/2)) %>% dplyr::select(seqnames, start, end), ignore.strand=T, keep.extra.columns=T)
        } else {
          x_ranges  = GenomicRanges::makeGRangesFromDataFrame(x %>% dplyr::mutate(seqnames=Rname, start=Junction, end=Junction+1) %>% dplyr::select(seqnames, start, end), ignore.strand=T, keep.extra.columns=T)
        }
      }
    }

    cov_ranges = as(GenomicRanges::coverage(x_ranges), "GRanges")
    ret_df = as.data.frame(cov_ranges) %>%
      dplyr::rename(tlxcov_chrom="seqnames", tlxcov_start="start", tlxcov_end="end", tlxcov_pileup="score") %>%
      dplyr::inner_join(x %>% dplyr::distinct(Rname, tlxcov_is_bait_chrom=tlx_is_bait_chrom), by=c("tlxcov_chrom"="Rname")) %>%
      dplyr::select(matches("tlxcov_"))
    ret_df
  }

  group_cols = tlx_get_group_cols(group, ignore.strand=ignore.strand)

  # Calculate coverage for each sample
  writeLines("Calculating each sample coverage...")
  tlxcov_df = tlx_df %>%
    dplyr::group_by(tlx_group, tlx_group_i, tlx_sample, tlx_control, tlx_path, tlx_strand) %>%
    dplyr::do(tlx_coverage_(., extsize=extsize, exttype=exttype)) %>%
    dplyr::ungroup()
  tlxcov_df = tlxcov_df %>%
    dplyr::left_join(libfactors_df %>% dplyr::select(tlx_sample, library_factor), by="tlx_sample") %>%
    dplyr::mutate(tlxcov_pileup.norm=tlxcov_pileup*library_factor)

  # Summarize group coverage by summing all samples in the group with each sample having a weight decided by library size
  writeLines("Adding up coverages from sample(s)...")
  tlxcov_df %>%
    dplyr::group_by_at(group_cols) %>%
    dplyr::do((function(z){
      z
      z_ranges = GenomicRanges::makeGRangesFromDataFrame(z %>% dplyr::mutate(seqnames=tlxcov_chrom, start=tlxcov_start, end=tlxcov_end), ignore.strand=T, keep.extra.columns=T)
      cov_ranges = as(GenomicRanges::coverage(z_ranges, weight=z$tlxcov_pileup.norm), "GRanges")
      ret_df = as.data.frame(cov_ranges) %>%
        dplyr::rename(tlxcov_chrom="seqnames", tlxcov_start="start", tlxcov_end="end", tlxcov_pileup="score") %>%
        dplyr::mutate(tlxcov_end=tlxcov_end+1) %>%
        dplyr::select(matches("tlxcov_"))
      ret_df
    })(.)) %>%
    dplyr::ungroup()
}

#' @title tlx_remove_rand_chromosomes
#' @export
#'
#' @param tlx_df Data frame with information from multiple TLX files
#' @description Remove non-conventional WGS sequence scaffolds (like Chr4_JH584294_random)
#' @seealso tlx_mark_rand_chromosomes
#'
#' @return Data frame with rand scaffolds removed
tlx_remove_rand_chromosomes = function(tlx_df) {
  tlx_df %>%
    dplyr::filter(Rname %in% paste0("chr", c(1:40, "X", "Y")))
}

#' @title tlx_mark_rand_chromosomes
#' @export
#'
#' @description Mark non-conventional WGS sequence scaffolds (like Chr4_JH584294_random). Junctions detected in these
#' chromosomes are marked appropriately in \code{tlx_is_rand_chrom} column
#'
#' @param tlx_df Data frame with information from multiple TLX files
#' @seealso tlx_remove_rand_chromosomes
#' @return Data frame with rand scaffolds removed
tlx_mark_rand_chromosomes = function(tlx_df) {
  tlx_df %>%
    dplyr::mutate(tlx_is_rand_chrom = !(Rname %in% paste0("chr", c(1:40, "X", "Y"))))
}


#' @title tlx_calc_copynumber
#' @export
#' @description Calculate copy number for each pray sequence in reference genome and save this value in \code{tlx_copynumber} column. The results are cached.
#'
#' @param tlx_df Data frame with information from multiple TLX files
#' @param bowtie2_index Path to bowtie2 index
#' @param max_hits Maximum number of matches to search for using bowtie2
#' @param threads Number of threads used in bowtie2 when aligning pray sequences
#' @param tmp_dir Path to folder with cached results
#'
#' @return Data frame with copy number in \code{tlx_copynumber} column
tlx_calc_copynumber = function(tlx_df, bowtie2_index, max_hits=500, threads=8, tmp_dir="tmp") {
  if(!bedr::check.binary("bowtie2")) {
    stop("dustmasker executible not found")
  }

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

      writeLines("Using bowtie2 to find number of copies in reference fasta file...")
      cmd = paste0("bowtie2 -f -x ", bowtie2_index, " -U ", qseq_fasta ," -k ", max_hits, " --threads ", threads, " -S ", qseq_count)
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


#' @title tlx_mark_dust
#' @export
#' @description Predict "dust" using dustmasker
#'
#' @param tlx_df Data frame with information from multiple TLX files
#' @param tmp_dir Path to folder with cached results
#'
#' @return Data frame with dust marked using \code{tlx_has_dust} column
tlx_mark_dust = function(tlx_df, tmp_dir="tmp") {
  if(!bedr::check.binary("dustmasker")) {
    stop("dustmasker executible not found")
  }

  dir.create(tmp_dir, recursive=T, showWarnings=F)
  qnames_hash = openssl::md5(paste0(tlx_df$Qname, collapse=""))
  qnames_fasta = file.path(tmp_dir, paste0(qnames_hash, ".fa"))
  qnames_dust = file.path(tmp_dir, paste0(qnames_hash, ".dust"))

  if(!file.exists(qnames_dust)) {
    if(!file.exists(qnames_fasta)) {
      writeLines("Writing temporary fasta file with sequences...")
      writeLines(with(tlx_df, paste0(">", Qname, "\n", QSeq)), con=qnames_fasta)
    }
    writeLines("Using dustmasker to calculate low complexity regions...")
    system(paste0("dustmasker -in ", qnames_fasta, " -out ", qnames_dust, " -outfmt acclist"))
  }

  dust_cols = readr::cols(Qname=readr::col_character(), dust_start=readr::col_double(), dust_end=readr::col_double())
  tlx_dust_df = readr::read_tsv(qnames_dust, col_names=names(dust_cols$cols), col_types=dust_cols) %>%
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
tlx_identify_baits = function(tlx_df, genome_fasta="") {
  if(is.null(tlx_df) || nrow(tlx_df)==0) {
    return(data.frame(bait_sample=NA, bait_chrom=NA, bait_strand=NA, bait_start=NA, bait_end=NA) %>% dplyr::slice(0))
  }

  baits_df = tlx_df %>%
    dplyr::group_by(tlx_group, tlx_sample, tlx_bait_chrom, tlx_bait_strand, tlx_bait_start, tlx_bait_end) %>%
    dplyr::summarize(n=dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(dplyr::desc(n)) %>%
    dplyr::select(bait_group=tlx_group, bait_sample=tlx_sample, bait_chrom=tlx_bait_chrom, bait_strand=tlx_bait_strand, bait_start=tlx_bait_start, bait_end=tlx_bait_end)

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
tlx_extract_bait = function(tlx_df, bait_size, bait_region) {
  tlx_df %>%
    dplyr::select(-dplyr::matches("^(tlx_bait_chrom|tlx_bait_start|tlx_bait_end|tlx_bait_strand|tlx_is_bait_junction)$")) %>%
    dplyr::group_by(tlx_group, tlx_sample, B_Rname, B_Strand) %>%
    dplyr::mutate(misprimed_max=max(misprimed-uncut)) %>%
    dplyr::mutate(tlx_bait_start=ifelse(B_Strand<0, B_Rstart + misprimed - misprimed_max-2, B_Rend - misprimed + misprimed_max)) %>%
    dplyr::mutate(tlx_bait_end=tlx_bait_start+bait_size - 1) %>%
    dplyr::mutate(tlx_bait_chrom=B_Rname, tlx_bait_strand=ifelse(B_Strand<0, "-", "+")) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(tlx_is_bait_chrom=B_Rname==Rname) %>%
    dplyr::mutate(tlx_is_bait_junction=B_Rname==Rname & (abs(tlx_bait_start-Junction)<=bait_region/2 | abs(tlx_bait_end-Junction)<=bait_region/2))
}

#' @export
tlx_mark_offtargets = function(tlx_df, offtargets_df, offtarget_region=1000, bait_region=1000) {
  tlx_df$tlx_id = 1:nrow(tlx_df)
  tlx_bait_ranges = tlx_df %>% df2ranges(tlx_bait_chrom, tlx_bait_start, tlx_bait_end)
  tlx_junc_ranges = tlx_df %>% df2ranges(Rname, Junction, Junction)
  offtargets_bait_ranges = offtargets_df %>% dplyr::mutate(bait_region=bait_region) %>% df2ranges(offtarget_bait_chrom, offtarget_bait_start-bait_region/2, offtarget_bait_end+bait_region/2)
  offtargets_offt_ranges = offtargets_df %>% dplyr::mutate(offtarget_region=offtarget_region) %>% df2ranges(offtarget_chrom, offtarget_start-offtarget_region/2, offtarget_end+offtarget_region/2)

  tlx_offtarget_ids = as.data.frame(IRanges::findOverlaps(tlx_bait_ranges, offtargets_bait_ranges)) %>%
    dplyr::rename(tlx_id="queryHits", o2b_id="subjectHits") %>%
    dplyr::inner_join(as.data.frame(IRanges::findOverlaps(tlx_junc_ranges, offtargets_offt_ranges)), by=c(tlx_id="queryHits", o2b_id="subjectHits")) %>%
    dplyr::distinct(tlx_id) %>%
    .$tlx_id

  tlx_df$tlx_is_offtarget = tlx_df$tlx_id %in% tlx_offtarget_ids | tlx_df$tlx_is_bait_junction

  tlx_df
}


#' @export
tlx_strand_crosscorrelation = function(z, step=1000, min_points=5, negative_correlation=F)
{
  zz<<-z
  z.range = c(min(z$tlxcov_start), max(z$tlxcov_end-1))
  z.sense = z %>%
    dplyr::filter(tlx_strand=="+") %>%
    dplyr::mutate(start=tlxcov_start, end=tlxcov_end-1) %>%
    reshape2::melt(measure.vars=c("start", "end"), value.name="position") %>%
    dplyr::arrange(position)
  z.anti = z %>%
    dplyr::filter(tlx_strand=="-") %>%
    dplyr::mutate(start=tlxcov_start, end=tlxcov_end-1) %>%
    reshape2::melt(measure.vars=c("start", "end"), value.name="position") %>%
    dplyr::arrange(position)

  if(nrow(z.sense)>=min_points & nrow(z.anti)>=min_points) {
    p.sense = approx(x=z.sense$position, y=z.sense$tlxcov_pileup, xout=seq(z.range[1], z.range[2], by=step), yleft=0, yright=0)$y
    p.anti = approx(x=z.anti$position, y=z.anti$tlxcov_pileup, xout=seq(z.range[1], z.range[2], by=step), yleft=0, yright=0)$y

    ccf.sense = ccf(p.sense, p.anti, lag.max=floor(length(p.anti)/2), plot=F)
    if(negative_correlation) {
      f = which.max(abs(ccf.sense$acf))
    } else {
      f = which.max(ccf.sense$acf)
    }

    res = data.frame(crosscorrelation_lag=ccf.sense$lag[f]*step, crosscorrelation_rellag=ccf.sense$lag[f]/length(p.sense), crosscorrelation_value=ccf.sense$acf[f])
  } else {
    res = data.frame(crosscorrelation_lag=NA_real_, crosscorrelation_rellag=NA_real_, crosscorrelation_value=NA_real_)
  }
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
geom_tlxcov = function(x, scale=1) {
  x_area = x %>%
    dplyr::group_by(rdc_chrom, rdc_cluster, rdc_cluster_display)  %>%
    dplyr::summarize(tlxcov_prevend=c(0, tlxcov_end[-dplyr::n()]), tlx_strand, tlxcov_start, tlxcov_end, tlxcov_pileup) %>%
    dplyr::rowwise() %>%
    dplyr::do((function(z) {
      if(z$tlxcov_prevend==0) {
        d = data.frame(tlxcov_pos=c(z$tlxcov_start, z$tlxcov_start, z$tlxcov_end), tlxcov_pileup=c(0,z$tlxcov_pileup,z$tlxcov_pileup), tlx_strand=z$tlx_strand)
      } else {
        if(z$tlxcov_start!=z$tlxcov_prevend) {
          d = data.frame(tlxcov_pos=c(z$tlxcov_prevend, z$tlxcov_start, z$tlxcov_start, z$tlxcov_end), tlxcov_pileup=c(0,0,z$tlxcov_pileup,z$tlxcov_pileup), tlx_strand=z$tlx_strand)
        } else {
         d = data.frame(tlxcov_pos=c(z$tlxcov_start, z$tlxcov_end), tlxcov_pileup=z$tlxcov_pileup, tlx_strand=z$tlx_strand)
        }
      }
      d = cbind(d, rdc_chrom=z$rdc_chrom, rdc_cluster=z$rdc_cluster, rdc_cluster_display=z$rdc_cluster_display)
      d
    })(.)) %>%
    dplyr::ungroup()


    geom_ribbon(aes(x=tlxcov_pos, ymin=0, ymax=tlxcov_pileup*scale, fill=tlx_strand), alpha=0.7, data=x_area)
}


#' @export
tlxcov_macs2 = function(tlxcov_df, group, params) {
  validate_group(group)
  group_cols = tlx_get_group_cols(group, ignore.strand=T, ignore.control=T)

  results_df = tlxcov_df %>%
    dplyr::group_by_at(group_cols) %>%
    dplyr::do((function(z){
      zz<<-z
      tlxcov_ranges = z %>% dplyr::mutate(score=tlxcov_pileup) %>% df2ranges(tlxcov_chrom, tlxcov_start, tlxcov_end)
      sample_ranges = tlxcov_ranges[!tlxcov_ranges$tlx_control]
      control_ranges = tlxcov_ranges[tlxcov_ranges$tlx_control]

      group_desc = paste0(group_cols, rep("=", length(group_cols)), sapply(distinct(z[,group_cols]), as.character), collapse=",")
      writeLines(paste0("Running MACS for group {", group_desc, "}"))
      results = macs2_coverage(sample_ranges=sample_ranges, control_ranges=control_ranges, params=params)
      results_df = dplyr::bind_cols(z[1,group_cols], dplyr::bind_rows(results[["islands"]], results[["qvalues"]])) %>%
        dplyr::relocate(dplyr::matches("qvalue_|island_"), .after=tidyselect::last_col())
      results_df
    })(.)) %>%
    dplyr::ungroup()

  islands_df = results_df %>%
    dplyr::filter_at(dplyr::vars(dplyr::starts_with("island_")), dplyr::all_vars(!is.na(.))) %>%
    dplyr::select(-dplyr::starts_with("qvalue_")) %>%
    dplyr::mutate(island_name=paste0("MACS3_", stringr::str_pad(1:dplyr::n(), 3, pad="0"))) %>%
    dplyr::relocate(island_name)
  qvalues_df = results_df %>%
    dplyr::filter_at(dplyr::vars(dplyr::starts_with("qvalue_")), dplyr::all_vars(!is.na(.))) %>%
    dplyr::select(-dplyr::starts_with("island_"))
  list(qvalues=qvalues_df, islands=islands_df)
}

#' @export
tlx_macs2 = function(tlx_df, effective_size, maxgap=NULL, qvalue=0.01, pileup=1, extsize=2000, slocal=200000, llocal=10000000, exclude_bait_region=F, exclude_repeats=F, exclude_offtargets=F, exttype, group, tmp_dir="tmp", tmp_name=NULL) {
  if(exclude_bait_region && !("tlx_is_bait_junction" %in% colnames(tlx_df))) {
    stop("tlx_is_bait_junction is not found in tlx data frame")
  }

  if(is.null(tmp_name)) tmp_name=basename(tempfile())

  validate_group(group)
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

  if(group=="all") {
    macs2_tlx_df$group = "all"
  }
  if(group=="sample") {
    macs2_tlx_df$group = paste(macs2_tlx_df$tlx_group, macs2_tlx_df$tlx_group_i)
  }
  if(group=="group") {
    macs2_tlx_df$group = macs2_tlx_df$tlx_group
  }
  dir.create(tmp_dir, recursive=T, showWarnings=F)

  macs_df.all = data.frame()
  for(gr in unique(macs2_tlx_df$group)) {
    tlx_df.gr = macs2_tlx_df %>% dplyr::filter(group==gr)

    f_input_bed = file.path(tmp_dir, paste0(tmp_name, "_", gr, "_input.bed"))
    f_control_bed = file.path(tmp_dir, paste0(tmp_name, "_", gr, "_control.bed"))

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
      macs_df = macs2(name=basename(f_input_bed), sample=f_input_bed, control=f_control_bed, maxgap=maxgap, effective_size=length(unique(tlx_df.gr$tlx_path))*effective_size, extsize=extsize, qvalue=qvalue, slocal=slocal, llocal=llocal, output_dir=dirname(f_input_bed))
    } else {
      log("Running MACS without control")
      macs_df = macs2(name=basename(f_input_bed), sample=f_input_bed, maxgap=maxgap, effective_size=effective_size, extsize=extsize, qvalue=qvalue, slocal=slocal, llocal=llocal, output_dir=dirname(f_input_bed))
    }

    if(group %in% c("sample", "group")) {
      macs_df$macs_group = tlx_df.gr$tlx_group[1]
    }
    if(group=="sample") {
      macs_df$tlx_group_i = tlx_df.gr$tlx_group_i[1]
    }
    if(group=="all") {
      macs_df$macs_group = "all"
    }
    macs_df = macs_df %>% dplyr::mutate(macs_group=gr) %>%
      dplyr::select(dplyr::matches("^(macs_group|macs_group_i)$"), dplyr::matches(".*"))

    macs_df.all = rbind(macs_df.all, macs_df)
  }

  macs_df.all = macs_df.all %>% dplyr::filter(macs_pileup>=pileup)

  macs_df.all
}
