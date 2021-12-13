

check_integer <- function(x){
  if (is.numeric(x)){
    return(mean(x%%1==0, na.rm =T)==1)
  } else{
    return(FALSE)
  }
}



peak_meta_data <- function(x){
  x <- dplyr::select(x,
                     any_of(
                       c(
                       "name", "X1",
                       "mzmed","mzmin", "mzmax", "mz",
                       "rtmed","rtmin","rtmax", "rt", "RT", "ESI", "ID",
                       "isotopes","adduct"
                       )
                     )
  )
  return(x)
}



intensity_data <- function(x){
  remove <- purrr::map_lgl(x, check_integer)
  x      <- x[,!remove]
  x <- dplyr::select(x,
                     -any_of(c("fold", "X1", "tstat",
                                         "anova", "pvalue", "isotopes",
                                         "adduct", "mzmed", "mzmin",
                                         "mzmax", "rtmed","rtmin", "rt", "mz", "RT", "ESI", "ID",
                                         "rtmax", "mean", "maxint")))
  colnames(x)<- gsub("/|[.]", "_", colnames(x))
  x <- x[!is.na(x$name),]
  new_cols <- x$name
  observation <- colnames(x)[2:ncol(x)]

  x <- dplyr::select(x, -name)
  x <- tibble::as_tibble(t(x))
  colnames(x) <- new_cols
  x <- dplyr::mutate(x, sample = observation)
  x <- dplyr::select(x, sample, everything())

  return(x)
}



#' Load data returned by XCMS getPeakTable function
#'
#' This function loads a csv as returned by getPeakTable (xcms). It wrangles
#' the data and returns two dataframes: a peak meta data data frame
#' describing mz, rt, isotopes, etc. and a peaks table with intensities and
#' observation ids.
#' The output is used in other ezmz functions.
#'
#' @param infile Path to the input file
#' @return A list of two dataframes
#' @export
load_data <- function(path){
  if (grepl("tsv|txt", path)){
    df   <- readr::read_tsv(path)
  } else {
    df   <- readr::read_csv(path)
  }
  if (!"name" %in% colnames(df)){
    df <- df %>%  mutate(name = paste("sample", 1:nrow(df), sep="")) %>%
      dplyr::select(name, everything())
  }
  meta <- peak_meta_data(df)
  peak <- intensity_data(df)
  data <- list("meta_data" = meta, "intensity_data" = peak)
  rm(df, meta, peak)
  return(data)
}
# x = df
# path = "C:/Users/johan/OneDrive - Aarhus Universitet/Universitetsnoter/PhD/GenomeDK/tox_age_prediction/xcms_preprocessed_uniform.csv"
# load_data(path)
# path = "wp2_wp3_n1_data.tsv"
# df   <- readr::read_tsv(path)
# meta <- peak_meta_data(df)
# peak <- intensity_data(df)
# df
