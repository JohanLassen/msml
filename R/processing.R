


#' Fourth root transform data
#'
#' @param ms
#'
#' @return ms
#' @export
#' @importFrom magrittr %>%
#' @import dplyr
#' @examples
#' fourth_root_transform(ms)
fourth_root_transform <- function(ms){

  ms$values     <- ms$values^0.25
  return(ms)

}


#' Remove features inflated with zeros (batch-wise)
#'
#' This might help avoiding batch effect and identification of batches through signals with binary (on/off) behaviour
#'
#' @param ms
#'
#' @return ms
#' @export
#' @importFrom magrittr %>%
#' @import dplyr
#' @examples
#' rm_feature_batch_inflated_zeros(ms)
rm_feature_batch_inflated_zeros <- function(ms){

  raw     <- ms$values
  rowinfo <- ms$rowinfo

  tmp1 <- rowinfo %>% select(BATCH) %>%
    bind_cols(raw) %>%
    pivot_longer(cols = starts_with("M")) %>%
    group_by(BATCH, name) %>%
    summarise(zeros = mean(value==0))

  bad_features <- tmp1[tmp1$zeros>0.2,] %>% pull(name) %>% table %>%  names() %>%  unique()

  raw <- raw %>% select(-all_of(bad_features))

  ms$values <- raw
  ms$rowinfo <- ms$rowinfo
}


#' Removes observations if their batches are below a given threshold
#'
#' Is usable if preceeded by batch normalization
#'
#' @param ms
#'
#' @return ms
#' @export
#' @importFrom magrittr %>%
#' @import dplyr
#' @examples
#' rm_sample_min_batch_occurence(ms)
rm_sample_min_batch_occurence <- function(ms, min_occurence = 5){

  tmp1 <- ms$rowinfo
  tmp2 <- ms$values0

  small_batches <- table(tmp1$BATCH)[table(tmp1$BATCH) < min_occurence] %>% names
  tmp1 <- tmp1 %>% filter(!BATCH %in% small_batches)
  tmp2 <- tmp2[tmp1$rowid,]
  tmp1 <- tmp1 %>% mutate(rowid = row_number())

  ms$rowinfo <- tmp1
  ms$values <- tmp2

}


#' Removes outliers by PCA
#'
#' Samples are removed if: abs(x - median(x)) > (1.5 * quantile(x,0.95)-quantile(x,0.05)), where x is the PCA scores.
#' Scores from the first 12 components are used.
#'
#' May be applied before row normalization but after removal of extreme features.
#'
#' @param ms
#'
#' @return ms
#' @export
#' @importFrom magrittr %>%
#' @import dplyr
#' @examples
#' rm_sample_pca_outliers(ms)
rm_sample_pca_outliers <- function(ms, plot = FALSE) {

  raw     <-  ms$values
  rowinfo <- ms$rowinfo

  tmp1 <- ms$rowinfo %>% mutate(rowid2 = row_number())
  tmp2 <- raw[tmp1$rowid,]

  r  <- prcomp(x = tmp2, retx = T, center=T, scale. = T, rank. = 12)

  bad_rows <- tibble(rowid2=apply(r$x, 2, function(x) {
    which(abs(x - median(x)) > (1.5 * quantile(x,0.95)-quantile(x,0.05)))
  }) %>%
    unlist() %>%
    as.vector()) %>%
    count(rowid2)

  tmp1 <- tmp1 %>%
    left_join(bad_rows) %>%
    mutate(n=ifelse(is.na(n), 0,n)) %>%
    mutate(label=case_when(
      n>0 ~ sample,
      n==0 ~ "")) %>%
    {.}

  pd <- r$x %>%
    as_tibble() %>%
    bind_cols(tmp1) %>%
    {.}

  pd <- pd %>%
    mutate(response = ifelse(n>0,"Outlier", "Not outlier")) %>%
    mutate(response = factor(response))

  if (plot) {
    plotlist <- list()

    for(i in 1:(ncol(r$x)/2)) {
      xvar <- names(pd)[2*i-1]
      yvar <- names(pd)[2*i]
      p1 <- ggplot(pd,aes(x=!!ensym(xvar), y=!!ensym(yvar),
                          fill=response, label=label))+
        geom_point(shape=21, color="#FFFFFFFF", size=3) +
        scale_fill_manual(values = c("#D0D0D0", "#D04040")) +
        theme(legend.position="none") +
        NULL

      plotlist[[length(plotlist)+1]] <- p1
      rm(p1)
    }

    cowplot::plot_grid(plotlist = plotlist,nrow=1)
    rm(plotlist)
  }

  bad_rows <- tmp1 %>% filter(n>0)

  if (nrow(bad_rows) > 0) {
    ms$values  <- raw[-bad_rows$rowid2,]
    ms$rowinfo <- rowinfo[-bad_rows$rowid2,]
  }

  ms$rowinfo <- ms$rowinfo %>%
    mutate(rowid = row_number())
}



#' Removes outliers by PCA
#'
#' Samples are removed if: abs(x - median(x)) > (1.5 * quantile(x,0.95)-quantile(x,0.05)), where x is the PCA scores.
#' Scores from the first 12 components are used.
#'
#' May be applied before row normalization but after removal of extreme features.
#'
#' @param ms
#'
#' @return ms
#' @export
#' @importFrom magrittr %>%
#' @import dplyr
#' @examples
#' rm_sample_pca_outliers(ms)
rm_feature_extreme_values <- function(ms) {

  raw     <- ms$values
  rowinfo <- ms$rowinfo

  tmp1 <- tibble(rowid = rowinfo$rowid) %>%
    bind_cols(as_tibble(raw))

  tmp1 <- tmp1 %>%
    pivot_longer(names_to = "compound", values_to = "value",  cols= c(-rowid))

  tmp2 <- tmp1 %>%
    group_by(compound) %>%
    summarise(n_bad = sum(value > median(value)+2*quantile(value,0.95))) %>%
    {.}

  bad_features <- tmp2 %>%
    ungroup() %>%
    filter(n_bad > 0) %>%
    select(compound) %>%
    distinct()

  ms$values <- raw %>% select(-any_of(bad_features$compound))
  ms$rowinfo <- rowinfo
}



#' Robust row normalization
#'
#' Only stable features used to normalize the rows (each sample)
#' Corrects matrix effects
#'
#' @param ms
#'
#' @return ms
#' @export
#' @importFrom magrittr %>%
#' @import dplyr
#' @examples
#' robust_row_normalize(ms)
robust_row_normalize <- function(ms) {

  target_info   <- ms$rowinfo
  target_values <- ms$values %>% as_tibble()

  tmp1 <- target_info
  tmp2 <- target_values[tmp1$rowid,]

  stable_features <- tmp1 %>%
    bind_cols(tmp2) %>%
    pivot_longer(starts_with("M")) %>%
    group_by(sample) %>%
    mutate(rank=rank(value)) %>%
    ungroup() %>%
    group_by(name) %>%
    summarise(median = median(rank),
              range = max(rank)-min(rank)) %>%
    ungroup() %>%
    slice_min(order_by = median, prop = 0.8) %>%
    slice_max(order_by = median, prop = 0.8) %>%
    slice_min(order_by = range, prop = 0.8)

  raw    <- target_values
  data.x <- raw
  tmp    <- rowSums(target_values %>% select(any_of(stable_features$name)))
  raw    <- max(raw)*raw / tmp

  ms$values  <- raw
  ms$rowinfo <- target_info
}



#' Select top n most variable features
#'
#' Effective for model screening to avoid heavy work loads. Often top 500 contains most of the signal.
#'
#' @param ms
#'
#' @return ms
#' @export
#' @importFrom magrittr %>%
#' @import dplyr
#' @examples
#' select_most_variable_features(ms)
select_most_variable_features <- function(ms, n=500) {
  raw <- ms$values
  rowinfo <- ms$rowinfo

  tmp1 <- rowinfo
  tmp2 <- raw[tmp1$rowid,]

  good_features <- tmp2 %>%
    bind_cols(tmp1) %>%
    pivot_longer(starts_with("M")) %>% # Convert to tidy format
    group_by(name) %>%
    summarise(mean = mean(value),
              variation = var(value)/mean(value)) %>%
    arrange(-variation) %>%
    slice(1:n) %>%
    pull(name)

  ms$rowinfo <- rowinfo
  ms$values  <- raw %>% select(any_of(good_features))
}



#' Z-standardization of batches
#'
#'
#' @param ms
#'
#' @return ms
#' @export
#' @importFrom magrittr %>%
#' @import dplyr
#' @examples
#' normalize_batch(ms)
normalize_batch <- function(ms, batch_column = "BATCH") {
  raw <- ms$values
  target_info <- ms$rowinfo

  tmp1 <- tibble(batch = target_info[[batch_column]], rowid = target_info$rowid) %>%
    bind_cols(as_tibble(raw)) %>%
    pivot_longer(names_to = "compound", values_to = "value",  cols= c(-batch, -rowid))

  tmp2 <- tmp1 %>%
    group_by(batch, compound) %>%
    summarise(value.batch_mean = mean(value),
              value.batch_sd   = sd(value))


  tmp1 <- full_join(x = tmp1, y = tmp2, all.x = TRUE) %>%
    mutate(value2 = (value- value.batch_mean)/(value.batch_sd)) %>% #
    mutate(value2 = ifelse((
      value.batch_mean==0 &
        value.batch_sd==0) |
        value == 0, 0, value2
    )) %>%
    select(batch, rowid, compound,value2) %>%
    pivot_wider(names_from = compound, values_from = value2)

  raw <- tmp1 %>%
    arrange(rowid) %>%
    select(-batch, -rowid)

  ms$rowinfo <- target_info
  ms$values  <- raw
}


