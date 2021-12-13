library(tidyverse)


#' PCA plot
#'
#' @param tmp1
#' @param tmp2
#' @param color_labels
#'
#' @return
#' @export
#'
#' @import tidyverse
#'
#' @examples
make_pca <- function(tmp1,tmp2,color_labels=c("AGE_YEARS", "BATCH")) {

  r  <- prcomp(x = tmp2, retx = TRUE, center = T, scale = T, rank. = 12)

  pd <- r$x %>%
    as_tibble() %>%
    bind_cols(tmp1 %>% select(all_of(color_labels))) %>%
    {.}

  plotlist <- list()
  titles <- list()

  for(i in 1:(ncol(r$x)/2)) {
    mini_plotlist <- list()
    #titles <- list()
    for (fill in color_labels) {
      #i <- 1
      xvar <- names(pd)[2*i-1]
      yvar <- names(pd)[2*i]
      p1 <-
        ggplot(pd, aes_string(x=xvar, y=yvar, fill=fill))+
        geom_point(shape=21, color="#FFFFFFFF", size=2, show.legend = F) +
        labs(fill = fill)+
        theme() +
        NULL
      if (i == 1) {
        title <- ggdraw() + draw_label(fill)
        titles[[length(titles)+1]] <- title
      }

      mini_plotlist[[length(mini_plotlist)+1]] <- p1
      rm(p1)
    }

    aggrgte <- plot_grid(plotlist = mini_plotlist, ncol = 1)
    plotlist[[length(plotlist)+1]] <- aggrgte
  }

  p1 <- cowplot::plot_grid(plotlist = plotlist, nrow=1)
  title <- ggdraw() + draw_label(paste(ncol(tmp2), "features"))
  p2 <- plot_grid(title, p1, ncol = 1, rel_heights=c(0.1, 1))

  titles <- plot_grid(plotlist = titles, ncol=1)
  final <- cowplot::plot_grid(titles, p2, nrow=1, rel_widths = c(2, 10))


  return(final)
}





#' Title
#'
#' @param tmp1
#' @param tmp2
#' @param color_label
#'
#' @return
#' @export
#' @import tidyverse
#' @import umap
#' @examples
make_umap <- function(tmp1, tmp2, color_label) {
  r <- umap::umap(tmp2)
  o <- tibble("umap1"=r$layout[,1],
              "umap2"=r$layout[,2])
  pd <- o %>%
    bind_cols(tmp1) %>%
    {.}

  p <- ggplot(pd,aes(x=umap1, y=umap2,
                     fill=.data[[color_label]]))+
    geom_point(shape=21, color="#FFFFFFFF", size=3, show.legend = F) +
    labs(fill = color_label)+
    theme() +
    NULL
  return(p)
}


#' probability histogram of caret model
#'
#' @param model
#'
#' @return
#' @export
#' @import tidyverse
#' @examples
probs_hist <- function(model){
  classes <- unique(model_cv_r$obs)
  ggplot(data = model$pred %>% as_tibble()) +
  geom_histogram(
    breaks = seq(0,1,length.out = 100),
    color = "gray20",
    alpha = 0.4,
    aes_string(x=as.character(classes[2]), fill="obs"),
    position = "identity"
  ) +
  theme_minimal()+
  scale_fill_manual(values = c("#d53e4f", "#3288bd"))+
  labs(title = best_model)
}


#' ROC of caret model
#'
#' @param model
#'
#' @return
#' @export
#'
#' @examples
plot_roc <- function(model){
  # ROC CURVE - ALL MODELS
  classes <- unique(model_cv_r$obs)
  g <- ggplot(model$pred,
              aes_string(m=classes[1],
                  d=factor("obs", levels = c(as.character(classes[1]), as.character(classes[2]))),
                  color = "model")
  ) +
    geom_roc(n.cuts=0) +
    coord_equal() +
    style_roc()

  g +
    annotate("text", label=paste("ROC", " =", round((calc_auc(g))$AUC, 4)))+
    scale_color_brewer(palette = "Spectral")
}

#' Importance from caret model
#'
#' @param model
#'
#' @return
#' @export
#' @import tidyverse
#' @import caret
#' @examples
plot_importance <- function(model){

  # FEATURE IMPORTANCE - BEST MODEL
  importance <- varImp(model) %>%
    {
      tibble(Feature = rownames(.$importance),
             Importance = .$importance %>% unlist %>% round(2)) %>%
        arrange(desc(Importance))
    }


  importance[1:10,] %>%
    ggplot(aes(y=factor(Feature, levels = rev(importance$Feature)), x=Importance)) +
    geom_col()+
    labs(y = "Features")+
    theme_minimal()

  important <- importance %>% filter(Importance>5) %>% arrange(-Importance) %>% slice(1:25)



  outcome <- model %>% select(obs, rowIndex) %>% distinct() %>% arrange(rowIndex)

  outcome %>%
    as_tibble() %>%
    bind_cols(model$trainingData) %>%
    select(obs, all_of(important$Feature)) %>%
    filter(!is.na(obs)) %>%
    pivot_longer(names_to = "Feature", cols=-obs) %>%
    ggplot() +
    geom_boxplot(aes(x=obs, y=value), fill = "#6DA398", show.legend = F)+
    facet_wrap(~Feature, scales = "free_y")+
    theme_bw()
}


#' Scatter plot of probabilities
#'
#' Used when classifying over a continuous variable
#'
#' @param model
#' @param target_info
#' @param var_name
#'
#' @return
#' @export
#'
#' @examples
prob_vs_continuous <- function(model, target_info, var_name){
  pdata <-
    model$pred %>%
    arrange(rowIndex) %>%
    mutate(rowid = rowIndex) %>%
    inner_join(target_info %>% drop_na() %>% mutate(rowid = row_number()), by = "rowid")

  classes <- unique(model_cv_r$obs)

  ggplot(pdata, aes_string(y=classes[1], x = var_name, fill = "obs")) +
    geom_point(shape = 21, color = "white") +
    scale_fill_manual(values = c("#d53e4f", "#3288bd"))+
    geom_segment(aes(x = min(age), xend = max(age), y=0.5, yend = 0.5), color = "black", linetype = "dashed")
}

