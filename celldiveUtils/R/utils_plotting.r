## It's a nice theme!
ggplot2::theme_set(ggthemes::theme_tufte(base_size = 14, base_family = "Helvetica"))

## Tufte is nice but removes boundaries from plots. Add it back here.
ggplot2::theme_update(panel.background = ggplot2::element_rect(fill = NA, color = "black"))

ggplot2::theme_update(
    axis.title = ggplot2::element_text(size = 12),
    axis.text = ggplot2::element_text(size = 10)
)

opts <- options()  # save old options
# options(ggplot2.continuous.colour="viridis")
# options(ggplot2.continuous.fill = "viridis")
# scale_colour_discrete <- function(...) {
#   scale_colour_manual(..., values = viridis::viridis(2))
# }





################# FOREST PLOTS ###########################
forest_multi <- function(data_df, sigma_max = Inf, reorder_levels=TRUE, grid.switch = NULL) {
    if (reorder_levels) {
        cluster_levels <- data.table::data.table(data_df)[, .(Z = mean(beta)), by = Cluster][order(Z), as.character(Cluster)]
        data_df <- data_df %>%
            dplyr::mutate(Cluster = factor(Cluster, cluster_levels))
    }

    data_df %>%
        ## shrink sigma if it's not plotable (make a note?)
        dplyr::mutate(sigma = pmin(sigma, sigma_max)) %>%
        ggplot2::ggplot(ggplot2::aes(Tissue, beta, color = Tissue)) +
            ggplot2::geom_point() +
            ggplot2::geom_errorbar(ggplot2::aes(ymin = beta - 1.96 * sigma, ymax = beta + 1.96 * sigma), width = 0) +
            ggplot2::coord_flip() +
#             theme_test(base_size = 14) +
#             scale_color_tableau() +
            ggplot2::facet_grid(Cluster~., space = 'free', scales = 'free', switch = grid.switch) +
            ggplot2::geom_hline(yintercept = c(0), linetype = 2) +
#             geom_hline(yintercept = c(-1, 0, 1), linetype = 2) +
            ggplot2::labs(x = '', y = 'Log2 fold change') +
            ggplot2::geom_hline(yintercept = c(0), linetype = 2) +
            ggplot2::theme(
                axis.text.y = ggplot2::element_blank(),
                axis.ticks.y = ggplot2::element_blank(),
                legend.position = 'bottom',
                strip.text.y = ggplot2::element_text(angle = 0, hjust = 0),
                panel.background = ggplot2::element_blank()
            ) +
#             guides(color = FALSE) +
            NULL

}

forest_uni <- function(data_df, sigma_max = Inf, fdr_max=0.05, beta_min=-Inf) {
    ## order clusters by effect size and direction
    cluster_levels <- data.table::data.table(data_df)[, .(Z = mean(beta)), by = Cluster][order(Z), as.character(Cluster)]

    data_df %>%
        dplyr::mutate(Cluster = factor(Cluster, cluster_levels)) %>%
        ## shrink sigma if it's not plotable (make a note?)
        dplyr::mutate(sigma = pmin(sigma, sigma_max)) %>%
        ggplot2::ggplot(ggplot2::aes(Cluster, beta, color = fdr < fdr_max & beta > beta_min)) +
            ggplot2::geom_point() +
            ggplot2::geom_errorbar(ggplot2::aes(ymin = beta - 1.96 * sigma, ymax = beta + 1.96 * sigma), width = 0) +
            ggplot2::coord_flip() +
#             theme_test(base_size = 14) +
            ggplot2::scale_color_manual(values = c('black', 'red')) +
            ggplot2::geom_hline(yintercept = c(0), linetype = 2) +
#             geom_hline(yintercept = c(-1, 0, 1), linetype = 2) +
            ggplot2::labs(x = '', y = 'Log2 fold change') +
            ggplot2::guides(color = FALSE) +
            NULL

}


################# FOREST PLOTS ###########################

remove_strip_box <- function(plt) {
    pg <- ggplot2::ggplotGrob(plt)

    for(i in which(grepl("^strip", pg$layout$name))){
      pg$grobs[[i]]$layout$clip <- "off"
    }

    return(pg)
#     return(wrap_elements(full = pg))

}
