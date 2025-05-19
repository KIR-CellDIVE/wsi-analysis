colors_types <- list(DP='lightgrey', Immuno=scales::muted('red'), Perivascular=scales::muted('blue'), Other='lightgrey')

plot_niche_expression <- function(niches, gene) {
    spots$metadata %>%
        ggplot2::ggplot(ggplot2::aes(y, -x)) +
            ggplot2::geom_point(data = . %>% subset(Niche %in% niches), size = .2, alpha = .4, shape = 21, fill='white') +
            ggplot2::geom_point(data = cbind(cells$metadata, as.matrix(spots$z)), ggplot2::aes_string(color = gene), shape = '.', alpha = .4) +
            ggplot2::scale_color_gradient2(low = 'white', high = scales::muted('red')) +
            ggplot2::labs(title = paste(niches, collapse = ' ')) +
            NULL
}


pal <- list(
    Other = 'lightgrey',
    Perivascular = scales::muted('red'),
    Vascular = scales::muted('blue'),
    Lymphoid = 'yellow',
    Lympho_vascular = 'green'
)

#' @export
make_spots <- function(cells, pow=1) {
    snn <- Matrix::Diagonal(n = nrow(cells$metadata))  ## Identity matrix

    spots <- list()
    spots$metadata <- tibble::tibble(
        dplyr::select(cells$metadata, SpotID = CellID, LibraryID, x, y, x_img, y_img),
        ncells = Matrix::colSums(snn),
        area = cells$metadata$area  ## each cell's own area (no pooling)
    )
    spots$intensity <- cells$intensity  ## directly copy intensity

    return(spots)
}


safe_make_spots <- function(cells, pow = 1, batch_size = 5000) {
    ## Split cells into batches
    n <- nrow(cells$metadata)
    split_indices <- split(seq_len(n), ceiling(seq_len(n) / batch_size))

    spots_list <- list()

    for (i in seq_along(split_indices)) {
        idx <- split_indices[[i]]

        cells_batch <- list(
            metadata = cells$metadata[idx, , drop = FALSE],
            intensity = cells$intensity[idx, , drop = FALSE]
        )

        ## Now call plain make_spots
        spots_batch <- make_spots(cells_batch, pow = pow)

        spots_list[[i]] <- spots_batch

        gc()
    }

    ## Combine results
    all_spots <- list(
        metadata = dplyr::bind_rows(lapply(spots_list, `[[`, "metadata")),
        intensity = do.call(rbind, lapply(spots_list, `[[`, "intensity"))
    )

    return(all_spots)
}

plotFeatures_split <- function(
    data_mat, dim_df, features, split_by, nrow = 1, qlo = 0.05, qhi = 1,
    order_by_expression = FALSE, pt_shape = 16, pt_size = 0.5,
    no_guide = FALSE, .xlim = c(NA, NA), .ylim = c(NA, NA), color_high = scales::muted("blue")
) {
    split(seq_len(length(split_by)), split_by) %>% purrr::map(function(idx) {
        plotFeatures(data_mat[, idx], dim_df[idx, ], features, nrow = 1, no_guide = TRUE)
    })
}


#' @export
do_norm <- function(obj) {
    obj[['intensity_norm']] <- sweep(obj$intensity, 1, obj$metadata$area, '/')
    return(obj)
}


#' @export
do_scale <- function(obj, z_thresh, within_batch=FALSE) {
    message("SCALE")
    if (within_batch) {
        .res <- split(seq_len(nrow(obj$metadata)), obj$metadata$LibraryID) %>% purrr::map(function(idx) {
            list(
                idx=idx,
                z=obj$intensity[idx, ] %>% log1p %>% scale %>% pmax(-z_thresh) %>% pmin(z_thresh)
            )
        })
        obj$z <- matrix(
            0,
            nrow=nrow(obj$intensity),
            ncol=ncol(obj$intensity),
            dimnames=dimnames(obj$intensity)
        )

        for (x in .res) {
            obj$z[x$idx, ] <- x$z
        }
    } else {
        obj$z <- obj$intensity %>% log1p %>% scale %>% pmax(-z_thresh) %>%
            pmin(z_thresh) %>% data.frame()
    }
    return(obj)
}


#' @export
do_pca <- function(obj) {
    message('PCA')
    svd_res <- svd(t(as.matrix(obj$z)))
    V <- with(svd_res, sweep(v, 2, d, '*'))
    V <- data.frame(V)
    rownames(V) <- rownames(obj$z)
    colnames(V) <- paste0('PC', 1:ncol(V))
    U <- data.frame(svd_res$u)
    colnames(U) <- paste0('PC', 1:ncol(U))
    rownames(U) <- colnames(obj$intensity)

    obj$V <- V
    obj$loadings <- U
    return(obj)
}


#' @export
do_umap <- function(obj, embedding, resname, ...) {
    message('UMAP')
    obj[[resname]] <- uwot::umap(obj[[embedding]], ret_extra = 'fgraph', ...)
    colnames(obj[[resname]]$embedding) <- paste0('UMAP', 1:2)
    return(obj)
}


#' @export
harmonize <- function(obj, var='LibraryID', theta=1) {
    obj$H <- harmony::RunHarmony(
        obj$V, obj$metadata, var,
        theta=theta,
        plot_convergence=TRUE,
        max_iter = 10,
        .options = harmony::harmony_options(epsilon.cluster = -Inf, epsilon.harmony = -Inf, max.iter.cluster = 10)
    ,)
    return(obj)
}


#' @export
do_louvain <- function(obj, embedding, resolutions, pow=1) {
    adjmat <- obj[[embedding]]$fgraph
    diag(adjmat) <- 1
    for (i in seq_len(pow-1)) {
        adjmat <- adjmat %*% adjmat
    }
    obj$Clusters <- RunModularityClustering(adjmat, 1, resolutions)
    # obj$Clusters <- paste0('C', RunModularityClustering(adjmat, 1, resolutions))
    return(obj)
}


#' @export
plot_heatmap <- function(markers, features, .scale=FALSE) {
    X <- markers %>%
        # subset(feature %in% c('ASMA', 'CD146', 'CD3', 'CD31', 'CD45', 'CD90', 'PDPN', 'CD68')) %>%
        # subset(!group %in% c('C22')) %>%
        dplyr::select(SCORE=auc, feature, group) %>%
        tidyr::spread(group, SCORE) %>%
        tibble::column_to_rownames('feature') %>%
        methods::as('matrix')

    if (.scale) X <- t(scale(t(X)))

    X <- X[features, ]
    ComplexHeatmap::Heatmap(X, column_names_rot = 45)
}


#' @export
do_markers <- function(obj) {
    obj$Markers <- apply(obj$Clusters, 2, function(y) {
        presto::wilcoxauc(t(obj$z), y)
    })

    names(obj$Markers) <- colnames(obj$Clusters)
    return(obj)
}


do_subcluster <- function(obj, embedding, resolution, clusters, clusters_zoom, pow=1) {
    idx <- which(clusters %in% clusters_zoom)
    adjmat <- obj[[embedding]]$fgraph
    adjmat <- adjmat[idx, idx]
    diag(adjmat) <- 1
    for (i in seq_len(pow - 1)) {
        adjmat <- adjmat %*% adjmat
    }
    # obj$Clusters <- RunModularityClustering(adjmat, 1, resolutions)
    new_clusters <- RunModularityClustering(adjmat, 1, resolution)
    new_clusters <- unlist(new_clusters)
    new_clusters <- paste0(clusters[idx], '_', new_clusters)

    res <- clusters
    res[idx] <- new_clusters
    return(res)
}

nice_plot_coloc <- function(spots, obj, niches, types, color='red', alpha=.1, show_libname=TRUE) {
    plt <- spots$metadata %>%
        ggplot2::ggplot(ggplot2::aes(y, -x)) +
            ggplot2::geom_point(
                data = . %>% dplyr::arrange(Niche_nice %in% niches),
                shape = '.', ggplot2::aes(color = Niche_nice %in% niches), fill = NA
            ) +
            ggplot2::theme_void() +
            ggplot2::scale_fill_manual(values = c(color)) +
            ggplot2::scale_color_manual(values = c('lightgrey', 'black')) +
            ggplot2::guides(
                color = 'none',
                fill = 'none',
                alpha = 'none'
            ) +
            ggplot2::facet_wrap(~LibraryID, scales='free') +
            ggplot2::theme(strip.text = ggplot2::element_text(size = 12)) +
            NULL

    if (!show_libname) {
        plt <- plt + ggplot2::theme(strip.text = ggplot2::element_blank())
    }

    plt +
        ggplot2::geom_point(
            data = subset(obj$metadata, Subtype %in% types),
            shape = 21,
            size = 1,
            ggplot2::aes(fill = Subtype),
            alpha = alpha
        ) +
        NULL

}

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

plot_biaxial <- function(dat, x, y, x0, y0, facet_by=NULL) {
    dat <- data.frame(dat)
    dat$density <- get_density(dat[, x], dat[, y], n = 100)

    if (!is.null(facet_by))
        dat$VAR <- facet_by

    plt <- ggplot2::ggplot(dat, ggplot2::aes_string(x, y, color = 'density')) +
        ggplot2::geom_point(shape = '.', alpha = .1, position = ggplot2::position_jitter(width = .02, height = .02)) +
        ggplot2::geom_hline(yintercept = y0, linetype = 2, color = 'red') +
        ggplot2::geom_vline(xintercept = x0, linetype = 2, color = 'red') +
        # scale_color_gradient2_tableau() +
        ggplot2::labs(x = x, y = y) +
        viridis::scale_color_viridis() +
        ggplot2::guides(color = 'none') +
        NULL

    if (!is.null(facet_by))
        plt <- plt + ggplot2::facet_wrap(~VAR, scales='free')

    return(plt)

}


coloc_heatmap <- function(stats) {
    stats %>%
        dplyr::mutate(estimate = exp(estimate)) %>%
        dplyr::select(Niche_nice, Subtype, estimate) %>%
        tidyr::spread(Niche_nice, estimate) %>%
        tibble::column_to_rownames('Subtype') %>%
        as.matrix() %>%
        ComplexHeatmap::Heatmap(column_names_rot = 45)
}


get_coloc_stats <- function (obj) {

    split(seq_len(nrow(obj$metadata)), obj$metadata$LibraryID) %>%
        purrr::imap(function(.idx, .libname) {
            X <- with(obj$metadata[.idx, ], table(Niche_nice, Subtype)) %>% data.frame()

            purrr::map(unique(X$Subtype), function(Subtype_level) {
                model <- lme4::glmer(
                    Subtype_test ~ 1 + (1|Niche_nice),
                    X %>% dplyr::mutate(Subtype_test = dplyr::case_when(
                        Subtype == Subtype_level ~ 1,
                        TRUE ~ 0
                    )),
                    stats::binomial,
                    weights = Freq
                )
                res <- data.frame(lme4::ranef(model)$Niche_nice)
                colnames(res) <- c('beta')
                res <- tibble::rownames_to_column(res, 'Niche_nice')
                res$Subtype <- Subtype_level

                ## compute p values
                .sim <- arm::sim(model, 100)
                res$sigma <- apply(data.frame(.sim@ranef$Niche_nice), 2, stats::sd)
                res$pval <- with(res, exp(pnorm(-beta/sigma, log.p = TRUE, lower.tail = TRUE))) ## 1-tailed
                return(res)
            }) %>%
                dplyr::bind_rows() %>%
                cbind(LibraryID = .libname) %>%
                dplyr::select(LibraryID, Subtype, Niche_nice, beta, sigma, pval) %>%
                dplyr::arrange(-beta/sigma)

        }) %>%
        dplyr::bind_rows() %>%
        dplyr::arrange(pval)

    # split(seq_len(nrow(obj$metadata)), obj$metadata$LibraryID) %>%
    # imap(function(.idx, .libname) {
    #     X <- with(obj$metadata[.idx, ], table(Niche_nice, Subtype)) %>% data.frame()
    #     expand.grid(Niche = unique(X$Niche_nice), Type = unique(X$Subtype)) %>%
    #         apply(1, function(vals) {
    #             glm(y ~ 1 + x, family = poisson, X %>% dplyr::mutate(y = Subtype ==
    #                 vals[["Type"]], x = Niche_nice == vals[["Niche"]]),
    #                 weights = Freq) %>% broom::tidy() %>% subset(term ==
    #                 "xTRUE") %>% dplyr::mutate(Niche_nice = vals[["Niche"]],
    #                 Subtype = vals[["Type"]]) %>% dplyr::select(-term) %>%
    #                 dplyr::select(Niche_nice, Subtype, everything())
    #         }) %>%
    #     bind_rows() %>%
    #     dplyr::mutate(LibraryID = .libname)
    # }) %>%
    #     bind_rows()
}

