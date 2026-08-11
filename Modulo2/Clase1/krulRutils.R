# krulRutils: A collection of utility functions for
# data manipulation and visualization

#' Construct simple outline polygons for clustered 2D data
#'
#' Builds per-group outline polygons over a 2D projection using a
#' convex hull when groups have 3+ points, a padded rectangle for 2 points,
#' and a small circle for singletons. Intended for layering with ggplot2
#' to visualize cluster boundaries across grouping variables (e.g., k, cluster).
#'
#' @param tbl A tibble or data.frame containing at least `xvar`, `yvar`, and
#'   the grouping columns passed in `...`.
#' @param xvar Bare column name for the x coordinate.
#' @param yvar Bare column name for the y coordinate.
#' @param ... Bare column name(s) used to group rows into clusters
#'   (e.g., `k, cluster`).
#'
#' @return A tibble with columns `xvar`, `yvar`, and the grouping columns
#'   passed in `...`, where each row is a vertex of an outline polygon for a
#'   group. Suitable for use with `geom_path()` and
#'   `group = interaction(...)`.
#'
#' @examples
#' library(dplyr)
#' iris_tbl <- as_tibble(iris) %>% mutate(.row_id = row_number())
#' # Example grouping by Species only
#' poly <- cluster_polygons(iris_tbl, Sepal.Length, Sepal.Width, Species)
#'
#' @export
cluster_polygons <- function(tbl, xvar, yvar, ...) {
  id_quos <- rlang::enquos(...)

  xrange <- tbl %>%
    dplyr::pull({{ xvar }}) %>%
    range(na.rm = TRUE)

  yrange <- tbl %>%
    dplyr::pull({{ yvar }}) %>%
    range(na.rm = TRUE)

  xpad <- diff(xrange) * 0.03
  ypad <- diff(yrange) * 0.03

  if (length(id_quos) == 0) {
    grouped <- dplyr::group_by(tbl, .dummy = 1L, .drop = TRUE)
  } else {
    grouped <- dplyr::group_by(tbl, !!!id_quos, .drop = TRUE)
  }

  groups <- dplyr::group_split(grouped)
  keys <- dplyr::group_keys(grouped)

  polys <- lapply(seq_along(groups), function(i) {
    d <- groups[[i]]
    xs <- dplyr::pull(d, {{ xvar }})
    ys <- dplyr::pull(d, {{ yvar }})
    n <- length(xs)
    if (n >= 3) {
      idx <- grDevices::chull(xs, ys)
      xy <- cbind(xs[idx], ys[idx])
      xy <- rbind(xy, xy[1, , drop = FALSE])
    } else if (n == 2) {
      xmin <- min(xs) - xpad
      xmax <- max(xs) + xpad
      ymin <- min(ys) - ypad
      ymax <- max(ys) + ypad
      xy <- rbind(
        c(xmin, ymin),
        c(xmin, ymax),
        c(xmax, ymax),
        c(xmax, ymin),
        c(xmin, ymin)
      )
    } else if (n == 1) {
      cx <- xs[1]
      cy <- ys[1]
      ang <- seq(0, 2 * pi, length.out = 40)
      r <- max(xpad, ypad, 1e-3)
      xy <- cbind(cx + r * cos(ang), cy + r * sin(ang))
    } else {
      return(NULL)
    }

    # Assign column names before converting to tibble to avoid name-repair warning
    colnames(xy) <- c(
      rlang::as_name(rlang::enquo(xvar)),
      rlang::as_name(rlang::enquo(yvar))
    )
    poly <- tibble::as_tibble(xy)

    meta <- keys[i, , drop = FALSE]
    meta_rep <- meta[rep(1, nrow(poly)), , drop = FALSE]
    dplyr::bind_cols(poly, meta_rep)
  })

  polys <- polys[!vapply(polys, is.null, logical(1))]
  if (!length(polys)) return(NULL)
  dplyr::bind_rows(polys)
}


#' Compute a dendrogram cut height that yields k clusters
#'
#' Given an `hclust` object, returns a height value `h` such that
#' cutting the dendrogram at that height with `cutree(hc, h = h)` produces
#' exactly `k` clusters (for 1 <= k <= n). The value is chosen midway between
#' the merge heights that transition from `k+1` to `k` clusters.
#'
#' @param hc An object of class `hclust`.
#' @param k Integer, desired number of clusters (1 to n).
#'
#' @return A numeric height at which to cut the dendrogram.
#' @examples
#' # hc <- hclust(dist(iris[, 1:4]))
#' # h  <- cut_dendrogram_height_for_k(hc, k = 3)
#' @export
cut_dendrogram_height_for_k <- function(hc, k) {
  h <- hc$height
  n <- length(h) + 1
  m <- n - k
  if (m <= 0) return(min(h, na.rm = TRUE) - 1e-6)     # effectively n clusters
  if (m >= length(h)) return(max(h, na.rm = TRUE) + 1e-6) # effectively 1 cluster
  (h[m] + h[m + 1]) / 2
}

#' Choose DBSCAN parameters targeting a clustered solution
#'
#' Tries a grid of `eps` values (with fixed `minPts`) and returns the
#' `dbscan::dbscan()` result whose number of non-noise clusters is closest to a
#' target and with the lowest noise proportion. If an exact match to the target
#' number of clusters is found with acceptable noise proportion, it is returned
#' immediately. The selected result carries an attribute `"eps"` with the
#' chosen epsilon.
#'
#' This helper is useful for quick, pedagogical exploration and small/medium
#' datasets. It performs a simple grid search and does not guarantee a globally
#' optimal solution.
#'
#' @param data_matrix A numeric matrix of observations (rows) by features (cols).
#' @param min_pts Integer; DBSCAN `minPts` parameter (default 5).
#' @param eps_grid Numeric vector of candidate `eps` values to try.
#'   Defaults to `seq(0.15, 0.8, by = 0.01)`.
#' @param max_noise_prop Maximum acceptable noise fraction for an exact target
#'   match (default `0.1`).
#' @param target_k Desired number of clusters to aim for (default `3`).
#'
#' @return An object returned by `dbscan::dbscan()` for the chosen `eps`, with
#'   an additional attribute `"eps"` storing the selected epsilon.
#'
#' @examples
#' \dontrun{
#'   set.seed(123)
#'   X <- scale(as.matrix(iris[, 1:2]))
#'   res <- choose_dbscan(X, min_pts = 5, target_k = 3)
#'   attr(res, "eps")
#'   table(res$cluster)
#' }
#'
#' @export
choose_dbscan <- function(
  data_matrix,
  min_pts = 5,
  eps_grid = seq(0.15, 0.8, by = 0.01),
  max_noise_prop = 0.1,
  target_k = 3
) {
  # fast path: try to find exact target_k with acceptable noise
  for (eps in eps_grid) {
    res <- dbscan::dbscan(data_matrix, eps = eps, minPts = min_pts)
    k <- length(unique(res$cluster[res$cluster > 0]))
    noise_prop <- mean(res$cluster == 0)
    if (k == target_k && noise_prop <= max_noise_prop) {
      attr(res, "eps") <- eps
      return(res)
    }
  }

  # fallback: rank by closeness to target_k, then by noise proportion
  cand <- lapply(eps_grid, function(eps) {
    res <- dbscan::dbscan(data_matrix, eps = eps, minPts = min_pts)
    k <- length(unique(res$cluster[res$cluster > 0]))
    list(eps = eps, res = res, k = k, noise = mean(res$cluster == 0))
  })
  ord <- order(
    vapply(cand, function(x) abs(x$k - target_k), numeric(1)),
    vapply(cand, function(x) x$noise, numeric(1))
  )
  best <- cand[[ord[1]]]
  attr(best$res, "eps") <- best$eps
  best$res
}

#' Reassign DBSCAN noise to the nearest cluster centroid
#'
#' Given a matrix of observations and a vector of cluster labels where
#' label 0 denotes noise (as returned by `dbscan::dbscan()`), this helper
#' reassigns each noise point to the nearest cluster based on squared
#' Euclidean distance to per-cluster centroids. If there is no noise, the
#' original labels are returned unchanged.
#'
#' This is a pragmatic post-processing step for visualization/teaching to
#' avoid orphaned points in faceted plots. It is not part of the DBSCAN
#' algorithm and should be used with care in analytical pipelines.
#'
#' @param data_matrix Numeric matrix with rows as observations and columns as
#'   features.
#' @param cluster_labels Integer vector of cluster assignments where 0 marks
#'   noise and positive integers mark clusters. Length must match
#'   `nrow(data_matrix)`.
#'
#' @return Integer vector of the same length as `cluster_labels` with noise
#'   points reassigned to the nearest centroid of existing clusters.
#'
#' @examples
#' \dontrun{
#'   set.seed(123)
#'   X <- scale(as.matrix(iris[, 1:2]))
#'   db <- choose_dbscan(X, min_pts = 5)
#'   labels_no_noise <- assign_noise_to_nearest(X, db$cluster)
#' }
#'
#' @export
assign_noise_to_nearest <- function(data_matrix, cluster_labels) {
  output_labels <- cluster_labels
  noise_idx <- which(output_labels == 0)
  if (length(noise_idx) == 0) return(output_labels)

  keep <- output_labels > 0
  if (!any(keep)) return(output_labels)

  centers <- t(sapply(
    sort(unique(output_labels[keep])),
    function(k) colMeans(data_matrix[output_labels == k, , drop = FALSE])
  ))
  rownames(centers) <- as.character(sort(unique(output_labels[keep])))

  for (i in noise_idx) {
    d <- apply(centers, 1, function(cn) sum((data_matrix[i, ] - cn)^2))
    output_labels[i] <- as.integer(names(which.min(d)))
  }
  output_labels
}

#' Convert codes in a tibble column to a labeled factor using a lookup table
#'
#' @param data_tbl A tibble containing the data with codes to convert.
#' @param code_col Unquoted name of the column in `data_tbl` with the codes.
#' @param lookup_tbl A tibble with code-label pairs for lookup.
#' @param lookup_code_col Unquoted name of the code column in `lookup_tbl`.
#' @param lookup_label_col Unquoted name of the label column in `lookup_tbl`.
#' @param factor_levels Optional character vector specifying the levels order
#' for the output factor. Defaults to all labels in `lookup_tbl`, in order.
#' @param new_factor_col_name Optional name for the new factor column.
#'
#' @return A tibble with a new factor vector with labels corresponding to codes,
#' same length/order as `data_tbl`.
#' Codes not found in the lookup table result in NA and trigger a warning.
#' @export
convert_codes_to_factor <- function(
    data_tbl,
    code_col, # Column in data_tbl to match on
    lookup_tbl, # The lookup table
    lookup_code_col, # Column in lookup_tbl that holds the codes
    lookup_label_col, # Column in lookup_tbl that holds the labels
    factor_levels = NULL, # Optional: explicit order of factor levels
    new_factor_col_name = NULL) {
  # Capture unquoted column names for tidy evaluation
  code_col_sym <- rlang::ensym(code_col)
  lookup_code_col_sym <- rlang::ensym(lookup_code_col)
  lookup_label_col_sym <- rlang::ensym(lookup_label_col)

  new_factor_col_name_input <- rlang::enquo(new_factor_col_name)
  if (rlang::quo_is_null(new_factor_col_name_input)) {
    new_factor_col_name_sym <- rlang::sym(
      paste0(rlang::as_label(code_col_sym), "_factor")
    )
  } else {
    new_factor_col_name_sym <- rlang::sym(
      rlang::as_label(new_factor_col_name_input)
    )
  }
  # Store original column names to ensure they are kept
  original_col_names <- names(data_tbl)

  # Perform the join to get the labels
  join_by_args <- stats::setNames(
    rlang::as_label(lookup_code_col_sym),
    rlang::as_label(code_col_sym)
  )

  joined_tbl <- data_tbl %>%
    left_join(
      lookup_tbl,
      by = join_by_args
    )

  # Check for missing codes and issue a warning
  missing_codes <- joined_tbl %>%
    filter(is.na(!!lookup_label_col_sym)) %>%
    distinct(!!code_col_sym) %>%
    pull(!!code_col_sym)

  if (length(missing_codes) > 0) {
    warning(
      "The following codes were not found in the lookup table ",
      "and will be set to NA in the factor: ",
      paste(unique(missing_codes), collapse = ", ")
    )
  }

  # Determine factor levels if not provided
  if (is.null(factor_levels)) {
    factor_levels <- lookup_tbl %>%
      pull(!!lookup_label_col_sym) %>%
      unique()
  }

  # Add the labeled factor column
  final_tbl <- joined_tbl %>%
    mutate(
      !!new_factor_col_name_sym := factor(
        !!lookup_label_col_sym,
        levels = factor_levels
      )
    )

  # Now, select the desired columns: original columns + the new factor column.
  # We need to explicitly pick the *original* columns that were in data_tbl,
  # and then add the new factor column.
  # The label column from the lookup_tbl that was joined is temporary.
  final_tbl <- final_tbl %>%
    select(
      all_of(original_col_names), # Selects all columns from the original data
      !!new_factor_col_name_sym # Add the newly created factor column
    )

  return(final_tbl)
}





#' Summarise numeric variables in a tidy format
#'
#' Computes min, quartiles, median, mean, and max for all numeric columns,
#' returning a tibble with statistics as rows and variables as columns.
#'
#' @param df A data frame or tibble.
#' @param na.rm Logical, whether to remove NA values. Default TRUE.
#'
#' @return A tibble with rows = statistics and columns = variables.
#' @export
summarise_numeric_tidy <- function(df, na.rm = TRUE) {
  # Check input
  if (!is.data.frame(df)) {
    stop("Input must be a data.frame or tibble.")
  }

  num_cols <- df %>%
    select(where(is.numeric)) %>%
    names()

  if (length(num_cols) == 0) {
    stop("No numeric columns found in input.")
  }

  # Define summary functions with safe wrappers

  # Function for finding the minimal value
  safe_min <- function(x) {
    if (all(is.na(x))) {
      NA_real_
    } else {
      min(x, na.rm = na.rm)
    }
  }

  # Function for finding the first quartile
  safe_q1 <- function(x) {
    if (all(is.na(x))) {
      NA_real_
    } else {
      quantile(x, 0.25, na.rm = na.rm)
    }
  }

  # Function for finding the median value
  safe_median <- function(x) {
    if (all(is.na(x))) {
      NA_real_
    } else {
      median(x, na.rm = na.rm)
    }
  }

  # Function for finding the third quartile
  safe_q3 <- function(x) {
    if (all(is.na(x))) {
      NA_real_
    } else {
      quantile(x, 0.75, na.rm = na.rm)
    }
  }

  # Function for finding the maximal value
  safe_max <- function(x) {
    if (all(is.na(x))) {
      NA_real_
    } else {
      max(x, na.rm = na.rm)
    }
  }

  # Function for finding the mean value
  safe_mean <- function(x) {
    if (all(is.na(x))) {
      NA_real_
    } else {
      mean(x, na.rm = na.rm)
    }
  }

  # Function for finding the mean value
  safe_var <- function(x) {
    if (all(is.na(x))) {
      NA_real_
    } else {
      var(x, na.rm = na.rm)
    }
  }


  # Compute summaries wide, then reshape long and wide as requested
  df %>%
    summarise(
      across(
        all_of(num_cols),
        list(
          min = safe_min,
          q1 = safe_q1,
          median = safe_median,
          q3 = safe_q3,
          max = safe_max,
          mean = safe_mean,
          var = safe_var
        ),
        .names = "{.col}_{.fn}"
      )
    ) %>%
    pivot_longer(
      everything(),
      names_to = c("variable", "statistic"),
      names_sep = "_",
      values_to = "value"
    ) %>%
    pivot_wider(
      names_from = variable,
      values_from = value
    ) %>%
    mutate(
      statistic = factor(
        statistic,
        levels = c("min", "q1", "median", "q3", "max", "mean", "var")
      )
    ) %>%
    arrange(statistic)
}


#' Modifies the color of the grid lines in ggplot2 plots
#' and sets the title size and face.
#'
#' @return A ggplot2 theme with modified grid line colors.
#' @export
theme_krul <- function() {
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    panel.grid.major = element_line(color = "gray80"),
    panel.grid.minor = element_line(color = "gray80")
  )
}

#' Applies the 'krul' theme to a GGally matrix
#'
#' @param ggpairs_obj A GGally ggpairs object.
#' @return The modified ggpairs object with the 'krul' theme applied.
#' @export
theme_krul_ggpairs <- function(ggpairs_obj) {
  n <- length(ggpairs_obj$xAxisLabels)

  for (i in 1:n) {
    for (j in 1:i) {
      ggpairs_obj[i, j] <- ggpairs_obj[i, j] +
        theme_krul()
    }
  }

  ggpairs_obj
}





#' Creates a boxplot legend explaining the components of a boxplot
#'
#' @param family Font family for the text in the legend.
#' @return A ggplot object with the boxplot legend.
#' @export
ggplot_box_legend <- function(family = "serif") {
  # Create data to use in the boxplot legend:
  set.seed(100)

  sample_df <- data.frame(
    parameter = "test",
    values = sample(500)
  )

  # Extend the top whisker a bit:
  sample_df$values[1:100] <- 701:800
  # Make sure there's only 1 lower outlier:
  sample_df$values[1] <- -350

  # Function to calculate important values:
  ggplot2_boxplot <- function(x) {
    quartiles <- as.numeric(quantile(x,
      probs = c(0.25, 0.5, 0.75)
    ))

    names(quartiles) <- c(
      "25th percentile",
      "50th percentile\n(median)",
      "75th percentile"
    )

    iqr <- diff(quartiles[c(1, 3)])

    upper_whisker <- max(x[x < (quartiles[3] + 1.5 * iqr)])
    lower_whisker <- min(x[x > (quartiles[1] - 1.5 * iqr)])

    upper_dots <- x[x > (quartiles[3] + 1.5 * iqr)]
    lower_dots <- x[x < (quartiles[1] - 1.5 * iqr)]

    list(
      "quartiles" = quartiles,
      "25th percentile" = as.numeric(quartiles[1]),
      "50th percentile\n(median)" = as.numeric(quartiles[2]),
      "75th percentile" = as.numeric(quartiles[3]),
      "IQR" = iqr,
      "upper_whisker" = upper_whisker,
      "lower_whisker" = lower_whisker,
      "upper_dots" = upper_dots,
      "lower_dots" = lower_dots
    )
  }

  # Get those values:
  ggplot_output <- ggplot2_boxplot(sample_df$values)

  # Lots of text in the legend, make it smaller and consistent font:
  update_geom_defaults(
    "text",
    list(
      size = 3,
      hjust = 0,
      family = family
    )
  )
  # Labels don't inherit text:
  update_geom_defaults(
    "label",
    list(
      size = 3,
      hjust = 0,
      family = family
    )
  )

  # Create the legend:
  # The main elements of the plot (the boxplot, error bars, and count)
  # are the easy part.
  # The text describing each of those takes a lot of fiddling to
  # get the location and style just right:
  explain_plot <- ggplot() +
    stat_boxplot(
      data = sample_df,
      aes(x = parameter, y = values),
      geom = "errorbar", width = 0.3
    ) +
    geom_boxplot(
      data = sample_df,
      aes(x = parameter, y = values),
      width = 0.3, fill = "lightgrey"
    ) +
    geom_text(aes(x = 1, y = 950, label = "500"), hjust = 0.5) +
    geom_text(
      aes(
        x = 1.17, y = 950,
        label = "Number of values"
      ),
      fontface = "bold", vjust = 0.4
    ) +
    theme_minimal(base_size = 5, base_family = family) +
    geom_segment(aes(
      x = 2.3, xend = 2.3,
      y = ggplot_output[["25th percentile"]],
      yend = ggplot_output[["75th percentile"]]
    )) +
    geom_segment(aes(
      x = 1.2, xend = 2.3,
      y = ggplot_output[["25th percentile"]],
      yend = ggplot_output[["25th percentile"]]
    )) +
    geom_segment(aes(
      x = 1.2, xend = 2.3,
      y = ggplot_output[["75th percentile"]],
      yend = ggplot_output[["75th percentile"]]
    )) +
    geom_text(aes(x = 2.4, y = ggplot_output[["50th percentile\n(median)"]]),
      label = "Interquartile\nrange", fontface = "bold",
      vjust = 0.4
    ) +
    geom_text(
      aes(
        x = c(1.17, 1.17),
        y = c(
          ggplot_output[["upper_whisker"]],
          ggplot_output[["lower_whisker"]]
        ),
        label = c(
          paste0(
            "Largest value within 1.5 times\n",
            "interquartile range above\n",
            "75th percentile"
          ),
          paste0(
            "Smallest value within 1.5 times\n",
            "interquartile range below\n",
            "25th percentile"
          )
        ),
      ),
      fontface = "bold", vjust = 0.9
    ) +
    geom_text(
      aes(
        x = c(1.17),
        y = ggplot_output[["lower_dots"]],
        label = "Outside value"
      ),
      vjust = 0.5, fontface = "bold"
    ) +
    geom_text(
      aes(
        x = c(1.17),
        y = ggplot_output[["lower_dots"]],
        label = paste0(
          "-Value is >1.5 times and \n",
          "<3 times the interquartile range\n",
          "beyond either end of the box"
        )
      ),
      vjust = 1.3
    ) +
    geom_label(
      aes(
        x = 1.17, y = ggplot_output[["quartiles"]],
        label = names(ggplot_output[["quartiles"]])
      ),
      vjust = c(0.4, 0.85, 0.4),
      fill = "white", label.size = 0
    ) +
    ylab("") +
    xlab("") +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      aspect.ratio = 4 / 3,
      plot.title = element_text(hjust = 0.5, size = 10)
    ) +
    coord_cartesian(xlim = c(1.4, 3.1), ylim = c(-600, 900)) +
    labs(title = "EXPLANATION")

  return(explain_plot)
}

#' Color palette from the 'Carnelian' color theme
#'
#' @export
c_palette <- c(
  `C red`          = "#B31B1B",
  `C orange`       = "#67411B",
  `C yellow`       = "#4F4F1B",
  `C chartreuse`   = "#3A591B",
  `C green`        = "#1B681B",
  `C spring green` = "#1B623E",
  `C cyan`         = "#1B5C5C",
  `C azure`        = "#1B548C",
  `C blue`         = "#3B3BBB",
  `C violet`       = "#6D21BB",
  `C magenta`      = "#8A1B8A",
  `C rose`         = "#9C1B5B",
  `C grey`         = "#494949"
)

#' Custom color palette for C colors
#'
#' Returns a named vector of hexadecimal color values.
#' @param ... Optional names of specific colors to return
#' @export
c_pal <- function(...) {
  cols <- c(...)
  if (length(cols) == 0) {
    return(c_palette)
  }
  c_palette[cols]
}

#' Integration with ggplot2 'scale_fill_manual()' function
#'
#' @param ... Optional colors to be chosen from the palette
#' @export
c_scale_fill <- function(...) {
  scale_fill_manual(values = unname(c_pal(...)))
}



#' Integration with ggplot2 'scale_color_manual()' function
#'
#' @param ... Optional colors to be chosen from the palette
#' @export
c_scale_color <- function(...) {
  scale_color_manual(values = unname(c_pal(...)))
}


#' Compute the Statistical Mode
#'
#' Returns the mode (most frequent value) of a numeric or character vector.
#' If there are multiple modes, the first one encountered is returned.
#'
#' @param x A vector of numeric or character values.
#'
#' @return A single value representing the mode of the input vector.
#'
#' @examples
#' kmode(c(1, 2, 2, 3, 3, 3, 4))
#' # Returns: 3
#'
#' kmode(c("apple", "banana", "apple", "orange"))
#' # Returns: "apple"
#'
#' @export
kmode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


#' Compute groupwise means and confidence intervals using t-tests
#'
#' @param data A data frame or tibble.
#' @param group_col A bare (unquoted) grouping variable.
#' @param value_col A bare (unquoted) numeric variable.
#' @param conf_level Confidence level for the t-test (default = 0.95).
#'
#' @return A tibble with group means and confidence intervals.
#' @export
compute_ci <- function(
    data,
    ...,
    value_col,
    conf_level = 0.95,
    alternative = "two.sided") {
  data %>%
    group_by(...) %>%
    summarise(
      broom::tidy(stats::t.test(
        {{ value_col }},
        conf.level = conf_level,
        alternative = alternative
      )),
      .groups = "drop"
    ) %>%
    select(
      ...,
      ci_lower = conf.low,
      x_bar = estimate,
      ci_upper = conf.high
    )
}

#' Perform a t-test for group means
#'
#' @param data A data frame or tibble.
#' @param group_col A bare (unquoted) grouping variable.
#' @param value_col A bare (unquoted) numeric variable.
#' @param alpha Significance level for the t-test (default = 0.05).
#'
#' @return A tibble with t-test results
#' @export
t_test <- function(
    data,
    ...,
    value_col,
    mu_0 = 0,
    alpha = 0.05,
    alternative = "two.sided") {
  conf_level <- 1 - alpha

  result_lookup_tbl <- tibble::tibble(
    code = c("reject", "fail_reject"),
    label = c("Reject null hypothesis", "Fail to reject null hypothesis")
  )

  data %>%
    group_by(...) %>%
    summarise(
      broom::tidy(stats::t.test(
        {{ value_col }},
        mu = mu_0,
        conf.level = conf_level,
        alternative = alternative
      )),
      .groups = "drop"
    ) %>%
    mutate(
      alpha = alpha,
      test_result = if_else(
        p.value < alpha,
        "reject",
        "fail_reject"
      )
    ) %>%
    convert_codes_to_factor(
      code_col = test_result,
      lookup_tbl = result_lookup_tbl,
      lookup_code_col = code,
      lookup_label_col = label,
      new_factor_col_name = result
    ) %>%
    select(
      ...,
      p_value = p.value,
      alpha,
      result
    )
}

#' Density function for location-scale t-distribution
#'
#' @param x Vector of quantiles.
#' @param df Degrees of freedom.
#' @param mu Location parameter.
#' @param sigma Scale parameter (σ > 0).
#'
#' @return Vector of density values.
#' @export
dlst <- function(mu, df, x_bar = 0, s_n = 1) {
  dt((mu - x_bar) / s_n, df) / s_n
}

#' Cumulative distribution function for location-scale t-distribution
#'
#' @param q Vector of quantiles.
#' @param df Degrees of freedom.
#' @param mu Location parameter.
#' @param sigma Scale parameter (σ > 0).
#'
#' @return Vector of probabilities.
#' @export
plst <- function(q, df, x_bar = 0, s_n = 1) {
  pt((q - x_bar) / s_n, df)
}

#' Quantile function for location-scale t-distribution
#'
#' @param p Vector of probabilities.
#' @param df Degrees of freedom.
#' @param mu Location parameter.
#' @param sigma Scale parameter (σ > 0).
#'
#' @return Vector of quantiles.
#' @export
qlst <- function(p, df, x_bar = 0, s_n = 1) {
  x_bar + s_n * qt(p, df)
}

#' Random generation from location-scale t-distribution
#'
#' @param n Number of observations.
#' @param df Degrees of freedom.
#' @param mu Location parameter.
#' @param sigma Scale parameter (σ > 0).
#'
#' @return Vector of random deviates.
#' @export
rlst <- function(m, df, x_bar = 0, s_n = 1) {
  x_bar + s_n * rt(m, df)
}


#' Standard Scaling of a Numeric vector
#'
#' @param x A numeric vector to be scaled.
#' @param na.rm Logical, whether to remove NA values before scaling. Default TRUE.
#'
#' @return A numeric vector with the same length as `x`, containing scaled values.
#' @export
standard_scale <- function(x, na.rm = TRUE) {
  if (na.rm) {
    x <- na.omit(x)
  }

  mean_x <- mean(x, na.rm = TRUE)
  sd_x <- sd(x, na.rm = TRUE)

  if (sd_x == 0) {
    warning("Standard deviation is zero; returning NA for all values.")
    return(rep(NA, length(x)))
  }

  (x - mean_x) / sd_x
}

#' Min-Max Scaling of a Numeric Vector
#'
#' @param x A numeric vector to be scaled.
#' @param na.rm Logical, whether to remove NA values before scaling. Default TRUE.
#'
#' @return A numeric vector with the same length as `x`, containing scaled values.
#' @export
min_max_scale <- function(x, na.rm = TRUE) {
  if (na.rm) {
    x <- na.omit(x)
  }

  min_x <- min(x, na.rm = TRUE)
  max_x <- max(x, na.rm = TRUE)

  if (max_x == min_x) {
    warning("Max and min are equal; returning NA for all values.")
    return(rep(NA, length(x)))
  }

  (x - min_x) / (max_x - min_x)
}

#' Robust Scaling of a Numeric Vector
#'
#' @param x A numeric vector to be scaled.
#' @param na.rm Logical, whether to remove NA values before scaling. Default TRUE.
#'
#' @return A numeric vector with the same length as `x`, containing scaled values.
#' @export
robust_scale <- function(x, na.rm = TRUE) {
  if (na.rm) {
    x <- na.omit(x)
  }

  median_x <- median(x, na.rm = TRUE)
  iqr_x <- IQR(x, na.rm = TRUE)

  if (iqr_x == 0) {
    warning("Interquartile range is zero; returning NA for all values.")
    return(rep(NA, length(x)))
  }

  (x - median_x) / iqr_x
}

#' Normalizing Scaling of a Numeric Vector
#'
#' @param x A numeric vector to be scaled.
#' @param na.rm Logical, whether to remove NA values before scaling. Default TRUE.
#'
#' @return A numeric vector with the same length as `x`, of norm one.
#' @export
normal_scale <- function(x, na.rm = TRUE) {
  if (na.rm) {
    x <- na.omit(x)
  }

  norm_x <- sqrt(sum(x^2, na.rm = TRUE))

  if (norm_x == 0) {
    warning("Norm is zero; returning NA for all values.")
    return(rep(NA, length(x)))
  }

  x / norm_x
}



#' Perform Shapiro-Wilk normality test for groups
#'
#' @param data A data frame or tibble.
#' @param group_col A bare (unquoted) grouping variable.
#' @param value_col A bare (unquoted) numeric variable.
#'
#' @return A tibble with Shapiro-Wilk test results
#' @export
tidy_shapiro_test <- function(
    data,
    ...,
    value_col) {
  data %>%
    group_by(...) %>%
    summarise(
      broom::tidy(stats::shapiro.test(
        {{ value_col }}
      )),
      .groups = "drop"
    ) %>%
    select(
      ...,
      p_value = p.value,
      statistic = statistic,
      method = method
    )
}

#' Perform Kolgomorov-Smirnov normality test for groups
#'
#' @param data A data frame or tibble.
#' @param group_col A bare (unquoted) grouping variable.
#' @param value_col A bare (unquoted) numeric variable.
#'
#' @return A tibble with Kolgomorov-Smirnov test results
#' @export
tidy_ks_test <- function(
    data,
    ...,
    value_col) {
  data %>%
    group_by(...) %>%
    summarise(
      broom::tidy(stats::ks.test(
        {{ value_col }},
        "pnorm",
        mean = mean({{ value_col }}, na.rm = TRUE),
        sd = sd({{ value_col }}, na.rm = TRUE)
      )),
      .groups = "drop"
    ) %>%
    select(
      ...,
      p_value = p.value,
      statistic = statistic,
      method = method
    )
}


#' Perform Anderson-Darling normality test for groups
#'
#' @param data A data frame or tibble.
#' @param group_col A bare (unquoted) grouping variable.
#' @param value_col A bare (unquoted) numeric variable.
#'
#' @return A tibble with Anderson-Darling test results
#' @export
tidy_ad_test <- function(
    data,
    ...,
    value_col) {
  data %>%
    group_by(...) %>%
    summarise(
      broom::tidy(nortest::ad.test(
        {{ value_col }}
      )),
      .groups = "drop"
    ) %>%
    select(
      ...,
      p_value = p.value,
      statistic = statistic,
      method = method
    )
}

#' Perform Jarque-Bera normality test for groups
#'
#' @param data A data frame or tibble.
#' @param group_col A bare (unquoted) grouping variable.
#' @param value_col A bare (unquoted) numeric variable.
#'
#' @return A tibble with Jarque-Bera test results
#' @export
tidy_jb_test <- function(
    data,
    ...,
    value_col) {
  data %>%
    group_by(...) %>%
    summarise(
      broom::tidy(tseries::jarque.bera.test(
        {{ value_col }}
      )),
      .groups = "drop"
    ) %>%
    select(
      ...,
      p_value = p.value,
      statistic = statistic,
      method = method
    )
}



#' Compute p-values from multiple normality tests for a sample
#'
#' Applies several standard normality tests (Shapiro–Wilk,
#' Kolmogorov–Smirnov, Anderson–Darling, and Jarque–Bera) to a numeric
#' sample contained in a data frame, and returns the resulting p-values
#' in long (tidy) format.
#'
#' @param tbl A data frame or tibble containing the sample data.
#'
#' @param value_col A numeric column in \code{tbl} containing the values
#'   to be tested for normality. This argument uses tidy evaluation and
#'   must be supplied unquoted.
#'
#' @param sample_name A character string giving the name of the output
#'   column that will store the p-values (e.g. \code{"Gamma"},
#'   \code{"Normal"}, \code{"Lognormal"}).
#'
#' @return A tibble with two columns:
#' \describe{
#'   \item{test}{A character vector identifying the normality test
#'   (\code{"shapiro"}, \code{"ks"}, \code{"ad"}, \code{"jb"}).}
#'   \item{<sample_name>}{A numeric column containing the corresponding
#'   p-values for each test.}
#' }
#'
#' @details
#' The following tests are applied:
#' \itemize{
#'   \item Shapiro–Wilk test
#'   \item Kolmogorov–Smirnov test
#'   \item Anderson–Darling test
#'   \item Jarque–Bera test
#' }
#'
#' The function returns results in long format to facilitate reshaping
#' (e.g. via \code{pivot_wider()}) when comparing multiple samples.
#'
#' @examples
#' \dontrun{
#' compute_normality_pvalues(
#'   gamma_sample_tbl,
#'   sample,
#'   "Gamma"
#' )
#' }
#'
#' @seealso
#' \code{\link{tidy_shapiro_test}},
#' \code{\link{tidy_ks_test}},
#' \code{\link{tidy_ad_test}},
#' \code{\link{tidy_jb_test}}
#'
#' @export
compute_normality_pvalues <- function(data, value_col, sample_name) {

  # Generate the lookup table for test labels
  normality_tests_lookup_tbl <- tibble(
    code = c("shapiro", "ks", "ad", "jb"),
    label = c(
      "Shapiro-Wilk",
      "Kolmogorov-Smirnov",
      "Anderson-Darling",
      "Jarque-Bera"
    ),
  )


  # each function must return a list/data object with $p_value
  tests <- list(
    shapiro = tidy_shapiro_test,
    ks      = tidy_ks_test,
    ad      = tidy_ad_test,
    jb      = tidy_jb_test
  )

  # extract the one-column tibble for testing
  one_col_tbl <- dplyr::select(data, {{ value_col }})


  # compute p-values for each test
  pvals <- purrr::imap_dbl(
    tests,
    ~ .x(data = one_col_tbl, value_col = {{ value_col }})$p_value
  )

  # return as tidy tibble
  tibble::tibble(
    test = names(pvals),
    !!sample_name := unname(pvals)
  ) |> convert_codes_to_factor(
    code_col = test,
    lookup_tbl = normality_tests_lookup_tbl,
    lookup_code_col = code,
    lookup_label_col = label,
    new_factor_col_name = Prueba
  ) |>
    select(Prueba, everything(), -test)
}
