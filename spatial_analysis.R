library(tidyverse)   
library(ggplot2)
library(spatstat)
library(deldir)
library(viridis)
library(gridExtra)

setwd("D:/Data/Nucleospira")   

import_distance_csv <- function(file_path) {
  raw <- read.csv(file_path, check.names = FALSE)

  # coordinate file: at least two columns whose names contain “x” and “y”
  has_xy <- any(grepl("^x$", colnames(raw), ignore.case = TRUE)) &&
            any(grepl("^y$", colnames(raw), ignore.case = TRUE))

  if (has_xy) {
    pts <- raw
  } else {
    message("Distance matrix detected, reconstructing coordinates with classical MDS …")
    m <- as.matrix(raw[,-1])
    rownames(m) <- raw[,1]
    mds <- cmdscale(m, k = 2)
    pts <- data.frame(
      id = rownames(mds),
      x  = mds[,1],
      y  = mds[,2]
    )
  }

  if (!"id" %in% names(pts)) pts$id <- paste0("P", seq_len(nrow(pts)))
  names(pts) <- sub("^x$", "x", names(pts), ignore.case = TRUE)
  names(pts) <- sub("^y$", "y", names(pts), ignore.case = TRUE)
  pts[,c("id","x","y")]
}

distance_matrix <- function(pts) {
  as.matrix(dist(pts[,c("x","y")], upper = TRUE, diag = TRUE))
}

nearest_distances <- function(pts) {
  d <- distance_matrix(pts)
  apply(d + diag(Inf, nrow(d)), 1, min)
}

nearest_neighbor_analysis <- function(pts, study_area = NULL) {
  nd     <- nearest_distances(pts)
  rA     <- mean(nd)

  if (is.null(study_area)) {
    bbox <- range(pts$x); width  <- diff(bbox)
    bboy <- range(pts$y); height <- diff(bboy)
    study_area <- width * height
    message(sprintf("Study area inferred from bounding box: %.2f", study_area))
  }

  n       <- nrow(pts)
  density <- n / study_area
  rE      <- 1 / (2 * sqrt(density))
  R       <- rA / rE
  SE      <- 0.26136 / sqrt(n * density)
  Z       <- (rA - rE) / SE
  type    <- case_when(R < 0.9 ~ "clustered",
                       R > 1.1 ~ "regular",
                       TRUE    ~ "random")
  signif  <- case_when(abs(Z) > 2.58 ~ "p < 0.01",
                       abs(Z) > 1.96 ~ "p < 0.05",
                       TRUE          ~ "not significant")
  list(n = n, area = study_area, density = density, rA = rA, rE = rE,
       R = R, Z = Z, CV = sd(nd)/rA, type = type, signif = signif,
       nd = nd)
}

plot_points <- function(pts) {
  ggplot(pts, aes(x, y)) +
    geom_point(color = "red", size = 3) +
    geom_text(aes(label = id), vjust = -1, size = 3) +
    theme_minimal() +
    labs(title = "Point distribution", x = "X", y = "Y") +
    coord_equal()
}

plot_voronoi <- function(pts) {
  xr <- range(pts$x); yr <- range(pts$y)
  m  <- 0.1 * max(diff(xr), diff(yr))
  v  <- deldir(pts$x, pts$y, rw = c(xr[1]-m, xr[2]+m, yr[1]-m, yr[2]+m))
  tiles <- tile.list(v)
  polygons <- lapply(seq_along(tiles), function(i) {
    data.frame(id = i, x = tiles[[i]]$x, y = tiles[[i]]$y)
  })
  vor <- do.call(rbind, polygons)
  ggplot() +
    geom_polygon(data = vor, aes(x, y, group = id),
                 fill = NA, colour = "#2196F3", size = 0.7) +
    geom_point(data = pts, aes(x, y), colour = "red", size = 3) +
    theme_minimal() +
    labs(title = "Voronoi diagram", x = "X", y = "Y") +
    coord_equal()
}

plot_nd_hist <- function(nd) {
  ggplot(data.frame(d = nd), aes(d)) +
    geom_histogram(bins = 10, fill = "#4CAF50", colour = "white") +
    theme_minimal() +
    labs(title = "Nearest‑neighbor distances", x = "Distance", y = "Count")
}

print_results <- function(res) {
  cat("\n=== Nearest‑neighbor analysis ===\n")
  cat(sprintf("Points: %d\nArea: %.2f\nDensity: %e\n", res$n, res$area, res$density))
  cat(sprintf("Mean nn distance (rA): %.2f\nExpected (rE): %.2f\n", res$rA, res$rE))
  cat(sprintf("R index: %.4f\nZ score: %.4f\nCV: %.4f\n", res$R, res$Z, res$CV))
  cat(sprintf("Pattern: %s (%s)\n", res$type, res$signif))
  cat(sprintf("Min: %.2f  Max: %.2f  Median: %.2f  SD: %.2f\n",
              min(res$nd), max(res$nd), median(res$nd), sd(res$nd)))
}

analyze_spatial_distribution <- function(file_path, study_area = NULL,
                                         output_prefix = "spatial_analysis") {

  message("Loading data …")
  pts <- import_distance_csv(file_path)
  message(sprintf("Loaded %d points.", nrow(pts)))

  message("Running nearest‑neighbor analysis …")
  res <- nearest_neighbor_analysis(pts, study_area)
  print_results(res)

  message("Generating plots …")
  p1 <- plot_points(pts)
  p2 <- plot_voronoi(pts)
  p3 <- plot_nd_hist(res$nd)

  message("Saving outputs …")
  ggsave(paste0(output_prefix, "_points.png"),   p1, width = 8, height = 6)
  ggsave(paste0(output_prefix, "_voronoi.png"),  p2, width = 8, height = 6)
  ggsave(paste0(output_prefix, "_nn_hist.png"),  p3, width = 8, height = 6)

  summary_csv <- data.frame(
    Parameter = c("Points", "Area", "Density", "rA", "rE",
                  "R", "Z", "CV", "Pattern", "Significance"),
    Value     = c(res$n, round(res$area,2),
                  format(res$density, scientific = TRUE, digits = 4),
                  round(res$rA,2), round(res$rE,2),
                  round(res$R,4),  round(res$Z,4),
                  round(res$CV,4), res$type, res$signif)
  )
  write.csv(summary_csv, paste0(output_prefix, "_summary.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")

  nd_csv <- data.frame(ID = pts$id, X = pts$x, Y = pts$y,
                       Nearest_Distance = res$nd)
  write.csv(nd_csv, paste0(output_prefix, "_nearest_distances.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")

  message("Done.")
  invisible(list(result = res, points = pts,
                 plots = list(points = p1, voronoi = p2, histogram = p3)))
}

# USAGE:
# analyze_spatial_distribution("Measurements.csv")