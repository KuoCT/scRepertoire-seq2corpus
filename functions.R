# seurat_workflow1 ----
seurat_workflow1 <- function(
    obj, 
    out_path = out2,
    cutoff = c(200, 7500, 10), 
    lim_nFeature_RNA = c(0, 10000),
    lim_nCount_RNA = c(0, 45000),
    lim_percent_mt = c(0, 35)
  ) {
  # Create out path
  if (!dir.exists(out_path)) dir.create(out_path)
  
  DefaultAssay(obj) <- "RNA"
  
  # Calculate the percentage of reads that map to the mitochondrial genome
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  
  # Extract project name
  obj_name <- Project(obj)
  
  # QC - violin plot
  features <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
  limits <- list(lim_nFeature_RNA, lim_nCount_RNA, lim_percent_mt)
  update_geom_defaults("violin", aes(linewidth = NA, alpha = 0))
  ps <- lapply(1:length(features), function(i) {
    p <- VlnPlot(
      object = obj, 
      features = features[i], 
      pt.size = 0, 
      adjust = 2
    ) + 
      geom_jitter(width = 0.45, height = 0, alpha = 0.03) +
      geom_violin(adjust = 2, linewidth = 1, alpha = 0.8) +
      ylim(limits[[i]]) +
      theme_prism(base_size = 12) +
      NoLegend()
    return(p)
  })
  
  p <- wrap_plots(ps, ncol = 3, widths = c(1, 1.05, 1.08))
  
  file_name <- paste0(obj_name, "_QC_vlnplot.svg")
  file_path <- file.path(out_path, file_name)
  svglite(file_path, width = 8, height = 5)
  print(p)
  dev.off()
  
  # This is a patch for remove the duplicate violin object in .svg file
  lines <- readLines(file_path)
  pattern <- "polygon\\s+points=.*style=['\"]?[^'\"]*stroke-width:\\s*nan;\\s*stroke:\\s*none;\\s*stroke-linecap:\\s*butt[^'\"]*['\"]?"
  lines <- lines[!grepl(pattern, lines)]
  writeLines(lines, file_path)
  
  # QC - scatter plot
  p1 <- FeatureScatter(
    obj, 
    feature1 = "nCount_RNA", 
    feature2 = "percent.mt", 
    shuffle = T,
    pt.size = 1.5
  ) +
    coord_cartesian(xlim = limits[[2]], ylim = limits[[3]]) +
    theme_prism(base_size = 10)
  
  p2 <- FeatureScatter(
    obj, 
    feature1 = "nCount_RNA", 
    feature2 = "nFeature_RNA",
    shuffle = T,
    pt.size = 1.5
  ) +
    coord_cartesian(xlim = limits[[2]], ylim = limits[[1]]) +
    theme_prism(base_size = 10)
  
  file_path <- paste0(out_path, "/", obj_name, "_QC_Scatterplot.svg")
  svglite(file_path, width = 9, height = 4)
  print(p1 + p2)
  dev.off()
  
  # QC - subset
  obj <- subset(obj, subset = nFeature_RNA > cutoff[1] & nFeature_RNA < cutoff[2] & percent.mt < cutoff[3])
  
  # Normalize
  obj <- NormalizeData(obj)
  
  # Find variable features
  obj <- FindVariableFeatures(obj, nfeatures = 3000)
  
  # Scale data
  obj <- ScaleData(obj)
}

# preset functions ----
UMAP_preset <- function(
    plot_title = "UMAP",
    legend_point_size = 2.8, 
    nrow_legend = 2,
    base_size = 10
) {
  list(
    theme_prism(base_size = base_size),
    labs(title = plot_title),
    theme(
      plot.title = element_text(hjust = 0),
      legend.position = "bottom",
      legend.box = "horizontal"
    ),
    guides(
      color = guide_legend(
        nrow = nrow_legend,
        byrow = TRUE,
        override.aes = list(size = legend_point_size)
      )
    )
  )
}

# plot_umap_density ----
plot_umap_density <- function(obj, Embedding = "umap") {
  umap_df <- as.data.frame(Embeddings(obj, Embedding)) %>% sample_frac(1) # mimic "shuffle" behavior
  print(head(umap_df))
  xcol <- colnames(umap_df)[1]
  ycol <- colnames(umap_df)[2]
  
  dp <- DimPlot(obj, reduction = Embedding)
  xlim <- ggplot_build(dp)$layout$panel_params[[1]]$x.range
  ylim <- ggplot_build(dp)$layout$panel_params[[1]]$y.range
  xbreaks <- ggplot_build(dp)$layout$panel_params[[1]]$x.sec$breaks
  ybreaks <- ggplot_build(dp)$layout$panel_params[[1]]$y.sec$breaks
  
  p <- ggplot(umap_df, aes(x = .data[[xcol]], y = .data[[ycol]])) +
    geom_point(size = 1) +
    stat_density_2d(aes(fill = after_stat(level)), geom = "polygon", bins = 15) +
    scale_x_continuous(limits = xlim * 1.05, breaks = xbreaks) + 
    scale_y_continuous(limits = ylim * 1.05, breaks = ybreaks) +
    scale_fill_viridis_c(option = "B") +
    theme_prism(base_size = 10) +
    labs(title = "UMAP density")
  return(p)
}

# FeaturePlot2 ----
FeaturePlot2 <- function(
    obj, 
    features, 
    min.cutoff = NA,
    max.cutoff = NA, 
    outfolder, 
    colours = replace(viridis_pal(option = "B")(10), 1, "gray90"),
    dim_x_range = NULL,
    dim_y_range = NULL,
    x_limits = waiver(),
    x_breaks = waiver(),
    y_limits = waiver(),
    y_breaks = waiver(),
    plot_width = 4,
    plot_height = 4.68
) {
  if (!dir.exists(outfolder)) dir.create(outfolder, recursive = T)
  for (feature in features) {
    FeaturePlot(
      obj,
      features = feature,
      order = T,
      pt.size = 1,
      min.cutoff = min.cutoff,
      max.cutoff = max.cutoff
    ) + 
      scale_color_gradientn(colours = colours, name = "Scaled Exp.") +
      coord_cartesian(xlim = dim_x_range, ylim = dim_y_range) +
      scale_x_continuous(breaks = x_breaks, limits = x_limits) +
      scale_y_continuous(breaks = y_breaks, limits = y_limits) + 
      theme_prism(base_size = 10) +
      theme(plot.title = element_text(hjust = 0),
            plot.subtitle = element_text(hjust = 0),
            legend.title = element_text(vjust = 0.75),
            legend.ticks = element_line(colour = "black"),
            legend.frame = element_rect(colour = "black"),
            legend.position = "bottom") -> p
    
    file_name <- paste0("FeaturePlot_", feature, ".svg")
    file_path <- file.path(outfolder, file_name)
    svglite(file_path, width = plot_width, height = plot_height)
    print(p)
    dev.off()
  }
}

# FeaturePlot3 ----
FeaturePlot3 <- function(
    obj, 
    features, 
    min.cutoff = NA,
    max.cutoff = NA, 
    colours = replace(viridis_pal(option = "B")(10), 1, "gray90"),
    dim_x_range = NULL,
    dim_y_range = NULL
) {
  plots <- list()
  
  for (feature in features) {
    p <- FeaturePlot(
      obj,
      features = feature,
      order = TRUE,
      pt.size = 1,
      min.cutoff = min.cutoff,
      max.cutoff = max.cutoff
    ) +
      scale_color_gradientn(colours = colours, name = "Scaled Exp.") +
      coord_cartesian(xlim = dim_x_range, ylim = dim_y_range) +
      theme_prism(base_size = 10) +
      theme(
        plot.title = element_text(hjust = 0),
        plot.subtitle = element_text(hjust = 0),
        legend.title = element_text(vjust = 0.75),
        legend.ticks = element_line(colour = "black"),
        legend.frame = element_rect(colour = "black"),
        legend.position = "bottom"
      )
    
    plots[[feature]] <- p
  }
  
  if (length(features) == 1) {
    return(plots[[1]])  # ggplot object
  } else {
    return(plots)       # list of ggplot objects
  }
}

# AddModuleScore2 ----
AddModuleScore2 <- function(obj, markers, ScoreName, assay) {
  existing_scores <- colnames(obj@meta.data)
  if (ScoreName %in% existing_scores) {
    obj@meta.data <- obj@meta.data[, !(colnames(obj@meta.data) %in% ScoreName)]
  }
  
  obj <- AddModuleScore(
    obj,
    features = list(markers),
    name = ScoreName,
    assay = assay
  )
  
  colnames(obj@meta.data) <- gsub(
    pattern = paste0(ScoreName, 1),
    replacement = ScoreName,
    colnames(obj@meta.data)
  )
  
  return(obj)
}

# RenameSeuratClusters ----
# Function to rename and reorder clusters
RenameSeuratClusters <- function(
    seu_object, 
    cluster_names, 
    order_levels = NULL, 
    new_column_name = "renamed_clusters"
) {
  # Ensure the number of provided cluster names matches the number of clusters
  if (length(cluster_names) != length(unique(Idents(seu_object)))) {
    stop("The number of cluster names must match the number of clusters in the object.")
  }
  
  # Create a named list for renaming clusters
  cluster_mapping <- setNames(cluster_names, levels(Idents(seu_object)))
  
  # Rename the clusters
  seu_object <- RenameIdents(seu_object, cluster_mapping)
  
  # Reorder the clusters if `order_levels` is provided
  if (!is.null(order_levels)) {
    if (!all(order_levels %in% cluster_names)) {
      stop("All order levels must be present in the cluster names.")
    }
    Idents(seu_object) <- factor(Idents(seu_object), levels = order_levels)
  }
  
  # Save the renamed clusters in a new metadata column
  seu_object[[new_column_name]] <- Idents(seu_object)
  
  return(seu_object)
}

# AssembleTrainingData ----
AssembleTrainingData <- function(seu_obj, scFv_dict, meatdata, genes) {
  expr <- FetchData(seu_obj, vars = c(meatdata, genes)) %>% rownames_to_column("seu_barcode")
  TrainingData <- left_join(scFv_dict, expr, by = "seu_barcode") %>% 
    select(
      any_of(c(
        "seu_barcode", 
        "scFv_type", 
        "linker_type", 
        "linker_type_protein", 
        meatdata,
        "scFv_length",
        "Full_V_domain_nt", 
        "Full_V_domain", 
        "V_left", 
        "linker", 
        "V_right",
        genes
      ))
    )
  return(TrainingData)
}


# add_GSEA_heatmap ----
add_GSEA_heatmap <- function(
    svg_path, 
    score_path, 
    output_folder,
    y_pos = 300, 
    height = 25,
    x_start = 100.73, 
    x_interval = (174.47 - 100.73) / 1000,
    palette = c("navy", "white", "firebrick3")
) {
  # Read SVG source code
  svg_lines <- readLines(svg_path, warn = FALSE)
  
  # Read score table
  df <- read.csv(score_path, row.names = NULL, sep = "\t")
  scores <- df$score
  
  # Order score
  df <- df[order(-df$score), ]
  scores <- df$score
  n <- length(scores)
  # print(n)
  
  # Color mapping
  lower_q <- quantile(scores, 0.10, na.rm = TRUE)
  upper_q <- quantile(scores, 0.90, na.rm = TRUE)
  mid <- 0
  scores_clipped <- squish(scores, range = c(lower_q, upper_q)) # Clipping
  
  neg_fn <- scales::col_numeric(palette = palette[1:2], domain = c(lower_q, mid))
  pos_fn <- scales::col_numeric(palette = palette[2:3], domain = c(mid, upper_q))
  
  mapped_colors <- vapply(scores_clipped, function(s) {
    if (s == mid) {
      palette[2]  # white
    } else if (s < mid) {
      neg_fn(s)
    } else {
      pos_fn(s)
    }
  }, character(1))
  
  # Coordinate mapping
  x_end <- x_start + x_interval * (n - 1)
  x_coords <- seq(x_start, x_end, length.out = n)
  rect_width <- x_interval
  
  # Build rect elements
  heatmap_rects <- paste0(
    "<rect x='", round(x_coords, 2), "' y='", y_pos, "' width='", round(rect_width, 2),
    "' height='", height, "' fill='", mapped_colors, "' stroke='none' />"
  )
  heatmap_svg <- c("<g id='heatmap-bar'>", heatmap_rects, "</g>")
  
  # Encode into SVG
  # insert_at <- grep("</svg>", svg_lines)[1] - 1  # Insert before the end
  insert_at <- grep("<svg[^>]*>", svg_lines)[1]  # Insert at top
  new_svg <- append(svg_lines, heatmap_svg, after = insert_at)
  
  # Ensure output folder exists
  if (!dir.exists(output_folder)) dir.create(output_folder, recursive = T)
  
  # Extract file name
  file_name <- basename(svg_path)
  file_base <- sub("\\.svg$", "", file_name)
  output_file <- file.path(output_folder, paste0(file_base, "_color.svg"))
  
  # Write SVG
  writeLines(new_svg, output_file)
  message("Heatmap bar added and saved to: ", output_file)
}

# read_GSEA ----
read_GSEA <- function(out_folder, cluser) {
  fp <- file.path(out_folder, glue::glue("GSEA/Project_{cluser}/enrichment_results_{cluser}.txt"))
  # Use read.delim for tab-separated
  df <- read.delim(fp, header = TRUE, sep = "\t", check.names = FALSE)
  df$cluster <- cluser
  df
}

# filter_GSEA_by_selection ----
filter_GSEA_by_selection <- function(df, selection_list) {
  keep <- imap(selection_list, ~ tibble(cluster = .y, description = .x)) %>%
    list_rbind()
  if (nrow(keep) == 0) {
    df[0, ]  # empty tibble if no selection specified
  } else {
    df %>% semi_join(keep, by = c("cluster", "description"))
  }
}

# add_GSEA_spacers ----
add_GSEA_spacers <- function(df) {
  last_cluster <- last(unique(df$cluster))
  df %>%
    group_by(cluster) %>%
    group_modify(function(.x, .y) {
      spacer <- tibble(
        label = paste0("<<SPACER_", .y$cluster, ">>"),
        NES = NA_real_,
        FDR = NA_real_,
        leading_edge_size = NA_real_
      )
      bind_rows(.x, spacer)
    }) %>%
    ungroup() %>%
    filter(label != paste0("<<SPACER_", last_cluster, ">>")) %>% # remove last row 
    mutate(label = forcats::fct_inorder(label))
}

# ExtractSeuratDataFrame ----
ExtractSeuratDataFrame <- function(
    seu_object, 
    reduction = NULL, 
    metadata = NULL, 
    genes = NULL
){
  # Initialize data.frame
  result_df <- data.frame(seu_barcode = colnames(seu_object), stringsAsFactors = F)
  
  # Extract embeddings
  if (!is.null(reduction)) {
    embedding_df <- Embeddings(seu_object, reduction = reduction) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("seu_barcode")
    result_df <- dplyr::left_join(result_df, embedding_df, by = "seu_barcode")
  }
  
  # Extract specified metadata columns
  if (!is.null(metadata) && length(metadata) > 0) {
    meta_available <- intersect(metadata, colnames(seu_object@meta.data))
    meta_df <- seu_object@meta.data[, meta_available, drop = FALSE] %>%
      tibble::rownames_to_column("seu_barcode")
    result_df <- dplyr::left_join(result_df, meta_df, by = "seu_barcode")
  }
  
  # Extract genes
  if (!is.null(genes) && length(genes) > 0) {
    expr_mat <- Seurat::GetAssayData(seu_object, assay = "RNA", layer = "data")
    genes_valid <- intersect(genes, rownames(expr_mat))
    
    if (length(genes_valid) > 0) {
      gene_expr_df <- expr_mat[genes_valid, , drop = FALSE] %>%
        as.matrix() %>%
        t() %>%
        as.data.frame() %>%
        tibble::rownames_to_column("seu_barcode")
      result_df <- dplyr::left_join(result_df, gene_expr_df, by = "seu_barcode")
    } else {
      warning("No valid genes found in the RNA assay.")
    }
  }
  
  return(result_df)
}

# AppendSeuratMetadata ----
AppendSeuratMetadata <- function(seu_object, df_to_Append, barcode_col = "seu_barcode") {
  metadata <- seu_object[[]] %>% 
    rownames_to_column(barcode_col) %>%
    left_join(df_to_Append, by = barcode_col) %>%
    column_to_rownames(barcode_col)
  seu_object <- AddMetaData(seu_object, metadata = metadata)
  return(seu_object)
  
}

# BubblePlot ----
BubblePlot <- function(
    seu_object, 
    used_features,
    row_hclust = T,
    row_expand = c(0.02, 0.02),
    col_hclust = T,
    col_expand = c(0.1, 0.1),
    rel_widths = c(1, 5), 
    rel_heights = c(1, 25)
) {
  # Get wide format for clustering
  df <- DotPlot(seu_object, features = used_features)$data
  df_wide <- df %>%
    select(id, features.plot, avg.exp.scaled) %>%
    pivot_wider(names_from = id, values_from = avg.exp.scaled)
  mat <- as.matrix(df_wide[,-1])
  rownames(mat) <- df_wide$features.plot
  
  # Hierarchical clustering
  if (row_hclust){
    hc_rows <- hclust(dist(mat))
    row_order <- rownames(mat)[hc_rows$order]
    
    # Dendrogram plot (rows)
    dendro_data_rows <- ggdendro::dendro_data(hc_rows)
    dendro_plot_rows <- ggplot(dendro_data_rows$segments) +
      geom_segment(aes(x = y, y = x, xend = yend, yend = xend)) +
      scale_y_continuous(
        breaks = seq_along(row_order),
        labels = row_order,
        expand = row_expand
      ) +
      scale_x_continuous(trans = "reverse") +
      theme_void()
    
    df$features.plot <- factor(df$features.plot, levels = row_order)
  }
  
  if (col_hclust){
    hc_columns <- hclust(dist(t(mat)))
    col_order <- colnames(mat)[hc_columns$order]
    
    # Dendrogram plot (columns)
    dendro_data_columns <- ggdendro::dendro_data(hc_columns)
    dendro_plot_columns <- ggplot(dendro_data_columns$segments) +
      geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
      scale_x_continuous(
        breaks = seq_along(col_order),
        labels = col_order,
        expand = col_expand
      ) + 
      theme_void()
    
    df$id <- factor(df$id, levels = col_order)
  }
  
  # Plot
  dotplot <- ggplot(df, aes(x = id, y = features.plot)) +
    geom_point(aes(size = pct.exp, fill = avg.exp.scaled), shape = 21, color = "black", stroke = 0.4) +
    scale_fill_viridis_c() +
    scale_y_discrete(expand = row_expand, position = "right") +
    scale_x_discrete(expand = col_expand) +
    theme_prism(base_size = 10, border = T) +
    theme(
      panel.grid.major = element_line(color = "grey80", size = 0.3, linetype = "dashed"),
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.ticks.y.right = element_blank(),
      axis.title = element_blank(),
      axis.ticks.length = unit(2, "pt"),
      legend.title = element_text(),
      legend.ticks = element_line(colour = "black"),
      legend.frame = element_rect(colour = "black"),
      panel.border = element_rect(size = 1)
    )
  
  if (row_hclust & col_hclust) {
    dotplot <- dotplot + theme(plot.margin = margin(l = 0, t = 0))
    p <- dendro_plot_columns + dendro_plot_rows + dotplot + plot_layout(design = "\n#A\nBC\n", widths = rel_widths, heights = rel_heights)
  } else if (row_hclust) {
    dotplot <- dotplot + theme(plot.margin = margin(l = 0))
    p <- dendro_plot_rows + dotplot + plot_layout(widths = rel_widths)
  } else if (col_hclust) {
    dotplot <- dotplot + theme(plot.margin = margin(t = 0))
    p <- dendro_plot_columns + dotplot + plot_layout(design = "\nA\nB\n", heights = rel_heights)
  } else {
    p <- dotplot
  }
  
  return(p)
}