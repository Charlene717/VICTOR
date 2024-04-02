
# source("Set_plot_color.R")  # source("Set_pbmc3k_plot.R")

## 函數：進行ROC分析
roc_analysis <- function(seurat_obj, Set_ACT = "Actual Cell Type",
                         Set_Anno = "Annotation", Set_Col = " Score",
                         DefaultThr = 0.5) {

  ## Load Packages
  if(!require("pROC")) install.packages("pROC"); library(pROC)
  if(!require("ggplot2")) install.packages("ggplot2"); library(ggplot2)
  if(!require("cowplot")) install.packages("cowplot"); library(cowplot)


  cell_types <- unique(seurat_obj@meta.data[[Set_Anno]])
  results_ROC_list <- list()
  metadata_sub <- seurat_obj@meta.data
  # metadata_sub <- seurat_obj@meta.data[!is.na(seurat_obj@meta.data$ScoreLR), ]
  for (cell_type in cell_types) {
    value_col <- paste0(cell_type, Set_Col)
    if (value_col %in% names(metadata_sub)) {
      actual_binary <- ifelse(metadata_sub[[Set_ACT]] == cell_type, 1, 0)
      if (!any(actual_binary == 1)) {
        cat(paste0("Skipping ROC analysis for ", cell_type, " due to lack of case observations.\n"))
        next
      }
      roc_obj <- roc(actual_binary, metadata_sub[[value_col]], direction = "<")
      optimal_threshold <- mean(coords(roc_obj, "best")$threshold)
      if (any(is.infinite(optimal_threshold) | is.nan(optimal_threshold))) {
        optimal_threshold <- DefaultThr
      }
      results_ROC_list[[cell_type]] <- list(roc_obj = roc_obj, optimal_threshold = optimal_threshold)
    }
  }
  return(results_ROC_list)
}


## 函數：繪製ROC曲線
plot_ROC <- function(ROC_list, combined = FALSE, color_vector = NULL,
                     title_prefix = "", DefaultThr = 0.5) {
  combined_data_ROC <- data.frame()

  for (cell_type in names(ROC_list)) {
    roc_obj <- ROC_list[[cell_type]]$roc_obj
    optimal_threshold <- ROC_list[[cell_type]][["optimal_threshold"]]
    # 如果optimal_threshold包含任何Inf, -Inf或NaN
    if (any(is.infinite(optimal_threshold) | is.nan(optimal_threshold))) {
      optimal_threshold <- DefaultThr
    }

    curve_data <- data.frame(
      specificity = 1 - roc_obj$specificities,
      sensitivity = roc_obj$sensitivities,
      cell_type = cell_type,
      AUC = round(auc(roc_obj), 3),
      OptimalThreshold = round(optimal_threshold, 3)
    )
    combined_data_ROC <- rbind(combined_data_ROC, curve_data)
  }

  if (!combined) {
    # Individual ROC plots
    roc_plots <- list()
    for (cell_type in names(ROC_list)) {
      curve_data <- combined_data_ROC[combined_data_ROC$cell_type == cell_type,]

      # line_color <- "black"
      line_color <- scales::hue_pal()(length(names(ROC_list)) )[which(cell_type == names(ROC_list))]
      line_width <- 1.5

      # Check if there's an existing color_vector and if current cell type exists in it
      if (!is.null(color_vector) && cell_type %in% names(color_vector)) {
        line_color <- color_vector[[cell_type]]
      }

      roc_plot <- ggplot(curve_data, aes(x = specificity, y = sensitivity)) +
        geom_line(color = line_color, size = line_width) +
        theme_minimal() +
        labs(title = paste0(title_prefix, cell_type),
             subtitle = paste("AUC:", round(curve_data$AUC, 3),
                              ", Thr: ", round(curve_data$OptimalThreshold, 3)),
             x = "1 - Specificity", y = "Sensitivity") +
        theme(legend.position = "none")
      roc_plots[[cell_type]] <- roc_plot
    }

    combined_plot <- plot_grid(plotlist = roc_plots, labels = "AUTO", ncol = 3)

  } else {
    # Combined ROC plot

    combined_data_ROC$legend_label <- paste(combined_data_ROC$cell_type,
                                            " (AUC: ", combined_data_ROC$AUC,
                                            ", Thr: ", combined_data_ROC$OptimalThreshold, ")")

    combined_data_ROC$legend_label <- as.character(combined_data_ROC$legend_label)
    combined_data_ROC$extracted_cell_type <- trimws(gsub(" \\(AUC:.*$", "", combined_data_ROC$legend_label))

    # 使用extracted_cell_type来获取颜色
    if (!is.null(color_vector) && all(combined_data_ROC$extracted_cell_type %in% names(color_vector))) {
      combined_data_ROC$color <- color_vector[combined_data_ROC$extracted_cell_type]
      # 確保每一行都有正確的顏色映射，否則使用預設顏色
      combined_data_ROC$color[is.na(combined_data_ROC$color)] <-scales::hue_pal()(length(names(ROC_list)) )[which(cell_type == names(ROC_list))]
    } else {
      combined_data_ROC$color <- scales::hue_pal()(nrow(combined_data_ROC))
    }

    # 繪圖
    # Before plotting, create a named color vector that ensures that each unique legend label gets its own color.
    named_color_vector <- setNames(as.character(combined_data_ROC$color), combined_data_ROC$legend_label)

    combined_plot <- ggplot(combined_data_ROC, aes(x = specificity, y = sensitivity, color = legend_label, group = legend_label)) +
      geom_line(size = 1.5) +
      theme_minimal() +
      labs(title = paste("Combined ROC Curves", title_prefix), x = "1 - Specificity", y = "Sensitivity") +
      scale_color_manual(values = named_color_vector, name = "Cell Type") +
      theme(
        axis.text.x = element_text(angle = 0, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16, angle = 0),
        title = element_text(size = 18, face = "bold"),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 11),
        panel.border = element_rect(color = "black", fill = NA, size = 1.5),
        aspect.ratio = 1
      )

  }

  return(list(plots = combined_plot, df = combined_data_ROC))
}


## 函數：生成結果的DataFrame
generate_roc_data <- function(ROC_list, title_prefix = "") {
  results_ROC_df <- data.frame()
  for (cell_type in names(ROC_list)) {
    roc_obj <- ROC_list[[cell_type]]$roc_obj
    optimal_threshold <- ROC_list[[cell_type]]$optimal_threshold
    temp_df <- data.frame(
      CellType = cell_type,
      AUC = auc(roc_obj),
      OptimalThreshold = optimal_threshold,
      Source = title_prefix
    )
    results_ROC_df <- rbind(results_ROC_df, temp_df)
  }
  return(results_ROC_df)
}

# ## 使用示例：
# ## 假設你的Seurat對象是 seuratObject_Ref 和 seuratObject_Sample
# ROC_list_Sample <- roc_analysis(seuratObject_Sample)
# ROC_list_Ref <- roc_analysis(seuratObject_Ref, Set_Anno = "Actual Cell Type")
#
# results_ROC_df_Sample <- generate_roc_data(ROC_list_Sample, title_prefix = "Sample")
# results_ROC_df_Ref <- generate_roc_data(ROC_list_Ref, title_prefix = "Reference")
#
# ROC_LR_results_sample_combined <- plot_ROC(ROC_list_Sample, combined = TRUE, color_vector=color_vector, title_prefix = " - Sample")
# ROC_LR_results_sample_combined$plots
# ROC_LR_results_sample_combined <- plot_ROC(ROC_list_Sample, combined = TRUE, title_prefix = " - Sample")
# ROC_LR_results_sample_combined$plots
# ROC_LR_results_sample <- plot_ROC(ROC_list_Sample, combined = FALSE, color_vector=color_vector, title_prefix = "Sample ")
# ROC_LR_results_sample$plots
# ROC_LR_results_sample <- plot_ROC(ROC_list_Sample, combined = FALSE, title_prefix = "Sample ")
# ROC_LR_results_sample$plots
#
# ROC_LR_results_ref_combined <- plot_ROC(ROC_list_Ref, combined = TRUE, color_vector=color_vector, title_prefix = " - Reference")
# ROC_LR_results_sample_combined$plots
# ROC_LR_results_ref <- plot_ROC(ROC_list_Ref, combined = FALSE, color_vector=color_vector, title_prefix = "Reference ")
#
