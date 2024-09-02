############ Format of function files  ############
############--------------------------------------------------------------------#############' A Self-made Median Function
#'
#' This function allows you to diagnose cell annotation labels.
#' @param seuratObj_Query Seurat Object of Query.
#' @param seuratObj_Ref Seurat Object of Reference.
#' @param ActualCellTypeColumn Column name of the actual cell type in meta.data.
#' @param AnnotCellTypeColumn Column name of the cell type annotation in meta.data.
#' @keywords scRNA-seq
#' @export
#' @examples
#' VICTOR(seuratObj_Query, seuratObj_Ref, ActualCellTypeColumn = "Actual_Cell_Type", AnnotCellTypeColumn = "Annotation", Add_nROC = FALSE)

VICTOR <- function(seuratObj_Query, seuratObj_Ref,
                   ActualCellTypeColumn = "Actual_Cell_Type",
                   AnnotCellTypeColumn = "Annotation",
                   Add_nROC = FALSE,
                   seurat_version = "V5",...) {

  if(!require("pbapply")) install.packages("pbapply"); library(pbapply)
  if(!require("caret")) install.packages("caret"); library(caret)

  if (!requireNamespace("remotes", quietly = TRUE)) {install.packages("remotes")}
  if(!require("harmony")) remotes::install_github("immunogenomics/harmony"); library(harmony)

  # source("FUN_VICTORPrep.R")

  #### Unified column name ####
  seuratObj_Ref[["Actual_Cell_Type"]] <- seuratObj_Ref[[ActualCellTypeColumn]]
  # seuratObj_Query[["Annotation"]] <- seuratObj_Query[[AnnotCellTypeColumn]]

  #### VICTORPrep score ####
  if(is.null(seuratObj_Ref@misc[["VICTOR"]])) {
    seuratObj_Ref <- getFeatureSpace(seuratObj_Ref, "Actual_Cell_Type", seurat_version = seurat_version,...)  ## Get the feature space to train the classifiers
  }

  if(length(seuratObj_Ref@misc[["VICTOR"]]@train) == 0) {
    seuratObj_Ref <- trainModel(seuratObj_Ref,model = "glmnet") # Train the model
  }

  ## Ref
  if (!any(grepl("VICTORScore$", colnames(seuratObj_Ref@meta.data)))) { # if(!"VICTOR_max" %in% colnames(seuratObj_Ref@meta.data)) {
    seuratObj_Ref <- VICTORPrep(seuratObj_Ref, seuratObj_Ref,seurat_version = seurat_version, ...) #, threshold = Set_scPredict_Thr)
  }

  ## Query
  if (!any(grepl("VICTORScore$", colnames(seuratObj_Query@meta.data)))) {  # if(!"VICTOR_max" %in% colnames(seuratObj_Query@meta.data)) {
    seuratObj_Query <- VICTORPrep(seuratObj_Query, seuratObj_Ref,seurat_version = seurat_version,...) #, threshold = Set_scPredict_Thr)
  }

  #### Modify column names ####
  modify_metadata_cols_VICTOR <- function(seuratObj) {
    colnames(seuratObj@meta.data) <- sapply(colnames(seuratObj@meta.data), function(name) {
      if(grepl("^VICTOR_", name) && !name %in% c("VICTOR_max", "VICTOR_prediction", "VICTOR_no_rejection")) {
        paste0(gsub("_", " ", gsub("_plus", "+", gsub("^VICTOR_", "", name))), " VICTORScore")
      } else { name }
    })
    return(seuratObj) # Return the modified Seurat object
  }

  seuratObj_Ref = modify_metadata_cols_VICTOR(seuratObj_Ref)
  seuratObj_Query = modify_metadata_cols_VICTOR(seuratObj_Query)

  #### Ref: Optimization threshold ####
  # source("ROC_Analysis.R")
  ROC_Ref.lt <-   roc_analysis(seuratObj_Ref, Set_ACT = "Actual_Cell_Type",
                               Set_Anno = "Actual_Cell_Type", Set_Col = " VICTORScore",
                               DefaultThr = 0.5)
  ROC_Ref.df <- generate_roc_data(ROC_Ref.lt, title_prefix = paste("Reference"))
  seuratObj_Ref@misc[["ROC"]][["ROC_Ref.lt"]] <- ROC_Ref.lt
  seuratObj_Ref@misc[["ROC"]][["ROC_Ref.df"]] <- ROC_Ref.df


  #### VICTOR Diagnosis ####
  score_type <- "VICTOR"
  metadata <- seuratObj_Query@meta.data
  metadata_col <- paste0(AnnotCellTypeColumn, "_", score_type, "Score")

  ## Calculate the score and update metadata
  metadata[[metadata_col]] <- apply(metadata, 1, function(row) {
    score_column <- paste0(row[[AnnotCellTypeColumn]], " ", score_type, "Score")
    if (score_column %in% colnames(metadata)) row[[score_column]] else NA
  })

  metadata[[metadata_col]] <- metadata[[metadata_col]] %>% as.numeric()

  if(Add_nROC){
    ## State
    ScoreStat_Thr <- 0.5
    metadata[[paste0("Diag_", score_type,"_", AnnotCellTypeColumn, "_Stat")]] <- ifelse(metadata[[metadata_col]] < ScoreStat_Thr, "F", "T")
  }

  ## State ROC
  thresholds <- sapply(metadata[[AnnotCellTypeColumn]], function(cell_type) {
    if (cell_type %in% ROC_Ref.df$CellType) {
      ROC_Ref.df[ROC_Ref.df$CellType == cell_type, "OptimalThreshold"]
    } else {
      NA
    }
  })

  # metadata[[paste0("Diag_", score_type,"_", AnnotCellTypeColumn,  "_StatROC")]] <- ifelse(metadata[[metadata_col]] < thresholds, "F", "T")
  metadata[[paste0("Diag_", score_type,"_", AnnotCellTypeColumn)]] <- ifelse(metadata[[metadata_col]] < thresholds, "F", "T")

  seuratObj_Query@meta.data <- metadata

  seuratObj_Query@meta.data <- seuratObj_Query@meta.data[, !colnames(seuratObj_Query@meta.data) %in% c("VICTOR_max", "VICTOR_prediction", "VICTOR_no_rejection")]
  seuratObj_Ref@meta.data <- seuratObj_Ref@meta.data[, !colnames(seuratObj_Ref@meta.data) %in% c("VICTOR_max", "VICTOR_prediction", "VICTOR_no_rejection")]


  #### Export ####
  Output.lt <- list(Query = seuratObj_Query, Reference = seuratObj_Ref)
  return(Output.lt)
}

############--------------------------------------------------------------------#############
