if(!require("pbapply")) install.packages("pbapply"); library(pbapply)
if(!require("caret")) install.packages("caret"); library(caret)

if (!requireNamespace("remotes", quietly = TRUE)) {install.packages("remotes")}
if(!require("harmony")) remotes::install_github("immunogenomics/harmony"); library(harmony)

#################################################################################
## ClassR
setClass("VICTOR", representation(pvar = "character",
                                  metadata = "data.frame",
                                  features = "list",
                                  cell_embeddings = "matrix",
                                  feature_loadings = "matrix",
                                  scaling = "data.frame",
                                  reduction = "character",
                                  reduction_key = "character",
                                  train = "list",
                                  misc = "list"),
         prototype(pvar = character(),
                   metadata = data.frame(),
                   features = list(),
                   cell_embeddings = matrix(),
                   feature_loadings = matrix(),
                   scaling = data.frame(),
                   reduction = character(),
                   reduction_key = character(),
                   train = list(),
                   misc = list())
)

#################################################################################
## generics.R
setGeneric("get_VICTOR", def = function(object) {
  standardGeneric("get_VICTOR")
})

setGeneric("get_metadata", def = function(object) {
  standardGeneric("get_metadata")
})


################################################################################
## getFeatureSpace.R

getFeatureSpace <- function(object, pvar, correction = "fdr", sig = 1, reduction = "pca"){


  # Validations -------------------------------------------------------------

  if(!is(object, "Seurat")){
    stop("Invalid class for object: must be 'Seurat'")
  }

  if(!any(correction %in% stats::p.adjust.methods)){
    stop("Invalid multiple testing correction method. See ?p.adjust function")
  }

  if(!pvar %in% names(object[[]])){
    stop("Prediction variable is not stored in metadata slot")
  }

  classes <- object[[pvar, drop = TRUE]]

  if(!is.factor(classes)){
    classes <- as.factor(classes)
  }

  assay <- DefaultAssay(object)
  # Filter principal components by variance ---------------------------------

  # Check if a PCA has been computed
  if(!(reduction %in% Reductions(object))){
    stop("No ",reduction, " reduction has been computet yet. See RunPCA() function?")
  }else{
    reduction_data <- Reductions(object, slot = reduction)
    if(reduction_data@assay.used != assay)
      stop("No ",
           reduction,
           " reduction associated with default assay: ",
           assay, "\nChange default assay or compute a new reduction")
  }

  # Check if available was normalized

  cellEmbeddings <- Embeddings(reduction_data)
  loadings <- Loadings(reduction_data)
  reduction_key <- reduction_data@key

  # Store original labels in metadata slot

  spmodel <- new("VICTOR", metadata = data.frame(pvar = classes))

  # Validate response variable values
  original_classes <- classes
  uniqueClasses <- unique(classes)
  isValidName <- uniqueClasses == make.names(uniqueClasses)

  if(!all(isValidName)){
    classes <- .make_names(classes)
    classes <- factor(classes, levels = unique(classes))
    names(classes) <- Cells(object)
  }

  spmodel@metadata$response <- classes

  # Get means and sds -------------------------------------------------------

  features <- rownames(loadings)

  data <- GetAssayData(object, "data", assay = assay)[features,]
  means <- Matrix::rowMeans(data)

  rowVar <- function(x, ...) {
    sqrt(Matrix::rowSums((x - means)^2, ...)/(ncol(x) - 1))
  }

  stdevs  <- rowVar(data)

  i <- stdevs == 0

  if(any(i)){
    warning(paste0(sum(i), " genes have zero variance but are present in the gene loadings. \nDid you subset or integrated this data before?"))
    cat(crayon::yellow("Removing zero-variance genes from loadings\n"))

    loadings <- loadings[!i,]
    means <- means[!i]
    stdevs <- stdevs[!i]
  }


  spmodel@scaling <- data.frame(means, stdevs)


  # Select informative principal components
  # If only 2 classes are present in prediction variable, train one model for the positive class
  # The positive class will be the first level of the factor variable

  cat(crayon::green(cli::symbol$record, " Extracting feature space for each cell type...\n"))
  if(length(levels(classes)) == 2){

    message("First factor level in '", pvar, "' metadata column considered as positive class:")
    message(levels(original_classes)[1])
    res <- .getFeatures(.make_names(levels(original_classes)[1]), classes, cellEmbeddings, correction, sig)
    res <- list(res)
    names(res) <- levels(original_classes)[1]

  }else{

    res <- pblapply(levels(classes), .getFeatures, classes, cellEmbeddings, correction, sig)
    dict <- data.frame(classes, original_classes)
    dict <- unique(dict)

    i <- match(levels(classes), dict$classes)
    names(res) <- as.character(dict$original_classes[i])

  }


  nFeatures <- unlist(lapply(res, nrow))

  noFeatures <- nFeatures == 0

  if(any(noFeatures)){

    warning("\nWarning: No features were found for classes:\n",
            paste0(names(res)[noFeatures], collapse = "\n"), "\n")
    res[[names(res)[noFeatures]]] <- NULL

  }

  # Create VICTOR object
  spmodel@pvar <- pvar
  spmodel@features <- res
  spmodel@cell_embeddings <- cellEmbeddings
  spmodel@feature_loadings <- loadings
  spmodel@reduction <- reduction
  spmodel@reduction_key <- reduction_key

  object@misc$VICTOR <- spmodel

  cat(crayon::green("DONE!\n"))

  object


}

.getFeatures <- function(positiveClass, classes, cellEmbeddings, correction, sig){

  # Set non-positive classes to "other"
  i <- classes != positiveClass
  newClasses <- as.character(classes)
  newClasses[i] <- "other"
  newClasses <- factor(newClasses, levels = c(positiveClass, "other"))

  # Get indices for positive and negative class cells
  positiveCells <- newClasses == positiveClass
  negativeCells <- newClasses == "other"


  # Get informative features
  wt_res <- apply(cellEmbeddings, 2, function(d) stats::wilcox.test(d[positiveCells], d[negativeCells]))
  wt_res <- lapply(wt_res, '[[', "p.value")
  wt_res <- data.frame(feature = names(wt_res), pValue = as.numeric(wt_res))
  wt_res$pValueAdj <- stats::p.adjust(wt_res$pValue, method = correction, n = nrow(wt_res))
  wt_res <- wt_res[order(wt_res$pValueAdj), ]
  wt_res <- wt_res[wt_res$pValueAdj < sig, ]
  wt_res
}


################################################################################
## methods.R
setMethod("get_VICTOR",
          signature("Seurat"),
          function(object){
            VICTOR <- object@misc$VICTOR
            if(is.null(VICTOR))
              stop("No VICTOR model is stored in this Seurat object")
            else
              VICTOR
          })


setMethod("get_metadata",
          signature("VICTOR"),
          function(object){
            object@metadata
          })


setMethod("get_metadata",
          signature("Seurat"),
          function(object){
            get_metadata(get_VICTOR(object))
          })

################################################################################
## project_query.R
project_query <- function(new,
                          reference,
                          max.iter.harmony = 20,
                          recompute_alignment = TRUE,
                          seed = 66,
                          ...){


  # Validate if provided object is an VICTOR object
  if(!(is(reference, "Seurat") | is(reference, "VICTOR"))) stop("'object' must be of class 'VICTOR' or 'Seurat'")

  if(is(reference, "Seurat")){
    spmodel <- reference@misc$VICTOR
  }else{
    spmodel <- reference
  }

  if(is.null(spmodel)) stop("No feature space has been determined!")
  if(!is(new, "Seurat")) stop("New data must be a Seurat object")


  # Dataset alignment -------------------------------------------------------

  if("VICTOR" %in% names(new@reductions)){
    if(recompute_alignment){
      alignment <- TRUE
      cat(crayon::yellow(cli::symbol$figure_dash, "Data has already been aligned to a reference.\n"), sep = "")
      cat(crayon::yellow(cli::symbol$sup_plus, "Skip data alignment using `recompute.alignment = FALSE`.\n"),  sep = "")

    }
    else {
      alignment <- FALSE
    }

  }else{
    alignment <- TRUE
  }

  if(alignment){

    cat(crayon::green(cli::symbol$record, " Matching reference with new dataset...\n"))

    # Subset data
    ref_loadings <- spmodel@feature_loadings
    ref_embeddings <- spmodel@cell_embeddings
    new_features <- rownames(new)

    # Get genes
    reference_features <- rownames(ref_loadings)


    # Get intersection between reference and new datasets
    shared_features <- intersect(reference_features, new_features)
    cat(crayon::cyan("\t", cli::symbol$line, paste(length(reference_features), "features present in reference loadings\n")))
    cat(crayon::cyan("\t", cli::symbol$line, paste(length(shared_features), "features shared between reference and new dataset\n")))
    cat(crayon::cyan("\t", cli::symbol$line, paste0(round(length(shared_features)/length(reference_features) * 100, 2),
                                                    "% of features in the reference are present in new dataset\n")))


    # Subset shared genes from reference
    ref_loadings <- ref_loadings[shared_features, ]



    new_data <- GetAssayData(new, "data")[shared_features,]
    means <- spmodel@scaling$means
    stdevs  <- spmodel@scaling$stdevs
    new_data <- Matrix::t(new_data)
    names(means) <- rownames(spmodel@scaling) -> names(stdevs)

    # Subset means and standard deviations
    means <- means[shared_features]
    stdevs <- stdevs[shared_features]
    #all(colnames(new_data) == names(means))

    i <- stdevs == 0

    if(any(i)){
      warning(paste0(sum(i), " features have zero variance but are present in the feature loadings. \nDid you subset or integrated this data before?"))
      cat(crayon::yellow("Removing zero-variance genes from projection\n"))

      new_data <- new_data[,!i]
      ref_loadings <- ref_loadings[!i,]
      means <- means[!i]
      stdevs <- stdevs[!i]
      #all(colnames(new_data) == rownames(ref_loadings))
    }


    scaled_data <- scale(new_data, means, stdevs)


    new_embeddings <- scaled_data %*% ref_loadings



    dataset <- factor(c(rep("reference", nrow(ref_embeddings)), rep("new", nrow(new_embeddings))),
                      levels = c("reference", "new"))


    rownames(ref_embeddings) <- paste0("ref_", rownames(ref_embeddings))
    rownames(new_embeddings) <- paste0("new_", rownames(new_embeddings))


    eigenspace <- as.data.frame(rbind(ref_embeddings, new_embeddings))
    meta_data <- data.frame(rownames(eigenspace), dataset = dataset)

    cat(crayon::green(cli::symbol$record, " Aligning new data to reference...\n"))

    set.seed(seed)
    harmony_embeddings <- HarmonyMatrix(eigenspace,
                                        meta_data,
                                        'dataset',
                                        do_pca = FALSE,
                                        reference_values = "reference",
                                        max.iter.harmony = max.iter.harmony,
                                        ...)

    new_embeddings_aligned <- harmony_embeddings[dataset == "new", , drop = FALSE]

  }else{
    new_embeddings_aligned <- Embeddings(new, reduction = "VICTOR")
    colnames(new_embeddings_aligned) <- gsub("VICTOR_", spmodel@reduction_key, colnames(new_embeddings_aligned))
  }


  rownames(new_embeddings_aligned) <- gsub("^new_", "", rownames(new_embeddings_aligned))
  new@reductions[["VICTOR"]] <- CreateDimReducObject(embeddings = new_embeddings_aligned,
                                                     key = "VICTOR_",
                                                     assay = DefaultAssay(object = new))
  if(recompute_alignment){

    rownames(new_embeddings) <- gsub("^new_", "", rownames(new_embeddings))
    new@reductions[["VICTOR_projection"]] <- CreateDimReducObject(embeddings = new_embeddings,
                                                                  key = "Projection_",
                                                                  assay = DefaultAssay(object = new))
  }

  new

}

################################################################################
## VICTOR.R
VICTORPrep <- function(new,
                    reference,
                    threshold = 0.55,
                    max.iter.harmony = 20,
                    recompute_alignment = TRUE,
                    seed = 66){

  # Function validations ----------------------------------------------------

  # Validate if provided object is an VICTOR object
  if(!(is(reference, "Seurat") | is(reference, "VICTOR"))) stop("'object' must be of class 'VICTOR' or 'Seurat'")

  if(is(reference, "Seurat")){
    spmodel <- reference@misc$VICTOR
  }else{
    spmodel <- reference
  }

  if(is.null(spmodel)) stop("No feature space has been determined!")
  if(!length(spmodel@train)) stop("No models have been trained!")
  if(!is(new, "Seurat")) stop("New data must be a Seurat object")


  # Project query data ------------------------------------------------------

  new <- project_query(new,
                       reference = spmodel,
                       max.iter.harmony = max.iter.harmony,
                       recompute_alignment = recompute_alignment,
                       seed = seed)

  new_embeddings_aligned <- Embeddings(new[["VICTOR"]])
  colnames(new_embeddings_aligned) <- colnames(spmodel@cell_embeddings)



  # Classify cells using all trained models
  cellTypeModelNames <- names(spmodel@features)
  .predictCellClass <-  function(cellType, spmodel, testEmbeddings){

    # Extract features for a given cell type
    as.character(spmodel@features[[cellType]]$feature) -> features

    # Extract cell type model
    model <- spmodel@train[[cellType]]

    # Perform predictions based on informative PCs
    prediction <- predict(model,
                          newdata = subsetMatrix(testEmbeddings, features),
                          type = "prob")

    # Add cell names to results
    rownames(prediction) <- rownames(testEmbeddings)

    # Return positive-class probability
    prediction[,1, drop = FALSE]

  }


  cat(crayon::green(cli::symbol$record, " Classifying cells...\n"))
  res <- sapply(cellTypeModelNames, .predictCellClass, spmodel, new_embeddings_aligned)

  # Gather results
  res <- as.data.frame(res)
  colnames(res) <- cellTypeModelNames
  rownames(res) <- colnames(new)

  classes <- cellTypeModelNames
  #plot(res$Lymphoid, col = as.factor(test$CellType))
  # If there is only 2 classes, compute complementary probability for negative class
  if(length(cellTypeModelNames) == 1){
    metadata <- get_metadata(spmodel)
    cellClasses <- levels(metadata$pvar)
    res_comp <- 1 - res[,1]
    negClass <- cellClasses[cellClasses != names(res)]
    res[[negClass]] <- res_comp

  }

  # Extract maximum probability for each class
  max_props <- as.data.frame(t(apply(res, 1, function(x) c(index = which.max(x),  max = x[which.max(x)]))))
  names(max_props) <- c("index", "max")


  # Store classification based on maximum probability
  max_props$generic_class <- names(res)[max_props$index]
  res <- cbind(res, max_props)

  # Classify cells according to probability threshold

  pred <- ifelse(res$max > threshold, res$generic_class, "unassigned")

  names(pred) <- colnames(new)

  # Format results
  res$prediction <- pred
  res$index <- NULL
  res$no_rejection <- res$generic_class
  res$generic_class <- NULL

  names(res) <- .make_names(paste0("VICTOR_", names(res)))


  # Return results
  new <- AddMetaData(new, res)

  cat(crayon::green("DONE!\n"))

  new


}

################################################################################
## trainModel.R
trainModel <- function(object,
                       model = c("glmnet","svmRadial"),
                       preProcess = c("center", "scale"),
                       resampleMethod = "cv",
                       number = 5,
                       seed = 66,
                       tuneLength = 3,
                       metric = c("ROC", "PR", "Accuracy", "Kappa"),
                       returnData = FALSE,
                       savePredictions = "final",
                       allowParallel = FALSE,
                       reclassify = NULL
){


  # Validations -------------------------------------------------------------

  # Check class
  if(!is(object, "Seurat") | is(object, "VICTOR")){
    stop("object must be 'Seurat' or 'VICTOR'")
  }


  if(is(object, "Seurat")){
    seurat_object <- object
    object <- get_VICTOR(object)

    if(is.null(object))
      stop("No features have been determined. Use 'getFeatureSpace()' function")

    object_class <- "Seurat"

  }else{
    object_class <- "VICTOR"
  }


  if(is.null(reclassify)){
    classes <- names(object@features)
  }else{
    classes <- reclassify
  }
  metric <- match.arg(metric)
  reduction <- object@reduction

  # Train a prediction model for each class
  cat(crayon::green(cli::symbol$record, " Training models for each cell type...\n"))



  if(length(classes) == 1){
    modelsRes <-  .trainModel(classes[1],
                              object,
                              model,
                              reduction,
                              preProcess,
                              resampleMethod,
                              tuneLength,
                              seed,
                              metric,
                              number,
                              returnData,
                              savePredictions,
                              allowParallel)
    modelsRes <- list(modelsRes)
    names(modelsRes) <- classes[1]


  }else{
    modelsRes <- pblapply(classes, .trainModel,
                          object,
                          model,
                          reduction,
                          preProcess,
                          resampleMethod,
                          tuneLength,
                          seed,
                          metric,
                          number,
                          returnData,
                          savePredictions,
                          allowParallel)
    names(modelsRes) <- classes
  }

  cat(crayon::green("DONE!\n"))

  if(is.null(reclassify)){
    object@train <- modelsRes
  }else{
    object@train[names(modelsRes)] <- modelsRes
  }

  if(object_class == "Seurat"){
    seurat_object@misc$VICTOR <- object
    seurat_object

  }else{
    object
  }
}


.trainModel <- function(positiveClass,
                        spmodel,
                        model,
                        reduction,
                        preProcess,
                        resampleMethod,
                        tuneLength,
                        seed,
                        metric,
                        number,
                        returnData,
                        savePredictions,
                        allowParallel){


  if(nrow(spmodel@features[[positiveClass]]) == 0){
    message("No informative principal components were identified for class: ", positiveClass)
  }

  names_features <- as.character(spmodel@features[[positiveClass]]$feature)
  features <- subsetMatrix(spmodel@cell_embeddings, names_features)
  response <- as.character(spmodel@metadata$response)


  i <- response != .make_names(positiveClass)
  response[i] <- "other"
  response <- factor(response, levels = c(.make_names(positiveClass), "other"))


  if(!is.null(seed)) set.seed(seed)

  if(metric == "ROC"){
    trCtrl <- trainControl(classProbs = TRUE,
                           method = resampleMethod,
                           number = number,
                           summaryFunction = twoClassSummary,
                           returnData = returnData,
                           savePredictions = savePredictions,
                           allowParallel = allowParallel)

  }else if(metric == "PR"){
    trCtrl <- trainControl(classProbs = TRUE,
                           method = resampleMethod,
                           number = number,
                           summaryFunction = prSummary,
                           returnData = returnData,
                           savePredictions = savePredictions,
                           allowParallel = allowParallel)
    # metric <- "AUC"
  }else{
    trCtrl <- trainControl(classProbs = TRUE,
                           method = resampleMethod,
                           number = number,
                           returnData = returnData,
                           savePredictions = savePredictions,
                           allowParallel = allowParallel)
  }

  ## ***
  # lassoGrid <- expand.grid(alpha = 0.5) #error: : The tuning parameter grid should have columns alpha, lambda#
  lassoGrid <- expand.grid(alpha = 0.5,
                           lambda = 10^seq(-6, 6, length=100))
  # cv <- cv.glmnet(x=X, y=y, family = "binomial", alpha = 0.5) # alpha = 0.5 不要精簡Var  ## 當係數值只有截距項則減少alpha
  # cv$lambda.min; cv$lambda.1se

  if(model== "glmnet"){
    fit <- train(x = features,
                 y = response,
                 method = model,
                 family = "binomial",
                 preProcess = preProcess,
                 metric = metric,
                 trControl = trCtrl,
                 tuneGrid = lassoGrid)


  }else if(model== "svmRadial"){
    fit <- train(x = features,
                 y = response,
                 method = model,
                 preProcess = preProcess,
                 metric = metric,
                 trControl = trCtrl)
  }


  fit
}

################################################################################
## utils.R
.make_names <- function(x){
  x <- gsub("\\+", "_plus", x)
  x <- gsub("\\-", "_minus", x)
  x <- gsub("\\ ", "_", x)
  x <- make.names(x)
}

subsetMatrix <- function(x, s, by.col = TRUE, drop = FALSE, verbose = FALSE, ...){

  if(by.col){
    ids <- colnames(x, ...)
    if(is.null(ids) & verbose) stop("No colnames were found!")
    i <- ids %in% s
    if(!any(i) & verbose){message("No matches were found")}
    x <- x[,i, drop = drop]
  }else{
    ids <- rownames(x, ...)
    if(is.null(ids) & verbose) stop("No rownames were found!")
    i <- ids %in% s
    if(!any(i) & verbose){message("No matches were found")}
    x <- x[i, , drop = drop]
  }

  x

}
