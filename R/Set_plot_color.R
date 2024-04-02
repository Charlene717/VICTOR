# -[] 設定當Set_Reference和Set_Sample不在的情況
if(!exists("Set_Reference")) {Set_Reference <- ""}
if(!exists("Set_Sample")) {Set_Sample <- ""}

if(Set_Reference == "HumanPrimaryCellAtlasData"){


# }else if(Set_Reference == "pbmc3k"){
}else if(grepl("pbmc3k", Set_Reference)){
  ## Set color
  color_CellType_Ref <- list(
    "Memory CD4 T" = "#1f77b4",  # blue
    "B" = "#ff7f0e",             # orange
    "CD14+ Mono" = "#2ca02c",    # green
    "NK" = "#d62728",            # red
    "CD8 T" = "#9467bd",         # purple
    "Naive CD4 T" = "#8c564b",   # brown
    "FCGR3A+ Mono" = "#e377c2",  # pink
    "Unknown" = "#7f7f7f",       # grey
    "Unassign" = "#84e8d4",
    "pruned" = "#f0ad8d",  #"#4a1919",       # dark grey
    "DC" = "#bcbd22",            # olive
    "Platelet" = "#17becf"       # light blue
  )

}else if(Set_Reference == "ZilionisLungData"){
  ## Set color
  color_CellType_Ref <- list(
    "T_cells_CD4_memory_activated" = "#1f77b4",  # blue
    "T_cells_CD4_memory_resting" = "#0f598c",
    "T_cells_regulatory_(Tregs)" = "#3283ba",
    "T_cells_follicular_helper" = "#416985",
    "B_cells_memory" = "#ff7f0e",             # orange
    "B_cells_naive" = "#e67b1e",
    "Monocytes" = "#2ca02c",    # green
    "NK_cells_activated" = "#d62728",            # red
    "NK_cells_resting"  = "#b01718",
    "T_cells_CD8" = "#9467bd",         # purple
    "T_cells_CD4_naive" = "#8c564b",   # brown
    "Macrophages_M0" = "#e377c2",  # pink
    "Macrophages_M1" = "#a84f8d",
    "Macrophages_M2" = "#7a205f",
    "Eosinophils" = "#7bcbd4",
    "Mast_cells_activated" = "#3b99a3",
    "Mast_cells_resting" = "#337a82",
    "Neutrophils" = "#56d166",
    "Unknown" = "#7f7f7f",       # grey
    "Unassign" = "#84e8d4",
    "pruned" = "#f0ad8d",  #"#4a1919",       # dark grey
    "Dendritic_cells_activated" = "#bcbd22",            # olive
    "Dendritic_cells_resting" = "#7f8016",
    "Plasma_cells" = "#17becf"       # light blue
  )

}else{
  ## Set color
  color_CellType_Ref <- list(
    "CD4+ T cell" = "#1f77b4",  # blue
    "B cell" = "#ff7f0e",             # orange
    "CD14+ monocyte" = "#2ca02c",    # green
    "Natural killer cell" = "#d62728",            # red
    "Cytotoxic T cell" = "#9467bd",         # purple
    "Megakaryocyte" = "#8c564b",   # brown
    "CD16+ monocyte" = "#e377c2",  # pink
    "Unknown" = "#7f7f7f",       # grey
    "unknown" = "#7f7f7f",       # grey
    "None" ="#7f7f7f",
    "Unassigned" = "#84e8d4",
    "Unassign" = "#84e8d4",
    "pruned" = "#f0ad8d",  #"#4a1919",       # dark grey
    "Dendritic cell" = "#bcbd22",            # olive
    "Plasmacytoid dendritic cell" = "#17becf" ,      # light blue


    ## PRJCA001063
    "Acinar cell" = "#e377c2",
    "Ductal cell type 1" = "#9467bd",
    "Ductal cell type 2" = "#62348c",
    "Endocrine cell" = "#34b399",
    "Endothelial cell" = "#83d663",
    "Fibroblast cell" = "#478047",
    "Macrophage cell" = "#8c564b",
    "Stellate cell" = "#5bb6f5",
    "T cell" = "#3f6b8a",

    ## BaronPanData
    # "Acinar cell" = "#e377c2",
    "Beta cell" = "#a64187",
    "Delta cell" = "#b37b84",
    "Stellate cell" = "",
    "Ductal cell" = "#9467bd",
    "Alpha cell" = "#17becf",
    "Epsilon cell" = "#609499",
    "Gamma cell" = "#bcbd22",
    # "Endothelial cell" = "#83d663",
    # "Stellate cell" = "#5bb6f5",
    # "Macrophage cell" = "#8c564b",
    "Schwann cell" = "#b05d7d",
    "Mast cell" = "#ad8832",
    # "T cell" = "#3f6b8a",

    ## MuraroPanData
    # "Alpha cell" = "#17becf",
    # "Endothelial cell" = "#83d663",
    # "Delta cell" = "#b37b84",
    # "Beta cell" = "#a64187",
    # "Unknown" = "#7f7f7f",       # grey
    # "Ductal cell" = "#9467bd",
    # "Acinar cell" = "#e377c2",
    "Pancreatic polypeptide cell" = "#1f585e",
    "Mesenchymal cell" = "#157f8a",
    # "Epsilon cell" = "#609499",

    ## SegerstolpePanData
    # "Gamma cell" = "#bcbd22",
    # "Alpha cell" = "#17becf",
    # "Beta cell" = "#a64187",
    # "Acinar cell" = "#e377c2",
    # "Epsilon cell" = "#609499",
    # "Unknown" = "#7f7f7f",       # grey
    # "Ductal cell" = "#9467bd",
    # "Delta cell" = "#b37b84",
    # "Endocrine cell" = "#34b399",
    # "Endothelial cell" = "#83d663",
    "co-expression cell" = "#663eab",
    # "Stellate cell" = "#5bb6f5",
    # "Mast cell" = "#ad8832",
    "MHC class II cell" = "#a87620",  ##


    ## PMID37027478
    "Type I pneumocyte" = "#9467bd",
    "Type II pneumocyte"  = "#62348c",
    # "B cell" = "#ff7f0e",
    "Natural killer T cell" = "#d95d5e",
    # "Mast cell" = "#ad8832",
    "Myeloid cell" = "#e377c2",
    "Proliferative cell"= "#8c564b",
    "Fibroblast" = "#478047",
    # "Endothelial cell" = "#83d663",
    "Epithelial cell" = "#83d663",
    "Lymphoid endothelial cell" = "#34b399",
    "Smooth muscle Cell" = "#ba5d8d",
    "Mito_high" =  "#84e8d4",


    "A549" = "#e377c2",
    "H1975" = "#62348c",
    "H2228" = "#83d663",
    "H838" = "#8c564b",
    "HCC827" = "#5bb6f5"
  )
}





################################################################################
## Set color for Sample
if(Set_Sample == "HumanPrimaryCellAtlasData"){

# }else if(Set_Sample == "pbmc3k"){
}else if(grepl("pbmc3k", Set_Sample)){
  color_CellType <- list(
    "Memory CD4 T" = "#1f77b4",  # blue
    "B" = "#ff7f0e",             # orange
    "CD14+ Mono" = "#2ca02c",    # green
    "NK" = "#d62728",            # red
    "CD8 T" = "#9467bd",         # purple
    "Naive CD4 T" = "#8c564b",   # brown
    "FCGR3A+ Mono" = "#e377c2",  # pink
    "Unknown" = "#7f7f7f",       # grey
    "Unassign" = "#84e8d4",
    "pruned" = "#f0ad8d",  #"#4a1919",       # dark grey
    "DC" = "#bcbd22",            # olive
    "Platelet" = "#17becf"       # light blue
  )

}else{
  ## Set color
  color_CellType <- list(
    "CD4+ T cell" = "#1f77b4",  # blue
    "B cell" = "#ff7f0e",             # orange
    "CD14+ monocyte" = "#2ca02c",    # green
    "Natural killer cell" = "#d62728",            # red
    "Cytotoxic T cell" = "#9467bd",         # purple
    "Megakaryocyte" = "#8c564b",   # brown
    "CD16+ monocyte" = "#e377c2",  # pink
    "Unknown" = "#7f7f7f",       # grey
    "unknown" = "#7f7f7f",       # grey
    "None" ="#7f7f7f",
    "Unassigned" = "#84e8d4",
    "Unassign" = "#84e8d4",
    "pruned" = "#f0ad8d",  #"#4a1919",       # dark grey
    "Dendritic cell" = "#bcbd22",            # olive
    "Plasmacytoid dendritic cell" = "#17becf" ,      # light blue

    ## PRJCA001063
    "Acinar cell" = "#e377c2",
    "Ductal cell type 1" = "#9467bd",
    "Ductal cell type 2" = "#62348c",
    "Endocrine cell" = "#34b399",
    "Endothelial cell" = "#83d663",
    "Fibroblast cell" = "#478047",
    "Macrophage cell" = "#8c564b",
    "Stellate cell" = "#5bb6f5",
    "T cell" = "#3f6b8a",

    ## BaronPanData
    # "Acinar cell" = "#e377c2",
    "Beta cell" = "#a64187",
    "Delta cell" = "#b37b84",
    "Stellate cell" = "",
    "Ductal cell" = "#9467bd",
    "Alpha cell" = "#17becf",
    "Epsilon cell" = "#609499",
    "Gamma cell" = "#bcbd22",
    # "Endothelial cell" = "#83d663",
    # "Stellate cell" = "#5bb6f5",
    # "Macrophage cell" = "#8c564b",
    "Schwann cell" = "#b05d7d",
    "Mast cell" = "#ad8832",
    # "T cell" = "#3f6b8a",

    ## MuraroPanData
    # "Alpha cell" = "#17becf",
    # "Endothelial cell" = "#83d663",
    # "Delta cell" = "#b37b84",
    # "Beta cell" = "#a64187",
    # "Unknown" = "#7f7f7f",       # grey
    # "Ductal cell" = "#9467bd",
    # "Acinar cell" = "#e377c2",
    "Pancreatic polypeptide cell" = "#1f585e",
    "Mesenchymal cell" = "#157f8a",
    # "Epsilon cell" = "#609499",

    ## SegerstolpePanData
    # "Gamma cell" = "#bcbd22",
    # "Alpha cell" = "#17becf",
    # "Beta cell" = "#a64187",
    # "Acinar cell" = "#e377c2",
    # "Epsilon cell" = "#609499",
    # "Unknown" = "#7f7f7f",       # grey
    # "Ductal cell" = "#9467bd",
    # "Delta cell" = "#b37b84",
    # "Endocrine cell" = "#34b399",
    # "Endothelial cell" = "#83d663",
    "co-expression cell" = "#663eab",
    # "Stellate cell" = "#5bb6f5",
    # "Mast cell" = "#ad8832",
    "MHC class II cell" = "#a87620",  ##


    ## PMID37027478
    "Type I pneumocyte" = "#9467bd",
    "Type II pneumocyte"  = "#62348c",
    # "B cell" = "#ff7f0e",
    "Natural killer T cell" = "#d95d5e",
    # "Mast cell" = "#ad8832",
    "Myeloid cell" = "#e377c2",
    "Proliferative cell"= "#8c564b",
    "Fibroblast" = "#478047",
    # "Endothelial cell" = "#83d663",
    "Epithelial cell" = "#83d663",
    "Lymphoid endothelial cell" = "#34b399",
    "Smooth muscle Cell" = "#ba5d8d",
    "Mito_high" =  "#84e8d4",

    "A549" = "#e377c2",
    "H1975" = "#62348c",
    "H2228" = "#83d663",
    "H838" = "#8c564b",
    "HCC827" = "#5bb6f5"

  )
}

################################################################################
## Set color for Classification
color_Class <- list(
  "Class 1" = "#b58b2a", # "#2c7856",
  "TP" = "#b58b2a", # "#2c7856",
  "Class 2" = "#e8bc56", # "#13c274",
  "TN" = "#e8bc56", # "#13c274",
  "Class 3" = "#368a5b", # "#8f174d",
  "FN" = "#368a5b", # "#8f174d",
  "Class 4" = "#73bd94", # "#696866" # "#d94185"
  "FP" = "#73bd94", # "#696866" # "#d94185"
  "Other" = "#919191",
  "NA" = "#919191"
)


color_State <- list(
  "T" = "#798ebd", # "#2c7856",
  "F" = "#eb8fb4", # "#2c7856",
  "NA" = "#9e9e9e"
)


## Set color for platform
color_platform <-list(
  "GSE132044_10xV2" = "#8AC0C0",
  "GSE132044_10xV2A" = "#AED4D4", # "#13c274",
  "GSE132044_10xV2B" = "#D2E8E8", # "#8f174d",
  "GSE132044_10xV3" = "#7ED6DF", # "#696866" # "#d94185"
  "GSE132044_CELSeq2" = "#FFF1E0", # "#2c7856",
  "GSE132044_DropSeq" = "#F7A072", # "#13c274",
  "GSE132044_inDrops" = "#FFC48C", # "#8f174d",
  "GSE132044_SeqWell" = "#D2B4DE", # "#696866" # "#d94185"
  "GSE132044_Smartseq2" = "#FFCF48" # "#8f174d",
)
