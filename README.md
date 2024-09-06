# VICTOR: Validation and Inspection of Cell Type Annotation through Optimal Regression
Chang, Chia-Jung, Chih-Yuan Hsu, Qi Liu, and Yu Shyr. "VICTOR: Validation and Inspection of Cell Type Annotation through Optimal Regression." Computational and Structural Biotechnology Journal (2024).
https://www.sciencedirect.com/science/article/pii/S2001037024002873


# Installation Instructions
To ensure that the necessary packages are installed and loaded into your R environment, follow these steps:

## Install and Load Devtools: 
Devtools is essential for installing packages directly from GitHub repositories, among other utilities. If it's not already installed, the script below will handle its installation and then load it for you.

```{r, eval = FALSE}
if(!require("devtools")) install.packages("devtools"); library(devtools)
```
## Install and Load VICTOR: 
VICTOR is a package designed for specific analyses. If it's not present in your R library, the following command will install it from GitHub and then load it. Ensure you have an internet connection for this step.
```{r, eval = FALSE}
if(!require("VICTOR")) install_github("Charlene717/VICTOR"); library(VICTOR)
```

<br>

# Operational Procedure
After installing the necessary packages, you will need to load additional packages, load a demo file, and perform operations with VICTOR as follows:

## Load Additional Packages:
For data manipulation and analysis, we'll be using tidyverse and Seurat. The script below checks for their presence and installs them if they're missing, before loading them into the session.

```{r, eval = FALSE}
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
if(!require("Seurat")) install.packages("Seurat"); library(Seurat)
```

## Load the Demo File:
Load your demo dataset by specifying the path to the .RData file. Replace Path/to/Demo.RData with the actual path to your file.
```{r, eval = FALSE}
load("Path/to/Demo.RData")
```

## Set Query and Reference Variables:
Configure the VICTOR function by specifying your query and reference Seurat objects, as well as the columns for the actual cell type and annotation.

```{r, eval = FALSE}
VICTOR.lt <- VICTOR(seuratObject_Query, seuratObject_Ref,
                      ActualCellTypeColumn = "Actual_Cell_Type",
                      AnnotCellTypeColumn = "Annotation")

```

## Update Query and Reference:
After running VICTOR, please update your Seurat objects as follows:

### Query Object: 
Results are in the query's meta.data. Update with
```{r, eval = FALSE}
seuratObject_Query <- VICTOR.lt$Query
```
### Reference Object: 
Model stored in the reference's misc. Update with 
```{r, eval = FALSE}
seuratObject_Ref <- VICTOR.lt$Reference
```



