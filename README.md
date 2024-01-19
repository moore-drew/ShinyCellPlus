# ShinyCellPlus
Expansion of ShinyCell: https://github.com/SGDDNB/ShinyCell

To install:

```
devtools::install_github("bioinformaticsMUSC/ShinyCellPlus")
```

Multiple tabs of ShinyCellPlus require specific data preparation with external libraries.  For full usage, please install Presto:

```
devtools::install_github("immunogenomics/presto")
```

Libra:

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("edgeR", "DESeq2", "limma"))

devtools::install_github("neurorestore/Libra")
```

AUCell:

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("AUCell")
```

In order to maintain Seurat meta.data column names and data types to go along with this tutorial, it is recommended to follow general R pipeline scripts, which are provided (https://github.com/BioinformaticsMUSC/biocm-scripts/blob/main/r_tools/interactive_scripts/, or https://github.com/BioinformaticsMUSC/biocm-scripts/blob/main/r_tools/auto_scripts/).

## Presto Markers and AUCell Gene Ranks

Presto markers data needs to be stored as such:

```
seurat@misc$markers$presto$all
seurat@misc$markers$presto$top_20
```

and AUCell gene ranks as:

```
seurat@misc$gene_ranks$aucell$all
```

The associated commands can be ran after standard processing and the removal of doublets:

```
# UPDATE THESE TO YOUR PERSONAL DATA DIRECTORIES, FILES, AND SELECTED RESOLUTION (EX: "integrated_snn_res.0.6")
opt = list(seurat_file_path = "~/example/file/path/seurat.rds",
           seurat_save_name = "example_output_seurat.rds",
	   seurat_output_path = "~/exmaple/file/output_path/",
	   resolution = "integrated_snn_res.0.6")

##############################
### LOAD/INSTALL LIBRARIES ###
##############################
cat("\nLoading/installing libraries ...\n")
suppressMessages(suppressWarnings({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(presto)
  library(AUCell)
}))

####################################
### LOAD DATA INTO SEURAT OBJECT ###
####################################
cat('\nLoading seurat object and converting to SCE...')
seurat <- readRDS(opt$seurat_file_path)


##############################
### RUN PRESTO ON CLUSTERS ###
##############################

seurat$seurat_clusters <- seurat[[opt$resolution]]

presto_markers <- wilcoxauc(seurat, "seurat_clusters")
top <- top_markers(presto_markers, n = 20)

seurat@misc$markers$presto$all <- presto_markers
seurat@misc$markers$presto$top_20 <- top

seurat@misc$gene_ranks$aucell$all <- AUCell_buildRankings(seurat@assays$SCT@counts)
```

Alternatively to Presto, we can use the Seurat library's FindAllMarkers() to calculate all markers and top 20:

```
library(tidyverse)
seurat@misc$markers$seurat$all <- FindAllMarkers(seurat)

auROC<-FindAllMarkers(seurat, test.use="roc")
seurat@misc$markers$seurat$top_20 <- auROC %>% group_by(.data$cluster) %>% top_n(n=20, wt=.data$myAUC) %>% mutate(rank=rank(-.data$myAUC, ties.method="random")) %>% ungroup() %>% select(.data$gene, .data$cluster, .data$rank) %>% spread(.data$cluster, .data$gene, fill=NA)
```

## Libra Differentially Expressed Genes

Preparing this can be a bit more complex, as depending on your data you can have varying amounts of differential expressions that you wish to compare.  The fields that you wish to differentiate may include two separate columns of meta.data.  For example, we will be using an example meta.data column 'Genotype', with set c("IC1", "IC2", "IC3", "IC4"):

```
genotypes = c("IC1", "IC2", "IC3", "IC4")
genotype_combs <- t(combn(genotypes, m=2)) |> as.data.frame()
```

We will define several functions to allow for the calculation of several sets of differential expressions and their accumulation into an overall data.frame along with a column defining which two genotypes were used for given set:

```
library(Libra)

# define function to run libra
run_libra <- function(seurat_obj, cell_type_col, replicate_col, label_col,
                      de_family, de_method, de_type){
  DefaultAssay(seurat_obj) <- "SCT"
  tmp <- seurat_obj
  
  tmp@meta.data$cell_type <- tmp@meta.data[[cell_type_col]]
  tmp@meta.data$replicate <- tmp@meta.data[[replicate_col]]
  tmp@meta.data$label <- tmp@meta.data[[label_col]]
  
  DE_output <- run_de(tmp, de_family = de_family, de_method = de_method,
                       de_type = de_type)
  
  return(DE_output) 
}

# define function to consolidate specific DE data tables into one 'overall'
add_libra_DE_table_to_seurat <- function(seurat, diff_exp_table_var, diff_exp_name) {

	if( !all(class(diff_exp_table_var) == c('tbl_df', 'tbl', 'data.frame')) ) {
		e<-simpleError("\"tbl_df, tbl, data.frame\" type expected (standard outcome class from Libra's run_de() function: https://github.com/neurorestore/Libra/blob/main/R/run_de.R).")
		stop(e)
	}

	diff_exp_table_var['de_name'] <- rep(diff_exp_name, times=nrow(diff_exp_table_var))

	if(is.null(seurat@misc$DE_genes$libra$overall)) {
		seurat@misc$DE_genes$libra$overall <- diff_exp_table_var
	}
	else {
		if( !all(class(seurat@misc$DE_genes$libra$overall) == class(diff_exp_table_var)) ) {
			e<-simpleError("seurat@misc$DE_genes$libra$overall class type is not the same as data, \"tbl_df, tbl, data.frame\" expected (standard outcome class from Libra's run_de() function: https://github.com/neurorestore/Libra/blob/main/R/run_de.R).")
			stop(e)
		}
		else{
			seurat@misc$DE_genes$libra$overall <- rbind(seurat@misc$DE_genes$libra$overall, diff_exp_table_var)
		}
	}

	return(seurat)
}

```

Now we will iterate through all defined combinations within 'genotype_combs', compute their differential expressions, combine them into a single data.frame called 'overall', and store them within the Seurat object (seurat@misc$DE_genes$libra$overall):

```
for (i in rownames(genotype_combs)) {
  title <- stringr::str_glue("{genotype_combs[i,1]}vs{genotype_combs[i,2]}")
  print(title)
  sub_seurat <- subset(seurat, Genotype %in% c(genotype_combs[i,1], genotype_combs[i,2]))
 
  #libra function on sub_seurat
  diff_exp <- run_libra(seurat_obj = sub_seurat,
                        cell_type_col = 'Cell',
                        replicate_col = 'Genotype',
                        label_col = 'Genotype',
                        de_family = 'singlecell',
                        de_method = 'wilcox',
                        de_type = NULL)

  seurat <- add_libra_DE_table_to_seurat(seurat, diff_exp, title)
}
```

## Creating the Shiny App from Prepared Seurat Object

IMPORTANT: It is recommended to save your Seurat object as this time, as sometimes interacting with local Shiny apps or closing them before their data is fully loaded and presented can cause RStudio to crash.

```
saveRDS(seurat, file=paste0(opt$seurat_output_path, '/SCP_prepped_', opt$seurat_save_name))
```

With a prepared Seurat object, the final script for Shiny App creation is as follows:

```
# UPDATE THESE TO YOUR PERSONAL DATA DIRECTORIES, FILES, AND SEURAT META.DATA COLUMNS
opt = list(output_dir = "~/tutorial/output/directory/",
           seurat_object = "~/tutorial/seurat.rds",
           app_name = "tutorial_app_name",
           scconf_defaults = c('Cell', 'seurat_clusters'),
           seurat_columns = "default")


##############################
### LOAD/INSTALL LIBRARIES ###
##############################
cat("\nLoading/installing libraries ...\n")
suppressMessages(suppressWarnings({
  library(Seurat)
  library(devtools)
  library(ShinyCellPlus)
  library(dplyr)
}))


##########################
### LOAD SEURAT OBJECT ###
##########################
seurat <- readRDS(opt$seurat_object)
setwd(opt$output_dir)


# filter the meta data
seurat@meta.data <- seurat@meta.data |>
  dplyr::select(Genotype, pMito, nCount_SCT, nFeature_SCT, seurat_clusters, Cell) # choose essential columns to visualize within the Shiny App's interactive graphs


# create app config
scConf = createConfig(seurat)

#modify defaults
scConf <- modDefault(scConf, default1 = opt$scconf_defaults[1],
                     default2 = opt$scconf_defaults[2])


makeShinyApp(seurat, scConf, gene.mapping = TRUE, shiny.title = opt$app_name,
 	     ### NEW, OPTIONAL PAGE ARGUMENTS; DEFAULT TO 'FALSE' ###
             markers.all=TRUE, markers.top20=TRUE, de.genes=TRUE,
             gene.ranks=TRUE, volc.plot=TRUE)
```

The Shiny App should run upon completion of this script.

## Deploying to shinyapps.io

ShinyCellPlus is written in Shiny and can be deployed to Posit's web service, https://shinyapps.io.  After setting up an account with shinyapps.io, copy your 'name'/'token'/'secret' pair through the "Tokens" page; it should be structured like this:

```
rsconnect::setAccountInfo(name='<ACCOUNT>',
			  token='<TOKEN>',
			  secret='<SECRET>')
```

Using this copied string, run the following script to deploy to shinyapps.io:

```
rsconnect::setAccountInfo(name='<ACCOUNT>',
                          token='<TOKEN>',
                          secret='<SECRET>')

options(repos = BiocManager::repositories())

rsconnect::deployApp("shinyApp/",
                     appName = opt$app_name,
                     account = '<ACCOUNT>', # REPLACE WITH YOUR ACCOUNT NAME
                     server = 'shinyapps.io')
```

Note that there are size restrictions depending on your account's subscription to the service, including the size of the data uploaded to shinyapps.io as well as runtime memory.  You can adjust the memory allocated to the site for larger applications through the shinyapps.io deployed applications settings (which should be a commmon  need with the data ShinyCellPlus is designed for).  You can also adjust the size of the data to be uploaded if you receive an associated error while deploying:

```
options(rsconnect.max.bundle.size=5368709000)
```

Which should adjust a global RStudio option to a few kilobytes below the maximum allowed shinyapps.io upload size (5 GiB).  At least a 'Basic' subscription plan is required for this adjustment.
