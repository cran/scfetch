## ----setup, echo=FALSE, warning=FALSE-----------------------------------------
library(knitr)
htmltools::tagList(rmarkdown::html_dependency_font_awesome())

# set dpi
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  dpi = 60
)

## ----install, eval=FALSE------------------------------------------------------
#  # install via CRAN (v0.5.0) # old version, it's better to install via Github
#  install.packages("scfetch")
#  # if you install from CRAN, you should install the following packages
#  # install.packages("devtools") #In case you have not installed it.
#  devtools::install_github("alexvpickering/GEOfastq") # download fastq
#  devtools::install_github("cellgeni/sceasy") # format conversion
#  devtools::install_github("mojaveazure/seurat-disk") # format conversion
#  devtools::install_github("satijalab/seurat-wrappers") # format conversion
#  
#  # install via Github (v0.5.0)
#  devtools::install_github("showteeth/scfetch")

## ----library, message=FALSE, warning=FALSE------------------------------------
library("scfetch")

## ----prepare_run, eval=FALSE--------------------------------------------------
#  GSE130636.runs <- ExtractRun(acce = "GSE130636", platform = "GPL20301")

## ----dwonload_sra, eval=FALSE-------------------------------------------------
#  # a small test
#  GSE130636.runs <- GSE130636.runs[GSE130636.runs$run %in% c("SRR9004346", "SRR9004351"), ]
#  # download, you may need to set prefetch.path
#  out.folder <- tempdir()
#  GSE130636.down <- DownloadSRA(
#    gsm.df = GSE130636.runs,
#    out.folder = out.folder
#  )
#  # GSE130636.down is null or dataframe contains failed runs

## ----split_sra, eval=FALSE----------------------------------------------------
#  # parallel-fastq-dump requires sratools.path
#  # you may need to set split.cmd.path and sratools.path
#  sra.folder <- tempdir()
#  GSE130636.split <- SplitSRA(
#    sra.folder = sra.folder,
#    fastq.type = "10x", split.cmd.threads = 4
#  )

## ----prepare_run_bam, eval=FALSE----------------------------------------------
#  GSE138266.runs <- ExtractRun(acce = "GSE138266", platform = "GPL18573")

## ----dwonload_bam, eval=FALSE-------------------------------------------------
#  # a small test
#  GSE138266.runs <- GSE138266.runs[GSE138266.runs$run %in% c("SRR10211566"), ]
#  # download, you may need to set prefetch.path
#  out.folder <- tempdir()
#  GSE138266.down <- DownloadBam(
#    gsm.df = GSE138266.runs,
#    out.folder = out.folder
#  )
#  # GSE138266.down is null or dataframe contains failed runs

## ----convert_bam_fastq, eval=FALSE--------------------------------------------
#  bam.folder <- tempdir()
#  # you may need to set bamtofastq.path and bamtofastq.paras
#  GSE138266.convert <- Bam2Fastq(
#    bam.folder = bam.folder
#  )

## ----geo_meta, eval=FALSE-----------------------------------------------------
#  # extract metadata of specified platform
#  GSE200257.meta <- ExtractGEOMeta(acce = "GSE200257", platform = "GPL24676")
#  # set VROOM_CONNECTION_SIZE to avoid error: Error: The size of the connection buffer (786432) was not large enough
#  Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 60)
#  # extract metadata of all platforms
#  GSE94820.meta <- ExtractGEOMeta(acce = "GSE94820", platform = NULL)

## ----geo_parse, eval=FALSE----------------------------------------------------
#  # for cellranger output
#  out.folder <- tempdir()
#  GSE200257.seu <- ParseGEO(
#    acce = "GSE200257", platform = NULL, supp.idx = 1, down.supp = TRUE, supp.type = "10x",
#    out.folder = out.folder
#  )
#  # for count matrix, no need to specify out.folder, download count matrix to tmp folder
#  GSE94820.seu <- ParseGEO(acce = "GSE94820", platform = NULL, supp.idx = 1, down.supp = TRUE, supp.type = "count")

## ----panglaodb_summary, eval=FALSE--------------------------------------------
#  # use cached metadata
#  StatDBAttribute(df = PanglaoDBMeta, filter = c("species", "protocol"), database = "PanglaoDB")

## ----panglaodb_meta, eval=FALSE-----------------------------------------------
#  hsa.meta <- ExtractPanglaoDBMeta(
#    species = "Homo sapiens", protocol = c("Smart-seq2", "10x chromium"),
#    show.cell.type = TRUE, cell.num = c(1000, 2000)
#  )

## ----panglaodb_celltype, eval=FALSE-------------------------------------------
#  hsa.composition <- ExtractPanglaoDBComposition(species = "Homo sapiens", protocol = c("Smart-seq2", "10x chromium"))

## ----panglaodb_parse, eval=FALSE----------------------------------------------
#  # small test
#  hsa.seu <- ParsePanglaoDB(hsa.meta[1:3, ], merge = TRUE)

## ----cb_show, eval=FALSE------------------------------------------------------
#  json.folder <- tempdir()
#  # first time run, the json files are stored under json.folder
#  # ucsc.cb.samples = ShowCBDatasets(lazy = TRUE, json.folder = json.folder, update = TRUE)
#  
#  # second time run, load the downloaded json files
#  ucsc.cb.samples <- ShowCBDatasets(lazy = TRUE, json.folder = json.folder, update = FALSE)
#  
#  # always read online
#  # ucsc.cb.samples = ShowCBDatasets(lazy = FALSE)

## ----cb_show_detail, eval=FALSE-----------------------------------------------
#  # the number of datasets
#  nrow(ucsc.cb.samples)
#  
#  # available species
#  unique(unlist(sapply(unique(gsub(pattern = "\\|parent", replacement = "", x = ucsc.cb.samples$organisms)), function(x) {
#    unlist(strsplit(x = x, split = ", "))
#  })))

## ----cb_summary, eval=FALSE---------------------------------------------------
#  StatDBAttribute(df = ucsc.cb.samples, filter = c("organism", "organ"), database = "UCSC")

## ----cb_extract, eval=FALSE---------------------------------------------------
#  hbb.sample.df <- ExtractCBDatasets(all.samples.df = ucsc.cb.samples, organ = c("brain", "blood"), organism = "Human (H. sapiens)", cell.num = c(1000, 2000))

## ----cb_celltype, eval=FALSE--------------------------------------------------
#  hbb.sample.ct <- ExtractCBComposition(json.folder = json.folder, sample.df = hbb.sample.df)

## ----cb_parse, eval=FALSE-----------------------------------------------------
#  hbb.sample.seu <- ParseCBDatasets(sample.df = hbb.sample.df)

## ----zenodo_meta, eval=FALSE--------------------------------------------------
#  # single doi
#  zebrafish.df <- ExtractZenodoMeta(doi = "10.5281/zenodo.7243603")
#  
#  # vector dois
#  multi.dois <- ExtractZenodoMeta(doi = c("1111", "10.5281/zenodo.7243603", "10.5281/zenodo.7244441"))

## ----zenodo_parse, eval=FALSE-------------------------------------------------
#  out.folder <- tempdir()
#  multi.dois.parse <- ParseZenodo(
#    doi = c("1111", "10.5281/zenodo.7243603", "10.5281/zenodo.7244441"),
#    file.ext = c("rdata", "rds"), out.folder = out.folder
#  )

## ----cellxgene_all, eval=FALSE------------------------------------------------
#  # all available datasets
#  all.cellxgene.datasets <- ShowCELLxGENEDatasets()

## ----cellxgene_summary, eval=FALSE--------------------------------------------
#  StatDBAttribute(df = all.cellxgene.datasets, filter = c("organism", "sex"), database = "CELLxGENE")

## ----cellxgene_meta, eval=FALSE-----------------------------------------------
#  # human 10x v2 and v3 datasets
#  human.10x.cellxgene.meta <- ExtractCELLxGENEMeta(
#    all.samples.df = all.cellxgene.datasets,
#    assay = c("10x 3' v2", "10x 3' v3"), organism = "Homo sapiens"
#  )

## ----cellxgene_parse, eval=FALSE----------------------------------------------
#  out.folder <- tempdir()
#  ParseCELLxGENE(
#    meta = human.10x.cellxgene.meta[1:5, ], file.ext = "rds",
#    out.folder = out.folder
#  )

## ----test_data, eval=FALSE----------------------------------------------------
#  # library
#  library(Seurat) # pbmc_small
#  library(scRNAseq) # seger

## ----test_seurat, eval=FALSE--------------------------------------------------
#  # object
#  pbmc_small

## ----testsce, eval=FALSE------------------------------------------------------
#  seger <- scRNAseq::SegerstolpePancreasData()

## ----seu2sce, eval=FALSE------------------------------------------------------
#  sce.obj <- ExportSeurat(seu.obj = pbmc_small, assay = "RNA", to = "SCE")

## ----seu2cds1, eval=FALSE-----------------------------------------------------
#  # BiocManager::install("monocle") # reuqire monocle
#  cds.obj <- ExportSeurat(seu.obj = pbmc_small, assay = "RNA", reduction = "tsne", to = "CellDataSet")

## ----seu2cds2, eval=FALSE-----------------------------------------------------
#  # remotes::install_github('cole-trapnell-lab/monocle3') # reuqire monocle3
#  cds3.obj <- ExportSeurat(seu.obj = pbmc_small, assay = "RNA", to = "cell_data_set")

## ----seu2anndata, eval=FALSE--------------------------------------------------
#  # remove pbmc_small.h5ad first
#  anndata.file <- tempfile(pattern = "pbmc_small_", fileext = ".h5ad")
#  # you may need to set conda.path
#  ExportSeurat(
#    seu.obj = pbmc_small, assay = "RNA", to = "AnnData",
#    anndata.file = anndata.file
#  )

## ----seu2loom, eval=FALSE-----------------------------------------------------
#  loom.file <- tempfile(pattern = "pbmc_small_", fileext = ".loom")
#  ExportSeurat(
#    seu.obj = pbmc_small, assay = "RNA", to = "loom",
#    loom.file = loom.file
#  )

## ----sce2seu, eval=FALSE------------------------------------------------------
#  seu.obj.sce <- ImportSeurat(obj = sce.obj, from = "SCE", count.assay = "counts", data.assay = "logcounts", assay = "RNA")

## ----cds2seu1, eval=FALSE-----------------------------------------------------
#  seu.obj.cds <- ImportSeurat(obj = cds.obj, from = "CellDataSet", count.assay = "counts", assay = "RNA")

## ----cds2seu2, eval=FALSE-----------------------------------------------------
#  seu.obj.cds3 <- ImportSeurat(obj = cds3.obj, from = "cell_data_set", count.assay = "counts", data.assay = "logcounts", assay = "RNA")

## ----anndata2seu, eval=FALSE--------------------------------------------------
#  # you may need to set conda.path
#  seu.obj.h5ad <- ImportSeurat(
#    anndata.file = anndata.file, from = "AnnData", assay = "RNA"
#  )

## ----loom2seu, eval=FALSE-----------------------------------------------------
#  # loom will lose reduction
#  seu.obj.loom <- ImportSeurat(loom.file = loom.file, from = "loom")

## ----sce2anndata, eval=FALSE--------------------------------------------------
#  # remove seger.h5ad first
#  seger.anndata.file <- tempfile(pattern = "seger_", fileext = ".h5ad")
#  SCEAnnData(
#    from = "SingleCellExperiment", to = "AnnData", sce = seger, X_name = "counts",
#    anndata.file = seger.anndata.file
#  )

## ----anndata2sce, eval=FALSE--------------------------------------------------
#  seger.anndata <- SCEAnnData(
#    from = "AnnData", to = "SingleCellExperiment",
#    anndata.file = seger.anndata.file
#  )

## ----sce2loom, eval=FALSE-----------------------------------------------------
#  # remove seger.loom first
#  seger.loom.file <- tempfile(pattern = "seger_", fileext = ".loom")
#  SCELoom(
#    from = "SingleCellExperiment", to = "loom", sce = seger,
#    loom.file = seger.loom.file
#  )

## ----loom2sce, eval=FALSE-----------------------------------------------------
#  seger.loom <- SCELoom(
#    from = "loom", to = "SingleCellExperiment",
#    loom.file = seger.loom.file
#  )

