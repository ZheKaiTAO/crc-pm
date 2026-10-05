### ==========
#   Seurat2Py: Convert RData into .h5ad format
### ==========

# ----------
# DATA LIST == single cell transcriptomic analysis
# HRA003293: tumor
# GSE234804: GSM747xxxx_xx_new.h5seurat
# GSE183916: GSE183916_merged_objects.rds
# ----------

# ----------
# Install Packages
# ----------
renv::init()
CranPacList <- c("Seurat","tidyverse","BiocManager","remotes")
library(BiocManager)
BiocPacList <- c("anndataR","SingleCellExperiment","rhdf5")
for (pac in CranPacList){
  if (! require(pac,character.only=T)){
    install.packages(pac,ask=F)
  }
}
for (pac in BiocPacList){
  if (! require(pac,character.only=T)){
    BiocManager::install(pac,ask = F)
  }
}
library(remotes)
RemotePacList <- c('SeuratDisk','SeuratData')
remotes::install_github("mojaveazure/seurat-disk")
remotes::install_github("satijalab/seurat-data")


# verify installation
for (pac in c(CranPacList,BiocPacList,RemotePacList)){
  require(pac, character.only=T)
}

# call packages
for (pac in c(CranPacList,BiocPacList,RemotePacList)){
  library(pac, character.only=T)
}

# ----------
# Manually Load Data, then Transform to .h5ad file
# ----------

# HRA003293
write_h5ad(tumor, path = "./data/Li_2025_NatCa/HRA003293.h5ad",compression = 'gzip')

# GSE234804
##########
# if Seurat version > 5.0, there is a compatibility issue between Seurat and SeuratDisk
# please refer to https://github.com/mojaveazure/seurat-disk/issues/192
# the following is a monkey patch generated via GPT
ns <- asNamespace("SeuratObject")
old_get <- getS3method("GetAssayData","Assay",envir = ns)
old_set <- getS3method("SetAssayData","Assay",envir = ns)
compat_get <- function(object,layer = c("data", "scale.data", "counts"),slot,...) {
  if (!missing(slot)) {layer <- slot}
  old_get(object = object,layer = layer,...)
}
compat_set <- function(object,layer = "data",new.data,slot,...) {
  if (!missing(slot)) {layer <- slot}
  old_set(object = object,layer = layer,new.data = new.data,...)
}
registerS3method("GetAssayData","Assay",compat_get,envir = ns)
registerS3method("SetAssayData","Assay",compat_set,envir = ns)
##########

# only fetch crc and pm, cast away all hepatic metastasis

DataRootDir = './seq-data/GSE234804'
CRCPattern = "^GSM747\\d+_CRC\\d+_new\\.h5seurat$"
PCPattern = "^GSM747\\d+_PC\\d+_new\\.h5seurat$"
CRCFileName <- list.files(
  path = DataRootDir,
  pattern = CRCPattern,
  full.names = T
) 
PCFileName <- list.files(
  path = DataRootDir,
  pattern = PCPattern,
  full.names = T
) 
for (dir in c(CRCFileName,PCFileName)){
  dest = paste0('./data/',str_extract(dir,pattern = "crc.+h5" ),'ad')
  h5seurat <- LoadH5Seurat(dir, images = FALSE)
  write_h5ad(h5seurat, path = dest,compression = 'gzip')
  rm('h5seurat')
}

# GSE183916
write_h5ad(GSE183916_merged_objects, path = "./data/Lenos_2022_NatCom/GSE183916/GSE183916.h5ad",compression = 'gzip')