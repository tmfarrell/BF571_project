#
# utils.R
# 
# Tim Farrell, tmf@bu.edu 
# Project, BF571, BU
# 20160415
# 
library(Biobase)
library(GEOquery)

basedir = '/Volumes/Data/work/courses/bf571/project/'

## Functions 
swap <- function(ss) {
  return(paste(ss[[1]], ss[[3]], ss[[2]]))
}

get_gsmlist <- function(gse_str) { 
  setwd(paste0(basedir, 'data/')) 
  gse = getGEO(gse_str, GSEMatrix=FALSE)
  gsm_list = GSMList(gse)
  return(gsm_list)
}

get_matched_gsms <- function(gsm_list, pattern) { 
  matched = c() 
  for (i in 1:length(gsm_list)) { 
    if (length(grep(pattern, gsm_list[[i]]@header$title)))
      matched = c(matched, gsm_list[[i]])
  }
  return(matched)
}

get_gene_expr <- function(gsm, prot_ref, id_ref, include_id_col=FALSE) {
  #prot_ref = read.csv('data/srep10775-s2-protein2gene.csv')
  #ilmn_id_ref = Table(getGEO('GPL10558'))
  table = Table(gsm)[,c("ID_REF","VALUE")]
  title = gsm@header$source_name_ch1
  names(table) = c("ID", title)
  idx = c()
  for (id in as.character(prot_ref[,"Entrez.Gene.id"])) { 
    idx = c(idx, which(id_ref["Entrez_Gene_ID"] == id))  
  }
  ids = id_ref[idx, "ID"]
  idx = c()
  for (id in unique(ids)) {
    idx = c(idx, which(table[,"ID"] == as.character(id)))
  }
  if(include_id_col) {
    result = as.data.frame(table[idx,])
    for (id in as.character(result[,"ID"]))
      result[which(as.character(result[,"ID"]) == id), "GENE"] = 
                            nid_ref[which(as.character(id_ref[,"ID"]) == id), "Symbol"]
  } else {
    result = as.data.frame(table[idx,c(title)])
    names(result) = title
  } 
  return(result)
}

build_gene_expr_matrix <- function() {
  gsm_list = get_gsmlist('GSE49577')
  prot_ref = read.csv('srep10775-s2-protein2gene.csv')
  ilmn_id_ref = Table(getGEO('GPL10558'))
  gene_expr_matrix = get_gene_expr(gsm_list[[1]], prot_ref, ilmn_id_ref,
                                   include_id_col=TRUE) 
  for (i in 2:length(gsm_list)) { 
    gene_expr_matrix = cbind(gene_expr_matrix, 
                             get_gene_expr(gsm_list[[i]], prot_ref, ilmn_id_ref))
  }
  return(gene_expr_matrix)
} 

