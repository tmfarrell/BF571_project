#
# utils.R
# 
# Tim Farrell, tmf@bu.edu 
# Project, BF571, BU
# 20160415
# 
## Prelims
library(lumi)
library(limma)
library(gtools)
library(ggplot2)
library(reshape)
library(Biobase)
library(GEOquery)

setwd(paste0(basedir, 'data/')) 

if (!exists("gsm_list"))
  gsm_list = GSMList(getGEO('GSE49577', GSEMatrix=FALSE))

if (!exists("ilmn_id_ref"))
  ilmn_id_ref = Table(getGEO('GPL10558'))

if (!exists("types")) {
  types = c("OV1002 Carbo", "OV1002 CarboTax", 
            "HOX424 Carbo", "HOX424 CarboTax")
} 

## Functions
nums2join <- function(vec, sep) { 
  return(paste(as.character(vec), collapse=sep))
}

split2nums <- function(str, sep=" ") { 
  return(as.numeric(strsplit(as.character(str), sep)[[1]]))
}

swap <- function(ss) {
  return(paste(ss[[1]], ss[[3]], ss[[2]]))
}

normalize <- function(arr) { 
  return((arr - mean(arr))/sd(arr))
}

normalize_df <- function(df, axis=2) { 
  mols = as.character(df[,"Molecule"])
  df[,"Molecule"] = NULL
  df = as.data.frame(apply(df, axis, FUN=normalize))
  colnames(df) = gsub("\\.", " ", colnames(df))
  df[,"Molecule"] = mols
  return(df)
}

remove_set <- function(names, sep="[.]") { 
  ns = c()
  for (n in strsplit(names, sep))
    ns = c(ns, paste(n[1], n[2], n[3])) 
  return(ns)
}

remove_set_single <- function(name, sep=" ") { 
  ns = strsplit(name, sep)[[1]]
  return(paste(ns[1],ns[2],ns[3]))
}

make_numeric <- function(expr_df) { 
  mols = as.character(expr_df[,"Molecule"])
  expr_df[,"Molecule"] = NULL
  cols = colnames(expr_df)
  expr_df = as.data.frame(lapply(expr_df,
                          FUN=function(x) as.numeric(as.character(x))))
  colnames(expr_df) = cols
  expr_df[,"Molecule"] = mols
  return(expr_df)
}

reduce <- function(expr_df) { 
  mols = as.character(expr_df[,"Molecule"])
  expr_df[,"Molecule"] = NULL
  cols = unique(colnames(expr_df))
  expr_df = reduce_df_by_mean(expr_df)
  colnames(expr_df) = cols
  expr_df[,"Molecule"] = mols
  return(expr_df)
}

reduce_df_by_mean <- function(df) { 
  return(as.data.frame(lapply(split(as.list(df), 
                                    f=colnames(df)),
                              function(x) Reduce(`+`,x) / length(x))))
}

get_gsmlist <- function(gse_str) { 
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

get_gene_prot_dict <- function() {
  prot_ref = read.csv(paste0(basedir,'data/srep10775-s2-protein2gene.csv')) 
  gene_prot_dict= list() 
  for (g in as.character(prot_ref[,"Gene.name"])) {
    p = as.character(prot_ref[which(prot_ref[,"Gene.name"] == g), 
                              "Protein.name"])
    gene_prot_dict[[g]] = p
    gene_prot_dict[[p]] = g
  } 
  return(gene_prot_dict)
}

get_gene_expr <- function(gsm, prot_ref, id_ref, include_id_col=FALSE) {
  #prot_ref = read.csv(paste0(basedir,'data/srep10775-s2-protein2gene.csv'))
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
    result = as.data.frame(table[idx,], stringAsFactors=FALSE)
    for (id in as.character(result[,"ID"]))
      result[which(as.character(result[,"ID"]) == id), "Molecule"] = 
                      id_ref[which(as.character(id_ref[,"ID"]) == id), "Symbol"]
  } else {
    result = as.data.frame(table[idx,c(title)], stringsAsFactors=FALSE)
    names(result) = title
  } 
  return(result)
}

get_gene_exprs <- function(normalize_preprocess=FALSE, lumi=TRUE) {
  prot_ref = read.csv(paste0(basedir,'data/srep10775-s2-protein2gene.csv')) 
  #ilmn_id_ref = Table(getGEO('GPL10558'))
  gene_expr = get_gene_expr(gsm_list[[1]], 
                            prot_ref, ilmn_id_ref, include_id_col=TRUE) 
  for (i in 2:length(gsm_list)) { 
    gene_expr = cbind(gene_expr,
                      get_gene_expr(gsm_list[[i]], prot_ref, ilmn_id_ref))
  }
  gene_expr["ID"] = NULL
  row.names(gene_expr) = 1:nrow(gene_expr)
  if (normalize_preprocess) {
    gene_expr = normalize_df(make_numeric(gene_expr), axis=2)
    comment(gene_expr) = "normalized"
  } else { 
    gene_expr = make_numeric(gene_expr)
    colnames(gene_expr) = gsub("[.]", " ", colnames(gene_expr))
    comment(gene_expr) = ""
  } 
  if (lumi) {
    mols = gene_expr[,"Molecule"]
    gene_expr[,"Molecule"] = NULL
    gene_expr = new("ExpressionSet", exprs=data.matrix(gene_expr))
    gene_expr = as.data.frame(as.matrix(lumiExpresso(gene_expr, 
                                                     normalize.param=list(method='rsn'))))
    gene_expr[,"Molecule"] = mols
    comment(gene_expr) = paste(comment(gene_expr),"lumi-preprocessed")
  }  
  return(gene_expr)
} 

get_lfc <- function(expr_df, type, day) {
  lfc_df = as.data.frame(matrix(ncol=0, nrow=nrow(expr_df)))
  cols = colnames(expr_df)
  expr_df = expr_df[intersect(grep(type, cols), 
                              grep(day, cols))] 
  cols = colnames(expr_df)
  control_idx = grep("Control", cols)
  conditions = c("Carbo", "CarboTax")
  for (c in conditions) { 
    idx = grep(paste(c,""), cols)
    fit = eBayes(lmFit(expr_df[,c(control_idx,idx)], 
                       design=c(rep(0, length(control_idx)), 
                                rep(1, length(idx)))))
    lfc = topTable(fit, adjust.method="fdr", number=100)[,c("logFC","adj.P.Val")]
    lfc = lfc[order(as.numeric(row.names(lfc))),]
    lfc_df[,paste(type,c,day)] = as.numeric(lfc[,"logFC"])
    lfc_df[,paste(type,c,day,"adj.p.val")] = as.numeric(lfc[,"adj.P.Val"])
  }
  return(lfc_df)
}

get_lfcs <- function(expr_df) { 
  types = c("HOX424", "OV1002")
  days = c("Day1", "Day2", "Day4", "Day7", "Day14")
  lfcs = as.data.frame(expr_df[,"Molecule"])
  names(lfcs) = c("Molecule")
  for (t in types) { 
    for (d in days) 
      lfcs = cbind(lfcs, get_lfc(expr_df, t, d))
  }
  return(lfcs)
}

get_gene_diff_expr <- function(gene_expr)  { 
  return(get_lfcs(gene_expr))
}

get_prot_diff_expr <- function() { 
  prot_expr = read.csv(paste0(basedir,'data/srep10775-s2.csv'))
  cs = prot_expr[,"condition"]
  cs = gsub("_", " ", sub("D", "Day", sub("/", "", cs)))
  cs = lapply(strsplit(cs, " "), FUN=swap)
  prot_expr["condition"] = NULL
  prot_expr = as.data.frame(t(prot_expr))
  colnames(prot_expr) = cs
  prot_expr[,"Molecule"] = row.names(prot_expr)
  row.names(prot_expr) = 1:nrow(prot_expr)
  return(prot_expr)
}

compute_differenced <- function(expr, sets=c(""), controls=FALSE) { 
  pval_idx = grepl("adj.p.val", names(expr))
  ns = names(expr)[!pval_idx]
  pval_idx = which(pval_idx == TRUE)[2:length(which(pval_idx == TRUE))]
  differenced_expr = as.data.frame(expr[])
  differenced_expr[,"Molecule"] = expr[,"Molecule"]
  if (controls)
    ts = c(types, "HOX424 Control", "OV1002 Control")
  else
    ts = types
  for (t in ts) { 
    for (s in sets) {
      ts = ns[grep(paste(t, ""), ns)]
      for (i in 2:length(ts))
        differenced_expr[,ts[i]] = expr[,ts[i]] - expr[,ts[i-1]]
    }
  }
  return(differenced_expr)
}

get_type_day <- function(expr_df, type, day, set="") { 
  idx = c() 
  cols = colnames(expr_df)
  for (i in 1:length(cols)) { 
    if (! set == "") {
      if (paste(type,day,set) == cols[i])
        idx = c(idx, i)
    } else if (paste(type,day) == cols[i])
      idx = c(idx, i)
  }
  res = as.data.frame(expr_df[,idx])
  colnames(res) = colnames(expr_df)[idx]
  return(res)
}

get_type <- function(expr_df, type,
                     days=c("Day1", "Day2", "Day4", "Day7", "Day14"), set="") { 
  mols = as.character(expr_df[,"Molecule"])
  df = as.data.frame(mols)
  names(df) = c("Molecule")
  for (d in days) { 
    if (!any(grepl(d, names(df))))
      df = cbind(df, get_type_day(expr_df, type, d, set=set))
  }
  df["Molecule"] = NULL
  df = df[, mixedsort(names(df))]
  df["Molecule"] = mols
  return(df)
}

plot_type <- function(expr_df, type, gene_or_prot="gene", ylab="expression", set="",
                      days=c("Day1", "Day2", "Day4", "Day7", "Day14")) {
  df = get_type(expr_df, type, set=set)
  mols = df[,"Molecule"]
  df_return = df
  df[,"Molecule"] = NULL
  names(df) = 1:length(names(df))
  if (gene_or_prot == "gene") {
    cts = as.data.frame(table(mols))
    unique_mols = c() 
    for (u in as.character(unique(mols)))
      unique_mols = c(unique_mols, paste(u, 1:(cts[cts["mols"] == u, "Freq"])))
    df[,"Group.1"] = unique_mols
  } else { 
    df[,"Group.1"] = mols
  }
  df_m = melt(df)
  p <- ggplot(df_m, aes(x=as.numeric(variable), y=value, colour=Group.1)) +
    geom_line() + ggtitle(paste(type,set)) + xlab("Day") + 
    ylab(ylab) + theme(legend.position="none")
  ggsave(paste0('/Volumes/Data/work/courses/bf571/project/results/', ylab, '-',
                sub(" ","_",type), set, lumi,'.png'), plot=p)
  return(df_return)
}
