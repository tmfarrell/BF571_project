# 
# causal_analysis.R
# 
# Tim Farrell, tmf@bu.edu 
# 20160426
# 
library(lmtest)
library(tseries)
source(paste0(basedir,'/code/utils.R'))

lumi=""
compute_gcause=FALSE

# compute adfs
ADFs = data.frame(expr=character(), type=character(), stat=double(), pval=double(), 
                  molecule=character(), set=character())
for (expr in c("gene_expression", "differential_gene_expression", "differenced_differential_gene_expression",
               "differential_protein_expression", "differenced_differential_protein_expression")) {
    for (t in types) {
      if (!grepl("differential", expr))
        sets = c("", "1", "2", "3")
      else
        sets = c("")
      for (s in sets) {
        e = read.csv(paste0(expr,'-',t,s,lumi,'.csv'), sep=',')
        mols = as.character(e[,"Molecule"])
        e[,"Molecule"] = NULL
        colnames(e) = remove_set(colnames(e)[!grepl("Molecule", colnames(e))])
        e[,"Molecule"] = mols
        for (i in 1:length(mols)) {
          result = tryCatch({
          adf = adf.test(as.numeric(e[i, grepl(t, colnames(e))]), 
                         alternative="stationary", k=0)
          ADFs = rbind(ADFs, data.frame(list(expr=expr, type=t, molecule=mols[[i]], 
                                             stat=adf$statistic, pval=adf$p.value, set=s)))
          }, error = function(e) { 
            print(paste("Error on:", expr, t, s, lumi))
          }, finally={})
        } 
      } 
    } 
} 

# get those ts that can be taken as stationary
stationary = ADFs[which(ADFs[,"pval"] <= 0.01),]
stationary_gene = stationary[grepl("gene", stationary[,"expr"]),]
stationary_prot = stationary[grepl("prot", stationary[,"expr"]),]

# compute granger causality b/t those ts
if(compute_gcause) {
g_cause = data.frame(gene_expr=character(), prot_expr=character(), type=character(),
                     gene=character(), prot=character(), set=character(), 
                     stat=double(), pval=double())
for (t in as.character(unique(stationary_gene[,"type"]))) { 
  g_idx = which(stationary_gene[,"type"] == t)
  for (i in g_idx) { 
    set = stationary_gene[i,"set"]
    g_expr = stationary_gene[i,"expr"]
    gene = as.character(stationary_gene[i, "molecule"])
    gene_expr = read.csv(paste0(g_expr, '-', t, s, lumi, '.csv'), sep=',')
    gene_expr = gene_expr[which(gene_expr[,"Molecule"] == gene), 
                                   !grepl("p.val", colnames(gene_expr)) &
                                   !grepl("Molecule", colnames(gene_expr))]
    for (ix in 1:nrow(gene_expr)) {
      gene_ts = as.numeric(gene_expr[ix,])
      p_idx = which(stationary_prot[,"type"] == t)
      for (j in p_idx) { 
        p_expr = stationary_prot[j, "expr"]
        prot = as.character(stationary_prot[j, "molecule"])
        prot_expr = read.csv(paste0(p_expr, '-', t, lumi, '.csv'), sep=',')
        prot_ts = as.numeric(prot_expr[which(prot_expr[,"Molecule"] == prot), 
                                       !grepl("Molecule", colnames(prot_expr))])
        result = tryCatch({
          gt = grangertest(prot_ts ~ gene_ts)
          g_cause = rbind(g_cause, data.frame(list(gene_expr=g_expr, prot_expr=p_expr, 
                                                   gene=gene, prot=prot, type=t, set=s,
                                                   stat=gt$F[2], pval=gt$"Pr(>F)"[2])))
        }, error = function(e) {
          print(paste("Error on:",g_expr,p_expr,t,e))
        }, finally = {}) 
      }
    }
  }
}
}

