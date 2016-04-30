# 
# correlation_analysis.R
# 
# Tim Farrell, tmf@bu.edu 
# 20160427
# 
library(Hmisc)
library(graphics)
source(paste0(basedir,'/code/utils.R'))

lumi = comment(gene_expr)

# recompute correlation df?  
compute_corr_df = TRUE
# split up significant vs nonsignificant logFCs in correlation calcs? 
splitup_corr_calcs = FALSE
# ecdfs analysis? 
ecdfs = TRUE

#
# Build correlation dataframe 
#
# first get correlations for each gene/protein expression for each condition
if (compute_corr_df) {
print("Building correlations dataframe...")
genes = as.character(gene_expr[,"Molecule"])
gene_prot_dict = get_gene_prot_dict()
corr_df = data.frame(gene_expr=character(), prot_expr=character(), type=character(), 
                     set=character(), gene=character(), sig_corr=double(), nonsig_corr=double(), 
                     gene_diff_expr_pvals=character(), stringsAsFactors=FALSE)
for (gene_exp in c("gene_expression", "differential_gene_expression",
                   "differenced_differential_gene_expression")) {
  for (prot_exp in c("differential_protein_expression",
                     "differenced_differential_protein_expression")) {
    for (t in types) {
      prot = read.csv(paste0(prot_exp,'-',t,lumi,'.csv'), sep=',')
      colnames(prot) = gsub("[.]", " ", colnames(prot))
      if (!grepl("differential", gene_exp))  {
        for (s in c("","1","2","3")) {
          gene = read.csv(paste0(gene_exp,'-',t,s,lumi,'.csv'), sep=',')
          colnames(gene) = c(remove_set(colnames(gene)[!grepl("p.val", colnames(gene)) & 
                                                         !grepl("Molecule", colnames(gene))]), 
                             "Molecule", colnames(gene)[grepl("p.val", colnames(gene))])
          cols = intersect(colnames(gene), colnames(prot))
          cols = cols[!grepl("Molecule", cols)]
          for (i in 1:length(genes)) {
            # compute cor that gene-type protein diff expr
            g = gene[i, cols]
            p = prot[which(prot[,"Molecule"] == gene_prot_dict[[genes[[i]]]]), cols]
            df = as.data.frame(list(gene_expr=gene_exp, prot_expr=prot_exp, type=t, set=s, 
                                    gene=genes[[i]], sig_corr=cor(as.numeric(g), as.numeric(p), method="pearson"),
                                    nonsig_corr=0.0, gene_diff_expr_pval=""))
            corr_df = rbind(corr_df, df)
          }
        }
      } else {
        # get gene exp
        if (splitup_corr_calcs) {
          gene = read.csv(paste0(gene_exp,'-',t,lumi,'.csv'), sep=',')
          colnames(gene) = c(remove_set(colnames(gene)[!grepl("p.val", colnames(gene)) & 
                                                       !grepl("Molecule", colnames(gene))]), 
                             "Molecule", colnames(gene)[grepl("p.val", colnames(gene))])
          sig_cols = colnames(gene)[1:(grep("Molecule", colnames(gene))-1)][
            as.numeric(gene[1, colnames(gene)[grepl('p.val', colnames(gene))]]) < 0.05]
          nonsig_cols = colnames(gene)[1:(grep("Molecule", colnames(gene))-1)][
            as.numeric(gene[1, colnames(gene)[grepl('p.val', colnames(gene))]]) > 0.05]
          sig_cols = intersect(sig_cols, colnames(prot))
          print(sig_cols)
          nonsig_cols = intersect(nonsig_cols, colnames(prot))
          print(nonsig_cols)
          sig_cols = sig_cols[!grepl("Molecule", sig_cols)]
          nonsig_cols = nonsig_cols[!grepl("Molecule", nonsig_cols)]
          print(sig_cols)
          print(nonsig_cols)
          for (i in 1:length(genes)) {
            # compute cor that gene-type protein diff expr
            diff_pvals = nums2join(gene[i, grepl("p.val", colnames(gene))], sep=" ")
            sig_g = gene[i, sig_cols]
            sig_p = prot[which(prot[,"Molecule"] == gene_prot_dict[[genes[[i]]]]), sig_cols]
            nonsig_g = gene[i, nonsig_cols]
            nonsig_p = prot[which(prot[,"Molecule"] == gene_prot_dict[[genes[[i]]]]), nonsig_cols]
            df = as.data.frame(list(gene_expr=gene_exp, prot_expr=prot_exp, type=t, set=s, 
                               gene=genes[[i]], sig_corr=cor(as.numeric(sig_g), as.numeric(sig_p), method="pearson"),
                               nonsig_corr=cor(as.numeric(nonsig_g), as.numeric(nonsig_p), method="pearson"),
                               gene_diff_expr_pval=diff_pvals))
            corr_df = rbind(corr_df, df)
          }
        } else {
          gene = read.csv(paste0(gene_exp,'-',t,lumi,'.csv'), sep=',')
          colnames(gene) = c(remove_set(colnames(gene)[!grepl("p.val", colnames(gene)) & 
                                                         !grepl("Molecule", colnames(gene))]), 
                             "Molecule", colnames(gene)[grepl("p.val", colnames(gene))])
          cols = intersect(colnames(gene), colnames(prot))
          cols = cols[!grepl("Molecule", cols)]
          for (i in 1:length(genes)) {
            # compute cor that gene-type protein diff expr
            diff_pvals = nums2join(gene[i, grepl("p.val", colnames(gene))], sep=" ")
            g = gene[i, cols]
            p = prot[which(prot[,"Molecule"] == gene_prot_dict[[genes[[i]]]]), cols]
            df = as.data.frame(list(gene_expr=gene_exp, prot_expr=prot_exp, type=t, set=s, 
                                    gene=genes[[i]], sig_corr=cor(as.numeric(g), as.numeric(p), method="pearson"),
                                    nonsig_corr=0.0, gene_diff_expr_pval=diff_pvals))
            corr_df = rbind(corr_df, df)
          }
        }
      }
    }
  }
} 
}


# 
# Empirical CDFs analysis
# 
if(ecdfs) {
print("Obtaining empirical CDFs for the DE vs. nonDE correlations...")
if(grepl("lumi", lumi)) {
corr_df[grepl("differential",corr_df[,"gene_expr"]), "gene_diff_expr_sig"] = 
  as.data.frame(list(gene_diff_expr_sig=c(sapply(corr_df[grepl("differential",corr_df[,"gene_expr"]), "gene_diff_expr_pval"],
                                                 FUN=function(x) all(split2nums(x) < 0.025)))))
} else { 
  corr_df[grepl("differential",corr_df[,"gene_expr"]), "gene_diff_expr_sig"] = 
    as.data.frame(list(gene_diff_expr_sig=c(sapply(corr_df[grepl("differential",corr_df[,"gene_expr"]), "gene_diff_expr_pval"],
                                                   FUN=function(x) all(split2nums(x) < 0.05)))))
}

# differential_gene and differential_protein
diff_expr = corr_df[intersect(which(corr_df[,"gene_expr"] == "differential_gene_expression"), 
                              which(corr_df[,"prot_expr"] == "differential_protein_expression")),]
# test
sig = diff_expr[which(!is.na(diff_expr[,"sig_corr"]) & diff_expr["gene_diff_expr_sig"] == TRUE), "sig_corr"]
nonsig = diff_expr[which(!is.na(diff_expr[,"sig_corr"]) & !diff_expr["gene_diff_expr_sig"] == TRUE), "sig_corr"]
ks = ks.test(sig, nonsig)
w = wilcox.test(sig, nonsig)
write("differential_gene & differential_protein", file=paste0("test_results_",lumi,".txt"))
write("\nKS-test:", file="test_results.txt", append=TRUE)
write(paste("stat:",ks$statistic,"\npval:",ks$p.value), file=paste0("test_results_",lumi,".txt"), append=TRUE)
write("\nWilcoxon-test:", file="test_results.txt", append=TRUE)
write(paste("stat:",w$statistic,"\npval:",w$p.value), file=paste0("test_results_",lumi,".txt"), append=TRUE)
# plot
if (grepl("lumi-preprocessed", lumi))
  pdf(file=paste0(basedir,"Empirical_Correlation_CDFs_DE_vs_nonDE-",lumi,".pdf"))
else
  pdf(file=paste0(basedir,"results/Empirical_Correlation_CDFs_DE_vs_nonDE.pdf"))
Ecdf(sig, main="Empirical Correlation CDFs (DE=black, nonDE=red)\ncor(differential_gene_expr, differential_prot_expr)",
     xlab="Correlation", ylab="CDF", xlim=c(-1.0, 1.0))
Ecdf(nonsig, add=TRUE, col="red")
text(0.5, 0.2, labels=paste("KS test p-value:",ks$p.value))
text(0.5, 0.1, labels=paste("Wilcoxon test p-value:",w$p.value))

# differenced_differential_gene and differenced_differential_protein
diff_expr = corr_df[intersect(which(corr_df[,"gene_expr"] == "differenced_differential_gene_expression"), 
                              which(corr_df[,"prot_expr"] == "differenced_differential_protein_expression")),]
# test 
sig = diff_expr[which(!is.na(diff_expr[,"sig_corr"]) & diff_expr["gene_diff_expr_sig"] == TRUE), "sig_corr"]
nonsig = diff_expr[which(!is.na(diff_expr[,"sig_corr"]) & !diff_expr["gene_diff_expr_sig"] == TRUE), "sig_corr"]
ks = ks.test(sig, nonsig)
w = wilcox.test(sig, nonsig)
write("\n\ndifferenced_differential_gene & differenced_differential_protein", file=paste0("test_results_",lumi,".txt"), append=TRUE)
write("\nKS-test:", file="test_results.txt", append=TRUE)
write(paste("stat:",ks$statistic,"\npval:",ks$p.value), file=paste0("test_results_",lumi,".txt"), append=TRUE)
write("\nWilcoxon-test:", file="test_results.txt", append=TRUE)
write(paste("stat:",w$statistic,"\npval:",w$p.value), file=paste0("test_results_",lumi,".txt"), append=TRUE)
# plot
Ecdf(sig, main="Empirical Correlation CDFs (DE=black, nonDE=red)\ncor(differenced_differential_gene_expr, differenced_differential_prot_expr)",
     xlab="Correlation", ylab="CDF", xlim=c(-1.0, 1.0))
Ecdf(nonsig, add=TRUE, col="red")
text(0.5, 0.2, labels=paste("KS test p-value:",ks$p.value))
text(0.5, 0.1, labels=paste("Wilcoxon test p-value:",w$p.value))

# differenced_differential_gene and differential_protein
diff_expr = corr_df[intersect(which(corr_df[,"gene_expr"] == "differenced_differential_gene_expression"), 
                              which(corr_df[,"prot_expr"] == "differential_protein_expression")),]
# test
sig = diff_expr[which(!is.na(diff_expr[,"sig_corr"]) & diff_expr["gene_diff_expr_sig"] == TRUE), "sig_corr"]
nonsig = diff_expr[which(!is.na(diff_expr[,"sig_corr"]) & !diff_expr["gene_diff_expr_sig"] == TRUE), "sig_corr"]
ks = ks.test(sig, nonsig)
w = wilcox.test(sig, nonsig)
write("\n\ndifferenced_differential_gene & differential_protein", file=paste0("test_results_",lumi,".txt"), append=TRUE)
write("\nKS-test:", file="test_results.txt", append=TRUE)
write(paste("stat:",ks$statistic,"\npval:",ks$p.value), file=paste0("test_results_",lumi,".txt"), append=TRUE)
write("\nWilcoxon-test:", file="test_results.txt", append=TRUE)
write(paste("stat:",w$statistic,"\npval:",w$p.value), file=paste0("test_results_",lumi,".txt"), append=TRUE)
# plot
Ecdf(sig, main="Empirical Correlation CDFs (DE=black, nonDE=red)\ncor(differenced_differential_gene_expr, differential_prot_expr)",
     xlab="Correlation", ylab="CDF", xlim=c(-1.0, 1.0))
Ecdf(nonsig, add=TRUE, col="red")
text(0.5, 0.2, labels=paste("KS test p-value:",ks$p.value))
text(0.5, 0.1, labels=paste("Wilcoxon test p-value:",w$p.value))

# differenced_differential_gene and differential_protein
diff_expr = corr_df[intersect(which(corr_df[,"gene_expr"] == "differential_gene_expression"), 
                              which(corr_df[,"prot_expr"] == "differenced_differential_protein_expression")),]
# test 
sig = diff_expr[which(!is.na(diff_expr[,"sig_corr"]) & diff_expr["gene_diff_expr_sig"] == TRUE), "sig_corr"]
nonsig = diff_expr[which(!is.na(diff_expr[,"sig_corr"]) & !diff_expr["gene_diff_expr_sig"] == TRUE), "sig_corr"]
ks = ks.test(sig, nonsig)
w = wilcox.test(sig, nonsig)
write("\n\ndifferential_gene & differenced_differential_protein", file=paste0("test_results_",lumi,".txt"), append=TRUE)
write("\nKS-test:", file="test_results.txt", append=TRUE)
write(paste("stat:",ks$statistic,"\npval:",ks$p.value), file=paste0("test_results_",lumi,".txt"), append=TRUE)
write("\nWilcoxon-test:", file="test_results.txt", append=TRUE)
write(paste("stat:",w$statistic,"\npval:",w$p.value), file=paste0("test_results_",lumi,".txt"), append=TRUE)
# plot
Ecdf(sig, main="Empirical Correlation CDFs (DE=black, nonDE=red)\ncor(differential_gene_expr, differenced_differential_prot_expr)",
     xlab="Correlation", ylab="CDF", xlim=c(-1.0, 1.0))
Ecdf(nonsig, add=TRUE, col="red")
text(0.5, 0.2, labels=paste("KS test p-value:",ks$p.value))
text(0.5, 0.1, labels=paste("Wilcoxon test p-value:",w$p.value))
#dev.off()
}


# 
# Plot logFC of gene-protein pairs
# 
print("Plotting logFCs of gene-protein pairs...")
ge = c() 
pe = c() 
diff_gene_expr_norm = normalize_df(diff_gene_expr[,colnames(diff_gene_expr)[
                                      !grepl("p.val", colnames(diff_gene_expr))]])
for (c in colnames(diff_gene_expr_norm)[!grepl("Molecule", colnames(diff_gene_expr_norm))]) {
  for (i in 1:length(genes)) {
    result = tryCatch({
      g = diff_gene_expr_norm[i, c]
      p = as.numeric(diff_prot_expr[which(prot_expr[,"Molecule"] == gene_prot_dict[[genes[[i]]]]), c])
      if (!is.na(g) && !is.na(p)) {
        ge = c(ge, g)
        pe = c(pe, p)
      } 
    }, error = function(e) {
      print(paste0("Error on ",c,g))
    }, finally = {})
  }
}
plot(ge, pe, ylab="protein logFC", xlab="mRNA logFC", main="All Measurements")
text(2, -1.1, labels=paste("r2:",summary(lm(pe ~ ge))$adj.r.squared))

ge = c() 
pe = c() 
for (c in colnames(diff_gene_expr_norm)[!grepl("Molecule", colnames(diff_gene_expr_norm))]) {
  for (i in 1:length(genes)) {
    result = tryCatch({
      gexp = mean(as.numeric(diff_gene_expr_norm[i, c]))
      pexp = as.numeric(diff_prot_expr[which(prot_expr[,"Molecule"] == gene_prot_dict[[genes[[i]]]]), c])
      if (!is.na(gexp) && !is.na(pexp)) {
        ge = c(ge, gexp)
        pe = c(pe, pexp)
      } 
    }, error = function(e) {
      print(paste("Error on",c,g))
    }, finally = {})
  }
}
plot(ge, pe, ylab="protein logFC", xlab="mRNA logFC", main="Averaging over genes")
text(2, -1.1, labels=paste("r2:",summary(lm(pe ~ ge))$adj.r.squared))
dev.off() 