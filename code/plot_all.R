# 
# plot_all.R 
# 
# Tim Farrell, tmf@bu.edu
# 20160425

source(paste0(basedir,'code/utils.R'))

sets = c("","1","2","3")
days = c("Day0","Day1","Day2","Day4","Day7","Day14")
lumi = comment(gene_expr)

# plot gene expr
lab = "gene_expression"
for (t in c(types, "HOX424 Control", "OV1002 Control")) {
  for (s in sets)
    result = tryCatch({
      df = plot_type(gene_expr, t, days=days, set=s, ylab=lab)
      write.table(df, file=paste0(lab,'-',t,s,lumi,".csv"), sep=",", row.names=FALSE)
    }, error = function(e) {
      print(paste0("Error on ",t,s))
    }, finally = {})
} 

# compute differenced gene expr; plot
differenced_gene_expr = compute_differenced(gene_expr, sets=sets, controls=TRUE)
lab = "differenced_gene_expression"
for (t in c(types, "HOX424 Control", "OV1002 Control")) { 
  for (s in sets)
    result = tryCatch({
      df = plot_type(differenced_gene_expr, t, days=days, set=s, ylab=lab)
      write.table(df, file=paste0(lab,'-',t,s,lumi,".csv"), sep=",", row.names=FALSE)
    }, error = function(e) {
      print(paste("Error on",t,s))
    }, finally = {})
}

# get differential gene expr; plot
diff_gene_expr = get_gene_diff_expr(gene_expr)
lab = "differential_gene_expression"
for (t in types) { 
  result = tryCatch({
    df = plot_type(diff_gene_expr, t, days=days, ylab=lab)
    df = cbind(df, diff_gene_expr[,grepl(paste(t,""), colnames(diff_gene_expr)) & 
                                   grepl("p.val", colnames(diff_gene_expr))])
    write.table(df, file=paste0(lab,'-',t,lumi,".csv"), sep=",", row.names=FALSE)
  }, error = function(e) {
    print(paste("Error on",t))
  }, finally = {})
}

# compute differenced differential gene expr; plot
differenced_diff_gene_expr = compute_differenced(diff_gene_expr)
lab = "differenced_differential_gene_expression"
for (t in types) { 
  result = tryCatch({
    df = plot_type(differenced_diff_gene_expr, t, days=days, ylab=lab)
    df = cbind(df, differenced_diff_gene_expr[,grepl(paste(t,""), colnames(differenced_diff_gene_expr)) & 
                                              grepl("p.val", colnames(differenced_diff_gene_expr))])
    write.table(df, file=paste0(lab,'-',t,lumi,".csv"), sep=",", row.names=FALSE)
  }, error = function(e) {
    print(paste("Error on",t,s))
  }, finally = {})
}

# get differential prot expr; plot
diff_prot_expr = get_prot_diff_expr()
lab = "differential_protein_expression"
for (t in types) { 
  result = tryCatch({
    df = plot_type(diff_prot_expr, t, days=days, gene_or_prot="protein", ylab=lab)
    write.table(df, file=paste0(lab,'-',t,lumi,".csv"), sep=",", row.names=FALSE)
  }, error = function(e) {
    print(paste("Error on",t))
  }, finally = {})
}

# compute differenced differential gene expr; plot
differenced_diff_prot_expr = compute_differenced(diff_prot_expr)
lab = "differenced_differential_protein_expression"
for (t in types) { 
  result = tryCatch({
    df = plot_type(differenced_diff_prot_expr, t, days=days, 
                   gene_or_prot="protein", ylab=lab)
    write.table(df, file=paste0(lab,'-',t,lumi,".csv"), sep=",", row.names=FALSE)
  }, error = function(e) {
    print(paste("Error on",t))
  }, finally = {})
}

closeAllConnections()