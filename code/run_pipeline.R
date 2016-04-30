# 
# run_pipeline.R
# 
# Tim Farrell, tmf@bu.edu 
# 20160426

# change this to run on your local machine
basedir = '/Volumes/Data/work/courses/bf571/BF571_project/'
source(paste0(basedir,'code/utils.R'))

# run the pipeline both with and without 
# preprocessing with lumi
for (lumi_preprocess in c(TRUE, FALSE)) {
  for (normalize_preprocess in c(FALSE, TRUE)) {
    if (!(lumi_preprocess && !normalize_preprocess)){
    print(paste0("Getting gene expression dataframe lumi=",as.character(lumi_preprocess),
                 " normalize=",as.character(normalize_preprocess),"..."))
    gene_expr = get_gene_exprs(normalize_preprocess=normalize_preprocess, lumi=lumi_preprocess)
    # plot
    print("Running plot_all.R...")
    source(paste0(basedir,'code/plot_all.R'))
    # run correlation analysis
    print("Running correlation_analysis.R...")
    source(paste0(basedir,'code/correlation_analysis.R'))
    # run casual analysis
    print("Running casual_analysis.R...")
    #source(paste0(basedir,'code/casual_analysis.R'))
    print("Done.")
    }
  } 
} 
