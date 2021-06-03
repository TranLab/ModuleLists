# ModuleLists

List of gene set modules for downstream applications like GSEA.

# Gene2ModuleExpressionScores function

Description
This function collapses a gene expression matrix to module-level expression scores

Usage
Gene2ModuleExpressionScores(gene_expression_dat, module_list = c("lowBTMs", "highBTMs", "BloodGen3Module","MonacoModules"), summary_stat = c(mean, median)) 

Arguments
gene_expression_dat    gene expression matrix of normalized logCPM values (not counts) or ExpressionSet object containing such a matrix with rownames as HUGO symbols or Ensembl IDs
module_list            name of module set to use. Can be "lowBTMs", "highBTMs", "BloodGen3Module" or "MonacoModules".
summary_stat           mean or median

Value
Output is a dataframe with rownames as module names and column values representing the mean or median or member genes within a module
