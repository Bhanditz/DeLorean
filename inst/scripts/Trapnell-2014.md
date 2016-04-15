




# Investigate Trapnell et al. (2014) single-cell RNA-seq data

Load data from GEO file.

```r
data.dir <- '../../data'
truseq <- 'GSE52529_truseq_fpkm_matrix.txt.gz'
trapnell <- read.csv(paste(data.dir, 'GSE52529_fpkm_matrix.txt.gz', sep='/'),
                     sep="\t")
```

```
## Warning in file(file, "rt"): cannot open file '../../data/
## GSE52529_fpkm_matrix.txt.gz': No such file or directory
```

```
## Error in file(file, "rt"): cannot open the connection
```

```r
trapnell$gene <- factor(rownames(trapnell), levels=rownames(trapnell))
```

```
## Error in rownames(trapnell): object 'trapnell' not found
```

```r
trapnell.l <- (
    melt(trapnell, variable.name="cell", value.name="fpkm")
    %>% mutate(fpkm.log10=log10(fpkm)))
```

```
## Error in melt(trapnell, variable.name = "cell", value.name = "fpkm"): object 'trapnell' not found
```

Parse the cell meta data.

```r
cell.meta <- (
    data.frame(cell=factor(levels(trapnell.l$cell)))
    %>% mutate(capture=str_match(levels(cell),
                                 "^T([0-9]+)")[,2],
               obstime=as.integer(capture),
               capture=factor(capture)))
```

```
## Error in levels(trapnell.l$cell): object 'trapnell.l' not found
```

```r
min.expr <- 0.1  # Minimum expression value for Tobit model
expr.threshold <- 1  # The value above which genes are considered expressed
(ggplot(trapnell.l
        %>% filter(fpkm > min.expr)
        %>% sample_n(10000)
        %>% left_join(cell.meta),
        aes(x=fpkm, color=factor(obstime)))
    + geom_density()
    + geom_rug(alpha=.01)
    + scale_x_log10()
)
```

```
## Error in eval(expr, envir, enclos): object 'trapnell.l' not found
```


Filter genes as in Trapnell et al.

```r
gene.meta <- (
    trapnell.l
    %>% group_by(gene)
    %>% dplyr::summarise(num.expr=sum(fpkm > expr.threshold),
                  lpe.sd=sd(fpkm.log10[fpkm > min.expr]))
    %>% filter(num.expr >= 50, lpe.sd > .7))
```

```
## Error in eval(expr, envir, enclos): object 'trapnell.l' not found
```

```r
qplot(gene.meta$num.expr, binwidth=5)
```

```
## Error in eval(expr, envir, enclos): object 'gene.meta' not found
```

```r
qplot(gene.meta$lpe.sd)
```

```
## Error in eval(expr, envir, enclos): object 'gene.meta' not found
```

```r
trapnell.f <- (
    trapnell.l
    %>% filter(gene %in% gene.meta$gene)
    %>% left_join(cell.meta)
)
```

```
## Error in eval(expr, envir, enclos): object 'trapnell.l' not found
```

```r
length(unique(trapnell.f$gene))
```

```
## Error in unique(trapnell.f$gene): object 'trapnell.f' not found
```


Use Tobit model to filter genes on a differential expression test.

```r
fit.model <- function(form) {
    vgam(form, family=tobit(Lower=log10(min.expr), Upper=Inf))
}
Tobit.Diff.Expr.LRT <- function(fpkm.log10, obstime) {
    FM.fit <- fit.model(fpkm.log10 ~ bs(obstime, df=3))
    # summary(FM.fit)
    RM.fit <- fit.model(fpkm.log10 ~ 1)
    # summary(RM.fit)
    if (is.null(FM.fit) == FALSE && is.null(RM.fit) == FALSE) {
        lrt <- lrtest(FM.fit, RM.fit)
        return(lrt@Body["Pr(>Chisq)"][2,])
    } else { return(1) }
}
genes.high.sd <- filter(gene.meta, num.expr >= 50, lpe.sd > 1)$gene
```

```
## Error in filter_(.data, .dots = lazyeval::lazy_dots(...)): object 'gene.meta' not found
```

```r
gene.diff.expr <- (trapnell.f
    %>% filter(gene %in% genes.high.sd)
    %>% group_by(gene)
    %>% dplyr::summarise(p.val=Tobit.Diff.Expr.LRT(fpkm.log10, obstime))
    %>% arrange(p.val)
)
```

```
## Error in eval(expr, envir, enclos): object 'trapnell.f' not found
```

```r
gene.diff.expr
```

```
## Error in eval(expr, envir, enclos): object 'gene.diff.expr' not found
```

```r
qplot(gene.diff.expr$p.val, binwidth=1) + scale_x_log10()
```

```
## Error in eval(expr, envir, enclos): object 'gene.diff.expr' not found
```

```r
genes.diff.expr <- filter(gene.diff.expr, p.val < 1e-2)$gene
```

```
## Error in filter_(.data, .dots = lazyeval::lazy_dots(...)): object 'gene.diff.expr' not found
```

```r
length(genes.diff.expr)
```

```
## Error in eval(expr, envir, enclos): object 'genes.diff.expr' not found
```


Get gene names.

```r
library("biomaRt")
```

```
## Error in library("biomaRt"): there is no package called 'biomaRt'
```

```r
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
```

```
## Error in eval(expr, envir, enclos): could not find function "useMart"
```

```r
filters <- listFilters(ensembl)
```

```
## Error in eval(expr, envir, enclos): could not find function "listFilters"
```

```r
attributes = listAttributes(ensembl)
```

```
## Error in eval(expr, envir, enclos): could not find function "listAttributes"
```

```r
gene.meta$ensembl_gene_id <- str_match(gene.meta$gene, '(ENSG[0-9]+)')[,1]
```

```
## Error in stri_match_first_regex(string, pattern, opts_regex = attr(pattern, : object 'gene.meta' not found
```

```r
gene.names <- getBM(attributes=c('ensembl_gene_id',
                                 'hgnc_symbol'),
                    mart = ensembl)
```

```
## Error in eval(expr, envir, enclos): could not find function "getBM"
```

```r
gene.meta <- (
    gene.meta
    %>% left_join(gene.names)  # Join the name to the ID
    %>% dplyr::select(-ensembl_gene_id)
    %>% group_by(gene)  # Just take first name for each gene
    %>% do(head(., 1))
    %>% ungroup()
)
```

```
## Error in eval(expr, envir, enclos): object 'gene.meta' not found
```


Reshape data into general format for model.

```r
fpkm <- trapnell.f %>% dcast(gene ~ cell, value.var="fpkm")
```

```
## Error in eval(expr, envir, enclos): object 'trapnell.f' not found
```

```r
fpkm.m <- as.matrix(fpkm %>% dplyr::select(-gene))
```

```
## Error in eval(expr, envir, enclos): object 'fpkm' not found
```

```r
rownames(fpkm.m) <- fpkm$gene
```

```
## Error in eval(expr, envir, enclos): object 'fpkm' not found
```

```r
expr <- log10(fpkm.m + 1)
```

```
## Error in eval(expr, envir, enclos): object 'fpkm.m' not found
```

Save data.

```r
rda.path <- paste(data.dir, "TrapnellDeLorean.rda", sep="/")
message('Saving expression data and meta-data to: ', rda.path)
```

```
## Saving expression data and meta-data to: ../../data/TrapnellDeLorean.rda
```

```r
rename.and.save <- function(..., file) {
    x <- list(...)
    save(list=names(x), file=file, envir=list2env(x))
}
rename.and.save(
    trapnell.expr=expr,
    trapnell.gene.meta=gene.meta,
    trapnell.cell.meta=cell.meta,
    file=rda.path)
```

```
## Error in rename.and.save(trapnell.expr = expr, trapnell.gene.meta = gene.meta, : object 'expr' not found
```

R version and packages used:

```r
sessionInfo()
```

```
## R version 3.2.3 (2015-12-10)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 15.10
## 
## locale:
##  [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
##  [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
##  [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] splines   stats4    stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] stringr_1.0.0       reshape2_1.4.1      dplyr_0.4.3        
##  [4] VGAM_1.0-0          ggplot2_2.0.0       DeLorean_1.0       
##  [7] rmarkdown_0.8.1     knitcitations_1.0.7 knitr_1.12.3       
## [10] devtools_1.9.1      vimcom_1.2-7        setwidth_1.0-4     
## [13] colorout_1.1-1     
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_0.12.2       formatR_1.2.1     plyr_1.8.3       
##  [4] bitops_1.0-6      tools_3.2.3       extrafont_0.17   
##  [7] digest_0.6.9      evaluate_0.8      lattice_0.20-33  
## [10] lubridate_1.5.0   memoise_0.2.1     gtable_0.1.2     
## [13] bibtex_0.4.0      DBI_0.3.1         yaml_2.1.13      
## [16] parallel_3.2.3    Rttf2pt1_1.3.3    coda_0.18-1      
## [19] gridExtra_2.0.0   RefManageR_0.8.63 httr_1.0.0       
## [22] roxygen2_5.0.1    grid_3.2.3        inline_0.3.14    
## [25] R6_2.1.1          XML_3.98-1.3      rstan_2.9.0-3    
## [28] RJSONIO_1.3-0     extrafontdb_1.0   magrittr_1.5     
## [31] MASS_7.3-44       scales_0.3.0      htmltools_0.2.6  
## [34] assertthat_0.1    colorspace_1.2-6  stringi_1.0-1    
## [37] functional_0.6    RCurl_1.95-4.7    munsell_0.4.2
```
