
#' KEGG analysis via clusterProfiler package
#'
#' @description FastKEGG is a fast way to do KEGG analysis via series funtion of clusterProfiler package and automatically save related files or plots
#' @param organism Species.Default is "hsa"(homo species).
#' @param pvalueCutoff a list of pvalue adjust method for \code{\link[clusterProfiler]{enrichKEGG}} and \code{\link[clusterProfiler]{enrichMKEGG}}
#' @inheritParams FastGO
#' @importFrom clusterProfiler enrichKEGG enrichMKEGG
#' @importFrom stringr str_c
#' @importFrom enrichplot dotplot emapplot
#' @importFrom ggplot2 ggtitle
#' @examples
#'Plus.library("clusterProfiler")
#'data(geneList, package='DOSE')
#'genes <- names(geneList);
#'genes <- genes[1:2000]
#'kegg <- FastKEGG(genes = genes,
#'                 organism = 'hsa',
#'                 pvalueCutoff = 0.05,
#'                 qvalueCutoff = 0.2,
#'                 save.path = "KEGG",
#'                 names = "love")
#' View(kegg)
#' @export
FastKEGG <- function(genes,
                     geneList,
                     default.universe = F,
                     organism = 'hsa',
                     pvalueCutoff = list(enrichKEGG = 0.05,
                                         enrichMKEGG = 0.05),
                     qvalueCutoff = 0.1,
                     verbose = T,
                     save.path = "KEGG",
                     names = "love"){
  ## 产生储存文件夹
  old <- options()
  dir <- paste0("./",names,"/",save.path,"/")
  options(warn = -1)
  dir.create(dir,recursive = T)
  options(warn = old$warn)

  ## target genes
  if(length(grep("ENSG",genes)) == 0){
    LuckyVerbose("Gene name is ENTREZ ID(May be)")
  } else {
    LuckyVerbose("Gene name is ENSEMBL ID and would be converted to ENTREZID.")
    ## genes
    entrezID <- convert(genes,totype = "ENTREZID")
    genes <- genes[!is.na(entrezID)]
    report_1 <- length(entrezID[is.na(entrezID)])
    LuckyVerbose(paste0("There are ",report_1," genes with NO ENTREZ ID."),levels = 2)
    genes <- convert(genes,totype = "ENTREZID")
    genes <- unique(genes)

    ## geneList
    entrezID <- convert(names(geneList),totype = "ENTREZID")
    geneList <- geneList[!is.na(entrezID)]
    report_1 <- length(entrezID[is.na(entrezID)])
    LuckyVerbose(paste0("There are ",report_1," geneList members with NO ENTREZ ID."),levels = 2)
    names(geneList) <- convert(names(geneList),totype = "ENTREZID")
  }

  ## repeat test
  #如果转换成entrezid后有重复者，则对genList取mean值
  repeat.test <- table(duplicated(names(geneList)))
  if(T %in% names(repeat.test)){
    r.n <- repeat.test[names(repeat.test) %in% T]
    LuckyVerbose(paste0("There are ",r.n," repeats in geneList and merge them by mean stragegy."),levels = 2)
    geneList2 <- tapply(geneList,names(geneList),mean)
    geneList <- as.numeric(geneList2)
    names(geneList) <- names(geneList2)
    geneList <- sort(geneList,decreasing = T)
  }

  ###============KEGG over-representation test===============###
  LuckyVerbose("Step1: KEGG over-representation test...")
  if(default.universe==T){
    kk1 <- enrichKEGG(gene = genes,
                      organism = organism,
                      keyType = "kegg",
                      pvalueCutoff = pvalueCutoff$enrichKEGG,
                      minGSSize = 10,
                      qvalueCutoff = qvalueCutoff)
  } else {
    ## 使用geneList中的数据作为背景。
    kk1 <- enrichKEGG(gene = genes,
                      universe = names(geneList),
                      organism = organism,
                      keyType = "kegg",
                      pvalueCutoff = pvalueCutoff$enrichKEGG,
                      minGSSize = 10,
                      qvalueCutoff = qvalueCutoff)
  }

  # head(kk)
  d_kk1 <- as.data.frame(kk1)
  write.csv(d_kk1,paste0(dir,names,"_KEGG_enrichment.csv"))

  if(nrow(d_kk1) == 0 & verbose){
    LuckyVerbose(paste0(i,": No KEGG enrichment"),levels = 2)
}

  ###==========KEGG Module over-representation test=========###
  LuckyVerbose("Step2: KEGG Module over-representation test...")
  #KEGG Module is a collection of manually defined function units. In some situation, KEGG Modules have a more straightforward interpretation.

  if(default.universe){
    mkk1 <- enrichMKEGG(gene = genes,
                        organism = 'hsa',
                        pvalueCutoff = pvalueCutoff$enrichMKEGG,
                        minGSSize = 10,
                        qvalueCutoff = 0.2)
  } else {
    ## 使用geneList中的数据作为背景。
    mkk1 <- enrichMKEGG(gene = genes,
                        universe = names(geneList),
                        organism = 'hsa',
                        pvalueCutoff = pvalueCutoff$enrichMKEGG,
                        minGSSize = 10,
                        qvalueCutoff = 0.2)
  }

  d_mkk1 <- as.data.frame(mkk1)
  write.csv(d_mkk1,paste0(dir,names,"_KEGG_Module over-representation test.csv"))
  if(nrow(d_mkk1) == 0 & verbose){
      LuckyVerbose("No enrichment in KEGG Module over-representation test...",levels = 2)
  }

  ## 输出结果
  kegg <- list(
    Repeat = list(
      genes = genes,
      geneList = geneList,
      organism = organism,
      pvalueCutoff = pvalueCutoff,
      qvalueCutoff = qvalueCutoff,
      verbose = verbose,
      save.path = save.path,
      names = names
    ),
    enrichKEGG = kk1,
    enrichMKEGG = mkk1
  )
  save(kegg,file = paste0(dir,"kegg.rda"))
  LuckyVerbose("All done!")
  return(kegg)
}



