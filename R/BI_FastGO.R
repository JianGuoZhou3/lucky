

#' @title GO analysis via clusterProfiler package
#' @description FastGOplot give fast way to draw GO series plot like barplot/dotplot/emapplot,Cnetplot and goplot for a enrichGO object based on clusterProfiler package.
#' @param genes a character of gene id.
#' @param geneList a numeric with id names.For example,logFoldChang with ENSEMBL id names.Note that the names of \code{geneList} must be the same type of \code{genes}
#' @param default.universe whether to use default universe.If \code{default.universe = F},the background genes would be provide by \code{geneList}
#' @param classlevel Default is 2.I had tested that 2:7 was still available.If many levels is set,the process is time-consuming
#' @param OrgDb OrgDb dataset.If NULL,use "org.Hs.eg.db"
#' @param keyType the type of gene.Support automatically test
#' @param pAdjustMethod the mathod of pvalue adjustment
#' @param pvalueCutoff  a list of pvalue adjust method for \code{\link[clusterProfiler]{enrichGO}}
#' @param qvalueCutoff  cutoff of q value
#' @param cnet.showCategory the number of showed cluster at cnetplot
#' @param verbose LuckyVerbose gseplot running message or not
#' @param save.path the sub path of saved files
#' @param names the main path and part names of saved files
#' @importFrom clusterProfiler groupGO enrichGO
#' @importFrom stringr str_c
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' library(lucky)
#' need <- c("clusterProfiler","org.Hs.eg.db");Plus.library(need)
#' data(geneList, package = "DOSE")
#' data(geneList, package='DOSE')
#' universe <- unique(as.character(common.annot$ENTREZID))
#' system.time(
#' l <- FastGO(genes,
#'             geneList,
#'             classlevel = 2:2,
#'             OrgDb  = NULL,
#'             keyType = NULL,
#'             pAdjustMethod = "BH",
#'             pvalueCutoff = list(
#'                gseGO = 0.05,
#'                enrichGO = 0.05),
#'             qvalueCutoff  = 0.05,
#'             cnet.showCategory = 5,
#'             save.path = "GO",
#'             names = "love")
#' )
#' ## enhanced plot
#' g <- enrichplot::cnetplot(goList$GOEnrichment$BP,
#'                           showCategory=10,
#'                           colorEdge=T,
#'                           node_label=T,
#'                           circular = F)
#' g1 <- g + theme(legend.text = element_text(face = "bold",size = 12))
#'
#' @export
FastGO <- function(genes,
                   geneList,
                   default.universe = F,
                   classlevel = 1:2,
                   OrgDb  = NULL,
                   keyType = NULL,
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff  = 0.05,
                   verbose = TRUE,
                   save.path = "GO",
                   names = "love"){
  ### 差异性基因的多个Categery的enrichGO analysis
  Go.categery <- c("MF","BP","CC")

  ### reference database
  if(is.null(OrgDb)){
    OrgDb <- "org.Hs.eg.db"
    LuckyVerbose("You don't specify the 'OrgDb' parameter.Here we use 'org.Hs.eg.db'(homo spacies)...")
  }

  ### whether is ensembl id
  if(is.null(keyType)){
    lg1 <- length(grep("ENSG",names(geneList))) == 0
    if(lg1){
      #说明不是ensemblid
      LuckyVerbose("The names of geneList is ENTREZ ID...")
      keyType <-  "ENTREZID"
    } else {
      LuckyVerbose("The names of geneList is ENSEMBL ID...")
      keyType <-  "ENSEMBL"
    }
  }

  ## 产生储存文件夹
  old <- options()
  dir <- paste0("./",names,"/",save.path,"/")
  options(warn = -1)
  dir.create(dir,recursive = T)
  options(warn = old$warn)

  ## geneList 排序
  LuckyVerbose("Step1: order 'geneList' from upper to lower...")
  geneList <- sort(geneList,decreasing = T)

  ###===================GO classification=================###
  go.classification <- NULL
  if(F){
  LuckyVerbose("Step2: GO classification...")
  #for(ont in Go.categery){
  #  a <- NULL
  #  for(level in classlevel){
  #    ggo <- groupGO(gene     = genes,
  #                   keyType  = keyType,
  #                   OrgDb    = OrgDb,
  #                   ont      = ont,
  #                   level    = 7,
  #                   readable = TRUE);
  #    go.classification <- c(go.classification,list(ggo))
  #    names(go.classification)[grep(ont,Go.categery)] <- ont
  #    a.i <- as.data.frame(ggo)
  #    if(verbose == TRUE){
  #      LuckyVerbose(paste0(ont,": GO classification_level ",level," was collected!"),levels = 2)
  #    }
  #    a.i <- a.i[a.i$Count != 0,]
  #    a <- rbind(a,a.i)
  #  }
  #   a <- a[!duplicated(a$ID),]
  #  a$categery <- ont
  #
  #}
  #save(go.classification,file = paste0(dir,names,"_go.classification.rda"))
  }

  ###==================GO over-representation test========###
  if(T){
  LuckyVerbose("GO over-representation test...")
  go.enrich <- NULL#记录egoList的向量。
  for (i in 1:length(Go.categery)) {# i=1
    ont <- Go.categery[i]
    if(verbose){
      LuckyVerbose(paste0(ont,": enrichGO List is on establishment..."),levels = 2)
    }
    if(default.universe==T){
      LuckyVerbose("You try to use a default background provided by ClusterProfilter.In factor,a self-defined background is more proper.")
      ego.sig.tv <- enrichGO(
        gene = genes,
        OrgDb  = OrgDb,
        ont   = ont,
        keyType = keyType,
        pAdjustMethod = pAdjustMethod,
        pvalueCutoff  = pvalueCutoff,
        qvalueCutoff  = qvalueCutoff,
        readable = TRUE)
    } else {
      ego.sig.tv <- enrichGO(
        gene = genes,
        universe  = names(geneList),
        OrgDb  = OrgDb,
        ont   = ont,
        keyType = keyType,
        pAdjustMethod = pAdjustMethod,
        pvalueCutoff  = pvalueCutoff,
        qvalueCutoff  = qvalueCutoff,
        readable = TRUE)
    }
    go.enrich <- c(go.enrich,list(ego.sig.tv))
    names(go.enrich)[i] <- ont
    if(verbose){
      LuckyVerbose(paste0(ont,": enrichGO List completed!"),levels = 2)
    }
    ego.i <- as.data.frame(ego.sig.tv)
    ## save result
    ego.sig.tv.savepath <- str_c(dir,names,"_",ont,"_enrichGO analysis.csv");
    write.csv(ego.i,ego.sig.tv.savepath)

    if(nrow(ego.i) == 0 & verbose){
      LuckyVerbose(paste0(i,": No GO enrichment"),levels = 2)
    }
  }
}

  ###==================Output the result===================###
  GO <- list(
   Repeat = list(
     genes = genes,
     geneList = geneList,
     classlevel = classlevel,
     OrgDb  = OrgDb,
     keyType = keyType,
     pAdjustMethod = pAdjustMethod,
     pvalueCutoff =  pvalueCutoff,
     qvalueCutoff  = qvalueCutoff,
     verbose = verbose,
     save.path = save.path,
     names = names,
     date = Sys.Date()
   ),
   GOClassification = go.classification,
   GOEnrichment = go.enrich
  )
  save(GO,file = paste0(dir,names,"_GO-enrichment.rda"))
  LuckyVerbose("All done!")
  return(GO)
}





