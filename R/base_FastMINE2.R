

#' @title Plus version of MINE correlation
#' @description Plus version of MINE correlation.We recommnand \code{FastMINE2} intead of old function \code{FastMINE}.
#' @param data a data frame or a matrix
#' @param transpositon if you have a data frame like sample on cols and markers on rows(etc. a conventional gene expression matrix), you should set \code{transpositon=T}
#' @param control.markers control markers
#' @param target.markers target markers
#' @param seed seed for random number generation reproducibility
#' @param nperm integer, number of permutation to perform
#' @param p.adjust.method method for pvalue adjustment, see \code{\link{p.adjust}} for available methods.
#' @param n.cores the core of parallel calculation
#' @param alpha In the original article, the authors state that the default value α=0.6 (which is the exponent of the search-grid size B(n)=n^{α}) has been empirically chosen. It is worthwhile noting that alpha and C are defined to obtain an heuristic approximation in a reasonable amount of time. In case of small sample size (etc. n<=63) it is preferable to increase alpha to 1 to obtain a solution closer to the theoretical one. Note: if you set α=1 in a big sample size data set, it would be time-consuming.
#' @importFrom minerva mine mictools
#' @importFrom reshape2 melt
#' @importFrom plyr llply desc
#' @importFrom dplyr left_join filter arrange
#' @details It seems that even if you set \code{control.markers=NULL} and {target.markers=NULL} at the same time to calculate all pairs, it would not take much time than hours, which is recommanded to do deeper exploration for your data.
#' @seealso \code{\link[minerva]{mine}}.
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' library(minerva)
#' data("Spellman")
#' data <- Spellman[1:100] # dim(Spellman)
#' colnames(data)
#'
#' dat.cor <- FastMINE2(data)
#'
#' dat.cor1 <- FastMINE2(data,control.markers = c("time"))
#'
#' dat.cor2 <- FastMINE2(data,control.markers = c("time","YAL001C"))
#' @export
FastMINE2 <- function(data,
                      transposition = F,
                      control.markers=NULL,
                      target.markers=NULL,
                      alpha=0.6,
                      nperm=1000,
                      seed=1234,
                      p.adjust.method="BH",
                      n.cores = 1){

  ## transposition
  if(transposition){
    data <- as.data.frame(t(data),stringsAsFactors = F)
  }

  if(is.null(control.markers) & is.null(target.markers)){
    ## all pairs
    dat.cor <- mine_all(data = data,
                        n.cores = n.cores,
                        alpha=alpha,
                        nperm=nperm,
                        seed=seed,
                        p.adjust.method=p.adjust.method)

  } else {
    ## Specified pairs
    dat.cor <- mine_part(data=data,
                         control.markers=control.markers,
                         target.markers=target.markers,
                         alpha=alpha,
                         nperm=nperm,
                         seed=seed,
                         p.adjust.method=p.adjust.method,
                         n.cores = n.cores)
  }

  return(dat.cor)

}

####===============Assistant Functions===================####
mine_all <- function(data,
                     n.cores = 1,
                     alpha=0.6,
                     nperm=1000,
                     seed=1234,
                     p.adjust.method="BH"){

  ## Calulate MINE results
  LuckyVerbose("Calculate all-pair MINE...")
  res.mine <- mine(x=data,master = NULL,n.cores = n.cores,alpha=alpha)

  ## Annotate mine result
  LuckyVerbose("Shape MINE result into data frame...")
  res.mine.df <- annotate_mine(res.mine)

  ## Calulate MIC P value and adjusted P value
  mt.test <- as.matrix(data)
  LuckyVerbose("Calculate all-pair P value for MIC...")
  ticenull <- mictools(mt.test,nperm=nperm,seed=seed,p.adjust.method=p.adjust.method)
  data.p1 <- ticenull$pval # colnames(data.p1)
  data.p <- data.p1[c("Var1","Var2","pval","adj.P.Val")]

  ## merge data
  data.p$f.repeat <- apply(data.p,1,function(x)find.repeat(x,var=c("Var1","Var2")))
  res.mine.df$f.repeat <- apply(res.mine.df,1,function(x)find.repeat(x,var=c("Var1","Var2")))
  newdata <- left_join(data.p,res.mine.df,by="f.repeat")
  newdata <- newdata[c("Var1.x","Var2.x","MIC","MAS","MEV","MCN","MICR2","GMIC","TIC","pval","adj.P.Val")]
  colnames(newdata)[c(1,2,10)] <- c("Var1","Var2","P.Val")
  newdata <- arrange(newdata,adj.P.Val,desc(MIC))

  ## output
  l <- list(
    long = newdata,
    cormatrix = res.mine
  )
  LuckyVerbose("Done!")
  return(l)
}


mine_part <- function(data,
                      control.markers="time",
                      target.markers=NULL,
                      alpha=0.6,
                      nperm=1000,
                      seed=1234,
                      p.adjust.method="BH",
                      n.cores = 1){

  ## Create pair matrix
  LuckyVerbose("Create pair matrix...")
  pm <- NULL
  for(i in 1:length(control.markers)){ # i=1
    c.i <- control.markers[i]
    o.i <- setdiff(colnames(data),c.i)
    pm.i <- data.frame(Var1=c.i,Var2=o.i)
    pm <- rbind(pm,pm.i)
  }
  test.r <- apply(pm,1,find.repeat)
  pm2 <- pm[!duplicated(test.r),]
  pm2 <- as.matrix(pm2)

  # old code
  if(F){
    pm <- as.data.frame(t(combn(colnames(data),2)),stringsAsFactors = F) # str(pm)
    colnames(pm) <- c("Var1","Var2")
    if(is.null(control.markers)){
      pm2 <- pm
    } else {
      if(is.null(target.markers)){
        tm <- setdiff(colnames(data),control.markers)
      } else {
        tm <- target.markers
      }

      pm2 <- filter(pm,Var1 %in% control.markers,Var2 %in% tm)
    }
    pm2 <- as.matrix(pm2)
  }

  ## MINE Value
  get1 <- function(x)one_mine(x[1],x[2],data = data,alpha = alpha,nperm = nperm,seed = seed,p.adjust.method = p.adjust.method)
  LuckyVerbose("Calcualte Specified pairs MINE...")
  dat.cor <- adply(pm2,1,get1)

  if(F){
    if(n.cores == 1){
      LuckyVerbose("Calcualte Specified pairs MINE...")
      dat.cor <- adply(pm2,1,get1)
    } else {
      cl <-  multidplyr::default_cluster()
      multidplyr::cluster_assign_partition(cl)
      multidplyr::cluster_call(cl,get1)
    }
  }

  dat.cor <- dat.cor[-1]
  dat.cor$adj.P.Val <- p.adjust(dat.cor$P.Val,method = p.adjust.method)
  newdata <- arrange(dat.cor,adj.P.Val,desc(MIC))
  LuckyVerbose("Done!")
  return(newdata)
}

## repeat-finder
find.repeat <- function(vt,var=c("Var1","Var2")){
  vt1 <- sort(vt[var])
  vt2 <- paste(vt1,collapse = "_")
  return(vt2)
}

## shape mine result into data frame
annotate_mine <- function(res.mine){

  ## Get Long Data
  get_long <- function(mt){
    mt1 <- melt(mt) # reshape2::
    mt1$f.repeat <- apply(mt1,1,find.repeat)
    mt2 <- mt1[!duplicated(mt1$f.repeat),]
    mt3 <- mt2[-4] # delete repeat label
    mt3 <- as.data.frame(as.matrix(mt3),stringsAsFactors = F) # delete factor level
    mt4 <- filter(mt3,`Var1` != `Var2`) # dplyr::

    ## Output
    return(mt4)
  }
  res.mine <- llply(res.mine,get_long)

  ## Annotation
  annot.markers <- names(res.mine)
  df <- res.mine[[1]]
  colnames(df)[3] <- annot.markers[1]
  for(i in 2:length(res.mine)){ # i=2
    df.i <- res.mine[[i]]
    colnames(df.i)[3] <- annot.markers[i]
    df <- left_join(df,df.i,by=c("Var1","Var2"))
  }
  return(df)

}

## assistant function: calculate one-pair MINE
one_mine <- function(marker.a,
                     marker.b,
                     data,
                     alpha=0.6,
                     nperm=1000,
                     seed=1234,
                     p.adjust.method = "BH"){
  ## Calculate one MINE
  x <- data[,c(marker.a,marker.b)]
  res.mine <- mine(x=x,master = match(marker.a,colnames(x)),alpha=alpha)
  res.mine.df <- annotate_mine(res.mine)

  ## Calculate P value
  ticenull.i <- mictools(as.matrix(x),nperm=nperm,seed=seed,p.adjust.method=p.adjust.method)
  p.i <- ticenull.i$pval$pval

  ## Merge data
  dat.mine <- cbind(res.mine.df,P.Val = p.i)

  ## Output
  return(dat.mine)
}


