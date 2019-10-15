








#' @title Fast way to calculate correlation between variables
#' @description Fast way to calculate correlation between variables
#' @inheritParams FastMINE2
#' @param method a character string indicating which correlation coefficient is to be used for the test. One of "pearson", "kendall", or "spearman", can be abbreviated.
#' @importFrom dplyr filter
#' @importFrom plyr adply
#' @return a data frame with p value and adjusted p value.
#' @seealso \code{\link{cor.test}}.
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' data <- state.x77;colnames(data)
#' a <- FastCorrelation(data,
#'                      transposition = F,
#'                      control.markers=NULL,
#'                      target.markers=NULL,
#'                      method = "pearson")
#' a1 <- dplyr::filter(a,adj.P.Val < 0.05)
#'
#' library(ggplot2)
#' ggplot(as.data.frame(data),aes(x=`Life Exp`,y=`Murder`)) + geom_point() + geom_smooth()
#' @export
FastCorrelation <- function(data,
                            transposition = F,
                            control.markers=NULL,
                            target.markers=NULL,
                            method = "pearson"){
  ## matrix preperation
  if(transposition){
    data <- as.data.frame(t(data))
  } else {
    data <- as.data.frame(data,stringsAsFactors = F)
  }

  ## Create pair matrix
  if(is.null(control.markers)){
    LuckyVerbose("Calculate all pairs...")
    pm2 <- as.data.frame(t(combn(colnames(data),2)),stringsAsFactors = F) # str(pm)
    colnames(pm2) <- c("Var1","Var2")
  } else {
    pm <- NULL
    for(i in 1:length(control.markers)){ # i=1
      c.i <- control.markers[i]
      o.i <- setdiff(colnames(data),c.i)
      pm.i <- data.frame(Var1=c.i,Var2=o.i)
      pm <- rbind(pm,pm.i)
    }
    test.r <- apply(pm,1,find.repeat)
    pm2 <- pm[!duplicated(test.r),]
  }
  pm2 <- as.matrix(pm2)

  ## assistant function: calculate one-pair correlation
  one_cor <- function(marker.a,
                      marker.b,
                      data,
                      method = "pearson",
                      alternative = "two.sided"){

    # select a subset of data
    x <- as.numeric(data[,marker.a])
    y <- as.numeric(data[,marker.b])

    ## calculate correlation and apply statistic analysis
    res <- cor.test(x, y, method = method, alternative = alternative)

    ## shape to data frame
    res2 <- data.frame(Var1 = marker.a,
                       Var2 = marker.b,
                       Cor = res$estimate,
                       P.Val = res$p.value)
    return(res2)
  }

  ## Multiple pairs of correlation
  LuckyVerbose("Multiple pairs of correlation...")
  dat.cor <- adply(pm2,1,function(x)one_cor(x[1],x[2],data= data,method = method))
  dat.cor <- dat.cor[-1] # delete redundancy
  dat.cor$adj.P.Val <- p.adjust(dat.cor$P.Val,method = "BH")
  dat.cor <- arrange(dat.cor,adj.P.Val)

  ## Out put
  LuckyVerbose("Done!")
  return(dat.cor)
}





















