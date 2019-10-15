


#' @title conversion between correlation matrix and related long data
#' @description conversion between correlation matrix and related long data
#' @param data matrix or data.frame
#' @param source.col source colnames
#' @param target.col target colnames
#' @param value value colnames
#' @param report whether to do report
#' @importFrom tidyr spread
#' @return a long data or a correlation matrix
#' @seealso \code{\link[tidyr]{spread}}.
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @export
makeCormatrix <- function(data,
                          source.col=NULL,
                          target.col=NULL,
                          value=NULL,
                          report = T){
  ## 判断是相关矩阵还是长型数据
  if(is.null(source.col)|is.null(target.col)|is.null(value)){
    if(report == T) {print("data是相关矩阵，将相关矩阵转化为长型矩阵。")}
    l1 <- list()
    for(i in 1:ncol(data)){
      l.i <- data[,i]
      l1 <- c(l1,list(l.i))
    }
    names(l1) <- colnames(data)
    cm1 <- stack(l1)
    cm1$target <- rep(rownames(data),ncol(data))
    colnames(cm1)[1:2] <- c("value","source")
    # 将值为1和重复的结果去除
    cm2 <- cm1[cm1$value !=1,]
    cm2 <- cm2[order(cm2$value,decreasing = T),]
    select <- seq(2,nrow(cm2),by=2) #选择偶数行
    cm3 <- cm2[select,];rownames(cm3) <- 1:nrow(cm3)
    cm3 <- subset(cm3,select = c("source","target","value"))
    logic1 <- length(unique(as.character(cm3$value))) == length(as.character(cm3$value))
    if(logic1){print("value值并无重复,结果可信。")} else {print("value值有重复，结果可疑，请检查算法。")}
    return(cm3)
  } else {
    if(report == T){print("data是长型矩阵，将长型矩阵转化为相关矩阵。")}
    colnames(data)[Fastmatch(c(source.col,target.col,value),colnames(data))] <- c("source","target","value")
    data1=data.frame(source = data[,"target"],
                     target = data[,"source"],
                     value = data[,"value"])
    co.genes <- unique(c(as.character(data[,"target"]),as.character(data[,"source"])))
    # reverse data
    data1.r <- data1
    colnames(data1.r) <- c("target","source","value")
    # data for self
    data2 = data.frame(source = co.genes,target = co.genes,value = 1)
    # all data
    data3 <- rbind(data1,data1.r,data2)
    # spread
    a <- tidyr::spread(data3,
                       key = source,
                       value = value)
    rownames(a) <- a$target;a1 <- a[,-1];
    a1 <- a1[co.genes,co.genes]
    return(a1)
  }
}

