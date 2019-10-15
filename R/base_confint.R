




#' @title Confidence interval calulation for sample
#' @description Confidence interval calulation for sample
#' @param x a numeric vector
#' @param sigma the standard deviation.If you know the sd value, please set it to \code{sigma}. Default is -1, which means that you didn't know the sd value.
#' @param alpha significance level
#' @return a data frame of basic statistics information
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' set.seed(2019);x <- sample(1:1000,500,replace = T)
#' confint_vector(x)
#' confint_vector(x,sigma = 2) # Not be used as usual.It's because we always don't know the value of standard deviation.
#' @export
confint_vector <- function(x,
                           sigma=-1,
                           alpha=0.05){
  n<-length(x)
  xb<-mean(x)
  if(sigma>=0)
  {
    tmp<-sigma/sqrt(n)*qnorm(1-alpha/2);df<-n
  }
  else{
    tmp<-sd(x)/sqrt(n)*qt(1-alpha/2,n-1);df<- n-1
  }
  data.frame(mean=xb,df=df,a=xb-tmp,b=xb+tmp)
}

