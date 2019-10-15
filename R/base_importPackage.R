
#' @title Fast way to use import for a new package
#' @description Fast way to use import for a new package
#' @param path.namespace the path of NAMESPACE file of the new package
#' @importFrom usethis use_package
#' @details 1. The package must be managed via \code{\link[roxygen2]{roxygen2-package}}; \cr 2. The package import relationship must be updated via \code{\link[roxygen2]{roxygenise}}
#' @return a complete import relationship in DESCRIPTION file of a new package
#' @seealso \code{\link[roxygen2]{roxygenise}}; \code{\link[devtools]{use_package}}
#' @note before use \code{importPackage()}, you have to run \code{roxygen2::roxygenise()} before in order to gather correct package usage from this new package.
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' # If you are in the space of a new package, just
#' importPackage()
#' @export
importPackage <- function(path.namespace = "."){

  ## get path of NAMESPACE of the package
  path.namespace <- list.files(path=path.namespace,
                               pattern = "NAMESPACE",
                               full.names = T)

  ## get unique dependant packages
  n1 <- read.table(path.namespace);n1 <- as.character(n1$V1)
  n2 <- n1[grep("importFrom",n1)]
  n3 <- Fastextra(n2,"[(]",2);n3 <- unique(Fastextra(n3,",",1))
  n3 <- sort(n3,decreasing = F)

  ## use usethis::use_package to import dependant packages
  for(i in n3){
    use_package(i,"Imports")
  }

}


