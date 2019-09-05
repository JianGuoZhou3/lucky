


#' @title Build markdown citation report
#' @description For .md files,when given specified reference form(.csl), citation library(.bib), \code{mardownCitation} product a .bat file to annotate the citations of .md by pandoc.
#' @param md.path the path of markdown doc
#' @param csl.path the path of .csl file. Default is "american-medical-association.csl" provided by system
#' @param bib.path the .bib citation library from cication editor(like Endnote). when you use endnote, you can select "BibTex Export" output form to produce a .txt
#' @param save.names the names of saved files
#' @param save.type the type of saved files. Default is "docx"
#' @param save.path the path of output files
#' @param sub.bat whether output indivitual .bat file. Default if \code{FALSE}
#' @details 1. The path should be a absolute path instead of relative one. \cr 2. You have to install pandoc before. More detailed installation please go to \url{https://github.com/jgm/pandoc/releases/tag/2.7.3} \cr 3.It seems that only "docx" is available by pandoc-citeproc filter
#' @return a list
#' @seealso csl library \url{https://editor.citationstyles.org/about/}
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' cit <- mdCitation()
#' cit <- mdCitation(csl.path = "american-medical-association")
#' @export
mdCitation <- function(md.path = NULL,
                       csl.path = NULL,
                       bib.path = NULL,
                       save.names = NULL,
                       save.type = "docx",
                       sub.bat = F,
                       save.path = NULL){
  ## Get .md path
  if(is.null(md.path)){
    mp <- getwd()
    # mp <- "E:/iProjects/exosome/learn/exosome_1A/hnRNPA2B1相关知识"
  } else {
    mp <- md.path
  }
  md.files <- list.files(mp,pattern = ".md",full.names = T)
  md.names <- Fastextra(list.files(mp,pattern = ".md",full.names = F),"[.]",1)

  ## Get csl path
  csl_dir <- system.file("extdata/CSL", package = "lucky")
  csl_path <-list.files(csl_dir,pattern = ".csl",full.names = T)
  csl_name <- list.files(csl_dir,pattern = ".csl",full.names = F)
  # report csl
  LuckyVerbose("You can select:")
  for(i in csl_name){
    LuckyVerbose(i,levels = 2)
  }

  # select a .csl
  if(is.null(csl.path)){
    csl_file <- csl_path[grep("american-medical-association",csl_path)]

  } else {
    if(length(grep(csl.path,csl_path)) != 0){
      # use specifed system .csl
      csl_file <- csl_path[grep(csl.path,csl_path)]
    } else {
      # user-defined .cls
      csl_file <- csl.path
    }
  }

  ## Get .bib path
  if(is.null(bib.path)){
    bib.path <- mp
    bib.files <- list.files(bib.path,pattern = ".bib",full.names = T)
    bib.file <- bib.files[1]
    if(is.na(bib.file)){
      LuckyVerbose("Notice! You have to give a .bid library!")
    }
  }

  ## Get save path
  if(is.null(save.path)){
    save.path <- mp
  }

  ## save.names
  if(is.null(save.names)){
    save.names <- md.names
  }

  ## pandoc expression
  PE <- NULL
  LuckyVerbose("Build pandoc expression:")
  for(i in 1:length(md.files)){ #i=1
    mf <- md.files[i]
    pe <- paste("pandoc",
                "--filter pandoc-citeproc",
                paste0("--bibliography=",bib.file),
                paste0("--csl=",csl_file),
                mf,
                "-o",
                paste0(save.path,"/",save.names[i],".",save.type)
    );

    ## save script as .bat
    if(sub.bat){
      write.table(pe,paste0(mp,"/",save.names[i],"_",save.type,".bat"),
                  sep = "\t",quote = F,
                  col.names = F,row.names = F)
      LuckyVerbose("Complete",save.names[i],levels = 2)
    }
    PE <- c(PE,pe)
  }

  ## save all scripts as .bat
  write.table(PE,paste0(mp,"/","all-fast","_",save.type,".bat"),
              sep = "\t",quote = F,
              col.names = F,row.names = F)
  LuckyVerbose("Complete all-md-fast",levels = 2)


  ## Out put
  l <- list(
    script = PE,
    md = md.files,
    csl = csl_file,
    bib = bib.file,
    save.path = save.path
  )
  LuckyVerbose("Done!")
  return(l)
}


