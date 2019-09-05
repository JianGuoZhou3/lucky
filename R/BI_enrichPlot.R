







####=============Assistant Functions==================####


enrich_goplot <- function(GO,
                          showCategory = list(CC=20,MF=10,BP=12),
                          save.path = "GO",
                          names = "test.enrichment"){
  
  ## 产生储存文件夹
  old <- options()
  dir <- paste0("./",names,"/",save.path,"/")
  options(warn = -1)
  dir.create(dir,recursive = T)
  options(warn = old$warn)
  
  ## GO category
  Go.category <- c("CC","MF","BP")

  ##goplot
  pdf.pathway <- str_c(dir,names,"_Goplot.pdf");
  pdf(pdf.pathway,width = 12,height = 8)
  for(i in 1:3){ # i=1
    ## Data preparation
    c.i <-  Go.category[i]
    go.enrich.i <- GO[["GOEnrichment"]][[c.i]]
    nego <- nrow(as.data.frame(go.enrich.i))
    sc <- showCategory[[c.i]]
    sc2 <- ifelse(sc >= nego,nego,sc)
    
    ## goplot
    if(sc2 > 2){
      g1 <- goplot(go.enrich.i,
                   showCategory = sc2) + 
        ggtitle(c.i) + 
        theme(plot.title = element_text(hjust = 0.5))
    }
    print(g1)
  }
  dev.off()
  LuckyVerbose(paste0("GOplot have been saved!"))
}
enrich_go_dotplot <- function(GO,
                              x = "GeneRatio",
                              color = "p.adjust", 
                              showCategory = list(CC=20,MF=10,BP=12),
                              font.size = 12,
                              title = "",
                              save.path = "GO",
                              names = "test.enrichment"){
  
  ## 产生储存文件夹
  old <- options()
  dir <- paste0("./",names,"/",save.path,"/")
  options(warn = -1)
  dir.create(dir,recursive = T)
  options(warn = old$warn)
  
  ## GO category
  Go.category <- c("CC","MF","BP")
  
  ##goplot
  pdf.pathway <- str_c(dir,names,"_dotplot_go enrichment.pdf");
  pdf(pdf.pathway,width = 12,height = 8)
  for(i in 1:3){ # i=1
    ## Data preparation
    c.i <-  Go.category[i]
    go.enrich.i <- GO[["GOEnrichment"]][[c.i]]
    nego <- nrow(as.data.frame(go.enrich.i))
    sc <- showCategory[[c.i]]
    sc2 <- ifelse(sc >= nego,nego,sc)
    
    ## goplot
    if(sc2 > 2){
      g1 <- dotplot(go.enrich.i, 
                    x = x,
                    color = color, 
                    showCategory = sc2,
                    font.size = font.size,
                    title = title)
      print(g1)
    }
    
  }
  dev.off()
  LuckyVerbose(paste0("Doplot enriched terms have been saved!"))
  
}













