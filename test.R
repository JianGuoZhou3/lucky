

library(lucky)
load("E:/RCloud/database/DataDownload/TCGA-TCGABiolinks/STAD_Download/design_STAD.rda")
load("E:/RCloud/database/DataDownload/TCGA-TCGABiolinks/STAD_Download/stad_FPKM.rda")


expr.matrix = STAD.FPKM
design = design[design$condition %in% "tumor",]
select = convert1("HNRNPA2B1")
event.status = "OS.status"
event.time = "OS.time"
event.lower = 0
mode = c("median","maxstat")[2]
pval.position =NULL
size = 12
save.file = T
width = 12
height = 12
names = "love"
print = T



p <- Survplot_genedata(expr.matrix = STAD.FPKM,
                       design= design[design$condition %in% "tumor",],
                       select=convert1("HNRNPA2B1"),
                       event.status= "OS.status",
                       event.time="OS.time",
                       event.lower = 0,
                       mode = c("median", "maxstat")[2],
                       pval.position = c(3500,1),
                       size = 12)


