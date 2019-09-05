## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----results = F---------------------------------------------------------
library(lucky)
library(ggpubr)
# Data
df <- data.frame(dose=c("D0.5", "D1", "D2"),
                 len=c(4.2, 10, 29.5))
print(df)

## ----fig.show = "hold"---------------------------------------------------
# draw a barplot
FastBar(data = df,
        x = "dose",y = "len",
        fill = "dose",
        title = "A test from ggpubr",
        legend.position = "right",
        size = 15)

## ---- fig.show='hold'----------------------------------------------------
plot(1:10)
plot(10:1)

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(head(mtcars, 10))

