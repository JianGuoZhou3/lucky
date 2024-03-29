% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BI_FastCeRNA.R
\name{FastCeRNA}
\alias{FastCeRNA}
\title{easy ceRNA network via GDCRNATools::gdcCEAnalysis}
\usage{
FastCeRNA(gene.counts, miRNA.counts, matrix.type = "counts",
  select = NULL, design, deMIR = NULL, lnc.targets = "starBase",
  pc.targets = "starBase", save.path = "ceRNA", names = "love")
}
\arguments{
\item{gene.counts}{the raw counts matrix of mRNA and lncRNA}

\item{miRNA.counts}{the raw counts matrix of miRNA}

\item{matrix.type}{the type of expression matrix.Default is "counts".If you have a normalized matrix,please set it "normalized".}

\item{design}{a design object}

\item{deMIR}{differencial expression miRNAs}

\item{lnc.targets}{a list with lcnRNA ENSEMBL names and miRNAs symbol elements.Or "starBase"}

\item{pc.targets}{a list with protein-coding RNA(mRNA) ENSEMBL names and miRNAs symbol elements.}

\item{save.path}{the space of the save file.Default is "WGCNA"}

\item{names}{part of saved files name}
}
\description{
easy ceRNA network via GDCRNATools::gdcCEAnalysis
}
\seealso{
\code{\link[GDCRNATools]{gdcCEAnalysis}}.
}
\author{
Weibin Huang<\email{654751191@qq.com}>
}
