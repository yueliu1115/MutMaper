% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vcf_snpindex.R
\name{vcf_snpindex}
\alias{vcf_snpindex}
\title{Title Compute SNP-index and ED value in Mutmap}
\usage{
vcf_snpindex(vcf, snpindex = TRUE, ED = c(FALSE, TRUE))
}
\arguments{
\item{vcf}{Input a vcf object (dataframe)}

\item{snpindex}{Whether to calculate SNP-index, default value is TRUE}

\item{ED}{Whether to calculate ED}
}
\value{
return a vcf dataframe with SNP-index | ED value
}
\description{
Title Compute SNP-index and ED value in Mutmap
}
