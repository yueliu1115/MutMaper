#' Title Read vcf file in R
#'
#' @param vcf vcf file is a data frame from GATK output
#' @param header The header of vcf file, the default value is TRUE
#' @param sep The delimter character of vcf file, the default value is tab
#' @param comment The comment character of vcf file, the default value is "##"
#'
#' @return The output is the raw vcf file
#' @export
#'
#' @examples \dontrun{vcf <- read_vcf(vcf = "./all.filter.snp.vcf")}
read_vcf <- function(vcf, header = TRUE, sep = "\t", comment = "##"){
  if (is.null(vcf)) {
    stop("The vcf file is missing! Please provide the vcf file base on GATK output")
  }
  vcffile <- readr::read_delim(file = vcf, col_names = header,
                               delim = sep, comment = comment) %>%
    dplyr::rename("CHROM" = "#CHROM") %>%
    purrr::set_names(c(colnames(.)[1:9], "Wild_Type", "Mutant"))
  return(vcffile)
}
