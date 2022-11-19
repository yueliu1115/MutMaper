#' Title Compute SNP-index and ED value in Mutmap
#'
#' @param vcf Input a vcf object (dataframe)
#' @param snpindex Whether to calculate SNP-index, default value is TRUE
#' @param ED Whether to calculate ED
#'
#' @return return a vcf dataframe with SNP-index | ED value
#' @export
#'
#' @examples
vcf_snpindex <- function(vcf, snpindex = TRUE, ED = c(FALSE, TRUE)){
  message("Please Make sure that the input vcf is clean vcf object!")
  if (snpindex == TRUE & isFALSE(ED)) {
    snpindex_out <- vcf %>%
      dplyr::mutate(mutant_snpindex = mutant_AD_alt / mutant_DP) %>%
      dplyr::mutate(control_snpindex = control_AD_alt / control_DP) %>%
      dplyr::mutate(snpindex = mutant_snpindex - control_snpindex) %>%
      dplyr::filter(snpindex > 0.2)
    write.csv(snpindex_out,
              file = paste0("snpindex_",Sys.Date(),".csv"))
    return(snpindex_out)
  }

  if (snpindex == TRUE & isTRUE(ED)) {
    snpindex_out <- vcf %>%
      dplyr::mutate(mutant_snpindex = mutant_AD_alt / mutant_DP) %>%
      dplyr::mutate(control_snpindex = control_AD_alt / control_DP) %>%
      dplyr::mutate(snpindex = mutant_snpindex - control_snpindex) %>%
      dplyr::filter(snpindex > 0.2) %>%
      dplyr::mutate(ED = ((2*(mutant_snpindex - control_snpindex)^2)^(1/2))) %>%
      dplyr::mutate(ED2 = ((2*(mutant_snpindex - control_snpindex)^2)^(1/2))^2) %>%
      dplyr::mutate(ED6 = ((2*(mutant_snpindex - control_snpindex)^2)^(1/2))^6)
    write.csv(snpindex_out,
              file = paste0("snpindex_",Sys.Date(),".csv"))
    return(snpindex_out)
  }
}
