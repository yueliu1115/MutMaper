#' Title Read vcf file and translate into clean vcf
#'
#' @param vcf Input a vcf dataframe or a vcf file path
#' @param filter  FALSE is not to filter vcf file, "EMS" is to select G->A and C->T
#'
#' @return retun a clean vcf dataframe
#' @export
#'
#' @examples
clean_vcf <- function(vcf, filter = c(FALSE, "EMS")){
  if (is.data.frame(vcf) & filter == FALSE) {
    vcffile <- vcf %>%
      dplyr::select(1,2,4,5,9,10,11) %>%
      tidyr::separate(Wild_Type,into = paste("control",c("GT", "AD","DP"),sep = "_"),
                      sep=":",remove= TRUE) %>%
      tidyr::separate(Mutant,into = paste("mutant",c("GT","AD","DP"),sep = "_"),
                      sep=":",remove= TRUE) %>%
      dplyr::mutate(mutant_GT = str_replace(mutant_GT,pattern = "/",replacement = "|")) %>%
      dplyr::mutate(control_GT = str_replace(control_GT,pattern = "/",replacement = "|")) %>%
      dplyr::filter(!(control_GT == "1|1" & mutant_GT == "1|1")) %>%
      dplyr::mutate(control_DP = as.numeric(control_DP)) %>%
      dplyr::mutate(mutant_DP = as.numeric(mutant_DP)) %>%
      dplyr::filter(mutant_DP > 10,mutant_DP < 300, control_DP > 10,control_DP < 300) %>%
      tidyr::separate(mutant_AD,into = c("mutant_AD_ref","mutant_AD_alt")) %>%
      tidyr::separate(control_AD,into = c("control_AD_ref",'control_AD_alt')) %>%
      dplyr::mutate(mutant_AD_ref = as.numeric(mutant_AD_ref)) %>%
      dplyr::mutate(mutant_AD_alt = as.numeric(mutant_AD_alt)) %>%
      dplyr::mutate(control_AD_ref = as.numeric(control_AD_ref)) %>%
      dplyr::mutate(control_AD_alt = as.numeric(control_AD_alt))
    write.csv(vcffile,
              file = paste0("clean_vcf_",Sys.Date(),".csv"))
    return(vcffile)
  }

  if (is.data.frame(vcf) & filter == "EMS") {
    vcffile <- vcf %>%
      dplyr::select(1,2,4,5,9,10,11) %>%
      tidyr::separate(Wild_Type,into = paste("control",c("GT", "AD","DP"),sep = "_"),
                      sep=":",remove= TRUE) %>%
      tidyr::separate(Mutant,into = paste("mutant",c("GT","AD","DP"),sep = "_"),
                      sep=":",remove= TRUE) %>%
      dplyr::mutate(mutant_GT = stringr::str_replace(mutant_GT,pattern = "/",replacement = "|")) %>%
      dplyr::mutate(control_GT = stringr::str_replace(control_GT,pattern = "/",replacement = "|")) %>%
      dplyr::filter(!(control_GT == "1|1" & mutant_GT == "1|1")) %>%
      dplyr::mutate(control_DP = as.numeric(control_DP)) %>%
      dplyr::mutate(mutant_DP = as.numeric(mutant_DP)) %>%
      dplyr::filter(mutant_DP > 10,mutant_DP < 300, control_DP > 10,control_DP < 300) %>%
      dplyr::filter((REF == "G" & ALT == "A") | (REF == "C" & ALT == "T")) %>%
      tidyr::separate(mutant_AD,into = c("mutant_AD_ref","mutant_AD_alt")) %>%
      tidyr::separate(control_AD,into = c("control_AD_ref",'control_AD_alt')) %>%
      dplyr::mutate(mutant_AD_ref = as.numeric(mutant_AD_ref)) %>%
      dplyr::mutate(mutant_AD_alt = as.numeric(mutant_AD_alt)) %>%
      dplyr::mutate(control_AD_ref = as.numeric(control_AD_ref)) %>%
      dplyr::mutate(control_AD_alt = as.numeric(control_AD_alt))
    write.csv(vcffile,
              file = paste0("clean_vcf_",Sys.Date(),".csv"))
    return(vcffile)
  }

  if (is.character(vcf) & filter == FALSE) {
    vcffile <- readr::read_delim(file = vcf, col_names = TRUE,
                                 delim = "\t", comment = "##") %>%
      dplyr::rename("CHROM" = "#CHROM") %>%
      dplyr::select(1,2,4,5,9,10,11) %>%
      tidyr::separate(Wild_Type,into = paste("control",c("GT", "AD","DP"),sep = "_"),
                      sep=":",remove= TRUE) %>%
      tidyr::separate(Mutant,into = paste("mutant",c("GT","AD","DP"),sep = "_"),
                      sep=":",remove= TRUE) %>%
      dplyr::mutate(mutant_GT = stringr::str_replace(mutant_GT,pattern = "/",replacement = "|")) %>%
      dplyr::mutate(control_GT = stringr::str_replace(control_GT,pattern = "/",replacement = "|")) %>%
      dplyr::filter(!(control_GT == "1|1" & mutant_GT == "1|1")) %>%
      dplyr::mutate(control_DP = as.numeric(control_DP)) %>%
      dplyr::mutate(mutant_DP = as.numeric(mutant_DP)) %>%
      dplyr::filter(mutant_DP > 10,mutant_DP < 300, control_DP > 10,control_DP < 300) %>%
      tidyr::separate(mutant_AD,into = c("mutant_AD_ref","mutant_AD_alt")) %>%
      tidyr::separate(control_AD,into = c("control_AD_ref",'control_AD_alt')) %>%
      dplyr::mutate(mutant_AD_ref = as.numeric(mutant_AD_ref)) %>%
      dplyr::mutate(mutant_AD_alt = as.numeric(mutant_AD_alt)) %>%
      dplyr::mutate(control_AD_ref = as.numeric(control_AD_ref)) %>%
      dplyr::mutate(control_AD_alt = as.numeric(control_AD_alt))
    write.csv(vcffile,
              file = paste0("clean_vcf_",Sys.Date(),".csv"))
    return(vcffile)
  }

  if (is.character(vcf) & filter == "EMS") {
    vcffile <- readr::read_delim(file = vcf, col_names = TRUE,
                                 delim = "\t", comment = "##") %>%
      dplyr::rename("CHROM" = "#CHROM") %>%
      dplyr::select(1,2,4,5,9,10,11) %>%
      tidyr::separate(Wild_Type,into = paste("control",c("GT", "AD","DP"),sep = "_"),
                      sep=":",remove= TRUE) %>%
      tidyr::separate(Mutant,into = paste("mutant",c("GT","AD","DP"),sep = "_"),
                      sep=":",remove= TRUE) %>%
      dplyr::mutate(mutant_GT = stringr::str_replace(mutant_GT,pattern = "/",replacement = "|")) %>%
      dplyr::mutate(control_GT = stringr::str_replace(control_GT,pattern = "/",replacement = "|")) %>%
      dplyr::filter(!(control_GT == "1|1" & mutant_GT == "1|1")) %>%
      dplyr::mutate(control_DP = as.numeric(control_DP)) %>%
      dplyr::mutate(mutant_DP = as.numeric(mutant_DP)) %>%
      dplyr::filter(mutant_DP > 10,mutant_DP < 300, control_DP > 10,control_DP < 300) %>%
      dplyr::filter((REF == "G" & ALT == "A") | (REF == "C" & ALT == "T")) %>%
      tidyr::separate(mutant_AD,into = c("mutant_AD_ref","mutant_AD_alt")) %>%
      tidyr::separate(control_AD,into = c("control_AD_ref",'control_AD_alt')) %>%
      dplyr::mutate(mutant_AD_ref = as.numeric(mutant_AD_ref)) %>%
      dplyr::mutate(mutant_AD_alt = as.numeric(mutant_AD_alt)) %>%
      dplyr::mutate(control_AD_ref = as.numeric(control_AD_ref)) %>%
      dplyr::mutate(control_AD_alt = as.numeric(control_AD_alt))
    write.csv(vcffile,
              file = paste0("clean_vcf_",Sys.Date(),".csv"))
    return(vcffile)
  }
}
