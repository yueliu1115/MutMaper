#' Title Onestep analysis mutmap based on vcf file and fixed windows
#'
#' @param vcf Input a vcf file path
#' @param windows Set fixed Windows in Chromosome (bp)
#'
#' @return A plot out SNP-index and ED based on Mutmap
#' @export
#'
#' @importFrom ggplot2 ggplot aes facet_grid ylab geom_smooth scale_x_continuous
#' @examples
onestep_mutmap <- function(vcf, windows){
  vcffile <- readr::read_delim(file = vcf, col_names = TRUE,
                               delim = "\t", comment = "##") %>%
    dplyr::rename("CHROM" = "#CHROM") %>%
    purrr::set_names(c(colnames(.)[1:9], "Wild_Type", "Mutant"))

  vcffile <- vcffile %>%
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

  vcf_snpindex <- vcffile %>%
    dplyr::mutate(mutant_snpindex = mutant_AD_alt / mutant_DP) %>%
    dplyr::mutate(control_snpindex = control_AD_alt / control_DP) %>%
    dplyr::mutate(snpindex = mutant_snpindex - control_snpindex) %>%
    dplyr::filter(snpindex > 0.2) %>%
    dplyr::mutate(ED = ((2*(mutant_snpindex - control_snpindex)^2)^(1/2))) %>%
    dplyr::mutate(ED2 = ((2*(mutant_snpindex - control_snpindex)^2)^(1/2))^2) %>%
    dplyr::mutate(ED6 = ((2*(mutant_snpindex - control_snpindex)^2)^(1/2))^6)
  write.csv(vcf_snpindex,
            file = paste0("snpindex_",Sys.Date(),".csv"))

  #stat chrom length
  chrom_length <- c()
  for (i in 1:length(table(vcf_snpindex$CHROM))) {
    length_chr <- max(vcf_snpindex[vcf_snpindex$CHROM == i,]$POS)
    chrom_length <- c(chrom_length, length_chr)
  }
  chrom_length <- chrom_length[!is.infinite(chrom_length)]
  #create windows
  win_start <- c()
  win_end <- c()
  win_chr <- c()
  win_skip <- windows

  for (i in 1:length(chrom_length)) {
    win_start_tmp <- seq(0,chrom_length[i],by=win_skip)
    win_end_tmp <- win_start_tmp + win_skip
    chr_tmp <- rep(i, length(win_start_tmp))
    win_start <- c(win_start, win_start_tmp)
    win_end <- c(win_end, win_end_tmp)
    win_chr <- c(win_chr,chr_tmp)
  }

  #make windows
  window_mutmap <- data.frame(
    chr = win_chr,
    start = win_start,
    end = win_end
  )

  #stat number of snp
  window_mutmap$nsnp <- apply(window_mutmap,1,function(x){
    tmp <- vcf_snpindex[vcf_snpindex$CHROM == x[1] &
                          vcf_snpindex$POS > x[2] &
                          vcf_snpindex$POS < x[3],]
    return(nrow(tmp))
  })

  #stat snp_index
  window_mutmap$snp_index <- apply(window_mutmap,1,function(x){
    tmp <- vcf_snpindex[vcf_snpindex$CHROM == x[1] &
                          vcf_snpindex$POS > x[2] &
                          vcf_snpindex$POS < x[3],]
    return(mean(tmp$snpindex))
  })

  #stat ED
  window_mutmap$ED <- apply(window_mutmap,1,function(x){
    tmp <- vcf_snpindex[vcf_snpindex$CHROM == x[1] &
                          vcf_snpindex$POS > x[2] &
                          vcf_snpindex$POS < x[3],]
    return(mean(tmp$ED))
  })

  #stat ED^2
  window_mutmap$ED2 <- apply(window_mutmap,1,function(x){
    tmp <- vcf_snpindex[vcf_snpindex$CHROM == x[1] &
                          vcf_snpindex$POS > x[2] &
                          vcf_snpindex$POS < x[3],]
    return(mean(tmp$ED2))
  })

  #stat ED^6
  window_mutmap$ED6 <- apply(window_mutmap,1,function(x){
    tmp <- vcf_snpindex[vcf_snpindex$CHROM == x[1] &
                          vcf_snpindex$POS > x[2] &
                          vcf_snpindex$POS < x[3],]
    return(mean(tmp$ED6))
  })


  windows_state <- na.omit(window_mutmap)

  write.csv(windows_state,
            file = paste0("windows_state_",Sys.Date(),".csv"))

  p_snp <- ggplot2::ggplot(data = windows_state) +
    #geom_point(aes(x=start, y=nsnp)) +
    ggplot2::geom_smooth(ggplot2::aes(x = start, y=nsnp,color = factor(chr)),stat = "identity") +
    ggplot2::scale_x_continuous(seq(0,max(windows_state$start),
                           by = 10^(floor(log10(max(windows_state$start))))),
                       labels = format_label(),
                       name = "Genomic Position(Mb)") +
    ggplot2::ylab("nSNPs/Mb") + ggplot2::theme(legend.position="none")+
    ggplot2::facet_grid(~chr,scales = "free_x",space = "free_x")

  p_snp_index <- ggplot2::ggplot(data = windows_state) +
    #geom_point(aes(x=start, y=snp_index)) +
    ggplot2::geom_smooth(ggplot2::aes(x = start, y=snp_index,color = factor(chr)),stat = "identity") +
    ggplot2::facet_grid(~chr, scales = "free") +
    ggplot2::scale_x_continuous(seq(0,max(windows_state$start),
                           by = 10^(floor(log10(max(windows_state$start))))),
                       labels = format_label(),
                       name = "Genomic Position(Mb)") +
    ggplot2::ylab("SNP-index") +
    ggplot2::ylim(c(0,1))+
    ggplot2::theme_bw()+
    ggplot2::geom_hline(ggplot2::aes(yintercept=0.5),color = 'red',size=1,linetype="dashed")+
    ggplot2::theme(legend.position="none") +
    ggplot2::facet_grid(~chr,scales = "free_x",space = "free_x")

  p_ed <- ggplot2::ggplot(data = windows_state) +
    #geom_point(aes(x=start, y=ED6)) +
    ggplot2::geom_smooth(ggplot2::aes(x = start, y=ED,color = factor(chr)),stat = "identity")+
    ggplot2::scale_x_continuous(seq(0,max(windows_state$start),
                           by = 10^(floor(log10(max(windows_state$start))))),
                       labels = format_label(),
                       name = "Genomic Position(Mb)") +
    ggplot2::ylab("ED")+
    #geom_hline(aes(yintercept=3),color = 'red',size=1,linetype="dashed") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position="none") +
    ggplot2::facet_grid(~chr,scales = "free_x",space = "free_x")

  p_ed2 <- ggplot2::ggplot(data = windows_state) +
    #geom_point(aes(x=start, y=ED6)) +
    ggplot2::geom_smooth(ggplot2::aes(x = start, y=ED2,color = factor(chr)),stat = "identity")+
    ggplot2::scale_x_continuous(seq(0,max(windows_state$start),
                           by = 10^(floor(log10(max(windows_state$start))))),
                       labels = format_label(),
                       name = "Genomic Position(Mb)") +
    ggplot2::ylab("ED^2")+
    #geom_hline(aes(yintercept=3),color = 'red',size=1,linetype="dashed") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position="none") +
    ggplot2::facet_grid(~chr,scales = "free_x",space = "free_x")

  p_ed6 <- ggplot2::ggplot(data = window_mutmap) +
    #geom_point(aes(x=start, y=ED6)) +
    ggplot2::geom_smooth(ggplot2::aes(x = start, y=ED6,color = factor(chr)),stat = "identity")+
    ggplot2::scale_x_continuous(seq(0,max(window_mutmap$start),
                           by = 10^(floor(log10(max(window_mutmap$start))))),
                       labels = format_label(),
                       name = "Genomic Position(Mb)") +
    ggplot2::ylab("ED^6")+
    #geom_hline(aes(yintercept=3),color = 'red',size=1,linetype="dashed") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position="none") +
    ggplot2::facet_grid(~chr,scales = "free_x",space = "free_x")

  p <- cowplot::plot_grid(
    p_snp,
    p_snp_index,
    p_ed,
    p_ed2,
    p_ed6,
    nrow = 5,
    align = "v"
  )
  return(p)
}
