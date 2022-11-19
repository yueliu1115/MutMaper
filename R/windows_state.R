#' Title Compute SNP-index and ED value in based on fixed Windows in Chromosome
#'
#' @param vcf_snpindex Input a vcf object (dataframe)
#' @param ED Whether to calculate ED
#' @param windows Set fixed Windows in Chromosome (bp)
#'
#' @return return a dataframe per fixed windows in Chromosome
#' @export
#'
#' @examples
windows_state <- function(vcf_snpindex, ED = c(FALSE, TRUE),windows){
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

  if (isTRUE(ED)) {
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
  }

  window_mutmap <- na.omit(window_mutmap)

  write.csv(window_mutmap,
            file = paste0("windows_state_",Sys.Date(),".csv"))

  return(window_mutmap)
}
