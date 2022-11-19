#' Title Plot SNP-index in Mutmap
#'
#' @param windows_state Input a dataframe with snpindex based on windows
#' @param ED Whether the windows_state contain ED value
#'
#' @return a Plot with SNP-index
#' @export
#'
#' @importFrom ggplot2 ggplot aes facet_grid ylab geom_smooth scale_x_continuous
#' @examples
plot_mutmap <- function(windows_state, ED = c(FALSE, TRUE)){
  if (is.vector(windows_state)) {
    windows_state <- utils::read.csv(file = windows_state,header = T)
    p_snp <- ggplot2::ggplot(data = windows_state) +
      #geom_point(aes(x=start, y=nsnp)) +
      ggplot2::geom_smooth(ggplot2::aes(x = start, y=nsnp,color = factor(chr)),stat = "identity") +
      ggplot2::scale_x_continuous(seq(0,max(windows_state$start),
                             by = 10^(floor(log10(max(windows_state$start))))),
                         labels = format_label(),
                         name = "Genomic Position(Mb)") +
      ggplot2::ylab("nSNPs/Mb") + theme(legend.position="none")+
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

    if (isTRUE(ED)) {
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
    }else{
      p <- cowplot::plot_grid(
        p_snp,
        p_snp_index,
        nrow = 2,
        align = "v"
      )
      return(p)
    }
  }

  if (is.data.frame(windows_state)) {
    p_snp <- ggplot2::ggplot(data = windows_state) +
      #geom_point(aes(x=start, y=nsnp)) +
      ggplot2::geom_smooth(ggplot2::aes(x = start, y=nsnp,color = factor(chr)),stat = "identity") +
      ggplot2::scale_x_continuous(seq(0,max(windows_state$start),
                             by = 10^(floor(log10(max(windows_state$start))))),
                         labels = format_label(),
                         name = "Genomic Position(Mb)") +
      ggplot2::ylab("nSNPs/Mb") +
      ggplot2::theme(legend.position="none")+
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

    if (isTRUE(ED)) {
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

      p_ed6 <- ggplot2::ggplot(data = windows_state) +
        #geom_point(aes(x=start, y=ED6)) +
        ggplot2:: geom_smooth(ggplot2::aes(x = start, y=ED6,color = factor(chr)),stat = "identity")+
        ggplot2::scale_x_continuous(seq(0,max(windows_state$start),
                               by = 10^(floor(log10(max(windows_state$start))))),
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
    }else{
      p <- cowplot::plot_grid(
        p_snp,
        p_snp_index,
        nrow = 2,
        align = "v"
      )
      return(p)
    }
  }
}
