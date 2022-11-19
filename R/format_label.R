format_label <- function(...) {
  function(x) {
    limits <- c(1e0,   1e3, 1e6)

    i <- findInterval(abs(x), limits)


    i <- ifelse(i==0, which(limits == 1e0), i)

    paste(format(round(x/limits[i], 1),
                 trim=TRUE, scientific=FALSE, ...)
    )
  }
}
