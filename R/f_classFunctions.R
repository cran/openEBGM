#' Print an openEBGM object
#'
#' @param x An openEBGM object constructed by \code{ebScores()}
#' @param threshold A numeric value indicating the minimum threshold
#'   for QUANT or EBGM values to display.
#' @param ... Arguments to be passed to other methods
#' @export
print.openEBGM <- function(x, threshold = 2, ...) {
  if(!is.numeric(threshold)) {
    stop("argument 'threshold' must be numeric")
  }
  if(any(grepl("QUANT", names(x$data)))) {
    quant_nums <- grep("QUANT", names(x$data), value = TRUE)
    quant_nums <- substring(quant_nums, first = 7, last = 1000000L)
    quant_min <- min(as.numeric(quant_nums))
    if(nchar(quant_min) == 1) {
      quant_min <- paste("0", quant_min, sep = "")
    }
    quant_sub <- x$data[,grep(quant_min, names(x$data))]
    cat("\n")
    cat("There were "); cat(length(quant_sub[quant_sub > threshold]))
    cat(" var1-var2 pairs with a ")
    cat(paste("QUANT", quant_min, sep = "_"))
    cat(" greater than ")
    cat(as.character(threshold))
    cat("\n")
    cat("\nTop 5 Highest ")
    cat(paste("QUANT", quant_min, sep = "_"))
    cat(" Scores\n")
    tmp <- x$data
    tmp <- head(tmp[order(quant_sub, decreasing = TRUE),], 5)
    if(any(nchar(as.character(tmp$var1)) > 50)) {
      tmp$var1 <- strtrim(as.character(tmp$var1), 50)
    }
    print(tmp[,c("var1", "var2", "N", "E", grep(quant_min, names(tmp),
                                                value = TRUE))])
    cat("\n")
  } else {
    cat("\n")
    cat("There were "); cat(length(x$data$EBGM[x$data$EBGM > threshold]))
    cat(" var1-var2 pairs with an EBGM score greater than ")
    cat(as.character(threshold))
    cat("\n")
    cat("\nTop 5 Highest EBGM Scores\n")
    tmp <- x$data
    tmp <- head(tmp[order(tmp$EBGM, decreasing = TRUE),], 5)
    if(any(nchar(as.character(tmp$var1)) > 50)) {
      tmp$var1 <- strtrim(as.character(tmp$var1), 50)
    }
    print(tmp[,c("var1", "var2", "N", "E", "EBGM")])
    cat("\n")
  }
}

#' Plot an openEBGM object
#'
#' @param x An openEBGM object constructed by \code{ebScores()}
#' @param y Unused parameter to satisfy generic function requirement
#' @param event An (optional) specification of an event to subset the data by.
#' @param plot.type A character vector specifying which type of plot should be
#'                  output. See details.
#' @param ... Arguments to be passed to methods
#'
#' @details There are three different types of plots that the plot function may
#'          produce when called on an openEBGM object. These are
#'          \itemize{
#'            \item bar
#'            \item shrinkage
#'            \item histogram
#' }
#'
#' @details A bar chart displays the top ten product-symptom EBGM scores, as
#'   well as error bars which display the highest and lowest of the quantiles
#'   chosen at the time of instantiating the openEBGM object. A shrinkage plot
#'   plots EBGM score on the y axis, and the natural log of the RR on the x
#'   axis. This plot is also called a squid plot and was conceived by Stuart
#'   Chirtel. Finally, a histogram simply displays a histogram of the EBGM
#'   scores.
#' @import ggplot2
#' @export
plot.openEBGM <- function(x, y = NULL, event = NULL, plot.type = "bar", ...) {
  tmp <- x$data
  if(plot.type == "bar") {
    if(!is.null(event)) {
      if(!any(grepl(event, tmp$var2, ignore.case = TRUE))) {
        stop("'event' does not match any events in var2. Are you sure this is an
             event in the data?")
      }
      tmp <- tmp[grep(event, tmp$var2, ignore.case = TRUE),]
      if(length(unique(tmp$var2)) > 1)
        warning("2 or more matches found for event specified")
      }
    tmp$EBGM_fuzzy <- ifelse(tmp$EBGM > 8, "8+",
                      ifelse(tmp$EBGM <=8 & tmp$EBGM > 4, "4-8",
                      ifelse(tmp$EBGM <= 4 & tmp$EBGM > 2, "2-4",
                      ifelse(tmp$EBGM <= 2 & tmp$EBGM > 1, "1-2", "<1"))))
    tmp$fac <- ifelse(tmp$EBGM > 8, "red4",
               ifelse(tmp$EBGM <= 8 & tmp$EBGM > 4, "orangered1",
               ifelse(tmp$EBGM <= 4 & tmp$EBGM > 2, "orange",
               ifelse(tmp$EBGM <= 2 & tmp$EBGM > 1, "gold", "yellow"))))
    tmp <- head(tmp[order(tmp$EBGM, decreasing = TRUE),], 15)
    tmp$var1var2 <- paste(tmp$var1, tmp$var2, sep = "_")
    g <- ggplot(aes(x = var1var2, y = EBGM, fill = EBGM_fuzzy), data = tmp)
    if(sum(grepl("QUANT", names(tmp))) >= 2) {
      quant_nums <- grep("QUANT", names(tmp), value = TRUE)
      quant_nums <- substring(quant_nums, first = 7, last = 100000000L)
      quant_min <- min(as.numeric(quant_nums))
      if(nchar(quant_min) == 1) {
        quant_min <- paste("0", quant_min, sep = "")
      }
      quant_max <- max(as.numeric(quant_nums))
      if(nchar(quant_max) == 1) {
        quant_max <- paste("0", quant_max, sep = "")
      }
      quant_min_name <- paste("QUANT", quant_min, sep = "_")
      quant_max_name <- paste("QUANT", quant_max, sep = "_")
      g + geom_bar(stat = "identity") +
        coord_flip() +
        scale_x_discrete(limits = tmp[order(tmp$EBGM, decreasing = FALSE), "var1var2"],
                         labels = strtrim(as.character(tmp[order(tmp$EBGM, decreasing = FALSE), "var1"]), 30)) +
        scale_y_continuous(limits = c(0, max(tmp[, grep(quant_max,
                                                        names(tmp))]) + 5),
                           expand = c(0, 0)) +
        geom_errorbar(aes(ymin = tmp[order(tmp$EBGM, decreasing = TRUE),
                                     grep(quant_min, names(tmp))],
                          ymax = tmp[order(tmp$EBGM, decreasing = TRUE),
                                     grep(quant_max, names(tmp))]),
                      position = "dodge",
                      width = 0.25) +
        geom_text(aes(label = paste("N = ", tmp[order(tmp$EBGM, decreasing = TRUE), "N"],
                                    sep  = ""),
                      y = tmp[order(tmp$EBGM, decreasing = TRUE),
                              grep(quant_max, names(tmp))]),
                  hjust = -0.25) +
        scale_fill_manual(values = c("8+" = "red4", "4-8" = "orangered1",
                                     "2-4" = "orange", "1-2" = "gold",
                                     "<1"= "yellow"),
                          name = "EBGM") +
        theme_bw() +
        ylab(paste(quant_min_name, quant_max_name, sep = " - EBGM - ")) +
        xlab("var1 observation") +
        ggtitle(ifelse(is.null(event), "EBGM Barplot", paste("EBGM Barplot with Event=", event, sep = ""))) +
        #Note that themes are applied *after* the coord_flip()
        theme(axis.title.x = element_text(face = "bold"),
              axis.title.y = element_text(face = "bold"),
              panel.grid = element_blank(),
              panel.border = element_blank(),
              axis.line = element_line(color = "black"),
              panel.spacing = unit(0, units = "cm"),
              plot.title = element_text(face = "bold"),
              legend.position = "bottom",
              legend.box = "horizontal")
    } else {
      g + geom_bar(stat = "identity") +
        coord_flip() +
        scale_x_discrete(limits = tmp[order(tmp$EBGM, decreasing = FALSE),"var1var2"],
                         labels = strtrim(as.character(tmp[order(tmp$EBGM, decreasing = FALSE), "var1"]), 30)) +
        scale_y_continuous(limits = c(0, max(tmp$EBGM) + 5), expand = c(0, 0)) +
        geom_text(aes(label = paste("N = ", tmp[order(tmp$EBGM, decreasing = TRUE), "N"], sep  = ""),
                      y = tmp[order(tmp$EBGM, decreasing = TRUE), "EBGM"]), hjust = -0.25) +
        scale_fill_manual(values = c("8+" = "red4", "4-8" = "orangered1",
                                     "2-4" = "orange", "1-2" = "gold",
                                     "<1"= "yellow"),
                          name = "EBGM") +
        theme_bw() +
        ylab("EBGM") + xlab("var1 observation") +
        ggtitle(ifelse(is.null(event), "EBGM Barplot", paste("EBGM Barplot with Event=", event, sep = ""))) +
        #Note that themes are applied *after* the coord_flip()
        theme(axis.title.x = element_text(face = "bold"),
              axis.title.y = element_text(face = "bold"),
              panel.grid = element_blank(),
              panel.border = element_blank(),
              axis.line = element_line(color = "black"),
              panel.spacing = unit(0, units = "cm"),
              plot.title = element_text(face = "bold"),
              legend.position = "bottom",
              legend.box = "horizontal")
    }
  } else if(plot.type == "shrinkage") {
    if(!is.null(event)) {
      tmp <- tmp[tmp$var2 == event,]
    }
    #Need to get the N column into colors
    tmp$N_col <- ifelse(tmp$N == 1, 1,
                 ifelse(tmp$N == 2, 2,
                 ifelse(tmp$N == 3, 3,
                 ifelse(tmp$N == 4, 4,
                 ifelse(tmp$N == 5, 5, "6+")))))
    g <- ggplot(aes(x = log(tmp$RR), y = tmp$EBGM, color = as.factor(N_col)),
                data = tmp)
    g + geom_point() +
      ggtitle("Chirtel Squid Shrinkage Plot") +
      xlab("ln(RR)") + ylab("EBGM") +
      theme_bw() +
      guides(color = guide_legend("N")) +
      theme(axis.title.x = element_text(face = "bold"),
            axis.title.y = element_text(face = "bold"),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            axis.line = element_line(color = "black"),
            panel.spacing = unit(0, units = "cm"),
            plot.title = element_text(face = "bold"))
  } else if(plot.type == "histogram") {
    if(!is.null(event)) {
      tmp <- tmp[tmp$var2 == event,]
    }
    g <- ggplot(aes(x = EBGM), data = tmp)
    g + geom_histogram(fill = "black") +
      ggtitle("Histogram of EBGM Scores") +
      xlab("EBGM") +
      theme_bw() +
      theme(axis.title.x = element_text(face = "bold"),
            axis.title.y = element_text(face = "bold"),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            axis.line = element_line(color = "black"),
            panel.spacing = unit(0, units = "cm"),
            plot.title = element_text(face = "bold"))
  }
}
#' Summarize an openEBGM object
#'
#' @param object An openEBGM object constructed by \code{\link{ebScores}}
#' @param plot.out A logical value indicating whether or not a histogram of the
#'                 EBGM scores should be displayed
#' @param log.trans A logical value indicating whether or not the data should be
#'                  log-transformed.
#' @param ... Additional arguments affecting the summary produced
#'
#' @details This function provides a brief summary of the results of the
#'          calculations performed in the \code{\link{ebScores}} function. In
#'          particular, it provides the numerical summary of the EBGM and
#'          QUANT_* vectors.
#'
#' @details Additionally, calling \code{\link[base]{summary}} on an openEBGM
#'          object will produce a histogram of the EBGM scores. By setting the
#'          log.trans parameter to TRUE, one can convert the EBGM score to
#'          EBlog2, which is a Bayesian version of the information criterion
#'          (DuMouchel).
#'
#' @examples
#' theta_init <- data.frame(alpha1 = c(0.2, 0.1),
#'                          beta1  = c(0.1, 0.1),
#'                          alpha2 = c(2,   10),
#'                          beta2  = c(4,   10),
#'                          p      = c(1/3, 0.2)
#'                          )
#' data(caers)
#' proc <- processRaw(caers)
#' squashed <- squashData(proc, bin_size = 100, keep_pts = 100)
#' squashed <- squashData(squashed, count = 2, bin_size = 10, keep_pts = 20)
#' suppressWarnings(
#'   hypers <- autoHyper(data = squashed, theta_init = theta_init)
#' )
#' ebout <- ebScores(processed = proc, hyper_estimate = hypers)
#' summary(ebout)
#' summary(ebout, plot.out = FALSE)
#' summary(ebout, log.trans = TRUE)
#'
#' @references DuMouchel W (1999). "Bayesian Data Mining in Large Frequency
#'   Tables, With an Application to the FDA Spontaneous Reporting System."
#'   \emph{The American Statistician}, 53(3), 177-190.
#' @keywords openEBGM
#' @export
summary.openEBGM <- function(object, plot.out = TRUE, log.trans = FALSE, ...) {
  if(any(grepl("QUANT", names(object$data)))) {
    tmp <- object$data[,grep("EB|QUANT", names(object$data))]
  } else {
    #tmp <- as.data.frame(object$data$EBGM)
    tmp <- as.data.frame(object$data$EBGM, stringsAsFactors = TRUE)
    names(tmp) <- "EBGM"
  }
  if(log.trans == TRUE) {
    tmp <- log(tmp, base = 2)
  }
  ebsummary <- summary(tmp)
  if(plot.out == TRUE) {
    hist(tmp$EBGM, xlab = "EBGM", main = "Histogram of Raw EBGM Scores")
  }
  cat("\nSummary of the EB-Metrics\n")
  print(ebsummary)
}

#Hack to trick 'R CMD check'
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("cor", "head", "hist", "var1var2", "EBGM", "EBGM_fuzzy",
                           "N_col"))
}
