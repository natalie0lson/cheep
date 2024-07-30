

plot.differentialTest <- function(x, level = NULL, data_only = FALSE, ...) {
  signif_taxa <- x$significant_taxa
  if ("phyloseq" %in% class(x$data)) {
    if (!is.null(x$data@tax_table)) {
      signif_taxa <- otu_to_taxonomy(signif_taxa, x$data, level = level)
      if (length(unique(signif_taxa)) != length(unique(x$significant_taxa))) {
        # Make sure if repeated taxa add unique otu identifiers
        signif_taxa <- paste0(signif_taxa, " (", x$significant_taxa, ")")
      }
    }
  }
  if (length(x$significant_models) != 0) {
    var_per_mod <- length(x$restrictions_DA) + length(x$restrictions_DV)
    total_var_count <- length(signif_taxa) * var_per_mod
    df <- as.data.frame(matrix(NA, nrow = total_var_count, ncol = 6))  # Increased column count to 6
    colnames(df) <- c("x", "xmin", "xmax", "taxa", "variable", "p_fdr")  # Added "p_fdr" column
    qval <- stats::qnorm(.975)
    restricts_mu <- attr(x$restrictions_DA, "index")
    restricts_phi <- attr(x$restrictions_DV, "index")
    
    count <- 1
    for (i in 1:length(x$significant_models)) {
      
      # Below from print_summary_bbdml, just to get coefficient names
      tmp <- x$significant_models[[i]]
      coefs.mu <- tmp$coefficients[1:tmp$np.mu,, drop = FALSE]
      rownames(coefs.mu) <- paste0(substring(rownames(coefs.mu), 4), "\nDifferential Abundance")
      coefs.mu <- coefs.mu[restricts_mu,, drop = FALSE]
      
      coefs.phi <- tmp$coefficients[(tmp$np.mu + 1):nrow(tmp$coefficients),, drop = FALSE]
      rownames(coefs.phi) <- paste0(substring(rownames(coefs.phi), 5), "\nDifferential Variability")
      coefs.phi <- coefs.phi[restricts_phi - tmp$np.mu,, drop = FALSE]
      
      coefs <- rbind(coefs.mu, coefs.phi)
      for (j in 1:var_per_mod) {
        df[count, 1:3] <- c(coefs[j, 1], coefs[j, 1] - qval * coefs[j, 2],
                            coefs[j, 1] + qval * coefs[j, 2])
        df[count, 4:6] <- c(signif_taxa[i], rownames(coefs)[j], x$p_fdr[i])
        count <- count + 1
      }
    }
    
    if (!(data_only)) {
      taxa <- xmin <- xmax <- NULL
      
      ggplot2::ggplot(df, ggplot2::aes(x = x, y = taxa)) +
        ggplot2::geom_vline(xintercept = 0, color = "gray50", lty = "dashed",
                            alpha = 0.75, lwd = 1) +
        ggplot2::geom_point() +
        ggplot2::geom_errorbarh(ggplot2::aes(xmin = xmin, xmax = xmax, colour = xmin <= 0 & xmax <= 0 | xmin >= 0 & xmax >= 0), height = .3) +
        ggplot2::geom_text(ggplot2::aes(label = sprintf("p-adj = %.4f", as.numeric(p_fdr))), size = 14, hjust = 0.5, vjust = 1.5) +  # Adjusted format specifier
        ggplot2::theme_bw() +
    
        ggplot2::facet_wrap(~variable, scales = "free_x", nrow = 1, labeller = ggplot2::labeller(variable = label_value)) +  # Specify facet labels
        ggplot2::labs(title = "", x = "", y = "") +
        ggplot2::scale_y_discrete(limits = rev(df$taxa)) +
        ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
        ggplot2::theme(
          legend.position = "none",
          axis.text.y = ggplot2::element_text(size = 50, margin = ggplot2::margin(r = 5)),
          axis.text.x = ggplot2::element_text(hjust = 1, size = 50),
          strip.text = ggplot2::element_text(size = 50),
          axis.title = ggplot2::element_text(size =50),
          plot.title = ggplot2::element_text(size = 50),
          strip.background = ggplot2::element_blank()
        )
    } else {
      return(df)
    }
    
  } else {
    message("No taxa were found to be significantly different using your model specification. \nPlease verify that your formulas are correctly specified.")
  }
}