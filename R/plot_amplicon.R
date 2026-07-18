#' Plot the amplicon map with primers and restriction site
#'
#' Generates a graphical representation of the amplicon showing the
#' position of the forward/reverse primers and the cut site of the
#' selected restriction enzyme.
#'
#' @param pcr_result `"rflp_amplicon"` object returned by [simulate_pcr()].
#' @param fwd_primer,rev_primer Primer sequences (character strings), the
#'   same ones used in [simulate_pcr()].
#' @param restriction_result List returned by [find_restriction_site()].
#' @param enzyme_name Enzyme name, for labeling purposes only (default
#'   `"Enzyme"`).
#' @param site_index Integer. If [find_restriction_site()] found more than
#'   one site, indicates which one to use for the plot (default `1`, the
#'   first one).
#'
#' @return A `ggplot` object. Can be saved with `ggplot2::ggsave()`.
#' @export
plot_amplicon_map <- function(pcr_result, fwd_primer, rev_primer,
                               restriction_result, enzyme_name = "Enzyme",
                               site_index = 1) {
  if (restriction_result$n_sites == 0) {
    stop("There are no restriction sites to plot. Check the result of find_restriction_site().")
  }
  if (site_index > restriction_result$n_sites) {
    stop(sprintf(
      "site_index = %d is out of range; only %d site(s) were found.",
      site_index, restriction_result$n_sites
    ))
  }

  size_amp <- pcr_result$size
  start_amp <- pcr_result$start

  rel_fwd_start <- 0
  rel_fwd_end <- nchar(fwd_primer)
  rel_rev_start <- size_amp
  rel_rev_end <- size_amp - nchar(rev_primer)
  rel_cut <- restriction_result$sites$cut_pos[site_index] - start_amp + 1

  fragment_1 <- rel_cut
  fragment_2 <- size_amp - rel_cut

  ggplot2::ggplot() +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = size_amp, y = 2, yend = 2),
                           color = "grey50", linewidth = 2) +
    ggplot2::geom_segment(ggplot2::aes(x = rel_fwd_start, xend = rel_fwd_end, y = 3, yend = 3),
                           arrow = ggplot2::arrow(length = ggplot2::unit(0.3, "cm")),
                           color = "#1f77b4", linewidth = 1.5) +
    ggplot2::geom_segment(ggplot2::aes(x = rel_rev_start, xend = rel_rev_end, y = 3, yend = 3),
                           arrow = ggplot2::arrow(length = ggplot2::unit(0.3, "cm")),
                           color = "#d62728", linewidth = 1.5) +
    ggplot2::geom_vline(xintercept = rel_cut, color = "#2ca02c",
                         linetype = "dashed", linewidth = 1.2) +
    ggplot2::annotate("text", x = rel_fwd_end / 2, y = 3.3, label = "Fwd primer",
                       color = "#1f77b4", fontface = "bold") +
    ggplot2::annotate("text", x = size_amp - (size_amp - rel_rev_end) / 2, y = 3.3,
                       label = "Rev primer", color = "#d62728", fontface = "bold") +
    ggplot2::annotate("text", x = rel_cut, y = 0.8,
                       label = paste0(enzyme_name, " site\n(bp ", rel_cut, ")"),
                       color = "#2ca02c", fontface = "bold") +
    ggplot2::annotate("text", x = size_amp / 2, y = 1.7,
                       label = paste("Total amplicon:", size_amp, "bp"), color = "grey30") +
    ggplot2::labs(
      title = "In silico PCR map: primers and restriction site",
      subtitle = paste0("Fragment 1: ", fragment_1, " bp | Fragment 2: ", fragment_2, " bp"),
      x = "Relative position in the amplicon (bp)", y = ""
    ) +
    ggplot2::xlim(-10, size_amp + 10) +
    ggplot2::ylim(0.5, 4) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      plot.subtitle = ggplot2::element_text(color = "grey40")
    )
}

#' Plot the amplicon sequence highlighting the restriction site
#'
#' Generates a letter-by-letter view of the amplicon, arranged in rows,
#' highlighting the bases that are part of the selected restriction
#' enzyme's cut site.
#'
#' @param pcr_result `"rflp_amplicon"` object returned by [simulate_pcr()].
#' @param restriction_result List returned by [find_restriction_site()].
#' @param enzyme_name Enzyme name, for labeling purposes only.
#' @param bases_per_row Number of nucleotides shown per row (default `50`).
#'
#' @return A `ggplot` object.
#' @export
plot_sequence_map <- function(pcr_result, restriction_result,
                               enzyme_name = "Enzyme", bases_per_row = 50) {
  amp_seq <- pcr_result$amplicon
  start_amp <- pcr_result$start
  end_amp <- pcr_result$end

  bases_vector <- strsplit(as.character(amp_seq), "")[[1]]
  genomic_pos <- start_amp:end_amp
  amplicon_pos <- seq_along(bases_vector)

  site_type <- paste0(enzyme_name, " site")

  seq_df <- data.frame(
    Base = bases_vector, GenPos = genomic_pos, AmpPos = amplicon_pos,
    Type = "Normal", stringsAsFactors = FALSE
  )

  if (restriction_result$n_sites > 0) {
    for (i in seq_len(nrow(restriction_result$sites))) {
      s <- restriction_result$sites$start[i]
      e <- restriction_result$sites$end[i]
      seq_df$Type[s:e] <- site_type
    }
  }

  seq_df$Row <- (seq_df$AmpPos - 1) %/% bases_per_row + 1
  seq_df$Column <- (seq_df$AmpPos - 1) %% bases_per_row + 1

  row_labels <- stats::aggregate(GenPos ~ Row, data = seq_df, FUN = min)
  row_labels$Label <- paste0("Base: ", row_labels$GenPos, " ->")

  fill_colors <- stats::setNames(c("#f7f9fa", "#ffcca3"), c("Normal", site_type))
  text_colors <- stats::setNames(c("#2c3e50", "#e65c00"), c("Normal", site_type))

  ggplot2::ggplot(seq_df, ggplot2::aes(x = Column, y = -Row)) +
    ggplot2::geom_tile(ggplot2::aes(fill = Type), color = "white",
                        linewidth = 0.5, width = 0.9, height = 0.9) +
    ggplot2::geom_text(ggplot2::aes(label = Base, color = Type),
                        fontface = "bold", size = 4.5) +
    ggplot2::geom_text(
      data = row_labels,
      ggplot2::aes(x = -2, y = -Row, label = Label),
      hjust = 1, color = "grey40", size = 3.5, fontface = "italic",
      inherit.aes = FALSE
    ) +
    ggplot2::scale_fill_manual(values = fill_colors) +
    ggplot2::scale_color_manual(values = text_colors) +
    ggplot2::labs(
      title = "Full amplicon sequence",
      subtitle = paste0(enzyme_name, " site highlighted (total length: ", length(amp_seq), " bp)"),
      x = "Relative position within the row (bp)", y = ""
    ) +
    ggplot2::xlim(-8, bases_per_row + 1) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 14, color = "#2c3e50"),
      plot.subtitle = ggplot2::element_text(size = 10, color = "grey40"),
      legend.position = "bottom",
      legend.title = ggplot2::element_blank()
    )
}
