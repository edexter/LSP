# This is a collection of functions that are needed for the dotplot script to
# run. Some of these are custom and some are forked from the pafr R package and
# subsequently modified for my purposes.

################################################################################
# Function to make sure that all of the alignment files end with an empty line
# Not all tools do this, but the PAFR reader package requires it
################################################################################

appendNewlineIfNeeded <- function(filePath) {
  # Read the last byte of the file
  con <- file(filePath, "rb")
  seek(con, where = -1, origin = "end")
  lastChar <- readBin(con, "raw", 1)
  close(con)
  
  # Check if the last character is not a newline (ASCII 10 in decimal)
  if (as.integer(lastChar) != 10) {
    # Append a newline to the file
    con <- file(filePath, "ab")
    writeBin(as.raw(10), con)
    close(con)
    cat("Newline appended to file", filePath, "\n")
  } else {
    cat("File", filePath, "already ends with a newline\n")
  }
}

################################################################################
# Modified version of the plot function from pafr to allow more customization
################################################################################

customPlotFun <- function (ali, order_by = c("size", "qstart", "provided"), 
                           label_seqs = FALSE, dashes = TRUE, ordering = list(), alignment_colour = "black", 
                           xlab = x_labels, ylab = x_labels, line_size = 2,
                           contigLineSize = 1, contigLineColor = "gray") 
{
  by <- match.arg(order_by)
  if (by == "provided") {
    check_ordering(ali, ordering)
    ali <- ali[ali$qname %in% ordering[[1]] & ali$tname %in% 
                 ordering[[2]], ]
  }
  seq_maps <- order_seqs(ali, by, ordering)
  ali <- add_pos_in_concatentaed_genome(ali, seq_maps)
  p <- ggplot() + geom_segment(data = ali[ali$strand == "+", 
  ], aes_string(x = "concat_qstart", xend = "concat_qend", 
                y = "concat_tstart", yend = "concat_tend"), 
  size = line_size, colour = alignment_colour) + geom_segment(data = ali[ali$strand == 
                                                                           "-", ], aes_string(x = "concat_qend", xend = "concat_qstart", 
                                                                                              y = "concat_tstart", yend = "concat_tend"), 
                                                              size = line_size, colour = alignment_colour) + coord_equal() + 
    scale_x_continuous(xlab, labels = Mb_lab_noUnit) + scale_y_continuous(ylab, 
                                                                   labels = Mb_lab_noUnit)
  if (!is.null(xlab)) {
    p <- p + scale_x_continuous(name = xlab, labels = Mb_lab_noUnit)
  }
  if (!is.null(ylab)) {
    p <- p + scale_y_continuous(name = ylab, labels = Mb_lab_noUnit)
  }
  
  if (dashes) {
    p <- p + geom_hline(yintercept = c(seq_maps[["tmap"]], 
                                       sum(unique(ali$tlen))), linetype = 3, size = contigLineSize, color = contigLineColor) + geom_vline(xintercept = c(seq_maps[["qmap"]], 
                                                                                                         sum(unique(ali$qlen))), linetype = 3, size = contigLineSize, color = contigLineColor)
  }
  if (label_seqs) {
    qname_df <- dotplot_name_df(seq_maps[["qmap"]], 
                                seq_maps[["qsum"]])
    tname_df <- dotplot_name_df(seq_maps[["tmap"]], 
                                seq_maps[["tsum"]])
    p <- p + geom_text(data = qname_df, aes_string(label = "seq_name", 
                                                   x = "centre", y = "0"), vjust = 1, check_overlap = TRUE) + 
      geom_text(data = tname_df, aes_string(label = "seq_name", 
                                            x = "0", y = "centre"), angle = 90, 
                vjust = 0, check_overlap = TRUE)
  }
  p$seq_map_fxn <- function(bed, query = TRUE, ...) {
    map_n <- if (query) 
      1
    else 3
    seq_map <- seq_maps[[map_n]]
    check_chroms <- bed[["chrom"]] %in% names(seq_map)
    if (!(all(check_chroms))) {
      if (!(any(check_chroms))) {
        stop("None of the chromosomes represented this bed file are part of the dotplot")
      }
      else {
        bed <- bed[bed$chrom %in% names(seq_map), ]
        missing <- unique(bed[["chrom"]][!check_chroms])
        msg <- paste(length(missing), "of the chromosomes in this bed file are not part of the dotplot:\n  ", 
                     paste(missing, collapse = ", "))
        warning(msg, call. = FALSE)
      }
    }
    data.frame(istart = bed[["start"]] + seq_map[bed[["chrom"]]], 
               iend = bed[["end"]] + seq_map[bed[["chrom"]]], 
               len = seq_maps[[map_n + 1]])
  }
  p
}

################################################################################
# "Check ordering" sub-function from pafr package for plotting (unmodified)
################################################################################
check_ordering <- function (ali, ordering) {
  q_in_order <- unique(ali[["qname"]]) %in% ordering[[1]]
  t_in_order <- unique(ali[["tname"]]) %in% ordering[[2]]
  if (any(!q_in_order)) {
    msg <- paste("Dropping data from sequences absent from ordering:\n", 
                 paste(unique(ali[["qname"]])[!q_in_order], 
                       collapse = ","))
    warning(msg, call. = FALSE)
  }
  if (any(!t_in_order)) {
    msg <- paste("Dropping data from sequences absent from ordering:\n", 
                 paste(unique(ali[["tname"]])[!t_in_order], 
                       collapse = ","))
    warning(msg, call. = FALSE)
  }
  return(invisible())
}

################################################################################
# "add_pos_in_concatentaed_genome" sub-function from pafr package for 
# plotting (unmodified)
################################################################################
add_pos_in_concatentaed_genome <- function (ali, maps) {
  ali$concat_qstart <- ali$qstart + maps[["qmap"]][ali$qname]
  ali$concat_qend <- ali$qend + maps[["qmap"]][ali$qname]
  ali$concat_tstart <- ali$tstart + maps[["tmap"]][ali$tname]
  ali$concat_tend <- ali$tend + maps[["tmap"]][ali$tname]
  ali
}

################################################################################
# "order_seqs" sub-function from pafr package for plotting (unmodified)
################################################################################

order_seqs <- function (ali, by, ordering = list()) {
  chrom_lens <- chrom_sizes(ali)
  qsum <- sum(chrom_lens[["qlens"]][, 2])
  tsum <- sum(chrom_lens[["tlens"]][, 2])
  if (by == "size") {
    q_idx <- order(chrom_lens[["qlens"]][, 2], decreasing = TRUE)
    t_idx <- order(chrom_lens[["tlens"]][, 2], decreasing = TRUE)
    qmap <- structure(.Names = chrom_lens[["qlens"]][q_idx, 
                                                     1], c(0, head(cumsum(chrom_lens[["qlens"]][q_idx, 
                                                                                                2]), -1)))
    tmap <- structure(.Names = chrom_lens[["tlens"]][t_idx, 
                                                     1], c(0, head(cumsum(chrom_lens[["tlens"]][t_idx, 
                                                                                                2]), -1)))
  }
  else if (by == "qstart") {
    q_idx <- order(chrom_lens[["qlens"]][, 2], decreasing = TRUE)
    qmap <- structure(.Names = chrom_lens[["qlens"]][q_idx, 
                                                     1], c(0, head(cumsum(chrom_lens[["qlens"]][q_idx, 
                                                                                                2]), -1)))
    longest_by_target <- slice_max(group_by(ali, .data[["tname"]]), 
                                   .data[["alen"]])
    t_idx <- order(qmap[longest_by_target$qname] + longest_by_target$qstart)
    tmap <- sort(structure(.Names = longest_by_target$tname[t_idx], 
                           c(0, head(cumsum(longest_by_target$tlen[t_idx]), 
                                     -1))))
  }
  else if (by == "provided") {
    if (length(ordering) != 2) {
      stop("To use 'provided' sequence ordering argument 'ordering' must by a list with two character vectors")
    }
    qord <- ordering[[1]]
    tord <- ordering[[2]]
    q_idx <- match(qord, chrom_lens[["qlens"]][["qname"]])
    if (any(is.na(q_idx))) {
      msg <- paste("Sequence(s) provided for ordering are not present in alignment:\n", 
                   qord[is.na(q_idx)], "\n")
      stop(msg)
    }
    t_idx <- match(tord, chrom_lens[["tlens"]][["tname"]])
    if (any(is.na(t_idx))) {
      msg <- paste("Sequence(s) provided for ordering are not present in alignment:\n", 
                   tord[is.na(t_idx)], "\n")
      stop(msg)
    }
    qmap <- structure(.Names = chrom_lens[["qlens"]][q_idx, 
                                                     1], c(0, head(cumsum(chrom_lens[["qlens"]][q_idx, 
                                                                                                2]), -1)))
    tmap <- structure(.Names = chrom_lens[["tlens"]][t_idx, 
                                                     1], c(0, head(cumsum(chrom_lens[["tlens"]][t_idx, 
                                                                                                2]), -1)))
  }
  list(qmap = qmap, qsum = qsum, tmap = tmap, tsum = tsum)
}

################################################################################
# "Mb_lab" sub-function from pafr package for plotting (modified to remove units)

Mb_lab_noUnit <- function (x) {
  paste(x/1e+06)
}
