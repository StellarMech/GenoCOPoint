#' Detect CO Points
#'
#' This function reads a genotype file, processes the data, detects CO points and generates a PDF visualization and a result file..
#'
#' @param geno_file The path to the genotype file with three columns of information 
#' (name of the chromosome where the SNP is located,
#' absolute position of the SNP,
#' code of the genotype corresponding to the SNP 0/1,
#' read counts of ref genotype,
#' read counts of alt genotype) .
#' @param truth_file The path to the truth file with CO intervals.
#' The file should contain the chromosome name, start and end positions of the CO intervals.
#' @param output_pdf Path to the output PDF file for visualization.
#' @param output_txt Path to the output text file for results.
#' @param window_size Size of the sliding window for smoothing.
#' @param smooth_threshold Threshold for smoothing.
#' @param min_support Minimum number of consecutive SNPs for noise filtering.
#' @param consistency Consistency ratio for filtering.
#' @return A list where each element contains the CO positions for a chromosome.
#' @examples
#' # Example usage:
#' # Create a temporary input file
#' geno_file <- "GenoCOPoint/examples/genofiles/AAACGCTCATATAGCC_smt_genotypes_withname_reads.txt"
#' truth_file <- "GenoCOPoint/examples/truthfiles/AAACGCTCATATAGCC_allele_cnts_at_markers_sorted_co_pred.txt"
#' # Run the function
#' result <- GenoCOPoint(geno_file,truth_file, output_pdf = "example_output.pdf", output_txt = "example_output.txt", window_size = 5, smooth_threshold = 0.6, min_support = 7, consistency = 0.85)
#' print(result)
#'
#' # Clean up
#' unlink(geno_file)
#' unlink(truth_file)
#' unlink("example_output.pdf")
#' unlink("example_output.txt")
#' @export
# 主函数：检测变点并生成结果
GenoCOPoint <- function(
  geno_file,
  truth_file,
  output_pdf = "CO_points.pdf",
  output_txt = "COresult.txt",
  window_size = 5,
  smooth_threshold = 0.6,
  min_support = 7,
  consistency = 0.85
) {
  library(changepoint)
  # 读取五列数据
  all_data <- read.table(geno_file,
                         header = FALSE,
                         col.names = c("chr", "position", "genotype", "ref_count", "alt_count"))
  # 读取真实CO区间
  truth_data <- read.table(truth_file, header = FALSE, col.names = c("chr", "start", "end", "x1", "x2", "x3"))
  # 按染色体分组
  chr_groups <- split(all_data, all_data$chr)

  # 平滑函数
  smooth_genotypes <- function(x, window_size = 5, smooth_threshold = 0.6) {
    n <- length(x)
    smoothed <- numeric(n)
    half_window <- floor(window_size / 2)
    for (i in 1:n) {
      start <- max(1, i - half_window)
      end <- min(n, i + half_window)
      window <- x[start:end]
      count0 <- sum(window == 0, na.rm = TRUE)
      count1 <- sum(window == 1, na.rm = TRUE)
      if (count0 > count1 && count0 / (count0 + count1) >= smooth_threshold) {
        smoothed[i] <- 0
      } else if (count1 > count0 && count1 / (count0 + count1) >= smooth_threshold) {
        smoothed[i] <- 1
      } else {
        smoothed[i] <- x[i]
      }
    }
    return(smoothed)
  }
  # 秩变换
  compute_centered_rank <- function(x) {
    rank_x <- rank(x, ties.method = "average", na.last = "keep")
    centered_rank <- rank_x - (length(na.omit(x)) + 1) / 2
    return(centered_rank)
  }
  # 颜色映射
  get_segment_color <- function(genotype) {
    ifelse(genotype == 0, "red", "blue")
  }

  pdf(output_pdf, width = 12, height = 9)
  all_change_points <- numeric()
  COresult <- data.frame(Chromosome = character(), Position = numeric(), stringsAsFactors = FALSE)

  for (chr_name in names(chr_groups)) {
    chr_data <- chr_groups[[chr_name]]
    chr_truth <- truth_data[truth_data$chr == chr_name, ]
    raw_positions <- chr_data$position
    raw_genotypes <- chr_data$genotype
    ref_counts <- chr_data$ref_count
    alt_counts <- chr_data$alt_count
    encoded_genotypes <- ifelse(raw_genotypes == 0.5, NA, raw_genotypes)
    encoded_genotypes <- smooth_genotypes(encoded_genotypes, window_size = window_size, smooth_threshold = smooth_threshold)
    valid_indices <- which(!is.na(encoded_genotypes))
    filtered_genotypes <- encoded_genotypes[valid_indices]
    filtered_positions <- raw_positions[valid_indices]
    if (length(filtered_genotypes) == 0) {
      message("染色体 ", chr_name, " 无有效数据，跳过分析")
      next
    }
    rank_data <- compute_centered_rank(filtered_genotypes)
    x_limits <- range(filtered_positions)
    tryCatch({
      result <- cpt.meanvar(rank_data, method = "PELT", penalty = "MBIC", pen.value = "2*log(n)")
      change_points <- cpts(result)
      if (length(change_points) > 0) {
        abs_points <- filtered_positions[change_points]
        sorted_order <- order(abs_points)
        sorted_cp <- change_points[sorted_order]
        sorted_abs <- abs_points[sorted_order]
        valid_cp <- logical(length(sorted_cp))
        for (i in seq_along(sorted_cp)) {
          cp_idx <- sorted_cp[i]
          left_start <- ifelse(i == 1, 1, sorted_cp[i - 1] + 1)
          left_end <- cp_idx
          right_start <- cp_idx + 1
          right_end <- ifelse(i == length(sorted_cp), length(filtered_genotypes), sorted_cp[i + 1] - 1)
          left_segment <- filtered_genotypes[left_start:left_end]
          right_segment <- filtered_genotypes[right_start:right_end]
          left_valid <- if (length(left_segment) == 0) {
            FALSE
          } else {
            left_type <- left_segment[1]
            sum(left_segment == left_type, na.rm = TRUE)/length(left_segment) >= consistency &&
              length(left_segment) >= min_support
          }
          right_valid <- if (length(right_segment) == 0) {
            FALSE
          } else {
            right_type <- right_segment[1]
            sum(right_segment == right_type, na.rm = TRUE)/length(right_segment) >= consistency &&
              length(right_segment) >= min_support
          }
          valid_cp[i] <- left_valid && right_valid
        }
        final_abs <- sorted_abs[valid_cp]
        if (length(final_abs) > 0) {
          cat("染色体", chr_name, "的CO位置在：",
              paste(format(final_abs, scientific = FALSE, big.mark = ","), collapse = ", "), "\n")
          COresult <- rbind(COresult, data.frame(Chromosome = chr_name, Position = final_abs))
          all_change_points <- c(all_change_points, final_abs)
          par(mfrow = c(2, 1), mar = c(1, 4, 1, 1), oma = c(3, 1, 2, 1), family = "serif", cex.lab = 0.88, cex.axis = 0.88, cex.main = 0.88, cex.sub = 0.88)
          # 第一图：Reads柱状图
          plot(NULL, xlim = x_limits, ylim = c(-max(alt_counts), max(ref_counts)),
               main = "", xlab = "", ylab = "Reads Count", xaxt = "n", bty = "n")
          segments(x0 = raw_positions, y0 = -alt_counts, y1 = ref_counts,
                   col = ifelse(filtered_genotypes == 0, "red", "blue"), lwd = 1.5)
          if (nrow(chr_truth) > 0) {
            for (i in 1:nrow(chr_truth)) {
              rect(chr_truth$start[i], -max(alt_counts), chr_truth$end[i], max(ref_counts),
                   col = rgb(1, 0.85, 0.2, 0.3), border = rgb(1, 0.6, 0, 0.8), lwd = 2)
            }
          }
          abline(v = final_abs, col = "black", lty = 2)
          # 第二图：基因型图（带颜色条背景）
          plot(filtered_positions, filtered_genotypes,
               type = "n", main = "", xaxt = "n",
               xlab = "", ylab = "Genotypes(0=Ref, 1=Alt)", cex.lab = 1.5,
               ylim = c(-0.1, 1.1), font.lab = 2, bty = "n")
          # 画颜色条背景
          final_abs_sorted <- sort(final_abs)
          breaks <- c(min(filtered_positions), final_abs_sorted, max(filtered_positions))
          cp_genotypes <- filtered_genotypes[match(final_abs_sorted, filtered_positions)]
          colors <- character(length(breaks) - 1)
          if (length(cp_genotypes) >= 1) {
            colors[1] <- get_segment_color(cp_genotypes[1])
          }
          if (length(cp_genotypes) >= 2) {
            for (i in 2:length(cp_genotypes)) {
              colors[i] <- ifelse(colors[i - 1] == "red", "blue", "red")
            }
          }
          if (length(colors) > length(cp_genotypes)) {
            colors[length(colors)] <- ifelse(colors[length(cp_genotypes)] == "red", "blue", "red")
          }
          for (i in 1:(length(breaks) - 1)) {
            rect(breaks[i], -0.1, breaks[i + 1], 1.1,
                 col = adjustcolor(colors[i], alpha.f = 0.15), border = NA)
          }
          points(filtered_positions, filtered_genotypes,
                 pch = 20, cex = 0.6,
                 col = ifelse(filtered_genotypes == 0, "red", "blue"))
          points(x = final_abs,
                 y = filtered_genotypes[match(final_abs, filtered_positions)],
                 pch = 19, cex = 1.2, col = "orange")
          x_range <- diff(range(filtered_positions))
          text(x = final_abs + 0.02 * x_range,
               y = rep(1.05, length(final_abs)),
               labels = format(final_abs, scientific = FALSE, big.mark = ","),
               col = "darkred", cex = 0.7, srt = 45, pos = 4, xpd = TRUE)
          abline(v = final_abs, col = "black", lty = 2)
          mtext(paste(chr_name, "Analysis"), side = 3, outer = TRUE, line = 0, cex = 1.5)
        } else {
          cat("染色体", chr_name, "所有CO位置均被过滤，不存在CO位点\n")
          par(mfrow = c(1, 1), mar = c(4, 4, 4, 2), family = "serif", cex.lab = 0.88, cex.axis = 0.88, cex.main = 0.88, cex.sub = 0.88)
          plot(NULL, xlim = x_limits, ylim = c(-max(alt_counts), max(ref_counts)),
               main = paste(chr_name, "Analysis"), xlab = "Position", ylab = "Reads Count", cex.lab = 1.5, font.lab = 2,  bty = "n")
          segments(x0 = raw_positions, y0 = -alt_counts, y1 = ref_counts,
                   col = ifelse(filtered_genotypes == 0, "red", "blue"), lwd = 1.5)
          if (nrow(chr_truth) > 0) {
            for (i in 1:nrow(chr_truth)) {
              rect(chr_truth$start[i], -max(alt_counts), chr_truth$end[i], max(ref_counts),
                   col = rgb(1, 0.85, 0.2, 0.3), border = rgb(1, 0.6, 0, 0.8), lwd = 2)
            }
          }
          mtext(paste(chr_name, "No valid CO detected"), side = 3, line = 1, cex = 1.2)
        }
      } else {
        cat("染色体", chr_name, "未检测到CO位置\n")
        par(mfrow = c(1, 1), mar = c(4, 4, 4, 2), family = "serif", cex.lab = 0.88, cex.axis = 0.88, cex.main = 0.88, cex.sub = 0.88)
        plot(NULL, xlim = x_limits, ylim = c(-max(alt_counts), max(ref_counts)),
             main = paste(chr_name, "Analysis"), xlab = "Position", ylab = "Reads Count", cex.lab = 1.5, font.lab = 2, bty = "n")
        segments(x0 = raw_positions, y0 = -alt_counts, y1 = ref_counts,
                 col = ifelse(filtered_genotypes == 0, "red", "blue"), lwd = 1.5)
        if (nrow(chr_truth) > 0) {
          for (i in 1:nrow(chr_truth)) {
            rect(chr_truth$start[i], -max(alt_counts), chr_truth$end[i], max(ref_counts),
                 col = rgb(1, 0.85, 0.2, 0.3), border = rgb(1, 0.6, 0, 0.8), lwd = 2)
          }
        }
        mtext(paste(chr_name, "No CO detected"), side = 3, line = 1, cex = 1.2)
      }
    }, error = function(e) {
      message("染色体 ", chr_name, " 分析出错：", e$message)
    })
  }
  dev.off()
  write.table(COresult, file = output_txt, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("该细胞的CO位置数量共计：", length(all_change_points), "个\n")
  return(COresult)
}