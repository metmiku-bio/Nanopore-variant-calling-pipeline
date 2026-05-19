args <- commandArgs(trailingOnly=TRUE)

library(tidyverse)
library(ggplot2)
library(readr)

all_data <- lapply(args, function(f) {
  df <- read.table(f, header=FALSE)

  colnames(df) <- c("chr", "start", "end","gene", "mean_cov")

  df$sample <- tools::file_path_sans_ext(basename(f))
  return(df)
}) %>% bind_rows()

# readr::write_tsv(all_data,"all_merged_coverage.tsv")
# Summary plot: coverage per sample
p <- ggplot(all_data, aes(x = gene, y = mean_cov,fill=gene)) +

  # boxplot with fill color
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +

  # jittered points inside box
  geom_jitter(
    width = 0.2,
    size = 1.5,
    alpha = 0.6,
    color = "black"
  ) +

    scale_y_sqrt(breaks = c(0, 50, 100,500, 1000, 2000, 4000)) +
    geom_hline(yintercept = 50, linetype = "dashed")+
    theme(
      axis.text.x = element_text(angle =45 , hjust = 1),
      panel.border = element_rect(color = "black", fill = NA, size = 1)
    )

ggsave("coverage_plot.pdf", p, width = 8, height = 5)

nc_data <- all_data %>%
  filter(grepl("NC", sample))
p_nc <- ggplot(nc_data, aes(x = sample, y = mean_cov)) +
  geom_boxplot(fill = "tomato", alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.7) +
    scale_y_sqrt(breaks = c(0, 10,20,30,40,50, 100,500)) +
    geom_hline(yintercept = 30, linetype = "dashed")+
    theme(
      axis.text.x = element_text(angle =45 , hjust = 1),
      panel.border = element_rect(color = "black", fill = NA, size = 1)
    )
ggsave("NC_coverage_plot.pdf", p_nc, width = 8, height = 5)

