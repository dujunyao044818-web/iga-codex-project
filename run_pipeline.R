dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("output/models", recursive = TRUE, showWarnings = FALSE)
dir.create("output/logs", recursive = TRUE, showWarnings = FALSE)

writeLines("IgA pipeline started", "output/logs/pipeline.log")

suppressPackageStartupMessages({
  library(GEOquery)
  library(limma)
  library(ggplot2)
})

gset <- GEOquery::getGEO("GSE104948", GSEMatrix = TRUE)
g <- gset[[1]]

expr <- Biobase::exprs(g)
pdata <- Biobase::pData(g)

if (max(expr, na.rm = TRUE) > 50) {
  expr <- log2(expr + 1)
}

expr <- limma::normalizeBetweenArrays(expr)

write.csv(expr, "output/tables/GSE104948_expression_matrix.csv")
write.csv(pdata, "output/tables/GSE104948_sample_metadata.csv")

pca <- prcomp(t(expr), scale. = TRUE)
pca_df <- data.frame(
  sample = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
)

p <- ggplot(pca_df, aes(PC1, PC2)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(title = "GSE104948 PCA")

ggsave("output/figures/Figure1_PCA.pdf", p, width = 7, height = 5)

writeLines(capture.output(sessionInfo()), "output/logs/sessionInfo.txt")
writeLines("IgA pipeline completed", "output/logs/pipeline.log")
