log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(DESeq2)
library(IHW)
library(ggplot2)
library(dplyr)
library(stringr)

if (snakemake@threads > 1) {
    library("BiocParallel")
    parallel <- TRUE
    register(MulticoreParam(snakemake@threads))
} else {
    parallel <- FALSE
}

t2g_fp <- snakemake@input[[2]]
t2g <- readr::read_table2(t2g_fp, col_names = c("tx", "ensgene", "symbol"))
t2g <- t2g %>%
    dplyr::distinct(ensgene, symbol)

dds <- readRDS(snakemake@input[[1]])

print(dds)

# Grab last variable at end of formula for contrasts
design_formula <- snakemake@params[['formula']]
s <- str_remove_all(design_formula, " |~")
s <- str_split(s, "\\+")
vars <- s[[1]]
var <- vars[length(vars)]

print(snakemake@params[['contrast']])
print("Creating results for the following variable:")
print(var)

print(resultsNames(dds))

contrast_coef <- paste(c(var, snakemake@params[['contrast']][1], "vs", snakemake@params[['contrast']][2]), collapse="_")
de_contrast <- c(var, snakemake@params[['contrast']][1], snakemake@params[['contrast']][2])

mle_res <- results(dds, contrast=de_contrast, filterFun=ihw, alpha = .05, parallel = parallel)
map_res <- lfcShrink(dds, coef=contrast_coef, type = "apeglm", parallel = parallel)

print("MLE LFC:")
print(mle_res)
print(summary(mle_res))

print("MAP LFC:")
print(map_res)
print(summary(map_res))

mle_df <- mle_res %>%
  data.frame() %>%
  tibble::rownames_to_column(var = "ensgene") %>%
  as_tibble() %>%
  left_join(t2g) %>%
  dplyr::select(symbol, ensgene, everything()) %>%
  dplyr::arrange(padj) %>%
  dplyr::distinct(symbol, .keep_all = TRUE)

map_df <- map_res %>%
    data.frame() %>%
    tibble::rownames_to_column(var = "ensgene") %>%
    as_tibble() %>%
    left_join(t2g) %>%
    dplyr::select(symbol, ensgene, everything()) %>%
    dplyr::arrange(padj) %>%
    dplyr::distinct(symbol, .keep_all = TRUE)

print("MLE dataframe")
print(mle_df)
print("MAP dataframe")
print(map_df)


mleplot <- mle_df %>%
    dplyr::mutate(
        significant = ifelse(padj < .05, "padj < .05", "padj >= .05"),
        direction = ifelse(log2FoldChange > 0, "Upregulated", "Downregulated")
    ) %>%
    ggplot(aes(log10(baseMean),log2FoldChange)) +
    geom_point(aes(color = significant, shape = direction)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    labs(x = "log10 Expression", y = "MLE Log2FoldChange") +
    theme_bw()

mapplot <- map_df %>%
    dplyr::mutate(
        significant = ifelse(padj < .05, "padj < .05", "padj >= .05"),
        direction = ifelse(log2FoldChange > 0, "Upregulated", "Downregulated")
    ) %>%
    ggplot(aes(log10(baseMean),log2FoldChange)) +
    geom_point(aes(color = significant, shape = direction)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    labs(x = "log10 Expression", y = "MAP Log2FoldChange") +
    theme_bw()

readr::write_csv(mle_df, snakemake@output[['mleres']])
readr::write_csv(map_df, snakemake@output[['mapres']])
ggsave(snakemake@output[['mlema']], plot = mleplot, width = 10, height = 7)
ggsave(snakemake@output[['mapma']], plot = mapplot, width = 10, height = 7)
