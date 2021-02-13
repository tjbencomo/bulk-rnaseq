log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(DESeq2)
library(IHW)
library(ggplot2)
library(dplyr)
library(stringr)
library(annotables)

if (snakemake@threads > 1) {
    library("BiocParallel")
    parallel <- TRUE
    register(MulticoreParam(snakemake@threads))
} else {
    parallel <- FALSE
}

dds <- readRDS(snakemake@input[[1]])

print(dds)

design_formula <- snakemake@params[['formula']]
s <- str_split(design_formula, '\\+', simplify=T)
var <- str_trim(s[1, dim(s)[2]])

print(snakemake@params[['contrast']])
# eb_coef <- str_c(c(var, snakemake@params[['contrast']][1], "vs", snakemake@params[['contrast']][2]), collapse="_")
# print(eb_coef)
mle_res <- results(dds, contrast=c(var, snakemake@params[['contrast']]),
                   filterFun=ihw, alpha = .05, parallel = parallel)
map_res <- lfcShrink(dds, contrast=c(var, snakemake@params[['contrast']]), type = "ashr",
                     parallel = parallel)
# map_res <- lfcShrink(dds, coef=eb_coef, type = "apeglm")

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
  dplyr::mutate(ensgene = str_split(ensgene, '\\.', simplify=T)[, 1]) %>%
  left_join(grch38) %>%
  dplyr::select(-chr, -start, -end, -strand, -biotype, -description) %>%
  dplyr::select(symbol, ensgene, everything()) %>%
  dplyr::arrange(padj) %>%
  dplyr::distinct(symbol, .keep_all = TRUE)

map_df <- map_res %>%
    data.frame() %>%
    tibble::rownames_to_column(var = "ensgene") %>%
    as_tibble() %>%
    dplyr::mutate(ensgene = str_split(ensgene, '\\.', simplify=T)[, 1]) %>%
    left_join(grch38) %>%
    dplyr::select(-chr, -start, -end, -strand, -biotype, -description) %>%
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
