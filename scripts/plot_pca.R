log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(DESeq2)
library(ggplot2)
library(stringr)

vsd <- readRDS(snakemake@input[[1]])

print(vsd)

label_vars <- str_split(snakemake@params[['label_vars']], ',', simplify=T)[1, ]

pcaplot <- plotPCA(vsd, intgroup=label_vars)

ggsave(snakemake@output[[1]], pcaplot)
