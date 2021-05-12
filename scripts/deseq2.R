log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(tximport)
library(DESeq2)
library(readr)
library(stringr)

res_dirs <- snakemake@input[['cts']]
tx2g_fp <- snakemake@input[['tx2g']]
samples_fp <- snakemake@input[['samples']]
quant_program <- snakemake@params[['aligner']]

design_formula <- snakemake@params[['formula']]

if (snakemake@threads > 1) {
    library("BiocParallel")
    parallel <- TRUE
    register(MulticoreParam(snakemake@threads))
} else {
    parallel <- FALSE
}

if (quant_program == 'kallisto') {
    files <- file.path(res_dirs, "abundance.h5")
} else {
    files <- file.path(res_dirs, "quant.sf")
}
names(files) <- basename(dirname(files))
tx2g <- read_table2(tx2g_fp, col_names = c("tx", "ensgene", "symbol"))
txi <- tximport(files, type = quant_program, txOut = FALSE, tx2g = tx2g[, 1:2])
# txi <- tximport(files, type = "kallisto", txOut = FALSE, tx2g = tx2g[, 1:2])

samples <- read.csv(samples_fp)
print("First")
print(samples)
samples$id <- paste(samples$patient, "-", samples$condition, sep = "")
print("Second")
print(samples)
# Reorder rows so they match files order
samples <- samples[match(names(files), samples$id),]
print("Third")
print(files)
print(samples)

## Ensure factor ordering based on config specifications
vars <- snakemake@params[['levels']]
var_levels <- str_split(vars, ';', simplify=T)
for (var in var_levels) {
    print(paste("Variable:", var))
    s <- str_split(var, '=|,', simplify=T)
    col <- s[1, 1]
    level_order = s[1, 2:dim(s)[2]]
    samples[, col] <- factor(samples[, col], level_order)
}

f <- as.formula(design_formula)
ddsTxi <- DESeqDataSetFromTximport(txi, colData = samples, design = f)

keep <- rowSums(counts(ddsTxi)) >= 1
ddsTxi <- ddsTxi[keep, ]

dds <- DESeq(ddsTxi, parallel=parallel)
print(dds)

vst_cts <- vst(dds, blind=FALSE)
print(vst_cts)

saveRDS(dds, file=snakemake@output[['deseq']])
saveRDS(vst_cts, file=snakemake@output[['cts']])
