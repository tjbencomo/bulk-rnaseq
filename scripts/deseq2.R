log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(tximeta)
library(DESeq2)
library(readr)
library(stringr)
library(org.Hs.eg.db)

res_dirs <- snakemake@input[['cts']]
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

samples$names <- samples$id
samples$files <- files
print(samples)

se <- tximeta(samples)
gse <- summarizeToGene(se)
gse <- addIds(gse, "SYMBOL", gene = T)

f <- as.formula(design_formula)
dds <- DESeqDataSet(gse, design = f)
print("Old ordering for condition")
print(colData(dds)$condition)

## Ensure factor ordering based on config specifications
vars <- snakemake@params[['levels']]
var_levels <- str_split(vars, ';', simplify=T)
for (var in var_levels) {
    print(paste("Variable:", var))
    s <- str_split(var, '=|,', simplify=T)
    col <- s[1, 1]
    level_order = s[1, 2:dim(s)[2]]
    colData(dds)[, col] <- factor(colData(dds)[, col], level_order)
    print(paste("Ordering for", col))
    print(levels(colData(dds)[, col]))
}

dds <- DESeq(dds, parallel=parallel)
print(dds)

vst_cts <- vst(dds, blind=FALSE)
print(vst_cts)

saveRDS(dds, file=snakemake@output[['deseq']])
saveRDS(vst_cts, file=snakemake@output[['cts']])
