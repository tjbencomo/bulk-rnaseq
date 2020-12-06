library(tximport)
library(DESea2)
library(readr)
library(stringr)

res_dirs <- snakemake@input[['cts']]
tx2g_fp <- snakemake@input[['tx2g']]
samples_fp <- snakemake@input[['samples']]

design_formula <- snakemake@params[['formula']]

if (snakemake@threads > 1) {
    library("BiocParallel")
    parallel <- TRUE
    register(MulticoreParam(snakemake@threads))
} else {
    parallel <- FALSE
}

print(res_dirs)
print(tx2g_fp)

files <- file.path(res_dirs, "abundance.h5")
names(files) <- basename(dirname(files))
tx2g <- read_table2(tx2g_fp, col_names = c("tx", "ensgene", "symbol"))
txi <- tximport(files, type = "kallisto", txOut = FALSE, tx2g = tx2g[, 1:2])
print(names(files))

samples <- read_csv(samples_fp)
samples$id <- paste(samples$patient, "-", samples$sample, sep = "")

## Ensure factor ordering based on config specifications
vars <- snakemake@params[['levels']]
var_levels <- str_split(s, ';')
for (var in var_levels) {
    s <- str_split(var, '=|,', simplify=T)
    col <- s[1, 1]
    level_order = s[1, 2:dim(s)[2]]
    samples[, col] <- factor(samples[, col], levels=level_order)
}
print(samples[match(names(files), samples$id),])

f <- as.formula(design_formula)
ddsTxi <- DESeqDataSetFromTximport(txi, colData = samples, design = f)

keep <- rowSums(counts(ddsTxi)) >= 1
ddsTxi <- ddsTxi[keep, ]

dds <- DESeq(ddsTxi, parallel=parallel)

vst_cts <- vst(dds, blind=FALSE)

saveRDS(dds, file=snakemake@output[['deseq']])
saveRDS(vst_cts, file=snakemake@output[['cts']])
