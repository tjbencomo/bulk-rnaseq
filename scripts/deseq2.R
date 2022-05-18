log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# Handle readr tzdata bug when using singularity mambaforge container
Sys.setenv("TZDIR"=paste0(Sys.getenv("CONDA_PREFIX"), "/share/zoneinfo"))

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

if (quant_program == 'salmon') {
    files <- file.path(res_dirs, "quant.sf")
} else{
    print("quant_program wasn't Salmon!")
    stop()
}

names(files) <- basename(dirname(files))

samples <- read.csv(samples_fp)
samples$id <- paste(samples$patient, "-", samples$condition, sep = "")
print("First")
print(samples)
# Reorder rows so they match files order
samples <- samples[match(names(files), samples$id),]
print("Second")
print(files)
print(samples)

# Check that matching worked
samples$names <- names(files)
samples$files <- files
stopifnot(all(samples$id == samples$names))
stopifnot(all(samples$names == basename(dirname(samples$files))))
print(samples)

# Save SummarizedExperiment with tx-level data for future use
se <- tximeta(samples)
se <- addIds(se, "SYMBOL", gene = T)
saveRDS(se, snakemake@output[['se_tx']])

# Collapse to gene-level for DE analysis
gse <- summarizeToGene(se)
gse <- addIds(gse, "SYMBOL", gene = T)
f <- as.formula(design_formula)

# No longer need se object - remove from memory
rm(se); gc()

## Gene-level
print("Building DESeq object for gene-level features")
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
