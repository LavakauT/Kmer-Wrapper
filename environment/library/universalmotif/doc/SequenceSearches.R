## ----setup, echo=FALSE--------------------------------------------------------
knitr::opts_chunk$set(collapse=TRUE, comment = "#>")
suppressPackageStartupMessages(library(universalmotif))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(ggbio))
data(ArabidopsisPromoters)
data(ArabidopsisMotif)
options(Biostrings.coloring = FALSE)

## -----------------------------------------------------------------------------
library(universalmotif)
library(Biostrings)

## Create some DNA sequences for use with an external program (default
## is DNA):

sequences.dna <- create_sequences(seqnum = 500,
                                  freqs = c(A=0.3, C=0.2, G=0.2, T=0.3))
## writeXStringSet(sequences.dna, "dna.fasta")
sequences.dna

## Amino acid:

create_sequences(alphabet = "AA")

## Any set of characters can be used

create_sequences(alphabet = paste0(letters, collapse = ""))

## -----------------------------------------------------------------------------
library(universalmotif)

## Background of DNA sequences:
dna <- create_sequences()
get_bkg(dna, k = 1:2)

## Background of non DNA/RNA sequences:
qwerty <- create_sequences("QWERTY")
get_bkg(qwerty, k = 1:2)

## ----fig.height=4.5,fig.width=4.5---------------------------------------------
library(universalmotif)

## Generate three random sets of sequences:
s1 <- create_sequences(seqnum = 20,
  freqs = c(A = 0.3, C = 0.2, G = 0.2, T = 0.3))
s2 <- create_sequences(seqnum = 20,
  freqs = c(A = 0.4, C = 0.4, G = 0.1, T = 0.1))
s3 <- create_sequences(seqnum = 20,
  freqs = c(A = 0.2, C = 0.3, G = 0.3, T = 0.2))

## Create a function to get properly formatted k-let counts:
get_klet_matrix <- function(seqs, k, groupName) {
  bkg <- get_bkg(seqs, k = k, merge.res = FALSE)
  bkg <- bkg[, c("sequence", "klet", "count")]
  bkg <- reshape(bkg, idvar = "sequence", timevar = "klet",
    direction = "wide")
  as.data.frame(cbind(Group = groupName, bkg))
}

## Calculate k-let content (up to you what size k you want!):
s1 <- get_klet_matrix(s1, 4, 1)
s2 <- get_klet_matrix(s2, 4, 2)
s3 <- get_klet_matrix(s3, 4, 3)

# Combine everything into a single object:
sAll <- rbind(s1, s2, s3)

## Do the PCA:
sPCA <- prcomp(sAll[, -(1:2)])

## Plot the PCA:
plot(sPCA$x, col = c("red", "forestgreen", "blue")[sAll$Group], pch = 19)

## -----------------------------------------------------------------------------
library(universalmotif)
library(Biostrings)
data(ArabidopsisPromoters)

## Potentially starting off with some external sequences:
# ArabidopsisPromoters <- readDNAStringSet("ArabidopsisPromoters.fasta")

euler <- shuffle_sequences(ArabidopsisPromoters, k = 2, method = "euler")
markov <- shuffle_sequences(ArabidopsisPromoters, k = 2, method = "markov")
linear <- shuffle_sequences(ArabidopsisPromoters, k = 2, method = "linear")
k1 <- shuffle_sequences(ArabidopsisPromoters, k = 1)

## -----------------------------------------------------------------------------
o.letter <- get_bkg(ArabidopsisPromoters, 1)
e.letter <- get_bkg(euler, 1)
m.letter <- get_bkg(markov, 1)
l.letter <- get_bkg(linear, 1)

data.frame(original=o.letter$count, euler=e.letter$count,
  markov=m.letter$count, linear=l.letter$count, row.names = DNA_BASES)

o.counts <- get_bkg(ArabidopsisPromoters, 2)
e.counts <- get_bkg(euler, 2)
m.counts <- get_bkg(markov, 2)
l.counts <- get_bkg(linear, 2)

data.frame(original=o.counts$count, euler=e.counts$count,
  markov=m.counts$count, linear=l.counts$count,
  row.names = get_klets(DNA_BASES, 2))

## -----------------------------------------------------------------------------
library(Biostrings)
library(universalmotif)
library(ggplot2)

myseq <- DNAStringSet(paste0(
  create_sequences(seqlen = 500, freqs = c(A=0.4, T=0.4, C=0.1, G=0.1)),
  create_sequences(seqlen = 500)
))

myseq_shuf <- shuffle_sequences(myseq)
myseq_shuf_local <- shuffle_sequences(myseq, window = TRUE)

myseq_bkg <- get_bkg(myseq, k = 1, window = TRUE)
myseq_shuf_bkg <- get_bkg(myseq_shuf, k = 1, window = TRUE)
myseq_shuf_local_bkg <- get_bkg(myseq_shuf_local, k = 1, window = TRUE)

myseq_bkg$group <- "original"
myseq_shuf_bkg$group <- "shuffled"
myseq_shuf_local_bkg$group <- "shuffled-local"

myseq_all <- as.data.frame(
  rbind(myseq_bkg, myseq_shuf_bkg, myseq_shuf_local_bkg)
)

ggplot(myseq_all, aes(x = start, y = probability, colour = klet)) +
  geom_line() +
  theme_minimal() +
  scale_colour_manual(values = universalmotif:::DNA_COLOURS) +
  xlab(element_blank()) +
  facet_wrap(~group, ncol = 1)


## -----------------------------------------------------------------------------
library(universalmotif)
m1 <- create_motif("TATATATATA", nsites = 50, type = "PWM", pseudocount = 1)
m2 <- matrix(c(0.10,0.27,0.23,0.19,0.29,0.28,0.51,0.12,0.34,0.26,
               0.36,0.29,0.51,0.38,0.23,0.16,0.17,0.21,0.23,0.36,
               0.45,0.05,0.02,0.13,0.27,0.38,0.26,0.38,0.12,0.31,
               0.09,0.40,0.24,0.30,0.21,0.19,0.05,0.30,0.31,0.08),
             byrow = TRUE, nrow = 4)
m2 <- create_motif(m2, alphabet = "DNA", type = "PWM")
m1["motif"]
m2["motif"]

## -----------------------------------------------------------------------------
motif_pvalue(m2, pvalue = 0.001)

## -----------------------------------------------------------------------------
library(universalmotif)
data(ArabidopsisMotif)
data(ArabidopsisPromoters)

(Example.Score <- score_match(ArabidopsisMotif, "TTCTCTTTTTTTTTT"))
(Example.Pvalue <- motif_pvalue(ArabidopsisMotif, Example.Score))

(Max.Possible.Hits <- sum(width(ArabidopsisPromoters) - ncol(ArabidopsisMotif) + 1))

## -----------------------------------------------------------------------------
(Example.bonferroni <- Example.Pvalue * Max.Possible.Hits)

## -----------------------------------------------------------------------------
(Scan.Results <- scan_sequences(ArabidopsisMotif, ArabidopsisPromoters,
  threshold = 0.8, threshold.type = "logodds", calc.qvals = FALSE))

## -----------------------------------------------------------------------------
Pvalues <- Scan.Results$pvalue
Pvalues.Ranks <- (rank(Pvalues) / Max.Possible.Hits) * 100
Qvalues.BH <- Pvalues / Pvalues.Ranks
(Example.BH <- Qvalues.BH[Scan.Results$match == "TTCTCTTTTTTTTTT"][1])

## -----------------------------------------------------------------------------
Scan.Results <- Scan.Results[order(Scan.Results$score, decreasing = TRUE), ]
Observed.Hits <- 1:nrow(Scan.Results)
Null.Hits <- Max.Possible.Hits * Scan.Results$pvalue
Qvalues.fdr <- Null.Hits / Observed.Hits
Qvalues.fdr <- rev(cummin(rev(Qvalues.fdr)))
(Example.fdr <- Qvalues.fdr[Scan.Results$match == "TTCTCTTTTTTTTTT"][1])

## -----------------------------------------------------------------------------
knitr::kable(
  data.frame(
    What = c("Score", "P-value", "bonferroni", "BH", "fdr"),
    Value = format(
      c(Example.Score, Example.Pvalue, Example.bonferroni, Example.BH, Example.fdr),
      scientific = FALSE
    )
  ),
  format = "markdown", caption = "Comparing P-value correction methods"
)

## -----------------------------------------------------------------------------
library(universalmotif)
library(Biostrings)
data(ArabidopsisPromoters)

## A 2-letter example:

motif.k2 <- create_motif("CWWWWCC", nsites = 6)
sequences.k2 <- DNAStringSet(rep(c("CAAAACC", "CTTTTCC"), 3))
motif.k2 <- add_multifreq(motif.k2, sequences.k2)

## -----------------------------------------------------------------------------
scan_sequences(motif.k2, ArabidopsisPromoters, RC = TRUE,
               threshold = 0.9, threshold.type = "logodds")

## -----------------------------------------------------------------------------
scan_sequences(motif.k2, ArabidopsisPromoters, use.freq = 2, RC = TRUE,
               threshold = 0.9, threshold.type = "logodds")

## -----------------------------------------------------------------------------
motif <- create_motif("AAAAAA")

## Leave in overlapping hits:

scan_sequences(motif, ArabidopsisPromoters, RC = TRUE, threshold = 0.9,
               threshold.type = "logodds")

## Only keep the highest scoring hit amongst overlapping hits:

scan_sequences(motif, ArabidopsisPromoters, RC = TRUE, threshold = 0.9,
               threshold.type = "logodds", no.overlaps = TRUE)

## -----------------------------------------------------------------------------
scan_sequences(motif.k2, ArabidopsisPromoters, RC = TRUE,
               threshold = 0.9, threshold.type = "logodds",
               return.granges = TRUE)

## ----fig.height = 2.5, fig.width = 6------------------------------------------
library(universalmotif)
library(GenomicRanges)
library(ggbio)

data(ArabidopsisPromoters)

motif1 <- create_motif("AAAAAA", name = "Motif A")
motif2 <- create_motif("CWWWWCC", name = "Motif B")

res <- scan_sequences(c(motif1, motif2), ArabidopsisPromoters[1:10],
  return.granges = TRUE, calc.pvals = TRUE, no.overlaps = TRUE,
  threshold = 0.2, threshold.type = "logodds")

## Just plot the motif hits:
autoplot(res, layout = "karyogram", aes(fill = motif, color = motif)) +
  theme(
    strip.background = element_rect(fill = NA, colour = NA),
    panel.background = element_rect(fill = NA, colour = NA)
  )

## Plot Motif A hits by P-value:
autoplot(res[res$motif.i == 1, ], layout = "karyogram",
  aes(fill = log10(pvalue), colour = log10(pvalue))) +
  scale_fill_gradient(low = "black", high = "grey75") +
  scale_colour_gradient(low = "black", high = "grey75") +
  theme(
    strip.background = element_rect(fill = NA, colour = NA),
    panel.background = element_rect(fill = NA, colour = NA)
  )

## ----fig.height = 4.5, fig.width=6.5------------------------------------------
library(universalmotif)
library(ggplot2)

data(ArabidopsisMotif)
data(ArabidopsisPromoters)

res <- scan_sequences(ArabidopsisMotif, ArabidopsisPromoters,
  threshold = 0, threshold.type = "logodds.abs")
res <- as.data.frame(res)
res$x <- mapply(function(x, y) mean(c(x, y)), res$start, res$stop)

ggplot(res, aes(x, sequence, fill = score)) +
  scale_fill_viridis_c() +
  scale_x_continuous(expand = c(0, 0)) +
  xlim(0, 1000) +
  xlab(element_blank()) +
  ylab(element_blank()) +
  geom_tile(width = ncol(ArabidopsisMotif)) +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text.y = element_text(size = 6)) 

## ----fig.height=5,fig.width=6.5-----------------------------------------------
library(universalmotif)
library(ggplot2)

data(ArabidopsisMotif)
data(ArabidopsisPromoters)

res <- scan_sequences(ArabidopsisMotif, ArabidopsisPromoters[1:5],
  threshold = -Inf, threshold.type = "logodds.abs")
res <- as.data.frame(res)
res$position <- mapply(function(x, y) mean(c(x, y)), res$start, res$stop)

ggplot(res, aes(position, score, colour = score)) +
  geom_line() +
  geom_hline(yintercept = 0, colour = "red", alpha = 0.3) +
  theme_bw() +
  scale_colour_viridis_c() +
  facet_wrap(~sequence, ncol = 1) +
  xlab(element_blank()) +
  ylab(element_blank()) +
  theme(strip.background = element_blank())

## -----------------------------------------------------------------------------
library(universalmotif)
data(ArabidopsisMotif)
data(ArabidopsisPromoters)

enrich_motifs(ArabidopsisMotif, ArabidopsisPromoters, shuffle.k = 3,
              threshold = 0.001, RC = TRUE)

## -----------------------------------------------------------------------------
library(universalmotif)
data(ArabidopsisPromoters)

m1 <- create_motif("TTTATAT", name = "PartA")
m2 <- create_motif("GGTTCGA", name = "PartB")

## -----------------------------------------------------------------------------
m <- cbind(m1, m2)
m <- add_gap(m, gaploc = ncol(m1), mingap = 4, maxgap = 6)
m

## -----------------------------------------------------------------------------
scan_sequences(m, ArabidopsisPromoters, threshold = 0.4, threshold.type = "logodds")

## -----------------------------------------------------------------------------
library(universalmotif)
library(Biostrings)

Ex.seq <- DNAStringSet(c(
  A = "GTTGAAAAAAAAAAAAAAAACAGACGT",
  B = "TTAGATGGCCCATAGCTTATACGGCAA",
  C = "AATAAAATGCTTAGGAAATCGATTGCC"
))

## -----------------------------------------------------------------------------
mask_seqs(Ex.seq, "AAAAAAAA")

## -----------------------------------------------------------------------------
(Ex.DUST <- sequence_complexity(Ex.seq, window.size = 10, method = "DUST",
    return.granges = TRUE))

## -----------------------------------------------------------------------------
(Ex.DUST <- Ex.DUST[Ex.DUST$complexity >= 3])
mask_ranges(Ex.seq, Ex.DUST)

## ----eval=FALSE---------------------------------------------------------------
#  library(universalmotif)
#  
#  ## 1. Once per session: via `options()`
#  
#  options(meme.bin = "/path/to/meme/bin/meme")
#  
#  run_meme(...)
#  
#  ## 2. Once per run: via `run_meme()`
#  
#  run_meme(..., bin = "/path/to/meme/bin/meme")

## ----eval=FALSE---------------------------------------------------------------
#  library(universalmotif)
#  data(ArabidopsisPromoters)
#  
#  ## 1. Read sequences from disk (in fasta format):
#  
#  library(Biostrings)
#  
#  # The following `read*()` functions are available in Biostrings:
#  # DNA: readDNAStringSet
#  # DNA with quality scores: readQualityScaledDNAStringSet
#  # RNA: readRNAStringSet
#  # Amino acid: readAAStringSet
#  # Any: readBStringSet
#  
#  sequences <- readDNAStringSet("/path/to/sequences.fasta")
#  
#  run_meme(sequences, ...)
#  
#  ## 2. Extract from a `BSgenome` object:
#  
#  library(GenomicFeatures)
#  library(TxDb.Athaliana.BioMart.plantsmart28)
#  library(BSgenome.Athaliana.TAIR.TAIR9)
#  
#  # Let us retrieve the same promoter sequences from ArabidopsisPromoters:
#  gene.names <- names(ArabidopsisPromoters)
#  
#  # First get the transcript coordinates from the relevant `TxDb` object:
#  transcripts <- transcriptsBy(TxDb.Athaliana.BioMart.plantsmart28,
#                               by = "gene")[gene.names]
#  
#  # There are multiple transcripts per gene, we only care for the first one
#  # in each:
#  
#  transcripts <- lapply(transcripts, function(x) x[1])
#  transcripts <- unlist(GRangesList(transcripts))
#  
#  # Then the actual sequences:
#  
#  # Unfortunately this is a case where the chromosome names do not match
#  # between the two databases
#  
#  seqlevels(TxDb.Athaliana.BioMart.plantsmart28)
#  #> [1] "1"  "2"  "3"  "4"  "5"  "Mt" "Pt"
#  seqlevels(BSgenome.Athaliana.TAIR.TAIR9)
#  #> [1] "Chr1" "Chr2" "Chr3" "Chr4" "Chr5" "ChrM" "ChrC"
#  
#  # So we must first rename the chromosomes in `transcripts`:
#  seqlevels(transcripts) <- seqlevels(BSgenome.Athaliana.TAIR.TAIR9)
#  
#  # Finally we can extract the sequences
#  promoters <- getPromoterSeq(transcripts,
#                              BSgenome.Athaliana.TAIR.TAIR9,
#                              upstream = 1000, downstream = 0)
#  
#  run_meme(promoters, ...)

## ----eval=FALSE---------------------------------------------------------------
#  run_meme(sequences, output = "/path/to/desired/output/folder")

## ----eval=FALSE---------------------------------------------------------------
#  motifs <- run_meme(sequences)
#  motifs.k23 <- mapply(add_multifreq, motifs$motifs, motifs$sites)

## -----------------------------------------------------------------------------
library(universalmotif)

string <- "DASDSDDSASDSSA"

calc_complexity(string)

count_klets(string, 2)

shuffle_string(string, 2)

## -----------------------------------------------------------------------------
library(universalmotif)

calc_windows(n = 12, window = 4, overlap = 2)

get_klets(c("A", "S", "D"), 2)

slide_fun("ABCDEFGH", charToRaw, raw(2), window = 2, overlap = 1)

window_string("ABCDEFGH", window = 2, overlap = 1)

## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()

