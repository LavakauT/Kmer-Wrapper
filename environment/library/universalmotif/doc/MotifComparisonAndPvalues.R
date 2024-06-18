## ----setup, echo=FALSE--------------------------------------------------------
knitr::opts_chunk$set(collapse=TRUE, comment = "#>")
suppressPackageStartupMessages(library(universalmotif))
suppressPackageStartupMessages(library(Biostrings))
suppressMessages(suppressPackageStartupMessages(library(MotifDb)))
suppressMessages(suppressPackageStartupMessages(library(ggplot2)))
suppressMessages(suppressPackageStartupMessages(library(ggtree)))
suppressMessages(suppressPackageStartupMessages(library(dplyr)))
data(ArabidopsisPromoters)
data(ArabidopsisMotif)
motdb <- convert_motifs(MotifDb)

## -----------------------------------------------------------------------------
EUCL <- function(c1, c2) {
  sqrt( sum( (c1 - c2)^2 ) )
}

WEUCL <- function(c1, c2, bkg1, bkg2) {
  sqrt( sum( (bkg1 + bkg2) * (c1 - c2)^2 ) )
}

KL <- function(c1, c2) {
  ( sum(c1 * log(c1 / c2)) + sum(c2 * log(c2 / c1)) ) / 2
}

HELL <- function(c1, c2) {
  sqrt( sum( ( sqrt(c1) - sqrt(c2) )^2 ) ) / sqrt(2)
}

SEUCL <- function(c1, c2) {
  sum( (c1 - c2)^2 )
}

MAN <- function(c1, c2) {
  sum ( abs(c1 - c2) )
}

PCC <- function(c1, c2) {
  n <- length(c1)
  top <- n * sum(c1 * c2) - sum(c1) * sum(c2)
  bot <- sqrt( ( n * sum(c1^2) - sum(c1)^2 ) * ( n * sum(c2^2) - sum(c2)^2 ) )
  top / bot
}

WPCC <- function(c1, c2, bkg1, bkg2) {
  weights <- bkg1 + bkg2
  mean1 <- sum(weights * c1)
  mean2 <- sum(weights * c2)
  var1 <- sum(weights * (c1 - mean1)^2)
  var2 <- sum(weights * (c2 - mean2)^2)
  cov <- sum(weights * (c1 - mean1) * (c2 - mean2))
  cov / sqrt(var1 * var2)
}

SW <- function(c1, c2) {
  2 - sum( (c1 - c2)^2 )
}

ALLR <- function(c1, c2, bkg1, bkg2, nsites1, nsites2) {
  left <- sum( c2 * nsites2 * log(c1 / bkg1) )
  right <- sum( c1 * nsites1 * log(c2 / bkg2) )
  ( left + right ) / ( nsites1 + nsites2 )
}

BHAT <- function(c1, c2) {
  sum( sqrt(c1 * c2) )
}

## -----------------------------------------------------------------------------
c1 <- c(0.7, 0.1, 0.1, 0.1)
c2 <- c(0.5, 0.0, 0.2, 0.3)

compare_columns(c1, c2, "PCC")
compare_columns(c1, c2, "EUCL")

## ----echo=FALSE,fig.wide=TRUE,fig.asp=1,fig.cap="\\label{fig:fig1}Distributions of scores from approximately 500 random motif and individual column comparisons"----
atm <- filter_motifs(motdb, organism = "Athaliana")
pool <- do.call(cbind, atm)@motif
pool <- pool + 0.01
metrics <- universalmotif:::COMPARE_METRICS

res <- vector("list", length(metrics))
names(res) <- metrics

res2 <- vector("list", length(metrics))
names(res2) <- metrics
res3 <- vector("list", length(metrics))
names(res3) <- metrics
res4 <- vector("list", length(metrics))
names(res4) <- metrics
res4 <- res4[!names(res4) %in% c("PCC", "ALLR", "ALLR_LL")]
res5 <- vector("list", length(metrics))
names(res5) <- metrics
res9 <- vector("list", length(metrics))
names(res9) <- metrics
res8 <- vector("list", length(metrics))
names(res8) <- metrics
res8 <- res8[!names(res8) %in% c("PCC", "ALLR", "ALLR_LL")]
res99 <- vector("list", length(metrics))
names(res99) <- metrics

y1 <- atm[sample(1:length(atm), 33)]
x1 <- pool[, sample(1:ncol(pool), 528)]
x2 <- pool[, sample(1:ncol(pool), 528)]

for (m in metrics) {
  res[[m]] <- numeric(528)
  for (i in 1:528) {
    res[[m]][i] <- compare_columns(x1[, i], x2[, i], m)
  }
  res2[[m]] <- as.numeric(as.dist(suppressWarnings(compare_motifs(y1, method=m, score.strat="sum", min.mean.ic=0))))
  res3[[m]] <- as.numeric(as.dist(suppressWarnings(compare_motifs(y1, method=m, score.strat="a.mean", min.mean.ic=0))))
  if (!m %in% c("PCC", "ALLR", "ALLR_LL")) {
    res4[[m]] <- as.numeric(as.dist(suppressWarnings(compare_motifs(y1, method=m, score.strat="g.mean", min.mean.ic=0))))
    res8[[m]] <- as.numeric(as.dist(suppressWarnings(compare_motifs(y1, method=m, score.strat="wg.mean", min.mean.ic=0))))
  }
  res5[[m]] <- as.numeric(as.dist(suppressWarnings(compare_motifs(y1, method=m, score.strat="median", min.mean.ic=0))))
  res9[[m]] <- as.numeric(as.dist(suppressWarnings(compare_motifs(y1, method=m, score.strat="wa.mean", min.mean.ic=0))))
  res99[[m]] <- as.numeric(as.dist(suppressWarnings(compare_motifs(y1, method=m, score.strat="fzt", min.mean.ic=0))))
}

l_2_df <- function(x) {

  d <- data.frame(key = character(), scores = numeric(), stringsAsFactors = FALSE)

  for (i in seq_along(x)) {
    y <- x[[i]]
    y <- y[!is.na(y)]
    d <- rbind(d, data.frame(key = rep(names(x)[i], length(y)), scores = y,
                             stringsAsFactors = FALSE))
  }

  d

}

res <- l_2_df(res)
res2 <- l_2_df(res2)
res3 <- l_2_df(res3)
res4 <- l_2_df(res4)
res5 <- l_2_df(res5)
res9 <- l_2_df(res9)
res8 <- l_2_df(res8)
res99 <- l_2_df(res99)

res$type <- "RawColumnScores"
res2$type <- "Sum"
res3$type <- "ArithMean"
res4$type <- "GeoMean"
res5$type <- "Median"
res9$type <- "WeightedArithMean"
res99$type <- "FisherZTrans"
res8$type <- "WeightedGeoMean"
res6 <- rbind(res, res3, res4, res5, res9, res8, res99)

dres <- res6 %>%
          group_by(key, type) %>%
          summarise(mean = mean(scores, na.rm = TRUE))

ggplot(res6, aes(x = scores, fill = type)) +
  geom_density(alpha = 0.3, adjust = 2) +
  geom_vline(aes(xintercept = 0), colour = "black", linetype = "dashed") +
  # geom_vline(data=dres, aes(xintercept=mean, colour = type)) +
  facet_wrap(key ~ ., ncol = 3, scales = "free") +
  theme_minimal() +
  theme(text = element_text(family = "Times"))

## ----echo=FALSE,fig.cap="\\label{fig:fig2}Example scores from comparing two motifs",fig.height=4,fig.width=5----
library(universalmotif)
library(MotifDb)

motifs <- convert_motifs(MotifDb)
motifs <- filter_motifs(motifs, altname = c("M0003_1.02", "M0004_1.02"))
# summarise_motifs(motifs)
view_motifs(motifs) +
  theme(text = element_text(family = "Times"))
try_all <- function(motifs, ...) {
  scores <- vector("list", 4)
  methods <- c("PCC", "WPCC", "EUCL", "SW", "KL", "ALLR",
               "BHAT", "HELL", "WEUCL", "SEUCL", "MAN", "ALLR_LL")
  for (i in seq_along(methods)) {
    scores[[1]][i] <- compare_motifs(motifs, method = methods[i])[1, 2]
    scores[[2]][i] <- compare_motifs(motifs, method = methods[i], normalise.scores = TRUE)[1, 2]
    scores[[3]][i] <- compare_motifs(motifs, method = methods[i], min.overlap = 99)[1, 2]
    scores[[4]][i] <- compare_motifs(motifs, method = methods[i], min.position.ic=0.25)[1, 2]
  }
  res <- data.frame(type = c("similarity", "similarity", "distance", 
                             "similarity", "distance", "similarity",
                             "similarity", "distance",
                             "distance", "distance", 
                             "distance", "similarity"),
                    method = methods,
                    default = scores[[1]],
                    normalised = scores[[2]],
                    # noRC = scores[[3]],
                    checkIC = scores[[4]])
  knitr::kable(res, format = "markdown", caption = "Comparing two motifs with various settings")
}
try_all(motifs)

## -----------------------------------------------------------------------------
library(universalmotif)
library(MotifDb)
motifs <- filter_motifs(MotifDb, organism = "Athaliana")

# Compare the first motif with everything and return P-values
head(compare_motifs(motifs, 1))

## -----------------------------------------------------------------------------
library(universalmotif)
library(MotifDb)

motifs <- filter_motifs(MotifDb, family = c("AP2", "B3", "bHLH", "bZIP",
                                            "AT hook"))
motifs <- motifs[sample(seq_along(motifs), 100)]
tree <- motif_tree(motifs, layout = "daylight", linecol = "family")

## Make some changes to the tree in regular ggplot2 fashion:
# tree <- tree + ...

tree

## -----------------------------------------------------------------------------
library(universalmotif)
library(MotifDb)
library(ggtree)
library(ggplot2)

motifs <- convert_motifs(MotifDb)
motifs <- filter_motifs(motifs, organism = "Athaliana")
motifs <- motifs[sample(seq_along(motifs), 25)]

## Step 1: compare motifs

comparisons <- compare_motifs(motifs, method = "PCC", min.mean.ic = 0,
                              score.strat = "a.mean")

## Step 2: create a "dist" object

# The current metric, PCC, is a similarity metric
comparisons <- 1 - comparisons

comparisons <- as.dist(comparisons)

# We also want to extract names from the dist object to match annotations
labels <- attr(comparisons, "Labels")

## Step 3: get the comparisons ready for tree-building

# The R package "ape" provides the necessary "as.phylo" function
comparisons <- ape::as.phylo(hclust(comparisons))

## Step 4: incorporate annotation data to colour tree lines

family <- sapply(motifs, function(x) x["family"])
family.unique <- unique(family)

# We need to create a list with an entry for each family; within each entry
# are the names of the motifs belonging to that family
family.annotations <- list()
for (i in seq_along(family.unique)) {
  family.annotations <- c(family.annotations,
                          list(labels[family %in% family.unique[i]]))
}
names(family.annotations) <- family.unique

# Now add the annotation data:
comparisons <- ggtree::groupOTU(comparisons, family.annotations)

## Step 5: draw the tree

tree <- ggtree(comparisons, aes(colour = group), layout = "rectangular") +
          theme(legend.position = "bottom", legend.title = element_blank())

## Step 6: add additional annotations

# If we wish, we can additional annotations such as tip labelling and size

# Tip labels:
tree <- tree + geom_tiplab()

# Tip size:
tipsize <- data.frame(label = labels,
                      icscore = sapply(motifs, function(x) x["icscore"]))

tree <- tree %<+% tipsize + geom_tippoint(aes(size = icscore))


## -----------------------------------------------------------------------------
library(universalmotif)
library(MotifDb)
library(cowplot)

## Get our starting set of motifs:
motifs <- convert_motifs(MotifDb[1:10])

## Get the tree: make sure it's a horizontal type layout
tree <- motif_tree(motifs, layout = "rectangular", linecol = "none")

## Now, make sure we order our list of motifs to match the order of tips:
mot.names <- sapply(motifs, function(x) x["name"])
names(motifs) <- mot.names
new.order <- tree$data$label[tree$data$isTip]
new.order <- rev(new.order[order(tree$data$y[tree$data$isTip])])
motifs <- motifs[new.order]

## Plot the two together (finessing of margins and positions may be required):
plot_grid(nrow = 1, rel_widths = c(1, -0.15, 1),
  tree + xlab(""), NULL,
  view_motifs(motifs, names.pos = "right") +
    ylab(element_blank()) +
    theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      axis.text = element_text(colour = "white")
    )
)

## -----------------------------------------------------------------------------
library(universalmotif)

m <- matrix(c(0.10,0.27,0.23,0.19,0.29,0.28,0.51,0.12,0.34,0.26,
              0.36,0.29,0.51,0.38,0.23,0.16,0.17,0.21,0.23,0.36,
              0.45,0.05,0.02,0.13,0.27,0.38,0.26,0.38,0.12,0.31,
              0.09,0.40,0.24,0.30,0.21,0.19,0.05,0.30,0.31,0.08),
            byrow = TRUE, nrow = 4)
motif <- create_motif(m, alphabet = "DNA", type = "PWM")
motif

## -----------------------------------------------------------------------------
data(ArabidopsisPromoters)

res <- scan_sequences(motif, ArabidopsisPromoters, verbose = 0,
  calc.pvals = FALSE, threshold = 0.8, threshold.type = "logodds")
head(res)

## -----------------------------------------------------------------------------
## One of the matches was CTCTAGAGAC, with a score of 5.869 (max possible = 6.531)

bkg <- get_bkg(ArabidopsisPromoters, 1)
bkg <- structure(bkg$probability, names = bkg$klet)
bkg

## -----------------------------------------------------------------------------
hit.prob <- bkg["A"]^3 * bkg["C"]^3 * bkg["G"]^2 * bkg["T"]^2
hit.prob <- unname(hit.prob)
hit.prob

## -----------------------------------------------------------------------------
res <- res[1:6, ]
pvals <- motif_pvalue(motif, res$score, bkg.probs = bkg)
res2 <- data.frame(motif=res$motif,match=res$match,pval=pvals)[order(pvals), ]
knitr::kable(res2, digits = 22, row.names = FALSE, format = "markdown")

## -----------------------------------------------------------------------------
res$score
motif_pvalue(motif, pvalue = pvals, bkg.probs = bkg)

## -----------------------------------------------------------------------------
data(ArabidopsisMotif)
ArabidopsisMotif

## -----------------------------------------------------------------------------
motif_range(ArabidopsisMotif)

## -----------------------------------------------------------------------------
(pvals2 <- motif_pvalue(ArabidopsisMotif, score = motif_range(ArabidopsisMotif)))

## -----------------------------------------------------------------------------
motif_pvalue(ArabidopsisMotif, pvalue = pvals2)

## -----------------------------------------------------------------------------
motif_pvalue(ArabidopsisMotif, score = c(-200, 100))

## -----------------------------------------------------------------------------
data(examplemotif2)
examplemotif2["multifreq"]["2"]
motif_range(examplemotif2, use.freq = 2)
motif_pvalue(examplemotif2, score = 15, use.freq = 2)
motif_pvalue(examplemotif2, pvalue = 0.00001, use.freq = 2)

## -----------------------------------------------------------------------------
(m <- create_motif(alphabet = "QWERTY"))
motif_pvalue(m, pvalue = c(1, 0.1, 0.001, 0.0001, 0.00001))

## -----------------------------------------------------------------------------
res <- res[1:6, ]
pvals <- motif_pvalue(motif, res$score, bkg.probs = bkg, method = "e")
res2 <- data.frame(motif=res$motif,match=res$match,pval=pvals)[order(pvals), ]
knitr::kable(res2, digits = 22, row.names = FALSE, format = "markdown")

## -----------------------------------------------------------------------------
scores <- c(-6, -3, 0, 3, 6)
k <- c(2, 4, 6, 8, 10)
out <- data.frame(k = c(2, 4, 6, 8, 10),
                  score.minus6 = rep(0, 5),
                  score.minus3 = rep(0, 5),
                  score.0 = rep(0, 5),
                  score.3 = rep(0, 5),
                  score.6 = rep(0, 5))

for (i in seq_along(scores)) {
  for (j in seq_along(k)) {
    out[j, i + 1] <- motif_pvalue(motif, scores[i], k = k[j], bkg.probs = bkg,
      method = "e")
  }
}

knitr::kable(out, format = "markdown", digits = 10)

## -----------------------------------------------------------------------------
bkg <- c(A=0.25, C=0.25, G=0.25, T=0.25)
pvals <- c(0.1, 0.01, 0.001, 0.0001, 0.00001)
scores <- motif_pvalue(motif, pvalue = pvals, bkg.probs = bkg, k = 10,
  method = "e")

scores.approx6 <- motif_pvalue(motif, pvalue = pvals, bkg.probs = bkg, k = 6,
  method = "e")
scores.approx8 <- motif_pvalue(motif, pvalue = pvals, bkg.probs = bkg, k = 8,
  method = "e")

pvals.exact <- motif_pvalue(motif, score = scores, bkg.probs = bkg, k = 10,
  method = "e")

pvals.approx6 <- motif_pvalue(motif, score = scores, bkg.probs = bkg, k = 6,
  method = "e")
pvals.approx8 <- motif_pvalue(motif, score = scores, bkg.probs = bkg, k = 8,
  method = "e")

res <- data.frame(pvalue = pvals, score = scores,
                  pvalue.exact = pvals.exact,
                  pvalue.k6 = pvals.approx6,
                  pvalue.k8 = pvals.approx8,
                  score.k6 = scores.approx6,
                  score.k8 = scores.approx8)
knitr::kable(res, format = "markdown", digits = 22)

## -----------------------------------------------------------------------------
bkg <- c(A=0.25, C=0.25, G=0.25, T=0.25)
scores <- -2:6
pvals <- motif_pvalue(motif, score = scores, bkg.probs = bkg, k = 10,
  method = "e")

scores.exact <- motif_pvalue(motif, pvalue = pvals, bkg.probs = bkg, k = 10,
  method = "e")

scores.approx6 <- motif_pvalue(motif, pvalue = pvals, bkg.probs = bkg, k = 6,
  method = "e")
scores.approx8 <- motif_pvalue(motif, pvalue = pvals, bkg.probs = bkg, k = 8,
  method = "e")

pvals.approx6 <- motif_pvalue(motif, score = scores, bkg.probs = bkg, k = 6,
  method = "e")
pvals.approx8 <- motif_pvalue(motif, score = scores, bkg.probs = bkg, k = 8,
  method = "e")

res <- data.frame(score = scores, pvalue = pvals,
                  pvalue.k6 = pvals.approx6,
                  pvalue.k8 = pvals.approx8,
                  score.exact = scores.exact,
                  score.k6 = scores.approx6,
                  score.k8 = scores.approx8)
knitr::kable(res, format = "markdown", digits = 22)

## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()

