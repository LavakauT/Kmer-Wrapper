## ----setup, echo=FALSE--------------------------------------------------------
knitr::opts_chunk$set(collapse=TRUE,comment = "#>")
suppressPackageStartupMessages(library(universalmotif))
suppressMessages(suppressPackageStartupMessages(library(MotifDb)))
suppressMessages(suppressPackageStartupMessages(library(TFBSTools)))
data(examplemotif)
data(MA0003.2)

## -----------------------------------------------------------------------------
library(universalmotif)
data(examplemotif)
examplemotif

## -----------------------------------------------------------------------------
library(universalmotif)
data(examplemotif)

## The various slots can be accessed individually using `[`

examplemotif["consensus"]

## To change a slot, use `[<-`

examplemotif["family"] <- "My motif family"
examplemotif

## ----error=TRUE---------------------------------------------------------------
library(universalmotif)
data(examplemotif)

## The consensus slot is dependent on the motif matrix

examplemotif["consensus"]

## Changing this would mean it no longer matches the motif

examplemotif["consensus"] <- "GGGAGAG"

## Another example of trying to change a protected slot:

examplemotif["strand"] <- "x"

## -----------------------------------------------------------------------------
library(universalmotif)
library(MotifDb)
data(examplemotif)
data(MA0003.2)

## convert from a `universalmotif` motif to another class

convert_motifs(examplemotif, "TFBSTools-PWMatrix")

## convert to universalmotif

convert_motifs(MA0003.2)

## convert between two packages

convert_motifs(MotifDb[1], "TFBSTools-ICMatrix")

## -----------------------------------------------------------------------------
library(universalmotif)
data(examplemotif)

## This motif is currently a PPM:

examplemotif["type"]

## -----------------------------------------------------------------------------
convert_type(examplemotif, "PCM")

## -----------------------------------------------------------------------------
examplemotif["pseudocount"]
convert_type(examplemotif, "PWM")

## -----------------------------------------------------------------------------
convert_type(examplemotif, "PWM", pseudocount = 1)

## -----------------------------------------------------------------------------
examplemotif["nsites"] <- 10
convert_type(examplemotif, "ICM", nsize_correction = FALSE)

convert_type(examplemotif, "ICM", nsize_correction = TRUE)

examplemotif["bkg"] <- c(A = 0.4, C = 0.1, G = 0.1, T = 0.4)
convert_type(examplemotif, "ICM", relative_entropy = TRUE)

## ----fig.height=4, fig.width=5------------------------------------------------
library(universalmotif)

m1 <- create_motif("TTAAACCCC", name = "1")
m2 <- create_motif("AACC", name = "2")
m3 <- create_motif("AACCCCGG", name = "3")

view_motifs(c(m1, m2, m3),
  show.positions.once = FALSE, show.names = FALSE)

## ----fig.height=2,fig.width=5-------------------------------------------------
view_motifs(merge_motifs(c(m1, m2, m3), method = "PCC"))

## -----------------------------------------------------------------------------
library(universalmotif)
library(MotifDb)

motifs <- filter_motifs(MotifDb, family = "bHLH")[1:100]
length(motifs)

motifs <- merge_similar(motifs)
length(motifs)

## -----------------------------------------------------------------------------
library(universalmotif)
data(examplemotif)

## Quickly switch to the reverse complement of a motif

## Original:

examplemotif

## Reverse complement:

motif_rc(examplemotif)

## -----------------------------------------------------------------------------
library(universalmotif)
data(examplemotif)

## DNA --> RNA

switch_alph(examplemotif)

## RNA --> DNA

motif <- create_motif(alphabet = "RNA")
motif

switch_alph(motif)

## -----------------------------------------------------------------------------
library(universalmotif)

motif <- create_motif("NNGCSGCGGNN")
motif

trim_motifs(motif)
trim_motifs(motif, trim.from = "right")

## ----fig.height=3.5, fig.width=5----------------------------------------------
motif1 <- create_motif("ATCGATGC", pseudocount = 10, type = "PPM", nsites = 100)
motif2 <- round_motif(motif1)
view_motifs(c(motif1, motif2))

## -----------------------------------------------------------------------------
motif.matrix <- matrix(c(0.7, 0.1, 0.1, 0.1,
                         0.7, 0.1, 0.1, 0.1,
                         0.1, 0.7, 0.1, 0.1,
                         0.1, 0.7, 0.1, 0.1,
                         0.1, 0.1, 0.7, 0.1,
                         0.1, 0.1, 0.7, 0.1,
                         0.1, 0.1, 0.1, 0.7,
                         0.1, 0.1, 0.1, 0.7), nrow = 4)

motif <- create_motif(motif.matrix, alphabet = "RNA", name = "My motif",
                      pseudocount = 1, nsites = 20, strand = "+")

## The 'type', 'icscore' and 'consensus' slots will be filled for you

motif

## -----------------------------------------------------------------------------
motif <- create_motif("CCNSNGG", nsites = 50, pseudocount = 1)

## Now to disk:
## write_meme(motif, "meme_motif.txt")

motif

## -----------------------------------------------------------------------------
create_motif()

## -----------------------------------------------------------------------------
create_motif(bkg = c(A = 0.2, C = 0.4, G = 0.2, T = 0.2))

## -----------------------------------------------------------------------------
create_motif(alphabet = "QWERTY")

## ----fig.height=2, fig.width=5------------------------------------------------
library(universalmotif)
data(examplemotif)

## With the native `view_motifs` function:
view_motifs(examplemotif)

## ----fig.height=2.5, fig.width=5----------------------------------------------
view_motifs(create_motif(15)) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)
  )

## ----fig.height=2.5, fig.width=5----------------------------------------------
## For all the following examples, simply passing the functions a PPM is
## sufficient
motif <- convert_type(examplemotif, "PPM")
## Only need the matrix itself
motif <- motif["motif"]

## seqLogo:
seqLogo::seqLogo(motif)

## motifStack:
motifStack::plotMotifLogo(motif)

## ggseqlogo:
ggseqlogo::ggseqlogo(motif)

## ----fig.height=5, fig.width=5------------------------------------------------
library(universalmotif)
library(MotifDb)

motifs <- convert_motifs(MotifDb[50:54])
view_motifs(motifs, show.positions.once = FALSE, names.pos = "right")

## ----fig.height=2.5, fig.width=5----------------------------------------------
library(universalmotif)
data(examplemotif)

## Start from a numeric matrix:
toplot <- examplemotif["motif"]

# Adjust the character heights as you wish (negative values are possible):
toplot[4] <- 2
toplot[20] <- -0.5

# Mix and match the number of characters per letter/position:
rownames(toplot)[1] <- "AA"

toplot <- toplot[c(1, 4), ]

toplot

view_logo(toplot)

## -----------------------------------------------------------------------------
library(universalmotif)

motif <- create_motif("CWWWWCC", nsites = 6)
sequences <- DNAStringSet(rep(c("CAAAACC", "CTTTTCC"), 3))
motif.k2 <- add_multifreq(motif, sequences, add.k = 2)

## Alternatively:
# motif.k2 <- create_motif(sequences, add.multifreq = 2)

motif.k2

## ----fig.height=2.5, fig.width=5----------------------------------------------
view_motifs(motif.k2, use.freq = 2)

## -----------------------------------------------------------------------------
library(universalmotif)
library(MotifDb)

## Obtain a `universalmotif_df` object
motifs <- to_df(MotifDb)
head(motifs)

## -----------------------------------------------------------------------------
library(dplyr)

motifs <- motifs %>%
  mutate(bkg = case_when(
    organism == "Athaliana" ~ list(c(A = 0.32, C = 0.18, G = 0.18, T = 0.32)),
    TRUE ~ list(c(A = 0.25, C = 0.25, G = 0.25, T = 0.25))
  ))
head(filter(motifs, organism == "Athaliana"))

## -----------------------------------------------------------------------------
motifs <- motifs %>%
  mutate(MotifIndex = 1:n())
head(motifs)

to_list(motifs)[[1]]

## -----------------------------------------------------------------------------
to_list(motifs, extrainfo = FALSE)[[1]]

## -----------------------------------------------------------------------------
library(universalmotif)

get_consensus(c(A = 0.7, C = 0.1, G = 0.1, T = 0.1))

consensus_to_ppm("G")

## -----------------------------------------------------------------------------
library(universalmotif)
library(MotifDb)

## Let us extract all of the Arabidopsis and C. elegans motifs 

motifs <- filter_motifs(MotifDb, organism = c("Athaliana", "Celegans"))

## Only keeping motifs with sufficient information content and length:

motifs <- filter_motifs(motifs, icscore = 10, width = 10)

head(summarise_motifs(motifs))

## -----------------------------------------------------------------------------
library(universalmotif)
data(examplemotif)

sample_sites(examplemotif)

## -----------------------------------------------------------------------------
library(universalmotif)
library(MotifDb)

motifs <- convert_motifs(MotifDb[1:50])
head(summarise_motifs(motifs))

motifs.shuffled <- shuffle_motifs(motifs, k = 3)
head(summarise_motifs(motifs.shuffled))

## -----------------------------------------------------------------------------
library(universalmotif)
data(examplemotif)
examplemotif

## Get the min and max possible scores:
motif_score(examplemotif)

## Show matches above a score of 10:
get_matches(examplemotif, 10)

## Get the probability of a match:
prob_match(examplemotif, "TTTTTTT", allow.zero = FALSE)

## Score a specific sequence:
score_match(examplemotif, "TTTTTTT")

## Take a look at the distribution of scores:
plot(density(get_scores(examplemotif), bw = 5))

## -----------------------------------------------------------------------------
library(universalmotif)

m <- create_motif(type = "PCM")["motif"]
m

apply(m, 2, pcm_to_ppm)

## -----------------------------------------------------------------------------
library(universalmotif)

position_icscore(c(0.7, 0.1, 0.1, 0.1))

## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()

