## ----setup, echo=FALSE--------------------------------------------------------
knitr::opts_chunk$set(collapse=TRUE, comment = "#>")
suppressPackageStartupMessages(library(universalmotif))

## -----------------------------------------------------------------------------
PPM <- function(C) C / sum(C)

## -----------------------------------------------------------------------------
PPMp <- function(C, p) (C + p / length(C)) / (sum(C) + p)

## -----------------------------------------------------------------------------
S <- function(C, B) log2(PPM(C) / B)

## ----logo1,fig.cap="Sequence logo of a Position Probability Matrix",echo=FALSE,fig.height=2.5,fig.width=5----
motif <- create_motif(c("AAGAAT", "ATCATA", "AAGTAA", "AACAAA", "ATTAAA",
                        "AAGAAT"), type = "PPM", pseudocount = 0)
view_motifs(motif, use.type="PPM") 

## ----logo2,fig.cap="Sequence logo of an Information Content Matrix",echo=FALSE,fig.height=2.5,fig.width=5----
motif <- create_motif(c("AAGAAT", "ATCATA", "AAGTAA", "AACAAA", "ATTAAA",
                        "AAGAAT"), type = "PPM", pseudocount = 0)
view_motifs(motif) 

## -----------------------------------------------------------------------------
ICtotal <- function(C) log2(length(C))

## -----------------------------------------------------------------------------
U <- function(C) -sum(PPM(C) * log2(PPM(C)), na.rm = TRUE)

## -----------------------------------------------------------------------------
ICfinal <- function(C) ICtotal(C) - U(C)

## -----------------------------------------------------------------------------
IC <- function(C) PPM(C) * ICfinal(C)

## -----------------------------------------------------------------------------
IC <- function(C, B) pmax(PPM(C) * log2(PPM(C) / B), 0)

