### R code from vignette source 'DiffLogoBasics.Rnw'

###################################################
### code chunk number 1: ImportLibrary
###################################################
library(DiffLogo)


###################################################
### code chunk number 2: ImportMotifsFromMotifDb
###################################################
library(MotifDb)

## import motifs
hitIndeces = grep ('CTCF', values (MotifDb)$geneSymbol, ignore.case=TRUE)
list    = as.list(MotifDb[hitIndeces])
sequenceCounts = as.numeric(values (MotifDb)$sequenceCount[hitIndeces])
names(sequenceCounts) = names(list)

pwm1 = reverseComplement(list$"Hsapiens-JASPAR_CORE-CTCF-MA0139.1"[, 2:18]) ## trim and reverse complement
pwm2 = list$"Hsapiens-jolma2013-CTCF"
n1 = sequenceCounts["Hsapiens-JASPAR_CORE-CTCF-MA0139.1"]
n2 = sequenceCounts["Hsapiens-jolma2013-CTCF"]

## DiffLogo can also handle motifs of different length
pwm_long = reverseComplement(list$"Hsapiens-JASPAR_CORE-CTCF-MA0139.1") ## reverse complement
pwm_short = list$"Hsapiens-jolma2013-CTCF"


###################################################
### code chunk number 3: ImportMotifsFromDiffLogoPackage
###################################################
## import five DNA motifs from matrix
motif_folder = "extdata/pwm"
motif_names_dna = c("H1-hESC", "MCF7", "HeLa-S3", "HepG2", "HUVEC")
motifs_dna = list()
for (name in motif_names_dna) {
  fileName = paste(motif_folder,"/",name,".pwm",sep="")
  file = system.file(fileName, package = "DiffLogo")
  motifs_dna[[name]] = getPwmFromPwmFile(file)
}
sampleSizes_dna = c("H1-hESC"=100, "MCF7"=100, "HeLa-S3"=100, "HepG2"=100, "HUVEC"=100)
## import three DNA motifs from table
motif_folder = "extdata/alignments"
motif_names_dna2 = c("Mad", "Max", "Myc")
motifs_dna2 = list()
for (name in motif_names_dna2) {
  fileName = paste(motif_folder,"/",name,".txt",sep="")
  file = system.file(fileName, package = "DiffLogo")
  motifs_dna2[[name]] = getPwmFromAlignmentFile(file)
}
## import three ASN motifs from fasta files
motif_folder = "extdata/alignments"
motif_names_asn = c("F-box_fungi.seq", "F-box_metazoa.seq", "F-box_viridiplantae.seq")
motifs_asn = list()
for (name in motif_names_asn) {
  fileName = paste(motif_folder,"/",name,".fa",sep="")
  file = system.file(fileName, package = "DiffLogo")
  motifs_asn[[name]] = getPwmFromFastaFile(file, FULL_ALPHABET)
}


###################################################
### code chunk number 4: PlotSequenceLogo
###################################################
## plot classic sequence logo
library(seqLogo)
seqLogo::seqLogo(pwm = pwm1)


###################################################
### code chunk number 5: PlotCustomSequenceLogo
###################################################
## plot custom sequence logo
par(mfrow=c(2,1), pin=c(3, 1), mar = c(2, 4, 1, 1))
DiffLogo::seqLogo(pwm = pwm1)
DiffLogo::seqLogo(pwm = pwm2, stackHeight = sumProbabilities)
par(mfrow=c(1,1), pin=c(1, 1), mar=c(5.1, 4.1, 4.1, 2.1))


###################################################
### code chunk number 6: PlotDiffLogo
###################################################
## plot DiffLogo
diffLogoFromPwm(pwm1 = pwm1, pwm2 = pwm2)

## diffLogoFromPwm is a convenience function for
diffLogoObj = createDiffLogoObject(pwm1 = pwm1, pwm2 = pwm2)
diffLogo(diffLogoObj)

## mark symbol stacks with significant changes
diffLogoObj = enrichDiffLogoObjectWithPvalues(diffLogoObj, n1, n2)
diffLogo(diffLogoObj)

## plot DiffLogo for PWMs of different length
diffLogoFromPwm(pwm1 = pwm_long, pwm2 = pwm_short, align_pwms=TRUE)


###################################################
### code chunk number 7: PlotDiffLogoTable
###################################################
## plot table of difference logos for CTFC motifs (DNA)
diffLogoTable(PWMs = motifs_dna)
## plot table of difference logos for E-Box motifs (DNA)
diffLogoTable(PWMs = motifs_dna2)
## plot table of difference logos for F-Box motifs (ASN)
diffLogoTable(PWMs = motifs_asn, alphabet = FULL_ALPHABET)
## CTFC motifs (DNA) with asterisks to indicate significant differences
diffLogoTable(PWMs = motifs_dna, sampleSizes = sampleSizes_dna)


###################################################
### code chunk number 8: Export
###################################################
## parameters
widthToHeightRatio = 16/10;
size = length(motifs_dna) * 2
resolution = 300
width = size * widthToHeightRatio
height = size

## export single DiffLogo as pdf document
fileName = "Comparison_of_two_motifs.pdf"
pdf(file = fileName, width = width, height = height)
diffLogoFromPwm(pwm1 = pwm1, pwm2 = pwm2)
dev.off()

## export DiffLogo table as png image
fileName = "Comparison_of_multiple_motifs.png"
png(
  filename = fileName, res = resolution,
  width = width * resolution, height = height * resolution)
diffLogoTable(PWMs = motifs_dna)
dev.off()


###################################################
### code chunk number 9: Alignment
###################################################
## align pwms and than plot seqLogo of them
pwms = list(pwm_long, pwm_short)
multiple_pwms_alignment = multipleLocalPwmsAlignment(pwms)
aligned_pwms = extendPwmsFromAlignmentVector(pwms, multiple_pwms_alignment$alignment$vector)

## seqLogo of aligned pwms
layout(mat = matrix(data = 1:4, nrow = 2, ncol = 2))
DiffLogo::seqLogo(pwms[[1]], main = "Unaligned")
DiffLogo::seqLogo(pwms[[2]])
DiffLogo::seqLogo(aligned_pwms[[1]], main = "Aligned")
DiffLogo::seqLogo(aligned_pwms[[2]])


