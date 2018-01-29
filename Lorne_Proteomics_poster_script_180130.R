######################################
##         Madeleine J Otway        ##
##        motway@cmri.org.au        ##
##            30/01/2018            ##
##      Lorne Proteomics 2018       ##
##      Source code for poster      ##
##     Data available on GitHub     ##
##   R version 3.4.3 (2017-11-30)   ##
##     Bioconductor version 3.6     ##
######################################

#install.source()
library(dialects)

#install.packages("fBasics")
library(fBasics)

#install.packages("ggplot2")
library(ggplot2)

#source("http://bioconductor.org/biocLite.R")
#biocLite("plyr")
library(plyr)

##

## Rat SRL convert to Human
rat.srl <- import.srl("/Users/MadeleineOtway/Documents/PhD/R/dialects/Data/SRL/MJO_PGH_Rat_SRL_dialects.txt")
human.fasta <- import.fasta("/Users/MadeleineOtway/Documents/PhD/R/dialects/Data/Fasta/uniprot-human.fasta")
digest.human <- digest.fasta(human.fasta)
human.from.rat <- convert.species(rat.srl, digest.human)
export.srl(human.from.rat, "/Users/MadeleineOtway/Documents/PhD/R/dialects/Data/SRL/MJO_PGH_Human_from_Rat_SRL_dialects.txt")


##

## Import results from human SRL
pep.human <- read.table("/Users/MadeleineOtway/Documents/PhD/R/dialects/Data/PeakExtractions/Human/MJO_PGH_Human_PeakExtraction_dialects_areas - Peptides.txt",
                        sep = "\t", header = T)
fdr.human <- read.table("/Users/MadeleineOtway/Documents/PhD/R/dialects/Data/PeakExtractions/Human/MJO_PGH_Human_PeakExtraction_dialects_FDR.txt",
                        sep = "\t", header = T)
## Remove decoy and match FDR and peptide data frames
fdr.decoy <- fdr.human[!(fdr.human$Decoy=="True"),]
fdr.decoy$fdr.pep <- ifelse(fdr.decoy$Protein %in% pep.human$Protein &
                              fdr.decoy$Peptide %in% pep.human$Peptide &
                              fdr.decoy$Precursor.MZ %in% pep.human$Precursor.MZ,
                            1, 0)
fdr.match <- fdr.decoy[!(fdr.decoy$fdr.pep == "0"),]
## Remove all peptides with an FDR greater that 1%
sort.fdr <- arrange(fdr.match, Protein, Peptide, Precursor.MZ)
sort.pep <- arrange(pep.human, Protein, Peptide, Precursor.MZ)
sort.fdr$Label <- NULL
sort.fdr$Decoy <- NULL
human.01 <- sort.pep
for (i in 6:10) {
  human.01[,i] <- ifelse(sort.fdr[,i] <= 0.01, human.01[,i], NA)
}
## Remove Missing values
human.01 <- human.01[(rowSums(is.na(human.01[6:10])) <= 0) ,]

##
## Repeat for human from rat SRL results
pep.human.r <- read.table("/Users/MadeleineOtway/Documents/PhD/R/dialects/Data/PeakExtractions/Human/MJO_PGH_Human_from_Rat_PeakExtraction_dialects_areas - Peptides.txt",
                          sep = "\t", header = T)
fdr.human.r <- read.table("/Users/MadeleineOtway/Documents/PhD/R/dialects/Data/PeakExtractions/Human/MJO_PGH_Human_from_Rat_PeakExtraction_dialects_FDR.txt",
                          sep = "\t", header = T)
fdr.decoy.h <- fdr.human.r[!(fdr.human.r$Decoy=="True"),]
fdr.decoy.h$fdr.pep <- ifelse(fdr.decoy.h$Protein %in% pep.human.r$Protein &
                                fdr.decoy.h$Peptide %in% pep.human.r$Peptide &
                                fdr.decoy.h$Precursor.MZ %in% pep.human.r$Precursor.MZ,
                              1, 0)
fdr.match.h <- fdr.decoy.h[!(fdr.decoy.h$fdr.pep == "0"),]
sort.fdr.h <- arrange(fdr.match.h, Protein, Peptide, Precursor.MZ)
sort.pep.h <- arrange(pep.human.r, Protein, Peptide, Precursor.MZ)
sort.fdr.h$Label <- NULL
sort.fdr.h$Decoy <- NULL
human.r.01 <- sort.pep.h
for (i in 6:10) {
  human.r.01[,i] <- ifelse(sort.fdr.h[,i] <= 0.01, human.r.01[,i], NA)
}
human.r.01 <- human.r.01[(rowSums(is.na(human.r.01[6:10])) <= 0) ,]


human.01$human.match <- ifelse(human.01$Protein %in% human.r.01$Protein &
                                 human.01$Peptide %in% human.r.01$Peptide &
                                 human.01$Precursor.MZ %in% human.r.01$Precursor.MZ &
                                 human.01$Precursor.Charge %in% human.r.01$Precursor.Charge,
                               1, 0)
human.match <- human.01[!(human.01$human.match == "0"),]

human.r.01$human.match <- ifelse(human.r.01$Protein %in% human.01$Protein &
                                   human.r.01$Peptide %in% human.01$Peptide &
                                   human.r.01$Precursor.MZ %in% human.01$Precursor.MZ &
                                   human.r.01$Precursor.Charge %in% human.01$Precursor.Charge,
                                 1, 0)

human.r.match <- human.r.01[!(human.r.01$human.match == "0"),]


human.r.match$mean <- rowMeans(human.r.match[4:10])
human.r.match$sd <- rowStdevs(human.r.match[4:10])
human.match$mean <- rowMeans(human.match[4:10])
human.match$sd <- rowStdevs(human.match[4:10])

rat.srl <- import.srl("/Users/MadeleineOtway/Documents/PhD/R/dialects/Data/SRL/MJO_PGH_Rat_SRL_dialects.txt")
rat.srl$rat.match <- ifelse(rat.srl$stripped_sequence %in% human.match$Peptide &
                              rat.srl$prec_z %in% human.match$Precursor.Charge,
                            1, 0)
rat.srl.match <- rat.srl[!(rat.srl$rat.match == "0"),]
rat.srl.match$new.u <- rat.match$uniprot_id
rat.srl.match$new.u <- gsub("..\\|.*\\|", "", rat.match$new.u)
rat.srl.match$new.u <- gsub("_RAT", "", rat.match$new.u)

human.match$new.u <- human.match$Protein
human.match$new.u <- gsub("..\\|.*\\|", "", human.match$new.u)
human.match$new.u <- gsub("_HUMAN", "", human.match$new.u)
human.match$u.match <- ifelse(human.match$new.u %in% rat.srl.match$new.u, 
                              1, 0)
human.u.match <- human.match[!(human.match$u.match == "1"),]
human.u.match <- human.u.match[!duplicated(human.u.match$new.u),]


human.r.match$new.u <- human.r.match$Protein
human.r.match$new.u <- gsub("..\\|.*\\|", "", human.r.match$new.u)
human.r.match$new.u <- gsub("_HUMAN", "", human.r.match$new.u)
human.r.match$u.match <- ifelse(human.r.match$new.u %in% rat.srl.match$new.u, 
                                1, 0)
human.r.u.match <- human.r.match[!(human.r.match$u.match == "1"),]
human.r.u.match <- human.r.u.match[!duplicated(human.r.u.match$new.u),]


h <- data.frame("Peptide" = human.u.match[1:5,2])
h$Intensity <- human.u.match[1:5,12]
h$group <- "human"

hr <- data.frame("Peptide" = human.r.u.match[1:5,2])
hr$Intensity <- human.r.u.match[1:5,12]
hr$group <- "human.r"
h.hr <- rbind(h, hr)

ggplot(h.hr, aes(x = Peptide, y = Intensity, fill = group)) +
  geom_bar(stat = "identity", color = "black", position=position_dodge()) +
  labs(title = NULL, x = NULL, y = NULL) +
  theme_classic() +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 20),
        axis.text.y = element_text(hjust = 1, size = 20)) +
  guides(fill = F) +
  scale_fill_manual(values=c('#008000','#FF0080'))

##

