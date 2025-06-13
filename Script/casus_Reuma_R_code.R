#het humane genoom indexeren.
buildindex(
  basename = 'ref_human',
  reference = 'Homo_sapiens.GRCh38.dna.toplevel',
  memory = 8000,
  indexSplit = TRUE)

# Mapping van de test patienten, tegen het geindexde humane genoom
align.SRR4785819 <- align(index = "ref_human", readfile1 = "data_RA_raw/Data_RA_raw/SRR4785819_1_subset40k.fastq", readfile2 = "data_RA_raw/Data_RA_raw/SRR4785819_2_subset40k.fastq", output_file = "SRR4785819.BAM")
align.SRR4785820 <- align(index = "ref_human", readfile1 = "data_RA_raw/Data_RA_raw/SRR4785820_1_subset40k.fastq", readfile2 = "data_RA_raw/Data_RA_raw/SRR4785820_2_subset40k.fastq", output_file = "SRR4785820.BAM")
align.SRR4785828 <- align(index = "ref_human", readfile1 = "data_RA_raw/Data_RA_raw/SRR4785828_1_subset40k.fastq", readfile2 = "data_RA_raw/Data_RA_raw/SRR4785828_2_subset40k.fastq", output_file = "SRR4785828.BAM")
align.SRR4785831 <- align(index = "ref_human", readfile1 = "data_RA_raw/Data_RA_raw/SRR4785831_1_subset40k.fastq", readfile2 = "data_RA_raw/Data_RA_raw/SRR4785831_2_subset40k.fastq", output_file = "SRR4785831.BAM")
align.SRR4785979 <- align(index = "ref_human", readfile1 = "data_RA_raw/Data_RA_raw/SRR4785979_1_subset40k.fastq", readfile2 = "data_RA_raw/Data_RA_raw/SRR4785979_2_subset40k.fastq", output_file = "SRR4785979.BAM")
align.SRR4785980 <- align(index = "ref_human", readfile1 = "data_RA_raw/Data_RA_raw/SRR4785980_1_subset40k.fastq", readfile2 = "data_RA_raw/Data_RA_raw/SRR4785980_2_subset40k.fastq", output_file = "SRR4785980.BAM")
align.SRR4785986 <- align(index = "ref_human", readfile1 = "data_RA_raw/Data_RA_raw/SRR4785986_1_subset40k.fastq", readfile2 = "data_RA_raw/Data_RA_raw/SRR4785986_2_subset40k.fastq", output_file = "SRR4785986.BAM")

# Laad Rsamtools voor sorteren en indexeren
library(Rsamtools)

# Bestandsnamen van de monsters
samples <- c('SRR4785819', 'SRR4785820', 'SRR4785828', 'SRR4785831', 'SRR4785979', 'SRR4785980', 'SRR4785986', 'SRR4785988')

# Voor elk monster: sorteer en indexeer de BAM-file
# Sorteer BAM-bestanden
lapply(samples, function(s) {sortBam(file = paste0(s, '.BAM'), destination = paste0(s, '.sorted'))
})
#indexeer de BAM-bestanden
lapply(samples, function(s) {indexBam(file = paste0(s, '.sorted.bam'))})
  
#dag 2: count matrix
library(readr)
library(dplyr)
library(Rsamtools)
library(Rsubread)
# Je definieert een vector met namen van BAM-bestanden. Elke BAM bevat reads van een RNA-seq-experiment (bijv. behandeld vs. controle).

allsamples <- c("SRR4785819.BAM", "SRR4785820.BAM", "SRR4785828.BAM", "SRR4785831.BAM", "SRR4785979.BAM", "SRR4785980.BAM", "SRR4785986.BAM", "SRR4785988.BAM")

#feature counts uitvoeren, De BAM-bestanden vergelijken met het gedownloade "homo_sapiens_GRCH38.114.gtf.gz" bestand. hierop wordt baseerd hoeveel reads op elke gen valt.
count_matrix <- featureCounts(
  files = allsamples,
  annot.ext = "Homo_sapiens.GRCh38.114.gtf.gz",
  isPairedEnd = TRUE,
  isGTFAnnotationFile = TRUE,
  GTF.attrType = "gene_id",
  useMetaFeatures = TRUE
)
#uitgebreide count_matrix inladen, bevat meer genomen.
count_matrix <- read.table(file = "count_matrix.txt", header = TRUE)

#.csv bestand maken van het count matrix
write_csv2(count_matrix, "countmatrix.csv")

#dag 3 statistiek en analyse
#behandelings tabel maken
Diagnose <- c("Normal", "Normal", "Normal", "Normal", "RA", "RA", "RA", "RA")
Diagnose_table <- data.frame(Diagnose)
rownames(Diagnose_table) <- c('SRR4785819', 'SRR4785820', 'SRR4785828', 'SRR4785831', 'SRR4785979', 'SRR4785980', 'SRR4785986', 'SRR4785988')

#benodigde packages downloaden/inladen.
BiocManager::install("DESeq2")
BiocManager::install("KEGGREST")
library(DESeq2)
library(KEGGREST)

# Maak DESeqDataSet aan. met DESeqDataSetFromMatrix zetten we de gegevens binnen de data set "count_matrix" in de juiste format voor de statistiche analyse
dds <- DESeqDataSetFromMatrix(countData = round(count_matrix),
                              colData = Diagnose_table,
                              design = ~ Diagnose,)

# Voer voer de statistische analyse uit
dds <- DESeq(dds)
#Isoleer de resultaten
resultaten <- results(dds)

# Resultaten opslaan in een bestand
#Bij het opslaan van je tabel kan je opnieuw je pad instellen met `setwd()` of het gehele pad waar je de tabel wilt opslaan opgeven in de code.

write.table(resultaten, file = 'ResultatenDSS.csv', row.names = TRUE, col.names = TRUE)

#kijken hoeveel reads op of neer gereguleerd zijn
sum(resultaten$padj < 0.05 & resultaten$log2FoldChange > 1, na.rm = TRUE)
sum(resultaten$padj < 0.05 & resultaten$log2FoldChange < -1, na.rm = TRUE)

#Sorteren van de meest significante resultaten. gebaseerd op hoogste/laagste fold change en laagste p-waarde.
hoogste_fold_change <- resultaten[order(resultaten$log2FoldChange, decreasing = TRUE), ]
laagste_fold_change <- resultaten[order(resultaten$log2FoldChange, decreasing = FALSE), ]
laagste_p_waarde <- resultaten[order(resultaten$padj, decreasing = FALSE), ]

#GO analyse voorbereiden
#Rij namen van alle resultaten isoleren
all <- rownames(resultaten)

#filtreren op significante p waarde en een log2foldchange van boven de 2 en onder de -2.

resultaten <- as.data.frame(resultaten)
resultaten_filtered <- resultaten %>%
  filter(padj < 0.05, log2FoldChange < -2 | log2FoldChange > 2)

#filtreren op rownames
DEG <- rownames(resultaten_filtered)


#GO analyse uitvoeren
#benodigde packages instaleren/inladen
BiocManager::install("goseq")
BiocManager::install("geneLenDataBase")
BiocManager""install("org.Hs.eg.db")
library("goseq")
library("geneLenDataBase")
library("org.Hs.eg.db")


gene.vector=as.integer(all%in%DEG)
names(gene.vector)=all 
head(gene.vector)
tail(gene.vector)

pwf=nullp(gene.vector,"hg19","geneSymbol")

#find enriched Go terms
GO.wall=goseq(pwf,"hg19","geneSymbol")

#How many enriched GO terms do we have
class(GO.wall)
head(GO.wall)
nrow(GO.wall)

#sorteren op significantie
enriched.GO=GO.wall$category[GO.wall$over_represented_pvalue<.05]

#check hoeveel enriched.GO we hebben
class(enriched.GO)
head(enriched.GO)
length(enriched.GO)

#downloaden/inladen van GO.db
BiocManager::install("GO.db")
library(GO.db)

#Go resultaten exporten naar een nieuwe bestand genaamd "SigGo.txt"
capture.output(for(go in enriched.GO[1:258]) { print(GOTERM[[go]])
  cat("--------------------------------------\n")
}
, file="SigGo.txt")

#top 10 visualiseren
library("ggplot2")
top10_GO <- GO.wall %>%
  dplyr::filter(category %in% enriched.GO) %>%
  arrange(over_represented_pvalue) %>%
  slice_head(n = 10)

top10_GO$term <- sapply(top10_GO$category, function(go_id) {
       term <- Term(GOTERM[[go_id]])
       if (is.null(term)) return(NA)
       else return(term)
  })

ggplot(top10_GO, aes(x = reorder(term, -log10(over_represented_pvalue)), 
                     y = -log10(over_represented_pvalue))) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Top 10 Enriched GO Terms",
       x = "GO Term",
       y = "-log10(P-value)") +
  theme_minimal()

#volcano plot uitvoeren en visualiseren.

if (!requireNamespace("EnhancedVolcano", quietly = TRUE)) {
  BiocManager::install("EnhancedVolcano")
}
library(EnhancedVolcano)

EnhancedVolcano(resultaten,
                lab = rownames(resultaten),
                x = 'log2FoldChange',
                y = 'padj')

# Alternatieve plot zonder p-waarde cutoff (alle genen zichtbaar)
EnhancedVolcano(resultaten,
                lab = rownames(resultaten),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0)

# het figuur opslaan
dev.copy(png, 'VolcanoplotWC.png', 
         width = 8,
         height = 10,
         units = 'in',
         res = 500)
dev.off()

#KEGG analyse uitvoeren
if (!requireNamespace("pathview", quietly = TRUE)) {
  BiocManager::install("pathview")
}
library(pathview)

#binnen "resultaten" de foldchange isoleren.
resultaten[1] <- NULL
resultaten[2:5] <- NULL

#pathway visualiseren waar RA het meest impact op heeft. hierbij is de kegg pathway "hsa04672" gekozen.
pathview(
  gene.data = resultaten,
  pathway.id = "hsa04672",  # KEGG pathway ID
  species = "hsa",          # 'hsa' = homo sapiens in KEGG
  gene.idtype = "SYMBOL",   
  limit = list(gene = 5)    # Kleurbereik voor log2FC van -5 tot +5
)
