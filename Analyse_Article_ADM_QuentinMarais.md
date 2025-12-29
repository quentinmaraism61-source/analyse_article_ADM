Analyse_Article_ADM_QuentinMarais
================

# Installation du package DADA2

``` r
library(dada2)
```

    ## Loading required package: Rcpp

``` r
library(Rcpp)
```

# Chargement du jeu de (méta)données

``` r
 runs <- read.csv("SraRunTable (3).csv", header = TRUE, stringsAsFactors = FALSE)

head(runs)
```

    ##           Run Assay.Type AvgSpotLen     Bases  BioProject    BioSample
    ## 1 SRR25410637        WGS        551 276561504 PRJNA962314 SAMN34394469
    ## 2 SRR25410638        WGS        599 151988990 PRJNA962314 SAMN34394469
    ## 3 SRR25410639        WGS        587 344283431 PRJNA962314 SAMN34394472
    ## 4 SRR25410640        WGS        599 148616214 PRJNA962314 SAMN34394472
    ## 5 SRR25410641        WGS        588    144219 PRJNA962314 SAMN34394473
    ## 6 SRR25410642        WGS        599 155963390 PRJNA962314 SAMN34394473
    ##                BioSampleModel     Bytes              Center.Name
    ## 1 Metagenome or environmental 182230791 ASTROBIOLOGY CENTER INTA
    ## 2 Metagenome or environmental 100785655 ASTROBIOLOGY CENTER INTA
    ## 3 Metagenome or environmental 226164938 ASTROBIOLOGY CENTER INTA
    ## 4 Metagenome or environmental  97397792 ASTROBIOLOGY CENTER INTA
    ## 5 Metagenome or environmental    170776 ASTROBIOLOGY CENTER INTA
    ## 6 Metagenome or environmental 101723058 ASTROBIOLOGY CENTER INTA
    ##   Collection_Date Consent DATASTORE.filetype DATASTORE.provider
    ## 1         2019-02  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 2         2019-02  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 3         2019-02  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 4         2019-02  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 5         2019-02  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 6         2019-02  public   fastq,run.zq,sra         gs,ncbi,s3
    ##                       DATASTORE.region  Experiment geo_loc_name_country
    ## 1 gs.us-east1,ncbi.public,s3.us-east-1 SRX21146351         uncalculated
    ## 2 gs.us-east1,ncbi.public,s3.us-east-1 SRX21146350         uncalculated
    ## 3 gs.us-east1,ncbi.public,s3.us-east-1 SRX21146349         uncalculated
    ## 4 gs.us-east1,ncbi.public,s3.us-east-1 SRX21146348         uncalculated
    ## 5 gs.us-east1,ncbi.public,s3.us-east-1 SRX21146347         uncalculated
    ## 6 gs.us-east1,ncbi.public,s3.us-east-1 SRX21146346         uncalculated
    ##   geo_loc_name_country_continent geo_loc_name     Instrument
    ## 1                   uncalculated   Antarctica Illumina MiSeq
    ## 2                   uncalculated   Antarctica Illumina MiSeq
    ## 3                   uncalculated   Antarctica Illumina MiSeq
    ## 4                   uncalculated   Antarctica Illumina MiSeq
    ## 5                   uncalculated   Antarctica Illumina MiSeq
    ## 6                   uncalculated   Antarctica Illumina MiSeq
    ##           isolation_source                     lat_lon  Library.Name
    ## 1 Antarctic volcanic rocks         62.9302 S 60.6853 W 18S-DIVOL-12A
    ## 2 Antarctic volcanic rocks         62.9302 S 60.6853 W 16S-DIVOL-12A
    ## 3 Antarctic volcanic rocks 62.97740199 S 60.68060199 W    18S-NB1-S2
    ## 4 Antarctic volcanic rocks 62.97740199 S 60.68060199 W    16S-NB1-S2
    ## 5 Antarctic volcanic rocks 62.98130699 S 60.70235000 W 18S-DIPV36-S3
    ## 6 Antarctic volcanic rocks 62.98130699 S 60.70235000 W 16S-DIPV36-S3
    ##   LibraryLayout LibrarySelection LibrarySource        Organism Platform
    ## 1        SINGLE      unspecified   METAGENOMIC rock metagenome ILLUMINA
    ## 2        SINGLE      unspecified   METAGENOMIC rock metagenome ILLUMINA
    ## 3        SINGLE      unspecified   METAGENOMIC rock metagenome ILLUMINA
    ## 4        SINGLE      unspecified   METAGENOMIC rock metagenome ILLUMINA
    ## 5        SINGLE      unspecified   METAGENOMIC rock metagenome ILLUMINA
    ## 6        SINGLE      unspecified   METAGENOMIC rock metagenome ILLUMINA
    ##            ReleaseDate          create_date version Sample.Name SRA.Study
    ## 1 2025-09-18T00:00:00Z 2023-07-25T03:50:00Z       1   DIVOL-12A SRP451238
    ## 2 2025-09-18T00:00:00Z 2023-07-25T03:51:00Z       1   DIVOL-12A SRP451238
    ## 3 2025-09-18T00:00:00Z 2023-07-25T03:51:00Z       1      NB1-S2 SRP451238
    ## 4 2025-09-18T00:00:00Z 2023-07-25T03:51:00Z       1      NB1-S2 SRP451238
    ## 5 2025-09-18T00:00:00Z 2023-07-25T03:50:00Z       1   DIPV36-S3 SRP451238
    ## 6 2025-09-18T00:00:00Z 2023-07-25T03:51:00Z       1   DIPV36-S3 SRP451238

``` r
getwd()
```

    ## [1] "/home/rstudio/analyse_article_ADM"

``` r
path <- "/home/rstudio/analyse_article_ADM/data/" # On donne le chemin d'accès vers les fichiers fastq 
list.files(path = "/home/rstudio/analyse_article_ADM/data/", pattern = "\\.fastq\\.gz$")
```

    ##  [1] "SRR25410637_1.fastq.gz" "SRR25410637_2.fastq.gz" "SRR25410638_1.fastq.gz"
    ##  [4] "SRR25410638_2.fastq.gz" "SRR25410639_1.fastq.gz" "SRR25410639_2.fastq.gz"
    ##  [7] "SRR25410640_1.fastq.gz" "SRR25410640_2.fastq.gz" "SRR25410641_1.fastq.gz"
    ## [10] "SRR25410641_2.fastq.gz" "SRR25410642_1.fastq.gz" "SRR25410642_2.fastq.gz"
    ## [13] "SRR25410643_1.fastq.gz" "SRR25410643_2.fastq.gz" "SRR25410644_1.fastq.gz"
    ## [16] "SRR25410644_2.fastq.gz" "SRR25410645_1.fastq.gz" "SRR25410645_2.fastq.gz"
    ## [19] "SRR25410646_1.fastq.gz" "SRR25410646_2.fastq.gz"

``` r
list.files(path) # On s'assure que nos fichiers fastq apparaissent bien dans path. 
```

    ##  [1] "filtered"               "SRR25410637_1.fastq.gz" "SRR25410637_2.fastq.gz"
    ##  [4] "SRR25410638_1.fastq.gz" "SRR25410638_2.fastq.gz" "SRR25410639_1.fastq.gz"
    ##  [7] "SRR25410639_2.fastq.gz" "SRR25410640_1.fastq.gz" "SRR25410640_2.fastq.gz"
    ## [10] "SRR25410641_1.fastq.gz" "SRR25410641_2.fastq.gz" "SRR25410642_1.fastq.gz"
    ## [13] "SRR25410642_2.fastq.gz" "SRR25410643_1.fastq.gz" "SRR25410643_2.fastq.gz"
    ## [16] "SRR25410644_1.fastq.gz" "SRR25410644_2.fastq.gz" "SRR25410645_1.fastq.gz"
    ## [19] "SRR25410645_2.fastq.gz" "SRR25410646_1.fastq.gz" "SRR25410646_2.fastq.gz"

Les données ont correctement été importées, on peut à présent commencer
l’analyse à proprement parler.

# Création des listes pour les fichiers Forward (fnFs) et Reverse (fnRs)

``` r
fnFs <- sort(list.files(path, pattern="_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq.gz", full.names = TRUE))
```

# Inspection des profils de qualité des reads

## Qualité de la lecture Forward

``` r
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1) # extraction du nom des échantillons

plotQualityProfile(fnFs[1:2])
```

![](Analyse_Article_ADM_QuentinMarais_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

COMMENTAIRE SUR LA QUALITE DES READS

## Qualtié de la lecture Reverse

``` r
plotQualityProfile(fnRs[1:2])
```

![](Analyse_Article_ADM_QuentinMarais_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->
COMMENTAIRE QUALITE LECTURE REVERSE

# Avan-Filtration

``` r
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
```

On attribut les noms des échantillons aux noms des fichiers filtrés.

# Filtration et pré-traitement des reads

``` r
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE (only needed for filterAndTrim)
head(out)
```

    ##                        reads.in reads.out
    ## SRR25410637_1.fastq.gz   501870    384402
    ## SRR25410638_1.fastq.gz   253452    201720
    ## SRR25410639_1.fastq.gz   585662    487287
    ## SRR25410640_1.fastq.gz   247778    198357
    ## SRR25410641_1.fastq.gz      245        20
    ## SRR25410642_1.fastq.gz   260028    209676

Le tableau ci-dessus nous montre le nombre de reads avant (reads.in) et
après (reads.out) filtration.En premier lieu, on peut remarquer que le
nombre de reads par échantillon est hétérogène (allant de 245 à plus de
500 000).

Ensuite concernant la filtration, on observe que dans la majorité des
cas, une grande parties des reads passent le filtre bien que la quantité
apparentre de reads ne passant le filtre puisse être non négligeable. En
effet, il faudra prendre en compte des biais comme le fait que les
échantillons avec peu de reads peuvent apparaître comme pauvres en
diversité car la profondeure de séquençage est insuffisante. Les
échantillons très riches en reads vont dominer certaines analyses
(ordination, clustering…). Il faudra garder en tête l’hétérogénéité du
nombre de reads et la perte d’informations (reads) dans la suite des
analyses (problèmes potentiels dans la normalisation, distances
écologiques, ordination, taxa rares…).

# Taux d’erreur de séquençage

Il arrive régulièrement que des reads contiennent des erreurs de
séquençage. La pipeline DADA2 permet de discriminer les séquences réels
des séquences contenant des erreurs de séquençage.

``` r
errF <- learnErrors(filtFs, multithread=TRUE)
```

    ## 140669280 total bases in 586122 reads from 2 samples will be used for learning the error rates.

``` r
errR <- learnErrors(filtRs, multithread=TRUE)
```

    ## 171745440 total bases in 1073409 reads from 3 samples will be used for learning the error rates.

## Visualisation du taux d’erreurs de séquençage

``` r
plotErrors(errF, nominalQ=TRUE)
```

    ## Warning: Transformation introduced infinite values in continuous y-axis

![](Analyse_Article_ADM_QuentinMarais_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

# Inférence des échantillons : Identification des ASV (amplicon sequence variant) présentes dans chaque échantillon

## Algorithme de DADA2 sur les reads Forward

``` r
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
```

    ## Sample 1 - 384402 reads in 63774 unique sequences.
    ## Sample 2 - 201720 reads in 31065 unique sequences.
    ## Sample 3 - 487287 reads in 110253 unique sequences.
    ## Sample 4 - 198357 reads in 45001 unique sequences.
    ## Sample 5 - 20 reads in 20 unique sequences.
    ## Sample 6 - 209676 reads in 33348 unique sequences.
    ## Sample 7 - 192969 reads in 36202 unique sequences.
    ## Sample 8 - 131986 reads in 31295 unique sequences.
    ## Sample 9 - 202333 reads in 52926 unique sequences.
    ## Sample 10 - 97632 reads in 27386 unique sequences.

## Algorithme de DADA2 sur les reads Reverse

``` r
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
```

    ## Sample 1 - 384402 reads in 43763 unique sequences.
    ## Sample 2 - 201720 reads in 23876 unique sequences.
    ## Sample 3 - 487287 reads in 67315 unique sequences.
    ## Sample 4 - 198357 reads in 28947 unique sequences.
    ## Sample 5 - 20 reads in 19 unique sequences.
    ## Sample 6 - 209676 reads in 22924 unique sequences.
    ## Sample 7 - 192969 reads in 27128 unique sequences.
    ## Sample 8 - 131986 reads in 21175 unique sequences.
    ## Sample 9 - 202333 reads in 38588 unique sequences.
    ## Sample 10 - 97632 reads in 19264 unique sequences.

## Résultats : Visualisation d’un exemple

``` r
dadaFs[[1]]
```

    ## dada-class: object describing DADA2 denoising results
    ## 937 sequence variants were inferred from 63774 input unique sequences.
    ## Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

Pour le premier échantillon, l’algorithme a identifié 937 variants
biologiques au sein 63 774 séquences.

COMMENTAIRE

# Fusion des reads appariés

``` r
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
```

    ## 37 paired-reads (in 1 unique pairings) successfully merged out of 381533 (in 18323 pairings) input.

    ## 93 paired-reads (in 42 unique pairings) successfully merged out of 200350 (in 3840 pairings) input.

    ## 38 paired-reads (in 2 unique pairings) successfully merged out of 480722 (in 41776 pairings) input.

    ## 31 paired-reads (in 15 unique pairings) successfully merged out of 194185 (in 12668 pairings) input.

    ## No paired-reads (in ZERO unique pairings) successfully merged out of 20 pairings) input.

    ## 50 paired-reads (in 21 unique pairings) successfully merged out of 206930 (in 6843 pairings) input.

    ## 6 paired-reads (in 2 unique pairings) successfully merged out of 191668 (in 4786 pairings) input.

    ## 563 paired-reads (in 54 unique pairings) successfully merged out of 127598 (in 10696 pairings) input.

    ## 5 paired-reads (in 3 unique pairings) successfully merged out of 198860 (in 16300 pairings) input.

    ## 0 paired-reads (in 0 unique pairings) successfully merged out of 93570 (in 8506 pairings) input.

``` r
head(mergers[[1]])
```

    ##                                                                                                                                                                                                                                                                                sequence
    ## 1003 GCCAGCAGCTGCGCGAACATGACAATCGGTCGGGCGCTCCGCCTGCTGCTGACGCTCACCGGCGGCGCTCGTCCGGGAGGGCTCGATCGCTCCACCCTCGGGCATCCCGGCAAGCTCGCCACGTGCTTCGCGGAGAACGAGGAGGCGAGCCCGTGGGAGCCGCTTCACGTCGAGCGCGGATTCGGCCCGGAGACCTCCACGGTCACCCTCGTCGCCGGGGACGCGCCGTTATCGATCTCCGACCATCGCAGACGCACTCCGAAGAAATTGACGG
    ##      abundance forward reverse nmatch nmismatch nindel prefer accept
    ## 1003        37     316     252    126         0      0      2   TRUE

# Construction d’une table contenant les séquences

``` r
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
```

    ## [1]  10 134

COMMENTAIRE

# Distruibution des longueurs des reads

``` r
table(nchar(getSequences(seqtab)))
```

    ## 
    ## 271 273 274 287 291 293 295 298 299 311 312 328 333 343 362 
    ##   1   1   7  54   2   2   7  42   9   2   2   1   1   2   1

COMMENTAIRE

# Suppression des chimères

``` r
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
```

    ## Identified 63 bimeras out of 134 input sequences.

``` r
dim(seqtab.nochim)
```

    ## [1] 10 71

# Proportion de reads ayant passé la filtration des chimères

``` r
sum(seqtab.nochim)/sum(seqtab)
```

    ## [1] 0.671932

COMMENTAIRE

# Suivi des reads au sein de la pipeline

``` r
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
```

    ##              input filtered denoisedF denoisedR merged nonchim
    ## SRR25410637 501870   384402    382730    383125     37      37
    ## SRR25410638 253452   201720    201020    200993     93      66
    ## SRR25410639 585662   487287    482684    485198     38      38
    ## SRR25410640 247778   198357    195746    196633     31      30
    ## SRR25410641    245       20         1         3      0       0
    ## SRR25410642 260028   209676    208117    208406     50      38

COMMENTAIRE / Grande perte d’information !

# Assignation taxonomique

COMMENTAIRE

ATTENTION : APRES ASSIGNATION, faut dégager les eucaryotes (séquences
mal assignées ou NA).
