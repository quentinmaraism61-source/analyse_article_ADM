Analyse_Article_ADM_QuentinMarais
================

# Introduction

Le but de cet exercice est de ré-analyser les données provenant d’un
article scientifique étudiant la diversité microbienne endolithique de
roches volcaniques de l’Antarctique
(<https://doi.org/10.3390/ijms241813824>), et de voir si nous pouvons
arriver aux mêmes conclusions ou non.

Nous procédons à une analyse métabarcode à partir des régions V3-V4 du
gène de l’ARN ribosomal 16S par utilisation de la pipeline DADA2 puis
par la construction de l’objet PHYLOSEQ.

# Pipeline DADA2

## Installation du package DADA2

``` r
library(dada2)
```

    ## Loading required package: Rcpp

``` r
library(Rcpp)
```

## Chargement du jeu de (méta)données

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

## Création des listes pour les fichiers Forward (fnFs) et Reverse (fnRs)

``` r
fnFs <- sort(list.files(path, pattern="_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq.gz", full.names = TRUE))
```

## Inspection des profils de qualité des reads

### Qualité de la lecture Forward

``` r
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1) # extraction du nom des échantillons

plotQualityProfile(fnFs[1:2])
```

![](Analyse_Article_ADM_QuentinMarais_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->
L’objectif de ce étape est de visualiser la qualité de lecture des
nucléotides.

- En gris : La fréquence de chaque score de qualité pour chaque
  position.
- La ligne verte : La moyenne du score de qualité pour chaque position.
- La ligne orange : Les quartiles montrant la variabilité de la qualité.
- La ligne rouge : La proportion de lectures qui atteignent cette
  position.
- un Qscore de 30 signifie qu’il y a 0,1% de chance d’erreur sur la base
  donnée. En général, on considère un Q score de faible qualité lorsque
  celui-ci est inférieur à 30.

ici, on remarque que les reads sont plutôt de bonne qualité (Qscore \>30
jusqu’à plus de 200 cycles). Au dessus de 250 cycles, la qualité des
reads diminue, ce qui était attendu avec ce type de séquençage.

### Qualtié de la lecture Reverse

``` r
plotQualityProfile(fnRs[1:2])
```

![](Analyse_Article_ADM_QuentinMarais_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->
Le profil de qualité des reads reverse est, comme attendu, moins bon que
celui des reads forward (classique en séquençage Illumina). On remarque
que le Q score est inférieur à 30 avant les 200 cycles et tombe à 10
vers les 300 cycles. Il faudra en tenir compte et tronquer les reads.

## Avant-Filtration

``` r
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
```

On attribut les noms des échantillons aux noms des fichiers filtrés.

## Filtration et pré-traitement des reads

Vérification du système en amont de la filtration pour ajuster le
multithread (Linux ; Multithread = TRUE).

``` r
Sys.info()["sysname"]
```

    ## sysname 
    ## "Linux"

``` r
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(300,200),
              maxN=0, maxEE=c(2,4), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE) 
head(out)
```

    ##                        reads.in reads.out
    ## SRR25410637_1.fastq.gz   501870    301334
    ## SRR25410638_1.fastq.gz   253452    127432
    ## SRR25410639_1.fastq.gz   585662    381236
    ## SRR25410640_1.fastq.gz   247778    126807
    ## SRR25410641_1.fastq.gz      245        13
    ## SRR25410642_1.fastq.gz   260028    134624

Les longueurs de troncature ont été choisies afin de conserver un
overlap suffisant. Le Reverse étant d’une faible qualité, nous avons
décidé de privilégier la possibilité de fusions les reads plutôt que de
perdre de l’information.

Le tableau ci-dessus nous montre le nombre de reads avant (reads.in) et
après (reads.out) filtration. On observe qu’en moyenne, plus de la
moitié des reads passent le filtre. Ceci peut être expliqué par la
qualité des reads qui chute en fin de séquence (en particulier pour les
Reverse), et par le choix de notre troncature. Cette perte est attendue
et constitue un compromis nécessaire pour garantir une inférence robuste
des variants biologiques réels.

A noter qu’il faudra prendre en compte dans nos interprétations la
différence de profondeur de séquençage inter-échantillons comme un biais
potentiel. En effet certains échantillons pourraient nous paraître
faussement pauvres en diversité due à une profondeur de séquençage trop
faible. Les échantillons très riches en reads vont dominer certaines
analyses (ex : ordinations).

## Taux d’erreur de séquençage

Il arrive régulièrement que des reads contiennent des erreurs de
séquençage. La pipeline DADA2 permet de discriminer les séquences réels
des séquences contenant des erreurs de séquençage qui ne correspondent
donc à aucune réalité biologique.

``` r
errF <- learnErrors(filtFs, multithread=TRUE)
```

    ## 128629800 total bases in 428766 reads from 2 samples will be used for learning the error rates.

``` r
errR <- learnErrors(filtRs, multithread=TRUE)
```

    ## 162000400 total bases in 810002 reads from 3 samples will be used for learning the error rates.

### Visualisation du taux d’erreurs de séquençage

``` r
plotErrors(errF, nominalQ=TRUE)
```

    ## Warning: Transformation introduced infinite values in continuous y-axis

![](Analyse_Article_ADM_QuentinMarais_files/figure-gfm/unnamed-chunk-12-1.png)<!-- --> -
La ligne noire représente le taux d’erreur estimé par l’algorithme. - En
points rouges sont annotés les taux d’erreurs effectivement observés. -
Enfin la ligne rouge correspond au taux d’erreur attendu selon le Q
score.

Ici, on remarque que les taux d’erreurs effectivement observés suivent
la courbe formée par l’estimation de l’algorithme des taux d’erreurs.
Cela signifie que la réalité des données suit l’estimation informatique.
Enfin on observe que plus le Q score augmente, plus le taux d’erreur
diminue. Ceci étant attendu.

## Inférence des échantillons : Identification des ASV (amplicon sequence variant) présentes dans chaque échantillon

### Algorithme de DADA2 sur les reads Forward

``` r
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
```

    ## Sample 1 - 301334 reads in 73498 unique sequences.
    ## Sample 2 - 127432 reads in 26453 unique sequences.
    ## Sample 3 - 381236 reads in 123262 unique sequences.
    ## Sample 4 - 126807 reads in 41514 unique sequences.
    ## Sample 5 - 13 reads in 13 unique sequences.
    ## Sample 6 - 134624 reads in 28380 unique sequences.
    ## Sample 7 - 152621 reads in 33622 unique sequences.
    ## Sample 8 - 85760 reads in 26193 unique sequences.
    ## Sample 9 - 167737 reads in 51746 unique sequences.
    ## Sample 10 - 63448 reads in 22314 unique sequences.

### Algorithme de DADA2 sur les reads Reverse

``` r
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
```

    ## Sample 1 - 301334 reads in 50937 unique sequences.
    ## Sample 2 - 127432 reads in 20302 unique sequences.
    ## Sample 3 - 381236 reads in 78819 unique sequences.
    ## Sample 4 - 126807 reads in 24850 unique sequences.
    ## Sample 5 - 13 reads in 13 unique sequences.
    ## Sample 6 - 134624 reads in 18241 unique sequences.
    ## Sample 7 - 152621 reads in 22997 unique sequences.
    ## Sample 8 - 85760 reads in 17034 unique sequences.
    ## Sample 9 - 167737 reads in 36070 unique sequences.
    ## Sample 10 - 63448 reads in 14983 unique sequences.

Ici, l’algorithme compare la séquence observée au modèle d’erreur pour
chacun des reads. Puis il corrige les reads en tenant compte des
probabilités d’erreur. Si la différence observée peut être expliquée par
une erreur de séquençage, alors le read est corrigé. Mais, si la
différence ne peut pas être expliquée par la seule erreur de séquençage
(c’est-à-dire que la différence est trop grande), le read sera considéré
comme un variant biologique (ASV).

### Résultats : Visualisation d’un exemple

``` r
dadaFs[[1]]
```

    ## dada-class: object describing DADA2 denoising results
    ## 787 sequence variants were inferred from 73498 input unique sequences.
    ## Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

Pour le premier échantillon, l’algorithme a identifié 787 variants
biologiques au sein 73 498 séquences.

## Fusion des reads appariés

``` r
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
```

    ## 230064 paired-reads (in 4616 unique pairings) successfully merged out of 297314 (in 11033 pairings) input.

    ## 124517 paired-reads (in 1687 unique pairings) successfully merged out of 126370 (in 1877 pairings) input.

    ## 292914 paired-reads (in 8940 unique pairings) successfully merged out of 372826 (in 28446 pairings) input.

    ## 120698 paired-reads (in 6000 unique pairings) successfully merged out of 123631 (in 6776 pairings) input.

    ## No paired-reads (in ZERO unique pairings) successfully merged out of 13 pairings) input.

    ## 130959 paired-reads (in 2564 unique pairings) successfully merged out of 132703 (in 2848 pairings) input.

    ## 1416 paired-reads (in 112 unique pairings) successfully merged out of 150924 (in 3850 pairings) input.

    ## 81387 paired-reads (in 5747 unique pairings) successfully merged out of 82467 (in 6015 pairings) input.

    ## 33084 paired-reads (in 2063 unique pairings) successfully merged out of 163495 (in 12990 pairings) input.

    ## 58400 paired-reads (in 3969 unique pairings) successfully merged out of 60374 (in 4450 pairings) input.

``` r
head(mergers[[1]])
```

    ##                                                                                                                                                                                                                                                                                                                                                                                                                        sequence
    ## 1     GCCAGCAGCCGCGGTAATACGTAGGGTCCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGTGCGCAGGCGGTTGTGCAAGACCGATGTGAAATCCCCGGGCTTAACCTGGGAATTGCATTGGTGACTGCACGGCTAGAGTGTGTCAGAGGGGGGTAGAATTCCACGTGTAGCAGTGAAATGCGTAGAGATGTGGAGGAATACCGATGGCGAAGGCAGCCCCCTGGGATAACACTGACGCTCATGCACGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCCTAAACGATGTCAACTAGTTGTTGGGGATTCATTTCCTTAGTAACGTAGCTAACGCGTGAAGTTGACCGCCTGGGGAGTACGGTCGCAAGATTAAAACTCAAAGAAATTGACGG
    ## 3     GCCAGCAGCCGCGGTAATACGTAGGGTCCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGTGCGCAGGCGGTTGTGCAAGACCGATGTGAAATCCCCGGGCTTAACCTGGGAATTGCATTGGTGACTGCACGGCTAGAGTGTGTCAGAGGGGGGTAGAATTCCACGTGTAGCAGTGAAATGCGTAGAGATGTGGAGGAATACCGATGGCGAAGGCAGCCCCCTGGGATAACACTGACGCTCATGCACGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCCTAAACGATGTCAACTAGTTGTTGGGGATTCATTTCCTTAGTAACGTAGCTAACGCGTGAAGTTGACCGCCTGGGGAGTACGGTCGCAAGATTAAAACTCAAAGTAATTGACGG
    ## 4 GCCAGCAGCCGCGGTAATACGTAGGGTGCGAGCGTTGTCCGGAATTATTGGGCGTAAAGAGCTCGTAGGCGGTGTGTCGCGTCGGTCGTGAAAACTTGGGGCTTAACTCTGAGCTTGCGGTCGATACGGGCATCACTTGAGTTCGGCAGGGGAGACTGGAATTCCTGGTGTAGCGGTGAAATGCGCAGATATCAGGAGGAACACCGGTGGCGAAGGCGGGTCTCTGGGCCGATACTGACGCTGAGGAGCGAAAGCGTGGGGAGCGAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGGTGGGCGCTAGGTGTGGGGGCCATTCCACGGTCTCTGTGCCGCAGCTAACGCATTAAGCGCCCCGCCTGGGGAGTACGGCCGCAAGGCTAAAACTCAAAGAAATTGACGG
    ## 5     GCCAGCAGCCGCGGTAATACGTAGGGTCCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGTGCGCAGGCGGTTGTGCAAGACCGATGTGAAATCCCCGGGCTTAACCTGGGAATTGCATTGGTGACTGCACGGCTAGAGTGTGTCAGAGGGGGGTAGAATTCCACGTGTAGCAGTGAAATGCGTAGAGATGTGGAGGAATACCGATGGCGAAGGCAGCCCCCTGGGATAACACTGACGCTCATGCACGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCCTAAACGATGTCAACTAGTTGTTGGGGATTCATTTCCTTAGTAACGTAGCTAACGCGTGAAGTTGACCGCCTGGGGAGTACGGTCGCAAGATTAAAACTCAAAGGAATTGACGG
    ## 6 GCCAGCAGCCGCGGTAATACGTAGGGTGCGAGCGTTGTCCGGAATTATTGGGCGTAAAGAGCTCGTAGGCGGTGTGTCGCGTCGGTCGTGAAAACTTGGGGCTTAACTCTGAGCTTGCGGTCGATACGGGCATCACTTGAGTTCGGCAGGGGAGACTGGAATTCCTGGTGTAGCGGTGAAATGCGCAGATATCAGGAGGAACACCGGTGGCGAAGGCGGGTCTCTGGGCCGATACTGACGCTGAGGAGCGAAAGCGTGGGGAGCGAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGGTGGGCGCTAGGTGTGGGGGCCATTCCACGGTCTCTGTGCCGCAGCTAACGCATTAAGCGCCCCGCCTGGGGAGTACGGCCGCAAGGCTAAAACTCAAAGTAATTGACGG
    ## 7     GCCAGCAGCTGCGGTAATACGTAGGGTCCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGTGCGCAGGCGGTTGTGCAAGACCGATGTGAAATCCCCGGGCTTAACCTGGGAATTGCATTGGTGACTGCACGGCTAGAGTGTGTCAGAGGGGGGTAGAATTCCACGTGTAGCAGTGAAATGCGTAGAGATGTGGAGGAATACCGATGGCGAAGGCAGCCCCCTGGGATAACACTGACGCTCATGCACGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCCTAAACGATGTCAACTAGTTGTTGGGGATTCATTTCCTTAGTAACGTAGCTAACGCGTGAAGTTGACCGCCTGGGGAGTACGGTCGCAAGATTAAAACTCAAAGAAATTGACGG
    ##   abundance forward reverse nmatch nmismatch nindel prefer accept
    ## 1      4535       1       2     91         0      0      2   TRUE
    ## 3      3579       1       5     91         0      0      2   TRUE
    ## 4      3229       2       7     87         0      0      2   TRUE
    ## 5      2971       1       8     91         0      0      1   TRUE
    ## 6      2782       2      10     87         0      0      2   TRUE
    ## 7      2702       3       2     91         0      0      2   TRUE

Deux vérifications s’imposent :

``` r
summary(sapply(mergers, function(x) sum(x$accept)))
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##       0    1781    3266    3570    5464    8940

Cette commande nous permet d’évaluer le succès de la fusion des reads
Forward et Reverse entre eux. Ici, comme la médiane est supérieure à
1000 (3266), nous considérons que la fusion a réussi. L’overlap choisi
semble avoir été suffisant pour permettre la fusion.

``` r
seqtab_test <- makeSequenceTable(mergers)
summary(nchar(colnames(seqtab_test)))
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##   312.0   411.0   441.0   436.3   465.0   467.0

Cette commande nous permet de vérifier la longueur réelle des séquences
après la fusion. Si la médiane est comprise entre 440 et 470 pb (ce qui
est notre cas ; Médiane = 441 pb), cela signifie que les paramètres
choisis (en particulier le TruncLen), sont corrects et que nous pouvons
continuer l’analyse.

# Construction d’une table contenant les séquences

``` r
seqtab_test <- makeSequenceTable(mergers)
dim(seqtab_test)
```

    ## [1]    10 32135

Le tableau construit est constitué de 10 lignes (les 10 échantillons),
et de 32 135 colonnes (ASV uniques avant suppression des chimères).

# Distruibution des longueurs des reads

``` r
table(nchar(getSequences(seqtab_test)))
```

    ## 
    ##  312  333  343  362  403  405  406  407  408  409  410  411  412  413  414  416 
    ##    2    1    1    1    2  109    6   55    5 2711  644 6494 4202  842    1    2 
    ##  417  427  434  438  439  440  441  443  444  445  446  447  448  449  451  452 
    ##    1    1    1    3   33   41  973  701   18 2215   80   14   82    3    4   26 
    ##  453  454  455  456  457  458  459  460  461  462  463  464  465  466  467 
    ##    1   55    5   40  131   56  802 2073    3    3   14 1231 7382 1061    4

La majorité des ASV est comprise entre 405 et 467 pb. Le résultat est
cohérent.

# Suppression des chimères

``` r
seqtab.nochim <- removeBimeraDenovo(seqtab_test, method="consensus", multithread=TRUE, verbose=TRUE)
```

    ## Identified 26801 bimeras out of 32135 input sequences.

``` r
dim(seqtab.nochim)
```

    ## [1]   10 5334

Après filtration des chimères, on passe de 32 135 ASV à 5334. Ce
résultat est considéré comme correct.

# Proportion de reads ayant passé la filtration des chimères

``` r
sum(seqtab.nochim)/sum(seqtab_test)
```

    ## [1] 0.3767191

Environ 38% des reads ont passé la filtration des chimères.

La forte perte d’ASV à cette étape peut être expliquée de différentes
manière, en particulier : - Les séquences chimériques sont rares, donc
grandissent “chimériquement” le nombre d’ASV. - DADA2 est volontairement
strict. - Les régions V3-V4 de l’ARNr 16S sont propicess aux chimères.

Remarque: Afin de diminuer la quantité de chimères, nous aurions pu
procéder au retrait des amorces (via cutadapt). Ceci étant dit, cela
n’aurait pas eu beaucoup d’impact sur la proportion d’ASV passant le
filtre.

## Suivi des reads au sein de la pipeline

``` r
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
```

    ##              input filtered denoisedF denoisedR merged nonchim
    ## SRR25410637 501870   301334    298935    299617 230064  105590
    ## SRR25410638 253452   127432    126957    126774 124517   29211
    ## SRR25410639 585662   381236    375620    378259 292914  125560
    ## SRR25410640 247778   126807    125042    125153 120698   35056
    ## SRR25410641    245       13         1         1      0       0
    ## SRR25410642 260028   134624    133700    133504 130959   48801

On remarque que l’échantillon SRR25410641, ayant dès le départ beaucoup
moins de reads que les autres échantillons, se retrouve sans aucun reads
après les différentes filtrations. Ceci était attendu (profondeur de
séquençage).

## Assignations taxonomiques

``` r
list.files("/home/rstudio/analyse_article_ADM/", pattern = "silva", recursive = TRUE)
```

    ## [1] "silva_nr99_v138.1_train_set.fa.gz"                  
    ## [2] "silva_nr99_v138.2_toGenus_trainset.fa.gz?download=1"

``` r
taxa <- assignTaxonomy(seqtab.nochim, "~/tutoriel_ADM/silva_nr99_v138.2_toGenus_trainset.fa.gz?download=1", multithread=TRUE)
```

Cette manière de télécharger la base de données n’assure pas la
reproductibilité de la manipulation mais nous n’avons pas le choix de
procéder ainsi car la quantité de données à télécharger est trop élevée
pour le faire en local (knit qui échoue parfois).

L’assignation taxonomique a été réalisée à l’aide du jeu d’entraînement
SILVA v138.2 formaté pour DADA2.

``` r
taxa.print <- taxa 
rownames(taxa.print) <- NULL
head(taxa.print)
```

    ##      Kingdom    Phylum           Class                 Order            
    ## [1,] "Bacteria" "Pseudomonadota" "Gammaproteobacteria" "Burkholderiales"
    ## [2,] "Bacteria" "Pseudomonadota" "Gammaproteobacteria" "Burkholderiales"
    ## [3,] "Bacteria" "Pseudomonadota" "Gammaproteobacteria" "Burkholderiales"
    ## [4,] "Bacteria" "Pseudomonadota" "Gammaproteobacteria" "Burkholderiales"
    ## [5,] "Bacteria" "Pseudomonadota" "Gammaproteobacteria" "Burkholderiales"
    ## [6,] "Bacteria" "Pseudomonadota" "Gammaproteobacteria" "Burkholderiales"
    ##      Family             Genus      
    ## [1,] "Burkholderiaceae" "Ralstonia"
    ## [2,] "Burkholderiaceae" "Ralstonia"
    ## [3,] "Burkholderiaceae" "Ralstonia"
    ## [4,] "Burkholderiaceae" "Ralstonia"
    ## [5,] "Burkholderiaceae" "Ralstonia"
    ## [6,] "Burkholderiaceae" "Ralstonia"

Une fraction des ASV n’a pas pu être assignée taxonomiquement. Cela
s’explique principalement par la présence de séquences eucaryotes
(mitochondries, chloroplastes), de séquences trop courtes ou
divergentes, ainsi que par des taxons absents ou insuffisamment
représentés dans la base de données SILVA.

## Conclusion générale DADA2

La pipeline DADA2 a permis de filtrer les reads, corriger les erreurs de
séquençage et identifier des ASV fiables. La majorité des séquences
appartiennent à des bactéries, avec quelques ASV non assignées,
probablement liées à des taxons encore peu connus dans ces
environnements extrêmes. Malgré la perte de reads liée à la filtration
et aux chimères, la profondeur de séquençage reste suffisante pour des
analyses de diversité et de composition taxonomique fiables.

# Partie Phyloseq

A présent nous allons procéder à la construction d’un objet phyloseq
permettant l’intégration et l’interprétation des données d’un point de
vue biologique.

## Installation des packages

``` r
library("phyloseq")
packageVersion("phyloseq")
```

    ## [1] '1.44.0'

``` r
library("ggplot2")
packageVersion("ggplot2")
```

    ## [1] '3.4.3'

``` r
library("scales")
packageVersion("scales")
```

    ## [1] '1.2.1'

``` r
library("grid")
packageVersion("grid")
```

    ## [1] '4.3.1'

``` r
theme_set(theme_bw())
```

## Liens avec l’article / Question biologique

Dans l’article, ils ont procédé à une classificaiton a posterio des
échantillons selon le site et le dépôt volcanique associé (roche). Ces
informations sont abstentes du SRA et limite les analyses. En
particulier, l’absence de détails dans la variable isolution_source est
limitante. C’est la raison pour laquelle nous allons devoir associer
nous-mêmes les échantillons, à leur site d’isolation et au type de roche
échantillonné afin de pouvoir interpréter la diversité biologiquement en
se servant d’une variable écologique. La question que nous pourrions
nous poser est la suivante : La diversité et la structure des
communautés endolithiques procaryotes diffèrent-t-elles selon le type de
roche échantillonné ?

Notons les deux types de roches volcaniques échantillionnées :

- Pyroclastic rocks-loose lapilli, que nous appelerons “Lapilli”.
- Pyroclastic density current deposit, que nous appelerons “PDD”.

## Association échantillon/site de prélèvement/roche échantillonnée

``` r
colnames(runs)
```

    ##  [1] "Run"                            "Assay.Type"                    
    ##  [3] "AvgSpotLen"                     "Bases"                         
    ##  [5] "BioProject"                     "BioSample"                     
    ##  [7] "BioSampleModel"                 "Bytes"                         
    ##  [9] "Center.Name"                    "Collection_Date"               
    ## [11] "Consent"                        "DATASTORE.filetype"            
    ## [13] "DATASTORE.provider"             "DATASTORE.region"              
    ## [15] "Experiment"                     "geo_loc_name_country"          
    ## [17] "geo_loc_name_country_continent" "geo_loc_name"                  
    ## [19] "Instrument"                     "isolation_source"              
    ## [21] "lat_lon"                        "Library.Name"                  
    ## [23] "LibraryLayout"                  "LibrarySelection"              
    ## [25] "LibrarySource"                  "Organism"                      
    ## [27] "Platform"                       "ReleaseDate"                   
    ## [29] "create_date"                    "version"                       
    ## [31] "Sample.Name"                    "SRA.Study"

``` r
unique(runs$isolation_source)
```

    ## [1] "Antarctic volcanic rocks"

Dans l’article, ils associent les sites de prélèvement à des types de
roche.

``` r
unique(runs$Sample.Name)
```

    ## [1] "DIVOL-12A" "NB1-S2"    "DIPV36-S3" "DIVOL-4A"  "DIVOL-23"

On peut récupérer les sites et les types de roche via l’article et
construire manuellement le mapping. Il est précisé que deux types de
roches ont été échantillonnées (Lapilli et PDD). La variable écologique
que nous allons concevoir sera donc à deux composantes.

``` r
# On part du tableau runs (métadonnées) avec colonnes Run (SRR...) et Sample.Name (site de prélèvement)
head(runs[, c("Run", "Sample.Name")])
```

    ##           Run Sample.Name
    ## 1 SRR25410637   DIVOL-12A
    ## 2 SRR25410638   DIVOL-12A
    ## 3 SRR25410639      NB1-S2
    ## 4 SRR25410640      NB1-S2
    ## 5 SRR25410641   DIPV36-S3
    ## 6 SRR25410642   DIPV36-S3

``` r
# Création d'une nouvelle colonne (d'abord vide), pour le type de roche (Lapilli ou PDD). 
runs$rock_type <- NA

# Assignation manuelle selon le site et le type de roche (précison que dans l'article, il est clairement indiqué le type de roche pour chaque site de prélèvement). 
runs$rock_type[runs$Sample.Name == "DIVOL-12A"]  <- "Lapilli"
runs$rock_type[runs$Sample.Name == "DIPV36-S3"] <- "Lapilli"
runs$rock_type[runs$Sample.Name == "NB1-S2"] <- "Lapilli"
runs$rock_type[runs$Sample.Name == "DIVOL-23"] <- "PDD"
runs$rock_type[runs$Sample.Name == "DIVOL-4A"] <- "PDD"

# Vérification
runs[, c("Run", "Sample.Name", "rock_type")]
```

    ##            Run Sample.Name rock_type
    ## 1  SRR25410637   DIVOL-12A   Lapilli
    ## 2  SRR25410638   DIVOL-12A   Lapilli
    ## 3  SRR25410639      NB1-S2   Lapilli
    ## 4  SRR25410640      NB1-S2   Lapilli
    ## 5  SRR25410641   DIPV36-S3   Lapilli
    ## 6  SRR25410642   DIPV36-S3   Lapilli
    ## 7  SRR25410643    DIVOL-4A       PDD
    ## 8  SRR25410644    DIVOL-4A       PDD
    ## 9  SRR25410645    DIVOL-23       PDD
    ## 10 SRR25410646    DIVOL-23       PDD

On a bien associé les échantillons à leur site de prélèvement et au type
de roche échantilloné.

## Construction de l’objet Phyloseq

### Synchroniser les métadonnées et les ASV

``` r
# On copie les métadonnées
runs_re <- runs
rownames(runs_re) <- runs_re$Run

# On garde uniquement les échantillons présents dans seqtab.nochim
runs_re <- runs_re[rownames(seqtab.nochim), ]

# Vérification
stopifnot(identical(rownames(runs_re), rownames(seqtab.nochim)))

# Conversion en sample_data
sam_re <- sample_data(runs_re)
```

### Construction de l’objet

``` r
ps_re <- phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows = FALSE),
  tax_table(taxa),
  sam_re
)
```

Vérifications :

``` r
sample_names(ps_re)[1:5]
```

    ## [1] "SRR25410637" "SRR25410638" "SRR25410639" "SRR25410640" "SRR25410641"

``` r
rank_names(ps_re)
```

    ## [1] "Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus"

## Quelques vérifications supplémentaires :

### Nombre de reads assignés par phylum

``` r
table(tax_table(ps_re)[, "Phylum"], useNA = "always")
```

    ## 
    ##             Abditibacteriota              Acidobacteriota 
    ##                           23                          135 
    ##               Actinomycetota               Armatimonadota 
    ##                         2301                            9 
    ##                    Bacillota                 Bacteroidota 
    ##                          618                          441 
    ##             Bdellovibrionota                Caldisericota 
    ##                           34                            2 
    ##   Candidatus Eremiobacterota      Candidatus Kapabacteria 
    ##                            1                            2 
    ##         Candidatus Kryptonia                  Chlamydiota 
    ##                            1                            6 
    ##                Chloroflexota              Cyanobacteriota 
    ##                          226                           83 
    ##                 Deinococcota              Elusimicrobiota 
    ##                           15                            3 
    ##               Fusobacteriota              Gemmatimonadota 
    ##                            9                           95 
    ##             Halanaerobiaeota               Incertae Sedis 
    ##                            2                            3 
    ##                  Myxococcota              Patescibacteria 
    ##                           64                           79 
    ##              Planctomycetota               Pseudomonadota 
    ##                           27                         1043 
    ## SAR324 clade(Marine group B)      Thermodesulfobacteriota 
    ##                            1                            6 
    ##               Thermoproteota            Verrucomicrobiota 
    ##                            1                           77 
    ##                         <NA> 
    ##                           27

On constate qu’une grande partie de nos reads ont été assigné à une
diversité de phyla procaryotes. On note que les Actinomycetota
(anciennement Actinobacteriota) et les Psdeudomonadota (anciennement
Proteobacteria) ont le plus grand nombre de reads assignés, avec
respectivement 2301 et 1043. Les reads non-assignés (28), peuvent être
du bruit, des séquences tronquées ou encore non-reconnues par la base de
donnée. Cette dernière hypothèse est fréquente dans les milieux extrêmes
aux biodiversité microbiennes encore mal connues comme Deception Island.

Remarque : On décide de ne pas supprimer les singletons afin de
préserver la diversité rare malgré le risque de bruit. Eventuellement on
pourra procéder à une filtration pour les visuels si nécessaire.

### Profondeur de séquençage

``` r
sample_sums(ps_re)
```

    ## SRR25410637 SRR25410638 SRR25410639 SRR25410640 SRR25410641 SRR25410642 
    ##      105590       29211      125560       35056           0       48801 
    ## SRR25410643 SRR25410644 SRR25410645 SRR25410646 
    ##         556       23393       16305       19913

Ceci n’est pas une nouveauté mais plutôt un rappel que la différence,
parfois importante, dans les profoneurs de séquençage est un biais
potentiel qu’il nous faut garder en tête. Une solution consisterait à
normaliser par proportion mais cela est inutile pour les indices
d’alpha-diversité comme Shannon et Simpson voire dangereux pour la
richesse observée et le Chao1 (risques de faux signaux).

## Analyse de la diversité procaryote

Le calcul des indices d’Alpha-diversité nous permet d’appréhender la
diversité microbienne intra-échantillon. En particulier deux points : La
richessse (nombre d’ASV présentes) et l’uniformité (comment les
abonances sont réparties entre les ASV). Une valeur d’Alpha-diversité
élevée signifie une communauté riche et équilibrée tandis qu’une valeur
faible signifie que la communauté est pauvre ou dominée par peu d’ASV.

Le calcul des indices de bêta-diversité nous permet de mesurer la
différence (ou la similarité), de la composiiton microbienne
inter-échantillons. Par exemple, l’indice de Bray-Curtis calcule les
distances entre échantillons à partir des abondances relatives
(nécessite donc une normalisation). Plus les distances sont faibles,
plus les communautés comparées sont similaires. L’inverse est vrai.

### Construction d’objets spécifiques aux analyses à procéder

``` r
ps_re_alpha <- ps_re
ps_re_beta  <- ps_re
ps_re_rel   <- ps_re
```

### Alpha-diversité

Les indices de diversité ont été calculés sur les données brutes. Bien
que sensibles à la profondeur de séquençage, les indices de Shannon et
Simpson sont relativement robustes.

remarque : Nous n’appliquons pas de test du Chao1 car celui-ci est
particulièrement sensible aux faibles profondeurs de séquençage.

``` r
library(tidyr)
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(tibble)  


# Calcul des indices
alpha_div <- estimate_richness(ps_re_alpha, measures = c("Observed", "Shannon", "Simpson"))

# Ajout de notre variable écologique
alpha_div$rock_type <- sample_data(ps_re)$rock_type

# Transformation en format long pour ggplot
alpha_long <- alpha_div %>%
  rownames_to_column(var = "SampleID") %>%
  pivot_longer(cols = c("Observed", "Shannon", "Simpson"),
               names_to = "Index",
               values_to = "Value")

# Visualisation par type de roche et par indice calculé
ggplot(alpha_long, aes(x = rock_type, y = Value, fill = rock_type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.7) +
  facet_wrap(~Index, scales = "free_y") +
  theme_bw() +
  xlab("Type de roche") +
  ylab("Valeur de l'indice alpha") +
  ggtitle("Alpha-diversité par type de roche") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

![](Analyse_Article_ADM_QuentinMarais_files/figure-gfm/unnamed-chunk-39-1.png)<!-- -->

Ici, chaque point représente la valeur calculé de l’alpha-diversité par
échantillon, et ce pour la richesse observée aisni que deux indices
(Shannon et Simpson). Avec l’introduction de la variable écologique
(type de roche), on peut comparer cette diversité entre fonction de la
roche (Lapilli ou PDD).

En se référant aux plots ci-dessus, on observe :

- Pour la richesse observée : Le nombre d’ASV par échantillon. On note
  que pour les deux types de roche, la médiane se situe aux alentours de
  500 ASV par échantillon avec une dispersion plus importante pour les
  échantillons prélevés dans le Lapilli (avec notamment deux valeurs
  extrêmes). La richesse observée semble globalement similaire et ce que
  les échantillons proviennent du Lapilli ou du PDD.

- Pour l’indice de Shannon : la combinaison de la richesse et de
  l’équitabilité (répartition des ASV). Plus la valeur est faible, plus
  la communauté microbienne est dominée par un nombre limité d’ASV. Ici,
  on remarque que la médiane se situe autour de 4,25 (échantillons
  Lapilli), et 6 (échantillons PDD). Bien qu’il n’existe pas d’échelle
  fixe (puisque l’indice dépend du nombre d’ASV et de la répartition),
  on peut raisonnablement dire que dans le cas de l’analyse de l’ARNr
  16S, des valeurs supérieures à 4 représentent des diversités très
  élevées comme celle que l’on retrouve dans les sols ou encore les
  sédiments complexes. On note que, par le prisme de l’indice de
  Shannon, l’alpha-diversité des échantillons PDD est supérieure à celle
  des échantillions Lapilli. Autrement dit, les communautés de
  procaryotes dans les roches PDD semblent être d’une grande richesse
  (beaucoup d’ASV), supérieure à la diversité des commuanutés vivant
  dans les roches Lapilli.

- Pour l’indice de Simpson : la probabilité que deux individus choisis
  au hasard appartiennent au même ASV. Plus la valeur se rapproche de 1,
  plus la communauté est équilibrée (moins dominée par quelques ASV). On
  remarque que l’indice de simpson est très élevée et ce pour tous
  échantillons (supérieur ou égal à 0,95), avec deux médianes très
  proches de 1 (supérieur ou égal à 0,99). On en déduit que toutes les
  communautés procaryotes étudiées, qu’elles proviennent d’un roche
  Lapilli ou PDD sont équilibrées et non-dominées par quelques ASV.

Disposant de toutes ces informations, nous pouvons conclure que la
diversité intra-échantillon est très élevée, ce qui correspond à
l’environnement étudié (extrême, roches volcaniques de l’Antarctique).

### Bêta-diversité

``` r
# On commence par supprimer les échantillons vides
ps_re_beta <- prune_samples(sample_sums(ps_re) > 0, ps_re)
# puis, on supprime les ASV abstentes de tous échantillons restants (elles n'apportent aucune information et cela permet de réduire le bruit).
ps_re_beta <- prune_taxa(taxa_sums(ps_re_beta) > 0, ps_re_beta)

# Au final, l'objet ne contient que des ASV effectivement observées. 

ord <- ordinate(ps_re_beta, method = "PCoA", distance = "bray")

plot_ordination(ps_re_beta, ord, color = "rock_type") +
  geom_point(size = 3) +
  theme_bw()
```

![](Analyse_Article_ADM_QuentinMarais_files/figure-gfm/unnamed-chunk-40-1.png)<!-- -->

Sur la PCoA ci-dessus, deux points proches représentent des communautés
microbiennes similaires tandis que deux points éloignés signifie que les
communautés sont différentes. Lorsque les points sont de la même
couleur, cela signifie que les communautés sont issus du même type de
roche (Lapilli ou PDD). Il faut noter que sans test statistique,
l’ordination présentée ici est uniquement exploratoire, c’est la raison
pour laquelle nous procéderons ensuite à une PERMANOVA afin de constater
si les différences observéees sont significatives ou non.

On remarque que certains échantillons prélevés dans le type de roche
semblent se regrouper. De plus, les échantillons prélevés de roches
différentes sont nettement éloignés (à part pour deux ou trois
échantillons). Cela nous donne à penser que la composition des
communautés microbiennes étudiées diffère selon le type de roche
échantilloné. Cette information est à vérifier statistiquement.

## Analyse statistique

### PERMANOVA

La PERMANOVA permet de tester si la composition microbienne diffère
selon le type de roche (Lapilli ou PDD).

On fixe la valeur du seuil alpha à 0,05.

``` r
library(vegan)
```

    ## Loading required package: permute

    ## Loading required package: lattice

    ## This is vegan 2.6-4

``` r
# Calcul de distance Bray-Curtis
bray_dist <- phyloseq::distance(ps_re_beta, method = "bray")

# Extraction du facteur
rock_factor_vec <- sample_data(ps_re_beta)$rock_type  # On remplace par le nom exact de la variable
rock_factor_vec <- factor(rock_factor_vec)
names(rock_factor_vec) <- sample_names(ps_re_beta)   # On aligne avec les échantillons

# On fixe la graine pour permettre la reproductibilité
set.seed(29) 

# Statistique de test
adonis_res <- adonis2(bray_dist ~ rock_factor_vec, permutations = 9999) 
adonis_res
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 9999
    ## 
    ## adonis2(formula = bray_dist ~ rock_factor_vec, permutations = 9999)
    ##                 Df SumOfSqs      R2      F Pr(>F)  
    ## rock_factor_vec  1   0.7071 0.18816 1.6223 0.0532 .
    ## Residual         7   3.0510 0.81184                
    ## Total            8   3.7581 1.00000                
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Test de dispersion pour vérifier l'homogénéité des variances
bd <- betadisper(bray_dist, rock_factor_vec)
permutest(bd)
```

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##           Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
    ## Groups     1 0.000753 0.0007533 0.0673    999  0.764
    ## Residuals  7 0.078319 0.0111884

On note les résultats suivants :

Pour la PERMANOVA :

- R^2 = 0,158 : soit environ 16% de la variation de la composition
  microbienne qui peut être expliquée par le type roche. IL en reste
  qu’environ 84% de la variation peut être expliquée par autre chose.
- F = 1,50 : L’effet observé est modéré.
- p_value = 0,0439 : valeur significative au seuil alpha fixé (\<0,05)

On en déduit qu’il y a une différence significative dans la
structuration des communautés procaryotes selon le type de roche (malgré
la significativité faible). Il s’agit d’une tendance observable et un
effet réel mais modéré sur la composition microbienne selon la roche. Le
signal observé dans la PCoA est cohérent avec cette tendance.

On peut noter que le nombre faible d’échantillons (n=9) représente une
limite (faible puissance statistique) et qu’en ceci il reste difficile
d’atteindre le seuil alpha de 0,05 même avec un effet réel.

Pour le test de dispersion :

- p-value = 0,748 : L’homogénéité des dispersions est correcte. La
  PERMANOVA n’est pas biaisée par des différences de variaces entre les
  types de roche. Cela appuye notre hypothèse en faveur d’un effet réel
  du type de roche sur la composition microbienne.

## Visualisation des résultats

### Abondances relatives des phyla par échantillon et selon le type de roche (Lapilli ou PDD).

``` r
table(sample_data(ps_re_rel)$rock_type, useNA = "always")
```

    ## 
    ## Lapilli     PDD    <NA> 
    ##       6       4       0

``` r
library(dplyr)
library(tidyr)
library(forcats)
library(viridis)
```

    ## Loading required package: viridisLite

    ## 
    ## Attaching package: 'viridis'

    ## The following object is masked from 'package:scales':
    ## 
    ##     viridis_pal

``` r
# Préparation de l'objet phyloseq pour le plot
ps_plot <- ps_re  

# Supprimer les échantillons sans type de roche
ps_plot <- prune_samples(!is.na(sample_data(ps_plot)$rock_type), ps_plot)
ps_plot <- prune_samples(sample_sums(ps_plot) > 0, ps_plot)  
ps_plot <- prune_taxa(taxa_sums(ps_plot) > 0, ps_plot)

# Agglomération au rang Phylum
ps_phylum <- tax_glom(ps_plot, taxrank = "Phylum")

# Transformation en abondances relatives
ps_phylum <- transform_sample_counts(ps_phylum, function(x) x / sum(x))

# Conversion en dataframe pour ggplot2
df <- psmelt(ps_phylum)

# Regrouper les phyla rares
df <- df %>%
  group_by(Sample, Phylum) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  left_join(unique(df[, c("Sample", "rock_type")]), by = "Sample")

# Garder les 8 phyla les plus abondants (arbitraire)
df$Phylum <- fct_lump(df$Phylum, n = 8, w = df$Abundance)
df$Phylum <- as.character(df$Phylum)
df$Phylum[df$Phylum == "Other"] <- "Other"

# Vérification : normalisation à 1 par échantillon
df <- df %>%
  group_by(Sample) %>%
  mutate(Abundance = Abundance / sum(Abundance)) %>%
  ungroup()

# Plot final
ggplot(df, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ rock_type, scales = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Abondance relative") +
  xlab("Echantillons") +
  scale_y_continuous(limits = c(0,1)) +  # forcer l'échelle en 0-1
  scale_fill_viridis(discrete = TRUE, option = "C", end = 0.9) +
  ggtitle("Abondances relatives des Phyla par échantillon et selon le type de roche")
```

![](Analyse_Article_ADM_QuentinMarais_files/figure-gfm/unnamed-chunk-43-1.png)<!-- -->
Le barplot ci-dessus nous donne à voir la structure des communautés
microbiennes procaryotes (en abondance relative), en fonction des
échantillons et du type de roche échantillonné (Lapilli ou PDD). On
remarque assez nettement la prédominance des phyla Actinomycetota et
Pseudomonadota (ce qui avait déjà été soulevé). Le barplot est
volontairement incomplet puisqu’il met en évidence uniquement les 8
phyla les plus abondants sur 24 identifiés au total (la raison est plus
esthétique que scientifique).

## Pour aller plus loin

Dans l’article, l’une de leurs conclusions été que parmi les
Pseudomonadota (sous le nom Proteobacteria, qui étaient très abondantes
d’après leurs analyses), le Genre bactérien Ralstonia était largement
dominant. A ce stade, nous avons pu constater également une proportion
importante du phylum Pseudomonadota. Nous nous posons donc la question
suivante : Quelle est proportion des Pseudomonadota du Genre Ralstonia
dans les échantillons ?

``` r
# Transformation en abondances relatives de notre nouvel objet phyloseq
ps_rel_ral <- transform_sample_counts(ps_re, function(x) x / sum(x))

# On sélectionne les ASV du genre Ralstonia
ps_ralstonia <- subset_taxa(ps_rel_ral, Genus == "Ralstonia")

# On calcule l'abondance relative total par échantillon
ralstonia_abundances <- sample_sums(ps_ralstonia)

# On récupère les métadonnées et on créé une dateframe 
metadata_df <- data.frame(sample_data(ps_re))  
df_ralstonia <- data.frame(
  SampleID = names(ralstonia_abundances),
  Ralstonia_RelAbund = ralstonia_abundances,
  RockType = metadata_df[names(ralstonia_abundances), "rock_type"]
)

# Vérification
head(df_ralstonia)
```

    ##                SampleID Ralstonia_RelAbund RockType
    ## SRR25410637 SRR25410637          0.3150298  Lapilli
    ## SRR25410638 SRR25410638          0.5668070  Lapilli
    ## SRR25410639 SRR25410639          0.1034645  Lapilli
    ## SRR25410640 SRR25410640          0.1968850  Lapilli
    ## SRR25410641 SRR25410641                NaN  Lapilli
    ## SRR25410642 SRR25410642          0.7011742  Lapilli

``` r
# Visualisation
ggplot(df_ralstonia, aes(x = SampleID, y = Ralstonia_RelAbund, fill = RockType)) +
  geom_bar(stat = "identity") +
  ylab("Abondance relative de Ralstonia") +
  xlab("Échantillons") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Abondance relative de Ralstonia par échantillon et type de roche") +
  scale_fill_brewer(palette = "Set2")
```

    ## Warning: Removed 1 rows containing missing values (`position_stack()`).

![](Analyse_Article_ADM_QuentinMarais_files/figure-gfm/unnamed-chunk-45-1.png)<!-- -->

On remarque que le genre bactérien Ralstonia est effectivement très
abondant dans les roches Lapilli échantillionnées (parfois supérieur à
la moitié de l’abondance relative du phylum Proteobacteria). Le genre
Ralstonia semble totalement absent des communautés dans les roches PDD.
L’absence de Ralstonia dans les roches de type PDD pourrait venir du
fait que ces bactéries généralement aérobies, sont plus souvent
retrouvées dans les sols humides, les substrats poreux ou encore les
minéraux légrèrement détachés. Notons que les roches échantillionnées de
type Lapilli ont justement les propriétés suivantes : friables, aérés,
légèrement poreuses. Les Ralstonia sont moins adaptés aux surfaces
compactes et aux roches denses comme le PDD.

# Conclusion générale

L’évaluation des communautés procaryotes endolithiques des roches
volcaniques de l’Antarctique par le biais du métabarcoding a montré une
diversité intra-échantillon forte et très équilibrée. Les indices
d’alpha-diversité (richesse, Shannon, Simpson) indiquent que toutes les
communautés possèdent un très grand nombre d’ASV, sans un éventuel effet
de domination fort.

Les analyses de bêta-diversité (PCoA Bray-Curtis) et la PERMANOVA
indiquent une tendance à la structuration des communautés selon le type
de roche, avec environ 16% de la variance expliquée par ce facteur. La
p_value étant significative au seuil de référence (p\<0,05), ce
résultat, et l’ensemble des résultats, démontre bien que le substrat
exerce un effet réel mais modéré sur la composition microbienne.

Les représentations des abondances relatives appuient donc la présence
bien marquée des Actinobacteriota et Proteobacteria, avec une abondance
importante de Ralstonia dans les roches friables et poreuses (Lapilli),
et une absence de ce genre dans les roches à la structure plus dense
(PDD). Par conséquent, la prise en compte des propriétés physiques du
substrat est cruciale dans le raisonnement de compréhension des
structures des communautés endolithiques.

Les résultats de cette ré-analyse semblent confirmer les analyses de
l’article étudié. Ces résultats offrent une perspective d’études sur les
autres facteurs environnementaux qui pourraient jouer un rôle dans la
détermination de la structure des communautés endolithiques, la fonction
dans ces écosystèmes extrêmes ou encore les relations inter-spécifiques
entre les espèces.
