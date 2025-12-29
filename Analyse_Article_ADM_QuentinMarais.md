Analyse_Article_ADM_QuentinMarais
================

# Installation du package DADA2

``` r
library(dada2)
```

    ## Loading required package: Rcpp

# Chargement du jeu de données

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
srr_list <- runs$Run
writeLines(srr_list, "SRR_Acc_list.txt")
```

Ensuite je souhaite télécharger tous les runs listés dans mon nouveau
dossier .txt
