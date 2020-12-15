Generate annotation file for mapping between ENST ID, ENSG ID and gene symbol
================

``` r
library(AnnotationHub)
library(ensembldb)
library(tidyverse)
library(purrr)
# Connect to AnnotationHub
ah <- AnnotationHub()
```

``` r
ah
```

    ## AnnotationHub with 52750 records
    ## # snapshotDate(): 2020-10-27
    ## # $dataprovider: Ensembl, BroadInstitute, UCSC, ftp://ftp.ncbi.nlm.nih.gov/g...
    ## # $species: Homo sapiens, Mus musculus, Drosophila melanogaster, Bos taurus,...
    ## # $rdataclass: GRanges, TwoBitFile, BigWigFile, Rle, EnsDb, OrgDb, ChainFile...
    ## # additional mcols(): taxonomyid, genome, description,
    ## #   coordinate_1_based, maintainer, rdatadateadded, preparerclass, tags,
    ## #   rdatapath, sourceurl, sourcetype 
    ## # retrieve records with, e.g., 'object[["AH5012"]]' 
    ## 
    ##             title                                       
    ##   AH5012  | Chromosome Band                             
    ##   AH5013  | STS Markers                                 
    ##   AH5014  | FISH Clones                                 
    ##   AH5015  | Recomb Rate                                 
    ##   AH5016  | ENCODE Pilot                                
    ##   ...       ...                                         
    ##   AH87062 | org.Schizosaccharomyces_pombe.eg.sqlite     
    ##   AH87063 | org.Burkholderia_anthina.eg.sqlite          
    ##   AH87064 | org.Ascoidea_rubescens_DSM_1968.eg.sqlite   
    ##   AH87065 | org.Burkholderia_pseudomultivorans.eg.sqlite
    ##   AH87066 | org.Halogeometricum_borinquense.eg.sqlite

## Search for human data using EnsDb package.

``` r
human_ens <- query(ah, c("Homo sapiens", "EnsDb"))
human_ens
```

    ## AnnotationHub with 16 records
    ## # snapshotDate(): 2020-10-27
    ## # $dataprovider: Ensembl
    ## # $species: Homo sapiens
    ## # $rdataclass: EnsDb
    ## # additional mcols(): taxonomyid, genome, description,
    ## #   coordinate_1_based, maintainer, rdatadateadded, preparerclass, tags,
    ## #   rdatapath, sourceurl, sourcetype 
    ## # retrieve records with, e.g., 'object[["AH53211"]]' 
    ## 
    ##             title                             
    ##   AH53211 | Ensembl 87 EnsDb for Homo Sapiens 
    ##   AH53715 | Ensembl 88 EnsDb for Homo Sapiens 
    ##   AH56681 | Ensembl 89 EnsDb for Homo Sapiens 
    ##   AH57757 | Ensembl 90 EnsDb for Homo Sapiens 
    ##   AH60773 | Ensembl 91 EnsDb for Homo Sapiens 
    ##   ...       ...                               
    ##   AH73986 | Ensembl 79 EnsDb for Homo sapiens 
    ##   AH75011 | Ensembl 98 EnsDb for Homo sapiens 
    ##   AH78783 | Ensembl 99 EnsDb for Homo sapiens 
    ##   AH79689 | Ensembl 100 EnsDb for Homo sapiens
    ##   AH83216 | Ensembl 101 EnsDb for Homo sapiens

## Use the latest version of Ensembl Genomes

``` r
human_ens <- human_ens[["AH83216"]]
```

    ## loading from cache

``` r
# Extract gene-level information
hgene=genes(human_ens, return.type = "data.frame")
# Extract transcript-level information
htranscript=transcripts(human_ens, return.type = "data.frame")
# Extract exon-level information
#hexon=exons(human_ens, return.type = "data.frame")
```

## There are many Ensembl identifiers that map to more than one Entrez (NCBI) identifier.

``` r
class(hgene$entrezid)
```

    ## [1] "list"

``` r
length(which(map(hgene$entrezid,length) > 1))
```

    ## [1] 322

## Keep the first identifier for these multiple mapping cases.

``` r
hgene$entrezid <- map(hgene$entrezid,1) %>%  unlist()
```

## The records without gene name and genes corresponding to multiple records.

``` r
length(which(is.na(hgene$symbol)))
```

    ## [1] 0

``` r
sum(duplicated(hgene$symbol))
```

    ## [1] 6641

## Generate the mapping from transcript id to gene id and gene symbol

``` r
txdb=htranscript%>%dplyr::select(tx_id,gene_id)

#remove Locus Reference Genomic records, only keep ENST
txdb=txdb[grep("ENST", txdb$tx_id),]

genedb=hgene%>%dplyr::select(gene_id, symbol,entrezid)

annotations <- inner_join(txdb, genedb)
```

    ## Joining, by = "gene_id"

``` r
saveRDS(annotations,'tx2gene_grch38_ens101.rds')
```
