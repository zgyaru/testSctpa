# testSctpa

A portable toolkit for pathway activity score calculation

## Introduction
Single-cell RNA sequencing (scRNA-seq) analysis enables researchers to uncover more refined and novel cell clusters, which have greatly advanced our understanding of cellular states. Pathway activity score (PAS) analysis has been applied to transform the gene-level data into explainable gene sets representing biological processes or pathways to uncover the potential mechanism of cell heterogeneity. This package provide six portable interfaces for PAS calculation tools and abundant pathway databases in human and mouse.

## Install
```
devtools::install_github('zgyaru/testSctpa')
```

## Quite start

### Step 1. loading test data
```
library(testSctpa)
counts = load_counts()
```
### Step 2. selecting a pathway database used for biological annotation
```
kegg = getPathways(species='mouse', pathway='kegg')
```
### Step 3. Calculating pathway activity score
```
pas = calVision(counts,kegg)
#pas = calGsva(counts,kegg)
#pas = calAUC(counts,kegg)
#pas = calSSgsea(counts,kegg)
#pas = calPlage(counts,kegg)
#pas = calZscore(counts,kegg)
print(pas[1:5,1:5])
```



## Pathway Detials
#### `human`
| Name | Detials  | Number of gene sets |
| - | :-: | -: |
|kegg | KEGG pathway database | 186 |
|scSignature | cell type signature | 671 |
| hallmarker | Hallmark gene sets | 50 |
| CGP | genetic and chemical perturbations | 3297 |
|biocarta | BioCarta pathway database | 289 |
|PID | PID pathway database | 196 |
|reactome | Reactome pathway database | 1532 |
|TFT | transcriptional factor targets | 1137 |
|CGN | cancer gene neighborhoods | 427 |
|CM | cancer models | 431|
|GO.bp | GO biological process | 7530 |
|GO.cc | Co cellular Component | 999 |
|GO.mf | GO molecular fucntion | 1663|
|OncoG | oncogenic signatures | 189 |
|Immu | immunologic signatures | 4872 |
|panther | protein annotation through evolutionary relationship | 94 |
|humancyc | human metabonomics | 127 |
|ASCN2 | retinoblastoma protein | 102 |
|pharmgkb | pharmacogenomics | 37 |

#### `mouse`
|Name | Detials  | Number of gene sets|
|- | :-: | -: |
|kegg | KEGG pathway database | 259|
|panther | protein annotation through evolutionary relationship | 151|
|mousecyc | mouse metabonomics | 321|
|biocarta | BioCarta pathway database | 176|
|reactome | Reactome pathway database | 1396|
|TFT | transcriptional factor targets | 373|
|GO.bp | GO biological process | 8203|
|GO.cc | Co cellular Component | 1082|
|GO.mf | GO molecular fucntion | 3240|
|drug | drug related | 844|



## Reference
>[1] Hänzelmann S, Castelo R, Guinney J (2013). “GSVA: gene set variation analysis for microarray and RNA-Seq data.” BMC Bioinformatics.

>[2] Aibar S, Bravo Gonzalez-Blas C, Moerman T, Huynh-Thu V, Imrichova H, Hulselmans G, Rambow F, Marine J, Geurts P, Aerts J, van den Oord J, Kalender Atak Z, Wouters J, Aerts S (2017). “SCENIC: Single-Cell Regulatory Network Inference And Clustering.” Nature Methods.

>[3] DeTomaso D, Jones MG, Subramaniam M, Ashuach T, Ye CJ, Yosef N (2019). Functional interpretation of single cell similarity maps. Nat Communication.