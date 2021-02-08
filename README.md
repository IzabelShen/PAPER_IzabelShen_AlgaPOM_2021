# PAPER_IzabelShen_AlgaPOM_2021

# the effects of phytoplankton deposition on water-sediment interface bacterial communities

A repository for R code,  scripts from the manuscript "Quality of phytoplankton deposition structures bacterial communities at the water-sediment interface"

by Dandan Izabel-Shen, Séréna Albert, Monika Winder, Hanna Farnelid, Francisco J. A. Nascimento


Corresponding author: Dandan Izabel-Shen, Department of Ecology, Environment and Plant Sciences, Stockholm University; email: dand.shen@gmail.com

_This work is Under Review_ 


## Data
The FASTQ files and associated metadata are available in the European Nucleotide Archive under the accession number PRJEB39288. 

## Installation

The following packages (and their dependencies) are required to run the whole analysis

```
library("vegan"); packageVersvegan") # 2.5.6
library("phyloseq"); packageVersion("phyloseq")  # 1.30.0
library("doParallel"); packageVersion("doParallel") # 1.0.15
library("foreach"); packageVersion("foreach") # 1.5.0
library("DESeq2"); packageVersion("DESeq2") # 1.26.0
library("ggplot2"); packageVersion("ggplot2") # 3.3.2
library("BiocManager"); packageVersion("BiocManager") # 1.30.10"
library("gplots"); packageVersion("gplots") # 3. 1. 0
library("gridExtra");packageVersion("gridExtra") # 2.3
library("picante");packageVersion("picante") # 1.8.2
library("dplyr");packageVersion("dplyr") # 1.0.2
library("stats");packageVersion("stats") # 3.6.2
library("RColorBrewer");packageVersion("RColorBrewer") # 1.1.2
library("VennDiagram"); packageVersion("VennDiagram") # 1.6.20
library("ggpubr"); packageVersion("ggpubr") # 0.4.0
library("tidyverse"); packageVersion("tidyverse") # 1.3.0
library("car"); packageVersion("car") # 3.0.10
library("agricolae"); packageVersion("agricolae") # 1.3.3
library("scales"); packageVersion("scales") # 1.1.1.
library("graphics"); packageVersion("graphics") # 3.6.0

```

## How to run

The R scripts in the base folder should be run interactively (for instance with RStudio).

The folder 'InputFiles' contains all data files required to run the analyses. In the repository, the input files are compressed with gzip. Before you proceed with the analyses, the files should be uncompressed.

The folder 'R_analysis' contains all scripts written for the analyses. 
