# EDMS practical course

This is the tutorial about bulk CUT&Tag data analysis for EDMS practical course. 

Before starting this tutorial:

1. Ensure you have a working terminal environment on your system.
2. Verify you have admin permissions to install software packages.
3. (**optional**) Know your shell type (e.g., [bash or zsh](https://apple.stackexchange.com/questions/361870/what-are-the-practical-differences-between-bash-and-zsh)). Run `echo $SHELL` to check.

## Part 1 : Install deepTools and MACS3

If you will run this on your local computer, make sure you have [(mini)conda](https://docs.anaconda.com/miniconda/) installed. Then we can install [deepTools](https://deeptools.readthedocs.io/en/develop/content/installation.html) and [MACS3](https://macs3-project.github.io/MACS/docs/INSTALL.html):

```
# Install to a new environment specified by -n
conda create -n epign_practical -c bioconda -c conda-forge macs3 deeptools
```
Due to compatibility issues with the ARM architecture (e.g., M1 or M2 Macs), deepTools may not install correctly via conda, so pip installation is recommended.

```
# Install to MACS3 to a new environment first 
conda create -n epign_practical -c bioconda -c conda-forge macs3
conda activate epign_practical

# Then install deepTools via pip
pip install git+https://github.com/deeptools/deepTools.git@develop
```

**Tips:** there are some small bugs in the release version so we install the development version here.

Then we can check the package installation:

```
pip list

# Or check these 2 packages specifically 
pip show deeptools
pip show macs3
```

## Part 2 : Install bedtools2
Another thing we need is bedtools follow the instruction [here](https://bedtools.readthedocs.io/en/latest/content/installation.html). The easiest way, without relying on additional package managers, is to download the *.tar.gz* file from the [release](https://github.com/arq5x/bedtools2/releases) tab on github, then run this:

```
# At the download path
tar -zxvf bedtools-2.31.1.tar.gz # Change the name if neccesary 
cd bedtools2
make
```
Adding bedtools2 to your $PATH allows you to run the bedtools command from any directory in your terminal.
```
export PATH="/PATH/TO/bedtools2/bin:$PATH" 
```
Replace /PATH/TO/bedtools2/bin with the actual path to the bin directory of your bedtools2 installation. To check if we can use bedtools now, we can run:

```
bedtools -h
```

### Optional: Making the export PATH command permanent

The export PATH command described earlier is not permanent; this means you’ll need to re-run it every time you open a new terminal session. To avoid this, you can add the command to your shell configuration file (e.g., ~/.zshrc, ~/.bash_profile, or ~/.bashrc), depending on your shell type.

#### Step 1 - Identify your shell

To determine which shell your terminal is using, run the following command:
```
echo $SHELL
```
- If the output contains zsh (e.g., /bin/zsh), you are using **zsh**, which is the default for most macOS systems.
- If the output contains bash (e.g., /bin/bash), you are using **bash**.

#### Step 2 - Edit the configuration file

Open the appropriate configuration file for editing:
```
nano ~/.zshrc  # For zsh (most macOS users)
nano ~/.bash_profile  # For bash
```
Add the following line to the end of the file:

```
export PATH=/PATH/TO/bedtools2/bin:$PATH
```
Replace /PATH/TO/bedtools2/bin with the actual path to the bin directory of your bedtools2 installation.

#### Step 3 - Save and reload the configuration

After saving the changes, reload the configuration file to apply the updated PATH:

```
source ~/.zshrc  # For zsh
source ~/.bash_profile  # For bash
```

## Part 3: install R packages 

Before start this part, make sure you have [R and Rstudio](https://posit.co/download/rstudio-desktop/) installed on your local computer. 

Many genomic databases and tools provide R interfaces, making R a powerful environment for genomic data analysis. In this tutorial, we will primarily use:

1. [ChIPseeker](https://www.bioconductor.org/packages/release/bioc/html/ChIPseeker.html) for functional annotation.
2. [UpSetR](https://github.com/hms-dbmi/UpSetR) for visualization.

Additionally, we will use the genomic annotation database [TxDb.Hsapiens.UCSC.hg38.knownGene](https://www.bioconductor.org/packages/release/data/annotation/html/TxDb.Hsapiens.UCSC.hg38.knownGene.html), and [org.Hs.eg.db](https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html)

```
# ChIPseeker
if (!require("BiocManager", quietly = TRUE)) 
	install.packages("BiocManager") 
BiocManager::install("ChIPseeker")

# UpSetR
install.packages("UpSetR")

# Annotation database 
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
BiocManager::install("org.Hs.eg.db")
```
**(Optional)** Bioconductor provides a wide range of annotation databases for various species and genome versions. You can explore them here: [R Annotation Databases](https://bioconductor.org/packages/3.20/data/annotation/).

Load the installed libraries in R before using them:

```
library(ChIPseeker)
library(UpSetR)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
```