The data processing and analysis outline is here: [[Bulk CUT&Tag data analysis chart.canvas|Bulk CUT&Tag data analysis chart]]
## Step 0 - Quality control 

Assuming we already have the FASTQ files from the sequencing facility, what can we do with the data?

[FASTQ files](https://en.wikipedia.org/wiki/FASTQ_format) contain the nucleotide sequence of each read, along with metadata (such as read identifiers and indices) and quality scores for each base in the sequence. The kind staff at the sequencing facility might have already done these two steps:

1. Raw data assessment using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [MultiQC](https://multiqc.info/). 
2. Trimming the adapters using [Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) or [cutadapt](https://cutadapt.readthedocs.io/). *It is therefore important to know which adapter (pairs) were used in the library preparation.*

If you want to do the QC by yourself, you can use the pre-installed modules on the uni's [high-computing clusters](https://scitas-doc.epfl.ch) :) 

---
**Now we start with the sequencing results: the FASTQ files. The following steps can give you an overview of what we are doing and why we are doing them. These might seem like a lot, but [snakePipes](https://snakepipes.readthedocs.io/en/stable/index.html) nicely compiled them and we will learn how to use it in this [section](#step-0-4-snakepipes).** 
## Step 1 - Alignment 

By mapping these reads to a reference genome, we can determine where in the genome our library is derived from. Several things we need to decide for alignment:

1. Which reference genome?
	1. Most of the time by Googling the species we can find what we want (typically from [UCSC](https://hgdownload.soe.ucsc.edu/downloads.html), [Ensembl](https://www.ensembl.org/index.html), [NCBI](https://www.ncbi.nlm.nih.gov/refseq/), [GENECODE](https://www.gencodegenes.org/human/) etc.). Many sources also contain the pre-build indices that are required for aligner to work (check [why](https://www.biostars.org/p/212594/)). 
	2. Pay attention to the version of the reference genome and the *style differences*: UCSC style includes the *chr* prefix (e.g chr1, chrM) while Ensembl/NCBI Style does not (e.g. 1, MT) 
2. What if spike-ins are used?
	1. The reference genome wouldn't have the spike-ins. To have a *hybrid* reference genome and corresponding index, we can use [createIndices](https://snakepipes.readthedocs.io/en/stable/content/workflows/createIndices.html) pipeline from snakePipes.  
3. Which aligner?
	1. There are multiple popular aligners, such as [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml), [BWA](https://bio-bwa.sourceforge.net), [BWA-MEM2](https://bio-bwa.sourceforge.net). Note that the choice needs to correspond to the indices you are using (e.g., Bowtie2 indices won’t work with BWA). 
	2. As a side note, the original CUT&Tag paper ([Kaya-okur et al. 2019](https://www.nature.com/articles/s41467-019-09982-5#Sec8)) used Bowtie2 (v2.2.5) with specific parameters. 
## Step 2 - Conversion and filtering 

Most aligners typically generate [SAM](conda create -n snakePipes -c mpi-ie -c conda-forge -c bioconda snakePipes) (Sequence Alignment Map) files, which are text-based format that are human-readable but large in size. [`samtools`](https://www.htslib.org) can convert SAM to [BAM](https://en.wikipedia.org/wiki/Binary_Alignment_Map) (Binary Alignment Map), which is also compatible with downstream tools. 

`samtools` is a powerful suite of programs with very useful commands. We commonly use the ones to mark duplicate reads (`markdup`), filter duplicates and low mapping quality reads (`filter --mapq --dedup`) and calculate the mapping statistics (`flagstat`). Other tools that can do BAM file operations (e.g. [Sambamba](https://lomereiter.github.io/sambamba/)).

---
**The following steps can start once we have the BAM files ready.**
## Step 3.1 - Coverage files 

We can first load the BAM files in to [IGV](https://www.igv.org) tools (NOTE: select the correct reference genome first!). If we zoom in a lot, we can see the exact nucleotide sequences of each read, and the mapping quality etc. This detailed information makes the BAM file still too big. To make it even smaller, we will use the powerful `deepTools`. [`bamCoverage`](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html) command can convert BAM file to [BigWig file](https://genome.ucsc.edu/goldenpath/help/bigWig.html), which is a *coverage file* calculated as the number of reads per bin. 

```
bamCoverage -b reads.bam -o coverage.bw
```

**(optional) Spike-in normalisation**. There are different options for read coverage normalisation. If spike-ins were used, we can normalise by the number of reads mapped to spike-in genome.

```
# Get the scale to the spike-in
multiBamSummary bins --region spike_in_region.bed --bamfiles reads.bam --scalingFactors scale_factor.txt 
```

`py-deepTools` is also pre-compiled in the HPC :)
## Step 3.2 - Peak calling 

If we load the BigWig files to the IGV, we see that they are histograms that have somewhat of *wave* shapes. [**Peak calling**](https://en.wikipedia.org/wiki/Peak_calling) is a computational method used to identify regions with high signals (enriched read alignment). 

There are again many computational tools out there for peak calling, such as [MACS3](https://github.com/macs3-project/MACS) and [SEACR](https://seacr.fredhutch.org/). Even though many of them are initially developed for ChIP-seq, by adjusting parameters, they can be adapted to ATAC-seq/ CUT&Tag-seq. An example could be:

```
macs3 callpeak -t reads.bam -f BAM --broad -g hs \
	--outdir MACS3_broad_peaks -B --broad-cutoff 0.1 --nolambda
```
## Step 4 - Visualisation & QC 

Once we have the BigWigs, we can visualise the average signal over regions of interest using `deepTools`:

```
computeMatrix scale-regions -S coverage.bigwig -R regions.bed -o matrix.gz
plotHeatmap -m matrix.gz -o heatmap.png
```

Here the `regions.bed` can be the peaks we called from step 3.2, then this heat-map allows us to visualise how 'good' is our peak calling or the data. For example, regions that are enriched with H3K27me3 should probably have less H3K27ac signals. 

If you notice something is off, maybe try these following QC/troubleshoot:

1. **Fraction of reads in peaks ([FRiP](https://www.encodeproject.org/data-standards/terms/#enrichment))**
	If the data quality and peak calling both worked, then we would expect more reads in the peaks (as if the signal-to-noise ratio is high). There are many tools to calculate FRiP, one of which is [`featureCounts`](https://subread.sourceforge.net/featureCounts.html) from `subread` packages:
	
	```
	# covert BED (the peaks) to SAF
	awk 'BEGIN{FS=OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"}{print $4, $1, $2+1, $3, "."}' peaks.bed > peaks.saf
	
	# featureCounts
	featureCounts -p -a eaks.saf -F SAF -o readCountInPeaks.txt reads.bam
	```
	
	The good thing about `featureCounts` (or `subread`) is that it is pre-compiled in the HPC cluster, but the annoying thing is that it need SAF file format. SAF start position is 1-based but BED is 0-based (like the difference between R and python), thus we need the conversion using `awk`.

2. **Replicates correlation**
	If we have multiple replicates, we can check if there is any outlier. `deepTools` provides tools to calculate the correlation from either BigWig files or BAM files:
	```
	# from BAM files
	multiBamSummary bins --bamfiles file1.bam file2.bam -o results.npz
	# from BW files
	multiBigwigSummary bins -b file1.bw file2.bw -o results.npz

	# calculate spearman correlation
	plotCorrelation -in results.npz -c spearman -p heatmap -o cor_plot.png
	# or PCA 
	plotPCA -in results.npz -o pca.png
	```

	Essentially [multiBamSummary](https://deeptools.readthedocs.io/en/develop/content/tools/multiBamSummary.html) or [multiBigwigSummary](https://deeptools.readthedocs.io/en/develop/content/tools/multiBigwigSummary.html) calculates the coverage on a binned genome or a given region (in BED format), then [plotCorrelation](https://deeptools.readthedocs.io/en/develop/content/tools/plotCorrelation.html) can compute Pearson or [Spearman](https://en.wikipedia.org/wiki/Spearman%27s_rank_correlation_coefficient) correlation, and [plotPCA](https://deeptools.readthedocs.io/en/develop/content/tools/plotPCA.html) can compute PCA.  
2. **Consensus of peaks**
	Another thing we can check when we have replicates is the consensus of peaks. Do we choose the peaks that are present in any/half/all of the replicates? Do we only keep the peaks that are larger/smaller than certain base pairs (bp) to remove the noises? Do we merge the peaks that are close to each other?
	Peaks are normally in BED file format, and one of the most powerful BED file operator is [`bedtools`](https://bedtools.readthedocs.io/en/latest/). Check the manual to see which command is the most suitable for specific application.
3. **Pre-defined regions (from databases)**
	Alternatively, a quick to check if the data make senses is to `plotHeatmap` on a given region, like transcription start site (TSS), enhancer regions or any the regions we are interested in. This way we can check the data quality without relying on our own peak calling. **NOTE**: consider if the pre-defined regions will contain the markers of your experiment :)
## Step 0-4: snakePipes 

[snakePipes](https://snakepipes.readthedocs.io/en/stable/index.html) are pipelines built using [snakemake](https://snakepipes.readthedocs.io/en/stable/snakemake.readthedocs.io) and _python_ for the analysis of epigenomic datasets. To use the pipelines, we can set up as:

1. Setting up snakePipes following the [instructions](https://snakepipes.readthedocs.io/en/stable/content/setting_up.html) 
	**Tips:** Before running [`snakePipes createEnvs`](https://snakepipes.readthedocs.io/en/stable/content/setting_up.html#create-the-conda-environments), make sure to change the conda-prefix in the profile file.
2. Reference genome by one of these 2 options:
	1. Download the pre-made indices as [here](https://snakepipes.readthedocs.io/en/stable/content/setting_up.html#download-premade-indices). **Make sure to change the organism profiles in the directory** (we can get the file locations by running `snakePipes info`)
	2. Make your own indices by running [`createIndices`](https://snakepipes.readthedocs.io/en/stable/content/workflows/createIndices.html) pipelines. **This is often needed when spike-ins are used.**

Once we have the environment and reference genome ready, we are ready to run one of the pipelines. [DNAmapping](https://snakepipes.readthedocs.io/en/stable/content/workflows/DNAmapping.html) covers the first few steps of both ChIPseq and ATACseq pipelines (alignment, filtering and conversion to BigWig files), while the other two include additional peak calling and differential peak analysis. In practice, we often use DNAmapping to have a flexible downstream analysis. 

DNAmapping also offers some visualisation and QC when the number of samples are less than 10. Check previous step 4 for more of what we can do :)

---
**The following lists some common downstream analysis we can do. It is again case-specific: what questions do we want to answer?**
## Step 5.1 - Differential enrichment analysis

Often we want to compare the differential enrichment between datasets, i.e., H3K27me3 differences between condition 1 and condition 2, or H3K27me3 and H3K27ac differences in the same condition. 

The most straightforward way is use the (spike-in normalised) BigWig files, and again `deepTools`. 
1. Compare 2 BigWig files. 
	[`bigwigCompare`](https://deeptools.readthedocs.io/en/develop/content/tools/bigwigCompare.html) can compare 2 BigWig files in many different ways, and the output can be in either BigWig (useful for IGV visualisation and other downstream analysis) or BED (useful to know the exact values at binned regions).
	
	```
	bigwigCompare -b1 sample1.bw -b2 sample2.bw -o log2.bw
	```
	
1. Calculate multiple BigWig signals on binned genome/pre-defined regions.
	If we want to get the enrichment of multiple BigWigs, [`multiBigwigSummary`](https://deeptools.readthedocs.io/en/develop/content/tools/multiBigwigSummary.html) is very helpful.

	```
	# Signal of each BW on binned genome
	multiBigwigSummary bins -b file1.bw file2.bw \ # can be multiple bw
		-o results.npz --outRawCounts results.tab
		
	# Signal of each BW on pre-defined BED file
	multiBigwigSummary BED-file -b file1.bw file2.bw \
		-o results.npz --outRawCounts results.tab \
		--BED selection.bed
	```
	The .*npz* file can be load to python for more (statistical) analysis, and *.tab* gives the full information of scores/signals/enrichment per genomic bin/BED. 

Here using `deepTools` we can process the BAM/BW files into matrices (Dim $X_{samples}, Y_{regions}$), something any program language can deal with :) Of note, many tools build for ChIP-seq analysis can also be applied here (e.g. [DiffBind](https://bioconductor.org/packages/release/bioc/html/DiffBind.html) or `bdgdiff` from [MACS3](https://macs3-project.github.io/MACS/docs/bdgdiff.html)).
## Step 5.2 - Combinatorial patterns

As we know, combinatorial histone modifications can mark functional regions, like promoter/enhancers etc. How do we identify the combinatorial regions?

Now we know `multiBigwigSummary` can give us the matrix with scores per bin/region, we can surely do some clustering there. But even better, `deepTools` already have such a function in [`plotHeatmap`](https://deeptools.readthedocs.io/en/develop/content/tools/plotHeatmap.html):

```
computeMatrix scale-regions -S coverage.bigwig -R regions.bed -o matrix.gz
plotHeatmap -m matrix.gz \
	--kmeans 4 \ # other clustering algorithm possible 
	--outFileSortedRegions heatmap_sortedRegions.bed \ # the cluster results
	--outFileNameMatrix cluster_martrix.gz \
	-o heatmap.png
```

 [`computeMatrix`](https://deeptools.readthedocs.io/en/develop/content/tools/computeMatrix.html) also somewhat calculate the score per region, but it further bins the region to 10bp (a user-defined parameter), so the matrix is of dimension $X_{samples}, Y_{regions*bins}$. Again this matrix can also be handled by any programming language. 
## Step 5.3 - Functional annotation 

All the peak files (in BED format) can be functionally annotated. Many of the genomic databases have an R interface and there are tools to use them as well, such as [ChIPseeker](https://www.bioconductor.org/packages/release/bioc/html/ChIPseeker.html):

```
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# Specify the TxDb object correctly
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
```

**Note:** Here the TxDb uses UCSC style (the one that uses *chr1* instead of *1*). Check what format the peak/BED files are using -  change the BED or TxDB style if needed. 

Once we have the database ready, we can use the function `annotatePeak` following the instructions [here](https://www.bioconductor.org/packages/devel/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html#peak-annotation). 

**(optional)** Another popular yet more advanced way to do functional annotation is [ChromHMM](https://ernstlab.github.io/ChromHMM/), which can use various histone modifications to discovery the major re-occuring combinations of markers. Check the [manual](https://ernstlab.github.io/ChromHMM/ChromHMM_manual.pdf) to see what how the software can `LearnModel` with in-house dataset. Alternatively, we can use the command `OverlapEnrichment` to overlap our results to a pre-trained model (such as models provided on [NIH Roadmap Epigenomics Mapping Consortium](https://egg2.wustl.edu/roadmap/web_portal/chr_state_learning.html#core_15state) , and this [full stack model](https://github.com/ernstlab/full_stack_ChromHMM_annotations) trained with over 1000 datasets).
