# SDRS_Transcriptional_Dynamics

This repository contains the code for the analysis of bulk RNA-seq data for 3, 5, and 6 hours of sleep deprivation (SD3, SD5, SD6) and 2 and 6 hours of subsequent recovery sleep (RS2, RS6) in adult male wild-type mice. This dataset is an integration of new unpublished data with the publicly available RNA-seq datasets GSE140345 (Hor et al., 2019) and GSE113754 (Ingiosi et al., 2019). There were circadian home cage (HC) controls for each time point: HC3, HC5, HC6, HC7, and HC12. The samples for the RS2 and HC7 time points have not previously been published, but they were collected at the same time as those from GSE113754. Samples for the SD5 and HC5 time points were generated in GSE113754 and extracted from frontal cortex tissue. Samples for the SD3, HC3, SD6, HC6, RS6, and HC12 time points were generated in GSE140345 and extracted from cerebral cortex tissue.

## Authors

- Alex Popescu (alex.popescu@yale.edu)
- Katie Ford (kaitlyn.ford@wsu.edu)
- Caitlin Ottaway (caitlin.ottaway@wsu.edu)
- Christine Muheim (christine.muheim@wsu.edu)
- Stephanie Hicks (shicks19@jhu.edu)

## Data

Sequencing data have been deposited in the Gene Expression Omnibus Database (GEO) under accession number GSE237419. All samples besides those for the RS2 and HC7 time points were previously deposited in GEO under either GSE113754 (SD5, HC5) or GSE140345 (SD3, HC3, SD6, HC6, RS6, and HC12).

## FASTQ Files

We already had gzipped FASTQ files for the new samples and the samples from GSE113754, which were paired-end sequenced.
FASTQ files for the samples from GSE140345 are available on SRA: https://www.ncbi.nlm.nih.gov/sra?term=SRP229742

To obtain these files, we installed the SRA toolkit with the following code, as described [here](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit):

```
# The current version as of May 2024 is 3.1.1
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-centos_linux64.tar.gz
tar -vxzf sratoolkit.2.11.0-centos_linux64.tar.gz

# Append the path to the binaries to your PATH environment variable (needed to run fastq-dump)
export PATH="/data/peixoto/Programs/sratoolkit.2.11.0-centos_linux64/bin:$PATH"
```

Alternatively, SRA toolkit can be installed with Conda:
``` conda install -c bioconda sra-tools ```

We downloaded FASTQ files from SRA with the code below. Note that ```fastq-dump``` is being deprecated and has been replaced by ```fasterq-dump```, which has different options and defaults and also requires more space. See [here](https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump) for details on ```fasterq-dump```.

```
# 'files' is a file containing the list of the SRR accession numbers (one per line)
for fname in $(cat files)
do
  fastq-dump --outdir fastq --gzip --skip-technical --read-filter pass --dumpbase --clip $fname
done

wait
```

We suggest using `prefetch` and `fasterq-dump` as described at https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump and https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/.


## Transcript Quantification

Raw sequencing reads were quantified using Salmon v1.8.0. The current version as of May 2024 is 1.10.2 and can be installed according to https://salmon.readthedocs.io/en/latest/building.html#installation.

The ```Quantification/SalmonIndex.sh``` script contains the code for building the index as described [here](https://salmon.readthedocs.io/en/latest/salmon.html#preparing-transcriptome-indices-mapping-based-mode) using GENCODE release M28 and genome assembly version GRCm39. As of May 2024, files for release M35 are available at https://www.gencodegenes.org/mouse/.

The second step of Salmon's [mapping-based mode](https://salmon.readthedocs.io/en/latest/salmon.html#quantifying-in-mapping-based-mode) is ```salmon quant```. All FASTQ files are paired-end, except those from GSE140345, which are single-end.

Paired-end FASTQs were quantified using the ```Quantification/SalmonQuant_Paired.sh``` script. Single-end FASTQs were quantified using the ```Quantification/SalmonQuant_Single.sh``` script with the ```-r``` option. In both cases, we specify ```--libType A``` to allow Salmon to automatically infer the library type. With ```--numBootstraps 30```, we chose to have Salmon compute inferential replicates (bootstrap samples) to incorporate uncertainty regarding transcript abundance estimates. These are used by certain differential expression methods, like Swish, but we later dropped them since they were not needed for the gene-level analysis.

The ```tximeta``` package (v1.18.0) from the R/Bioconductor project was used to import the ```quant.sf``` files output by Salmon into R. The relevant R script is ```Quantification/AP_SDRS_Tximeta.R```. The purpose of the optional first step of creating a ```linkedTxome``` is to record the sources of the FASTA file of transcript sequences and the GTF annotation file, assisting in a reproducible analysis. The resulting JSON file is ```Quantification/gencode.vM28-salmon-index-v1.0.0-mouse-withdecoys.json```. The subsequent steps of constructing a sample table (specific to the location of the ```quant.sf``` files) and running ```tximeta``` are as described in the [vignette](https://bioconductor.org/packages/devel/bioc/vignettes/tximeta/inst/doc/tximeta.html).

Because we want ```tximeta``` to use the local annotation information and metadata and not reference data from AnnotationHub, we specify ```useHub=FALSE```. For the SummarizedExperiment object used for the gene-level analysis, we drop inferential replicates from Salmon with ```dropInfReps = TRUE```. We summarize transcript-level quantifications to the gene level using ```summarizeToGene``` on the SummarizedExperiment object from ```tximeta``` and add gene symbols. We finally extract the gene-level and transcript-level counts from the SummarizedExperiment object, and these can be found at ```Quantification/Outputs/SDRS_Gene.txt``` and ```Quantification/Outputs/SDRS_Transcript.txt```.


## RUVs and Differential Expression Analysis

The file ```Differential_Expression/Data/gse_HCSDRS_WT_SleepIntegration_salmon.txt``` is the gene-level counts matrix from ```tximeta``` but with its columns reordered by experimental condition and time of day. It serves as the input to the ```Differential_Expression/AP_SDRS_RUV.R``` script.

Based on the recovery of positive controls, we filtered out genes present less than 10x across more than 5 samples to remove genes expressed at a very low level. Selecting five samples made sense for our dataset because the majority of experimental conditions had five replicates. As a result, the original 54307 x 42 matrix is reduced to 18872 x 42. As explained in the next section, the list of 18872 genes expressed in the mouse cortex serves as a gene background for the downstream functional enrichment analysis and can be found at ```Functional_Enrichment/Data/DAVID_Gene_Background.xlsx```.

To correct for library size, we implemented upper-quartile normalization using the ```betweenLaneNormalization``` method from EDASeq (v2.34.0). The ```RUVSeq``` package (v1.34.0) from the R/Bioconductor project was used to correct for confounding biological and technical factors, such as the fact that the samples came from different labs with different sleep deprivation procedures. In particular, we used the ```RUVs``` method, specifying a set of empirical negative controls from Gerstner et al. (2016) (see ```Differential_Expression/Data/AdditionalFile4_BMC_Full.xlsx```, adj. p-value > 0.9) for the ```cIdx``` argument as described in Risso et al. (2014). For the ```scIdx``` argument, we used a matrix specifying the replicates with SD5/SD6 and HC5/HC6 grouped together to correct for batch effects and to obtain a single estimate of the effects of a long duration of SD.

By examining PCA plots following RUVs normalization and results from the downstream differential expression analysis (specifically the recovery of positive controls) for a range of k values, we found k = 3 to be optimal. The list of positive controls was also assembled in Gerstner et al. (2016) and can be found at ```Differential_Expression/Data/BMC_Genomics_AF2_SD_PosCtrls.xlsx``` for SD5-6 and at ```Differential_Expression/Data/AdditionalFile4_BMC_Full.xlsx``` (adj. p-value < 0.01) for RS2 and RS6. 

Note that for all control lists from Gerstner et al. (2016), Affymetrix probe IDs needed to be mapped to Ensembl IDs (Ensembl release 105, https://dec2021.archive.ensembl.org/index.html) and this was done using the ```Differential_Expression/AP_SDRS_Controls.R``` script. Following annotation with the ```Differential_Expression/Annotate.pl``` Perl script, we obtained the four positive and negative control lists available at ```Differential_Expression/Data/Supplemental Table S1.xlsx```.

We used edgeR (v3.38.1) for the differential expression analysis, FDR < 0.05. We applied a double thresholding approach, whereby an additional filter on expression level (log2CPM > 0) was used to maximize the ratio of positive controls recovered to empirical negative controls recovered. Another benefit was that this addressed the limitation on the number of gene IDs that can be input into DAVID for the functional enrichment analysis. However, we note that the number of DEGs we identified (over 8500 genes for the SD5-6 vs HC5-6 comparison) should be considered as the number before applying the log2CPM cutoff.


## Functional Enrichment Analysis

To examine transcriptional dynamics, different patterns of expression across the four SD and RS time points, we intersected upregulated and downregulated DEGs (FDR < 0.05, log2CPM > 0) separately across SD3, SD5-6, RS2, and RS6. All non-empty intersections were visualized using UpSet plots, created with UpSetR (v1.4.0) and ComplexUpset (v1.3.3), as seen in the ```Functional_Enrichment/AP_SDRS_UpSetPlot.R``` script. The input to the script is the output from the previous section (```Differential_Expression/Outputs/Differential_Expression_Full_Results.xlsx```) but with the double threshold on FDR and log2CPM applied and is available at ```Functional_Enrichment/Data/Supplemental_Table_S3_log2CPM_Cut.xlsx```.

To better understand how varying quantities of sleep deprivation and recovery sleep affect pathways, biological processes, and molecular functions, functional enrichment analysis was carried out using using the Database for Annotation, Visualization and Integrated Discovery v2021 (DAVID, Sherman et al., 2022 and Huang et al., 2009) with DAVID Knowledgebase v2023q3. Note that the Knowledgebase is updated quarterly; for example, the current version as of July 2024 is v2024q2. The following annotation categories were used: UniProt Biological Process, UniProt Molecular Function, and KEGG Pathways. Enrichment was determined relative to the list of 18872 expressed genes, available at ```Functional_Enrichment/Data/DAVID_Gene_Background.xlsx```. A p-value threshold (EASE Score, a modified Fisher Exact p-value) < 0.05 was used. A similarity threshold of 0.20 was used to allow for inclusive clustering.

The intesections output by the ```UpSetPlot.R``` script that were chosen for functional enrichment can be found at ```Functional_Enrichment/Data/Dec_2023_Functional_Annotation_Lists.xlsx```. From our experience, it is difficult to get informative results from functional enrichment with fewer than 50-60 genes so we excluded those intersections as well as additional intersections (namely SD5-6/RS6) that were not biologically relevant. For each of the three levels of recovery (recover within 2 hours; require 2-6 hours to recover; require more than 6 hours to recover), we also performed functional enrichment analysis on the union across the two SD time points.

The output from DAVID can be found at ```Functional_Enrichment/Data/Dec_2023_Functional_Annotation_Output.xlsx``` and serves as the input to the ```Functional_Enrichment/AP_SDRS_BubblePlot.R``` script for generating the visualizations of clustered and unclustered terms with ggplot2 (v3.5.1). Terms within a given cluster were condensed into a single informative cluster name by taking the geometric mean of their p-values and fold enrichment values with ```geometric.mean()``` from the ```psych``` package.

The complete results of functional enrichment that also include selected positive controls and hub genes can be found in Supplemental Table S5 of the associated manuscript.
