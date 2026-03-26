This repository contains supplemental material and all code necessary to reproduce the analyses reported in C.F. Webster, "Seasonal changes in the transcriptome of bat eyes reveal adaptive changes to exploit summer twilight."

1. Supplemental Tables (Supplemental Tables CWebster.xlsx)<br>
   This Excel file includes the supplemental data tables described in the supplemental material section of the manuscript:<br>
   - Supplemental Table S1: Metadata for biological replicates for both seasonal and daily datasets
   - Supplemental Table S2: Final list of target genes that were examined with functional categorization. Data includes gene ID, gene name, ensembl transcript ID as reported in species' genome, and functional category the gene was assigned to.<br>
   - Supplemental Table S3: Significant seasonal genes (DESeq2 padj < 0.05) with WGCNA and mFuzz test results.<br>
   - Supplemental Table S4: Significant daily genes (DESeq2 padj < 0.05) with WGCNA and mFuzz test results.<br>

2. HTSeq-Count Results (HTSeq Counts.zip)<br>
A zip file containing per-sample HTSeq-count text files derived from raw data files. Each file corresponds to a biological replicate for both the daily and seasonal datasets. A total of 45 text files are contained within the zip folder.

3. Code for all analyses<br>
- CWebster_eyepaper.R<br>
  Performs DESeq2 differential expression analysis, WGCNA analysis for identifying expression modules, mFuzz analysis for clustering expression data over time, and goprofiler2 functional enrichment analysis<br>
    - Input: HTSeq-count text files, sample metadata information.<br>
    - Output: DESeq2 results, WGCNA correlations and pvalues (module-gene and trait-gene), mFuzz cluster assignment and graphs, GO analysis results using custom GMT file<br>
    - R libraries required: DESeq2, sva, WGCNA, mFuzz, gprofiler2<br>
- CWebster_quantum_catch.R<br>
  Performs quantum catch modeling<br>
    - Input: Global to horizontal irradiance values (300-700 nm), expression data from DESeq2 results of opsins<br>
    - Output: Raw quantum catch models (opsin sensitivity and irradiance), total quantum catch model (raw QC weighted by opsin expression), relative quantum catch (relative contribution of each opsin)<br>
    - R libraries required: pavo<br>
- hprc_expression<br>
  Command-line batch script for trimming and QC check of raw fastq reads, aligning validated reads to species' genome, and using generated BAM file to get read counts.<br>
    - Utilizes Trim_Galore/0.6.7, STAR/2.7.10b, SAMtools/1.17, HTSeq/2.0.1<br>

Repository Structure<br>
├── Supplemental Tables CWebster.xlsx<br>
├── HTSeq Counts.zip<br>
├── CWebster_eyepaper.R<br>
├── CWebster_quantum_catch.R<br>
├── hprc_expression<br>
└── README.md<br>

For questions, please contact:
Cara Webster
cfjwebster17@tamu.edu
