<img src="materials/images/introduction-to-epigenetics-cover.png"/>

# **Introduction to Epigenetics**

`🕒 This module should take less than 1 hour to complete.`

`✍️ This notebook is written using Python.`

Epigenetics is a field of study focused on changes in DNA that do not involve alterations to the underlying sequence. In comparison to static nature of the genome, epigenetics is dynamic, and change substantially during development, aging, cancers and in response to environmental exposure.

Chemical modification of DNA and the protein that interact with DNA can change the degrees to which genes are switched **ON** or **OFF**, causing an effect on the phenotype. The most common epigenetic regulation include DNA methylation, histone modifications, and non-coding RNAs.

The collection of all epigenetic changes in a genome is called an epigenome. In this module, you will learn approaches to study epigenomes, and interpret them.

<img src="materials/images/epigenetic-regulation.png"/>

<div class="alert alert-block alert-info">
<h3>⌨️ Keyboard shortcut</h3>

These common shortcut could save your time going through this notebook:
- Run the current cell: **`Enter + Shift`**.
- Add a cell above the current cell: Press **`A`**.
- Add a cell below the current cell: Press **`B`**.
- Change a code cell to markdown cell: Select the cell, and then press **`M`**.
- Delete a cell: Press **`D`** twice.

Need more help with keyboard shortcut? Press **`H`** to look it up.
</div>

---

## **Epigenetic regulation in development and disease**

Here we have two examples around epigenetics regulation:

1. The first example shows a stem cell uses epigenetic factors in the cell differentiation process.

<img src="materials/images/epigenetic-regulation-development.png"/>

2. The second example indicates our lifestyle can switch ON or OFF genes. Studies show smoking and drinking turn on genes that are associated with the development of addiction.

<img src="materials/images/epigenetic-regulation-disease.png"/>

---

## Epigenetic assays

**DNA methylation analysis**

Cytosine methylation (5-methylcytosine, 5mC) is one of the main covalent base modifications in eukaryotic genomes, generally observed on CpG dinucleotides.  Once methyl compounds are present on a gene, the gene is muted, or turned off. As a result, no protein is generated. That is how DNA methylation regulates gene expression.

Genome-wide DNA methylation can be mapped using either Whole Genome Bisulphite Sequencing (WBGS), Reduced-Representation Bisulfite Sequencing (RRBS), Methylation-sensitive Restriction Enzyme (MRE) and immunoprecipitation based assays.

<img src="materials/images/dna-methylation.png"/>

**Histone modification analysis**

DNA is wrapped around a protein complex called the histone complex. Histones form a chain of beads along the DNA. Histone modifications at specific locations (e.g., lysine acetylation and methylation) can affect whether a gene is wrapped or unwrapped ("ON" versus "OFF"). This alters the expression of genes.

Proteins that read genes cannot reach DNA wrapped tightly around histones. Consequently, this mechanism switches some genes off because the DNA is wrapped around the histone proteins, and is inaccessible to the proteins that read DNA to make RNA, whereas other genes get expressed because they are not wrapped around histones.

Genome-wide histone modifications can be measured using ChIP-Sequencing. ChIP-Seq, a combination of chromatin immunoprecipitation (ChIP), and massively parallel sequencing, delves into the interactions between proteins, DNA, and RNA, revealing critical regulatory events in many biological processes and disease states. ChIP-Seq is used to identify transcription factor binding sites, track histone modifications across the genome, and constrain chromatin structure and function.

<img src="materials/images/histone-modification.png"/>

---

## The data

In this module we will take a look at the epigenome data (DNA methylation and histone modifications) from human reference epigenomes generated as part of **NIH Roadmap Epigenomics Program**. You can learn more at: https://egg2.wustl.edu/roadmap/web_portal/index.html


**Differentially Methylated Region (DMR) calls across reference epigenomes**

As a general resource for epigenomic comparisons across all epigenomes, Differentially Methylated Region (DMR) were defined using the Lister et al method (Lister et al., 2013), combining all Differentially Methylated Sites (DMSs) within 250bp of one another into a single DMR, and excluding any DMR with less than 3 DMSs.

For each DMR in each sample, average methylation level was computed, weighted by the number of reads overlapping it (Schultz et al., 2012). This resulted in a methylation level matrix with rows of DMRs and columns of samples.

**ChIP-seq peak calls**

For each gene, a region of 10.000bp around the transcription start site of the gene is extracted (5000bp upstream and 5000bp downstream). This region is binned in 100 bins of 100bp. For each bin, five core histone modification marks are counted.

**Whole Genome Bisulphite Sequencing data from 111 reference epigenomes**



```python
import sys
import os
# import pyathena
import pandas as pd
from pathlib import Path
```


```python
file_path = Path("data") / "WGBS_DMRs_v2.tsv"
wgbs = pd.read_csv(file_path, sep="\t")
```


```python
print(wgbs.head(10))
```

The data matrix shows rows of DMRs and columns of samples including chromosome number, location 'start' and 'end'.

1. `chr`: Chromosome.
2. `start`: Start of the location.
3. `end`: End of the location.

Metadata on the samples used in the analysis are available for downloading below: https://docs.google.com/spreadsheets/d/1yikGx4MsO9Ei36b64yOy9Vb6oPC5IBGlFbYEt-N6gOM/edit#gid=15

If the meta data link above doesn't work, you could download the file to your local:


```python
from IPython.display import FileLink, HTML

# absolute or relative path to the Excel file inside the notebook’s workspace
file_path = "data/Roadmap.metadata.qc.jul2013.xlsx"

# 1-click download link (renders as a hyperlink in the output cell)
FileLink(file_path)

# ――― optional: make it look like a button ―――
HTML(f"""
<a download href="{file_path}"
   style="font-family:system-ui;padding:8px 14px;border:1px solid #ccc;
          border-radius:6px;text-decoration:none;">
  ⬇️ Download Roadmap.metadata.qc.jul2013.xlsx
</a>
""")

```

**Querying DNA methylation and histone modification states of genes in UCSC genome Browser**

`FractionMethylation.tar.gz` files are used for visualization of DNA methylation states, and provides fractional methylation calls at CpG. It contains 25 files. One for each chromosome.

Format of each file: Tab separated table CpGs (rows) X epigenomes (columns)
Methylation calls are round to two decimal digits.

Each file has the same matrix format:
- The first column is a position of C or G in the CpG
- The rest of the columns are epigenomes.

Only CpG info is present (as it came from EDACC files); for those CpG where coverage was <=3 we - for both coverage and Methylation (as missing data)


For ChIP-Seq data visualization, negative log10 of the Poisson p-value of ChIP-seq counts relative to expected background counts local were used. These signal confidence scores provides a measure of statistical significance of the observed enrichment.


The NCBI RefSeq Genes composite track shows human protein-coding and non-protein-coding genes taken from the NCBI RNA reference sequences collection (RefSeq) [hg19 refGene]. You could learn more at the following:
- https://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/
- https://genome.ucsc.edu/cgi-bin/hgTables?db=hg19&hgta_group=genes&hgta_track=refSeqComposite&hgta_table=refGene&hgta_doSchema=describe+table+schema (Schema for NCBI RefSeq - RefSeq genes from NCBI)

Visualize DNA methylation and histone modification (an active marker-H3K4me3) here. http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr6%3A31128153%2D31142413&hgsid=1467202695_p1HPrLZa2tMzJREfrrGplnR1tNsc

Compare Pou5f1 promoter methylation and H3K4me3 levels between ESCs and neuronal progenitor cultured cells. You could query other crucial genes in ESC maintenance for e.g., Nanog, Sox2. TDGF1, LEFTY1, GDF3, FOXD3 and in neuronal progenitor cells, PAX6, SIX3, LHX2, OTX2, PLZF, SOX1, FOXG1.

<img src="materials/images/visualization.png"/>

---

**Reference**

- Roadmap Epigenomics Consortium., Kundaje, A., Meuleman, W. et al. Integrative analysis of 111 reference human epigenomes. Nature 518, 317–330 (2015). https://doi.org/10.1038/nature14248.
https://www.illumina.com/content/dam/illumina-marketing/documents/products/appnotes/appnote-methylseq-wgbs.pdf

---

# Contributions & acknowledgment

Thank the following team to work on this module:

- **Module Content:** Abtin Tondar, Mohan Babu
- **Engineering:** Amit Dixit
- **UX/UI Design & Illustration:** Kexin Cha
- **Video Production:** Francesca Goncalves
- **Project Management:** Amir Bahmani, Kexin Cha

---

Copyright (c) 2022 Stanford Data Ocean (SDO)

All rights reserved.
