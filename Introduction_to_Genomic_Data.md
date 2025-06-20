<img src="materials/images/introduction-to-genomic-data-cover.png"/>

# Introduction to Genomics

Genomics is the branch of molecular biology concerned with the structure, function, evolution, and mapping of genomes.

In this module, you will be learning about how to process and annotate **Variant Call Format (VCF)** file using Amazon Athena. The Variant Call Format specifies the format of a text file used in bioinformatics for storing gene sequence variations.

The format has been developed with the advent of large-scale genotyping and DNA sequencing projects, such as the ***1000 Genomes Project***.

`🕒 This module should take about 30 minutes to complete.`

`✍️ This notebook is written using Python.`

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



```python
# Tables you can query
# ['default.g1000vcf_csv_int', 'default.g1000vcf_csv', 'default.g1000vcf_parquet', 'default.g1000vcf_partitioned']
# COSMIC68 Annotation Dataset ['1000_genomes.hg19_cosmic68_int']
# UCSC RefGene Annotation Dataset ['1000_genomes.hg19_ucsc_refgene_int']
```

---

# Initial Setup

We'll use the PyAthena library to get access to a database stored in AWS S3. You can read more about PyAthena here:

• https://pypi.org/project/pyathena/

• https://aws.amazon.com/athena/?whats-new-cards.sort-by=item.additionalFields.postDateTime&whats-new-cards.sort-order=desc


```python
import sys
import pyathena
import pandas as pd

from IPython import display

conn = pyathena.connect(s3_staging_dir="s3://athena-output-351869726285/", region_name='us-east-1', encryption_option='SSE_S3')
```

---

# Query the 1000 Genomes Project Dataset

It's usually helpful to picture the data before performing any analysis, so we are going to import the database as the first step, and then view a few random rows.

<div class="alert alert-block alert-info">
<b>Tip:</b> If you are new to Jupyter Notebook, try run the code cell below using keyboard shortcut: "Shift" + "Enter". You could look up more keyboard shortcuts by pressing "H".
</div>


```python
# pd.set_option('display.max_colwidth', None) # This code expands the table horizontally so that all table cells are visible.
pd.read_sql('SELECT * FROM default.g1000vcf_csv_int LIMIT 10', conn).head(10)
```

We will go through a few important columns for you to understand in the table above. If you want to dive deeper, you may find this document helpful:https://samtools.github.io/hts-specs/VCFv4.2.pdf

1. `chrm` is the chromosome location.
2. `start_position` is the start position of the DNA variant.
3. `end_position` is the end position of the DNA variant.
4. `reference_base` is the allele at a specific location in the reference genomes, which are considered as the "approximated normal".
5. `alternate_base` is the allele that shows up at a specific location in the sample, but does not exist at the corresponding location in the reference genomes.
6. `rsid` is the identification number of the SNPs (Single Nucleotide Polymorphism), see here: https://www.snpedia.com/index.php/SNPedia
7. `qual` is the Phred-scaled quality score for the assertion made in ALT.
8. `filter` indicates PASS when the position has passed all filters.

The `info` column explains additional information about the data. In the example below, we can tell the percentage of East Asian population who has `rs559815820` SNP is only **0.0399361%**. That is 2 alleles count out of 5008 alleles from 2504 samples in the ***1000 Genomes Project***. The percentage of East Asian population who has `rs559815820` SNP is **0.1%**.

<img src="materials/images/example-1.png"/>

Now, try to reference the explanation below about the additional information section. Select a SNP from the data you queried: How does this SNP present in the ethnic group you belong to, in comparison to other ethnic groups?

1. `AC` stands for allele count in genotypes.
2. `AF` stands for allele frequency. Note: allele frequency (AF) = allele count (AC)/ allele number (AN)
3. `AN` stands for total number of alleles in genotypes.
4. `NS` stands for number of samples with data.
5. `EAS_AF` is the allele frequency of **East Asian population**.
6. `AMR_AF` is the allele frequency of **Ad Mixed American population**.
7. `EUR_AF` is the allele frequency of **European population**.
8. `AFR_AF` is the allele frequency of **African population**.
9. `SAS_AF` is the allele frequency of **South Asian population**.

---

# Search for SNPs (Single Nucleotide Polymorphism)

Next, pick a SNPs you are interested in investigating from this website: https://www.snpedia.com

The website has SNPs associated with a wide range of phenotypes. You could start exploring the popular SNPs on the homepage.

Here are a few examples:
1. `rs53576` is highly associated with the ability to empathize with others.
2. `rs72921001` is responsible for certain population's dislike towards the taste of cilantro.
3. `rs28936679` is associated with sleep disorder.
4. `rs1805009` is associated with skin cancer.

The example code below calls the data in the ***1000 Genomes Project*** that has `rs12913832`, see the code  ` WHERE rsid='rs12913832'`.

`rs12913832` is the SNP associated with blue or brown eye color.


```python
pd.set_option('display.max_colwidth', None) # This code expands the table horizontally so that all table cells are visible.
pd.read_sql("SELECT * FROM \"default\".g1000vcf_csv_int WHERE rsid='rs12913832'", conn).head()
```

Does the result above make sense to you?

We can see allele frequency of `rs12913832` among **Ad Mixed American** and **European** populations are 20.17% and 63.62% respectively, while only 0.2% among **East Asian** and 2.8% among **African** populations.

---

# Query COSMIC68 Annotation Dataset (hg19)

***COSMIC*** database is short for Catalogue of Somatic Mutations in Cancer [hg19 cosmic 68]. Learn more at
https://cancer.sanger.ac.uk/cosmic.

Now, let's take a look at few random rows of the ***COSMIC*** database:


```python
pd.read_sql('SELECT * FROM "1000_genomes".hg19_cosmic68 LIMIT 10', conn).head(10)
```

Like previous datasets, we see the chromosome location, start and end positions in DNA sequence, reference allele and alternate allele.

The `cosmic_info` column tells us the types of cancer and mutation occurence.

---

# Query UCSC RefGene Annotation Dataset (hg19)

The NCBI RefSeq Genes composite track shows human protein-coding and non-protein-coding genes taken from the NCBI RNA reference sequences collection (RefSeq) [hg19 refGene]. You could learn more at the following:
- https://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/
- https://genome.ucsc.edu/cgi-bin/hgTables?db=hg19&hgta_group=genes&hgta_track=refSeqComposite&hgta_table=refGene&hgta_doSchema=describe+table+schema (Schema for NCBI RefSeq - RefSeq genes from NCBI)

<img src="materials/images/genetic-mutations.png"/>

Now, let's take a look at random rows of the dataset:


```python
pd.read_sql('SELECT * FROM "1000_genomes".hg19_ucsc_refgene_int LIMIT 10', conn).head(10)
```

Here are what some of the columns mean according to Schema for NCBI RefSeq - RefSeq genes from NCBI:

1. `cdsstart`: Coding region start.
2. `cdsend`: Coding region end.
3. `exoncount`: Number of exons.
4. `strand`: + or - for strand


---

# Variant-Based Annotation

Variant-based annotation aims to look for **exact matches** between a query variant and a record in annotation datasets (i.e., two items have identical  chromosome, start position, end position, reference allele and alternative allele).

The code below uses the `JOIN` function to look for exact matches between the ***1000 Genomes Project*** dataset and the ***COSMIC*** dataset. It compares the start and end positions, reference allele, alternate allele at chromosome 2 between the two datasets.


```python
pd.read_sql("SELECT A.chrm, A.start_position, A.end_position, A.reference_bases, A.alternate_bases,B.cosmic_info, A.info "
+ " FROM (SELECT * FROM \"default\".g1000vcf_csv_int WHERE chrm='2') as A "
+ " JOIN "
+ " (SELECT * FROM \"1000_genomes\".hg19_cosmic68_int WHERE chrm='2') as B "
+ " ON A.start_position=B.start_position AND A.alternate_bases=B.alternate_bases "
+ " ORDER By  A.start_position", conn).head()
```

The returned table tells us information about genes associated with different types of cancers.

In the example below, the cancer related to central nervous system is at chromosome 2, where the reference allele is C and alternate allele is T. It appears to only happen to **East Asian** population with a frequency of 1.69% while absent among other populations.

<img src="materials/images/example-2.png"/>

---

# Interval-Based Annotation [TP53: chrm17]

The aim of interval-based annotation is to look for overlap of a query variant with a region (this region could be a single position) in annotation databases.

<img src="materials/images/overlapping-condition.png"/>

<div class="alert alert-block alert-info">
<b> Green.start_position[X]<=Blue.end_position[B] AND Blue.start_position[A]<=Green.end_position[Y]</b>

Source: https://stackoverflow.com/questions/20981783/how-to-sum-overlapping-interval-parts-in-sql.
</div>

The code below uses the `JOIN` function to compare two datasets using overlapping condition.

1. Rather than comparing whether the two datasets are exactly the same, this method focuses on overlapping regions. In the following example, we are running chromosome 17 from the ***1000 Genomes Project*** against gene TP53 from chromosome 17 in the ***UCSC RefGene Annotation*** dataset.
2. The `ON` condition is trying to see if there is any overlapping between the two datasets in the selected region. See the graph below:



```python
pd.read_sql("SELECT A.chrm, A.start_position, A.end_position, A.reference_bases, A.alternate_bases,B.name, B.name2, A.info "
+ " FROM (SELECT * FROM \"default\".g1000vcf_csv_int WHERE chrm='17') as A "
+ " JOIN "
+ " (SELECT * FROM \"1000_genomes\".hg19_ucsc_refgene_int WHERE chrm='17' and name2='TP53') as B "
+ " ON A.start_position<=B.end_position AND B.start_position<=A.end_position "
+ " ORDER By  A.start_position", conn).head()
```

<img src="materials/images/conductor.png"/>

For your information, TP53 gene provides instructions for making a protein called tumor protein p53 (or p53). This protein acts as a tumor suppressor, which means that it regulates cell division by keeping cells from growing and dividing (proliferating) too fast, or in an uncontrolled way. TP53 acts like a conductor in an orchestra.

The table gives us a view of all the mutations on gene TP53 on chromosome 2. Typically, mutations on TP53 are considered as rare mutations.

---

# Contributions & acknowledgement

Thank the following team for working on this module:

- **Module Content**: Amir Bahmani
- **Engineering**: Amit Dixit
- **UX/UI Design & Illustration**: Kexin Cha
- **Video Production**: Francesca Goncalves
- **Project Management**: Amir Bahmani, Kexin Cha

---

Copyright (c) 2022 Stanford Data Ocean (SDO)

All rights reserved.
