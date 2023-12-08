# DatasetEnrichmentAnalysis

The idea is similar to Gene Set Enrichment Analysis (GSEA), which actually is **(Reference) Gene Set Enrichment (in input ranks) Analysis**. Our proposed method is **(Reference) Dataset Enrichment (in input genes) Analysis**.

We will include following reference data:

1.  Expression patterns:

```         
a.   Tissue level (Gene): GTex; Human Protein Atlas;

b.   Tissue or sub-tissue level (Protein): Human Protein Atlas. Need some review/modification for \"protein\" data.

c.   Cell level: Single Cell Atlas
```

2.  Perturbation-driven gene expression changes: [DATA LIBRARY [clue.io]](https://clue.io/data/CMap2020#LINCS2020). Need some data re-format. Almost ready to go.

3.  Disease related:

```         
a.   TCGA: Gene expression. Consider other types of data except RNA-Seq?

b.   Consider many tissues or disease specific bulk or single cell RNA-Seq database?
```

4.  Genetics: Can be used but data cleaning works needed;

```         
a.   GWAS: [GWAS Catalog (ebi.ac.uk)](https://www.ebi.ac.uk/gwas/docs/file-downloads) (Need some data cleaning/aggregation)

b.   PheWAS: <https://pheweb.sph.umich.edu/>   (Seems not including everything. Need more work on this)

c.   COSMIC: [ABL1 Gene  - Somatic Mutations in Cancer (sanger.ac.uk)](https://cancer.sanger.ac.uk/cosmic/gene/analysis?ln=ABL1); \"Tissue distribution\"
```

5.  Epigenetic:

    a.  Methylation:  [A DNA methylation atlas of normal human cell types \| Nature](https://www.nature.com/articles/s41586-022-05580-6). Issue is no Intensity provided (Supplementary Dataset 1), only region and Methylation/Not (Figure 1 as an example). May be methylation regions can be counted and methylation count for each protein can be used.

6.  Spatial Transcriptome data:

    a\. Need find a comprehensive dataset. And then a "Spatial pattern" summary statistic value will be generated for each gene in each sample to be tested.

7.  User input:

    1.  We can also download count table from Bulk Rna-Seq databases for similar analysis, for example, [recount3: uniformly processed RNA-seq](https://rna.recount.bio/). This may be an add-on for users to select any data as reference to test.
