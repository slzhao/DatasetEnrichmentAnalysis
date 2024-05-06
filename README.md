# DatasetEnrichmentAnalysis

The idea is similar to Gene Set Enrichment Analysis (GSEA), which actually is **(Reference) Gene Set Enrichment (in input ranks) Analysis**. Our proposed method is **(Reference) Dataset Enrichment (in input genes) Analysis**.

We will include following reference data:

1.  **Expression patterns**:

    log(nTPM+1) was performed on gene expression data. sqrt(sum(x\^2)) was performed on gene level for normalization and then ranked by sample level.

```         
a.   Tissue level (Gene): GTex; Human Protein Atlas;

b.   Tissue or sub-tissue level (Protein): Human Protein Atlas. Need some review/modification for \"protein\" data.

c.   Cell level: Single Cell Atlas (*TODO*)
```

2.  Perturbation-driven gene expression changes: [DATA LIBRARY [clue.io]](https://clue.io/data/CMap2020#LINCS2020).

    **2.1** Some introductions from here: <https://clue.io/connectopedia/data_levels>. We used level 5 data.

    Level 3a - NORM - Gene expression (GEX, Level 2) are normalized to invariant gene set curves and quantile normalized across each plate. Here, the data from each perturbagen treatment is referred to as a **profile**, **experiment**, or **instance**.

    Level 3b - INF- Additional values for 11,350 additional genes not directly measured in the L10000 assay are inferred based on the normalized values for the 978 landmark genes.

    Level 4 - ZS - Z-scores for each gene based on Level 3 with respect to the entire plate population. This comparison of profiles to their appropriate population control generates a list of differentially expressed genes.

    Level 5 - MODZ - replicate-collapsed z-score vectors based on Level 4. Replicate collapse generates one differential expression vector, which we term a **signature**. Connectivity analyses are performed on signatures.

    **2.2** Genes were ranked directly (as they were already Z score across samples).

3.  Disease related:

```         
a.   TCGA: Gene expression. Consider other types of data except RNA-Seq?

b.   Consider many tissues or disease specific bulk or single cell RNA-Seq database?
```

4.  Genetics: Can be used but data cleaning works needed;

For PheWAS, variants were filtered by p value (from website, seems p\<10\^-6) and beta value was used. Variant level data were summarized to gene level by which.max(abs(beta)). Then sqrt(sum(x\^2)) was performed on gene level for normalization and then ranked by sample level.

For GWAS, variants were filtered by p value (p\<10\^-8) and "OR.or.BETA" value was used (this is in fact OR as there is no value \<0). Variant level data were summarized to gene level by which.max(abs(log(OR.or.BETA))). Then sqrt(sum(x\^2)) was performed on gene level for normalization and then ranked by sample level.

```         
a.   GWAS: [GWAS Catalog (ebi.ac.uk)]<https://www.ebi.ac.uk/gwas/docs/file-downloads>

b.   PheWAS: <https://pheweb.sph.umich.edu/>

c.   COSMIC: [ABL1 Gene  - Somatic Mutations in Cancer (sanger.ac.uk)](https://cancer.sanger.ac.uk/cosmic/gene/analysis?ln=ABL1); \"Tissue distribution\"
```

5.  Epigenetic:

    a.  Methylation:  [A DNA methylation atlas of normal human cell types \| Nature](https://www.nature.com/articles/s41586-022-05580-6). Issue is no Intensity provided (Supplementary Dataset 1), only region and Methylation/Not (Figure 1 as an example). May be methylation regions can be counted and methylation count for each protein can be used.

6.  Spatial Transcriptome data:

    a\. Need find a comprehensive dataset. And then a "Spatial pattern" summary statistic value will be generated for each gene in each sample to be tested.

7.  User input:

    1.  We can also download count table from Bulk Rna-Seq databases for similar analysis, for example, [recount3: uniformly processed RNA-seq](https://rna.recount.bio/). This may be an add-on for users to select any data as reference to test.

## Usage

### Reference Data

| Name                              | Source         | Notes                             |
|-----------------------------|------------------|-------------------------|
| HPA_rna_tissue_consensus          | HPA            | Tissue level gene expression data |
| CMAPLINCS_CRISPR                  | CMAP           | CMAP knockout                     |
| CMAPLINCS_OverExpression          | CMAP           | CMAP OverExpression               |
| gwas_catalog_beta                 | gwas catalog   | gwas catalog beta value           |
| UKB.PheWAS.top_hits.PhenoSub.beta | PheWAS catalog | PheWAS catalog beta value         |

### Example

``` r
library(DatasetEnrichmentAnalysis)

selectedGenes=c("CARMIL1" ,   "HLA-B" ,     "PLCD1",      "WNT7B" ,     "HFE" ,     
"F5" ,        "ABCG8"  ,  "PITX2"  ,    "HLA-C" ,     "HLA-DQB1",   "HLA-DRB1",  
"TCF7L2"  , "IRF4"  ,     "AL359922.1","APOE"   ,    "HLA-DQA1" ,  "HLA-DRB5" ,  
"PTPN22" ,    "SFRP4"   , "CASZ1"  )

testResult=dataRankTest(selectedGenes,referenceData="HPA_rna_tissue_consensus",nRep=1000)
#testResult=dataRankTest(selectedGenes,referenceData="CMAPLINCS_CRISPR",nRep=1000)
#testResult=dataRankTest(selectedGenes,referenceData="gwas_catalog_beta",nRep=1000)
#testResult=dataRankTest(selectedGenes,referenceData="UKB.PheWAS.top_hits.PhenoSub.beta",nRep=1000)

library(tidyverse)
head(testResult %>% arrange(pAdj))
```
