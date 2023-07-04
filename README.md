# MLIIC
Artificial intelligence learning landscape of triple-negative breast cancer uncovers new opportunities for enhancing outcomes and immunotherapy responses

## MLIIC signature construction

A MLIIC signature was determined by a comprehensive analysis of purified immune cells, TNBC cell lines, and TNBC solid tumor tissues using an innovative computational framework according to a series of sequential ML algorithms. 

(1)  The highest 15% expressed RNAs of each immune cell line were extracted as potential screening immune-related RNAs.

(2)  The tissue specificity index (TSI)  was employed to determine the specificity of the above candidate RNA expression per immune cell. The TSI varied from 0 to 1, where when the value was 0, the RNA was defined as immune cell-general RNA, and when the value was 1, the RNA was defined as immune cell-specific RNA. The highly expressed RNAs in all immune cell types were authenticated as immune-related generic RNAs (igRNAs).

(3)  Differentially expressed igRNAs that showed a pattern of both upregulations in multiple immune cell lines and downregulation in TNBC cell lines were determined as IIC-RNAs. 

(4)  ML algorithms for classification were further used for classification.  This step aimed to filter worthwhile IIC-RNAs by extracting the intersected IIC-RNAs identified by  ML algorithms for classification.

(5)  IIC-RNAs with prognostic potential were then screened in the TCGA TNBC dataset using univariate Cox regression analysis.

(6)  Next, ML algorithms for survival were subsequently applied to assess the significance of the prognostic IIC-RNAs and conduct the dimensionality reduction accordingly.

(7)  ML algorithms for scoring were used to determine the most reliable model.

(8)  The MLIIC signature was established according to the prognostic IIC-RNAs via performing the RSF algorithm. 

## Comparison of published signatures

Due to the rapid advancements in omics technologies, numerous studies have been reported to construct and analyze signatures based on specific gene combinations with promising predictive efficacy. We then aimed to systematically compare these relevant signatures with our MLIIC signature in the past decade. After a detailed investigation, we included a total of 148 signatures in terms of RNA signatures. It should be noted that the MLIIC signature exhibited superior performance with respect to C-index in the TCGA TNBC, METABRIC, GSE96058, and GSE103091 datasets, compared to nearly all of the previous models.
