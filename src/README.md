# Analysis scripts
Note that scripts are ordered based on the analysis workflow.

## Alignment based on 10x CellRanger
### 10x Alignments
02_Cellranger.sh

### Script to generate genome with RFP protein
02_make_custom_10x_ref.sh

### Aggregate samples into one
03_CellrangerAggregate.sh


## Statistical analysis in R
### Sets up R environment (called by each script)
00-init.R

### External datasets required prior to analysis
1. GEO datasets

EXT_01_GetGEO.sh

EXT_02_GEOAnalysis.R

2. Human Cell Atlas

EXT_03_HCA.R

3. Organ mixture experiments from 10x Genomics

EXT_04_OrganismMixtureExperiments.sh

### Get 10x sequencing statistics
04_SequencingStats.R

### Process reference spike in cells
05_PrepSpikeIns_Clean.R

### Processes in-sample spike-ins
06_SpikeIns_Islets.R

### Clean up data
07_03_CleanHuman_CleanSpikeIns.R

07_04_CleanHuman3_CleanSpikeIns.R

07_11_CleanMouse_CleanSpikeIns.R

### Seurat analyses
10_H_01_Seurat.R

10_H3_01_Seurat.R

10_M_01_Seurat.R

10_SCRIPT_01_Seurat.R

### Cell type assignment based on marker genes
11_H3_AssignCelltypes.R

11_H_AssignCelltypes.R

11_M_AssignCelltypes.R

### Script to identify differential genes (called later)
12_DiffGenes_SCRIPT.R

### Machine-learning-based annotation of cell types
14_02_Classifier_moreCelltypes_NoEndocrine.R

### Some simple analyses on the cells
15_Counts.R

### Identification of differential genes (based on script 12)
16_02_H3_DiffGenes_NoInsulin_FindMore.R

16_02_H_DiffGenes_NoInsulin_FindMore.R

16_02_M_DiffGenes_NoInsulin_FindMore.R

### Summarizing of differential genes
17_02_06_GenePlots_ByGroup_q0.1_noGABAII.R

### Identification of alpha and beta cell signatures
20_01_AlphaCellSignature2.R

20_02_AlphaCellSignature.R

### Calculation of correlation to alpha and beta cell signatures
21_01_CellCorrelation.R

### Additional clustering experiments
30_01_ClusteringExperiments.R

### Analysis of drop-seq dataset
DropSeq

### Export gene averages across cells
EXPORT_01_Averages.R

### Figures related to contamination
FIG_01_Contamination.R

### Demonstration that data look better after correction
FIG_01_Correlate_Replicates.R

### Analyses of single-cell dataset (clustering, celltypes,...)
FIG_02_SingleCellAnalyses.R

### Differential expression
FIG_03_DiffExpr.R

### Additional, detailed analyses of spike-in cells, contamination
FIG_04_SpikeInDetails.R

## R Functions
FUNC_Enrichr.R

FUNC_Seurat2_fast_simple.R

FUNC_Seurat_subcluster.R

FUNCTIONS_HELPERS.R

FUNCTIONS.R

FUNCTIONS_Synonyms_etc.R