# Artemether_scRNA
Single cell RNA-seq analysis of pancreatic islets with drug treatment. This is the code used to analyze this dataset.

# Structure
Scripts to analyze all data are in src/, see the readme file in this directory for more details.
Additional data (i.e. list of marker genes, etc) are in the folder metadata/.


# System requirements
The code provided requires the setting of two environmental variables $CODEBASE (where the code is) and $PROCESSED (where the data is and analysis results will be stored). For example:

```
export CODEBASE="$HOME/code/"
export PROCESSED="$HOME/projects/"
```

Genome alignments are based on CellRanger (version 2.1.0).
Statistical analysis is based on R (version 3.4.0), using the following packages:
  - SoupX (version 0.1.1)
  - Seurat (version 2.0.1)
  - cowplot (version 0.8.0)
  - Rtsne (version 0.13)
  - ROCR (version 1.0-7)
  - gplots (version 3.0.1.1)
  - rhdf5 (version 2.20.0)
  - pryr (version 0.1.4)
  - project.init (version 0.0.1)
  - pheatmap (version 1.0.12)
  - glmnet (version 2.0-18)
  - enrichR (version 1.0)
  - dplyr (version 0.7.2)
  - doMC (version 1.3.4)
  - iterators (version 1.0.12)
  - foreach (version 1.4.3)
  - data.table (version 1.12.2)
  - cellrangerRkit (version 2.0.0)
  - Rmisc (version 1.5)
  - plyr (version 1.8.4)
  - lattice (version 0.20-38)
  - ggplot2 (version 3.1.1)
  - RColorBrewer (version 1.1-2)
  - Matrix (version 1.2-17)

