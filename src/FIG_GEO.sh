cd /home/nfortelny/projects_shared/pathway_learning/results_analysis/GEO/FTP/

ncftp
set passive on
set so-bufsize 33554432
open ftp://geo:33%259uyj_fCh%3FM16H@ftp-private.ncbi.nlm.nih.gov
mkdir christoph_bock_artemether/
cd christoph_bock_artemether/

# PROCESSED DATA
put -z /scratch/lab_bock/shared/projects/artemether/analysis/GEO/FTP_processed//Human1_2_CellAnnotation_final.tsv Human1_2_CellAnnotation_final.tsv
put -z /scratch/lab_bock/shared/projects/artemether/analysis/GEO/FTP_processed//Human1_2_correctedTPM_Matrix_barcodes.txt Human1_2_correctedTPM_Matrix_barcodes.txt
put -z /scratch/lab_bock/shared/projects/artemether/analysis/GEO/FTP_processed//Human1_2_correctedTPM_Matrix_genes.txt Human1_2_correctedTPM_Matrix_genes.txt
put -z /scratch/lab_bock/shared/projects/artemether/analysis/GEO/FTP_processed//Human1_2_correctedTPM_Matrix.mtx Human1_2_correctedTPM_Matrix.mtx
put -z /scratch/lab_bock/shared/projects/artemether/analysis/GEO/FTP_processed//Human3_CellAnnotation_final.tsv Human3_CellAnnotation_final.tsv
put -z /scratch/lab_bock/shared/projects/artemether/analysis/GEO/FTP_processed//Human3_correctedTPM_Matrix_barcodes.txt Human3_correctedTPM_Matrix_barcodes.txt
put -z /scratch/lab_bock/shared/projects/artemether/analysis/GEO/FTP_processed//Human3_correctedTPM_Matrix_genes.txt Human3_correctedTPM_Matrix_genes.txt
put -z /scratch/lab_bock/shared/projects/artemether/analysis/GEO/FTP_processed//Human3_correctedTPM_Matrix.mtx Human3_correctedTPM_Matrix.mtx
put -z /scratch/lab_bock/shared/projects/artemether/analysis/GEO/FTP_processed//MF179_CellAnnotation_final.tsv MF179_CellAnnotation_final.tsv
put -z /scratch/lab_bock/shared/projects/artemether/analysis/GEO/FTP_processed//MouseLTI_CellAnnotation_final.tsv MouseLTI_CellAnnotation_final.tsv
put -z /scratch/lab_bock/shared/projects/artemether/analysis/GEO/FTP_processed//MouseLTI_correctedTPM_Matrix_barcodes.txt MouseLTI_correctedTPM_Matrix_barcodes.txt
put -z /scratch/lab_bock/shared/projects/artemether/analysis/GEO/FTP_processed//MouseLTI_correctedTPM_Matrix_genes.txt MouseLTI_correctedTPM_Matrix_genes.txt
put -z /scratch/lab_bock/shared/projects/artemether/analysis/GEO/FTP_processed//MouseLTI_correctedTPM_Matrix.mtx MouseLTI_correctedTPM_Matrix.mtx
put -z /scratch/lab_bock/shared/projects/artemether/results_pipeline/cellranger_count/MF179_hmG/outs/raw_gene_bc_matrices_h5.h5 MF179_raw_gene_bc_matrices_h5.h5


# RAW DATA
put -z /data/groups/lab_bsf/sequences/BSF_0409_HNCYFBBXX_l1_l2_l3_l4_l5_l6_l7_10x/fastq_path/HNCYFBBXX//hIslets_I_A10_1_S33933/hIslets_I_A10_1_S33933_S18_L006_I1_001.fastq.gz hIslets_I_A10_1_S33933_S18_L006_I1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0409_HNCYFBBXX_l1_l2_l3_l4_l5_l6_l7_10x/fastq_path/HNCYFBBXX//hIslets_I_A10_1_S33933/hIslets_I_A10_1_S33933_S18_L006_R1_001.fastq.gz hIslets_I_A10_1_S33933_S18_L006_R1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0409_HNCYFBBXX_l1_l2_l3_l4_l5_l6_l7_10x/fastq_path/HNCYFBBXX//hIslets_I_A10_1_S33933/hIslets_I_A10_1_S33933_S18_L006_R2_001.fastq.gz hIslets_I_A10_1_S33933_S18_L006_R2_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0409_HNCYFBBXX_l1_l2_l3_l4_l5_l6_l7_10x/fastq_path/HNCYFBBXX//hIslets_I_DMSO_1_S33929/hIslets_I_DMSO_1_S33929_S17_L006_I1_001.fastq.gz hIslets_I_DMSO_1_S33929_S17_L006_I1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0409_HNCYFBBXX_l1_l2_l3_l4_l5_l6_l7_10x/fastq_path/HNCYFBBXX//hIslets_I_DMSO_1_S33929/hIslets_I_DMSO_1_S33929_S17_L006_R1_001.fastq.gz hIslets_I_DMSO_1_S33929_S17_L006_R1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0409_HNCYFBBXX_l1_l2_l3_l4_l5_l6_l7_10x/fastq_path/HNCYFBBXX//hIslets_I_DMSO_1_S33929/hIslets_I_DMSO_1_S33929_S17_L006_R2_001.fastq.gz hIslets_I_DMSO_1_S33929_S17_L006_R2_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0409_HNCYFBBXX_l1_l2_l3_l4_l5_l6_l7_10x/fastq_path/HNCYFBBXX//mIslets_I_A1/mIslets_I_A1_S3_L001_I1_001.fastq.gz mIslets_I_A1_S3_L001_I1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0409_HNCYFBBXX_l1_l2_l3_l4_l5_l6_l7_10x/fastq_path/HNCYFBBXX//mIslets_I_A1/mIslets_I_A1_S3_L001_R1_001.fastq.gz mIslets_I_A1_S3_L001_R1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0409_HNCYFBBXX_l1_l2_l3_l4_l5_l6_l7_10x/fastq_path/HNCYFBBXX//mIslets_I_A1/mIslets_I_A1_S3_L001_R2_001.fastq.gz mIslets_I_A1_S3_L001_R2_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0409_HNCYFBBXX_l1_l2_l3_l4_l5_l6_l7_10x/fastq_path/HNCYFBBXX//mIslets_I_A10/mIslets_I_A10_S4_L002_I1_001.fastq.gz mIslets_I_A10_S4_L002_I1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0409_HNCYFBBXX_l1_l2_l3_l4_l5_l6_l7_10x/fastq_path/HNCYFBBXX//mIslets_I_A10/mIslets_I_A10_S4_L002_R1_001.fastq.gz mIslets_I_A10_S4_L002_R1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0409_HNCYFBBXX_l1_l2_l3_l4_l5_l6_l7_10x/fastq_path/HNCYFBBXX//mIslets_I_A10/mIslets_I_A10_S4_L002_R2_001.fastq.gz mIslets_I_A10_S4_L002_R2_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0409_HNCYFBBXX_l1_l2_l3_l4_l5_l6_l7_10x/fastq_path/HNCYFBBXX//mIslets_I_DMSO/mIslets_I_DMSO_S2_L001_I1_001.fastq.gz mIslets_I_DMSO_S2_L001_I1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0409_HNCYFBBXX_l1_l2_l3_l4_l5_l6_l7_10x/fastq_path/HNCYFBBXX//mIslets_I_DMSO/mIslets_I_DMSO_S2_L001_R1_001.fastq.gz mIslets_I_DMSO_S2_L001_R1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0409_HNCYFBBXX_l1_l2_l3_l4_l5_l6_l7_10x/fastq_path/HNCYFBBXX//mIslets_I_DMSO/mIslets_I_DMSO_S2_L001_R2_001.fastq.gz mIslets_I_DMSO_S2_L001_R2_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0409_HNCYFBBXX_l1_l2_l3_l4_l5_l6_l7_10x/fastq_path/HNCYFBBXX//mIslets_I_FoxO/mIslets_I_FoxO_S6_L002_I1_001.fastq.gz mIslets_I_FoxO_S6_L002_I1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0409_HNCYFBBXX_l1_l2_l3_l4_l5_l6_l7_10x/fastq_path/HNCYFBBXX//mIslets_I_FoxO/mIslets_I_FoxO_S6_L002_R1_001.fastq.gz mIslets_I_FoxO_S6_L002_R1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0409_HNCYFBBXX_l1_l2_l3_l4_l5_l6_l7_10x/fastq_path/HNCYFBBXX//mIslets_I_FoxO/mIslets_I_FoxO_S6_L002_R2_001.fastq.gz mIslets_I_FoxO_S6_L002_R2_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0409_HNCYFBBXX_l1_l2_l3_l4_l5_l6_l7_10x/fastq_path/HNCYFBBXX//mIslets_I_GABA/mIslets_I_GABA_S5_L002_I1_001.fastq.gz mIslets_I_GABA_S5_L002_I1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0409_HNCYFBBXX_l1_l2_l3_l4_l5_l6_l7_10x/fastq_path/HNCYFBBXX//mIslets_I_GABA/mIslets_I_GABA_S5_L002_R1_001.fastq.gz mIslets_I_GABA_S5_L002_R1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0409_HNCYFBBXX_l1_l2_l3_l4_l5_l6_l7_10x/fastq_path/HNCYFBBXX//mIslets_I_GABA/mIslets_I_GABA_S5_L002_R2_001.fastq.gz mIslets_I_GABA_S5_L002_R2_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0409_HNCYFBBXX_l1_l2_l3_l4_l5_l6_l7_10x/fastq_path/HNCYFBBXX//mIslets_II_A1/mIslets_II_A1_S8_L003_I1_001.fastq.gz mIslets_II_A1_S8_L003_I1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0409_HNCYFBBXX_l1_l2_l3_l4_l5_l6_l7_10x/fastq_path/HNCYFBBXX//mIslets_II_A1/mIslets_II_A1_S8_L003_R1_001.fastq.gz mIslets_II_A1_S8_L003_R1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0409_HNCYFBBXX_l1_l2_l3_l4_l5_l6_l7_10x/fastq_path/HNCYFBBXX//mIslets_II_A1/mIslets_II_A1_S8_L003_R2_001.fastq.gz mIslets_II_A1_S8_L003_R2_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0409_HNCYFBBXX_l1_l2_l3_l4_l5_l6_l7_10x/fastq_path/HNCYFBBXX//mIslets_II_A10/mIslets_II_A10_S9_L003_I1_001.fastq.gz mIslets_II_A10_S9_L003_I1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0409_HNCYFBBXX_l1_l2_l3_l4_l5_l6_l7_10x/fastq_path/HNCYFBBXX//mIslets_II_A10/mIslets_II_A10_S9_L003_R1_001.fastq.gz mIslets_II_A10_S9_L003_R1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0409_HNCYFBBXX_l1_l2_l3_l4_l5_l6_l7_10x/fastq_path/HNCYFBBXX//mIslets_II_A10/mIslets_II_A10_S9_L003_R2_001.fastq.gz mIslets_II_A10_S9_L003_R2_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0409_HNCYFBBXX_l1_l2_l3_l4_l5_l6_l7_10x/fastq_path/HNCYFBBXX//mIslets_II_DMSO/mIslets_II_DMSO_S7_L003_I1_001.fastq.gz mIslets_II_DMSO_S7_L003_I1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0409_HNCYFBBXX_l1_l2_l3_l4_l5_l6_l7_10x/fastq_path/HNCYFBBXX//mIslets_II_DMSO/mIslets_II_DMSO_S7_L003_R1_001.fastq.gz mIslets_II_DMSO_S7_L003_R1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0409_HNCYFBBXX_l1_l2_l3_l4_l5_l6_l7_10x/fastq_path/HNCYFBBXX//mIslets_II_DMSO/mIslets_II_DMSO_S7_L003_R2_001.fastq.gz mIslets_II_DMSO_S7_L003_R2_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0409_HNCYFBBXX_l1_l2_l3_l4_l5_l6_l7_10x/fastq_path/HNCYFBBXX//mIslets_II_FoxO_1_S33905/mIslets_II_FoxO_1_S33905_S11_L004_I1_001.fastq.gz mIslets_II_FoxO_1_S33905_S11_L004_I1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0409_HNCYFBBXX_l1_l2_l3_l4_l5_l6_l7_10x/fastq_path/HNCYFBBXX//mIslets_II_FoxO_1_S33905/mIslets_II_FoxO_1_S33905_S11_L004_R1_001.fastq.gz mIslets_II_FoxO_1_S33905_S11_L004_R1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0409_HNCYFBBXX_l1_l2_l3_l4_l5_l6_l7_10x/fastq_path/HNCYFBBXX//mIslets_II_FoxO_1_S33905/mIslets_II_FoxO_1_S33905_S11_L004_R2_001.fastq.gz mIslets_II_FoxO_1_S33905_S11_L004_R2_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0409_HNCYFBBXX_l1_l2_l3_l4_l5_l6_l7_10x/fastq_path/HNCYFBBXX//mIslets_II_GABA_1_S33901/mIslets_II_GABA_1_S33901_S10_L004_I1_001.fastq.gz mIslets_II_GABA_1_S33901_S10_L004_I1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0409_HNCYFBBXX_l1_l2_l3_l4_l5_l6_l7_10x/fastq_path/HNCYFBBXX//mIslets_II_GABA_1_S33901/mIslets_II_GABA_1_S33901_S10_L004_R1_001.fastq.gz mIslets_II_GABA_1_S33901_S10_L004_R1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0409_HNCYFBBXX_l1_l2_l3_l4_l5_l6_l7_10x/fastq_path/HNCYFBBXX//mIslets_II_GABA_1_S33901/mIslets_II_GABA_1_S33901_S10_L004_R2_001.fastq.gz mIslets_II_GABA_1_S33901_S10_L004_R2_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0409_HNCYFBBXX_l1_l2_l3_l4_l5_l6_l7_10x/fastq_path/HNCYFBBXX//mIslets_III_A1/mIslets_III_A1_S13_L005_I1_001.fastq.gz mIslets_III_A1_S13_L005_I1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0409_HNCYFBBXX_l1_l2_l3_l4_l5_l6_l7_10x/fastq_path/HNCYFBBXX//mIslets_III_A1/mIslets_III_A1_S13_L005_R1_001.fastq.gz mIslets_III_A1_S13_L005_R1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0409_HNCYFBBXX_l1_l2_l3_l4_l5_l6_l7_10x/fastq_path/HNCYFBBXX//mIslets_III_A1/mIslets_III_A1_S13_L005_R2_001.fastq.gz mIslets_III_A1_S13_L005_R2_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0409_HNCYFBBXX_l1_l2_l3_l4_l5_l6_l7_10x/fastq_path/HNCYFBBXX//mIslets_III_A10/mIslets_III_A10_S14_L005_I1_001.fastq.gz mIslets_III_A10_S14_L005_I1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0409_HNCYFBBXX_l1_l2_l3_l4_l5_l6_l7_10x/fastq_path/HNCYFBBXX//mIslets_III_A10/mIslets_III_A10_S14_L005_R1_001.fastq.gz mIslets_III_A10_S14_L005_R1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0409_HNCYFBBXX_l1_l2_l3_l4_l5_l6_l7_10x/fastq_path/HNCYFBBXX//mIslets_III_A10/mIslets_III_A10_S14_L005_R2_001.fastq.gz mIslets_III_A10_S14_L005_R2_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0409_HNCYFBBXX_l1_l2_l3_l4_l5_l6_l7_10x/fastq_path/HNCYFBBXX//mIslets_III_DMSO_1_S33909/mIslets_III_DMSO_1_S33909_S12_L004_I1_001.fastq.gz mIslets_III_DMSO_1_S33909_S12_L004_I1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0409_HNCYFBBXX_l1_l2_l3_l4_l5_l6_l7_10x/fastq_path/HNCYFBBXX//mIslets_III_DMSO_1_S33909/mIslets_III_DMSO_1_S33909_S12_L004_R1_001.fastq.gz mIslets_III_DMSO_1_S33909_S12_L004_R1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0409_HNCYFBBXX_l1_l2_l3_l4_l5_l6_l7_10x/fastq_path/HNCYFBBXX//mIslets_III_DMSO_1_S33909/mIslets_III_DMSO_1_S33909_S12_L004_R2_001.fastq.gz mIslets_III_DMSO_1_S33909_S12_L004_R2_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0409_HNCYFBBXX_l1_l2_l3_l4_l5_l6_l7_10x/fastq_path/HNCYFBBXX//mIslets_III_FoxO_1_S33925/mIslets_III_FoxO_1_S33925_S16_L006_I1_001.fastq.gz mIslets_III_FoxO_1_S33925_S16_L006_I1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0409_HNCYFBBXX_l1_l2_l3_l4_l5_l6_l7_10x/fastq_path/HNCYFBBXX//mIslets_III_FoxO_1_S33925/mIslets_III_FoxO_1_S33925_S16_L006_R1_001.fastq.gz mIslets_III_FoxO_1_S33925_S16_L006_R1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0409_HNCYFBBXX_l1_l2_l3_l4_l5_l6_l7_10x/fastq_path/HNCYFBBXX//mIslets_III_FoxO_1_S33925/mIslets_III_FoxO_1_S33925_S16_L006_R2_001.fastq.gz mIslets_III_FoxO_1_S33925_S16_L006_R2_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0409_HNCYFBBXX_l1_l2_l3_l4_l5_l6_l7_10x/fastq_path/HNCYFBBXX//mIslets_III_GABA/mIslets_III_GABA_S15_L005_I1_001.fastq.gz mIslets_III_GABA_S15_L005_I1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0409_HNCYFBBXX_l1_l2_l3_l4_l5_l6_l7_10x/fastq_path/HNCYFBBXX//mIslets_III_GABA/mIslets_III_GABA_S15_L005_R1_001.fastq.gz mIslets_III_GABA_S15_L005_R1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0409_HNCYFBBXX_l1_l2_l3_l4_l5_l6_l7_10x/fastq_path/HNCYFBBXX//mIslets_III_GABA/mIslets_III_GABA_S15_L005_R2_001.fastq.gz mIslets_III_GABA_S15_L005_R2_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0410_HNF7TBBXX_l1_10x/fastq_path/HNF7TBBXX//hIslets_I_FoxO_1_S33941/hIslets_I_FoxO_1_S33941_S2_L001_I1_001.fastq.gz hIslets_I_FoxO_1_S33941_S2_L001_I1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0410_HNF7TBBXX_l1_10x/fastq_path/HNF7TBBXX//hIslets_I_FoxO_1_S33941/hIslets_I_FoxO_1_S33941_S2_L001_R1_001.fastq.gz hIslets_I_FoxO_1_S33941_S2_L001_R1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0410_HNF7TBBXX_l1_10x/fastq_path/HNF7TBBXX//hIslets_I_FoxO_1_S33941/hIslets_I_FoxO_1_S33941_S2_L001_R2_001.fastq.gz hIslets_I_FoxO_1_S33941_S2_L001_R2_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0410_HNF7TBBXX_l1_10x/fastq_path/HNF7TBBXX//hIslets_I_GABA_1_S33937/hIslets_I_GABA_1_S33937_S1_L001_I1_001.fastq.gz hIslets_I_GABA_1_S33937_S1_L001_I1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0410_HNF7TBBXX_l1_10x/fastq_path/HNF7TBBXX//hIslets_I_GABA_1_S33937/hIslets_I_GABA_1_S33937_S1_L001_R1_001.fastq.gz hIslets_I_GABA_1_S33937_S1_L001_R1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0410_HNF7TBBXX_l1_10x/fastq_path/HNF7TBBXX//hIslets_I_GABA_1_S33937/hIslets_I_GABA_1_S33937_S1_L001_R2_001.fastq.gz hIslets_I_GABA_1_S33937_S1_L001_R2_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0429_HNTVGBBXX_l5_l6_l7_l8_10x/fastq_path/HNTVGBBXX//hIslets_II_A1/hIslets_II_A1_S8_L007_I1_001.fastq.gz hIslets_II_A1_S8_L007_I1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0429_HNTVGBBXX_l5_l6_l7_l8_10x/fastq_path/HNTVGBBXX//hIslets_II_A1/hIslets_II_A1_S8_L007_R1_001.fastq.gz hIslets_II_A1_S8_L007_R1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0429_HNTVGBBXX_l5_l6_l7_l8_10x/fastq_path/HNTVGBBXX//hIslets_II_A1/hIslets_II_A1_S8_L007_R2_001.fastq.gz hIslets_II_A1_S8_L007_R2_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0429_HNTVGBBXX_l5_l6_l7_l8_10x/fastq_path/HNTVGBBXX//hIslets_II_A2/hIslets_II_A2_S9_L007_I1_001.fastq.gz hIslets_II_A2_S9_L007_I1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0429_HNTVGBBXX_l5_l6_l7_l8_10x/fastq_path/HNTVGBBXX//hIslets_II_A2/hIslets_II_A2_S9_L007_R1_001.fastq.gz hIslets_II_A2_S9_L007_R1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0429_HNTVGBBXX_l5_l6_l7_l8_10x/fastq_path/HNTVGBBXX//hIslets_II_A2/hIslets_II_A2_S9_L007_R2_001.fastq.gz hIslets_II_A2_S9_L007_R2_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0429_HNTVGBBXX_l5_l6_l7_l8_10x/fastq_path/HNTVGBBXX//hIslets_II_DMSO/hIslets_II_DMSO_S10_L007_I1_001.fastq.gz hIslets_II_DMSO_S10_L007_I1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0429_HNTVGBBXX_l5_l6_l7_l8_10x/fastq_path/HNTVGBBXX//hIslets_II_DMSO/hIslets_II_DMSO_S10_L007_R1_001.fastq.gz hIslets_II_DMSO_S10_L007_R1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0429_HNTVGBBXX_l5_l6_l7_l8_10x/fastq_path/HNTVGBBXX//hIslets_II_DMSO/hIslets_II_DMSO_S10_L007_R2_001.fastq.gz hIslets_II_DMSO_S10_L007_R2_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0429_HNTVGBBXX_l5_l6_l7_l8_10x/fastq_path/HNTVGBBXX//hIslets_II_FoxO/hIslets_II_FoxO_S12_L008_I1_001.fastq.gz hIslets_II_FoxO_S12_L008_I1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0429_HNTVGBBXX_l5_l6_l7_l8_10x/fastq_path/HNTVGBBXX//hIslets_II_FoxO/hIslets_II_FoxO_S12_L008_R1_001.fastq.gz hIslets_II_FoxO_S12_L008_R1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0429_HNTVGBBXX_l5_l6_l7_l8_10x/fastq_path/HNTVGBBXX//hIslets_II_FoxO/hIslets_II_FoxO_S12_L008_R2_001.fastq.gz hIslets_II_FoxO_S12_L008_R2_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0429_HNTVGBBXX_l5_l6_l7_l8_10x/fastq_path/HNTVGBBXX//hIslets_II_GABA/hIslets_II_GABA_S13_L008_I1_001.fastq.gz hIslets_II_GABA_S13_L008_I1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0429_HNTVGBBXX_l5_l6_l7_l8_10x/fastq_path/HNTVGBBXX//hIslets_II_GABA/hIslets_II_GABA_S13_L008_R1_001.fastq.gz hIslets_II_GABA_S13_L008_R1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0429_HNTVGBBXX_l5_l6_l7_l8_10x/fastq_path/HNTVGBBXX//hIslets_II_GABA/hIslets_II_GABA_S13_L008_R2_001.fastq.gz hIslets_II_GABA_S13_L008_R2_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0460_HV235BBXX_l1_l7_l8_10x/fastq_path/HV235BBXX//MF179/MF179_S2_L007_I1_001.fastq.gz MF179_S2_L007_I1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0460_HV235BBXX_l1_l7_l8_10x/fastq_path/HV235BBXX//MF179/MF179_S2_L007_R1_001.fastq.gz MF179_S2_L007_R1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0460_HV235BBXX_l1_l7_l8_10x/fastq_path/HV235BBXX//MF179/MF179_S2_L007_R2_001.fastq.gz MF179_S2_L007_R2_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0460_HV235BBXX_l1_l7_l8_10x/fastq_path/HV235BBXX//MF179/MF179_S2_L008_I1_001.fastq.gz MF179_S2_L008_I1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0460_HV235BBXX_l1_l7_l8_10x/fastq_path/HV235BBXX//MF179/MF179_S2_L008_R1_001.fastq.gz MF179_S2_L008_R1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0460_HV235BBXX_l1_l7_l8_10x/fastq_path/HV235BBXX//MF179/MF179_S2_L008_R2_001.fastq.gz MF179_S2_L008_R2_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0525_H2K53BBXY_l2_l3_l4_l5_l6_l7_l8_10x/fastq_path/H2K53BBXY//hIslets_III_A1_36hr/hIslets_III_A1_36hr_S4_L004_I1_001.fastq.gz hIslets_III_A1_36hr_S4_L004_I1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0525_H2K53BBXY_l2_l3_l4_l5_l6_l7_l8_10x/fastq_path/H2K53BBXY//hIslets_III_A1_36hr/hIslets_III_A1_36hr_S4_L004_R1_001.fastq.gz hIslets_III_A1_36hr_S4_L004_R1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0525_H2K53BBXY_l2_l3_l4_l5_l6_l7_l8_10x/fastq_path/H2K53BBXY//hIslets_III_A1_36hr/hIslets_III_A1_36hr_S4_L004_R2_001.fastq.gz hIslets_III_A1_36hr_S4_L004_R2_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0525_H2K53BBXY_l2_l3_l4_l5_l6_l7_l8_10x/fastq_path/H2K53BBXY//hIslets_III_A1_72hr/hIslets_III_A1_72hr_S7_L005_I1_001.fastq.gz hIslets_III_A1_72hr_S7_L005_I1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0525_H2K53BBXY_l2_l3_l4_l5_l6_l7_l8_10x/fastq_path/H2K53BBXY//hIslets_III_A1_72hr/hIslets_III_A1_72hr_S7_L005_R1_001.fastq.gz hIslets_III_A1_72hr_S7_L005_R1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0525_H2K53BBXY_l2_l3_l4_l5_l6_l7_l8_10x/fastq_path/H2K53BBXY//hIslets_III_A1_72hr/hIslets_III_A1_72hr_S7_L005_R2_001.fastq.gz hIslets_III_A1_72hr_S7_L005_R2_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0525_H2K53BBXY_l2_l3_l4_l5_l6_l7_l8_10x/fastq_path/H2K53BBXY//hIslets_III_A2_36hr/hIslets_III_A2_36hr_S8_L005_I1_001.fastq.gz hIslets_III_A2_36hr_S8_L005_I1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0525_H2K53BBXY_l2_l3_l4_l5_l6_l7_l8_10x/fastq_path/H2K53BBXY//hIslets_III_A2_36hr/hIslets_III_A2_36hr_S8_L005_R1_001.fastq.gz hIslets_III_A2_36hr_S8_L005_R1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0525_H2K53BBXY_l2_l3_l4_l5_l6_l7_l8_10x/fastq_path/H2K53BBXY//hIslets_III_A2_36hr/hIslets_III_A2_36hr_S8_L005_R2_001.fastq.gz hIslets_III_A2_36hr_S8_L005_R2_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0525_H2K53BBXY_l2_l3_l4_l5_l6_l7_l8_10x/fastq_path/H2K53BBXY//hIslets_III_A2_72hr/hIslets_III_A2_72hr_S10_L006_I1_001.fastq.gz hIslets_III_A2_72hr_S10_L006_I1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0525_H2K53BBXY_l2_l3_l4_l5_l6_l7_l8_10x/fastq_path/H2K53BBXY//hIslets_III_A2_72hr/hIslets_III_A2_72hr_S10_L006_R1_001.fastq.gz hIslets_III_A2_72hr_S10_L006_R1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0525_H2K53BBXY_l2_l3_l4_l5_l6_l7_l8_10x/fastq_path/H2K53BBXY//hIslets_III_A2_72hr/hIslets_III_A2_72hr_S10_L006_R2_001.fastq.gz hIslets_III_A2_72hr_S10_L006_R2_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0525_H2K53BBXY_l2_l3_l4_l5_l6_l7_l8_10x/fastq_path/H2K53BBXY//hIslets_III_A3_72hr/hIslets_III_A3_72hr_S13_L007_I1_001.fastq.gz hIslets_III_A3_72hr_S13_L007_I1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0525_H2K53BBXY_l2_l3_l4_l5_l6_l7_l8_10x/fastq_path/H2K53BBXY//hIslets_III_A3_72hr/hIslets_III_A3_72hr_S13_L007_R1_001.fastq.gz hIslets_III_A3_72hr_S13_L007_R1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0525_H2K53BBXY_l2_l3_l4_l5_l6_l7_l8_10x/fastq_path/H2K53BBXY//hIslets_III_A3_72hr/hIslets_III_A3_72hr_S13_L007_R2_001.fastq.gz hIslets_III_A3_72hr_S13_L007_R2_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0525_H2K53BBXY_l2_l3_l4_l5_l6_l7_l8_10x/fastq_path/H2K53BBXY//hIslets_III_DMSO1_36hr/hIslets_III_DMSO1_36hr_S5_L004_I1_001.fastq.gz hIslets_III_DMSO1_36hr_S5_L004_I1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0525_H2K53BBXY_l2_l3_l4_l5_l6_l7_l8_10x/fastq_path/H2K53BBXY//hIslets_III_DMSO1_36hr/hIslets_III_DMSO1_36hr_S5_L004_R1_001.fastq.gz hIslets_III_DMSO1_36hr_S5_L004_R1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0525_H2K53BBXY_l2_l3_l4_l5_l6_l7_l8_10x/fastq_path/H2K53BBXY//hIslets_III_DMSO1_36hr/hIslets_III_DMSO1_36hr_S5_L004_R2_001.fastq.gz hIslets_III_DMSO1_36hr_S5_L004_R2_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0525_H2K53BBXY_l2_l3_l4_l5_l6_l7_l8_10x/fastq_path/H2K53BBXY//hIslets_III_DMSO1_72hr/hIslets_III_DMSO1_72hr_S9_L005_I1_001.fastq.gz hIslets_III_DMSO1_72hr_S9_L005_I1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0525_H2K53BBXY_l2_l3_l4_l5_l6_l7_l8_10x/fastq_path/H2K53BBXY//hIslets_III_DMSO1_72hr/hIslets_III_DMSO1_72hr_S9_L005_R1_001.fastq.gz hIslets_III_DMSO1_72hr_S9_L005_R1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0525_H2K53BBXY_l2_l3_l4_l5_l6_l7_l8_10x/fastq_path/H2K53BBXY//hIslets_III_DMSO1_72hr/hIslets_III_DMSO1_72hr_S9_L005_R2_001.fastq.gz hIslets_III_DMSO1_72hr_S9_L005_R2_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0525_H2K53BBXY_l2_l3_l4_l5_l6_l7_l8_10x/fastq_path/H2K53BBXY//hIslets_III_DMSO2_36hr/hIslets_III_DMSO2_36hr_S6_L004_I1_001.fastq.gz hIslets_III_DMSO2_36hr_S6_L004_I1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0525_H2K53BBXY_l2_l3_l4_l5_l6_l7_l8_10x/fastq_path/H2K53BBXY//hIslets_III_DMSO2_36hr/hIslets_III_DMSO2_36hr_S6_L004_R1_001.fastq.gz hIslets_III_DMSO2_36hr_S6_L004_R1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0525_H2K53BBXY_l2_l3_l4_l5_l6_l7_l8_10x/fastq_path/H2K53BBXY//hIslets_III_DMSO2_36hr/hIslets_III_DMSO2_36hr_S6_L004_R2_001.fastq.gz hIslets_III_DMSO2_36hr_S6_L004_R2_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0525_H2K53BBXY_l2_l3_l4_l5_l6_l7_l8_10x/fastq_path/H2K53BBXY//hIslets_III_DMSO2_72hr/hIslets_III_DMSO2_72hr_S11_L006_I1_001.fastq.gz hIslets_III_DMSO2_72hr_S11_L006_I1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0525_H2K53BBXY_l2_l3_l4_l5_l6_l7_l8_10x/fastq_path/H2K53BBXY//hIslets_III_DMSO2_72hr/hIslets_III_DMSO2_72hr_S11_L006_R1_001.fastq.gz hIslets_III_DMSO2_72hr_S11_L006_R1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0525_H2K53BBXY_l2_l3_l4_l5_l6_l7_l8_10x/fastq_path/H2K53BBXY//hIslets_III_DMSO2_72hr/hIslets_III_DMSO2_72hr_S11_L006_R2_001.fastq.gz hIslets_III_DMSO2_72hr_S11_L006_R2_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0525_H2K53BBXY_l2_l3_l4_l5_l6_l7_l8_10x/fastq_path/H2K53BBXY//hIslets_III_DMSO3_72hr/hIslets_III_DMSO3_72hr_S12_L006_I1_001.fastq.gz hIslets_III_DMSO3_72hr_S12_L006_I1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0525_H2K53BBXY_l2_l3_l4_l5_l6_l7_l8_10x/fastq_path/H2K53BBXY//hIslets_III_DMSO3_72hr/hIslets_III_DMSO3_72hr_S12_L006_R1_001.fastq.gz hIslets_III_DMSO3_72hr_S12_L006_R1_001.fastq.gz
put -z /data/groups/lab_bsf/sequences/BSF_0525_H2K53BBXY_l2_l3_l4_l5_l6_l7_l8_10x/fastq_path/H2K53BBXY//hIslets_III_DMSO3_72hr/hIslets_III_DMSO3_72hr_S12_L006_R2_001.fastq.gz hIslets_III_DMSO3_72hr_S12_L006_R2_001.fastq.gz
