cd $PROCESSED/artemether/analysis/EXT_04_SpeciesMixtures/

wget http://cf.10xgenomics.com/samples/cell-exp/3.0.0/hgmm_1k_v3/hgmm_1k_v3_filtered_feature_bc_matrix.tar.gz
wget http://cf.10xgenomics.com/samples/cell-exp/3.0.0/hgmm_1k_v2/hgmm_1k_v2_filtered_feature_bc_matrix.tar.gz
wget http://cf.10xgenomics.com/samples/cell-exp/3.0.0/hgmm_5k_v3/hgmm_5k_v3_filtered_feature_bc_matrix.tar.gz
wget http://cf.10xgenomics.com/samples/cell-exp/3.0.0/hgmm_10k_v3/hgmm_10k_v3_filtered_feature_bc_matrix.tar.gz

tar -xvzf hgmm_1k_v3_filtered_feature_bc_matrix.tar.gz
mv filtered_feature_bc_matrix hgmm_1k_v3_filtered_feature_bc_matrix
cd hgmm_1k_v3_filtered_feature_bc_matrix
gunzip barcodes.tsv.gz
gunzip features.tsv.gz
gunzip matrix.mtx.gz
mv features.tsv genes.tsv
cd ..

tar -xvzf hgmm_1k_v2_filtered_feature_bc_matrix.tar.gz
mv filtered_feature_bc_matrix hgmm_1k_v2_filtered_feature_bc_matrix
cd hgmm_1k_v2_filtered_feature_bc_matrix
gunzip barcodes.tsv.gz
gunzip features.tsv.gz
gunzip matrix.mtx.gz
mv features.tsv genes.tsv
cd ..

tar -xvzf hgmm_5k_v3_filtered_feature_bc_matrix.tar.gz
mv filtered_feature_bc_matrix hgmm_5k_v3_filtered_feature_bc_matrix
cd hgmm_5k_v3_filtered_feature_bc_matrix
gunzip barcodes.tsv.gz
gunzip features.tsv.gz
gunzip matrix.mtx.gz
mv features.tsv genes.tsv
cd ..

tar -xvzf hgmm_10k_v3_filtered_feature_bc_matrix.tar.gz
mv filtered_feature_bc_matrix hgmm_10k_v3_filtered_feature_bc_matrix
cd hgmm_10k_v3_filtered_feature_bc_matrix
gunzip barcodes.tsv.gz
gunzip features.tsv.gz
gunzip matrix.mtx.gz
mv features.tsv genes.tsv
cd ..