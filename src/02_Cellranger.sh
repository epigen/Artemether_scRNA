cr="$CODEBASE/cellranger-2.1.0/cellranger"
mouseGenome=/home/nfortelny/resources_nfortelny/10X_Genomics/refdata-cellranger-mm10-1.2.0/
humanGenome=/home/nfortelny/resources_nfortelny/10X_Genomics/refdata-cellranger-GRCh38-1.2.0/
combinedGenome=/home/nfortelny/resources_nfortelny/10X_Genomics/refdata-cellranger-hg19_and_mm10-1.2.0/
rfpGenome=$PROCESSED/artemether/custom_genomes/RFP_Genome/refdata-cellranger-mm10-1.2.0_ext/

outPath=$PROCESSED/artemether/results_pipeline/cellranger_count/



###
### DATA FROM BSF 0409 (Jan 24 2018)
###
inPath=/data/groups/lab_bsf/sequences/BSF_0409_HNCYFBBXX_l1_l2_l3_l4_l5_l6_l7_10x/fastq_path/HNCYFBBXX/
mkdir $outPath
cd $outPath

# MURINE DATA (No lineage tracing)
array=("mIslets_I_A1" "mIslets_I_A10" "mIslets_I_DMSO" "mIslets_I_FoxO" "mIslets_I_GABA" "mIslets_II_A1" "mIslets_II_A10" "mIslets_II_DMSO" "mIslets_II_FoxO_1_S33905" "mIslets_II_GABA_1_S33901" "mIslets_III_A1" "mIslets_III_A10" "mIslets_III_DMSO_1_S33909" "mIslets_III_FoxO_1_S33925" "mIslets_III_GABA")
for file in ${array[@]}
do            
    echo $file
    sbatch --job-name="10X $file" --cpus-per-task=16 --ntasks=1 --mem=90000 --partition=longq --time=5-00:00:00 \
        --wrap="$cr count --id=$file --transcriptome=$rfpGenome --fastqs=$inPath/${file} --expect-cells=10000 --localcores=16 --localmem=90" \
            --output="$outPath/$file.log"
done

# LINEAGE TRACING
array=("mIslets_LT1_FoxO_1_S33961")
for file in ${array[@]}
do            
    echo $file
    sbatch --job-name="10X $file" --cpus-per-task=16 --ntasks=1 --mem=90000 --partition=longq --time=5-00:00:00 \
        --wrap="$cr count --id=$file --transcriptome=$rfpGenome --fastqs=$inPath/${file} --expect-cells=10000 --localcores=16 --localmem=90" \
            --output="$outPath/$file.log"
done

# HUMAN DATA
array=("hIslets_I_A10_1_S33933" "hIslets_I_DMSO_1_S33929")
for file in ${array[@]}
do            
    echo $file
    sbatch --job-name="10X $file" --cpus-per-task=16 --ntasks=1 --mem=90000 --partition=longq --time=5-00:00:00 \
        --wrap="$cr count --id=$file --transcriptome=$humanGenome --fastqs=$inPath/${file} --expect-cells=10000 --localcores=16 --localmem=90" \
            --output="$outPath/$file.log"
done



###
### DATA FROM BSF 0410 (Jan 26 2018)
###
inPath=/data/groups/lab_bsf/sequences/BSF_0410_HNF7TBBXX_l1_10x/fastq_path/HNF7TBBXX/
mkdir $outPath
cd $outPath

# HUMAN DATA
array=("hIslets_I_FoxO_1_S33941" "hIslets_I_GABA_1_S33937" "hIslets_I_pla2g16_1_S33945")
for file in ${array[@]}
do            
    echo $file
    sbatch --job-name="10X $file" --cpus-per-task=16 --ntasks=1 --mem=90000 --partition=longq --time=5-00:00:00 \
        --wrap="$cr count --id=$file --transcriptome=$humanGenome --fastqs=$inPath/${file} --expect-cells=10000 --localcores=16 --localmem=90" \
            --output="$outPath/$file.log"
done

###
### DATA FROM BSF 0412 (Feb 5 2018)
###
inPath=/data/groups/lab_bsf/sequences/BSF_0412_HNJH5BBXX_l1_10x/fastq_path/HNJH5BBXX/
mkdir $outPath
cd $outPath

# LTI data
array=("mIslets_LT1_A10" "mIslets_LT1_DMSO" "mIslets_LT1_GABA")
for file in ${array[@]}
do            
    echo $file
    sbatch --job-name="10X $file" --cpus-per-task=16 --ntasks=1 --mem=90000 --partition=longq --time=5-00:00:00 \
        --wrap="$cr count --id=$file --transcriptome=$rfpGenome --fastqs=$inPath/${file} --expect-cells=10000 --localcores=16 --localmem=90 --chemistry=SC3Pv2" \
            --output="$outPath/$file.log"
done

###
### DATA FROM BSF 0429 (March 19, 2018), Human islets
###
inPath=/data/groups/lab_bsf/sequences/BSF_0429_HNTVGBBXX_l5_l6_l7_l8_10x/fastq_path/HNTVGBBXX/
mkdir $outPath
cd $outPath
array=("hIslets_II_A1" "hIslets_II_A2" "hIslets_II_DMSO" "hIslets_II_FoxO" "hIslets_II_GABA")
for file in ${array[@]}
do   
  echo $file
  sbatch --job-name="10X $file" --cpus-per-task=16 --ntasks=1 --mem=90000 --partition=longq --time=5-00:00:00 \
      --wrap="$cr count --id=$file --transcriptome=$humanGenome --fastqs=$inPath/${file} --expect-cells=10000 --localcores=16 --localmem=90" \
          --output="${file}.log"
done






################
################ ALIGN TO MOUSE AND HUMAN
################

###
### DATA FROM BSF 0409 (Jan 24 2018)
###
inPath=/data/groups/lab_bsf/sequences/BSF_0409_HNCYFBBXX_l1_l2_l3_l4_l5_l6_l7_10x/fastq_path/HNCYFBBXX/
mkdir $outPath
cd $outPath

# MURINE DATA (No lineage tracing)
array=("mIslets_I_A1" "mIslets_I_A10" "mIslets_I_DMSO" "mIslets_I_FoxO" "mIslets_I_GABA" "mIslets_II_A1" "mIslets_II_A10" "mIslets_II_DMSO" "mIslets_II_FoxO_1_S33905" "mIslets_II_GABA_1_S33901" "mIslets_III_A1" "mIslets_III_A10" "mIslets_III_DMSO_1_S33909" "mIslets_III_FoxO_1_S33925" "mIslets_III_GABA")
for file in ${array[@]}
do            
    echo $file
    sbatch --job-name="10X $file" --cpus-per-task=16 --ntasks=1 --mem=90000 --partition=longq --time=5-00:00:00 \
        --wrap="$cr count --id=${file}_hmG --transcriptome=$combinedGenome --fastqs=$inPath/${file} --expect-cells=10000 --localcores=16 --localmem=90" \
            --output="$outPath/${file}_hmG.log"
done

# LINEAGE TRACING
array=("mIslets_LT1_FoxO_1_S33961")
for file in ${array[@]}
do            
    echo $file
    sbatch --job-name="10X $file" --cpus-per-task=16 --ntasks=1 --mem=90000 --partition=longq --time=5-00:00:00 \
        --wrap="$cr count --id=${file}_hmG --transcriptome=$combinedGenome --fastqs=$inPath/${file} --expect-cells=10000 --localcores=16 --localmem=90" \
            --output="$outPath/${file}_hmG.log"
done

# HUMAN DATA
array=("hIslets_I_A10_1_S33933" "hIslets_I_DMSO_1_S33929")
for file in ${array[@]}
do            
    echo $file
    sbatch --job-name="10X $file" --cpus-per-task=16 --ntasks=1 --mem=90000 --partition=longq --time=5-00:00:00 \
        --wrap="$cr count --id=${file}_hmG --transcriptome=$combinedGenome --fastqs=$inPath/${file} --expect-cells=10000 --localcores=16 --localmem=90" \
            --output="$outPath/${file}_hmG.log"
done



###
### DATA FROM BSF 0410 (Jan 26 2018)
###
inPath=/data/groups/lab_bsf/sequences/BSF_0410_HNF7TBBXX_l1_10x/fastq_path/HNF7TBBXX/
mkdir $outPath
cd $outPath

# HUMAN DATA
array=("hIslets_I_FoxO_1_S33941" "hIslets_I_GABA_1_S33937" "hIslets_I_pla2g16_1_S33945")
for file in ${array[@]}
do            
    echo $file
    sbatch --job-name="10X $file" --cpus-per-task=16 --ntasks=1 --mem=90000 --partition=longq --time=5-00:00:00 \
        --wrap="$cr count --id=${file}_hmG --transcriptome=$combinedGenome --fastqs=$inPath/${file} --expect-cells=10000 --localcores=16 --localmem=90" \
            --output="$outPath/${file}_hmG.log"
done

###
### DATA FROM BSF 0412 (Feb 5 2018)
###
inPath=/data/groups/lab_bsf/sequences/BSF_0412_HNJH5BBXX_l1_10x/fastq_path/HNJH5BBXX/
mkdir $outPath
cd $outPath

# LTI data
array=("mIslets_LT1_A10" "mIslets_LT1_DMSO" "mIslets_LT1_GABA")
for file in ${array[@]}
do            
    echo $file
    sbatch --job-name="10X $file" --cpus-per-task=16 --ntasks=1 --mem=90000 --partition=longq --time=5-00:00:00 \
        --wrap="$cr count --id=${file}_hmG --transcriptome=$combinedGenome --fastqs=$inPath/${file} --expect-cells=10000 --localcores=16 --localmem=90 --chemistry=SC3Pv2" \
            --output="$outPath/${file}_hmG.log"
done

###
### DATA FROM BSF 0429 (March 19, 2018), Human islets
###
inPath=/data/groups/lab_bsf/sequences/BSF_0429_HNTVGBBXX_l5_l6_l7_l8_10x/fastq_path/HNTVGBBXX/
mkdir $outPath
cd $outPath
array=("hIslets_II_A1" "hIslets_II_A2" "hIslets_II_DMSO" "hIslets_II_FoxO" "hIslets_II_GABA" "Hypothalamus_traject_E15" "Hypothalamus_traject_E17" "Hypothalamus_traject_P0" "Hypothalamus_traject_P10" "Hypothalamus_traject_P2" "Hypothalamus_traject_P23" "Hypothalamus_traject_P2SN")
for file in ${array[@]}
do   
  echo $file
  sbatch --job-name="10X $file" --cpus-per-task=16 --ntasks=1 --mem=90000 --partition=longq --time=5-00:00:00 \
      --wrap="$cr count --id=${file}_hmG --transcriptome=$combinedGenome --fastqs=$inPath/${file} --expect-cells=10000 --localcores=16 --localmem=90" \
          --output="${file}_hmG.log"
done
array=("Hypothalamus_traject_E15" "Hypothalamus_traject_E17" "Hypothalamus_traject_P0" "Hypothalamus_traject_P10" "Hypothalamus_traject_P2" "Hypothalamus_traject_P23" "Hypothalamus_traject_P2SN")
for file in ${array[@]}
do   
  echo $file
  sbatch --job-name="10X $file human" --cpus-per-task=16 --ntasks=1 --mem=90000 --partition=longq --time=5-00:00:00 \
      --wrap="$cr count --id=${file}_human --transcriptome=$humanGenome --fastqs=$inPath/${file} --expect-cells=10000 --localcores=16 --localmem=90" \
          --output="${file}_human.log"
done
array=("Hypothalamus_traject_E15" "Hypothalamus_traject_E17" "Hypothalamus_traject_P0" "Hypothalamus_traject_P10" "Hypothalamus_traject_P2" "Hypothalamus_traject_P23" "Hypothalamus_traject_P2SN")
for file in ${array[@]}
do   
  echo $file
  sbatch --job-name="10X $file mouse" --cpus-per-task=16 --ntasks=1 --mem=90000 --partition=longq --time=5-00:00:00 \
      --wrap="$cr count --id=${file}_mouse --transcriptome=$mouseGenome --fastqs=$inPath/${file} --expect-cells=10000 --localcores=16 --localmem=90" \
          --output="${file}_mouse.log"
done



###
### DATA FROM BSF 0429 (March 19, 2018), References
###
inPath=/data/groups/lab_bsf/sequences/BSF_0460_HV235BBXX_l1_l7_l8_10x/fastq_path/HV235BBXX/
mkdir $outPath
cd $outPath
array=("MF179")
for file in ${array[@]}
do   
  echo $file
  sbatch --job-name="10X $file" --cpus-per-task=16 --ntasks=1 --mem=90000 --partition=longq --time=5-00:00:00 \
      --wrap="$cr count --id=${file}_hmG --transcriptome=$combinedGenome --fastqs=$inPath/${file} --expect-cells=10000 --localcores=16 --localmem=90" \
          --output="${file}_hmG.log"
done
for file in ${array[@]}
do   
  echo $file
  sbatch --job-name="10X $file human" --cpus-per-task=16 --ntasks=1 --mem=90000 --partition=longq --time=5-00:00:00 \
      --wrap="$cr count --id=${file}_human --transcriptome=$humanGenome --fastqs=$inPath/${file} --expect-cells=10000 --localcores=16 --localmem=90" \
          --output="${file}_human.log"
done
for file in ${array[@]}
do   
  echo $file
  sbatch --job-name="10X $file mouse" --cpus-per-task=16 --ntasks=1 --mem=90000 --partition=longq --time=5-00:00:00 \
      --wrap="$cr count --id=${file}_mouse --transcriptome=$mouseGenome --fastqs=$inPath/${file} --expect-cells=10000 --localcores=16 --localmem=90" \
          --output="${file}_mouse.log"
done




###
### DATA FROM BSF_0525_H2K53BBXY (October 2, 2018), Human islets III
###
inPath=/data/groups/lab_bsf/sequences/BSF_0525_H2K53BBXY_l2_l3_l4_l5_l6_l7_l8_10x/fastq_path/H2K53BBXY/
mkdir $outPath
cd $outPath
array=("hIslets_III_A1_36hr" "hIslets_III_A2_72hr" "hIslets_III_DMSO1_72hr" "hIslets_III_DMSO3_72hr" "hIslets_III_A1_72hr" "hIslets_III_A3_72hr" "hIslets_III_DMSO2_36hr" "hIslets_III_A2_36hr" "hIslets_III_DMSO1_36hr" "hIslets_III_DMSO2_72hr")
for file in ${array[@]}
do   
  echo $file
  sbatch --job-name="10X $file combined" --cpus-per-task=16 --ntasks=1 --mem=90000 --partition=longq --time=5-00:00:00 \
      --wrap="$cr count --id=${file}_hmG --transcriptome=$combinedGenome --fastqs=$inPath/${file} --expect-cells=10000 --localcores=16 --localmem=90" \
          --output="${file}_hmG.log"
done
for file in ${array[@]}
do   
  echo $file
  sbatch --job-name="10X $file human" --cpus-per-task=16 --ntasks=1 --mem=90000 --partition=longq --time=5-00:00:00 \
      --wrap="$cr count --id=${file}_human --transcriptome=$humanGenome --fastqs=$inPath/${file} --expect-cells=10000 --localcores=16 --localmem=90" \
          --output="${file}_human.log"
done