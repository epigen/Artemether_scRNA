#!/bin/bash

cr="$CODEBASE/cellranger-2.1.0/cellranger"

basePW=$PROCESSED/artemether/results_pipeline/cellranger_count/
cd ${basePW}

sbatch --job-name="10X  aggr Mouse samples" --cpus-per-task=5 --ntasks=1 --mem=30000 --partition=mediumq --time=1-00:00:00 \
    --wrap="$cr aggr --id=Mouse --csv=$CODEBASE/artemether/metadata/Aggregate_Mouse.csv --normalize=none" \
    --output="Mouse.log"

sbatch --job-name="10X  aggr MouseLTI samples" --cpus-per-task=5 --ntasks=1 --mem=30000 --partition=mediumq --time=1-12:00:00 \
    --wrap="$cr aggr --id=MouseLTI --csv=$CODEBASE/artemether/metadata/Aggregate_Mouse_withLTI.csv --normalize=none" \
    --output="MouseLTI.log"

file=Human1
sbatch --job-name="10X  aggr $file" --cpus-per-task=5 --ntasks=1 --mem=30000 --partition=mediumq --time=1-00:00:00 \
    --wrap="$cr aggr --id=$file --csv=$CODEBASE/artemether/metadata/Aggregate_$file.csv --normalize=none" \
    --output="$file.log"

file=Human1_2
sbatch --job-name="10X  aggr $file" --cpus-per-task=5 --ntasks=1 --mem=30000 --partition=mediumq --time=1-00:00:00 \
    --wrap="$cr aggr --id=$file --csv=$CODEBASE/artemether/metadata/Aggregate_$file.csv --normalize=none" \
    --output="$file.log"

file=Human3
sbatch --job-name="10X  aggr $file" --cpus-per-task=5 --ntasks=1 --mem=30000 --partition=mediumq --time=1-00:00:00 \
    --wrap="$cr aggr --id=$file --csv=$CODEBASE/artemether/metadata/Aggregate_$file.csv --normalize=none" \
    --output="$file.log"