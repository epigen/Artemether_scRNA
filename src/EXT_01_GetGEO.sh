cd $PROCESSED/artemether/analysis/
mkdir EXT_01_GEO
cd EXT_01_GEO


geo=GSE73727
geo=GSE81547
geo=GSE84133

array=("GSE73727" "GSE81547" "GSE84133" "GSE85241" "GSE83139" "GSE81608")
array=("GSE85241")
for geo in ${array[@]}
do            
  geo_short=${geo:0:5}

  echo $geo
  echo $geo_short


  files=$(curl -l ftp://ftp.ncbi.nlm.nih.gov/geo/series/${geo_short}nnn/$geo/suppl/)
  for ff in ${files[@]}
  do
    echo $ff
    curl -O ftp://ftp.ncbi.nlm.nih.gov/geo/series/${geo_short}nnn/$geo/suppl/$ff
  done
  
  files=$(curl -l ftp://ftp.ncbi.nlm.nih.gov/geo/series/${geo_short}nnn/$geo/matrix/)
  for ff in ${files[@]}
  do
    echo $ff
    curl -O curl -O ftp://ftp.ncbi.nlm.nih.gov/geo/series/${geo_short}nnn/$geo/matrix/$ff
  done

  mkdir $geo
  files=$(ls ${geo}*tar)
  for ff in ${files[@]}
  do
    echo $ff
    tar -xvf $ff -C $geo/
  done
done



curl -O https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-5061/E-MTAB-5061.processed.1.zip
curl -O https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-5061/E-MTAB-5061.sdrf.txt

unzip E-MTAB-5061.processed.1.zip
mv pancreas_refseq_rpkms_counts_3514sc.txt E-MTAB-5061.pancreas_refseq_rpkms_counts_3514sc.txt