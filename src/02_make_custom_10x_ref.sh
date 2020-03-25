#!/bin/bash
#
# execute:
# $ ./make_custom_10x_ref.sh <NAME-OF-GENOME-IN-RESOURCE-FOLDER> <OUTPUT-FOLDER> [<FASTA-FILE> <FASTA-FILE> ..ls]
# (each FASTA file *must* contain exactly (and only) one sequence entry including a header line)
#
# see also: https://github.com/epigen/crop-seq/blob/master/src/guides_to_ref.py

GENOME=$1
OUTROOT=$2
CELLRANGER=/home/nfortelny/code/cellranger-2.1.0/cellranger

GENOME_ROOT_10X=/scratch/lab_bock/nfortelny/resources/10X_Genomics/



mkdir -p ${OUTROOT}

# gather and filter genome sequence...
if [ ! -f ${OUTROOT}/${GENOME}-ext.fa ] ; then
	#SRC=${RESOURCES}/genomes/${GENOME}/${GENOME}.fa
	SRC=${GENOME_ROOT_10X}/${GENOME}/fasta/genome.fa
	echo "Get reference from $SRC"
	cp ${SRC} ${OUTROOT}/${GENOME}-ext.fa
fi
# ... and gene annotation:
if [ ! -f ${OUTROOT}/${GENOME}-ext.gtf ] ; then
	#SRC=${RESOURCES}/genomes/${GENOME}/${GENOME}.gtf
	SRC=${GENOME_ROOT_10X}/${GENOME}/genes/genes.gtf
	echo "Get annotations from $SRC"
	# N.B. this uses only protein-coding genes
	#${CELLRANGER} mkgtf ${SRC} ${OUTROOT}/${GENOME}-ext.gtf --attribute=gene_biotype:protein_coding
	cp ${SRC} ${OUTROOT}/${GENOME}-ext.gtf
fi

# append the sequence from each provided FASTA file to the genome (as a "chromosome") 
# and add a mock annotation entry for a "gene" spanning the entire pseudo-chromosome
for FASTA in "${@:3}"
do
	if [ -f $FASTA ] ; then
		NAME=$(basename $FASTA)
		NAME=${NAME/.fa/}
		echo "$NAME ($FASTA)"
		echo "assuming one sequence per file (EVERYTHING ELSE WILL CRASH AND BURN!)"
		
		echo "  * add to FASTA"
		echo ">${NAME}_chrom" >> ${OUTROOT}/${GENOME}-ext.fa
		tail -n +2 $FASTA >> ${OUTROOT}/${GENOME}-ext.fa

		echo "  * add to GTF"
		LEN=$(tail -n +2 $FASTA | wc -m)
		echo "${NAME}_chrom	ext	gene	1	${LEN}	.	+	.	gene_id \"${NAME}_gene\"; gene_name \"${NAME}_gene\"; gene_source \"ext\"; gene_biotype \"lincRNA\";
${NAME}_chrom	ext	transcript	1	${LEN}	.	+	.	gene_id \"${NAME}_gene\"; transcript_id \"${NAME}_transcript\"; gene_name \"${NAME}_gene\"; gene_source \"ext\"; gene_biotype \"lincRNA\"; transcript_name \"${NAME}_transcript\"; transcript_source \"ext\";
${NAME}_chrom	ext	exon	1	${LEN}	.	+	.	gene_id \"${NAME}_gene\"; transcript_id \"${NAME}_transcript\"; exon_number \"1\"; gene_name \"${NAME}_gene\"; gene_source \"ext\"; gene_biotype \"lincRNA\"; transcript_name \"${NAME}_transcript\"; transcript_source \"ext\"; exon_id \"${NAME}_exon\"" >> ${OUTROOT}/${GENOME}-ext.gtf
	else
		echo "File not found: $FASTA"
	fi
done

# then build the 10x-compatible reference:
cd ${OUTROOT}
${CELLRANGER} mkref --genome=${GENOME}_ext --fasta=${GENOME}-ext.fa --genes=${GENOME}-ext.gtf
