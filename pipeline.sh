#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=100G
#SBATCH -t 48:00:00
#SBATCH --ntasks=25
#SBATCH --partition cpu
#SBATCH --mem-per-cpu=200G
#SBATCH --cpus-per-task=1

conda activate longRnaScan
module load nextflow

cd /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/rnaDevelopment/Nanopore/RNA/nextflow_RNAseq/

#mkdir -p fastq_split
#mkdir -p 1-7
#mkdir -p 8-16
#mkdir -p 17-23

#while read p; do
#  echo "$p"
#  cat /scratch/prj/sgdp_nanopore/Projects/BR22_00024_ALS_Liverpool_RNA_seq/1-7/20230504_1318_1C_PAK60199_b6ada2b6/fastq_pass/$p/* > 1-7/$p.gz
#done < barcodes.txt
#
#while read p; do
#  echo "$p"
#  cat /scratch/prj/sgdp_nanopore/Projects/BR22_00024_ALS_Liverpool_RNA_seq/8-16/20230509_1236_1A_PAM38314_badd002d/fastq_pass/$p/* > 8-16/$p.gz
#done < barcodes.txt
#
#while read p; do
#  echo "$p"
#  cat /scratch/prj/sgdp_nanopore/Projects/BR22_00024_ALS_Liverpool_RNA_seq/17-23/20230517_1635_1C_PAM35116_e82c66d2/fastq_pass/$p/* > 17-23/$p.gz
#done < barcodes.txt

#while read p; do
#  echo "$p"
#  cat /scratch/prj/sgdp_nanopore/Projects/BR22_00024_ALS_Liverpool_RNA_seq/1-7/20230504_1318_1C_PAK60199_b6ada2b6/fastq_pass/$p/* \
#  /scratch/prj/sgdp_nanopore/Projects/BR22_00024_ALS_Liverpool_RNA_seq/8-16/20230509_1236_1A_PAM38314_badd002d/fastq_pass/$p/* \
#  /scratch/prj/sgdp_nanopore/Projects/BR22_00024_ALS_Liverpool_RNA_seq/17-23/20230517_1635_1C_PAM35116_e82c66d2/fastq_pass/$p/* \
# > fastq_split/$p.gz
#done < barcodes.txt

#while read p; do
#  echo "$p"
#  cat /scratch/prj/sgdp_nanopore/Projects/BR22_00024_ALS_Liverpool_RNA_seq/1-7/20230504_1318_1C_PAK60199_b6ada2b6/fastq_pass/$p/* \
#  > fastq_split/$p.gz
#done < barcodes.txt

nextflow run rnaseq.nf -resume