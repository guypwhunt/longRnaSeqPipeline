# Nanopore
Useful scripts for Nanopore + makings of a pipeline

## RNAseq
You will need to copy the files into your project directory (containing 'fastq_pass', 'fastq_fail' etc. folders).

- rnaseq.nf is the main Nextflow file. Here you will find all the processes. Feel free to fork and add more, if needed! 
- barcodes.txt is an array with barcodes. It tells 'start_nextflow.sh' how to concatenate the barcoded read files. 
- start_nextflow.sh is a bash script to launch Nextflow. It first concatenates all the files into one fastq.gz file and then start the pipeline
- resume_nextflow can be used to simply re-launch Nextflow with all cached analyses (from where you left off)
- nextflow.config is a file used to input configurations into the pipeline. It's currently empty so that the pipeline can remain generic.
- samples.tsv is a file needed for downstream anlayses (see below). It is used to differentiate between groups/conditions. 

### Note: RNAseq pipeline needs to be finished with Ballgown / DEseq2 / EDGER 
