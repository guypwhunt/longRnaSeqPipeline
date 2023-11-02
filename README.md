# RNA Seq Nanopore Pipeline
A nextflow pipeline enabling the conversion of nanopore fastqs to a transcript count matrix including several quality control steps.

## Pipeline
You will need to copy the files into your project directory (containing 'fastq_pass', 'fastq_fail' etc. folders).

- rnaseq.nf is the main Nextflow file. Here you will find all the processes.
- barcodes.txt is an array with barcodes. It tells 'start_nextflow.sh' how to concatenate the barcoded read files. 
- start_nextflow.sh is a bash script to launch Nextflow. It first concatenates all the files into one fastq.gz file and then start the pipeline
- resume_nextflow.sh can be used to simply re-launch Nextflow with all cached analyses (from where you left off)
- nextflow.config is a file used to input configurations into the pipeline. It's currently empty so that the pipeline can remain generic.