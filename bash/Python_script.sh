#!/bin/bash -l
#$ -P mRNA velocyto
#$ -N python script
#$ -j y
#$ -m a
#$ -M yah2014@med.cornell.edu
#$ -l h_rt=4:00:00
#$ -l athena=true
#$ -l h_vmem=50G
#$ -q *@@red
spack load -r python@3.6.0
spack load -r samtools@1.8
spack load -r zlib@1.2.11
spack load -r openssl@1

samtools --version

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date : $(date)"
echo "Job name : $JOB_NAME"
echo "Job ID : $JOB_ID  $SGE_TASK_ID"
echo "=========================================================="

file=/athena/elementolab/scratch/yah2014/Projects/scRNAseq-Glioma/bash/velocyto_single.py
echo $(ls -l $file)
python $file $SGE_TASK_ID
