#!/bin/bash -l
#$ -N scRNAseq-velocyto
#$ -j y
#$ -m a
#$ -M yah2014@med.cornell.edu
#$ -l h_vmem=60G
#$ -l athena=true
#$ -q *@@red

spack load -r bzip2@1.0.6
spack load -r xz@5.2.3

#---------------------Variables to be set-------------------------#
PROJECT_NAME="scRNAseq-Glioma"
path=/athena/elementolab/scratch/yah2014/Projects/${PROJECT_NAME}/data/bam
file_folder=$(ls ${path} | tail -n +${SGE_TASK_ID}| head -1) # Uses job array for each sample in the folder
file="${file_folder}.bam" # add .bam
rmsk_gtf=/athena/elementolab/scratch/yah2014/Indexed_genome/hg19_rmsk.gtf
genes_gtf=/athena/elementolab/scratch/yah2014/Indexed_genome/hg19-and-mm10/refdata-cellranger-hg19_and_mm10-1.2.0/genes/genes.gtf
echo "path="
echo "$path"
echo " "
echo $(ll -l $path/$file_folder/$file)
echo $(ll -l $rmsk_gtf)
echo $(ll -l $genes_gtf)

#----------------Files Transfer---------------------------#

cd $TMPDIR
echo "Start to transfer bam file."
rsync -v -a -z --exclude 'Summary' $path/$file_folder $TMPDIR
echo "total size"
echo $(ll -l $TMPDIR)
echo "files transferring accomplished."
echo " "

#----------------rename BAM File-------------------
echo "to sort by cellID"
mv $TMPDIR/$file_folder/$file $TMPDIR/$file_folder/possorted_genome_bam.bam
echo $(ll -l $TMPDIR/$file_folder/possorted_genome_bam.bam)
echo Pipestance completed successfully! > $TMPDIR/$file_folder/_log

#-----------velocyto Command--------------------------------#
echo "Processing velocyto run10x"
echo " "
echo "-------------------------------- "
echo "Processing $file_folder"
velocyto run10x -m $rmsk_gtf $TMPDIR/$file_folder $genes_gtf
echo "velocyto run10x Complished"
echo "velocyto output files:"
echo $(ls -l $TMPDIR/$file_folder/velocyto/)
echo " "

#---------------------------------------------------------------
rsync -rav $TMPDIR/$file_folder/velocyto $path/$file_folder
