#!/bin/bash -l
#$ -N scRNAseq-velocyto
#$ -j y
#$ -m a
#$ -M yah2014@med.cornell.edu
#$ -l h_vmem=60G
#$ -l athena=true
#$ -q *@@red

spack load -r python@3.6.0
spack load -r samtools@1.8
spack load -r zlib@1.2.11
spack load -r openssl@1

#---------------------Variables to be set-------------------------#
PROJECT_NAME="scRNAseq-Glioma"
path=/athena/elementolab/scratch/yah2014/Projects/${PROJECT_NAME}/data/bam
file_folder=$(ls ${path} | tail -n +${SGE_TASK_ID}| head -1) # Uses job array for each sample in the folder
file="${file_folder}.bam" # add .bam
rmsk_gtf=/athena/elementolab/scratch/yah2014/Indexed_genome/hg19_rmsk.gtf
genes_gtf=/athena/elementolab/scratch/yah2014/Indexed_genome/refdata-cellranger-hg19-3.0.0/genes/genes.gtf
echo "path="
echo "$path"
echo " "
echo $(ls -l $path/$file_folder/$file)
echo $(ls -l $rmsk_gtf)
echo $(ls -l $genes_gtf)

#----------------Files Transfer---------------------------#
cd $TMPDIR
echo "Start to transfer bam file."
rsync -v -a -z --exclude 'Summary' $path/$file_folder ./
echo "total size"
echo $(ls -l $TMPDIR/$file_folder)
echo "files transferring accomplished."
echo " "

#----------------rename BAM File-------------------
echo "to sort by cellID"
mv $TMPDIR/$file_folder/$file $TMPDIR/$file_folder/possorted_genome_bam.bam
echo $(ls -l $TMPDIR/$file_folder/possorted_genome_bam.bam)
echo Pipestance completed successfully! > $TMPDIR/$file_folder/_log

#-----------velocyto Command--------------------------------#
echo "Processing velocyto run10x"
echo " "
echo "-------------------------------- "
echo "Processing $file_folder"
velocyto run -b barcodes.tsv -o velocyto -e $file_folder -m $rmsk_gtf possorted_genome_bam.bam $genes_gtf
echo "velocyto run10x Complished"
echo "velocyto output files:"
echo $(ls -l $TMPDIR/$file_folder/velocyto/)
echo " "

#---------------------------------------------------------------
rsync -rav $TMPDIR/$file_folder/velocyto ${path%bam}
