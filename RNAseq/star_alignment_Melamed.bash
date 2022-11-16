#!/bin/bash
#SBATCH -c 16                               # Request four cores
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 0-3:00                        # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=50G                          # Memory total in MB (for all cores)
#SBATCH -o /home/ah390/logs/%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e /home/ah390/logs/%j.err                 # File to which STDERR will be written, including job ID

module load gcc/6.2.0
module load star/2.7.3a
module load samtools/1.13

cd /n/scratch3/users/a/ah390/melamed2019
mkdir star_out

for var in SRR8144907_1.fastq.bz2 SRR8144908_1.fastq.bz2 SRR8144909_1.fastq.bz2 SRR8144910_1.fastq.bz2 SRR8144911_1.fastq.bz2 SRR8144912_1.fastq.bz2   
do
	filename1=${var}
	filename3=${var:0:10}
	mkdir star_out/$filename3
	
	STAR --runThreadN 16 \
	--genomeDir /home/ah390/hg38_ref/assembly200 \
	--readFilesIn reads/$filename1 \
	--readFilesCommand bunzip2 -c \
	--outFileNamePrefix star_out/$filename3/$filename3 \
	--quantMode GeneCounts \
	--outSAMtype BAM SortedByCoordinate \
	--limitBAMsortRAM 10000000000 \
	--genomeLoad LoadAndKeep
	
	cd star_out/$filename3
	samtools index ${var:0:10}Aligned.sortedByCoord.out.bam
	cd ../..
done
