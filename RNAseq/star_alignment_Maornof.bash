#!/bin/bash
#SBATCH -c 16                               # Request four cores
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 0-4:00                        # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=50G                          # Memory total in MB (for all cores)
#SBATCH -o /home/ah390/logs/%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e /home/ah390/logs/%j.err                 # File to which STDERR will be written, including job ID

module load gcc/6.2.0
module load star/2.7.3a
module load samtools/1.13

cd /n/scratch3/users/a/ah390/maornof2020
mkdir star_out

for var in SRR13120267_1.fastq.bz2 SRR13120268_1.fastq.bz2 SRR13120269_1.fastq.bz2 SRR13120270_1.fastq.bz2 SRR13120271_1.fastq.bz2 SRR13120272_1.fastq.bz2 SRR13120275_1.fastq.bz2 SRR13120276_1.fastq.bz2 SRR13120277_1.fastq.bz2 SRR13120278_1.fastq.bz2 SRR13120279_1.fastq.bz2 SRR13120280_1.fastq.bz2 SRR13120283_1.fastq.bz2 SRR13120284_1.fastq.bz2 SRR13120285_1.fastq.bz2 SRR13120286_1.fastq.bz2 SRR13120287_1.fastq.bz2 SRR13120288_1.fastq.bz2
do
	filename1=${var}
	filename2=${var:0:11}_2.fastq.bz2
	filename3=${var:0:11}
	mkdir star_out/$filename3
	
	STAR --runThreadN 16 \
	--genomeDir /home/ah390/mm_genome/mm_assembly \
	--readFilesIn reads/$filename1 reads/$filename2 \
	--readFilesCommand bunzip2 -c \
	--outFileNamePrefix star_out/$filename3/$filename3 \
	--quantMode GeneCounts \
	--outSAMtype BAM SortedByCoordinate \
	--limitBAMsortRAM 10000000000 \
	--genomeLoad LoadAndKeep
	
	cd star_out/$filename3
	samtools index ${var:0:11}Aligned.sortedByCoord.out.bam
	cd ../..
done
