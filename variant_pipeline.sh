#!/bin/bash
if [ $# -lt 3 ]; then
	echo "USAGE:$0 <output work dir> <your input dir> <filename>"
	exit 1;
fi

filename = $2

STAR \
	--runThreadN 10 \
	--twopassMode Basic \
	--genomeDir /home/tchen/mouse_project/refer/star-ref \
	--readFilesIn "$2"/"$filename"_1.clean.fq "$2"/"$filename"_2.clean.fq \
	--outFilterType BySJout \
	--outFilterMultimapNmax 20 \
	--alignSJoverhangMin 8 \
	--alignSJDBoverhangMin 1 \
	--outFilterMismatchNmax 999 \
	--outFilterMismatchNoverReadLmax 0.04 \
	--alignIntronMin 20 \
	--alignIntronMax 1000000 \
	--alignMatesGapMax 1000000 \
	--outSAMstrandField intronMotif \
	--outSAMattrRGline ID:$filename SM:$filename PL:illumina LB:$filename PU:$filename \
	--outSAMtype BAM SortedByCoordinate \
	--outFileNamePrefix "$1"/"$filename".

export picard=/home/tchen/Chentao/picard.jar

java -jar $picard MarkDuplicates \
	I="$1"/"$filename".Aligned.sortedByCoord.out.bam \
	O="$1"/"$filename".md.bam \
	M="$1"/"$filename".md.metrics.txt \
	CREATE_INDEX=true

export gatk=/home/tchen/Chentao/software/gatk-4.0.6.0/gatk
ref=/home/tchen/mouse_project/refer/mouse_GRCm38.fasta

$gatk SplitNCigarReads \
	-R $ref \
	-I "$1"/"$filename".md.bam \
	-O "$1"/"$filename".md.split.bam

$gatk HaplotypeCaller \
	-R $ref \
	--native-pair-hmm-threads 6 \
	-I "$1"/"$filename".md.split.bam \
	--dont-use-soft-clipped-bases true \
	-stand-call-conf 20.0 \
	-O "$1"/"$filename".vcf

$gatk VariantFiltration \
	-R $ref \
	-V "$1"/"$filename".vcf \
	-window 35 \
	-cluster 3 \
	--filter-name FS \
	-filter "FS > 30.0" \
	--filter-name QD \
	-filter "QD < 2.0" \
	-O "$1"/vcfdir/"$filename".filter.vcf

rm "$1"/"$filename".Aligned.sortedByCoord.out.bam
rm "$1"/"$filename"_1.fastq "$1"/"$filename"_2.fastq
