workflow RNACallVariant {
	String inputdir
	String outputdir

	Array[String] inputfirname

	String? starpath
	String star = select_first([starpath, "/home/tchen/Chentao/software/STAR-2.5.4b/bin/Linux_x86_64/STAR"])

	String? picardpath
	String picard = select_first([picardpath, "/home/tchen/Chentao/picard.jar"])


	String? gatkpath
	String gatk = select_first([gatkpath, "/home/tchen/Chentao/software/gatk-4.0.6.0/gatk"])

	String starGenomeRefer
	String gatkGenomeFastaRefer


	Int? thread
	Int Thread = select_first([thread, 8])


	scatter (files in inputfirname) {
		call star_mappint {
			input:
				STARS        =    star,
				threads      =    Thread,
				STARrefer    =    starGenomeRefer,
				Inputdir     =    inputdir,
				Outputdir    =    outputdir,
				Filename     =    files
		}
	}

	Array[File] mapbams = star_mappint.outbam
	scatter (bam in mapbams) {
		call picard_markduplicated {
			input:
				picards      =    picard,
				bamfile      =    bam,
				Filename     =    basename(bam, ".Aligned.sortedByCoord.out.bam"),
				Outputdir    =    outputdir
		}

	}

	Array[File] mardupbams = picard_markduplicated.markdupbam
	scatter (splitbam in mardupbams) {
		call SplitNigarReads {
			input:
				gatks         =     gatk,
				ref           =     gatkGenomeFastaRefer,
				Inputfile     =     splitbam,
				Filename      =     basemane(splitbam, ".md.bam"),
				Outputdir     =     outputdir
		}
	}


	Array[File] splitbams = SplitNigarReads.split_bam
	scatter (vcf_bam in splitbams) {
		call HaplotypeCaller {
			input:
				gatks         =     gatk,
				ref           =     gatkGenomeFastaRefer,
				Inputfile     =     vcf_bam,
				Outputdir     =     outputdir,
				Filename      =     basename(vcf_bam, ".md.split.bam")
		}
	}

	Array[File] vcfs = HaplotypeCaller.vcffile
	scatter (filter_vcf in vcfs) {
		call VariantFilteration {
			input:
				gatks         =       gatk,
				ref           =       gatkGenomeFastaRefer,
				Inputfile     =       filter_vcf,
				Outputdir     =       outputdir,
				Filename      =       basename(filter_vcf, ".vcf")
		}
	}

}


task star_mappint {
	String STARS
	Int threads
	String STARrefer
	String Inputdir
	String Outputdir
	String Filename

	command <<<
		${STARS} \
			--runThreadN ${threads} \
			--twopassMode Basic \
			--genomeDir ${STARrefer} \
			--readFilesIn ${Inputdir}${Filename}_1.clean.fq ${Inputdir}${Filename}_1.clean.fq \
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
			--outSAMattrRGline ID:${Filename} SM:${Filename} PL:illumina LB:${Filename} PU:${Filename} \
			--outSAMtype BAM SortedByCoordinate \
			--outFileNamePrefix ${Outputdir}${Filename}.
	>>>

	output {
		File outbam = "*.out.bam"
	}
}


task picard_markduplicated {
	String picards
	File bamfile
	String Filename
	String Outputdir
	
	command <<<
		java -jar ${picards} MarkDuplicates \
			I=${bamfile} \
			O=${Outputdir}${Filename}.md.bam \
			M=${Outputdir}${Filename}.md.metrics.txt \
			CREATE_INDEX=true
	>>>

	output {
		File  markdupbam = "${Outputdir}${Filename}.md.bam"
	}
	
}

task SplitNigarReads {
	String gatks
	String ref
	File Inputfile
	String Outputdir
	String Filename	
	command <<<
		${gatks} SplitNCigarReads \
			-R ${ref} \
			-I ${Inputfile} \
			-O ${Outputdir}${Filename}.md.split.bam
	>>>

	output {
		File split_bam = "${Outputdir}${Filename}.md.split.bam"
	}
	
}

task HaplotypeCaller {
	String gatks
	String ref
	String Inputfile
	String Outputdir
	String Filename

	command <<<
		${gatks} HaplotypeCaller \
			-R ${ref} \
			--native-pair-hmm-threads 6 \
			-I ${Inputfile} \
			--dont-use-soft-clipped-bases true \
			-stand-call-conf 20.0 \
			-O ${Outputdir}${Filename}.vcf
	>>>

	output {
		File vcffile = "${Outputdir}${Filename}.vcf"
	}
	
}

task VariantFilteration {
	String gatks
	String ref
	File Inputfile
	String Outputdir
	String Filename

	command <<<
		${gatks} VariantFiltration \
			-R ${ref} \
			-V ${Inputfile} \
			-window 35 \
			-cluster 3 \
			--filter-name FS \
			-filter "FS > 30.0" \
			--filter-name QD \
			-filter "QD < 2.0" \
			-O ${Outputdir}${Filename}.filter.vcf
	>>>

	output {
		File filtervcf = "${Outputdir}${Filename}.filter.vcf"
	}
}
