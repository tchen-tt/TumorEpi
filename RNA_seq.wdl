workflow RNAseq {
	Array[String] inputfastq

	String outputdir
	String inputdir

	String refFastaIndex
	String? star_path_override
	String star_path = select_first([star_path_override, "/home/tchen/Chentao/software/STAR-2.5.4b/bin/Linux_x86_64/"])

	String? rsem_calculate_expression
	String rsem_calculate_path = select_first([rsem_calculate_expression, "/home/tchen/Chentao/software/RSEM-master/rsem-calculate-expression"])

	Int? threads
	Int Threads = select_first([threads, 10])

	scatter (fastqfiles in inputfastq) {
		call calculate_expression {
			input:
				fastq_file = fastqfiles,
				output_dir = outputdir,
				input_dir = inputdir,

				refFasta_Index = refFastaIndex,
				star_path_ = star_path,
				rsem_calculate_path_ = rsem_calculate_path,
				Threads_ = Threads
		}
	}
}


task calculate_expression {
	String fastq_file
	String output_dir
	String input_dir

	String refFasta_Index
	String star_path_
	String rsem_calculate_path_
	Int Threads_

	command <<< 
		${rsem_calculate_path_} --star \
			--star-path ${star_path_} \
			--paired-end \
			-p ${Threads_} \
			${input_dir}${fastq_file}_1.clean.fq \
			${input_dir}${fastq_file}_2.clean.fq \
			${refFasta_Index} \
			${output_dir}${fastq_file}
	>>>

	output {
		File out1 = "${output_dir}${fastq_file}.genes.results"
	}

}
