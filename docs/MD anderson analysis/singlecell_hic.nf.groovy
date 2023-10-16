// Declare syntax version
//https://github.com/danrlu/Nextflow_cheatsheet/blob/main/nextflow_cheatsheet.pdf
nextflow.enable.dsl=2
// Script parameters
params.flowcellDir = "/volumes/seq/flowcells/MDA/nextseq2000/2023/230306_VH00219_371_AACJJFWM5"
params.outname = "test"
params.outdir = "/volumes/seq/projects/gccACT/nextflow_test"
params.ref = "/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa"
params.chrsizes = "/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/hg38.chrom.sizes"
params.bins = "/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/1mb.bins"
params.bins_bed="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/1mb.bins.bed"
params.gc_bins = "/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/1mb.gc.bins"
params.gc_bins_bed = "/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/1mb.gc.bed"

log.info """\

		================================================
		             scWGS/HiC NF PIPELINE
		================================================
		Flowcell Dir : ${params.flowcellDir}
		Output Prefix : ${params.outname}
		NF Working Dir : ${workflow.launchDir}
		Reference Fasta: ${params.ref}
		Bins : ${params.bins}
		Output Directory : ${params.outdir}

""".stripIndent()

/* BCL TO FASTQ PIPELINE */
process BCL_TO_FASTQ { 
	//Generate Undetermined Fastq Files from BCL Files.
	conda "/volumes/USR2/Ryan/miniconda3/envs/r4.2"
	cpus 100

	input:
		path flowcellDir
	output:
		tuple path("Undetermined_S0_L001_R1_001.fastq.gz"),
		path("Undetermined_S0_L001_R2_001.fastq.gz"),
		path("Undetermined_S0_L001_I1_001.fastq.gz"),
		path("Undetermined_S0_L001_I2_001.fastq.gz")
	script:
		"""
		bcl2fastq -R $flowcellDir \\
		-o . \\
		-r 10 \\
		-p ${task.cpus} \\
		-w 10 \\
		--ignore-missing-bcls \\
		--ignore-missing-filter \\
		--ignore-missing-positions \\
		--ignore-missing-controls \\
		--create-fastq-for-index-reads

		"""
}

process CHUNK_FASTQ {
	//Chunk fastq files to 50 chunks reads to parallelize demultiplexing
	conda '/volumes/USR2/Ryan/miniconda3/envs/r4.2'
	cpus 150

	input:
		tuple path(read1), path(read2), path(ind1), path(ind2)
	output:
		tuple path("${read1.simpleName}.chunk*gz"), path("${read2.simpleName}.chunk*gz"),
			path("${ind1.simpleName}.chunk*gz"), path("${ind2.simpleName}.chunk*gz")
	script:
	"""
	seqkit split2 ${read1} -p 50 -j ${task.cpus} -O . --by-part-prefix ${read1.simpleName}.chunk -e .gz &
	seqkit split2 ${read2} -p 50 -j ${task.cpus} -O . --by-part-prefix ${read2.simpleName}.chunk -e .gz &
	seqkit split2 ${ind1} -p 50 -j ${task.cpus} -O . --by-part-prefix ${ind1.simpleName}.chunk -e .gz &
	seqkit split2 ${ind2} -p 50 -j ${task.cpus} -O . --by-part-prefix ${ind2.simpleName}.chunk -e .gz &
	"""
}

process DEMUX_FASTQ { 
	//Assign fastq files to determined index pairs.
	conda '/volumes/USR2/Ryan/miniconda3/envs/r4.2'
	cpus 50

	input:
		tuple path(read1), path(read2), path(idx1), path(idx2)
	output:
		path("*R1*.barc.fastq")
		path("*R2*.barc.fastq")
	script:
		"""
		python /volumes/USR2/Ryan/src/plate2_fastqsplitter.nf.py \\
		$read1 \\
		$read2 \\
		$idx1 \\
		$idx2

		"""
}

process UNCHUNK_FASTQ {
	//Chunk fastq files to 5M reads to parallelize demultiplexing
	conda '/volumes/USR2/Ryan/miniconda3/envs/r4.2'
	
	input:
		path(read1)
		path(read2)
	output:
		tuple path("test.R1.fq.gz"),
			path("test.R2.fq.gz")
	script:
	"""
	cat ${read1} | gzip > test.R1.fq.gz &
	cat ${read2} | gzip > test.R2.fq.gz &
	"""
}


/* ALIGNMENT AND SPLITTING CELLS */
process BWA_ALIGN {
	//Map reads with BWA mem
	conda '/volumes/USR2/Ryan/miniconda3/envs/r4.2'
	cpus 100

	input:
		tuple path(read1), path(read2)
		path bwa_index
		val outname
	output:
		tuple val("${outname}"), path("${outname}.bam")
	script:
		def idxbase = bwa_index[0].baseName
		"""
		bwa mem \\
		-t ${task.cpus} \\
		${idxbase} \\
		${read1} \\
		${read2} \\
		| samtools view -@ ${task.cpus} -b - > ${outname}.bam
		"""
}

process SPLIT_BAM_BY_READNAME {
	//Split bam file by read names
	conda '/volumes/USR2/Ryan/miniconda3/envs/r4.2'

	input:
		tuple val(outname), path(bam)
	output:
		path("*sam")

	script:
	"""
	samtools view ${bam} \\
	| awk \\
	\'OFS=\"\\t\" {split(\$1,a,\":\"); print \$0,\"XM:Z:\"a[1] \\
	> a[1]\".${bam.simpleName}.sam\"}\'
	"""
}


/* CONVERT TO BAMS AND RUN QC */
process BAM_CONVERT {
	//Convert sam to bam and headers on single cell bams
	conda '/volumes/USR2/Ryan/miniconda3/envs/r4.2'
	cpus 100

	input:
		path sc_sams
		path fasta_ref
	output:
		path("*.sc.bam")

	script:
	"""
	samtools view -bT ${fasta_ref} ${sc_sams} \\
	| samtools sort -o ${sc_sams.baseName}.sc.bam -
	"""
}

process BAM_MARKDUP {
	//Fix mates, sort and mark duplicates in bam
	conda '/volumes/USR2/Ryan/miniconda3/envs/r4.2'
	cpus 50
	publishDir "${params.outdir}/data/cells", mode: 'copy'

	input:
		path bam
	output:
		path("*.sc.bbrd.bam")
	script:
	"""
	samtools sort -n -o - ${bam} \\
	| samtools fixmate -m - - \\
	| samtools sort -T . -o - - \\
	| samtools markdup -s - ${bam.baseName}.bbrd.bam 2> ${bam.baseName}.rmdup.stats.txt
	"""
}

process FASTQC {
	//Generate FastQC per bam
	//Based on https://github.com/nf-core/rnaseq/blob/37f260d360e59df7166cfd60e2b3c9a3999adf75/main.nf#L473
	conda '/volumes/USR2/Ryan/miniconda3/envs/r4.2'
	cpus 50
	publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

	input:
		path bam
	output:
    	file "*_fastqc.{zip,html}" into fastqc_results
	script:
	"""
	fastqc ${bam}
	"""
}

process MULTIQC {
	//Run MultiQC
	conda '/volumes/USR2/Ryan/miniconda3/envs/r4.2'
	cpus 50	
	publishDir '/data/qc', mode: 'copy'

	input:
		path bam
	output:
		path("*.sc.bbrd.bam")	
	script:
	"""
	multiqc --outdir /data/qc/ -f .
	"""
}

/* CNV PROFILING */
process CNV_CLONES {
	conda '/volumes/USR2/Ryan/miniconda3/envs/r4.2'
	cpus 50
	publishDir "${params.outdir}/data", mode: 'copy'

	input:
		path sc_dedup_bams
	output:
		path("*.bam_list.txt")
	script:
		"""
		Rscript /volumes/USR2/Ryan/src/copykit_cnv_clones.nf.R \\
		${sc_dedup_bams[0]}.Parent \\
		${task.cpus}
		
		"""
}

/* MAKE CONTACT MATRICES */
process CLONE_CONTACT_MATRICES {
	conda '/volumes/USR2/Ryan/miniconda3/envs/EagleC'
	cpus 100
	publishDir '/data/clone_contacts', mode: 'copy'


	input:
		path bam_list
		path ref
	output:
		path("${bam_list.simpleName}.contacts.bam")
	script:
	"""	
	samtools merge -b ${bam_list} -O SAM -@ ${task.cpus} - \\
	| awk \'{if (sqrt((\$9^2))>=1000 || \$7 != \"=\") print \$0}\' \\
	| samtools view -bT ${ref} - > \\
	${bam_list.simpleName}.contacts.bam

	bam2pairs \\
	${bam_list.simpleName}.contacts.bam \\
	${bam_list.simpleName}
	"""
}

process CELL_CONTACT_MATRICES {
	conda '/volumes/USR2/Ryan/miniconda3/envs/EagleC'
	cpus 100
	publishDir '/data/sc_contacts', mode: 'copy'


	input: 
		path sc_dedup_bams
		path ref
	output:
		path("${sc_dedup_bams.simpleName}.contacts.bam")
	script:
	"""
	samtools view -F 1024 ${sc_dedup_bams} \\
	| awk \'{if (sqrt((\$9^2))>=1000 || \$7 != \"=\") print \$0}\' \\
	| samtools view -bT ${ref} - > \\
	${sc_dedup_bams.simpleName}.contacts.bam
	
	bam2pairs \\
	${sc_dedup_bams.simpleName}.contacts.bam \\
	${sc_dedup_bams.simpleName}
	"""
}

process COOLER_MATRICES {
	conda '/volumes/USR2/Ryan/miniconda3/envs/EagleC'
	cpus 100
	input:
		path pairix_in
		path bins_bed
	output:
		path("${pairix_in.simpleName}.${bins_bed.simpleName}.cool")
	script:
	"""
	cooler cload pairix \\
	-0 \\
	-p ${task.cpus} \\
	--assembly hg38 \\
	${bins_bed} \\
	${pairix_in} \\
	${pairix_in.simpleName}.${bins_bed.simpleName}.cool
	"""
}

process PLOTTING_TRACKS {
	conda '/volumes/USR2/Ryan/miniconda3/envs/r4.2'
	publishDir '/data', mode: 'copy'
	cpus 50

	input:
		path sc_dedup_bams
		path bins_bed
	output:
		path("*.segmented.{bigWig,bedgraph}")
	script:
		"""
		Rscript /volumes/USR2/Ryan/src/copykit_plotting_tracks.nf.R \\
		${sc_dedup_bams[0]}.Parent \\
		${bins_bed} \\
		${task.cpus}
		
		"""
}


workflow {

	/* SETTING UP VARIABLES */
	def fasta_ref = Channel.value(params.ref)
	def outname = Channel.value(params.outname)
	def chrsizes = Channel.value(params.chrsizes)
	def bins = Channel.value(params.bins)
	def bins_bed = Channel.value(params.bins_bed)
	def gc_bins = Channel.value(params.gc_bins)
	def gc_bins_bed = Channel.value(params.gc_bins_bed)
	bwa_index = file('/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa{,.amb,.ann,.bwt,.pac,.sa}' )

	/* BCL TO FASTQ PIPELINE */
	(fq1, fq2) =Channel.fromPath(params.flowcellDir) \
	| BCL_TO_FASTQ \
	| CHUNK_FASTQ \
	| transpose \
	| DEMUX_FASTQ 
	fq1 = fq1 | collect
	fq2 = fq2 | collect
	unchunked_fqs=UNCHUNK_FASTQ(fq1, fq2)
	
	/* ALIGNMENT AND SPLITTING CELLS */
	sc_sams = \
	BWA_ALIGN(unchunked_fqs,bwa_index,outname) \
	| SPLIT_BAM_BY_READNAME \
	| flatten 

	sc_dedup_bams = \
	BAM_CONVERT(sc_sams,fasta_ref) \
	| flatten \
	| BAM_MARKDUP

	/* QC ON CELLS
	| FASTQC(sc_dedup_bams.flatten()) \
	| MULTIQC */
	
	/* CALL CNVs AND GENERATE CLONE LISTS 
	clone_lists= \
	CNV_CLONES(sc_dedup_bams) */

	/* GENERATE COOLER FILES 
	CLONE_CONTACT_MATRICES(clone_lists,fasta_ref)
	CELL_CONTACT_MATRICES(sc_dedup_bams,fasta_ref) */

	/* GENERATE PLOTTING TRACKS 
	PLOTTING_TRACKS */


}
