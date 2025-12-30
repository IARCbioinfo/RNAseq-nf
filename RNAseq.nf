#!/usr/bin/env nextflow

// Copyright (C) 2025 IARC/WHO
// License: GNU GPL v3+

nextflow.enable.dsl=2

// ---------------------------
// DEFAULT PARAMETERS
// ---------------------------

params.input_folder      = null
params.input_file        = null
params.ref_folder        = null
params.gtf               = null
params.bed               = null
params.cpu               = 4
params.cpu_gatk          = 1
params.mem               = 50
params.mem_QC            = 2
params.fastq_ext         = "fq.gz"
params.suffix1           = "_1"
params.suffix2           = "_2"
params.output_folder     = "."
params.ref               = "ref.fa"
params.snp_vcf           = "dbsnp.vcf"
params.indel_vcf         = "Mills_1000G_indels.vcf"
params.RG                = "PL:ILLUMINA"
params.STAR_mapqUnique   = 255
params.stranded          = "no"
params.hisat2_idx        = "genome_tran"
params.cpu_trim          = 15
params.htseq_maxreads    = null
params.multiqc_config    = 'NO_FILE'
params.sjtrim            = null
params.recalibration     = null
params.cutadapt          = null
params.hisat2            = null
params.help              = null

log.info ""
log.info "--------------------------------------------------------"
log.info "  RNAseq-nf 2.4.0: alignment, QC, and reads counting workflow for RNA sequencing "
log.info "--------------------------------------------------------"
log.info "Copyright (C) IARC/WHO"
log.info ""

// ---------------------------
// HELP MESSAGE
// ---------------------------

if (params.help) {
    log.info '-------------------------------------------------------------'
    log.info ' USAGE  '
    log.info '-------------------------------------------------------------'
    log.info ''
    log.info 'nextflow run iarcbioinfo/RNAseq.nf [-with-docker] --input_folder input/ --ref_folder ref/ [OPTIONS]'
    log.info ''
    log.info 'Mandatory arguments:'
    log.info '    --input_folder   FOLDER                 Folder containing BAM or fastq files to be aligned.'
    log.info '    --input_file     STRING                 Input file (tab-separated values) with 4 columns: SM (sample name), RG (read group), pair1 (first fastq pair file), and pair2 (second fastq pair file)'
    log.info '    --ref_folder     FOLDER                 Folder with genome reference files (with index).'
    log.info '    --gtf            FILE                   Annotation file.'
    log.info '    --bed            STRING                 Bed file with interval list'
    log.info ""
    log.info 'Optional arguments:'
    log.info '    --ref            FILE                   Reference fasta file (with index) for splice junction trimming and base recalibration.'
    log.info '    --output_folder  STRING                 Output folder (default: .).'
    log.info '    --cpu            INTEGER                Number of cpu used by bwa mem and sambamba (default: 4).'
    log.info '    --cpu_gatk       INTEGER                Number of cpu used by gatk (default: 1).'
    log.info '    --cpu_trim       INTEGER                Number of cpu used by cutadapt (default: 15).'
    log.info '    --mem            INTEGER                Size of memory used for mapping (in GB) (default: 50).'
    log.info '    --mem_QC         INTEGER                Size of memory used for QC and cutadapt (in GB) (default: 2).'
    log.info '    --RG             STRING                 Samtools read group specification (default : PL:ILLUMINA).'
    log.info '    --STAR_mapqUnique INTEGER               STAR default mapping quality for unique mappers (default : 255).'
    log.info '    --fastq_ext      STRING                 Extension of fastq files (default : fq.gz)'
    log.info '    --suffix1        STRING                 Suffix of fastq files 1 (default : _1)'
    log.info '    --suffix2        STRING                 Suffix of fastq files 2 (default : _2)'
    log.info '    --htseq_maxreads INTEGER               Maximum number of reads taken into account by htseq-count (default: htseq-count default=30000000)'
    log.info '    --snp_vcf        STRING                 Path to SNP VCF from GATK bundle (default : dbsnp.vcf)'
    log.info '    --indel_vcf      STRING                 Path to indel VCF from GATK bundle (default : Mills_1000G_indels.vcf)'
    log.info '    --stranded       STRING                 Are reads stranded? (default : no; alternatives : yes, r)'
    log.info '    --hisat2_idx     STRING                 hisat2 index file prefix (default : genome_tran)'
    log.info '    --multiqc_config STRING                 Config yaml file for multiqc (default : none)'
    log.info ''
    log.info 'Flags:'
    log.info '    --sjtrim                                enable splice junction trimming'
    log.info '    --recalibration                         perform base quality score recalibration (GATK)'
    log.info '    --hisat2                                use hisat2 instead of STAR for reads mapping'
    log.info '    --cutadapt                              perform adapter sequence trimming'
    log.info ''
    exit 0
} else {
    log.info "input_folder   = ${params.input_folder}"
    log.info "input_file     = ${params.input_file}"
    log.info "ref            = ${params.ref}"
    log.info "cpu            = ${params.cpu}"
    log.info "cpu_gatk       = ${params.cpu_gatk}"
    log.info "cpu_trim       = ${params.cpu_trim}"
    log.info "mem            = ${params.mem}"
    log.info "mem_QC         = ${params.mem_QC}"
    log.info "fastq_ext      = ${params.fastq_ext}"
    log.info "suffix1        = ${params.suffix1}"
    log.info "suffix2        = ${params.suffix2}"
    log.info "output_folder  = ${params.output_folder}"
    log.info "bed            = ${params.bed}"
    log.info "ref_folder     = ${params.ref_folder}"
    log.info "gtf            = ${params.gtf}"
    log.info "RG             = ${params.RG}"
    log.info "STAR_mapqUnique = ${params.STAR_mapqUnique}"
    log.info "stranded       = ${params.stranded}"
    log.info "hisat2_idx     = ${params.hisat2_idx}"
    log.info "hisat2         = ${params.hisat2}"
    log.info "htseq_maxreads = ${params.htseq_maxreads}"
    log.info "multiqc_config = ${params.multiqc_config}"
    log.info "recalibration  = ${params.recalibration}"
    log.info "sjtrim         = ${params.sjtrim}"
    log.info "cutadapt       = ${params.cutadapt}"
    log.info "snp_vcf        = ${params.snp_vcf}"
    log.info "indel_vcf      = ${params.indel_vcf}"
    log.info "help           = ${params.help}"
}

// ---------------------------
// REF FILES DEFINITION
// ---------------------------

fasta_ref = file(params.ref)
fasta_ref_fai = file("${params.ref}.fai")
if ((params.sjtrim != null) || (params.recalibration != null)) {
    def fasta_ref_dictn = params.ref[0..<params.ref.lastIndexOf('.')]
    fasta_ref_dict = file("${fasta_ref_dictn}.dict")
}

bed = file(params.bed)
gtf = file(params.gtf)

multiqc = params.multiqc_config == 'NO_FILE'
    ? Channel.empty()
    : file(params.multiqc_config)

def aligner_ref
if (params.hisat2) {
    def pfx = "${params.ref_folder}/${params.hisat2_idx}"
    aligner_ref = Channel.fromPath("${pfx}.1.ht2")
            .concat(Channel.fromPath("${pfx}.2.ht2"),
                    Channel.fromPath("${pfx}.3.ht2"),
                    Channel.fromPath("${pfx}.4.ht2"),
                    Channel.fromPath("${pfx}.5.ht2"),
                    Channel.fromPath("${pfx}.6.ht2"),
                    Channel.fromPath("${pfx}.7.ht2"),
                    Channel.fromPath("${pfx}.8.ht2"))
} else { // STAR
    aligner_ref = Channel.fromPath("${params.ref_folder}/chrStart.txt")
            .concat(Channel.fromPath("${params.ref_folder}/chrNameLength.txt"),
                    Channel.fromPath("${params.ref_folder}/chrName.txt"),
                    Channel.fromPath("${params.ref_folder}/chrLength.txt"),
                    Channel.fromPath("${params.ref_folder}/exonGeTrInfo.tab"),
                    Channel.fromPath("${params.ref_folder}/exonInfo.tab"),
                    Channel.fromPath("${params.ref_folder}/geneInfo.tab"),
                    Channel.fromPath("${params.ref_folder}/Genome"),
                    Channel.fromPath("${params.ref_folder}/genomeParameters.txt"),
                    Channel.fromPath("${params.ref_folder}/SA"),
                    Channel.fromPath("${params.ref_folder}/SAindex"),
                    Channel.fromPath("${params.ref_folder}/sjdbInfo.txt"),
                    Channel.fromPath("${params.ref_folder}/transcriptInfo.tab"),
                    Channel.fromPath("${params.ref_folder}/sjdbList.fromGTF.out.tab"),
                    Channel.fromPath("${params.ref_folder}/sjdbList.out.tab"))
}

known_snps = file( params.snp_vcf )
known_snps_index = file( params.snp_vcf+'.tbi' )
known_indels = file( params.indel_vcf )
known_indels_index = file( params.indel_vcf+'.tbi' )

// ---------------------------
// INPUT CHECKS - either infile mode (tab) or fastq/bam scan
// ---------------------------

def mode = null
if (params.input_file) {
    mode = 'infile'
} else {
    if (file(params.input_folder).listFiles().findAll { it.name ==~ /.*${params.fastq_ext}/ }.size() > 0) {
		mode = 'fastq'
        println "fastq files found, proceed with alignment"
		
    } else {
        if (file(params.input_folder).listFiles().findAll { it.name ==~ /.*bam/ }.size() > 0) {
			mode = 'bam'
            println "BAM files found, proceed with realignment"
        } else {
            println "ERROR: input folder contains no fastq nor BAM files"
            System.exit(0)
        }
    }
}


// ---------------------------
// PROCESSES
// ---------------------------

    process BAM2FASTQ {
        tag { file_tag }
        cpus 1
        memory "${params.mem_QC}G"

        input:
        tuple val(file_tag), val(rg), path(infile) 

        output:
        tuple val(file_tag), val(rg), path("${file_tag}${params.suffix1}.${params.fastq_ext}"), path("${file_tag}${params.suffix2}.${params.fastq_ext}"), emit: readPairs0

        script:
        """
        set -euo pipefail
        samtools collate -u -O -n 128 ${infile} tmp_${file_tag} |
		samtools fastq -1 ${file_tag}${params.suffix1}.${params.fastq_ext} -2 ${file_tag}${params.suffix2}.${params.fastq_ext} -0 /dev/null -s /dev/null -n -
        """
    }

	process FASTQC_PRETRIM {
		tag { file_tag }
		cpus params.cpu
		memory "${params.mem_QC}GB"

		input:
		tuple val(file_tag), val(rg), path(pair1), path(pair2) 

		output:
		path "*_pretrim_fastqc.zip", emit: fastqc_pairs

		publishDir "${params.output_folder}/QC/fastq", mode: 'copy', pattern: '*fastqc.zip'

		script:
		"""
    	set -euo pipefail

    	ext='${params.fastq_ext}'
    	suffix1='${params.suffix1}'
    	suffix2='${params.suffix2}'

    	pair1_file='${pair1}'
    	pair2_file='${pair2}'
    	file_tag='${file_tag}'
    	rg='${rg}'

    	basename1=\$(basename "\$pair1_file" ".\$ext")

    	if [ "\$(basename "\$pair2_file")" != "NO_fastq2" ]; then
        	basename2=\$(basename "\$pair2_file" ".\$ext")
        	fastqc -t ${task.cpus} "\$pair1_file" "\$pair2_file"
       	 	mv "\${basename1}_fastqc.zip" "\${file_tag}\${suffix1}\${rg}_pretrim_fastqc.zip"
        	mv "\${basename2}_fastqc.zip" "\${file_tag}\${suffix2}\${rg}_pretrim_fastqc.zip"
    	else
        	fastqc -t ${task.cpus} "\$pair1_file"
        	mv "\${basename1}_fastqc.zip" "\${file_tag}\${suffix1}\${rg}_pretrim_fastqc.zip"
    	fi
    	"""
}

	process MULTIQC_PRETRIM {
		tag { "all" }
		cpus 1
		memory "${params.mem_QC}GB"

		input:
		path fastqc1 
		path multiqc_config 

		output:
		path "multiqc_pretrim_report.html" , emit: multiqc_pre
		path "multiqc_pretrim_report_data" , emit: multiqc_pre_data

		publishDir "${params.output_folder}/QC", mode: 'copy'

		script:
 		"""
    	set -euo pipefail

    	config_file='${multiqc_config}'

    	if [ "\$(basename "\$config_file")" = "NO_FILE" ]; then
        	opt=""
    	else
        	opt="--config \$config_file"
    	fi

 		while IFS= read -r -d '' f; do
        	cp --remove-destination "\$(readlink "\$f")" "\$f" || true
    		done < <(find . -name "*_pretrim_fastqc.zip" -type l -print0)
  		multiqc . -n multiqc_pretrim_report.html -m fastqc \$opt --comment "RNA-seq Pre-trimming QC report"
    		"""
		}

    process ADAPTER_TRIMMING {
        tag { file_tag + rg }
        cpus params.cpu_trim
        memory "${params.mem_QC}GB"

        input:
        tuple val(file_tag), val(rg), path(pair1), path(pair2) 
		// from readPairs

        output:
//        tuple val(file_tag), val(rg), path("${file_tag}${rg}*val_1.fq.gz"), path("${file_tag}${rg}*val_2.fq.gz") , emit: readPairs2
		tuple val(file_tag), val(rg), path("${file_tag}${rg}_val_1.fq.gz"), path("${file_tag}${rg}_val_2.fq.gz"), emit: readPairs2
        path "*_fastqc.zip" , emit: fastqc_postpairs
        path "*trimming_report.txt" , emit: trimming_reports

        publishDir "${params.output_folder}/QC/adapter_trimming", mode: 'copy', pattern: '*report.txt,*fastqc.zip'

        script:
        """
		set -euo pipefail

        cpu_tg=\$(( ${task.cpus} - 1 ))
        cpu_tg3=\$(python - <<PY
		import math
		print(max(1, int(math.ceil(${task.cpus} / 3.5))))
		PY
		)
		if [ "\$(basename ${pair2})" != "NO_fastq2" ]; then
        	opts="--paired"
        	trim_galore \$opts --fastqc --gzip \
            	--basename ${file_tag}${rg} \
            	-j \$cpu_tg3 \
            	${pair1} ${pair2}
    	else
        	trim_galore --fastqc --gzip \
            	--basename ${file_tag}${rg} \
            	-j \$cpu_tg3 \
            	${pair1}

        # Normalise to paired-like outputs
        	mv ${file_tag}${rg}_trimmed.fq.gz ${file_tag}${rg}_val_1.fq.gz
        	touch ${file_tag}${rg}_val_2.fq.gz
   		fi
        """
    }

	process ALIGNMENT {
		tag { file_tag }
		cpus params.cpu
		memory "${params.mem}G"

		input:
		tuple val(file_tag), val(rg), path(pair1), path(pair2) 
		path star_index
		// path ref 
		path gtf

		output:
		tuple val(file_tag), val(rg), path("${file_tag}.bam"), path("${file_tag}.bam.bai") , emit: bam_files
		path "*Log*" , emit: align_out
		tuple val(file_tag), path("*SJ.out.junction") , emit: SJ_out
		path "*SJ.out.tab" , emit: SJ_out_others

		script:
    	"""
    	set -euo pipefail

		align_threads=$(( ${params.cpu} / 2 ))
    	(( align_threads < 1 )) && align_threads=1

    	sort_threads=$(( ${params.cpu} / 2 - 1 ))
    	(( sort_threads < 1 )) && sort_threads=1

    	sort_mem=$(( ${params.mem} / 4 ))

		input_f1="${pair1}"
		rgline="ID:${file_tag} SM:${file_tag} ${params.RG}"
		if [ -n "${pair2}" ] && [ "$(basename ${pair2})" != "NO_fastq2" ]; then
    		pairs="${pair1} ${pair2}"
		else
    		pairs="${pair1}"
		fi

		STAR \
        --genomeDir ${star_index} \
        --sjdbGTFfile ${gtf} \
        --runThreadN ${align_threads} \
        --readFilesCommand zcat \
        --readFilesIn \$pairs \
        --outStd SAM \
        --outSAMattrRGline "${rgline}" \
        --outSAMmapqUnique ${params.STAR_mapqUnique} \
        --chimSegmentMin 12 \
        --chimJunctionOverhangMin 12 \
        --chimSegmentReadGapMax 3 \
        --alignSJDBoverhangMin 10 \
        --alignMatesGapMax 100000 \
        --alignIntronMax 100000 \
        --alignSJstitchMismatchNmax 5 -1 5 5 \
        --outSAMstrandField intronMotif \
        --chimMultimapScoreRange 10 \
        --chimMultimapNmax 10 \
        --chimNonchimScoreDropMin 10 \
        --peOverlapNbasesMin 12 \
        --peOverlapMMp 0.1 \
        --chimOutJunctionFormat 1 \
        --twopassMode Basic \
        --outReadsUnmapped None \
    	| samblaster --addMateTags \
    	| sambamba view -S -f bam -l 0 /dev/stdin \
    	| sambamba sort -t ${sort_threads} -m ${sort_mem}G --tmpdir=${file_tag}_tmp -o ${file_tag}.bam /dev/stdin
	  	
		sambamba index -t \$sort_threads ${file_tag}.bam

    	mv Chimeric.out.junction STAR.${file_tag}.Chimeric.SJ.out.junction || true
    	mv SJ.out.tab STAR.${file_tag}.SJ.out.tab || true
    	mv Log.final.out STAR.${file_tag}.Log.final.out || true
    	mv Log.out STAR.${file_tag}.Log.out || true
    	mv Log.progress.out STAR.${file_tag}.Log.progress.out || true
    	mv Log.std.out STAR.${file_tag}.Log.std.out || true
    	"""
}

    process SPLICE_JUNCT_TRIM {
        tag { file_tag }
        cpus params.cpu_gatk
        memory "${params.mem}G"

        input:
        tuple val(file_tag), val(rg), path(bam), path(bai)
        path fasta_ref
        path fasta_ref_fai
        path fasta_ref_dict

        output:
        tuple val("${file_tag}_split"), val(rg), path("${file_tag}_split.bam"), path("${file_tag}_split.bam.bai"), emit: bam_files2

        script:
        """
  		set -euo pipefail
        gatk SplitNCigarReads --java-options "-Xmx${task.memory.toGiga()}G" -R ${fasta_ref} -I ${bam} -O ${file_tag}_split.bam
        """
    }

    process BASE_QUALITY_SCORE_RECALIBRATION {
        tag { file_tag }
        cpus params.cpu_gatk
        memory "${params.mem}G"

        publishDir "${params.output_folder}/BAM/", mode: 'copy', pattern: "*bam*"
        publishDir "${params.output_folder}/QC/BAM/BQSR/", mode: 'copy',
                   saveAs: { filename ->
                       if (filename.indexOf("table") > 0) "${filename}"
                       else if (filename.indexOf("plots") > 0) "${filename}"
                       else null
                   }

        input:
        tuple val(file_tag), val(rg), path(bam), path(bai)
    	path known_snps
    	path known_snps_index
    	path known_indels
    	path known_indels_index
    	path fasta_ref
    	path fasta_ref_fai
    	path fasta_ref_dict

        output:
    	path("${file_tag}_recal.table"), emit: recal_table_files
    	path("${file_tag}_recalibration_plots.pdf"), emit: recal_plots_files
    tuple val("${file_tag}_BQSRecalibrated"), val(rg), path("${file_tag}_BQSRecalibrated.bam"), path("${file_tag}_BQSRecalibrated.bam.bai"), emit: bam_files3

        script:
		"""
    	set -euo pipefail
    	file_tag_new=${file_tag}_BQSRecalibrated
		gatk BaseRecalibrator --java-options "-Xmx${task.memory.toGiga()}G" -R ${fasta_ref} -I ${bam} --known-sites ${known_snps} --known-sites ${known_indels} -O ${file_tag}_recal.table
   		gatk ApplyBQSR --java-options "-Xmx${task.memory.toGiga()}G" -R ${fasta_ref} -I ${bam} --bqsr-recal-file ${file_tag}_recal.table -O ${file_tag_new}.bam
    	gatk BaseRecalibrator --java-options "-Xmx${task.memory.toGiga()}G" -R ${fasta_ref} -I ${file_tag_new}.bam --known-sites ${known_snps} --known-sites ${known_indels} -O ${file_tag_new}_recal.table
   		gatk AnalyzeCovariates --java-options "-Xmx${task.memory.toGiga()}G" -before ${file_tag}_recal.table -after ${file_tag_new}_recal.table -plots ${file_tag_new}_recalibration_plots.pdf
    	"""
    }

	process RSEQC { //(read distribution, clipping, junction saturation)
		tag { file_tag }
		cpus 1
		memory "${params.mem_QC}GB"

		input:
			tuple val(file_tag), val(rg), path(bam), path(bai)
			// from bam_files_for_quantif
			file bed
			// from bed

		output:
			path "${file_tag}_readdist.txt" , emit: rseqc_files
			path "*clipping*" , emit: rseqc_clip_files
			path "*jun_saturation*" , emit: rseqc_jsat_files

		publishDir "${params.output_folder}/QC/bam", mode: 'copy'

		script:
		'''
		read_distribution.py -i ${bam} -r ${bed} > ${file_tag}_readdist.txt
		clipping_profile.py  -i ${bam} -s "PE" -o ${file_tag}_clipping
		junction_saturation.py -i ${bam} -r ${bed} -o ${file_tag}_jun_saturation
		'''
}

	process RSEQCSPLIT {
		tag { file_tag }
		cpus 1
		memory "${params.mem_QC}GB"
	
		input:
		tuple val(file_tag), val(rg), path(bam), path(bai)
		// from bam_files_for_quantif
		file bed
		// from bed

		output:
		path "*readdist.txt" , emit: rseqc_files_split

		publishDir "${params.output_folder}/QC/bam", mode: 'copy'

		script:
		'''
		basename=$(basename ${bam})
		samtools split !{bam} -f "%*_%!.%."
		for f in ${basename}_*.bam; do
			read_distribution.py -i $f -r ${bed} > ${f%.bam}_readdist.txt
		done
		'''
}

	process QUANTIFICATION {
		tag { file_tag }
		cpus params.cpu
		memory { (params.sjtrim || params.recalibration) ? "${params.mem}G" : "${params.mem_QC}G" }

		input:
		tuple val(file_tag), val(rg), path(bam), path(bai)
		// from bam_files_for_quantif
		file gtf
		// from gtf

		output:
		path "${file_tag}_count.txt" , emit: htseq_files

		publishDir "${params.output_folder}/counts", mode: 'copy'

		script:
		'''
		buffer=""
		if [ -n "${params.htseq_maxreads}" ]; then
			buffer="--max-reads-in-buffer ${params.htseq_maxreads}"
		fi
		if [ -n "${params.sjtrim}" ] || [ -n "${params.recalibration}" ]; then
			mv ${file_tag}.bam ${file_tag}_coordinate_sorted.bam
			sambamba sort -n -t ${task.cpus} -m ${params.mem}G --tmpdir=${file_tag}_tmp -o ${file_tag}.bam ${file_tag}_coordinate_sorted.bam
			htseq-count -n ${params.cpu} -r name -s ${params.stranded} -f bam ${file_tag}.bam ${gtf} ${buffer} --additional-attr=gene_name > ${file_tag}_count.txt
		else
			htseq-count -n ${params.cpu} -r pos -s ${params.stranded} -f bam ${file_tag}.bam ${gtf} ${buffer} --additional-attr=gene_name > ${file_tag}_count.txt
		fi
		'''
}

	process MULTIQC_POSTTRIM {
		tag { "all" }
		cpus 1
		memory "${params.mem_QC}GB"

		input:
		file STAR
		// from align_out
		file htseq
		// from htseq_files
		file rseqc_clip
		// from rseqc_clip_files
		file rseqc
		// from rseqc_files
		file rseqc_jsat
		// from rseqc_jsat_files
		file trim
		// from trimming_reports.ifEmpty([])
		file fastqcpost
		// from fastqc_postpairs.ifEmpty([])
		file rseqc_split
		// from rseqc_files_split.ifEmpty([])
		file multiqc_config
		// from multiqc

		output:
		path "multiqc_posttrim_report.html" , emit: multiqc_post
		path "multiqc_posttrim_report_data" , emit: multiqc_post_data

		publishDir "${params.output_folder}/QC", mode: 'copy'

		script:
		'''
		if [ "$(basename ${multiqc_config})" == "NO_FILE" ]; then
			opt=""
		else
			opt="--config ${multiqc_config}"
		fi
		if compgen -G "*fastq.zip" > /dev/null; then
			for f in $(find . -name "*_fastqc.zip" -type l); do cp --remove-destination $(readlink $f) $f || true; done
		fi
		multiqc . -n multiqc_posttrim_report.html -m fastqc -m cutadapt -m star -m rseqc -m htseq ${opt} --comment "RNA-seq Post-trimming QC report"
		'''
}

// ========================================================================================================================================================
// WORKFLOW
// ========================================================================================================================================================

workflow {

     // ----------------------------------------
     // 0. INPUT NORMALISATION
     // ----------------------------------------

    def readPairs
    def readPairs2
    def files

	// If file as input
     if (mode == 'infile') {
   		readPairs = Channel.fromPath(params.input_file)
        .splitCsv(header: true, sep: '\t', strip: true)
        .map { row -> tuple(row.SM, row.RG, file(row.pair1), file(row.pair2)) }
		readPairs2 = readPairs
		}
	 
	 // If BAM as input : process bam -> fastq //
		else if (mode == 'bam') {
		files = Channel.fromPath("${params.input_folder}/*.bam")
                          .map { path -> tuple(path.baseName, '', path) }
		
		def bam2fq_out = BAM2FASTQ(files)
        readPairs = bam2fq_out.readPairs0
        readPairs2 = bam2fq_out.readPairs0
		} 
	
	// IF FASTQ as input: build readPairs/readPairs2 channels if not already filled /////
		else if (mode == 'fastq') {
    		if (params.suffix2) {
        		readPairs = Channel.fromFilePairs("${params.input_folder}/*{${params.suffix1},${params.suffix2}}.${params.fastq_ext}")
               .map { row -> tuple(row[0], '', row[1][0], row[1][1]) }
    			} else {
       	 				readPairs = Channel.fromPath("${params.input_folder}/*${params.suffix1}.${params.fastq_ext}")
               			.map { row -> tuple(row.name.replace("${params.suffix1}.${params.fastq_ext}", ''),'', row, file('NO_fastq2') ) }
   				 		}
						readPairs2 = readPairs
			}

     // ----------------------------------------
     // 1. FASTQC PRETRIM
     // ----------------------------------------

	def fastqc1 = FASTQC_PRETRIM(readPairs)

    // --------------------------------------------------------------
    // 2. MULTIQC PRETRIM
    // --------------------------------------------------------------

	fastqc_pretrim_all = fastqc1.collect()
	MULTIQC_PRETRIM(fastqc_pretrim_all,multiqc)

    // --------------------------------------------------------------
    // 3. OPTIONAL ADAPTER TRIMMING
    // --------------------------------------------------------------

	def readPairs_for_align
	def trim_reports_ch = Channel.empty()
	def fastqc_postpairs_ch = Channel.empty()
	if (params.cutadapt) {
 	def trim = ADAPTER_TRIMMING(readPairs)
	readPairs_for_align = readPairs2
	trim_reports_ch = trim.trimming_reports
	fastqc_postpairs_ch = trim.fastqc_postpairs
	} else {
			readPairs_for_align = readPairs
    		}

    // --------------------------------------------------------------
    // 4. ALIGNMENT
    // --------------------------------------------------------------
     
	def align = ALIGNMENT(readPairs_for_align,aligner_ref,gtf)

    // --------------------------------------------------------------
    // 5. OPTIONAL SPLICE JUNCTION TRIM
    // --------------------------------------------------------------


	def bam_files_for_bqsr
	if (params.sjtrim) {
        def sjt = SPLICE_JUNCT_TRIM(align.bam_files,fasta_ref,fasta_ref_fai,fasta_ref_dict)
		bam_files_for_bqsr = sjt.bam_files2
		} else {
				 bam_files_for_bqsr = align.bam_files
				}

    // --------------------------------------------------------------
    // 6. OPTIONAL BQSR
    // --------------------------------------------------------------

	def bam_files_for_quantif
	if (params.recalibration) {
        def bq = BASE_QUALITY_SCORE_RECALIBRATION(bam_files_for_bqsr,known_snps,known_snps_index,known_indels,known_indels_index,fasta_ref,fasta_ref_fai,fasta_ref_dict)
		bam_files_for_quantif = bq.bam_files3
		} else {
        		bam_files_for_quantif = bam_files_for_bqsr
    			} 
/*
    // --------------------------------------------------------------
    // 7. RSEQC
    // --------------------------------------------------------------

    def rs = RSEQC(bam_files_for_quantif,bed)

    // --------------------------------------------------------------
    // 8. RSEQCSPLIT
    // --------------------------------------------------------------

    def rss = RSEQCSPLIT(bam_files_for_quantif,bed)

    // --------------------------------------------------------------
    // 9. QUANTIFICATION (RNA only)
    // --------------------------------------------------------------
   
	def quant = QUANTIFICATION(bam_files_for_quantif,gtf)

    // --------------------------------------------------------------
    // 10. MULTIQC POSTRIM
    // --------------------------------------------------------------

	align_out_all        = align.align_out.collect()
	htseq_all            = quant.htseq_files.collect()
	rseqc_clip_all       = rs.rseqc_clip_files.collect()
	rseqc_all            = rs.rseqc_files.collect()
	rseqc_jsat_all       = rs.rseqc_jsat_files.collect()
	trim_reports_all     = trim_reports_ch.collect()
	fastqc_post_all      = fastqc_postpairs_ch.collect()
	rseqc_split_all      = rss.rseqc_files_split.collect()

	MULTIQC_POSTTRIM(align_out_all,htseq_all,rseqc_clip_all,rseqc_all,rseqc_jsat_all,trim_reports_all,fastqc_post_all,rseqc_split_all,multiqc)
*/
}
