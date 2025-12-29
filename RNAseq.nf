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

bed = file(params.bed) // OR    bed_ch = params.bed ? Channel.value(file(params.bed)) : Channel.empty()
gtf = file(params.gtf) // OR	gtf_ch = params.gtf ? Channel.value(file(params.gtf)) : Channel.empty()

    ref_ch = Channel.value(file(params.ref))



ch_config_for_multiqc = file(params.multiqc_config)

Def aligner_ref // STAR or HISAT2
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

// ---------------------------
// INPUT CHECKS - either infile mode (tab) or fastq/bam scan
// ---------------------------


Def readPairs = Channel.create()
Def readPairs2 = Channel.create()

def mode = null
if (params.input_file) {
    mode = 'infile'
    Channel.fromPath("${params.input_file}")
        .splitCsv(header: true, sep: '\t', strip: true)
        .map { row -> [ row.SM, row.RG, file(row.pair1), file(row.pair2) ] }
        .into(readPairs, readPairs2)
} else {
    if (file(params.input_folder).listFiles().findAll { it.name ==~ /.*${params.fastq_ext}/ }.size() > 0) {
        println "fastq files found, proceed with alignment"
		mode = 'fastq'
    } else {
        if (file(params.input_folder).listFiles().findAll { it.name ==~ /.*bam/ }.size() > 0) {
            println "BAM files found, proceed with realignment"
            mode = 'bam'
            Channel.fromPath("${params.input_folder}/*.bam")
                   .map { path -> [ path.name.replace(".bam", ""), "", path ] }
                   .into(readPairs)   // into single channel for downstream processing
            // we will set readPairs2 later in workflow when necessary
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
        set val(file_tag), val(rg), path(infile) from files

        output:
        set val(file_tag), val(file_tag), file("${file_tag}_1.fq.gz"), file("${file_tag}_2.fq.gz") into readPairs0

        script:
        '''
        set -o pipefail
        samtools collate -uOn 128 !{infile} tmp_!{infile.baseName} | samtools fastq -1 !{infile.baseName}!{params.suffix1}.!{params.fastq_ext} -2 !{infile.baseName}!{params.suffix2}.!{params.fastq_ext} -
        '''
    }

    process ADAPTER_TRIMMING {
        tag { file_tag + rg }
        cpus params.cpu_trim
        memory "${params.mem_QC}GB"

        input:
        set val(file_tag), val(rg), path(pair1), path(pair2) from readPairs2

        output:
        set val(file_tag), val(rg), file("${file_tag}${rg}*val_1.fq.gz"), file("${file_tag}${rg}*val_2.fq.gz") into readPairs3
        file("*_fastqc.zip") into fastqc_postpairs
        file("*trimming_report.txt") into trimming_reports

        publishDir "${params.output_folder}/QC/adapter_trimming", mode: 'copy', pattern: '{*report.txt,*fastqc.zip}'

        script:
        '''
        cpu_tg=$((params.cpu_trim - 1))
        cpu_tg2=$(echo "$cpu_tg/3.5" | bc -l)
        cpu_tg3=$(python - <<'PY'
		import math,sys
		v=float(sys.argv[1])
		print(int(math.ceil(v)))
		PY
		$cpu_tg2)
        if [ -n "${pair2}" ] && [ "$(basename ${pair2})" != "NO_fastq2" ]; then
            opts="--paired"
        else
            opts=""
        fi
        trim_galore ${opts} --fastqc --gzip --basename ${file_tag}${rg} -j ${cpu_tg3} ${pair1} ${pair2}
        if [ ! -L NO_fastq2 ]; then
            mv ${file_tag}${rg}_trimmed.fq.gz ${file_tag}${rg}_val_1.fq.gz
            touch ${file_tag}${rg}_val_2.fq.gz
        fi
        '''
    }


	process FASTQC_PRETRIM {
		tag { file_tag }
		cpus params.cpu
		memory "${params.mem_QC}GB"

		input:
		set val(file_tag), val(rg), path(pair1), path(pair2) from readPairs

		output:
		file("*_pretrim_fastqc.zip") into fastqc_pairs

		publishDir "${params.output_folder}/QC/fastq", mode: 'copy', pattern: '{*fastqc.zip}'

		script:
		'''
		basename1=$(basename ${pair1} .${params.fastq_ext})
		if [ -n "${pair2}" ] && [ "$(basename ${pair2})" != "NO_fastq2" ]; then
			fastqc -t ${task.cpus} ${pair1} ${pair2}
			mv ${basename1}_fastqc.zip ${file_tag}${params.suffix1}${rg}_pretrim_fastqc.zip
		else
			fastqc -t ${task.cpus} ${pair1}
			mv ${basename1}_fastqc.zip ${file_tag}${params.suffix1}${rg}_pretrim_fastqc.zip
		fi
    '''
}

	process ALIGNMENT {
		tag { file_tag }
		cpus params.cpu
		memory "${params.mem}G"

		input:
		set val(file_tag), val(rg), path(pair1), path(pair2) from readPairs_align
		file ref from ref.collect()
		file gtf from gtf

		output:
		set val(file_tag), val(rg), file("${file_tag}.bam"), file("${file_tag}.bam.bai") into bam_files
		file("*Log*") into align_out
		set val(file_tag), file("*SJ.out.junction") into SJ_out
		file("*SJ.out.tab") into SJ_out_others

		// conditional publish logic is achieved post-run by where outputs are copied in shell commands
		script:
		'''
		align_threads=$(( ${params.cpu} / 2 ))
		sort_threads=$(( ${params.cpu} / 2 - 1 ))
		sort_mem=$(( ${params.mem} / 4 ))
		input_f1="${pair1}"
		rgline="ID:${file_tag} SM:${file_tag} ${params.RG}"
		if [ -n "${pair2}" ] && [ "$(basename ${pair2})" != "NO_fastq2" ]; then
			pairs="${pair1} ${pair2}"
		else
			pairs="${pair1}"
		fi
		STAR --outSAMattrRGline "${rgline}" --outSAMmapqUnique ${params.STAR_mapqUnique} --chimSegmentMin 12 --chimJunctionOverhangMin 12 \
		--chimSegmentReadGapMax 3 --alignSJDBoverhangMin 10 --alignMatesGapMax 100000 --alignIntronMax 100000 \
		--alignSJstitchMismatchNmax 5 -1 5 5 --outSAMstrandField intronMotif --chimMultimapScoreRange 10 --chimMultimapNmax 10 \
		--chimNonchimScoreDropMin 10 --peOverlapNbasesMin 12 --peOverlapMMp 0.1 --chimOutJunctionFormat 1 --twopassMode Basic \
		--outReadsUnmapped None --runThreadN ${align_threads} --genomeDir . --sjdbGTFfile ${gtf} --readFilesCommand zcat \
		--readFilesIn ${pairs} --outStd SAM | samblaster --addMateTags | sambamba view -S -f bam -l 0 /dev/stdin | sambamba sort -t ${sort_threads} -m ${sort_mem}G --tmpdir=${file_tag}_tmp -o ${file_tag}.bam /dev/stdin
		mv Chimeric.out.junction STAR.${file_tag}.Chimeric.SJ.out.junction || true
		mv SJ.out.tab STAR.${file_tag}.SJ.out.tab || true
		mv Log.final.out STAR.${file_tag}.Log.final.out || true
		mv Log.out STAR.${file_tag}.Log.out || true
		mv Log.progress.out STAR.${file_tag}.Log.progress.out || true
		mv Log.std.out STAR.${file_tag}.Log.std.out || true
		'''
}

    process SPLICE_JUNCT_TRIM {
        tag { file_tag }
        cpus params.cpu_gatk
        memory "${params.mem}G"

        input:
        set val(file_tag), val(rg), path(bam), path(bai) from bam_files
        file fasta_ref from fasta_ref
        file fasta_ref_fai from fasta_ref_fai
        file fasta_ref_dict from fasta_ref_dict

        output:
        set val(file_tag_new), val(rg), file("${file_tag_new}.bam"), file("${file_tag_new}.bam.bai") into bam_files2

        script:
        '''
        file_tag_new=${file_tag}_split
        gatk SplitNCigarReads --java-options "-Xmx${params.mem}G" -R ${fasta_ref} -I ${bam} -O ${file_tag_new}.bam
        mv ${file_tag_new}.bai ${file_tag_new}.bam.bai || true
        '''
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
        set val(file_tag), val(rg), path("${file_tag}.bam"), path("${file_tag}.bam.bai") from bam_files2
        file known_snps from known_snps
        file known_snps_index from known_snps_index
        file known_indels from known_indels
        file known_indels_index from known_indels_index
        file fasta_ref from fasta_ref
        file fasta_ref_fai from fasta_ref_fai
        file fasta_ref_dict from fasta_ref_dict

        output:
        file("*_recal.table") into recal_table_files
        file("*plots.pdf") into recal_plots_files
        set val(file_tag_new), val(rg), file("${file_tag_new}.bam"), file("${file_tag_new}.bam.bai") into recal_bam_files

        script:
        '''
        file_tag_new=${file_tag}_BQSRecalibrated
        gatk BaseRecalibrator --java-options "-Xmx${params.mem}G" -R ${fasta_ref} -I ${file_tag}.bam --known-sites ${known_snps} --known-sites ${known_indels} -O ${file_tag}_recal.table
        gatk ApplyBQSR --java-options "-Xmx${params.mem}G" -R ${fasta_ref} -I ${file_tag}.bam --bqsr-recal-file ${file_tag}_recal.table -O ${file_tag_new}.bam
        gatk BaseRecalibrator --java-options "-Xmx${params.mem}G" -R ${fasta_ref} -I ${file_tag_new}.bam --known-sites ${known_snps} --known-sites ${known_indels} -O ${file_tag_new}_recal.table
        gatk AnalyzeCovariates --java-options "-Xmx${params.mem}G" -before ${file_tag}_recal.table -after ${file_tag_new}_recal.table -plots ${file_tag_new}_recalibration_plots.pdf
        mv ${file_tag_new}.bai ${file_tag_new}.bam.bai || true
        '''
    }



	process RSEQC { //(read distribution, clipping, junction saturation)
		tag { file_tag }
		cpus 1
		memory "${params.mem_QC}GB"

		input:
			set val(file_tag), val(rg), path(bam), path(bai) from recal_bam_files4QC
			file bed from bed

		output:
			file("${file_tag}_readdist.txt") into rseqc_files
			file("*clipping*") into rseqc_clip_files
			file("*jun_saturation*") into rseqc_jsat_files

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
		set val(file_tag), val(rg), path(bam), path(bai) from recal_bam_files4QCsplit
		file bed from bed

		output:
		file("*readdist.txt") into rseqc_files_split

		publishDir "${params.output_folder}/QC/bam", mode: 'copy'

		script:
		'''
		basename=$(basename ${bam})
		samtools split ${bam} -f "%*_%!.%."
		for f in ${basename}_*.bam; do
			read_distribution.py -i $f -r ${bed} > ${f%.bam}_readdist.txt
		done
		'''
}

	process QUANTIFICATION {
		tag { file_tag }
		cpus params.cpu
		memory { (params.sjtrim || params.recalibration) ? "${params.mem}G" : "${params.mem_QC}G" }()

		input:
		set val(file_tag), val(rg), path(bam), path(bai) from recal_bam_files4quant
		file gtf from gtf

		output:
		file("${file_tag}_count.txt") into htseq_files

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

	process MULTIQC_PRETRIM {
		tag { "all" }
		cpus 1
		memory "${params.mem_QC}GB"

		input:
		file fastqc1 from fastqc_pairs.collect()
		file multiqc_config from ch_config_for_multiqc

		output:
		file("multiqc_pretrim_report.html") into multiqc_pre
		file("multiqc_pretrim_report_data") into multiqc_pre_data

		publishDir "${params.output_folder}/QC", mode: 'copy'

		script:
		'''
		if [ "$(basename ${multiqc_config})" == "NO_FILE" ]; then
			opt=""
		else
			opt="--config ${multiqc_config}"
		fi
		for f in $(find . -name "*_pretrim_fastqc.zip" -type l); do cp --remove-destination $(readlink $f) $f || true; done
		multiqc . -n multiqc_pretrim_report.html -m fastqc ${opt} --comment "RNA-seq Pre-trimming QC report"
		'''
	}

	process MULTIQC_POSTTRIM {
		tag { "all" }
		cpus 1
		memory "${params.mem_QC}GB"

		input:
		file STAR from align_out.collect()
		file htseq from htseq_files.collect()
		file rseqc_clip from rseqc_clip_files.collect()
		file rseqc from rseqc_files.collect()
		file rseqc_jsat from rseqc_jsat_files.collect()
		file trim from trimming_reports.collect().ifEmpty([])
		file fastqcpost from fastqc_postpairs.collect().ifEmpty([])
		file rseqc_split from rseqc_files_split.collect().ifEmpty([])
		file multiqc_config from ch_config_for_multiqc

		output:
		file("multiqc_posttrim_report.html") into multiqc_post
		file("multiqc_posttrim_report_data") into multiqc_post_data

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

#TO-DO
