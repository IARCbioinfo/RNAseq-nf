#! /usr/bin/env nextflow
// usage : ./RNAseq.nf --input_folder input/ --cpu 8 --mem 32 --ref hg19.fasta 
/*
vim: syntax=groovy
-*- mode: groovy;-*- */

// requirement:
// - fastQC
// - STAR

//default values
params.help         = null
params.input_folder = '.'
params.cpu          = 8
params.mem          = 32
params.fastq_ext    = "fq.gz"
params.suffix1      = "_1"
params.suffix2      = "_2"
params.gendir       = "ref"
params.fasta_ref    = "ref.fa"
params.output_folder   = "."
params.annot_gtf    = "Homo_sapiens.GRCh38.79.gtf"
params.annot_gff    = "Homo_sapiens.GRCh38.79.gff"
params.GATK_folder  = "GATK"
params.GATK_bundle  = "GATK_bundle"
params.intervals    = "intervals.bed"
params.RG           = "PL:ILLUMINA"
params.sjtrim       = "false"
params.bqsr         = "false"

if (params.help) {
    log.info ''
    log.info '-------------------------------------------------------------'
    log.info 'NEXTFLOW RNASEQ ANALYSIS PIPELINE'
    log.info '-------------------------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info 'nextflow run RNAseq.nf --input_folder input/ --gendir ref/ [--cpu 8] [--mem 32] [--suffix1 _1] [--suffix2 _2] [--output_folder output/]'
    log.info ''
    log.info 'Mandatory arguments:'
    log.info '    --input_folder   FOLDER                  Folder containing BAM or fastq files to be aligned.'
    log.info 'Optional arguments:'
    log.info '    --cpu          INTEGER                 Number of cpu used by bwa mem and sambamba (default: 8).'
    log.info '    --mem          INTEGER                 Size of memory used by sambamba (in GB) (default: 32).'
    log.info '    --fastq_ext        STRING                Extension of fastq files (default : fq.gz)'
    log.info '    --suffix1        STRING                Suffix of fastq files 1 (default : _1)'
    log.info '    --suffix2        STRING                Suffix of fastq files 2 (default : _2)'
    log.info '    --gendir     STRING                Folder with reference genome and STAR index (default: ref).'
    log.info '    --output_folder     STRING                Output folder (default: results_RNAseq).'
    log.info ''
    exit 1
}

//read files
keys1 = file(params.input_folder).listFiles().findAll { it.name ==~ /.*${params.suffix1}.${params.fastq_ext}/ }.collect { it.getName() }
                                                                                                               .collect { it.replace("${params.suffix1}.${params.fastq_ext}",'') }
keys2 = file(params.input_folder).listFiles().findAll { it.name ==~ /.*${params.suffix2}.${params.fastq_ext}/ }.collect { it.getName() }
                                                                                                               .collect { it.replace("${params.suffix2}.${params.fastq_ext}",'') }

if ( !(keys1.containsAll(keys2)) || !(keys2.containsAll(keys1)) ) {println "\n ERROR : There is not at least one fastq without its mate, please check your fastq files."; System.exit(0)}

println keys1

// Gather files ending with _1 suffix
reads1 = Channel
    .fromPath( params.input_folder+'/*'+params.suffix1+'.'+params.fastq_ext )
    .map {  path -> [ path.name.replace("${params.suffix1}.${params.fastq_ext}",""), path ] }

// Gather files ending with _2 suffix
reads2 = Channel
    .fromPath( params.input_folder+'/*'+params.suffix2+'.'+params.fastq_ext )
    .map {  path -> [ path.name.replace("${params.suffix2}.${params.fastq_ext}",""), path ] }

// Match the pairs on two channels having the same 'key' (name) and emit a new pair containing the expected files
readPairs = reads1
    .phase(reads2)
    .map { pair1, pair2 -> [ pair1[1], pair2[1] ] }

println reads1
        
// pre-trimming QC
process fastqc_pretrim {
	cpus params.cpu
        memory params.mem+'GB'    
        tag { file_tag }
        
        input:
        file pairs from readPairs
        output:
	file pairs into readPairs2
	val(file_tag) into filetag
	file("${file_tag}${params.suffix1}_fastqc.zip") into fastqc_pair1
	file("${file_tag}${params.suffix2}_fastqc.zip") into fastqc_pair2

	shell:
        file_tag = pairs[0].name.replace("${params.suffix1}.${params.fastq_ext}","")
        '''
	fastqc -t !{task.cpus} !{pairs[0]} !{pairs[1]}
        '''
}

process multiqc_pretrim {
    cpus params.cpu
    memory '1G'
    tag { "multiqc pretrim"}
        
    input:
    val(file_tag) from filetag
    file pairs2 from readPairs2
    file fastqc1 from fastqc_pair1.toList()
    file fastqc2 from fastqc_pair2.toList()
    
    output:
    val(file_tag) into filetag2
    file("multiqc_pretrim_report.html") into multiqc_report
    file("multiqc_pretrim_report_data") into multiqc_data
    file pairs2 into readPairs3

    publishDir params.output_folder, mode: 'copy', pattern: 'multiqc_pretrim_report*'

    shell:
    '''
    for f in $(find *fastqc.zip -type l);do cp --remove-destination $(readlink $f) $f;done;
    multiqc . -n multiqc_pretrim_report.html
    '''
}


// adapter sequence trimming and post trimming QC
process adapter_trimming {
            cpus params.cpu
            memory params.mem+'G'
            tag { file_tag }
	    
            input:
	    val(file_tag) from filetag2
	    file pairs3 from readPairs3
            output:
            val(file_tag) into filetag3
	    file("${file_tag}*val*.fq.gz") into readPairs4
	    file("${file_tag}${params.suffix1}_val${params.suffix1}_fastqc.zip") into fastqc_postpair1
	    file("${file_tag}${params.suffix2}_val${params.suffix2}_fastqc.zip") into fastqc_postpair2
	    file("${file_tag}*trimming_report.txt") into trimming_reports
	    
	    publishDir params.output_folder, mode: 'copy', pattern: '{*report.txt}'
	    
            shell:
            '''
	    trim_galore --paired --fastqc !{pairs3[0]} !{pairs3[1]}
            '''
}


process multiqc_posttrim {
    cpus params.cpu
    memory '1G'
    tag { "multiqc posttrim"}
        
    input:
    val(file_tag) from filetag3
    file pairs4 from readPairs4
    file fastqc1 from fastqc_postpair1.toList()
    file fastqc2 from fastqc_postpair2.toList()
        
    output:
    val(file_tag) into filetag4
    file("multiqc_posttrim_report.html") into multiqc_post
    file("multiqc_posttrim_report_data") into multiqc_post_data
    file pairs4 into readPairs5

    publishDir params.output_folder, mode: 'copy', pattern: 'multiqc_posttrim*'

    shell:
    '''
    for f in $(find *fastqc.zip -type l);do cp --remove-destination $(readlink $f) $f;done;
    multiqc -n multiqc_posttrim_report.html $PWD/
    '''
}


//Mapping, mark duplicates and sorting
process alignment {
      cpus params.cpu
      memory params.mem+'G'
      tag { file_tag }
      
      input:
      val(file_tag) from filetag4
      file pairs5  from readPairs5
            
      output:
      val(file_tag) into filetag5
      file("${file_tag}.bam") into bam_files
      file("${file_tag}.bam.bai") into bai_files
      file("STAR.${file_tag}.Log.final.out") into STAR_out
      publishDir params.output_folder, mode: 'copy', pattern: "STAR.${file_tag}.Log.final.out"
            
      shell:
      STAR_threads = params.cpu.intdiv(2) - 1
      sort_threads = params.cpu.intdiv(2) - 1
      sort_mem     = params.mem.intdiv(4)
      '''
      STAR --outSAMattrRGline "ID:!{file_tag}\\tSM:!{file_tag}\\t!{params.RG}" --chimSegmentMin 12 --chimJunctionOverhangMin 12 --chimSegmentReadGapMax 3 --alignSJDBoverhangMin 10 --alignMatesGapMax 200000 --alignIntronMax 200000 --alignSJstitchMismatchNmax 5 -1 5 5 --twopassMode Basic --runThreadN !{STAR_threads} --genomeDir !{params.gendir} --sjdbGTFfile !{params.annot_gtf} --readFilesCommand zcat --readFilesIn !{pairs5[0]} !{pairs5[1]} --outStd SAM | samblaster --addMateTags | sambamba view -S -f bam -l 0 /dev/stdin | sambamba sort -t !{sort_threads} -m !{sort_mem}G --tmpdir=!{file_tag}_tmp -o !{file_tag}.bam /dev/stdin
      mv Log.final.out STAR.!{file_tag}.Log.final.out
      '''
}

//Splice junctions trimming
if(params.sjtrim != "false"){
   process splice_junct_trim {
      cpus params.cpu
      memory params.mem+'G'
      tag { file_tag }
      
      input:
      val(file_tag) from filetag5
      file bam  from bam_files
      file bai  from bai_files
            
      output:
      val(file_tag) into filetag6
      file("${file_tag}_split.bam") into bam_files2
      file("${file_tag}_split.bam.bai") into bai_files2
            
      shell:
      '''
      java -jar !{params.GATK_folder}/GenomeAnalysisTK.jar -T SplitNCigarReads -R !{params.fasta_ref} -I !{bam} -o !{file_tag}_split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS
      '''
   }
}else{
   process no_splice_junct_trim {
      cpus '1'
      memory '100M'
      tag { file_tag }
      
      input:
      val(file_tag) from filetag5
      file bam  from bam_files
      file bai  from bai_files
            
      output:
      val(file_tag) into filetag6
      file bam into bam_files2
      file bai into bai_files2
            
      shell:
      '''
      touch !{file_tag}.bam
      '''
   }
}


//BQSrecalibration
if(params.bqsr != "false"){
   process base_quality_score_recalibration {
    	cpus params.cpu
	memory params.mem+'G'
    	tag { file_tag }
        
    	input:
    	val(file_tag) from filetag6
	file bam from bam_files2
    	file bai from bai_files2
    	output:
	val(file_tag) into filetag7
    	file("${file_tag}_recal.table") into recal_table_files
    	file("${file_tag}_post_recal.table") into recal_table_post_files
    	file("${file_tag}_recalibration_plots.pdf") into recal_plots_files
    	file("${file_tag}.bam") into recal_bam_files
    	file("${file_tag}.bam.bai") into recal_bai_files
    	publishDir params.output_folder, mode: 'move'

    	shell:
    	'''
    	indelsvcf=(`ls !{params.GATK_bundle}/*indels*.vcf* | grep -v ".tbi" | grep -v ".idx"`)
    	dbsnpvcfs=(`ls !{params.GATK_bundle}/*dbsnp*.vcf* | grep -v ".tbi" | grep -v ".idx"`)
    	dbsnpvcf=${dbsnpvcfs[@]:(-1)}
    	knownSitescom=''
    	for ll in $indelsvcf; do knownSitescom=$knownSitescom' -knownSites '$ll; done
    	knownSitescom=$knownSitescom' -knownSites '$dbsnpvcf
    	java -jar !{params.GATK_folder}/GenomeAnalysisTK.jar -T BaseRecalibrator -nct !{params.cpu} -R !{params.fasta_ref} -I !{file_tag}.bam $knownSitescom -L !{params.intervals} -o !{file_tag}_recal.table
    	java -jar !{params.GATK_folder}/GenomeAnalysisTK.jar -T BaseRecalibrator -nct !{params.cpu} -R !{params.fasta_ref} -I !{file_tag}.bam $knownSitescom -BQSR !{file_tag}_recal.table -L !{params.intervals} -o !{file_tag}_post_recal.table		
    	java -jar !{params.GATK_folder}/GenomeAnalysisTK.jar -T AnalyzeCovariates -R !{params.fasta_ref} -before !{file_tag}_recal.table -after !{file_tag}_post_recal.table -plots !{file_tag}_recalibration_plots.pdf	
    	java -jar !{params.GATK_folder}/GenomeAnalysisTK.jar -T PrintReads -nct !{params.cpu} -R !{params.fasta_ref} -I !{file_tag}.bam -BQSR !{file_tag}_recal.table -L !{params.intervals} -o !{file_tag}.bam
    	mv !{file_tag}.bai !{file_tag}.bam.bai
    	'''
   }
}else{
 process no_BQSR {
      cpus '1'
      memory '100M'
      tag { file_tag }
      
      input:
      val(file_tag) from filetag6
      file bam  from bam_files2
      file bai  from bai_files2
            
      output:
      val(file_tag) into filetag7
      file bam into recal_bam_files
      file bai into recal_bai_files
            
      shell:
      '''
      touch !{file_tag}.bam
      '''
   }
}

//Quantification
process quantification{
    	cpus params.cpu
	memory params.mem+'G'
    	tag { file_tag }
        
    	input:
    	val(file_tag) from filetag7
	file bam from recal_bam_files
    	file bai from recal_bai_files
    	output:
	val(file_tag) into filetag8
	file bam into recal_bam_files2
    	file bai into recal_bai_files2
    	file("*.txt") into htseq_files
    	publishDir params.output_folder, mode: 'move'

    	shell:
    	'''
	htseq-count -r pos -s yes -f bam !{file_tag}.bam !{params.annot_gff}
    	'''
}
