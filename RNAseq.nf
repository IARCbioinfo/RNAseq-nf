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
params.ref          = 'hg19.fasta'
params.genparams    = 'genomeParameters.txt'
params.GTF          = 'hg19.gtf'
params.cpu          = 8
params.mem          = 32
params.fastq_ext    = "fq.gz"
params.suffix1      = "_1"
params.suffix2      = "_2"
params.gendir       = "ref"
params.output_folder   = "."

if (params.help) {
    log.info ''
    log.info '-------------------------------------------------------------'
    log.info 'NEXTFLOW RNASEQ ANALYSIS PIPELINE'
    log.info '-------------------------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info 'nextflow run RNAseq.nf --input_folder input/ --ref hg19.fasta [--cpu 8] [--mem 32] [--suffix1 _1] [--suffix2 _2] [--output_folder output/]'
    log.info ''
    log.info 'Mandatory arguments:'
    log.info '    --input_folder   FOLDER                  Folder containing BAM or fastq files to be aligned.'
    log.info '    --ref          FILE                    Reference fasta file (with index).'
    log.info 'Optional arguments:'
    log.info '    --cpu          INTEGER                 Number of cpu used by bwa mem and sambamba (default: 8).'
    log.info '    --mem          INTEGER                 Size of memory used by sambamba (in GB) (default: 32).'
    log.info '    --suffix1        STRING                Suffix of fastq files 1 (default : _1)'
    log.info '    --suffix2        STRING                Suffix of fastq files 2 (default : _2)'
    log.info '    --output_folder     STRING                Output folder (default: results_RNAseq).'
    log.info ''
    exit 1
}

//read files
ref     = file( params.ref )

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
	file('${file_tag}${params.suffix1}_fastqc.html') into fastqc_pair1
	file('${file_tag}${params.suffix2}_fastqc.html') into fastqc_pair2
	
	publishDir params.output_folder, mode: 'move', pattern: '*.html'

	shell:
        file_tag = pairs[0].name.replace("${params.suffix1}.${params.fastq_ext}","")
        '''
	fastqc -t !{task.cpus} !{pairs[0]} !{pairs[1]}
        '''
}

// adapter sequence trimming and post trimming QC
process trimming {
            cpus params.cpu
            memory params.mem+'G'
            tag { file_tag }
	    
            input:
	    val(file_tag) from filetag
	    file pairs2 from readPairs2
            output:
            val(file_tag) into filetag2
	    file("${file_tag}_*_trimmed.fq.gz") into readPairs3
	    file('${file_tag}*.html') into fastqc_posttrim
	
	    publishDir params.output_folder, mode: 'move', pattern: '*.html'
	    
            shell:
            '''
	    trim_galore --paired --fastqc !{pairs2[0]} !{pairs2[1]}
            '''
}

// alignment
process alignment {
      cpus params.cpu
      memory params.mem+'G'
      perJobMemLimit = true
      tag { file_tag }
      
      input:
      file pairs3  from readPairs3
      val(file_tag) from filetag2
      
      output:
      val(file_tag) into filetag3
      file("${file_tag}*.bam") into bam_files
      file("${file_tag}*.bai") into bai_files
      publishDir params.output_folder, mode: 'move'
      
      shell:
      '''
      STAR --runThreadN !{task.cpus} --genomeDir !{params.gendir} --readFilesCommand zcat --readFilesIn !{pairs3[0]} !{pairs3[1]} --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts
      '''
}

//2nd pass mapping (TODO)
//process 1:
//- merge SJ.out.tab files from all runs
//- Filter junctions by removing likely false positives, e.g.  junctions in the mitochondrion genome, or non-canonical junctions supported by a few reads.  If using annotations, only novel junctions need to be considered here, since annotated junctions will be re-used in the 2nd pass anyway.
//-  Use the filtered list of junctions from the 1st pass with --sjdbFileChrStartEnd option, together with annotations (via --sjdbGTFfile option) to generate the new genome indices for the 2nd pass mapping.  This needs to be done only once for all samples.

//process 2: Run the 2nd pass mapping for all samples with the new genome index


//Splice junctions trimming


//BQSrecalibration