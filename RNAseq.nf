#! /usr/bin/env nextflow

// vim: syntax=groovy -*- mode: groovy;-*-

// Copyright (C) 2017 IARC/WHO

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

params.input_folder = null
params.ref_folder   = null
params.starfusion_folder = null
params.gtf          = null
params.bed          = null

params.cpu          = 4
params.mem          = 50
params.mem_QC       = 2
params.fastq_ext    = "fq.gz"
params.suffix1      = "_1"
params.suffix2      = "_2"
params.output_folder= "."
params.ref          = "ref.fa"
params.GATK_jar     = "GenomeAnalysisTK.jar"
params.GATK_bundle  = "GATK_bundle"
params.RG           = "PL:ILLUMINA"
params.stranded     = "no"
params.hisat2_idx   = "genome_tran"
params.clustering_n = 500 
params.clustering_t = "vst"
params.clustering_c = "hc"
params.clustering_l = "complete"

params.sjtrim       = null
params.recalibration= null
params.hisat2       = null
params.clustering   = null

params.cutadapt     = null

params.htseq_maxreads= null //default value of htseq-count is 30000000
params.help         = null


log.info "" 
log.info "--------------------------------------------------------"
log.info "  RNAseq-nf 1.0.0: alignment, QC, and reads counting workflow for RNA sequencing "
log.info "--------------------------------------------------------"
log.info "Copyright (C) IARC/WHO"
log.info "This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE"
log.info "This is free software, and you are welcome to redistribute it"
log.info "under certain conditions; see LICENSE for details."
log.info "--------------------------------------------------------"
log.info ""


if (params.help) {
    log.info '-------------------------------------------------------------'
    log.info ' USAGE  '
    log.info '-------------------------------------------------------------'
    log.info ''
    log.info 'nextflow run iarcbioinfo/RNAseq.nf [-with-docker] --input_folder input/ --ref_folder ref/ [OPTIONS]'
    log.info ''
    log.info 'Mandatory arguments:'
    log.info '    --input_folder   FOLDER                  Folder containing BAM or fastq files to be aligned.'
    log.info '    --ref_folder          FOLDER                   Folder with genome reference files (with index).'
    log.info '    --starfusion_folder          FOLDER                   Folder with STAR-Fusion reference files.'
    log.info '    --gtf          FILE                    Annotation file.'
    log.info '    --bed        STRING                bed file with interval list'
    log.info ""
    log.info 'Optional arguments:'
    log.info '    --ref          FILE                    Reference fasta file (with index) for splice junction trimming and base recalibration.'
    log.info '    --output_folder     STRING                Output folder (default: results_alignment).'
    log.info '    --cpu          INTEGER                 Number of cpu used by bwa mem and sambamba (default: 8).'
    log.info '    --mem          INTEGER                 Size of memory used for mapping (in GB) (default: 32).'
    log.info '    --mem_QC     INTEGER                 Size of memory used for QC and cutadapt (in GB) (default: 32).'
    log.info '    --RG           STRING                  Samtools read group specification (default : PL:ILLUMINA).'
    log.info '    --fastq_ext        STRING                Extension of fastq files (default : fq.gz)'
    log.info '    --suffix1        STRING                Suffix of fastq files 1 (default : _1)'
    log.info '    --suffix2        STRING                Suffix of fastq files 2 (default : _2)'
    log.info '    --GATK_bundle        STRING                path to GATK bundle files (default : .)'
    log.info '    --GATK_jar        STRING                path to GATK GenomeAnalysisTK.jar file (default : .)'
    log.info '    --stranded        STRING                are reads stranded? (default : no; alternatives : yes, r)'
    log.info '    --hisat2_idx        STRING                hisat2 index file prefix (default : genome_tran)'
    log.info '    --clustering_n        INTEGER               number of genes to use for clustering  (default : 500)'
    log.info '    --clustering_t        STRING                count transformation method; rld, vst, or auto (default : vst)'
    log.info '    --clustering_c        STRING                clustering algorithm to be passed to ConsensusClusterPlus (default : hc)'
    log.info '    --clustering_l        STRING                method for hierarchical clustering to be passed to ConsensusClusterPlus (default : complete)'
    log.info ''
    log.info 'Flags:'
    log.info '    --sjtrim                    enable splice junction trimming'
    log.info '    --recalibration                    performs base quality score recalibration (GATK)'
    log.info '    --hisat2                    use hisat2 instead of STAR for reads mapping'
    log.info '    --clustering                    perform unsupervised analysis of read counts'
    log.info '    --cutadapt                  perform adapter sequence trimming'
    log.info ''
    exit 0

params.clustering_n = 500 
params.clustering_t = "vst"
params.clustering_c = "hc"
params.clustering_l = "complete"

params.sjtrim       = null
params.recalibration = null
params.hisat2       = null
params.clustering   = null

params.htseq_maxreads = null //default value of htseq-count is 30000000
params.help         = null

}else {
  /* Software information */
  log.info "input_folder = ${params.input_folder}"
  log.info "ref          = ${params.ref}"
  log.info "cpu          = ${params.cpu}"
  log.info "mem          = ${params.mem}"
  log.info "fastq_ext    = ${params.fastq_ext}"
  log.info "suffix1      = ${params.suffix1}"
  log.info "suffix2      = ${params.suffix2}"
  log.info "output_folder= ${params.output_folder}"
  log.info "bed          = ${params.bed}"
  log.info "GATK_bundle  = ${params.GATK_bundle}"
  log.info "GATK_jar     = ${params.GATK_jar}"
  log.info "mem_QC       = ${params.mem_QC}"
  log.info "ref_folder   = ${params.ref_folder}"
  log.info "starfusion_folder = ${params.starfusion_folder}"
  log.info "gtf          = ${params.gtf}"
  log.info "RG           = ${params.RG}"
  log.info "stranded     = ${params.stranded}"
  log.info "hisat2_idx   = ${params.hisat2_idx}"
  log.info "sjtrim       = ${params.sjtrim}"
  log.info "hisat2       = ${params.hisat2}"
  log.info "clustering   = ${params.clustering}"
  log.info "htseq_maxreads=${params.htseq_maxreads}"
  log.info "recalibration= ${params.recalibration}"

  log.info "clustering_n = ${params.clustering_n}"
  log.info "clustering_t = ${params.clustering_t}"
  log.info "clustering_c = ${params.clustering_c}"
  log.info "clustering_l = ${params.clustering_l}"
  log.info "help=${params.help}"
}

//read ref files
if(params.hisat2){
	ref_1  = Channel.fromPath(params.ref_folder + '/' + params.hisat2_idx + '.1.ht2')
	ref_2  = Channel.fromPath(params.ref_folder + '/' + params.hisat2_idx + '.2.ht2')
	ref_3  = Channel.fromPath(params.ref_folder + '/' + params.hisat2_idx + '.3.ht2')
	ref_4  = Channel.fromPath(params.ref_folder + '/' + params.hisat2_idx + '.4.ht2')
	ref_5  = Channel.fromPath(params.ref_folder + '/' + params.hisat2_idx + '.5.ht2')
	ref_6  = Channel.fromPath(params.ref_folder + '/' + params.hisat2_idx + '.6.ht2')
	ref_7  = Channel.fromPath(params.ref_folder + '/' + params.hisat2_idx + '.7.ht2')
	ref_8  = Channel.fromPath(params.ref_folder + '/' + params.hisat2_idx + '.8.ht2')
	ref    = ref_1.concat( ref_2,ref_3,ref_4,ref_5,ref_6,ref_7,ref_8)
}else{
	ref_1  = Channel.fromPath(params.ref_folder +'/chrStart.txt')
	ref_2  = Channel.fromPath(params.ref_folder +'/chrNameLength.txt')
	ref_3  = Channel.fromPath(params.ref_folder +'/chrName.txt')
	ref_4  = Channel.fromPath(params.ref_folder +'/chrLength.txt')
	ref_5  = Channel.fromPath(params.ref_folder +'/exonGeTrInfo.tab')
	ref_6  = Channel.fromPath(params.ref_folder +'/exonInfo.tab')
	ref_7  = Channel.fromPath(params.ref_folder +'/geneInfo.tab')
	ref_8  = Channel.fromPath(params.ref_folder +'/Genome')
	ref_9  = Channel.fromPath(params.ref_folder +'/genomeParameters.txt')
	ref_10 = Channel.fromPath(params.ref_folder +'/SA')
	ref_11 = Channel.fromPath(params.ref_folder +'/SAindex')
	ref_12 = Channel.fromPath(params.ref_folder +'/sjdbInfo.txt')
	ref_13 = Channel.fromPath(params.ref_folder +'/transcriptInfo.tab')
	ref_14 = Channel.fromPath(params.ref_folder +'/sjdbList.fromGTF.out.tab')
	ref_15 = Channel.fromPath(params.ref_folder +'/sjdbList.out.tab')
	ref    = ref_1.concat( ref_2,ref_3,ref_4,ref_5,ref_6,ref_7,ref_8,ref_9,ref_10,ref_11,ref_12,ref_13,ref_14,ref_15)
	ref.into { ref_align; ref_fusion }
	blast_ava = file("${params.starfusion_folder}/blast_pairs.idx")
	gtf_gs = file("${params.starfusion_folder}/ref_annot.gtf.gene_spans")
}

gtf    = file(params.gtf)
bed    = file(params.bed)

//read files
mode = 'fastq'
if (file(params.input_folder).listFiles().findAll { it.name ==~ /.*${params.fastq_ext}/ }.size() > 0){
    println "fastq files found, proceed with alignment"
}else{
    if (file(params.input_folder).listFiles().findAll { it.name ==~ /.*bam/ }.size() > 0){
        println "BAM files found, proceed with realignment"; mode ='bam'; files = Channel.fromPath( params.input_folder+'/*.bam' )
    }else{
        println "ERROR: input folder contains no fastq nor BAM files"; System.exit(0)
    }
}

if(mode=='bam'){
    process bam2fastq {
        cpus '1'
        memory params.mem_QC+'G'
        tag { file_tag }
        
        input:
        file infile from files
     
        output:
	set file("${file_tag}_1.fq.gz"), file("${file_tag}_2.fq.gz")  into readPairs

        shell:
	file_tag = infile.baseName
		
        '''
        set -o pipefail
        samtools collate -uOn 128 !{file_tag}.bam tmp_!{file_tag} | samtools fastq -1 !{file_tag}_1.fq -2 !{file_tag}_2.fq -
	gzip !{file_tag}_1.fq
	gzip !{file_tag}_2.fq
        '''
    }
}else{
if(mode=='fastq'){
    println "fastq mode"
    
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
}
}

// pre-trimming QC
process fastqc_pretrim {
	cpus params.cpu
        memory params.mem_QC+'GB'    
        tag { file_tag }
        
        input:
        file pairs from readPairs
	
        output:
	set file("${file_tag}${params.suffix1}_pretrim_fastqc.zip"), file("${file_tag}${params.suffix2}_pretrim_fastqc.zip") into fastqc_pairs
	set val(file_tag), file(pairs) into readPairs3
	
	publishDir "${params.output_folder}/QC/fastq", mode: 'copy', pattern: '{*fastqc.zip}'

	shell:
        file_tag = pairs[0].name.replace("${params.suffix1}.${params.fastq_ext}","")
        '''
	fastqc -t !{task.cpus} !{pairs[0]} !{pairs[1]}
	mv !{file_tag}!{params.suffix1}_fastqc.zip !{file_tag}!{params.suffix1}_pretrim_fastqc.zip
	mv !{file_tag}!{params.suffix2}_fastqc.zip !{file_tag}!{params.suffix2}_pretrim_fastqc.zip 
        '''
}

// adapter sequence trimming and post trimming QC
if(params.cutadapt!=null){
	process adapter_trimming {
            cpus '1'
            memory params.mem_QC+'GB'
            tag { file_tag }
	    
            input:
	    set val(file_tag), file(pairs3) from readPairs3
	    
            output:
            set val(file_tag), file("${file_tag}*val*.fq.gz") into readPairs4
	    set file("${file_tag}${params.suffix1}_val_1_fastqc.zip"), file("${file_tag}${params.suffix2}_val_2_fastqc.zip") into fastqc_postpairs
	    file("${file_tag}*trimming_report.txt") into trimming_reports
	    
	    publishDir "${params.output_folder}/QC/adapter_trimming", mode: 'copy', pattern: '{*report.txt,*fastqc.zip}'
	    
            shell:
            '''
	    trim_galore --paired --fastqc !{pairs3[0]} !{pairs3[1]}
            '''
	}
}else{
	readPairs4 = readPairs3
	fastqc_postpairs=null
	trimming_reports=null
}

readPairs4.into { readPairs_align; readPairs_fusion }


//Mapping, mark duplicates and sorting
process alignment {
      cpus params.cpu
      memory params.mem+'G'
      tag { file_tag }
      
      input:
      set val(file_tag), file(pairs5)  from readPairs_align
      file ref from ref_align.collect()
      file gtf
                  
      output:
      set val(file_tag), file("${file_tag}.bam"), file("${file_tag}.bam.bai") into bam_files
      file("*Log*") into align_out
      set val(file_tag), file("*Chimeric.out.junction") into SJ_out
      if( (params.sjtrim == null)&&(params.recalibration == null) ){
      	publishDir params.output_folder, mode: 'copy', saveAs: {filename ->
                 if (filename.indexOf(".bam") > 0)                      "BAM/$filename"
            else if (filename.indexOf("SJ") > 0)              "BAM/$filename"
	    else if (filename.indexOf("Chimeric") > 0)              "BAM/$filename"
            else if (filename.indexOf("Log") > 0)             "QC/alignment/$filename"
        }
      }else{
	publishDir params.output_folder, mode: 'copy', saveAs: {filename ->
            if (filename.indexOf("SJ") > 0)              "BAM/$filename"
	    else if (filename.indexOf("Chimeric") > 0)             "BAM/$filename"
            else if (filename.indexOf("Log") > 0)             "QC/alignment/$filename"
        }
      }
            
      shell:
      align_threads = params.cpu.intdiv(2)
      sort_threads = params.cpu.intdiv(2) - 1
      sort_mem     = params.mem.intdiv(4)

      if(params.hisat2){
            '''
            hisat2 --rg-id !{file_tag} --rg SM:!{file_tag} --rg !{params.RG} --met-file hisat2.!{file_tag}.Log.final.out -p !{align_threads} -x !{params.hisat2_idx} -1 !{pairs5[0]} -2 !{pairs5[1]} | samblaster --addMateTags | sambamba view -S -f bam -l 0 /dev/stdin | sambamba sort -t !{sort_threads} -m !{sort_mem}G --tmpdir=!{file_tag}_tmp -o !{file_tag}.bam /dev/stdin
	    '''
      }else{
      '''
      STAR --outSAMattrRGline ID:!{file_tag} SM:!{file_tag} !{params.RG} --chimSegmentMin 12 --chimJunctionOverhangMin 12 --chimSegmentReadGapMax 3 --alignSJDBoverhangMin 10 --alignMatesGapMax 200000 --alignIntronMax 200000 --alignSJstitchMismatchNmax 5 -1 5 5 --twopassMode Basic --runThreadN !{align_threads} --genomeDir . --sjdbGTFfile !{gtf} --readFilesCommand zcat --readFilesIn !{pairs5[0]} !{pairs5[1]} --outStd SAM | samblaster --addMateTags | sambamba view -S -f bam -l 0 /dev/stdin | sambamba sort -t !{sort_threads} -m !{sort_mem}G --tmpdir=!{file_tag}_tmp -o !{file_tag}.bam /dev/stdin
      mv Chimeric.out.junction STAR.!{file_tag}.Chimeric.out.junction 
      mv Chimeric.out.sam STAR.!{file_tag}.Chimeric.out.sam
      mv SJ.out.tab STAR.!{file_tag}.SJ.out.tab
      mv Log.final.out STAR.!{file_tag}.Log.final.out
      mv Log.out STAR.!{file_tag}.Log.out
      mv Log.progress.out    STAR.!{file_tag}.Log.progress.out
      mv Log.std.out STAR.!{file_tag}.Log.std.out
      '''
      }
}

fasta_ref       = file(params.ref)
fasta_ref_fai   = file(params.ref + '.fai')

// Fusion-genes detection
if(params.hisat2==null){
process fusion {
      cpus params.cpu
      memory params.mem_QC+'G'
      tag { file_tag }
      
      input:
      set val(file_tag), file(pairs)  from readPairs_fusion
      set val(file_tag), file(SJ) from SJ_out
      file ref from ref_fusion.collect()
      file gtf
      file gtf_gs
      file fasta_ref
      file blast_ava
            
      output:
      file("star_fusion_${file_tag}")
      publishDir "${params.output_folder}/fusion", mode: 'copy'
            
      shell:
      '''
      if [ ! -f ref_annot.gtf ]; then
	 ln -s !{gtf} ref_annot.gtf
      fi
      if [ ! -f ref_genome.fa ]; then
	 ln -s !{fasta_ref} ref_genome.fa
      	 ln -s !{fasta_ref_fai} ref_genome.fa.fai
      fi
      STAR-Fusion --genome_lib_dir . -J !{SJ} --left_fq !{pairs[0]} --right_fq !{pairs[1]} --output_dir star_fusion_!{file_tag} --FusionInspector validate --denovo_reconstruct --annotate --examine_coding_effect
      '''
   }

}

if( (params.sjtrim!=null)||(params.recalibration!=null) ){
    fasta_ref_dictn = params.ref[0..<params.ref.lastIndexOf('.')]
    fasta_ref_dict  = file(fasta_ref_dictn  + '.dict')
}

//Splice junctions trimming
if(params.sjtrim){
   GATK_jar=file(params.GATK_jar)
   
   process splice_junct_trim {
      cpus params.cpu
      memory params.mem+'G'
      tag { file_tag }
      
      input:
      set val(file_tag), file(bam), file(bai)  from bam_files
      file fasta_ref
      file fasta_ref_fai
      file fasta_ref_dict     
      file GATK_jar
            
      output:
      set val(file_tag_new), file("${file_tag_new}.bam"), file("${file_tag_new}.bam.bai") into bam_files2
      if(params.recalibration == null){
        publishDir "${params.output_folder}/BAM", mode: 'copy'
      }
            
      shell:
      file_tag_new = file_tag+'_split'
      '''
      java -Xmx!{params.mem}g -Djava.io.tmpdir=. -jar !{GATK_jar} -T SplitNCigarReads -R !{fasta_ref} -I !{bam} -o !{file_tag_new}.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS
      mv !{file_tag_new}.bai !{file_tag_new}.bam.bai
      '''
   }
}else{
      bam_files2=bam_files
}


//BQSrecalibration
if(params.recalibration){
   GATK_jar     = file(params.GATK_jar)
   bundle_indel = Channel.fromPath(params.GATK_bundle + '/*indels*.vcf')
   bundle_dbsnp = Channel.fromPath(params.GATK_bundle + '/*dbsnp*.vcf')

   process base_quality_score_recalibration {
    	cpus params.cpu
	memory params.mem+'G'
    	tag { file_tag }
        
    	input:
	set val(file_tag), file(bam), file(bai) from bam_files2
	file fasta_ref
      	file fasta_ref_fai
	file fasta_ref_dict
      	file bed
	file GATK_jar
    	file indel from bundle_indel.collect()
	file dbsnp from bundle_dbsnp.collect()
	
    	output:
	set val(file_tag_new), file("${file_tag_new}.bam"), file("${file_tag_new}.bam.bai") into recal_bam_files
    	file("${file_tag}_recal.table") into recal_table_files
    	file("${file_tag}_post_recal.table") into recal_table_post_files
    	file("${file_tag}_recalibration_plots.pdf") into recal_plots_files
    	publishDir params.output_folder, mode: 'copy', saveAs: {filename ->
                 if (filename.indexOf(".bam") > 0)                      "BAM/$filename"
            else "QC/BQSR/$filename"
        }

    	shell:
	file_tag_new = file_tag+'_recal'
    	'''
    	indelsvcf=(`ls *indels*.vcf`)
    	dbsnpvcfs=(`ls *dbsnp*.vcf`)
    	dbsnpvcf=${dbsnpvcfs[@]:(-1)}
    	knownSitescom=''
    	for ll in $indelsvcf; do knownSitescom=$knownSitescom' -knownSites:VCF '$ll; done
    	knownSitescom=$knownSitescom' -knownSites:VCF '$dbsnpvcf
    	java -Xmx!{params.mem}g -Djava.io.tmpdir=. -jar !{GATK_jar} -T BaseRecalibrator -filterRNC -nct !{params.cpu} -R !{fasta_ref} -I !{file_tag}.bam $knownSitescom -L !{bed} -o !{file_tag}_recal.table
    	java -Xmx!{params.mem}g -Djava.io.tmpdir=. -jar !{GATK_jar} -T BaseRecalibrator -filterRNC -nct !{params.cpu} -R !{fasta_ref} -I !{file_tag}.bam $knownSitescom -BQSR !{file_tag}_recal.table -L !{bed} -o !{file_tag}_post_recal.table		
    	java -Xmx!{params.mem}g -Djava.io.tmpdir=. -jar !{GATK_jar} -T AnalyzeCovariates -R !{fasta_ref} -before !{file_tag}_recal.table -after !{file_tag}_post_recal.table -plots !{file_tag}_recalibration_plots.pdf	
    	java -Xmx!{params.mem}g -Djava.io.tmpdir=. -jar !{GATK_jar} -T PrintReads -filterRNC -nct !{params.cpu} -R !{fasta_ref} -I !{file_tag}.bam -BQSR !{file_tag_new}.table -L !{bed} -o !{file_tag_new}.bam
    	mv !{file_tag_new}.bai !{file_tag_new}.bam.bai
    	'''
   }
}else{      
      recal_bam_files=bam_files2
}

recal_bam_files.into { recal_bam_files4QC; recal_bam_files4quant }

//RSEQC
process RSEQC{
    		cpus '1'
		memory params.mem_QC+'GB'
    		tag { file_tag }
        
		input:
    		set val(file_tag), file(bam), file(bai) from recal_bam_files4QC
		file bed
		
    		output:
		file("${file_tag}_readdist.txt") into rseqc_files
    		publishDir "${params.output_folder}/QC/bam", mode: 'copy'

    		shell:
    		'''
		read_distribution.py -i !{file_tag}".bam" -r !{bed} > !{file_tag}"_readdist.txt"
    		'''
}

//Quantification
process quantification{
    	if( (params.sjtrim)||(params.recalibration) ){
		cpus params.cpu
		memory params.mem+'GB'
	}else{
		cpus '1'
		memory params.mem_QC+'GB'
	}
	
    	tag { file_tag }
        
    	input:
    	set val(file_tag), file(bam), file(bai) from recal_bam_files4quant
	file gtf

    	output:
	file("${file_tag}_count.txt") into htseq_files
    	publishDir "${params.output_folder}/counts", mode: 'copy'

    	shell:
	buffer=''
	if(params.htseq_maxreads) buffer='--max-reads-in-buffer '+params.htseq_maxreads
	
	if( (params.sjtrim)||(params.recalibration) ){
	'''
	htseq-count -h
	mv !{file_tag}.bam !{file_tag}_coordinate_sorted.bam
	sambamba sort -n -t !{task.cpus} -m !{params.mem}G --tmpdir=!{file_tag}_tmp -o !{file_tag}.bam !{file_tag}_coordinate_sorted.bam
	htseq-count -r name -s !{params.stranded} -f bam !{file_tag}.bam !{gtf} !{buffer} --additional-attr=gene_name > !{file_tag}_count.txt 
	'''
	}else{
	 	'''
		htseq-count -r pos -s !{params.stranded} -f bam !{file_tag}.bam !{gtf} !{buffer} --additional-attr=gene_name > !{file_tag}_count.txt 
    		'''
	}
}


process multiqc_pretrim {
    cpus '1'
    memory params.mem_QC+'GB'
    tag { "all"}
        
    input:
    file fastqc1 from fastqc_pairs.collect()
        
    output:
    file("multiqc_pretrim_report.html") into multiqc_pre
    file("multiqc_pretrim_report_data") into multiqc_pre_data
    publishDir "${params.output_folder}/QC", mode: 'copy'

    shell:
    '''
    for f in $(find *fastqc.zip -type l);do cp --remove-destination $(readlink $f) $f;done;
    multiqc . -n multiqc_pretrim_report.html -m fastqc
    '''
}


if(params.clustering) htseq_files.into { htseq_files ; htseq_files4clust }


process multiqc_posttrim {
    cpus '1'
    memory params.mem_QC+'GB'
    tag { "all"}
        
    input:
    file STAR from align_out.collect()
    file htseq from htseq_files.collect()
    file rseqc from rseqc_files.collect()
    file trim from trimming_reports.collect()
    file fastqcpost from fastqc_postpairs.collect()
        
    output:
    file("multiqc_posttrim_report.html") into multiqc_post
    file("multiqc_posttrim_report_data") into multiqc_post_data

    publishDir "${params.output_folder}/QC", mode: 'copy'

    shell:
    '''
    for f in $(find *fastqc.zip -type l);do cp --remove-destination $(readlink $f) $f;done;
    multiqc . -n multiqc_posttrim_report.html -m fastqc -m cutadapt -m star -m rseqc -m htseq
    '''
}

if(params.clustering){
   process clustering {
    	cpus params.cpu
	memory params.mem_QC+'G'
    	tag { "all" }
        
    	input:
	file htseq from htseq_files4clust.collect()
	output:
    	file("unsupervised_analysis") into unsup_res

	publishDir params.output_folder, mode: 'move'

    	shell:
    	'''
	RNAseq_unsupervised.R -o unsupervised_analysis -n !{params.clustering_n} -t !{params.clustering_t} -c !{params.clustering_c} -l !{params.clustering_l}
    	'''
   }
}
