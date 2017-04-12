# RNAseq
RNAseq mapping, quality control, and reads counting pipeline

## Overview of pipeline workflow
![workflow](RNAseqpipeline.png?raw=true "Scheme of alignment/realignment Workflow")

## Prerequisites

### General prerequisites
The following programs need to be installed and in the PATH environment variable:
- [*fastqc*](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/INSTALL.txt)
- [*cutadapt*](http://cutadapt.readthedocs.io/en/stable/installation.html), which requires Python version > 2.7
- [*trim_galore*](https://github.com/FelixKrueger/TrimGalore)
- [*fastQC*](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [*multiQC*](http://multiqc.info/docs/)
- [*STAR*](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)
- [*htseq*](http://www-huber.embl.de/HTSeq/doc/install.html#install); the python script htseq-count must also be in the PATH
- [*nextflow*](https://www.nextflow.io/docs/latest/getstarted.html)

In addition, STAR requires genome indices that can be generated from a genome fasta file ref.fa and a splice junction annotation file ref.gtf using the following command:
```bash
STAR --runThreadN n --runMode genomeGenerate --genomeDir ref --genomeFastaFiles ref.fa --sjdbGTFfile ref.gtf --sjdbOverhang 99
```

### Prerequisites for reads trimming at splice junctions
In order to perform the optional reads trimming at splice junctions, GATK must be installed:
- GATK [*GenomeAnalysisTK.jar*](https://software.broadinstitute.org/gatk/guide/quickstart)
- [GATK bundle](https://software.broadinstitute.org/gatk/download/bundle) VCF files with lists of indels and SNVs (recommended: 1000 genomes indels, Mills gold standard indels VCFs, dbsnp VCF)

In addition, index *.fai* and dictionnary *.dict* must be generated from the fasta reference genome using the following commands:
```bash
samtools faidx ref.fa
java -jar picard.jar CreateSequenceDictionary R= ref.fa O= ref.dict
```

### Prerequisites for base quality score recalibration
- GATK [*GenomeAnalysisTK.jar*](https://software.broadinstitute.org/gatk/guide/quickstart)
- [GATK bundle](https://software.broadinstitute.org/gatk/download/bundle) VCF files with lists of indels and SNVs (recommended: 1000 genomes indels, Mills gold standard indels VCFs, dbsnp VCF)
- bed file with intervals to be considered

## Usage
To run the pipeline on a series of paired-end fastq files (with suffixes *_1* and *_2*) in folder *fastq*, and a reference genome with indexes in folder *ref_genome*, one can type:
```bash
nextflow run iarcbioinfo/RNAseq-nf --input_folder fastq --gendir ref_genome --suffix1 _1 --suffix2 _2
```
### Enable reads trimming at splice junctions
To use the reads trimming at splice junctions step, you must add the ***--sjtrim* option**, specify the path to the folder containing the GenomeAnalysisTK jar file, as well as satisfy the requirements above mentionned. For example:
```bash
nextflow run iarcbioinfo/RNAseq-nf --input_folder fastq --gendir ref_genome --suffix1 _1 --suffix2 _2 --sjtrim --GATK_folder /home/user/GATK 
```

### Enable Base Quality Score Recalibration
To use the base quality score recalibration step, you must add the ***--bqsr* option**, specify the path to the folder containing the GenomeAnalysisTK jar file, the path to the GATK bundle folder for your reference genome, specify the path to the bed file with intervals to be considered, as well as satisfy the requirements above mentionned. For example:
```bash
nextflow run iarcbioinfo/RNAseq-nf --input_folder fastq --gendir ref_genome --suffix1 _1 --suffix2 _2 --bqsr --GATK_folder /home/user/GATK --GATK_bundle /home/user/GATKbundle --intervals intervals.bed
```

## All parameters
| **PARAMETER** | **DEFAULT** | **DESCRIPTION** |
|-----------|--------------:|-------------| 
| *--help* | null | print usage and optional parameters |
*--input_folder* | . | input folder |
*--output_folder* |   . | output folder |
*--gendir* | ref | reference genome folder |
*--cpu*          | 8 | number of CPUs |
*--mem*         | 32 | memory|
*--fastq_ext*    | fq.gz | extension of fastq files|
*--suffix1*      | \_1 | suffix for second element of read files pair|
*--suffix2*      | \_2 | suffix for second element of read files pair|
*--output_folder*   | . | output folder for aligned BAMs|
*--fasta_ref* |    ref.fa | reference genome fasta file for GATK |
*--annot_gtf*   |  Homo_sapiens.GRCh38.79.gtf | annotation GTF file |
*--annot_gff*   |  Homo_sapiens.GRCh38.79.gff | annotation GFF file |
*--GATK_folder* |  GATK | folder with jar file GenomeAnalysisTK.jar |
*--GATK_bundle* |  GATK_bundle | folder with files for BQSR |
*--intervals*   |  intervals.bed | bed file with intervals for BQSR | 
*--RG*          |  PL:ILLUMINA | string to be added to read group information in BAM file |
*--sjtrim*      |  false | enable reads trimming at splice junctions | 
*--bqsr*        |  false | enable base quality score recalibration |
