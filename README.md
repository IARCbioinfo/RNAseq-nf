# RNAseq
RNAseq mapping, quality control, and analysis pipeline

## Prerequisites
The following programs need to be installed and in the PATH environment variable:
- [*fastqc*](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/INSTALL.txt)
- [*cutadapt*](http://cutadapt.readthedocs.io/en/stable/installation.html), which requires Python version > 2.7
- [*trim_galore*](https://github.com/FelixKrueger/TrimGalore)
- [*STAR*](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)
- [*nextflow*](https://www.nextflow.io/docs/latest/getstarted.html)

In addition, STAR requires genome indexes that can be generated from a genome fasta file ref.fa and a splice junction annotation file ref.gtf using the following command:
```bash
STAR --runThreadN n --runMode genomeGenerate --genomeDir ref --genomeFastaFiles ref.fa --sjdbGTFfile ref.gtf --sjdbOverhang 99
```

## Usage
To run the pipeline on a series of paired-end fastq files (with suffixes *_1* and *_2*) in folder *fastq*, and a reference genome with indexes in folder *ref_genome*, one can type:
```bash
nextflow run iarcbioinfo/RNAseq-nf --input_folder fastq --gendir ref_genome --suffix1 _1 --suffix2 _2
```
