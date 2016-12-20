# RNAseq
RNAseq analysis pipeline

## Prerequisites
STAR requires genome indexes that can be generated from a genome fasta file ref.fa and a splice junction annotation file ref.gtf using the following command:

```bash
STAR --runThreadN n --runMode genomeGenerate --genomeDir ref --genomeFastaFiles ref.fa --sjdbGTFfile ref.gtf --sjdbOverhang 99
```
