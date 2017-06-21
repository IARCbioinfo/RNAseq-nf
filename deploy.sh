#!/bin/bash
cd ~/RNAseq-nf/
git config --global user.email "alcalan@fellows.iarc.fr"
git add dag_STAR.png
git add dag_hisat2.png
git add dag_STAR.html
git add dag_hisat2.html
git commit -m "Generated DAG [skip ci]"
git push origin $CIRCLE_BRANCH

curl -H "Content-Type: application/json" --data "{\"source_type\": \"Branch\", \"source_name\": \"$CIRCLE_BRANCH\"}" -X POST https://registry.hub.docker.com/u/iarcbioinfo/rnaseq-nf/trigger/a3a4905f-e04d-4a2c-8b8d-acf592b667e4/


