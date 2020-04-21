#!/bin/bash
cd ~/project/
git config --global user.email "alcalan@fellows.iarc.fr"
git config --global user.name "Circle CI_$CIRCLE_PROJECT_REPONAME_$CIRCLE_BRANCH"
git add .
git status
git commit -m "Generated DAG [skip ci]"
git push origin $CIRCLE_BRANCH

curl -H "Content-Type: application/json" --data "{\"source_type\": \"Branch\", \"source_name\": \"$CIRCLE_BRANCH\"}" -X POST https://registry.hub.docker.com/u/iarcbioinfo/rnaseq-nf/trigger/a3a4905f-e04d-4a2c-8b8d-acf592b667e4/


