#!/bin/bash
cd ~/alignment-nf/
git config --global user.email "alcalan@fellows.iarc.fr"
git add dag.png
git add dag.html
git commit -m "Generated DAG [skip ci]"
git push origin $CIRCLE_BRANCH

curl -H "Content-Type: application/json" --data "{\"source_type\": \"Branch\", \"source_name\": \"$CIRCLE_BRANCH\"}" -X POST https://registry.hub.docker.com/u/nalcala/rnaseq-nf/trigger/a07f34d1-b3ce-4675-a1b2-3fee9362eae9/


