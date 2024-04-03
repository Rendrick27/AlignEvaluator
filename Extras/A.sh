#!/bin/bash

# Run Snakemake
snakemake --use-conda all --cores all

#iqtree -s Main.fasta --trees all.trees --test-weight --test-au -n 0 --test 10000 -pre IQtree

