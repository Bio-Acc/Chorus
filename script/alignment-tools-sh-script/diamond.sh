#!/bin/bash
QUERY="/home/wangziyuan/fxw/work/TASK03/query2/Q_$1.fa"
DIAMOND="/home/wangziyuan/bio/diamond/diamond"
RESULT="/home/wangziyuan/fxw/work/TASK03/diamond/result"

# TARGET="/home/wangziyuan/bio/benchmark/mmseqs2-benchmark-pub/db/targetdb.fasta"
# REFERENCE="/data-nvme/wangziyuan/fxw/db/diamond/diamond"

TARGET="/pentapool/program/alphafold/data_update/download/uniprot/2023_03/uniref90.fasta"
REFERENCE="/nethome/wangziyuan/fxw/db/diamond/diamond"


# mkdir -p "${REFERENCE}"

# "${DIAMOND}" makedb --in "${TARGET}" --db "${REFERENCE}"

"${DIAMOND}" blastp --query "${QUERY}"  --db "${REFERENCE}" --out "${RESULT}"/output.m8 --fast
