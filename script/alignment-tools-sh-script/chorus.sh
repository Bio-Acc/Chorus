#!/bin/bash
CHROUS="/home/wangziyuan/bio/chorus/build"
QUERY="/home/wangziyuan/fxw/work/TASK03/query2/Q_$1.fa"
RESULT="/home/wangziyuan/fxw/work/TASK03/chorus/result"

# TARGET="/home/wangziyuan/bio/benchmark/mmseqs2-benchmark-pub/db/targetdb.fasta"
# DB="/data-nvme/wangziyuan/fxw//db/chorus"

TARGET="/pentapool/program/alphafold/data_update/download/uniprot/2023_03/uniref90.fasta"
DB="/nethome/wangziyuan/fxw/db/chorus"

mkdir -p "${DB}"
mkdir -p "${RESULT}"
# 创建数据库
# sh ../other/cache.sh

# 匹配
# echo "createdb..."
# "${CHROUS}/createDB" "${TARGET}" "${DB}/nr" 32

echo "query..."
"${CHROUS}/query" -q "${QUERY}" -d "${DB}/nr" -o "${RESULT}"/output.m8


