#!/bin/bash
# mmseqs2跑mmseqs2的benchmark
MMSEQS="/home/wangziyuan/bio/mmseqs/mmseqs/bin/mmseqs"
RESULT="/home/wangziyuan/fxw/work/TASK03/mmseqs2/result"
QUERY="/home/wangziyuan/fxw/work/TASK03/query2/Q_$1.fa"

TARGET="/home/wangziyuan/bio/benchmark/mmseqs2-benchmark-pub/db/targetdb.fasta"
DB="/nethome/wangziyuan/fxw/db/mmseqs2"
TMP="/nethome/wangziyuan/fxw/db/mmseqs2/tmp"

mkdir -p "${DB}"
mkdir -p "${TMP}"

# "${MMSEQS}" createdb "${TARGET}" "${DB}"/targetDB
# "${MMSEQS}" createindex "${DB}"/targetDB "${TMP}"

# sh /home/wangziyuan/cache.sh
"${MMSEQS}" easy-search "${QUERY}" "${DB}"/targetDB "${RESULT}"/alnRes.m8 "${TMP}" -s 1

# "${MMSEQS}" easy-search "${QUERY}" "${TARGET}" "${RESULT}"/alnRes.m8 "${TMP}"