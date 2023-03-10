# Chorus: Heterogeneous GPU+CPU Multiple Protein Sequences Alignment Search for Large Database

Chorus is an efficient protein-protein local alignment search tool for multiple query sequences and large databases. Our method is up to 300x faster than NCBI-BLASTP when inputting multiple query seqs at a time, and can stay ahead of the state-of-the-art database-index methods when there are less than 5000 input query seqs (2.2-27x faster than fastest MMseqs2, 1.1-9x than fastest DIAMOND), while maintaining a low memory footprint.

## Requirements

Hardware:

A linux (ubuntu) server with a NVIDIA GPU.

Software:

    g++ 9
    CUDA 11.7

CUDA Environment (add to ~/.bashrc):

    export CUDA_HOME=/usr/local/cuda
    export PATH=$CUDA_HOME/bin:$PATH
    export CPATH=$CUDA_HOME/include:$CPATH
    export LD_LIBRARY_PATH=$CUDA_HOME/lib64:$LD_LIBRARY_PATH

\*Don't forget to source ~/.bashrc

## Quick Start

### 1. Building from source code

    mkdir build
    cd build
    cmake ..
    make

### 2. Create database (process only once for each database)

Input: A protein database file (in fasta format).

    ./createDB <protein_db.fasta> <output db> <batch size (GB)>

Output: This step will create database files for each part, including "dbnameX.seq", "dbnameX.name", "dbnameX.sofs" and "dbnameX.nofs", to the same directory with the original db file.

E.g.

    ./createDB ../db/nr.fasta ../db/nr 4

### 3. Query the sequences

Input: The pre-processed database in step 2. And a query sequences file (in fasta format).

    ./query -q <query.fasta> -d <db1> [db2 ... (optional)] -o <output file>

Chorus will automatically detect all parts of the database created in the same directory, iteratively process them, and finally show all results in the output file.

E.g.

    ./query -d ../db/nr -q ../example_query -o res.out

## More Arguments for Seaching

Please use

    ./query --help

to get help message.

| Arg                | Description                                                                                                                                                                   | Default                   |
| ------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------- |
| --outfmt           | Output format: 0: m8 tabular (like blast+ outfmt 6); 1: detailed alignment; 2: tabular and ref seqs (fasta); 3: Opfi format; 4: a3m format; 5: Diamond benchmark format.      | 0                         |
| -l, --filter-level | The filter level to pass a database sequence to smith waterman align. The smaller the value, the more sensitive the result.                                                   | 1                         |
| --min-score        | The minimum score for each alignment to display.                                                                                                                              | 0                         |
| -e, --max-evalue   | The maximum expectation value for each alignment to display                                                                                                                   | 1e1                       |
| --max-output-align | The maximum number of alignments to display for each query. (0 means no limit)                                                                                                | 0                         |
| --num-threads      | Number of CPU threads.                                                                                                                                                        | Maximum threads available |
| -k, --seed-length  | Length of k-mer. (3, 4, or 5)                                                                                                                                                 | 5                         |
| -w, --qit-width    | Width of query index table. (Number of indices for each row) The larger the value, the less likely the overflow occurs.                                                       | 4                         |
| -h, --hash-size    | Voting hash table size ratio. The larger the value, the larger the GPU voting hash table (doubly for each increment). (0: maybe OK; 1: OK for most query seqs; 2: redundancy) | 2                         |
| --band-width       | The bandwidth in banded smith-waterman. The larger the value, the more sensitive the result.                                                                                  | 8                         |
| --must-include     | Requires the ref sequence must match a regular expression n times. (e.g. --must-include "R[A-Z]{4,7}H" 2)                                                                     | No limit                  |
