# Algorithm Engineering project — Genome Mapping 

## Overview

I have implemented a mapper that aligns all reads to a genome using the FM-index method. The mapper supports only exact matching with wildcards 'N' handling.

Key characteristics of the approach:

- 'N' characters in reads are treated as wildcards (match any base in the genome)
- Each read is searched for exact matches in both the original genome and its reverse complement
- Results are written to `results.txt` in the following format:

```
<id> | <cnt_occur_gen> <index1,index2,...> | <cnt_occur_comp> <index1,index2,...>
```

where each index is a 0-based position in the original genome sequence.

## Algorithm

The algorithm proceeds as follows:

1. Build a suffix array for the genome in O(|s| log |s|) time
2. Construct the BWT from the suffix array in O(|s|)
3. Precompute the letter-count array and occurrence table for FM-index queries

For each read t:

1. Find the longest substring t[l:r] that contains no 'N' characters
2. Use FM-index backward search to find all occurrences of this substring in O(|t|)
3. For each candidate position, verify that the full read matches in O(|t|)
4. Repeat the same process for the reverse complement of the read

The worst-case complexity is O(|t| * num_candidates) per read, but in practice the maximum number of candidates observed was only 46.

Since reads are processed independently, the search is parallelized using OpenMP.

## Program report

```txt
Genome size: 4641652 chars
Total reads: 22720100
Unmapped: 3100549 (13.6467%)
Mapped: 19619551 (86.3533%)
  Unique (1 occurrence): 19257169 (84.7583%)
  Double (2 occurrences): 68069 (0.299598%)
  Multi (3+ occurrences): 294313 (1.29539%)
Alignment quality: 100%, exact matching only, 'N' is a wildcard
Genome coverage: 98.3908%
Total time: 42.8 s
```

## Requirements:

- `ERR022075_1.fastq` and `GCF_000005845.2_ASM584v2_genomic.fna` must be in the same directory as `main.cpp`
- OpenMP support required