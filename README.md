## <a name="start"></a>Getting Started
```sh
# Compile
git clone https://github.com/lh3/ropebwt3
cd ropebwt3
make  # use "make omp=0" if your compiler doesn't suport OpenMP

# Toy examples
echo -e 'AGG\nAGC' | ./ropebwt3 build -LR -
echo TGAACTCTACACAACATATTTTGTCACCAAG | ./ropebwt3 build -Lbo idx.fmr -
echo ACTCTACACAAgATATTTTGTC | ./ropebwt3 match -Ll10 idx.fmr -

# Acquire the BWT of a human pangenome
wget -O human100.fmr.gz https://zenodo.org/records/11533211/files/human100.fmr.gz?download=1
gzip -d human100.fmr.gz  # decompress
./ropebwt build -i human100.fmr -do human100.fmd  # not required but recommended

# Count super-maximal exact matches (no locations)
echo CTCCAGTTGACACAAAATAGtCTACGAAAGTGGCTTTAACAT|./ropebwt3 match -L human100.fmd -l20 -

# Retrieve chrM of CHM13
./ropebwt3 get human100.fmd 48 > CHM13-chrM.fa
```

## Table of Contents

- [Getting Started](#start)
- [Introduction](#intro)
- [Usage](#use)
  - [Constructing a BWT](#build)
  - [Binary BWT formats](#format)
  - [Counting exact matches](#match)
- [Performance](#perf)
- [Limitations](#limit)

## <a name="intro"></a>Introduction

Ropebwt3 constructs the BWT of a large DNA sequence set and finds super-maximal
exact matches of a query sequence against the BWT. It is optimized for
repetitive sequences such as a pangenome or sequence reads at high coverage. It
can incrementally add new sequences to an existing BWT and is one of the few
methods that can construct the double-strand BWT of 100 human genomes using
reasonable resources (see [Performance](perf) below).

Ropebwt3 has most of the functionality of [ropebwt2][rb2] but works better for
long sequences such as chromsomes and assembled contigs. It additionally
includes some functionality in [fermi2][fm2].

## <a name="use"></a>Usage

### <a name="build"></a>Constructing a BWT

Ropebwt3 implements three algorithms for BWT construction. For the best
performance, you need to choose an algorithm based on the input date types.

1. If you are not sure, use the general command line
   ```sh
   ropebwt3 build -t24 -bo bwt.fmr file1.fa file2.fa filen.fa
   ```
   You can also append another file to an existing index
   ```sh
   ropebwt3 build -t24 -bo bwt-new.fmr -i bwt-old.fmr filex.fa
   ```

2. For a set of large genomes (e.g. a human pangenome), you may generate the
   BWT of each individual genome on a cluster and merge them togather. This
   parallelizes sub-BWT construction and speeds up the overall process.
   ```sh
   ropebwt3 build -t8 -bo genome1.fmr genome1.fa.gz
   ropebwt3 build -t8 -bo genome2.fmr genome2.fa.gz
   ropebwt3 build -t8 -bo genomen.fmr genomen.fa.gz
   ropebwt3 merge -t24 -bo bwt.fmr genome1.fmr genome2.fmr genomen.fmr
   ```

3. For a set of small genomes, it is better to concatenate them together:
   ```sh
   cat file1.fa file2.fa filen.fa | ropebwt3 build -t24 -bo bwt.fmr -
   ```

4. For short reads, use the ropebwt2 algorithm and enable the RCLO sorting:
   ```sh
   ropebwt3 build -r -bo bwt.fmr reads.fq.gz
   ```

5. Use [grlBWT][grlbwt], which is faster than ropebwt3 for large pangenomes:
   ```sh
   ropebwt3 fa2line genome1.fa genome2.fa genomen.fa > all.txt
   grlbwt-cli all.txt -t 32 -T . -o bwt.grl
   grl2plain bwt.rl_bwt bwt.txt
   ropebwt3 plain2fmd -o bwt.fmd bwt.txt
   ```

These command lines construct a BWT for both strands of the input sequences.
You can skip the reverse strand by adding option `-R`.

When you use the `build` command for one genome, the peak memory by default is
$`B+17\cdot\min\{2S,7{\rm g}\}`$ where $S$ is the input file size and $B$ is
the final BWT size in run-length encoding. You can reduce the peak memory by
reducing the batch size via option `-m`.

The peak memory for the `merge` command is
$`B+\max\{B_1,\ldots,B_n\}+8\max\{L_1,\ldots,L_n\}`$, where $B$ is the final
BWT size in run-length encoding, $`B_i`$ is the size of the $i$-th input BWT
and $`L_i`$ is the number of symbols in the $i$-th BWT to be merged.

If you provide multiple files on a `build` command line, ropebwt3 internally
will run `build` on each input file and then incrementally merge each
individual BWT to the final BWT. The peak memory will be the higher one between
`build` and `merge`.

### <a name="format"></a>Binary BWT file formats

Ropebwt3 uses two binary formats to store run-length encoded BWTs: the ropebwt2
FMR format and the fermi FMD format. The FMR format is dynamic in that you can
add new sequences or merge BWTs to an existing FMR file. The same BWT does not
necessarily lead to the same FMR. The FMD format is simpler in structure,
faster to load, smaller in memory and can be memory-mapped.

At present, the FMR-to-FMD conversion works but the FMD-to-FMR conversion is
buggy. It is recommended to use FMR for BWT construction. You may convert FMR
to FMD at the end for matching.

### <a name="match"></a>Counting exact matches

A maximal exact match (MEM) is an exact alignment between the index and a query
that cannot be extended in either direction. A super MEM (SMEM) is a MEM that
is not contained in any other MEM. You can find SMEM of a query if your BWT is
constructed from both strands of sequences.
```sh
ropebwt3 match bwt.fmd query.fa > matches.bed
```
In the output, the first three columns give the query sequence name, start and
end of a match and the fourth column gives the number of hits. As of now,
**ropebwt3 does not report the locations of matches**.

If the BWT only contains one strand, you can use the `suffix` command to find
the longest matching suffix of a query sequence:
```sh
ropebwt3 suffix bwt.fmd query.fa > suffixes.bed
```

## <a name="perf"></a>Performance

The following table shows the time to construct the BWT for 100 human haplotype
assemblies on both strands (~600 Gbp sequences in input). The following methods
were evaluated:

* `rb3 build`: construct BWT from input FASTA files with `ropebwt3 build`. This
  is the only method here that does not use working disk space.

* `rb3 merge`: merge 100 BWTs constructed from 100 FASTA files, respectively.
  Constructing the BWT for one human genome takes 10-12 minutes, which is not
  counted in the table.

* `grlBWT`: construct BWT using [grlBWT][grlbwt]. We need to concatenate all
  input FASTA files and convert them to the one-sequence-per-line format with
  `ropebwt3 fa2line`. Conversion time is not counted.

* `pfp-thresholds`: launched via the [Movi][movi] indexing script. The
  `newscanNT.x` step (from [Big-BWT][bigbwt]) finished in 15.9 hours but the
  `pfp-thresholds` command failed apparently due to insufficient memory.

|                 |rb3 bulid|rb3 merge|grlBWT|pfp-thresholds|
|:----------------|--------:|--------:|-----:|:-------------|
|Elaspsed time (h)|     49.2|     24.2|   8.3|>15.9 (unfinished)|
|CPU time (h)     |    792.6|    757.2|  29.6|>15.7 (unfinished)|
|Peak memory (GB) |    114.2|     70.7|  84.8|>300 (out-of-memory)|

grlBWT is clearly the winner for BWT construction and it also works for non-DNA
alphabet. Ropebwt3 has acceptable performance and its support of incremental
build may be helpful for large datasets.

## <a name="limit"></a>Limitations

* The "match" command of ropebwt3 only counts the number of hits but does not
  report the locations of the hits. [Fermi2][fm2] already supports such
  functionality using standard sampled suffix array but it needs to be
  reworked.

* The "merge" command can be accelerated by 10-30% with a more efficient data
  structure but grlBWT will be faster anyway.

[grlbwt]: https://github.com/ddiazdom/grlBWT
[movi]: https://github.com/mohsenzakeri/Movi
[bigbwt]: https://gitlab.com/manzai/Big-BWT
[fm2]: https://github.com/lh3/fermi2
[rb2]: https://github.com/lh3/ropebwt2
