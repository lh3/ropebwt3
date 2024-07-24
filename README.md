## <a name="start"></a>Getting Started
```sh
# Compile
git clone https://github.com/lh3/ropebwt3
cd ropebwt3
make  # use "make omp=0" if your compiler doesn't suport OpenMP

# Toy examples
echo -e 'AGG\nAGC' | ./ropebwt3 build -LR -
echo TGAACTCTACACAACATATTTTGTCACCAAG | ./ropebwt3 build -Lbo idx.fmr -
echo ACTCTACACAAgATATTTTGTC | ./ropebwt3 search -Ll10 idx.fmr -

# Download the prebuilt BWT of 152 M. tuberculosis genomes
wget -O- https://zenodo.org/records/12803206/files/mtb152.tar.gz?download=1 | tar -zxf -

# Count super-maximal exact matches (no contig positions)
echo ACCTACAACACCGGTGGCTACAACGTGG  | ./ropebwt3 mem -L mtb152.fmd -

# Local alignment
echo ACCTACAACACCGGTaGGCTACAACGTGG | ./ropebwt3 sw -Lm20 mtb152.fmd -

# Retrieve R15311, the 46th genome in the collection, where 90=(46-1)*2
./ropebwt3 get mtb152.fmd 90 > R15311.fa
```

## Table of Contents

- [Getting Started](#start)
- [Introduction](#intro)
- [Usage](#use)
  - [Constructing a BWT](#build)
  - [Binary BWT formats](#format)
  - [Counting maximal exact matches](#mem)
  - [Local alignment](#bwasw)
- [For developers](#dev)
- [Algorithms](#algo)
  - [BWT construction](#algo-build)
  - [Searching](#algo-match)
- [Performance](#perf)
- [Limitations](#limit)

## <a name="intro"></a>Introduction

Ropebwt3 constructs the BWT of a large DNA sequence set and searches for
matches of a query sequence against the BWT. It is optimized for
repetitive sequences such as a pangenome or sequence reads at high coverage.
On BWT construction, ropebwt3 is one of the few methods that can construct the
double-strand BWT of 100 human genomes using reasonable resources (see
[Performance](#perf) below). It also indexes 7.3Tb of common bacterial
genomes into a 30GB run-length encoded BWT file. These BWTs can be downloaded
[from Zenodo][zenodo]. On alignment, ropebwt3 finds supermaximal exact matches
(SMEMs) efficiently or inexact hits with mismatches and gaps at
slower speed. It does not report the locations of these hits for now.
<!--
Ropebwt3 largely replaces [ropebwt2][rb2] and works better for long sequences
such as chromsomes and assembled contigs. It additionally includes
functionality in [fermi2][fm2] such as counting and searching.

Ropebwt3 is mostly a research project I use to understand the performance of
BWT construction. It may also be useful if you want to get the occurrence of a
string of arbitrary length, or count long matches against a large pangenome.
-->
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
   cat file1.fa file2.fa filen.fa | ropebwt3 build -t24 -m2g -bo bwt.fmr -
   ```

4. For short reads, use the ropebwt2 algorithm and enable the RCLO sorting:
   ```sh
   ropebwt3 build -r -bo bwt.fmr reads.fq.gz
   ```

5. Use [grlBWT][grlbwt], which is [faster](#perf) than ropebwt3 for large pangenomes:
   ```sh
   ropebwt3 fa2line genome1.fa genome2.fa genomen.fa > all.txt
   grlbwt-cli all.txt -t 32 -T . -o bwt.grl
   grl2plain bwt.rl_bwt bwt.txt
   ropebwt3 plain2fmd -o bwt.fmd bwt.txt
   ```

These command lines construct a BWT for both strands of the input sequences.
You can skip the reverse strand by adding option `-R`.

When you use the `build` command for one genome, the peak memory by default is
$`B+11\cdot\min\{2S,7{\rm g}\}`$ where $S$ is the input file size and $B$ is
the final BWT size in run-length encoding. If you have more than 65536
sequences in a batch, factor 11 will be increased to 17 due to the use of a
different algorithm. You can reduce the peak memory by reducing the batch size
via option `-m`.

The peak memory for the `merge` command is
$`B+\max\{B_1,\ldots,B_n\}+8\max\{L_1,\ldots,L_n\}`$, where $B$ is the final
BWT size in run-length encoding, $`B_i`$ is the size of the $i$-th input BWT
to be merged and $`L_i`$ is the number of symbols in the $i$-th BWT.

If you provide multiple files on a `build` command line, ropebwt3 internally
will run `build` on each input file and then incrementally merge each
individual BWT to the final BWT. The peak memory will be the higher one between
the `build` step and the `merge` step.

### <a name="format"></a>Binary BWT file formats

Ropebwt3 uses two binary formats to store run-length encoded BWTs: the ropebwt2
FMR format and the fermi FMD format. The FMR format is dynamic in that you can
add new sequences or merge BWTs to an existing FMR file. The same BWT does not
necessarily lead to the same FMR. The FMD format is simpler in structure,
faster to load, smaller in memory and can be memory-mapped. The two formats can
be used interchangeably in ropebwt3, but it is recommended to use FMR for BWT
construction and FMD for finding exact matches. You can explicitly convert
between the two formats with:
```sh
ropebwt3 build -i in.fmd -bo out.fmr
ropebwt3 build -i in.fmr -do out.fmd
```

### <a name="mem"></a>Counting maximal exact matches

A maximal exact match (MEM) is an exact alignment between the index and a query
that cannot be extended in either direction. A super MEM (SMEM) is a MEM that
is not contained in any other MEM on the query sequence. You can find the SMEMs
of a query **provided that your BWT is constructed from both strands of sequences**.
```sh
ropebwt3 mem -t4 bwt.fmd query.fa > matches.bed
```
In the output, the first three columns give the query sequence name, start and
end of a match and the fourth column gives the number of hits.
<!--
If searching for SMEMs is slow, you may add option `-g` to look for greedy MEMs
which are found by a forward search followed by a backward search from the
furthest forward search. Similar to SMEMs, greedy MEMs are MEMs and are not
contained in each other on the query sequence. However, greedy MEMs are not
always SMEMs. They are approximate.

If the BWT only contains one strand, you can use the `suffix` command to find
the longest matching suffix of a query sequence:
```sh
ropebwt3 suffix bwt.fmd query.fa > suffixes.bed
```
-->
### <a name="bwasw"></a>Local alignment
```sh
ropebwt3 sw -t4 bwt.fmd query.fa > aln.paf
```

## <a name="dev"></a>For Developers

You can encode and decode a FMD file with [rld0.h](rld0.h) and
[rld0.c](rld0.c). The two-file library also supports the rank() operator. Here
is a small program to convert FMD to plain text:
```c
// compile with "gcc -O3 rld0.c this.c"; run with "./a.out idx.fmd > out.txt"
#include <stdio.h>
#include "rld0.h"
int main(int argc, char *argv[]) {
  if (argc < 2) return 1;
  rld_t *e = rld_restore(argv[1]);
  rlditr_t ei; // iterator
  rld_itr_init(e, &ei, 0);
  int c;
  int64_t i, l;
  while ((l = rld_dec(e, &ei, &c, 0)) > 0)
    for (i = 0; i < l; ++i) putchar("\nACGTN"[c]);
  rld_destroy(e);
  return 0;
}
```
and to count a string in an FMD file:
```c
// compile with "gcc -O3 rld0.c this.c"; run with "./a.out idx.fmd AGCATAG"
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include "rld0.h"
int main(int argc, char *argv[]) {
  if (argc < 3) return 1;
  rld_t *e = rld_restore(argv[1]);
  uint64_t k = 0, l = e->cnt[6], ok[6], ol[6];
  const char *s = argv[2];
  int i, len = strlen(s);
  for (i = len - 1; i >= 0; --i) { // backward search
    int c = s[i];
    c = c=='A'?1:c=='C'?2:c=='G'?3:c=='T'?4:5;
    rld_rank2a(e, k, l, ok, ol);
    k = e->cnt[c] + ok[c];
    l = e->cnt[c] + ol[c];
    if (k == l) break;
  }
  printf("%ld\n", (long)(l - k));
  rld_destroy(e);
  return 0;
}
```

## <a name="algo"></a>Algorithms

### <a name="algo-build"></a>BWT construction

Ropebwt3 effectively appends a distinct sentinel to each string in the string
set such that we never need to compare suffixes beyond sentinels. This is an
essential assumption behind ropebwt3. Ropebwt3 would become much slower if you
concatenate all strings without sentinels.

Like ropebwt2, ropebwt3 uses a B+-tree to store a run-length encoded BWT.
It can either insert sequences or merge BWTs into the B+-tree. The sequence
insertion algorithm is identical to ropebwt2. Please see [its paper][rb2-paper]
for details.

BWT merging can be done in several equivalent ways and has been
described in multiple papers. More exactly in ropebwt3, suppose we want to
append the $`i_0`$ sequence in BWT $`B'`$ into $`B`$. We start with $`k\gets
C({\rm A})`$ and $`i\gets i_0`$, calculate the position of $`B'[i]`$ in the
final merged BWT with
```math
{\rm rank}(B'[i],k)+{\rm rank}'(B'[i],i)
```
and update $k$ and $i$ by
```math
k\gets C(B'[i])+{\rm rank}(B'[i],k),\, i\gets C'(B'[i])+{\rm rank}'(B'[i],i)
```
until $`B'[i]`$ is a seninel. When we have the final position of each symbol
$`B'[i]`$, we insert them into $B$ to generate the merged BWT.

### <a name="algo-match"></a>Searching

A classical BWT only supports backward search but if a BWT contains both
strands of DNA sequences, it will support both forward and backward searches.
And with search in both directions, it is possible to find SMEMs. Please see
the [fermi paper][fm-paper] for details.

## <a name="perf"></a>Performance

The following table shows the time to construct the BWT for three datasets:

1. human100 (300Gb): 100 human genomes assembled with long reads from the pangene paper
2. ecoli315k (1.6Tb): 315k *E. coli* genomes from [AllTheBacteria v0.2][atb02]
3. CommonBacteria (7.3Tb): genomes from AllTheBacteria excluding those in the "dustbin" and "unknown" categories

BWTs are constructed from both strands, so the size of each BWT doubles the
number of input bases.

|Dataset        | Metric          |rb3 bulid|rb3 merge|grlBWT|pfp-thresholds|
|:--------------|:----------------|--------:|--------:|-----:|-------------:|
|human100       |Elapsed time (h) |     33.7|     24.2|   8.3| 51.7 |
|               |CPU time (h)     |    803.6|    757.2|  29.6| 51.5 |
|               |Peak memory (GB) |     82.3|     70.7|  84.8| 788.1 |
|ecoli315k      |Elapsed time (h) |    128.7|
|               |CPU time (h)     |   3826.8|
|               |Peak memory (GB) |     20.5|
|CommonBacteria |Elapsed time (d) |     26.5|
|               |CPU time (d)     |    830.3|
|               |Peak memory (GB) |     67.3|

For human100, the following methods were evaluated:

* `rb3 build`: construct BWT from input FASTA files with `ropebwt3 build -t48`
  (using up to 48 threads). This is the only method here that does not use
  working disk space.

* `rb3 merge`: merge 100 BWTs constructed from 100 FASTA files, respectively.
  Constructing the BWT for one human genome takes around 10 minutes, which is not
  counted in the table.

* `grlBWT`: construct BWT using [grlBWT][grlbwt]. We need to concatenate all
  input FASTA files and convert them to the one-sequence-per-line format with
  `ropebwt3 fa2line`. Conversion time is not counted.

* `pfp-thresholds`: launched via the [Movi][movi] indexing script. It was run
  on a slower machine with more RAM. The time for `prepare_ref` is not counted,
  either.

**grlBWT is clearly the winner for BWT construction** and it also works for non-DNA
alphabet. Ropebwt3 has acceptable performance and its support of incremental
build may be helpful for large datasets.

## <a name="limit"></a>Limitations

* Ropebwt3 only supports sampled suffix array for retrieving positions. An
  r-index will probably do better.

* The "merge" command can be accelerated by 10-30% with a more efficient data
  structure but grlBWT will be faster anyway.

[grlbwt]: https://github.com/ddiazdom/grlBWT
[movi]: https://github.com/mohsenzakeri/Movi
[bigbwt]: https://gitlab.com/manzai/Big-BWT
[fm2]: https://github.com/lh3/fermi2
[rb2]: https://github.com/lh3/ropebwt2
[zenodo]: https://zenodo.org/records/11533210
[rb2-paper]: https://academic.oup.com/bioinformatics/article/30/22/3274/2391324
[fm-paper]: https://academic.oup.com/bioinformatics/article/28/14/1838/218887
[atb02]: https://ftp.ebi.ac.uk/pub/databases/AllTheBacteria/Releases/0.2/
