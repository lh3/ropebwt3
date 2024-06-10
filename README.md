## Getting Started
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
./ropebwt build -i human100.fmr -do human100.fmd  # not required by recommended

# Count super-maximal exact matches (no locations)
echo CTCCAGTTGACACAAAATAGtCTACGAAAGTGGCTTTAACAT|./ropebwt3 match -L human100.fmd -l20 -
```

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

TODO

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
  functionality using standard suffix array subsampling but it needs to be
  reworked.

* The "merge" command can be accelerated by 10-30% with a more efficient data
  structure but grlBWT will be faster anyway.

[grlbwt]: https://github.com/ddiazdom/grlBWT
[movi]: https://github.com/mohsenzakeri/Movi
[bigbwt]: https://gitlab.com/manzai/Big-BWT
[fm2]: https://github.com/lh3/fermi2
[rb2]: https://github.com/lh3/ropebwt2
