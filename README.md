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
```

## Introduction

## Performance

The following table shows the time to construct the BWT for 100 human haplotype
assemblies on both strands (~600 Gbp sequences in input). The following methods
were evaluated:

* `rb3 build`: construct BWT from input FASTA files without using disk space

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

[grlbwt]: https://github.com/ddiazdom/grlBWT
[movi]: https://github.com/mohsenzakeri/Movi
[bigbwt]: https://gitlab.com/manzai/Big-BWT
