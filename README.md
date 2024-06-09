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

|                 |ropebwt3 bulid|ropebwt3 merge|grlBWT|pfp-thresholds|
|:----------------|-------------:|-------------:|-----:|:-------------|
|Elaspsed time (h)|          49.2|          24.2|   8.3|>15.9 (unfinished)|
|CPU time (h)     |         792.6|         757.2|  29.6|>15.7 (unfinished)|
|Peak memory (GB) |         114.2|          70.7|  84.8|>300 (out-of-memory)|
