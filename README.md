## Getting Started
```sh
# compile
git clone https://github.com/lh3/ropebwt3
cd ropebwt3
make  # add "omp=0" if your compiler doesn't suport OpenMP

# construct BWT
echo -e 'AGG\nAGC'|./ropebwt3 build -LR -
```
