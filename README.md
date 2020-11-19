# cuba

Ideally, a collection of modules, written in c++, to perform indexing/searching (and other stuff that are yet to come) efficiently on sequencing data.
For the time being, this is written on the top of [SeqAn3](https://github.com/seqan/seqan3.git).

Available features:

- index. Build a mono/bi-directional FM-index from FASTA/FASTQ files (optionally gzipped).
- find. Search for a given string in a mono/bi-directional FM-index. Approximate search is implemented.


## Installation

``` bash
apt-get update
apt-get install -y git build-essential g++ cmake libz-dev
git clone --recursive https://github.com/davidebolo1993/cuba
mkdir cuba/build
cd cuba/build
cmake ../src
make
./cuba --help
```
