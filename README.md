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

## Usage

``` bash
#simple FM-index of a FASTA/FASTQ file
./cuba index -o test.fmi ../test/test.fa
#bidirectional FM-index of a FASTA/FASTQ file
./cuba index -b -o test.bifmi ../test/test.fa #this is nearly double the size of a monodirectional FM-index
#exact search of a string in the FM-index
./cuba find -f test.fmi GGGGGGGGGGGG #returns one hit in reference_id:1 (second sequence)
#approximate match of a string in the FM-index (bidirectional FM-indexes allow for faster approximate search). Allow up-to 2 errors (max 1 mismatch and 1 insertion).
./cuba find -b -f test.bifmi -e 2 ACGAT #returns one perfect match in reference_id:0 (first squence) and one approximate match (2 dels) in reference_id:1 (second sequence)
```
