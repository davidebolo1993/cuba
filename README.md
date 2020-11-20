# cuba

Ideally, a collection of modules, written in c++, to perform indexing/searching (and other stuff that are yet to come) efficiently on sequencing data.
For the time being, this is written on the top of [SeqAn3](https://github.com/seqan/seqan3.git).

Available features:

- index. Build a mono/bi-directional FM-index from FASTA/FASTQ files (optionally gzipped).
- find. Search for a given string in a mono/bi-directional FM-index. Approximate search is implemented.

This is a work-in-progress.

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

### index

``` bash
#simple FM-index of a FASTA/FASTQ file. Multiple FASTA/FASTQ can be specified as positional arguments
./cuba index -f test.fmi ../test/test.fa
#bidirectional FM-index of a FASTA/FASTQ file. This is nearly double the size of a monodirectional FM-index
./cuba index -b -f test.bifmi ../test/test.fa
#additionally store a vector of FASTA/FASTQ sequences to file. Can be used for printing alignments
./cuba index -b -f test.bifmi -v seqvec.obj ../test/test.fa
```

### find

``` bash
#exact search of a string in the FM-index
./cuba find -f test.fmi GGGGGGGGGGGG #returns one hit in the second sequence (starting at base 12)
#approximate match of a string in the FM-index (bidirectional FM-indexes allow for faster approximate search). Allow 1 error
./cuba find -b -f test.bifmi -e 1 ATTTAT #return multiple hits in the first sequence (and one in the second)
#print alignment of query string against original FASTA/FASTQ sequence
./cuba find -b -f test.bifmi -e 1 -v seqvec.obj ATTTAT #print alignment. Gaps are insertions in the query (or deletions in the original FASTA/FASTQ sequence).
```
