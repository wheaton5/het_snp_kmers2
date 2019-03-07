# het_snp_kmers
memory efficient de novo detection of het snp kmers using counting bloom filters.

Install requirements: rust ver 1.3 or later, clang
```
curl https://sh.rustup.rs -sSf | sh
echo 'export PATH=~/.cargo/bin:$PATH' >> ~/.bashrc
source ~/.bashrc
which cargo
```
If the build fails on the htslib dependency you might need xz. You will need xz for the htslib dependency. 
```
export CFLAGS='-I/path/to/xz/<version>/include'
or add that to your .bashrc and source it
```
Then you should be able to clone and install the project.
```
git clone git@github.com:wheaton5/het_snp_kmers2.git
cd het_snp_kmers2
cargo build --release
```
To get all het kmers and run distributed on lsf run the python script named het_kmer_distributor. If you have a different cluster setup it should be relatively easy to edit this script accordingly. 
requires python 3 but just change print statements if you prefer python 2
```
./het_kmer_distributer -h
usage: het_kmer_distributer [-h] -n NUM_JOBS -i INPUTS [INPUTS ...] -o OUTPUT
                            [-k KMER_SIZE] -m MIN_COUNT -H OUTPUT_HIST -u
                            ESTIMATED_UNIQUE_KMERS [-j JOB_OUTPUT] -M
                            MEM_PER_JOB

distribute het kmer detection

optional arguments:
  -h, --help            show this help message and exit
  -n NUM_JOBS, --num_jobs NUM_JOBS
                        number of jobs
  -i INPUTS [INPUTS ...], --inputs INPUTS [INPUTS ...]
                        input files, takes fasta/fastq (optionally gzipped),
                        sam, bam
  -o OUTPUT, --output OUTPUT
                        output filename for het kmers
  -k KMER_SIZE, --kmer_size KMER_SIZE
                        kmer size to use. default = 21, supports up to 31.
                        must be odd
  -m MIN_COUNT, --min_count MIN_COUNT
                        mininum kmer count to output.
  -H OUTPUT_HIST, --output_hist OUTPUT_HIST
                        output hist file name
  -u ESTIMATED_UNIQUE_KMERS, --estimated_unique_kmers ESTIMATED_UNIQUE_KMERS
                        estimated number of unique kmers in your dataset, rule
                        of thumb is 2-3x genome size
  -j JOB_OUTPUT, --job_output JOB_OUTPUT
                        file to output bsub messaging
  -M MEM_PER_JOB, --mem_per_job MEM_PER_JOB
                        mem per job in megabytes

```

To get a subset of the het snps you might want to run the rust executable directly with an appropriate modimizer. It will only count kmers whose 64bit representation in 2bit encoding is kmer % modimizer == mod_remainder
```
./target/release/het_snp_kmers -h
het_snp_kmers 1.0
Haynes Heaton <whheaton@gmail.com>
Finds kmer pairs that are different in the middle base and each have roughly haploid coverage. Meant for illumina data
as an initial step for de novo phasing.

USAGE:
    het_snp_kmers [OPTIONS] --estimated_kmers <estimated_kmers> --inputs <inputs>... --min_coverage <min_coverage> --output_full_hist <output_full_hist>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
        --estimated_kmers <estimated_kmers>
            estimated total unique kmers. good rule of thumb is roughly 2 * genome size

    -i, --inputs <inputs>...
            input sequence files (fastq,fasta can be gzipped, sam, bam) from which to find het snp kmers

    -k, --kmer_size <kmer_size>                  kmer size to use, defaults to 21
        --min_coverage <min_coverage>            min coverage for each kmer of the pair
        --mod_remainder <mod_remainder>          this is a test
        --modimizer <modimizer>                  this is a test
    -o, --output <output>
        --output_full_hist <output_full_hist>    file name for full kmer histogram
```


Example on made up small test data
```
./target/release/het_snp_kmers --inputs test/data/test.fastq.gz --estimated_kmers 200 --min_coverage 4 --output_full_hist hist.tsv
AAAAAAGGGGACCCCCTTTTT   7       AAAAAAGGGGGCCCCCTTTTT   4       0
```


If using this in conjunction with haplovalidate (https://github.com/wheaton5/haplovalidate) I suggest you run the binary with a reasonable modimizer such as 11 or so and with a min coverage that removes most of the error kmers (10, maybe 20 if high coverage). Then you can look can look at the full histogram and decide on thresholds for filtering.
Then there is a filter.py script to filter your het kmers.
```
./filter.py -h
usage: filter.py [-h] [--min MIN] [--max MAX] [--max_offsite MAX_OFFSITE]
                 [--max_sum MAX_SUM] [--min_sum MIN_SUM]
                 input

filter het kmer set on coverage thresholds

positional arguments:
  input                 input kmer set file as output by het_snp_kmers

optional arguments:
  -h, --help            show this help message and exit
  --min MIN             min coverage
  --max MAX             max coverage
  --max_offsite MAX_OFFSITE
                        max value of 3rd+4th middle base kmers
  --max_sum MAX_SUM     max sum of coverage
  --min_sum MIN_SUM     min sum of coverage
```
which will dump the filtered het kmers to stdout. This filtered het kmer set can then be fed into haplovalidate.
