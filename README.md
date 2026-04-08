# Ancient DNA C++ Tools

A collection of high-performance C++ tools for analyzing and processing ancient DNA sequencing data in BAM and VCF formats.

## Tools

This toolkit includes seven tools:

- **analyzeBAM**: Analyze BAM files to generate comprehensive statistics on alignments with filtering and optional duplicate removal
- **dedupBAM**: Multi-threaded PCR duplicate removal tool optimized for paired-end sequencing data
- **filterBAM**: Filter BAM files based on deamination patterns (C→T substitutions)
- **deamBAM**: Analyze damage patterns in ancient DNA, particularly terminal substitution frequencies
- **analyzeBAMCoverage**: Analyze coverage in BAM files with optional GC content stratification and read-length binning
- **analyzeVCF**: Analyze VCF files to generate statistics on non-reference alleles
- **splitBAM**: Split BAM files into multiple chunks for parallel processing

## Installation

### Prerequisites

- C++17 compiler (GCC recommended)
- HTSlib (1.4 or newer)
- Required libraries: libz, libbz2, liblzma, libpthread

### Project Structure

```
ancient_dna_cpp_tools/
├── src/           # Source files (.cpp)
├── include/       # Header files (.h)
├── obj/           # Object files (generated during build)
├── Makefile       # Build configuration
└── README.md      # Documentation
```

### Building from source

```bash
git clone [repository-url]
cd ancient_dna_cpp_tools
make
```

The build system will automatically:
- Create the `obj/` directory for object files
- Compile source files from `src/` using headers from `include/`
- Generate executables in the root directory

To clean the build:
```bash
make clean
```

## Usage

### analyzeBAM

Analyze one or more BAM files and generate statistics with optional filtering and duplicate removal.

```bash
./analyzeBAM [options] <BAM files>
```

Options:
- `-out_folder <path>`: Output folder (default: current directory)
- `-min_len <length>`: Minimum alignment length (default: 35)
- `-min_map_qual <quality>`: Minimum mapping quality (default: 0)
- `-targetfile <path>`: BED file with target regions (tab-separated, 0-based)
- `-paired`: Do not disregard paired reads
- `-count_f`: Ignore filter (QC failed) flag
- `-remove_dups`: Remove PCR duplicates during analysis
- `-ignore_read_groups`: Ignore read group information when identifying duplicates
- `-help`: Display help message

### dedupBAM

A multi-threaded PCR duplicate removal tool that can handle paired-end sequencing data.

```bash
./dedupBAM [options] <BAM file>
```

Options:
- `-out_folder <path>`: Output folder (default: current directory)
- `-min_len <length>`: Minimum alignment length (default: 0)
- `-min_map_qual <quality>`: Minimum mapping quality (default: 0)
- `-ignore_read_groups`: Ignore read group information when identifying duplicates
- `-ignore_qc_fail`: Process reads that failed QC
- `-threads <num>`: Number of threads to use (default: 4)
- `-keep_chr_bams`: Keep chromosome-level BAM files
- `-help`: Display help message

#### Output Format

The output BAM includes:
- `XP:i` tag: Number of PCR duplicates for this read (1 = unique)
- `XY:Z` tag: Coordinates of duplicate reads (for non-unique reads)

### filterBAM

Filter BAM files based on deamination patterns at specific positions.

```bash
./filterBAM [options] <input.bam>
```

Options:
- `-out_folder <output_path>`: Output folder (default: current directory)
- `-suffix <string>`: String to add to output filename
- `-p3`: 3' C→T differences at specified positions (e.g., 0,-1,-2)
- `-p5`: 5' C→T differences at specified positions (e.g., 0,1,2)
- `-h`: Display help message

### deamBAM

Analyze substitution patterns in ancient DNA, particularly focusing on damage.

```bash
./deamBAM [options] <BAM file>
```

Options:
- `-threads <num>`: Number of threads (default: 40)
- `-out_folder <path>`: Output folder (default: current directory)
- `-min_len <num>`: Minimum read length (default: 35)
- `-max_len <num>`: Maximum read length (default: 300)
- `-help`: Show help message

### analyzeBAMCoverage

Generate coverage depth histograms across mappable regions, optionally stratified by read-length bins and GC bins.

```bash
./analyzeBAMCoverage [options]
```

Required arguments:
- `-b <path>`: Input BAM file (coordinate-sorted, indexed)
- `-o <path>`: Output directory
- `-ref <str>`: Reference genome name: `t2t | hg38 | hg19`
- `-r <csv>`: Read-length bin edges CSV (edges are inclusive/exclusive as `[lo, hi)`), e.g. `35,60,100,142`
- `-m <path>`: Mappability BED file (0-based, half-open)

Optional arguments:
- `-gc <path>`: GC bin directory (enables GC-stratified coverage)
  - Expected per-chrom files: `<dir>/chr<N>.tab.gz`
  - Each line is parsed as: `<ignored_token> <1-based_pos> <gc_percent_int>`
- `-chr <c1,c2,...>`: Chromosomes to process (e.g. `1,2,3,5`)  
  Default: `1-22`
- `-t <int>`: Threads per chromosome (default: `40`)
- `-h`: Print help and exit

Notes:
- Coverage is computed only on positions covered by the mappability BED mask.
- Read-length binning uses the alignment query length (`l_qseq`) and bins are `[edge[i], edge[i+1])`, with the last bin extending to infinity.

### splitBAM

Split a BAM file into multiple chunks containing roughly equal numbers of reads for parallel processing.

```bash
./splitBAM <input_bam> <output_prefix> <num_chunks> [num_threads]
```

Parameters:
- `input_bam`: Path to the input BAM file
- `output_prefix`: Prefix for the output BAM files
- `num_chunks`: Number of chunks to split the BAM file into
- `num_threads`: Number of threads to use for writing (default: 4)

### analyzeVCF

Analyze VCF files to generate statistics on non-reference alleles across chromosomes.

```bash
./analyzeVCF [options]
```

Options:
- `-sample <sample_name>`: Sample name to analyze (required)
- `-out_folder <output_folder>`: Output folder for results (default: current directory)
- `-vcf_pattern <pattern>`: Path pattern to VCF files with {chr} placeholder (required)
- `-filter_pattern <pattern>`: Path pattern to filter BED files with {chr} placeholder (required)
- `-reference <genome>`: Genome reference (hg19, hg38, or t2t) (required)
- `-block_size <size>`: Size in base pairs for heterozygosity distribution blocks (default: 10000)
- `-output_distributions`: Output additional distribution files for analysis
- `-help`: Display help message

## Output Files

### analyzeBAM

When used without `-remove_dups`:
- `summary_stats.[ontarget_<bed_name>].L<min_len>MQ<min_map_qual>.txt`: Summary statistics
- `read_length_distribution.<bam_name>.[ontarget_<bed_name>].L<min_len>MQ<min_map_qual>.tsv`: Length distribution
- `<bam_name>.[ontarget].L<min_len>MQ<min_map_qual>.bam`: Filtered BAM
- `<bam_name>.[ontarget].L<min_len>MQ<min_map_qual>.bam.bai`: BAM index

When used with `-remove_dups`:
- `summary_stats.[ontarget_<bed_name>].uniq.L<min_len>MQ<min_map_qual>.txt`: Summary with deduplication stats
- `read_length_distribution.<bam_name>.[ontarget_<bed_name>].uniq.L<min_len>MQ<min_map_qual>.tsv`: Length distribution including unique reads
- `<bam_name>.[ontarget].uniq.L<min_len>MQ<min_map_qual>.bam`: Filtered and deduplicated BAM
- `<bam_name>.[ontarget].uniq.L<min_len>MQ<min_map_qual>.bam.bai`: BAM index

### dedupBAM

- `<bam_name>.uniq.L<min_len>MQ<min_map_qual>.bam`: Deduplicated BAM file with XP and XY tags
- `<bam_name>.uniq.L<min_len>MQ<min_map_qual>.bam.bai`: BAM index

### filterBAM

- `<bam_name>.<suffix>.bam`: Filtered BAM file with specified deamination patterns
- `<bam_name>.<suffix>.bam.bai`: BAM index

### deamBAM

- `substitution_patterns.<bam_name>.txt`: Overall substitution pattern statistics
- `substitution_patterns_<read_length>.<bam_name>.txt`: Detailed substitution patterns by read length
- `read_length_distribution.<bam_name>.txt`: Read length distribution

### analyzeBAMCoverage

For each processed chromosome `<chr>` (e.g. `1`) the tool writes:

- `coverage.<bam_name>.chr<chr>.tab`  
  Two columns: `Cov<TAB>Total` (depth histogram for all reads)

- `coverage.<bam_name>.chr<chr>.<lo>_<hi>.tab`  
  Depth histogram restricted to reads whose read length falls in that bin.

If `-gc` is provided, it also writes GC-stratified depth histograms:

- `gc.<bam_name>.chr<chr>.tab`  
  Table with header `Cov` then GC bins `0..10` plus `GC_missing`.
  Each row gives counts for a depth at each GC bin.

- `gc.<bam_name>.chr<chr>.<lo>_<hi>.tab`  
  Same as above, stratified by read-length bin.

Additionally, “global” (all processed chromosomes combined) files are written:

- `coverage.<bam_name>.tab`
- `coverage.<bam_name>.<lo>_<hi>.tab`
- (if GC enabled) `gc.<bam_name>.tab` and `gc.<bam_name>.<lo>_<hi>.tab`

### splitBAM

- `<output_prefix>.0.bam`, `<output_prefix>.1.bam`, ...: Split BAM files

### analyzeVCF

- `non_ref_results.txt`: Comprehensive statistics per chromosome
- `het_gq_dist.txt`: Genotype quality distribution for heterozygous sites (with `-output_distributions`)
- `alt_gq_dist.txt`: Genotype quality distribution for homozygous alternative sites (with `-output_distributions`)
- `het_dist_<chr>_<block_size>.tab`: Heterozygosity distribution per chromosome (with `-output_distributions`)

## Notes

- Input BAM files must be coordinate-sorted
- Output BAM files are automatically indexed
