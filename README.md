# Ancient DNA C++ Tools

A collection of high-performance C++ tools for analyzing and processing ancient DNA sequencing data in BAM and VCF formats.

## Tools

This toolkit includes six specialized tools:

- **analyzeBAM**: Analyze BAM files to generate comprehensive statistics on alignments with filtering and optional duplicate removal
- **dedupBAM**: Multi-threaded PCR duplicate removal tool optimized for paired-end sequencing data
- **filterBAM**: Filter BAM files based on deamination patterns (C→T substitutions)
- **deamBAM**: Analyze damage patterns in ancient DNA, particularly terminal substitution frequencies
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
- `-help`: Display help message

### dedupBAM

A multi-threaded PCR duplicate removal tool specifically designed for paired-end sequencing data.

```bash
./dedupBAM [options] <BAM file>
```

Options:
- `-out_folder <path>`: Output folder (default: current directory)
- `-min_len <length>`: Minimum alignment length (default: 0)
- `-min_map_qual <quality>`: Minimum mapping quality (default: 0)
- `-ignore_read_groups`: Ignore read group information when identifying duplicates
- `-ignore_length`: Ignore alignment length when identifying duplicates
- `-ignore_qc_fail`: Process reads that failed QC
- `-threads <num>`: Number of threads to use (default: 4)
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

### splitBAM

- `<output_prefix>.0.bam`, `<output_prefix>.1.bam`, ...: Split BAM files

### analyzeVCF

- `non_ref_results.txt`: Comprehensive statistics per chromosome
- `het_gq_dist.txt`: Genotype quality distribution for heterozygous sites (with `-output_distributions`)
- `alt_gq_dist.txt`: Genotype quality distribution for homozygous alternative sites (with `-output_distributions`)
- `het_dist_<chr>_<block_size>.tab`: Heterozygosity distribution per chromosome (with `-output_distributions`)

## Examples

### Basic Analysis and Filtering

Analyze BAM files with minimum length 35 and mapping quality 25:
```bash
./analyzeBAM -min_len 35 -min_map_qual 25 sample1.bam sample2.bam
```

### On-Target Analysis with Duplicate Removal

Analyze reads on target regions and remove duplicates:
```bash
./analyzeBAM -targetfile targets.bed -remove_dups -min_len 35 -min_map_qual 25 sample.bam
```

### Standalone Duplicate Removal

Remove duplicates using 8 threads:
```bash
./dedupBAM -threads 8 -min_len 35 -min_map_qual 25 sample.bam
```

Remove duplicates ignoring read groups:
```bash
./dedupBAM -ignore_read_groups -threads 16 sample.bam
```

### Damage Pattern Analysis

Filter for reads with 5' C→T deamination at positions 0, 1, and 2:
```bash
./filterBAM -p5 0,1,2 -p3 0,-1,-2 -suffix deam3_or_5 sample.bam
```

Analyze damage patterns using 16 threads:
```bash
./deamBAM -threads 16 -min_len 35 -max_len 150 sample.bam
```

### VCF Analysis

Analyze VCF files for sample NA12878 with hg38 reference:
```bash
./analyzeVCF -sample NA12878 \
  -vcf_pattern /path/to/sample.{chr}.vcf.gz \
  -filter_pattern /path/to/filter.{chr}.bed \
  -reference hg38 \
  -out_folder results/
```

### Parallel Processing

Split a BAM file into 10 chunks using 8 threads:
```bash
./splitBAM sample.bam output/sample_chunk 10 8
```

## Notes

- Input BAM files must be coordinate-sorted
- Output BAM files are automatically indexed
