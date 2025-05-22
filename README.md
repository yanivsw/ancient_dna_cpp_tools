# Ancient DNA C++ Tools

A collection of high-performance C++ tools for analyzing and processing ancient DNA sequencing data in BAM format.

## Tools

This toolkit includes three specialized tools:

- **analyzeBAM**: Analyze BAM files to generate comprehensive statistics on alignments with filtering options
- **filterBAM**: Filter BAM files based on deamination patterns (C→T substitutions)
- **deamBAM**: Analyze damage patterns in ancient DNA, particularly terminal substitution frequencies

## Installation

### Prerequisites

- C++17 compiler (GCC recommended)
- HTSlib (1.4 or newer)
- Required libraries: libz, libbz2, liblzma, libpthread

### Building from source

```bash
git clone [repository-url]
cd ancient_dna_cpp_tools
make
```

## Usage

### analyzeBAM

Analyze one or more BAM files and generate statistics with optional filtering.

```bash
./analyzeBAM [options] <BAM files>
```

Options:
- `-out_folder <path>`: Output folder (default: current directory)
- `-min_len <length>`: Minimum alignment length (default: 35)
- `-min_map_qual <quality>`: Minimum mapping quality (default: 0)
- `-targetfile <path>`: BED file with target regions (tab-separated, 0-based)
- `-paired`: Do not disregard paired reads
- `-check_dups`: Generate duplication stats for post bam-rmdup BAM files
- `-count_f`: Ignore filter (QC failed) flag
- `-help`: Display help message

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

## Output Files

### analyzeBAM

- `summary_stats.L{min_len}MQ{min_map_qual}.txt`: Summary statistics for all processed BAM files
- `summary_stats.{bam_name}.L{min_len}MQ{min_map_qual}.txt`: Summary statistics for individual BAM files
- `read_length_distribution.{bam_name}.L{min_len}MQ{min_map_qual}.tsv`: Read length distribution
- Filtered BAM files: `{bam_name}.L{min_len}MQ{min_map_qual}.bam`

### filterBAM

- Filtered BAM files with specified suffix

### deamBAM

- `substitution_patterns.{bam_name}.txt`: Overall substitution pattern statistics
- `substitution_patterns_{read_length}.{bam_name}.txt`: Detailed substitution patterns by read length
- `read_length_distribution.{bam_name}.txt`: Read length distribution

## Examples

Analyze BAM files with minimum length 35 and mapping quality 30:
```bash
./analyzeBAM -min_len 35 -min_map_qual 30 sample1.bam sample2.bam
```

Filter for reads with 5' C→T deamination at positions 0, 1, and 2:
```bash
./filterBAM -p5 0,1,2 -suffix deam5 sample.bam
```

Analyze damage patterns using 16 threads:
```bash
./deamBAM -threads 16 -min_len 35 -max_len 150 sample.bam
```
