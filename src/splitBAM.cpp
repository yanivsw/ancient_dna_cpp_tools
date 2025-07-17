#include <htslib/sam.h>
#include <iostream>
#include <vector>
#include <string>
#include <thread>
#include <future>
#include <fstream>

#include "bam_processing.h"

#define ALN_BUFFER_SIZE 1000000

int write_chunk(
    uint64_t start_read,
    uint64_t chunk_size,
    std::string out_file_name,
    const std::string &input_bam)
{
    bam_file_config_t input_bam_config = {};
    bam_constructor(input_bam, &input_bam_config, "r");

    bam_file_config_t output_bam_config = {};
    bam_constructor(out_file_name, &output_bam_config, "wb");
    output_bam_config.header = bam_hdr_dup(input_bam_config.header);

    // Write the header
    if (sam_hdr_write(output_bam_config.bam_file, output_bam_config.header) < 0)
    {
        std::cerr << "Failed to write BAM header to output file: " << out_file_name << std::endl;
        bam_destructor(&input_bam_config);
        bam_destructor(&output_bam_config);
        return 1;
    }

    // Allocate memory for the alignment
    bam1_t *alignment = bam_init1();
    std::vector<bam1_t *> alignment_buffer;

    uint64_t current_read = 0;
    while (sam_read1(input_bam_config.bam_file, input_bam_config.header, alignment) >= 0)
    {
        if (current_read < start_read)
        {
            current_read++;
            continue;
        }
        else if (current_read >= start_read && current_read < start_read + chunk_size)
        {
            current_read++;
            if (sam_write1(output_bam_config.bam_file, output_bam_config.header, alignment) < 0)
            {
                std::cerr << "Error writing alignment to: " << out_file_name << std::endl;
                bam_destroy1(alignment);
                bam_destructor(&input_bam_config);
                bam_destructor(&output_bam_config);
                return 1;
            }
        }
        else if (current_read >= start_read + chunk_size)
        {
            break;
        }
    }

    bam_destroy1(alignment);
    bam_destructor(&input_bam_config);
    bam_destructor(&output_bam_config);
    return 0;
}

uint64_t count_reads_in_bam(const std::string &input_bam)
{
    bam_file_config_t input_bam_config = {};
    if (!bam_constructor(input_bam, &input_bam_config, "r"))
    {
        std::cerr << "Error opening input BAM file: " << input_bam << std::endl;
        return 0;
    }

    // Count the reads
    uint64_t total_reads = 0;
    bam1_t *alignment = bam_init1();

    while (sam_read1(input_bam_config.bam_file, input_bam_config.header, alignment) >= 0)
    {
        total_reads++;
    }

    // Clean up
    bam_destroy1(alignment);
    bam_destructor(&input_bam_config);

    return total_reads;
}

void split_bam_file(
    const std::string &input_bam,
    const std::string &output_prefix,
    int num_chunks,
    int num_threads)
{
    // Check if BAM index exists, create if not
    if (!bam_index_exists(input_bam))
    {
        if (create_bam_index(input_bam) != 0)
        {
            std::cerr << "Failed to create BAM index. Exiting." << std::endl;
            std::exit(1);
        }
    }

    // Get the total number of reads
    uint64_t total_reads = count_reads_in_bam(input_bam);
    if (total_reads == 0)
    {
        std::cerr << "No reads found in BAM file: " << input_bam << std::endl;
        std::exit(1);
    }

    std::cout << "Total reads in BAM file: " << total_reads << std::endl;

    // Calculate the chunk size
    uint64_t chunk_size = total_reads / num_chunks;

    // Write the chunks
    std::vector<std::future<int>> futures;

    for (int i = 0; i < num_chunks; ++i)
    {
        std::string out_file_name = output_prefix + "." + std::to_string(i) + ".bam";
        uint64_t start_read = i * chunk_size;
        uint64_t end_read = (i == num_chunks - 1) ? total_reads : (i + 1) * chunk_size;

        // Wait if we've reached the thread limit
        if ((int)futures.size() >= num_threads)
        {
            // Wait for the first future to complete
            int result = futures[0].get();
            if (result != 0)
            {
                std::cerr << "Error: write_chunk failed. Terminating program." << std::endl;
                std::exit(1);
            }
            futures.erase(futures.begin());
        }

        futures.push_back(std::async(std::launch::async, write_chunk, start_read, end_read - start_read, out_file_name, input_bam));
    }

    // Wait for remaining futures
    for (auto &future : futures)
    {
        int result = future.get();
        if (result != 0)
        {
            std::cerr << "Error: write_chunk failed. Terminating program." << std::endl;
            std::exit(1);
        }
    }
}

int print_help()
{
    std::cout << "Usage: splitBAM <input_bam> <output_prefix> <num_chunks> [num_threads]" << std::endl;
    std::cout << "Splits a BAM file into multiple chunks." << std::endl;
    std::cout << "Parameters:" << std::endl;
    std::cout << "  input_bam: Path to the input BAM file." << std::endl;
    std::cout << "  output_prefix: Prefix for the output BAM files." << std::endl;
    std::cout << "  num_chunks: Number of chunks to split the BAM file into." << std::endl;
    std::cout << "  num_threads: Number of threads to use for writing (default: 4)." << std::endl;
    return 0;
}

int main(int argc, char *argv[])
{
    if (argc < 4)
    {
        std::cerr << "Error: Not enough arguments provided." << std::endl;
        print_help();
        return 1;
    }

    std::string input_bam = argv[1];
    std::string output_prefix = argv[2];
    uint32_t num_chunks = std::stoi(argv[3]);
    uint32_t num_threads = (argc >= 5) ? std::stoi(argv[4]) : 4;  // Default to 4 threads

    if (num_chunks <= 0)
    {
        std::cerr << "Number of chunks must be greater than 0." << std::endl;
        return 1;
    }

    if (num_threads <= 0)
    {
        std::cerr << "Number of threads must be greater than 0." << std::endl;
        return 1;
    }

    split_bam_file(input_bam, output_prefix, num_chunks, num_threads);

    std::cout << "BAM file split into " << num_chunks << " chunks." << std::endl;
    return 0;
}
