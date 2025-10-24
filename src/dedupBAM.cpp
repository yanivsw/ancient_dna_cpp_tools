#include <htslib/sam.h>
#include <htslib/kstring.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <unordered_map>
#include <sstream>
#include <random>
#include <map>
#include <thread>
#include <mutex>
#include <queue>
#include <condition_variable>

#include "bam_processing.h"
#include "duplicate_handler.h"
#include "types.h"

#define ALN_BUFFER_SIZE 1000000

#define READ_LEN_DIST_SIZE 512
#define VERSION_NUMBER 0.1

void process_chromosome(
    int tid,
    std::string input_bam_name,
    std::string chr_output_bam,
    bool ignore_length,
    bool ignore_read_groups,
    bool ignore_qc_fail,
    uint32_t min_length,
    uint32_t min_map_quality,
    duplication_stats_struct& dup_stats,
    std::map<int, bool>& chromosome_ready,
    std::mutex& ready_mutex,
    std::condition_variable& ready_cv)
{
    bam_file_config_t input_bam_config = {};
    bam_constructor(input_bam_name, &input_bam_config, "r");

    // Create iterator for this chromosome
    hts_itr_t* iter = sam_itr_queryi(input_bam_config.index, tid, 0, input_bam_config.header->target_len[tid]);
    if (!iter)
    {
        std::cerr << "Failed to create iterator for chromosome " << tid << std::endl;
        return;
    }

    bam_file_config_t output_bam_config = {};
    bam_constructor(chr_output_bam, &output_bam_config, "wb");
    output_bam_config.header = bam_hdr_dup(input_bam_config.header);
    sam_hdr_write(output_bam_config.bam_file, output_bam_config.header);

    std::vector<bam1_t*> alignment_buffer;
    std::vector<bam1_t*> deduped_alignment_buffer;
    int32_t current_position = -1;

    bam1_t* alignment = bam_init1();
    while (sam_itr_next(input_bam_config.bam_file, iter, alignment) >= 0)
    {
        // Skip reads that don't meet quality/length criteria
        if ((alignment->core.l_qseq < static_cast<int32_t>(min_length)) ||
            (alignment->core.qual < min_map_quality) ||
            (alignment->core.flag & BAM_FMUNMAP) ||
            ((alignment->core.flag & BAM_FPAIRED) && !(alignment->core.flag & BAM_FPROPER_PAIR)) ||
            ((alignment->core.flag & BAM_FPAIRED) && (alignment->core.mtid != tid)) ||
            ((alignment->core.flag & BAM_FQCFAIL) && (ignore_qc_fail == false)))
        {
            continue;
        }

        // Process buffer when it gets big enough and we've moved past the current position
        if ((current_position != -1) && (alignment->core.pos > current_position) && (alignment_buffer.size() >= ALN_BUFFER_SIZE))
        {
            remove_duplicates(alignment_buffer, deduped_alignment_buffer, ignore_length, ignore_read_groups, &dup_stats);
            if (write_alignment_buffer_to_bam(&output_bam_config, deduped_alignment_buffer) != 0)
            {
                std::cerr << "Failed to write alignment buffer for chromosome " << tid << std::endl;
            }
            deduped_alignment_buffer.clear();
        }

        current_position = alignment->core.pos;

        bam1_t* new_alignment = bam_init1();
        bam_copy1(new_alignment, alignment);
        alignment_buffer.push_back(new_alignment);
    }

    // Process any remaining alignments
    if (!alignment_buffer.empty())
    {
        remove_duplicates(alignment_buffer, deduped_alignment_buffer, ignore_length, ignore_read_groups, &dup_stats);
        if (write_alignment_buffer_to_bam(&output_bam_config, deduped_alignment_buffer) != 0)
        {
            std::cerr << "Failed to write alignment buffer for chromosome " << tid << std::endl;
        }
    }

    bam_destroy1(alignment);
    hts_itr_destroy(iter);

    if (hts_close(output_bam_config.bam_file) < 0)
    {
        std::cerr << "Error closing chromosome BAM file" << std::endl;
    }

    bam_hdr_destroy(output_bam_config.header);

    // Signal that this chromosome is ready for merging
    std::lock_guard<std::mutex> lock(ready_mutex);
    chromosome_ready[tid] = true;
    // std::cout << "Chromosome " << tid << " ready for merging" << std::endl;

    ready_cv.notify_all(); // Notify merger thread
}

int merge_chromosome_bams(
    const std::map<int, std::string>& chr_bam_files,
    bam_file_config_t& output_bam_config,
    std::map<int, bool>& chromosome_ready,
    std::mutex& ready_mutex,
    std::condition_variable& ready_cv,
    bool& all_processing_complete)
{
    // std::cout << "Starting concurrent merge thread..." << std::endl;

    // Process chromosomes in order
    for (const auto& [tid, file_path] : chr_bam_files)
    {
        // Wait until this chromosome is ready
        {
            std::unique_lock<std::mutex> lock(ready_mutex);
            ready_cv.wait(lock, [&]()
                {
                    return chromosome_ready[tid] || all_processing_complete;
                });

            // If all processing is complete but this chromosome isn't ready,
            // there was an error and we should skip it
            if (!chromosome_ready[tid] && all_processing_complete)
            {
                std::cerr << "Warning: Chromosome " << tid << " never became available" << std::endl;
                continue;
            }
        }

        // std::cout << "Merging chromosome " << tid << " from " << file_path << std::endl;

        // Verify the file exists and is readable
        std::ifstream file_check(file_path);
        if (!file_check.good())
        {
            std::cerr << "Warning: Chromosome file not found: " << file_path << std::endl;
            continue;
        }
        file_check.close();

        // Open and merge this chromosome file
        bam_file_config_t chr_config = {};
        bam_constructor(file_path, &chr_config, "r");

        // Copy alignments from this chromosome to output
        bam1_t* aln = bam_init1();
        uint64_t alignment_count = 0;
        std::vector<bam1_t*> alignment_buffer;

        while (sam_read1(chr_config.bam_file, chr_config.header, aln) >= 0)
        {
            bam1_t* new_alignment = bam_init1();
            bam_copy1(new_alignment, aln);
            alignment_buffer.push_back(new_alignment);

            // If buffer is full, process it
            if (alignment_buffer.size() >= ALN_BUFFER_SIZE)
            {
                if (write_alignment_buffer_to_bam(&output_bam_config, alignment_buffer) != 0)
                {
                    std::cerr << "Failed to write alignment buffer" << std::endl;
                    return 1;
                }
            }
            alignment_count++;
        }

        // Process any remaining alignments in the buffer
        if (!alignment_buffer.empty())
        {
            if (write_alignment_buffer_to_bam(&output_bam_config, alignment_buffer) != 0)
            {
                std::cerr << "Failed to write final alignment buffer" << std::endl;
            }

            alignment_count += alignment_buffer.size();
        }

        // std::cout << "  Merged " << alignment_count << " alignments from chromosome " << tid << std::endl;

        // Clean up
        bam_destroy1(aln);
        bam_destructor(&chr_config);
    }

    // std::cout << "Merge thread complete" << std::endl;
    return 0;
}

void print_help()
{
    std::cout << "Usage: dedupBAM [options] <BAM file>\n"
        << "Options:\n"
        << "  -out_folder <path>       Specify the output folder, otherwise files will be written to the current folder\n"
        << "  -min_len <length>        Specify the minimum alignment length [default 0]\n"
        << "  -min_map_qual <quality>  Specify the minimum mapping quality [default 0]\n"
        << "  -ignore_read_groups      Ignore read group information when identifying duplicates\n"
        << "  -ignore_length           Ignore alignment length when identifying duplicates\n"
        << "  -ignore_qc_fail          Ignore QC fail flag when identifying duplicates\n"
        << "  -threads <num>           Specify the number of threads to use [default 4]\n"
        << "  -help                    Display this help message\n";
}

int main(int argc, char* argv[])
{
    std::cout << "dedupBAM.cpp v" << VERSION_NUMBER << " - " << get_date_time() << std::endl;

    bool ignore_length = false;
    bool ignore_read_groups = false;
    bool ignore_qc_fail = false;
    uint32_t min_length = 0;
    uint32_t min_map_quality = 0;
    uint32_t thread_count = 4;

    std::string bam_file = "";
    std::string out_folder = "";
    std::string input_bam_name = "";
    std::string min_len_str = std::to_string(min_length);
    std::string min_map_qual_str = std::to_string(min_map_quality);

    // Save the command line arguments to a string
    std::ostringstream command_line;
    for (int i = 0; i < argc; i++)
    {
        command_line << argv[i] << " ";
    }
    std::string command_line_str = command_line.str();

    for (int i = 0; i < argc; i++)
    {
        if (std::string(argv[i]) == "-out_folder")
        {
            out_folder = argv[i + 1];
            if (out_folder.back() != '/')
            {
                out_folder += "/";
            }
        }
        if (std::string(argv[i]) == "-min_len")
        {
            min_length = std::stoi(argv[i + 1]);
            min_len_str = std::to_string(min_length);
        }
        if (std::string(argv[i]) == "-min_map_qual")
        {
            min_map_quality = std::stoi(argv[i + 1]);
            min_map_qual_str = std::to_string(min_map_quality);
        }
        if (std::string(argv[i]) == "-ignore_read_groups")
        {
            ignore_read_groups = true;
        }
        if (std::string(argv[i]) == "-ignore_length")
        {
            ignore_length = true;
        }
        if (std::string(argv[i]) == "-ignore_qc_fail")
        {
            ignore_qc_fail = true;
        }
        if (std::string(argv[i]) == "-threads")
        {
            thread_count = std::stoi(argv[i + 1]);
        }
        if (std::string(argv[i]).find(".bam") != std::string::npos)
        {
            bam_file = argv[i];
        }
        if (std::string(argv[i]) == "-help" || argc == 1)
        {
            print_help();
            return 0;
        }
    }

    if (bam_file.empty())
    {
        std::cerr << "No BAM file provided" << std::endl;
        return 1;
    }
    else
    {
        std::cout << "Deduplicating " << bam_file << std::endl;
    }

    input_bam_name = extract_file_name(bam_file);

    std::string suffix = std::string(".uniq.L") + min_len_str + "MQ" + min_map_qual_str;
    std::string output_file_location = out_folder + input_bam_name + suffix + ".bam";

    bam_file_config_t input_bam_config = {};
    bam_constructor(bam_file, &input_bam_config, "r");

    bam_file_config_t output_bam_config = {};
    bam_constructor(output_file_location, &output_bam_config, "wb");
    output_bam_config.header = bam_hdr_dup(input_bam_config.header);

    // Update header with command line 
    std::string version_str = std::to_string(VERSION_NUMBER);
    version_str.erase(version_str.find_last_not_of('0') + 1, std::string::npos);
    append_line_to_bam_header(output_bam_config.header, "@PG\tID:dedupBAM\tVN:" + std::string(version_str) + "\tCL:" + command_line_str + "\n");
    if (sam_hdr_write(output_bam_config.bam_file, output_bam_config.header) < 0)
    {
        std::cerr << "Failed to write BAM header to output file: " << output_file_location << std::endl;
        bam_destructor(&input_bam_config);
        bam_destructor(&output_bam_config);
        return 1;
    }

    // Create thread pool
    if (thread_count == 0)
    {
        thread_count = 4;
    }

    std::cout << "Using " << thread_count << " threads for processing" << std::endl;
    std::vector<std::thread> threads;

    duplication_stats_struct dup_stats_total(READ_LEN_DIST_SIZE);
    std::vector<duplication_stats_struct> dup_stats(input_bam_config.header->n_targets, duplication_stats_struct(READ_LEN_DIST_SIZE));

    // Merge chromosome BAMs
    std::map<int, std::string> chr_bam_files;
    for (int tid = 0; tid < input_bam_config.header->n_targets; ++tid)
    {
        std::string chr_name = input_bam_config.header->target_name[tid];
        std::string chr_output_bam = out_folder + input_bam_name + ".chr" + chr_name + ".bam";
        chr_bam_files[tid] = chr_output_bam;
    }

    // Map to track which chromosomes are ready for merging
    std::map<int, bool> chromosome_ready;
    std::mutex ready_mutex;
    std::condition_variable ready_cv;
    bool all_processing_complete = false;

    // Start merger thread
    std::thread merger_thread(merge_chromosome_bams,
        std::ref(chr_bam_files),
        std::ref(output_bam_config),
        std::ref(chromosome_ready),
        std::ref(ready_mutex),
        std::ref(ready_cv),
        std::ref(all_processing_complete));

    // Process chromosomes in batches based on thread count
    for (int tid = 0; tid < input_bam_config.header->n_targets; tid += thread_count)
    {
        for (int i = 0; i < static_cast<int>(thread_count) && (tid + i) < input_bam_config.header->n_targets; ++i)
        {
            int current_tid = tid + i;
            std::string chr_name = input_bam_config.header->target_name[current_tid];
            std::string chr_output_bam = chr_bam_files[current_tid];

            // Launch thread for this chromosome
            threads.emplace_back(process_chromosome,
                current_tid,
                bam_file,
                chr_output_bam,
                ignore_length,
                ignore_read_groups,
                ignore_qc_fail,
                min_length,
                min_map_quality,
                std::ref(dup_stats[current_tid]),
                std::ref(chromosome_ready),
                std::ref(ready_mutex),
                std::ref(ready_cv));
        }

        // Clear previous threads
        for (auto& t : threads)
        {
            if (t.joinable())
            {
                t.join();
            }
        }
        threads.clear();
    }

    // Signal that all processing is complete
    {
        std::lock_guard<std::mutex> lock(ready_mutex);
        all_processing_complete = true;
    }
    ready_cv.notify_all();

    // Wait for remaining threads
    for (auto& t : threads)
    {
        if (t.joinable())
        {
            t.join();
        }
    }

    // // Wait for merger thread to complete
    merger_thread.join();

    for (const auto& stats : dup_stats)
    {
        // Aggregate stats if needed
        dup_stats_total.single += stats.single;
        dup_stats_total.uniq += stats.uniq;
        dup_stats_total.total += stats.total;
    }

    // Delete temporary chromosome BAM files
    for (const auto& [tid, file_path] : chr_bam_files)
    {
        if (std::remove(file_path.c_str()) != 0)
        {
            std::cerr << "Warning: Failed to delete temporary file " << file_path << std::endl;
        }
    }

    // Clean up
    bam_destructor(&input_bam_config);
    bam_destructor(&output_bam_config);

    if (sam_index_build(output_file_location.c_str(), 0) < 0)
    {
        std::cerr << "Failed to generate BAM index" << std::endl;
        return 1;
    }

    std::cout << "\nAlignments processed: " << dup_stats_total.total << std::endl;
    std::cout << "Alignments written (deduplicated): " << dup_stats_total.uniq << std::endl;
    std::cout << "Singleton alignments: " << dup_stats_total.single << std::endl;
    std::cout << "QC failed alignments: " << dup_stats_total.qc_failed << std::endl;
    std::cout << "Deduplication complete. Output written to " << output_file_location << std::endl;

    return 0;
}
