#include <iostream>
#include <variant>
#include <htslib/sam.h>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <cstring>

#include "bed_processing.h"
#include "bam_processing.h"
#include "duplicate_handler.h"
#include "pair_handler.h"
#include "utils.h"
#include "types.h"

#define VERSION_NUMBER 0.8

#define ALN_BUFFER_SIZE 1000000
#define READ_LEN_DIST_SIZE 512

std::unordered_set<std::string> off_genome_chromosomes = {
    "phiX", "NC_007605", "hs37d5",
    "GL000207.1", "GL000226.1", "GL000229.1", "GL000231.1", "GL000210.1", "GL000239.1",
    "GL000235.1", "GL000201.1", "GL000247.1", "GL000245.1", "GL000197.1", "GL000203.1",
    "GL000246.1", "GL000249.1", "GL000196.1", "GL000248.1", "GL000244.1", "GL000238.1",
    "GL000202.1", "GL000234.1", "GL000232.1", "GL000206.1", "GL000240.1", "GL000236.1",
    "GL000241.1", "GL000243.1", "GL000242.1", "GL000230.1", "GL000237.1", "GL000233.1",
    "GL000204.1", "GL000198.1", "GL000208.1", "GL000191.1", "GL000227.1", "GL000228.1",
    "GL000214.1", "GL000221.1", "GL000209.1", "GL000218.1", "GL000220.1", "GL000213.1",
    "GL000211.1", "GL000199.1", "GL000217.1", "GL000216.1", "GL000215.1", "GL000205.1",
    "GL000219.1", "GL000224.1", "GL000223.1", "GL000195.1", "GL000212.1", "GL000222.1",
    "GL000200.1", "GL000193.1", "GL000194.1", "GL000225.1", "GL000192.1" };

void write_length_dist_file(
    const std::string& len_dist_file_name,
    bool is_target_file,
    const summary_stats_struct& summary_stats,
    const duplication_stats_struct& dup_stats,
    bool remove_dups = false)
{
    std::ofstream len_dist_file(len_dist_file_name, std::ios::app);
    for (size_t j = 0; j < summary_stats.length_distribution_matrix.size(); j++)
    {
        const auto& row = summary_stats.length_distribution_matrix[j];
        len_dist_file << j << "\t";
        for (size_t i = 0; i < DIST_FILTERED_LEN_MAPPED_QUAL_ONTARGET; i++)
        {
            len_dist_file << row[i] << "\t";

            if (i == DIST_FILTERED_LEN_MAPPED_QUAL)
            {
                if (row[DIST_FILTERED_LEN] > 0)
                {
                    len_dist_file << (static_cast<double>(row[DIST_FILTERED_LEN_MAPPED_QUAL]) / row[DIST_FILTERED_LEN]) * 100;
                }
                else
                {
                    len_dist_file << "NA";
                }
                if (is_target_file)
                {
                    len_dist_file << "\t" << row[i + 1];
                }
                if (remove_dups)
                {
                    if (dup_stats.length_distribution_uniq[j] > 0)
                    {
                        len_dist_file << "\t" << dup_stats.length_distribution_uniq[j]
                            << "\t" << (static_cast<double>(dup_stats.length_distribution_total[j]) / dup_stats.length_distribution_uniq[j])
                            << "\t" << dup_stats.length_distribution_single[j];
                    }
                    else
                    {
                        len_dist_file << "\tNA\tNA\tNA";
                    }
                }
                len_dist_file << "\n";
            }
        }
    }
    len_dist_file.close();
}

void append_summary_stats_to_file(
    const std::string& input_bam_name,
    const std::string& summary_stats_file_name,
    const summary_stats_struct& summary_stats,
    const duplication_stats_struct& dup_stats,
    bool is_target_file,
    bool remove_dups = false)
{
    std::ofstream summary_stats_file(summary_stats_file_name, std::ios::app);
    summary_stats_file << input_bam_name << "\t"
        << summary_stats.raw_alignments << "\t"
        << summary_stats.merged_alignments << "\t"
        << summary_stats.filtered_alignments << "\t"
        << summary_stats.filtered_len_alignments << "\t"
        << summary_stats.filtered_len_mapped_alignments << "\t"
        << summary_stats.filtered_len_mapped_qual_alignments << "\t";

    if (summary_stats.filtered_len_alignments > 0)
    {
        summary_stats_file << (static_cast<double>(summary_stats.filtered_len_mapped_qual_alignments) / summary_stats.filtered_len_alignments) * 100;
    }
    else
    {
        summary_stats_file << "NA";
    }

    if (is_target_file)
    {
        summary_stats_file << "\t" << summary_stats.filtered_len_mapped_qual_target_alignments;
    }

    if (remove_dups)
    {
        summary_stats_file << "\t" << dup_stats.uniq << "\t"
            << (static_cast<double>(dup_stats.total) / dup_stats.uniq) << "\t"
            << dup_stats.single;
    }

    summary_stats_file << "\n";
    summary_stats_file.close();
}

bool check_on_target(
    bam1_t* alignment,
    const std::string& current_chromosome,
    const bed_file_t& bed_region_map,
    uint64_t& bed_region_ptr,
    uint32_t& target_end)
{
    if (bed_region_map.find(current_chromosome) == bed_region_map.end())
    {
        return false;
    }

    int aln_len = get_alignment_length(alignment);
    int aln_start = alignment->core.pos;
    int aln_end = aln_start + aln_len;

    // Increment the target region pointer until we find one that could overlap with the alignment
    while ((aln_start > static_cast<int32_t>(target_end)) && 
           (bed_region_ptr < bed_region_map.at(current_chromosome).size()))
    {
        bed_region_ptr++;
        if (bed_region_ptr < bed_region_map.at(current_chromosome).size())
        {
            target_end = bed_region_map.at(current_chromosome)[bed_region_ptr].end;
        }
    }

    // Check if the alignment overlaps with any target region
    for (size_t i = bed_region_ptr; i < bed_region_map.at(current_chromosome).size(); i++)
    {
        int bed_region_start = bed_region_map.at(current_chromosome)[i].start;
        int bed_region_end = bed_region_map.at(current_chromosome)[i].end;

        // The alignment is on target
        if ((aln_end > bed_region_start) && (aln_start < bed_region_end))
        {
            return true;
        }

        // The alignment doesn't overlap with any target regions
        if (bed_region_start > aln_end)
        {
            break;
        }
    }

    return false;
}

void print_help()
{
    std::cout << "Usage: analyzeBAM [options] <BAM files>\n"
        << "Options:\n"
        << "  -out_folder <path>          Specify the output folder, otherwise files will be written to the current folder\n"
        << "  -min_len <length>           Specify the minimum alignment length [default 35]\n"
        << "  -min_map_qual <quality>     Specify the minimum mapping quality [default 0]\n"
        << "  -targetfile <path>          Specify the target file path. ***NB*** This expects a bed file (i.e. tab separated, 0-based)\n"
        << "  -paired                     Do not disregard paired reads\n"
        << "  -count_f                    Ignore filter (QC failed) flag\n"
        << "  -remove_dups                Remove duplicate reads\n"
        << "  -help                       Display this help message\n";
}

int main(int argc, char* argv[])
{
    int32_t min_len = 35;
    int32_t min_map_qual = 0;
    std::vector<std::string> bam_files;
    std::string target_file_path = "";
    std::string target_file_name = "";
    bool is_target_file = false;
    bool paired = false;
    bool count_f = false;
    std::string min_len_str;
    std::string min_map_qual_str;
    bed_file_t bed_region_map;
    std::string out_folder = "";

    bool remove_dups = false;

    std::cout << "analyzeBAM.cpp v" << VERSION_NUMBER << " - " << get_date_time() << std::endl;
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
            min_len = std::stoi(argv[i + 1]);
            min_len_str = std::to_string(min_len);
        }
        if (std::string(argv[i]) == "-min_map_qual")
        {
            min_map_qual = std::stoi(argv[i + 1]);
            min_map_qual_str = std::to_string(min_map_qual);
        }
        if (std::string(argv[i]) == "-targetfile")
        {
            target_file_path = argv[i + 1];
            target_file_name = extract_file_name(target_file_path);
            is_target_file = true;
        }
        if (std::string(argv[i]) == "-paired")
        {
            paired = true;
        }
        if (std::string(argv[i]).find(".bam") != std::string::npos)
        {
            bam_files.push_back(argv[i]);
        }
        if (std::string(argv[i]) == "-remove_dups")
        {
            remove_dups = true;
        }
        if (std::string(argv[i]) == "-count_f")
        {
            count_f = true;
        }
        if ((std::string(argv[i]) == "-help") || (argc == 1))
        {
            print_help();
            return 0;
        }
    }

    if (bam_files.size() == 0)
    {
        std::cerr << "No BAM files provided" << std::endl;
        return 1;
    }

    std::string summary_stats_file_name = out_folder + "summary_stats" +
        (is_target_file ? ".ontarget_" + target_file_name : "") +
        (remove_dups ? ".uniq" : "") +
        ".L" + min_len_str + "MQ" + min_map_qual_str + ".txt";

    std::ofstream summary_stats_file(summary_stats_file_name);
    summary_stats_file << "#file\traw\tmerged\tfilter_passed\tL" << min_len_str
        << "\tmappedL" << min_len_str
        << "\tmappedL" << min_len_str << "MQ" << min_map_qual_str
        << "\t%mappedL" << min_len_str << "MQ" << min_map_qual_str;
    if (is_target_file)
    {
        summary_stats_file << "\ttargetL" << min_len_str << "MQ" << min_map_qual_str;
    }
    if (remove_dups)
    {
        summary_stats_file << "\tunique\taverage_dups\tsingletons";
    }
    summary_stats_file << "\n";
    summary_stats_file.close();

    // Read target file to an unordered_map if it exists
    if (is_target_file)
    {
        read_bed_file(target_file_path, bed_region_map);
    }

    // Loop over bam files
    for (const auto& bam_file : bam_files)
    {
        // Extract the file name without the extension and path
        std::string input_bam_name;
        input_bam_name = extract_file_name(bam_file);

        // Output file name
        std::string suffix = std::string(is_target_file ? ".ontarget" : "") + std::string(remove_dups ? ".uniq" : "") + ".L" + min_len_str + "MQ" + min_map_qual_str;
        std::string output_file_location = out_folder + input_bam_name + suffix + ".bam";

        // Summary stats and length distribution file names
        std::string target_suffix = std::string(is_target_file ? ".ontarget_" + target_file_name : "") + std::string(remove_dups ? ".uniq" : "") + ".L" + min_len_str + "MQ" + min_map_qual_str;
        std::string summary_stats_file_name_2 = out_folder + "summary_stats." + input_bam_name + target_suffix + ".txt";
        std::string len_dist_file_name = out_folder + "read_length_distribution." + input_bam_name + target_suffix + ".tsv";

        std::ofstream summary_stats_file_2(summary_stats_file_name_2);
        summary_stats_file_2 << "#file\traw\tmerged\tfilter_passed\tL" << min_len_str
            << "\tmappedL" << min_len_str
            << "\tmappedL" << min_len_str << "MQ" << min_map_qual_str
            << "\t%mappedL" << min_len_str << "MQ" << min_map_qual_str
            << "\ttargetL" << min_len_str << "MQ" << min_map_qual_str;
        if (remove_dups)
        {
            summary_stats_file_2 << "\tunique\taverage_dups\tsingletons";
        }
        summary_stats_file_2 << "\n";
        summary_stats_file_2.close();

        std::ofstream len_dist_file(len_dist_file_name);
        len_dist_file << "Length\traw\tmerged\tfilter_passed\tL" << min_len_str
            << "\tmappedL" << min_len_str
            << "\tmappedL" << min_len_str << "MQ" << min_map_qual_str
            << "\t%mappedL" << min_len_str << "MQ" << min_map_qual_str
            << "\ttargetL" << min_len_str << "MQ" << min_map_qual_str;
        if (remove_dups)
        {
            len_dist_file << "\tunique\taverage_dups\tsingletons";
        }
        len_dist_file << "\n";
        len_dist_file.close();

        std::cout << "\nProcessing " << bam_file << "\n";

        bam_file_config_t input_bam_config = {};
        bam_constructor(bam_file, &input_bam_config, "r");

        bam_file_config_t output_bam_config = {};
        bam_constructor(output_file_location, &output_bam_config, "wb");
        output_bam_config.header = bam_hdr_dup(input_bam_config.header);

        // Get the chromosome names from the BAM header
        std::vector<std::string> chromosomes;
        get_chromosomes_from_bam(chromosomes, bam_file);

        std::string version_str = std::to_string(VERSION_NUMBER);
        version_str.erase(version_str.find_last_not_of('0') + 1, std::string::npos);
        append_line_to_bam_header(output_bam_config.header, "@PG\tID:analyzeBAM\tVN:" + std::string(version_str) + "\tCL:" + command_line_str + "\n");

        // Write the header
        if (sam_hdr_write(output_bam_config.bam_file, output_bam_config.header) < 0)
        {
            std::cerr << "Failed to write BAM header to output file: " << output_file_location << std::endl;
            bam_destructor(&input_bam_config);
            bam_destructor(&output_bam_config);
            return 1;
        }

        summary_stats_struct summary_stats(READ_LEN_DIST_SIZE);
        duplication_stats_struct dup_stats(READ_LEN_DIST_SIZE);

        uint64_t counter = 0;
        uint64_t qc_fail_reads = 0;
        uint64_t unmapped_reads = 0;
        uint64_t off_genome_reads = 0;

        uint64_t bed_region_ptr = 0;
        uint32_t target_end = 0;
        bool on_target = false;
        int32_t current_tid = -1;
        int32_t current_position = -1;
        std::string current_chromosome = "";

        std::vector<bam1_t*> alignment_buffer;
        std::vector<bam1_t*> deduped_alignment_buffer;

        paired_read_tracker_t pair_tracker;

        // bam flags: "qdfs21RrUuPp"
        // Read each alignment
        bam1_t* alignment = bam_init1();
        while (sam_read1(input_bam_config.bam_file, input_bam_config.header, alignment) >= 0)
        {
            on_target = false;
            counter++;
            if ((counter % 100000) == 0)
            {
                std::cout << "\r" << counter << " alignments processed";
                std::cout.flush();
            }

            // If this is a new chromosome, update current chromosome and reset target region pointer
            // If removing duplicates, process the buffered alignments for the previous chromosome
            if (current_tid != alignment->core.tid)
            {
                current_chromosome = chromosomes[alignment->core.tid];

                if (remove_dups && !alignment_buffer.empty())
                {
                    if (paired)
                    {
                        rescue_failed_mates(alignment_buffer, pair_tracker);
                    }

                    remove_duplicates(alignment_buffer, deduped_alignment_buffer, false, true, &dup_stats);

                    if (write_alignment_buffer_to_bam(&output_bam_config, deduped_alignment_buffer) != 0)
                    {
                        std::cerr << "Failed to write alignment buffer" << std::endl;
                        return 1;
                    }

                    if (paired)
                    {
                        cleanup_pair_tracker(pair_tracker);
                    }
                }

                if (is_target_file)
                {
                    if (bed_region_map[current_chromosome].size() > 0)
                    {
                        bed_region_ptr = 0;
                        target_end = bed_region_map[current_chromosome][0].end;
                    }
                    else if (current_chromosome == "MT" || current_chromosome == "chrM")
                    {
                        // Skip processing for mitochondria if not in the target file
                        continue;
                    }
                }

                current_tid = alignment->core.tid;
                current_position = -1; // Reset position when chromosome changes
            }

            // Skip if the alignment is the second read in a pair and we are not doing paired read processing
            if (((alignment->core.flag & BAM_FREAD2) != 0) && (paired == false))
            {
                continue;
            }
            summary_stats.raw_alignments++;
            summary_stats.length_distribution_matrix[alignment->core.l_qseq][DIST_RAW]++;

            // Skip if not merged or not paired
            if ((((alignment->core.flag & BAM_FREAD1) != 0) && (paired == false)) ||
                (((alignment->core.flag & BAM_FPAIRED) == 0) && (paired == true)))
            {
                continue;
            }
            summary_stats.merged_alignments++;
            summary_stats.length_distribution_matrix[alignment->core.l_qseq][DIST_MERGED]++;

            // Skip if the alignment failed quality control
            if (((alignment->core.flag & BAM_FQCFAIL) != 0) && (count_f == false))
            {
                qc_fail_reads++;
                continue;
            }
            summary_stats.filtered_alignments++;
            summary_stats.length_distribution_matrix[alignment->core.l_qseq][DIST_FILTERED]++;

            // Skip if alignment length < min_len
            if (alignment->core.l_qseq < min_len)
            {
                continue;
            }
            summary_stats.filtered_len_alignments++;
            summary_stats.length_distribution_matrix[alignment->core.l_qseq][DIST_FILTERED_LEN]++;

            // Skip if unmapped
            if ((alignment->core.flag & BAM_FUNMAP) != 0)
            {
                unmapped_reads++;
                continue;
            }

            // Skip if alignment is mapped to a contig in the off_genome_chromosomes set
            if (off_genome_chromosomes.count(current_chromosome) > 0)
            {
                off_genome_reads++;
                continue;
            }
            summary_stats.filtered_len_mapped_alignments++;
            summary_stats.length_distribution_matrix[alignment->core.l_qseq][DIST_FILTERED_LEN_MAPPED]++;

            bool passes_mapq = alignment->core.qual >= min_map_qual;
            bool is_paired = paired && (alignment->core.flag & BAM_FPAIRED);

            // Write the alignment buffer to the output BAM file if it is full or we have moved past the current position
            if ((current_position != -1) && (alignment->core.pos > current_position) && (alignment_buffer.size() >= ALN_BUFFER_SIZE))
            {
                if (paired)
                {
                    rescue_failed_mates(alignment_buffer, pair_tracker);
                }

                if (remove_dups)
                {
                    remove_duplicates(alignment_buffer, deduped_alignment_buffer, false, true, &dup_stats);
                    if (write_alignment_buffer_to_bam(&output_bam_config, deduped_alignment_buffer) != 0)
                    {
                        std::cerr << "Failed to write alignment buffer to " << output_bam_config.bam_file_location << std::endl;
                        bam_destructor(&input_bam_config);
                        bam_destructor(&output_bam_config);
                        return 1;
                    }
                }
                else
                {
                    if (write_alignment_buffer_to_bam(&output_bam_config, alignment_buffer) != 0)
                    {
                        std::cerr << "Failed to write alignment buffer to " << output_bam_config.bam_file_location << std::endl;
                        bam_destructor(&input_bam_config);
                        bam_destructor(&output_bam_config);
                        return 1;
                    }
                }
            }

            current_position = alignment->core.pos;

            if (is_target_file)
            {
                on_target = check_on_target(alignment, current_chromosome, bed_region_map, 
                                           bed_region_ptr, target_end);
            }
            else
            {
                on_target = true;  // Consider all reads "on-target" if no target file
            }

            // Handle paired or unpaired reads
            if (is_paired)
            {
                // Check if this is part of a valid pair
                if ((alignment->core.flag & BAM_FMUNMAP) ||
                    (alignment->core.mtid != alignment->core.tid) ||
                    !(alignment->core.flag & BAM_FPROPER_PAIR))
                {
                    continue;
                }

                bool pair_passed = handle_paired_read(alignment, passes_mapq, 
                                                      alignment_buffer, pair_tracker,
                                                      on_target);
                
                // Update stats only if this read originally passed (not a rescued mate)
                if (pair_passed)
                {
                    summary_stats.filtered_len_mapped_qual_alignments++;
                    summary_stats.length_distribution_matrix[alignment->core.l_qseq][DIST_FILTERED_LEN_MAPPED_QUAL]++;
                    
                    if (on_target && is_target_file)
                    {
                        summary_stats.filtered_len_mapped_qual_target_alignments++;
                        summary_stats.length_distribution_matrix[alignment->core.l_qseq][DIST_FILTERED_LEN_MAPPED_QUAL_ONTARGET]++;
                    }
                }
            }
            else
            {
                // Unpaired read
                if (!passes_mapq)
                {
                    continue;
                }
                
                summary_stats.filtered_len_mapped_qual_alignments++;
                summary_stats.length_distribution_matrix[alignment->core.l_qseq][DIST_FILTERED_LEN_MAPPED_QUAL]++;

                if (on_target)
                {
                    if (is_target_file)
                    {
                        summary_stats.filtered_len_mapped_qual_target_alignments++;
                        summary_stats.length_distribution_matrix[alignment->core.l_qseq][DIST_FILTERED_LEN_MAPPED_QUAL_ONTARGET]++;
                    }
                    
                    bam1_t* new_alignment = bam_init1();
                    bam_copy1(new_alignment, alignment);
                    alignment_buffer.push_back(new_alignment);
                }
            }
        }

        if (paired)
        {
            rescue_failed_mates(alignment_buffer, pair_tracker);
        }

        // Write any remaining buffered alignments to the output BAM file
        if (remove_dups)
        {
            remove_duplicates(alignment_buffer, deduped_alignment_buffer, false, true, &dup_stats);
            if (write_alignment_buffer_to_bam(&output_bam_config, deduped_alignment_buffer) != 0)
            {
                std::cerr << "Failed to write alignment buffer to " << output_bam_config.bam_file_location << std::endl;
                bam_destructor(&input_bam_config);
                bam_destructor(&output_bam_config);
                return 1;
            }
        }
        else
        {
            if (write_alignment_buffer_to_bam(&output_bam_config, alignment_buffer) != 0)
            {
                std::cerr << "Failed to write alignment buffer to " << output_bam_config.bam_file_location << std::endl;
                bam_destructor(&input_bam_config);
                bam_destructor(&output_bam_config);
                return 1;
            }
        }

        // Write summary stats
        write_length_dist_file(len_dist_file_name, is_target_file, summary_stats, dup_stats, remove_dups);

        append_summary_stats_to_file(input_bam_name, summary_stats_file_name, summary_stats, dup_stats, is_target_file, remove_dups);
        append_summary_stats_to_file(input_bam_name, summary_stats_file_name_2, summary_stats, dup_stats, is_target_file, remove_dups);

        // Clean up
        bam_destructor(&input_bam_config);
        bam_destructor(&output_bam_config);

        // Generate BAM index
        if (sam_index_build(output_file_location.c_str(), 0) < 0)
        {
            std::cerr << "Failed to generate BAM index" << std::endl;
            return 1;
        }
    }

    return 0;
}