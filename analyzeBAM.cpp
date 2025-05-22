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
#include "utils.h"

#define VERSION_NUMBER 0.6

#define ALN_BUFFER_SIZE 1000000
#define READ_LEN_DIST_SIZE 512

enum filter_stage_t
{
    DIST_RAW = 0,
    DIST_MERGED,
    DIST_FILTERED,
    DIST_FILTERED_LEN,
    DIST_FILTERED_LEN_MAPPED,
    DIST_FILTERED_LEN_MAPPED_QUAL,
    DIST_FILTERED_LEN_MAPPED_QUAL_ONTARGET
};

struct summary_stats_struct
{
    uint64_t raw_alignments;
    uint64_t merged_alignments;
    uint64_t filtered_alignments;
    uint64_t filtered_len_alignments;
    uint64_t filtered_len_mapped_alignments;
    uint64_t filtered_len_mapped_qual_alignments;
    uint64_t filtered_len_mapped_qual_target_alignments;

    std::vector<std::vector<int>> length_distribution_matrix;

    summary_stats_struct(size_t length)
    {
        raw_alignments = 0;
        merged_alignments = 0;
        filtered_alignments = 0;
        filtered_len_alignments = 0;
        filtered_len_mapped_alignments = 0;
        filtered_len_mapped_qual_alignments = 0;
        filtered_len_mapped_qual_target_alignments = 0;

        initialize_vector(length);
    }

    void initialize_vector(size_t length)
    {
        length_distribution_matrix = std::vector<std::vector<int>>(length, std::vector<int>(DIST_FILTERED_LEN_MAPPED_QUAL_ONTARGET + 1, 0));
    }
};

struct duplication_stats_struct
{
    uint64_t uniq;
    uint64_t total;
    uint64_t single;
    std::vector<uint64_t> length_distribution_uniq;
    std::vector<uint64_t> length_distribution_single;
    std::vector<uint64_t> length_distribution_total;

    duplication_stats_struct(size_t length)
    {
        uniq = 0;
        total = 0;
        single = 0;

        initialize_vector(length);
    }

    void initialize_vector(size_t length)
    {
        length_distribution_uniq = std::vector<uint64_t>(length, 0);
        length_distribution_single = std::vector<uint64_t>(length, 0);
        length_distribution_total = std::vector<uint64_t>(length, 0);
    }
};

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
    "GL000200.1", "GL000193.1", "GL000194.1", "GL000225.1", "GL000192.1"};

void write_length_dist_file(
    const std::string &len_dist_file_name,
    bool is_target_file,
    const summary_stats_struct &summary_stats)
{
    std::ofstream len_dist_file(len_dist_file_name, std::ios::app);
    for (size_t j = 0; j < summary_stats.length_distribution_matrix.size(); j++)
    {
        const auto &row = summary_stats.length_distribution_matrix[j];
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
                len_dist_file << "\n";
            }
        }
    }
    len_dist_file.close();
}

void append_summary_stats_to_file(
    const std::string &input_bam_name,
    const std::string &summary_stats_file_name,
    const summary_stats_struct &summary_stats,
    bool is_target_file)
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
    summary_stats_file << "\n";
    summary_stats_file.close();
}

int check_dups(
    std::string bam_file,
    duplication_stats_struct &dup_stats)
{
    bam_file_config_t bam_config = {};
    bam_constructor(bam_file, &bam_config, "r");

    bam1_t *alignment = bam_init1();
    while (sam_read1(bam_config.bam_file, bam_config.header, alignment) >= 0)
    {
        uint8_t *xp_ptr = bam_aux_get(alignment, "XP");
        if (xp_ptr)
        {
            auto xp_flag = bam_aux2i(xp_ptr);
            dup_stats.total += xp_flag;
            dup_stats.length_distribution_total[alignment->core.l_qseq] += xp_flag;
        }
        else
        {
            dup_stats.single++;
            dup_stats.total++;

            dup_stats.length_distribution_total[alignment->core.l_qseq]++;
            dup_stats.length_distribution_single[alignment->core.l_qseq]++;
        }

        dup_stats.length_distribution_uniq[alignment->core.l_qseq]++;
        dup_stats.uniq++;
    }

    bam_destructor(&bam_config);
    bam_destroy1(alignment);
    return 0;
}

void update_summary_stats_for_dups(
    const std::string &in_name,
    const std::string &out_name,
    const std::string &input_bam_name,
    const duplication_stats_struct &dup_stats)
{
    std::ifstream in(in_name);
    std::ofstream out(out_name, std::ofstream::trunc);
    std::string line;
    while (std::getline(in, line))
    {
        if ((line.find("#file") != std::string::npos) && (line.find("unique") == std::string::npos))
        {
            out << line << "\tunique\taverage_dups\tsingletons\n";
        }
        else if ((line.find(input_bam_name) != std::string::npos) && (count_tabs(line) <= 8))
        {
            out << line << "\t" << dup_stats.uniq << "\t"
                << (static_cast<double>(dup_stats.total) / dup_stats.uniq) << "\t"
                << dup_stats.single << "\n";
        }
        else
        {
            out << line << "\n";
        }
    }
}

void handle_deduped_bam(
    int32_t min_len,
    int32_t min_map_qual,
    const std::string& output_file_location,
    const std::string& input_bam_name,
    const std::string& summary_stats_file_name,
    const std::string& len_dist_file_name)
{
    std::string bam_file_deduped = replace_substring(output_file_location, 
        ".L" + std::to_string(min_len) + "MQ" + std::to_string(min_map_qual), 
        ".uniq.L" + std::to_string(min_len) + "MQ" + std::to_string(min_map_qual));

    duplication_stats_struct dup_stats(READ_LEN_DIST_SIZE);
    check_dups(bam_file_deduped, dup_stats);

    // Update summary stats files
    update_summary_stats_for_dups(summary_stats_file_name, "temp.txt", input_bam_name, dup_stats);
    std::rename("temp.txt", summary_stats_file_name.c_str());

    std::string summary_stats_file_name_2 = replace_substring(summary_stats_file_name, "summary_stats", "summary_stats." + input_bam_name);
    update_summary_stats_for_dups(summary_stats_file_name_2, "temp.txt", input_bam_name, dup_stats);
    std::rename("temp.txt", summary_stats_file_name_2.c_str());

    // Update length distribution file
    std::ifstream len_dist_file_in(len_dist_file_name);
    std::ofstream len_dist_file_out("temp.txt", std::ofstream::trunc);
    std::string line;
    int idx = 0;
    while (std::getline(len_dist_file_in, line))
    {
        if ((line.find("Length") != std::string::npos) &&
            (line.find("unique") == std::string::npos))
        {
            len_dist_file_out << line << "\tunique\taverage_dups\tsingletons\n";
        }
        else
        {
            if (dup_stats.length_distribution_uniq[idx] > 0)
            {
                len_dist_file_out << line << "\t" << dup_stats.length_distribution_uniq[idx]
                << "\t" << (static_cast<double>(dup_stats.length_distribution_total[idx]) / dup_stats.length_distribution_uniq[idx])
                << "\t" << dup_stats.length_distribution_single[idx] << "\n";
            }
            else
            {
                len_dist_file_out << line << "\tNA\tNA\tNA\n";
            }
            idx++;
        }
    }
    len_dist_file_out.close();
    len_dist_file_in.close();
    std::rename("temp.txt", len_dist_file_name.c_str());
}

void print_help()
{
    std::cout << "Usage: analyzeBAM [options] <BAM files>\n"
              << "Options:\n"
              << "  -out_folder <path>       Specify the output folder, otherwise files will be written to the current folder\n"
              << "  -min_len <length>        Specify the minimum alignment length [default 35]\n"
              << "  -min_map_qual <quality>  Specify the minimum mapping quality [default 0]\n"
              << "  -targetfile <path>       Specify the target file path. ***NB*** This expects a bed file (i.e. tab separated, 0-based)\n"
              << "  -paired                  Do not disregard paired reads\n"
              << "  -check_dups              Generate duplication stats for post bam-rmdup BAM files. All other options should be the same as the ones used to generate the pre bam-rmdup summary statistics and BAM files. \n"
              << "  -count_f                 Ignore filter (QC failed) flag\n"
              << "  -help                    Display this help message\n";
}

int main(int argc, char *argv[])
{
    int32_t min_len = 35;
    int32_t min_map_qual = 0;
    std::vector<std::string> bam_files;
    std::string target_file_path = "";
    std::string target_file_name = "";
    bool is_target_file = false;
    bool paired = false;
    bool check_dups = false;
    bool count_f = false;
    std::string min_len_str;
    std::string min_map_qual_str;
    bed_file_t bed_hash;
    std::string out_folder = "";

    std::cout << "analyzeBAM.cpp v" << VERSION_NUMBER << " - " << get_date_time() << std::endl;

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
        if (std::string(argv[i]) == "-check_dups")
        {
            check_dups = true;
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

    std::string summary_stats_file_name = "summary_stats" +
                                          (is_target_file ? ".ontarget_" + target_file_name : "") +
                                          ".L" + min_len_str + "MQ" + min_map_qual_str + ".txt";

    if (!check_dups)
    {
        std::ofstream summary_stats_file(summary_stats_file_name);
        summary_stats_file << "#file\traw\tmerged\tfilter_passed\tL" << min_len_str
                           << "\tmappedL" << min_len_str
                           << "\tmappedL" << min_len_str << "MQ" << min_map_qual_str
                           << "\t%mappedL" << min_len_str << "MQ" << min_map_qual_str;
        if (is_target_file)
        {
            summary_stats_file << "\ttargetL" << min_len_str << "MQ" << min_map_qual_str;
        }
        summary_stats_file << "\n";
    }

    // Read target file to an unordered_map if it exists
    if ((is_target_file) && (check_dups == false))
    {
        read_bed_file(target_file_path, bed_hash);
    }

    // Loop over bam files
    for (const auto &bam_file : bam_files)
    {
        // Extract the file name without the extension and path
        std::string input_bam_name;
        input_bam_name = extract_file_name(bam_file);

        // Output file name
        std::string suffix = (std::string(is_target_file ? ".ontarget" : "") + ".L" + min_len_str + "MQ" + min_map_qual_str);
        std::string output_file_location = out_folder + input_bam_name + suffix + ".bam";

        std::string target_suffix = (is_target_file ? ".ontarget_" + target_file_name : "") + ".L" + min_len_str + "MQ" + min_map_qual_str;
        std::string summary_stats_file_name_2 = "summary_stats." + input_bam_name + target_suffix + ".txt";
        std::string len_dist_file_name = "read_length_distribution." + input_bam_name + target_suffix + ".tsv";

        if (!check_dups)
        {
            std::ofstream summary_stats_file_2(summary_stats_file_name_2);
            summary_stats_file_2 << "#file\traw\tmerged\tfilter_passed\tL" << min_len_str
                                 << "\tmappedL" << min_len_str
                                 << "\tmappedL" << min_len_str << "MQ" << min_map_qual_str
                                 << "\t%mappedL" << min_len_str << "MQ" << min_map_qual_str
                                 << "\ttargetL" << min_len_str << "MQ" << min_map_qual_str << "\n";

            std::ofstream len_dist_file(len_dist_file_name);
            len_dist_file << "Length\traw\tmerged\tfilter_passed\tL" << min_len_str
                          << "\tmappedL" << min_len_str
                          << "\tmappedL" << min_len_str << "MQ" << min_map_qual_str
                          << "\t%mappedL" << min_len_str << "MQ" << min_map_qual_str
                          << "\ttargetL" << min_len_str << "MQ" << min_map_qual_str << "\n";
        }
        else
        {
            handle_deduped_bam(min_len, min_map_qual, output_file_location, input_bam_name, summary_stats_file_name, len_dist_file_name);
            continue;
        }

        summary_stats_struct summary_stats(READ_LEN_DIST_SIZE);

        uint64_t counter = 0;
        uint64_t qc_fail_reads = 0;
        uint64_t unmapped_reads = 0;
        uint64_t off_genome_reads = 0;

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

        uint64_t bed_hash_ptr = 0;
        int target_end = 0;
        bool on_target = false;
        std::string last_seen_chr = "";

        // bam flags: "qdfs21RrUuPp"
        // Read each alignment
        std::vector<bam1_t *> alignment_buffer;
        bam1_t *alignment = bam_init1();
        while (sam_read1(input_bam_config.bam_file, input_bam_config.header, alignment) >= 0)
        {
            on_target = false;
            counter++;
            if ((counter % 100000) == 0)
            {
                std::cout << "\r" << counter << " alignments processed";
                std::cout.flush();
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

            // Skip if alignment length <= min_len
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
            std::string chromosome = chromosomes[alignment->core.tid];
            if (off_genome_chromosomes.count(chromosome) > 0)
            {
                off_genome_reads++;
                continue;
            }
            summary_stats.filtered_len_mapped_alignments++;
            summary_stats.length_distribution_matrix[alignment->core.l_qseq][DIST_FILTERED_LEN_MAPPED]++;

            // Skip if the mapping quality <= min_map_qual
            if (alignment->core.qual < min_map_qual)
            {
                continue;
            }
            summary_stats.filtered_len_mapped_qual_alignments++;
            summary_stats.length_distribution_matrix[alignment->core.l_qseq][DIST_FILTERED_LEN_MAPPED_QUAL]++;

            // If we have a target file
            if (is_target_file)
            {
                // Check if we have not seen this chromosome before (we are assuming the BAM file is sorted by chromosome)
                if (last_seen_chr != chromosome)
                {
                    if (bed_hash[chromosome].size() > 0)
                    {
                        bed_hash_ptr = 0;
                        target_end = bed_hash[chromosome][0].end;
                        last_seen_chr = chromosome;
                    }
                    else if (chromosome == "MT" || chromosome == "chrM") // Handle mitochondria explicitly
                    {
                        // Skip processing for mitochondria if not in the target file
                        continue;
                    }
                }

                int aln_len = get_alignment_length(alignment);
                int aln_start = alignment->core.pos;
                int aln_end = aln_start + aln_len;

                // Increment the target region pointer until we find one that could overlap with the alignment
                while ((aln_start > target_end) && (bed_hash_ptr < bed_hash[chromosome].size()))
                {
                    bed_hash_ptr++;
                    target_end = bed_hash[chromosome][bed_hash_ptr].end;
                }

                // Check for overlap between the alignment and the target region
                for (size_t i = bed_hash_ptr; i < bed_hash[chromosome].size(); i++)
                {
                    int region_start = bed_hash[chromosome][i].start;
                    int region_end = bed_hash[chromosome][i].end;

                    // The alignment is on target
                    if ((aln_end > region_start) && (aln_start < region_end))
                    {
                        on_target = true;
                        break;
                    }

                    // The alignment doesn't overlap with any target regions
                    if (region_start > aln_end)
                    {
                        break;
                    }
                }

                if (on_target)
                {
                    summary_stats.filtered_len_mapped_qual_target_alignments++;
                    summary_stats.length_distribution_matrix[alignment->core.l_qseq][DIST_FILTERED_LEN_MAPPED_QUAL_ONTARGET]++;
                }
            }

            // Buffer the alignment before writing
            if ((is_target_file && on_target) || !is_target_file)
            {
                bam1_t *new_alignment = bam_init1();
                bam_copy1(new_alignment, alignment);
                alignment_buffer.push_back(new_alignment);
            }

            // Write the alignment buffer to the output BAM
            if (alignment_buffer.size() >= ALN_BUFFER_SIZE)
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

        // Write any remaining alignments in the buffer
        if (write_alignment_buffer_to_bam(&output_bam_config, alignment_buffer) != 0)
        {
            std::cerr << "Failed to write alignment buffer to " << output_bam_config.bam_file_location << std::endl;
            bam_destructor(&input_bam_config);
            bam_destructor(&output_bam_config);
            return 1;
        }

        // Write summary stats
        write_length_dist_file(len_dist_file_name, is_target_file, summary_stats);

        append_summary_stats_to_file(input_bam_name, summary_stats_file_name, summary_stats, is_target_file);
        append_summary_stats_to_file(input_bam_name, summary_stats_file_name_2, summary_stats, is_target_file);

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