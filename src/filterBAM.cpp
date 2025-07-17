#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <random>
#include <thread>
#include <string>
#include <cstring>
#include <sstream>

#include "utils.h"
#include "bam_processing.h"

#define VERSION_NUMBER 0.3
#define ALN_BUFFER_SIZE 1000000

void print_help()
{
    std::cout << "Usage: filterBAM [options] <input.bam>" << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "  -out_folder <output_path> specify the output folder, otherwise files will be written to the current folder" << std::endl;
    std::cout << "  -suffix                   string to be added to output file name before '.bam'" << std::endl;
    std::cout << "  -p3                       require 3' C->T differences at the indicated positions (e.g. 0,-1,-2)" << std::endl;
    std::cout << "  -p5                       require 5' C->T differences at the indicated positions (e.g. 0,1,2)" << std::endl;
    std::cout << "  -h                        display this help message" << std::endl;
}

int main(int argc, char *argv[])
{
    std::vector<std::string> bam_files;
    std::string out_folder = "";
    std::string suffix = "";
    std::vector<int16_t> p3;
    std::vector<int16_t> p5;

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
        if (std::string(argv[i]).find(".bam") != std::string::npos)
        {
            bam_files.push_back(argv[i]);
        }
        if (std::string(argv[i]) == "-suffix")
        {
            suffix = argv[i + 1];
        }
        if (std::string(argv[i]) == "-p3")
        {
            std::string p3_arg = argv[i + 1];
            size_t pos = 0;
            while ((pos = p3_arg.find(',')) != std::string::npos)
            {
                p3.push_back(std::stoi(p3_arg.substr(0, pos)));
                p3_arg.erase(0, pos + 1);
            }
            if (!p3_arg.empty())
            {
                p3.push_back(std::stoi(p3_arg));
            }
        }
        if (std::string(argv[i]) == "-p5")
        {
            std::string p5_arg = argv[i + 1];
            size_t pos = 0;
            while ((pos = p5_arg.find(',')) != std::string::npos)
            {
                p5.push_back(std::stoi(p5_arg.substr(0, pos)));
                p5_arg.erase(0, pos + 1);
            }
            if (!p5_arg.empty())
            {
                p5.push_back(std::stoi(p5_arg));
            }
        }
        if ((std::string(argv[i]) == "-h") || (argc == 1))
        {
            print_help();
            return 0;
        }
    }

    std::cout << "filterBAM.cpp v" << VERSION_NUMBER << " - " << get_date_time() << std::endl;

    for (const auto &bam_file : bam_files)
    {
        std::cout << "\nProcessing " << bam_file << "\n";

        std::string input_bam_name = extract_file_name(bam_file);
        std::string output_file_location = out_folder + input_bam_name + "." + suffix + ".bam";

        bam_file_config_t input_bam_config = {};
        bam_constructor(bam_file, &input_bam_config, "r");

        bam_file_config_t output_bam_config = {};
        bam_constructor(output_file_location, &output_bam_config, "wb");
        output_bam_config.header = bam_hdr_dup(input_bam_config.header);

        std::string version_str = std::to_string(VERSION_NUMBER);
        version_str.erase(version_str.find_last_not_of('0') + 1, std::string::npos);
        append_line_to_bam_header(output_bam_config.header, "@PG\tID:filterBAM\tVN:" + std::string(version_str) + "\tCL:" + command_line_str + "\n");

        // Write the header
        if (sam_hdr_write(output_bam_config.bam_file, output_bam_config.header) < 0)
        {
            std::cerr << "Failed to write BAM header to output file: " << output_file_location << std::endl;
            bam_destructor(&input_bam_config);
            bam_destructor(&output_bam_config);
            return 1;
        }

        int64_t counter = 0;
        int64_t deaminated_reads = 0;
        bam1_t *alignment = bam_init1();
        std::vector<bam1_t *> alignment_buffer;

        while (sam_read1(input_bam_config.bam_file, input_bam_config.header, alignment) >= 0)
        {
            counter++;
            if ((counter % 100000) == 0)
            {
                std::cout << "\r" << counter << " alignments processed";
                std::cout.flush();
            }

            if (alignment->core.flag & BAM_FUNMAP)
            {
                continue;
            }

            int nm = bam_aux2i(bam_aux_get(alignment, "NM"));
            if (nm == 0)
            {
                continue;
            }
            bool is_reverse = alignment->core.flag & BAM_FREVERSE;

            std::string reconstructed_alignment = reconstruct_alignment(alignment);
            std::string reconstructed_reference = reconstruct_reference(alignment, reconstructed_alignment);

            // Check for deamination at ends
            std::string alignment_ends;
            std::string reference_ends;

            for (auto p5_val : p5)
            {
                alignment_ends += reconstructed_alignment[p5_val];
                reference_ends += reconstructed_reference[p5_val];
            }
            for (auto p3_val : p3)
            {
                alignment_ends += reconstructed_alignment[reconstructed_alignment.size() + p3_val - 1];
                reference_ends += reconstructed_reference[reconstructed_reference.size() + p3_val - 1];
            }

            for (size_t i = 0; i < alignment_ends.size(); i++)
            {
                if ((reference_ends[i] == 'C' && alignment_ends[i] == 'T' && !is_reverse) ||
                    (reference_ends[i] == 'G' && alignment_ends[i] == 'A' && is_reverse))
                {
                    deaminated_reads++;

                    bam1_t *new_alignment = bam_init1();
                    bam_copy1(new_alignment, alignment);
                    alignment_buffer.push_back(new_alignment);
                    break;
                }
            }

            // Write the alignment buffer to the output BAM file
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

        bam_destructor(&input_bam_config);
        bam_destructor(&output_bam_config);

        // Generate BAM index
        if (sam_index_build(output_file_location.c_str(), 0) < 0)
        {
            std::cerr << "Failed to generate BAM index" << std::endl;
            return 1;
        }

        std::cout << "\n" << deaminated_reads << " deaminated reads found" << std::endl;
    }

    return 0;
}
