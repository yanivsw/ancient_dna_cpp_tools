#ifndef BAM_PROCESSING_H
#define BAM_PROCESSING_H

#include <string>
#include <vector>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <iostream>
#include <cstring>

#include "utils.h"

struct bam_file_config_t
{
    std::string bam_file_location;
    samFile *bam_file;
    bam_hdr_t *header;
    hts_idx_t *index;
};

bool bam_constructor(
    std::string bam_file_location,
    bam_file_config_t* bam_config,
    const char* mode = "r")
{
    bam_config->bam_file_location = bam_file_location;
    bam_config->bam_file = nullptr;
    bam_config->header = nullptr;
    bam_config->index = nullptr;

    // Open BAM file with the given mode
    if (!(bam_config->bam_file = sam_open(bam_config->bam_file_location.c_str(), mode)))
    {
        std::cerr << "Failed to open BAM file: " << bam_config->bam_file_location << std::endl;
        return false;
    }

    // Read or write header as appropriate
    if (mode[0] == 'r')
    {
        if (!(bam_config->header = sam_hdr_read(bam_config->bam_file)))
        {
            std::cerr << "Failed to open BAM file header" << std::endl;
            return false;
        }
        // Only load index for reading
        if (!(bam_config->index = sam_index_load(bam_config->bam_file, bam_config->bam_file_location.c_str())))
        {
            std::cerr << "Failed to open index for BAM file: " << bam_config->bam_file_location << std::endl;
            return false;
        }
    }

    return true;
}

void bam_destructor(bam_file_config_t* bam_config)
{
    bam_hdr_destroy(bam_config->header);
    sam_close(bam_config->bam_file);

    if (bam_config->index)
    {
        hts_idx_destroy(bam_config->index);
    }
}

std::pair<int, int> check_bam(
    std::string bam_file_location,
    std::string chr,
    int num_threads)
{
    bam_file_config_t bam_config = {};
    if (!bam_constructor(bam_file_location, &bam_config))
    {
        return {-1, -1};
    }

    int tid;
    int chunk_size;

    if (chr[0] == '\0')
    {
        tid = -1;
        chunk_size = 1;
    }
    else
    {
        tid = bam_name2id(bam_config.header, chr.c_str());
        chunk_size = bam_config.header->target_len[tid] / num_threads;
    }

    bam_destructor(&bam_config);

    return {tid, chunk_size};
}

void get_chromosomes_from_bam(
    std::vector<std::string> &chromosomes,
    const std::string &bam_file_location)
{
    bam_file_config_t bam_config = {};
    if (!bam_constructor(bam_file_location, &bam_config))
    {
        std::cerr << "Failed to open BAM file: " << bam_file_location << std::endl;
        return;
    }

    // Get the chromosome names from the BAM header
    for (int i = 0; i < bam_config.header->n_targets; i++)
    {
        chromosomes.push_back(bam_config.header->target_name[i]);
    }
    bam_destructor(&bam_config);
}

int write_alignment_buffer_to_bam(
    bam_file_config_t* out_bam_config,
    std::vector<bam1_t *> &alignment_buffer)
{
    for (bam1_t *aln : alignment_buffer)
    {
        if (sam_write1(out_bam_config->bam_file, out_bam_config->header, aln) < 0)
        {
            return 1;
        }
        bam_destroy1(aln);
    }
    alignment_buffer.clear();
    return 0;
}

void append_line_to_bam_header(
    bam_hdr_t *header,
    const std::string &header_line)
{
    if (header == nullptr)
    {
        std::cerr << "Header is null" << std::endl;
        return;
    }

    // Create a new header text with the additional line
    std::string new_header_text(header->text, header->l_text);
    new_header_text += header_line;

    // Free the old header text and assign the new one
    free(header->text);
    header->text = strdup(new_header_text.c_str());
    header->l_text = new_header_text.size();
}

int get_alignment_length(const bam1_t *alignment)
{
    int length = 0;
    uint32_t *cigar = bam_get_cigar(alignment);
    for (uint32_t i = 0; i < alignment->core.n_cigar; i++)
    {
        int op = bam_cigar_op(cigar[i]);
        int oplen = bam_cigar_oplen(cigar[i]);
        if ((op == BAM_CMATCH) || (op == BAM_CEQUAL) || (op == BAM_CDIFF) || (op == BAM_CDEL))
        {
            length += oplen;
        }
    }
    return length;
}

std::string reconstruct_alignment(const bam1_t *alignment)
{
    std::string reconstructed_sequence;
    uint8_t *seq = bam_get_seq(alignment);
    uint32_t *cigar = bam_get_cigar(alignment);
    int32_t ref_len = 0;

    for (uint32_t i = 0; i < alignment->core.n_cigar; ++i)
    {
        uint32_t op_len = bam_cigar_oplen(cigar[i]);
        char op = bam_cigar_opchr(cigar[i]);

        switch (op)
        {
        case 'M': // Match or mismatch
        case '=': // Match
        case 'X': // Mismatch
            for (uint32_t j = 0; j < op_len; ++j)
            {
                reconstructed_sequence += seq_nt16_str[bam_seqi(seq, ref_len + j)];
            }
            ref_len += op_len;
            break;
        case 'I': // Insertion to the reference
        case 'S': // Soft clipping (clipped sequences present in SEQ)
            ref_len += op_len;
            break;
        case 'D': // Deletion from the reference
        case 'N': // Skipped region from the reference
            for (uint32_t j = 0; j < op_len; ++j)
            {
                reconstructed_sequence += '-';
            }
            break;
        case 'H': // Hard clipping (clipped sequences NOT present in SEQ)
            break;
        case 'P': // Padding (silent deletion from padded reference)
            break;
        default:
            std::cerr << "Unknown CIGAR operation: " << op << std::endl;
            break;
        }
    }

    return reconstructed_sequence;
}

std::string reconstruct_reference(
    const bam1_t *alignment,
    std::string reconstructed_alignment)
{
    std::string reconstructed_reference;

    // Get the MD tag
    uint8_t *md_ptr = bam_aux_get(alignment, "MD");
    if (md_ptr == nullptr)
    {
        std::cerr << "MD tag not found" << std::endl;
        return "";
    }

    const char *md_str = bam_aux2Z(md_ptr);
    std::string md_string(md_str);
    std::vector<std::string> md_segments = split_by_numeric(md_string);

    // Parse the MD tag and reconstruct the reference sequence
    size_t seq_index = 0;
    for (const std::string &segment : md_segments)
    {
        if (segment.empty() || segment == "0") continue;

        if (isdigit(segment[0]))
        {
            int match_length = std::stoi(segment);
            for (int i = 0; i < match_length; i++)
            {
                reconstructed_reference += reconstructed_alignment[seq_index];
                seq_index++;
            }
        }
        else if (segment[0] != '^')
        {
            // Handle single base pair changes in the reference
            for (size_t i = 0; i < segment.length(); i++)
            {
                reconstructed_reference += segment[i];
                seq_index++;
            }
        }
        else if (segment[0] == '^')
        {
            // Handle deletions in the reference
            for (size_t i = 1; i < segment.length(); i++)
            {
                reconstructed_reference += segment[i];
                seq_index++;
            }
        }
    }

    if (seq_index != reconstructed_alignment.size())
    {
        std::cerr << "Error: MD tag and reconstructed alignment lengths do not match" << std::endl;
    }

    return reconstructed_reference;
}

#endif // BAM_PROCESSING_H