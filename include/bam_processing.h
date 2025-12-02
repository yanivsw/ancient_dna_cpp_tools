#ifndef BAM_PROCESSING_H
#define BAM_PROCESSING_H

#include <string>
#include <vector>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/kstring.h>
#include <iostream>
#include <cstring>

#include "utils.h"

struct bam_file_config_t
{
    std::string bam_file_location;
    samFile* bam_file;
    bam_hdr_t* header;
    hts_idx_t* index;
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

    // Read header and index if in read mode
    if (mode[0] == 'r')
    {
        if (!(bam_config->header = sam_hdr_read(bam_config->bam_file)))
        {
            std::cerr << "Failed to open BAM file header" << std::endl;
            return false;
        }
        if (!(bam_config->index = sam_index_load(bam_config->bam_file, bam_config->bam_file_location.c_str())))
        {
            // std::cerr << "Failed to open index for BAM file: " << bam_config->bam_file_location << std::endl;
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

bool bam_index_exists(const std::string& bam_file)
{
    std::string index_file = bam_file + ".bai";
    std::ifstream file(index_file);
    return file.good();
}

int create_bam_index(const std::string& bam_file)
{
    std::cout << "Creating BAM index for: " << bam_file << std::endl;

    if (sam_index_build(bam_file.c_str(), 0) != 0)
    {
        std::cerr << "Failed to create BAM index for: " << bam_file << std::endl;
        return 1;
    }

    // std::cout << "BAM index created successfully." << std::endl;
    return 0;
}

std::pair<int, int> check_bam(
    std::string bam_file_location,
    std::string chr,
    int num_threads)
{
    bam_file_config_t bam_config = {};
    if (!bam_constructor(bam_file_location, &bam_config))
    {
        return { -1, -1 };
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

    return { tid, chunk_size };
}

void get_chromosomes_from_bam(
    std::vector<std::string>& chromosomes,
    const std::string& bam_file_location)
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
    std::vector<bam1_t*>& alignment_buffer)
{
    for (bam1_t* aln : alignment_buffer)
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
    bam_hdr_t* header,
    const std::string& header_line)
{
    if (header == nullptr)
    {
        std::cerr << "Header is null" << std::endl;
        return;
    }

    // Create a string from the current header text
    std::string current_header_text(header->text, header->l_text);

    // Make sure the header line ends with a newline
    std::string new_header_line = header_line;
    if (!new_header_line.empty() && new_header_line.back() != '\n')
    {
        new_header_line += '\n';
    }

    // Find the position of the first @PG line
    size_t pg_pos = std::string::npos;
    if (current_header_text.compare(0, 3, "@PG") == 0)
    {
        pg_pos = 0;
    }
    else
    {
        pg_pos = current_header_text.find("\n@PG");
        if (pg_pos != std::string::npos)
        {
            // Skip the newline
            pg_pos += 1;
        }
    }

    // If we found a @PG line and the new line is also a @PG line
    if (pg_pos != std::string::npos && new_header_line.compare(0, 3, "@PG") == 0)
    {
        // Get the ID of the previous @PG line
        std::string pg_line;
        size_t line_end = current_header_text.find('\n', pg_pos);
        if (line_end != std::string::npos)
        {
            pg_line = current_header_text.substr(pg_pos, line_end - pg_pos);
        }
        else
        {
            pg_line = current_header_text.substr(pg_pos);
        }

        // Find the ID field
        size_t id_pos = pg_line.find("\tID:");
        if (id_pos != std::string::npos)
        {
            // Skip "\tID:"
            id_pos += 4;
            // Find the end of the ID value (next tab or end of line)
            size_t id_end = pg_line.find('\t', id_pos);
            if (id_end == std::string::npos)
            {
                id_end = pg_line.length();
            }

            std::string prev_id = pg_line.substr(id_pos, id_end - id_pos);

            // Add PP tag to the new @PG line
            if (!new_header_line.empty() && new_header_line.back() == '\n')
            {
                new_header_line.pop_back();
            }
            new_header_line += "\tPP:" + prev_id + "\n";
        }
    }

    // Create the new header text
    std::string new_header_text;
    if (pg_pos != std::string::npos)
    {
        // Insert before the @PG line
        new_header_text = current_header_text.substr(0, pg_pos) +
            new_header_line +
            current_header_text.substr(pg_pos);
    }
    else
    {
        // No @PG line found, append to the end
        new_header_text = current_header_text + new_header_line;
    }

    // Update the header
    free(header->text);
    header->text = strdup(new_header_text.c_str());
    header->l_text = new_header_text.size();
}

int get_alignment_length(const bam1_t* alignment)
{
    int length = 0;
    uint32_t* cigar = bam_get_cigar(alignment);
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

std::string reconstruct_alignment(const bam1_t* alignment)
{
    std::string reconstructed_sequence;
    uint8_t* seq = bam_get_seq(alignment);
    uint32_t* cigar = bam_get_cigar(alignment);
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
    const bam1_t* alignment,
    std::string reconstructed_alignment)
{
    std::string reconstructed_reference;

    // Get the MD tag
    uint8_t* md_ptr = bam_aux_get(alignment, "MD");
    if (md_ptr == nullptr)
    {
        std::cerr << "MD tag not found" << std::endl;
        return "";
    }

    const char* md_str = bam_aux2Z(md_ptr);
    std::string md_string(md_str);
    std::vector<std::string> md_segments = split_by_numeric(md_string);

    // Parse the MD tag and reconstruct the reference sequence
    size_t seq_index = 0;
    for (const std::string& segment : md_segments)
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

std::vector<uint8_t> get_base_quality(const bam1_t* alignment)
{
    std::vector<uint8_t> quality;
    uint8_t* qual = bam_get_qual(alignment);
    for (int i = 0; i < alignment->core.l_qseq; ++i)
    {
        quality.push_back(qual[i]);
    }
    return quality;
}

uint8_t get_minimum_base_quality(const bam1_t* alignment)
{
    std::vector<uint8_t> quality = get_base_quality(alignment);

    return *std::min_element(quality.begin(), quality.end());
}

// aux_type hack
static inline char aux_type(uint8_t* s)
{
    return *(s);
}

void update_xp_tag(
    bam1_t* alignment,
    uint32_t group_size)
{
    // Update XP tag in consensus alignment
    uint8_t* old_tag = bam_aux_get(alignment, "XP");
    if (old_tag)
    {
        bam_aux_del(alignment, old_tag);
    }

    bam_aux_append(alignment, "XP", 'i', 4, (uint8_t*) &group_size);
}

bool parse_coordinates_from_read_name(
    const char* qname,
    int32_t& x_coord,
    int32_t& y_coord)
{
    std::string name(qname);
    std::vector<std::string> parts;
    std::stringstream ss(name);
    std::string part;

    while (std::getline(ss, part, ':'))
    {
        parts.push_back(part);
    }

    // Format: MACHINE:RUN:FLOWCELL:LANE:TILE:X:Y or MACHINE_LANE:TILE:X:Y
    if (parts.size() >= 4)
    {
        try
        {
            // For newer format with 7 parts: MACHINE:RUN:FLOWCELL:LANE:TILE:X:Y
            if (parts.size() >= 7)
            {
                x_coord = std::stoi(parts[5]);

                // Handle Y coordinate (may include read number after #)
                std::string y_part = parts[6];
                size_t hash_pos = y_part.find('#');
                if (hash_pos != std::string::npos)
                {
                    y_part = y_part.substr(0, hash_pos);
                }
                y_coord = std::stoi(y_part);
                return true;
            }
            // For older format with 4+ parts: MACHINE_LANE:TILE:X:Y
            else
            {
                x_coord = std::stoi(parts[parts.size() - 2]);

                // Handle Y coordinate (may include read number after #)
                std::string y_part = parts[parts.size() - 1];
                size_t hash_pos = y_part.find('#');
                if (hash_pos != std::string::npos)
                {
                    y_part = y_part.substr(0, hash_pos);
                }
                y_coord = std::stoi(y_part);
                return true;
            }
        }
        catch (const std::exception& e)
        {
            // Failed to parse coordinates
            return false;
        }
    }
    return false;
}

void store_xy_coordinates(
    bam1_t* alignment,
    const std::vector<bam1_t*>& alignment_group)
{
    // Vectors to store all X and Y coordinates from the alignment group
    std::vector<int32_t> x_coords;
    std::vector<int32_t> y_coords;

    // Extract coordinates from all reads in the group
    for (auto& aln : alignment_group)
    {
        int32_t x = -1, y = -1;
        if (parse_coordinates_from_read_name(bam_get_qname(aln), x, y))
        {
            // Check if this coordinate is not already in our vectors
            if (std::find(x_coords.begin(), x_coords.end(), x) == x_coords.end())
            {
                x_coords.push_back(x);
            }
            if (std::find(y_coords.begin(), y_coords.end(), y) == y_coords.end())
            {
                y_coords.push_back(y);
            }
        }
    }

    // Skip if no coordinates found
    if (x_coords.empty())
    {
        return;
    }

    // Remove any existing X/Y tags
    uint8_t* old_x = bam_aux_get(alignment, "XC");
    if (old_x)
    {
        bam_aux_del(alignment, old_x);
    }

    uint8_t* old_y = bam_aux_get(alignment, "YC");
    if (old_y)
    {
        bam_aux_del(alignment, old_y);
    }

    // Method 1: Using HTSlib's array helpers (preferred)
    kstring_t str = { 0, 0, NULL };
    int32_t i;
    // For X coordinates
    if (!x_coords.empty())
    {
        str.l = 0;
        ksprintf(&str, "B:i,%d", (int32_t) x_coords.size());
        for (i = 0; i < static_cast<int32_t>(x_coords.size()); ++i)
        {
            ksprintf(&str, ",%d", x_coords[i]);
        }
        bam_aux_append(alignment, "XC", 'Z', str.l + 1, (uint8_t*) str.s);
    }

    // For Y coordinates
    if (!y_coords.empty())
    {
        str.l = 0;
        ksprintf(&str, "B:i,%d", (int32_t) y_coords.size());
        for (i = 0; i < static_cast<int32_t>(y_coords.size()); ++i)
        {
            ksprintf(&str, ",%d", y_coords[i]);
        }
        bam_aux_append(alignment, "YC", 'Z', str.l + 1, (uint8_t*) str.s);
    }

    if (str.s)
    {
        free(str.s);
    }
}

void update_read_name(
    bam1_t* aln,
    const std::string& new_name,
    const std::string& old_qname)
{
    // Validate inputs
    if (!aln || !aln->data || new_name.empty() || new_name.length() >= 254)
    {
        return;
    }

    size_t old_qname_len = old_qname.length() + 1; // Include null terminator
    size_t new_qname_len = new_name.length() + 1; // Include null terminator

    // Calculate padded lengths to multiples of 4
    size_t old_padded_len = (old_qname_len + 3) & ~3; // Round up to multiple of 4
    int extranul = (new_qname_len & 3) ? (4 - (new_qname_len & 3)) : 0;
    size_t new_padded_len = new_qname_len + extranul;

    // Calculate size difference and new data length
    int64_t length_diff = new_padded_len - old_padded_len;
    size_t new_data_len = aln->l_data + length_diff;

    // Reallocate if needed
    if (new_data_len > aln->m_data)
    {
        size_t new_m_data = new_data_len;
        kroundup32(new_m_data);
        uint8_t* new_data = (uint8_t*) realloc(aln->data, new_m_data);
        if (!new_data)
        {
            fprintf(stderr, "Failed to reallocate memory for read name\n");
            return;
        }
        aln->data = new_data;
        aln->m_data = new_m_data;
    }

    // Move existing data if name length changed
    if (length_diff != 0)
    {
        memmove(aln->data + new_padded_len,
            aln->data + old_padded_len,
            aln->l_data - old_padded_len);
    }

    // Copy new name
    memcpy(aln->data, new_name.c_str(), new_name.length());

    // Null terminate the new name
    aln->data[new_name.length()] = '\0';

    // Add padding nulls
    for (int i = 0; i < extranul; i++)
    {
        aln->data[new_qname_len + i] = '\0';
    }

    // Update BAM record fields
    aln->l_data = new_data_len;
    aln->core.l_qname = new_padded_len;
    aln->core.l_extranul = extranul;
}

#endif // BAM_PROCESSING_H