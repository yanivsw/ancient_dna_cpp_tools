#ifndef TYPES_H
#define TYPES_H

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

// Filters for analyzeBAM
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

// Summary statistics for analyzeBAM
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

// Duplication statistics for analyzeBAM and dedupBAM
struct duplication_stats_struct
{
    uint64_t uniq;
    uint64_t total;
    uint64_t single;
    uint64_t qc_failed;
    std::vector<uint64_t> length_distribution_uniq;
    std::vector<uint64_t> length_distribution_single;
    std::vector<uint64_t> length_distribution_total;

    duplication_stats_struct(size_t length)
    {
        uniq = 0;
        total = 0;
        single = 0;
        qc_failed = 0;

        initialize_vector(length);
    }

    void initialize_vector(size_t length)
    {
        length_distribution_uniq = std::vector<uint64_t>(length, 0);
        length_distribution_single = std::vector<uint64_t>(length, 0);
        length_distribution_total = std::vector<uint64_t>(length, 0);
    }
};

// Key structure for grouping alignments for duplicate removal
struct AlignmentGroupKey
{
    int32_t position;
    uint8_t orientation;
    uint32_t read_group_hash;
    bool is_paired;
    uint32_t length;

    bool operator==(const AlignmentGroupKey& other) const
    {
        return position == other.position &&
            orientation == other.orientation &&
            read_group_hash == other.read_group_hash &&
            is_paired == other.is_paired &&
            length == other.length;
    }
};

// Fast hash function from https://stackoverflow.com/questions/17016175/c-unordered-map-with-a-custom-class-type-as-key
struct AlignmentGroupKeyHash
{
    std::size_t operator()(const AlignmentGroupKey& key) const
    {
        // FNV-1a hash combination
        std::size_t h = 0x811c9dc5;
        h ^= key.position;        h *= 0x01000193;
        h ^= key.orientation;     h *= 0x01000193;
        h ^= key.read_group_hash; h *= 0x01000193;
        h ^= key.is_paired;       h *= 0x01000193;
        h ^= key.length;          h *= 0x01000193;
        return h;
    }
};

#endif // TYPES_H