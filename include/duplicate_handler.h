#ifndef DUPLICATE_HANDLER_H
#define DUPLICATE_HANDLER_H

#include <htslib/sam.h>
#include <htslib/kstring.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <unordered_map>
#include <sstream>
#include <random>
#include <map>
#include <set>
#include <string_view>

#include "bam_processing.h"
#include "types.h"

// Global RNG for random selection of a representative alignment among duplicates
std::mt19937 gen(std::random_device{}());

// Set to track which read names should be kept (first seen mates that passed deduplication)
std::unordered_set<std::string> retained_read_names;

void group_alignments(
    std::vector<bam1_t*>& alignment_buffer,
    std::unordered_map<AlignmentGroupKey, std::vector<bam1_t*>, AlignmentGroupKeyHash>& alignment_groups,
    bool ignore_length = false,
    bool ignore_read_groups = true)
{
    // First group by read name to find pairs
    std::unordered_map<std::string_view, std::vector<bam1_t*>> read_pairs;
    read_pairs.reserve(alignment_buffer.size());
    for (const auto& aln : alignment_buffer)
    {
        std::string_view qname(bam_get_qname(aln));
        read_pairs[qname].push_back(aln);
    }

    alignment_groups.reserve(read_pairs.size());

    std::unordered_map<std::string, uint32_t> rg_hash_cache;
    uint32_t next_rg_hash = 1;
    for (const auto& [qname, reads] : read_pairs)
    {
        // Get read group hash (only once per read name)
        uint32_t read_group_hash = 42069;
        if (!ignore_read_groups && !reads.empty())
        {
            uint8_t* rg_aux = bam_aux_get(reads[0], "RG");
            if (rg_aux)
            {
                std::string rg(bam_aux2Z(rg_aux));
                auto it = rg_hash_cache.find(rg);
                if (it == rg_hash_cache.end())
                {
                    read_group_hash = next_rg_hash++;
                    rg_hash_cache[rg] = read_group_hash;
                }
                else
                {
                    read_group_hash = it->second;
                }
            }
        }

        // Handle paired reads where both mates are present
        if (reads.size() == 2 && (reads[0]->core.flag & BAM_FPAIRED) && (reads[1]->core.flag & BAM_FPAIRED))
        {
            bam1_t* left_read = (reads[0]->core.pos <= reads[1]->core.pos) ? reads[0] : reads[1];
            bam1_t* right_read = (left_read == reads[0]) ? reads[1] : reads[0];

            bool is_reverse = left_read->core.flag & BAM_FREVERSE;
            uint32_t length = 0;
            if (!ignore_length)
            {
                length = right_read->core.pos;
            }

            AlignmentGroupKey key1 = {
                left_read->core.pos,
                is_reverse,
                read_group_hash,
                true,
                length
            };

            AlignmentGroupKey key2 = {
                right_read->core.pos,
                255, // special orientation value for retained second mates
                0,
                true,
                0
            };

            alignment_groups[key1].push_back(left_read);
            alignment_groups[key2].push_back(right_read);
        }
        // Handle paired reads where only one mate is present
        else if ((reads.size() == 1) && (reads[0]->core.flag & BAM_FPAIRED))
        {
            bam1_t* aln = reads[0];

            // If this read is the second mate, check if the first mate was retained
            if (retained_read_names.count(bam_get_qname(aln)))
            {
                AlignmentGroupKey key = {
                    aln->core.pos,
                    255, // special orientation value for retained second mates
                    0,
                    true,
                    0
                };
                alignment_groups[key].push_back(aln);
            }
            // If this is the first mate, process normally
            else
            {
                bool is_reverse = aln->core.flag & BAM_FREVERSE;
                uint32_t length = 0;
                if (!ignore_length)
                {
                    length = aln->core.mpos;
                }

                AlignmentGroupKey key = {
                    aln->core.pos,
                    is_reverse,
                    read_group_hash,
                    true,
                    length
                };

                alignment_groups[key].push_back(aln);
            }
        }
        // Handle merged reads
        else
        {
            bam1_t* aln = reads[0];

            uint32_t length = 0;
            if (!ignore_length)
            {
                length = get_alignment_length(aln);
            }

            bool is_reverse = aln->core.flag & BAM_FREVERSE;

            AlignmentGroupKey key = {
                aln->core.pos,
                is_reverse,
                read_group_hash,
                0,
                length
            };

            alignment_groups[key].push_back(aln);
        }
    }
}

void remove_duplicates(
    std::vector<bam1_t*>& alignment_buffer,
    std::vector<bam1_t*>& deduped_alignment_buffer,
    bool ignore_length = false,
    bool ignore_read_groups = true,
    duplication_stats_struct* dup_stats = nullptr)
{
    std::unordered_map<AlignmentGroupKey, std::vector<bam1_t*>, AlignmentGroupKeyHash> alignment_groups;
    group_alignments(alignment_buffer, alignment_groups, ignore_length, ignore_read_groups);

    // Convert to vector and sort by position
    std::vector<std::pair<AlignmentGroupKey, std::vector<bam1_t*>>> sorted_groups;
    sorted_groups.reserve(alignment_groups.size());
    for (auto& [key, group] : alignment_groups)
    {
        sorted_groups.push_back({ key, group });
    }
    std::sort(sorted_groups.begin(), sorted_groups.end(),
        [](const auto& a, const auto& b)
        {
            return a.first.position < b.first.position;
        });

    // Process each group
    for (auto& [key, alignment_group] : sorted_groups)
    {
        // If the mate position (key.length) is less than the position of the first seen mate
        // skip because we have not seen the mate yet which means it did not pass filters
        if (key.is_paired && (key.length < static_cast<uint32_t>(key.position)) && (key.orientation != 255))
        {
            continue;
        }

        // 255 orientation means this is the second mate of a pair - retain these alignments unconditionally
        if (key.orientation == 255)
        {
            for (auto* aln : alignment_group)
            {
                if (retained_read_names.count(bam_get_qname(aln)))
                {
                    bam1_t* deduped_alignment = bam_init1();
                    bam_copy1(deduped_alignment, aln);
                    update_xp_tag(deduped_alignment, alignment_group.size());
                    deduped_alignment_buffer.push_back(deduped_alignment);
                    break;
                }
            }
        }
        // If only one alignment in this group, retain it
        else if (alignment_group.size() == 1)
        {
            bam1_t* aln = alignment_group[0];
            bam1_t* deduped_alignment = bam_init1();
            bam_copy1(deduped_alignment, aln);

            update_xp_tag(deduped_alignment, 1);
            deduped_alignment_buffer.push_back(deduped_alignment);

            // If this is paired, record its name so we can look for the mate later
            if (aln->core.flag & BAM_FPAIRED)
            {
                std::string qname = bam_get_qname(aln);
                retained_read_names.insert(qname);
            }

            if (dup_stats)
            {
                dup_stats->single++;
                dup_stats->total++;
                dup_stats->uniq++;

                dup_stats->length_distribution_total[aln->core.l_qseq]++;
                dup_stats->length_distribution_single[aln->core.l_qseq]++;
                dup_stats->length_distribution_uniq[aln->core.l_qseq]++;
            }
        }
        // Multiple alignments in this group - need to select one using the highest minimum base quality
        else
        {
            std::vector<bam1_t*> qualified_alignments;
            uint8_t highest_quality = 0;
            bam1_t* highest_quality_alignment = nullptr;

            for (auto* aln : alignment_group)
            {
                uint8_t min_base_quality = get_minimum_base_quality(aln);
                if (min_base_quality > highest_quality)
                {
                    highest_quality = min_base_quality;
                    qualified_alignments.clear();
                    qualified_alignments.push_back(aln);
                    highest_quality_alignment = aln;
                }
                else if (min_base_quality == highest_quality)
                {
                    qualified_alignments.push_back(aln);
                }
            }

            // Select random alignment from qualified alignments
            bam1_t* selected_alignment = nullptr;
            if (qualified_alignments.size() > 1)
            {
                std::uniform_int_distribution<size_t> dist(0, qualified_alignments.size() - 1);
                selected_alignment = qualified_alignments[dist(gen)];
            }
            else
            {
                selected_alignment = highest_quality_alignment;
            }

            // Add selected alignment to output
            bam1_t* deduped_alignment = bam_init1();
            bam_copy1(deduped_alignment, selected_alignment);
            update_xp_tag(deduped_alignment, alignment_group.size());
            store_xy_coordinates(deduped_alignment, alignment_group);
            deduped_alignment_buffer.push_back(deduped_alignment);

            // If this is paired, record its name so we can look for the mate later
            if (selected_alignment->core.flag & BAM_FPAIRED)
            {
                std::string qname = bam_get_qname(deduped_alignment);
                retained_read_names.insert(qname);
            }

            if (dup_stats)
            {
                dup_stats->uniq++;
                dup_stats->length_distribution_uniq[selected_alignment->core.l_qseq]++;

                dup_stats->total += alignment_group.size();
                for (auto* aln : alignment_group)
                {
                    dup_stats->length_distribution_total[aln->core.l_qseq]++;
                }
            }
        }
    }

    // Delete all alignments in alignment_buffer to prevent memory leaks
    for (auto* aln : alignment_buffer)
    {
        bam_destroy1(aln);
    }

    alignment_buffer.clear();
    alignment_groups.clear();
}

#endif // DUPLICATE_HANDLER_H