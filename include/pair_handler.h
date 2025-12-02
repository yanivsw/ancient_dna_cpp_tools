#ifndef PAIRED_READ_HANDLER_H
#define PAIRED_READ_HANDLER_H

#include <htslib/sam.h>
#include <unordered_set>
#include <unordered_map>
#include <string>
#include <vector>

struct paired_read_tracker_t
{
    std::unordered_set<std::string> passed_filter_read_names;
    std::unordered_map<std::string, bam1_t*> failed_mates_buffer;
    
    paired_read_tracker_t() = default;
    
    // Destructor to clean up any remaining failed mates
    ~paired_read_tracker_t()
    {
        clear_failed_mates();
    }
    
    void clear_failed_mates()
    {
        for (auto& [qname, failed_mate] : failed_mates_buffer)
        {
            bam_destroy1(failed_mate);
        }
        failed_mates_buffer.clear();
        passed_filter_read_names.clear();
    }
};

bool handle_paired_read(
    bam1_t* alignment,
    bool passes_mapq,
    std::vector<bam1_t*>& alignment_buffer,
    paired_read_tracker_t& tracker,
    bool passes_on_target = true)
{
    bool pair_passed = false;

    std::string qname = bam_get_qname(alignment);
    bool mate_passed = tracker.passed_filter_read_names.count(qname) > 0;
    
    if (passes_mapq && passes_on_target)
    {
        // This read passes - add it to buffer
        bam1_t* new_alignment = bam_init1();
        bam_copy1(new_alignment, alignment);
        alignment_buffer.push_back(new_alignment);
        tracker.passed_filter_read_names.insert(qname);

        // Check if we have a failed mate waiting in this buffer
        auto it = tracker.failed_mates_buffer.find(qname);
        if (it != tracker.failed_mates_buffer.end())
        {
            // Add the failed mate to buffer (will be rescued)
            alignment_buffer.push_back(it->second);
            tracker.failed_mates_buffer.erase(it);
        }

        pair_passed = true;
    }
    else if (mate_passed)
    {
        // This read failed but mate passed in this buffer - rescue it
        bam1_t* new_alignment = bam_init1();
        bam_copy1(new_alignment, alignment);
        alignment_buffer.push_back(new_alignment);
        tracker.passed_filter_read_names.erase(qname); // Clean up tracking
    }
    else
    {
        // This read failed and mate not seen yet in this buffer - cache it
        bam1_t* failed_read = bam_init1();
        bam_copy1(failed_read, alignment);
        
        // Replace if already exists (shouldn't happen normally)
        auto it = tracker.failed_mates_buffer.find(qname);
        if (it != tracker.failed_mates_buffer.end())
        {
            bam_destroy1(it->second);
        }
        tracker.failed_mates_buffer[qname] = failed_read;
    }

    return pair_passed;
}

void rescue_failed_mates(
    std::vector<bam1_t*>& alignment_buffer,
    paired_read_tracker_t& tracker)
{
    // Add any rescued mates from failed_mates_buffer
    for (const auto& qname_in_buffer : tracker.passed_filter_read_names)
    {
        auto it = tracker.failed_mates_buffer.find(qname_in_buffer);
        if (it != tracker.failed_mates_buffer.end())
        {
            alignment_buffer.push_back(it->second);
            tracker.failed_mates_buffer.erase(it);
        }
    }
    
    // Discard any failed mates that were never rescued
    for (auto& [qname, failed_mate] : tracker.failed_mates_buffer)
    {
        bam_destroy1(failed_mate);
    }
    tracker.failed_mates_buffer.clear();
    tracker.passed_filter_read_names.clear();
}

void cleanup_pair_tracker(paired_read_tracker_t& tracker)
{
    tracker.clear_failed_mates();
}

#endif // PAIRED_READ_HANDLER_H