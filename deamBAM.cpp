#include <iostream>
#include <htslib/sam.h>
#include <vector>
#include <thread>
#include <numeric>
#include <fstream>
#include <map>
#include <cstring>
#include <sstream>
#include <string>
#include <iomanip>

#include "bam_processing.h"

#define END_LENGTH 21
#define VERSION_NUMBER 0.1

struct substitution_patterns_t
{
    uint64_t AA = 0;
    uint64_t AC = 0;
    uint64_t AG = 0;
    uint64_t AT = 0;
    uint64_t AN = 0;

    uint64_t GG = 0;
    uint64_t GA = 0;
    uint64_t GC = 0;
    uint64_t GT = 0;
    uint64_t GN = 0;

    uint64_t CC = 0;
    uint64_t CA = 0;
    uint64_t CG = 0;
    uint64_t CT = 0;
    uint64_t CN = 0;

    uint64_t TT = 0;
    uint64_t TA = 0;
    uint64_t TC = 0;
    uint64_t TG = 0;
    uint64_t TN = 0;

    bool is_non_zero() const
    {
        return AA || AC || AG || AT || AN ||
               GG || GA || GC || GT || GN ||
               CC || CA || CG || CT || CN ||
               TT || TA || TC || TG || TN;
    }
};

void increment_subst_pattern(
    substitution_patterns_t& pattern,
    char ref_base,
    char aln_base)
{
    if (ref_base == 'A')
    {
        if (aln_base == 'A') pattern.AA++;
        else if (aln_base == 'C') pattern.AC++;
        else if (aln_base == 'G') pattern.AG++;
        else if (aln_base == 'T') pattern.AT++;
        else if (aln_base == 'N') pattern.AN++;
    }
    else if (ref_base == 'C') 
    {
        if (aln_base == 'A') pattern.CA++;
        else if (aln_base == 'C') pattern.CC++;
        else if (aln_base == 'G') pattern.CG++;
        else if (aln_base == 'T') pattern.CT++;
        else if (aln_base == 'N') pattern.CN++;
    }
    else if (ref_base == 'G') 
    {
        if (aln_base == 'A') pattern.GA++;
        else if (aln_base == 'C') pattern.GC++;
        else if (aln_base == 'G') pattern.GG++;
        else if (aln_base == 'T') pattern.GT++;
        else if (aln_base == 'N') pattern.GN++;
    }
    else if (ref_base == 'T')
    {
        if (aln_base == 'A') pattern.TA++;
        else if (aln_base == 'C') pattern.TC++;
        else if (aln_base == 'G') pattern.TG++;
        else if (aln_base == 'T') pattern.TT++;
        else if (aln_base == 'N') pattern.TN++;
    }
}

void process_chunk(
    int chunk_num,
    std::string bam_file_location,
    int tid,
    int start,
    int end,
    std::map<int, std::vector<std::vector<substitution_patterns_t>>>& subst_patterns,
    std::map<int, std::vector<std::vector<substitution_patterns_t>>>& subst_patterns_deam_5,
    std::map<int, std::vector<std::vector<substitution_patterns_t>>>& subst_patterns_deam_3,
    std::vector<std::vector<int>>& read_length_distribution_per_thread)
{
    bam_file_config_t bam_config = {};
    bam_constructor(bam_file_location, &bam_config);
    hts_itr_t* iter = sam_itr_queryi(bam_config.index, tid, start, end);

    bam1_t* alignment = bam_init1();
    while (sam_itr_next(bam_config.bam_file, iter, alignment) >= 0)
    {
        int read_length = alignment->core.l_qseq;

        read_length_distribution_per_thread[chunk_num][read_length]++;

        bool is_reverse = alignment->core.flag & BAM_FREVERSE;

        std::string reconstructed_alignment = reconstruct_alignment(alignment);
        std::string reconstructed_reference = reconstruct_reference(alignment, reconstructed_alignment);

        // bool valid_bases = true;
        // for (char base : {reconstructed_alignment.front(), reconstructed_alignment.back(), reconstructed_reference.front(), reconstructed_reference.back()})
        // {
        //     if (base != 'A' && base != 'C' && base != 'G' && base != 'T')
        //     {
        //         valid_bases = false;
        //         break;
        //     }
        // }

        // if (!valid_bases)
        // {
        //     continue;
        // }

        if (is_reverse)
        {
            reconstructed_alignment = reverse_complement(reconstructed_alignment);
            reconstructed_reference = reverse_complement(reconstructed_reference);
        }

        std::string alignment_ends;
        std::string reference_ends;

        for (int i = 0; i < END_LENGTH; i++)
        {
            alignment_ends += reconstructed_alignment[i];
            reference_ends += reconstructed_reference[i];
        }

        for (int i = reconstructed_alignment.size() - END_LENGTH; i < reconstructed_alignment.size(); i++)
        {
            alignment_ends += reconstructed_alignment[i];
            reference_ends += reconstructed_reference[i];
        }

        bool deam_5 = false;
        bool deam_3 = false;

        if (reference_ends[0] == 'C' && alignment_ends[0] == 'T')
        {
            deam_5 = true;
        }

        if (reference_ends.back() == 'C' && alignment_ends.back() == 'T')
        {
            deam_3 = true;
        }

        for (int i = 0; i < reference_ends.size(); i++)
        {
            increment_subst_pattern(subst_patterns[read_length][i][chunk_num], reference_ends[i], alignment_ends[i]);
            if (deam_5)
            {
                increment_subst_pattern(subst_patterns_deam_5[read_length][i][chunk_num], reference_ends[i], alignment_ends[i]);
            }
            if (deam_3)
            {
                increment_subst_pattern(subst_patterns_deam_3[read_length][i][chunk_num], reference_ends[i], alignment_ends[i]);
            }
        }
    }

    bam_destructor(&bam_config);
}

void update_subst_patterns(
    std::map<int, std::vector<substitution_patterns_t>>& subst_patterns_total,
    const std::map<int, std::vector<substitution_patterns_t>>& patterns)
{
    for (const auto& [read_length, pattern_vectors] : patterns)
    {
        for (int i = 0; i < END_LENGTH * 2; i++)
        {
            subst_patterns_total[read_length][i].AA += pattern_vectors[i].AA;
            subst_patterns_total[read_length][i].AC += pattern_vectors[i].AC;
            subst_patterns_total[read_length][i].AG += pattern_vectors[i].AG;
            subst_patterns_total[read_length][i].AT += pattern_vectors[i].AT;
            subst_patterns_total[read_length][i].AN += pattern_vectors[i].AN;

            subst_patterns_total[read_length][i].GG += pattern_vectors[i].GG;
            subst_patterns_total[read_length][i].GA += pattern_vectors[i].GA;
            subst_patterns_total[read_length][i].GC += pattern_vectors[i].GC;
            subst_patterns_total[read_length][i].GT += pattern_vectors[i].GT;
            subst_patterns_total[read_length][i].GN += pattern_vectors[i].GN;

            subst_patterns_total[read_length][i].CC += pattern_vectors[i].CC;
            subst_patterns_total[read_length][i].CA += pattern_vectors[i].CA;
            subst_patterns_total[read_length][i].CG += pattern_vectors[i].CG;
            subst_patterns_total[read_length][i].CT += pattern_vectors[i].CT;
            subst_patterns_total[read_length][i].CN += pattern_vectors[i].CN;

            subst_patterns_total[read_length][i].TT += pattern_vectors[i].TT;
            subst_patterns_total[read_length][i].TA += pattern_vectors[i].TA;
            subst_patterns_total[read_length][i].TC += pattern_vectors[i].TC;
            subst_patterns_total[read_length][i].TG += pattern_vectors[i].TG;
            subst_patterns_total[read_length][i].TN += pattern_vectors[i].TN;
        }
    }
}

void update_subst_patterns_threads(
    uint64_t num_threads,
    std::map<int, std::vector<substitution_patterns_t>>& subst_patterns_total,
    const std::map<int, std::vector<std::vector<substitution_patterns_t>>>& patterns)
{
    for (const auto& [read_length, pattern_vectors] : patterns)
    {
        for (int i = 0; i < num_threads; i++)
        {
            for (int j = 0; j < END_LENGTH * 2; j++)
            {
                subst_patterns_total[read_length][j].AA += pattern_vectors[j][i].AA;
                subst_patterns_total[read_length][j].AC += pattern_vectors[j][i].AC;
                subst_patterns_total[read_length][j].AG += pattern_vectors[j][i].AG;
                subst_patterns_total[read_length][j].AT += pattern_vectors[j][i].AT;
                subst_patterns_total[read_length][j].AN += pattern_vectors[j][i].AN;

                subst_patterns_total[read_length][j].GG += pattern_vectors[j][i].GG;
                subst_patterns_total[read_length][j].GA += pattern_vectors[j][i].GA;
                subst_patterns_total[read_length][j].GC += pattern_vectors[j][i].GC;
                subst_patterns_total[read_length][j].GT += pattern_vectors[j][i].GT;
                subst_patterns_total[read_length][j].GN += pattern_vectors[j][i].GN;

                subst_patterns_total[read_length][j].CC += pattern_vectors[j][i].CC;
                subst_patterns_total[read_length][j].CA += pattern_vectors[j][i].CA;
                subst_patterns_total[read_length][j].CG += pattern_vectors[j][i].CG;
                subst_patterns_total[read_length][j].CT += pattern_vectors[j][i].CT;
                subst_patterns_total[read_length][j].CN += pattern_vectors[j][i].CN;

                subst_patterns_total[read_length][j].TT += pattern_vectors[j][i].TT;
                subst_patterns_total[read_length][j].TA += pattern_vectors[j][i].TA;
                subst_patterns_total[read_length][j].TC += pattern_vectors[j][i].TC;
                subst_patterns_total[read_length][j].TG += pattern_vectors[j][i].TG;
                subst_patterns_total[read_length][j].TN += pattern_vectors[j][i].TN;
            }
        }
    }
}

bool get_substitution_patterns_for_chr(
    const int num_threads,
    const int min_len,
    const int max_len,
    const std::string bam_file_location,
    const std::string chr,
    std::map<int, std::vector<substitution_patterns_t>>& subst_patterns_total,
    std::map<int, std::vector<substitution_patterns_t>>& subst_patterns_deam_5_total,
    std::map<int, std::vector<substitution_patterns_t>>& subst_patterns_deam_3_total,
    std::map<std::string, std::vector<int>>& read_length_distribution_per_chromosome)
{
    std::map<int, std::vector<std::vector<substitution_patterns_t>>> subst_patterns;
    std::map<int, std::vector<std::vector<substitution_patterns_t>>> subst_patterns_deam_5;
    std::map<int, std::vector<std::vector<substitution_patterns_t>>> subst_patterns_deam_3;
    for (int i = min_len; i <= max_len; i++)
    {
        subst_patterns[i] = std::vector<std::vector<substitution_patterns_t>>(END_LENGTH * 2, std::vector<substitution_patterns_t>(num_threads));
        subst_patterns_deam_5[i] = std::vector<std::vector<substitution_patterns_t>>(END_LENGTH * 2, std::vector<substitution_patterns_t>(num_threads));
        subst_patterns_deam_3[i] = std::vector<std::vector<substitution_patterns_t>>(END_LENGTH * 2, std::vector<substitution_patterns_t>(num_threads));
    }

    // Create a read length distribution vector for each thread
    std::vector<std::vector<int>> read_length_distribution_per_thread(num_threads, std::vector<int>(max_len, 0));

    std::pair<int, int> tid_chunk_size = check_bam(bam_file_location, chr, num_threads);
    if (tid_chunk_size.first < 0)
    {
        return false;
    }
    int tid = tid_chunk_size.first;
    int chunk_size = tid_chunk_size.second;
    if (tid < 0)
    {
        return false;
    }

    std::vector<std::thread> threads;
    for (int i = 0; i < num_threads; i++)
    {
        int start = i * chunk_size;
        int end = start + chunk_size;

        // process_chunk(i, bam_file_location, tid, start, end, subst_patterns, subst_patterns_deam_5, subst_patterns_deam_3);
        threads.push_back(std::thread(process_chunk, i, bam_file_location, tid, start, end, std::ref(subst_patterns), std::ref(subst_patterns_deam_5), std::ref(subst_patterns_deam_3), std::ref(read_length_distribution_per_thread)));
    }

    for (auto& thread : threads)
    {
        thread.join();
    }

    // Combine the read length distribution from all threads
    for (int i = 0; i < num_threads; i++)
    {
        for (int j = 0; j < max_len; j++)
        {
            read_length_distribution_per_chromosome[chr][j] += read_length_distribution_per_thread[i][j];
        }
    }

    // Combine the substitution patterns from all threads
    update_subst_patterns_threads(num_threads, subst_patterns_total, subst_patterns);
    update_subst_patterns_threads(num_threads, subst_patterns_deam_5_total, subst_patterns_deam_5);
    update_subst_patterns_threads(num_threads, subst_patterns_deam_3_total, subst_patterns_deam_3);

    return true;
}

void write_ct_stats(
    std::string output_file_name,
    const int min_len,
    const int max_len,
    std::map<int, std::vector<substitution_patterns_t>>& subst_patterns_total,
    std::map<int, std::vector<substitution_patterns_t>>& subst_patterns_deam_5_total,
    std::map<int, std::vector<substitution_patterns_t>>& subst_patterns_deam_3_total
)
{
    std::ofstream outfile(output_file_name);

    outfile << "len\t5'CT\t5'CT_95CI\t5'#refC\tcond5'CT\tcond5'CT_95CI\tcond5'#refC\t3'CT\t3'CT_95CI\t3'#refC\tcond3'CT\tcond3'CT_95CI\tcond3'#refC\n";

    uint64_t cond_subst_ct_3_total = 0;
    uint64_t cond_subst_c_3_total = 0;

    uint64_t cond_subst_ct_5_total = 0;
    uint64_t cond_subst_c_5_total = 0;

    uint64_t subst_ct_3_total = 0;
    uint64_t subst_c_3_total = 0;

    uint64_t subst_ct_5_total = 0;
    uint64_t subst_c_5_total = 0;

    for (int i = min_len; i <= max_len; i++)
    {
        uint64_t cond_subst_c_3 = subst_patterns_deam_5_total[i][END_LENGTH * 2 - 1].CC + subst_patterns_deam_5_total[i][END_LENGTH * 2 - 1].CT + subst_patterns_deam_5_total[i][END_LENGTH * 2 - 1].CG + subst_patterns_deam_5_total[i][END_LENGTH * 2 - 1].CA + subst_patterns_deam_5_total[i][END_LENGTH * 2 - 1].CN;
        uint64_t cond_subst_ct_3 = subst_patterns_deam_5_total[i][END_LENGTH * 2 - 1].CT;

        cond_subst_ct_3_total += cond_subst_ct_3;
        cond_subst_c_3_total += cond_subst_c_3;

        uint64_t cond_subst_c_5 = subst_patterns_deam_3_total[i][0].CC + subst_patterns_deam_3_total[i][0].CT + subst_patterns_deam_3_total[i][0].CG + subst_patterns_deam_3_total[i][0].CA + subst_patterns_deam_3_total[i][0].CN;
        uint64_t cond_subst_ct_5 = subst_patterns_deam_3_total[i][0].CT;

        cond_subst_ct_5_total += cond_subst_ct_5;
        cond_subst_c_5_total += cond_subst_c_5;

        uint64_t subst_c_3 = subst_patterns_total[i][END_LENGTH * 2 - 1].CC + subst_patterns_total[i][END_LENGTH * 2 - 1].CT + subst_patterns_total[i][END_LENGTH * 2 - 1].CG + subst_patterns_total[i][END_LENGTH * 2 - 1].CA + subst_patterns_total[i][END_LENGTH * 2 - 1].CN;
        uint64_t subst_ct_3 = subst_patterns_total[i][END_LENGTH * 2 - 1].CT;

        subst_ct_3_total += subst_ct_3;
        subst_c_3_total += subst_c_3;

        uint64_t subst_c_5 = subst_patterns_total[i][0].CC + subst_patterns_total[i][0].CT + subst_patterns_total[i][0].CG + subst_patterns_total[i][0].CA + subst_patterns_total[i][0].CN;
        uint64_t subst_ct_5 = subst_patterns_total[i][0].CT;

        subst_ct_5_total += subst_ct_5;
        subst_c_5_total += subst_c_5;

        auto cond_ci_3 = binomial_confidence_interval(cond_subst_ct_3, cond_subst_c_3);
        auto cond_ci_3_str = binomial_ci_string(cond_ci_3);
    
        auto cond_ci_5 = binomial_confidence_interval(cond_subst_ct_5, cond_subst_c_5);
        auto cond_ci_5_str = binomial_ci_string(cond_ci_5);
    
        auto ci_3 = binomial_confidence_interval(subst_ct_3, subst_c_3);
        auto ci_3_str = binomial_ci_string(ci_3);

        auto ci_5 = binomial_confidence_interval(subst_ct_5, subst_c_5);
        auto ci_5_str = binomial_ci_string(ci_5);

        if (subst_c_5 == 0 || subst_c_3 == 0)
        {
            continue;
        }

        outfile << i << "\t"
                << std::fixed << std::setprecision(2) << ((double)subst_ct_5 / subst_c_5) * 100 << "\t" << ci_5_str << "\t" << subst_c_5
                << "\t" << ((double)cond_subst_ct_5 / cond_subst_c_5) * 100 << "\t" << cond_ci_5_str << "\t" << cond_subst_c_5
                << "\t" << ((double)subst_ct_3 / subst_c_3) * 100 << "\t" << ci_3_str << "\t" << subst_c_3
                << "\t" << ((double)cond_subst_ct_3 / cond_subst_c_3) * 100 << "\t" << cond_ci_3_str << "\t" << cond_subst_c_3 << "\n";
    }

    auto ci_3_total = binomial_confidence_interval(subst_ct_3_total, subst_c_3_total);
    auto ci_3_total_str = binomial_ci_string(ci_3_total);

    auto ci_5_total = binomial_confidence_interval(subst_ct_5_total, subst_c_5_total);
    auto ci_5_total_str = binomial_ci_string(ci_5_total);

    auto cond_ci_3_total = binomial_confidence_interval(cond_subst_ct_3_total, cond_subst_c_3_total);
    auto cond_ci_3_total_str = binomial_ci_string(cond_ci_3_total);

    auto cond_ci_5_total = binomial_confidence_interval(cond_subst_ct_5_total, cond_subst_c_5_total);
    auto cond_ci_5_total_str = binomial_ci_string(cond_ci_5_total);
    outfile << "total\t" << std::fixed << std::setprecision(2) << ((double)subst_ct_5_total / subst_c_5_total) * 100 << "\t" << ci_5_total_str << "\t" << subst_c_5_total
            << "\t" << ((double)cond_subst_ct_5_total / cond_subst_c_5_total) * 100 << "\t" << cond_ci_5_total_str << "\t" << cond_subst_c_5_total
            << "\t" << ((double)subst_ct_3_total / subst_c_3_total) * 100 << "\t" << ci_3_total_str << "\t" << subst_c_3_total
            << "\t" << ((double)cond_subst_ct_3_total / cond_subst_c_3_total) * 100 << "\t" << cond_ci_3_total_str << "\t" << cond_subst_c_3_total << "\n";

    outfile.close();
}

void write_subst_patterns(
    const std::string input_bam_name,
    const std::string out_folder,
    const std::map<int, std::vector<substitution_patterns_t>>& subst_patterns_total
)
{
    for (const auto& [read_length, patterns] : subst_patterns_total)
    {
        // Check if all patterns are zero
        bool all_zero = true;
        for (const auto &pattern : patterns)
        {
            if (pattern.is_non_zero())
            {
                all_zero = false;
                break;
            }
        }
        if (all_zero)
        {
            continue;
        }

        std::string filename = out_folder + "substitution_patterns_" + std::to_string(read_length) + "." + input_bam_name + ".txt";
        std::ofstream outfile(filename);

        if (!outfile.is_open())
        {
            std::cerr << "Error opening file: " << filename << std::endl;
            continue;
        }

        outfile << "Position\tAA\tAC\tAG\tAT\tAN\tGG\tGA\tGC\tGT\tGN\tCC\tCA\tCG\tCT\tCN\tTT\tTA\tTC\tTG\tTN\n";
        for (int i = 0; i < END_LENGTH; i++)
        {
            outfile << i << "\t" << patterns[i].AA << "\t" << patterns[i].AC << "\t" << patterns[i].AG << "\t" << patterns[i].AT << "\t" << patterns[i].AN << "\t" << patterns[i].GG << "\t" << patterns[i].GA << "\t" << patterns[i].GC << "\t" << patterns[i].GT << "\t" << patterns[i].GN << "\t" << patterns[i].CC << "\t" << patterns[i].CA << "\t" << patterns[i].CG << "\t" << patterns[i].CT << "\t" << patterns[i].CN << "\t" << patterns[i].TT << "\t" << patterns[i].TA << "\t" << patterns[i].TC << "\t" << patterns[i].TG << "\t" << patterns[i].TN << "\n";
        }

        for (int i = END_LENGTH; i < END_LENGTH * 2; i++)
        {
            outfile << i - END_LENGTH * 2 + 1 << "\t" << patterns[i].AA << "\t" << patterns[i].AC << "\t" << patterns[i].AG << "\t" << patterns[i].AT << "\t" << patterns[i].AN << "\t" << patterns[i].GG << "\t" << patterns[i].GA << "\t" << patterns[i].GC << "\t" << patterns[i].GT << "\t" << patterns[i].GN << "\t" << patterns[i].CC << "\t" << patterns[i].CA << "\t" << patterns[i].CG << "\t" << patterns[i].CT << "\t" << patterns[i].CN << "\t" << patterns[i].TT << "\t" << patterns[i].TA << "\t" << patterns[i].TC << "\t" << patterns[i].TG << "\t" << patterns[i].TN << "\n";
        }

        outfile.close();
    }
}

void write_read_length_dist(
    const int max_len,
    const std::string input_bam_name,
    const std::string out_folder,
    std::vector<std::string>& chromosomes,
    std::map<std::string, std::vector<int>>& read_length_distribution_per_chromosome
)
{
    std::string read_length_distribution_file_name = out_folder + "read_length_distribution." + input_bam_name + ".txt";
    std::ofstream read_length_distribution_file(read_length_distribution_file_name);

    // Write header
    read_length_distribution_file << "Length\tTotal";

    // Only include chromosomes with at least one read
    std::vector<std::string> nonzero_chromosomes;
    for (const auto& chr : chromosomes)
    {
        bool has_reads = false;
        for (int i = 0; i < max_len; i++)
        {
            if (read_length_distribution_per_chromosome[chr][i] > 0)
            {
                has_reads = true;
                break;
            }
        }
        if (has_reads)
        {
            nonzero_chromosomes.push_back(chr);
            read_length_distribution_file << "\t" << chr;
        }
    }
    read_length_distribution_file << "\n";

    for (int i = 0; i < max_len; i++)
    {
        int total = 0;
        for (const auto& chr : nonzero_chromosomes)
        {
            total += read_length_distribution_per_chromosome[chr][i];
        }
        read_length_distribution_file << i << "\t" << total;
        for (const auto& chr : nonzero_chromosomes)
        {
            read_length_distribution_file << "\t" << read_length_distribution_per_chromosome[chr][i];
        }
        read_length_distribution_file << "\n";
    }
    read_length_distribution_file.close();
}

void print_help()
{
    std::cout << "Usage: deamBAM [options] <BAM file>\n";
    std::cout << "Options:\n";
    std::cout << "  -threads <num>        Number of threads to use (default: 40)\n";
    std::cout << "  -out_folder <path>    Output folder (default: current directory)\n";
    std::cout << "  -min_len <num>        Minimum read length (default: 35)\n";
    std::cout << "  -max_len <num>        Maximum read length (default: 300)\n";
    std::cout << "  -help                 Show this help message\n";
}

int main(int argc, char *argv[])
{
    int32_t min_len = 35;
    int32_t max_len = 300;
    int32_t num_threads = 40;
    std::string out_folder = "";
    std::string bam_file_location = "";

    std::cout << "deamBAM.cpp v" << VERSION_NUMBER << " - " << get_date_time() << std::endl;

    for (int i = 0; i < argc; i++)
    {
        if (std::string(argv[i]) == "-threads")
        {
            num_threads = std::stoi(argv[i + 1]);
        }
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
        }
        if (std::string(argv[i]) == "-max_len")
        {
            max_len = std::stoi(argv[i + 1]);
        }
        if (std::string(argv[i]).find(".bam") != std::string::npos)
        {
            bam_file_location = argv[i];
        }
        if ((std::string(argv[i]) == "-help") || (argc == 1))
        {
            print_help();
            return 0;
        }
    }

    std::vector<std::string> chromosomes;
    get_chromosomes_from_bam(chromosomes, bam_file_location);
    if (chromosomes.empty())
    {
        std::cerr << "Error: No chromosomes found in the BAM file" << std::endl;
        return 1;
    }

    std::map<std::string, std::vector<int>> read_length_distribution_per_chromosome;
    for (const auto& chr : chromosomes) {
        read_length_distribution_per_chromosome[chr] = std::vector<int>(max_len, 0);
    }

    std::map<int, std::vector<substitution_patterns_t>> subst_patterns_total;
    std::map<int, std::vector<substitution_patterns_t>> subst_patterns_deam_5_total;
    std::map<int, std::vector<substitution_patterns_t>> subst_patterns_deam_3_total;
    for (int i = min_len; i <= max_len; i++)
    {
        subst_patterns_total[i] = std::vector<substitution_patterns_t>(END_LENGTH * 2);
        subst_patterns_deam_5_total[i] = std::vector<substitution_patterns_t>(END_LENGTH * 2);
        subst_patterns_deam_3_total[i] = std::vector<substitution_patterns_t>(END_LENGTH * 2);
    }

    for (const std::string& chr : chromosomes)
    {
        std::cout << chr << "\n";

        std::map<int, std::vector<substitution_patterns_t>> subst_patterns;
        std::map<int, std::vector<substitution_patterns_t>> subst_patterns_deam_5;
        std::map<int, std::vector<substitution_patterns_t>> subst_patterns_deam_3;
        for (int i = min_len; i <= max_len; i++)
        {
            subst_patterns[i] = std::vector<substitution_patterns_t>(END_LENGTH * 2);
            subst_patterns_deam_5[i] = std::vector<substitution_patterns_t>(END_LENGTH * 2);
            subst_patterns_deam_3[i] = std::vector<substitution_patterns_t>(END_LENGTH * 2);
        }

        if (!get_substitution_patterns_for_chr(num_threads, min_len, max_len, bam_file_location, chr, subst_patterns, subst_patterns_deam_5, subst_patterns_deam_3, read_length_distribution_per_chromosome))
        {
            return 1;
        }
        
        update_subst_patterns(subst_patterns_total, subst_patterns);
        update_subst_patterns(subst_patterns_deam_5_total, subst_patterns_deam_5);
        update_subst_patterns(subst_patterns_deam_3_total, subst_patterns_deam_3);
    }

    // Extract the file name without the extension and path
    std::string input_bam_name;
    input_bam_name = extract_file_name(bam_file_location);

    std::string output_file_name = out_folder + "substitution_patterns." + input_bam_name + ".txt";

    write_ct_stats(output_file_name, min_len, max_len, subst_patterns_total, subst_patterns_deam_5_total, subst_patterns_deam_3_total);
    write_subst_patterns(input_bam_name, out_folder, subst_patterns_total);
    write_read_length_dist(max_len, input_bam_name, out_folder, chromosomes, read_length_distribution_per_chromosome);

    return 0;
}
