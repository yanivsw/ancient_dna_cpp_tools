#include <iostream>
#include <htslib/vcf.h>
#include <vector>
#include <numeric>
#include <thread>
#include <cmath>
#include <fstream>
#include <unordered_map>
#include <filesystem>

#include "bed_processing.h"
#include "utils.h"
#include "chromosome_sizes.h"

#define VERSION_NUMBER 0.1

struct non_ref_stats_t
{
    std::string chr;
    uint64_t n_bp;
    uint64_t n_hets;
    uint64_t n_hets_transversions;
    float heterozygosity;
    uint64_t ref_count;
    uint64_t alt_count;
    float alt_allele_prop;
    float het_coverage;
    std::vector<uint32_t> het_dist;
    std::vector<float> het_gq_dist;
    std::unordered_map<uint64_t, uint64_t> het_gq_hist;
    float het_gq_mean;
    float het_gq_std;
    uint64_t n_ts;
    uint64_t n_tv;
    float ts_tv_ratio;
    std::vector<float> alt_gq_dist;
    std::unordered_map<uint64_t, uint64_t> alt_gq_hist;
    float alt_gq_mean;
    float alt_gq_std;
};

void get_coverage(bcf_hdr_t* header,
                  bcf1_t* record,
                  int sample_index,
                  int chr_index,
                  std::vector<std::vector<uint64_t>>& coverage_dist)
{
    // Get the coverage from the FORMAT field
    int32_t* coverage = nullptr;
    int len_coverage = 0;
    if (bcf_get_format_int32(header, record, "DP", &coverage, &len_coverage) > 0)
    {
        auto cov = coverage[sample_index];
        coverage_dist[cov][chr_index]++;
    }

    // Free the memory allocated by bcf_get_format_int32
    if (coverage != nullptr)
    {
        free(coverage);
    }
}

float get_gq(bcf_hdr_t* header,
                bcf1_t* record,
                int sample_index,
                int chr_index)
{
    // Get the GQ from the FORMAT field
    float* gq = nullptr;
    int len_gq = 0;
    uint64_t gq_value = 0;
    if (bcf_get_format_float(header, record, "GQ", &gq, &len_gq) > 0)
    {
        gq_value = gq[sample_index];
    }

    // Free the memory allocated by bcf_get_format_int32
    if (gq != nullptr)
    {
        free(gq);
    }

    return gq_value;
}

bool is_transition(const std::string& base1, const std::string& base2)
{
    return (base1 == "A" && base2 == "G") || (base1 == "G" && base2 == "A") ||
           (base1 == "C" && base2 == "T") || (base1 == "T" && base2 == "C");
}

int get_non_ref_stats(std::string vcf_file,
                      std::string filter_file,
                      std::string sample_name,
                      std::string chr,
                      int chr_index,
                      std::vector<non_ref_stats_t>& non_ref_stats,
                      std::unordered_map<std::string, long>& len_map,
                      uint64_t het_block_size = 10000)
{
    uint64_t het_gq_total = 0;
    uint64_t alt_gq_total = 0;

    non_ref_stats[chr_index].chr = chr;

    htsFile* vcf = bcf_open(vcf_file.c_str(), "r");
    if (vcf == nullptr)
    {
        std::cerr << "Failed to open VCF file: " << vcf_file << std::endl;
        return 1;
    }

    // Read the VCF header
    bcf_hdr_t* header = bcf_hdr_read(vcf);
    if (header == nullptr)
    {
        std::cerr << "Failed to read VCF header" << std::endl;
        bcf_close(vcf);
        return 1;
    }

    int sample_index = bcf_hdr_id2int(header, BCF_DT_SAMPLE, sample_name.c_str());
    if (sample_index < 0)
    {
        std::cerr << "Sample not found: " << sample_name << std::endl;
        bcf_hdr_destroy(header);
        bcf_close(vcf);
        return 1;
    }

    non_ref_stats[chr_index].n_bp = get_bed_size(filter_file, chr);
    non_ref_stats[chr_index].het_dist.resize(std::ceil(static_cast<double>(len_map[chr]) / het_block_size), 0);

    // Iterate over the VCF records
    bcf1_t* record = bcf_init();
    while (bcf_read(vcf, header, record) == 0)
    {
        int* genotype = nullptr;
        int ngt = 0;
        ngt = bcf_get_genotypes(header, record, &genotype, &ngt);

        int allele1 = bcf_gt_allele(genotype[0]);
        int allele2 = bcf_gt_allele(genotype[1]);

        bcf_unpack(record, BCF_UN_STR);

        std::string base1 = record->d.allele[allele1];
        std::string base2 = record->d.allele[allele2];

        // If this is homozygous alt
        if ((allele1 == allele2) && (allele1 == 1))
        {
            std::string ref = record->d.allele[0];
            std::string alt = record->d.allele[1];

            if (is_transition(ref, alt))
            {
                non_ref_stats[chr_index].n_ts++;
            }
            else
            {
                non_ref_stats[chr_index].n_tv++;
            }

            auto gq = get_gq(header, record, sample_index, chr_index);
            non_ref_stats[chr_index].alt_gq_dist.push_back(gq);
        }

        // If this is heterozygous
        if ((allele1 != allele2) && (allele1 == 0))
        {
            // Check if the alleles are transversions
            if (!is_transition(base1, base2))
            {
                int32_t* allele1_cov = NULL;
                int len_allele1_cov = 0;
                len_allele1_cov = bcf_get_format_int32(header, record, base1.c_str(), &allele1_cov, &len_allele1_cov);
                int cov_allele1 = allele1_cov[0] + allele1_cov[1];

                int32_t* allele2_cov = NULL;
                int len_allele2_cov = 0;
                len_allele2_cov = bcf_get_format_int32(header, record, base2.c_str(), &allele2_cov, &len_allele2_cov);
                int cov_allele2 = allele2_cov[0] + allele2_cov[1];

                non_ref_stats[chr_index].ref_count += cov_allele1;
                non_ref_stats[chr_index].alt_count += cov_allele2;
                non_ref_stats[chr_index].n_hets_transversions++;

                free(allele1_cov);
                free(allele2_cov);
            }

            non_ref_stats[chr_index].n_hets++;

            auto block_index = record->pos / het_block_size;
            non_ref_stats[chr_index].het_dist[block_index]++;

            auto gq = get_gq(header, record, sample_index, chr_index);
            non_ref_stats[chr_index].het_gq_dist.push_back(gq);
        }
        // If this is heterozygous with non-reference alleles
        else if ((allele1 != allele2) && (allele1 != 0))
        {
            non_ref_stats[chr_index].n_hets++;

            auto block_index = record->pos / het_block_size;
            non_ref_stats[chr_index].het_dist[block_index]++;

            auto gq = get_gq(header, record, sample_index, chr_index);
            non_ref_stats[chr_index].het_gq_dist.push_back(gq);
        }

        // Free the memory allocated by bcf_get_genotypes
        if (genotype != nullptr)
        {
            free(genotype);
        }
    }

    non_ref_stats[chr_index].alt_allele_prop = (float)non_ref_stats[chr_index].alt_count / (non_ref_stats[chr_index].alt_count + non_ref_stats[chr_index].ref_count);
    non_ref_stats[chr_index].heterozygosity = ( (float)non_ref_stats[chr_index].n_hets / non_ref_stats[chr_index].n_bp ) * 10000;
    non_ref_stats[chr_index].het_coverage = (float)(non_ref_stats[chr_index].alt_count + non_ref_stats[chr_index].ref_count) / non_ref_stats[chr_index].n_hets_transversions;
    non_ref_stats[chr_index].ts_tv_ratio = (float)non_ref_stats[chr_index].n_ts / non_ref_stats[chr_index].n_tv;
    non_ref_stats[chr_index].het_gq_mean = (float)std::accumulate(non_ref_stats[chr_index].het_gq_dist.begin(), non_ref_stats[chr_index].het_gq_dist.end(), 0) / non_ref_stats[chr_index].n_hets;
    non_ref_stats[chr_index].alt_gq_mean = (float)std::accumulate(non_ref_stats[chr_index].alt_gq_dist.begin(), non_ref_stats[chr_index].alt_gq_dist.end(), 0) / (non_ref_stats[chr_index].n_ts + non_ref_stats[chr_index].n_tv);

    for (const auto& gq : non_ref_stats[chr_index].het_gq_dist)
    {
        int bin = static_cast<int>(std::floor(gq));
        non_ref_stats[chr_index].het_gq_hist[bin]++;
    }

    for (const auto& gq : non_ref_stats[chr_index].alt_gq_dist)
    {
        int bin = static_cast<int>(std::floor(gq));
        non_ref_stats[chr_index].alt_gq_hist[bin]++;
    }

    // Clean up
    bcf_destroy(record);
    bcf_hdr_destroy(header);
    bcf_close(vcf);
}

void print_help()
{
    std::cout << "Usage: analyzeVCF -sample <sample_name> -out_folder <output_folder> -vcf_pattern <vcf_file_pattern> -filter_pattern <filter_file_pattern> -reference <reference_genome> [-block_size <size>] [-output_distributions]\n";
    std::cout << "Version: " << VERSION_NUMBER << "\n\n";
    std::cout << "Parameters:\n";
    std::cout << "  -sample: Sample name to analyze\n";
    std::cout << "  -out_folder: Output folder for results. If not specified, results will be saved in the current directory\n";
    std::cout << "  -vcf_pattern: Path pattern to VCF files. Use {chr} placeholder which will be replaced with chromosome numbers (1-22)\n";
    std::cout << "                Example: /path/to/sample.{chr}.vcf.gz\n";
    std::cout << "  -filter_pattern: Path pattern to filter files. Use {chr} placeholder which will be replaced with chromosome numbers\n";
    std::cout << "                      Example: /path/to/filter.{chr}.bed\n";
    std::cout << "  -reference: Genome reference (hg19, hg38, or t2t)\n";
    std::cout << "  -block_size: Size in base pairs for heterozygosity distribution blocks (default: 10000)\n";
    std::cout << "  -output_distributions: Output additional distribution files (het_gq_dist.txt, alt_gq_dist.txt and het_dist_*.tab) for analysis\n";
}

int main(int argc, char *argv[])
{
    std::string sample_name = "";
    std::string out_folder = "";
    std::string vcf_file_pattern = "";
    std::string filter_file_pattern = "";
    std::string reference_genome = "";
    uint64_t het_block_size = 10000;
    bool output_dist = false;

    for (int i = 0; i < argc; i++)
    {
        if (std::string(argv[i]) == "-sample")
        {
            sample_name = argv[i + 1];
            if (sample_name.empty())
            {
                std::cerr << "Error: Sample name cannot be empty" << std::endl;
                return 1;
            }
        }
        else if (std::string(argv[i]) == "-out_folder")
        {
            out_folder = argv[i + 1];
            if (out_folder.back() != '/')
            {
                out_folder += "/";
            }
            if (!out_folder.empty() && !std::filesystem::exists(out_folder))
            {
                std::cerr << "Error: Output folder does not exist: " << out_folder << std::endl;
                return 1;
            }
        }
        else if (std::string(argv[i]) == "-vcf_pattern")
        {
            vcf_file_pattern = argv[i + 1];
            if (vcf_file_pattern.find("{chr}") == std::string::npos)
            {
                std::cerr << "Error: VCF file pattern must contain {chr} placeholder" << std::endl;
                return 1;
            }
        }
        else if (std::string(argv[i]) == "-filter_pattern")
        {
            filter_file_pattern = argv[i + 1];
            if (filter_file_pattern.find("{chr}") == std::string::npos)
            {
                std::cerr << "Error: Filter file pattern must contain {chr} placeholder" << std::endl;
                return 1;
            }
        }
        else if (std::string(argv[i]) == "-reference")
        {
            reference_genome = argv[i + 1];
            if (reference_genome != "hg19" && reference_genome != "hg38" && reference_genome != "t2t")
            {
                std::cerr << "Error: Reference genome must be one of: hg19, hg38, t2t" << std::endl;
                return 1;
            }
        }
        else if (std::string(argv[i]) == "-block_size")
        {
            het_block_size = std::stoull(argv[i + 1]);
        }
        else if (std::string(argv[i]) == "-output_distributions")
        {
            output_dist = true;
        }
        if ((std::string(argv[i]) == "-help") || (argc == 1))
        {
            print_help();
            return 0;
        }
    }

    if (sample_name.empty() || vcf_file_pattern.empty() || filter_file_pattern.empty() || reference_genome.empty())
    {
        std::cerr << "Error: Missing required parameters." << std::endl;
        print_help();
        return 1;
    }

    std::vector<std::string> chromosomes = {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"};
    std::vector<non_ref_stats_t> non_ref_stats(chromosomes.size());

    std::unordered_map<std::string, long> len_map;
    if (reference_genome == "hg19")
    {
        len_map = hg19_len_map;
    }
    else if (reference_genome == "hg38")
    {
        len_map = hg38_len_map;
    }
    else if (reference_genome == "t2t")
    {
        len_map = t2t_len_map;
    }

    std::unordered_map<uint64_t, uint64_t> het_gq_hist;
    std::unordered_map<uint64_t, uint64_t> alt_gq_hist;

    std::vector<std::thread> threads;
    for (int i = 0; i < chromosomes.size(); ++i)
    {
        auto chr = chromosomes[i];
        std::string vcf_file_temp = vcf_file_pattern;
        std::string filter_file_temp = filter_file_pattern;

        std::string vcf_file = vcf_file_temp.replace(vcf_file_temp.find("{chr}"), 5, chr);
        std::string filter_file = filter_file_temp.replace(filter_file_temp.find("{chr}"), 5, chr);
        threads.push_back(std::thread(get_non_ref_stats, vcf_file, filter_file, sample_name, chr, i, std::ref(non_ref_stats), std::ref(len_map), het_block_size));
    }

    for (auto& thread : threads)
    {
        thread.join();
    }

    std::ofstream res_file;
    std::ofstream het_gq_dist_file;
    std::ofstream alt_gq_dist_file;
    if (out_folder.empty())
    {
        res_file.open("non_ref_results.txt");
        if (output_dist)
        {
            het_gq_dist_file.open("het_gq_dist.txt");
            alt_gq_dist_file.open("alt_gq_dist.txt");
        }
    }
    else
    {
        res_file.open(out_folder + "/non_ref_results.txt");
        if (output_dist)
        {
            het_gq_dist_file.open(out_folder + "/het_gq_dist.txt");
            alt_gq_dist_file.open(out_folder + "/alt_gq_dist.txt");
        }
    }
    res_file << "chr\tn_bp\tn_hets\tn_hets_transversions\theterozygosity\thet_coverage\thet_gq_mean\tref_count\talt_count\talt_allele_prop\talt_gq_mean\tn_ts\tn_tv\tts_tv_ratio\n";
    for (const auto& result : non_ref_stats)
    {
        res_file << result.chr << "\t" << result.n_bp << "\t" << result.n_hets << "\t" << result.n_hets_transversions << "\t" << result.heterozygosity << "\t" << result.het_coverage << "\t" << result.het_gq_mean << "\t" << result.ref_count << "\t" << result.alt_count << "\t" << result.alt_allele_prop << "\t" << result.alt_gq_mean << "\t" << result.n_ts << "\t" << result.n_tv << "\t" << result.ts_tv_ratio << "\n";

        if (output_dist)
        {
            std::ofstream het_dist_file;
            if (out_folder.empty())
            {
                het_dist_file.open("het_dist_" + result.chr + "_" + std::to_string(het_block_size) + ".tab");
            }
            else
            {
                het_dist_file.open(out_folder + "/het_dist_" + result.chr + "_" + std::to_string(het_block_size) + ".tab");
            }
            for (size_t i = 0; i < result.het_dist.size(); ++i)
            {
                het_dist_file << i << "\t" << result.het_dist[i] << std::endl;
            }
            het_dist_file.close();
        }
        
        for (const auto& gq : result.het_gq_dist)
        {
            int bin = static_cast<int>(std::floor(gq));
            het_gq_hist[bin]++;
        }

        for (const auto& gq : result.alt_gq_dist)
        {
            int bin = static_cast<int>(std::floor(gq));
            alt_gq_hist[bin]++;
        }
    }
    res_file.close();

    if (output_dist)
    {
        het_gq_dist_file << "gq\tcount\n";
        for (const auto& [gq, count] : het_gq_hist)
        {
            het_gq_dist_file << gq << "\t" << count << "\n";
        }
        het_gq_dist_file.close();

        alt_gq_dist_file << "gq\tcount\n";
        for (const auto& [gq, count] : alt_gq_hist)
        {
            alt_gq_dist_file << gq << "\t" << count << "\n";
        }
        alt_gq_dist_file.close();
    }

    return 0;
}
