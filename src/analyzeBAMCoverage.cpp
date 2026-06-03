#include <htslib/sam.h>
#include <iostream>
#include <vector>
#include <numeric>
#include <thread>
#include <limits>
#include <fstream>
#include <unordered_map>
#include <sstream>
#include <zlib.h>
#include <mutex>
#include <atomic>
#include <algorithm>
#include <cctype>
#include <array>

#include "bam_processing.h"
#include "bed_processing.h"
#include "chromosome_sizes.h"
#include "types.h"
#include "utils.h"

#define COVERAGE_DISTRIBUTION_SIZE 100000

static constexpr uint8_t GC_BIN_MISSING = 255;
static constexpr size_t GC_BINS = 12;
static constexpr uint32_t GC_BIN_MISSING_OUT = 11;

struct readlen_bins_t
{
    std::vector<int> edges; // size = nbins+1; last edge may be INT_MAX for inf

    size_t nbins() const
    {
        return edges.size() > 1 ? edges.size() - 1 : 0;
    }

    int bin_of(int readlen) const
    {
        if (edges.size() < 2) return -1;

        // bins are [edges[i], edges[i+1])
        for (size_t i = 0; i + 1 < edges.size(); ++i)
        {
            if (readlen >= edges[i] && readlen < edges[i + 1]) return static_cast<int>(i);
        }
        return -1;
    }

    std::string label(size_t i) const
    {
        if (edges.size() < 2 || i >= nbins())
        {
            return "NA";
        }

        const int lo = edges[i];
        const int hi = edges[i + 1];
        if (hi == std::numeric_limits<int>::max())
        {
            return std::to_string(lo) + "_max";
        }
        return std::to_string(lo) + "_" + std::to_string(hi);
    }
};

static inline readlen_bins_t parse_readlen_bins(const std::string& csv)
{
    // Examples:
    //  "35,60,100,141"      -> edges {35,60,100,141,INT_MAX}
    //  "35,60,100,141,inf"  -> edges {35,60,100,141,INT_MAX}
    readlen_bins_t out;
    std::vector<int> edges;

    size_t start = 0;
    while (start < csv.size())
    {
        size_t comma = csv.find(',', start);
        if (comma == std::string::npos)
        {
            comma = csv.size();
        }
        std::string tok = trim_copy(csv.substr(start, comma - start));
        start = comma + 1;

        if (tok.empty())
        {
            continue;
        }

        std::string lowtok;
        lowtok.reserve(tok.size());
        for (char c : tok)
        {
            lowtok.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
        }

        if (lowtok == "inf" || lowtok == "infinity")
        {
            edges.push_back(std::numeric_limits<int>::max());
            continue;
        }

        char* endp = nullptr;
        long v = std::strtol(tok.c_str(), &endp, 10);
        if (endp == tok.c_str() || *endp != '\0')
        {
            throw std::runtime_error("Bad readlen bin token: '" + tok + "'");
        }
        if (v < 0 || v > std::numeric_limits<int>::max())
        {
            throw std::runtime_error("Out-of-range readlen bin edge: '" + tok + "'");
        }
        edges.push_back(static_cast<int>(v));
    }

    if (edges.size() < 2)
    {
        throw std::runtime_error("Need at least 2 edges for readlen bins (e.g. 35,60)");
    }

    // Extend last edge to infinity automatically
    if (edges.back() != std::numeric_limits<int>::max())
    {
        edges.push_back(std::numeric_limits<int>::max());
    }

    // Validate strictly increasing
    for (size_t i = 1; i < edges.size(); ++i)
    {
        if (!(edges[i] > edges[i - 1]))
        {
            throw std::runtime_error("Readlen bin edges must be strictly increasing");
        }
    }

    out.edges = std::move(edges);
    return out;
}

static readlen_bins_t readlen_bins_global;

static inline int readlen_to_bin(int readlen)
{
    return readlen_bins_global.bin_of(readlen);
}

static inline std::string readlen_bin_label(size_t bin)
{
    return readlen_bins_global.label(bin);
}

static inline size_t gc_flat_index(size_t cov, size_t gc_bin)
{
    return cov * GC_BINS + gc_bin;
}

struct chrom_accum_t
{
    std::vector<uint64_t> cov_total;
    std::vector<std::vector<uint64_t>> cov_total_by_bin;

    std::vector<uint64_t> cov_gc_total_flat;
    std::vector<std::vector<uint64_t>> cov_gc_total_by_bin_flat;

    std::vector<uint64_t> base_cov_by_lenbin;

    void init(size_t cov_size, size_t nbins, bool compute_gc)
    {
        cov_total.assign(cov_size, 0);
        cov_total_by_bin.assign(nbins, std::vector<uint64_t>(cov_size, 0));

        if (compute_gc)
        {
            cov_gc_total_flat.assign(cov_size * GC_BINS, 0);
            cov_gc_total_by_bin_flat.assign(nbins, std::vector<uint64_t>(cov_size * GC_BINS, 0));
        }

        base_cov_by_lenbin.assign(nbins, 0);
    }
};

static inline void add_chrom_accum(chrom_accum_t& dst, const chrom_accum_t& src)
{
    const size_t cov_size = dst.cov_total.size();
    const size_t nbins    = dst.cov_total_by_bin.size();
    const bool   has_gc   = !dst.cov_gc_total_flat.empty();

    if (src.cov_total.size() != cov_size || src.cov_total_by_bin.size() != nbins)
    {
        std::cerr << "add_chrom_accum: size mismatch\n";
        std::abort();
    }

    // non-GC totals
    for (size_t i = 0; i < cov_size; ++i)
    {
        dst.cov_total[i] += src.cov_total[i];
    }

    // GC totals
    if (has_gc)
    {
        if (src.cov_gc_total_flat.size() != cov_size * GC_BINS)
        {
            std::cerr << "add_chrom_accum: GC flat size mismatch\n";
            std::abort();
        }
        for (size_t i = 0; i < cov_size * GC_BINS; ++i)
        {
            dst.cov_gc_total_flat[i] += src.cov_gc_total_flat[i];
        }
    }

    for (size_t length_bin = 0; length_bin < nbins; ++length_bin)
    {
        if (dst.cov_total_by_bin[length_bin].size() != cov_size || src.cov_total_by_bin[length_bin].size() != cov_size)
        {
            std::cerr << "add_chrom_accum: bin size mismatch at b=" << length_bin << "\n";
            std::abort();
        }

        for (size_t i = 0; i < cov_size; ++i)
        {
            dst.cov_total_by_bin[length_bin][i] += src.cov_total_by_bin[length_bin][i];
        }

        if (has_gc)
        {
            for (size_t i = 0; i < cov_size * GC_BINS; ++i)
            {
                dst.cov_gc_total_by_bin_flat[length_bin][i] += src.cov_gc_total_by_bin_flat[length_bin][i];
            }
        }

        dst.base_cov_by_lenbin[length_bin] += src.base_cov_by_lenbin[length_bin];
    }
}

static inline void add_zero_coverage_into_histograms(
    const std::vector<bool>& bed_mask,
    const std::vector<uint8_t>& gc_bin_by_pos,
    chrom_accum_t& out)
{
    const size_t cov_size = out.cov_total.size();
    const bool   has_gc   = !out.cov_gc_total_flat.empty();
    if (cov_size == 0)
    {
        return;
    }

    uint64_t N_bed = 0;
    std::array<uint64_t, GC_BINS> N_by_gc{};
    N_by_gc.fill(0);

    for (size_t pos0 = 0; pos0 < bed_mask.size(); ++pos0)
    {
        if (!bed_mask[pos0])
        {
            continue;
        }
        ++N_bed;

        if (has_gc)
        {
            const uint8_t gc_u8 = gc_bin_by_pos[pos0];
            uint32_t gc_bin = (gc_u8 == GC_BIN_MISSING) ? GC_BIN_MISSING_OUT : std::min<uint32_t>(gc_u8, 10);
            ++N_by_gc[gc_bin];
        }
    }

    // total cov==0
    uint64_t nonzero_total = 0;
    for (size_t d = 1; d < cov_size; ++d)
    {
        nonzero_total += out.cov_total[d];
    }
    out.cov_total[0] = (N_bed >= nonzero_total) ? (N_bed - nonzero_total) : 0;

    // GC cov==0
    if (has_gc)
    {
        for (size_t gc = 0; gc < GC_BINS; ++gc)
        {
            uint64_t nonzero_gc = 0;
            for (size_t d = 1; d < cov_size; ++d)
            {
                nonzero_gc += out.cov_gc_total_flat[gc_flat_index(d, gc)];
            }
            out.cov_gc_total_flat[gc_flat_index(0, gc)] = (N_by_gc[gc] >= nonzero_gc) ? (N_by_gc[gc] - nonzero_gc) : 0;
        }
    }

    for (int length_bin = 0; length_bin < (int)readlen_bins_global.nbins(); ++length_bin)
    {
        uint64_t nonzero_b = 0;
        for (size_t d = 1; d < cov_size; ++d)
        {
            nonzero_b += out.cov_total_by_bin[length_bin][d];
        }
        out.cov_total_by_bin[length_bin][0] = (N_bed >= nonzero_b) ? (N_bed - nonzero_b) : 0;

        if (has_gc)
        {
            for (size_t gc = 0; gc < GC_BINS; ++gc)
            {
                uint64_t nonzero_gc_b = 0;
                for (size_t d = 1; d < cov_size; ++d)
                {
                    nonzero_gc_b += out.cov_gc_total_by_bin_flat[length_bin][gc_flat_index(d, gc)];
                }
                out.cov_gc_total_by_bin_flat[length_bin][gc_flat_index(0, gc)] = (N_by_gc[gc] >= nonzero_gc_b) ? (N_by_gc[gc] - nonzero_gc_b) : 0;
            }
        }
    }
}

static int get_coverage_histograms_region(
    bam_file_config_t& bam_config,
    int tid,
    int64_t beg,
    int64_t end,
    const std::vector<bool>& bed_mask,
    const std::vector<uint8_t>& gc_bin_by_pos,
    chrom_accum_t& out)
{
    const int64_t chr_len = static_cast<int64_t>(bam_config.header->target_len[tid]);
    if (beg < 0)
    {
        beg = 0;
    }
    if (end > chr_len)
    {
        end = chr_len;
    }
    if (end <= beg)
    {
        return 0;
    }

    const int ibeg = static_cast<int>(beg);
    const int iend = static_cast<int>(end);

    hts_itr_t* iter = sam_itr_queryi(bam_config.index, tid, ibeg, iend);
    if (!iter)
    {
        std::cerr << "Failed to create iterator for tid=" << tid << " beg=" << beg << " end=" << end << "\n";
        return 1;
    }

    bam_plp_t plp = bam_plp_init(nullptr, nullptr); // manual push mode
    if (!plp)
    {
        hts_itr_destroy(iter);
        std::cerr << "Failed to init pileup\n";
        return 1;
    }
    bam_plp_set_maxcnt(plp, std::numeric_limits<int>::max());

    bam1_t* b = bam_init1();
    if (!b)
    {
        bam_plp_destroy(plp);
        hts_itr_destroy(iter);
        return 1;
    }

    const bam_pileup1_t* pileup = nullptr;
    int plp_tid = -1, refpos = -1, n = 0;
    const size_t cov_size = out.cov_total.size();

    auto drain = [&]() {
        while ((pileup = bam_plp_next(plp, &plp_tid, &refpos, &n)) != nullptr)
        {
            if (plp_tid != tid)
            {
                continue;
            }
            if (refpos < ibeg || refpos >= iend)
            {
                continue;
            }

            const size_t pos0 = static_cast<size_t>(refpos);
            if (pos0 >= bed_mask.size() || !bed_mask[pos0])
            {
                continue;
            }

            uint32_t gc_bin = 0;
            const uint8_t gc_u8 = gc_bin_by_pos[pos0];
            if (gc_u8 == GC_BIN_MISSING)
            {
                gc_bin = GC_BIN_MISSING_OUT;
            }
            else
            {
                gc_bin = static_cast<uint32_t>(gc_u8);
                if (gc_bin > 10)
                {
                    gc_bin = 10;
                }
            }

            uint32_t depth_total = 0;
            uint32_t nbins = readlen_bins_global.nbins();
            std::vector<uint32_t> depth_by_bin(nbins, 0);

            for (int i = 0; i < n; ++i)
            {
                if (pileup[i].is_del || pileup[i].is_refskip)
                {
                    continue;
                }
                const bam1_t* alignment = pileup[i].b;
                if (!alignment)
                {
                    continue;
                }
                // Don't count pairs for now (may be added in the future as an option)
                if (alignment->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FPAIRED))
                {
                    continue;
                }
                depth_total++;

                const int bin = readlen_to_bin(alignment->core.l_qseq);
                if (bin >= 0)
                {
                    depth_by_bin[bin]++;
                    out.base_cov_by_lenbin[bin]++;
                }
            }

            if (depth_total < cov_size)
            {
                out.cov_total[depth_total]++;
                if (!out.cov_gc_total_flat.empty())
                {
                    out.cov_gc_total_flat[gc_flat_index(depth_total, gc_bin)]++;
                }
            }

            for (int length_bin = 0; length_bin < (int)readlen_bins_global.nbins(); ++length_bin)
            {
                const uint32_t d = depth_by_bin[length_bin];
                if (d == 0)
                {
                    continue;
                }
                if (d < cov_size)
                {
                    out.cov_total_by_bin[length_bin][d]++;
                    if (!out.cov_gc_total_by_bin_flat.empty())
                    {
                        out.cov_gc_total_by_bin_flat[length_bin][gc_flat_index(d, gc_bin)]++;
                    }
                }
            }
        }
    };

    // Stream alignments into pileup and drain as we go
    while (sam_itr_next(bam_config.bam_file, iter, b) >= 0)
    {
        bam_plp_push(plp, b);
        drain();
    }

    bam_plp_push(plp, nullptr); // signal end
    drain();

    bam_destroy1(b);
    bam_plp_destroy(plp);
    hts_itr_destroy(iter);
    return 0;
}

static int build_bed_mask_and_gc_bins(
    const std::string& chr,
    long chromosome_length,
    const std::string& gc51_file,
    bed_file_t& map35,
    std::vector<bool>& bed_mask,
    std::vector<uint8_t>& gc_bin_by_pos,
    bool compute_gc)
{
    bed_mask.assign(static_cast<size_t>(chromosome_length), false);
    gc_bin_by_pos.assign(static_cast<size_t>(chromosome_length), GC_BIN_MISSING);

    // mark mask from BED
    const auto& regions = map35[chr];
    for (const auto& region : regions)
    {
        const uint64_t start = std::min<uint64_t>(region.start, bed_mask.size());
        const uint64_t end = std::min<uint64_t>(region.end, bed_mask.size());
        for (uint64_t i = start; i < end; ++i)
        {
            bed_mask[i] = true;
        }
    }

    // read GC file and store GC-bin for masked positions
    if (!compute_gc)
    {
        return 0;
    }
    gzFile file = gzopen(gc51_file.c_str(), "rb");
    if (!file)
    {
        std::cerr << "Failed to open GC file: " << gc51_file << "\n";
        return 1;
    }

    char buffer[1024];
    while (gzgets(file, buffer, sizeof(buffer)) != NULL)
    {
        std::istringstream iss(buffer);
        std::string temp;
        int key = 0, value = 0;

        // expected: <something> <pos> <gc>
        while (iss >> temp >> key >> value)
        {
            if (key <= 0)
            {
                continue;
            }
            const size_t pos0 = static_cast<size_t>(key - 1);
            if (pos0 >= bed_mask.size())
            {
                continue;
            }
            if (!bed_mask[pos0])
            {
                continue;
            }

            uint32_t gc_bin = static_cast<uint32_t>(value / 5);
            if (gc_bin > 10)
            {
                gc_bin = 10;
            }
            gc_bin_by_pos[pos0] = static_cast<uint8_t>(gc_bin);
        }
    }
    gzclose(file);

    return 0;
}

static int process_one_chromosome_in_threads(
    const std::string& bam_file,
    const std::string& gc51_file,
    const std::string& chr,
    long chromosome_length,
    bed_file_t& map35,
    size_t n_threads,
    bool compute_gc,
    chrom_accum_t& chrom_out)
{
    std::vector<bool> bed_mask;
    std::vector<uint8_t> gc_bin_by_pos;
    if (build_bed_mask_and_gc_bins(chr, chromosome_length, gc51_file, map35, bed_mask, gc_bin_by_pos, compute_gc) != 0)
    {
        return 1;
    }

    const int64_t chunk_size = 1000000;
    std::vector<std::pair<int64_t, int64_t>> regions;
    for (int64_t chunk_start = 0; chunk_start < chromosome_length; chunk_start += chunk_size)
    {
        const int64_t chunk_end = std::min<int64_t>(chunk_start + chunk_size, chromosome_length);
        regions.emplace_back(chunk_start, chunk_end);
    }

    std::cout << "Processing " << chr << " in " << regions.size() << " regions with " << n_threads << " threads...\n";
    std::atomic<size_t> next_region{0};
    std::vector<chrom_accum_t> local_accumulators(n_threads);
    for (auto& accumulator : local_accumulators)
    {
        accumulator.init(chrom_out.cov_total.size(), readlen_bins_global.nbins(), compute_gc);
    }

    std::vector<std::thread> threads;
    threads.reserve(n_threads);

    for (size_t t = 0; t < n_threads; ++t)
    {
        threads.emplace_back([&, t]() {
            bam_file_config_t bam_config = {};
            bam_constructor(bam_file, &bam_config);

            const int tid = bam_name2id(bam_config.header, chr.c_str());
            if (tid < 0)
            {
                bam_destructor(&bam_config);
                return;
            }

            while (true)
            {
                const size_t r = next_region.fetch_add(1);
                if (r >= regions.size())
                {
                    break;
                }

                const auto [beg, end] = regions[r];
                (void)get_coverage_histograms_region(bam_config, tid,beg, end,
                                                     bed_mask, gc_bin_by_pos,
                                                     local_accumulators[t]);
                
                std::cout << "Thread " << t << " done with region " << r << " (" << beg << "-" << end << ")\n";
            }

            bam_destructor(&bam_config);
        });
    }

    for (auto& th : threads)
    {
        th.join();
    }

    // merge thread-local accumulators into chrom_out
    for (const auto& accumulator : local_accumulators)
    {
        add_chrom_accum(chrom_out, accumulator);
    }

    add_zero_coverage_into_histograms(bed_mask, gc_bin_by_pos, chrom_out);
    return 0;
}

void print_help(const char* prog)
{
    std::cout << "Usage: " << prog << " [options]\n"
              << "Required:\n"
              << "  -b   <path>        Input BAM file\n"
              << "  -o   <path>        Output directory\n"
              << "  -m   <path>        Mappability BED file\n"
              << "  -r   <str>         Read length bins CSV e.g. 35,60,100,142\n"
              << "  -ref <str>         Reference genome: t2t | hg38 | hg19\n"
              << "Optional:\n"
              << "  -gc  <path>        GC bin directory (enables GC computation)\n"
              << "                     Expected files: <dir>/chr<N>.tab.gz\n"
              << "  -chr <c1,c2,...>   Chromosomes to process [default: 1-22]\n"
              << "  -t   <int>         Threads per chromosome [default: 40]\n"
              << "  -h                 Print this help\n";
}

static void write_cov_table(
    const std::string& path,
    const std::vector<uint64_t>& cov_total)
{
    std::ofstream f(path);
    f << "Cov\tTotal\n";
    for (size_t cov = 0; cov < cov_total.size(); ++cov)
    {
        if (cov_total[cov] > 0) f << cov << "\t" << cov_total[cov] << "\n";
    }
}

static void write_gc_cov_table_flat(
    const std::string& path,
    const std::vector<uint64_t>& flat,
    size_t cov_size)
{
    std::ofstream f(path);

    f << "Cov";
    for (size_t gc = 0; gc < GC_BINS; ++gc)
    {
        if (gc == GC_BIN_MISSING_OUT)
        {
            f << "\tGC_missing";
        }
        else
        {
            f << "\t" << gc;
        }
    }
    f << "\n";

    for (size_t cov = 0; cov < cov_size; ++cov) // include cov==0
    {
        uint64_t row_sum = 0;
        for (size_t gc = 0; gc < GC_BINS; ++gc)
        {
            row_sum += flat[gc_flat_index(cov, gc)];
        }
        if (row_sum == 0)
        {
            continue;
        }

        f << cov;
        for (size_t gc = 0; gc < GC_BINS; ++gc)
        {
            f << "\t" << flat[gc_flat_index(cov, gc)];
        }
        f << "\n";
    }
}

int main(int argc, char* argv[])
{
    std::string bam_file;
    std::string output_file_location;
    std::string bed_file;
    std::string bins_csv;
    std::string gc_dir;
    std::string chromosomes_str;
    std::string ref_name;
    size_t n_threads_per_chr = 40;

    if (argc == 1)
    {
        print_help(argv[0]);
        return 0;
    }

    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];
        if (arg == "-b" && i+1 < argc)
        {
            bam_file = argv[++i];
        }
        else if (arg == "-o" && i+1 < argc)
        {
            output_file_location = argv[++i];
        }
        else if (arg == "-m" && i+1 < argc)
        {
            bed_file = argv[++i];
        }
        else if (arg == "-r" && i+1 < argc)
        {
            bins_csv = argv[++i];
        }
        else if (arg == "-gc" && i+1 < argc)
        {
            gc_dir = argv[++i];
        }
        else if (arg == "-chr" && i+1 < argc)
        {
            chromosomes_str = argv[++i];
        }
        else if (arg == "-t" && i+1 < argc)
        {
            n_threads_per_chr = std::stoul(argv[++i]);
        }
        else if (arg == "-ref" && i+1 < argc)
        {
            ref_name = argv[++i];
        }
        else if (arg == "-h")
        {
            print_help(argv[0]);
            return 0;
        }
        else
        {
            std::cerr << "Unknown argument: " << arg << "\n";
            print_help(argv[0]);
            return 1;
        }
    }

    // Select length map
    std::unordered_map<std::string, long> len_map;
    if (ref_name == "t2t"  || ref_name == "T2T")
    {
        len_map = t2t_len_map;
    }
    else if (ref_name == "hg38" || ref_name == "HG38")
    {
        len_map = hg38_len_map;
    }
    else if (ref_name == "hg19" || ref_name == "HG19")
    {
        len_map = hg19_len_map;
    }
    else
    {
        std::cerr << "Error: unknown reference '" << ref_name << "'. " << "Choose from: t2t, hg38, hg19\n";
        return 1;
    }

    const bool compute_gc = !gc_dir.empty();

    std::vector<std::string> chromosomes;
    if (chromosomes_str.empty())
    {
        chromosomes = {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"};
    }
    else
    {
        std::stringstream ss(chromosomes_str);
        std::string tok;
        while (std::getline(ss, tok, ','))
        {
            if (!tok.empty())
            {
                chromosomes.push_back(tok);
            }
        }
    }

    readlen_bins_global = parse_readlen_bins(bins_csv);
    const size_t nbins = readlen_bins_global.nbins();

    chrom_accum_t global_accumulator;
    std::vector<chrom_accum_t> chromosome_accumulators(chromosomes.size());

    global_accumulator.init(COVERAGE_DISTRIBUTION_SIZE, nbins, compute_gc);
    for (auto& accumulator : chromosome_accumulators)
    {
        accumulator.init(COVERAGE_DISTRIBUTION_SIZE, nbins, compute_gc);
    }
    std::string input_bam_name = extract_file_name(bam_file);

    bed_file_t map35;
    read_bed_file(bed_file, map35);
    
    for (size_t chr_i = 0; chr_i < chromosomes.size(); ++chr_i)
    {
        const std::string& chr = chromosomes[chr_i];
        const std::string gc51_file = compute_gc ? gc_dir + "/chr" + chr + ".tab.gz" : "";
        const long chr_len = len_map[chr];

        const int rc = process_one_chromosome_in_threads(bam_file,
                                                         gc51_file,
                                                         chr,
                                                         chr_len,
                                                         map35,
                                                         n_threads_per_chr,
                                                         compute_gc,
                                                         chromosome_accumulators[chr_i]);
        if (rc != 0)
        {
            std::cerr << "Failed chromosome " << chr << "\n";
            return 1;
        }
    }

    for (const auto& accumulator : chromosome_accumulators)
    {
        add_chrom_accum(global_accumulator, accumulator);
    }
    
    for (size_t chr_i = 0; chr_i < chromosomes.size(); ++chr_i)
    {
        const auto& c = chromosome_accumulators[chr_i];
        const std::string chr_label = chromosomes[chr_i];

        std::string output_prefix_cov = output_file_location + "/coverage." + input_bam_name + ".chr" + chr_label;
        std::string output_prefix_gc = output_file_location + "/gc." + input_bam_name + ".chr" + chr_label;

        write_cov_table(output_prefix_cov + ".tab", c.cov_total);
        if (compute_gc)
        {
            write_gc_cov_table_flat(output_prefix_gc + ".tab", c.cov_gc_total_flat, COVERAGE_DISTRIBUTION_SIZE);
        }

        for (int bin = 0; bin < (int)readlen_bins_global.nbins(); ++bin)
        {
            write_cov_table(output_prefix_cov + "." + readlen_bin_label(bin) + ".tab", c.cov_total_by_bin[bin]);
            if (compute_gc)
            {
                write_gc_cov_table_flat(output_prefix_gc + "." + readlen_bin_label(bin) + ".tab", c.cov_gc_total_by_bin_flat[bin], COVERAGE_DISTRIBUTION_SIZE);
            }
        }
    }

    // ---- global outputs ----
    std::string output_prefix_cov = output_file_location + "/coverage." + input_bam_name;
    std::string output_prefix_gc = output_file_location + "/gc." + input_bam_name;
    write_cov_table(output_prefix_cov + ".tab", global_accumulator.cov_total);
    if (compute_gc)
    {
        write_gc_cov_table_flat(output_prefix_gc + ".tab", global_accumulator.cov_gc_total_flat, COVERAGE_DISTRIBUTION_SIZE);
    }
    for (int bin = 0; bin < (int)readlen_bins_global.nbins(); ++bin)
    {
        write_cov_table(output_prefix_cov + "." + readlen_bin_label(bin) + ".tab", global_accumulator.cov_total_by_bin[bin]);
        if (compute_gc)
        {
            write_gc_cov_table_flat(output_prefix_gc + "." + readlen_bin_label(bin) + ".tab", global_accumulator.cov_gc_total_by_bin_flat[bin], COVERAGE_DISTRIBUTION_SIZE);
        }
    }

    return 0;
}
