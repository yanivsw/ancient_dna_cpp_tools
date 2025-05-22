#ifndef BED_PROCESSING_H
#define BED_PROCESSING_H

#include <iostream>
#include <variant>
#include <htslib/sam.h>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include <unordered_map>
#include <zlib.h>

struct bed_region_t
{
    uint64_t start;
    uint64_t end;
    uint64_t size = 0;
};

using bed_file_t = std::unordered_map<std::string, std::vector<bed_region_t>>;

void process_bed_line(
    std::string line,
    bed_file_t &bed_file_map,
    bool get_region_size)
{
    std::istringstream iss(line);
    std::string chromosome;
    int start, end;
    if (!(iss >> chromosome >> start >> end))
    {
        std::cerr << "Failed to parse line: " << line << std::endl;
        // continue;
    }

    bed_region_t region;
    region.start = start;
    region.end = end;
    if (get_region_size)
    {
        region.size = end - start;
    }

    bed_file_map[chromosome].push_back(region);
}

void read_gzipped_file(
    std::string bed_file_location,
    bed_file_t &bed_file_map,
    bool get_region_size)
{
    std::string line;
    gzFile bed_file = gzopen(bed_file_location.c_str(), "rb");

    char gz_line[1024];
    while (gzgets(bed_file, gz_line, sizeof(gz_line)))
    {
        line = gz_line;
        process_bed_line(line, bed_file_map, get_region_size);
    }

    gzclose(bed_file);
}

void read_text_file(
    std::string bed_file_location,
    bed_file_t &bed_file_map,
    bool get_region_size)
{
    std::string line;
    std::ifstream bed_file(bed_file_location);

    while (std::getline(bed_file, line))
    {
        process_bed_line(line, bed_file_map, get_region_size);
    }

    bed_file.close();
}

int read_bed_file(
    std::string bed_file_location,
    bed_file_t &bed_file_map,
    bool get_region_size = false)
{
    bool is_gzipped = false;

    std::cerr << "Loading target file: " << bed_file_location << std::endl;

    std::ifstream bed_file(bed_file_location, std::ios::binary);
    if (!bed_file)
    {
        std::cerr << "Failed to open BED file: " << bed_file_location << std::endl;
        return 1;
    }

    // Check if the file is gzipped
    char buffer[2];
    bed_file.read(buffer, 2);
    bed_file.close();

    if (buffer[0] == '\x1F' && buffer[1] == '\x8B')
    {
        is_gzipped = true;
    }

    if (is_gzipped)
    {
        read_gzipped_file(bed_file_location, bed_file_map, get_region_size);
    }
    else
    {
        read_text_file(bed_file_location, bed_file_map, get_region_size);
    }

    return 0;
}

uint64_t get_bed_size(std::string bed_file_location,
                      std::string chr = "")
{
    bed_file_t bed_file_map;
    read_bed_file(bed_file_location, bed_file_map, true);

    uint64_t bed_size = 0;
    if (!chr.empty())
    {
        std::vector<bed_region_t> regions = bed_file_map[chr];
        for (const auto &region : regions)
        {
            bed_size += region.size;
        }
    }
    else
    {
        for (const auto &chr_map : bed_file_map)
        {
            const std::vector<bed_region_t> &regions = chr_map.second;
            for (const auto &region : regions)
            {
                bed_size += region.size;
            }
        }
    }

    return bed_size;
}

#endif // BED_PROCESSING_H