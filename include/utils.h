#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <chrono>
#include <ctime>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <algorithm>
#include <unordered_map>
#include <cmath>
#include <utility>
#include <iomanip>

std::pair<double, double> binomial_confidence_interval(
    uint64_t successes,
    uint64_t trials)
{
    if (trials == 0)
        return {0.0, 0.0};

    double p = static_cast<double>(successes) / trials;
    double z = 1.96;

    double se = std::sqrt(p * (1 - p) / trials);
    double lower = p - z * se;
    double upper = p + z * se;

    if (lower < 0)
        lower = 0;
    if (upper > 1)
        upper = 1;

    return {lower, upper};
}

std::string binomial_ci_string(std::pair<double, double> ci)
{
    std::ostringstream ci_s;
    ci_s << std::fixed << std::setprecision(2) << ci.first * 100 << "-" << ci.second * 100;

    return ci_s.str();
}

std::string extract_file_name(const std::string &path)
{
    size_t pos_start = path.find_last_of("/") + 1;
    size_t pos_end = path.find_last_of(".");

    if (pos_end == std::string::npos || pos_end < pos_start)
    {
        pos_end = path.length();
    }

    return path.substr(pos_start, pos_end - pos_start);
}

std::string get_date_time()
{
    auto now = std::chrono::system_clock::now();
    std::time_t time = std::chrono::system_clock::to_time_t(now);
    std::string date_time = std::ctime(&time);
    date_time.pop_back();
    return date_time;
}

std::string replace_substring(
    const std::string &input,
    const std::string &old_str,
    const std::string &new_str)
{
    std::string result = input;
    size_t pos = 0;
    while ((pos = result.find(old_str, pos)) != std::string::npos)
    {
        result.replace(pos, old_str.length(), new_str);
        pos += new_str.length();
    }
    return result;
}

auto count_tabs = [](const std::string &s)
{
    return std::count(s.begin(), s.end(), '\t');
};

std::string reverse_complement(const std::string &sequence)
{
    std::unordered_map<char, char> complement = {
        {'A', 'T'}, {'T', 'A'},
        {'C', 'G'}, {'G', 'C'},
        {'a', 't'}, {'t', 'a'},
        {'c', 'g'}, {'g', 'c'},
        {'N', 'N'}, {'n', 'n'}};

    std::string rev_comp(sequence.rbegin(), sequence.rend());
    std::transform(rev_comp.begin(), rev_comp.end(), rev_comp.begin(),
                   [&complement](char nucleotide)
                   { return complement[nucleotide]; });

    return rev_comp;
}

std::vector<std::string> split_by_numeric(const std::string &str)
{
    std::vector<std::string> result;
    std::string current_segment;
    bool is_current_segment_numeric = isdigit(str[0]);

    for (char ch : str)
    {
        if (isdigit(ch) == is_current_segment_numeric)
        {
            current_segment += ch;
        }
        else
        {
            result.push_back(current_segment);
            current_segment = ch;
            is_current_segment_numeric = !is_current_segment_numeric;
        }
    }

    if (!current_segment.empty())
    {
        result.push_back(current_segment);
    }

    return result;
}

#endif // UTILS_H