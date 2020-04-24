#include <string>
#include <string_view>
#include <regex>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <vector>

#include "lattice.h"

static void parse_row(std::string_view line, lattice::rmat_t::array_view<1>::type row) {
    using svregex_iterator = std::regex_iterator<std::string_view::const_iterator>;

    std::regex entry_regex("[-+]?([0-9]+)(/([0-9]+))?");

    svregex_iterator begin(line.begin(), line.end(), entry_regex);
    svregex_iterator end;

    std::int64_t i = 0;
    svregex_iterator j = begin;

    for (; i != row.shape()[0] && j != end; ++i, ++j) {
        row[i] = lattice::rat_t(j->str());
    }

    if (i != row.shape()[0] || j != end) {
        throw std::invalid_argument("row has incorrect size");
    }
}

void lattice::parse(std::string_view name, std::int64_t& size, lattice::rmat_t& basis, lattice::rmat_t& lower, lattice::rmat_t& upper) {
    using range = lattice::rmat_t::index_range;

    std::ifstream fstream;
    std::istream& stream = name == "" ? std::cin : fstream;

    if (name != "") {
        fstream.open(name);
    }

    if (!stream) {
        throw std::runtime_error("could not open file");
    }

    std::string line;
    std::getline(stream, line);

    size = std::strtoll(line.data(), nullptr, 10);

    if (size < 1) {
        throw std::runtime_error("invalid size");
    }

    basis.resize(boost::extents[size][size]);
    lower.resize(boost::extents[size][1]);
    upper.resize(boost::extents[size][1]);

    std::getline(stream, line);

    for (std::int64_t i = 0; i < size; ++i) {
        std::getline(stream, line);
        ::parse_row(line, basis[boost::indices[range()][i]]);
    }

    std::getline(stream, line);

    std::getline(stream, line);
    ::parse_row(line, lower[boost::indices[range()][0]]);

    std::getline(stream, line);
    ::parse_row(line, upper[boost::indices[range()][0]]);
}
