#pragma once

#include <boost/multi_array.hpp>
#include <boost/multiprecision/gmp.hpp>
#include <string_view>
#include <vector>

#include "fmt/format.h"

namespace lattice {
    typedef boost::multiprecision::mpz_int int_t;
    typedef boost::multiprecision::mpq_rational rat_t;

    typedef boost::multi_array<rat_t, 2> rmat_t;
    typedef boost::multi_array<rat_t, 2>::array_view<2>::type rmat_view_t;

    void parse(std::string_view name, std::int64_t& size, rmat_t& basis, rmat_t& lower, rmat_t& upper);
    void lu(rmat_view_t m, std::vector<std::int64_t>& p);
    void solve(rmat_view_t x, const rmat_view_t& m, const std::vector<std::int64_t>& p);
    void search(const rmat_t& basis, const rmat_t& lower, const rmat_t& upper, std::vector<rmat_t>& results, std::int64_t thread_max);
    void optimize(rmat_view_t x, const rmat_view_t& initial_table, std::int64_t dimensions, std::int64_t depth);

    inline rmat_view_t view(rmat_t& m) {
        return m[boost::indices[rmat_t::index_range()][rmat_t::index_range()]];
    }

    inline rmat_view_t view(rmat_t& m, std::int64_t row0, std::int64_t row1, std::int64_t col0, std::int64_t col1) {
        return m[boost::indices[rmat_t::index_range(row0, row1)][rmat_t::index_range(col0, col1)]];
    }

    inline int_t floor(const rat_t& x) {
        const int_t& n = boost::multiprecision::numerator(x);
        const int_t& m = boost::multiprecision::denominator(x);

        return n / m - (x.sign() == -1);
    }

    inline int_t ceil(const rat_t& x) {
        return -lattice::floor(-x);
    }

    inline void print(const rmat_view_t& view) {
        for (std::int64_t i = 0; i < view.shape()[0]; ++i) {
            for (std::int64_t j = 0; j < view.shape()[1]; ++j) {
                fmt::print("{} ", view[i][j].str());
            }

            fmt::print("\n");
        }
    }

    inline void print(const rmat_t& mat) {
        for (std::int64_t i = 0; i < mat.shape()[0]; ++i) {
            for (std::int64_t j = 0; j < mat.shape()[1]; ++j) {
                fmt::print("{} ", mat[i][j].str());
            }

            fmt::print("\n");
        }
    }
}
