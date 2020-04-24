#include <boost/core/swap.hpp>
#include <cassert>
#include <stdexcept>

#include "lattice.h"

void lattice::lu(lattice::rmat_view_t m, std::vector<std::int64_t>& p) {
    std::int64_t size = m.shape()[0];

    assert(size == m.shape()[0]);
    assert(size == m.shape()[1]);

    p.resize(size);

    for (std::int64_t i = 0; i < size; ++i) {
        std::int64_t pivot = -1;

        for (std::int64_t row = i; row < size; ++row) {
            if (m[row][i].sign() != 0) {
                pivot = row;
                break;
            }
        }

        if (pivot == -1) {
            throw std::overflow_error("matrix is not invertible");
        }

        p[i] = pivot;

        if (pivot != i) {
            for (std::int64_t col = 0; col < size; ++col) {
                boost::multiprecision::swap(m[i][col], m[pivot][col]);
            }
        }

        for (std::int64_t row = i + 1; row < size; ++row) {
            m[row][i] /= m[i][i];
        }

        for (std::int64_t row = i + 1; row < size; ++row) {
            for (std::int64_t col = i + 1; col < size; ++col) {
                m[row][col] -= m[row][i] * m[i][col];
            }
        }
    }
}

void lattice::solve(rmat_view_t x, const lattice::rmat_view_t& m, const std::vector<std::int64_t>& p) {
    const std::int64_t size = m.shape()[0];

    assert(m.shape()[0] == size);
    assert(m.shape()[1] == size);
    assert(x.shape()[0] == size);
    assert(p.size() == size);

    for (std::int64_t dcol = 0; dcol < x.shape()[1]; ++dcol) {
        for (std::int64_t row = 0; row < size; ++row) {
            boost::multiprecision::swap(x[row][dcol], x[p[row]][dcol]);
        }
    }

    for (std::int64_t dcol = 0; dcol < x.shape()[1]; ++dcol) {
        for (std::int64_t row = 0; row < size; ++row) {
            for (std::int64_t col =0; col < row; ++col) {
                x[row][dcol] -= m[row][col] * x[col][dcol];
            }
        }
    }

    for (std::int64_t dcol = 0; dcol < x.shape()[1]; ++dcol) {
        for (std::int64_t row = size - 1; row >= 0; --row) {
            for (std::int64_t col = size - 1; col > row; --col) {
                x[row][dcol] -= m[row][col] * x[col][dcol];
            }

            x[row][dcol] /= m[row][row];
        }
    }
}
