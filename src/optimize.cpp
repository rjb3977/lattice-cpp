#include <cassert>
#include <condition_variable>
#include <thread>
#include <mutex>
#include <vector>

#include "lattice.h"
#include "fmt/format.h"

static void pivot(lattice::rmat_view_t table, std::vector<std::int64_t>& B, std::vector<std::int64_t>& N, std::int64_t variables, std::int64_t constraints, std::int64_t entering, std::int64_t exiting) {
    std::int64_t a = table.shape()[0] - constraints;
    std::int64_t b = table.shape()[1] - 1;

    assert(0 <= entering && entering < b);
    assert(0 <= exiting && exiting < constraints);

    for (std::int64_t col = 0; col < table.shape()[1]; ++col) {
        if (col == entering) {
            continue;
        }

        table[a + exiting][col] /= table[a + exiting][entering];
    }

    for (std::int64_t row = 0; row < table.shape()[0]; ++row) {
        if (row == a + exiting) {
            continue;
        }

        for (std::int64_t col = 0; col < table.shape()[1]; ++col) {
            if (col == entering) {
                continue;
            }

            table[row][col] -= table[row][entering] * table[a + exiting][col];
        }

        table[row][entering] /= -table[a + exiting][entering];
    }

    table[a + exiting][entering] = 1 / table[a + exiting][entering];

    std::swap(N[entering], B[exiting]);
}

static bool step(lattice::rmat_view_t table, std::vector<std::int64_t>& B, std::vector<std::int64_t>& N, std::int64_t variables, std::int64_t constraints) {
    std::int64_t a = table.shape()[0] - constraints;
    std::int64_t b = table.shape()[1] - 1;

    bool bland = false;

    std::int64_t entering = -1;
    std::int64_t exiting = -1;

    lattice::rat_t t0, t1;

    for (std::int64_t row = 0; row < constraints; ++row) {
        if (table[a + row][b].is_zero()) {
            bland = true;
            break;
        }
    }

    for (std::int64_t col = 0; col < b; ++col) {
        if (table[0][col].sign() > 0 && (entering == -1 || table[0][col] > t0)) {
            entering = col;
            t0 = table[0][col];

            if (bland) {
                break;
            }
        }
    }

    if (entering == -1) {
        return true;
    }

    for (std::int64_t row = 0; row < constraints; ++row) {
        if (table[a + row][entering].sign() > 0) {
            t1 = table[a + row][b] / table[a + row][entering];

            if (exiting == -1 || t1 < t0) {
                exiting = row;
                t0 = t1;
            }
        }
    }

    assert(exiting != -1);

    ::pivot(table, B, N, variables, constraints, entering, exiting);

    return false;
}

void lattice::optimize(lattice::rmat_view_t x, const lattice::rmat_view_t& initial_table, std::int64_t size, std::int64_t depth) {
    assert(x.shape()[0] == size);
    assert(x.shape()[1] == 1);

    lattice::rmat_t table(boost::extents[2 + size + depth][size + 1]);

    for (std::int64_t row = 1; row < 2 + size + depth; ++row) {
        for (std::int64_t col = 0; col < size + 1; ++col) {
            table[row][col] = initial_table[row - 1][col];
        }
    }

    for (std::int64_t row = 2 + size; row < 2 + size + depth; ++row) {
        for (std::int64_t col = 0; col < size + 1; ++col) {
            table[0][col] += table[row][col];
        }
    }

    std::vector<std::int64_t> B(size + depth, 0);
    std::vector<std::int64_t> N(size, 0);

    for (std::int64_t i = 0; i < size; ++i) {
        N[i] = i;
    }

    for (std::int64_t i = 0; i < size + depth; ++i) {
        B[i] = size + i;
    }

    lattice::rmat_view_t view0 = lattice::view(table, 0, 2 + size + depth, 0, size + 1);

    while (!::step(view0, B, N, 2 * size + depth, size + depth)) {
        //
    }

    for (std::int64_t row = 0; row < size + depth; ++row) {
        if (B[row] >= 2 * size) {
            assert(table[2 + row][size].is_zero());

            for (std::int64_t col = 0; col < size; ++col) {
                if (N[col] < 2 * size && !table[2 + row][col].is_zero()) {
                    ::pivot(view0, B, N, 2 * size + depth, size + depth, col, row);
                    break;
                }
            }
        }
    }

    for (std::int64_t c0 = 0, c1 = size - 1; c0 < size - depth; ++c0) {
        if (N[c0] >= 2 * size) {
            for (;; --c1) {
                if (N[c1] < 2 * size) {
                    for (std::int64_t row = 0; row < table.shape()[0]; ++row) {
                        boost::multiprecision::swap(table[row][c0], table[row][c1]);
                    }

                    std::swap(N[c0], N[c1]);
                    break;
                }
            }
        }
    }

    for (std::int64_t row = 0; row < 2  + size + depth; ++row) {
        boost::multiprecision::swap(table[row][size - depth], table[row][size]);
    }

    lattice::rmat_view_t view1 = lattice::view(table, 1, 2 + size + depth, 0, size - depth + 1);

    while (!::step(view1, B, N, 2 * size, size + depth)) {
        //
    }

    for (std::int64_t i = 0; i < size; ++i) {
        x[i][0] = 0;
    }

    for (std::int64_t i = 0; i < size + depth; ++i) {
        if (B[i] < size) {
            x[B[i]][0] = table[2 + i][size - depth];
        }
    }
}
