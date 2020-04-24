#include <cassert>
#include <condition_variable>
#include <thread>
#include <mutex>
#include <vector>

#include "fmt/format.h"
#include "lattice.h"

struct search_info_t {
    std::int64_t size;
    std::int64_t depth;

    const lattice::rmat_view_t& transform;
    const lattice::rmat_view_t& offset;
    lattice::rmat_t fixed;

    lattice::rmat_t table;
    lattice::rmat_t x;

    std::vector<lattice::rmat_t>& results;

    std::mutex& mutex;
    std::condition_variable& finished;
    std::int64_t& available;
};

static void search_thread(search_info_t info);

static void search(search_info_t& info) {
    // if (info.depth > 0) {
    //     return;
    // }

    if (info.depth == info.size) {
        std::unique_lock lock(info.mutex);

        info.results.push_back(info.fixed);

        fmt::print("found: ");

        for (std::int64_t i = 0; i < info.size; ++i) {
            fmt::print("{} ", info.fixed[i][0].str());
        }

        fmt::print("\n");
    } else {
        lattice::rat_t t0;
        lattice::int_t min, max;

        // lower bound

        for (std::int64_t col = 0; col < info.size; ++col) {
            info.table[0][col] = -info.transform[info.depth][col];
        }

        lattice::optimize(lattice::view(info.x), lattice::view(info.table), info.size, info.depth);

        t0 = info.offset[info.depth][0];

        for (std::int64_t col = 0; col < info.size; ++col) {
            t0 += info.transform[info.depth][col] * info.x[col][0];
        }

        min = lattice::ceil(t0);

        // upper bound

        for (std::int64_t col = 0; col < info.size; ++col) {
            info.table[0][col] = info.transform[info.depth][col];
        }

        lattice::optimize(lattice::view(info.x), lattice::view(info.table), info.size, info.depth);

        t0 = info.offset[info.depth][0];

        for (std::int64_t col = 0; col < info.size; ++col) {
            t0 += info.transform[info.depth][col] * info.x[col][0];
        }

        max = lattice::floor(t0);

        //fmt::print("{}: {} -> {}\n", info.depth, min.str(), max.str());

        lattice::rat_t offset = info.offset[info.depth][0];

        for (lattice::int_t x = min; x <= max; ++x) {
            if (x <= offset) {
                for (std::int64_t col = 0; col < info.size; ++col) {
                    info.table[1 + info.size + info.depth][col] = -info.transform[info.depth][col];
                }

                info.table[1 + info.size + info.depth][info.size] = offset - x;
            } else {
                for (std::int64_t col = 0; col < info.size; ++col) {
                    info.table[1 + info.size + info.depth][col] = info.transform[info.depth][col];
                }

                info.table[1 + info.size + info.depth][info.size] = x - offset;
            }

            info.fixed[info.depth][0] = x;
            info.depth += 1;

            std::unique_lock lock(info.mutex);

            if (info.available > 0) {
                --info.available;
                lock.unlock();
                std::thread(search_thread, info).detach();
            } else {
                lock.unlock();
                search(info);
            }

            info.depth -= 1;
        }
    }
}

static void search_thread(search_info_t info) {
    ::search(info);

    std::unique_lock lock(info.mutex);

    ++info.available;
    info.finished.notify_all();
}

static void search_parallel(const search_info_t& root) {
    std::unique_lock lock(root.mutex);

    std::int64_t threads = root.available--;
    std::thread(search_thread, root).detach();

    while (root.available < threads) {
        root.finished.wait(lock);
    }
}

void lattice::search(const lattice::rmat_t& basis, const lattice::rmat_t& lower, const lattice::rmat_t& upper, std::vector<lattice::rmat_t>& results, std::int64_t thread_max) {
    using boost::extents;
    using boost::indices;
    using range = lattice::rmat_t::index_range;

    std::int64_t size = basis.shape()[0];

    assert(basis.shape()[0] == size);
    assert(basis.shape()[1] == size);
    assert(lower.shape()[0] == size);
    assert(lower.shape()[1] == 1);
    assert(upper.shape()[0] == size);
    assert(upper.shape()[1] == 1);

    lattice::rmat_t transform(extents[size][size]);
    lattice::rmat_t offset(extents[size][1]);

    for (std::int64_t i = 0; i < size; ++i) {
        transform[i][i] = 1;
        offset[i][0] = lower[i][0];
    }

    lattice::rmat_t lu(extents[size][size]);
    std::vector<std::int64_t> p;

    lu = basis;

    lattice::lu(lu[indices[range()][range()]], p);
    lattice::solve(transform[indices[range()][range()]], lu[indices[range()][range()]], p);
    lattice::solve(offset[indices[range()][range()]], lu[indices[range()][range()]], p);

    std::mutex mutex;
    std::condition_variable finished;

    ::search_info_t root {
        .size       = size,
        .depth      = 0,
        .transform  = transform[indices[range()][range()]],
        .offset     = offset[indices[range()][range()]],
        .fixed      = lattice::rmat_t(extents[size][1]),
        .table      = lattice::rmat_t(extents[2 * size + 1][size + 1]),
        .x          = lattice::rmat_t(extents[size][1]),
        .results    = results,
        .mutex      = mutex,
        .finished   = finished,
        .available  = thread_max
    };

    for (std::int64_t i = 0; i < size; ++i) {
        root.table[i + 1][i] = 1;
        root.table[i + 1][size] = upper[i][0] - lower[i][0];
    }

    ::search_parallel(root);
}
