#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multi_array.hpp>
#include <chrono>
#include <thread>

#include "lattice.h"
#include "fmt/format.h"

int main(int argc, char** argv) {
    std::int64_t size;
    lattice::rmat_t basis;
    lattice::rmat_t lower;
    lattice::rmat_t upper;

    std::vector<lattice::rmat_t> results;

    using clock = std::chrono::high_resolution_clock;

    clock::time_point start = clock::now();

    lattice::parse(argc < 2 ? "" : argv[1], size, basis, lower, upper);
    lattice::search(basis, lower, upper, results, std::thread::hardware_concurrency());
    // lattice::search(basis, lower, upper, results, 1);

    clock::time_point end = clock::now();

    std::chrono::duration<float> elapsed = end - start;

    for (const lattice::rmat_t& result : results) {
        for (std::int64_t i = 0; i < size; ++i) {
            fmt::print("{} ", result[i][0].str());
        }

        fmt::print("\n");
    }

    fmt::print("elapsed: {:7.3f}\n", elapsed.count());
    fmt::print("count:   {:7d}\n", results.size());
}
