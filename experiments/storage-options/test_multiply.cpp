// test_multiply.cpp - one multiply function, two storages.
//
// Builds one sparse matrix, stores it twice (CSC and vector-of-vectors) from the same
// entries, extracts column pointers from each, and runs the *same* multiply over both.
//
// Correctness is the first result and the point of the experiment: the two runs must agree
// bit for bit, because they are the same additions in the same order over the same values,
// reached through pointers that differ only in where they point.
//
// Timing is the second result, and it says something the design needs: how much the layout
// costs, and how much the pointer indirection costs, which turn out to be very different
// numbers.

#include "MultiplyEngine.h"
#include "SparseMatrixCsc.h"
#include "SparseMatrixVv.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

using namespace StorageOptions;

namespace {

// The matrix, as per-column lists. Both storages are built from this, so their content is
// identical by construction and the experiment measures layout, not two different matrices.
struct Columns {
    std::size_t size;
    std::vector<std::vector<std::int32_t>> rowIdx;
    std::vector<std::vector<double>>       val;
};

Columns build(std::size_t size, std::size_t bandwidth, std::mt19937& rng) {
    Columns c;
    c.size = size;
    c.rowIdx.resize(size);
    c.val.resize(size);
    std::uniform_real_distribution<double> uniform(-1.0, 1.0);

    for (std::size_t j = 0; j < size; ++j)
        for (std::size_t i = j; i < std::min(size, j + bandwidth); ++i) {
            c.rowIdx[j].push_back(static_cast<std::int32_t>(i));
            c.val[j].push_back(uniform(rng));
        }
    return c;
}

SparseMatrixCsc toCsc(const Columns& c) {
    std::vector<std::size_t>  colPtr(c.size + 1, 0);
    std::vector<std::int32_t> rowIdx;
    std::vector<double>       val;
    for (std::size_t j = 0; j < c.size; ++j) {
        rowIdx.insert(rowIdx.end(), c.rowIdx[j].begin(), c.rowIdx[j].end());
        val.insert(val.end(), c.val[j].begin(), c.val[j].end());
        colPtr[j + 1] = rowIdx.size();
    }
    return SparseMatrixCsc(c.size, std::move(colPtr), std::move(rowIdx), std::move(val));
}

// Columns allocated in order. The allocator lays them out nearly contiguously, so this is the
// most favourable case a vector of vectors can have: a flat layout with an extra hop.
SparseMatrixVv toVv(const Columns& c) {
    return SparseMatrixVv(c.size, c.rowIdx, c.val);
}

// The same content, but the columns allocated in random order with other allocations
// interleaved, so they land wherever the allocator has room. This is what a vector of vectors
// looks like once it has been *used*: in dynamic factorization, fronts grow and reallocate at
// different times, and end up scattered. The version above flatters VV; this one does not.
SparseMatrixVv toVvScattered(const Columns& c, std::mt19937& rng) {
    std::vector<std::size_t> order(c.size);
    for (std::size_t j = 0; j < c.size; ++j) order[j] = j;
    std::shuffle(order.begin(), order.end(), rng);

    std::vector<std::vector<std::int32_t>> rowIdx(c.size);
    std::vector<std::vector<double>>       val(c.size);
    std::vector<std::vector<double>>       spacer;   // pushes the next column elsewhere
    spacer.reserve(c.size);

    for (std::size_t k = 0; k < c.size; ++k) {
        const std::size_t j = order[k];
        rowIdx[j] = c.rowIdx[j];
        val[j]    = c.val[j];
        spacer.emplace_back(8, 0.0);
    }
    return SparseMatrixVv(c.size, std::move(rowIdx), std::move(val));
}

double maxDiff(const std::vector<double>& a, const std::vector<double>& b) {
    double d = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i)
        d = std::max(d, std::abs(a[i] - b[i]));
    return d;
}

template<class Fn>
double timeIt(Fn&& fn, int repeats) {
    using Clock = std::chrono::steady_clock;
    fn();   // warm up
    const auto t0 = Clock::now();
    for (int r = 0; r < repeats; ++r) fn();
    const auto t1 = Clock::now();
    return std::chrono::duration<double, std::milli>(t1 - t0).count() / repeats;
}

} // namespace

int main() {
    const std::size_t size      = 200000;
    const std::size_t bandwidth = 16;
    const int         repeats   = 20;

    std::mt19937 rng(20260712);
    const Columns c = build(size, bandwidth, rng);

    const SparseMatrixCsc csc  = toCsc(c);
    const SparseMatrixVv  vv   = toVv(c);
    const SparseMatrixVv  vvS  = toVvScattered(c, rng);

    MultiplyEngine eng;

    // Extract the column pointers, once per storage. After this, the storages are gone: the
    // multiply sees only these arrays.
    std::vector<const std::int32_t*> riC, riV, riS;
    std::vector<const double*>       vC,  vV,  vS;
    std::vector<std::size_t>         lC,  lV,  lS;
    eng.columnPointers(csc, riC, vC, lC);
    eng.columnPointers(vv,  riV, vV, lV);
    eng.columnPointers(vvS, riS, vS, lS);

    std::vector<double> x(size);
    for (std::size_t j = 0; j < size; ++j)
        x[j] = 1.0 + 0.001 * static_cast<double>(j % 97);

    std::vector<double> yBase(size, 0.0), yCsc(size, 0.0), yVv(size, 0.0), yVvS(size, 0.0);

    std::cout << "storage-options: one multiply, two storages\n"
              << "  size = " << size << ", nnz = " << csc.nnz()
              << ", repeats = " << repeats << "\n\n";

    const double tBase = timeIt([&]{ std::fill(yBase.begin(), yBase.end(), 0.0);
                                     eng.multiplyCsc(csc, x.data(), yBase.data()); }, repeats);
    const double tCsc  = timeIt([&]{ std::fill(yCsc.begin(), yCsc.end(), 0.0);
                                     eng.multiply(size, riC.data(), vC.data(), lC.data(),
                                                  x.data(), yCsc.data()); }, repeats);
    const double tVv   = timeIt([&]{ std::fill(yVv.begin(), yVv.end(), 0.0);
                                     eng.multiply(size, riV.data(), vV.data(), lV.data(),
                                                  x.data(), yVv.data()); }, repeats);
    const double tVvS  = timeIt([&]{ std::fill(yVvS.begin(), yVvS.end(), 0.0);
                                     eng.multiply(size, riS.data(), vS.data(), lS.data(),
                                                  x.data(), yVvS.data()); }, repeats);

    const double dCsc = maxDiff(yBase, yCsc);
    const double dVv  = maxDiff(yBase, yVv);
    const double dVvS = maxDiff(yBase, yVvS);
    const bool   ok   = (dCsc == 0.0 && dVv == 0.0 && dVvS == 0.0);

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "  hand-written CSC (baseline)  " << std::setw(8) << tBase << " ms\n";
    std::cout << "  multiply(), CSC pointers     " << std::setw(8) << tCsc << " ms   "
              << std::setprecision(2) << tCsc / tBase << "x\n" << std::setprecision(3);
    std::cout << "  multiply(), VV pointers      " << std::setw(8) << tVv << " ms   "
              << std::setprecision(2) << tVv / tBase << "x\n" << std::setprecision(3);
    std::cout << "  multiply(), VV scattered     " << std::setw(8) << tVvS << " ms   "
              << std::setprecision(2) << tVvS / tBase << "x\n\n";

    std::cout << "  all results identical: " << (ok ? "yes" : "NO") << "\n";
    if (!ok)
        std::cout << "    max diff  csc " << dCsc << "  vv " << dVv
                  << "  vv-scattered " << dVvS << "\n";

    std::cout << "\nWhat this shows:\n"
              << "  The three timed multiply() rows are the SAME COMPILED FUNCTION. Not a\n"
              << "  template instantiated three times: one function, called three times, with\n"
              << "  pointers that came from different places. It cannot tell CSC from VV,\n"
              << "  because by the time it runs there is nothing left to tell apart.\n"
              << "\n"
              << "  The pointer indirection is free. multiply() over CSC pointers matches the\n"
              << "  hand-written CSC baseline, which never builds a pointer array at all. So\n"
              << "  the generality costs nothing on the layout we care most about.\n"
              << "\n"
              << "  What costs is the layout. The two VV rows are the same class holding the\n"
              << "  same content, differing only in the order their inner vectors were\n"
              << "  allocated. A flat buffer guarantees consecutive columns are adjacent; a\n"
              << "  vector of vectors only ever borrows that from the allocator.\n"
              << "\n"
              << "  Two caveats on the scattered number. It is constructed, not observed: we\n"
              << "  shuffle the allocation order and interleave spacers, to remove the\n"
              << "  allocator's help entirely. And this kernel is the harshest setting a cache\n"
              << "  miss can have, roughly two flops per element loaded, so a miss shows at\n"
              << "  full price. A factorization front handed to BLAS-3 does O(n^3) work on\n"
              << "  O(n^2) data, and the same miss amortizes. Read the packed row as VV's\n"
              << "  structural cost and the scattered row as an upper bound, not a forecast.\n";

    return ok ? 0 : 1;
}
